//! Reproducción de Gómez-Gardeñes, Gómez, Arenas, Moreno,
//! "Explosive Synchronization Transitions in Scale-free Networks",
//! Phys. Rev. Lett. 106, 128701 (2011). arXiv:1102.4823
//!
//! Kuramoto: dθ_i/dt = ω_i + λ Σ_j A_ij sin(θ_j − θ_i), con ω_i = k_i.
//!
//! Uso:
//!   explosive_sync sweep <alpha> <n> <k_medio> <n_redes> <lmin> <lmax> <dl> <salida.csv>
//!   explosive_sync star  <K> <dl> <salida.csv>
//!   explosive_sync plot  <entrada.csv> <salida.svg> "<titulo>"
//!
//! Variables de entorno opcionales: TTRANS (40), TMEAS (10), DTF (1.0, factor del paso).
//! CSV: lambda, r_fwd, std_fwd, r_bwd, std_bwd, wHub_fwd, wHub_bwd
//! (media sobre redes; wHub = frecuencia efectiva del nodo de mayor grado).

use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use std::fmt::Write as _;
use std::io::Write;

/// Grafo en formato CSR (listas de vecinos contiguas en memoria).
struct Graph {
    offsets: Vec<usize>,
    neighbors: Vec<u32>,
}

impl Graph {
    fn from_edges(n: usize, edges: &[(usize, usize)]) -> Self {
        let mut deg = vec![0usize; n];
        for &(a, b) in edges {
            deg[a] += 1;
            deg[b] += 1;
        }
        let mut offsets = vec![0usize; n + 1];
        for i in 0..n {
            offsets[i + 1] = offsets[i] + deg[i];
        }
        let mut fill = offsets.clone();
        let mut neighbors = vec![0u32; offsets[n]];
        for &(a, b) in edges {
            neighbors[fill[a]] = b as u32;
            fill[a] += 1;
            neighbors[fill[b]] = a as u32;
            fill[b] += 1;
        }
        Graph { offsets, neighbors }
    }

    fn n(&self) -> usize {
        self.offsets.len() - 1
    }

    fn degree(&self, i: usize) -> usize {
        self.offsets[i + 1] - self.offsets[i]
    }

    fn hub(&self) -> usize {
        (0..self.n()).max_by_key(|&i| self.degree(i)).unwrap()
    }
}

/// Modelo BA–ER de Gómez-Gardeñes & Moreno (2006): cada enlace nuevo se une,
/// con probabilidad `alpha`, a un nodo elegido al azar (ER) y, con probabilidad
/// `1 - alpha`, a un nodo elegido con probabilidad proporcional a su grado (BA).
/// Cada nodo nuevo añade `m` enlaces; <k> ≈ 2m.
fn ba_er(n: usize, m: usize, alpha: f64, rng: &mut ChaCha8Rng) -> Graph {
    let mut edges: Vec<(usize, usize)> = Vec::new();
    if alpha < 0.0 {
        // Erdős–Rényi G(N, M) con M = N·m enlaces (⟨k⟩ = 2m), degrees de Poisson.
        let mut seen = std::collections::HashSet::new();
        while edges.len() < n * m {
            let (a, b) = (rng.gen_range(0..n), rng.gen_range(0..n));
            if a != b && seen.insert((a.min(b), a.max(b))) {
                edges.push((a, b));
            }
        }
        return Graph::from_edges(n, &edges);
    }
    // `stubs` contiene cada nodo tantas veces como su grado: muestrear de aquí es preferencial.
    let mut stubs: Vec<usize> = Vec::new();
    for i in 0..=m {
        for j in 0..i {
            edges.push((i, j));
            stubs.push(i);
            stubs.push(j);
        }
    }
    for new in (m + 1)..n {
        let mut targets: Vec<usize> = Vec::with_capacity(m);
        while targets.len() < m {
            let t = if rng.gen::<f64>() < alpha {
                rng.gen_range(0..new)
            } else {
                stubs[rng.gen_range(0..stubs.len())]
            };
            if !targets.contains(&t) {
                targets.push(t);
            }
        }
        for t in targets {
            edges.push((new, t));
            stubs.push(new);
            stubs.push(t);
        }
    }
    Graph::from_edges(n, &edges)
}

fn star(k: usize) -> Graph {
    let edges: Vec<(usize, usize)> = (1..=k).map(|j| (0, j)).collect();
    Graph::from_edges(k + 1, &edges)
}

/// Buffers de trabajo reutilizados en cada paso de integración.
struct Work {
    k: [Vec<f64>; 4],
    tmp: Vec<f64>,
    sin: Vec<f64>,
    cos: Vec<f64>,
}

impl Work {
    fn new(n: usize) -> Self {
        Work {
            k: std::array::from_fn(|_| vec![0.0; n]),
            tmp: vec![0.0; n],
            sin: vec![0.0; n],
            cos: vec![0.0; n],
        }
    }
}

/// dθ_i/dt. Se usa sin(θ_j − θ_i) = sinθ_j cosθ_i − cosθ_j sinθ_i, de modo que
/// solo se calculan N senos/cosenos por evaluación en vez de 2E senos.
fn deriv(g: &Graph, omega: &[f64], lambda: f64, th: &[f64], s: &mut [f64], c: &mut [f64], out: &mut [f64]) {
    for i in 0..th.len() {
        let (si, ci) = th[i].sin_cos();
        s[i] = si;
        c[i] = ci;
    }
    for i in 0..th.len() {
        let (a, b) = (g.offsets[i], g.offsets[i + 1]);
        let (mut ss, mut cc) = (0.0, 0.0);
        for &j in &g.neighbors[a..b] {
            ss += s[j as usize];
            cc += c[j as usize];
        }
        out[i] = omega[i] + lambda * (ss * c[i] - cc * s[i]);
    }
}

/// Un paso de Runge-Kutta de 4º orden.
fn rk4(g: &Graph, omega: &[f64], lambda: f64, th: &mut [f64], dt: f64, w: &mut Work) {
    let n = th.len();
    let Work { k, tmp, sin, cos } = w;
    let [k1, k2, k3, k4] = k;
    deriv(g, omega, lambda, th, sin, cos, k1);
    for i in 0..n {
        tmp[i] = th[i] + 0.5 * dt * k1[i];
    }
    deriv(g, omega, lambda, tmp, sin, cos, k2);
    for i in 0..n {
        tmp[i] = th[i] + 0.5 * dt * k2[i];
    }
    deriv(g, omega, lambda, tmp, sin, cos, k3);
    for i in 0..n {
        tmp[i] = th[i] + dt * k3[i];
    }
    deriv(g, omega, lambda, tmp, sin, cos, k4);
    for i in 0..n {
        th[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
}

fn order_parameter(th: &[f64]) -> f64 {
    let (mut c, mut s) = (0.0, 0.0);
    for &t in th {
        let (si, ci) = t.sin_cos();
        c += ci;
        s += si;
    }
    (c * c + s * s).sqrt() / th.len() as f64
}

#[derive(Clone, Copy, Default)]
struct Point {
    r: f64,
    w_hub: f64,
}

#[derive(Clone, Copy)]
struct Params {
    t_trans: f64,
    t_meas: f64,
    dtf: f64,
}

impl Params {
    fn from_env() -> Self {
        let get = |k: &str, d: f64| std::env::var(k).ok().and_then(|v| v.parse().ok()).unwrap_or(d);
        Params { t_trans: get("TTRANS", 40.0), t_meas: get("TMEAS", 10.0), dtf: get("DTF", 1.0) }
    }
}

/// Continuación adiabática: se mantiene el estado entre valores consecutivos de λ.
/// Para cada λ: transitorio `t_trans`, luego media de r y frecuencia efectiva del hub
/// (Δθ_hub/T, con la fase sin envolver) durante `t_meas`.
fn continuation(g: &Graph, omega: &[f64], lambdas: &[f64], th: &mut [f64], p: Params) -> Vec<Point> {
    let kmax = g.degree(g.hub()) as f64;
    let hub = g.hub();
    let mut w = Work::new(g.n());
    let mut out = Vec::with_capacity(lambdas.len());
    for &lam in lambdas {
        // Estabilidad de RK4: dt·λ·k_max del orden de 1 (el hub es el nodo más rápido).
        let dt = (0.01_f64).min(p.dtf / (lam.max(0.1) * kmax));
        let steps_trans = (p.t_trans / dt) as usize;
        let steps_meas = (p.t_meas / dt) as usize;
        for _ in 0..steps_trans {
            rk4(g, omega, lam, th, dt, &mut w);
        }
        let th0 = th[hub];
        let mut acc = 0.0;
        for _ in 0..steps_meas {
            rk4(g, omega, lam, th, dt, &mut w);
            acc += order_parameter(th);
        }
        out.push(Point {
            r: acc / steps_meas as f64,
            w_hub: (th[hub] - th0) / (steps_meas as f64 * dt),
        });
    }
    out
}

fn arange(a: f64, b: f64, d: f64) -> Vec<f64> {
    let n = ((b - a) / d).round() as usize;
    (0..=n).map(|i| a + i as f64 * d).collect()
}

/// Barrido hacia adelante y hacia atrás sobre una red, con ω_i = k_i.
fn run(g: &Graph, lambdas: &[f64], seed: u64, p: Params) -> (Vec<Point>, Vec<Point>) {
    let omega: Vec<f64> = (0..g.n()).map(|i| g.degree(i) as f64).collect();
    let mut rng = ChaCha8Rng::seed_from_u64(seed ^ 0x9E37_79B9);
    let mut th: Vec<f64> = (0..g.n()).map(|_| rng.gen_range(0.0..std::f64::consts::TAU)).collect();
    let fwd = continuation(g, &omega, lambdas, &mut th, p);
    let rev: Vec<f64> = lambdas.iter().rev().copied().collect();
    let mut bwd = continuation(g, &omega, &rev, &mut th, p);
    bwd.reverse();
    (fwd, bwd)
}

fn mean_std(xs: &[f64]) -> (f64, f64) {
    let m = xs.iter().sum::<f64>() / xs.len() as f64;
    let v = xs.iter().map(|x| (x - m).powi(2)).sum::<f64>() / xs.len() as f64;
    (m, v.sqrt())
}

/// Promedia varias realizaciones y escribe el CSV.
fn write_csv(path: &str, lambdas: &[f64], runs: &[(Vec<Point>, Vec<Point>)]) {
    let mut f = std::fs::File::create(path).expect("no se pudo crear el CSV");
    writeln!(f, "lambda,r_fwd,std_fwd,r_bwd,std_bwd,whub_fwd,whub_bwd").unwrap();
    for i in 0..lambdas.len() {
        let col = |sel: &dyn Fn(&(Vec<Point>, Vec<Point>)) -> f64| -> Vec<f64> { runs.iter().map(sel).collect() };
        let (rf, sf) = mean_std(&col(&|r| r.0[i].r));
        let (rb, sb) = mean_std(&col(&|r| r.1[i].r));
        let (wf, _) = mean_std(&col(&|r| r.0[i].w_hub));
        let (wb, _) = mean_std(&col(&|r| r.1[i].w_hub));
        writeln!(f, "{:.4},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6}", lambdas[i], rf, sf, rb, sb, wf, wb).unwrap();
    }
}

/// Gráfica SVG de r(λ) con ramas hacia adelante y hacia atrás (± desviación entre redes).
fn plot_svg(csv: &str, out: &str, title: &str) {
    let text = std::fs::read_to_string(csv).expect("no se pudo leer el CSV");
    let rows: Vec<Vec<f64>> = text
        .lines()
        .skip(1)
        .map(|l| l.split(',').map(|x| x.parse().unwrap()).collect())
        .collect();
    let (w, h, ml, mr, mt, mb) = (680.0, 440.0, 60.0, 20.0, 40.0, 50.0);
    let (x0, x1) = (rows[0][0], rows[rows.len() - 1][0]);
    let px = |x: f64| ml + (x - x0) / (x1 - x0) * (w - ml - mr);
    let py = |y: f64| h - mb - y * (h - mt - mb);
    let mut s = String::new();
    write!(s, r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" font-family="sans-serif" font-size="12"><rect width="100%" height="100%" fill="white"/>"##).unwrap();
    write!(s, r##"<text x="{}" y="22" font-size="14" text-anchor="middle">{}</text>"##, w / 2.0, title).unwrap();
    for i in 0..=5 {
        let y = i as f64 / 5.0;
        let yy = py(y);
        write!(s, r##"<line x1="{ml}" x2="{}" y1="{yy:.1}" y2="{yy:.1}" stroke="#ddd"/><text x="{}" y="{:.1}" text-anchor="end">{:.1}</text>"##, w - mr, ml - 6.0, yy + 4.0, y).unwrap();
    }
    for i in 0..=6 {
        let x = x0 + (x1 - x0) * i as f64 / 6.0;
        write!(s, r##"<text x="{:.1}" y="{}" text-anchor="middle">{:.2}</text>"##, px(x), h - mb + 18.0, x).unwrap();
    }
    write!(s, r##"<text x="{}" y="{}" text-anchor="middle">λ</text><text x="16" y="{}" text-anchor="middle" transform="rotate(-90 16 {})">r</text>"##, w / 2.0, h - 10.0, h / 2.0, h / 2.0).unwrap();
    for (ci, si, color, label, lx) in [(1usize, 2usize, "#1f77b4", "adelante", 90.0), (3, 4, "#d62728", "atrás", 190.0)] {
        let upper: Vec<String> = rows.iter().map(|r| format!("{:.1},{:.1}", px(r[0]), py((r[ci] + r[si]).min(1.0)))).collect();
        let lower: Vec<String> = rows.iter().rev().map(|r| format!("{:.1},{:.1}", px(r[0]), py((r[ci] - r[si]).max(0.0)))).collect();
        write!(s, r##"<polygon points="{} {}" fill="{color}" opacity="0.15"/>"##, upper.join(" "), lower.join(" ")).unwrap();
        let line: Vec<String> = rows.iter().map(|r| format!("{:.1},{:.1}", px(r[0]), py(r[ci]))).collect();
        write!(s, r##"<polyline points="{}" fill="none" stroke="{color}" stroke-width="1.6"/>"##, line.join(" ")).unwrap();
        write!(s, r##"<circle cx="{lx}" cy="{}" r="4" fill="{color}"/><text x="{}" y="{}">{label}</text>"##, mt + 8.0, lx + 8.0, mt + 12.0).unwrap();
    }
    write!(s, r##"<rect x="{ml}" y="{mt}" width="{}" height="{}" fill="none" stroke="black"/></svg>"##, w - ml - mr, h - mt - mb).unwrap();
    std::fs::write(out, s).expect("no se pudo escribir el SVG");
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let p = Params::from_env();
    match a.get(1).map(String::as_str).unwrap_or("") {
        "sweep" if a.len() == 10 => {
            let alpha: f64 = a[2].parse().unwrap();
            let n: usize = a[3].parse().unwrap();
            let kmean: f64 = a[4].parse().unwrap();
            let nets: u64 = a[5].parse().unwrap();
            let (lmin, lmax, dl): (f64, f64, f64) =
                (a[6].parse().unwrap(), a[7].parse().unwrap(), a[8].parse().unwrap());
            let lambdas = arange(lmin, lmax, dl);
            eprintln!("N={n} <k>≈{kmean} alpha={alpha} redes={nets} puntos={} {:?}", lambdas.len(), (p.t_trans, p.t_meas, p.dtf));
            // Una red por hilo (hilos de std: sin dependencias con script de compilación).
            let next = std::sync::atomic::AtomicU64::new(1);
            let workers = std::thread::available_parallelism().map_or(1, |x| x.get());
            let mut runs: Vec<(u64, (Vec<Point>, Vec<Point>))> = std::thread::scope(|sc| {
                let handles: Vec<_> = (0..workers)
                    .map(|_| {
                        sc.spawn(|| {
                            let mut local = Vec::new();
                            loop {
                                let seed = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                                if seed > nets {
                                    break;
                                }
                                let mut rng = ChaCha8Rng::seed_from_u64(seed);
                                let g = ba_er(n, (kmean / 2.0).round() as usize, alpha, &mut rng);
                                local.push((seed, run(&g, &lambdas, seed, p)));
                            }
                            local
                        })
                    })
                    .collect();
                handles.into_iter().flat_map(|h| h.join().unwrap()).collect()
            });
            runs.sort_by_key(|r| r.0);
            let runs: Vec<_> = runs.into_iter().map(|r| r.1).collect();
            write_csv(&a[9], &lambdas, &runs);
        }
        "star" if a.len() == 5 => {
            let k: usize = a[2].parse().unwrap();
            let dl: f64 = a[3].parse().unwrap();
            let g = star(k);
            eprintln!(
                "Estrella K={k}: teoría λc={:.4}, r_c={:.4}, Ω={:.4}",
                (k as f64 - 1.0) / (k as f64 + 1.0),
                k as f64 / (k as f64 + 1.0),
                2.0 * k as f64 / (k as f64 + 1.0)
            );
            let lambdas = arange(0.0, 2.0, dl);
            let r = run(&g, &lambdas, 1, Params { t_trans: 200.0, t_meas: 50.0, ..p });
            write_csv(&a[4], &lambdas, &[r]);
        }
        "starscan" if a.len() == 4 => {
            // Compara la estrella con la teoría para varios K (K = 5, 10, 20, ..., 70).
            // Rama enganchada exacta: sin δ = (K−1)/((K+1)λ), r = |K + e^{iδ}|/(K+1), cos δ > 0.
            let dl: f64 = a[2].parse().unwrap();
            let lambdas = arange(0.0, 3.0, dl);
            let mut f = std::fs::File::create(&a[3]).expect("no se pudo crear el CSV");
            writeln!(f, "K,lc_teoria,lc_sim,err_lc,rc_paper,rc_modulo,r_sim_en_lc,rmse_rama,lambda_fwd").unwrap();
            for k in [5usize, 10, 20, 30, 40, 50, 60, 70] {
                let kf = k as f64;
                let g = star(k);
                let omega_eff = 2.0 * kf / (kf + 1.0);
                let pp = Params { t_trans: 200.0, t_meas: 50.0, ..p };
                let omega: Vec<f64> = (0..g.n()).map(|i| g.degree(i) as f64).collect();
                let locked = |pt: &Point| (pt.w_hub - omega_eff).abs() < 1e-3 * omega_eff;
                // Rama hacia atrás desde un estado sincronizado (la rama hacia adelante solo se
                // engancha a λ grande cuando K crece, así que no se usa para arrancarla).
                let mut th = vec![0.0; g.n()];
                let rev: Vec<f64> = lambdas.iter().rev().copied().collect();
                let mut bwd = continuation(&g, &omega, &rev, &mut th, pp);
                bwd.reverse();
                // Menor λ a partir del cual toda la rama hacia atrás está enganchada.
                let i_b = (0..lambdas.len()).rev().take_while(|&i| locked(&bwd[i])).last().expect("rama no enganchada");
                // Rama hacia adelante (desde fases aleatorias): solo para K ≤ 20, con λ ≤ 3.
                let lam_fwd = if k <= 20 {
                    let (fwd, _) = run(&g, &lambdas, 1, pp);
                    (0..lambdas.len()).find(|&i| locked(&fwd[i])).map_or(f64::NAN, |i| lambdas[i])
                } else {
                    f64::NAN
                };
                let exact_r = |lam: f64| {
                    let s = (kf - 1.0) / ((kf + 1.0) * lam);
                    (kf * kf + 2.0 * kf * (1.0 - s * s).max(0.0).sqrt() + 1.0).sqrt() / (kf + 1.0)
                };
                let lc = (kf - 1.0) / (kf + 1.0);
                let sq: Vec<f64> = (i_b..lambdas.len()).map(|i| (bwd[i].r - exact_r(lambdas[i])).powi(2)).collect();
                let rmse = (sq.iter().sum::<f64>() / sq.len() as f64).sqrt();
                writeln!(
                    f,
                    "{k},{lc:.4},{:.4},{:.4},{:.4},{:.4},{:.4},{:.5},{:.4}",
                    lambdas[i_b],
                    lambdas[i_b] - lc,
                    kf / (kf + 1.0),
                    (kf * kf + 1.0).sqrt() / (kf + 1.0),
                    bwd[i_b].r,
                    rmse,
                    lam_fwd
                )
                .unwrap();
            }
        }
        "plot" if a.len() == 5 => plot_svg(&a[2], &a[3], &a[4]),
        _ => {
            eprintln!("Uso:\n  explosive_sync sweep <alpha> <n> <k_medio> <n_redes> <lmin> <lmax> <dl> <salida.csv>\n  explosive_sync star <K> <dl> <salida.csv>\n  explosive_sync plot <entrada.csv> <salida.svg> \"<titulo>\"");
            std::process::exit(2);
        }
    }
}
