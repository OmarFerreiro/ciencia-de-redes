//! Reproducción de Gómez-Gardeñes, Gómez, Arenas, Moreno,
//! "Explosive Synchronization Transitions in Scale-free Networks",
//! Phys. Rev. Lett. 106, 128701 (2011). arXiv:1102.4823
//!
//! Kuramoto: dθ_i/dt = ω_i + λ Σ_j A_ij sin(θ_j − θ_i), con ω_i = k_i.
//!
//! Uso:
//!   explosive_sync sweep <alpha> <n> <k_medio> <seed> <lmin> <lmax> <dl> <salida.csv>
//!   explosive_sync star  <K> <dl> <salida.csv>

use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
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
}

/// Modelo BA–ER de Gómez-Gardeñes & Moreno (2006): cada enlace nuevo se une,
/// con probabilidad `alpha`, a un nodo elegido al azar (ER) y, con probabilidad
/// `1 - alpha`, a un nodo elegido con probabilidad proporcional a su grado (BA).
/// Cada nodo nuevo añade `m` enlaces; <k> ≈ 2m.
fn ba_er(n: usize, m: usize, alpha: f64, rng: &mut ChaCha8Rng) -> Graph {
    let mut edges: Vec<(usize, usize)> = Vec::new();
    // `stubs` contiene cada nodo tantas veces como su grado: muestrear de aquí es preferencial.
    let mut stubs: Vec<usize> = Vec::new();
    // Núcleo inicial: grafo completo de m+1 nodos.
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

struct Kuramoto<'a> {
    g: &'a Graph,
    omega: Vec<f64>,
    lambda: f64,
}

impl Kuramoto<'_> {
    fn deriv(&self, th: &[f64], out: &mut [f64]) {
        for i in 0..self.g.n() {
            let (a, b) = (self.g.offsets[i], self.g.offsets[i + 1]);
            let mut s = 0.0;
            for &j in &self.g.neighbors[a..b] {
                s += (th[j as usize] - th[i]).sin();
            }
            out[i] = self.omega[i] + self.lambda * s;
        }
    }

    /// Un paso de Runge-Kutta de 4º orden con buffers reutilizados.
    fn rk4(&self, th: &mut [f64], dt: f64, buf: &mut [Vec<f64>; 5]) {
        let n = th.len();
        let [k1, k2, k3, k4, tmp] = buf;
        self.deriv(th, k1);
        for i in 0..n {
            tmp[i] = th[i] + 0.5 * dt * k1[i];
        }
        self.deriv(tmp, k2);
        for i in 0..n {
            tmp[i] = th[i] + 0.5 * dt * k2[i];
        }
        self.deriv(tmp, k3);
        for i in 0..n {
            tmp[i] = th[i] + dt * k3[i];
        }
        self.deriv(tmp, k4);
        for i in 0..n {
            th[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
    }
}

fn order_parameter(th: &[f64]) -> f64 {
    let (mut c, mut s) = (0.0, 0.0);
    for &t in th {
        c += t.cos();
        s += t.sin();
    }
    (c * c + s * s).sqrt() / th.len() as f64
}

/// Continuación adiabática: se mantiene el estado entre valores consecutivos de λ.
/// Para cada λ: transitorio `t_trans`, luego media de r durante `t_meas`.
fn continuation(
    g: &Graph,
    omega: &[f64],
    lambdas: &[f64],
    th: &mut Vec<f64>,
    t_trans: f64,
    t_meas: f64,
) -> Vec<f64> {
    let kmax = (0..g.n()).map(|i| g.degree(i)).max().unwrap() as f64;
    let n = g.n();
    let mut buf: [Vec<f64>; 5] = std::array::from_fn(|_| vec![0.0; n]);
    let mut result = Vec::with_capacity(lambdas.len());
    for &lam in lambdas {
        let sys = Kuramoto { g, omega: omega.to_vec(), lambda: lam };
        // Paso estable para RK4: dt·λ·k_max ≲ 1 (el nodo más conectado es el más rápido).
        let dt = (0.01_f64).min(0.5 / (lam.max(0.1) * kmax));
        let steps_trans = (t_trans / dt) as usize;
        let steps_meas = (t_meas / dt) as usize;
        for _ in 0..steps_trans {
            sys.rk4(th, dt, &mut buf);
        }
        let mut acc = 0.0;
        for _ in 0..steps_meas {
            sys.rk4(th, dt, &mut buf);
            acc += order_parameter(th);
        }
        result.push(acc / steps_meas as f64);
    }
    result
}

fn arange(a: f64, b: f64, d: f64) -> Vec<f64> {
    let n = ((b - a) / d).round() as usize;
    (0..=n).map(|i| a + i as f64 * d).collect()
}

fn write_csv(path: &str, lambdas: &[f64], fwd: &[f64], bwd: &[f64]) {
    let mut f = std::fs::File::create(path).expect("no se pudo crear el CSV");
    writeln!(f, "lambda,r_forward,r_backward").unwrap();
    for i in 0..lambdas.len() {
        writeln!(f, "{:.4},{:.6},{:.6}", lambdas[i], fwd[i], bwd[i]).unwrap();
    }
}

fn run(g: &Graph, lambdas: &[f64], seed: u64, t_trans: f64, t_meas: f64) -> (Vec<f64>, Vec<f64>) {
    let omega: Vec<f64> = (0..g.n()).map(|i| g.degree(i) as f64).collect();
    let mut rng = ChaCha8Rng::seed_from_u64(seed ^ 0x9E37_79B9);
    let mut th: Vec<f64> = (0..g.n()).map(|_| rng.gen_range(0.0..std::f64::consts::TAU)).collect();
    let fwd = continuation(g, &omega, lambdas, &mut th, t_trans, t_meas);
    let rev: Vec<f64> = lambdas.iter().rev().copied().collect();
    let mut bwd = continuation(g, &omega, &rev, &mut th, t_trans, t_meas);
    bwd.reverse();
    (fwd, bwd)
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let mode = a.get(1).map(String::as_str).unwrap_or("");
    match mode {
        "sweep" if a.len() == 10 => {
            let alpha: f64 = a[2].parse().unwrap();
            let n: usize = a[3].parse().unwrap();
            let kmean: f64 = a[4].parse().unwrap();
            let seed: u64 = a[5].parse().unwrap();
            let (lmin, lmax, dl): (f64, f64, f64) =
                (a[6].parse().unwrap(), a[7].parse().unwrap(), a[8].parse().unwrap());
            let mut rng = ChaCha8Rng::seed_from_u64(seed);
            let g = ba_er(n, (kmean / 2.0).round() as usize, alpha, &mut rng);
            let kmax = (0..n).map(|i| g.degree(i)).max().unwrap();
            let mean = g.neighbors.len() as f64 / n as f64;
            eprintln!("N={n} <k>={mean:.3} k_max={kmax} alpha={alpha} seed={seed}");
            let lambdas = arange(lmin, lmax, dl);
            let (fwd, bwd) = run(&g, &lambdas, seed, 40.0, 10.0);
            write_csv(&a[9], &lambdas, &fwd, &bwd);
        }
        "star" if a.len() == 5 => {
            let k: usize = a[2].parse().unwrap();
            let dl: f64 = a[3].parse().unwrap();
            let g = star(k);
            let lc = (k as f64 - 1.0) / (k as f64 + 1.0);
            let rc = k as f64 / (k as f64 + 1.0);
            eprintln!("Estrella K={k}: teoría λc={lc:.4}, r_c={rc:.4}");
            // ω_hub = K, ω_hoja = 1 (grado)
            let lambdas = arange(0.0, 2.0, dl);
            let (fwd, bwd) = run(&g, &lambdas, 1, 200.0, 50.0);
            write_csv(&a[4], &lambdas, &fwd, &bwd);
        }
        _ => {
            eprintln!("Uso:\n  explosive_sync sweep <alpha> <n> <k_medio> <seed> <lmin> <lmax> <dl> <salida.csv>\n  explosive_sync star <K> <dl> <salida.csv>");
            std::process::exit(2);
        }
    }
}
