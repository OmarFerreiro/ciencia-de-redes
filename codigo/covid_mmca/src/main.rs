//! Núcleo del modelo MMCA de Arenas et al., Phys. Rev. X 10, 041055 (2020), Ecs. 2-27.
//! Reimplementación independiente con datos sintéticos; no usa datos reales de España.
//!
//! Uso:
//!   covid_mmca kappa             barrido de R(κ) y κ crítico (parche único) frente a la Ec. 51
//!   covid_mmca demo <salida.csv> escenario sintético de 4 parches y 3 edades

use std::fmt::Write as _;

const S: usize = 0;
const E: usize = 1;
const A: usize = 2;
const I: usize = 3;
const PH: usize = 4;
const PD: usize = 5;
const HR: usize = 6;
const HD: usize = 7;
const R: usize = 8;
const D: usize = 9;
const CH: usize = 10;
const NM: usize = 11;

/// Tasas por paso (1 paso = 1 día) o probabilidades, por estrato de edad.
#[derive(Clone)]
struct Epi {
    eta: Vec<f64>,
    alpha: Vec<f64>,
    mu: Vec<f64>,
    theta: Vec<f64>,
    gamma: Vec<f64>,
    zeta: Vec<f64>,
    lambda: Vec<f64>,
    omega: Vec<f64>,
    psi: Vec<f64>,
    chi: Vec<f64>,
}

#[derive(Clone)]
struct Params {
    ng: usize,
    np: usize,
    pop: Vec<f64>,       // n^g_i, índice g*np + i
    area: Vec<f64>,      // s_i
    mob: Vec<Vec<f64>>,  // R^g_ij por estrato, índice i*np + j
    p: Vec<f64>,         // grado de movilidad p_g
    k_tot: Vec<f64>,     // <k_g>
    k_h: Vec<f64>,       // contactos en casa
    k_hw: Vec<f64>,      // contactos en casa + trabajo
    elderly: Option<usize>,
    c: Vec<f64>,         // C_gh, índice g*ng + h
    xi: f64,
    sigma: f64,
    beta_a: f64,
    beta_i: f64,
    epi: Epi,
    phi: f64,
    delta: f64,
    kappa: Vec<f64>,     // κ0(t) por día desde t = 0 (se aplica a partir de tc)
    tc: usize,
}

struct Ctx {
    kap: f64,
    conf: bool,
    mm: Vec<Vec<f64>>, // mm[g][j*np + i] = M^g_ji
    nt: Vec<f64>,
    f: Vec<f64>,
    kg: Vec<f64>,
    z: Vec<f64>,
}

fn idx(ng: usize, np: usize, m: usize, g: usize, i: usize) -> usize {
    (m * ng + g) * np + i
}

impl Params {
    /// Confinamiento κ0(t): 0 antes de tc y, después, la serie diaria (el último valor se mantiene).
    fn kap_at(&self, t: usize) -> f64 {
        if t < self.tc {
            0.0
        } else {
            self.kappa[t.min(self.kappa.len() - 1)]
        }
    }

    /// Cantidades que dependen de t pero no de las probabilidades ρ: matrices M, ñ, f, <k_g>(t), z_g(t).
    fn ctx(&self, t: usize) -> Ctx {
        let (ng, np) = (self.ng, self.np);
        let conf = t >= self.tc;
        let kap = self.kap_at(t);

        // Matrices efectivas M^g_ji = (1-p) δ_ij + p R_ji y población efectiva
        let mut mm = vec![vec![0.0; np * np]; ng]; // mm[g][j*np + i]
        let mut nt = vec![0.0; ng * np]; // ñ^g_i
        for g in 0..ng {
            let pg = (1.0 - kap) * self.p[g];
            for j in 0..np {
                for i in 0..np {
                    let d = if i == j { 1.0 } else { 0.0 };
                    mm[g][j * np + i] = (1.0 - pg) * d + pg * self.mob[g][j * np + i];
                }
            }
            for i in 0..np {
                nt[g * np + i] = (0..np).map(|j| mm[g][j * np + i] * self.pop[g * np + j]).sum();
            }
        }
        let f: Vec<f64> = (0..np)
            .map(|i| {
                let n: f64 = (0..ng).map(|g| nt[g * np + i]).sum();
                1.0 + (1.0 - (-self.xi * n / self.area[i]).exp())
            })
            .collect();

        // contactos por estrato y normalización z_g
        let kg: Vec<f64> = (0..ng)
            .map(|g| {
                if !conf {
                    self.k_tot[g]
                } else if self.elderly == Some(g) {
                    self.k_h[g]
                } else {
                    kap * self.k_h[g] + (1.0 - kap) * (1.0 - self.delta) * self.k_hw[g]
                }
            })
            .collect();
        let z: Vec<f64> = (0..ng)
            .map(|g| {
                let ng_tot: f64 = (0..np).map(|i| self.pop[g * np + i]).sum();
                let den: f64 = (0..np).map(|i| f[i] * nt[g * np + i]).sum();
                ng_tot / den
            })
            .collect();
        Ctx { kap, conf, mm, nt, f, kg, z }
    }

    fn step(&self, rho: &mut [f64], t: usize) {
        let (ng, np) = (self.ng, self.np);
        let Ctx { kap, conf, mm, nt, f, kg, z } = self.ctx(t);

        // infecciosos presentes por parche: ih[m][h][i] = Σ_j n^h_j ρ^{m,h}_j M^h_ji / ñ^h_i
        let mut ih = vec![vec![0.0; ng * np]; 2];
        for (mi, m) in [A, I].into_iter().enumerate() {
            for h in 0..ng {
                for i in 0..np {
                    let s: f64 = (0..np)
                        .map(|j| self.pop[h * np + j] * rho[idx(ng, np, m, h, j)] * mm[h][j * np + i])
                        .sum();
                    ih[mi][h * np + i] = if nt[h * np + i] > 0.0 { s / nt[h * np + i] } else { 0.0 };
                }
            }
        }
        let (la, li) = ((1.0 - self.beta_a).ln(), (1.0 - self.beta_i).ln());

        // P^g_i y Π^g_i (Ecs. 12-13 y 23)
        let mut pp = vec![0.0; ng * np];
        for g in 0..ng {
            for i in 0..np {
                let mut lq = 0.0;
                for h in 0..ng {
                    lq += self.c[g * ng + h] * (la * ih[0][h * np + i] + li * ih[1][h * np + i]);
                }
                pp[g * np + i] = -(z[g] * kg[g] * f[i] * lq).exp_m1();
            }
        }
        let mut pi = vec![0.0; ng * np];
        for g in 0..ng {
            let pg = (1.0 - kap) * self.p[g];
            for i in 0..np {
                let rem: f64 = (0..np).map(|j| self.mob[g][i * np + j] * pp[g * np + j]).sum();
                pi[g * np + i] = (1.0 - pg) * pp[g * np + i] + pg * rem;
            }
        }

        // fracción de hogares sin infecciosos (Ec. 24)
        let chh: Vec<f64> = (0..np)
            .map(|i| {
                let ni: f64 = (0..ng).map(|g| self.pop[g * np + i]).sum();
                let s: f64 = (0..ng)
                    .map(|g| {
                        let x = 1.0
                            - rho[idx(ng, np, E, g, i)]
                            - rho[idx(ng, np, A, g, i)]
                            - rho[idx(ng, np, I, g, i)];
                        x * self.pop[g * np + i]
                    })
                    .sum();
                (s / ni).powf(self.sigma)
            })
            .collect();

        let ep = &self.epi;
        let old = rho.to_vec();
        for g in 0..ng {
            for i in 0..np {
                let o = |m: usize| old[idx(ng, np, m, g, i)];
                let x = if conf { (1.0 - self.phi) * kap * chh[i] } else { 0.0 };
                let pool = o(S) + o(CH);
                let pig = pi[g * np + i];
                let mut set = |m: usize, v: f64| rho[idx(ng, np, m, g, i)] = v;
                set(S, pool * (1.0 - pig) * (1.0 - x));
                set(CH, pool * x);
                set(E, pool * pig * (1.0 - x) + (1.0 - ep.eta[g]) * o(E));
                set(A, ep.eta[g] * o(E) + (1.0 - ep.alpha[g]) * o(A));
                set(I, ep.alpha[g] * o(A) + (1.0 - ep.mu[g]) * o(I));
                set(PD, ep.mu[g] * ep.theta[g] * o(I) + (1.0 - ep.zeta[g]) * o(PD));
                set(PH, ep.mu[g] * (1.0 - ep.theta[g]) * ep.gamma[g] * o(I) + (1.0 - ep.lambda[g]) * o(PH));
                set(
                    R,
                    ep.mu[g] * (1.0 - ep.theta[g]) * (1.0 - ep.gamma[g]) * o(I) + ep.chi[g] * o(HR) + o(R),
                );
                set(HD, ep.lambda[g] * ep.omega[g] * o(PH) + (1.0 - ep.psi[g]) * o(HD));
                set(HR, ep.lambda[g] * (1.0 - ep.omega[g]) * o(PH) + (1.0 - ep.chi[g]) * o(HR));
                set(D, ep.zeta[g] * o(PD) + ep.psi[g] * o(HD) + o(D));
            }
        }
    }
}

/// Valores de las Tablas I-IV del artículo. Las "tasas" de las tablas son duraciones en días.
fn epi_spain() -> Epi {
    let inv = |d: [f64; 3]| d.iter().map(|x| 1.0 / x).collect::<Vec<_>>();
    Epi {
        eta: inv([2.444; 3]),
        alpha: inv([5.671, 2.756, 2.756]),
        mu: inv([1.0, 3.915, 3.915]),
        theta: vec![0.0, 0.008, 0.047],
        gamma: vec![0.0003, 0.003, 0.026],
        zeta: inv([7.084; 3]),
        lambda: inv([4.084; 3]),
        omega: vec![0.3; 3],
        psi: inv([7.0; 3]),
        chi: inv([21.0; 3]),
    }
}

const C_SPAIN: [f64; 9] = [0.5980, 0.3849, 0.0171, 0.2440, 0.7210, 0.0350, 0.1919, 0.5705, 0.2376];

/// Parche único, un estrato (adultos): el caso de la Ec. 45.
fn single_patch(kappa: f64) -> Params {
    let full = epi_spain();
    let pick = |v: &Vec<f64>| vec![v[1]];
    Params {
        ng: 1,
        np: 1,
        pop: vec![1.0e6],
        area: vec![100.0],
        mob: vec![vec![1.0]],
        p: vec![0.0],
        k_tot: vec![13.3],
        k_h: vec![3.2],
        k_hw: vec![3.2 + 5.2],
        elderly: None,
        c: vec![1.0],
        xi: 0.01,
        sigma: 2.5,
        beta_a: 0.5 * 0.075,
        beta_i: 0.075,
        epi: Epi {
            eta: pick(&full.eta),
            alpha: pick(&full.alpha),
            mu: pick(&full.mu),
            theta: pick(&full.theta),
            gamma: pick(&full.gamma),
            zeta: pick(&full.zeta),
            lambda: pick(&full.lambda),
            omega: pick(&full.omega),
            psi: pick(&full.psi),
            chi: pick(&full.chi),
        },
        phi: 0.174,
        delta: 0.207,
        kappa: vec![kappa],
        tc: 0,
    }
}

/// Escenario sintético: 4 parches de tamaños distintos, 3 edades, adultos con movilidad.
fn synthetic(kappa: f64, tc: usize) -> Params {
    let np = 4;
    let ng = 3;
    let sizes = [400_000.0, 250_000.0, 150_000.0, 60_000.0];
    let frac = [0.25, 0.55, 0.20];
    let mut pop = vec![0.0; ng * np];
    for g in 0..ng {
        for i in 0..np {
            pop[g * np + i] = sizes[i] * frac[g];
        }
    }
    let ident: Vec<f64> = (0..np * np).map(|k| if k / np == k % np { 1.0 } else { 0.0 }).collect();
    let mut adult = vec![0.0; np * np];
    for i in 0..np {
        for j in 0..np {
            adult[i * np + j] = if i == j { 0.7 } else { 0.3 / (np as f64 - 1.0) };
        }
    }
    Params {
        ng,
        np,
        pop,
        area: vec![300.0, 200.0, 150.0, 80.0],
        mob: vec![ident.clone(), adult, ident],
        p: vec![1.0; 3],
        k_tot: vec![11.8, 13.3, 6.8],
        k_h: vec![3.1, 3.2, 3.3],
        k_hw: vec![4.9, 8.4, 3.3],
        elderly: Some(2),
        c: C_SPAIN.to_vec(),
        xi: 0.01,
        sigma: 2.5,
        beta_a: 0.5 * 0.075,
        beta_i: 0.075,
        epi: epi_spain(),
        phi: 0.174,
        delta: 0.207,
        kappa: vec![kappa],
        tc,
    }
}

fn init(p: &Params, seed: f64, seed_patch: usize, seed_g: usize) -> Vec<f64> {
    let (ng, np) = (p.ng, p.np);
    let mut rho = vec![0.0; NM * ng * np];
    for g in 0..ng {
        for i in 0..np {
            rho[idx(ng, np, S, g, i)] = 1.0;
        }
    }
    rho[idx(ng, np, S, seed_g, seed_patch)] -= seed;
    rho[idx(ng, np, E, seed_g, seed_patch)] += seed;
    rho
}

/// Pendiente diaria de ln(E+A+I) entre los días 60 y 120 (parche único, siembra 1e-40 para no agotar susceptibles).
fn growth_rate(kappa: f64) -> f64 {
    let p = single_patch(kappa);
    let mut rho = init(&p, 1e-40, 0, 0);
    let inf = |r: &[f64]| r[idx(1, 1, E, 0, 0)] + r[idx(1, 1, A, 0, 0)] + r[idx(1, 1, I, 0, 0)];
    let mut l200 = 0.0;
    for t in 0..120 {
        if t == 60 {
            l200 = inf(&rho).ln();
        }
        p.step(&mut rho, t);
    }
    (inf(&rho).ln() - l200) / 60.0
}

/// R estático de un parche (Ec. 45) con ρ_S = <CH> = 1 y la infectividad dada.
fn r_static(kappa: f64, beta_a: f64, beta_i: f64) -> f64 {
    let p = single_patch(kappa);
    let (k_h, k_hw) = (p.k_h[0], p.k_hw[0]);
    let dur = beta_a / p.epi.alpha[0] + beta_i / p.epi.mu[0];
    dur * (1.0 - kappa * (1.0 - p.phi)) * (kappa * k_h + (1.0 - p.delta) * (1.0 - kappa) * k_hw)
}

/// κ crítico de la Ec. 51 (raíz en [0,1]) con ρ_S = <CH> = 1.
fn kappa_c_analytic(beta_a: f64, beta_i: f64) -> f64 {
    let p = single_patch(0.0);
    let (k_h, k_hw, phi, delta) = (p.k_h[0], p.k_hw[0], p.phi, p.delta);
    let dur = beta_a / p.epi.alpha[0] + beta_i / p.epi.mu[0];
    let a = (1.0 - phi) * (k_hw * (1.0 - delta) - k_h);
    let b = k_hw * (1.0 - delta) * (1.0 + (1.0 - phi)) - k_h;
    let c = k_hw * (1.0 - delta) - 1.0 / dur;
    let disc = (b * b - 4.0 * a * c).sqrt();
    let r1 = (b - disc) / (2.0 * a);
    let r2 = (b + disc) / (2.0 * a);
    [r1, r2].into_iter().find(|r| (0.0..=1.0).contains(r)).unwrap_or(f64::NAN)
}

fn kappa_c_numeric() -> f64 {
    let (mut lo, mut hi) = (0.2, 0.95);
    for _ in 0..40 {
        let mid = 0.5 * (lo + hi);
        if growth_rate(mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

fn simulate(p: &Params, rho0: Vec<f64>, steps: usize) -> Vec<Vec<f64>> {
    let mut traj = vec![rho0];
    for t in 0..steps {
        let mut r = traj[t].clone();
        p.step(&mut r, t);
        traj.push(r);
    }
    traj
}

/// Peso de infectividad τ pasos después de la exposición: ζ_A(τ) β_A + ζ_I(τ) β_I (Ecs. 34-36).
fn infectivity_profile(p: &Params, g: usize, len: usize, beta_a: f64, beta_i: f64) -> Vec<f64> {
    let ep = &p.epi;
    let (mut ze, mut za, mut zi) = (1.0, 0.0, 0.0);
    let mut w = Vec::with_capacity(len);
    for _ in 0..len {
        w.push(za * beta_a + zi * beta_i);
        let (e2, a2, i2) = (
            (1.0 - ep.eta[g]) * ze,
            ep.eta[g] * ze + (1.0 - ep.alpha[g]) * za,
            ep.alpha[g] * za + (1.0 - ep.mu[g]) * zi,
        );
        ze = e2;
        za = a2;
        zi = i2;
    }
    w
}

/// Q^g_i(s): contactos susceptibles de un infectado de (g, i) en el instante s (Ec. 32).
fn kernel(p: &Params, c: &Ctx, rho: &[f64]) -> Vec<f64> {
    let (ng, np) = (p.ng, p.np);
    let mut rt = vec![0.0; ng * np]; // ρ̃^{S,h}_j, Ec. 31
    for h in 0..ng {
        for j in 0..np {
            let s: f64 = (0..np)
                .map(|k| p.pop[h * np + k] * rho[idx(ng, np, S, h, k)] * c.mm[h][k * np + j])
                .sum();
            rt[h * np + j] = if c.nt[h * np + j] > 0.0 { s / c.nt[h * np + j] } else { 0.0 };
        }
    }
    let mut q = vec![0.0; ng * np];
    for g in 0..ng {
        for i in 0..np {
            let mut s = 0.0;
            for j in 0..np {
                for h in 0..ng {
                    s += c.z[g] * c.kg[g] * c.f[j] * p.c[g * ng + h] * c.mm[g][i * np + j] * rt[h * np + j];
                }
            }
            q[g * np + i] = s;
        }
    }
    q
}

struct Repro {
    r: Vec<f64>,         // R(t), Ec. 42
    rg: Vec<Vec<f64>>,   // R_g(t), Ec. 43, índice [g][t]
    expo: Vec<f64>,      // nuevos expuestos totales en t (individuos)
}

/// R^g_i(t) (Ec. 40) para t = 0..=tmax. `traj` debe llegar al menos hasta tmax + horizonte.
fn reproduction(p: &Params, traj: &[Vec<f64>], beta_a: f64, beta_i: f64, tmax: usize) -> Repro {
    let (ng, np) = (p.ng, p.np);
    let q: Vec<Vec<f64>> = (0..traj.len()).map(|s| kernel(p, &p.ctx(s), &traj[s])).collect();
    let w: Vec<Vec<f64>> = (0..ng).map(|g| infectivity_profile(p, g, traj.len(), beta_a, beta_i)).collect();
    let mut rgi = vec![vec![0.0; ng * np]; tmax + 1];
    for t in 0..=tmax {
        for g in 0..ng {
            for i in 0..np {
                rgi[t][g * np + i] = (0..traj.len() - t).map(|tau| w[g][tau] * q[t + tau][g * np + i]).sum();
            }
        }
    }
    let pool = |rho: &[f64], g: usize, i: usize| rho[idx(ng, np, S, g, i)] + rho[idx(ng, np, CH, g, i)];
    let mut r = vec![0.0; tmax + 1];
    let mut rg = vec![vec![0.0; tmax + 1]; ng];
    let mut expo = vec![0.0; tmax + 1];
    for t in 0..=tmax {
        let (mut num, mut den) = (0.0, 0.0);
        let mut numg = vec![0.0; ng];
        let mut deng = vec![0.0; ng];
        for g in 0..ng {
            for i in 0..np {
                let e = if t == 0 {
                    traj[0][idx(ng, np, E, g, i)]
                } else {
                    pool(&traj[t - 1], g, i) - pool(&traj[t], g, i)
                } * p.pop[g * np + i];
                num += e * rgi[t][g * np + i];
                den += e;
                numg[g] += e * rgi[t][g * np + i];
                deng[g] += e;
            }
        }
        r[t] = if den > 0.0 { num / den } else { f64::NAN };
        for g in 0..ng {
            rg[g][t] = if deng[g] > 0.0 { numg[g] / deng[g] } else { f64::NAN };
        }
        expo[t] = den;
    }
    Repro { r, rg, expo }
}

/// Serie κ0(t) sintética con la forma de un confinamiento: subida, meseta y relajación parcial.
fn kappa_series() -> Vec<f64> {
    (0..200)
        .map(|t| {
            let t = t as f64;
            if t < 14.0 {
                0.0
            } else if t < 24.0 {
                0.85 * (t - 14.0) / 10.0
            } else if t < 60.0 {
                0.85
            } else if t < 90.0 {
                0.85 - 0.20 * (t - 60.0) / 30.0
            } else {
                0.65
            }
        })
        .collect()
}

/// Día (primer índice) en que la serie baja de 1 tras haber estado por encima de 1.
fn first_drop_below_one(v: &[f64]) -> Option<usize> {
    let up = v.iter().position(|&x| x > 1.0)?;
    (up..v.len()).find(|&t| v[t] < 1.0)
}

/// Tiempo de generación medio de un perfil de infectividad w(τ).
fn mean_generation_time(w: &[f64]) -> f64 {
    let tot: f64 = w.iter().sum();
    w.iter().enumerate().map(|(t, x)| t as f64 * x).sum::<f64>() / tot
}

/// Estimador de renovación instantáneo (tipo EpiEstim): R̂(t) = e(t) / Σ_{s≥1} e(t-s) w(s).
fn renewal_estimate(expo: &[f64], w: &[f64]) -> Vec<f64> {
    let tot: f64 = w.iter().skip(1).sum();
    (0..expo.len())
        .map(|t| {
            let lam: f64 = (1..=t.min(w.len() - 1)).map(|s| expo[t - s] * w[s] / tot).sum();
            if lam > 1e-12 { expo[t] / lam } else { f64::NAN }
        })
        .collect()
}

fn svg_lines(title: &str, x: &[f64], series: &[(&str, &str, &Vec<f64>)], ymax: f64) -> String {
    let (w, h, l, r, t, b) = (680.0, 440.0, 60.0, 20.0, 40.0, 50.0);
    let xmax = x[x.len() - 1];
    let px = |v: f64| l + (v / xmax) * (w - l - r);
    let py = |v: f64| h - b - (v / ymax) * (h - t - b);
    let mut s = format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {w} {h}\" font-family=\"sans-serif\" font-size=\"12\"><rect width=\"100%\" height=\"100%\" fill=\"#fff\"/>\
         <text x=\"{}\" y=\"22\" text-anchor=\"middle\" font-size=\"14\">{title}</text>",
        w / 2.0
    );
    for k in 0..=4 {
        let v = ymax * k as f64 / 4.0;
        write!(s, "<line x1=\"{l}\" x2=\"{}\" y1=\"{y}\" y2=\"{y}\" stroke=\"#ddd\"/><text x=\"{}\" y=\"{}\" text-anchor=\"end\">{v:.1}</text>", w - r, l - 6.0, py(v) + 4.0, y = py(v)).unwrap();
    }
    for k in 0..=5 {
        let v = xmax * k as f64 / 5.0;
        write!(s, "<text x=\"{}\" y=\"{}\" text-anchor=\"middle\">{v:.0}</text>", px(v), h - b + 18.0).unwrap();
    }
    write!(s, "<line x1=\"{l}\" x2=\"{}\" y1=\"{y}\" y2=\"{y}\" stroke=\"#444\" stroke-dasharray=\"4 3\"/>", w - r, y = py(1.0)).unwrap();
    write!(s, "<text x=\"{}\" y=\"{}\" text-anchor=\"middle\">día</text><text x=\"16\" y=\"{}\" text-anchor=\"middle\" transform=\"rotate(-90 16 {})\">R(t)</text>", w / 2.0, h - 10.0, h / 2.0, h / 2.0).unwrap();
    for (k, (name, color, ys)) in series.iter().enumerate() {
        let pts: Vec<String> = x.iter().zip(ys.iter()).filter(|(_, y)| y.is_finite()).map(|(a, y)| format!("{:.1},{:.1}", px(*a), py(y.min(ymax)))).collect();
        write!(s, "<polyline fill=\"none\" stroke=\"{color}\" stroke-width=\"2\" points=\"{}\"/>", pts.join(" ")).unwrap();
        write!(s, "<rect x=\"{}\" y=\"{}\" width=\"12\" height=\"3\" fill=\"{color}\"/><text x=\"{}\" y=\"{}\">{name}</text>", w - 150.0, t + 6.0 + 16.0 * k as f64, w - 134.0, t + 11.0 + 16.0 * k as f64).unwrap();
    }
    s.push_str("</svg>");
    s
}

/// Escenario con κ constante desde t = 0. Con `recip`, la matriz de contactos se simetriza para que
/// n_g <k_g> C_gh = n_h <k_h> C_hg (reciprocidad), condición bajo la cual la Ec. 40 es exacta.
fn const_scenario(recip: bool) -> Params {
    let mut p = synthetic(0.6, 0);
    if recip {
        let kg = p.ctx(0).kg;
        let ng = p.ng;
        let n: Vec<f64> = (0..ng).map(|g| (0..p.np).map(|i| p.pop[g * p.np + i]).sum()).collect();
        let fl = |g: usize, h: usize| n[g] * kg[g] * p.c[g * ng + h];
        let mut c = vec![0.0; ng * ng];
        for g in 0..ng {
            for h in 0..ng {
                c[g * ng + h] = 0.5 * (fl(g, h) + fl(h, g)) / (n[g] * kg[g]);
            }
        }
        p.c = c;
    }
    p
}

/// Σ_t e(t) R(t) / Σ_{t≥1} e(t): cada nueva exposición tiene un único infectante, así que el cociente vale 1.
fn identity_ratio(p: &Params, traj: &[Vec<f64>], beta_a: f64, beta_i: f64) -> f64 {
    let tmax = traj.len() - 150;
    let rep = reproduction(p, traj, beta_a, beta_i, tmax);
    let lhs: f64 = (0..=tmax).filter(|&t| rep.expo[t] > 0.0).map(|t| rep.expo[t] * rep.r[t]).sum();
    let rhs: f64 = (1..=tmax).map(|t| rep.expo[t]).sum();
    lhs / rhs
}

fn fig3_scenario() -> (Params, Vec<Vec<f64>>) {
    let mut p = synthetic(0.0, 14);
    p.kappa = kappa_series();
    let rho0 = init(&p, 1e-3, 0, 1);
    let traj = simulate(&p, rho0, 200);
    (p, traj)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    match args.get(1).map(String::as_str) {
        Some("kappa") => {
            let beta_a = 0.5 * 0.075;
            let beta_i = 0.075;
            let (ea, ei) = (-(1.0f64 - beta_a).ln(), -(1.0f64 - beta_i).ln());
            println!("kappa  R(Ec.45)  R(β efectiva)  crecimiento/día (simulado)");
            for k in 0..=9 {
                let kap = k as f64 * 0.1;
                println!(
                    "{:.1}    {:.4}    {:.4}         {:+.5}",
                    kap,
                    r_static(kap, beta_a, beta_i),
                    r_static(kap, ea, ei),
                    growth_rate(kap)
                );
            }
            println!("\nκ crítico:");
            println!("  Ec. 51 con β del artículo   : {:.4}", kappa_c_analytic(beta_a, beta_i));
            println!("  Ec. 51 con β efectiva -ln(1-β): {:.4}", kappa_c_analytic(ea, ei));
            println!("  simulación (crecimiento = 0)  : {:.4}", kappa_c_numeric());
        }
        Some("demo") => {
            let out = args.get(2).expect("falta <salida.csv>");
            let p = synthetic(0.6, 30);
            let mut rho = init(&p, 1e-4, 0, 1);
            let mut csv = String::from("dia,expuestos_nuevos,infectados,fallecidos,suma_max_error\n");
            let tot = |rho: &[f64], m: usize| -> f64 {
                let mut s = 0.0;
                for g in 0..p.ng {
                    for i in 0..p.np {
                        s += rho[idx(p.ng, p.np, m, g, i)] * p.pop[g * p.np + i];
                    }
                }
                s
            };
            let mut prev_s = tot(&rho, S) + tot(&rho, CH);
            for t in 0..200 {
                p.step(&mut rho, t);
                let s_now = tot(&rho, S) + tot(&rho, CH);
                let mut err: f64 = 0.0;
                for g in 0..p.ng {
                    for i in 0..p.np {
                        let sum: f64 = (0..NM).map(|m| rho[idx(p.ng, p.np, m, g, i)]).sum();
                        err = err.max((sum - 1.0).abs());
                    }
                }
                writeln!(
                    csv,
                    "{},{:.3},{:.3},{:.3},{:.2e}",
                    t + 1,
                    prev_s - s_now,
                    tot(&rho, I),
                    tot(&rho, D),
                    err
                )
                .unwrap();
                prev_s = s_now;
            }
            std::fs::write(out, csv).unwrap();
        }
        Some("fig3") => {
            let (csv_path, svg_path) = (args.get(2).expect("falta <salida.csv>"), args.get(3).expect("falta <salida.svg>"));
            let (p, traj) = fig3_scenario();
            let (ba, bi) = (p.beta_a, p.beta_i);
            let rep = reproduction(&p, &traj, ba, bi, 100);
            let mut csv = String::from("dia,kappa,R,R_jovenes,R_adultos,R_mayores\n");
            for t in 0..=100 {
                writeln!(csv, "{t},{:.4},{:.4},{:.4},{:.4},{:.4}", p.kap_at(t), rep.r[t], rep.rg[0][t], rep.rg[1][t], rep.rg[2][t]).unwrap();
            }
            std::fs::write(csv_path, csv).unwrap();
            let x: Vec<f64> = (0..=100).map(|t| t as f64).collect();
            let svg = svg_lines(
                "R(t) en un escenario sintético (κ0(t) inventado, no el de España)",
                &x,
                &[("R global", "#111", &rep.r), ("jóvenes", "#1f77b4", &rep.rg[0]), ("adultos", "#d62728", &rep.rg[1]), ("mayores", "#2ca02c", &rep.rg[2])],
                rep.rg.iter().flatten().chain(rep.r.iter()).cloned().filter(|v| v.is_finite()).fold(1.0, f64::max).ceil(),
            );
            std::fs::write(svg_path, svg).unwrap();
        }
        Some("validar") => {
            println!("Identidad de conteo: Σ_t e(t)·R(t) / Σ_(t≥1) e(t), confinamiento constante κ = 0.6, régimen subcrítico");
            let (ba, bi) = (0.5 * 0.075, 0.075);
            let (ea, ei) = (-(1.0f64 - ba).ln(), -(1.0f64 - bi).ln());
            for (nombre, recip) in [("matriz C del artículo", false), ("C simetrizada (recíproca)", true)] {
                let p = const_scenario(recip);
                let traj = simulate(&p, init(&p, 1e-3, 0, 1), 400);
                for (bn, a, b) in [("β del artículo", ba, bi), ("β efectiva", ea, ei)] {
                    println!("  {nombre}, {bn}: {:.4}", identity_ratio(&p, &traj, a, b));
                }
            }
            let (p, traj) = fig3_scenario();
            let rep = reproduction(&p, &traj, ba, bi, 100);
            let w = infectivity_profile(&p, 1, 100, ba, bi);
            let ren = renewal_estimate(&rep.expo, &w);
            println!(
                "\nEscenario de la Fig. 3 sintética: R(t) analítico baja de 1 el día {:?}; estimador de renovación el día {:?}; tiempo de generación medio (adultos) {:.1} d",
                first_drop_below_one(&rep.r),
                first_drop_below_one(&ren),
                mean_generation_time(&w)
            );
        }
        _ => eprintln!("uso: covid_mmca kappa | demo <salida.csv> | fig3 <salida.csv> <salida.svg> | validar"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn conserva_probabilidad_y_no_negativos() {
        let p = synthetic(0.5, 30);
        let mut rho = init(&p, 1e-3, 1, 1);
        for t in 0..120 {
            p.step(&mut rho, t);
            assert!(rho.iter().all(|&x| x >= -1e-15), "negativo en t={t}");
            for g in 0..p.ng {
                for i in 0..p.np {
                    let s: f64 = (0..NM).map(|m| rho[idx(p.ng, p.np, m, g, i)]).sum();
                    assert!((s - 1.0).abs() < 1e-12, "suma {s} en t={t}");
                }
            }
        }
    }

    #[test]
    fn sin_siembra_no_hay_epidemia() {
        let p = synthetic(0.0, usize::MAX);
        let mut rho = init(&p, 0.0, 0, 0);
        for t in 0..50 {
            p.step(&mut rho, t);
        }
        assert!((0..p.ng * p.np).all(|k| (rho[idx(p.ng, p.np, S, 0, 0) + k] - 1.0).abs() < 1e-15));
    }

    fn pool_total(p: &Params, rho: &[f64]) -> f64 {
        let (ng, np) = (p.ng, p.np);
        let mut s = 0.0;
        for g in 0..ng {
            for i in 0..np {
                s += (rho[idx(ng, np, S, g, i)] + rho[idx(ng, np, CH, g, i)]) * p.pop[g * np + i];
            }
        }
        s / p.pop.iter().sum::<f64>()
    }

    fn permutar(p: &Params, perm: &[usize]) -> Params {
        let (ng, np) = (p.ng, p.np);
        let mut q = p.clone();
        for g in 0..ng {
            for (new, &old) in perm.iter().enumerate() {
                q.pop[g * np + new] = p.pop[g * np + old];
                for (new2, &old2) in perm.iter().enumerate() {
                    q.mob[g][new * np + new2] = p.mob[g][old * np + old2];
                }
            }
        }
        q.area = perm.iter().map(|&o| p.area[o]).collect();
        q
    }

    #[test]
    fn invariante_a_reetiquetar_parches() {
        let p = synthetic(0.6, 14);
        let perm = [2usize, 0, 3, 1];
        let q = permutar(&p, &perm);
        let new_seed = perm.iter().position(|&o| o == 0).unwrap();
        let a = simulate(&p, init(&p, 1e-3, 0, 1), 150);
        let b = simulate(&q, init(&q, 1e-3, new_seed, 1), 150);
        let da = pool_total(&p, &a[150]);
        let db = pool_total(&q, &b[150]);
        println!("reetiquetado: fracción susceptible final {da:.12} frente a {db:.12}");
        assert!((da - db).abs() < 1e-12);
    }

    #[test]
    fn invariante_a_escalar_poblacion_y_superficie() {
        let p = synthetic(0.6, 14);
        let mut q = p.clone();
        q.pop.iter_mut().for_each(|x| *x *= 10.0);
        q.area.iter_mut().for_each(|x| *x *= 10.0);
        let a = simulate(&p, init(&p, 1e-3, 0, 1), 150);
        let b = simulate(&q, init(&q, 1e-3, 0, 1), 150);
        let d = a[150].iter().zip(&b[150]).map(|(x, y)| (x - y).abs()).fold(0.0, f64::max);
        println!("escala x10: diferencia máxima de probabilidades {d:.2e}");
        assert!(d < 1e-12);
    }

    #[test]
    fn tamano_final_frente_a_la_ecuacion_clasica() {
        let (ba, bi) = (0.5 * 0.075, 0.075);
        let (ea, ei) = (-(1.0f64 - ba).ln(), -(1.0f64 - bi).ln());
        let mut p = single_patch(0.0);
        p.tc = usize::MAX;
        let r0 = (ea / p.epi.alpha[0] + ei / p.epi.mu[0]) * p.k_tot[0];
        let traj = simulate(&p, init(&p, 1e-6, 0, 0), 500);
        let sim = 1.0 - pool_total(&p, &traj[500]);
        let mut z = 0.9;
        for _ in 0..200 {
            z = 1.0 - (-r0 * z).exp();
        }
        println!("tamaño final: R0 = {r0:.4}, simulado {sim:.4}, ecuación 1-exp(-R0 z) {z:.4}");
        assert!((sim - z).abs() < 0.01);
    }

    #[test]
    fn mas_confinamiento_menos_casos() {
        let mut prev = f64::INFINITY;
        for kap in [0.0, 0.2, 0.4, 0.6, 0.8] {
            let p = synthetic(kap, 14);
            let traj = simulate(&p, init(&p, 1e-3, 0, 1), 300);
            let ar = 1.0 - pool_total(&p, &traj[300]);
            println!("κ = {kap:.1}: tasa de ataque {ar:.5}");
            assert!(ar <= prev + 1e-12);
            prev = ar;
        }
    }

    #[test]
    fn r_no_depende_del_tamano_de_la_siembra() {
        let (ba, bi) = (0.5 * 0.075, 0.075);
        let p = single_patch(0.4);
        let r = |seed: f64| reproduction(&p, &simulate(&p, init(&p, seed, 0, 0), 200), ba, bi, 30).r[30];
        let (a, b) = (r(1e-12), r(1e-5));
        println!("R(30) con siembra 1e-12: {a:.6}; con 1e-5: {b:.6}");
        assert!((a - b).abs() / a < 1e-3);
    }

    #[test]
    fn r_estatico_coincide_con_la_ec_45() {
        let (ba, bi) = (0.5 * 0.075, 0.075);
        let mut p = single_patch(0.4);
        p.pop = vec![1.0e6];
        let traj = simulate(&p, init(&p, 1e-12, 0, 0), 200);
        let rep = reproduction(&p, &traj, ba, bi, 30);
        let est = r_static(0.4, ba, bi);
        assert!((rep.r[30] - est).abs() / est < 1e-4, "R(30) = {} frente a Ec. 45 = {est}", rep.r[30]);
    }

    #[test]
    fn identidad_de_conteo_con_contactos_reciprocos() {
        let (ba, bi) = (0.5 * 0.075, 0.075);
        let (ea, ei) = (-(1.0f64 - ba).ln(), -(1.0f64 - bi).ln());
        let p = const_scenario(true);
        let traj = simulate(&p, init(&p, 1e-3, 0, 1), 400);
        let q = identity_ratio(&p, &traj, ea, ei);
        assert!((q - 1.0).abs() < 0.01, "cociente {q}");
    }

    #[test]
    fn cruce_de_r_igual_a_uno_frente_a_estimador_de_renovacion() {
        let (p, traj) = fig3_scenario();
        let (ba, bi) = (p.beta_a, p.beta_i);
        let rep = reproduction(&p, &traj, ba, bi, 100);
        let w = infectivity_profile(&p, 1, 100, ba, bi);
        let ren = renewal_estimate(&rep.expo, &w);
        let (a, b) = (first_drop_below_one(&rep.r).unwrap() as f64, first_drop_below_one(&ren).unwrap() as f64);
        assert!((a - (b - mean_generation_time(&w))).abs() <= 2.0, "analítico {a}, renovación {b}");
    }

    #[test]
    fn kappa_critico_coincide_con_la_ecuacion_51() {
        let (ba, bi) = (0.5 * 0.075, 0.075);
        let (ea, ei) = (-(1.0f64 - ba).ln(), -(1.0f64 - bi).ln());
        let num = kappa_c_numeric();
        let ana = kappa_c_analytic(ea, ei);
        assert!((num - ana).abs() < 0.01, "numérico {num} frente a analítico {ana}");
    }
}
