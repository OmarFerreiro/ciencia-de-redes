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
    kappa: f64,          // confinamiento constante a partir de tc
    tc: usize,
}

fn idx(ng: usize, np: usize, m: usize, g: usize, i: usize) -> usize {
    (m * ng + g) * np + i
}

impl Params {
    fn step(&self, rho: &mut [f64], t: usize) {
        let (ng, np) = (self.ng, self.np);
        let conf = t >= self.tc;
        let kap = if conf { self.kappa } else { 0.0 };

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
        kappa,
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
        kappa,
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
        _ => eprintln!("uso: covid_mmca kappa | demo <salida.csv>"),
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

    #[test]
    fn kappa_critico_coincide_con_la_ecuacion_51() {
        let (ba, bi) = (0.5 * 0.075, 0.075);
        let (ea, ei) = (-(1.0f64 - ba).ln(), -(1.0f64 - bi).ln());
        let num = kappa_c_numeric();
        let ana = kappa_c_analytic(ea, ei);
        assert!((num - ana).abs() < 0.01, "numérico {num} frente a analítico {ana}");
    }
}
