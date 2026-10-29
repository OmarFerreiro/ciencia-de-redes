# covid_mmca: núcleo del modelo MMCA (Arenas et al., Phys. Rev. X 10, 041055, 2020)

Reimplementación independiente en Rust de las Ecs. 2-27 del artículo, con **datos sintéticos** (no usa datos reales de España). El código oficial de los autores está en Julia (`MMCAcovid19.jl`, AGPL-3); este crate no deriva de él ni se ha comparado con él.

**Alcance actual:** dinámica de diez compartimentos más el compartimento de hogares confinados (CH), tres edades, movilidad solo en adultos, confinamiento κ constante desde un día t_c, distanciamiento δ y permeabilidad φ. **No incluye:** datos del INE, calibración ABC, cálculo de R(t) general (Ec. 42), ni el material suplementario.

## Comprobaciones (`cargo test --release`)

1. En cada paso, la suma de los once compartimentos de cada parche y edad es 1 (error máximo 2.4·10⁻¹⁵ en 200 días) y no hay valores negativos.
2. Sin siembra, la población permanece susceptible.
3. **κ crítico:** parche único, un estrato (adultos), ⟨k_h⟩ = 3.2, ⟨k_h+w⟩ = 8.4, δ = 0.207, φ = 0.174. La Ec. 51 da κ₀ᶜ = 0.563 con la β del artículo. La simulación (pendiente de ln(E+A+I) igual a cero) da 0.578. La diferencia se explica por la forma de la probabilidad de contagio: el modelo usa 1 − (1 − β)^x, que para x pequeño es x·(−ln(1 − β)) y no x·β (un 4 % mayor con β = 0.075). Con β efectiva = −ln(1 − β) la Ec. 51 da 0.578, igual que la simulación (diferencia < 0.001). Es una precisión sobre la fórmula del artículo, que es una aproximación, no un error.

Resultados en `resultados/kappa.txt` (tabla R(κ) y crecimiento simulado) y `resultados/demo.csv` (escenario sintético de 4 parches, κ = 0.6 desde el día 30).

## Notas de implementación

- Las Tablas II y III del artículo rotulan como "days⁻¹" valores que son duraciones en días; aquí la tasa por paso es 1/duración (ver ficha 02).
- `P = 1 − exp(x)` se calcula con `exp_m1`: con infecciones del orden de 10⁻²⁰, `1 − exp(x)` se redondea a cero y falsea el crecimiento inicial.

## Uso

```bash
CARGO_TARGET_DIR=$HOME/rt/cm cargo test --release
CARGO_TARGET_DIR=$HOME/rt/cm cargo run --release -- kappa
CARGO_TARGET_DIR=$HOME/rt/cm cargo run --release -- demo resultados/demo.csv
```
