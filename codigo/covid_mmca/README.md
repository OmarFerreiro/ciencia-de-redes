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

## R(t) (Ecs. 29-43) y figura tipo Fig. 3

`reproduction()` calcula R^g_i(t) (Ec. 40, con los perfiles ζ de las Ecs. 34-36), el global R(t) (Ec. 42) y el de cada edad R_g(t) (Ec. 43). Comandos: `validar` (comprobaciones) y `fig3 <csv> <svg>` (figura). Resultados en `resultados/validacion_R.txt`, `resultados/fig3_sintetica.csv` y `.svg`.

**La figura NO es la Fig. 3 del artículo:** usa un κ₀(t) inventado (subida, meseta en 0.85, relajación a 0.65) sobre 4 parches sintéticos, no la serie del INE. Solo comprueba que la forma cualitativa sale (R baja de 1 al confinar y sube al relajar; los adultos tienen la R más alta). Los valores de R al inicio quedan vacíos (NaN) en los días sin nuevos expuestos.

Comprobaciones nuevas (`cargo test --release`, 6 en total):
4. **R estático:** en un parche con κ constante, R(30) de la Ec. 40 coincide con la Ec. 45 (error relativo < 10⁻⁴).
5. **Identidad de conteo:** cada nueva exposición tiene un único infectante, por tanto Σ_t e(t)·R(t) = Σ_(t≥1) e(t). Con una matriz de contactos **recíproca** (n_g⟨k_g⟩C_gh = n_h⟨k_h⟩C_hg) y β efectiva = −ln(1−β), el cociente es 0.997. Con la β del artículo es 0.965 (el 3.5 % es la misma diferencia entre β y −ln(1−β) vista en κ crítico).
6. **Estimador de renovación** (tipo EpiEstim, sobre la incidencia simulada): R(t) baja de 1 el día 15; el estimador, el día 21. La diferencia, 6 días, es el tiempo de generación medio (7.1 d en adultos): el R(t) del artículo mira hacia delante (infecciones que causará quien se expone en t) y el estimador instantáneo mira hacia atrás. Corrigiendo ese desfase la diferencia es de 1.1 días (el artículo da unos 2 días con datos reales).

**Hallazgo (a verificar con los autores o el código en Julia):** la Ec. 40/41 usa los contactos del *infectante* (z_g⟨k_g⟩C_gh), mientras que la dinámica (Ec. 13) usa los del *susceptible* (z_h⟨k_h⟩C_hg). Solo coinciden si la matriz de contactos cumple reciprocidad, y la matriz de la Tabla I no la cumple (por ejemplo, n_Y⟨k_Y⟩C_YM ≈ 1.14 frente a n_M⟨k_M⟩C_MY ≈ 1.79 con las fracciones sintéticas). Con esa matriz, la identidad de conteo da 0.993 (con la β del artículo) en el escenario global, pero el R por edades queda sesgado. No sé si el código oficial lo trata distinto.

**Duda abierta:** con los parámetros de las tablas tal como los leo, R al inicio sale en torno a 5 en adultos (0.397 días de infectividad × 13.3 contactos), 3.4 en jóvenes y 2.7 en mayores, por encima del R₀ ≈ 3 que menciona el resumen del artículo. No sé si es una diferencia de definición o de lectura de parámetros por mi parte.
