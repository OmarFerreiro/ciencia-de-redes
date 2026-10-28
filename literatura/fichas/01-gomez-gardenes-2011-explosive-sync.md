# Ficha 01: Explosive Synchronization Transitions in Scale-free Networks

- **Referencia:** J. Gómez-Gardeñes, S. Gómez, A. Arenas, Y. Moreno, Phys. Rev. Lett. 106, 128701 (2011). DOI 10.1103/PhysRevLett.106.128701 · arXiv:1102.4823
- **Estado:** leído (versión ar5iv) · reproducción parcial en Rust
- **Fecha:** 2026-10-28

## Resumen

Kuramoto sobre redes cuya frecuencia natural es el grado del nodo (ω_i = k_i). En redes libres de escala (BA) la sincronización aparece de forma abrupta, con histéresis, es decir, una transición de primer orden. Con frecuencias asignadas al azar, o en redes menos heterogéneas, la transición es continua.

## Modelo

- Dinámica (Ec. 1): θ̇_i = ω_i + λ Σ_j A_ij sin(θ_j − θ_i). El acoplamiento λ no está normalizado.
- Parámetro de orden (Ec. 2): r e^{iΨ} = (1/N) Σ_j e^{iθ_j}.
- Redes: interpolación BA–ER (α = 1 es ER, α = 0 es BA, P(k) ~ k⁻³), con N = 10³ y ⟨k⟩ = 6. También redes tipo configuración con P(k) ~ k⁻γ, γ = 2.4–3.3.
- Protocolo: continuación adiabática en λ con δλ = 0.02, hacia adelante y hacia atrás. El artículo **no especifica** el integrador, dt, el transitorio ni la ventana de promedio.

## Resultado analítico (grafo estrella, hub de grado K con ω = K, K hojas con ω = 1)

- λ_c = (K − 1)/(K + 1) y r_c = K/(K + 1), con salto finito de r (primer orden).
- Para K = 10: λ_c = 0.818 y r_c = 0.909.

## Reproducción en Rust

Código: `codigo/explosive_sync/` (RK4, grafos CSR, ChaCha8 con semilla, hilos de std). Datos y gráficas en `codigo/explosive_sync/resultados/`. Se ejecuta en WSL2 (Ubuntu), porque Windows bloqueaba los ejecutables nuevos con una política de control de aplicaciones.

**Condiciones:** N = 1000, ⟨k⟩ ≈ 6 (m = 3), 8 redes por caso (semillas 1 a 8), δλ = 0.02 entre λ = 1.0 y 2.2, transitorio 40 y promedio 10 unidades de tiempo, dt ≈ 1/(λ·k_max) (límite 0.01).

| Caso | Resultado |
|---|---|
| **Estrella K = 10** (λ_c = 0.8182 y r_c = 0.9091 según el artículo) | Rama hacia atrás: el hub queda enganchado (frecuencia efectiva = Ω = 1.8182) para λ ≥ 0.82 y no lo está en λ = 0.80. Coincide con λ_c = 0.818. En λ = 0.82, r = 0.920 (r_c = 0.909 en λ_c exacto). Rama hacia adelante: el hub se engancha en λ ≈ 1.60, y hay histéresis. |
| **BA (α = 0)**, ω = k | Rama hacia atrás: r = 0.84 ya en λ = 1.3, y cae a 0.06 en λ = 1.2. Rama hacia adelante: sube entre λ = 1.4 y 1.7, y el punto de salto varía de red a red (desviación entre redes ≈ 0.27–0.37 en el salto). Hay **histéresis** clara, de unos 0.3 en λ. Ver `ba_alpha0_fine.svg`. |
| **ER (α = 1)**, ω = k | Subida rápida entre λ = 1.3 y 1.5. Con δλ = 0.02 aparece una histéresis pequeña (λ = 1.3: r = 0.19 hacia adelante y 0.48 hacia atrás; ambas ramas coinciden desde λ = 1.5). Ver `er_alpha1_fine.svg`. |

**Lo que se concluye:**
- La estrella reproduce el resultado analítico: el desenganche del hub ocurre en λ_c con la resolución usada (δλ = 0.02).
- En BA hay histéresis y el salto es abrupto, como en el artículo.
- **Diferencia con el artículo:** el artículo describe ER como de segundo orden (curvas que coinciden). En nuestras simulaciones ER sube de forma muy empinada y muestra una histéresis pequeña en una ventana estrecha. No sabemos si es un efecto del tamaño finito, del tiempo de transitorio corto (cerca de la transición el sistema tarda más en relajar) o de otra causa. Hay que comprobarlo antes de afirmar nada.

**Pendiente:**
- Convergencia numérica: repetir con TTRANS mayor y con DTF menor, y comprobar que no cambia nada.
- Más redes y N mayor para ver la dependencia con el tamaño.
- Comparar las curvas con la Fig. 1 del artículo en valores concretos de λ (aún no se ha hecho con los datos de la figura).
- Barrer γ y las correlaciones parciales (Fig. 3 y notas 24 y 27).

## Comparación con la teoría de la estrella (K = 5 a 70)

Datos: `codigo/explosive_sync/resultados/starscan.csv` (generados con `run_starscan.sh`, δλ = 0.01, transitorio 200, promedio 50). La rama hacia atrás arranca de un estado sincronizado. La rama enganchada exacta es sin δ = (K−1)/((K+1)λ) y r = |K + e^{iδ}|/(K+1), con cos δ > 0 (derivada aquí, no tomada del artículo).

| K | λ_c teoría | λ_c simulado | diferencia | r_c del artículo K/(K+1) | r_c módulo √(K²+1)/(K+1) |
|---|---|---|---|---|---|
| 5 | 0.6667 | 0.670 | +0.003 | 0.833 | 0.850 |
| 10 | 0.8182 | 0.820 | +0.002 | 0.909 | 0.914 |
| 20 | 0.9048 | 0.910 | +0.005 | 0.952 | 0.954 |
| 30 | 0.9355 | 0.940 | +0.005 | 0.968 | 0.968 |
| 40 | 0.9512 | 0.960 | +0.009 | 0.976 | 0.976 |
| 50 | 0.9608 | 0.970 | +0.009 | 0.980 | 0.981 |
| 60 | 0.9672 | 0.970 | +0.003 | 0.984 | 0.984 |
| 70 | 0.9718 | 0.980 | +0.008 | 0.986 | 0.986 |

- **λ_c:** la simulación queda siempre por encima de la teoría, entre 0.002 y 0.009. Con una malla de 0.01 eso es el orden de la resolución (el valor simulado es el primero de la malla ya enganchado).
- **Rama enganchada:** el error cuadrático medio entre r simulado y la fórmula exacta es menor que 5·10⁻⁶ en todos los K (se muestra como 0.00000). El integrador reproduce la curva analítica.
- **r en λ_c:** el artículo da K/(K+1) (Ec. 10, que calcula r como promedio de cosenos). El módulo exacto es √(K²+1)/(K+1), que es algo mayor (0.5 % para K = 10, ya despreciable para K ≥ 30). La simulación mide el módulo, así que queda por encima de ambos. No es un error del artículo: es una cuestión de definición. Se puede mencionar como matiz, no como corrección.
- **Histéresis:** la rama hacia adelante se engancha en λ = 0.85 (K = 5), 1.56 (K = 10) y 2.34 (K = 20), mientras que λ_c solo sube de 0.67 a 0.90. La histéresis crece con K, de acuerdo con la descripción de la Fig. 4b. Para K ≥ 30 la rama hacia adelante no se engancha antes de λ = 3 (no se ha medido).

## Posibles aportes (a valorar, sin afirmar novedad)

- Código abierto y reproducible con semilla del protocolo completo, que el artículo no detalla.
- Barrido fino con muchas realizaciones de red, en paralelo.
- Estudio de convergencia numérica: integrador, dt, tiempo de transitorio.
- Variantes citadas en el artículo (ω ~ k^a, correlación parcial p).

## Notas técnicas

- El paso dt se limita por la rigidez del nodo hub (dt·λ·k_max ≈ 1). Con senos y cosenos precalculados el barrido grueso de BA bajó de unos 8 minutos a 2. El barrido fino (61 puntos, 8 redes, 3 casos) tarda unos 15–20 minutos en 4 núcleos.
- Las rutas largas de Windows rompen el enlazador y `git init` en la carpeta temporal. Se compila en WSL con `CARGO_TARGET_DIR=$HOME/rt/es` y el código se lee desde `/mnt/c/...`.
- `run_fine.sh` reproduce todos los datos y las gráficas de esta ficha.
