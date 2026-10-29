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

## Control de convergencia (ER y BA)

Datos: `codigo/explosive_sync/resultados/conv/` (generados con `run_conv.sh`, 8 redes por caso, δλ = 0.02). Valores de r (media entre redes), rama hacia adelante / rama hacia atrás:

| λ | ER base (N=1000) | ER transitorio ×5 | ER dt/2 | ER N=4000 |
|---|---|---|---|---|
| 1.28 | 0.109 / 0.245 | 0.096 / 0.170 | 0.109 / 0.245 | 0.048 / 0.055 |
| 1.32 | 0.201 / 0.580 | 0.290 / 0.536 | 0.201 / 0.580 | 0.066 / 0.571 |
| 1.36 | 0.537 / 0.672 | 0.671 / 0.674 | 0.537 / 0.672 | 0.092 / 0.675 |
| 1.40 | 0.765 / 0.765 | 0.765 / 0.766 | 0.765 / 0.765 | 0.168 / 0.735 |
| 1.44 | 0.806 / 0.806 | 0.806 / 0.806 | 0.806 / 0.806 | 0.533 / 0.783 |
| 1.48 | 0.828 / 0.828 | 0.828 / 0.828 | 0.828 / 0.828 | 0.816 / 0.816 |

- **Paso de integración:** con dt a la mitad los resultados son iguales en los tres decimales mostrados. Para ER, el paso no es la causa.
- **Transitorio:** con transitorio ×5 la rama hacia adelante se acerca a la de vuelta (en λ = 1.36 pasa de 0.54 a 0.67) y la histéresis se estrecha, pero no desaparece del todo: en λ = 1.32 sigue habiendo 0.29 frente a 0.54. Con transitorio corto la histéresis de ER estaba **sobrestimada**.
- **Tamaño del sistema:** con N = 4000 y el transitorio corto (40), la histéresis es **más ancha** (λ ≈ 1.32–1.48) que con N = 1000. Esa ejecución no es comparable a la del transitorio largo, porque cerca de la transición el sistema relaja lentamente y 40 unidades pueden no bastar. Falta N = 4000 con transitorio largo (unos 45 min estimados).
- **BA:** con transitorio ×5 la rama hacia atrás no cambia y la de adelante salta en λ = 1.64 en todas las redes (desviación entre redes ≈ 0 a partir de ahí). La histéresis de BA (λ ≈ 1.3–1.64) es robusta.

**Conclusión provisional:** en BA la histéresis es robusta frente a dt y transitorio. En ER, parte de la histéresis era un efecto del transitorio corto, pero queda un resto que no se explica con lo medido. No se puede afirmar todavía ni que ER sea continuo ni que no lo sea. Falta N mayor con transitorio largo, y comparar con la Fig. 1a del artículo.

## Comparación con la Fig. 1 del artículo

Valores del artículo leídos a ojo de la imagen de ar5iv (error de lectura ≈ ±0.03 en λ y ±0.03 en r). La imagen no se copia al repositorio. Datos propios en `codigo/explosive_sync/resultados/fig1/` (`run_fig1.sh`; N = 1000, ⟨k⟩ = 6, 4 redes, δλ = 0.04, transitorio 40).

**Generador.** El artículo remite al modelo BA–ER de Gómez-Gardeñes y Moreno, Phys. Rev. E 73, 056124 (2006), [arXiv:cond-mat/0604100](https://arxiv.org/pdf/cond-mat/0604100). Lo reimplementé según su descripción: N fijo; cada nodo es una vez fuente de m enlaces; destino uniforme entre todos los nodos con probabilidad α y preferencial con 1 − α. **Supuestos míos:** clique inicial de m + 1 nodos, preferencial del modelo A con A = m, sin enlaces repetidos ni autoenlaces. No sé qué variante del modelo usó el artículo de 2011.

Una primera versión mía del generador (red que crece eligiendo destinos solo entre nodos ya presentes) daba hubs mayores y resultados distintos de la Fig. 1. Los resultados de las secciones de histéresis de ER y de convergencia de ER corresponden a ese generador, no al del artículo.

| Panel | Red | Artículo (lectura aproximada) | Mis simulaciones (generador de 2006) |
|---|---|---|---|
| (a) | α = 1 | Continua, sube hacia λ ≈ 0.62; r ≈ 0.9 en λ = 1.2 | Continua, sin histéresis; r = 0.44 en λ = 0.60; r = 0.92 en λ = 1.16 |
| (b) | α = 0.6 | Continua, hacia λ ≈ 0.77 | Continua; r = 0.60 en λ = 0.80; r = 0.92 en λ = 1.36 |
| (c) | α = 0.2 | Casi abrupta, λ ≈ 1.22; r ≈ 0.94 en λ = 1.8 | Muy empinada, r = 0.65 en λ = 1.20 y 0.84 en 1.28; r = 0.94 en λ = 1.76; sin histéresis visible |
| (d) | α = 0, BA | Caída hacia atrás en λ ≈ 1.28 y salto hacia adelante en λ ≈ 1.42 | Caída hacia atrás entre λ = 1.20 y 1.28 (r = 0.64 en 1.28); salto hacia adelante entre 1.5 y más de 1.76 según la red |

**Coincidencia:** en (a), (b) y (c), tanto la posición de la transición como los valores de r a λ grande coinciden con el artículo dentro del error de lectura. En (d), la rama hacia atrás coincide.

**Diferencia que queda:** en (d) mi rama hacia adelante salta más tarde que la del artículo (1.42) y con gran dispersión entre redes (desviación ≈ 0.38 en λ = 1.76: algunas redes aún no habían saltado). El artículo muestra una sola red. Con transitorio corto, el estado incoherente metaestable tarda en desestabilizarse. Falta comprobar con transitorio mayor y con un barrido más fino.

**Pendiente:** Fig. 3 (γ = 2.4–3.3 y frecuencias aleatorias), Fig. 4b (hasta K = 70, ya medida de forma parcial) y la comparación de la rama hacia adelante de (d) con mayor transitorio.

## Posibles aportes (a valorar, sin afirmar novedad)

- Código abierto y reproducible con semilla del protocolo completo, que el artículo no detalla.
- Barrido fino con muchas realizaciones de red, en paralelo.
- Estudio de convergencia numérica: integrador, dt, tiempo de transitorio.
- Variantes citadas en el artículo (ω ~ k^a, correlación parcial p).

## Notas técnicas

- El paso dt se limita por la rigidez del nodo hub (dt·λ·k_max ≈ 1). Con senos y cosenos precalculados el barrido grueso de BA bajó de unos 8 minutos a 2. El barrido fino (61 puntos, 8 redes, 3 casos) tarda unos 15–20 minutos en 4 núcleos.
- Las rutas largas de Windows rompen el enlazador y `git init` en la carpeta temporal. Se compila en WSL con `CARGO_TARGET_DIR=$HOME/rt/es` y el código se lee desde `/mnt/c/...`.
- `run_fine.sh` reproduce todos los datos y las gráficas de esta ficha.
