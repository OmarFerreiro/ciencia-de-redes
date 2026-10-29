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

Código: `codigo/explosive_sync/` (RK4, grafos CSR, ChaCha8 con semilla, hilos de la biblioteca estándar). Datos y gráficas en `codigo/explosive_sync/resultados/`.

**Condiciones comunes:** N = 1000, ⟨k⟩ = 6 (m = 3), ω_i = k_i, λ sin normalizar, continuación adiabática hacia adelante y hacia atrás. Transitorio de 40 y promedio de 10 unidades de tiempo (salvo que se indique), dt = min(0.01, 1/(λ·k_max)).

**Qué se comparó y dónde se detalla:**
- La estrella frente a la teoría exacta: sección siguiente.
- Los cuatro paneles de la Fig. 1: sección "Comparación con la Fig. 1".
- El control de dt y del transitorio: sección "Control de convergencia con el generador final".

**Pendiente:**
- Guardar el λ de salto de cada red en BA y comparar su distribución con el valor publicado.
- Fig. 3 (γ = 2.4–3.3 y frecuencias aleatorias), Fig. 4b (rama hacia adelante para K = 20–70) y las notas 24 y 27 (ω ~ k^a, correlación parcial).
- N mayor con transitorio largo.

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

## Comparación con la Fig. 1 del artículo

Valores del artículo leídos a ojo de la imagen de ar5iv (error de lectura ≈ ±0.03 en λ y ±0.03 en r). La imagen no se copia al repositorio. Datos propios en `codigo/explosive_sync/resultados/fig1/` (`run_fig1.sh`; N = 1000, ⟨k⟩ = 6, 4 redes, δλ = 0.04, transitorio 40).

**Generador.** El artículo remite al modelo BA–ER de Gómez-Gardeñes y Moreno, Phys. Rev. E 73, 056124 (2006), [arXiv:cond-mat/0604100](https://arxiv.org/pdf/cond-mat/0604100). Lo reimplementé según su descripción: N fijo; cada nodo es una vez fuente de m enlaces; destino uniforme entre todos los nodos con probabilidad α y preferencial con 1 − α. **Supuestos míos:** clique inicial de m + 1 nodos, preferencial del modelo A con A = m, sin enlaces repetidos ni autoenlaces. No sé qué variante del modelo usó el artículo de 2011.

Una primera versión mía del generador (red que crece eligiendo destinos solo entre nodos ya presentes) daba hubs mayores y resultados distintos de la Fig. 1: una red "aleatoria" con transición abrupta y al doble de λ. Los resultados de esa primera versión no están en este repositorio, porque no corresponden al modelo del artículo.

| Panel | Red | Artículo (lectura aproximada) | Mis simulaciones (generador de 2006) |
|---|---|---|---|
| (a) | α = 1 | Continua, sube hacia λ ≈ 0.62; r ≈ 0.9 en λ = 1.2 | Continua, sin histéresis; r = 0.44 en λ = 0.60; r = 0.92 en λ = 1.16 |
| (b) | α = 0.6 | Continua, hacia λ ≈ 0.77 | Continua; r = 0.60 en λ = 0.80; r = 0.92 en λ = 1.36 |
| (c) | α = 0.2 | Casi abrupta, λ ≈ 1.22; r ≈ 0.94 en λ = 1.8 | Muy empinada, r = 0.65 en λ = 1.20 y 0.84 en 1.28; r = 0.94 en λ = 1.76; sin histéresis visible |
| (d) | α = 0, BA | Caída hacia atrás en λ ≈ 1.28 y salto hacia adelante en λ ≈ 1.42 | Caída hacia atrás entre λ = 1.20 y 1.28 (r = 0.64 en 1.28); salto hacia adelante entre 1.5 y más de 1.76 según la red |

**Coincidencia:** en (a), (b) y (c), tanto la posición de la transición como los valores de r a λ grande coinciden con el artículo dentro del error de lectura. En (d), la rama hacia atrás coincide.

**Diferencia que queda (matizada tras el control de convergencia de abajo):** en (d) el salto hacia adelante no ocurre en el mismo λ en todas mis redes: según las medias y desviaciones, unas redes saltan hacia λ ≈ 1.44 y otras hasta λ ≈ 1.64. El valor del artículo (1.42) corresponde a una sola red y queda en el borde inferior de ese intervalo. No guardé el salto de cada red por separado, así que esto es una inferencia a partir de las medias y desviaciones, no una medida directa.

**Pendiente:** Fig. 3 (γ = 2.4–3.3 y frecuencias aleatorias), Fig. 4b (hasta K = 70, ya medida de forma parcial) y guardar el λ de salto de cada red por separado en BA, para comparar su distribución con el valor publicado.

## Control de convergencia con el generador final

Datos: `codigo/explosive_sync/resultados/conv_final/` (generados con `run_conv_final.sh`; N = 1000, ⟨k⟩ = 6, 8 redes por caso, δλ = 0.02). Variantes: base (transitorio 40, dt estándar), transitorio ×5 (200) y dt/2.

- **α = 1 (panel a), λ de 0.4 a 0.9:** dt/2 da las mismas curvas que la base (iguales en los tres decimales mostrados). Con transitorio ×5 las diferencias son de 0.02 como mucho (r = 0.558 frente a 0.555 en λ = 0.62), por debajo de la dispersión entre redes. Las ramas hacia adelante y hacia atrás coinciden: no hay histéresis.
- **BA (α = 0), λ de 1.1 a 1.9:** la rama hacia atrás es la misma en las tres variantes (r ≈ 0.85 en λ = 1.32, 0.89 en 1.40). Su caída se produce entre λ = 1.24 y 1.32 en todas; el valor de r en λ = 1.28 varía entre 0.30 y 0.60 según la variante, con desviación grande entre redes. La rama hacia adelante tiene saltos en λ ≈ 1.44–1.64 según la red; con transitorio ×5 la sincronización completa se alcanza algo antes (λ ≈ 1.60 frente a 1.68), lo que concuerda con un estado incoherente metaestable que se desestabiliza con el tiempo.
- **Conclusión:** para α = 1 los resultados no dependen del paso ni del transitorio dentro de lo medido. Para BA, la rama hacia atrás es robusta, y la posición del salto hacia adelante depende de la red y del transitorio.

## Posibles aportes (a valorar, sin afirmar novedad)

- Código abierto y reproducible con semilla del protocolo completo, que el artículo no detalla.
- Barrido fino con muchas realizaciones de red, en paralelo.
- Estudio de convergencia numérica: integrador, dt, tiempo de transitorio.
- Variantes citadas en el artículo (ω ~ k^a, correlación parcial p).

## Notas técnicas

- El paso dt se limita por la rigidez del nodo hub (dt·λ·k_max ≈ 1). Con senos y cosenos precalculados, un barrido grueso de BA (25 puntos, una red) bajó de unos 8 minutos a 2. Los barridos de convergencia (8 redes, 41 puntos, 3 variantes) tardan de 2 a 26 minutos cada uno en 4 núcleos.
- Se ejecuta en WSL2 (Ubuntu): Windows bloqueaba los ejecutables recién compilados con su política de control de aplicaciones. Se compila con `CARGO_TARGET_DIR=$HOME/rt/es cargo build --release`.
- Los scripts `run_fig1.sh`, `run_conv_final.sh` y `run_starscan.sh` (en `codigo/explosive_sync/`) reproducen los datos de esta ficha.
