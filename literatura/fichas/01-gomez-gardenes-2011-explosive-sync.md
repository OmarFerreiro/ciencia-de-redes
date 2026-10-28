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

Código: `codigo/explosive_sync/` (RK4, grafos CSR, generador ChaCha8 con semilla). Datos en `codigo/explosive_sync/resultados/`.

**Condiciones:** N = 1000, ⟨k⟩ ≈ 6 (m = 3), semilla 1, δλ = 0.25 (más grueso que el 0.02 del artículo), transitorio 40 y promedio 10 unidades de tiempo.

| Caso | Resultado |
|---|---|
| BA (α = 0), k_max = 88 | Rama hacia adelante: r ≈ 0.07 en λ = 1.5 y r = 0.946 en λ = 1.75. Rama hacia atrás: r = 0.918 en λ = 1.5, y r ≈ 0.07 en λ = 1.25. **Hay histéresis**, con una ventana entre λ ≈ 1.25 y 1.75. |
| ER (α = 1), k_max = 25 | Salto de r entre λ = 1.25 y 1.5. Las dos ramas coinciden con esta resolución (sin histéresis visible). |
| Estrella K = 10 | En la rama hacia atrás, el hub pierde el enganche entre λ = 0.80 y 0.85, coherente con λ_c = 0.818. |

**Lo que se concluye y lo que no:**
- Cualitativamente se reproduce lo esencial: la histéresis aparece en BA con ω = k.
- Con δλ = 0.25 no se puede medir λ_c con precisión, ni afirmar que ER es de segundo orden. Falta refinar el barrido (δλ = 0.02), promediar sobre varias redes y comprobar la convergencia en dt y en el tiempo de transitorio.
- En la estrella, las hojas tienen todas la misma frecuencia, así que permanecen en fase entre sí aunque el hub se desenganche. Por eso r de la rama hacia atrás no cae a 0 y no basta para ver el salto. Para validar la teoría habrá que medir también la frecuencia efectiva del hub (Ec. 3). Eso queda pendiente.
- No se ha hecho ninguna figura (no hay matplotlib instalado).

## Posibles aportes (a valorar, sin afirmar novedad)

- Código abierto y reproducible con semilla del protocolo completo, que el artículo no detalla.
- Barrido fino con muchas realizaciones de red, en paralelo (rayon).
- Estudio de convergencia numérica: integrador, dt, tiempo de transitorio.
- Variantes citadas en el artículo (ω ~ k^a, correlación parcial p).

## Notas técnicas

- El paso dt se limita por la rigidez del nodo hub (dt·λ·k_max ≲ 1), lo que encarece los barridos en BA: unos 8 min para 25 puntos, hacia adelante y hacia atrás. Hay que acelerar con un integrador adaptativo o con paralelismo antes de afinar el barrido.
- En Windows, la ruta de la carpeta es demasiado larga para el enlazador. Se compila con `CARGO_TARGET_DIR` en una ruta corta.
