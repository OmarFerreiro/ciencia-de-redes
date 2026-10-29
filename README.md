# Ciencia de Redes

Estudio gradual de los trabajos publicados por el grupo de Biocomputación y Física de Sistemas Complejos (BIFI / COSNET, Universidad de Zaragoza), con el objetivo de:

1. Entender cada trabajo a fondo (modelo, método, resultados).
2. Reproducir sus cálculos y análisis, automatizándolos con **Rust**.
3. Valorar qué aporta la reproducción: rendimiento, verificación independiente, extensiones de parámetros o casos nuevos.
4. Publicar artículos breves en LinkedIn, siempre con referencia al trabajo original y enlace al repositorio con el código.

## Flujo de trabajo por cada paper

1. **Lectura** → ficha en `literatura/fichas/`.
2. **Reproducción** → crate de Rust en `codigo/` (un repositorio por trabajo, o un workspace).
3. **Valoración** → ¿qué se reproduce, qué difiere, qué se puede aportar? (en la ficha).
4. **Artículo corto** → `articulos/AAAA-MM-slug.md`, a partir de `articulos/PLANTILLA.md`.
5. **Publicación** → LinkedIn, con enlace al repo y cita completa del paper.

## Estructura

- `literatura/registro.md`: lista de trabajos y estado de cada uno.
- `literatura/fichas/`: una ficha por trabajo.
- `codigo/`: simulaciones y análisis en Rust (con scripts de gráficos en Python si hacen falta).
- `articulos/`: borradores de los artículos para LinkedIn.

## Reglas de honestidad científica

- Toda referencia se verifica en la fuente (DOI/arXiv) antes de citarla. Nada se cita de memoria.
- Los artículos dicen claramente qué es del grupo original y qué es reproducción o aportación nuestra.
- Los resultados se publican con semilla aleatoria, versión del código y parámetros, para que se puedan reproducir.
- No se afirma que algo sea "inédito" sin una revisión de literatura que lo respalde.

## Reproducir los resultados

Requisitos: Rust (edición 2021) en Linux o WSL2. Desde `codigo/explosive_sync/`:

```bash
CARGO_TARGET_DIR=$HOME/rt/es cargo build --release
bash run_starscan.sh      # grafo estrella, K = 5..70
bash run_fig1.sh          # paneles de la Fig. 1
bash run_conv_final.sh    # control de convergencia
bash run_plots.sh         # gráficas SVG
```

Los benchmarks (Rust, C++, Python) están en `codigo/benchmark/` (ver su README). Los resultados y sus condiciones se detallan en `literatura/fichas/01-gomez-gardenes-2011-explosive-sync.md`.

## Licencia

MIT (ver `LICENSE`).
