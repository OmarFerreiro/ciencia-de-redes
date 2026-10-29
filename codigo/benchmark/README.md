# Benchmark: Rust, C++ y Python (Kuramoto con ω = k, RK4)

**Problema (idéntico en las tres versiones):** red BA (modelo de 2006, α = 0), N = 1000, ⟨k⟩ = 6, semilla 1; λ = 1.5; dt = 0.002; fases iniciales deterministas θ_i = 2π·frac(i·0.618…). Cada versión implementa el mismo algoritmo y el mismo orden de operaciones (senos y cosenos por nodo, suma sobre vecinos CSR, RK4).

**Entorno:** Intel Core i5-7200U (4 hilos, 2.50 GHz), WSL2 Ubuntu 26.04, rustc 1.99.1 (`--release`), g++ 15.2.0 (`-O3`), Python 3.14.4 sin NumPy. Un solo hilo en todas las versiones.

## Verificación

Tras 1000 pasos, las tres versiones dan **los mismos valores a 12 cifras**: r = 0.121070072935 y θ₀ = 121.641836800533. Después de 20000 pasos coinciden Rust y C++ (r = 0.028276645414). Python solo se midió con 300 pasos.

## Tiempos (µs por paso RK4; mediana de las repeticiones)

| Versión | µs/paso | Repeticiones (min–máx) |
|---|---|---|
| Rust (`sin_cos`) | ≈ 166 | 162–188 (8 medidas) |
| C++ `-O3` con `sin` y `cos` por separado | ≈ 209 | 200–228 (8 medidas) |
| C++ `-O3` con `sincos` | ≈ 177 | 169–203 (5 medidas) |
| Python puro (sin NumPy) | ≈ 3260 | 3251–3318 (3 medidas) |

## Qué se puede concluir y qué no

- **Rust y C++ rinden lo mismo.** La diferencia inicial (≈ 20 %) venía de cómo se calculan senos y cosenos: Rust llama a `sin_cos` (una llamada) y mi primera versión de C++ llamaba a `sin` y `cos` por separado. Con `sincos` en C++ la diferencia entra en el ruido de las medidas (la variación entre repeticiones es del 10–15 %).
- **Python puro es unas 20 veces más lento** que Rust o C++ en este problema. Es la comparación con un código sin vectorizar.
- **Falta NumPy/SciPy.** No estaba instalado en WSL y no lo he instalado. Un Python con matriz dispersa puede quedar bastante más cerca, así que la cifra de 20× no vale como "Python frente a Rust" en general.
- **No se midió el paralelismo** (los barridos reales usan varios hilos para repartir redes distintas).
- **Una sola máquina y una sola red.** Los valores sirven como orden de magnitud en este equipo.
- Rust no gana en velocidad a C++. Lo que ofrece es seguridad de memoria y un paralelismo más sencillo de escribir sin errores, y eso no lo ha medido este benchmark.

## Reproducir

```bash
bash run_bench.sh          # construye, verifica y mide; escribe resultados.txt
bash run_cpp_sincos.sh     # variante de C++ con sincos; resultados_sincos.txt
```

Los resultados en bruto están en `resultados.txt` y `resultados_sincos.txt`.
