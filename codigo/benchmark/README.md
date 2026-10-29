# Benchmark: Rust, C++ y Python (Kuramoto con ω = k, RK4)

**Problema (idéntico en las tres versiones):** red BA (modelo de 2006, α = 0), N = 1000, ⟨k⟩ = 6, semilla 1; λ = 1.5; dt = 0.002; fases iniciales deterministas θ_i = 2π·frac(i·0.618…). Cada versión implementa el mismo algoritmo y el mismo orden de operaciones (senos y cosenos por nodo, suma sobre vecinos CSR, RK4).

**Entorno:** Intel Core i5-7200U (4 hilos, 2.50 GHz), WSL2 Ubuntu 26.04, rustc 1.99.1 (`--release`), g++ 15.2.0 (`-O3`), Python 3.14.4 (puro) y Python con NumPy 2.5.4 + SciPy 1.18.1. Un solo hilo en todas las versiones.

## Verificación

Tras 1000 pasos, las tres versiones dan **los mismos valores a 12 cifras**: r = 0.121070072935 y θ₀ = 121.641836800533. Después de 20000 pasos coinciden Rust y C++ (r = 0.028276645414). Python solo se midió con 300 pasos.
**Verificación.** Tras 1000 pasos, las cuatro versiones dan **los mismos valores a 12 cifras**: r = 0.121070072935 y θ₀ = 121.641836800533. Tras 20000 pasos coinciden Rust y C++ (r = 0.028276645414). NumPy da r = 0.0688 y θ₀ = 2444.6 (frente a 2450.9): las sumas por vecinos se hacen en otro orden y la dinámica, cerca del estado incoherente, amplifica esas diferencias de redondeo. Por tanto, solo hay verificación a 12 cifras hasta los 1000 pasos. Python puro se midió con 300 pasos.

## Tiempos (µs por paso RK4; mediana de las repeticiones)

| Versión | µs/paso | Repeticiones (min–máx) |
|---|---|---|
| Rust (`sin_cos`) | ≈ 170 | 162–188 (13 medidas) |
| C++ `-O3` con `sin` y `cos` por separado | ≈ 209 | 200–228 (8 medidas) |
| C++ `-O3` con `sincos` | ≈ 177 | 169–203 (5 medidas) |
| Python con NumPy + SciPy (matriz dispersa CSR) | ≈ 403 | 403–425 (5 medidas) |
| Python puro (sin NumPy) | ≈ 3260 | 3251–3318 (3 medidas) |

## Qué se puede concluir y qué no

- **Rust y C++ rinden lo mismo.** La diferencia inicial (≈ 20 %) venía de cómo se calculan senos y cosenos: Rust llama a `sin_cos` (una llamada) y mi primera versión de C++ llamaba a `sin` y `cos` por separado. Con `sincos` en C++ la diferencia entra en el ruido de las medidas (la variación entre repeticiones es del 10–15 %).
- **NumPy con matriz dispersa es unas 2,4 veces más lento que Rust** en este problema, y **Python puro unas 20 veces más lento.** Con N = 1000 el coste de NumPy está dominado por la sobrecarga de cada llamada (cada paso hace decenas de operaciones sobre vectores pequeños). Con redes mucho mayores la diferencia podría reducirse. No se ha medido.
- **No se midió el paralelismo** (los barridos reales usan varios hilos para repartir redes distintas). NumPy y SciPy tampoco se probaron con varios hilos.
- **Una sola máquina y una sola red.** Los valores sirven como orden de magnitud en este equipo.
- Rust no gana en velocidad a C++. Lo que ofrece es seguridad de memoria y un paralelismo más sencillo de escribir sin errores, y eso no lo ha medido este benchmark.

## Reproducir

```bash
bash run_bench.sh          # construye, verifica y mide Rust, C++ y Python puro -> resultados.txt
bash run_cpp_sincos.sh     # variante de C++ con sincos -> resultados_sincos.txt
bash setup_venv.sh         # crea el entorno virtual con NumPy y SciPy (en ~/bench/venv)
bash run_numpy.sh          # NumPy + SciPy intercalado con Rust -> resultados_numpy.txt
```

`setup_venv.sh` necesita el paquete `python3-venv` de Ubuntu (`sudo apt install python3-venv`). Los resultados en bruto están en `resultados.txt`, `resultados_sincos.txt` y `resultados_numpy.txt`.
