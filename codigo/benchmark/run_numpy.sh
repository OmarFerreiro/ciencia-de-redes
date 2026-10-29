cd "$(dirname "$0")"
P=$HOME/bench/venv/bin/python
G=$HOME/bench/red.txt
echo "# NumPy+SciPy: $($P -c 'import numpy,scipy;print("numpy",numpy.__version__,"scipy",scipy.__version__)')"
echo "# Verificación (1000 pasos)"
$HOME/rt/es/release/explosive_sync bench $G 1.5 0.002 1000
$P bench_numpy.py $G 1.5 0.002 1000
echo "# Tiempos (20000 pasos, 5 repeticiones, intercalado con Rust)"
for i in 1 2 3 4 5; do $P bench_numpy.py $G 1.5 0.002 20000; $HOME/rt/es/release/explosive_sync bench $G 1.5 0.002 20000; done
