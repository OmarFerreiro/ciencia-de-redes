# Benchmark Rust vs C++ vs Python puro: mismo problema (BA, N=1000, <k>=6, lambda=1.5, dt=0.002).
source ~/.cargo/env
cd /mnt/c/Users/jorge/proyectos/ciencia-de-redes/codigo/benchmark
(cd ../explosive_sync && CARGO_TARGET_DIR=$HOME/rt/es cargo build --release 2>&1 | tail -1)
RS=$HOME/rt/es/release/explosive_sync
mkdir -p $HOME/bench
$RS export 0 1000 6 1 $HOME/bench/red.txt
g++ -O3 -std=c++17 bench.cpp -o $HOME/bench/bench_cpp
G=$HOME/bench/red.txt; L=1.5; DT=0.002
{
echo "# Entorno"
lscpu | grep -E "Model name|^CPU\(s\)"
echo "rustc: $(rustc --version)"; echo "g++: $(g++ --version | head -1)"; echo "python: $(python3 --version)"
echo; echo "# Verificación (1000 pasos): r y theta_0 deben coincidir"
$RS bench $G $L $DT 1000; $HOME/bench/bench_cpp $G $L $DT 1000; python3 bench.py $G $L $DT 1000
echo; echo "# Tiempos (µs por paso RK4; 5 repeticiones Rust/C++, 3 Python)"
for i in 1 2 3 4 5; do $RS bench $G $L $DT 20000; done
for i in 1 2 3 4 5; do $HOME/bench/bench_cpp $G $L $DT 20000; done
for i in 1 2 3; do python3 bench.py $G $L $DT 300; done
} > resultados.txt 2>&1
cat resultados.txt
