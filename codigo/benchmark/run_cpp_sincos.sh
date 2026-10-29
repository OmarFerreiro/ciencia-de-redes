cd /mnt/c/Users/jorge/proyectos/ciencia-de-redes/codigo/benchmark
g++ -O3 -std=c++17 -DUSE_SINCOS bench.cpp -o $HOME/bench/bench_cpp_sincos
G=$HOME/bench/red.txt
echo "# C++ con sincos (una llamada para seno y coseno, como hace Rust con sin_cos)"
$HOME/bench/bench_cpp_sincos $G 1.5 0.002 1000
for i in 1 2 3 4 5; do $HOME/bench/bench_cpp_sincos $G 1.5 0.002 20000; done
echo "# Rust y C++ originales, de nuevo, intercalados"
for i in 1 2 3; do $HOME/rt/es/release/explosive_sync bench $G 1.5 0.002 20000; $HOME/bench/bench_cpp $G 1.5 0.002 20000; done
