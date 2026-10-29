// Benchmark: Kuramoto con omega_i = k_i, RK4. Mismo algoritmo y mismo orden de operaciones que la
// versión en Rust (explosive_sync bench). Compilar: g++ -O3 -std=c++17 bench.cpp -o bench_cpp
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

struct Graph {
    std::vector<size_t> off;
    std::vector<unsigned> nb;
};

static void deriv(const Graph& g, const std::vector<double>& omega, double lam,
                  const std::vector<double>& th, std::vector<double>& s, std::vector<double>& c,
                  std::vector<double>& out) {
    const size_t n = th.size();
    for (size_t i = 0; i < n; i++) {
#ifdef USE_SINCOS
        sincos(th[i], &s[i], &c[i]);
#else
        s[i] = std::sin(th[i]);
        c[i] = std::cos(th[i]);
#endif
    }
    for (size_t i = 0; i < n; i++) {
        double ss = 0.0, cc = 0.0;
        for (size_t k = g.off[i]; k < g.off[i + 1]; k++) {
            ss += s[g.nb[k]];
            cc += c[g.nb[k]];
        }
        out[i] = omega[i] + lam * (ss * c[i] - cc * s[i]);
    }
}

int main(int argc, char** argv) {
    if (argc != 5) {
        std::fprintf(stderr, "uso: bench_cpp <enlaces.txt> <lambda> <dt> <pasos>\n");
        return 2;
    }
    FILE* f = std::fopen(argv[1], "r");
    size_t n, e;
    if (!f || std::fscanf(f, "%zu %zu", &n, &e) != 2) return 1;
    std::vector<std::pair<unsigned, unsigned>> edges(e);
    std::vector<size_t> deg(n, 0);
    for (size_t k = 0; k < e; k++) {
        if (std::fscanf(f, "%u %u", &edges[k].first, &edges[k].second) != 2) return 1;
        deg[edges[k].first]++;
        deg[edges[k].second]++;
    }
    std::fclose(f);
    Graph g;
    g.off.assign(n + 1, 0);
    for (size_t i = 0; i < n; i++) g.off[i + 1] = g.off[i] + deg[i];
    g.nb.assign(g.off[n], 0);
    std::vector<size_t> fill(g.off.begin(), g.off.end() - 1);
    for (auto& ed : edges) {
        g.nb[fill[ed.first]++] = ed.second;
        g.nb[fill[ed.second]++] = ed.first;
    }
    const double lam = std::atof(argv[2]), dt = std::atof(argv[3]);
    const long steps = std::atol(argv[4]);
    const double TAU = 6.283185307179586;
    std::vector<double> omega(n), th(n), k1(n), k2(n), k3(n), k4(n), tmp(n), s(n), c(n);
    for (size_t i = 0; i < n; i++) {
        omega[i] = (double)deg[i];
        double x = (double)i * 0.6180339887498949;
        th[i] = TAU * (x - std::floor(x));
    }
    auto t0 = std::chrono::steady_clock::now();
    for (long st = 0; st < steps; st++) {
        deriv(g, omega, lam, th, s, c, k1);
        for (size_t i = 0; i < n; i++) tmp[i] = th[i] + 0.5 * dt * k1[i];
        deriv(g, omega, lam, tmp, s, c, k2);
        for (size_t i = 0; i < n; i++) tmp[i] = th[i] + 0.5 * dt * k2[i];
        deriv(g, omega, lam, tmp, s, c, k3);
        for (size_t i = 0; i < n; i++) tmp[i] = th[i] + dt * k3[i];
        deriv(g, omega, lam, tmp, s, c, k4);
        for (size_t i = 0; i < n; i++) th[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
    double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    double cs = 0, sn = 0;
    for (size_t i = 0; i < n; i++) {
        cs += std::cos(th[i]);
        sn += std::sin(th[i]);
    }
    std::printf("cpp r=%.12f th0=%.12f us_por_paso=%.3f\n", std::sqrt(cs * cs + sn * sn) / (double)n, th[0],
                sec / (double)steps * 1e6);
    return 0;
}
