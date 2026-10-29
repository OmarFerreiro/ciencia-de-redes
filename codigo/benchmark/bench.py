"""Benchmark: Kuramoto con omega_i = k_i, RK4, en Python puro (sin NumPy).

Mismo algoritmo y mismo orden de operaciones que las versiones de Rust y C++.
Uso: python3 bench.py <enlaces.txt> <lambda> <dt> <pasos>
"""
import math
import sys
import time


def main() -> None:
    path, lam, dt, steps = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
    with open(path) as f:
        n, e = map(int, f.readline().split())
        nb = [[] for _ in range(n)]
        for _ in range(e):
            a, b = map(int, f.readline().split())
            nb[a].append(b)
            nb[b].append(a)
    omega = [float(len(v)) for v in nb]
    tau = 6.283185307179586
    th = [tau * ((i * 0.6180339887498949) % 1.0) for i in range(n)]
    sin, cos = math.sin, math.cos
    rng = range(n)

    def deriv(t):
        s = [sin(x) for x in t]
        c = [cos(x) for x in t]
        out = [0.0] * n
        for i in rng:
            ss = 0.0
            cc = 0.0
            for j in nb[i]:
                ss += s[j]
                cc += c[j]
            out[i] = omega[i] + lam * (ss * c[i] - cc * s[i])
        return out

    t0 = time.perf_counter()
    for _ in range(steps):
        k1 = deriv(th)
        k2 = deriv([th[i] + 0.5 * dt * k1[i] for i in rng])
        k3 = deriv([th[i] + 0.5 * dt * k2[i] for i in rng])
        k4 = deriv([th[i] + dt * k3[i] for i in rng])
        th = [th[i] + dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) for i in rng]
    sec = time.perf_counter() - t0
    cs = sum(cos(x) for x in th)
    sn = sum(sin(x) for x in th)
    print(f"python r={math.hypot(cs, sn) / n:.12f} th0={th[0]:.12f} us_por_paso={sec / steps * 1e6:.3f}")


main()
