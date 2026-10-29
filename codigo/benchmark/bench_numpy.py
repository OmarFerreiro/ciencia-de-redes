"""Benchmark: Kuramoto con omega_i = k_i, RK4, con NumPy y matriz dispersa SciPy (CSR).

Mismo problema que bench.cpp / bench.py / `explosive_sync bench`. La suma sobre vecinos se hace con
un producto matriz dispersa por vector (el orden de las sumas difiere del de Rust y C++, así que los
resultados pueden diferir en las últimas cifras).
Uso: python bench_numpy.py <enlaces.txt> <lambda> <dt> <pasos>
"""
import sys
import time

import numpy as np
import scipy.sparse as sp


def main() -> None:
    path, lam, dt, steps = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
    with open(path) as f:
        n, e = map(int, f.readline().split())
        ed = np.loadtxt(f, dtype=np.int64, max_rows=e)
    rows = np.concatenate([ed[:, 0], ed[:, 1]])
    cols = np.concatenate([ed[:, 1], ed[:, 0]])
    a = sp.csr_matrix((np.ones(rows.size), (rows, cols)), shape=(n, n))
    omega = np.asarray(a.sum(axis=1)).ravel()
    th = 2.0 * np.pi * ((np.arange(n) * 0.6180339887498949) % 1.0)

    def deriv(t):
        s, c = np.sin(t), np.cos(t)
        return omega + lam * ((a @ s) * c - (a @ c) * s)

    t0 = time.perf_counter()
    for _ in range(steps):
        k1 = deriv(th)
        k2 = deriv(th + 0.5 * dt * k1)
        k3 = deriv(th + 0.5 * dt * k2)
        k4 = deriv(th + dt * k3)
        th = th + dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    sec = time.perf_counter() - t0
    r = abs(np.exp(1j * th).sum()) / n
    print(f"numpy r={r:.12f} th0={th[0]:.12f} us_por_paso={sec / steps * 1e6:.3f}")


main()
