#!/usr/bin/env python3
"""
benchmark_interp_weights.py
---------------------------
Compares three strategies for obtaining RBF-GA weights inside _interp_table:

  (A) Table  — precomputed weight table, linear or cubic lookup   O(N·M)  runtime
  (B) Exact  — weights_batch() solves the M×M system per step     O(N·M³) runtime
  (C) Poly   — precomputed polynomial fit of w(δ), Horner eval    O(N·M²) runtime
               (analytical, zero table-interpolation error)

Complexity summary
------------------
                Build                  Per-call          Table error
  Table deg=1   O(n_tab · M³)         O(N · M)          O(Δδ²)
  Table deg=3   O(n_tab · M³)         O(N · M)          O(Δδ⁴)
  Exact         —                     O(N · M³)         0
  Poly          O(K · M³ + M · P²)   O(N · M · P)      ~machine ε
                  K = oversampled nodes (~4M)
                  P = poly degree (~2r+2)

For r=1,2,3: M = 3,5,7 and P ≈ 4,6,8.  Poly ≈ Exact in cost at small N
but much faster at large N since M·P << M³.
"""

import os, sys, time
import numpy as np

current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(current_dir)

import importlib
_mod = importlib.import_module("1D_Euler_RBFGA_ImplicitSL")
LPTSolver  = _mod.LPTSolver
from rbf_ga_weights_1d import RBFGAStencil

# ── Polynomial weight representation ─────────────────────────────────────────

def build_poly_weights(dx, r, eps=1e-3, n_fit=None, poly_deg=None):
    """
    Fit a polynomial p_k(δ) for each of the M weight components.

    Strategy
    --------
    1. Sample w(δ) at n_fit Chebyshev nodes on [-0.5, 0.5] (exact solves).
    2. Fit each component with a polynomial of degree poly_deg via least squares.

    Returns
    -------
    coeffs : (M, poly_deg+1) array of Horner coefficients (highest power first)
    """
    M = 2 * r + 1
    if poly_deg is None:
        poly_deg = 2 * r + 2          # one degree above the theoretical minimum
    if n_fit is None:
        n_fit = max(4 * (poly_deg + 1), 50)   # oversample for stable LS fit

    stencil     = RBFGAStencil(eps)
    x_nodes_loc = dx * np.arange(-r, r + 1, dtype=float)

    # Chebyshev nodes on [-0.5, 0.5]
    k_idx  = np.arange(1, n_fit + 1)
    delta_nodes = -0.5 * np.cos(np.pi * (2 * k_idx - 1) / (2 * n_fit))

    W_samples = np.zeros((n_fit, M))
    for i, d in enumerate(delta_nodes):
        W_samples[i] = stencil.weights(x_nodes_loc, -d * dx, deriv=0)

    # Least-squares polynomial fit for each weight component
    coeffs = np.zeros((M, poly_deg + 1))
    for m in range(M):
        coeffs[m] = np.polyfit(delta_nodes, W_samples[:, m], poly_deg)

    return coeffs


def poly_interp(f_src, x_dep, x_ref, dx, r, coeffs):
    """
    Evaluate weights via Horner's method on precomputed polynomial coefficients.
    O(N · M · poly_deg) per call.
    """
    N  = len(x_ref)
    x0 = x_ref[0]
    M  = 2 * r + 1

    frac_idx = (x_dep - x0) / dx
    j = np.round(frac_idx).astype(int)
    j = np.clip(j, r, N - r - 1)
    delta = j.astype(float) - frac_idx          # δ ∈ [−0.5, 0.5]

    # Horner evaluation for each weight component: W shape (N, M)
    W = np.zeros((len(delta), M))
    for m in range(M):
        W[:, m] = np.polyval(coeffs[m], delta)

    proto = np.arange(-r, r + 1, dtype=int)
    si    = np.clip(j[:, None] + proto[None, :], 0, N - 1)
    return (W * f_src[si]).sum(axis=1)


# ── Benchmark helpers ─────────────────────────────────────────────────────────

def _make_solver(Nx, r, xmin, xmax, exact):
    dx  = (xmax - xmin) / Nx
    x   = np.linspace(xmin + dx / 2, xmax - dx / 2, Nx)
    rho_fn = lambda xi: 1.5 + np.exp(-xi**2 / 4)
    u_fn   = lambda xi: -0.2 * np.exp(-xi**2 / 4)
    return LPTSolver(x, rho_fn, u_fn, r=r, ck=0,
                     n_table=5_000, table_deg=3, exact_weights=exact)


def _time(fn, repeats=5):
    fn()                               # warm-up
    t0 = time.perf_counter()
    for _ in range(repeats):
        fn()
    return (time.perf_counter() - t0) / repeats


# ── Main benchmark ────────────────────────────────────────────────────────────

def run_benchmark():
    xmin, xmax = -6.0, 6.0
    r_vals  = [1, 2, 3]
    Nx_vals = [128, 512, 2048]
    n_table = 5_000
    repeats = 10

    print("=" * 78)
    print(f"  Weight-table interpolation benchmark")
    print(f"  n_table={n_table}  repeats={repeats}")
    print("=" * 78)
    print(f"  {'Nx':>6}  {'r':>2}  {'M':>2}  "
          f"{'Build(tab)':>12}  {'Build(poly)':>12}  "
          f"{'Table-1':>10}  {'Table-3':>10}  "
          f"{'Exact':>10}  {'Poly':>10}")
    print("-" * 78)

    for Nx in Nx_vals:
        for r in r_vals:
            M   = 2 * r + 1
            dx  = (xmax - xmin) / Nx
            x   = np.linspace(xmin + dx / 2, xmax - dx / 2, Nx)
            rho_fn = lambda xi: 1.5 + np.exp(-xi**2 / 4)
            u_fn   = lambda xi: -0.2 * np.exp(-xi**2 / 4)

            # ── Build times ───────────────────────────────────────────────────
            t_build_tab = _time(
                lambda: LPTSolver(x, rho_fn, u_fn, r=r, ck=0,
                                  n_table=n_table, exact_weights=False),
                repeats=3)

            t_build_poly = _time(
                lambda: build_poly_weights(dx, r),
                repeats=3)

            # ── Solvers for call timing ────────────────────────────────────────
            sol_tab1  = LPTSolver(x, rho_fn, u_fn, r=r, ck=0,
                                  n_table=n_table, table_deg=1, exact_weights=False)
            sol_tab3  = LPTSolver(x, rho_fn, u_fn, r=r, ck=0,
                                  n_table=n_table, table_deg=3, exact_weights=False)
            sol_exact = LPTSolver(x, rho_fn, u_fn, r=r, ck=0,
                                  exact_weights=True)
            poly_coeffs = build_poly_weights(dx, r)

            # Fake departure points (CFL=0.05)
            lam_max = 1.4
            dt      = 0.05 * dx / lam_max
            x_dep   = x - lam_max * dt + 0.01 * dx * np.random.randn(Nx)
            f_src   = rho_fn(x)

            t_tab1  = _time(lambda: sol_tab1._interp_table(f_src, x_dep),  repeats)
            t_tab3  = _time(lambda: sol_tab3._interp_table(f_src, x_dep),  repeats)
            t_exact = _time(lambda: sol_exact._interp_table(f_src, x_dep), repeats)
            t_poly  = _time(
                lambda: poly_interp(f_src, x_dep, x, dx, r, poly_coeffs),
                repeats)

            def _ms(t): return f"{t*1e3:8.3f}ms"

            print(f"  {Nx:>6}  {r:>2}  {M:>2}  "
                  f"  {_ms(t_build_tab)}    {_ms(t_build_poly)}  "
                  f"{_ms(t_tab1)}  {_ms(t_tab3)}  "
                  f"{_ms(t_exact)}  {_ms(t_poly)}")

        print()

    # ── Accuracy check at a fixed point ──────────────────────────────────────
    print("=" * 78)
    print("  Accuracy check: |w_table − w_exact| and |w_poly − w_exact|")
    print(f"  (Nx=256, r=3, δ≈0.05, n_table={n_table})")
    print("-" * 78)

    Nx  = 256
    r   = 3
    M   = 2 * r + 1
    dx  = (xmax - xmin) / Nx
    x   = np.linspace(xmin + dx / 2, xmax - dx / 2, Nx)
    f_src = 1.5 + np.exp(-x**2 / 4)

    lam_max = 1.4
    dt      = 0.05 * dx / lam_max
    x_dep   = x - lam_max * dt       # uniform CFL departure

    sol_tab1  = LPTSolver(x, rho_fn, u_fn, r=r, ck=0,
                          n_table=n_table, table_deg=1, exact_weights=False)
    sol_tab3  = LPTSolver(x, rho_fn, u_fn, r=r, ck=0,
                          n_table=n_table, table_deg=3, exact_weights=False)
    sol_exact = LPTSolver(x, rho_fn, u_fn, r=r, ck=0, exact_weights=True)
    poly_coeffs = build_poly_weights(dx, r)

    f_tab1  = sol_tab1._interp_table(f_src, x_dep)
    f_tab3  = sol_tab3._interp_table(f_src, x_dep)
    f_exact = sol_exact._interp_table(f_src, x_dep)
    f_poly  = poly_interp(f_src, x_dep, x, dx, r, poly_coeffs)

    print(f"  max|tab1  − exact| = {np.max(np.abs(f_tab1  - f_exact)):.3e}")
    print(f"  max|tab3  − exact| = {np.max(np.abs(f_tab3  - f_exact)):.3e}")
    print(f"  max|poly  − exact| = {np.max(np.abs(f_poly  - f_exact)):.3e}")
    print("=" * 78)


if __name__ == "__main__":
    np.random.seed(0)
    run_benchmark()
