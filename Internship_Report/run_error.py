#!/usr/bin/env python3
"""
run_error.py
------------
Single run of the BFECC solver on the simple-wave IC and comparison
against the exact solution.

Euler (ck=0) and BFECC (ck=1) are shown side by side.
"""

import os, sys, importlib
import numpy as np
import matplotlib.pyplot as plt

current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(current_dir)

from exact_simple_wave import (
    make_simple_wave_ic, SimpleWaveExact, _sound_speed,
)

_isl      = importlib.import_module("1D_Euler_RBFGA_ImplicitSL")
LPTSolver = _isl.LPTSolver

# ── Parameters ────────────────────────────────────────────────────────────────

K, gamma   = 1.0, 1.4
rho_bg     = 1.2
xmin, xmax = -6.0, 6.0
Nx         = 256
CFL        = 0.5
r          = 2
T_max      = 0.5

def perturbation(x):
    return -0.5 * np.exp(-x**2 / (2 * 0.8**2))

# ── Build IC ──────────────────────────────────────────────────────────────────

dx  = (xmax - xmin) / Nx
x   = np.linspace(xmin + dx / 2, xmax - dx / 2, Nx)

rho0, u0, _ = make_simple_wave_ic(x, rho_bg, K, gamma, perturbation)
rho0_fn     = lambda xi: np.interp(xi, x, rho0)
u0_fn       = lambda xi: np.interp(xi, x, u0)

exact  = SimpleWaveExact(x, rho0, u0, K=K, gamma=gamma,
                         rho_bg=rho_bg, perturbation=perturbation)
t_star = exact.shock_time()
print(f"Shock formation time: t* = {t_star:.4f}  (running to T = {T_max})")
if T_max >= t_star:
    raise ValueError(f"T_max={T_max} >= t*={t_star:.4f}. Reduce T_max.")

# ── Run both schemes ──────────────────────────────────────────────────────────

cases = [(0, 'Euler (ck=0)', 'C0'), (1, 'BFECC (ck=1)', 'C1')]
results = {}

for ck_val, label, col in cases:
    sol = LPTSolver(x, rho0_fn, u0_fn, K=K, gamma=gamma, r=r, ck=ck_val)
    t   = 0.0
    while t < T_max:
        lam  = np.abs(sol.u) + _sound_speed(sol.rho, K, gamma)
        dt   = min(CFL * dx / float(np.max(lam)), T_max - t)
        sol.step(dt)
        t   += dt
    err_rho, err_u = exact.l2_error(x, sol.rho, sol.u, T_max)
    results[label] = (err_rho, err_u, sol, col)
    print(f"[{label}]  err_rho = {err_rho:.4e}   err_u = {err_u:.4e}")

# ── Plot ──────────────────────────────────────────────────────────────────────

rho_ex, u_ex = exact.solve(x, T_max)

fig, axes = plt.subplots(1, 2, figsize=(13, 5))
fig.suptitle(f"BFECC vs exact simple wave  (Nx={Nx}, CFL={CFL}, r={r}, T={T_max})")

for label, (_, _, sol, col) in results.items():
    axes[0].plot(x, sol.rho, marker='.', ms=3, label=label, color=col)
    axes[1].plot(x, sol.u,   marker='.', ms=3, label=label, color=col)

for ax, exact_arr, title in zip(axes, [rho_ex, u_ex], ['ρ', 'u']):
    ax.plot(x, exact_arr, 'k--', lw=1.5, label='exact')
    ax.set_title(title)
    ax.set_xlabel('x')
    ax.legend()
    ax.grid(True)

plt.tight_layout()
plt.show()
