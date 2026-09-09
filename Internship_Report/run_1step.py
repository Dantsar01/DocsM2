#!/usr/bin/env python3
"""
run_1step.py
------------
Single-step per-step error diagnostic.

One BFECC step from the exact simple-wave IC; error compared to the exact
solution at t = dt.  dt scales as CFL * dx / lam_max so the per-step error
should scale as O(dx^{2r+1}) by the RBF-GA interpolation theory.

This test isolates the spatial interpolation order independently of any
accumulation or BFECC correction.
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

K, gamma   = 1.0, 1.1
rho_bg     = 1.2
xmin, xmax = -6.0, 6.0
CFL_diag   = 1.
ck         = 1

Nx_list = [32, 64, 128, 256, 512, 1024, 2048]

def perturbation(x):
    return -0.5 * np.exp(-x**2 / (2 * 0.8**2))

# (r_val, label)
cases = [
    (1, 'r=1  (M=3,  theory O(h^3))'),
    (2, 'r=2  (M=5,  theory O(h^5))'),
    (3, 'r=3  (M=7,  theory O(h^7))'),
    (4, 'r=4  (M=9,  theory O(h^9))'),
]

# ── Diagnostic loop ───────────────────────────────────────────────────────────

diag_results = {}

for r_val, label in cases:
    theory_p = 2 * r_val + 1
    print(f"\n[{label}]  — expected O(h^{theory_p})")
    errs, dxs = [], []

    for Nx in Nx_list:
        dx   = (xmax - xmin) / Nx
        x    = np.linspace(xmin + dx / 2, xmax - dx / 2, Nx)

        rho0, u0, _ = make_simple_wave_ic(x, rho_bg, K, gamma, perturbation)
        rho0_fn = lambda xi, _r=rho0, _x=x: np.interp(xi, _x, _r)
        u0_fn   = lambda xi, _u=u0,   _x=x: np.interp(xi, _x, _u)

        sol = LPTSolver(x, rho0_fn, u0_fn, K=K, gamma=gamma, r=r_val, ck=ck)
        lam = np.abs(sol.u) + _sound_speed(sol.rho, K, gamma)
        dt  = CFL_diag * min(dx **((2 * r_val + 1)/(ck + 2)), dx / float(np.max(lam)))

        sol.step(dt)

        exact = SimpleWaveExact(x, rho0, u0, K=K, gamma=gamma,
                                rho_bg=rho_bg, perturbation=perturbation)
        err_rho, _ = exact.l2_error(x, sol.rho, sol.u, dt)
        errs.append(err_rho)
        dxs.append(dx)

        if len(errs) > 1:
            p_obs = np.log(errs[-2] / err_rho) / np.log(dxs[-2] / dx)
            print(f"  Nx={Nx:5d}  dx={dx:.5f}  dt={dt:.2e}  err={err_rho:.3e}  p={p_obs:.2f}")
        else:
            print(f"  Nx={Nx:5d}  dx={dx:.5f}  dt={dt:.2e}  err={err_rho:.3e}  (  -  )")

    diag_results[label] = (np.array(dxs), np.array(errs))

# ── Plot ──────────────────────────────────────────────────────────────────────

fig, ax = plt.subplots(figsize=(8, 6))
ax.set_title("Single-step error vs dx  (CFL = {:.2f})".format(CFL_diag))
colors = ['C0', 'C1', 'C2', 'C3']

for (r_val, label), col in zip(cases, colors):
    dxs_d, errs_d = diag_results[label]
    ax.loglog(dxs_d, errs_d, 'o-', label=label, color=col)

dx_ref = np.array([(xmax - xmin) / Nx for Nx in Nx_list])
err0   = diag_results[cases[0][1]][1][0]
for p_ref, ls in [(3, '--'), (5, ':'), (7, '-.'), (9, (0, (3, 1, 1, 1)))]:
    ax.loglog(dx_ref, err0 * (dx_ref / dx_ref[0]) ** p_ref,
              linestyle=ls, color='gray', alpha=0.5, label=f'O(h^{p_ref})')

ax.set_xlabel('dx')
ax.set_ylabel('L2 error (ρ) — one step only')
ax.legend(fontsize=9)
ax.grid(True, which='both', alpha=0.3)
plt.tight_layout()
plt.show()
