#!/usr/bin/env python3
"""
run_temporal_convergence.py
---------------------------
Temporal convergence study: fix Nx, vary CFL (→ dt).

Toggle wave_dir to select the IC and reference:
  'right' / 'left'  →  simple-wave IC with exact analytical reference.
                        Clean temporal order is visible down to the
                        spatial floor C_s * dx^p / dt.
  None              →  two-wave IC with a pseudo-spectral (DOP853) reference.

For each tested CFL the total error is fitted to
    E(dt) = C_t * dt^q + C_sp / dt
to locate the crossover CFL where reducing dt starts increasing the error.
"""

import os, sys, importlib
import numpy as np
import matplotlib.pyplot as plt

current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(current_dir)

from exact_simple_wave import (
    make_simple_wave_ic, SimpleWaveExact,
    temporal_convergence_study, _sound_speed,
)
from exact_two_wave import build_two_wave_reference

_isl      = importlib.import_module("1D_Euler_RBFGA_ImplicitSL")
LPTSolver = _isl.LPTSolver

# ── Parameters ────────────────────────────────────────────────────────────────

K, gamma   = 1.0, 3.
rho_bg     = 1.2
xmin, xmax = -6.0, 6.0
kernel     = 'se'
r          = 2

Nx_t  = 512
T_t   = 1.

def perturbation(x):
    return -0.5 * np.exp(-x**2 / (2 * 0.8**2))

# 'right', 'left', or None (two-wave + fine BFECC reference)
wave_dir = None

# ── Build IC and reference ────────────────────────────────────────────────────

dx_t = (xmax - xmin) / Nx_t
x_t  = np.linspace(xmin + dx_t / 2, xmax - dx_t / 2, Nx_t)

if wave_dir in ('right', 'left'):
    rho0_t, u0_t, _ = make_simple_wave_ic(x_t, rho_bg, K, gamma,
                                           perturbation, direction=wave_dir)
    exact_t = SimpleWaveExact(x_t, rho0_t, u0_t, K=K, gamma=gamma,
                              rho_bg=rho_bg, perturbation=perturbation,
                              direction=wave_dir)
    t_star = exact_t.shock_time()
    print(f"\n── Temporal study ({wave_dir}-going IC, exact ref)  Nx={Nx_t} ──")
    print(f"   Shock time t* = {t_star:.4f}  (T = {T_t})")

    def rho_0_t(xi):
        r0, _, _ = make_simple_wave_ic(xi, rho_bg, K, gamma,
                                        perturbation, direction=wave_dir)
        return r0

    def u_0_t(xi):
        _, u0, _ = make_simple_wave_ic(xi, rho_bg, K, gamma,
                                        perturbation, direction=wave_dir)
        return u0

    def exact_fn_t(xi, rho, u, t):
        rho_ex, u_ex = exact_t.solve(xi, t)
        dxi = xi[1] - xi[0]
        return (float(np.sqrt(np.sum((rho - rho_ex)**2) * dxi)),
                float(np.sqrt(np.sum((u   - u_ex  )**2) * dxi)))

    study_title = f"Temporal convergence — {wave_dir}-going IC, exact ref  (Nx={Nx_t})"

else:
    def rho_0_t(xi): return 1.5 + np.exp(-xi**2 / 4)
    def u_0_t(xi):   return -0.2 * np.exp(-xi**2 / 4)

    print(f"\n── Temporal study (two-wave IC, spectral ref) ──")
    print(f"   Building spectral reference ...")
    _ref_tw = build_two_wave_reference(rho_0_t, u_0_t,
                  K=K, gamma=gamma, xmin=xmin, xmax=xmax, T_target=T_t)
    print("   Reference done.\n")

    def exact_fn_t(xi, rho, u, _t):
        return _ref_tw.l2_error(xi, rho, u)

    study_title = f"Temporal convergence — two-wave IC, spectral ref  (Nx={Nx_t})"

# ── CFL / dt list ─────────────────────────────────────────────────────────────

lam_max_t  = float(np.max(np.abs(u_0_t(x_t)) + _sound_speed(rho_0_t(x_t), K, gamma)))
CFL_list   = np.array([5, 4, 3, 2, 1], dtype=float)
dt_list_t  = [cfl * dx_t / lam_max_t for cfl in CFL_list]

print(f"   dx={dx_t:.5f}  lam_max={lam_max_t:.4f}")
print(f"   CFL_list = {np.round(CFL_list, 3)}")
print(f"   dt_list  = {[f'{d:.5f}' for d in dt_list_t]}\n")

# ── Run ───────────────────────────────────────────────────────────────────────

ck_cases = [(0, 'Euler', 'C0'), (1, 'BFECC', 'C1')]
temp_results = {}

for ck_val, ck_label, _ in ck_cases:
    print(f"[{ck_label}]")
    temp_results[ck_label] = temporal_convergence_study(
        xmin, xmax, Nx_t, T_t,
        rho_0_t, u_0_t, K, gamma,
        dt_list=dt_list_t,
        solver_cls=LPTSolver,
        solver_kwargs={'kernel': kernel, 'r': r, 'ck': ck_val},
        exact_fn=exact_fn_t,
    )
    res = temp_results[ck_label]
    print(f"  orders (rho): {np.round(res['order_rho'], 2)}")
    print(f"  orders (u)  : {np.round(res['order_u'],   2)}\n")

# ── Crossover CFL fit:  E(dt) = C_t * dt^q + C_sp / dt ──────────────────────
# Fit is linear in [dt^q, 1/dt] → [C_t, C_sp] via least squares.
# Crossover: C_t * dt_cross^q = C_sp / dt_cross  → dt_cross = (C_sp/C_t)^{1/(q+1)}

floor_info = {}
for ck_val, ck_label, _ in ck_cases:
    q     = ck_val + 1    # global order: Euler→1, BFECC→2
    dts   = np.array(temp_results[ck_label]['dt'])
    info  = {}
    for key in ['err_rho', 'err_u']:
        errs  = np.array(temp_results[ck_label][key])
        A_mat = np.column_stack([dts**q, 1.0 / dts])
        coeffs, *_ = np.linalg.lstsq(A_mat, errs, rcond=None)
        C_t, C_sp  = float(coeffs[0]), float(coeffs[1])
        if C_t > 0 and C_sp > 0:
            dt_cross  = (C_sp / C_t) ** (1.0 / (q + 1))
            cfl_cross = dt_cross * lam_max_t / dx_t
            E_cross   = C_t * dt_cross**q + C_sp / dt_cross
        else:
            dt_cross = cfl_cross = E_cross = float('nan')
        info[key] = (C_t, C_sp, dt_cross, cfl_cross, E_cross)
        print(f"  [{key}] {ck_label} (q={q}): C_t={C_t:.3e}  C_sp={C_sp:.3e}  "
              f"CFL_cross={cfl_cross:.3f}")
    floor_info[ck_label] = info

# ── Plot ──────────────────────────────────────────────────────────────────────

fig, axes = plt.subplots(1, 2, figsize=(13, 5))
fig.suptitle(study_title)

for ax, key, ylbl in zip(axes, ['err_rho', 'err_u'], ['ρ', 'u']):
    for ck_val, ck_label, col in ck_cases:
        r_res = temp_results[ck_label]
        ax.loglog(r_res['dt'], r_res[key], 'o-', color=col, label=ck_label)

    for (ck_val, ck_label, col), (p, ls) in zip(ck_cases, [(1,'--'),(2,':')]):
        r_res = temp_results[ck_label]
        dts   = np.array(r_res['dt'])
        e0    = r_res[key][0]
        ax.loglog(dts, e0 * (dts / dts[0])**p,
                  linestyle=ls, color=col, alpha=0.35, label=f'O(dt^{p})')
    '''
    for ck_val, ck_label, col in ck_cases:
        C_t, C_sp, dt_cross, cfl_cross, E_cross = floor_info[ck_label][key]
        if not np.isnan(dt_cross):
            ax.axvline(dt_cross, color=col, linestyle=':', linewidth=1.2, alpha=0.7,
                       label=f'CFL_cross≈{cfl_cross:.2f} ({ck_label})')
            ax.axhline(E_cross,  color=col, linestyle='--', linewidth=0.8, alpha=0.4)
    '''
    ax.set_xlabel('dt')
    ax.set_ylabel(f'L2 error ({ylbl})')
    ax.set_title(ylbl)
    ax.legend(fontsize=8)
    ax.grid(True, which='both', alpha=0.3)

plt.tight_layout()
plt.show()
