#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_burgers_convergence.py
--------------------------
Spatial convergence study for BurgersSolver (exact GP, 1D_Burgers_BFECC.py).

IC : smooth Gaussian bump  u_0(x) = A * exp(-x²/σ²)
     (shock forms at T_shock ≈ σ√e / (A√2) ≈ 2.33 for A=1, σ=2)

Reference : BurgersExact — method of characteristics, machine-precision.

dt scaling : dt = C_cfl · dx^{p/(q+1)}  to balance temporal O(dt^q) and
             spatial O(dx^p) errors, where p = 2r+1 and q = ck+1.
             Resulting combined order : q·p/(q+1).

Two groups of cases:
  I.  Euler SL  (ck=0, q=1) : expected O(h^{p/2})
  II. Midpoint  (ck=1, q=2) : expected O(h^{2p/3})
"""

import os, sys, importlib, time
import numpy as np
import matplotlib.pyplot as plt

current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(current_dir)


# ── Load solver and exact reference from 1D_Burgers_BFECC ────────────────────

def _load(name):
    if name in sys.modules:
        importlib.reload(sys.modules[name])
    return importlib.import_module(name)

_mod         = _load("1D_Burgers_BFECC")
BurgersSolver = _mod.BurgersSolver
build_burgers_exact = _mod.build_burgers_exact


# ── Parameters ───────────────────────────────────────────────────────────────

A, sigma = 1.0, 2.0
xmin, xmax = -10.0, 10.0

T_shock = sigma * np.sqrt(np.e) / (A * np.sqrt(2.0))   # ≈ 2.33

T_max   = 1. #* T_shock   # well before shock; ≈ 1.28

eps          = 1e-3    # RBF-GA shape parameter
newton_iter  = 1       # Newton iterations per departure solve
weight_mode  = 'poly'  # 'poly' | 'exact' | 'table'
C_cfl        = 0.9     # prefactor in  dt = C · dx^{p/(q+1)}

Nx_list = [50, 100, 200, 300, 400]

def u_0(x):
    return A * np.exp(-x**2 / sigma**2)


# ── Helpers ───────────────────────────────────────────────────────────────────

_markers = ['o', 's', '^', 'D', 'v', 'p']

def _orders(errs, dxs):
    return ([np.log(errs[i] / errs[i+1]) / np.log(dxs[i] / dxs[i+1])
             for i in range(len(errs) - 1)] + [np.nan])

def _ref_slopes(ax, dx_arr, err_arr, slopes):
    styles = [('--','silver'), (':','gray'), ('-.','darkgray'),
              (':','black'), ('--','dimgray')]
    e0, d0 = float(err_arr[0]), float(dx_arr[0])
    for (p_ref, (ls, col)) in zip(slopes, styles):
        ax.loglog(dx_arr, e0 * (dx_arr / d0)**p_ref,
                  ls, color=col, alpha=0.5, label=f'O(h^{p_ref:.2f})')


# ── Build exact reference once (extended IC so all departure points stay inside)

print(f"T_shock ≈ {T_shock:.4f},  running to T = {T_max:.4f}")
print("Building exact characteristics reference ...")
exact_ref = build_burgers_exact(u_0, xmin, xmax, T_max, N_ic=32768)
print("  done.\n")


# ════════════════════════════════════════════════════════════════════════════════
# Convergence loop
# ════════════════════════════════════════════════════════════════════════════════

cases = [
    # (ck,  r,  label)
    (0, 1, 'Euler  r=1'),
    (0, 2, 'Euler  r=2'),
    (0, 3, 'Euler  r=3'),
    #(1, 1, 'trapz  r=1'),
    #(1, 2, 'trapz  r=2'), 
    #(1, 3, 'BFECC  r=3'),
]

conv = {}

for ck_val, r_val, label in cases:
    q        = ck_val + 1           # temporal order
    p        = 2 * r_val + 1        # spatial order (stencil exactness)
    exp_ord  = q * p / (q + 1)

    print("=" * 64)
    print(f"[{label}]  ck={ck_val}  r={r_val}  p={p}  q={q}  "
          f"dt∝dx^{p/(q+1):.3f}   expected O(h^{exp_ord:.3f})")
    print("=" * 64)

    errs, dxs, T_reached = [], [], []

    for Nx in Nx_list:
        dx    = (xmax - xmin) / Nx
        x_ref = np.linspace(xmin + dx / 2, xmax - dx / 2, Nx)

        sol = BurgersSolver(
            x_ref, u_0,
            r=r_val, eps=eps, ck=ck_val, newton_iter=newton_iter,
            weight_mode=weight_mode,
        )

        dt_ref = C_cfl * dx ** (p / (q + 1))

        t = 0.0
        while t < T_max:
            # CFL-safe cap: u_max * dt / dx ≤ 1
            u_max  = float(np.max(np.abs(sol.u))) + 1e-14
            dt_cfl = dx / u_max
            dt     = min(dt_ref, dt_cfl, T_max - t)
            sol.step(dt)
            t += dt

        err = exact_ref.l2_error(x_ref, sol.u, t)
        errs.append(err); dxs.append(dx); T_reached.append(t)

        if len(errs) > 1:
            p_emp = np.log(errs[-2] / err) / np.log(dxs[-2] / dx)
            print(f"  Nx={Nx:4d}  dx={dx:.5f}  dt_ref={dt_ref:.3e}  "
                  f"err={err:.4e}  order={p_emp:.2f}")
        else:
            print(f"  Nx={Nx:4d}  dx={dx:.5f}  dt_ref={dt_ref:.3e}  "
                  f"err={err:.4e}")

    conv[label] = {
        'ck'      : ck_val,
        'r'       : r_val,
        'Nx'      : np.array(Nx_list),
        'dx'      : np.array(dxs),
        'err'     : np.array(errs),
        'T'       : np.array(T_reached),
        'orders'  : np.array(_orders(errs, dxs)),
        'exp_ord' : exp_ord,
    }

# ── Summary table ────────────────────────────────────────────────────────────

print("\n" + "="*64)
print("Summary of empirical convergence orders")
print("="*64)
for ck_val, r_val, label in cases:
    d = conv[label]
    ords = np.round(d['orders'][:-1], 2)
    print(f"  [{label}]  expected O(h^{d['exp_ord']:.3f})  orders: {ords}")


# ── Plot ─────────────────────────────────────────────────────────────────────

colors  = {cl: f'C{i}' for i, (*_, cl) in enumerate(cases)}
markers = {cl: _markers[i % len(_markers)] for i, (*_, cl) in enumerate(cases)}

fig, axes = plt.subplots(1, 2, figsize=(14, 5))
fig.suptitle(
    f"Spatial convergence — Burgers,  u₀ = exp(−x²/σ²),  T={T_max:.3f}\n"
    f"weight_mode='{weight_mode}',  newton_iter={newton_iter},  "
    f"dt = {C_cfl}·dx^{{p/(q+1)}}",
    fontsize=9,
)

# ── Left panel: all cases ─────────────────────────────────────────────────────
ax = axes[0]
for *_, label in cases:
    d = conv[label]
    ax.loglog(d['dx'], d['err'],
              marker=markers[label], linestyle='-',
              color=colors[label], label=label)

# Reference slopes spanning the full dx range
all_dx  = np.concatenate([d['dx'] for d in conv.values()])
dx_min  = all_dx.min();  dx_max = all_dx.max()
dx_line = np.array([dx_max, dx_min])
# anchor: first case, largest dx
first_label = cases[0][-1]
err0 = conv[first_label]['err'][0]
dx0  = conv[first_label]['dx'][0]
for slope, ls, col in [(1.5,'--','silver'), (2.0,':','gray'),
                       (10/3,'-.','darkgray'), (14/3,'--','black')]:
    ax.loglog(dx_line, err0 * (dx_line / dx0)**slope,
              ls, color=col, alpha=0.5, label=f'O(h^{slope:.2f})')

ax.set_xlabel('dx');  ax.set_ylabel('L² error (u)')
ax.set_title('All cases');  ax.legend(fontsize=7);  ax.grid(True, which='both', alpha=0.3)

# ── Right panel: highlight ck=1 cases and show pre-asymptotic regime ──────────
ax = axes[1]
for *_, label in cases:
    d = conv[label]
    ls = '-' if d['ck'] == 1 else '--'
    ax.loglog(d['dx'], d['err'],
              marker=markers[label], linestyle=ls,
              color=colors[label], label=f"{label} (exp {d['exp_ord']:.2f})")

# Reference slopes
for slope, ls, col in [(2.0,'--','silver'), (10/3,':','gray'), (14/3,'-.','darkgray')]:
    ax.loglog(dx_line, err0 * (dx_line / dx0)**slope,
              ls, color=col, alpha=0.5, label=f'O(h^{slope:.2f})')

ax.set_xlabel('dx');  ax.set_ylabel('L² error (u)')
ax.set_title('Midpoint (ck=1) vs Euler (ck=0, dashed)')
ax.legend(fontsize=7);  ax.grid(True, which='both', alpha=0.3)

plt.tight_layout()
'''
plt.savefig(os.path.join(current_dir, 'burgers_convergence.png'),
            dpi=150, bbox_inches='tight')

print("\nSaved burgers_convergence.png")
'''
plt.show()
