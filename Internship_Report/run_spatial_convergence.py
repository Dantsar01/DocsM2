#!/usr/bin/env python3
"""
run_spatial_convergence.py
--------------------------
Spatial convergence study for the BFECC solver.

Two sections:
  A) Simple-wave IC  — exact solution available.
     dt ∝ dx^(p/(q+1)) per case to balance temporal and spatial errors
     (Euler q=1: exponent p/2; BFECC q=2: exponent p/3).
     Expected order: q*p/(q+1).

  B) Two-wave IC     — no exact solution; uses a fine-grid BFECC run as
     reference.  Same dt scaling as section A.
"""

import os, sys, importlib
import numpy as np
import matplotlib.pyplot as plt

current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(current_dir)

from exact_simple_wave import (
    make_simple_wave_ic, SimpleWaveExact,
    _sound_speed,
)



def load(name):
    if name in sys.modules:
        importlib.reload(sys.modules[name])
    return importlib.import_module(name)

_isl      = load("1D_Euler_RBFGA_ImplicitSL")
LPTSolver = _isl.LPTSolver


# ── Shared parameters ─────────────────────────────────────────────────────────

K, gamma   = 1.0, 1.4
rho_bg     = 1.2
xmin, xmax = -6.0, 6.0
kernel     = 'se'
CFL_cap    = 0.9

weight_mode      = 'cheb'
cheb_pad         = 0
newton_iter      = 1

solver_kw = dict(weight_mode=weight_mode, cheb_pad=cheb_pad, newton_iter=newton_iter,
                 monotone_limiter=False)

def _fmtv(v):
    return f"'{v}'" if isinstance(v, str) else str(v)

param_str = "  ".join(f"{k}={_fmtv(v)}" for k, v in solver_kw.items())

def perturbation(x):
    return -0.5 * np.exp(-x**2 / (2 * 0.8**2))

_markers = ['o', 's', '^', 'D', 'v', 'p']

def _orders_from_list(errs, dxs):
    return [np.log(errs[i] / errs[i+1]) / np.log(dxs[i] / dxs[i+1])
            for i in range(len(errs) - 1)] + [np.nan]

# %%



# ════════════════════════════════════════════════════════════════════════════
# A) Simple-wave IC — exact reference, dt ∝ dx^(p/(q+1))
# ════════════════════════════════════════════════════════════════════════════

print("\n" + "="*60)
print("A) Simple-wave IC  (dt∝dx^(p/(q+1)), exact reference)")
print("="*60)

T_sw   = .5
Nx_sw  = 100 * np.arange(1,6)


cases_sw = [
    (0, 1, 'Euler  r=1'),
    (0, 2, 'Euler  r=2'),
    (1, 1, 'trapz  r=1'),
    (1, 2, 'trapz  r=2'),
    (1, 3, 'trapz  r=3'),
]

conv_sw = {}
for ck_val, r_val, label in cases_sw:
    q        = ck_val + 1          # Euler→1, BFECC→2
    p        = 2 * r_val + 1
    expected = q * p / (q + 1)
    print(f"\n[{label}]  p={p}  q={q}  dt∝dx^{p/(q+1):.3f}  expected O(h^{expected:.2f})")
    errs_rho, errs_u, T_list, dxs = [], [], [], []

    for Nx in Nx_sw:
        dx  = (xmax - xmin) / Nx
        x   = np.linspace(xmin + dx / 2, xmax - dx / 2, Nx)

        rho0, u0, _ = make_simple_wave_ic(x, rho_bg, K, gamma, perturbation)
        rho0_fn = lambda xi, _r=rho0, _x=x: np.interp(xi, _x, _r)
        u0_fn   = lambda xi, _u=u0,   _x=x: np.interp(xi, _x, _u)

        sol = LPTSolver(x, rho0_fn, u0_fn, K=K, gamma=gamma,
                        r=r_val, ck=ck_val, kernel=kernel, **solver_kw)

        dt_n  = dx ** (p / (q + 1))

        exact = SimpleWaveExact(x, rho0, u0, K=K, gamma=gamma,
                                rho_bg=rho_bg, perturbation=perturbation)

        t = 0.0
        while t < T_sw:
            sol.step(min(dt_n, T_sw - t))
            t += min(dt_n, T_sw - t)

        err_rho, err_u = exact.l2_error(x, sol.rho, sol.u, T_sw)
        errs_rho.append(err_rho); errs_u.append(err_u)
        T_list.append(t); dxs.append(dx)

        if len(errs_rho) > 1:
            p_rho = np.log(errs_rho[-2] / err_rho) / np.log(dxs[-2] / dx)
            p_u   = np.log(errs_u[-2]   / err_u)   / np.log(dxs[-2] / dx)
            print(f"  Nx={Nx:5d}  dx={dx:.5f}  dt={dt_n:.3e}  "
                  f"err_rho={err_rho:.3e} (p={p_rho:.2f})  err_u={err_u:.3e} (p={p_u:.2f})")
        else:
            print(f"  Nx={Nx:5d}  dx={dx:.5f}  dt={dt_n:.3e}  "
                  f"err_rho={err_rho:.3e}  err_u={err_u:.3e}")

    dxs_arr = np.array(dxs)
    conv_sw[label] = {
        'Nx'       : np.array(Nx_sw),
        'dx'       : dxs_arr,
        'T'        : np.array(T_list),
        'err_rho'  : np.array(errs_rho),
        'err_u'    : np.array(errs_u),
        'order_rho': np.array(_orders_from_list(errs_rho, dxs)),
        'order_u'  : np.array(_orders_from_list(errs_u,   dxs)),
    }

print("\n── Empirical orders (simple-wave, dt∝dx^(p/(q+1))) ──")
for ck_val, r_val, label in cases_sw:
    q        = ck_val + 1
    p        = 2 * r_val + 1
    expected = q * p / (q + 1)
    r = conv_sw[label]
    print(f"  [{label}]  expected O(h^{expected:.2f})  "
          f"rho: {np.round(r['order_rho'][:-1], 2)}  u: {np.round(r['order_u'][:-1], 2)}")

colors_sw  = {cl: f'C{i}' for i, (*_, cl) in enumerate(cases_sw)}
markers_sw = {cl: _markers[i % len(_markers)] for i, (*_, cl) in enumerate(cases_sw)}
anchor_sw  = cases_sw[0][-1]


fig, axes = plt.subplots(1, 2, figsize=(14, 5))
fig.suptitle(f"Spatial convergence — simple-wave IC, dt∝dx^(p/(q+1))  (T={T_sw})\n"
             f"{param_str}", fontsize=9)
for ax, key, ylbl in zip(axes, ['err_rho', 'err_u'], ['ρ', 'u']):
    for *_, label in cases_sw:
        r = conv_sw[label]
        ax.loglog(r['dx'], r[key], marker=markers_sw[label], linestyle='-',
                  color=colors_sw[label], label=label)
    dx_ref  = conv_sw[anchor_sw]['dx']
    err_ref = conv_sw[anchor_sw][key]
    for p_ref, ls, col in [(1.5,'--','silver'),(2.5,':','gray'),(3.5,'-.','darkgray'),
                           (10/3,':','black'),(14/3,'--','dimgray')]:
        ax.loglog(dx_ref, float(err_ref[0]) * (dx_ref / float(dx_ref[0]))**p_ref,
                  linestyle=ls, color=col, alpha=0.45, label=f'O(h^{p_ref:.2f})')
    ax.set_xlabel('dx'); ax.set_ylabel(f'L2 error ({ylbl})')
    ax.set_title(ylbl); ax.legend(fontsize=8); ax.grid(True, which='both', alpha=0.3)
plt.tight_layout()
plt.show()

# %%

weight_mode      = 'cheb'
cheb_pad         = 0
newton_iter      = 1
interp_backend   = 'gp'
#adaptive_newton  = False

solver_kw = dict(weight_mode=weight_mode, cheb_pad=cheb_pad, newton_iter=newton_iter)
param_str = "  ".join(f"{k}={_fmtv(v)}" for k, v in solver_kw.items())

def load(name):
    if name in sys.modules:
        importlib.reload(sys.modules[name])
    return importlib.import_module(name)

__isl      = load("1D_Euler_RBFGA_ImplicitSL")
LPTSolver = __isl.LPTSolver

# ════════════════════════════════════════════════════════════════════════════
# B) Two-wave IC — fine BFECC/CK reference, dt = CFL·dx/lam
# ════════════════════════════════════════════════════════════════════════════

print("\n" + "="*60)
print("B) Two-wave IC  (dt=CFL·dx/lam, pseudo-spectral reference)")
print("="*60)

T_tw          = .8
Nx_tw         = [50, 100, 150, 200, 400]
gamma         = 1.4
xmin_tw, xmax_tw = -12.0, 12.0   # wider domain: IC tail at ±10 is ~1e-11, floor < 1e-9

def rho_0_tw(xi): return 1.5 + np.exp(-xi**2 / 4)
def u_0_tw(xi):   return -0.2 * np.exp(-xi**2 / 4)

from exact_two_wave import build_two_wave_reference, build_gamma3_exact

if gamma == 3.0:
    print("  gamma=3: using exact Burgers reference (machine precision).")
    _ref_tw     = build_gamma3_exact(rho_0_tw, u_0_tw, K=K,
                                     xmin=xmin_tw, xmax=xmax_tw, T_max=T_tw)
    exact_fn_tw = lambda x, rho, u, t: _ref_tw.l2_error(x, rho, u, t)
else:
    print("  Building spectral reference ...")
    _ref_tw     = build_two_wave_reference(rho_0_tw, u_0_tw,
                      K=K, gamma=gamma, xmin=xmin_tw, xmax=xmax_tw, T_target=T_tw)
    exact_fn_tw = _ref_tw.l2_error

cases_tw = [
    (0, 1, 'Euler ISL r=1'),
    (0, 2, 'Euler ISL r=2'),
    (1, 1, '2nd order ISL  r=1'),
    (1, 2, '2nd order ISL r=2'),
]

conv_tw = {}
for ck_val, r_val, label in cases_tw:
    q        = ck_val + 1          # temporal order: Euler→1, BFECC→2
    p        = 2 * r_val + 1
    # With CFL dt (dt∝dx), the global error is O(dt^q)=O(dx^q).
    # At γ=3, BFECC reduces to implicit Euler, so both give q_eff=1.
    expected = p*q / (q+1) 
    print(f"\n[{label}]  p={p}  q={q}  dt=CFL·dx/lam  expected O(h^{expected:.2f})")
    errs_rho, errs_u, T_list, dxs = [], [], [], []

    for Nx in Nx_tw:
        dx  = (xmax_tw - xmin_tw) / Nx
        x   = np.linspace(xmin_tw + dx / 2, xmax_tw - dx / 2, Nx)
        sol = LPTSolver(x, rho_0_tw, u_0_tw, K=K, gamma=gamma,
                        r=r_val, ck=ck_val, kernel=kernel, **solver_kw)


        t = 0.0
        while t < T_tw:
            lam  = np.abs(sol.u) + _sound_speed(sol.rho, K, gamma)
            dt_n =  dx**(p/(q+1)) / float(np.max(lam))
            dt_actual = min(dt_n, T_tw - t)
            sol.step(dt_actual)
            t += dt_actual
        #print(f'SL done')

        err_rho, err_u = exact_fn_tw(x, sol.rho, sol.u, t)
        #print(f'spectral done')
        errs_rho.append(err_rho); errs_u.append(err_u)
        T_list.append(t); dxs.append(dx)

        if len(errs_rho) > 1:
            p_rho = np.log(errs_rho[-2] / err_rho) / np.log(dxs[-2] / dx)
            p_u   = np.log(errs_u[-2]   / err_u)   / np.log(dxs[-2] / dx)
            print(f"  Nx={Nx:4d}  dx={dx:.4f}  dt={dt_n:.3e}  "
                  f"err_rho={err_rho:.3e} (p={p_rho:.2f})  err_u={err_u:.3e} (p={p_u:.2f})")
        else:
            print(f"  Nx={Nx:4d}  dx={dx:.4f}  dt={dt_n:.3e}  "
                  f"err_rho={err_rho:.3e}  err_u={err_u:.3e}")

    dxs_arr = np.array(dxs)
    conv_tw[label] = {
        'Nx'       : np.array(Nx_tw),
        'dx'       : dxs_arr,
        'T'        : np.array(T_list),
        'err_rho'  : np.array(errs_rho),
        'err_u'    : np.array(errs_u),
        'order_rho': np.array(_orders_from_list(errs_rho, dxs)),
        'order_u'  : np.array(_orders_from_list(errs_u,   dxs)),
    }
# %%


print("\n── Empirical orders (two-wave, dt=dx**(p/(q+1))/lam) ──")
for ck_val, r_val, label in cases_tw:
    q        = ck_val + 1          # temporal order
    p        = 2 * r_val + 1
    expected = p * q /(q + 1)                  
    r = conv_tw[label]
    print(f"  [{label}]  expected O(h^{expected:.0f})  "
          f"rho: {np.round(r['order_rho'][:-1], 2)}  u: {np.round(r['order_u'][:-1], 2)}")


colors_tw  = {cl: f'C{i}' for i, (*_, cl) in enumerate(cases_tw)}
markers_tw = {cl: _markers[i % len(_markers)] for i, (*_, cl) in enumerate(cases_tw)}
anchor_tw  = cases_tw[0][-1]


fig, axes = plt.subplots(1, 2, figsize=(14, 5))
fig.suptitle(f"Spatial convergence — two-wave IC, dt=CFL·dx/lam  (T={T_tw}, γ={gamma}, spectral ref)\n"
             f"{param_str}", fontsize=9)
for ax, key, ylbl in zip(axes, ['err_rho', 'err_u'], ['ρ', 'u']):
    for *_, label in cases_tw:
        r = conv_tw[label]
        ax.loglog(r['dx'], r[key], marker=markers_tw[label], linestyle='-',
                  color=colors_tw[label], label=label)
    dx_ref2  = conv_tw[anchor_tw]['dx']
    err_ref2 = conv_tw[anchor_tw][key]
    for slope, ls, col in [(2,'--','silver'),(3.33 ,':','gray'),(4.67,'-.','darkgray')]:
        ax.loglog(dx_ref2, float(err_ref2[0]) * (dx_ref2 / float(dx_ref2[0]))**slope,
                  linestyle=ls, color=col, alpha=0.45, label=f'O(h^{slope})')
    ax.set_xlabel('dx'); ax.set_ylabel(f'L2 error ({ylbl})')
    ax.set_title(ylbl); ax.legend(fontsize=8); ax.grid(True, which='both', alpha=0.3)
plt.tight_layout()
plt.show()
