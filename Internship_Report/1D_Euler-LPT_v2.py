#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
1D Euler – shock-adaptive LPT  (v2 rewrite)
============================================

Two Lagrangian grids track the two families of characteristics:

  x-grid (R+ chars) :  x_i(t) advects at  λ+ = u + c,  carries  R+ = u + 2c
  y-grid (R- chars) :  y_j(t) advects at  λ- = u - c,  carries  R- = u - 2c

Algorithm per time step
-----------------------
1.  Reconstruct (ρ, u) on each grid by cross-interpolating the complementary
    Riemann invariant from the other grid.
2.  Compute the adaptive dt = min crossing time across both grids
    (i.e. step to the very next characteristic crossing).
3.  Advance x by λ+·dt,  y by λ-·dt.
4.  Detect collision blocks on each grid independently:
      - Interpolate the complementary invariant at the outermost points of
        each block to recover (ρ_L, u_L) and (ρ_R, u_R).
      - Solve HLL → U* = (ρ*, u*).
      - Replace the k-point block by a single point carrying R*;
        the grid loses k−1 points per collision.
    The other grid is never touched during this step.
5.  No re-anchoring to a fixed Eulerian mesh.  The grids coarsen naturally
    as shocks accumulate.
"""

import numpy as np
import matplotlib.pyplot as plt


# ─────────────────────────────────────────────────────────────────────────────
# EOS  /  Riemann invariants
# ─────────────────────────────────────────────────────────────────────────────

def sound_speed(rho, K=1.0, gamma=2.0):
    return np.sqrt(np.maximum(K * gamma * rho**(gamma - 1.0), 0.0))


def to_riemann(rho, u, K=1.0, gamma=2.0):
    """(ρ, u)  →  (R+, R-, c)"""
    c = sound_speed(rho, K, gamma)
    return u + 2.0*c/(gamma-1),  u - 2.0*c/(gamma-1),  c


def from_riemann(Rp, Rm, K=1.0, gamma=2.0):
    """(R+, R-)  →  (ρ, u)   for a γ-law gas"""
    u   = 0.5  * (Rp + Rm)
    c   = 0.25 * (gamma - 1.0) * (Rp - Rm)
    c   = np.maximum(c, 0.0)          # safeguard against negative sound speed
    rho = (c**2 / (K * gamma))**(1.0 / (gamma - 1.0))
    return rho, u


def hll_state(rhoL, uL, rhoR, uR, K=1.0, gamma=2.0):
    """
    HLL intermediate conserved state  U* = (ρ*, u*).

    Vectorised: works on scalars or equal-length arrays.

    U  = [ρ,  ρu]
    F  = [ρu, ρu² + p]

    U* = (S_R · U_R  −  S_L · U_L  −  F_R  +  F_L) / (S_R − S_L)
    """
    pL = K * rhoL**gamma;   pR = K * rhoR**gamma
    cL = sound_speed(rhoL, K, gamma)
    cR = sound_speed(rhoR, K, gamma)

    SL = np.minimum(uL - cL,  uR - cR)
    SR = np.maximum(uL + cL,  uR + cR)
    dS = SR - SL                               # > 0 always

    rho_star  = (SR*rhoR        - SL*rhoL        - (rhoR*uR        - rhoL*uL       )) / dS
    rhoU_star = (SR*rhoR*uR     - SL*rhoL*uL     - (rhoR*uR**2+pR  - rhoL*uL**2-pL)) / dS
    u_star    = rhoU_star / rho_star

    return rho_star, u_star


# ─────────────────────────────────────────────────────────────────────────────
# HLLC Eulerian solver  (reference scheme)
# ─────────────────────────────────────────────────────────────────────────────

def _euler_flux(U, K=1.0, gamma=2.0):
    rho = U[0];  u = U[1] / rho
    p   = K * rho**gamma
    return np.array([rho*u, rho*u**2 + p])

def _hllc_flux(UL, UR, K=1.0, gamma=2.0):
    rhoL, uL = UL[0], UL[1]/UL[0]
    rhoR, uR = UR[0], UR[1]/UR[0]
    pL = K * rhoL**gamma
    pR = K * rhoR**gamma
    cL = sound_speed(rhoL, K, gamma)
    cR = sound_speed(rhoR, K, gamma)

    SL = min(uL - cL, uR - cR)
    SR = max(uL + cL, uR + cR)

    FL = _euler_flux(UL, K, gamma)
    FR = _euler_flux(UR, K, gamma)

    if SL >= 0:  return FL
    if SR <= 0:  return FR

    # Contact wave speed S* from momentum jump condition
    num = pR - pL + rhoL*uL*(SL - uL) - rhoR*uR*(SR - uR)
    den = rhoL*(SL - uL) - rhoR*(SR - uR)
    Ss  = num / den

    # Intermediate state U*_K = rho_K*(S_K - u_K)/(S_K - S*) * [1, S*]
    if Ss >= 0:
        fac    = rhoL * (SL - uL) / (SL - Ss)
        UL_s   = np.array([fac, fac * Ss])
        return FL + SL * (UL_s - UL)
    else:
        fac    = rhoR * (SR - uR) / (SR - Ss)
        UR_s   = np.array([fac, fac * Ss])
        return FR + SR * (UR_s - UR)

def compute_dt_hll(U, dx, CFL=0.5, K=1.0, gamma=2.0):
    rho = U[0];  u = U[1] / rho
    c   = sound_speed(rho, K, gamma)
    return CFL * dx / np.max(np.abs(u) + c)

def godunov_step(U, dx, dt, K=1.0, gamma=2.0):
    nx    = U.shape[1]
    U_new = U.copy()
    F     = np.zeros((2, nx + 1))
    for i in range(1, nx):
        F[:, i] = _hllc_flux(U[:, i-1], U[:, i], K, gamma)
    F[:, 0]  = F[:, 1]
    F[:, -1] = F[:, -2]
    for i in range(nx):
        U_new[:, i] -= dt / dx * (F[:, i+1] - F[:, i])
    return U_new


# ─────────────────────────────────────────────────────────────────────────────
# LPT core utilities
# ─────────────────────────────────────────────────────────────────────────────

def next_crossing_dt(pos, speed, dt_fallback=1e10):
    """
    Smallest positive time before any adjacent pair in sorted `pos`
    collides, given their current advection speeds.

    Two adjacent points i and i+1 (pos[i] < pos[i+1]) converge when
    speed[i] > speed[i+1], with crossing time:

        dt = (pos[i+1] - pos[i]) / (speed[i] - speed[i+1])

    Parameters
    ----------
    pos          : 1-D array, sorted positions
    speed        : 1-D array, advection speed at each position
    dt_fallback  : returned when no convergent pairs exist

    Returns
    -------
    dt : float
    """
    gap    = pos[1:] - pos[:-1]           # > 0 because pos is sorted
    dspeed = speed[:-1] - speed[1:]       # > 0  ↔  left char overtakes right

    converging = (gap > 0.0) & (dspeed > 0.0)
    if not np.any(converging):
        return dt_fallback

    return float(np.min(gap[converging] / dspeed[converging]))


def merge_crossings(pos, R_own, pos_other, R_other,
                    is_Rplus, K=1.0, gamma=2.0, tol=1e-10):
    """
    After advancing `pos` by one time step, collapse every group of
    co-located (crossed) points into a single HLL star-state point.

    For each collision block [i … j]  (k = j−i+1 points):
      1. Interpolate the complementary Riemann invariant from `pos_other`
         at positions pos[i] (left) and pos[j] (right).
      2. Reconstruct (ρ_L, u_L) and (ρ_R, u_R).
      3. Compute U* via HLL.
      4. Store a single merged point with R_merged = u* ± 2c*
         at the midpoint of the block.
      → k points become 1;  k−1 points are removed.

    The other grid (`pos_other`, `R_other`) is NOT modified here.

    Parameters
    ----------
    pos       : 1-D array, positions of this grid (may contain crossings)
    R_own     : Riemann invariant carried by this grid (R+ or R-)
    pos_other : 1-D array, positions of the other grid (sorted)
    R_other   : Riemann invariant carried by the other grid
    is_Rplus  : True  if R_own = R+ = u+2c  (forward / x-grid)
                False if R_own = R- = u-2c  (backward / y-grid)
    tol       : absolute gap threshold to declare a collision

    Returns
    -------
    new_pos   : merged positions  (len ≤ len(pos))
    new_R     : merged R_own values
    n_lost    : total number of points removed
    Rc_new    : complementary R at every new position
                  surviving points → interpolated from (pos_other, R_other)
                  merged blocks    → R_comp_star from HLL (NOT re-interpolated)
    """
    pos       = np.asarray(pos,       dtype=float)
    R_own     = np.asarray(R_own,     dtype=float)
    pos_other = np.asarray(pos_other, dtype=float)
    R_other   = np.asarray(R_other,   dtype=float)

    # Sort other grid for np.interp (should already be sorted, just in case)
    sort_idx  = np.argsort(pos_other)
    pos_other = pos_other[sort_idx]
    R_other   = R_other[sort_idx]

    gap       = pos[1:] - pos[:-1]
    colliding = gap < tol              # True where adjacent pts have crossed

    new_pos = []
    new_R   = []
    Rc_new  = []
    n_lost  = 0
    i = 0
    n = len(pos)

    while i < n:
        # Extend block as far as collision flags are True
        j = i
        while j < n - 1 and colliding[j]:
            j += 1

        if j == i:
            # Isolated (non-colliding) point – keep as-is
            new_pos.append(pos[i])
            new_R.append(R_own[i])
            Rc_new.append(float(np.interp(pos[i], pos_other, R_other)))
            i += 1
            continue

        # ── Collision block  pos[i … j] ──────────────────────────────────────
        x_L = pos[i];   x_R = pos[j]

        # Complementary invariant at outermost edges by interpolation
        R_comp_L = float(np.interp(x_L, pos_other, R_other))
        R_comp_R = float(np.interp(x_R, pos_other, R_other))

        if is_Rplus:                   # R_own = R+  →  complement is R-
            Rp_L, Rm_L = R_own[i], R_comp_L
            Rp_R, Rm_R = R_own[j], R_comp_R
        else:                          # R_own = R-  →  complement is R+
            Rp_L, Rm_L = R_comp_L, R_own[i]
            Rp_R, Rm_R = R_comp_R, R_own[j]

        rhoL, uL = from_riemann(Rp_L, Rm_L, K, gamma)
        rhoR, uR = from_riemann(Rp_R, Rm_R, K, gamma)

        rho_s, u_s = hll_state(rhoL, uL, rhoR, uR, K, gamma)
        x_merged   = 0.5 * (x_L + x_R)
        c_s        = float(sound_speed(rho_s, K, gamma))

        if is_Rplus:
            R_merged    = u_s + 2.0 * c_s / (gamma - 1.0)   # R+*
            R_comp_star = u_s - 2.0 * c_s / (gamma - 1.0)   # R-*  from HLL
        else:
            R_merged    = u_s - 2.0 * c_s / (gamma - 1.0)   # R-*
            R_comp_star = u_s + 2.0 * c_s / (gamma - 1.0)   # R+*  from HLL

        new_pos.append(x_merged)
        new_R.append(float(R_merged))
        Rc_new.append(float(R_comp_star))
        n_lost += j - i                # k points → 1,  k-1 removed
        i = j + 1

    return np.array(new_pos), np.array(new_R), n_lost, np.array(Rc_new)


# ─────────────────────────────────────────────────────────────────────────────
# Plotting helper
# ─────────────────────────────────────────────────────────────────────────────

def fig_show(x, u, rho, x_hll, u_hll, rho_hll, t, plot_hll=True):
    plt.figure(figsize=(8, 6))

    # Density (rho)
    plt.plot(x, rho, label='rho (LPT)', marker='+')
    if plot_hll:
        plt.plot(x_hll, rho_hll, label='rho (HLL)', linestyle='--')

    # Velocity (u)
    plt.plot(x, u, label='u (LPT)', marker='o')
    if plot_hll:
        plt.plot(x_hll, u_hll, label='u (HLL)', linestyle='--', alpha = 0.7)

    plt.title(f"Solution at time t = {t:.4f}")
    plt.xlabel("x")
    plt.ylabel("Values")

    # Optional: adjust limits to cover both variables
    plt.ylim([-1, 2])

    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.show()


# ─────────────────────────────────────────────────────────────────────────────
# Initial conditions
# ─────────────────────────────────────────────────────────────────────────────

xmin, xmax = -1.,4 
Nx    = 256
dx    = 1.0 / Nx
K     = 1.0
gamma = 1.1

t      = 0.0
T_max  = 0.2

freq   = 10          # plot every `freq` steps
k_step = 0

'''
def rho_0(xi):
    return 0.125 + (1.0 - 0.125) * 0.5 * (1.0 - np.tanh((xi - 0.5) / 0.02))

def u_0(xi):
    return np.zeros_like(xi)
'''

def rho_0(x):
    epsilon = 2
    return 1.0 + epsilon * np.exp(-x**2)

def u_0(x):
    epsilon = .8
    return epsilon * np.exp(-x**2)

# LPT grids – both start co-located, then drift freely (no re-anchoring)
x = np.linspace(xmin + dx/2, xmax - dx/2, Nx)
y = x.copy()

rho = rho_0(x)
u   = u_0(x)
R_plus, R_minus, _ = to_riemann(rho, u, K, gamma)

# Initial state saved as fallback for out-of-range cross-interpolation.
# HLL reference grid (fixed Eulerian)
x_hll   = np.linspace(xmin + dx/2, xmax - dx/2, Nx)
rho_hll = rho_0(x_hll)
U_hll   = np.vstack((rho_hll, rho_hll * u_0(x_hll)))

# Carried complementary invariants: each grid particle remembers the last
# accurate R_comp it received.  At t=0 both grids are co-located so the
# values are exact; after that, merged blocks receive R_comp_star from HLL
# instead of the wrong linear interpolation across the initial discontinuity.
Rm_at_x = R_minus.copy()   # R- carried by x-grid particles
Rp_at_y = R_plus.copy()    # R+ carried by y-grid particles


# ─────────────────────────────────────────────────────────────────────────────
# Time loop
# ─────────────────────────────────────────────────────────────────────────────

%matplotlib inline

while t < T_max:

    # ── 1. Reconstruct (ρ, u) and characteristic speeds on both grids ─────────
    #
    #  Use the CARRIED complementary invariants (set at end of previous step).
    #  Merged shock points have R_comp from the HLL star state, not interpolation.
    rho, u  = from_riemann(R_plus, Rm_at_x, K, gamma)
    c       = sound_speed(rho, K, gamma)
    lam_p   = u + c                                      # λ+  for x-grid

    rho_y, u_y = from_riemann(Rp_at_y, R_minus, K, gamma)
    c_y        = sound_speed(rho_y, K, gamma)
    lam_m      = u_y - c_y                               # λ-  for y-grid

    # ── 2. Adaptive dt: march to the very first crossing ──────────────────────
    dt_p   = next_crossing_dt(x, lam_p)
    dt_m   = next_crossing_dt(y, lam_m)
    dt_lpt = min(dt_p, dt_m)
    dt_hll = compute_dt_hll(U_hll, dx, K=K, gamma=gamma)
    dt     = min(dt_lpt , dt_hll, T_max - t)
    
    #dt = dt_lpt * 0.99

    if dt <= 0.0:
        break

    # ── 3. Advance both Lagrangian grids ──────────────────────────────────────
    x = x + lam_p * dt
    y = y + lam_m * dt
    # Rm_at_x / Rp_at_y are carried unchanged (Riemann invariants along chars)

    # ── 4. Detect and resolve crossings on each grid independently ────────────
    #
    #  merge_crossings now returns Rc_new: complementary R at every new point.
    #  For collision blocks Rc_new = R_comp_star from HLL (exact star state).
    #  For surviving points Rc_new = np.interp from other grid (smooth region).
    x, R_plus,  n_x, Rm_at_x = merge_crossings(
        x, R_plus,  y, R_minus, is_Rplus=True,  K=K, gamma=gamma, tol = 1e-6)

    y, R_minus, n_y, Rp_at_y = merge_crossings(
        y, R_minus, x, R_plus,  is_Rplus=False, K=K, gamma=gamma, tol = 1e-6)

    # Clip to domain (transmissive BC)
    x = np.clip(x, xmin, xmax)
    y = np.clip(y, xmin, xmax)

    # Boundary injection when edge gap > dx; new points get R_comp by interp
    if x[0] - xmin > dx:
        x       = np.concatenate([[xmin], x])
        R_plus  = np.concatenate([[R_plus[0]], R_plus])
        Rm_at_x = np.concatenate([[float(np.interp(xmin, y, R_minus))], Rm_at_x])
    if xmax - x[-1] > dx:
        x       = np.concatenate([x, [xmax]])
        R_plus  = np.concatenate([R_plus, [R_plus[-1]]])
        Rm_at_x = np.concatenate([Rm_at_x, [float(np.interp(xmax, y, R_minus))]])

    if y[0] - xmin > dx:
        y       = np.concatenate([[xmin], y])
        R_minus = np.concatenate([[R_minus[0]], R_minus])
        Rp_at_y = np.concatenate([[float(np.interp(xmin, x, R_plus))], Rp_at_y])
    if xmax - y[-1] > dx:
        y       = np.concatenate([y, [xmax]])
        R_minus = np.concatenate([R_minus, [R_minus[-1]]])
        Rp_at_y = np.concatenate([Rp_at_y, [float(np.interp(xmax, x, R_plus))]])

    # ── 5. HLL reference step ─────────────────────────────────────────────────
    U_hll   = godunov_step(U_hll, dx, dt, K, gamma)
    rho_hll = U_hll[0]
    m_hll   = U_hll[1]

    # ── 6. Plot ───────────────────────────────────────────────────────────────
    if k_step % freq == 0:
        rho_p, u_p = from_riemann(R_plus, Rm_at_x, K, gamma)
        fig_show(x, u_p, rho_p, x_hll, m_hll / rho_hll, rho_hll, t)

    t      += dt
    k_step += 1

# ── Final plot ────────────────────────────────────────────────────────────────
%matplotlib qt
rho_p, u_p = from_riemann(R_plus, Rm_at_x, K, gamma)
fig_show(x, u_p, rho_p, x_hll, m_hll / rho_hll, rho_hll, t, plot_hll=True)
print(f"\nFinal:  t = {t:.4f},  N_x = {len(x)},  N_y = {len(y)},  steps = {k_step}")
