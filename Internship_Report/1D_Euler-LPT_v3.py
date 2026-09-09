#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  3 14:01:50 2026

@author: dantsar
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
1D Euler – single-grid shock-adaptive LPT  (v3)
================================================

Key differences from v2:
  - ONE grid  {pos, R+, R-}  instead of two separate grids.
  - When advancing along λ+, R+ is exact (frozen).  R- is updated with
    the first-order coupling correction:
        ΔR- = 2c · (∂R-/∂x) · dt        (derived from ∂_t R- + λ- ∂_x R- = 0)
  - At a forward-shock collision: HLL gives U*; BOTH R+* and R-* are set
    from U*, and the shock point is assigned the Rankine-Hugoniot speed
        s = [ρu] / [ρ]   (from mass RH)
    instead of the characteristic speed u*+c*.  This fixes the wrong
    shock speed produced by v2.
  - Symmetric treatment for backward shocks (λ- crossing).
"""

import numpy as np
import matplotlib.pyplot as plt


# ─────────────────────────────────────────────────────────────────────────────
# EOS / Riemann invariants
# ─────────────────────────────────────────────────────────────────────────────

def sound_speed(rho, K=1.0, gamma=2.0):
    return np.sqrt(np.maximum(K * gamma * rho**(gamma - 1.0), 0.0))

def to_riemann(rho, u, K=1.0, gamma=2.0):
    c = sound_speed(rho, K, gamma)
    return u + 2.0*c, u - 2.0*c

def from_riemann(Rp, Rm, K=1.0, gamma=2.0):
    u   = 0.5  * (Rp + Rm)
    c   = np.maximum(0.25 * (Rp - Rm), 0.0)
    rho = (c**2 / (K * gamma))**(1.0 / (gamma - 1.0))
    return rho, u


# ─────────────────────────────────────────────────────────────────────────────
# Riemann solvers
# ─────────────────────────────────────────────────────────────────────────────

def hll_state(rhoL, uL, rhoR, uR, K=1.0, gamma=2.0):
    """HLL intermediate conserved state U* = (ρ*, u*)."""
    pL = K * rhoL**gamma;  pR = K * rhoR**gamma
    cL = sound_speed(rhoL, K, gamma)
    cR = sound_speed(rhoR, K, gamma)
    SL = np.minimum(uL - cL, uR - cR)
    SR = np.maximum(uL + cL, uR + cR)
    dS = SR - SL
    rho_s  = (SR*rhoR       - SL*rhoL       - (rhoR*uR       - rhoL*uL       )) / dS
    rhoU_s = (SR*rhoR*uR    - SL*rhoL*uL    - (rhoR*uR**2+pR - rhoL*uL**2-pL)) / dS
    u_s    = rhoU_s / rho_s
    return float(rho_s), float(u_s)


def rh_shock_speed(rhoL, uL, rhoR, uR):
    """
    Rankine-Hugoniot shock speed from mass conservation:
        s = [ρu] / [ρ]

    This is what the shock ACTUALLY moves at, which is generally different
    from any characteristic speed u* ± c*.  Using u* + c* for the merged
    point (as in v2) is the primary cause of the wrong shock speed.
    """
    drho = rhoR - rhoL
    if abs(drho) < 1e-14:
        return 0.5 * (uL + uR)
    return (rhoR * uR - rhoL * uL) / drho


# ─────────────────────────────────────────────────────────────────────────────
# HLL Eulerian reference solver  (unchanged)
# ─────────────────────────────────────────────────────────────────────────────

def _euler_flux(U, K=1.0, gamma=2.0):
    rho = U[0];  u = U[1] / rho
    return np.array([rho*u, rho*u**2 + K*rho**gamma])

def _hll_flux(UL, UR, K=1.0, gamma=2.0):
    rhoL, uL = UL[0], UL[1]/UL[0]
    rhoR, uR = UR[0], UR[1]/UR[0]
    cL = sound_speed(rhoL, K, gamma)
    cR = sound_speed(rhoR, K, gamma)
    SL = min(uL-cL, uR-cR);  SR = max(uL+cL, uR+cR)
    FL = _euler_flux(UL, K, gamma);  FR = _euler_flux(UR, K, gamma)
    if   SL >= 0: return FL
    elif SR <= 0: return FR
    return (SR*FL - SL*FR + SL*SR*(UR - UL)) / (SR - SL)

def compute_dt_hll(U, dx, CFL=0.5, K=1.0, gamma=2.0):
    rho = U[0];  u = U[1]/rho
    return CFL * dx / np.max(np.abs(u) + sound_speed(rho, K, gamma))

def godunov_step(U, dx, dt, K=1.0, gamma=2.0):
    nx = U.shape[1];  U_new = U.copy()
    F  = np.zeros((2, nx+1))
    for i in range(1, nx):
        F[:, i] = _hll_flux(U[:, i-1], U[:, i], K, gamma)
    F[:, 0] = F[:, 1];  F[:, -1] = F[:, -2]
    for i in range(nx):
        U_new[:, i] -= dt/dx * (F[:, i+1] - F[:, i])
    return U_new


# ─────────────────────────────────────────────────────────────────────────────
# LPT helpers
# ─────────────────────────────────────────────────────────────────────────────

def next_crossing_dt(pos, speed, dt_fallback=1e10):
    """Min positive time before any adjacent converging pair in sorted pos crosses."""
    gap    = pos[1:] - pos[:-1]
    dspeed = speed[:-1] - speed[1:]
    conv   = (gap > 0) & (dspeed > 0)
    if not np.any(conv):
        return dt_fallback
    return float(np.min(gap[conv] / dspeed[conv]))


def lpt_step(pos, Rp, Rm, vel, dt, K=1.0, gamma=2.0, tol=1e-10):
    """
    Advance the single LPT grid by dt.

    Parameters
    ----------
    pos : sorted positions
    Rp  : R+ at each position
    Rm  : R- at each position
    vel : advection velocity at each position
          (characteristic speed  OR  RH shock speed for shock points)
    dt  : time step

    Returns
    -------
    new_pos, new_Rp, new_Rm, new_vel : updated arrays (possibly shorter)

    Algorithm
    ---------
    The grid is advanced along whichever λ± family has the next crossing.
    Say it is λ+ (forward shock):
      - R+ is exact along λ+  →  unchanged at non-shock points.
      - R- evolves along λ+ as:
            dR-/dt |_{along λ+} = (λ+ − λ−) ∂_x R- = 2c ∂_x R-
        so  ΔR- ≈ 2c (∂R-/∂x) dt   at each non-shock point.
      - At the shock pair (j, j+1):
          * Solve HLL  →  U* = (ρ*, u*).
          * Set  R+* = u* + 2c*,  R-* = u* − 2c*  from U*.
          * Compute the true RH shock speed:
                s = [ρu] / [ρ]   (mass conservation)
            and assign it as vel for the merged point.
            This replaces the WRONG choice  vel = u* + c*  from v2,
            which caused the shock to propagate at the wrong speed.
    """
    rho, u = from_riemann(Rp, Rm, K, gamma)
    c      = sound_speed(rho, K, gamma)
    lam_p  = u + c
    lam_m  = u - c

    # Determine dominant family for this step
    forward_shock = (next_crossing_dt(pos, lam_p) <= next_crossing_dt(pos, lam_m))

    # Advance all positions using the per-point velocity (respects shock pts)
    new_pos = pos + vel * dt

    # ── Update invariants ────────────────────────────────────────────────────
    new_Rp = Rp.copy()
    new_Rm = Rm.copy()

    if forward_shock:
        # R+ exact (frozen along λ+).
        # R- coupling: ΔR- = +2c (∂R-/∂x) dt
        dRm_dx  = np.gradient(Rm, pos)
        new_Rm += 2.0 * c * dRm_dx * dt
    else:
        # R- exact (frozen along λ-).
        # R+ coupling: ΔR+ = -2c (∂R+/∂x) dt
        dRp_dx  = np.gradient(Rp, pos)
        new_Rp -= 2.0 * c * dRp_dx * dt

    # ── Resolve collisions ───────────────────────────────────────────────────
    gap       = new_pos[1:] - new_pos[:-1]
    colliding = gap < tol

    out_pos = [];  out_Rp = [];  out_Rm = [];  out_vel = []
    i = 0;  n = len(new_pos)

    while i < n:
        # Find end of collision block
        j = i
        while j < n - 1 and colliding[j]:
            j += 1

        if j == i:
            # Isolated (smooth) point
            # Velocity for next step: characteristic speed (re-evaluated below)
            out_pos.append(new_pos[i])
            out_Rp.append(new_Rp[i])
            out_Rm.append(new_Rm[i])
            out_vel.append(None)      # placeholder; set after loop
            i += 1
            continue

        # ── Collision block  [i … j] → 1 shock point ─────────────────────
        rhoL, uL = from_riemann(new_Rp[i], new_Rm[i], K, gamma)
        rhoR, uR = from_riemann(new_Rp[j], new_Rm[j], K, gamma)

        rho_s, u_s = hll_state(rhoL, uL, rhoR, uR, K, gamma)
        c_s        = sound_speed(rho_s, K, gamma)

        Rp_s = u_s + 2.0 * c_s
        Rm_s = u_s - 2.0 * c_s

        # ── Rankine-Hugoniot shock speed ─────────────────────────────────
        #
        # This is the fix for the wrong shock speed in v2.
        # In v2, the merged point was advected at u*+c* (characteristic
        # speed at U*), but the shock physically moves at:
        #
        #     s = [ρu] / [ρ] = (ρR uR − ρL uL) / (ρR − ρL)
        #
        # which is generally ≠ u*±c*.  For a strong shock the difference
        # is O(1) in wave-speed units.
        s_rh = rh_shock_speed(rhoL, uL, rhoR, uR)

        out_pos.append(0.5 * (new_pos[i] + new_pos[j]))
        out_Rp.append(Rp_s)
        out_Rm.append(Rm_s)
        out_vel.append(s_rh)          # shock tracked at RH speed, not λ±

        i = j + 1

    # ── Fill in characteristic speeds for non-shock points ───────────────────
    out_pos = np.array(out_pos)
    out_Rp  = np.array(out_Rp)
    out_Rm  = np.array(out_Rm)

    rho_new, u_new = from_riemann(out_Rp, out_Rm, K, gamma)
    c_new          = sound_speed(rho_new, K, gamma)
    lam_new        = (u_new + c_new) if forward_shock else (u_new - c_new)

    out_vel_arr = np.where(
        [v is None for v in out_vel],   # None → smooth point
        lam_new,                         # characteristic speed
        [v if v is not None else 0.0 for v in out_vel]  # RH speed for shocks
    )

    return out_pos, out_Rp, out_Rm, out_vel_arr


# ─────────────────────────────────────────────────────────────────────────────
# Plotting
# ─────────────────────────────────────────────────────────────────────────────

xmin, xmax = 0., 1.0

def fig_show(pos, u, rho, x_hll, u_hll, rho_hll, t, plot_hll=True):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    axes[0].plot(pos, rho, 'o-', ms=2, lw=1.0, label=f'LPT (N={len(pos)})')
    if plot_hll:
        axes[0].plot(x_hll, rho_hll, lw=1.5, label='HLL')
    axes[0].set_title(f"Density   t = {t:.4f}")
    axes[0].set_xlabel("x"); axes[0].set_ylabel("ρ")
    axes[0].set_ylim([0, 2]); axes[0].legend()

    axes[1].plot(pos, u, 'o-', ms=2, lw=1.0, label=f'LPT (N={len(pos)})')
    if plot_hll:
        axes[1].plot(x_hll, u_hll, lw=1.5, label='HLL')
    axes[1].set_title(f"Velocity  t = {t:.4f}")
    axes[1].set_xlabel("x"); axes[1].set_ylabel("u")
    axes[1].set_ylim([-1, 2]); axes[1].legend()

    plt.tight_layout()
    plt.show()


# ─────────────────────────────────────────────────────────────────────────────
# Initial conditions
# ─────────────────────────────────────────────────────────────────────────────

Nx    = 100
dx    = 1.0 / Nx
K     = 1.0
gamma = 2.0

def rho_0(xi):
    return 0.125 + (1.0 - 0.125) * 0.5 * (1.0 - np.tanh((xi - 0.5) / 0.02))

def u_0(xi):
    return np.zeros_like(xi)


# Single LPT grid
pos = np.linspace(xmin + dx/2, xmax - dx/2, Nx)
rho = rho_0(pos)
u   = u_0(pos)
Rp, Rm = to_riemann(rho, u, K, gamma)
c       = sound_speed(rho, K, gamma)
lam_p   = u + c                      # initial velocity = λ+ (no shock points yet)
vel     = lam_p.copy()

# HLL reference grid
x_hll   = np.linspace(xmin + dx/2, xmax - dx/2, Nx)
rho_hll = rho_0(x_hll)
U_hll   = np.vstack((rho_hll, rho_hll * u_0(x_hll)))


# ─────────────────────────────────────────────────────────────────────────────
# Time loop
# ─────────────────────────────────────────────────────────────────────────────

t      = 0.0
T_max  = 0.1
freq   = 1
k_step = 0

%matplotlib inline

while t < T_max:

    # ── Adaptive dt: step to the next crossing ─────────────────────────────
    rho_cur, u_cur = from_riemann(Rp, Rm, K, gamma)
    c_cur          = sound_speed(rho_cur, K, gamma)
    lam_p_cur      = u_cur + c_cur
    lam_m_cur      = u_cur - c_cur

    dt_lpt = min(next_crossing_dt(pos, lam_p_cur),
                 next_crossing_dt(pos, lam_m_cur))
    dt_hll = compute_dt_hll(U_hll, dx, K=K, gamma=gamma)
    dt     = min(dt_lpt, dt_hll, T_max - t)

    if dt <= 0.0:
        break

    # ── LPT step (single grid) ─────────────────────────────────────────────
    pos, Rp, Rm, vel = lpt_step(pos, Rp, Rm, vel, dt, K=K, gamma=gamma)
    pos = np.clip(pos, xmin, xmax)

    # ── HLL step ───────────────────────────────────────────────────────────
    U_hll   = godunov_step(U_hll, dx, dt, K, gamma)
    rho_hll = U_hll[0]
    m_hll   = U_hll[1]

    # ── Plot ───────────────────────────────────────────────────────────────
    if k_step % freq == 0:
        rho_p, u_p = from_riemann(Rp, Rm, K, gamma)
        fig_show(pos, u_p, rho_p, x_hll, m_hll/rho_hll, rho_hll, t, plot_hll=False)

    t      += dt
    k_step += 1

# ── Final plot ────────────────────────────────────────────────────────────────
rho_p, u_p = from_riemann(Rp, Rm, K, gamma)
fig_show(pos, u_p, rho_p, x_hll, m_hll/rho_hll, rho_hll, t, plot_hll=True)
print(f"\nFinal:  t = {t:.4f},  N_lpt = {len(pos)},  steps = {k_step}")