#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  8 12:20:51 2026

@author: dantsar
"""


#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
1D Full Euler equations solved with HLLC (Godunov) and LPT (Lagrangian-Particle Tracking),
organized as classes.

Three characteristic families:
  λ- = u - c  →  R- = u - 2c/(γ-1)          (y-grid)
  λ0 = u      →  s  = p / ρ^γ  (entropy)     (z-grid)
  λ+ = u + c  →  R+ = u + 2c/(γ-1)          (x-grid)

HLLC runs on a finer grid (N_factor * Nx) and serves as the reference solution.
LPT runs on Nx cells with its own time step.
L2 error of LPT vs HLLC is reported at T_max.
"""

import numpy as np
import matplotlib.pyplot as plt
import time

import os, sys
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(current_dir)
from gp_interp import gp_interpolate
from CIP import cip_step



# %%


# ── Equation of state (full Euler: p = (γ-1)(E - ρu²/2)) ────────────────────

def pressure_from_cons(rho, rhou, E, gamma=1.4):
    u = rhou / rho
    return (gamma - 1.0) * (E - 0.5 * rho * u**2)

def sound_speed_prim(rho, p, gamma=1.4):
    return np.sqrt(gamma * p / rho)

def entropy_var(rho, p, gamma=1.4):
    """s = p / rho^gamma  (conserved along particle paths in smooth flow)."""
    return p / rho**gamma

def prim_from_chars(Rp, Rm, s, gamma=1.4):
    """
    Recover (rho, u, p) from the three characteristic variables:
      R+ = u + 2c/(γ-1),  R- = u - 2c/(γ-1),  s = p/ρ^γ

    Inversion:
      u = (R+ + R-) / 2
      c = (R+ - R-) * (γ-1) / 4
      ρ = (c² / (γ s))^(1/(γ-1))     [from c² = γ s ρ^(γ-1)]
      p = s ρ^γ
    """
    u   = 0.5 * (Rp + Rm)
    c   = np.maximum(0.25 * (gamma - 1.0) * (Rp - Rm), 1e-12)
    rho = (c**2 / (gamma * s))**(1.0 / (gamma - 1.0))
    p   = s * rho**gamma
    return rho, u, p


# ── HLLC Riemann solver helpers (full 3-equation Euler) ──────────────────────

def _cons_to_prim_full(U, gamma=1.4):
    rho = U[0]
    u   = U[1] / rho
    E   = U[2]
    p   = pressure_from_cons(rho, U[1], E, gamma)
    return rho, u, p, E

def _flux_full(U, gamma=1.4):
    rho, u, p, E = _cons_to_prim_full(U, gamma)
    return np.array([rho * u,
                     rho * u**2 + p,
                     u * (E + p)])

def _hllc_flux_full(UL, UR, gamma=1.4):
    rhoL, uL, pL, EL = _cons_to_prim_full(UL, gamma)
    rhoR, uR, pR, ER = _cons_to_prim_full(UR, gamma)
    cL = sound_speed_prim(rhoL, pL, gamma)
    cR = sound_speed_prim(rhoR, pR, gamma)

    SL = min(uL - cL, uR - cR)
    SR = max(uL + cL, uR + cR)

    num    = pR - pL + rhoL * uL * (SL - uL) - rhoR * uR * (SR - uR)
    denom  = rhoL * (SL - uL) - rhoR * (SR - uR)
    Sstar  = num / denom

    FL = _flux_full(UL, gamma)
    FR = _flux_full(UR, gamma)

    def _star_state(rho_k, u_k, p_k, E_k, S_k):
        fac = rho_k * (S_k - u_k) / (S_k - Sstar)
        return fac * np.array([
            1.0,
            Sstar,
            E_k / rho_k + (Sstar - u_k) * (Sstar + p_k / (rho_k * (S_k - u_k))),
        ])

    if SL >= 0:
        return FL
    elif SR <= 0:
        return FR
    elif Sstar >= 0:
        return FL + SL * (_star_state(rhoL, uL, pL, EL, SL) - UL)
    else:
        return FR + SR * (_star_state(rhoR, uR, pR, ER, SR) - UR)


# ── HLLCSolver ────────────────────────────────────────────────────────────────

class HLLCSolver:
    """
    Fixed-grid Godunov solver with HLLC Riemann fluxes for 1D full Euler.

    State: conserved vector U = [rho, rho*u, E], shape (3, Nx).
    """

    def __init__(self, x, rho_0, u_0, p_0, gamma=1.4, CFL=0.9):
        self.x     = x.copy()
        self.dx    = x[1] - x[0]
        self.gamma = gamma
        self.CFL   = CFL

        rho = rho_0(x)
        u   = u_0(x)
        p   = p_0(x)
        E   = p / (gamma - 1.0) + 0.5 * rho * u**2
        self.U = np.vstack((rho, rho * u, E))   # shape (3, Nx)

    @property
    def rho(self):
        return self.U[0].copy()

    @property
    def u(self):
        return (self.U[1] / self.U[0]).copy()

    @property
    def p(self):
        return pressure_from_cons(self.U[0], self.U[1], self.U[2], self.gamma)

    def compute_dt(self):
        rho = self.U[0]
        u   = self.U[1] / rho
        p   = pressure_from_cons(rho, self.U[1], self.U[2], self.gamma)
        c   = sound_speed_prim(rho, p, self.gamma)
        return self.CFL * self.dx / np.max(np.abs(u) + c)

    def step(self, dt):
        U  = self.U
        nx = U.shape[1]
        F  = np.zeros((3, nx + 1))

        for i in range(1, nx):
            F[:, i] = _hllc_flux_full(U[:, i - 1], U[:, i], self.gamma)

        # Transmissive boundary conditions
        F[:, 0]  = F[:, 1]
        F[:, -1] = F[:, -2]

        self.U = U - (dt / self.dx) * (F[:, 1:] - F[:, :-1])

    def run_to(self, T_max):
        t = 0.0
        while t < T_max:
            dt = min(self.compute_dt(), T_max - t)
            self.step(dt)
            t += dt


# ── LPTSolver ─────────────────────────────────────────────────────────────────

class LPTSolver:
    """
    Lagrangian-Particle Tracking solver for 1D full Euler.

    Three moving grids:
      x  → tracks λ+ = u + c  characteristics  (carries R+)
      y  → tracks λ- = u - c  characteristics  (carries R-)
      z  → tracks λ0 = u      characteristics  (carries entropy s = p/ρ^γ)

    After each step all grids are re-anchored to a fixed reference grid
    (semi-Lagrangian approach).
    """

    def __init__(self, x, rho_0, u_0, p_0, gamma=1.4):
        dx = x[1] - x[0]
        self.xmin  = x[0]  - dx / 2
        self.xmax  = x[-1] + dx / 2
        self.Nx    = len(x)
        self.gamma = gamma

        rho = rho_0(x)
        u   = u_0(x)
        p   = p_0(x)

        c = sound_speed_prim(rho, p, gamma)
        self.R_plus  = u + 2.0 * c / (gamma - 1.0)
        self.R_minus = u - 2.0 * c / (gamma - 1.0)
        self.S       = entropy_var(rho, p, gamma)

        self.rho = rho.copy()
        self.u   = u.copy()
        self.p   = p.copy()

        self.x = x.copy()   # λ+ grid
        self.y = x.copy()   # λ- grid
        self.z = x.copy()   # λ0 grid

    def _char_speeds(self):
        c     = sound_speed_prim(self.rho, self.p, self.gamma)
        lam_p = self.u + c
        lam_m = self.u - c
        lam_0 = self.u.copy()
        return lam_p, lam_m, lam_0

    def compute_shock_dt(self):
        lam_p, lam_m, lam_0 = self._char_speeds()
        dp = (self.x[1:] - self.x[:-1]) / np.abs(lam_p[1:] - lam_p[:-1] + 1e-10)
        dm = (self.y[1:] - self.y[:-1]) / np.abs(lam_m[1:] - lam_m[:-1] + 1e-10)
        dz = (self.z[1:] - self.z[:-1]) / np.abs(lam_0[1:] - lam_0[:-1] + 1e-10)
        return min(np.min(dp), np.min(dm), np.min(dz))

    def step(self, dt, interp_method='numpy', monotone=False):
        lam_p, lam_m, lam_0 = self._char_speeds()

        # Advect all three grids
        self.x += lam_p * dt
        self.y += lam_m * dt
        self.z += lam_0 * dt

        # Clip to domain and sort each grid independently
        x_s = np.clip(self.x, self.xmin, self.xmax)
        y_s = np.clip(self.y, self.xmin, self.xmax)
        z_s = np.clip(self.z, self.xmin, self.xmax)

        sx = np.argsort(x_s);  x_s = x_s[sx];  Rp_s = self.R_plus[sx]
        sy = np.argsort(y_s);  y_s = y_s[sy];  Rm_s = self.R_minus[sy]
        sz = np.argsort(z_s);  z_s = z_s[sz];  S_s  = self.S[sz]

        x_ref = np.linspace(
            self.xmin + (self.xmax - self.xmin) / self.Nx / 2,
            self.xmax - (self.xmax - self.xmin) / self.Nx / 2,
            self.Nx,
        )
        
        if interp_method == 'numpy':
            Rp_ref = np.interp(x_ref, x_s, Rp_s)
            Rm_ref = np.interp(x_ref, y_s, Rm_s)
            S_ref  = np.interp(x_ref, z_s, S_s)

        elif interp_method == 'GP':
            r, nugget = 12, 1e-10
            Rp_ref = gp_interpolate(x_s, Rp_s, x_ref, r, nugget=nugget)
            Rm_ref = gp_interpolate(y_s, Rm_s, x_ref, r, nugget=nugget)
            S_ref  = gp_interpolate(z_s, S_s,  x_ref, r, nugget=nugget)

        elif interp_method == 'GP_CIP':
            r, nugget = 12, 1e-10
            Rp_ref = cip_step(x_s, Rp_s, lam_p, dt, r, monotone=monotone)
            Rm_ref = cip_step(y_s, Rm_s, lam_m, dt, r, monotone=monotone)
            S_ref  = cip_step(z_s, S_s,  lam_0, dt, r, monotone=monotone)

        # Enforce positivity before inversion
        S_ref = np.maximum(S_ref, 1e-12)

        rho_new, u_new, p_new = prim_from_chars(Rp_ref, Rm_ref, S_ref, self.gamma)

        rho_new = np.maximum(rho_new, 1e-12)
        p_new   = np.maximum(p_new,   1e-12)

        c_new = sound_speed_prim(rho_new, p_new, self.gamma)
        self.R_plus  = u_new + 2.0 * c_new / (self.gamma - 1.0)
        self.R_minus = u_new - 2.0 * c_new / (self.gamma - 1.0)
        self.S       = entropy_var(rho_new, p_new, self.gamma)

        self.rho = rho_new
        self.u   = u_new
        self.p   = p_new

        # Re-anchor all grids to x_ref
        self.x = x_ref.copy()
        self.y = x_ref.copy()
        self.z = x_ref.copy()


# ── EulerSimulation ───────────────────────────────────────────────────────────

class EulerSimulation:
    """
    Runs HLLCSolver (fine grid) and LPTSolver (coarse grid) independently,
    then reports the L2 error of LPT vs HLLC at T_max.
    """

    def __init__(
        self,
        xmin, xmax, Nx,
        rho_0, u_0, p_0,
        gamma=1.4,
        CFL_hllc=0.9, CFL_lpt=0.99,
        T_max=0.4,
        N_factor=8,
        interp_method='numpy',
        monotone=False
    ):
        self.T_max    = T_max
        self.CFL_hllc = CFL_hllc
        self.CFL_lpt  = CFL_lpt

        dx_lpt  = abs(xmax - xmin) / Nx
        x_lpt   = np.linspace(xmin + dx_lpt / 2, xmax - dx_lpt / 2, Nx)

        Nx_hllc = N_factor * Nx
        dx_hllc = abs(xmax - xmin) / Nx_hllc
        x_hllc  = np.linspace(xmin + dx_hllc / 2, xmax - dx_hllc / 2, Nx_hllc)

        self.hllc = HLLCSolver(x_hllc, rho_0, u_0, p_0, gamma, CFL_hllc)
        self.lpt  = LPTSolver(x_lpt,  rho_0, u_0, p_0, gamma)

        self.x_lpt = x_lpt

    # ── Plotting ──────────────────────────────────────────────────────────────

    def _plot_final(self):
        lpt  = self.lpt
        hllc = self.hllc

        rho_hllc = np.interp(self.x_lpt, hllc.x, hllc.rho)
        u_hllc   = np.interp(self.x_lpt, hllc.x, hllc.u)
        p_hllc   = np.interp(self.x_lpt, hllc.x, hllc.p)

        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        fig.suptitle(f"Solutions at T_max = {self.T_max:.4f}")

        axes[0].plot(self.x_lpt, lpt.rho, label="LPT",  marker="+")
        axes[0].plot(self.x_lpt, rho_hllc, label="HLLC (ref)", linestyle="--")
        axes[0].set_title("Density ρ")
        axes[0].set_xlabel("x"); axes[0].legend(); axes[0].grid(True)

        axes[1].plot(self.x_lpt, lpt.u, label="LPT",  marker="+")
        axes[1].plot(self.x_lpt, u_hllc, label="HLLC (ref)", linestyle="--")
        axes[1].set_title("Velocity u")
        axes[1].set_xlabel("x"); axes[1].legend(); axes[1].grid(True)

        axes[2].plot(self.x_lpt, lpt.p, label="LPT",  marker="+")
        axes[2].plot(self.x_lpt, p_hllc, label="HLLC (ref)", linestyle="--")
        axes[2].set_title("Pressure p")
        axes[2].set_xlabel("x"); axes[2].legend(); axes[2].grid(True)

        plt.tight_layout()
        plt.show()

    # ── L2 error ──────────────────────────────────────────────────────────────

    def compute_l2_error(self):
        hllc = self.hllc
        lpt  = self.lpt

        rho_ref = np.interp(self.x_lpt, hllc.x, hllc.rho)
        u_ref   = np.interp(self.x_lpt, hllc.x, hllc.u)
        p_ref   = np.interp(self.x_lpt, hllc.x, hllc.p)

        dx = self.x_lpt[1] - self.x_lpt[0]
        err_rho = np.sqrt(np.sum((lpt.rho - rho_ref)**2) * dx)
        err_u   = np.sqrt(np.sum((lpt.u   - u_ref  )**2) * dx)
        err_p   = np.sqrt(np.sum((lpt.p   - p_ref  )**2) * dx)
        return err_rho, err_u, err_p

    # ── Main run ──────────────────────────────────────────────────────────────

    def run(self, freq=10, interp_method='numpy', monotone=False):
        """
        Run HLLC and LPT independently to T_max, then report L2 error.

        Parameters
        ----------
        freq : plot LPT every `freq` steps during its time loop
        """
        print("Running HLLC (reference) ...")
        start_hllc = time.time()
        self.hllc.run_to(self.T_max)
        end_hllc = time.time()
        print(f"  HLLC done in {end_hllc - start_hllc:.6f} seconds ")

        print("Running LPT ...")
        t = 0.0

        start_lpt = time.time()
        while t < self.T_max:
            dt_shock = self.lpt.compute_shock_dt()

            x  = self.lpt.x
            dx = min(x[1:] - x[:-1])
            dt = min(dt_shock * self.CFL_lpt, self.T_max - t, dx)
            self.lpt.step(dt, interp_method=interp_method, monotone=monotone)
            t += dt

        end_lpt = time.time()
        print(f"  LPT done in {end_lpt - start_lpt:.6f} seconds with {interp_method}")

        self._plot_final()

        err_rho, err_u, err_p = self.compute_l2_error()
        print(f"\nL2 error at T_max = {self.T_max}:")
        print(f"  rho : {err_rho:.6e}")
        print(f"  u   : {err_u:.6e}")
        print(f"  p   : {err_p:.6e}")
        return err_rho, err_u, err_p


# ── Initial conditions ────────────────────────────────────────────────────────

def rho_0(x):
    epsilon = 1
    return 1.5 + epsilon * np.exp(-x**2)

def u_0(x):
    epsilon = -0.2
    return epsilon * np.exp(-x**2)

def p_0(x):
    K, gamma = 1.0, 1.4
    return K * rho_0(x)**gamma   # isentropic initial pressure (s = K everywhere)


# ── Run ───────────────────────────────────────────────────────────────────────

interp_method = 'GP_CIP'

if __name__ == "__main__":
    sim = EulerSimulation(
        xmin=-6.0, xmax=6.0, Nx=512,
        rho_0=rho_0, u_0=u_0, p_0=p_0,
        gamma=1.4,
        CFL_hllc=0.9, CFL_lpt=1,
        T_max=1.,
        N_factor=1,
        monotone=False
    )
    sim.run(interp_method=interp_method)
