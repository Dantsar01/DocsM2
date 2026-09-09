
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Wed Apr 29 16:25:18 2026

@author: dantsar

1D Euler equations solved with HLLC (Godunov) and LPT (Lagrangian-Particle Tracking),
organized as classes.

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
from rbf_ga_weights_1d import RBFGAUniformGrid
from scipy.interpolate import CubicSpline, make_interp_spline


# ── Equation of state (module-level, pure functions) ──────────────────────────

def pressure(rho, K=1.0, gamma=2.0):
    return K * rho**gamma

def sound_speed(rho, K=1.0, gamma=2.0):
    return np.sqrt(K * gamma * rho**(gamma - 1.0))


# ── HLLC Riemann solver helpers ───────────────────────────────────────────────

def _conserved_to_primitive(U):
    rho = U[0]
    u   = U[1] / rho
    return rho, u

def _flux(U, K=1.0, gamma=2.0):
    rho, u = _conserved_to_primitive(U)
    p = pressure(rho, K, gamma)
    return np.array([rho * u, rho * u**2 + p])

def _hllc_flux(UL, UR, K=1.0, gamma=2.0):
    rhoL, uL = _conserved_to_primitive(UL)
    rhoR, uR = _conserved_to_primitive(UR)
    cL = sound_speed(rhoL, K, gamma)
    cR = sound_speed(rhoR, K, gamma)

    SL = min(uL - cL, uR - cR)
    SR = max(uL + cL, uR + cR)

    pL = pressure(rhoL, K, gamma)
    pR = pressure(rhoR, K, gamma)
    num   = pR - pL + rhoL * uL * (SL - uL) - rhoR * uR * (SR - uR)
    denom = rhoL * (SL - uL) - rhoR * (SR - uR)
    Sstar = num / denom

    FL = _flux(UL, K, gamma)
    FR = _flux(UR, K, gamma)

    if SL >= 0:
        return FL
    elif SR <= 0:
        return FR
    elif Sstar >= 0:
        coeff  = rhoL * (SL - uL) / (SL - Sstar)
        UstarL = coeff * np.array([1.0, Sstar])
        return FL + SL * (UstarL - UL)
    else:
        coeff  = rhoR * (SR - uR) / (SR - Sstar)
        UstarR = coeff * np.array([1.0, Sstar])
        return FR + SR * (UstarR - UR)


# ── HLLCSolver ────────────────────────────────────────────────────────────────

class HLLCSolver:
    """
    Fixed-grid Godunov solver with HLLC Riemann fluxes for 1D isentropic Euler.

    State: conserved vector U = [rho, rho*u], shape (2, Nx).
    """

    def __init__(self, x, rho_0, u_0, K=1.0, gamma=2.0, CFL=0.9):
        self.x     = x.copy()
        self.dx    = x[1] - x[0]
        self.K     = K
        self.gamma = gamma
        self.CFL   = CFL

        rho = rho_0(x)
        u   = u_0(x)
        self.U = np.vstack((rho, rho * u))   # shape (2, Nx)

    @property
    def rho(self):
        return self.U[0].copy()

    @property
    def u(self):
        return (self.U[1] / self.U[0]).copy()

    def compute_dt(self):
        rho = self.U[0]
        u   = self.U[1] / rho
        c   = sound_speed(rho, self.K, self.gamma)
        return self.CFL * self.dx / np.max(np.abs(u) + c)

    def step(self, dt):
        U  = self.U
        nx = U.shape[1]
        F  = np.zeros((2, nx + 1))

        for i in range(1, nx):
            F[:, i] = _hllc_flux(U[:, i - 1], U[:, i], self.K, self.gamma)

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
    Lagrangian-Particle Tracking solver for 1D isentropic Euler.

    Two moving grids:
      x  → tracks λ+ = u + c  characteristics  (carries R+)
      y  → tracks λ- = u - c  characteristics  (carries R-)

    After each step both grids are re-anchored to a fixed reference grid
    (semi-Lagrangian approach).
    """

    def __init__(self, x, rho_0, u_0, K=1.0, gamma=2.0, kernel='se', r=2, rbfga_eps=None, k_spline=5):
        dx = x[1] - x[0]
        self.xmin   = x[0]  - dx / 2
        self.xmax   = x[-1] + dx / 2
        self.Nx     = len(x)
        self.K      = K
        self.gamma  = gamma
        self.kernel   = kernel
        self.r        = r
        self.k_spline = k_spline

        eps = rbfga_eps if rbfga_eps is not None else 1e-3
        
        M = 2 * r + 1
        M_der = M + 2
        
        self.uniform_stencil_d1  = RBFGAUniformGrid(dx=dx, M=M, eps=eps, deriv=1)
        self.uniform_stencil_d2  = RBFGAUniformGrid(dx=dx, M=M_der, eps=eps, deriv=2)

        self.rho = rho_0(x)
        self.u   = u_0(x)
        self.x   = x.copy()   # λ+ grid
        self.y   = x.copy()   # λ- grid

        c = sound_speed(self.rho, K, gamma)
        self.R_plus  = self.u + 2 * c / (gamma - 1)
        self.R_minus = self.u - 2 * c / (gamma - 1)

        M  = 2 * r + 1
        w1 = self.uniform_stencil_d1.w
        Rp_p = np.pad(self.R_plus,  r, mode='edge')
        Rm_p = np.pad(self.R_minus, r, mode='edge')
        self.G_plus  = np.lib.stride_tricks.sliding_window_view(Rp_p, M) @ w1
        self.G_minus = np.lib.stride_tricks.sliding_window_view(Rm_p, M) @ w1

    def _char_speeds(self):
        c     = sound_speed(self.rho, self.K, self.gamma)
        lam_p = self.u + c
        lam_m = self.u - c
        return lam_p, lam_m

    def compute_shock_dt(self):
        lam_p, lam_m = self._char_speeds()
        dp = (self.x[1:] - self.x[:-1]) / np.abs(lam_p[1:] - lam_p[:-1] + 1e-6)
        dm = (self.y[1:] - self.y[:-1]) / np.abs(lam_m[1:] - lam_m[:-1] + 1e-6)
        return min(np.min(dp), np.min(dm))

    def step(self, dt, interp_method='numpy', monotone=False, midpoint_tracing=False):
        lam_p, lam_m = self._char_speeds()


        x       = self.x
        y       = self.y
        R_plus  = self.R_plus
        R_minus = self.R_minus
        eps     = self.uniform_stencil_d1.eps
        r       = self.r

        M  = 2 * r + 1
        w1 = self.uniform_stencil_d1.w                         # (M,)

        # ── Use carried G instead of recomputing from R ───────────────────────
        dR_plus  = self.G_plus
        dR_minus = self.G_minus

        # ── ddR = d/dx(G) with the same w1 stencil (one derivative of G) ─────
        Gp_pad    = np.pad(dR_plus,  r, mode='edge')
        Gm_pad    = np.pad(dR_minus, r, mode='edge')
        ddR_plus  = np.lib.stride_tricks.sliding_window_view(Gp_pad, M) @ w1
        ddR_minus = np.lib.stride_tricks.sliding_window_view(Gm_pad, M) @ w1


        # ── Cauchy–Kovalevskaya corrections ──────────────────────────────────
        g = self.gamma
        c = (g - 1) / 4 * (R_plus - R_minus)
        '''
        # ── Second derivatives for CK2 (finite differences, 4th-order) ───────
        dx  = self.x[1] - self.x[0]
        N   = self.Nx
        idx = np.arange(N)
        im2 = np.clip(idx - 2, 0, N - 1)
        im1 = np.clip(idx - 1, 0, N - 1)
        ip1 = np.clip(idx + 1, 0, N - 1)
        ip2 = np.clip(idx + 2, 0, N - 1)
        dx2 = dx * dx
        ddR_plus  = (-R_plus[ip2]  + 16*R_plus[ip1]  - 30*R_plus[idx]
                     + 16*R_plus[im1]  - R_plus[im2])  / (12 * dx2)
        ddR_minus = (-R_minus[ip2] + 16*R_minus[ip1] - 30*R_minus[idx]
                     + 16*R_minus[im1] - R_minus[im2]) / (12 * dx2)
        '''


        # ── Auxiliary spatial quantities ──────────────────────────────────────
        c_x     = (g - 1) / 4 * (dR_plus - dR_minus)        # ∂_x c
        lam_p_x = (g + 1) / 4 * dR_plus + (3 - g) / 4 * dR_minus   # ∂_x λ+
        lam_m_x = (3 - g) / 4 * dR_plus + (g + 1) / 4 * dR_minus   # ∂_x λ-

        # ── CK1: A± = dλ±/dt (material derivative of characteristic speed) ──
        # Algebraic identity: A+ = (3-g)/2·c·∂R-/∂x,  A- = -(3-g)/2·c·∂R+/∂x
        # (λ λ_x cancels in the compact form — do NOT add it separately).
        A_p =  (3 - g) / 2 * c * dR_minus
        A_m = -(3 - g) / 2 * c * dR_plus

        # ── CK2: B± = -∂_t A± + A± ∂_x λ± - λ± ∂_x A±  (doc eq., recursion P3) ──
        # ∂_t c   from ∂_t R± = -λ± R±_x (PDE):
        c_t = (g - 1) / 4 * (-lam_p * dR_plus + lam_m * dR_minus)

        # ∂_t(∂_x R±) = ∂_x(∂_t R±) = ∂_x(-λ± R±_x) = -λ±_x R±_x - λ± R±_xx
        dtdxRp = -lam_p_x * dR_plus  - lam_p * ddR_plus
        dtdxRm = -lam_m_x * dR_minus - lam_m * ddR_minus

        # ∂_t A± = (3-g)/2 · [∂_t c · R±_x  +  c · ∂_t(∂_x R±)] (chain rule)
        A_p_t = (3 - g) / 2 * (c_t * dR_minus + c * dtdxRm)
        A_m_t = -(3 - g) / 2 * (c_t * dR_plus  + c * dtdxRp)

        # ∂_x A± = (3-g)/2 · [∂_x c · R±_x  +  c · R±_xx]
        A_p_x = (3 - g) / 2 * (c_x * dR_minus + c * ddR_minus)
        A_m_x = -(3 - g) / 2 * (c_x * dR_plus  + c * ddR_plus)

        B_p = -A_p_t + A_p * lam_p_x - lam_p * A_p_x
        B_m = -A_m_t + A_m * lam_m_x - lam_m * A_m_x

        # ── Mesh advection ── Euler [active] ─────────────────────────────────
        # 1st-order trajectory: departure point = x + λ± · dt
        x += lam_p * dt
        y += lam_m * dt


        # ── Apply CK corrections: R̃± = R± - A± dt²/2 R±_x + B± dt³/6 R±_x ──
        ck2_p = (A_p * dR_plus)  * (dt**2 / 2) - (B_p * dR_plus)  * (dt**3 / 6)
        ck2_m = (A_m * dR_minus) * (dt**2 / 2) - (B_m * dR_minus) * (dt**3 / 6)

        R_plus  -= ck2_p
        R_minus -= ck2_m

        '''
        sort_x  = np.argsort(x)
        x_s     = x[sort_x]
        Rp_s    = R_plus[sort_x]
        lam_p_s = lam_p[sort_x]


        sort_y  = np.argsort(self.y)
        y_s     = np.clip(y, self.xmin, self.xmax)[sort_y]
        Rm_s    = R_minus[sort_y]
        lam_m_s = lam_m[sort_y]
        '''
        x_ref  = np.linspace(
            self.xmin + (self.xmax - self.xmin) / self.Nx / 2,
            self.xmax - (self.xmax - self.xmin) / self.Nx / 2,
            self.Nx,
        )

        nugget = 0
        if interp_method == 'numpy' :
            Rp_ref = np.interp(x_ref, x, R_plus)
            Rm_ref = np.interp(x_ref, y, R_minus)
            Gp_ref = np.interp(x_ref, x, dR_plus)
            Gm_ref = np.interp(x_ref, y, dR_minus)

        elif interp_method == 'cubic':
            Rp_ref = make_interp_spline(x, R_plus,  k=self.k_spline)(x_ref)
            Rm_ref = make_interp_spline(y, R_minus, k=self.k_spline)(x_ref)
            Gp_ref = make_interp_spline(x, dR_plus,  k=self.k_spline)(x_ref)
            Gm_ref = make_interp_spline(y, dR_minus, k=self.k_spline)(x_ref)

        elif interp_method == 'GP':
            Rp_ref = gp_interpolate(x, R_plus, x_ref, r, nugget=nugget,
                                    kernel=self.kernel, rbfga_eps = eps)
            Rm_ref = gp_interpolate(y, R_minus, x_ref, r, nugget=nugget,
                                    kernel=self.kernel, rbfga_eps = eps)
            Gp_ref = gp_interpolate(x, dR_plus, x_ref, r, nugget=nugget,
                                    kernel=self.kernel, rbfga_eps = eps)
            Gm_ref = gp_interpolate(y, dR_minus, x_ref, r, nugget=nugget,
                                    kernel=self.kernel, rbfga_eps = eps)
        '''
        elif interp_method == 'GP_CIP':
            Rp_ref = cip_step(x, Rp_s, lam_p_s, dt, r, monotone=monotone)
            Rm_ref = cip_step(y_s, Rm_s, lam_m_s, dt, r, monotone=monotone)
        '''

        u_ref = 0.5 * (Rp_ref + Rm_ref)
        # c is already known exactly from the Riemann invariants — no need to go
        # through rho.  With gamma=1.01 the exponent 1/(gamma-1)=100 amplifies any
        # error in c catastrophically, so avoid the roundtrip c→rho→c for R±.
        c_ref = 0.25 * (self.gamma - 1) * (Rp_ref - Rm_ref)

        # rho is only needed for _char_speeds() next step; clamp to prevent overflow
        # from the ^100 power (float64 overflows when c_ref > ~35).
        rho_new = (c_ref**2 / (self.K * self.gamma))**(1 / (self.gamma - 1))

        # Rebuild R± directly from c_ref — NOT from sound_speed(rho_new)
        # which would reintroduce the rho^100 round-trip error.
        Rp_ref = u_ref + 2 * c_ref / (self.gamma - 1)
        Rm_ref = u_ref - 2 * c_ref / (self.gamma - 1)

        self.u       = u_ref
        self.rho     = rho_new
        self.R_plus  = Rp_ref
        self.R_minus = Rm_ref
        self.G_plus  = Gp_ref
        self.G_minus = Gm_ref

        self.x = x_ref.copy()
        self.y = x_ref.copy()


# ── EulerSimulation ───────────────────────────────────────────────────────────

class EulerSimulation:
    """
    Runs HLLCSolver (fine grid) and LPTSolver (coarse grid) independently,
    then reports the L2 error of LPT vs HLLC at T_max.
    """

    def __init__(
        self,
        xmin, xmax, Nx,
        rho_0, u_0,
        K=1.0, gamma=2.0,
        CFL_hllc=0.9, CFL_lpt=0.99,
        T_max=0.4,
        N_factor=8,
        interp_method = 'numpy',
        monotone = False,
        kernel="se",
    ):
        self.T_max    = T_max
        self.CFL_hllc = CFL_hllc
        self.CFL_lpt  = CFL_lpt

        dx_lpt  = abs(xmax - xmin) / Nx
        x_lpt   = np.linspace(xmin + dx_lpt / 2, xmax - dx_lpt / 2, Nx)

        Nx_hllc = N_factor * Nx
        dx_hllc = abs(xmax - xmin) / Nx_hllc
        x_hllc  = np.linspace(xmin + dx_hllc / 2, xmax - dx_hllc / 2, Nx_hllc)

        self.hllc = HLLCSolver(x_hllc, rho_0, u_0, K, gamma, CFL_hllc)
        self.lpt  = LPTSolver(x_lpt,  rho_0, u_0, K, gamma, kernel=kernel)

        self.x_lpt = x_lpt

    # ── Plotting ──────────────────────────────────────────────────────────────
    
    def _plot_final(self):
        lpt  = self.lpt
        hllc = self.hllc

        rho_hllc_on_lpt = np.interp(self.x_lpt, hllc.x, hllc.rho)
        u_hllc_on_lpt   = np.interp(self.x_lpt, hllc.x, hllc.u)

        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        fig.suptitle(f"Solutions at T_max = {self.T_max:.4f}")

        axes[0].plot(self.x_lpt, lpt.rho, label="LPT",  marker="+")
        axes[0].plot(self.x_lpt, rho_hllc_on_lpt, label="HLLC (ref)", linestyle="--")
        axes[0].set_title("Density ρ")
        axes[0].set_xlabel("x")
        axes[0].legend()
        axes[0].grid(True)

        axes[1].plot(self.x_lpt, lpt.u, label="LPT",  marker="+")
        axes[1].plot(self.x_lpt, u_hllc_on_lpt, label="HLLC (ref)", linestyle="--")
        axes[1].set_title("Velocity u")
        axes[1].set_xlabel("x")
        axes[1].legend()
        axes[1].grid(True)

        plt.tight_layout()
        plt.show()
    
    # ── L2 error ──────────────────────────────────────────────────────────────

    def compute_l2_error(self):
        hllc = self.hllc
        lpt  = self.lpt

        rho_ref = np.interp(self.x_lpt, hllc.x, hllc.rho)
        u_ref   = np.interp(self.x_lpt, hllc.x, hllc.u)

        dx = self.x_lpt[1] - self.x_lpt[0]
        err_rho = np.sqrt(np.sum((lpt.rho - rho_ref)**2) * dx)
        err_u   = np.sqrt(np.sum((lpt.u   - u_ref  )**2) * dx)
        return err_rho, err_u

    # ── Main run ──────────────────────────────────────────────────────────────

    def run(self, freq=20, interp_method='numpy', monotone=False, kernel='se'):
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
        #k = 0
        
        start_lpt = time.time()
        while t < self.T_max:
            #dt_shock = self.lpt.compute_shock_dt()
            
            u = self.lpt.u
            rho = self.lpt.rho
            
            x = self.lpt.x
            
            dx = min(x[1:] - x[:-1])
            dt_ref = min(dx / max(abs((rho))), dx / max(abs((u))))
            #dt = min(dt_shock * self.CFL_lpt, self.T_max - t, dt_ref)
            dt = min(self.T_max - t, dt_ref)
            self.lpt.step(dt, interp_method=interp_method, monotone=monotone)
            t += dt
            '''
            k = 0
            #%matplotlib inline
            while k % freq == 0:
                plt.figure(figsize=(8, 6))
                plt.plot(self.lpt.x, self.lpt.rho, label="rho (LPT)", marker="+")
                plt.plot(self.lpt.x, self.lpt.u,   label="u (LPT)",   marker="o")
                plt.title(f"LPT solution at t = {t:.4f}")
                plt.xlabel("x")
                plt.legend()
                plt.grid(True)
                plt.tight_layout()
                plt.pause(0.01)   # renders the frame and keeps the loop going
                plt.close()       # optional: prevents figure accumulation
                
                k+=1
            '''
        end_lpt = time.time()
        print(f"  LPT done in {end_lpt - start_lpt:.6f} seconds with {interp_method}_{kernel}")
        
        self._plot_final()
        
        err_rho, err_u = self.compute_l2_error()
        print(f"\nL2 error at T_max = {self.T_max}:")
        print(f"  rho : {err_rho:.6e}")
        print(f"  u   : {err_u:.6e}")
        #print(f' order of convergence : , {np.log()/np.log():.6f}')
        return err_rho, err_u


# ── Initial conditions ────────────────────────────────────────────────────────

def rho_0(x):
    epsilon = 1
    return 1.5 + epsilon * np.exp(-x**2)

def u_0(x):
    epsilon = -0.2
    return epsilon * np.exp(-x**2)


'''
def rho_0(x):
    return 0.125 + (1.0 - 0.125) * 0.5 * (1 - np.tanh((x - 0.5) / 0.02))

def u_0(x):
    return np.zeros_like(x)
'''


# ── Run ───────────────────────────────────────────────────────────────────────

interp_method = 'GP'
kernel = 'se'

if __name__ == "__main__":
    sim = EulerSimulation(
        xmin=-6.0, xmax=6.0, Nx=256,
        rho_0=rho_0, u_0=u_0,
        K=1., gamma=1.4,
        CFL_hllc=0.9, CFL_lpt=1,
        T_max=2.,
        N_factor=1,
        monotone = False,
        kernel = kernel
    )
    sim.run(interp_method=interp_method, kernel = kernel)
