
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
1D_Euler_RBFGA_FD.py
--------------------
Same physics as 1D_Euler_LPTpy.py, but the semi-Lagrangian step uses
**backward characteristic tracing** instead of forward advection.

Forward scheme (1D_Euler_LPTpy.py)
    1. Advect the Lagrangian grid forward:  x_adv = x_ref + λ± · dt  (+ CK value corrections)
    2. Interpolate R± from non-uniform x_adv  →  uniform x_ref.

Backward scheme (this file)
    1. Compute departure points on the uniform grid:
           x_dep = x_ref − λ± · dt − A± · dt²/2 + B± · dt³/6   (CK shifts the departure point)
    2. Interpolate R± from uniform x_ref  →  non-uniform x_dep.

Because the SOURCE is always the fixed uniform grid, the RBF-GA stencil
geometry (uniform spacing h) is the same for every target point and every
time step.  The weight vector depends only on the scalar fractional offset
    δᵢ = round(frac_idx_i) − frac_idx_i  ∈ [−0.5, 0.5]
so it can be precomputed once as a 1-D table w(δ) and looked up in O(N·M)
per step, vs O(N·M³) for the forward scheme.
"""

import numpy as np
import matplotlib.pyplot as plt
import time

import os, sys
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(current_dir)
from rbf_ga_weights_1d import RBFGAUniformGrid, RBFGAStencil
from scipy.interpolate import make_interp_spline


# ── Equation of state ─────────────────────────────────────────────────────────

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
    """Fixed-grid Godunov solver with HLLC Riemann fluxes for 1D isentropic Euler."""

    def __init__(self, x, rho_0, u_0, K=1.0, gamma=2.0, CFL=0.9):
        self.x     = x.copy()
        self.dx    = x[1] - x[0]
        self.K     = K
        self.gamma = gamma
        self.CFL   = CFL
        rho = rho_0(x)
        u   = u_0(x)
        self.U = np.vstack((rho, rho * u))

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
        F[:, 0]  = F[:, 1]
        F[:, -1] = F[:, -2]
        self.U = U - (dt / self.dx) * (F[:, 1:] - F[:, :-1])

    def run_to(self, T_max):
        t = 0.0
        while t < T_max:
            dt = min(self.compute_dt(), T_max - t)
            self.step(dt)
            t += dt


# ── LPTSolver (backward tracing) ─────────────────────────────────────────────

class LPTSolver:
    """
    Backward semi-Lagrangian solver for 1D isentropic Euler.

    At each step the characteristic departure points are computed on the
    fixed uniform grid x_ref and R± are interpolated FROM x_ref TO x_dep.
    The grid x_ref never moves.

    GP interpolation uses a precomputed weight table w(δ) where δ is the
    fractional offset of x_dep from the nearest source grid point.  This
    reduces the per-step interpolation cost from O(N·M³) (non-uniform source)
    to O(N·M) (table lookup + stencil dot-product).

    Parameters
    ----------
    n_table : int
        Number of δ samples in the precomputed weight table.  500 is more
        than enough for linear interpolation error < 1e-10 relative to the
        stencil accuracy.
    ck : int (0, 1, or 2)
        CK correction order applied to the departure point.
        0 → pure Euler departure  (1st-order in time)
        1 → subtract A±·dt²/2    (2nd-order)
        2 → also add  B±·dt³/6   (3rd-order)
    """

    def __init__(self, x, rho_0, u_0, K=1.0,
                 gamma=2.0, kernel='se', r=4, rbfga_eps=None,
                 k_spline=5, ck=2, n_table=8_000, table_deg=3,
                 weight_mode='cheb', cheb_pad=10):
        """
        weight_mode : str
            'cheb'  — Chebyshev series per weight component, Clenshaw eval.
                      O(N·M·deg) per call, machine ε error, no Runge.  Default.
            'exact' — solve the M×M RBF-GA system at every departure point.
                      O(N·M³) per call, zero error.
            'table' — precomputed uniform table with cubic Catmull-Rom lookup.
                      O(N·M) per call.
        cheb_pad : int
            Extra grid cells near each boundary that fall back to cubic table
            in 'cheb' mode (default 10).
        """

        dx = x[1] - x[0]
        self.xmin        = x[0]  - dx / 2
        self.xmax        = x[-1] + dx / 2
        self.Nx          = len(x)
        self.dx          = dx
        self.K           = K
        self.gamma       = gamma
        self.kernel      = kernel
        self.r           = r
        self.k_spline    = k_spline
        self.table_deg   = table_deg
        self.weight_mode = weight_mode
        self._cheb_pad   = int(cheb_pad)
        self.ck          = ck   # 0 = no CK, 1 = CK1 (dt²), 2 = CK1+CK2 (dt²+dt³)

        eps       = rbfga_eps if rbfga_eps is not None else 1e-3
        self._eps = eps

        self.M     = 2 * r + 1
        self.stencil_der_add = 2
        self.M_der = self.M + 2 * self.stencil_der_add

        self.uniform_stencil_d1 = RBFGAUniformGrid(dx=dx, M=self.M,     eps=eps, deriv=1)
        self.uniform_stencil_d2 = RBFGAUniformGrid(dx=dx, M=self.M_der, eps=eps, deriv=2)

        # Fixed reference grid — never changes
        self.x_ref = np.linspace(
            self.xmin + dx / 2,
            self.xmax - dx / 2,
            self.Nx,
        )
        self.x = self.x_ref.copy()   # kept for API compatibility
        self.y = self.x_ref.copy()

        self.rho = rho_0(self.x_ref)
        self.u   = u_0(self.x_ref)

        c = sound_speed(self.rho, K, gamma)
        self.R_plus  = self.u + 2 * c / (gamma - 1)
        self.R_minus = self.u - 2 * c / (gamma - 1)

        self._stencil = RBFGAStencil(eps)
        self._build_weight_table(n_table)   # always built; used as boundary fallback
        if weight_mode == 'cheb':
            self._build_cheb_weights()
        elif weight_mode not in ('table', 'exact'):
            raise ValueError(f"weight_mode must be 'cheb', 'exact' or 'table', got {weight_mode!r}")

    # ── Weight table ──────────────────────────────────────────────────────────

    def _build_weight_table(self, n_table):
        dx  = self.dx
        r   = self.r
        M   = 2 * r + 1

        x_nodes_loc = dx * np.arange(-r, r + 1, dtype=float)
        delta_tab   = np.linspace(-0.5, 0.5, n_table)
        W_tab       = np.zeros((n_table, M))

        stencil = RBFGAStencil(self._eps)
        for i, d in enumerate(delta_tab):
            W_tab[i] = stencil.weights(x_nodes_loc, -d * dx, deriv=0)

        self._delta_tab = delta_tab
        self._W_tab     = W_tab
        self._n_table   = n_table
        self._d_lo      = delta_tab[0]
        self._d_step    = delta_tab[1] - delta_tab[0]

    # ── Chebyshev weight representation ──────────────────────────────────────

    def _build_cheb_weights(self, cheb_deg=None):
        from numpy.polynomial.chebyshev import chebfit

        dx  = self.dx
        r   = self.r
        M   = 2 * r + 1
        deg = cheb_deg if cheb_deg is not None else 2 * r + 4
        n_fit = 20

        x_nodes_loc = dx * np.arange(-r, r + 1, dtype=float)

        k_idx       = np.arange(1, n_fit + 1)
        delta_nodes = -0.5 * np.cos(np.pi * (2 * k_idx - 1) / (2 * n_fit))
        s_nodes     = 2.0 * delta_nodes

        W_samples = np.zeros((n_fit, M))
        for i, d in enumerate(delta_nodes):
            W_samples[i] = self._stencil.weights(x_nodes_loc, -d * dx, deriv=0)

        self._cheb_coeffs = np.array([
            chebfit(s_nodes, W_samples[:, m], deg)
            for m in range(M)
        ])   # (M, deg+1)

    @staticmethod
    def _chebval_batch(s, coeffs):
        """Evaluate M Chebyshev series at N points via vectorised Clenshaw recursion."""
        M, K = coeffs.shape
        N    = len(s)
        s2   = 2.0 * s
        b1   = np.zeros((N, M))
        b2   = np.zeros((N, M))
        for j in range(K - 1, 0, -1):
            b_new = s2[:, None] * b1 - b2 + coeffs[:, j][None, :]
            b2    = b1
            b1    = b_new
        return s[:, None] * b1 - b2 + coeffs[:, 0][None, :]

    # ── Fast backward GP interpolation ───────────────────────────────────────

    def _interp_table(self, f_src, x_dep):
        """Interpolate f_src (on uniform self.x_ref) at non-uniform x_dep."""
        r  = self.r
        N  = self.Nx
        x0 = self.x_ref[0]

        xmin      = x0 - 0.5 * self.dx
        xmax      = self.x_ref[-1] + 0.5 * self.dx
        out_left  = x_dep < xmin
        out_right = x_dep > xmax
        if out_left.any() or out_right.any():
            result               = np.empty(len(x_dep))
            result[out_left]     = f_src[0]
            result[out_right]    = f_src[-1]
            mask                 = ~(out_left | out_right)
            if mask.any():
                result[mask] = self._interp_table(f_src, x_dep[mask])
            return result

        frac_idx = (x_dep - x0) / self.dx
        j        = np.round(frac_idx).astype(int)
        j        = np.clip(j, r, N - r - 1)

        proto = np.arange(-r, r + 1, dtype=int)
        si    = np.clip(j[:, None] + proto[None, :], 0, N - 1)

        delta = j.astype(float) - frac_idx   # δ ∈ [−0.5, 0.5]

        if self.weight_mode == 'cheb':
            pad      = self._cheb_pad
            W        = np.empty((len(delta), self.M))
            interior = (np.abs(delta) <= 0.5) & (j > r + pad) & (j < N - r - 1 - pad)
            if interior.any():
                s = 2.0 * delta[interior]
                W[interior] = self._chebval_batch(s, self._cheb_coeffs)
            bnd = ~interior
            if bnd.any():
                d_bnd = np.clip(delta[bnd], -0.5, 0.5)
                t_bnd = (d_bnd - self._d_lo) / self._d_step
                t_bnd = np.clip(t_bnd, 1.0, self._n_table - 3.0)
                idx   = t_bnd.astype(int)
                f_t   = (t_bnd - idx)[:, None]
                f2, f3 = f_t * f_t, f_t * f_t * f_t
                c0 = -0.5*f3 + 1.0*f2 - 0.5*f_t
                c1 =  1.5*f3 - 2.5*f2 + 1.0
                c2 = -1.5*f3 + 2.0*f2 + 0.5*f_t
                c3 =  0.5*f3 - 0.5*f2
                W[bnd] = (c0 * self._W_tab[idx - 1] + c1 * self._W_tab[idx]
                        + c2 * self._W_tab[idx + 1] + c3 * self._W_tab[idx + 2])

        elif self.weight_mode == 'exact':
            W = self._stencil.weights_batch(self.x_ref[si], x_dep)

        else:  # 'table'
            t   = (delta - self._d_lo) / self._d_step
            t   = np.clip(t, 1.0, self._n_table - 3.0)
            idx = t.astype(int)
            f_t = (t - idx)[:, None]
            f2, f3 = f_t * f_t, f_t * f_t * f_t
            c0 = -0.5*f3 + 1.0*f2 - 0.5*f_t
            c1 =  1.5*f3 - 2.5*f2 + 1.0
            c2 = -1.5*f3 + 2.0*f2 + 0.5*f_t
            c3 =  0.5*f3 - 0.5*f2
            W  = (c0 * self._W_tab[idx - 1] + c1 * self._W_tab[idx]
                + c2 * self._W_tab[idx + 1] + c3 * self._W_tab[idx + 2])

        return (W * f_src[si]).sum(axis=1)

    # ── Helpers ───────────────────────────────────────────────────────────────

    def _char_speeds(self):
        c     = sound_speed(self.rho, self.K, self.gamma)
        lam_p = self.u + c
        lam_m = self.u - c
        return lam_p, lam_m

    def compute_shock_dt(self):
        lam_p, lam_m = self._char_speeds()
        dx = self.dx
        dp = dx / np.abs(lam_p[1:] - lam_p[:-1] + 1e-14)
        dm = dx / np.abs(lam_m[1:] - lam_m[:-1] + 1e-14)
        return min(float(np.min(dp)), float(np.min(dm)))

    # ── Step ─────────────────────────────────────────────────────────────────

    def step(self, dt, interp_method='numpy'):
        lam_p, lam_m = self._char_speeds()

        x_ref   = self.x_ref
        R_plus  = self.R_plus
        R_minus = self.R_minus
        r       = self.r

        M               = self.M
        M_der           = self.M_der

        g = self.gamma
        c = (g - 1) / 4 * (R_plus - R_minus)

        # ── Strided RBF-GA derivatives on the fixed uniform grid ─────────────
        w1 = self.uniform_stencil_d1.w
        w2 = self.uniform_stencil_d2.w

        Rp_pad  = np.pad(R_plus,  r, mode='edge')
        Rm_pad  = np.pad(R_minus, r, mode='edge')
        wins_Rp = np.lib.stride_tricks.sliding_window_view(Rp_pad, M)
        wins_Rm = np.lib.stride_tricks.sliding_window_view(Rm_pad, M)
        dR_plus  = wins_Rp @ w1
        dR_minus = wins_Rm @ w1

        Rp_pad  = np.pad(R_plus,  r + self.stencil_der_add, mode='edge')
        Rm_pad  = np.pad(R_minus, r + self.stencil_der_add, mode='edge')
        wins_Rp = np.lib.stride_tricks.sliding_window_view(Rp_pad, M_der)
        wins_Rm = np.lib.stride_tricks.sliding_window_view(Rm_pad, M_der)
        
        ddR_plus  = wins_Rp @ w2
        ddR_minus = wins_Rm @ w2

        # ── Auxiliary spatial quantities on uniform grid ──────────────────────
        c_x     = (g - 1) / 4 * (dR_plus - dR_minus)
        lam_p_x = (g + 1) / 4 * dR_plus + (3 - g) / 4 * dR_minus
        lam_m_x = (3 - g) / 4 * dR_plus + (g + 1) / 4 * dR_minus
        
        # ── Backward characteristic tracing ───────────────────────────────────
        # Fixed-point iteration on:
        #   x_dep = x_ref − λ(x_dep)·dt − A(x_dep)·dt²/2 + B(x_dep)·dt³/6
        #
        # Iteration scheme and departure accuracy:
        #   x⁽⁰⁾ = x_ref − λ(x_ref)·dt                              O(dt²)  CK0
        #   x⁽¹⁾ = x_ref − λ(x⁽⁰⁾)·dt − A(x⁽⁰⁾)·dt²/2             O(dt³)  CK1
        #   x⁽²⁾ = x_ref − λ(x⁽¹⁾)·dt − A(x⁽¹⁾)·dt²/2
        #                              + B(x⁽⁰⁾)·dt³/6               O(dt⁴)  CK2
        #
        # B is kept at x⁽⁰⁾: moving it to x⁽¹⁾ only saves O(dt⁵).
        # Evaluating λ and A at x⁽⁰⁾ in the CK2 step (old scheme) leaves an
        # O(dt³) residual λ_x·(A/2−λ_t)·dt³, giving only 2nd-order global — same
        # as CK1.  The second iteration (λ,A at x⁽¹⁾) cancels that term → 3rd order.
        x_dep_p0 = x_ref - lam_p * dt   # x⁽⁰⁾ Euler departure, + char
        x_dep_m0 = x_ref - lam_m * dt   # x⁽⁰⁾ Euler departure, - char
        x_dep_p  = x_dep_p0             # default CK0
        x_dep_m  = x_dep_m0

        if self.ck >= 1:
            A_p_grid =  (3 - g) / 2 * c * dR_minus
            A_m_grid = -(3 - g) / 2 * c * dR_plus

            # Evaluate λ, A at x⁽⁰⁾ → x⁽¹⁾  (O(dt³) departure accuracy)
            lam_p_d = self._interp_table(lam_p,    x_dep_p0)
            A_p_dep = self._interp_table(A_p_grid, x_dep_p0)

            lam_m_d = self._interp_table(lam_m,    x_dep_m0)
            A_m_dep = self._interp_table(A_m_grid, x_dep_m0)

            x_dep_p1 = x_ref - lam_p_d * dt - A_p_dep * (dt**2 / 2)   # x⁽¹⁾
            x_dep_m1 = x_ref - lam_m_d * dt - A_m_dep * (dt**2 / 2)
            x_dep_p  = x_dep_p1
            x_dep_m  = x_dep_m1
        
        if self.ck == 2:
            # Precompute B± grids
            
            c_t = (g - 1) / 4 * (-lam_p * dR_plus + lam_m * dR_minus)

            dtRm_grid = -lam_m_x * dR_minus - lam_m * ddR_minus
            dtRp_grid = -lam_p_x * dR_plus  - lam_p * ddR_plus
            
            Ap_t_grid = (3 - g) / 2 * (c_t * dR_minus + c * dtRm_grid)
            Ap_x_grid = (3 - g) / 2 * (c_x * dR_minus + c * ddR_minus)
            Am_t_grid = -(3 - g) / 2 * (c_t * dR_plus  + c * dtRp_grid)
            Am_x_grid = -(3 - g) / 2 * (c_x * dR_plus  + c * ddR_plus)
            
            # B evaluated at x⁽⁰⁾ (moving to x⁽¹⁾ only saves O(dt⁵))
            Ap_t_p  = self._interp_table(Ap_t_grid, x_dep_p0)
            Ap_x_p  = self._interp_table(Ap_x_grid, x_dep_p0)
            B_p_dep = -Ap_t_p - lam_p_d * Ap_x_p

            Am_t_m  = self._interp_table(Am_t_grid, x_dep_m0)
            Am_x_m  = self._interp_table(Am_x_grid, x_dep_m0)
            B_m_dep = -Am_t_m - lam_m_d * Am_x_m
            
            # Second iteration: re-evaluate λ, A at x⁽¹⁾ → x⁽²⁾  (O(dt⁴) departure accuracy)
            lam_p_d2 = self._interp_table(lam_p,    x_dep_p1)
            A_p_dep2 = self._interp_table(A_p_grid, x_dep_p1)

            lam_m_d2 = self._interp_table(lam_m,    x_dep_m1)
            A_m_dep2 = self._interp_table(A_m_grid, x_dep_m1)

            #print(max(abs(A_p_dep2)))
            x_dep_p = x_ref - lam_p_d2 * dt - A_p_dep2 * (dt**2 / 2) + B_p_dep * (dt**3 / 6)
            x_dep_m = x_ref - lam_m_d2 * dt - A_m_dep2 * (dt**2 / 2) + B_m_dep * (dt**3 / 6)
        
        # ── Interpolate R± from uniform x_ref to departure points ────────────
        if interp_method == 'numpy':
            Rp_new = np.interp(x_dep_p, x_ref, R_plus)
            Rm_new = np.interp(x_dep_m, x_ref, R_minus)

        elif interp_method == 'cubic':
            Rp_new = make_interp_spline(x_ref, R_plus,  k=self.k_spline)(x_dep_p)
            Rm_new = make_interp_spline(x_ref, R_minus, k=self.k_spline)(x_dep_m)

        elif interp_method == 'GP':
            Rp_new = self._interp_table(R_plus,  x_dep_p)
            Rm_new = self._interp_table(R_minus, x_dep_m)

        # ── Recover primitive variables ───────────────────────────────────────
        u_new   = 0.5 * (Rp_new + Rm_new)
        c_new   = 0.25 * (self.gamma - 1) * (Rp_new - Rm_new)
        rho_new = (c_new**2 / (self.K * self.gamma))**(1 / (self.gamma - 1))

        # Rebuild R± from c_new directly to avoid rho^{1/(γ-1)} round-trip error
        self.R_plus  = u_new + 2 * c_new / (self.gamma - 1)
        self.R_minus = u_new - 2 * c_new / (self.gamma - 1)
        self.u       = u_new
        self.rho     = rho_new
        # x_ref, self.x, self.y are unchanged (always the fixed uniform grid)


# ── EulerSimulation ───────────────────────────────────────────────────────────

class EulerSimulation:
    """Runs HLLCSolver (fine grid) and LPTSolver (coarse grid) and reports L2 error."""

    def __init__(
        self,
        xmin, xmax, Nx,
        rho_0, u_0,
        K=1.0, gamma=2.0,
        CFL_hllc=0.9,
        T_max=0.4,
        N_factor=8,
        kernel="se",
        cfl=1,
    ):
        self.T_max    = T_max
        self.CFL_hllc = CFL_hllc
        self.cfl      = cfl

        dx_lpt  = abs(xmax - xmin) / Nx
        x_lpt   = np.linspace(xmin + dx_lpt / 2, xmax - dx_lpt / 2, Nx)

        Nx_hllc = N_factor * Nx
        dx_hllc = abs(xmax - xmin) / Nx_hllc
        x_hllc  = np.linspace(xmin + dx_hllc / 2, xmax - dx_hllc / 2, Nx_hllc)

        self.hllc  = HLLCSolver(x_hllc, rho_0, u_0, K, gamma, CFL_hllc)
        self.lpt   = LPTSolver(x_lpt,  rho_0, u_0, K, gamma, kernel=kernel)
        self.x_lpt = x_lpt

    def _plot_final(self):
        lpt  = self.lpt
        hllc = self.hllc
        rho_hllc_on_lpt = np.interp(self.x_lpt, hllc.x, hllc.rho)
        u_hllc_on_lpt   = np.interp(self.x_lpt, hllc.x, hllc.u)
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        fig.suptitle(f"Solutions at T_max = {self.T_max:.4f}")
        axes[0].plot(self.x_lpt, lpt.rho, label="LPT (bwd)",  marker="+")
        axes[0].plot(self.x_lpt, rho_hllc_on_lpt, label="HLLC (ref)", linestyle="--")
        axes[0].set_title("Density ρ"); axes[0].set_xlabel("x")
        axes[0].legend(); axes[0].grid(True)
        axes[1].plot(self.x_lpt, lpt.u, label="LPT (bwd)",  marker="+")
        axes[1].plot(self.x_lpt, u_hllc_on_lpt, label="HLLC (ref)", linestyle="--")
        axes[1].set_title("Velocity u"); axes[1].set_xlabel("x")
        axes[1].legend(); axes[1].grid(True)
        plt.tight_layout(); plt.show()

    def compute_l2_error(self):
        hllc = self.hllc
        lpt  = self.lpt
        rho_ref = np.interp(self.x_lpt, hllc.x, hllc.rho)
        u_ref   = np.interp(self.x_lpt, hllc.x, hllc.u)
        dx = self.x_lpt[1] - self.x_lpt[0]
        err_rho = np.sqrt(np.sum((lpt.rho - rho_ref)**2) * dx)
        err_u   = np.sqrt(np.sum((lpt.u   - u_ref  )**2) * dx)
        return err_rho, err_u

    def run(self, interp_method='numpy', kernel='se'):
        print("Running HLLC (reference) ...")
        t0 = time.time()
        self.hllc.run_to(self.T_max)
        print(f"  HLLC done in {time.time()-t0:.4f} s")

        print("Running LPT (backward tracing) ...")
        t = 0.0
        t0 = time.time()
        while t < self.T_max:
            dx    = self.lpt.dx
            u     = self.lpt.u
            rho   = self.lpt.rho
            lam   = np.abs(u) + sound_speed(rho, self.lpt.K, self.lpt.gamma)
            dt_cfl = self.cfl * dx / float(np.max(lam))
            dt    = min(self.T_max - t, dt_cfl)
            self.lpt.step(dt, interp_method=interp_method)
            t += dt
        print(f"  LPT done in {time.time()-t0:.4f} s with {interp_method}")

        self._plot_final()
        err_rho, err_u = self.compute_l2_error()
        print(f"\nL2 error at T_max = {self.T_max}:")
        print(f"  rho : {err_rho:.6e}")
        print(f"  u   : {err_u:.6e}")
        return err_rho, err_u


# ── Initial conditions ────────────────────────────────────────────────────────

def rho_0(x):
    return 1.5 + np.exp(-x**2 / 4)

def u_0(x):
    return -0.2 * np.exp(-x**2 / 4)


# ── Run ───────────────────────────────────────────────────────────────────────

interp_method = 'GP'

if __name__ == "__main__":
    sim = EulerSimulation(
        xmin=-10.0, xmax=10.0, Nx=256,
        rho_0=rho_0, u_0=u_0,
        K=1., gamma=1.4,
        CFL_hllc=0.1,
        T_max=1.,
        N_factor=10,
        kernel='se',
        cfl=1,
    )
    sim.run(interp_method=interp_method, kernel='se')
