
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
1D_Euler_RBFGA_ImplicitSL.py
-----------------------------
Implicit semi-Lagrangian (ISL) solver for 1D isentropic Euler in Riemann-invariant form.

RBF-GA (or barycentric Lagrange) interpolation is used to evaluate R± at
departure points.  Two temporal schemes are available via the `ck` flag:

ck=0  — implicit Euler departure.  Newton on (R+ⁿ⁺¹, R-ⁿ⁺¹):
           x_D± = x − λ±(Rⁿ⁺¹(x)) · dt
           R±ⁿ⁺¹(x) = R±ⁿ(x_D±)
        Per-step departure error O(dt²) → global O(dt).

ck=1  — implicit trapezoidal departure.  Newton on departure points (x_D+, x_D-):
           x_D± = x − dt/2 · (λ±ⁿ(x_D±)  +  λ±ⁿ⁺¹(x))
        where the arrival speed uses Rⁿ⁺¹(x) = (R+ⁿ(x_D+), R-ⁿ(x_D-)).
        Expanding into the coupled system solved by Newton_trap:
           F+ = x_D+ − x + dt/2·(2a·R+ⁿ(x_D+) + b·R-ⁿ(x_D+) + b·R-ⁿ(x_D-)) = 0
           F- = x_D- − x + dt/2·(b·R+ⁿ(x_D+) + 2a·R-ⁿ(x_D-) + b·R+ⁿ(x_D-)) = 0
        with a=(γ+1)/4, b=(3−γ)/4.
        Per-step departure error O(dt³) → global O(dt²).
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


# ── LPTSolver (implicit semi-Lagrangian) ─────────────────────────────────────

class LPTSolver:
    """
    Implicit semi-Lagrangian solver for 1D isentropic Euler in Riemann-invariant form.

    The grid x_ref is fixed and uniform.  RBF-GA (or barycentric Lagrange)
    interpolation weights are precomputed once as a 1-D table w(δ),
    δ ∈ [−0.5, 0.5], and looked up in O(N·M) per interpolation call.

    Two temporal schemes:
        ck=0 : implicit Euler departure — Newton on (R+ⁿ⁺¹, R-ⁿ⁺¹). O(dt) global.
        ck=1 : implicit trapezoidal departure — Newton on (x_D+, x_D-). O(dt²) global.

    Parameters
    ----------
    n_table : int
        Number of δ samples in the precomputed weight table.
    """

    def __init__(self, x, rho_0, u_0, K=1.0,
                 gamma=2.0, kernel='se', r=2, rbfga_eps=None,
                 k_spline=5, ck=1, n_table=8_000, table_deg=3,
                 weight_mode='cheb', cheb_pad=0,
                 newton_iter=1, interp_backend='gp',
                 adaptive_newton=False, adaptive_newton_tol=1e-6,
                 monotone_limiter=True, limiter_threshold=0.5):
        """
        weight_mode : str
            'poly'  — monomial fit at uniform nodes; O(N·M·deg) per call,
                      machine-ε accuracy.  Default and recommended.
            'cheb'  — Chebyshev series fit (Clenshaw eval); kept for
                      comparison with 'poly'.  O(N·M·deg) per call.
            'exact' — solve the M×M RBF-GA system at every departure point.
                      O(N·M³) per call, zero error.  Use as reference only.
            'table' — precomputed uniform table with linear or Catmull-Rom cubic
                      lookup (controlled by table_deg).  O(N·M) per call.
            table_deg : 1 → linear, 3 → cubic Catmull-Rom  (only used when
                        weight_mode='table')
        cheb_pad : int
            Extra boundary cells that fall back to table in 'cheb' mode.
            Default 0 (only j-clipped cells fall back).
        cross_corr : bool
            If True, add a product correction after the ½-step targeting the
            cross-derivative coupling:
                δR± = cross_beta · (R̃⁺ − R⁺) · (R̃⁻ − R⁻) · dt
            O(dt³) per step, so does not change convergence order but reduces
            the coupling constant C.  Default False.
        cross_beta : float
            Coefficient for the cross-coupling correction.  Default 20.
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
        self.weight_mode    = weight_mode
        self._cheb_pad      = int(cheb_pad)
        self.interp_backend = interp_backend   # 'gp' or 'lagrange'

        eps       = rbfga_eps if rbfga_eps is not None else dx
        self._eps = eps
        self.M    = 2 * r + 1
        self.r_d  = r + 1            # derivative stencil half-width: one wider
        self.M_d  = 2 * (r + 1) + 1 # gives O(h^{M}) derivative vs O(h^{M}) interp
        self.ck                  = ck   # 0 → implicit Euler, 1 → implicit trapezoidal
        self.newton_iter         = newton_iter
        self.adaptive_newton     = adaptive_newton
        self.adaptive_newton_tol = adaptive_newton_tol
        self.monotone_limiter    = monotone_limiter
        self.limiter_threshold   = limiter_threshold
        self._shock_indicator       = np.zeros(len(x))
        self._newton_residual       = np.zeros(len(x))
        self._newton_residual_iter1 = np.zeros(len(x))

        # Fixed reference grid — never changes
        self.x_ref = np.linspace(
            self.xmin + dx / 2,
            self.xmax - dx / 2,
            self.Nx,
        )
        self.x = self.x_ref.copy()
        self.y = self.x_ref.copy()

        self.rho = rho_0(self.x_ref)
        self.u   = u_0(self.x_ref)

        if float(gamma) == 1.0:
            a = np.sqrt(K)
            self.R_plus  = self.u + a * np.log(self.rho)
            self.R_minus = self.u - a * np.log(self.rho)
        else:
            c = sound_speed(self.rho, K, gamma)
            self.R_plus  = self.u + 2 * c / (gamma - 1)
            self.R_minus = self.u - 2 * c / (gamma - 1)

        self._stencil = RBFGAStencil(eps)
        self._build_weight_table(n_table)   # always built; used as boundary fallback
        if weight_mode == 'poly':
            self._build_poly_weights()
        elif weight_mode == 'cheb':
            self._build_cheb_weights()
        elif weight_mode not in ('table', 'exact'):
            raise ValueError(f"weight_mode must be 'poly', 'exact' or 'table', got {weight_mode!r}")

    # ── Barycentric Lagrange weights (machine-ε, no RBF-GA ill-conditioning) ──

    @staticmethod
    def _lagrange_weights(r, delta):
        """
        Barycentric Lagrange interpolation weights for the prototype uniform
        stencil with nodes at {-r, ..., r} evaluated at fractional offset delta.
        Returns (M,) weights w s.t.  w @ f_stencil ≈ f(delta).
        Exact for polynomials of degree ≤ 2r; machine-ε accurate for any dx.
        """
        M     = 2 * r + 1
        nodes = np.arange(M, dtype=float) - r          # [-r, ..., r]
        diff  = delta - nodes                           # δ − n_k

        at_node = np.abs(diff) < 1e-14
        if np.any(at_node):
            w = np.zeros(M)
            w[np.where(at_node)[0][0]] = 1.0
            return w

        from math import comb
        bary = np.array([(-1) ** k * comb(M - 1, k) for k in range(M)], dtype=float)
        val  = bary / diff
        return val / val.sum()

    @staticmethod
    def _lagrange_deriv_weights(r, delta, dx):
        """
        Weights for d/dx of the Lagrange interpolant at fractional offset delta
        for the prototype stencil {-r, ..., r} with physical spacing dx.
        Returns (M,) weights w s.t.  w @ f_stencil ≈ f'(delta·dx).
        """
        M     = 2 * r + 1
        nodes = np.arange(M, dtype=float) - r
        diff  = delta - nodes

        from math import comb
        bary = np.array([(-1) ** k * comb(M - 1, k) for k in range(M)], dtype=float)

        at_node = np.abs(diff) < 1e-14
        if np.any(at_node):
            j = np.where(at_node)[0][0]
            k = np.arange(M, dtype=float)
            w = np.where(k != j,
                         (bary / bary[j]) / (nodes[j] - nodes),
                         0.0)
            w[j] = -w.sum()
            return w / dx

        val    = bary / diff
        L      = val / val.sum()         # interpolation weights
        S1     = (1.0 / diff).sum()      # Σ_k 1/(δ−n_k)
        return L * (S1 - 1.0 / diff) / dx

    # ── Weight table ──────────────────────────────────────────────────────────

    def _build_weight_table(self, n_table):
        dx    = self.dx
        r     = self.r
        eps   = self._eps
        M     = 2 * r + 1

        r_d  = self.r_d
        M_d  = self.M_d

        stencil        = RBFGAStencil(eps)
        x_nodes_loc    = dx * np.arange(-r,   r   + 1, dtype=float)
        x_nodes_loc_d1 = dx * np.arange(-r_d, r_d + 1, dtype=float)
        delta_tab   = np.linspace(-0.5, 0.5, n_table)
        W_tab    = np.zeros((n_table, M))
        W_tab_d1 = np.zeros((n_table, M_d))

        for i, d in enumerate(delta_tab):
            if self.interp_backend == 'gp':
                W_tab[i]    = stencil.weights(x_nodes_loc,    -d * dx, deriv=0)
                W_tab_d1[i] = stencil.weights(x_nodes_loc_d1, -d * dx, deriv=1)
            else:  # 'lagrange'
                W_tab[i]    = self._lagrange_weights(r,   -d)
                W_tab_d1[i] = self._lagrange_deriv_weights(r_d, -d, dx)

        self._delta_tab = delta_tab
        self._W_tab     = W_tab
        self._W_tab_d1  = W_tab_d1
        self._n_table   = n_table
        self._d_lo      = delta_tab[0]
        self._d_step    = delta_tab[1] - delta_tab[0]

    # ── Polynomial weight representation ─────────────────────

    def _build_poly_weights(self):
        """
        Fit w_k(δ) as a polynomial in δ ∈ [−0.5, 0.5] for each of the M weight
        components.  Samples at uniform nodes, least-squares fit.
        Degree 3r+3 (not 2r+2) — the RBF-GA weights are analytic but not
        polynomial, and at fine resolution each per-call fit residual accumulates
        over O(1/dx) steps; higher degree pushes that residual to machine ε.
        """
        dx  = self.dx
        r   = self.r;   r_d = self.r_d
        M   = self.M;   M_d = self.M_d
        deg   = min(3 * r   + 3, 16)
        deg_d = min(3 * r_d + 3, 16)
        K   = 2 * max(deg, deg_d) + 4

        x_nodes_loc    = dx * np.arange(-r,   r   + 1, dtype=float)
        x_nodes_loc_d1 = dx * np.arange(-r_d, r_d + 1, dtype=float)
        delta_nodes = np.linspace(-0.5, 0.5, K)

        W_samples    = np.zeros((K, M))
        W_samples_d1 = np.zeros((K, M_d))
        for i, d in enumerate(delta_nodes):
            if self.interp_backend == 'gp':
                W_samples[i]    = self._stencil.weights(x_nodes_loc,    -d * dx, deriv=0)
                W_samples_d1[i] = self._stencil.weights(x_nodes_loc_d1, -d * dx, deriv=1)
            else:  # 'lagrange'
                W_samples[i]    = self._lagrange_weights(r,   -d)
                W_samples_d1[i] = self._lagrange_deriv_weights(r_d, -d, dx)

        self._poly_coeffs    = np.array([
            np.polyfit(delta_nodes, W_samples[:, m], deg)
            for m in range(M)
        ])                                   # shape (M, deg+1)
        self._poly_coeffs_d1 = np.array([
            np.polyfit(delta_nodes, W_samples_d1[:, m], deg_d)
            for m in range(M_d)
        ])

    # ── Chebyshev weight representation ────────────────

    def _build_cheb_weights(self):
        from numpy.polynomial.chebyshev import chebfit
        dx    = self.dx
        r     = self.r;   r_d = self.r_d
        M     = self.M;   M_d = self.M_d
        deg   = min(3 * r   + 3, 20)
        deg_d = min(3 * r_d + 3, 20)
        n_fit = 2 * max(deg, deg_d) + 4
        x_nodes_loc    = dx * np.arange(-r,   r   + 1, dtype=float)
        x_nodes_loc_d1 = dx * np.arange(-r_d, r_d + 1, dtype=float)
        k_idx       = np.arange(1, n_fit + 1)
        delta_nodes = -0.5 * np.cos(np.pi * (2 * k_idx - 1) / (2 * n_fit))
        s_nodes      = 2.0 * delta_nodes
        W_samples    = np.zeros((n_fit, M))
        W_samples_d1 = np.zeros((n_fit, M_d))
        for i, d in enumerate(delta_nodes):
            if self.interp_backend == 'gp':
                W_samples[i]    = self._stencil.weights(x_nodes_loc,    -d * dx, deriv=0)
                W_samples_d1[i] = self._stencil.weights(x_nodes_loc_d1, -d * dx, deriv=1)
            else:  # 'lagrange'
                W_samples[i]    = self._lagrange_weights(r,   -d)
                W_samples_d1[i] = self._lagrange_deriv_weights(r_d, -d, dx)
        self._cheb_coeffs    = np.array([
            chebfit(s_nodes, W_samples[:, m], deg) for m in range(M)
        ])   # (M, deg+1)
        self._cheb_coeffs_d1 = np.array([
            chebfit(s_nodes, W_samples_d1[:, m], deg_d) for m in range(M_d)
        ])

    @staticmethod
    def _chebval_batch(s, coeffs):
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

    # ── Interpolation (deriv=0) or derivative (deriv=1) ─────────────────────

    def _interp_table(self, f_src, x_dep, deriv=0):
        """
        Interpolate f_src (on uniform x_ref) at non-uniform x_dep.
        deriv=0 : value interpolation  (w @ f ≈ f(x_dep))
        deriv=1 : derivative via GP    (w @ f ≈ f'(x_dep)), uses precomputed _d1 weights
        """
        r  = self.r  if deriv == 0 else self.r_d
        M  = self.M  if deriv == 0 else self.M_d
        N  = self.Nx
        x0 = self.x_ref[0]

        # Constant-BC: backward leg can push boundary cells outside the physical domain.
        # Return the nearest boundary cell value (deriv=0) or 0 (deriv=1).
        xmin      = x0 - 0.5 * self.dx
        xmax      = self.x_ref[-1] + 0.5 * self.dx
        out_left  = x_dep < xmin
        out_right = x_dep > xmax
        if out_left.any() or out_right.any():
            result               = np.empty(len(x_dep))
            result[out_left]     = f_src[0]  if deriv == 0 else 0.0
            result[out_right]    = f_src[-1] if deriv == 0 else 0.0
            mask                 = ~(out_left | out_right)
            if mask.any():
                result[mask] = self._interp_table(f_src, x_dep[mask], deriv)
            return result

        frac_idx = (x_dep - x0) / self.dx
        j = np.round(frac_idx).astype(int)
        j = np.clip(j, r, N - r - 1)

        proto = np.arange(-r, r + 1, dtype=int)
        si    = np.clip(j[:, None] + proto[None, :], 0, N - 1)

        delta = j.astype(float) - frac_idx   # δ ∈ [−0.5, 0.5]

        W_tab    = self._W_tab    if deriv == 0 else self._W_tab_d1

        if self.weight_mode == 'cheb':
            cheb_c   = self._cheb_coeffs    if deriv == 0 else self._cheb_coeffs_d1
            pad      = self._cheb_pad
            W        = np.empty((len(delta), M))
            interior = (np.abs(delta) <= 0.5) & (j > r + pad) & (j < N - r - 1 - pad)
            if interior.any():
                s = 2.0 * delta[interior]
                W[interior] = self._chebval_batch(s, cheb_c)
            bnd = ~interior
            if bnd.any():
                d_bnd = np.clip(delta[bnd], -0.5, 0.5)
                t_bnd = np.clip((d_bnd - self._d_lo) / self._d_step, 1.0, self._n_table - 3.0)
                idx   = t_bnd.astype(int)
                f_t   = (t_bnd - idx)[:, None]
                f2, f3 = f_t * f_t, f_t * f_t * f_t
                c0 = -0.5*f3 + 1.0*f2 - 0.5*f_t
                c1 =  1.5*f3 - 2.5*f2 + 1.0
                c2 = -1.5*f3 + 2.0*f2 + 0.5*f_t
                c3 =  0.5*f3 - 0.5*f2
                W[bnd] = (c0 * W_tab[idx - 1] + c1 * W_tab[idx]
                        + c2 * W_tab[idx + 1] + c3 * W_tab[idx + 2])

        elif self.weight_mode == 'exact':
            if self.interp_backend == 'gp':
                # O(N·M³): solve the RBF-GA system exactly at each departure point
                x_stencil = self.x_ref[si]
                W = self._stencil.weights_batch(x_stencil, x_dep, deriv=deriv)
            else:  # 'lagrange'
                r_w = self.r if deriv == 0 else self.r_d
                W   = np.stack([
                    self._lagrange_weights(r_w, -delta[i]) if deriv == 0
                    else self._lagrange_deriv_weights(r_w, -delta[i], self.dx)
                    for i in range(len(delta))
                ])

        elif self.weight_mode == 'poly':
            poly_c = self._poly_coeffs if deriv == 0 else self._poly_coeffs_d1

            W = np.zeros((len(delta), M))

            # Within ?? stencil positions of either boundary → table
            # (Chebyshev node spacing causes Runge-like drift very near boundaries)
            # five nodes spacing near boundaries, need more
            near_boundary = (j <= r) | (j >= N - r - 1)

            interior = (np.abs(delta) <= 0.5) & ~near_boundary

            if interior.any():
                d_int = delta[interior]
                for m in range(M):
                    W[interior, m] = np.polyval(poly_c[m], d_int)

            # Near-boundary cells with |δ|≤0.5 → table lookup
            tbl = near_boundary & (np.abs(delta) <= 0.5)
            if tbl.any():
                t_tbl = (delta[tbl] - self._d_lo) / self._d_step
                if self.table_deg == 1:
                    t_tbl = np.clip(t_tbl, 0.0, self._n_table - 1 - 1e-10)
                    idx_t = t_tbl.astype(int)
                    f_t   = (t_tbl - idx_t)[:, None]
                    W[tbl] = (1.0 - f_t) * W_tab[idx_t] + f_t * W_tab[idx_t + 1]
                else:  # cubic Catmull-Rom
                    t_tbl = np.clip(t_tbl, 1.0, self._n_table - 3.0)
                    idx_t = t_tbl.astype(int)
                    f_t   = (t_tbl - idx_t)[:, None]
                    f2t, f3t = f_t * f_t, f_t * f_t * f_t
                    c0 = -0.5*f3t + 1.0*f2t - 0.5*f_t
                    c1 =  1.5*f3t - 2.5*f2t + 1.0
                    c2 = -1.5*f3t + 2.0*f2t + 0.5*f_t
                    c3 =  0.5*f3t - 0.5*f2t
                    W[tbl] = (c0 * W_tab[idx_t - 1]
                            + c1 * W_tab[idx_t]
                            + c2 * W_tab[idx_t + 1]
                            + c3 * W_tab[idx_t + 2])

            # Clipped departure (|δ|>0.5) → exact
            clipped = ~interior & ~tbl
            if clipped.any():
                x_stencil = self.x_ref[si[clipped]]
                W[clipped] = self._stencil.weights_batch(
                    x_stencil, x_dep[clipped], deriv=deriv)

        else:  # 'table'
            t = (delta - self._d_lo) / self._d_step
            if self.table_deg == 1:
                t   = np.clip(t, 0.0, self._n_table - 1 - 1e-10)
                idx = t.astype(int)
                f   = (t - idx)[:, None]
                W   = (1.0 - f) * W_tab[idx] + f * W_tab[idx + 1]
            elif self.table_deg == 3:
                # Catmull-Rom cubic: 4 nodes at idx-1, idx, idx+1, idx+2
                t   = np.clip(t, 1.0, self._n_table - 3.0)
                idx = t.astype(int)
                f   = (t - idx)[:, None]
                f2, f3 = f * f, f * f * f
                c0 = -0.5*f3 + 1.0*f2 - 0.5*f
                c1 =  1.5*f3 - 2.5*f2 + 1.0
                c2 = -1.5*f3 + 2.0*f2 + 0.5*f
                c3 =  0.5*f3 - 0.5*f2
                W  = (c0 * W_tab[idx - 1]
                    + c1 * W_tab[idx]
                    + c2 * W_tab[idx + 1]
                    + c3 * W_tab[idx + 2])
            else:
                raise ValueError(f"table_deg must be 1 or 3, got {self.table_deg}")

        return (W * f_src[si]).sum(axis=1)

    # ── Helpers ───────────────────────────────────────────────────────────────

    def _char_speeds(self, R_plus, R_minus):
        u = 0.5 * (R_plus + R_minus)
        if float(self.gamma) == 1.0:
            c = np.full_like(u, np.sqrt(self.K))
        else:
            c = (self.gamma - 1) / 4 * (R_plus - R_minus)
        lam_p = u + c
        lam_m = u - c
        return lam_p, lam_m

    def compute_shock_dt(self):
        lam_p, lam_m = self._char_speeds(self.R_plus, self.R_minus)
        dx = self.dx
        # Only converging pairs can form a shock: lam[i] > lam[i+1]
        conv_p = lam_p[:-1] - lam_p[1:]   # positive = converging right-going chars
        conv_m = lam_m[:-1] - lam_m[1:]   # positive = converging left-going chars
        candidates = []
        if np.any(conv_p > 0):
            candidates.append(float(np.min(dx / conv_p[conv_p > 0])))
        if np.any(conv_m > 0):
            candidates.append(float(np.min(dx / conv_m[conv_m > 0])))
        return min(candidates) if candidates else np.inf

    # ── Newton solver for the departure equation ──────────────────────────────

    def Newton(self, x_ref, dt, newton_iter):
        """
        Implicit Euler departure for isentropic Euler in Riemann-invariant form.

        Base step (always): explicit Euler departure with old speeds.
            0 iterations → explicit Euler SL (1st order temporal).
            n iterations → n Newton corrections toward implicit Euler.

        Residual (implicit Euler equation):
            r± = R±ⁿ⁺¹ − R±ⁿ(x − λ±(Rⁿ⁺¹) · dt) = 0

        Jacobian:  J = I + dt · [[a·p_p, b·p_p], [b·p_m, a·p_m]]
            p_p = ∂R+ⁿ/∂x(x_D+),   p_m = ∂R-ⁿ/∂x(x_D-)
        """
        R_plus_n  = self.R_plus
        R_minus_n = self.R_minus
        a = (self.gamma + 1) / 4
        b = (3 - self.gamma) / 4

        lam_p0, lam_m0 = self._char_speeds(R_plus_n, R_minus_n)
        x_dep_p = x_ref - lam_p0 * dt
        x_dep_m = x_ref - lam_m0 * dt
        Rp_new  = self._interp_table(R_plus_n,  x_dep_p)
        Rm_new  = self._interp_table(R_minus_n, x_dep_m)

        p_p = np.zeros(len(x_ref))
        p_m = np.zeros(len(x_ref))
        for _ in range(newton_iter):
            lam_p_new, lam_m_new = self._char_speeds(Rp_new, Rm_new)
            x_dep_p = x_ref - lam_p_new * dt
            x_dep_m = x_ref - lam_m_new * dt

            r_p = Rp_new - self._interp_table(R_plus_n,  x_dep_p)
            r_m = Rm_new - self._interp_table(R_minus_n, x_dep_m)

            p_p = self._interp_table(R_plus_n,  x_dep_p, deriv=1)
            p_m = self._interp_table(R_minus_n, x_dep_m, deriv=1)

            det_J = (1 + a*dt*p_p) * (1 + a*dt*p_m) - (b*dt)**2 * p_p * p_m

            Rp_new = Rp_new - ((1 + a*dt*p_m) * r_p - b*dt*p_p * r_m) / det_J
            Rm_new = Rm_new - (-b*dt*p_m * r_p + (1 + a*dt*p_p) * r_m) / det_J

        p_p_full, p_m_full = p_p, p_m   # save before adaptive block clobbers p_p/p_m

        # Post-correction residual: re-evaluate at updated departure points
        lam_p_f, lam_m_f = self._char_speeds(Rp_new, Rm_new)
        x_dep_p_f = x_ref - lam_p_f * dt
        x_dep_m_f = x_ref - lam_m_f * dt
        r_p_post  = Rp_new - self._interp_table(R_plus_n,  x_dep_p_f)
        r_m_post  = Rm_new - self._interp_table(R_minus_n, x_dep_m_f)
        self._newton_residual       = np.sqrt(r_p_post**2 + r_m_post**2)
        self._newton_residual_iter1 = self._newton_residual.copy()

        # Adaptive extra iteration for cells where residual exceeds tolerance
        if self.adaptive_newton:
            mask = self._newton_residual > self.adaptive_newton_tol
            if mask.any():
                Rp_m = Rp_new[mask];  Rm_m = Rm_new[mask];  xrf = x_ref[mask]
                lam_p_m, lam_m_m = self._char_speeds(Rp_m, Rm_m)
                xdp = xrf - lam_p_m * dt
                xdm = xrf - lam_m_m * dt

                r_p = Rp_m - self._interp_table(R_plus_n,  xdp)
                r_m = Rm_m - self._interp_table(R_minus_n, xdm)
                p_p = self._interp_table(R_plus_n,  xdp, deriv=1)
                p_m = self._interp_table(R_minus_n, xdm, deriv=1)

                det = (1+a*dt*p_p)*(1+a*dt*p_m) - (b*dt)**2 * p_p*p_m
                Rp_new[mask] = Rp_m - ((1+a*dt*p_m)*r_p - b*dt*p_p*r_m) / det
                Rm_new[mask] = Rm_m - (-b*dt*p_m*r_p + (1+a*dt*p_p)*r_m) / det

                lam_pf, lam_mf = self._char_speeds(Rp_new[mask], Rm_new[mask])
                r_p2 = Rp_new[mask] - self._interp_table(R_plus_n,  xrf - lam_pf*dt)
                r_m2 = Rm_new[mask] - self._interp_table(R_minus_n, xrf - lam_mf*dt)
                self._newton_residual[mask] = np.sqrt(r_p2**2 + r_m2**2)

        return Rp_new, Rm_new, p_p_full, p_m_full, x_dep_p_f, x_dep_m_f

    def Newton_trap(self, x_ref, dt, newton_iter):
        """
        Newton solver for the coupled implicit trapezoidal departure equations.

        Unknowns: departure points (x_D+, x_D-) for both characteristics.
        Residuals (derived by substituting R±ⁿ⁺¹(x) = R±ⁿ(x_D±) into the
        trapezoidal departure formula x_D± = x − dt/2·(λ±ⁿ(x_D±)+λ±ⁿ⁺¹(x))):

            F+ = x_D+ − x + dt/2·(2a·R+ⁿ(x_D+) + b·R-ⁿ(x_D+) + b·R-ⁿ(x_D-)) = 0
            F- = x_D- − x + dt/2·(b·R+ⁿ(x_D+) + 2a·R-ⁿ(x_D-) + b·R+ⁿ(x_D-)) = 0

        Approximate 2×2 Jacobian (2 derivative calls, same structure as Newton):
            p_p = ∂R+ⁿ/∂x(x_D+),   p_m = ∂R-ⁿ/∂x(x_D-)

            J11 = 1 + a·dt·p_p          J12 = dt/2·b·p_m
            J21 = dt/2·b·p_p            J22 = 1 + a·dt·p_m

        The exact J also contains b·∂R-/∂x(x_D+) in J11 and b·∂R+/∂x(x_D-) in J22
        (same-departure cross-field terms).  These vanish for γ=3 (b=0) and are
        O(b) small elsewhere, so dropping them costs nothing at the primary
        debugging case and saves 2 derivative calls in general.

        Initial guess: explicit Euler departure with old characteristic speeds.
        One Newton iteration gives O(dt³) departure accuracy → O(dt²) globally.
        """
        R_plus_n  = self.R_plus
        R_minus_n = self.R_minus
        a = (self.gamma + 1) / 4
        b = (3 - self.gamma) / 4


        lam_p0, lam_m0 = self._char_speeds(R_plus_n, R_minus_n)
        x_dep_p = x_ref - dt * lam_p0
        x_dep_m = x_ref - dt * lam_m0

        p_p = np.zeros(len(x_ref))
        p_m = np.zeros(len(x_ref))
        for _ in range(newton_iter):
            Rp_at_p = self._interp_table(R_plus_n,  x_dep_p)
            Rm_at_p = self._interp_table(R_minus_n, x_dep_p)
            Rp_at_m = self._interp_table(R_plus_n,  x_dep_m)
            Rm_at_m = self._interp_table(R_minus_n, x_dep_m)

            F_p = x_dep_p - x_ref + 0.5*dt*(2*a*Rp_at_p + b*Rm_at_p + b*Rm_at_m)
            F_m = x_dep_m - x_ref + 0.5*dt*(b*Rp_at_p  + 2*a*Rm_at_m + b*Rp_at_m)

            p_p = self._interp_table(R_plus_n,  x_dep_p, deriv=1)
            p_m = self._interp_table(R_minus_n, x_dep_m, deriv=1)

            J11 = 1.0 + a*dt*p_p
            J12 = 0.5*dt*b*p_m
            J21 = 0.5*dt*b*p_p
            J22 = 1.0 + a*dt*p_m

            det_J = J11*J22 - J12*J21
            det_J = np.where(np.abs(det_J) > 1e-15, det_J, np.sign(det_J)*1e-15)


            x_dep_p -= ( J22*F_p - J12*F_m) / det_J
            x_dep_m -= (-J21*F_p + J11*F_m) / det_J

        Rp_new = self._interp_table(R_plus_n,  x_dep_p)
        Rm_new = self._interp_table(R_minus_n, x_dep_m)

        p_p_full, p_m_full = p_p, p_m   # save before adaptive block clobbers p_p/p_m

        # Post-correction residual: re-evaluate F at the final departure points
        Rp_p = self._interp_table(R_plus_n,  x_dep_p)
        Rm_p = self._interp_table(R_minus_n, x_dep_p)
        Rp_m = self._interp_table(R_plus_n,  x_dep_m)
        Rm_m = self._interp_table(R_minus_n, x_dep_m)
        F_p_post = x_dep_p - x_ref + 0.5*dt*(2*a*Rp_p + b*Rm_p + b*Rm_m)
        F_m_post = x_dep_m - x_ref + 0.5*dt*(b*Rp_p  + 2*a*Rm_m + b*Rp_m)
        self._newton_residual       = np.sqrt(F_p_post**2 + F_m_post**2)
        self._newton_residual_iter1 = self._newton_residual.copy()

        # Adaptive extra iteration for cells where residual exceeds tolerance
        if self.adaptive_newton:
            mask = self._newton_residual > self.adaptive_newton_tol
            if mask.any():
                xdp = x_dep_p[mask].copy()
                xdm = x_dep_m[mask].copy()
                xrf = x_ref[mask]

                Rp_at_p = self._interp_table(R_plus_n,  xdp)
                Rm_at_p = self._interp_table(R_minus_n, xdp)
                Rp_at_m = self._interp_table(R_plus_n,  xdm)
                Rm_at_m = self._interp_table(R_minus_n, xdm)

                F_p = xdp - xrf + 0.5*dt*(2*a*Rp_at_p + b*Rm_at_p + b*Rm_at_m)
                F_m = xdm - xrf + 0.5*dt*(b*Rp_at_p   + 2*a*Rm_at_m + b*Rp_at_m)

                p_p = self._interp_table(R_plus_n,  xdp, deriv=1)
                p_m = self._interp_table(R_minus_n, xdm, deriv=1)

                J11 = 1.0 + a*dt*p_p;  J12 = 0.5*dt*b*p_m
                J21 = 0.5*dt*b*p_p;    J22 = 1.0 + a*dt*p_m
                det = J11*J22 - J12*J21
                det = np.where(np.abs(det) > 1e-15, det, np.sign(det)*1e-15)

                xdp -= ( J22*F_p - J12*F_m) / det
                xdm -= (-J21*F_p + J11*F_m) / det

                Rp_new[mask] = self._interp_table(R_plus_n,  xdp)
                Rm_new[mask] = self._interp_table(R_minus_n, xdm)

                Rm_p2 = self._interp_table(R_minus_n, xdp)
                Rp_m2 = self._interp_table(R_plus_n,  xdm)
                F_p2 = xdp - xrf + 0.5*dt*(2*a*Rp_new[mask] + b*Rm_p2 + b*Rm_new[mask])
                F_m2 = xdm - xrf + 0.5*dt*(b*Rp_new[mask]   + 2*a*Rm_new[mask] + b*Rp_m2)
                self._newton_residual[mask] = np.sqrt(F_p2**2 + F_m2**2)

        return Rp_new, Rm_new, p_p_full, p_m_full, x_dep_p, x_dep_m

    # ── Step ─────────────────────────────────────────────────────────────────

    def step(self, dt):
        """
        ck=0 : implicit Euler departure (Newton on R±ⁿ⁺¹).  O(dt) global.
        ck=1 : implicit trapezoidal departure (Newton_trap on x_D±).  O(dt²) global.
        """
        R_plus  = self.R_plus
        R_minus = self.R_minus
        x_ref   = self.x_ref
        
        dt_star = self.compute_shock_dt()

        if self.ck == 0:
            Rp_new, Rm_new, p_p, p_m, x_dep_p, x_dep_m = self.Newton(x_ref, dt, self.newton_iter)
        else:
            Rp_new, Rm_new, p_p, p_m, x_dep_p, x_dep_m = self.Newton_trap(x_ref, dt, self.newton_iter)

        u_new = 0.5 * (Rp_new + Rm_new)
        if float(self.gamma) == 1.0:
            a       = np.sqrt(self.K)
            rho_new = np.exp((Rp_new - Rm_new) / (2.0 * a))
        else:
            c_new   = 0.25 * (self.gamma - 1) * (Rp_new - Rm_new)
            rho_new = (c_new**2 / (self.K * self.gamma))**(1 / (self.gamma - 1))

        self.R_plus  = Rp_new
        self.R_minus = Rm_new
        self.u       = u_new
        self.rho     = rho_new

        # Shock indicator: dimensionless gradient max(|∂R+/∂x|, |∂R-/∂x|)·dx.
        # Reuses p_p/p_m already computed at departure points during Newton — no extra interpolation.
        self._shock_indicator = np.maximum(np.abs(p_p), np.abs(p_m)) * self.dx
        
        

    def print_residual_diagnostics(self):
        """Compare post-Newton residual at shock cells vs smooth cells."""
        r1   = self._newton_residual_iter1
        res  = self._newton_residual
        detj = self._shock_indicator
        threshold = np.percentile(detj, 90)
        shock  = detj >= threshold
        smooth = detj <= np.percentile(detj, 10)
        adapted = (res < r1)   # cells where adaptive iter fired and helped

        print(f"  shock_indicator range : [{detj.min():.4f}, {detj.max():.4f}]  "
              f"shock threshold (p90) = {threshold:.4f}")
        print(f"  Residual @ shock cells  (n={shock.sum():3d}):  "
              f"iter1 max={r1[shock].max():.3e}  →  final max={res[shock].max():.3e}")
        print(f"  Residual @ smooth cells (n={smooth.sum():3d}):  "
              f"iter1 max={r1[smooth].max():.3e}  →  final max={res[smooth].max():.3e}")
        if adapted.any():
            print(f"  Adaptive fired on {adapted.sum()} cells: "
                  f"mean reduction {r1[adapted].mean():.3e} → {res[adapted].mean():.3e}  "
                  f"({r1[adapted].mean()/res[adapted].mean():.1f}x)")


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
        r=2,
        ck=1,
        newton_iter=1,
        interp_backend='lagrange',
        adaptive_newton=False,
        adaptive_newton_tol=1e-6,
        monotone_limiter=True,
        limiter_threshold=0.5,
        safety=0.99,
        shock_tol=2.0,
    ):
        self.T_max    = T_max
        self.CFL_hllc = CFL_hllc
        self.cfl      = cfl
        self.safety   = safety
        self.shock_tol = shock_tol

        dx_lpt  = abs(xmax - xmin) / Nx
        x_lpt   = np.linspace(xmin + dx_lpt / 2, xmax - dx_lpt / 2, Nx)

        Nx_hllc = N_factor * Nx
        dx_hllc = abs(xmax - xmin) / Nx_hllc
        x_hllc  = np.linspace(xmin + dx_hllc / 2, xmax - dx_hllc / 2, Nx_hllc)

        self.hllc  = HLLCSolver(x_hllc, rho_0, u_0, K, gamma, CFL_hllc)
        self.lpt   = LPTSolver(x_lpt, rho_0, u_0, K, gamma,
                               kernel=kernel, r=r,
                               ck=ck, newton_iter=newton_iter,
                               interp_backend=interp_backend,
                               adaptive_newton=adaptive_newton,
                               adaptive_newton_tol=adaptive_newton_tol,
                               monotone_limiter=monotone_limiter,
                               limiter_threshold=limiter_threshold)
        self.x_lpt = x_lpt

    def _plot_final(self):
        lpt  = self.lpt
        hllc = self.hllc
        rho_hllc_on_lpt = np.interp(self.x_lpt, hllc.x, hllc.rho)
        u_hllc_on_lpt   = np.interp(self.x_lpt, hllc.x, hllc.u)

        si   = lpt._shock_indicator
        vmin = 0.0
        vmax = float(np.max(si)) if si.max() > 0 else 1.0
        norm = plt.Normalize(vmin=vmin, vmax=vmax)
        cmap = plt.cm.RdYlGn_r   # green=smooth(0), red=shock(large)

        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        fig.suptitle(
            f"Solutions at T_max = {self.T_max:.4f}  —  dots colored by "
            f"max(|∂R+/∂x|, |∂R-/∂x|)·dx  (0=green, {vmax:.3f}=red)"
        )

        for ax, y_lpt, y_hllc, title in [
            (axes[0], lpt.rho, rho_hllc_on_lpt, "Density ρ"),
            (axes[1], lpt.u,   u_hllc_on_lpt,   "Velocity u"),
        ]:
            ax.plot(self.x_lpt, y_hllc, label="HLLC (ref)", linestyle="--", zorder=1)
            sc = ax.scatter(self.x_lpt, y_lpt, c=si, cmap=cmap, norm=norm,
                            marker="+", s=60, label="LPT ImplicitSL", zorder=2)
            fig.colorbar(sc, ax=ax, label="shock indicator")
            ax.set_title(title); ax.set_xlabel("x")
            ax.legend(); ax.grid(True)

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

    def run(self):
        print("Running HLLC (reference) ...")
        t0 = time.time()
        self.hllc.run_to(self.T_max)
        print(f"  HLLC done in {time.time()-t0:.4f} s")

        label = "LPT ImplicitSL trapezoidal" if self.lpt.ck >= 1 else "LPT ImplicitSL Euler"
        print(f"Running {label} ...")
        t = 0.0
        t0 = time.time()
        shock_detected = False
        while t < self.T_max:
            dx      = self.lpt.dx
            lam     = np.abs(self.lpt.u) + sound_speed(self.lpt.rho, self.lpt.K, self.lpt.gamma)
            dt_cfl  = self.cfl * dx / float(np.max(lam))
            dt_star = self.lpt.compute_shock_dt()

            # Stop just before characteristics cross: hand off to MGFM here.
            '''
            if dt_star <= self.shock_tol * dx:
                si = self.lpt._shock_indicator.max()
                print(f"  --> shock detected at t={t:.6f}  "
                      f"dt*={dt_star:.4e}  shock_ind={si:.3f}")
                shock_detected = True
                break

            '''
            
            plt.figure()
            plt.plot(self.x, self.rho)
            plt.show()
            
            # CFL step for Newton accuracy, capped at safety*dt_star so we
            # approach the shock landing precisely without overshooting it.
            dt = min(dt_cfl, self.safety * dt_star, self.T_max - t)

            si = self.lpt._shock_indicator.max()
            if si > 0.5:
                print(f"  t={t:.4f}  dt={dt:.4e}  dt*={dt_star:.4e}  shock_ind={si:.3f}")

            self.lpt.step(dt)
            t += dt

        elapsed = time.time() - t0
        if shock_detected:
            print(f"  SL stopped at shock formation in {elapsed:.4f} s  (t={t:.6f})")
        else:
            print(f"  SL done in {elapsed:.4f} s")

        self._plot_final()
        err_rho, err_u = self.compute_l2_error()
        print(f"\nL2 error at t={t:.4f}:")
        print(f"  rho : {err_rho:.6e}")
        print(f"  u   : {err_u:.6e}")
        return err_rho, err_u


# ── Initial conditions ────────────────────────────────────────────────────────

def rho_0(x):
    return 2. + np.exp(-x**2 / 2)

def u_0(x):
    return -0.2 * np.exp(-x**2 / 2)


# ── Run ───────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="1D isentropic Euler – RBF-GA implicit semi-Lagrangian solver")
    parser.add_argument("--T_max",        type=float, default=1.8)
    parser.add_argument("--Nx",           type=int,   default=256)
    parser.add_argument("--xmin",         type=float, default=-18.0)
    parser.add_argument("--xmax",         type=float, default=18.0)
    parser.add_argument("--K",            type=float, default=1.0)
    parser.add_argument("--gamma",        type=float, default=2.)
    parser.add_argument("--CFL_hllc",     type=float, default=0.9)
    parser.add_argument("--N_factor",     type=int,   default=1)
    parser.add_argument("--cfl",          type=float, default=1.)
    parser.add_argument("--kernel",       type=str,   default="se")
    parser.add_argument("--radius",       type=int,   default=2)
    parser.add_argument("--ck",           type=int,   default=1,
                        help="0=implicit Euler departure, 1=implicit trapezoidal departure")
    parser.add_argument("--newton_iter",  type=int,   default=1)
    parser.add_argument("--interp_backend", type=str, default='gp',
                        help="'lagrange' or 'gp'")
    parser.add_argument("--adaptive_newton", action="store_true", default=False,
                        help="retry Newton on cells where post-step residual > tol")
    parser.add_argument("--adaptive_newton_tol", type=float, default=1e-6)
    parser.add_argument("--safety",     type=float, default=0.99,
                        help="dt = safety * dt_star  (default 0.99)")
    parser.add_argument("--shock_tol",  type=float, default=2.,
                        help="stop when dt* <= shock_tol * dx  (default 2.0)")
    args = parser.parse_args()

    sim = EulerSimulation(
        xmin=args.xmin, xmax=args.xmax, Nx=args.Nx,
        rho_0=rho_0, u_0=u_0,
        K=args.K, gamma=args.gamma,
        CFL_hllc=args.CFL_hllc,
        T_max=args.T_max,
        N_factor=args.N_factor,
        kernel=args.kernel,
        cfl=args.cfl,
        r=args.radius,
        ck=args.ck,
        newton_iter=args.newton_iter,
        interp_backend=args.interp_backend,
        adaptive_newton=args.adaptive_newton,
        adaptive_newton_tol=args.adaptive_newton_tol,
        safety=args.safety,
        shock_tol=args.shock_tol,
    )
    sim.run()
