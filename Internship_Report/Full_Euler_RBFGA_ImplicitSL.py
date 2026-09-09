
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Full_Euler_RBFGA_ImplicitSL.py
-----------------------------
Implicit semi-Lagrangian (ISL) solver for 1D isentropic Euler in Riemann-invariant form.
Draft base for the full Euler extension (third field: entropy s).

RBF-GA interpolation is used to evaluate R± at departure points, with an
implicit trapezoidal departure.  Newton on departure points (x_D+, x_D-):
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
# Full-Euler EOS: p = rho^gamma * exp(sigma), sigma := ln K a per-point field
# (reduces to the isentropic p = K*rho^gamma when sigma is spatially constant).

def pressure(rho, sigma, gamma=2.0):
    return rho**gamma * np.exp(sigma)

def sound_speed(rho, sigma, gamma=2.0):
    return np.sqrt(gamma * rho**(gamma - 1.0) * np.exp(sigma))


# ── HLLC Riemann solver helpers ───────────────────────────────────────────────

def _conserved_to_primitive(U):
    rho = U[0]
    u   = U[1] / rho
    return rho, u

def _flux(U, sigma, gamma=2.0):
    rho, u = _conserved_to_primitive(U)
    p = pressure(rho, sigma, gamma)
    return np.array([rho * u, rho * u**2 + p])

def _hllc_flux(UL, UR, sigmaL, sigmaR, gamma=2.0):
    rhoL, uL = _conserved_to_primitive(UL)
    rhoR, uR = _conserved_to_primitive(UR)
    cL = sound_speed(rhoL, sigmaL, gamma)
    cR = sound_speed(rhoR, sigmaR, gamma)

    SL = min(uL - cL, uR - cR)
    SR = max(uL + cL, uR + cR)

    pL = pressure(rhoL, sigmaL, gamma)
    pR = pressure(rhoR, sigmaR, gamma)
    num   = pR - pL + rhoL * uL * (SL - uL) - rhoR * uR * (SR - uR)
    denom = rhoL * (SL - uL) - rhoR * (SR - uR)
    Sstar = num / denom

    FL = _flux(UL, sigmaL, gamma)
    FR = _flux(UR, sigmaR, gamma)

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
    """Fixed-grid Godunov solver with HLLC Riemann fluxes for 1D full Euler
    (per-point reduced entropy field sigma; reduces to isentropic when
    sigma is spatially constant)."""

    def __init__(self, x, rho_0, u_0, sigma_0=None, gamma=2.0, CFL=0.9):
        self.x     = x.copy()
        self.dx    = x[1] - x[0]
        self.gamma = gamma
        self.CFL   = CFL
        if sigma_0 is None:
            sigma_0 = lambda x: np.zeros_like(x)
        self.sigma = sigma_0(x)
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
        c   = sound_speed(rho, self.sigma, self.gamma)
        return self.CFL * self.dx / np.max(np.abs(u) + c)

    def step(self, dt):
        U  = self.U
        nx = U.shape[1]
        F  = np.zeros((2, nx + 1))
        for i in range(1, nx):
            F[:, i] = _hllc_flux(U[:, i - 1], U[:, i],
                                  self.sigma[i - 1], self.sigma[i], self.gamma)
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
    Implicit semi-Lagrangian solver for 1D full Euler in Riemann-invariant
    form, with a per-point reduced entropy field sigma := ln K
    (p = rho^gamma * exp(sigma)).  Reduces to the isentropic solver when
    sigma is spatially constant.  sigma is currently held fixed at its IC
    value (not yet advected by its own decoupled characteristic -- see
    full_euler_riemann_invariants_source_terms.tex, Section 7.1).

    The grid x_ref is fixed and uniform.  RBF-GA interpolation weights are
    precomputed once as a 1-D table w(δ), δ ∈ [−0.5, 0.5], and looked up in
    O(N·M) per interpolation call.

    Temporal scheme: implicit trapezoidal departure — Newton on (x_D+, x_D-).
    O(dt²) global.

    Parameters
    ----------
    n_table : int
        Number of δ samples in the precomputed weight table.
    """

    def __init__(self, x, rho_0, u_0, sigma_0=None,
                 gamma=2.0, kernel='se', r=2, rbfga_eps=None,
                 k_spline=5, n_table=8_000, table_deg=3,
                 weight_mode='cheb', cheb_pad=0,
                 newton_iter=1,
                 adaptive_newton=False, adaptive_newton_tol=1e-6,
                 monotone_limiter=True, limiter_threshold=0.5):
        """
        weight_mode : str
            'cheb'  — Chebyshev series fit (Clenshaw eval).  Default.
                      O(N·M·deg) per call, falls back to a precomputed
                      table near the domain boundary (see cheb_pad).
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
        self.gamma       = gamma
        self.kernel      = kernel
        self.r           = r
        self.k_spline    = k_spline
        self.table_deg   = table_deg
        self.weight_mode    = weight_mode
        self._cheb_pad      = int(cheb_pad)

        eps       = rbfga_eps if rbfga_eps is not None else dx
        self._eps = eps
        self.M    = 2 * r + 1
        self.r_d  = r + 1            # derivative stencil half-width: one wider
        self.M_d  = 2 * (r + 1) + 1 # gives O(h^{M}) derivative vs O(h^{M}) interp
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

        if sigma_0 is None:
            sigma_0 = lambda x: np.zeros_like(x)
        self.sigma = sigma_0(self.x_ref)

        self.rho = rho_0(self.x_ref)
        self.u   = u_0(self.x_ref)

        if float(gamma) == 1.0:
            a = np.exp(0.5 * self.sigma)
            self.R_plus  = self.u + a * np.log(self.rho)
            self.R_minus = self.u - a * np.log(self.rho)
        else:
            c = sound_speed(self.rho, self.sigma, gamma)
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
            W_tab[i]    = stencil.weights(x_nodes_loc,    -d * dx, deriv=0)
            W_tab_d1[i] = stencil.weights(x_nodes_loc_d1, -d * dx, deriv=1)

        self._delta_tab = delta_tab
        self._W_tab     = W_tab
        self._W_tab_d1  = W_tab_d1
        self._n_table   = n_table
        self._d_lo      = delta_tab[0]
        self._d_step    = delta_tab[1] - delta_tab[0]

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
            W_samples[i]    = self._stencil.weights(x_nodes_loc,    -d * dx, deriv=0)
            W_samples_d1[i] = self._stencil.weights(x_nodes_loc_d1, -d * dx, deriv=1)
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
                # j is pinned to [r+pad, N-r-1-pad] near the domain edges, so
                # |delta|=|j-frac_idx| can exceed 0.5 by up to ~r there (the
                # true query point is several nodes away from the pinned
                # stencil center).  The ±0.5 table only ever represents an
                # offset of at most half a cell, so clipping delta into it
                # silently swaps in the wrong stencil weights -- an error
                # that grows with r (more cells get pinned) rather than
                # shrinking, exactly backwards from what higher r should give.
                # Fix: solve the RBF-GA system directly for these few
                # boundary-adjacent points at their true (unclipped) offset,
                # reusing the same exact solve used by weight_mode='exact'.
                x_stencil_bnd = self.x_ref[si[bnd]]
                W[bnd] = self._stencil.weights_batch(x_stencil_bnd, x_dep[bnd], deriv=deriv)

        elif self.weight_mode == 'exact':
            # O(N·M³): solve the RBF-GA system exactly at each departure point
            x_stencil = self.x_ref[si]
            W = self._stencil.weights_batch(x_stencil, x_dep, deriv=deriv)

        else:  # 'table'
            # Same ±0.5-only validity issue as the 'cheb' boundary fallback
            # above: j is pinned near the domain edges, so |delta| can exceed
            # 0.5 there.  Route those few points through the exact solve;
            # everywhere else (interior) delta is always in [-0.5, 0.5] by
            # construction (unclipped round-to-nearest), so the table applies.
            interior = np.abs(delta) <= 0.5
            W        = np.empty((len(delta), M))
            bnd      = ~interior
            if bnd.any():
                x_stencil_bnd = self.x_ref[si[bnd]]
                W[bnd] = self._stencil.weights_batch(x_stencil_bnd, x_dep[bnd], deriv=deriv)
            if interior.any():
                t = (delta[interior] - self._d_lo) / self._d_step
                if self.table_deg == 1:
                    t   = np.clip(t, 0.0, self._n_table - 1 - 1e-10)
                    idx = t.astype(int)
                    f   = (t - idx)[:, None]
                    W[interior] = (1.0 - f) * W_tab[idx] + f * W_tab[idx + 1]
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
                    W[interior] = (c0 * W_tab[idx - 1]
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
            c = np.exp(0.5 * self.sigma)
        else:
            c = (self.gamma - 1) / 4 * (R_plus - R_minus)
        lam_p = u + c
        lam_m = u - c
        return lam_p, lam_m

    def _Lambda(self, R_plus, R_minus):
        """Lambda := c^2/(gamma*(gamma-1)), c recovered algebraically from
        (R_plus, R_minus) via eq. (1) -- the coefficient of the entropy-gradient
        source term in D+/-R+/-/Dt = Lambda * dsigma/dx (full_euler_riemann_
        invariants_source_terms.tex, (R+)/(R-))."""
        gamma = self.gamma
        return (gamma - 1.0) * (R_plus - R_minus)**2 / (16.0 * gamma)

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

    def Newton_trap(self, x_ref, dt, newton_iter):
        """
        Newton solver for the coupled implicit trapezoidal departure equations,
        now including the entropy-gradient source term of full_euler_riemann_
        invariants_source_terms.tex (Section 5, eqs. (R+)/(R-)) for non-constant
        sigma.

        Unknowns: departure points (x_D+, x_D-) for both characteristics.
        Residuals (derived by substituting
            R±ⁿ⁺¹(x) = R±ⁿ(x_D±) + src±(x)
        -- the trapezoidal quadrature of the source integral, Section 6.2 --
        into the trapezoidal departure formula
        x_D± = x − dt/2·(λ±ⁿ(x_D±)+λ±ⁿ⁺¹(x))):

            F+ = x_D+ − x + dt/2·(2a·R+ⁿ(x_D+) + b·R-ⁿ(x_D+) + b·R-ⁿ(x_D-))
                          + dt/2·(a·src+ + b·src-) = 0
            F- = x_D- − x + dt/2·(b·R+ⁿ(x_D+) + 2a·R-ⁿ(x_D-) + b·R+ⁿ(x_D-))
                          + dt/2·(b·src+ + a·src-) = 0

        with a=(γ+1)/4, b=(3−γ)/4, and the trapezoidal source quadrature
            src+ = dt/2·[Λ∂ₓσ|_{x_D+} + Λ∂ₓσ|_x]   (Λ = c²/(γ(γ−1)))
            src- = dt/2·[Λ∂ₓσ|_{x_D-} + Λ∂ₓσ|_x]
        the arrival-side term Λ∂ₓσ|_x uses the lagged (already-known) sigmaⁿ
        field, per the "Cost of each endpoint" discussion (Section 6.2). The
        source is folded into the residual (not added post-hoc) to preserve
        the scheme's O(dt²) global accuracy (Section 5.4); its own derivative
        is dropped from the Jacobian (O(dt²) smaller, costs only convergence
        rate, never accuracy).

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
        sigma_n   = self.sigma
        gamma     = self.gamma
        has_source = (float(gamma) != 1.0)
        a = (self.gamma + 1) / 4
        b = (3 - self.gamma) / 4

        # Arrival-side source contribution Λ∂ₓσ|_x: fixed for the whole step,
        # uses only already-known (lagged) fields at the grid nodes.
        if has_source:
            Lambda_x = self._Lambda(R_plus_n, R_minus_n)
            dsigma_x = self._interp_table(sigma_n, x_ref, deriv=1)
            src_x    = Lambda_x * dsigma_x
        else:
            src_x = np.zeros(len(x_ref))

        def _src_terms(x_dep_p_, x_dep_m_, Rp_at_p_, Rm_at_p_, Rp_at_m_, Rm_at_m_):
            if not has_source:
                z = np.zeros(len(x_dep_p_))
                return z, z
            Lambda_p = self._Lambda(Rp_at_p_, Rm_at_p_)
            Lambda_m = self._Lambda(Rp_at_m_, Rm_at_m_)
            dsigma_p = self._interp_table(sigma_n, x_dep_p_, deriv=1)
            dsigma_m = self._interp_table(sigma_n, x_dep_m_, deriv=1)
            src_p = 0.5*dt*(Lambda_p*dsigma_p + src_x)
            src_m = 0.5*dt*(Lambda_m*dsigma_m + src_x)
            return src_p, src_m

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

            src_p, src_m = _src_terms(x_dep_p, x_dep_m, Rp_at_p, Rm_at_p, Rp_at_m, Rm_at_m)

            F_p = (x_dep_p - x_ref + 0.5*dt*(2*a*Rp_at_p + b*Rm_at_p + b*Rm_at_m)
                   + 0.5*dt*(a*src_p + b*src_m))
            F_m = (x_dep_m - x_ref + 0.5*dt*(b*Rp_at_p  + 2*a*Rm_at_m + b*Rp_at_m)
                   + 0.5*dt*(b*src_p + a*src_m))

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

        p_p_full, p_m_full = p_p, p_m   # save before adaptive block clobbers p_p/p_m

        # Final source-corrected values and post-correction residual, both
        # evaluated fresh at the converged departure points.
        Rp_p = self._interp_table(R_plus_n,  x_dep_p)
        Rm_p = self._interp_table(R_minus_n, x_dep_p)
        Rp_m = self._interp_table(R_plus_n,  x_dep_m)
        Rm_m = self._interp_table(R_minus_n, x_dep_m)
        src_p, src_m = _src_terms(x_dep_p, x_dep_m, Rp_p, Rm_p, Rp_m, Rm_m)

        Rp_new = Rp_p + src_p
        Rm_new = Rm_m + src_m

        F_p_post = (x_dep_p - x_ref + 0.5*dt*(2*a*Rp_p + b*Rm_p + b*Rm_m)
                    + 0.5*dt*(a*src_p + b*src_m))
        F_m_post = (x_dep_m - x_ref + 0.5*dt*(b*Rp_p  + 2*a*Rm_m + b*Rp_m)
                    + 0.5*dt*(b*src_p + a*src_m))
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

                src_p_a, src_m_a = _src_terms(xdp, xdm, Rp_at_p, Rm_at_p, Rp_at_m, Rm_at_m)

                F_p = (xdp - xrf + 0.5*dt*(2*a*Rp_at_p + b*Rm_at_p + b*Rm_at_m)
                       + 0.5*dt*(a*src_p_a + b*src_m_a))
                F_m = (xdm - xrf + 0.5*dt*(b*Rp_at_p   + 2*a*Rm_at_m + b*Rp_at_m)
                       + 0.5*dt*(b*src_p_a + a*src_m_a))

                p_p = self._interp_table(R_plus_n,  xdp, deriv=1)
                p_m = self._interp_table(R_minus_n, xdm, deriv=1)

                J11 = 1.0 + a*dt*p_p;  J12 = 0.5*dt*b*p_m
                J21 = 0.5*dt*b*p_p;    J22 = 1.0 + a*dt*p_m
                det = J11*J22 - J12*J21
                det = np.where(np.abs(det) > 1e-15, det, np.sign(det)*1e-15)

                xdp -= ( J22*F_p - J12*F_m) / det
                xdm -= (-J21*F_p + J11*F_m) / det

                Rp_p2 = self._interp_table(R_plus_n,  xdp)
                Rm_p2 = self._interp_table(R_minus_n, xdp)
                Rp_m2 = self._interp_table(R_plus_n,  xdm)
                Rm_m2 = self._interp_table(R_minus_n, xdm)
                src_p2, src_m2 = _src_terms(xdp, xdm, Rp_p2, Rm_p2, Rp_m2, Rm_m2)

                Rp_new[mask] = Rp_p2 + src_p2
                Rm_new[mask] = Rm_m2 + src_m2

                F_p2 = (xdp - xrf + 0.5*dt*(2*a*Rp_p2 + b*Rm_p2 + b*Rm_m2)
                        + 0.5*dt*(a*src_p2 + b*src_m2))
                F_m2 = (xdm - xrf + 0.5*dt*(b*Rp_p2   + 2*a*Rm_m2 + b*Rp_m2)
                        + 0.5*dt*(b*src_p2 + a*src_m2))
                self._newton_residual[mask] = np.sqrt(F_p2**2 + F_m2**2)

                x_dep_p[mask] = xdp
                x_dep_m[mask] = xdm

        return Rp_new, Rm_new, p_p_full, p_m_full, x_dep_p, x_dep_m

    def _solve_sigma_departure(self, x_ref, dt, u_new, newton_iter):
        """
        Decoupled scalar departure solve for sigma's own characteristic,
        dx/dt = u (Section 6.3 of full_euler_riemann_invariants_source_terms.tex):

            x_D^sigma = x − dt/2·(uⁿ(x_D^sigma) + uⁿ⁺¹(x))

        Source-free (Dsigma/Dt = 0 exactly), so the update is pure advection:
        sigmaⁿ⁺¹(x) = interp(sigmaⁿ, x_D^sigma).  uⁿ⁺¹(x) = u_new is already
        known at this point (block-triangular structure: the (R+,R-) block
        does not depend on sigma, sigma's row depends only on the already-
        resolved u).  Newton on the single scalar unknown x_D^sigma, same
        trapezoidal philosophy as Newton_trap.
        """
        u_n     = self.u        # old (tⁿ) field -- must be read before self.u is overwritten
        sigma_n = self.sigma

        x_dep_s = x_ref - dt * u_n   # explicit Euler initial guess

        for _ in range(newton_iter):
            u_at_s = self._interp_table(u_n, x_dep_s)
            F      = x_dep_s - x_ref + 0.5*dt*(u_at_s + u_new)
            du_at_s = self._interp_table(u_n, x_dep_s, deriv=1)
            J       = 1.0 + 0.5*dt*du_at_s
            J       = np.where(np.abs(J) > 1e-15, J, np.sign(J)*1e-15)
            x_dep_s -= F / J

        sigma_new = self._interp_table(sigma_n, x_dep_s)
        return sigma_new, x_dep_s

    # ── Step ─────────────────────────────────────────────────────────────────

    def step(self, dt):
        """
        Implicit trapezoidal departure (Newton_trap on x_D±).  O(dt²) global.
        Solve order (Section 7.1): (i) 2×2 system for R+,R- (residual includes
        the sigma source, using old sigma); (ii) recover uⁿ⁺¹; (iii) decoupled
        scalar departure solve for x_D^sigma using that uⁿ⁺¹; (iv) update
        sigmaⁿ⁺¹, then reconstruct rho using the *new* sigma (consistent EOS
        evaluation point).
        """
        R_plus  = self.R_plus
        R_minus = self.R_minus
        x_ref   = self.x_ref

        dt_star = self.compute_shock_dt()

        Rp_new, Rm_new, p_p, p_m, x_dep_p, x_dep_m = self.Newton_trap(x_ref, dt, self.newton_iter)

        u_new = 0.5 * (Rp_new + Rm_new)
        if float(self.gamma) == 1.0:
            a         = np.exp(0.5 * self.sigma)
            rho_new   = np.exp((Rp_new - Rm_new) / (2.0 * a))
            sigma_new = self.sigma   # isothermal branch: sigma not advected here
        else:
            c_new     = 0.25 * (self.gamma - 1) * (Rp_new - Rm_new)
            sigma_new, x_dep_s = self._solve_sigma_departure(x_ref, dt, u_new, self.newton_iter)
            rho_new   = (c_new**2 / (self.gamma * np.exp(sigma_new)))**(1 / (self.gamma - 1))

        self.R_plus  = Rp_new
        self.R_minus = Rm_new
        self.u       = u_new
        self.rho     = rho_new
        self.sigma   = sigma_new

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
        sigma_0=None, gamma=2.0,
        CFL_hllc=0.9,
        T_max=0.4,
        N_factor=8,
        kernel="se",
        cfl=1,
        r=2,
        newton_iter=1,
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

        self.hllc  = HLLCSolver(x_hllc, rho_0, u_0, sigma_0, gamma, CFL_hllc)
        self.lpt   = LPTSolver(x_lpt, rho_0, u_0, sigma_0, gamma,
                               kernel=kernel, r=r,
                               newton_iter=newton_iter,
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

        p_lpt           = pressure(lpt.rho,  lpt.sigma,  lpt.gamma)
        p_hllc          = pressure(hllc.rho, hllc.sigma, hllc.gamma)
        p_hllc_on_lpt   = np.interp(self.x_lpt, hllc.x, p_hllc)

        si   = lpt._shock_indicator
        vmin = 0.0
        vmax = float(np.max(si)) if si.max() > 0 else 1.0
        norm = plt.Normalize(vmin=vmin, vmax=vmax)
        cmap = plt.cm.RdYlGn_r   # green=smooth(0), red=shock(large)

        fig = plt.figure(figsize=(13, 9))
        gs = fig.add_gridspec(2, 2)
        ax_rho = fig.add_subplot(gs[0, 0])
        ax_u   = fig.add_subplot(gs[0, 1])
        ax_p   = fig.add_subplot(gs[1, :])
        fig.suptitle(
            f"Solutions at T_max = {self.T_max:.4f}  —  dots colored by "
            f"max(|∂R+/∂x|, |∂R-/∂x|)·dx  (0=green, {vmax:.3f}=red)"
        )

        for ax, y_lpt, y_hllc, title in [
            (ax_rho, lpt.rho, rho_hllc_on_lpt, "Density ρ"),
            (ax_u,   lpt.u,   u_hllc_on_lpt,   "Velocity u"),
            (ax_p,   p_lpt,   p_hllc_on_lpt,   "Pressure p"),
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

        print("Running LPT ImplicitSL trapezoidal ...")
        t = 0.0
        t0 = time.time()
        shock_detected = False
        while t < self.T_max:
            dx      = self.lpt.dx
            lam     = np.abs(self.lpt.u) + sound_speed(self.lpt.rho, self.lpt.sigma, self.lpt.gamma)
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
            '''
            plt.figure()
            plt.plot(pressure(self.lpt.rho,self.lpt.sigma,self.lpt.gamma))
            plt.show()
            '''
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

def sigma_0(x):
    # Spatially varying reduced entropy profile; amplitude set by --sigma0
    # (sigma0=0 reduces to the isentropic sigma=0 special case).  Full-period
    # cosine spanning the domain (wavenumber inferred from x itself, so it
    # stays resolvable regardless of domain size), with zero slope at BOTH
    # the domain center and the domain boundaries.  Zero slope at the center
    # matters because the IC pulses (rho_0, u_0) sit there -- a domain-
    # spanning sine instead puts its steepest gradient at that same point
    # and drastically front-loads shock formation.  Zero slope at the
    # boundaries matters separately: a half-period bump (steepest gradient
    # at the edges) was tried first and made the LPT/HLLC boundary
    # treatments disagree there, since the entropy-gradient source term is
    # largest exactly where the LPT departure-point derivative stencil is
    # least accurate (domain edge).
    L  = x.max() - x.min()
    xc = 0.5 * (x.max() + x.min())
    return np.cos(2 * np.pi * (x - xc) / L)


# ── Run ───────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="1D isentropic Euler – RBF-GA implicit semi-Lagrangian solver")
    parser.add_argument("--T_max",        type=float, default=2.5)
    parser.add_argument("--Nx",           type=int,   default=256)
    parser.add_argument("--xmin",         type=float, default=-10.0)
    parser.add_argument("--xmax",         type=float, default=10.0)
    parser.add_argument("--sigma0",        type=float, default=0.1,
                        help="amplitude of the spatially varying sigma=ln K profile "
                             "cos(2 pi (x-xc)/L); sigma0=0 <=> uniform K=1 (isentropic default)")
    parser.add_argument("--gamma",        type=float, default=1.4)
    parser.add_argument("--CFL_hllc",     type=float, default=0.9)
    parser.add_argument("--N_factor",     type=int,   default=1)
    parser.add_argument("--cfl",          type=float, default=1.)
    parser.add_argument("--kernel",       type=str,   default="se")
    parser.add_argument("--radius",       type=int,   default=1)
    parser.add_argument("--newton_iter",  type=int,   default=1)
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
        sigma_0=lambda x: args.sigma0 * sigma_0(x), gamma=args.gamma,
        CFL_hllc=args.CFL_hllc,
        T_max=args.T_max,
        N_factor=args.N_factor,
        kernel=args.kernel,
        cfl=args.cfl,
        r=args.radius,
        newton_iter=args.newton_iter,
        adaptive_newton=args.adaptive_newton,
        adaptive_newton_tol=args.adaptive_newton_tol,
        safety=args.safety,
        shock_tol=args.shock_tol,
    )
    sim.run()
