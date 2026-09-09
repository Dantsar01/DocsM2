#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
1D_Burgers_BFECC.py
-------------------
Backward semi-Lagrangian GP solver for the inviscid Burgers equation

    ∂_t u + u ∂_x u = 0

with smooth (Gaussian) initial data.  Purpose: debug the gamma=3 isentropic
Euler code (1D_Euler_RBFGA_ImplicitSL.py), because at gamma=3 both Riemann
invariants R± satisfy independent Burgers equations.

Only exact GP interpolation is used (weight_mode='exact'); no polynomial /
Chebyshev / table approximations.

Three solvers are compared:
  1. BurgersExact     – method-of-characteristics, machine-precision reference.
  2. BurgersSpectral  – pseudo-spectral (Fourier + DOP853), ~1e-10 accuracy.
  3. BurgersSolver    – backward semi-Lagrangian + Newton, exact GP interp.

Step modes (ck parameter):
  ck=0  — plain implicit Euler departure:
            x_D = x - u^{n+1}(x_D)·dt   (Newton)
  ck=1  — midpoint / leapfrog (same pattern as Euler code):
            1. half-step Newton from x_ref → u^{n+1/2}
            2. full-step Newton from x_ref - dt/2·u^{n+1/2}
            Global order O(dt²).
"""

import numpy as np
import matplotlib.pyplot as plt
import time

import os, sys
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(current_dir)
from rbf_ga_weights_1d import RBFGAStencil

from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline, make_interp_spline
from scipy.optimize  import brentq


# ─────────────────────────────────────────────────────────────────────────────
# Exact characteristics solution
# ─────────────────────────────────────────────────────────────────────────────

class BurgersExact:
    """
    Exact solution of  ∂_t u + u ∂_x u = 0  via the method of characteristics.

    u(x,t) = u_0(ξ)  where  ξ + u_0(ξ)·t = x.
    ξ is found using Newton's method with CubicSpline for u_0.
    Brent fallback for robustness near shocks.

    Parameters
    ----------
    x_ic  : 1-D array of IC sample points (ascending, covers extended domain)
    u0    : u_0(x_ic)
    """

    def __init__(self, x_ic, u0, newton_tol=1e-13, newton_max=40):
        self._x_ic  = np.asarray(x_ic, dtype=float)
        self._cs    = CubicSpline(x_ic, u0)
        self._tol   = newton_tol
        self._max   = newton_max

    def _solve_xi_pointwise(self, x_i, t):
        cs = self._cs
        # Newton: f(ξ) = ξ + cs(ξ)·t - x_i
        xi = x_i  # initial guess ξ = x
        for _ in range(self._max):
            f  = xi + cs(xi) * t - x_i
            fp = 1.0 + cs(xi, 1) * t
            dxi = f / (fp if abs(fp) > 1e-14 else 1e-14)
            xi -= dxi
            if abs(dxi) < self._tol:
                return float(xi)
        # Brent fallback
        u_max = float(np.abs(self._cs(self._x_ic)).max()) + 1.0
        lo = x_i - u_max * t - 1.0
        hi = x_i + u_max * abs(t) + 1.0
        try:
            xi = brentq(lambda z: z + cs(z) * t - x_i,
                        lo, hi, xtol=self._tol, rtol=self._tol)
        except ValueError:
            xi = x_i  # outside wave support: u ≈ u_0
        return float(xi)

    def solve(self, x_eval, t):
        """Return u(x_eval, t).  Machine-precision for smooth (pre-shock) data."""
        x_eval = np.asarray(x_eval, dtype=float)
        if t == 0.0:
            return self._cs(x_eval)
        u_out = np.empty_like(x_eval)
        for i, x_i in enumerate(x_eval):
            xi       = self._solve_xi_pointwise(x_i, t)
            u_out[i] = float(self._cs(xi))
        return u_out

    def l2_error(self, x_eval, u_num, t):
        u_ex = self.solve(x_eval, t)
        dx   = float(x_eval[1] - x_eval[0])
        return float(np.sqrt(np.sum((u_num - u_ex)**2) * dx))


def build_burgers_exact(u_0_fn, xmin, xmax, T_max=1.0, N_ic=16384):
    """
    Build a BurgersExact reference from a callable IC.

    The IC is sampled on an extended domain so that departure points at
    t ≤ T_max always fall inside the sample region.
    """
    x_phys  = np.linspace(xmin, xmax, 4096)
    u_phys  = u_0_fn(x_phys)
    u_max   = float(np.max(np.abs(u_phys))) + 1.0
    margin  = u_max * T_max
    x_ic    = np.linspace(xmin - margin, xmax + margin, N_ic)
    return BurgersExact(x_ic, u_0_fn(x_ic))


# ─────────────────────────────────────────────────────────────────────────────
# Pseudo-spectral reference
# ─────────────────────────────────────────────────────────────────────────────

class BurgersSpectral:
    """
    Pre-computed pseudo-spectral reference for inviscid Burgers.

    Conservative form:  ∂_t u + ∂_x(u²/2) = 0
    Fourier pseudo-spectral + 2/3 dealiasing, DOP853 time integration.
    Combined accuracy: ~1e-10 for smooth solutions.
    """

    def __init__(self, x_out, u_f, t_target):
        self.t_target = t_target
        self._cs      = make_interp_spline(x_out, u_f, k=5)

    def solve(self, x_eval, t=None):
        import warnings
        if t is not None and not np.isclose(float(t), self.t_target, rtol=1e-5):
            warnings.warn(
                f"Spectral ref built at t={self.t_target:.6g}; "
                f"requested t={float(t):.6g}.  Returning built solution.",
                stacklevel=2,
            )
        return self._cs(x_eval)

    def l2_error(self, x_eval, u_num, t=None):
        u_ex = self.solve(x_eval, t)
        dx   = float(x_eval[1] - x_eval[0])
        return float(np.sqrt(np.sum((u_num - u_ex)**2) * dx))


def build_burgers_spectral(
    u_0_fn,
    xmin, xmax,
    T_target,
    N_spec   = 4096,
    D_factor = 4.0,
    rtol     = 1e-11,
    atol     = 1e-13,
    verbose  = True,
):
    """
    Build a BurgersSpectral reference from a callable IC.

    Parameters
    ----------
    u_0_fn   : callable  x → array
    xmin, xmax : physical domain
    T_target : target time
    N_spec   : Fourier modes on the extended domain
    D_factor : extended domain length = D_factor * (xmax - xmin)
    """
    import warnings

    L_phys = xmax - xmin
    cx     = 0.5 * (xmin + xmax)
    L_ext  = D_factor * L_phys
    x0_ext = cx - 0.5 * L_ext
    x_spec = x0_ext + np.arange(N_spec) * (L_ext / N_spec)

    u0     = u_0_fn(x_spec)

    v_max  = float(np.max(np.abs(u0))) + 1e-14
    reach  = v_max * T_target
    margin = 0.5 * (L_ext - L_phys)
    if reach >= margin:
        warnings.warn(
            f"Waves may reach periodic boundary before T={T_target}: "
            f"reach={reach:.3f} >= margin={margin:.3f}.  Increase D_factor.",
            stacklevel=2,
        )

    if verbose:
        print(f"[spectral ref]  N_spec={N_spec},  L_ext={L_ext:.1f},  "
              f"v_max≈{v_max:.3f},  reach≈{reach:.3f},  margin={margin:.3f}")

    k_arr  = 2.0 * np.pi / L_ext * np.arange(N_spec // 2 + 1, dtype=float)
    cutoff = N_spec // 3

    def rhs(t, state):
        u_hat          = np.fft.rfft(state)
        u_hat[cutoff:] = 0.0
        u              = np.fft.irfft(u_hat, n=N_spec)
        flux_hat          = np.fft.rfft(0.5 * u**2)
        flux_hat[cutoff:] = 0.0
        return np.fft.irfft(-1j * k_arr * flux_hat, n=N_spec)

    if verbose:
        print(f"  Integrating Burgers t ∈ [0, {T_target}] with DOP853 ...")

    sol = solve_ivp(rhs, [0.0, T_target], u0,
                    method='DOP853', rtol=rtol, atol=atol, dense_output=False)

    if not sol.success:
        raise RuntimeError(f"DOP853 failed: {sol.message}")

    u_f = sol.y[:, -1]

    if verbose:
        u_hat    = np.fft.rfft(u_f)
        amp_max  = float(np.abs(u_hat).max())
        amp_tail = float(np.abs(u_hat[cutoff - 10 : cutoff]).mean())
        ratio    = amp_tail / amp_max if amp_max > 0 else np.nan
        print(f"  nfev={sol.nfev},  spectral tail |û_tail|/|û_max| = {ratio:.2e} "
              f"({'OK' if ratio < 1e-8 else 'WARN: increase N_spec'})")

    buf   = 0.1 * L_phys
    mask  = (x_spec >= xmin - buf) & (x_spec <= xmax + buf)
    return BurgersSpectral(x_spec[mask], u_f[mask], float(sol.t[-1]))


# ─────────────────────────────────────────────────────────────────────────────
# GP semi-Lagrangian Burgers solver
# ─────────────────────────────────────────────────────────────────────────────

class BurgersSolver:
    """
    Backward semi-Lagrangian solver for  ∂_t u + u ∂_x u = 0.

    RBF-GA (GP) interpolation with selectable weight mode (same infrastructure
    as LPTSolver in 1D_Euler_RBFGA_ImplicitSL.py).

    Parameters
    ----------
    x_ref       : uniform reference grid (cell-centred)
    u_0_fn      : callable, initial condition
    r           : stencil half-width (default 2 → 5-point stencil)
    eps         : RBF-GA shape parameter
    ck          : 0 = plain Euler SL,  1 = midpoint (2nd order)
    newton_iter : Newton iterations per departure solve
    weight_mode : 'exact' | 'poly' | 'table'
                  'exact' — solve the M×M RBF-GA system at every departure point.
                            O(N·M³) per call, zero weight error.  Reference only.
                  'poly'  — polynomial fit to w_k(δ) ∈ [−0.5,0.5].
                            O(N·M·deg) per call, machine-ε accuracy.  Default.
                  'table' — uniform lookup table with Catmull-Rom cubic interp.
                            O(N·M) per call, slight approximation error.
    n_table     : number of δ samples in the precomputed weight table (used by
                  'table' mode and as boundary fallback in 'poly' mode).
    table_deg   : 1 = linear, 3 = Catmull-Rom cubic (only for 'table' mode).
    """

    def __init__(self, x_ref, u_0_fn,
                 r=2, eps=1e-3, ck=1, newton_iter=3,
                 weight_mode='poly', n_table=8_000, table_deg=3):
        self.x_ref      = np.asarray(x_ref, dtype=float)
        self.Nx         = len(x_ref)
        self.dx         = float(x_ref[1] - x_ref[0])
        self.r          = r
        self.r_d        = r + 1
        self.M          = 2 * r + 1
        self.M_d        = 2 * (r + 1) + 1
        self.ck         = ck
        self.newton_iter    = newton_iter
        self.weight_mode    = weight_mode
        self.table_deg      = table_deg
        self._eps           = eps
        self._stencil       = RBFGAStencil(eps)

        self._build_weight_table(n_table)
        if weight_mode == 'poly':
            self._build_poly_weights()
        elif weight_mode not in ('exact', 'table'):
            raise ValueError(f"weight_mode must be 'exact', 'poly', or 'table', got {weight_mode!r}")

        self._build_beta_precomp()
        self.u = u_0_fn(self.x_ref).copy()

    # ── GP-WENO smoothness indicator precomputation ───────────────────────

    def _build_beta_precomp(self):
        """
        Precompute K^{-1} for the (r+1)-node sub-stencils used in GP-WENO.

        On a uniform grid all R+1 sub-stencils of size m=r+1 share the same
        inter-node spacing dx, so K is identical for every sub-stencil and
        every position.  K^{-1} is computed once via direct inversion (safe
        for fixed eps and small m).

        Sub-stencil layout (r+1 stencils, each of size r+1):
            S_s : node offsets  [s-r, s-r+1, ..., s]  for s = 0 ... r
        """
        dx  = self.dx
        #eps = self._eps
        eps = 3 * self.dx
        m   = self.r + 1                              # sub-stencil size
        x_sub = np.arange(m, dtype=float) * dx        # relative positions [0, dx, ..., (m-1)*dx]
        K     = np.exp(-eps**2 * (x_sub[:, None] - x_sub[None, :])**2)
        self._beta_K_inv = np.linalg.inv(K)           # (m, m) — precomputed once
        self._beta_m     = m
        self._beta_R1    = self.r + 1                 # number of sub-stencils

    def compute_beta(self):
        """
        Compute GP-WENO smoothness indicators beta_m at all N grid nodes.

        Returns
        -------
        beta : (N, R+1) array
            beta[i, s] = f_{S_s}^T K^{-1} f_{S_s}   (quadratic form, >= 0)

        High beta → data inconsistent with smooth GP prior → near a shock.
        Low  beta → smooth region.
        """
        N     = self.Nx
        u     = self.u
        r     = self.r
        m     = self._beta_m
        K_inv = self._beta_K_inv
        R1    = self._beta_R1

        beta = np.zeros((N, R1))
        for s in range(R1):
            offsets = np.arange(s - r, s + 1)                    # (m,) node offsets
            idx     = np.clip(
                np.arange(N)[:, None] + offsets[None, :],
                0, N - 1,
            )                                                     # (N, m)
            F     = u[idx]                                        # (N, m) stencil values
            alpha = F @ K_inv                                     # (N, m) = K^{-1} f_s per row
            beta[:, s] = (F * alpha).sum(axis=1)                  # (N,) quadratic form
        return beta

    # ── Weight table (always built; used as boundary fallback) ────────────

    def _build_weight_table(self, n_table):
        dx   = self.dx
        r    = self.r;    r_d = self.r_d
        M    = self.M;    M_d = self.M_d
        eps  = self._eps
        stencil        = RBFGAStencil(eps)
        x_nodes_loc    = dx * np.arange(-r,   r   + 1, dtype=float)
        x_nodes_loc_d1 = dx * np.arange(-r_d, r_d + 1, dtype=float)
        delta_tab      = np.linspace(-0.5, 0.5, n_table)
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

    # ── Polynomial weight representation ──────────────────────────────────

    def _build_poly_weights(self):
        dx    = self.dx
        r     = self.r;    r_d = self.r_d
        M     = self.M;    M_d = self.M_d
        deg   = min(3 * r   + 3, 16)
        deg_d = min(3 * r_d + 3, 16)
        K     = 2 * max(deg, deg_d) + 4
        x_nodes_loc    = dx * np.arange(-r,   r   + 1, dtype=float)
        x_nodes_loc_d1 = dx * np.arange(-r_d, r_d + 1, dtype=float)
        delta_nodes    = np.linspace(-0.5, 0.5, K)
        W_samples    = np.zeros((K, M))
        W_samples_d1 = np.zeros((K, M_d))
        for i, d in enumerate(delta_nodes):
            W_samples[i]    = self._stencil.weights(x_nodes_loc,    -d * dx, deriv=0)
            W_samples_d1[i] = self._stencil.weights(x_nodes_loc_d1, -d * dx, deriv=1)
        self._poly_coeffs    = np.array([
            np.polyfit(delta_nodes, W_samples[:, m],    deg)   for m in range(M)
        ])
        self._poly_coeffs_d1 = np.array([
            np.polyfit(delta_nodes, W_samples_d1[:, m], deg_d) for m in range(M_d)
        ])

    # ── Interpolation / differentiation ───────────────────────────────────

    def _interp(self, f_src, x_dep, deriv=0):
        """
        Interpolate (deriv=0) or differentiate (deriv=1) f_src at x_dep.
        Dispatches to the chosen weight_mode; boundary points use constant BC.
        """
        r  = self.r   if deriv == 0 else self.r_d
        M  = self.M   if deriv == 0 else self.M_d
        N  = self.Nx
        dx = self.dx
        x0 = self.x_ref[0]

        xmin_  = x0 - 0.5 * dx
        xmax_  = self.x_ref[-1] + 0.5 * dx
        out_L  = x_dep < xmin_
        out_R  = x_dep > xmax_
        if out_L.any() or out_R.any():
            result          = np.empty(len(x_dep))
            result[out_L]   = f_src[0]  if deriv == 0 else 0.0
            result[out_R]   = f_src[-1] if deriv == 0 else 0.0
            mask = ~(out_L | out_R)
            if mask.any():
                result[mask] = self._interp(f_src, x_dep[mask], deriv)
            return result

        frac_idx = (x_dep - x0) / dx
        j  = np.round(frac_idx).astype(int)
        j  = np.clip(j, r, N - r - 1)
        proto = np.arange(-r, r + 1, dtype=int)
        si    = np.clip(j[:, None] + proto[None, :], 0, N - 1)
        delta = j.astype(float) - frac_idx   # δ ∈ [−0.5, 0.5]

        W_tab = self._W_tab if deriv == 0 else self._W_tab_d1

        if self.weight_mode == 'exact':
            x_stencil = self.x_ref[si]
            W = self._stencil.weights_batch(x_stencil, x_dep, deriv=deriv)

        elif self.weight_mode == 'poly':
            poly_c       = self._poly_coeffs if deriv == 0 else self._poly_coeffs_d1
            W            = np.zeros((len(delta), M))
            near_bnd     = (j <= r) | (j >= N - r - 1)
            interior     = (np.abs(delta) <= 0.5) & ~near_bnd
            if interior.any():
                d_int = delta[interior]
                for m in range(M):
                    W[interior, m] = np.polyval(poly_c[m], d_int)
            tbl = near_bnd & (np.abs(delta) <= 0.5)
            if tbl.any():
                t_tbl = np.clip((delta[tbl] - self._d_lo) / self._d_step,
                                1.0, self._n_table - 3.0)
                idx_t = t_tbl.astype(int)
                f_t   = (t_tbl - idx_t)[:, None]
                f2t, f3t = f_t * f_t, f_t * f_t * f_t
                c0 = -0.5*f3t + 1.0*f2t - 0.5*f_t
                c1 =  1.5*f3t - 2.5*f2t + 1.0
                c2 = -1.5*f3t + 2.0*f2t + 0.5*f_t
                c3 =  0.5*f3t - 0.5*f2t
                W[tbl] = (c0 * W_tab[idx_t - 1] + c1 * W_tab[idx_t]
                        + c2 * W_tab[idx_t + 1] + c3 * W_tab[idx_t + 2])
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
            else:  # Catmull-Rom cubic
                t   = np.clip(t, 1.0, self._n_table - 3.0)
                idx = t.astype(int)
                f   = (t - idx)[:, None]
                f2, f3 = f * f, f * f * f
                c0 = -0.5*f3 + 1.0*f2 - 0.5*f
                c1 =  1.5*f3 - 2.5*f2 + 1.0
                c2 = -1.5*f3 + 2.0*f2 + 0.5*f
                c3 =  0.5*f3 - 0.5*f2
                W  = (c0 * W_tab[idx - 1] + c1 * W_tab[idx]
                    + c2 * W_tab[idx + 1] + c3 * W_tab[idx + 2])

        return (W * f_src[si]).sum(axis=1)

    # ── Newton solver for the implicit departure equation ─────────────────

    def _newton_euler(self, u_n, x_start, dt):
        """
        Implicit Euler departure:  u_new = u_n(x_start − u_new · dt)

        Residual : r = u_new − interp(u_n, x_D),  x_D = x_start − u_new·dt
        Jacobian : 1 + u_n'(x_D)·dt

        For scalar Burgers this departure is EXACT: u is conserved along
        characteristics, so ξ + u_n(ξ)·dt = x is exactly the Burgers
        characteristic equation.  There is zero temporal discretisation error
        regardless of dt.  The total error is purely interpolation accumulated
        over T/dt steps: O(T · dx^p / dt).
        """
        u_new = self._interp(u_n, x_start)
        for _ in range(self.newton_iter):
            x_D   = x_start - u_new * dt
            res   = u_new - self._interp(u_n, x_D)
            du_dx = self._interp(u_n, x_D, deriv=1)
            J     = 1.0 + du_dx * dt
            J     = np.where(np.abs(J) > 1e-10, J, np.sign(J) * 1e-10)
            u_new = u_new - res / J
        return u_new

    # ── Time step ─────────────────────────────────────────────────────────

    def step(self, dt):
        """
        Both ck=0 and ck=1 use the implicit Euler Newton departure, which is
        the EXACT Burgers characteristic (zero temporal error for all dt).

        The convergence order difference between ck=0 and ck=1 comes entirely
        from the dt scaling used in the convergence study:
          ck=0 : dt ∝ dx^{p/2}   → total error O(dx^{p/2})
          ck=1 : dt ∝ dx^{p/3}   → total error O(dx^{2p/3})  [same as BFECC]

        Why a single Newton suffices for both: the Euler ck=1 scheme gains its
        second-order temporal accuracy through the off-diagonal coupling term
        b·p in the 2×2 Jacobian (R+ couples R−).  For scalar Burgers b=0 so
        that correction vanishes; the implicit Euler departure is already exact,
        and no separate trapezoidal formulation is needed.
        """
        u_new = self._newton_euler(self.u, self.x_ref, dt)
        self.u = u_new

    def compute_dt(self, cfl=1.0):
        return cfl * self.dx / (float(np.max(np.abs(self.u))) + 1e-14)


# ─────────────────────────────────────────────────────────────────────────────
# BurgersSimulation: tie everything together
# ─────────────────────────────────────────────────────────────────────────────

class BurgersSimulation:
    """
    Runs BurgersSolver (GP SL) against exact and spectral references.

    Parameters
    ----------
    xmin, xmax : physical domain
    Nx         : number of cells for the GP solver
    u_0_fn     : callable, initial condition
    T_max      : final time (must be < shock formation time for smooth solution)
    r          : stencil half-width for GP solver
    eps        : RBF-GA shape parameter
    ck         : 0 = Euler SL,  1 = midpoint
    cfl        : CFL number for the GP solver time-step
    newton_iter : Newton iterations per step
    use_spectral : if False, skip building/using the pseudo-spectral reference
                   (saves the DOP853 integration cost; useful for discontinuous
                   or non-periodic-friendly ICs such as the Riemann problem).
    N_spec     : Fourier modes for spectral reference
    D_factor   : extended domain factor for spectral / exact references
    """

    def __init__(
        self,
        xmin, xmax, Nx,
        u_0_fn,
        T_max       = 0.5,
        r           = 2,
        eps         = 1e-3,
        ck          = 1,
        cfl         = 1.0,
        newton_iter = 3,
        weight_mode = 'poly',
        use_spectral = True,
        N_spec      = 4096,
        D_factor    = 4.0,
        verbose     = True,
    ):
        self.xmin    = xmin
        self.xmax    = xmax
        self.Nx      = Nx
        self.T_max   = T_max
        self.cfl     = cfl
        self.verbose = verbose
        self.use_spectral = use_spectral

        dx    = (xmax - xmin) / Nx
        x_ref = np.linspace(xmin + dx / 2, xmax - dx / 2, Nx)
        self.x_ref = x_ref

        # GP semi-Lagrangian solver
        self.solver = BurgersSolver(
            x_ref, u_0_fn,
            r=r, eps=eps, ck=ck, newton_iter=newton_iter,
            weight_mode=weight_mode,
        )

        # Exact characteristics reference
        if verbose:
            print("Building exact characteristics reference ...")
        self.exact = build_burgers_exact(u_0_fn, xmin, xmax, T_max)

        # Spectral reference (optional)
        if use_spectral:
            self.spectral = build_burgers_spectral(
                u_0_fn, xmin, xmax, T_max,
                N_spec=N_spec, D_factor=D_factor, verbose=verbose,
            )
        else:
            self.spectral = None
            if verbose:
                print("Spectral reference disabled (use_spectral=False).")

    def run(self):
        t  = 0.0
        t0 = time.time()
        label = f"BurgersSolver (ck={self.solver.ck})"
        if self.verbose:
            print(f"Running {label} ...")

        while t < self.T_max:
            dt = min(self.solver.compute_dt(self.cfl), self.T_max - t)
            self.solver.step(dt)
            t += dt

        if self.verbose:
            print(f"  Done in {time.time()-t0:.3f} s")

        self._plot_final(t)
        err_ex, err_sp = self._report_errors(t)
        return err_ex, err_sp

    def _report_errors(self, t):
        x      = self.x_ref
        u_num  = self.solver.u
        u_ex   = self.exact.solve(x, t)
        dx     = float(x[1] - x[0])
        err_ex = float(np.sqrt(np.sum((u_num - u_ex)**2) * dx))

        print(f"\nL2 errors at T = {t:.6g}:")
        print(f"  vs exact      : {err_ex:.4e}")

        err_sp = None
        if self.spectral is not None:
            u_sp   = self.spectral.solve(x)
            err_sp = float(np.sqrt(np.sum((u_num - u_sp)**2) * dx))
            print(f"  vs spectral   : {err_sp:.4e}")
        return err_ex, err_sp

    def _plot_final(self, t):
        x      = self.x_ref
        u_num  = self.solver.u
        u_ex   = self.exact.solve(x, t)
        u_sp   = self.spectral.solve(x) if self.spectral is not None else None
        beta   = self.solver.compute_beta()          # (N, R+1)

        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        fig.suptitle(f"Burgers  ∂_t u + u ∂_x u = 0   "
                     f"(T={t:.4g},  Nx={self.Nx},  ck={self.solver.ck})")

        ax = axes[0]
        ax.plot(x, u_ex, 'k-',   lw=1.5, label='Exact')
        if u_sp is not None:
            ax.plot(x, u_sp, 'b--',  lw=1.2, label='Spectral')
        ax.plot(x, u_num,'r+',   ms=4,   label=f'GP SL (ck={self.solver.ck})')
        ax.set_xlabel('x');  ax.set_title('Solution  u(x, T)')
        ax.legend();         ax.grid(True, alpha=0.4)

        ax = axes[1]
        ax.semilogy(x, np.abs(u_num - u_ex) + 1e-17, 'r-', label='|GP − exact|')
        if u_sp is not None:
            ax.semilogy(x, np.abs(u_num - u_sp) + 1e-17, 'b--', label='|GP − spectral|')
        ax.set_xlabel('x');  ax.set_title('Pointwise error')
        ax.legend();         ax.grid(True, alpha=0.4)

        ax = axes[2]
        R1     = beta.shape[1]
        colors = plt.cm.tab10(np.linspace(0, 1, R1))
        for s in range(R1):
            lbl = f'S{s}  (offsets [{s - self.solver.r}..{s}])'
            ax.semilogy(x, beta[:, s] + 1e-17, color=colors[s], lw=1.2, label=lbl)
        ax.semilogy(x, beta.min(axis=1) + 1e-17, 'k--', lw=1.5, label='min β (WENO selector)')
        ax.set_xlabel('x');  ax.set_title('GP-WENO smoothness indicator  β_m')
        ax.legend(fontsize=8); ax.grid(True, alpha=0.4)

        plt.tight_layout()
        plt.show()


# ─────────────────────────────────────────────────────────────────────────────
# Initial condition
# ─────────────────────────────────────────────────────────────────────────────

def u_0(x):
    """Smooth Gaussian bump.  Shock forms at T_shock ≈ sqrt(e)·sigma/(A·sqrt(2))."""
    A, sigma = 1.0, 2.0
    return A - A * np.exp(-x**2 / sigma**2)

# Shock formation time: T_shock = -1 / min(∂_x u_0)
# ∂_x u_0 = -2x/sigma^2 * A * exp(-x^2/sigma^2)
# min at x = sigma/sqrt(2): ∂_x u_0_min = -A*sqrt(2/e)/sigma
# T_shock = sigma / (A*sqrt(2/e)) = sigma*sqrt(e) / (A*sqrt(2))
_T_shock = 2.0 * np.sqrt(np.e) / np.sqrt(2.0)   # ≈ 2.33 for A=1, sigma=2


def make_riemann_ic(u_L=-1.0, u_R=1.0, x0=0.0, width=0.15):
    """
    Build a Riemann-problem initial condition:

        u_0(x) = u_L   for x << x0
        u_0(x) = u_R   for x >> x0

    smoothed over a transition of the given `width` with a tanh profile
    (BurgersExact/BurgersSpectral both assume a smooth IC, so a genuine jump
    discontinuity is regularised rather than used verbatim).

    u_L < u_R  →  rarefaction wave (characteristics diverge, no shock).
    u_L > u_R  →  shock forms immediately at t=0+.

    Returns a callable u_0_fn(x).  Because the far-field states u_L ≠ u_R
    never decay to a common value, the pseudo-spectral reference (which
    assumes periodicity on its extended domain) is not meaningful here —
    run with `use_spectral=False` in BurgersSimulation.
    """
    def u_0_fn(x):
        #return u_L + 0.5 * (u_R - u_L) * (1.0 + np.tanh((x - x0) / width))
        return np.where(x < x0, u_L, u_R)
    return u_0_fn


# Rarefaction wave: u_L < u_R  →  characteristics fan out, stays smooth for all t>0.
u_0_rarefaction = make_riemann_ic(u_L=1., u_R=1.1, x0=0.0, width=0.15)


# ─────────────────────────────────────────────────────────────────────────────
# Run
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    # ── Toggles ────────────────────────────────────────────────────────────
    USE_SPECTRAL = True      # set False to skip the pseudo-spectral reference
    IC           = "rarefaction"  # "gaussian" | "rarefaction"

    if IC == "gaussian":
        print(f"Shock forms at T_shock ≈ {_T_shock:.4f}  (running to T < T_shock)")
        u_0_fn = u_0
        T_max  = 0.8 * _T_shock   # stay in smooth regime
    elif IC == "rarefaction":
        u_0_fn = u_0_rarefaction
        T_max  = 3.
        # Rarefaction: far-field states differ, so the spectral reference's
        # periodicity assumption doesn't apply — force it off.
        USE_SPECTRAL = False
    else:
        raise ValueError(f"Unknown IC {IC!r}")

    sim = BurgersSimulation(
        xmin        = -8.0,
        xmax        =  8.0,
        Nx          = 128,
        u_0_fn      = u_0_fn,
        T_max       = T_max,
        r           = 2,
        eps         = 1e-3,
        ck          = 1,
        cfl         = 1.,
        newton_iter = 1,
        use_spectral = USE_SPECTRAL,
        N_spec      = 4096,
        D_factor    = 4.0,
        verbose     = True,
    )
    sim.run()
