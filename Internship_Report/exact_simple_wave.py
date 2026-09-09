#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
exact_simple_wave.py
--------------------
Exact solution for 1D isentropic Euler via the *right-going simple wave*.

Physics
-------
A right-going simple wave has R- = const everywhere, so only R+ varies.
Each R+ value rides its λ+ = u+c characteristic:

    x(t) = x0 + λ+(x0) * t

Inverting this map gives R+(x,t) = R+(x0), and hence ρ and u everywhere.

Public API
----------
make_simple_wave_ic(x, rho_bg, K, gamma, perturbation)
    Build (rho0, u0) that is a pure right-going simple wave.

SimpleWaveExact(x0, rho0, u0, K, gamma)
    Exact solver.
      .shock_time()          → estimated shock-formation time t*
      .solve(x_eval, t)      → (rho, u) at time t

check_simple_wave(rho0, u0, K, gamma, tol)
    Check whether an IC satisfies R- ≈ const.

How to plug into 1D_Euler_LPTpy.py
------------------------------------
    from exact_simple_wave import SimpleWaveExact, make_simple_wave_ic

    x      = sim.x_lpt
    rho0, u0, _ = make_simple_wave_ic(x, rho_bg=1.5, K=1.0, gamma=1.1,
                                       perturbation=lambda x: 0.3*np.exp(-x**2))
    exact  = SimpleWaveExact(x, rho0, u0, K=1.0, gamma=1.1)

    # after sim.run():
    rho_ex, u_ex = exact.solve(sim.x_lpt, sim.T_max)
    dx = sim.x_lpt[1] - sim.x_lpt[0]
    err_rho = np.sqrt(np.sum((sim.lpt.rho - rho_ex)**2) * dx)
    err_u   = np.sqrt(np.sum((sim.lpt.u   - u_ex  )**2) * dx)
"""

import warnings
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline, make_interp_spline



# ── EOS helpers ───────────────────────────────────────────────────────────────

def _sound_speed(rho, K, gamma):
    return np.sqrt(K * gamma * rho ** (gamma - 1.0))

def _riemann_invariants(rho, u, K, gamma):
    c  = _sound_speed(rho, K, gamma)
    Rp = u + 2.0 * c / (gamma - 1.0)
    Rm = u - 2.0 * c / (gamma - 1.0)
    return Rp, Rm

def _rho_from_char(Rp, Rm_const, K, gamma):
    c = 0.25 * (gamma - 1.0) * (Rp - Rm_const)
    return (c**2 / (K * gamma)) ** (1.0 / (gamma - 1.0))


# ── Public helpers ────────────────────────────────────────────────────────────

def check_simple_wave(rho0, u0, K=1.0, gamma=2.0, tol=1e-4):
    """
    Check whether (rho0, u0) is approximately a right-going simple wave.

    Returns
    -------
    ok        : bool   – True if R- variation is below tol (relative)
    variation : float  – max - min of R- across the grid
    """
    _, Rm = _riemann_invariants(rho0, u0, K, gamma)
    variation = float(np.max(Rm) - np.min(Rm))
    scale     = float(abs(np.mean(Rm))) + 1.0
    return variation < tol * scale, variation


def make_simple_wave_ic(x, rho_bg=0.5, K=1.0, gamma=2.0, perturbation=None,
                        direction='right'):
    """
    Build a simple-wave initial condition.

    direction='right'  (default)
        Right-going wave: R- = const.  dR-/dx = 0  →  A+ = 0.
        Perturbation δ is added to R+.

    direction='left'
        Left-going wave: R+ = const.  dR+/dx = 0  →  A- = 0.
        Perturbation δ is added to R-.
        Use this to isolate the A- CK term for manual verification.

    Parameters
    ----------
    x           : 1D array, spatial grid
    rho_bg      : float, background density
    K, gamma    : EOS constants
    perturbation: callable f(x) → δ(x), or None
    direction   : 'right' or 'left'

    Returns
    -------
    rho0, u0      : arrays, initial density and velocity
    R_const       : float, the constant Riemann invariant
                    (R- for direction='right', R+ for direction='left')
    """
    c_bg = float(_sound_speed(np.array([rho_bg]), K, gamma)[0])

    if direction == 'right':
        Rm_const = -2.0 * c_bg / (gamma - 1.0)
        Rp0 = np.full_like(x, 2.0 * c_bg / (gamma - 1.0))
        if perturbation is not None:
            Rp0 = Rp0 + perturbation(x)
        u0   = 0.5 * (Rp0 + Rm_const)
        c0   = 0.25 * (gamma - 1.0) * (Rp0 - Rm_const)
        rho0 = (c0**2 / (K * gamma)) ** (1.0 / (gamma - 1.0))
        return rho0, u0, Rm_const

    elif direction == 'left':
        Rp_const = 2.0 * c_bg / (gamma - 1.0)
        Rm0 = np.full_like(x, -2.0 * c_bg / (gamma - 1.0))
        if perturbation is not None:
            Rm0 = Rm0 + perturbation(x)
        u0   = 0.5 * (Rp_const + Rm0)
        c0   = 0.25 * (gamma - 1.0) * (Rp_const - Rm0)
        rho0 = (c0**2 / (K * gamma)) ** (1.0 / (gamma - 1.0))
        return rho0, u0, Rp_const

    else:
        raise ValueError(f"direction must be 'right' or 'left', got {direction!r}")


# ── Exact solver ──────────────────────────────────────────────────────────────

class SimpleWaveExact:
    """
    Exact solution for a 1D simple wave.

    direction='right'  (default): R- = const, traces λ+ = u+c characteristics.
    direction='left'             : R+ = const, traces λ- = u-c characteristics.
                                   Use to verify the A- CK term in isolation.

    Parameters
    ----------
    x0         : 1D array, initial positions (sorted)
    rho0, u0   : 1D arrays, initial density and velocity
    K, gamma   : EOS constants
    n_dense    : internal resolution for the characteristic inversion
    direction  : 'right' or 'left'

    Usage
    -----
    solver = SimpleWaveExact(x, rho0, u0, K=1.0, gamma=1.1)
    t_star = solver.shock_time()
    rho_ex, u_ex = solver.solve(x_eval, t=0.5)
    """

    def __init__(self, x0, rho0, u0, K=1.0, gamma=2.0, n_dense=100_000,
                 rho_bg=None, perturbation=None, direction='right'):
        self.K          = K
        self.gamma      = gamma
        self.direction  = direction

        # Extend the dense grid slightly beyond [x0[0], x0[-1]]
        span   = x0[-1] - x0[0]
        x_min  = x0[0]  - 0.05 * span
        x_max  = x0[-1] + 0.05 * span
        self._x0d = np.linspace(x_min, x_max, n_dense)

        if rho_bg is not None and perturbation is not None:
            # Exact IC on the dense grid — no interpolation error at all.
            rho0d, u0d, _ = make_simple_wave_ic(
                self._x0d, rho_bg, K, gamma, perturbation, direction=direction)
        else:
            # Fallback: cubic-spline interpolation from the coarse grid (O(h⁴)).
            rho0d = CubicSpline(x0, rho0)(self._x0d)
            u0d   = CubicSpline(x0, u0)(self._x0d)

        Rp0d, Rm0d = _riemann_invariants(rho0d, u0d, K, gamma)
        c0d        = _sound_speed(rho0d, K, gamma)

        if direction == 'right':
            # R- = const; characteristics are λ+ = u + c
            self._R_var_spline = CubicSpline(self._x0d, Rp0d)  # varying invariant
            self._R_const      = float(np.mean(Rm0d))
            self._lam0d        = u0d + c0d                      # λ+

            variation = float(np.max(Rm0d) - np.min(Rm0d))
            const_name = 'R-'
        else:
            # R+ = const; characteristics are λ- = u - c
            self._R_var_spline = CubicSpline(self._x0d, Rm0d)  # varying invariant
            self._R_const      = float(np.mean(Rp0d))
            self._lam0d        = u0d - c0d                      # λ-

            variation = float(np.max(Rp0d) - np.min(Rp0d))
            const_name = 'R+'

        if variation > 1e-3 * (abs(self._R_const) + 1.0):
            warnings.warn(
                f"{const_name} variation = {variation:.3e}  (not a clean simple wave). "
                "Consider using make_simple_wave_ic() to build consistent ICs. "
                "Exact solution accuracy may be reduced.",
                UserWarning,
                stacklevel=2,
            )

    # ── diagnostics ──────────────────────────────────────────────────────────

    def shock_time(self):
        """
        Estimate the shock-formation time t*.

        Right-going: shock when dλ+/dx0 < 0 (compressive).
        Left-going:  shock when dλ-/dx0 > 0 (compressive going left).
        """
        dlam = np.diff(self._lam0d) / np.diff(self._x0d)
        if self.direction == 'right':
            neg = dlam[dlam < 0.0]
            return np.inf if len(neg) == 0 else 1.0 / (-neg.min())
        else:
            pos = dlam[dlam > 0.0]
            return np.inf if len(pos) == 0 else 1.0 / pos.max()

    # ── main solve ────────────────────────────────────────────────────────────

    def solve(self, x_eval, t):
        """
        Return (rho, u) at positions x_eval and time t.

        Raises ValueError if t >= shock_time().
        """
        t_star = self.shock_time()
        if t >= t_star:
            raise ValueError(
                f"t = {t:.4f} >= shock formation time t* = {t_star:.4f}. "
                "The solution has developed a shock; the exact simple-wave "
                "formula is no longer valid."
            )

        # Forward map: x0 → x(t) along the active characteristic
        x_fwd = self._x0d + self._lam0d * t

        x0_at_x  = CubicSpline(x_fwd, self._x0d)(x_eval)
        Rvar_at_x = self._R_var_spline(x0_at_x)

        if self.direction == 'right':
            Rp_at_x = Rvar_at_x
            u_ex    = 0.5 * (Rp_at_x + self._R_const)
            rho_ex  = _rho_from_char(Rp_at_x, self._R_const, self.K, self.gamma)
        else:
            Rm_at_x = Rvar_at_x
            u_ex    = 0.5 * (self._R_const + Rm_at_x)
            # c = (gamma-1)/4 * (R+ - R-)
            rho_ex  = _rho_from_char(self._R_const, Rm_at_x, self.K, self.gamma)

        return rho_ex, u_ex

    def l2_error(self, x_eval, rho_num, u_num, t):
        """
        Compute discrete L2 errors ||rho_num - rho_exact|| and ||u_num - u_exact||.
        """
        rho_ex, u_ex = self.solve(x_eval, t)
        dx = x_eval[1] - x_eval[0]
        err_rho = float(np.sqrt(np.sum((rho_num - rho_ex)**2) * dx))
        err_u   = float(np.sqrt(np.sum((u_num   - u_ex  )**2) * dx))
        return err_rho, err_u


# ── Convergence study helper ──────────────────────────────────────────────────

def convergence_study(
    xmin, xmax, T_max,
    rho_bg, K, gamma, perturbation,
    Nx_list, solver_cls, solver_kwargs=None,
    interp_method='numpy',
    N_steps=None, CFL=0.9,
):
    """
    Run a spatial convergence study for any LPT-like solver.

    Parameters
    ----------
    solver_cls   : class with __init__(x, rho_0, u_0, K, gamma) and
                   .step(dt, interp_method) and .rho, .u, .x attributes,
                   plus a compute_shock_dt() method.
    solver_kwargs: extra kwargs forwarded to solver_cls.__init__
    Nx_list      : list of grid sizes, e.g. [32, 64, 128, 256]
    N_steps      : if given, fix the number of time steps (T_max per resolution
                   becomes N_steps * dt(Nx)). Avoids different error accumulation
                   across resolutions. T_max is ignored when N_steps is set.
    CFL          : CFL number used to set dt = CFL * dx / max|λ|

    Returns
    -------
    dict with keys 'Nx', 'dx', 'T', 'err_rho', 'err_u', 'order_rho', 'order_u'
    """
    if solver_kwargs is None:
        solver_kwargs = {}

    Nx_list  = list(Nx_list)
    err_rhos = []
    err_us   = []
    T_list   = []

    for Nx in Nx_list:
        dx  = (xmax - xmin) / Nx
        x   = np.linspace(xmin + dx / 2, xmax - dx / 2, Nx)

        rho0, u0, _ = make_simple_wave_ic(x, rho_bg, K, gamma, perturbation)

        rho0_fn = lambda xi, _r=rho0, _x=x: np.interp(xi, _x, _r)
        u0_fn   = lambda xi, _u=u0,   _x=x: np.interp(xi, _x, _u)

        solver = solver_cls(x, rho0_fn, u0_fn, K=K, gamma=gamma, **solver_kwargs)

        if N_steps is not None:
            # Fixed step count: T_max adapts to dx so all resolutions accumulate
            # the same number of interpolation errors.
            c0    = _sound_speed(solver.rho, K, gamma)
            dt0   = CFL * dx / float(np.max(np.abs(solver.u) + c0))
            T_run = N_steps * dt0
        else:
            T_run = T_max
            exact_check = SimpleWaveExact(x, rho0, u0, K=K, gamma=gamma,
                                          rho_bg=rho_bg, perturbation=perturbation)
            t_star = exact_check.shock_time()
            if T_run >= t_star:
                raise ValueError(
                    f"T_max={T_run:.3f} >= shock time t*={t_star:.3f}. "
                    "Reduce T_max or use a smaller perturbation."
                )

        t = 0.0
        while t < T_run:
            dt_cfl = CFL * dx / np.max(np.abs(solver.u) + _sound_speed(solver.rho, K, gamma))
            dt = min(solver.compute_shock_dt() * CFL, dt_cfl, T_run - t)
            solver.step(dt, interp_method=interp_method)
            t += dt

        T_list.append(t)

        # Reference: exact simple-wave solution at t=T_run
        exact = SimpleWaveExact(x, rho0, u0, K=K, gamma=gamma,
                                rho_bg=rho_bg, perturbation=perturbation)
        err_rho, err_u = exact.l2_error(x, solver.rho, solver.u, t)
        err_rhos.append(err_rho)
        err_us.append(err_u)
        if len(err_rhos) > 1:
            prev_dx  = (xmax - xmin) / Nx_list[len(err_rhos) - 2]
            ord_rho  = np.log(err_rhos[-2] / err_rhos[-1]) / np.log(prev_dx / dx)
            ord_u    = np.log(err_us[-2]   / err_us[-1])   / np.log(prev_dx / dx)
            print(f"  Nx={Nx:4d}  dx={dx:.4f}  T={t:.4f}  err_rho={err_rho:.3e} (p={ord_rho:.2f})  err_u={err_u:.3e} (p={ord_u:.2f})")
        else:
            print(f"  Nx={Nx:4d}  dx={dx:.4f}  T={t:.4f}  err_rho={err_rho:.3e}  (  -  )  err_u={err_u:.3e}  (  -  )")

    dxs = [(xmax - xmin) / Nx for Nx in Nx_list]

    def _orders(errs):
        return [
            np.log(errs[i] / errs[i+1]) / np.log(dxs[i] / dxs[i+1])
            for i in range(len(errs) - 1)
        ] + [np.nan]

    return {
        'Nx'       : np.array(Nx_list),
        'dx'       : np.array(dxs),
        'T'        : np.array(T_list),
        'err_rho'  : np.array(err_rhos),
        'err_u'    : np.array(err_us),
        'order_rho': np.array(_orders(err_rhos)),
        'order_u'  : np.array(_orders(err_us)),
    }


def spatial_convergence_study(
    xmin, xmax, T_max,
    rho_0_fn, u_0_fn, K, gamma,
    Nx_list, solver_cls, solver_kwargs=None,
    interp_method='GP', CFL=0.9,
    exact_fn=None,
):
    """
    Spatial convergence study for an arbitrary IC.

    Parameters
    ----------
    rho_0_fn, u_0_fn : callables (x) -> array  — initial condition
    exact_fn         : callable (xi, rho, u, t) -> (err_rho, err_u)
                       Reference solution. For two-wave IC pass a spline built
                       from a fine-grid reference run.
    CFL              : CFL number used to set dt = CFL * dx / max|λ|.
                       Use a small value (e.g. 0.05) so temporal errors stay
                       below the spatial floor across all Nx levels.

    Returns dict with keys 'Nx', 'dx', 'T', 'err_rho', 'err_u',
    'order_rho', 'order_u'.
    """
    if solver_kwargs is None:
        solver_kwargs = {}
    if exact_fn is None:
        raise ValueError("exact_fn is required.")

    Nx_list  = list(Nx_list)
    err_rhos, err_us, T_list = [], [], []

    for Nx in Nx_list:
        dx = (xmax - xmin) / Nx
        x  = np.linspace(xmin + dx / 2, xmax - dx / 2, Nx)

        solver = solver_cls(x, rho_0_fn, u_0_fn, K=K, gamma=gamma, **solver_kwargs)

        t = 0.0
        while t < T_max:
            dt_cfl = CFL * dx / np.max(np.abs(solver.u) + _sound_speed(solver.rho, K, gamma))
            dt = min(dt_cfl, T_max - t)
            solver.step(dt, interp_method=interp_method)
            t += dt
        T_list.append(t)

        err_rho, err_u = exact_fn(x, solver.rho, solver.u, t)
        err_rhos.append(err_rho)
        err_us.append(err_u)

        if len(err_rhos) > 1:
            prev_dx = (xmax - xmin) / Nx_list[len(err_rhos) - 2]
            p_rho   = np.log(err_rhos[-2] / err_rhos[-1]) / np.log(prev_dx / dx)
            p_u     = np.log(err_us[-2]   / err_us[-1])   / np.log(prev_dx / dx)
            print(f"  Nx={Nx:4d}  dx={dx:.4f}  T={t:.4f}  "
                  f"err_rho={err_rho:.3e} (p={p_rho:.2f})  err_u={err_u:.3e} (p={p_u:.2f})")
        else:
            print(f"  Nx={Nx:4d}  dx={dx:.4f}  T={t:.4f}  "
                  f"err_rho={err_rho:.3e}  (  -  )  err_u={err_u:.3e}  (  -  )")

    dxs = [(xmax - xmin) / Nx for Nx in Nx_list]

    def _orders(errs):
        return [
            np.log(errs[i] / errs[i+1]) / np.log(dxs[i] / dxs[i+1])
            for i in range(len(errs) - 1)
        ] + [np.nan]

    return {
        'Nx'       : np.array(Nx_list),
        'dx'       : np.array(dxs),
        'T'        : np.array(T_list),
        'err_rho'  : np.array(err_rhos),
        'err_u'    : np.array(err_us),
        'order_rho': np.array(_orders(err_rhos)),
        'order_u'  : np.array(_orders(err_us)),
    }


def temporal_convergence_study(
    xmin, xmax, Nx, T_max,
    rho_0_fn, u_0_fn, K, gamma,
    dt_list, solver_cls, solver_kwargs=None,
    interp_method='GP',
    dt_ref=None,
    exact_fn=None,
):
    """
    Temporal convergence study: fix Nx, vary dt directly.

    dt_list values are used as the literal time-step size (not CFL × dx/|λ|).

    Reference options (mutually exclusive):
      exact_fn : callable(x, rho, u, t) -> (err_rho, err_u)
                 If provided, errors are measured against this exact solution.
                 This avoids the self-reference bias from accumulated spatial
                 error when dt_ref is too small.
      dt_ref   : float — fallback self-reference; the same solver is run with
                 this dt and used as the baseline.  Choose dt_ref near dt_cross
                 = h^{p/(q+1)} (the minimum-error point), NOT much smaller,
                 otherwise spatial error accumulates over T/dt_ref steps and
                 the reference becomes less accurate than the tested solutions.

    Returns dict with keys 'dt', 'err_rho', 'err_u', 'order_rho', 'order_u'.
    """
    if solver_kwargs is None:
        solver_kwargs = {}
    if exact_fn is None and dt_ref is None:
        raise ValueError("Provide either exact_fn or dt_ref.")

    dx = (xmax - xmin) / Nx
    x  = np.linspace(xmin + dx/2, xmax - dx/2, Nx)

    def _run(dt_val):
        solver = solver_cls(x, rho_0_fn, u_0_fn, K=K, gamma=gamma, **solver_kwargs)
        t = 0.0
        step_dts = []
        while t < T_max:
            dt = min(dt_val, T_max - t)
            solver.step(dt, interp_method=interp_method)
            step_dts.append(dt)
            t += dt
        return solver.rho.copy(), solver.u.copy(), float(np.mean(step_dts))

    if exact_fn is not None:
        rho_ref = u_ref = None   # not used in exact mode
    else:
        print(f"  Computing reference (dt={dt_ref:.4g})...")
        rho_ref, u_ref, _ = _run(dt_ref)

    err_rhos, err_us, dts_used = [], [], []
    for dt_val in dt_list:
        rho_num, u_num, dt_avg = _run(dt_val)
        if exact_fn is not None:
            err_rho, err_u = exact_fn(x, rho_num, u_num, T_max)
        else:
            err_rho = float(np.sqrt(np.sum((rho_num - rho_ref)**2) * dx))
            err_u   = float(np.sqrt(np.sum((u_num   - u_ref  )**2) * dx))
        err_rhos.append(err_rho)
        err_us.append(err_u)
        dts_used.append(dt_avg)
        print(f"  dt={dt_val:.4g}  dt≈{dt_avg:.4f}  err_rho={err_rho:.3e}  err_u={err_u:.3e}")

    def _orders(errs, dts):
        return [
            np.log(errs[i] / errs[i+1]) / np.log(dts[i] / dts[i+1])
            for i in range(len(errs) - 1)
        ] + [np.nan]

    return {
        'dt'       : np.array(dts_used),
        'err_rho'  : np.array(err_rhos),
        'err_u'    : np.array(err_us),
        'order_rho': np.array(_orders(err_rhos, dts_used)),
        'order_u'  : np.array(_orders(err_us,   dts_used)),
    }


def plot_convergence(results, title="Convergence"):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(title)
    for ax, key, label in zip(axes, ['err_rho', 'err_u'], ['ρ', 'u']):
        dx  = results['dx']
        err = results[key]
        ax.loglog(dx, err, 'o-', label='error')
        # reference slopes
        for p, ls in [(1, '--'), (2, ':'), (3, '-.')]:
            ref = err[0] * (dx / dx[0]) ** p
            ax.loglog(dx, ref, ls, color='gray', alpha=0.6, label=f'O(h^{p})')
        ax.set_xlabel('dx')
        ax.set_ylabel(f'L2 error ({label})')
        ax.set_title(label)
        ax.legend()
        ax.grid(True, which='both', alpha=0.3)
    plt.tight_layout()
    plt.show()


