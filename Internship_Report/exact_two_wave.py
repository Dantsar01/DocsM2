#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
exact_two_wave.py
-----------------
Reference solutions for 1D isentropic Euler two-wave ICs.

Two backends:

1. BurgersExactGamma3  (γ=3 ONLY)
   At γ=3 the system decouples into two independent Burgers equations:
       ∂_t R± + R± ∂_x R± = 0
   The exact solution is R±(x,t) = R±_0(ξ) where ξ satisfies
       ξ + R±_0(ξ)·t = x    (solved pointwise by Newton's method).
   Machine-precision accuracy for smooth ICs before shock formation.

   API:
     ref = build_gamma3_exact(rho_0_fn, u_0_fn, K, xmin, xmax, N_ic)
         → BurgersExactGamma3
     ref.solve(x_eval, t)        → (rho, u)
     ref.l2_error(x_eval, rho_num, u_num, t) → (err_rho, err_u)

2. TwoWaveReference  (general γ, pseudo-spectral)
   For general γ the characteristics are coupled — no closed-form solution.
   This backend integrates the conservative Euler equations with a Fourier
   pseudo-spectral spatial discretisation + DOP853 ODE solver.

   Spatial:   Fourier pseudo-spectral + 2/3-dealiasing.
   Temporal:  DOP853 (8th-order), rtol=1e-11, atol=1e-13.
   Combined accuracy for smooth solutions: ~1e-10.

   API:
     ref = build_two_wave_reference(rho_0_fn, u_0_fn, K, gamma,
                                    xmin, xmax, T_target, ...)
         → TwoWaveReference
     ref.solve(x_eval, t=None)              → (rho, u)
     ref.l2_error(x_eval, rho_num, u_num, t) → (err_rho, err_u)
"""

import warnings
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import make_interp_spline, CubicSpline
from scipy.optimize import brentq


# ── Exact Burgers solution for γ=3 ───────────────────────────────────────────

class BurgersExactGamma3:
    """
    Machine-precision reference for γ=3 isentropic Euler.

    At γ=3 the system decouples:  ∂_t R± + R± ∂_x R± = 0  (two Burgers).
    Exact solution: R±(x,t) = R±_0(ξ)  where  ξ + R±_0(ξ)·t = x.
    ξ is found by Newton's method (typically converges in <20 iterations).

    Parameters
    ----------
    x_ic   : 1-D array of IC sample points (uniform, ascending).
    Rp0    : R+(x_ic, 0)
    Rm0    : R-(x_ic, 0)
    K      : EOS constant (p = K ρ^γ, γ=3).
    newton_tol, newton_max : Newton convergence controls.
    """

    gamma = 3.0

    def __init__(self, x_ic, Rp0, Rm0, K,
                 newton_tol=1e-13, newton_max=50):
        self._x_ic       = np.asarray(x_ic,  dtype=float)
        self._Rp0        = np.asarray(Rp0,   dtype=float)
        self._Rm0        = np.asarray(Rm0,   dtype=float)
        self._K          = float(K)
        self._tol        = newton_tol
        self._maxiter    = newton_max

        # Cubic splines for R0 and its derivative — used in Newton iterations.
        self._cs_Rp = CubicSpline(self._x_ic, self._Rp0)
        self._cs_Rm = CubicSpline(self._x_ic, self._Rm0)

    # -- internal: R_0 and its derivative at arbitrary ξ via cubic spline ------
    def _R0_and_dR0(self, xi, which):
        cs = self._cs_Rp if which == '+' else self._cs_Rm
        r  = cs(xi)
        dr = cs(xi, 1)   # first derivative
        return r, dr

    def _solve_xi(self, x_eval, t, which):
        """
        Solve  ξ + R0(ξ)·t = x  for ξ using Brent's method pointwise.

        Brackets:
          R+ (speed > 0, characteristics move right):  ξ ∈ [x - R+_max·t, x]
          R- (speed < 0, characteristics move left):   ξ ∈ [x, x - R-_min·t]
        """
        cs  = self._cs_Rp if which == '+' else self._cs_Rm

        if which == '+':
            R_max = float(cs(self._x_ic).max())
            bracket_lo = float(self._x_ic[0])     # IC left boundary
            def get_bracket(xi):
                lo = max(bracket_lo, xi - R_max * t - 1.0)
                hi = xi                            # f(hi) = R+(xi)·t > 0
                return lo, hi
        else:
            R_min = float(cs(self._x_ic).min())   # most-negative R-
            bracket_hi = float(self._x_ic[-1])    # IC right boundary
            def get_bracket(xi):
                lo = xi                            # f(lo) = R-(xi)·t < 0
                hi = min(bracket_hi, xi - R_min * t + 1.0)
                return lo, hi

        tol  = self._tol
        R0_out = np.empty_like(x_eval)
        for i, x_i in enumerate(x_eval):
            lo, hi = get_bracket(x_i)

            def f(xi):
                return xi + float(cs(xi)) * t - x_i

            # Make sure the bracket is valid (f(lo)<0<f(hi) or vice versa).
            # If f(lo) and f(hi) have the same sign the function is constant
            # in that region → ξ ≈ x (background, no wave).
            flo, fhi = f(lo), f(hi)
            if flo * fhi > 0:
                # No root in bracket: fall back to ξ = x (background region)
                xi_sol = x_i
            else:
                xi_sol = brentq(f, lo, hi, xtol=tol, rtol=tol)

            R0_out[i] = float(cs(xi_sol))
        return R0_out

    def solve(self, x_eval, t):
        """Return (rho, u) at x_eval at time t (machine-precision for smooth IC)."""
        x_eval = np.asarray(x_eval, dtype=float)
        t      = float(t)
        Rp = self._solve_xi(x_eval, t, '+')
        Rm = self._solve_xi(x_eval, t, '-')
        u   = 0.5 * (Rp + Rm)
        c   = 0.5 * (Rp - Rm)          # 2c/(γ-1) = 2c/2 = c  for γ=3
        rho = (c**2 / (self._K * self.gamma)) ** (1.0 / (self.gamma - 1.0))
        return rho, u

    def l2_error(self, x_eval, rho_num, u_num, t):
        """Discrete L2 errors (rho, u) against the exact solution at time t."""
        rho_ex, u_ex = self.solve(x_eval, t)
        dx      = float(x_eval[1] - x_eval[0])
        err_rho = float(np.sqrt(np.sum((rho_num - rho_ex)**2) * dx))
        err_u   = float(np.sqrt(np.sum((u_num   - u_ex  )**2) * dx))
        return err_rho, err_u


def build_gamma3_exact(rho_0_fn, u_0_fn, K, xmin, xmax, T_max=2.0, N_ic=16384):
    """
    Build a BurgersExactGamma3 reference from callable ICs.

    The IC is sampled on an extended domain [xmin - margin, xmax + margin]
    where margin = lambda_max * T_max so that departure points at t ≤ T_max
    always fall inside the IC sample region.

    Parameters
    ----------
    rho_0_fn, u_0_fn : callables  x → array
    K                : EOS constant
    xmin, xmax       : physical domain (reference is valid here)
    T_max            : maximum time at which solve() will be called (default 2.0)
    N_ic             : number of IC sample points on the extended domain

    Returns
    -------
    BurgersExactGamma3 instance
    """
    gamma = 3.0
    # Estimate lambda_max from the physical domain IC
    x_phys = np.linspace(xmin, xmax, 4096)
    rho_phys = rho_0_fn(x_phys)
    u_phys   = u_0_fn(x_phys)
    c_phys   = np.sqrt(K * gamma * rho_phys ** (gamma - 1.0))
    lam_max  = float(np.max(np.abs(u_phys) + c_phys)) + 1.0  # +1 safety margin

    margin   = lam_max * T_max
    x_ic  = np.linspace(xmin - margin, xmax + margin, N_ic)
    rho0  = rho_0_fn(x_ic)
    u0    = u_0_fn(x_ic)
    c0    = np.sqrt(K * gamma * rho0 ** (gamma - 1.0))
    Rp0   = u0 + 2.0 * c0 / (gamma - 1.0)   # = u + c  (γ=3)
    Rm0   = u0 - 2.0 * c0 / (gamma - 1.0)   # = u - c
    return BurgersExactGamma3(x_ic, Rp0, Rm0, K)


# ── Spectral RHS ──────────────────────────────────────────────────────────────

def _make_rhs(N, L, K, gamma):
    """
    Conservative isentropic Euler RHS via Fourier pseudo-spectral + 2/3 rule.

    State vector layout: [rho (N,), rhou (N,)]
    Equations:
        ∂_t ρ   = −∂_x(ρu)
        ∂_t(ρu) = −∂_x(ρu² + K ρ^γ)
    """
    k_arr  = 2.0 * np.pi / L * np.arange(N // 2 + 1, dtype=float)
    cutoff = N // 3   # 2/3 dealiasing: zero Fourier modes with index ≥ N//3

    def rhs(t, state):
        rho  = state[:N]
        rhou = state[N:]

        # Dealias the conserved fields
        rho_hat  = np.fft.rfft(rho);  rho_hat[cutoff:]  = 0.0
        rhou_hat = np.fft.rfft(rhou); rhou_hat[cutoff:] = 0.0
        rho  = np.fft.irfft(rho_hat,  n=N)
        rhou = np.fft.irfft(rhou_hat, n=N)

        u     = rhou / rho
        flux2 = rhou * u + K * rho**gamma   # ρu² + Kρ^γ

        def ddx(f):
            fhat          = np.fft.rfft(f)
            fhat[cutoff:] = 0.0
            return np.fft.irfft(1j * k_arr * fhat, n=N)

        return np.concatenate([-ddx(rhou), -ddx(flux2)])

    return rhs


# ── Reference object ──────────────────────────────────────────────────────────

class TwoWaveReference:
    """
    Pre-computed pseudo-spectral reference for the two-wave isentropic Euler IC.

    Not meant to be instantiated directly — use build_two_wave_reference().
    """

    def __init__(self, x_out, rho_f, u_f, t_target):
        self.t_target = t_target
        self._rho_cs  = make_interp_spline(x_out, rho_f, k=5)
        self._u_cs    = make_interp_spline(x_out, u_f,   k=5)

    def solve(self, x_eval, t=None):
        """
        Return (rho, u) at x_eval via B-spline.

        t is checked for consistency with the built time but is not used
        to re-integrate (the reference is pre-built at t_target).
        """
        if t is not None and not np.isclose(float(t), self.t_target, rtol=1e-5):
            warnings.warn(
                f"Reference was built at t={self.t_target:.6g}; "
                f"requested t={float(t):.6g}.  Returning built solution.",
                UserWarning, stacklevel=2,
            )
        return self._rho_cs(x_eval), self._u_cs(x_eval)

    def l2_error(self, x_eval, rho_num, u_num, t=None):
        """Discrete L2 errors ||rho_num - rho_ref||_2 and ||u_num - u_ref||_2."""
        rho_ex, u_ex = self.solve(x_eval, t)
        dx      = float(x_eval[1] - x_eval[0])
        err_rho = float(np.sqrt(np.sum((rho_num - rho_ex)**2) * dx))
        err_u   = float(np.sqrt(np.sum((u_num   - u_ex  )**2) * dx))
        return err_rho, err_u


# ── Builder ───────────────────────────────────────────────────────────────────

def build_two_wave_reference(
    rho_0_fn, u_0_fn,
    K, gamma,
    xmin, xmax,
    T_target,
    N_spec   = 2048*8,
    D_factor = 4.0,
    rtol     = 1e-11,
    atol     = 1e-13,
    verbose  = True,
):
    """
    Compute a spectral-accuracy reference for a two-wave isentropic Euler IC.

    Parameters
    ----------
    rho_0_fn, u_0_fn : callables  x → array  (initial density and velocity)
    K, gamma         : EOS constants  (p = K ρ^γ)
    xmin, xmax       : physical domain; reference is returned on this interval
    T_target         : target time
    N_spec           : Fourier modes on the extended domain.  Default 2048.
                       For γ close to 1 or large amplitude, increase to 4096.
    D_factor         : domain extension factor.  The periodic solver domain has
                       length D_factor*(xmax-xmin), centred on (xmin+xmax)/2.
                       Must be large enough that no wave reaches the periodic
                       boundary before T_target.  Default 4.0.
    rtol, atol       : DOP853 tolerances.  Defaults give ~1e-10 accuracy.
    verbose          : print progress and diagnostics

    Returns
    -------
    TwoWaveReference instance
    """
    L_phys = xmax - xmin
    cx     = 0.5 * (xmin + xmax)
    L_ext  = D_factor * L_phys
    x0_ext = cx - 0.5 * L_ext

    x_spec = x0_ext + np.arange(N_spec) * (L_ext / N_spec)

    rho0 = rho_0_fn(x_spec)
    u0   = u_0_fn(x_spec)

    # Safety: check that waves won't wrap around before T_target
    c0_max = float(np.sqrt(K * gamma * np.max(rho0)**(gamma - 1.0)))
    v_max  = float(np.max(np.abs(u0))) + c0_max
    reach  = v_max * T_target
    margin = 0.5 * (L_ext - L_phys)   # distance from physical edge to periodic edge
    if reach >= margin:
        warnings.warn(
            f"Waves may reach the periodic boundary before T={T_target}: "
            f"estimated reach={reach:.3f} >= margin={margin:.3f}.  "
            "Increase D_factor.",
            UserWarning, stacklevel=2,
        )

    if verbose:
        print(f"[two-wave ref]  N_spec={N_spec},  L_ext={L_ext:.1f},  "
              f"v_max≈{v_max:.3f},  reach≈{reach:.3f},  margin={margin:.3f}")

    rhs    = _make_rhs(N_spec, L_ext, K, gamma)
    state0 = np.concatenate([rho0, rho0 * u0])

    if verbose:
        print(f"  Integrating t ∈ [0, {T_target}] with DOP853 "
              f"(rtol={rtol:.0e}, atol={atol:.0e}) ...")

    sol = solve_ivp(
        rhs,
        [0.0, T_target],
        state0,
        method       = 'DOP853',
        rtol         = rtol,
        atol         = atol,
        dense_output = False,
    )

    if not sol.success:
        raise RuntimeError(f"DOP853 failed: {sol.message}")

    rho_f = sol.y[:N_spec, -1]
    u_f   = sol.y[N_spec:, -1] / rho_f

    # Spectral resolution check: Fourier amplitudes near the dealiasing cutoff
    rho_hat    = np.fft.rfft(rho_f)
    amp_max    = float(np.abs(rho_hat).max())
    cutoff     = N_spec // 3
    amp_tail   = float(np.abs(rho_hat[cutoff - 10 : cutoff]).mean())
    tail_ratio = amp_tail / amp_max if amp_max > 0 else np.nan

    if verbose:
        print(f"  nfev={sol.nfev},  t_final={sol.t[-1]:.10f}")
        print(f"  Spectral tail check: |ρ̂_tail| / |ρ̂_max| = {tail_ratio:.2e}  "
              f"({'OK' if tail_ratio < 1e-8 else 'WARN: increase N_spec'})")

    if tail_ratio > 1e-6:
        warnings.warn(
            f"Fourier tail is large (ratio={tail_ratio:.2e}): solution may not be "
            "fully resolved.  Increase N_spec.",
            UserWarning, stacklevel=2,
        )

    # Restrict to physical domain + small buffer for the B-spline
    buf   = 0.1 * L_phys
    mask  = (x_spec >= xmin - buf) & (x_spec <= xmax + buf)
    x_out = x_spec[mask]

    return TwoWaveReference(x_out, rho_f[mask], u_f[mask], float(sol.t[-1]))


# ── Standalone test and plot ──────────────────────────────────────────────────

if __name__ == '__main__':
    import matplotlib.pyplot as plt

    K, gamma   = 1.0, 1.4
    xmin, xmax = -6.0, 6.0
    T_target   = 1.0

    # IC from run_spatial_convergence.py  (section B)
    def rho_0(x): return 1.5 + np.exp(-x**2 / 4.0)
    def u_0(x):   return -0.2 * np.exp(-x**2 / 4.0)

    ref = build_two_wave_reference(
        rho_0, u_0,
        K=K, gamma=gamma,
        xmin=xmin, xmax=xmax,
        T_target=T_target,
        N_spec=2048,
        D_factor=4.0,
        verbose=True,
    )

    x_plot  = np.linspace(xmin, xmax, 4000)
    rho_ex, u_ex = ref.solve(x_plot)
    x0 = x_plot

    # ── Riemann invariants ────────────────────────────────────────────────────
    def riemann_invariants(rho, u, K, gamma):
        c  = np.sqrt(K * gamma * rho**(gamma - 1.0))
        Rp = u + 2.0 * c / (gamma - 1.0)
        Rm = u - 2.0 * c / (gamma - 1.0)
        return Rp, Rm

    Rp0, Rm0 = riemann_invariants(rho_0(x0), u_0(x0), K, gamma)
    Rp_ex, Rm_ex = riemann_invariants(rho_ex, u_ex, K, gamma)

    # ── Fourier spectrum of final ρ (resolution diagnostic) ──────────────────
    rho_hat = np.fft.rfft(ref._rho_cs(x_plot))
    freqs   = np.arange(len(rho_hat))

    fig, axes = plt.subplots(2, 2, figsize=(13, 8))
    fig.suptitle(
        f"Two-wave pseudo-spectral reference  "
        f"(γ={gamma}, K={K}, T={T_target})",
        fontsize=13,
    )

    ax = axes[0, 0]
    ax.plot(x0,     rho_0(x0), '--', color='gray', alpha=0.6, label='IC (t=0)')
    ax.plot(x_plot, rho_ex,    '-',  color='C0',              label=f't={T_target}')
    ax.set_title('Density ρ'); ax.set_xlabel('x')
    ax.legend(); ax.grid(True, alpha=0.4)

    ax = axes[0, 1]
    ax.plot(x0,     u_0(x0), '--', color='gray', alpha=0.6, label='IC (t=0)')
    ax.plot(x_plot, u_ex,    '-',  color='C1',              label=f't={T_target}')
    ax.set_title('Velocity u'); ax.set_xlabel('x')
    ax.legend(); ax.grid(True, alpha=0.4)

    ax = axes[1, 0]
    ax.plot(x0,     Rp0,    '--', color='gray',  alpha=0.5, label='R+ IC')
    ax.plot(x_plot, Rp_ex,  '-',  color='C2',               label='R+ final')
    ax.plot(x0,     Rm0,    ':',  color='gray',  alpha=0.5, label='R− IC')
    ax.plot(x_plot, Rm_ex,  '--', color='C3',               label='R− final')
    ax.set_title('Riemann invariants R±  (both vary → two-wave)'); ax.set_xlabel('x')
    ax.legend(); ax.grid(True, alpha=0.4)

    ax = axes[1, 1]
    ax.semilogy(freqs, np.abs(rho_hat) + 1e-16, color='C0')
    ax.axvline(4000 // 3, color='red', linestyle='--', alpha=0.6, label='2/3 cutoff')
    ax.set_title('Fourier spectrum |ρ̂_k|  (resolution check)')
    ax.set_xlabel('mode k'); ax.legend(); ax.grid(True, alpha=0.4)

    plt.tight_layout()
    plt.savefig('two_wave_reference.png', dpi=150, bbox_inches='tight')
    print("  Saved two_wave_reference.png")
    plt.show()

    # ── Usage hint ────────────────────────────────────────────────────────────
    print("""
Usage in run_spatial_convergence.py:

    from exact_two_wave import build_two_wave_reference

    ref_tw = build_two_wave_reference(
        rho_0_tw, u_0_tw,
        K=K, gamma=gamma,
        xmin=xmin, xmax=xmax,
        T_target=T_tw,
    )

    # drop-in replacement for exact_fn_tw
    exact_fn_tw = ref_tw.l2_error
""")
