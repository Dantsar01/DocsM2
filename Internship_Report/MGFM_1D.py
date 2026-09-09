#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Usage: python MGFM_1D.py [--Nx 256] [--T_max 3.0] [--cfl 0.9] ...
"""
MGFM_1D.py
==========
Skeleton for the Modified Ghost Fluid Method (MGFM) applied to the 1D
isentropic Euler equations with two fluids.

Reference: GFM_implementation.tex (Liu, Khoo & Yeo 2003 formulation).

Intended use
------------
Instantiate MGFMSimulation, which owns two LPTSolver instances (one per
fluid), a LevelSet, an MGFMInterfaceSolver, and a GhostCellConstructor.
At each time step the driver calls them in order:

    1. LevelSet.locate()            →  x_Gamma, i_Gamma
    2. MGFMInterfaceSolver.solve()  →  p_star, u_star
    3. GhostCellConstructor.build() →  ghost arrays for each fluid
    4. LPTSolver.step()             ×2 (each fluid sees the other's ghosts)
    5. LevelSet.advance()           →  updated phi

Author: dantsar
"""

import os, sys, importlib
import numpy as np
import matplotlib.pyplot as plt

current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(current_dir)

from exact_simple_wave import (
    make_simple_wave_ic, SimpleWaveExact,
    _sound_speed,
)



def load(name):
    if name in sys.modules:
        importlib.reload(sys.modules[name])
    return importlib.import_module(name)

_isl        = load("1D_Euler_RBFGA_ImplicitSL")
LPTSolver   = _isl.LPTSolver
HLLCSolver  = _isl.HLLCSolver


# ── Equation of state helpers ─────────────────────────────────────────────────

def pressure(rho, K, gamma):
    """Isentropic pressure: p = K * rho^gamma."""
    return K * rho**gamma


def sound_speed(rho, K, gamma):
    """Sound speed: c = sqrt(gamma * p / rho), with p = K * rho^gamma."""
    return np.sqrt(gamma * pressure(rho, K, gamma) / rho)


def to_riemann(u, rho, K, gamma):
    """
    Convert (u, rho) → (R+, R-).

    gamma == 1  (linear EOS p = K*rho, c = sqrt(K) = const):
        R± = u ± sqrt(K) * ln(rho)

    gamma != 1  (isentropic EOS p = K*rho^gamma):
        R± = u ± 2c/(gamma-1),  c = sound_speed(rho, K, gamma)
    """
    if np.ndim(gamma) == 0 and float(gamma) == 1.0:
        a = np.sqrt(K)
        dR = a * np.log(rho)
    else:
        c  = sound_speed(rho, K, gamma)
        dR = 2.0 * c / (gamma - 1.0)
    return u + dR, u - dR


def from_riemann(R_plus, R_minus, K, gamma):
    """
    Convert (R+, R-) → (u, rho, c).

    gamma == 1:
        u = (R+ + R-)/2,  rho = exp((R+ - R-)/(2*sqrt(K)))

    gamma != 1:
        u = (R+ + R-)/2,  c = (gamma-1)/4 * (R+ - R-)
        rho = (c^2 / (K*gamma))^(1/(gamma-1))
    """
    u = 0.5 * (R_plus + R_minus)
    if np.ndim(gamma) == 0 and float(gamma) == 1.0:
        a   = np.sqrt(K)
        rho = np.exp((R_plus - R_minus) / (2.0 * a))
        c   = np.full_like(u, a) if np.ndim(u) > 0 else float(a)
    else:
        c   = 0.25 * (gamma - 1.0) * (R_plus - R_minus)
        rho = (c**2 / (K * gamma)) ** (1.0 / (gamma - 1.0))
    return u, rho, c


def K_from_state(rho0, p0, gamma):
    """
    Compute the isentropic constant K = p0 / rho0^gamma from a reference
    state.  Call once during initialisation for each fluid.
    """
    return p0 / rho0**gamma


# ── Shock tracker (replaces level set in 1D) ─────────────────────────────────

class ShockTracker:
    """
    Tracks multiple shocks simultaneously, one per converging characteristic
    family (λ+ → right-going shock, λ− → left-going shock).

    Each shock is a dict  {'x': float, 'i': int}  where 'x' is the current
    position and 'i' is the index of the last cell strictly left of 'x'.

    Backward-compat properties x_Gamma / i_Gamma expose the first shock.
    """

    def __init__(self, x):
        self.x      = x
        self.dx     = x[1] - x[0]
        self.shocks = []          # list of {'x': float, 'i': int}  (shock/shock case)
        self.interfaces = []      # [left, right] dicts (mixed / wave-wave case)

    # ── backward compat ──────────────────────────────────────────────────────

    @property
    def x_Gamma(self):
        return self.shocks[0]['x'] if self.shocks else None

    @property
    def i_Gamma(self):
        return self.shocks[0]['i'] if self.shocks else None

    @property
    def shock_formed(self):
        return len(self.shocks) > 0

    # ── location ─────────────────────────────────────────────────────────────

    def locate_all(self, solver):
        """
        Find one shock per converging characteristic family.

        For each family (λ+, λ−), find the pair (i, i+1) with the minimum
        collision time dt*_i = dx / (λ_i − λ_{i+1}) and predict the shock
        formation point  x = x[i] + λ[i] * dt*_i.

        Replaces self.shocks with the newly found list.
        Returns dt_star (minimum across all found shocks).
        """
        x  = self.x
        dx = self.dx
        lam_p, lam_m = solver._char_speeds(solver.R_plus, solver.R_minus)

        candidates = []
        best_dt   = np.inf

        for lam in (lam_p, lam_m):
            diff = lam[:-1] - lam[1:]
            mask = diff > 0
            if not np.any(mask):
                continue
            idx  = np.where(mask)[0]
            dts  = dx / diff[mask]
            k    = int(np.argmin(dts))
            i    = idx[k]
            dt_k = dts[k]
            x_s  = float(x[i] + lam[i] * dt_k)
            i_s  = i
            candidates.append({'x': x_s, 'i': i_s, 'dt_k': dt_k, 'i_src': i})
            best_dt = min(best_dt, dt_k)

        # De-duplicate: when both families share the same source cell (or adjacent),
        # the smaller-dt_k one is the genuine shock; the other is spurious convergence
        # from the steep gradient of the dominant wave.
        unique = []
        for c in sorted(candidates, key=lambda d: d['dt_k']):
            if not any(abs(c['i_src'] - u['i_src']) <= 1 for u in unique):
                unique.append(c)

        self.shocks = [{'x': c['x'], 'i': c['i']} for c in unique]
        return best_dt

    def locate_from_characteristics(self, solver):
        return self.locate_all(solver)

    def add_new_shocks(self, solver, dt_cfl, shock_tol, min_sep_cells=8):
        """
        During the post-shock phase, scan for new shocks that are not yet tracked.

        Calls locate_all on a temporary copy to get candidate positions, then adds
        only those candidates that are far (> min_sep_cells*dx) from every
        already-tracked shock AND whose dt_k passes the detection threshold.

        Does NOT overwrite existing tracked shocks.
        Returns True if any new shock was added.
        """
        x = self.x; dx = self.dx
        lam_p, lam_m = solver._char_speeds(solver.R_plus, solver.R_minus)

        added = False
        candidates = []
        for lam in (lam_p, lam_m):
            diff = lam[:-1] - lam[1:]
            mask = diff > 0
            if not np.any(mask):
                continue
            idx  = np.where(mask)[0]
            dts  = dx / diff[mask]
            k    = int(np.argmin(dts))
            i    = idx[k]
            dt_k = dts[k]
            if dt_k > shock_tol * dt_cfl:
                continue          # not imminent
            x_s  = float(x[i] + 0.5 * dx)
            i_s  = i
            candidates.append({'x': x_s, 'i': i_s, 'dt_k': dt_k, 'i_src': i})

        # de-duplicate candidates among themselves
        unique_cands = []
        for c in sorted(candidates, key=lambda d: d['dt_k']):
            if not any(abs(c['i_src'] - u['i_src']) <= 1 for u in unique_cands):
                unique_cands.append(c)

        for c in unique_cands:
            already = any(abs(c['x'] - s['x']) < min_sep_cells * dx
                          for s in self.shocks)
            if not already:
                self.shocks.append({'x': c['x'], 'i': c['i']})
                print(f"  [new shock added] x={c['x']:.4f}  i={c['i']}")
                added = True
        return added

    # ── advancement ──────────────────────────────────────────────────────────

    def advance_all_rh(self, rho, u, dt):
        """
        Advance every shock by one step using s = [ρu]/[ρ] evaluated at the
        cells immediately adjacent to each shock (indices i and i+1 in the
        full solver arrays).

        This uses the cells right at the shock face — the SL-evolved state
        with ghost contamination removed — rather than the far-field sampling
        used for the MGFM solve.
        """
        for shock in self.shocks:
            i = shock['i']
            rho_L = float(rho[i])
            rho_R = float(rho[i + 1])
            u_L   = float(u  [i])
            u_R   = float(u  [i + 1])
            if abs(rho_R - rho_L) < 1e-14 * max(rho_L, rho_R):
                s = 0.5 * (u_L + u_R)          # degenerate / smooth region
            else:
                s = (rho_R * u_R - rho_L * u_L) / (rho_R - rho_L)
            shock['x'] += s * dt
            shock['i']  = int(np.searchsorted(self.x, shock['x']) - 1)

    # kept for single-shock compat
    def advance_rh(self, rho_L, u_L, rho_R, u_R, dt):
        if not self.shocks:
            return
        s = (rho_R * u_R - rho_L * u_L) / (rho_R - rho_L)
        self.shocks[0]['x'] += s * dt
        self.shocks[0]['i']  = int(np.searchsorted(self.x, self.shocks[0]['x']) - 1)

    def advance_speed(self, s, dt):
        if not self.shocks:
            return
        self.shocks[0]['x'] += s * dt
        self.shocks[0]['i']  = int(np.searchsorted(self.x, self.shocks[0]['x']) - 1)


# ── Wave (rarefaction fan) tracker ───────────────────────────────────────────

class WaveTracker:
    """
    Tracks rarefaction fans as (x_lo, x_hi) boundary pairs.

    Each fan is a dict:
        {'x_lo': float, 'x_hi': float, 'family': int}

    family=1  →  1-wave (left-going fan), both boundaries advance at λ⁻ = u − c
    family=2  →  2-wave (right-going fan), both boundaries advance at λ⁺ = u + c

    Convention (matching the LaTeX document):
        1-wave: x_lo is the head (undisturbed L state, travels at u_L − c_L),
                x_hi is the tail (plateau * state, travels at u* − c*).
                x_lo < x_hi always (fan spreads to the left and rightward).
        2-wave: x_lo is the tail (plateau * state, travels at u* + c*),
                x_hi is the head (undisturbed R state, travels at u_R + c_R).
                x_lo < x_hi always.
    """

    def __init__(self, x):
        self.x     = x
        self.dx    = x[1] - x[0]
        self.waves = []

    def add_wave(self, x_lo, x_hi, family, s_lo, s_hi):
        """
        Register a new rarefaction fan.

        s_lo, s_hi are the analytically known initial characteristic speeds
        for each boundary (e.g. u_L−c_L for the head, u*−c* for the tail of
        a 1-wave).  They are used on the first advance call so the boundaries
        spread correctly even though both start at the same x_iface at t=0.
        """
        self.waves.append({'x_lo': float(x_lo), 'x_hi': float(x_hi),
                           'family': int(family),
                           's_lo': float(s_lo), 's_hi': float(s_hi)})

    def advance(self, dt):
        """
        Move each fan boundary by one time step.

        For a centered RP (shock-tube IC) the characteristic speeds at the
        head and tail are CONSTANT: u_L−c_L and u*−c* for a 1-wave,
        u*+c* and u_R+c_R for a 2-wave.  The stored s_lo / s_hi are set
        analytically at initialisation and never change.

        Reading speeds back from the solver inside the fan gives wrong results:
        the SL kernel smooths the solution so both boundaries see an
        intermediate state rather than the correct undisturbed or plateau values.
        """
        for wave in self.waves:
            wave['x_lo'] += wave['s_lo'] * dt 
            wave['x_hi'] += wave['s_hi'] * dt
            


'''
# ── Level set (kept for 2D extension) ─────────────────────────────────────────

class LevelSet:
    """
    Tracks the interface position via a scalar field phi defined on the
    shared grid x.

    Sign convention (choose once and stick to it):
        phi < 0  →  fluid 1 (left fluid)
        phi > 0  →  fluid 2 (right fluid)

    The interface is located at the unique zero crossing (single interface
    assumption).

    Parameters
    ----------
    x   : (N,) array  — cell-centre positions
    phi : (N,) array  — initial level-set values

    Attributes set after locate()
    ------------------------------
    x_Gamma : float   — interface position (linear interpolation)
    i_Gamma : int     — index of the last cell with phi < 0, so the
                         interface lies between x[i_Gamma] and x[i_Gamma+1]
    """

    def __init__(self, x, phi):
        self.x   = x
        self.phi = phi.copy()
        self.dx  = x[1] - x[0]

        self.x_Gamma = None
        self.i_Gamma = None

    def locate(self):
        """
        Find x_Gamma and i_Gamma from the current phi.

        The interface is the linear zero of phi between two sign-changing
        neighbours (eq. (xGamma) in the document):

            x_Gamma = x[i] - phi[i] * dx / (phi[i+1] - phi[i])

        Update self.x_Gamma and self.i_Gamma.
        Raise ValueError if no crossing is found.
        """
        crossings = np.where(self.phi[:-1] * self.phi[1:] < 0)[0]
        if len(crossings) == 0:
            raise ValueError("LevelSet.locate: no zero crossing found in phi.")
        if len(crossings) > 1:
            raise ValueError(f"LevelSet.locate: {len(crossings)} crossings found; "
                             "single-interface assumption violated.")

        i = int(crossings[0])
        self.i_Gamma = i
        self.x_Gamma = self.x[i] - self.phi[i] * self.dx / (self.phi[i + 1] - self.phi[i])

    def advance(self, u, dt):
        """
        Advance the level set by one time step using the transport equation

            d phi / dt + u * d phi / dx = 0

        A simple upwind or WENO discretisation in x is appropriate.

        Parameters
        ----------
        u  : (N,) array — velocity field on the grid
        dt : float      — time step
        """
        # TODO: implement (upwind advection is fine for a first pass)
        raise NotImplementedError

    def real_mask(self, fluid_id):
        """
        Return a boolean mask of length N marking the real cells of fluid_id.

            fluid_id == 1  →  phi < 0
            fluid_id == 2  →  phi > 0

        Ghost cells are those just outside the real domain (a few cells wide,
        controlled by the stencil width of LPTSolver).
        """
        # TODO: implement
        raise NotImplementedError

'''
# ── MGFM interface solver ─────────────────────────────────────────────────────

class MGFMInterfaceSolver:
    """
    Solves the two-shock Riemann problem at the interface to predict the
    interface state (p*, u*).

    Exact isentropic Riemann solver (future_ideas.pdf, Section 3.2).

    The residual is

        F(p*) = f_L(p*) + f_R(p*) - (u_L - u_R) = 0

    where each wave function branches on the wave type:

        f_k(p*) = (p* - p_k) / W_k(p*)                          p* > p_k  [shock]
                  2*c_k/(gamma_k-1) * [(p*/p_k)^e_k - 1]        p* ≤ p_k  [rarefaction]

        e_k = (gamma_k - 1) / (2 * gamma_k),  c_k = sqrt(gamma_k * p_k / rho_k)

    The two-shock MGFM is the special case where the top branch is always used.
    Once p* is found, u* is recovered by u* = u_L - f_L(p*), which is exact
    for both wave types.

    Newton's method is used with the acoustic (PVRS) initial guess.

    Parameters
    ----------
    K1, gamma1 : EOS constants for fluid 1
    K2, gamma2 : EOS constants for fluid 2
    tol        : convergence tolerance for Newton (default 1e-12)
    max_iter   : safety cap on Newton iterations (default 20)
    """

    def __init__(self, K1, gamma1, K2, gamma2, tol=1e-12, max_iter=10):
        self.K1       = K1
        self.gamma1   = gamma1
        self.K2       = K2
        self.gamma2   = gamma2
        self.tol      = tol
        self.max_iter = max_iter

    # ── private helpers ──────────────────────────────────────────────────────

    def _rho_I(self, rho_k, p_k, p_star, gamma_k):
        """Isentropic closure (eq. 17): rho_I^k = rho_k * (p*/p_k)^{1/gamma_k}."""
        return rho_k * (p_star / p_k) ** (1.0 / gamma_k)

    def _W(self, rho_k, p_k, p_star, gamma_k):
        """
        RH mass flux magnitude (eq. 14):
            W_k = sqrt( rho_k * rho_I^k * (p* - p_k) / (rho_I^k - rho_k) )

        Acoustic limit (p* -> p_k): W_k -> rho_k * c_k.
        Guard against division by zero when rho_I ≈ rho_k.
        """
        rho_I  = self._rho_I(rho_k, p_k, p_star, gamma_k)
        d_rho  = rho_I - rho_k
        d_p    = p_star - p_k

        if abs(d_rho) < 1e-14 * rho_k:
            # acoustic limit
            K_k = p_k / rho_k ** gamma_k
            return rho_k * np.sqrt(K_k * gamma_k * rho_k ** (gamma_k - 1.0))

        W2 = rho_k * rho_I * d_p / d_rho
        if W2 <= 0.0:
            # guard: should not happen for a genuine compression, but if p* dips
            # below p_k during Newton fall back to acoustic value
            K_k = p_k / rho_k ** gamma_k
            return rho_k * np.sqrt(K_k * gamma_k * rho_k ** (gamma_k - 1.0))

        return np.sqrt(W2)

    def _f_k(self, rho_k, p_k, p_star, gamma_k):
        """
        Exact isentropic wave function f_k(p*).

        Shock branch  (p* > p_k):
            f_k = (p* - p_k) / W_k(p*)
            u*  = u_k ∓ f_k  (minus for left side, plus for right side)

        Rarefaction branch  (p* ≤ p_k):
            f_k = 2*c_k/(gamma_k-1) * [(p*/p_k)^e - 1],  e=(gamma_k-1)/(2*gamma_k)
            Derived from the isentropic Riemann invariant: u ± 2c/(gamma-1) = const.
            f_k ≤ 0 here (expansion lowers pressure and accelerates the fluid).
        """
        c_k = np.sqrt(gamma_k * p_k / rho_k)
        if p_star > p_k:
            return (p_star - p_k) / self._W(rho_k, p_k, p_star, gamma_k)
        else:
            if float(gamma_k) == 1.0:
                a = np.sqrt(p_k / rho_k)          # = sqrt(K), constant c
                return a * np.log(p_star / p_k)   # from R+ = u + a*ln(rho) = const
            e = (gamma_k - 1.0) / (2.0 * gamma_k)
            return 2.0 * c_k / (gamma_k - 1.0) * ((p_star / p_k) ** e - 1.0)

    def _df_k(self, rho_k, p_k, p_star, gamma_k):
        """
        Derivative f_k'(p*) for Newton.

        Shock  (p* > p_k):
            f_k' = 1/(2W) * [1 + rho_k*(p*-p_k) / (gamma_k*p* * (rho_I - rho_k))]
            Acoustic limit: f_k' → 1/(rho_k*c_k).

        Rarefaction  (p* ≤ p_k):
            f_k' = 1/(rho_k*c_k) * (p*/p_k)^(-(gamma_k+1)/(2*gamma_k))
        """
        c_k = np.sqrt(gamma_k * p_k / rho_k)
        if p_star > p_k:
            rho_I = self._rho_I(rho_k, p_k, p_star, gamma_k)
            d_rho = rho_I - rho_k
            W     = self._W(rho_k, p_k, p_star, gamma_k)
            if abs(d_rho) < 1e-14 * rho_k:
                return 1.0 / W
            return 1.0 / (2.0 * W) * (1.0 + rho_k * (p_star - p_k) / (gamma_k * p_star * d_rho))
        else:
            if float(gamma_k) == 1.0:
                a = np.sqrt(p_k / rho_k)   # = sqrt(K)
                return a / p_star           # d/dp* of a*ln(p*/p_k)
            e = -(gamma_k + 1.0) / (2.0 * gamma_k)
            return 1.0 / (rho_k * c_k) * (p_star / p_k) ** e

    def _F(self, p_star, rho_L, u_L, p_L, rho_R, u_R, p_R):
        """
        Exact Riemann residual: F = f_L(p*) + f_R(p*) - (u_L - u_R) = 0.

        F is monotone increasing in p* for both shock and rarefaction branches,
        so Newton converges globally from the PVRS guess.
        """
        fL = self._f_k(rho_L, p_L, p_star, self.gamma1)
        fR = self._f_k(rho_R, p_R, p_star, self.gamma2)
        return fL + fR - (u_L - u_R)

    def _dF(self, p_star, rho_L, u_L, p_L, rho_R, u_R, p_R):
        """F'(p*) = f_L'(p*) + f_R'(p*)."""
        return self._df_k(rho_L, p_L, p_star, self.gamma1) + \
               self._df_k(rho_R, p_R, p_star, self.gamma2)

    def _pvrs_guess(self, rho_L, u_L, p_L, c_L, rho_R, u_R, p_R, c_R):
        """
        Acoustic (PVRS) initial guess: exact solution of the linearised
        Riemann problem (W_k replaced by rho_k * c_k):

            p*_0 = (z_R*p_L + z_L*p_R + z_L*z_R*(u_L - u_R)) / (z_L + z_R)

        where z_k = rho_k * c_k is the acoustic impedance.
        Guarantees 2-3 Newton iterations for moderate pressure ratios.
        """
        z_L = rho_L * c_L
        z_R = rho_R * c_R
        return (z_R * p_L + z_L * p_R + z_L * z_R * (u_L - u_R)) / (z_L + z_R)

    # ── public interface ─────────────────────────────────────────────────────

    def solve(self, rho_L, u_L, p_L, rho_R, u_R, p_R):
        """
        Compute the interface state (p*, u*) via Newton on F(p*) = 0.

        Parameters
        ----------
        rho_L, u_L, p_L : left sampling state  (fluid 1, cell i_Gamma - 1)
        rho_R, u_R, p_R : right sampling state (fluid 2, cell i_Gamma + 2)

        Returns
        -------
        p_star : float
        u_star : float  (recovered from left shock relation, eq. 15)
        """
        K1, g1 = self.K1, self.gamma1
        K2, g2 = self.K2, self.gamma2

        c_L = np.sqrt(K1 * g1 * rho_L ** (g1 - 1.0))
        c_R = np.sqrt(K2 * g2 * rho_R ** (g2 - 1.0))

        p_star = self._pvrs_guess(rho_L, u_L, p_L, c_L, rho_R, u_R, p_R, c_R)
        p_floor = 1e-10 * min(p_L, p_R)
        p_star  = max(p_star, p_floor)

        for _ in range(self.max_iter):
            F  = self._F( p_star, rho_L, u_L, p_L, rho_R, u_R, p_R)
            dF = self._dF(p_star, rho_L, u_L, p_L, rho_R, u_R, p_R)
            dp = F / dF
            p_star -= dp
            p_star  = max(p_star, p_floor)

            if abs(dp) < self.tol * p_star:
                break

        # Recover u* from the left wave relation (exact for both shock and rarefaction):
        #   u* = u_L - f_L(p*)
        # Shock:       f_L = (p*-p_L)/W_L > 0  →  u* < u_L  (deceleration)
        # Rarefaction: f_L < 0                  →  u* > u_L  (acceleration)
        u_star = u_L - self._f_k(rho_L, p_L, p_star, self.gamma1)

        return float(p_star), float(u_star)


# ── Ghost cell constructor ────────────────────────────────────────────────────

class GhostCellConstructor:
    """
    Builds the ghost arrays that extend each fluid's real domain past the
    interface.

    After the MGFM solve we have (p*, u*).  Each fluid extrapolates (p, u)
    from the interface into its ghost region, then recovers the ghost density
    from its own isentropic EOS.

    Three extrapolation levels are available (section 6 of the document):

        'constant'  — p^g = p*,  u^g = u*
        'linear'    — p^g = p* + grad_p * d_g,  etc.
        'rbfga'     — high-order kernel extrapolation (future work)

    Parameters
    ----------
    x        : (N,) array  — grid
    n_ghost  : int         — number of ghost cells to fill on each side
    method   : str         — 'constant' | 'linear' (default 'constant')
    """

    def __init__(self, x, n_ghost, method='constant'):
        self.x       = x
        self.n_ghost = n_ghost
        self.method  = method

    def build(self, fluid_id, i_Gamma, p_star, u_star,
              K_k, gamma_k, rho_real, u_real, p_real):
        """
        Fill ghost arrays for fluid_id past the interface.

        Layout:
            fluid 1  →  real domain is x[: i_Gamma+1],
                         ghost cells are x[i_Gamma+1 : i_Gamma+1+n_ghost]
            fluid 2  →  real domain is x[i_Gamma+1 :],
                         ghost cells are x[i_Gamma+1-n_ghost : i_Gamma+1]

        Steps (see algorithm summary, section 7):
            a. Determine signed distance d_g from each ghost node to x_Gamma.
            b. Extrapolate (p^g, u^g) from (p*, u*) using self.method.
            c. Recover ghost density:
                   rho^g = (p^g / K_k)^{1/gamma_k}
            d. Return ghost arrays (rho_ghost, u_ghost, p_ghost).

        Parameters
        ----------
        fluid_id           : 1 or 2
        i_Gamma            : interface cell index (from LevelSet)
        p_star, u_star     : interface state
        K_k, gamma_k       : EOS constants of fluid_id
        rho_real, u_real, p_real : full real arrays (used for linear extrap)

        Returns
        -------
        rho_ghost, u_ghost, p_ghost : (n_ghost,) arrays
        """
        n = self.n_ghost
        N = len(rho_real)

        # Ghost cell indices and the first real cell outside the ghost zone.
        # fluid 1: ghost = [i+1, i+n],   real boundary = i+n+1
        # fluid 2: ghost = [i+1-n, i],   real boundary = i-n
        if fluid_id == 1:
            ghost_slice = slice(i_Gamma + 1, i_Gamma + 1 + n)
            i_bnd = min(i_Gamma + n + 1, N - 1)
        else:
            ghost_slice = slice(i_Gamma + 1 - n, i_Gamma + 1)
            i_bnd = max(i_Gamma - n, 0)       # same anchor → p_bnd = p* → uniform ghost

        n_cells = ghost_slice.stop - ghost_slice.start

        if self.method == 'constant' or n_cells <= 1:
            p_ghost = np.full(n_cells, p_star)
            u_ghost = np.full(n_cells, u_star)

        elif self.method == 'linear':
            # Linear transition from the real boundary value (no jump at the
            # ghost/real edge) to (p*, u*) at the interface cell.
            # This makes the interpolation kernel see a kink (C0) rather than
            # a step (C-1), eliminating the dominant source of Gibbs oscillations.
            p_bnd = float(p_real[i_bnd])
            u_bnd = float(u_real[i_bnd])
            if fluid_id == 1:
                # k=0 → i+1 (interface side, p*)  ;  k=n-1 → i+n (real boundary, p_bnd)
                t = np.linspace(0.0, 1.0, n_cells)
                p_ghost = p_star  + (p_bnd  - p_star)  * t
                u_ghost = u_star  + (u_bnd  - u_star)  * t
            else:
                # k=0 → i+1-n (real boundary, p_bnd)  ;  k=n-1 → i (interface, p*)
                t = np.linspace(0.0, 1.0, n_cells)
                p_ghost = p_bnd   + (p_star  - p_bnd)   * t
                u_ghost = u_bnd   + (u_star  - u_bnd)   * t

        else:
            raise ValueError(f"GhostCellConstructor: unknown method {self.method!r}")

        rho_ghost = (p_ghost / K_k) ** (1.0 / gamma_k)
        return rho_ghost, u_ghost, p_ghost


# ── Top-level simulation ──────────────────────────────────────────────────────

class MGFMSimulation:
    """
    Orchestrates the full two-fluid MGFM time integration.

    Owns:
        - two LPTSolver instances (fluid 1 and fluid 2)
        - a LevelSet
        - an MGFMInterfaceSolver
        - a GhostCellConstructor

    Time-stepping loop (one call to step()):
        1. Locate interface  (LevelSet.locate)
        2. Sample states one cell back from interface on each side
        3. Solve MGFM Riemann problem  (MGFMInterfaceSolver.solve)
        4. Build ghost cells for each fluid  (GhostCellConstructor.build)
        5. Merge real + ghost arrays and inject into each LPTSolver
        6. Advance each LPTSolver by dt
        7. Advance level set by dt  (LevelSet.advance)

    Parameters
    ----------
    x                       : (N,) array — shared uniform grid
    rho1_0, u1_0            : initial conditions for fluid 1
    rho2_0, u2_0            : initial conditions for fluid 2
    phi_0                   : initial level-set field
    K1, gamma1, K2, gamma2  : EOS constants
    T_max                   : final time
    cfl                     : CFL number (applied to both fluids)
    n_ghost                 : ghost-cell width (should match LPTSolver stencil)
    solver_kw               : keyword dict forwarded to both LPTSolver.__init__
    """

    def __init__(self, x,
                 rho1_0, u1_0,
                 rho2_0, u2_0,
                 K1=1.0, gamma1=2.0,
                 K2=1.0, gamma2=2.0,
                 T_max=1.0,
                 cfl=2.,
                 safety=0.99,
                 n_ghost=2,
                 solver_kw=None):

        solver_kw = solver_kw or {}

        self.K1      = K1;  self.gamma1 = gamma1
        self.K2      = K2;  self.gamma2 = gamma2
        self.T_max   = T_max
        self.cfl     = cfl
        self.safety  = safety
        self.n_ghost = n_ghost

        # Shock tracker (no level set in 1D)
        self.tracker = ShockTracker(x)

        # Interface solver
        self.iface = MGFMInterfaceSolver(K1, gamma1, K2, gamma2)

        # Ghost cell constructor
        self.ghosts = GhostCellConstructor(x, n_ghost, method='linear')

        # One LPTSolver per fluid
        # Note: each solver works on the *full* grid; ghost values overwrite
        # the cells outside the real domain before each step.
        p1_0 = pressure(rho1_0, K1, gamma1)
        p2_0 = pressure(rho2_0, K2, gamma2)

        self.solver1 = LPTSolver(x, rho1_0, u1_0, K=K1, gamma=gamma1,
                                 **solver_kw)
        self.solver2 = LPTSolver(x, rho2_0, u2_0, K=K2, gamma=gamma2,
                                 **solver_kw)

    # ── private helpers ──────────────────────────────────────────────────────

    def _inject_ghosts(self, solver, fluid_id,
                       rho_ghost, u_ghost, p_ghost, i_Gamma):
        """
        Overwrite the ghost region of solver.rho, solver.u (and the
        corresponding Riemann invariants R_plus, R_minus) with the ghost
        values computed by GhostCellConstructor.

        For fluid 1 the ghost region is to the right of i_Gamma.
        For fluid 2 the ghost region is to the left  of i_Gamma.

        Remember to recompute R_plus = u + 2c/(gamma-1) and
        R_minus = u - 2c/(gamma-1) from the new rho, u before the step.
        """
        # TODO: implement
        raise NotImplementedError

    def _compute_dt(self):
        """
        Time step restricted by three constraints:
            1. CFL over both fluids: dt_cfl = cfl * dx / max(|u| + c)
            2. Shock formation time for fluid 1: solver1.compute_shock_dt()
            3. Shock formation time for fluid 2: solver2.compute_shock_dt()

        The shock_dt terms prevent the SL departure point from crossing a
        forming shock within one step (same logic as in EulerSimulation.run).
        """
        dx = self.solver1.dx

        lam1 = np.abs(self.solver1.u) + sound_speed(self.solver1.rho,
                                                     self.K1, self.gamma1)
        lam2 = np.abs(self.solver2.u) + sound_speed(self.solver2.rho,
                                                     self.K2, self.gamma2)
        dt_cfl = self.cfl * dx / max(float(np.max(lam1)), float(np.max(lam2)))

        dt_star1 = self.solver1.compute_shock_dt()
        dt_star2 = self.solver2.compute_shock_dt()

        return min(dt_cfl, self.safety * dt_star1, self.safety * dt_star2)

    # ── public interface ─────────────────────────────────────────────────────

    def step(self, dt):
        """
        Advance the two-fluid system by one time step dt.

        Follow the seven-step algorithm described in the class docstring.
        """
        # 1. Locate interface
        self.ls.locate()
        i_Gamma = self.ls.i_Gamma

        # 2. Sample states (one cell back from interface on each side)
        #    fluid 1: cell i_Gamma - 1  (index i-1 in the document)
        #    fluid 2: cell i_Gamma + 2  (index i+2 in the document)
        # TODO

        # 3. Solve MGFM Riemann problem
        # p_star, u_star = self.iface.solve(...)
        # TODO

        # 4. Build ghost cells
        # rho_g1, u_g1, p_g1 = self.ghosts.build(1, i_Gamma, p_star, u_star,
        #                                         self.K1, self.gamma1, ...)
        # rho_g2, u_g2, p_g2 = self.ghosts.build(2, i_Gamma, p_star, u_star,
        #                                         self.K2, self.gamma2, ...)
        # TODO

        # 5. Inject ghost values into each solver
        # self._inject_ghosts(self.solver1, 1, rho_g1, u_g1, p_g1, i_Gamma)
        # self._inject_ghosts(self.solver2, 2, rho_g2, u_g2, p_g2, i_Gamma)
        # TODO

        # 6. Advance each LPTSolver
        # self.solver1.step(dt)
        # self.solver2.step(dt)
        # TODO

        # 7. Advance level set using the velocity of fluid 1 (or an average)
        # self.ls.advance(self.solver1.u, dt)
        # TODO

    def run(self):
        """
        Main loop.  Computes dt adaptively and calls step() until T_max.
        Returns a list of snapshots for post-processing.
        """
        t = 0.0
        snapshots = []

        while t < self.T_max:
            dt = min(self._compute_dt(), self.T_max - t)
            self.step(dt)
            t += dt
            snapshots.append(self._snapshot(t))

        return snapshots

    def _snapshot(self, t):
        """
        Collect current state into a dict for plotting / error analysis.
        """
        return dict(
            t      = t,
            rho1   = self.solver1.rho.copy(),
            u1     = self.solver1.u.copy(),
            rho2   = self.solver2.rho.copy(),
            u2     = self.solver2.u.copy(),
            x_Gam  = self.tracker.x_Gamma,
            i_Gam  = self.tracker.i_Gamma,
        )


# ── Single-fluid MGFM driver ──────────────────────────────────────────────────

class MGFMSolver1F:
    """
    Single-fluid MGFM driver for a smooth IC that develops a shock.

    Phases
    ------
    Pre-shock  (tracker.shock_formed = False):
        LPTSolver.step() runs normally with no ghost cells.
        ShockTracker.locate_from_characteristics() predicts x_Gamma each step
        but does NOT inject anything.
        The loop stops when dt* <= shock_tol * dx (shock about to form).

    Post-shock (tracker.shock_formed = True):
        Each step:
          1. Sample states one cell back from the shock on each side.
          2. MGFMInterfaceSolver.solve() → (p*, u*).
          3. GhostCellConstructor.build() for left and right ghost regions.
          4. _inject_ghosts() overwrites R±, rho, u in both ghost slices.
          5. LPTSolver.step(dt).
          6. ShockTracker.advance_rh() moves x_Gamma with the RH speed.

    Parameters
    ----------
    x          : (N,) array — uniform grid
    rho_0      : callable or (N,) array — initial density
    u_0        : callable or (N,) array — initial velocity
    K, gamma   : isentropic EOS constants (same fluid on both sides)
    T_max      : final time
    cfl        : CFL number
    safety     : multiplier on dt* (< 1) to avoid exactly hitting shock time
    shock_tol  : threshold: shock declared formed when dt* <= shock_tol * dx
    solver_kw  : extra kwargs forwarded to LPTSolver.__init__
                 n_ghost is set automatically to solver.r after construction
    """

    def __init__(self, x, rho_0, u_0,
                 K=1.0, gamma=2.0,
                 T_max=1.0, cfl=0.9, safety=0.99,
                 shock_tol=1.0,
                 ghost_method='constant',
                 minmax_limiter=False,
                 small_cell_fix=True,
                 solver_kw=None):

        solver_kw = solver_kw or {}

        self.K              = K
        self.gamma          = gamma
        self.T_max          = T_max
        self.cfl            = cfl
        self.safety         = safety
        self.shock_tol      = shock_tol
        self.minmax_limiter = minmax_limiter
        self.small_cell_fix = small_cell_fix

        # LPTSolver expects callables; wrap arrays if needed
        rho_fn = rho_0 if callable(rho_0) else lambda xi, _r=rho_0: _r
        u_fn   = u_0   if callable(u_0)   else lambda xi, _u=u_0:   _u
        self.solver1 = LPTSolver(x, rho_fn, u_fn, K=K, gamma=gamma, **solver_kw)
        self.solver2 = LPTSolver(x, rho_fn, u_fn, K=K, gamma=gamma, **solver_kw)
        self.tracker = ShockTracker(x)

        # n_ghost = 2*r: the outermost real cell (i ± n+1) has a departure point
        # that may land at i ± n (last ghost cell). The stencil of radius r there
        # spans entirely within the ghost zone, preventing the interpolation from
        # crossing the sharp ghost/real boundary and generating Gibbs oscillations.
        self.n_ghost = 2 * self.solver1.r

        # HLLC reference on a finer grid (8× cells, same domain)
        self.hllc = HLLCSolver(x, rho_fn, u_fn, K=K, gamma=gamma, CFL=0.9)

        # Same EOS on both sides (single fluid)
        self.iface  = MGFMInterfaceSolver(K, gamma, K, gamma)
        self.ghosts = GhostCellConstructor(x, self.n_ghost, method=ghost_method)

        # Wave (rarefaction fan) tracker
        self.wave_tracker = WaveTracker(x)

    # ── wave detection ───────────────────────────────────────────────────────

    def detect_waves_from_rp(self, x_iface, t=0.0):
        """
        Sample states on each side of x_iface, solve the RP, and register
        any rarefaction fans with self.wave_tracker.

        Rarefaction criterion: p* < p_k on side k.
            1-wave (left-going fan): p* < p_L
                x_lo = x_iface + (u_L  − c_L ) * t   [head, undisturbed L side]
                x_hi = x_iface + (u*   − c*  ) * t   [tail, plateau * side]
            2-wave (right-going fan): p* < p_R
                x_lo = x_iface + (u*   + c*  ) * t   [tail, plateau * side]
                x_hi = x_iface + (u_R  + c_R ) * t   [head, undisturbed R side]

        At t=0 (RP IC), both boundaries collapse to x_iface and spread
        immediately as the solver advances.
        """
        x = self.solver1.x
        n = self.n_ghost
        N = len(x)
        i   = int(np.searchsorted(x, x_iface) - 1)
        i_L = max(i - n, 0)
        i_R = min(i + n + 1, N - 1)

        rho_L = float(self.solver1.rho[i_L])
        u_L   = float(self.solver1.u[i_L])
        p_L   = float(pressure(rho_L, self.K, self.gamma))
        rho_R = float(self.solver2.rho[i_R])
        u_R   = float(self.solver2.u[i_R])
        p_R   = float(pressure(rho_R, self.K, self.gamma))

        p_star, u_star = self.iface.solve(rho_L, u_L, p_L, rho_R, u_R, p_R)
        rho_star = (p_star / self.K) ** (1.0 / self.gamma)
        print(f"  [IRP]  u*={u_star:.6f}  p*={p_star:.6f}  rho*={rho_star:.6f}")
        c_L    = float(sound_speed(rho_L,    self.K, self.gamma))
        c_R    = float(sound_speed(rho_R,    self.K, self.gamma))
        c_star = float(sound_speed(rho_star, self.K, self.gamma))

        if p_star < p_L:   # 1-rarefaction (left-going fan)
            x_lo = x_iface + (u_L    - c_L   ) * t
            x_hi = x_iface + (u_star - c_star) * t
            # s_lo = head speed (u_L − c_L), s_hi = tail speed (u* − c*)
            self.wave_tracker.add_wave(x_lo, x_hi, family=1,
                                       s_lo=u_L    - c_L,
                                       s_hi=u_star - c_star)
            print(f"  [1-rarefaction] x_lo={x_lo:.4f}  x_hi={x_hi:.4f}  "
                  f"s_lo={u_L-c_L:.4f}  s_hi={u_star-c_star:.4f}")

        if p_star < p_R:   # 2-rarefaction (right-going fan)
            x_lo = x_iface + (u_star + c_star) * t
            x_hi = x_iface + (u_R    + c_R   ) * t
            # s_lo = tail speed (u* + c*), s_hi = head speed (u_R + c_R)
            self.wave_tracker.add_wave(x_lo, x_hi, family=2,
                                       s_lo=u_star + c_star,
                                       s_hi=u_R    + c_R)
            print(f"  [2-rarefaction] x_lo={x_lo:.4f}  x_hi={x_hi:.4f}  "
                  f"s_lo={u_star+c_star:.4f}  s_hi={u_R+c_R:.4f}")

    # ── interface classification ─────────────────────────────────────────────

    def init_interfaces_from_rp(self, x_iface=0.0, t=0.0):
        """
        Solve the interface Riemann problem once and set up whatever tracking
        structure the resulting wave configuration needs:

          shock/shock  →  a single interface, RH-tracked as before
                           (self.tracker.shocks, unchanged path).
          otherwise    →  two interfaces bounding the constant plateau
                           (rho*, u*, p*) (self.tracker.interfaces):
                             - a shock side uses the one-sided RH jump
                               between its far state and the plateau;
                             - a wave side tracks ONLY the tail (the head is
                               not tracked) at the constant characteristic
                               speed u* -/+ c* (1-wave / 2-wave).

        Both speeds are constant because this is a self-similar Riemann
        problem: the star state (and hence every wave speed) does not change
        in time, so there is no need to re-solve the IRP on subsequent steps.

        Also registers the full rarefaction fan(s) (head-to-tail) with
        self.wave_tracker via detect_waves_from_rp, purely for diagnostics /
        plotting -- it plays no role in the ghost injection below.
        """
        x = self.solver1.x
        n = self.n_ghost
        N = len(x)
        i   = int(np.searchsorted(x, x_iface) - 1)
        i_L = max(i - n, 0)
        i_R = min(i + n + 1, N - 1)

        rho_L = float(self.solver1.rho[i_L])
        u_L   = float(self.solver1.u[i_L])
        p_L   = float(pressure(rho_L, self.K, self.gamma))
        rho_R = float(self.solver2.rho[i_R])
        u_R   = float(self.solver2.u[i_R])
        p_R   = float(pressure(rho_R, self.K, self.gamma))

        p_star, u_star = self.iface.solve(rho_L, u_L, p_L, rho_R, u_R, p_R)
        rho_star = (p_star / self.K) ** (1.0 / self.gamma)
        c_star   = float(sound_speed(rho_star, self.K, self.gamma))

        # Cache the plateau state -- it is used as the ghost value at every
        # interface for the whole run.
        self.p_star, self.u_star, self.rho_star, self.c_star = (
            p_star, u_star, rho_star, c_star)

        left_wave  = p_star < p_L    # 1-rarefaction
        right_wave = p_star < p_R    # 2-rarefaction
        i_s = int(np.searchsorted(x, x_iface) - 1)

        if not left_wave and not right_wave:
            print("  [config] shock/shock -> single interface (RH jump)")
            self.tracker.shocks = [{'x': x_iface, 'i': i_s,
                                     'rho_L': rho_L, 'u_L': u_L,
                                     'rho_R': rho_R, 'u_R': u_R}]
            self.tracker.interfaces = []
        else:
            kind_L = 'wave' if left_wave  else 'shock'
            kind_R = 'wave' if right_wave else 'shock'
            print(f"  [config] {kind_L}/{kind_R} -> two interfaces "
                  f"(plateau u*={u_star:.6f}  p*={p_star:.6f}  rho*={rho_star:.6f})")

            s_L = (u_star - c_star if kind_L == 'wave' else
                   (rho_star * u_star - rho_L * u_L) / (rho_star - rho_L))
            s_R = (u_star + c_star if kind_R == 'wave' else
                   (rho_R * u_R - rho_star * u_star) / (rho_R - rho_star))

            print(f"    left  ({kind_L}):  x0={x_iface:.4f}  speed={s_L:.6f}")
            print(f"    right ({kind_R}):  x0={x_iface:.4f}  speed={s_R:.6f}")

            self.tracker.shocks = []
            self.tracker.interfaces = [
                {'x': x_iface, 'i': i_s, 'kind': kind_L, 'side': 'L', 'speed': s_L},
                {'x': x_iface, 'i': i_s, 'kind': kind_R, 'side': 'R', 'speed': s_R},
            ]

        # Register the full fan boundaries (head + tail) for visualization.
        self.detect_waves_from_rp(x_iface, t)

    # ── private helpers ──────────────────────────────────────────────────────

    def _inject_wave_ghosts(self, solver):
        """
        Inject linear-extrapolated + clipped ghost cells at each tracked
        rarefaction fan boundary so the SL stencil doesn't cross the sharp
        plateau edge.  Ported from wave_ghost_injection.py.
        """
        if not self.wave_tracker.waves:
            return
        n   = self.n_ghost
        N   = len(solver.R_plus)
        Rp  = solver.R_plus
        Rm  = solver.R_minus
        x   = solver.x
        dx  = solver.dx

        for wave in self.wave_tracker.waves:
            i_lo = int(np.searchsorted(x, wave['x_lo']) - 1)
            i_hi = int(np.searchsorted(x, wave['x_hi']) - 1)
            i_lo = int(np.clip(i_lo, 1, N - 2))
            i_hi = int(np.clip(i_hi, 1, N - 2))
            fam  = wave['family']
            if i_hi <= i_lo:
                continue

            # lo-zone: [i_lo+1 : i_lo+n+1] — L-plateau extrapolated into fan
            p0 = max(i_lo - 1, 0)
            p1 = i_lo
            sdx  = (x[p1] - x[p0]) if p1 > p0 else dx
            slRp = (Rp[p1] - Rp[p0]) / sdx if sdx > 0 else 0.0
            slRm = (Rm[p1] - Rm[p0]) / sdx if sdx > 0 else 0.0
            g_lo = slice(i_lo + 1, min(i_lo + 1 + n, i_hi + 1))
            if g_lo.start < g_lo.stop:
                d = x[g_lo] - x[p1]
                Rp[g_lo] = Rp[p1] + slRp * d
                Rm[g_lo] = Rm[p1] + slRm * d
                if fam == 2:
                    # lo-zone is the tail for a 2-fan: clip R+ at plateau
                    np.minimum(Rp[g_lo], float(Rp[p1]), out=Rp[g_lo])
                    Rm[g_lo] = float(Rm[p1])

            # hi-zone: [i_hi+1 : i_hi+n+1] — R-plateau extrapolated backward
            q0 = min(i_hi + 1 + n, N - 2)
            q1 = min(i_hi + 2 + n, N - 1)
            sdx  = (x[q1] - x[q0]) if q1 > q0 else dx
            slRp = (Rp[q1] - Rp[q0]) / sdx if sdx > 0 else 0.0
            slRm = (Rm[q1] - Rm[q0]) / sdx if sdx > 0 else 0.0
            g_hi = slice(i_hi + 1, min(i_hi + 1 + n, N))
            if g_hi.start < g_hi.stop:
                d = x[g_hi] - x[q0]
                Rp[g_hi] = Rp[q0] + slRp * d
                Rm[g_hi] = Rm[q0] + slRm * d
                if fam == 1:
                    # hi-zone is the tail for a 1-fan: clip R- at plateau
                    np.minimum(Rm[g_hi], float(Rm[q0]), out=Rm[g_hi])
                    Rp[g_hi] = float(Rp[q0])

            # recover (u, rho) from R± in both ghost zones
            for sl in (g_lo, g_hi):
                if sl.start >= sl.stop:
                    continue
                u_g, rho_g, _ = from_riemann(Rp[sl].copy(), Rm[sl].copy(),
                                              self.K, self.gamma)
                solver.u  [sl] = u_g
                solver.rho[sl] = rho_g

    def _inject_ghosts(self, fluid_id, rho_ghost, u_ghost, i_Gamma):
        """
        Overwrite R±, rho, u in the ghost slice of the appropriate solver.

        fluid_id == 1: right ghost [i_Gamma+1 : i_Gamma+1+n_ghost] → solver1
        fluid_id == 2: left  ghost [i_Gamma+1-n_ghost : i_Gamma+1] → solver2

        R± are recomputed from (rho_ghost, u_ghost) using self.K, self.gamma
        so the solver sees thermodynamically consistent values.
        """
        n  = self.n_ghost
        sl = (slice(i_Gamma + 1, i_Gamma + 1 + n) if fluid_id == 1
              else slice(i_Gamma + 1 - n, i_Gamma + 1))
        s  = self.solver1 if fluid_id == 1 else self.solver2

        s.R_plus[sl], s.R_minus[sl] = to_riemann(u_ghost, rho_ghost,
                                                   self.K, self.gamma)
        s.rho    [sl] = rho_ghost
        s.u      [sl] = u_ghost

    def _compute_dt(self):
        """Return (dt_cfl, dt_star) for the current solver state."""
        dx     = self.solver1.dx
        lam1   = np.abs(self.solver1.u) + sound_speed(self.solver1.rho, self.K, self.gamma)
        lam2   = np.abs(self.solver2.u) + sound_speed(self.solver2.rho, self.K, self.gamma)
        dt_cfl = self.cfl * dx / max(float(np.max(lam1)), float(np.max(lam2)))
        dt_star = min(self.solver1.compute_shock_dt(), self.solver2.compute_shock_dt())
        return dt_cfl, dt_star

    def _sync_solvers(self):
        """Merge real domains after each step: solver1 owns 0..i, solver2 owns i+1..N-1."""
        if not self.tracker.shock_formed:
            np.copyto(self.solver2.rho,     self.solver1.rho)
            np.copyto(self.solver2.u,       self.solver1.u)
            np.copyto(self.solver2.R_plus,  self.solver1.R_plus)
            np.copyto(self.solver2.R_minus, self.solver1.R_minus)
            return
        i = self.tracker.i_Gamma
        self.solver2.rho    [:i+1] = self.solver1.rho    [:i+1]
        self.solver2.u      [:i+1] = self.solver1.u      [:i+1]
        self.solver2.R_plus [:i+1] = self.solver1.R_plus [:i+1]
        self.solver2.R_minus[:i+1] = self.solver1.R_minus[:i+1]
        self.solver1.rho    [i+1:] = self.solver2.rho    [i+1:]
        self.solver1.u      [i+1:] = self.solver2.u      [i+1:]
        self.solver1.R_plus [i+1:] = self.solver2.R_plus [i+1:]
        self.solver1.R_minus[i+1:] = self.solver2.R_minus[i+1:]

    # ── time-stepping ────────────────────────────────────────────────────────

    def step_pre_shock(self, dt):
        """Normal SL step, no ghost cells."""
        self.solver1.step(dt)
        self.solver2.step(dt)
        self._sync_solvers()

    def step_post_shock(self, dt):
        """
        One MGFM-augmented SL step after shock formation.

        Loops over ALL tracked shocks (one per converging characteristic family).
        For each shock:
          1. Sample real states n cells away from the shock on each side.
          2. MGFM Riemann solve → (p*, u*).
          3. Build and inject ghost cells on both sides.
        Then advance the solver once, then advance every shock position with RH.
        """
        n      = self.n_ghost
        N      = len(self.solver1.rho)
        p_full1 = pressure(self.solver1.rho, self.K, self.gamma)
        p_full2 = pressure(self.solver2.rho, self.K, self.gamma)

        # Collect pre-step states for RH advancement (must be done before step)
        rh_states = []
        for shock in self.tracker.shocks:
            i = shock['i']
            i_L = max(i - n, 0)
            i_R = min(i + n + 1, N - 1)

            # Use frozen far-field states if available (avoids startup-error
            # contamination propagating back into the Riemann solve).
            if 'rho_L' in shock:
                rho_L = shock['rho_L'];  u_L = shock['u_L']
                rho_R = shock['rho_R'];  u_R = shock['u_R']
            else:
                rho_L = float(self.solver1.rho[i_L])
                u_L   = float(self.solver1.u  [i_L])
                rho_R = float(self.solver2.rho[i_R])
                u_R   = float(self.solver2.u  [i_R])
            p_L   = pressure(rho_L, self.K, self.gamma)
            p_R   = pressure(rho_R, self.K, self.gamma)

            p_star, u_star = self.iface.solve(rho_L, u_L, p_L, rho_R, u_R, p_R)

            # Mass-flux RH formula: s = [ρu]/[ρ].  Direction-agnostic and exact
            # for any sampling distance.
            d_rho = rho_R - rho_L
            s_rh  = ((rho_R * u_R - rho_L * u_L) / d_rho
                     if abs(d_rho) > 1e-14 * max(rho_L, rho_R)
                     else 0.5 * (u_L + u_R))

            is_compressive = p_star > min(p_L, p_R)

            if is_compressive:
                rho_gL, u_gL, p_gL = self.ghosts.build(
                    1, i, p_star, u_star, self.K, self.gamma,
                    self.solver1.rho, self.solver1.u, p_full1)

                rho_gR, u_gR, p_gR = self.ghosts.build(
                    2, i, p_star, u_star, self.K, self.gamma,
                    self.solver2.rho, self.solver2.u, p_full2)

                self._inject_ghosts(1, rho_gL, u_gL, i)   # → solver1
                self._inject_ghosts(2, rho_gR, u_gR, i)   # → solver2

            rh_states.append((rho_L, u_L, rho_R, u_R, s_rh, p_star, u_star))

        # Step each fluid in its own solver independently.
        # solver1 real domain: 0..i  (post-shock, ghost zone i+1..i+n)
        # solver2 real domain: i+1..N-1 (pre-shock, ghost zone i+1-n..i)
        self.solver1.step(dt)
        self.solver2.step(dt)

        # Small-cell fix: the r cells nearest the interface in each solver have
        # Lagrange stencils that extend into the ghost zone, contaminating them.
        # Replace the last r real cells of solver1 (i-r+1..i) with cell i-r, and
        # the first r real cells of solver2 (i+1..i+r) with cell i+r+1 -- the
        # first genuinely-clean interior cells on each side.
        if self.small_cell_fix:
            _r = self.solver1.r
            for shock in self.tracker.shocks:
                _i = shock['i']
                _src1 = max(_i - _r, 0)
                for _j in range(max(_i - _r + 1, 0), _i + 1):
                    for _s in ('rho', 'u', 'R_plus', 'R_minus'):
                        getattr(self.solver1, _s)[_j] = getattr(self.solver1, _s)[_src1]
                _src2 = min(_i + _r + 1, N - 1)
                for _j in range(_i + 1, min(_i + _r + 1, N)):
                    for _s in ('rho', 'u', 'R_plus', 'R_minus'):
                        getattr(self.solver2, _s)[_j] = getattr(self.solver2, _s)[_src2]

        # Post-step minmax limiter applied to each solver's real interface cells.
        if self.minmax_limiter:
            gam = self.gamma
            for shock, (rho_L, u_L, rho_R, u_R, _, p_star, u_star) in zip(
                    self.tracker.shocks, rh_states):
                i      = shock['i']
                c_L    = sound_speed(rho_L, self.K, gam)
                c_R    = sound_speed(rho_R, self.K, gam)
                rho_st = (p_star / self.K) ** (1.0 / gam)
                c_st   = sound_speed(rho_st, self.K, gam)

                def _clamp(slv, sl, u_lo, u_hi, c_lo, c_hi):
                    Rp = slv.R_plus [sl]
                    Rm = slv.R_minus[sl]
                    u_g, rho_g, c_g = from_riemann(Rp, Rm, self.K, gam)
                    u_g   = np.clip(u_g, u_lo, u_hi)
                    c_g   = np.clip(c_g, c_lo, c_hi)
                    if float(gam) == 1.0:
                        a     = np.sqrt(self.K)
                        rho_g = np.exp((Rp - Rm) / (2.0 * a))
                    else:
                        rho_g = (c_g**2 / (self.K * gam)) ** (1.0 / (gam - 1))
                    slv.u      [sl] = u_g
                    slv.rho    [sl] = rho_g
                    slv.R_plus [sl], slv.R_minus[sl] = to_riemann(u_g, rho_g,
                                                                   self.K, gam)

                r       = self.solver1.r
                n_clamp = r

                # Soft clamp inner real cells near interface
                _clamp(self.solver1, slice(max(i+1-n_clamp, 0), i+1),
                       min(u_L, u_star), max(u_L, u_star),
                       min(c_L, c_st),   max(c_L, c_st))
                _clamp(self.solver2, slice(i+1, min(i+1+n_clamp, N)),
                       min(u_star, u_R), max(u_star, u_R),
                       min(c_st, c_R),   max(c_st, c_R))

        # Merge real domains: solver1 owns 0..i, solver2 owns i+1..N-1
        self._sync_solvers()

        # Advance each shock with the RH speed from its pre-step state.
        # Track old indices so newly-promoted cells can be repaired below.
        old_i_vals = [shock['i'] for shock in self.tracker.shocks]
        for shock, (rho_L, u_L, rho_R, u_R, s_rh, *_) in zip(self.tracker.shocks, rh_states):
            shock['x'] += s_rh * dt
            shock['i']  = int(np.searchsorted(self.tracker.x, shock['x']) - 1)

        # When the shock moves, cells that were just promoted from ghost to real
        # still hold the wrong fluid's state from _sync_solvers, and would
        # contaminate their neighbours on the very next SL step.
        #
        # Rightward (new_i > old_i): cells old_i+1..new_i enter solver1's real
        # domain holding rho_pre (from sync) — overwrite with post-shock state.
        #
        # Leftward  (new_i < old_i): cells new_i+1..old_i enter solver2's real
        # domain holding rho_post (from sync) — overwrite with pre-shock state.
        gam = self.gamma
        for shock, old_i, (rho_L, u_L, rho_R, u_R, *_) in zip(
                self.tracker.shocks, old_i_vals, rh_states):
            new_i = shock['i']
            if new_i > old_i:
                Rp, Rm = to_riemann(u_L, rho_L, self.K, gam)
                for j in range(old_i + 1, min(new_i + 1, N)):
                    self.solver1.rho    [j] = rho_L
                    self.solver1.u      [j] = u_L
                    self.solver1.R_plus [j] = Rp
                    self.solver1.R_minus[j] = Rm
                    self.solver2.rho    [j] = rho_L
                    self.solver2.u      [j] = u_L
                    self.solver2.R_plus [j] = Rp
                    self.solver2.R_minus[j] = Rm
            elif new_i < old_i:
                Rp, Rm = to_riemann(u_R, rho_R, self.K, gam)
                for j in range(max(new_i + 1, 0), old_i + 1):
                    self.solver2.rho    [j] = rho_R
                    self.solver2.u      [j] = u_R
                    self.solver2.R_plus [j] = Rp
                    self.solver2.R_minus[j] = Rm
                    self.solver1.rho    [j] = rho_R
                    self.solver1.u      [j] = u_R
                    self.solver1.R_plus [j] = Rp
                    self.solver1.R_minus[j] = Rm

    def step_post_shock_split(self, dt):
        """
        One step when the interface RP is NOT shock/shock: two interfaces
        (self.tracker.interfaces = [left, right]) bound the constant plateau
        (self.rho_star, self.u_star, self.p_star).

        Each interface is self-similar (constant speed, cached in
        init_interfaces_from_rp), so there is nothing to re-solve here:
          - solver1 gets a ghost zone beyond the LEFT interface filled with
            the plateau star state (exactly like the shock/shock ghost
            injection, just anchored at i_L instead of i_Gamma);
          - solver2 gets a ghost zone beyond the RIGHT interface filled with
            the same plateau star state, anchored at i_R;
          - the plateau strip between the two interfaces (if any cells wide)
            is pinned to the exact star state after the step, since neither
            solver's SL stencil should be trusted there (it may otherwise
            pick up rarefaction-fan or ghost contamination).
        """
        ifL, ifR = self.tracker.interfaces
        i_L, i_R = ifL['i'], ifR['i']
        N = len(self.solver1.rho)

        p_full1 = pressure(self.solver1.rho, self.K, self.gamma)
        p_full2 = pressure(self.solver2.rho, self.K, self.gamma)

        rho_g1, u_g1, _ = self.ghosts.build(
            1, i_L, self.p_star, self.u_star, self.K, self.gamma,
            self.solver1.rho, self.solver1.u, p_full1)
        self._inject_ghosts(1, rho_g1, u_g1, i_L)

        rho_g2, u_g2, _ = self.ghosts.build(
            2, i_R, self.p_star, self.u_star, self.K, self.gamma,
            self.solver2.rho, self.solver2.u, p_full2)
        self._inject_ghosts(2, rho_g2, u_g2, i_R)

        self.solver1.step(dt)
        self.solver2.step(dt)

        # Small-cell fix: only meaningful at a SHOCK interface, where the real
        # cells nearest the boundary have Lagrange stencils that reach into
        # the ghost zone (same contamination issue as step_post_shock's RH
        # shock). A wave interface tracks only the rarefaction tail -- there
        # is no ghost-injected jump there for the stencil to straddle, so it
        # is left untouched.
        if self.small_cell_fix:
            _r = self.solver1.r
            if ifL['kind'] == 'shock':
                _i    = ifL['i']
                _src1 = max(_i - _r, 0)
                for _j in range(max(_i - _r + 1, 0), _i + 1):
                    for _s in ('rho', 'u', 'R_plus', 'R_minus'):
                        getattr(self.solver1, _s)[_j] = getattr(self.solver1, _s)[_src1]
            if ifR['kind'] == 'shock':
                _i    = ifR['i']
                _src2 = min(_i + _r + 1, N - 1)
                for _j in range(_i + 1, min(_i + _r + 1, N)):
                    for _s in ('rho', 'u', 'R_plus', 'R_minus'):
                        getattr(self.solver2, _s)[_j] = getattr(self.solver2, _s)[_src2]

        # Pin the plateau strip to the exact star state in both solvers.
        if i_R > i_L:
            Rp_s, Rm_s = to_riemann(self.u_star, self.rho_star, self.K, self.gamma)
            sl = slice(i_L + 1, i_R + 1)
            for slv in (self.solver1, self.solver2):
                slv.rho[sl]      = self.rho_star
                slv.u[sl]        = self.u_star
                slv.R_plus[sl]   = Rp_s
                slv.R_minus[sl]  = Rm_s

        # Merge: solver1 owns [0..i_L], solver2 owns [i_R+1..N-1], the
        # plateau strip is already identical in both after the pin above.
        self.solver2.rho    [:i_L+1] = self.solver1.rho    [:i_L+1]
        self.solver2.u      [:i_L+1] = self.solver1.u      [:i_L+1]
        self.solver2.R_plus [:i_L+1] = self.solver1.R_plus [:i_L+1]
        self.solver2.R_minus[:i_L+1] = self.solver1.R_minus[:i_L+1]
        self.solver1.rho    [i_R+1:] = self.solver2.rho    [i_R+1:]
        self.solver1.u      [i_R+1:] = self.solver2.u      [i_R+1:]
        self.solver1.R_plus [i_R+1:] = self.solver2.R_plus [i_R+1:]
        self.solver1.R_minus[i_R+1:] = self.solver2.R_minus[i_R+1:]

        # Advance each interface at its own fixed self-similar speed.
        for iface in (ifL, ifR):
            iface['x'] += iface['speed'] * dt
            iface['i']  = int(np.searchsorted(self.tracker.x, iface['x']) - 1)

    # ── main loop ────────────────────────────────────────────────────────────

    def run(self):
        """
        Integrate from t=0 to T_max.  Returns a list of snapshot dicts.

        Single continuous loop.  Shock detection is event-driven: the moment
        dt* drops below shock_tol * dt_cfl, the shock is registered and GFM
        protection is applied on that same step (no unprotected approach step).
        The snap-to-gradient correction is removed; the RH-projected position
        from locate_from_characteristics is used directly.
        """
        t         = 0.0
        snapshots = []

        print(self._compute_dt())
        if self.tracker.interfaces:
            print(f"  [pre-shock skipped] two interfaces pre-initialised: "
                  f"{[(f['kind'], f['x']) for f in self.tracker.interfaces]}")
        elif self.tracker.shock_formed:
            print(f"  [pre-shock skipped] shock pre-initialised at "
                  f"x={[s['x'] for s in self.tracker.shocks]}")

        while t < self.T_max:
            '''
            plt.figure()
            plt.plot(self.solver1.rho)
            plt.show()
            '''
            dt_cfl, dt_star = self._compute_dt()
            imminent = dt_star <= self.shock_tol * dt_cfl

            # Register shock on first detection only; once formed, RH tracking
            # takes over — do NOT overwrite the tracked position.
            if imminent and not self.tracker.shock_formed and not self.tracker.interfaces:
                self.tracker.locate_from_characteristics(self.solver1)
                xs_str = [f"{s['x']:.4f}" for s in self.tracker.shocks]
                print(f"  shock detected at t={t:.6f}  "
                      f"dt*={dt_star:.4e}  x_shocks={xs_str}")
                # Freeze far-field states now (solver still clean at detection).
                _Ng = self.n_ghost
                _Nf = len(self.solver1.rho)
                for _sh in self.tracker.shocks:
                    if 'rho_L' not in _sh:
                        _iL = max(_sh['i'] - _Ng, 0)
                        _iR = min(_sh['i'] + _Ng + 1, _Nf - 1)
                        _sh['rho_L'] = float(self.solver1.rho[_iL])
                        _sh['u_L']   = float(self.solver1.u  [_iL])
                        _sh['rho_R'] = float(self.solver2.rho[_iR])
                        _sh['u_R']   = float(self.solver2.u  [_iR])

            if self.tracker.interfaces:
                # Mixed / wave-wave configuration: two self-similar interfaces.
                dt = min(dt_cfl, self.T_max - t)
                self.step_post_shock_split(dt)
            elif self.tracker.shock_formed:
                # Scan for additional new shocks during the post-shock phase.
                self.tracker.add_new_shocks(
                    self.solver1, dt_cfl, self.shock_tol,
                    min_sep_cells=4 * self.n_ghost)
                # Cap dt to safety*dt_star on formation steps so the shock
                # lands precisely; subsequent steps use dt_cfl freely.
                dt = min(dt_cfl, self.T_max - t)
                if imminent:
                    dt = min(dt, self.safety * dt_star)
                self.step_post_shock(dt)
            else:
                dt = min(dt_cfl, self.safety * dt_star, self.T_max - t)
                self.step_pre_shock(dt)
            
            '''
            plt.figure()
            plt.plot(self.solver1.u)
            plt.show()
            '''
    
            # Advance rarefaction fan boundaries at their characteristic speeds.
            self.wave_tracker.advance(dt)

            t += dt
            snapshots.append(self._snapshot(t))

        self.hllc.run_to(self.T_max)
        return snapshots

    def _snapshot(self, t):
        return dict(
            t            = t,
            x            = self.solver1.x.copy(),
            rho          = self.solver1.rho.copy(),
            u            = self.solver1.u.copy(),
            x_Gamma      = self.tracker.x_Gamma,
            i_Gamma      = self.tracker.i_Gamma,
            x_shocks     = [s['x'] for s in self.tracker.shocks],
            x_interfaces = [{'x': f['x'], 'kind': f['kind'], 'side': f['side']}
                            for f in self.tracker.interfaces],
            waves        = [{'x_lo': w['x_lo'], 'x_hi': w['x_hi'],
                              'family': w['family']}
                             for w in self.wave_tracker.waves],
        )

    def plot_final(self, t):
        """
        Plot rho and u at the final time.

        Shock/shock case: fluid 1 (left, post-shock) is blue, fluid 2 (right,
        pre-shock) is red, one dashed line marks the RH-tracked shock.

        Mixed / wave-wave case: fluid 1 (left of the left interface) is blue,
        the constant plateau strip is green, fluid 2 (right of the right
        interface) is red; two dashed lines mark the two self-similar
        interfaces, labelled by their kind (shock / wave).
        """

        #%matplotlib qt
        x   = self.solver1.x
        N   = len(x)

        shocks     = self.tracker.shocks
        interfaces = self.tracker.interfaces

        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        rho_hllc = np.interp(x, self.hllc.x, self.hllc.rho)
        u_hllc   = np.interp(x, self.hllc.x, self.hllc.u)

        if interfaces:
            ifL, ifR = interfaces
            i_L, i_R = ifL['i'], ifR['i']
            is_fluid1  = np.arange(N) <= i_L
            is_plateau = (np.arange(N) > i_L) & (np.arange(N) <= i_R)
            is_fluid2  = np.arange(N) > i_R

            x_str = (f"L({ifL['kind']})={ifL['x']:.4f}  "
                     f"R({ifR['kind']})={ifR['x']:.4f}")
            fig.suptitle(
                f"MGFM solution at t = {t:.4f}  —  "
                f"blue = fluid 1  |  green = plateau  |  red = fluid 2  |  "
                f"interfaces: {x_str}"
            )

            for ax, y_mgfm, y_hllc, title in [
                (axes[0], self.solver1.rho, rho_hllc, "Density ρ"),
                (axes[1], self.solver1.u,   u_hllc,   "Velocity u"),
            ]:
                ax.plot(x, y_hllc, color="gray", linewidth=1.5,
                        label="HLLC (ref)", zorder=1)
                ax.scatter(x[is_fluid1],  y_mgfm[is_fluid1],
                           color="blue",  marker="+", s=60, zorder=2, label="MGFM fluid 1")
                ax.scatter(x[is_plateau], y_mgfm[is_plateau],
                           color="green", marker="+", s=60, zorder=2, label="MGFM plateau")
                ax.scatter(x[is_fluid2],  y_mgfm[is_fluid2],
                           color="red",   marker="+", s=60, zorder=2, label="MGFM fluid 2")
                ax.axvline(ifL['x'], color="k", linestyle="--", linewidth=1,
                           label=f"L ({ifL['kind']})")
                ax.axvline(ifR['x'], color="k", linestyle=":", linewidth=1,
                           label=f"R ({ifR['kind']})")
                ax.set_title(title)
                ax.set_xlabel("x")
                ax.legend()
                ax.grid(True)
        else:
            x_gam_str = (", ".join(f"{s['x']:.4f}" for s in shocks)
                         if shocks else "not formed")

            # Build a per-cell fluid mask based on shock['i']:
            #   cells 0..i    → fluid 1 real domain (post-shock)
            #   cells i+1..N  → fluid 2 real domain (pre-shock)
            if shocks:
                i_bnd = shocks[0]['i']                # last cell of fluid 1 real domain
            else:
                i_bnd = N // 2 - 1
            is_fluid1 = np.arange(N) <= i_bnd

            fig.suptitle(
                f"MGFM solution at t = {t:.4f}  —  "
                f"blue = fluid 1 (post-shock)  |  red = fluid 2 (pre-shock)  |  "
                f"x_Γ(RH) = [{x_gam_str}]"
            )

            for ax, y_mgfm, y_hllc, title in [
                (axes[0], self.solver1.rho, rho_hllc, "Density ρ"),
                (axes[1], self.solver1.u,   u_hllc,   "Velocity u"),
            ]:
                ax.plot(x, y_hllc, color="gray", linewidth=1.5,
                        label="HLLC (ref)", zorder=1)
                ax.scatter(x[is_fluid1],  y_mgfm[is_fluid1],
                           color="blue",  marker="+", s=60, zorder=2, label="MGFM fluid 1")
                ax.scatter(x[~is_fluid1], y_mgfm[~is_fluid1],
                           color="red",   marker="+", s=60, zorder=2, label="MGFM fluid 2")
                for k, shock in enumerate(shocks):
                    ax.axvline(shock['x'], color="k", linestyle="--", linewidth=1,
                               label=f"x_Γ{k+1}(RH)" if k == 0 else f"x_Γ{k+1}(RH)")
                ax.set_title(title)
                ax.set_xlabel("x")
                ax.legend()
                ax.grid(True)

        plt.tight_layout()
        plt.show()

        #%matplotlib inline


# ── Initial conditions (same as SL code) ─────────────────────────────────────

def _rho_0(x):
    return 1. + np.exp(-x**2 / 4)

def _u_0(x):
    return -0.1 * np.exp(-x**2 / 4)


# ── Shock-tube IC builder ─────────────────────────────────────────────────────
# Right-going shock at x=0 moving at lab speed s = delta.
# Use delta < 0 for a left-going shock; delta = 0 for stationary.
#
# Free parameters: K, gamma, rho_pre, chi = rho_post/rho_pre, delta.
# Velocities follow from mass + momentum RH:
#
#   J       = sqrt( (p_post - p_pre)*rho_post*rho_pre / (rho_post - rho_pre) )
#   u_post  = delta - J / rho_post       (x < 0, compressed side)
#   u_pre   = delta - J / rho_pre        (x > 0, undisturbed side)

def make_shock_tube_ic(K, gamma, rho_pre=1.0, chi=4.0, delta=0.1):
    """
    Exact RH shock-tube states for isentropic Euler p = K*rho^gamma.

    Parameters
    ----------
    K, gamma  : EOS constants
    rho_pre   : undisturbed density (x > 0)
    chi       : compression ratio rho_post / rho_pre  (> 1)
    delta     : desired lab-frame shock speed  (s = delta)
                delta > 0 → right-going, delta < 0 → left-going, 0 → stationary

    Returns
    -------
    rho_post, u_post, rho_pre, u_pre, S
    """
    rho_post = chi * rho_pre
    p_pre    = K * rho_pre  ** gamma
    p_post   = K * rho_post ** gamma
    J        = np.sqrt((p_post - p_pre) * rho_post * rho_pre / (rho_post - rho_pre))
    u_post   = delta - J / rho_post
    u_pre    = delta - J / rho_pre
    c_pre    = np.sqrt(K * gamma * rho_pre ** (gamma - 1.0))
    print(f"[shock IC]  K={K}, gamma={gamma:.4f}, chi={chi}, delta={delta}")
    print(f"  rho_post={rho_post:.6f}  u_post={u_post:.6f}")
    print(f"  rho_pre ={rho_pre:.6f}  u_pre ={u_pre:.6f}")
    print(f"  J={J:.6f}  c_pre={c_pre:.4f}  "
          f"Mach_lab={delta/c_pre:.4f}  Mach_shock={(delta - u_pre)/c_pre:.4f}")
    return rho_post, u_post, rho_pre, u_pre, float(delta)


# ── Custom direct IC (toggle with --direct_ic) ───────────────────────────────
_K_ST, _GAM_ST = 1., 1.4


_S_CUSTOM   = .1
_RHO_L_CUSTOM = 4.0
_U_L_CUSTOM   = _S_CUSTOM - 3.      # = -0.45 with _S_CUSTOM - .5    
_RHO_R_CUSTOM = 1.0
_U_R_CUSTOM   = -2.0 + _S_CUSTOM      # = -1.95

cell_fix = True

# ── Active IC — edit only these three lines ───────────────────────────────────
#_K_ST, _GAM_ST = 1.0, 5/3
_CHI_ST        = 4   # compression ratio rho_post / rho_pre
_DELTA_ST      = .1  # lab-frame shock speed  (0 = stationary, < 0 = left-going)

(_RHO_ST_L, _U_ST_L,
 _RHO_ST_R, _U_ST_R, _S_ST) = make_shock_tube_ic(
    _K_ST, _GAM_ST, rho_pre=1, chi=_CHI_ST, delta=_DELTA_ST)


def _rho_0_st(x):
    """Shock-tube initial density: post-shock for x<0, pre-shock for x≥0."""
    return np.where(x < 0.0, float(_RHO_ST_L), float(_RHO_ST_R))

def _u_0_st(x):
    """Shock-tube initial velocity: exact RH states on each side."""
    return np.where(x < 0.0, float(_U_ST_L), float(_U_ST_R))


# ── CLI ───────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="1D isentropic Euler – MGFM + RBF-GA implicit semi-Lagrangian solver")

    # Grid / time
    parser.add_argument("--Nx",      type=int,   default=256)
    parser.add_argument("--xmin",    type=float, default=-2.0)
    parser.add_argument("--xmax",    type=float, default=2.0)
    parser.add_argument("--T_max",   type=float, default=.5)

    # EOS
    parser.add_argument("--K",       type=float, default=1.0)
    parser.add_argument("--gamma",   type=float, default=2.5)

    # Solver
    parser.add_argument("--cfl",     type=float, default=1.0)
    parser.add_argument("--kernel",  type=str,   default="se")
    parser.add_argument("--radius",  type=int,   default=2)
    parser.add_argument("--ck",      type=int,   default=1,
                        help="0=implicit Euler, 1=implicit trapezoidal departure")
    parser.add_argument("--newton_iter",   type=int,   default=1)
    parser.add_argument("--interp_backend", type=str,  default="gp",
                        help="'lagrange' or 'gp'")
    parser.add_argument("--adaptive_newton",     action="store_true", default=False)
    parser.add_argument("--adaptive_newton_tol", type=float, default=1e-7)

    # Shock detection
    parser.add_argument("--safety",    type=float, default=0.99,
                        help="dt = safety * dt_star before shock forms (default 0.99)")
    parser.add_argument("--shock_tol", type=float, default=2.0,
                        help="shock declared at dt* <= shock_tol * dx  (default 2.0)")

    # IC selection
    parser.add_argument("--shock_tube", action="store_true", default=False,
                        help="Use shock-tube IC (exact RH states, shock at x=0)")
    parser.add_argument("--delta", type=float, default=None,
                        help="Lab-frame shock speed (0=stationary, <0=left-going). "
                             "Overrides _DELTA_ST at the top of the file.")
    parser.add_argument("--chi",   type=float, default=None,
                        help="Compression ratio rho_post/rho_pre. "
                             "Overrides _CHI_ST at the top of the file.")
    parser.add_argument("--minmax", action="store_true", default=False,
                        help="Apply min-max limiter to ghost cells at each shock: "
                             "left ghosts clamped to [u_L, u*], right to [u*, u_R].")
    parser.add_argument("--direct_ic", action="store_true", default=True,
                        help="Use custom direct IC: W_L=(rho=4, u=-2+s), "
                             "W_R=(rho=2, u=s-0.5) with s=0.05, shock at x=0.")


    args = parser.parse_args()

    x = np.linspace(args.xmin, args.xmax, args.Nx)

    solver_kw = dict(
        kernel              = args.kernel,
        r                   = args.radius,
        ck                  = args.ck,
        newton_iter         = args.newton_iter,
        interp_backend      = args.interp_backend,
        adaptive_newton     = args.adaptive_newton,
        adaptive_newton_tol = args.adaptive_newton_tol,
    )

    if args.direct_ic:
        # Force EOS to match the hardcoded IC (linear EOS: p = K*rho, c = sqrt(K))
        if args.K != _K_ST or args.gamma != _GAM_ST:
            print(f"[direct IC]  overriding K={args.K}->{_K_ST}, gamma={args.gamma}->{_GAM_ST} "
                  f"to match IC EOS")
        args.K     = _K_ST
        args.gamma = _GAM_ST
        rho_post, u_post = _RHO_L_CUSTOM, _U_L_CUSTOM
        rho_pre,  u_pre  = _RHO_R_CUSTOM, _U_R_CUSTOM
        S = None
        rho_fn = lambda x, _rL=rho_post, _rR=rho_pre: np.where(x < 0.0, _rL, _rR)
        u_fn   = lambda x, _uL=u_post,   _uR=u_pre:   np.where(x < 0.0, _uL, _uR)
        print(f"[direct IC]  K={_K_ST}  gamma={_GAM_ST}  s={_S_CUSTOM}")
        print(f"  W_L: rho={rho_post}  u={u_post}")
        print(f"  W_R: rho={rho_pre}   u={u_pre}")
    elif args.shock_tube:
        # CLI --delta / --chi override the module-level defaults
        chi   = args.chi   if args.chi   is not None else _CHI_ST
        delta = args.delta if args.delta is not None else _DELTA_ST
        rho_post, u_post, rho_pre, u_pre, S = make_shock_tube_ic(
            args.K, args.gamma, rho_pre=1.0, chi=chi, delta=delta)
        rho_fn = lambda x, _rL=rho_post, _rR=rho_pre: np.where(x < 0.0, _rL, _rR)
        u_fn   = lambda x, _uL=u_post,   _uR=u_pre:   np.where(x < 0.0, _uL, _uR)
    else:
        S      = None
        rho_fn = _rho_0
        u_fn   = _u_0

    sim = MGFMSolver1F(
        x, rho_fn, u_fn,
        K               = args.K,
        gamma           = args.gamma,
        T_max           = args.T_max,
        cfl             = args.cfl,
        safety          = args.safety,
        shock_tol       = args.shock_tol,
        minmax_limiter  = args.minmax,
        #small cell fix only useful when dealing with a left-going transonic (close to 1)
        # Mach Number. Practically, the compression ratio is close to 1, making the MGFM fail.
        # Root cause : one cell contamination error causing intensifying oscillations
        small_cell_fix  = cell_fix,
        solver_kw       = solver_kw,
    )

    if args.direct_ic or args.shock_tube:
        # Solve the IRP once, classify shock/shock vs mixed vs wave/wave, and
        # pre-initialise whichever tracking structure (single shock or the
        # two self-similar interfaces) the configuration needs.
        sim.init_interfaces_from_rp(x_iface=0.0, t=0.0)

    sim.run()
    sim.plot_final(args.T_max)

    if args.shock_tube:
        x_exact = S * args.T_max
        mgfm_x  = sim.tracker.shocks[0]['x'] if sim.tracker.shocks else float('nan')
        hllc_g  = np.abs(np.diff(sim.hllc.rho))
        i_h     = int(np.argmax(hllc_g))
        x_hllc  = 0.5 * (sim.hllc.x[i_h] + sim.hllc.x[i_h + 1])
        print(f"\n── Shock position at T = {args.T_max} ──")
        print(f"  Exact  : {x_exact:.4f}")
        print(f"  MGFM   : {mgfm_x:.4f}   error = {mgfm_x - x_exact:.8f}")
        print(f"  HLLC   : {x_hllc:.4f}   error = {x_hllc - x_exact:.8f}")
