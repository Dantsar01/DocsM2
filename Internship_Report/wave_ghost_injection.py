"""
wave_ghost_injection.py
=======================
Standalone copy of _inject_wave_ghosts() developed in session 2026-06-15.

Injects ghost cells at both boundaries of each tracked rarefaction fan using
linear extrapolation from the corresponding real plateau, then clips at the
plateau value on the tail side to prevent overshoot.

Derivation: ghost_injection_wave_boundaries.pdf (2026-06-15)

Ghost zone layout (both zones to the RIGHT of their boundary):
    lo-zone: [i_lo+1 : i_lo+n+1]  — inside the fan, at x_lo
    hi-zone: [i_hi+1 : i_hi+n+1]  — in the far plateau past x_hi

Slope source (real plateau cells, NOT fan cells):
    lo-zone: cells [i_lo-1, i_lo]         (plateau left of x_lo)
    hi-zone: cells [i_hi+n+1, i_hi+n+2]   (plateau right of ghost zone)

Clipping (tail side only — see pdf for one-sided argument):
    1-fan: hi-zone is the tail → clip R- at R-[i_hi+n+1]; R+ held flat
    2-fan: lo-zone is the tail → clip R+ at R+[i_lo-1];   R- held flat

After setting R±_ghost, recover (u, rho) from EOS:
    u   = (R+ + R-) / 2
    c   = (gamma-1)/4 * (R+ - R-)
    rho = (c² / (K*gamma))^(1/(gamma-1))
"""

import numpy as np


def inject_wave_ghosts(solver, wave_tracker, n_ghost, K, gamma):
    """
    Parameters
    ----------
    solver       : LPTSolver instance — must expose R_plus, R_minus, x, dx, rho, u
    wave_tracker : WaveTracker instance — .waves is a list of dicts with keys
                   'x_lo', 'x_hi', 'family' (1 or 2)
    n_ghost      : int  — number of ghost cells on each side
    K, gamma     : float — isentropic EOS constants
    """
    if not wave_tracker.waves:
        return

    n  = n_ghost
    N  = len(solver.R_plus)
    Rp = solver.R_plus
    Rm = solver.R_minus
    x  = solver.x
    dx = solver.dx

    for wave in wave_tracker.waves:
        i_lo = int(np.searchsorted(x, wave['x_lo']) - 1)
        i_hi = int(np.searchsorted(x, wave['x_hi']) - 1)
        i_lo = int(np.clip(i_lo, 1, N - 2))
        i_hi = int(np.clip(i_hi, 1, N - 2))
        fam  = wave['family']
        if i_hi <= i_lo:
            continue

        # ── lo-zone: [i_lo+1 : i_lo+n+1], inside fan ────────────────────────
        # Slope from the two plateau cells immediately LEFT of x_lo.
        p0 = max(i_lo - 1, 0)
        p1 = i_lo
        sdx  = (x[p1] - x[p0]) if p1 > p0 else dx
        slRp = (Rp[p1] - Rp[p0]) / sdx if sdx > 0 else 0.0
        slRm = (Rm[p1] - Rm[p0]) / sdx if sdx > 0 else 0.0
        g_lo = slice(i_lo + 1, min(i_lo + 1 + n, i_hi + 1))
        if g_lo.start < g_lo.stop:
            d = x[g_lo] - x[p1]          # distances rightward from p1, > 0
            Rp[g_lo] = Rp[p1] + slRp * d
            Rm[g_lo] = Rm[p1] + slRm * d
            if fam == 2:
                # lo-zone is the TAIL for a 2-fan: clip R+ at plateau value
                np.minimum(Rp[g_lo], float(Rp[p1]), out=Rp[g_lo])
                Rm[g_lo] = float(Rm[p1])  # R- is flat in a 2-fan

        # ── hi-zone: [i_hi+1 : i_hi+n+1], far plateau ───────────────────────
        # Slope from the two plateau cells immediately RIGHT of the ghost zone.
        q0 = min(i_hi + 1 + n, N - 2)
        q1 = min(i_hi + 2 + n, N - 1)
        sdx  = (x[q1] - x[q0]) if q1 > q0 else dx
        slRp = (Rp[q1] - Rp[q0]) / sdx if sdx > 0 else 0.0
        slRm = (Rm[q1] - Rm[q0]) / sdx if sdx > 0 else 0.0
        g_hi = slice(i_hi + 1, min(i_hi + 1 + n, N))
        if g_hi.start < g_hi.stop:
            d = x[g_hi] - x[q0]          # distances leftward from q0, < 0
            Rp[g_hi] = Rp[q0] + slRp * d
            Rm[g_hi] = Rm[q0] + slRm * d
            if fam == 1:
                # hi-zone is the TAIL for a 1-fan: clip R- at plateau value
                np.minimum(Rm[g_hi], float(Rm[q0]), out=Rm[g_hi])
                Rp[g_hi] = float(Rp[q0])  # R+ is flat in a 1-fan

        # ── recover (u, rho) from R± in both ghost zones ─────────────────────
        for sl in (g_lo, g_hi):
            if sl.start >= sl.stop:
                continue
            rp = Rp[sl].copy()
            rm = Rm[sl].copy()
            u_g   = 0.5 * (rp + rm)
            c_g   = np.maximum(0.25 * (rp - rm) * (gamma - 1.0), 0.0)
            rho_g = np.where(
                c_g > 0.0,
                (c_g**2 / (K * gamma)) ** (1.0 / (gamma - 1.0)),
                solver.rho[sl],
            )
            solver.u  [sl] = u_g
            solver.rho[sl] = rho_g
