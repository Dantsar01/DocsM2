import numpy as np
from gp_interp import gp_interpolate


def fritsch_carlson_monotone(R, G, h_all):
    """
    Vectorized Fritsch-Carlson (1980) monotone derivative adjustment.

    Modifies derivative estimates G so that the cubic Hermite built from
    (R, G) is monotone on every interval.  Three steps:
      1. Flat intervals  → zero both endpoint derivatives.
      2. Sign check      → zero any derivative that points away from the
                           chord slope of a neighbouring interval.
      3. Circle clamp    → project (alpha, beta) = (G_j, G_{j+1}) / delta
                           onto the disk alpha^2 + beta^2 <= 9, which is a
                           sufficient condition for monotonicity.

    Parameters
    ----------
    R     : (N,)      function values at nodes
    G     : (N,)      derivative estimates (e.g. from GP)
    h_all : (N-1,) or scalar  cell widths  x[k+1] - x[k]

    Returns
    -------
    G_fc : (N,) adjusted derivatives
    """
    n    = len(R)
    G_fc = G.copy().astype(float)
    h    = np.broadcast_to(h_all, n - 1).copy() if np.ndim(h_all) == 0 \
           else np.asarray(h_all, dtype=float)

    delta = (R[1:] - R[:-1]) / h          # chord slopes, shape (n-1,)

    # --- Step 1: flat intervals ---
    flat = delta == 0.0
    G_fc[:-1] = np.where(flat, 0.0, G_fc[:-1])
    G_fc[1:]  = np.where(flat, 0.0, G_fc[1:])

    # --- Step 2: sign check ---
    # node 0: only right interval
    if G_fc[0] * delta[0] < 0.0:
        G_fc[0] = 0.0
    # interior nodes: must agree with BOTH neighbouring chord slopes
    if n > 2:
        dl  = delta[:-1]                   # left  chord for nodes 1..n-2
        dr  = delta[1:]                    # right chord for nodes 1..n-2
        bad = (G_fc[1:-1] * dl < 0.0) | (G_fc[1:-1] * dr < 0.0)
        G_fc[1:-1] = np.where(bad, 0.0, G_fc[1:-1])
    # node n-1: only left interval
    if G_fc[-1] * delta[-1] < 0.0:
        G_fc[-1] = 0.0

    # --- Step 3: monotone-region constraint (circle clamp) ---
    safe  = np.where(flat, 1.0, delta)     # avoid division by zero
    alpha = G_fc[:-1] / safe              # (n-1,)
    beta  = G_fc[1:]  / safe              # (n-1,)

    tau     = alpha**2 + beta**2
    outside = (~flat) & (tau > 9.0)
    phi_int = np.where(outside, 3.0 / np.sqrt(np.where(outside, tau, 1.0)), 1.0)

    # Each node belongs to up to two intervals; take the tightest scaling.
    phi_node = np.ones(n)
    np.minimum.at(phi_node, np.arange(n - 1), phi_int)   # left  endpoints
    np.minimum.at(phi_node, np.arange(1, n),  phi_int)   # right endpoints

    G_fc *= phi_node
    return G_fc


def cip_step(x, R, lam, dt, r, nugget=0, n_nearest=3,
             monotone=False, rbga_eps = 1e-3):
    """
    One CIP timestep: advect R along characteristics and return R*.

    Parameters
    ----------
    x               : (N,) node positions (non-uniform allowed)
    R               : (N,) Riemann invariant values at x
    lam             : (N,) characteristic speed at each node (e.g. u+c, u-c, or u)
    dt              : float - timestep
    r               : float - GP stencil radius multiplier
    nugget          : float - GP regularisation (default 0)
    n_nearest       : int   - GP fallback stencil size (default 3)
    sigma_scale     : float - sigma_k = sigma_scale * nearest-neighbour distance
    monotone        : bool  - if True, apply Fritsch-Carlson adjustment to GP
                              gradients before building the Hermite polynomial,
                              guaranteeing monotonicity on each cell (default False)
    midpoint_tracing: bool  - if True, use midpoint (2nd-order) characteristic
                              tracing instead of straight-line (1st-order).
                              Reduces departure-point error from O(dt²) to O(dt³),
                              unblocking the cubic interpolant for 3rd-order
                              convergence on nonlinear problems (default False).

    Returns
    -------
    R_star : (N,) interpolated values at departure points
    """
    N = len(x)

    # ------------------------------------------------------------------
    # 1. Gradient G = dR/dx at every node via GP derivative kernel
    # ------------------------------------------------------------------
    if rbga_eps == 0:
        G = gp_interpolate(x, R, x, r, derivative=True)
    else :
        G = gp_interpolate(x, R, x, r, derivative=True, rbfga_eps = rbga_eps)
    # ------------------------------------------------------------------
    # 1b. Optional Fritsch-Carlson monotone adjustment
    # ------------------------------------------------------------------
    if monotone:
        G = fritsch_carlson_monotone(R, G, x[1:] - x[:-1])

    # ------------------------------------------------------------------
    # 2. Departure points
    # ------------------------------------------------------------------

    x_dep = np.clip(x - lam * dt, x[0], x[-1])

    j   = np.searchsorted(x, x_dep, side='right') - 1
    j   = np.clip(j, 0, N - 2)

    # Local cell width — varies per cell when x is non-uniform
    h   = x[j + 1] - x[j]
    xi  = x_dep - x[j]                  # in [0, h]

    R_j   = R[j];   R_jp1 = R[j + 1]
    G_j   = G[j];   G_jp1 = G[j + 1]

    # ------------------------------------------------------------------
    # 3. Hermite cubic coefficients (vectorised)
    # ------------------------------------------------------------------
    a0 = R_j
    a1 = G_j
    a2 = (3.0 * (R_jp1 - R_j) - h * (2.0 * G_j + G_jp1)) / h**2
    a3 = (2.0 * (R_j - R_jp1) + h * (G_j + G_jp1))        / h**3

    # ------------------------------------------------------------------
    # 4. Evaluate + USCIP clamp
    # ------------------------------------------------------------------
    R_star = a0 + xi * (a1 + xi * (a2 + xi * a3))
    R_star = np.clip(R_star, np.minimum(R_j, R_jp1), np.maximum(R_j, R_jp1))

    return R_star