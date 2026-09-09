import numpy as np
import os, sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from kernels import get_kernel
from rbf_ga_weights_1d import RBFGAStencil

def gp_interpolate(x_star, u_star, x, r,
                   nugget=0, sigma_scale=3.,
                   derivative=False,
                   kernel="se",
                   rbfga_eps=None,
    ):
    """
    Interpolate u_star (defined on x_star) onto equally-spaced grid x
    using local GP with squared-exponential kernel and spatially variable
    length scales: sigma_k = sigma_scale * nearest-neighbor-distance_k.

    Following Fornberg & Zuev (2007): epsilon_k ∝ 1/d_k suppresses the Runge
    phenomenon that arises when the source grid is non-uniform (e.g. near shocks).

    Parameters
    ----------
    x_star      : (N*,) source grid points (may be non-uniform)
    u_star      : (N*,) values at source grid
    x           : (N,)  target equally-spaced grid
    r           : int   - stencil half-width in nodes; each target point uses
                          the (2*r + 1) nearest source nodes (constant stencil size)
    nugget      : float - regularization added to kernel diagonal (default 0)
    sigma_scale : float - sigma_k = sigma_scale * d_k  (single free constant;
                          replaces the old fixed sigma = 80*dx)
    derivative : int or bool
                  0 / False : interpolation
                  1 / True  : first derivative
                  2         : second derivative (RBF-GA path only)
    rbfga_eps  : float or None
        If not None, bypass the GP kernel and use RBF-GA weights (Fornberg et al.
        2013) with this Gaussian shape parameter ε.  Stable for all ε, including
        ε → 0.  When None (default), the regular GP kernel path is used.

    Returns
    -------
    y : (N,) interpolated values at x
    """
    # Sort source grid — caller may pass unsorted x_star (e.g. reordered by a different sort)
    sort_s  = np.argsort(x_star)
    x_star  = x_star[sort_s]
    u_star  = u_star[sort_s]

    N = len(x)
    dx = (x[-1] - x[0]) / (N - 1)

    # ── Sanitise source data ──────────────────────────────────────────────────
    # 1. Drop NaN/Inf values (can arise from crashed shock in Lagrangian schemes).
    finite  = np.isfinite(u_star)
    if not finite.all():
        x_star = x_star[finite]
        u_star = u_star[finite]

    # 2. Merge near-duplicate source positions (Lagrangian particle crossings).
    #    Two nodes closer than tol are collapsed to their centroid with averaged
    #    value; otherwise the SE kernel matrix is rank-deficient → singular solve.
    if len(x_star) > 1:
        tol  = dx * 1e-8
        gaps = np.diff(x_star)
        if (gaps < tol).any():
            keep   = np.concatenate([[True], gaps >= tol])   # first of each group
            labels = np.cumsum(keep) - 1                     # group index per node
            n_grp  = int(labels[-1]) + 1
            x_new  = np.zeros(n_grp)
            u_new  = np.zeros(n_grp)
            count  = np.zeros(n_grp, dtype=int)
            np.add.at(x_new, labels, x_star)
            np.add.at(u_new, labels, u_star)
            np.add.at(count, labels, 1)
            x_star = x_new / count
            u_star = u_new / count
    r = int(r)
    stencil_size = 2 * r + 1

    # ── Spatially variable sigma (Fornberg & Zuev 2007) ───────────────────────
    # sigma_k = sigma_scale * nearest-neighbor distance at each source point.
    # Where the source grid compresses (near shocks), sigma shrinks accordingly,
    # keeping the kernel width proportional to the local node spacing.
    if len(x_star) > 1:
        d_l = np.empty_like(x_star)
        d_r = np.empty_like(x_star)
        d_l[0],  d_l[1:]  = x_star[1] - x_star[0],   x_star[1:] - x_star[:-1]
        d_r[-1], d_r[:-1] = x_star[-1] - x_star[-2],  x_star[1:] - x_star[:-1]
        d_nn = np.minimum(np.abs(d_l), np.abs(d_r))
        sigma_star = sigma_scale * np.maximum(d_nn, dx * 1e-10)  # floor avoids /0 on duplicates
    else:
        sigma_star = np.array([sigma_scale * dx])

    # --- Constant BC ghost cells -------------------------------------------------
    # Enough ghost cells to ensure every target point has stencil_size neighbours.
    n_ghost = r + 1

    x_left  = x_star[0]  - dx * np.arange(n_ghost, 0, -1)
    x_right = x_star[-1] + dx * np.arange(1, n_ghost + 1)
    u_left  = np.full(n_ghost, u_star[0])
    u_right = np.full(n_ghost, u_star[-1])

    # Ghost sigmas: copy the edge value
    sigma_left  = np.full(n_ghost, sigma_star[0])
    sigma_right = np.full(n_ghost, sigma_star[-1])

    x_ext     = np.concatenate([x_left,       x_star,      x_right      ])
    u_ext     = np.concatenate([u_left,        u_star,      u_right      ])
    sigma_ext = np.concatenate([sigma_left,    sigma_star,  sigma_right  ])
    # -----------------------------------------------------------------------------

    # (N, N_ext) absolute distances to extended source grid
    D = np.abs(x[:, None] - x_ext[None, :])

    # Pick the stencil_size nearest nodes per target point (uniform stencil size)
    M     = stencil_size
    order = np.argpartition(D, min(M - 1, D.shape[1] - 1), axis=1)[:, :M]  # (N, M)

    x_stencil     = x_ext[order]        # (N, M)
    u_stencil     = u_ext[order]        # (N, M)
    sigma_stencil = sigma_ext[order]    # (N, M)

    # ── RBF-GA path ──────────────────────────────────────────────────────────
    if rbfga_eps is not None:
        stencil = RBFGAStencil(rbfga_eps)
        deriv   = int(derivative)   # accepts 0/False, 1/True, or 2

        # All rows have the same stencil size M — single vectorised call.
        W = stencil.weights_batch(x_stencil, x, deriv)           # (N, M)
        return (W * u_stencil).sum(axis=1)

    # ── GP kernel path (default) ──────────────────────────────────────────────
    kern = get_kernel(kernel) if isinstance(kernel, str) else kernel

    # Local kernel matrices (N, M, M)
    diff = x_stencil[:, :, None] - x_stencil[:, None, :]    # (N, M, M)
    K    = kern(diff, sigma_stencil[:, None, :])             # (N, M, M) — sigma broadcast over rows
    K   += nugget * np.eye(M)

    alpha = np.linalg.solve(K, u_stencil[..., None]).squeeze(-1)   # (N, M)

    # Prediction
    C      = x[:, None] - x_stencil                          # (N, M)
    k_star = kern(C, sigma_stencil)
    if derivative:
        k_star *= -(C / sigma_stencil**2)

    return (k_star * alpha).sum(axis=1)