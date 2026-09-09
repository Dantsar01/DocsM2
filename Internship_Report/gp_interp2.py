import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve


def gp_interpolate(x_star, u_star, x, m=6, nugget=1e-10, sigma_scale=10., derivative=False):
    """
    Remap u_star (defined on non-uniform x_star) onto the uniform grid x
    by solving the GP system  A @ U = u_star, with periodic boundary conditions.

    The GP is oriented FROM the uniform grid (source) TO the non-uniform points
    (target).  Because the source is uniform the local M×M kernel matrix is
    identical for every stencil — inverted once — reducing cost from O(N·M³)
    to O(N·M²) for assembly, plus a sparse banded solve for the N×N system.

    With periodic BC the derivative operator D (uniform → uniform) is circulant,
    so derivative=True costs O(N log N) via FFT.

    Parameters
    ----------
    x_star     : (N,)  non-uniform points  (e.g. Lagrangian particle positions)
    u_star     : (N,)  values at x_star
    x          : (N,)  uniform target grid (equally spaced, periodic domain)
    m          : int   stencil size (uniform grid points per local kernel)
    nugget     : float diagonal regularisation added to K_uniform
    sigma_scale: float sigma = sigma_scale * dx  (constant)
    derivative : bool  if True return dU/dx on the uniform grid

    Returns
    -------
    U : (N,) interpolated values (or their x-derivative) on the uniform grid

    Notes
    -----
    Assumes len(x_star) == len(x).  Once shocks cause particle clustering or
    rarefaction this assumption breaks; switch to a least-squares solver then.
    """
    N = len(x)
    assert len(x_star) == N, "x_star and x must have equal length (no-shock regime)"

    dx = (x[-1] - x[0]) / (N - 1)
    sigma = sigma_scale * dx

    # ── K_uniform: M×M kernel for m consecutive uniform points, inverted once ──
    p = np.arange(m)
    diff_local = (p[:, None] - p[None, :]) * dx          # (m, m)
    K_unif = np.exp(-diff_local**2 / (2.0 * sigma**2))
    K_unif += nugget * np.eye(m)
    K_inv = np.linalg.inv(K_unif)                         # single inversion

    # ── Stencil indices with periodic wrap ──────────────────────────────────────
    center   = np.round((x_star - x[0]) / dx).astype(int)  # (N,) nearest grid index
    raw_cols = center[:, None] - m // 2 + p[None, :]        # (N, m) unwrapped offsets
    col_idx  = raw_cols % N                                  # (N, m) periodic indices

    # ── Cross-covariance using unwrapped distances (kernel sees true geometry) ──
    x_stencil = x[0] + raw_cols * dx                        # (N, m) unwrapped positions
    delta     = x_star[:, None] - x_stencil                 # (N, m)
    k_star    = np.exp(-delta**2 / (2.0 * sigma**2))

    # ── Interpolation weights: one matmul, no per-row solve ─────────────────────
    weights = k_star @ K_inv                                 # (N, m)

    # ── Sparse forward operator A (N × N) ───────────────────────────────────────
    row_idx = np.repeat(np.arange(N), m)
    A = csr_matrix((weights.ravel(), (row_idx, col_idx.ravel())), shape=(N, N))

    # ── Recover U on uniform grid ────────────────────────────────────────────────
    U = spsolve(A, u_star)

    if not derivative:
        return U

    # ── Derivative: D is circulant (uniform → uniform) → FFT in O(N log N) ──────
    # D[i,j] depends only on (i-j) % N  because both source and target are uniform.
    # Weights for a stencil centred at any uniform grid point (same for all rows):
    delta_d0 = (m // 2 - p) * dx                            # (m,) distances to stencil
    k_d0     = np.exp(-delta_d0**2 / (2.0 * sigma**2)) * (-delta_d0 / sigma**2)
    w_d0     = k_d0 @ K_inv                                  # (m,) derivative weights

    # Build first column of D then apply as circulant matvec via FFT
    d_col              = np.zeros(N)
    d_col[(m // 2 - p) % N] = w_d0

    return np.fft.ifft(np.fft.fft(d_col) * np.fft.fft(U)).real
