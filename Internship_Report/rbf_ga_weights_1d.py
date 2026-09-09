"""
rbf_ga_weights_1d.py
Stable RBF-FD weight computation in 1D via the RBF-GA algorithm.

Reference: Fornberg, Lehto, Powell (2013)
           "Stable calculation of Gaussian-based RBF-FD stencils"
           Comput. Math. Appl. 65, 627-637.

Public API
----------
RBFGAStencil(eps)                        — main class; instantiate once per shape parameter
  .weights(x_nodes, x_c, deriv)          → (n,) FD weights  w  s.t.  w @ f ≈ L[f](x_c)
  .weights_batch(x_stencil, x_c, deriv)  → (N, M) weights for N stencils at once
  .predict(x_nodes, f, x_c, …)          → scalar prediction
  .weights_direct(…)                     → RBF-Direct (fast but unstable for small ε)
  .weights_auto(…)                       → auto-select GA vs Direct
"""

import numpy as np


# ─────────────────────────────────────────────────────────────────────────────
# Module-level helpers (private)
# ─────────────────────────────────────────────────────────────────────────────

def _build_Bk(xk, k):
    """
    Null vector of the (k × k+1) Vandermonde P[p,j] = xk[j]^p, p=0..k-1.
    The null space is 1-D (distinct nodes assumed). Returns a (k+1,) unit vector.
    """
    if k == 0:
        return np.array([1.0])
    P = xk[None, :] ** np.arange(k)[:, None]    # (k, k+1)
    _, _, Vh = np.linalg.svd(P, full_matrices=True)
    return Vh[-1]


def _build_Bpad(xn):
    """
    Precompute all B_k null vectors into a zero-padded (n, n) matrix.
    B_pad[k, :k+1] = null vector of the (k × k+1) Vandermonde of xn[:k+1].
    B_pad[k, k+1:] = 0  (padding).
    """
    n = len(xn)
    B_pad = np.zeros((n, n))
    for k in range(n):
        B_pad[k, :k + 1] = _build_Bk(xn[:k + 1], k)
    return B_pad


def _build_Gk_tensor(Z, n, extra=50):
    """
    Build the (n, n, n) tensor  GG[k, m, i] = G_k(Z[m, i]).

    Uses the Taylor tail  G_k(z) = Σ_{j≥k}  z^j / j!  via a **suffix cumsum**:
      1. Accumulate Taylor terms  all_terms[j] = Z^j / j!  for j = 0 … T.
      2. GG[k] = Σ_{j=k}^T  all_terms[j]  = suffix_cumsum[k].

    This is stable for any z (no cancellation: we never subtract exp(z)).
    The recurrence  all_terms[j] = all_terms[j-1] * Z / j  avoids overflow
    because multiplying by Z/j shrinks the term for j > |Z|.

    Parameters
    ----------
    Z     : (n, n) ndarray — Z[m, i] = 2ε² · xn[m] · xn[i]
    n     : int            — number of stages (= number of stencil nodes)
    extra : int            — extra Taylor terms beyond j = n−1 for accuracy
    """
    T = n - 1 + extra          # sum from j=0 to T  (T+1 terms)
    all_terms = np.empty((T + 1, n, n))
    all_terms[0] = 1.0
    for j in range(1, T + 1):
        all_terms[j] = all_terms[j - 1] * Z / j   # z^j / j! via recurrence

    # suffix_cumsum[k] = Σ_{j=k}^T  all_terms[j]  =  G_k(Z)
    suffix = np.cumsum(all_terms[::-1], axis=0)[::-1]
    return suffix[:n]          # (n, n, n): GG[k] = G_k(Z)


def _build_Bpad_batch(xn, M):
    """
    Batched null vectors for N stencils.

    Parameters
    ----------
    xn : (N, M) centred node positions
    M  : stencil width

    Returns
    -------
    B_pad : (N, M, M)  where B_pad[i, k, :k+1] = null vec of (k × k+1) Vandermonde of xn[i, :k+1]
    """
    N = xn.shape[0]
    B_pad = np.zeros((N, M, M))
    B_pad[:, 0, 0] = 1.0
    for k in range(1, M):
        # P[i, p, j] = xn[i, j]^p,  p = 0..k-1,  j = 0..k  →  (N, k, k+1)
        P = xn[:, :k + 1][:, None, :] ** np.arange(k)[None, :, None]
        _, _, Vh = np.linalg.svd(P, full_matrices=True)   # Vh: (N, k+1, k+1)
        B_pad[:, k, :k + 1] = Vh[:, -1, :]
    return B_pad


def _build_Gk_tensor_batch(Z, n, extra=20):
    """
    Batched version of _build_Gk_tensor for N stencils.

    Parameters
    ----------
    Z : (N, M, M) ndarray — Z[i, m, j] = 2ε² · xn[i, m] · xn[i, j]
    n : int               — number of stages (= stencil width M)

    Returns
    -------
    GG : (n, N, M, M)  — GG[k, i, m, j] = G_k(Z[i, m, j])
    """
    T = n - 1 + extra
    all_terms = np.empty((T + 1,) + Z.shape)   # (T+1, N, M, M)
    all_terms[0] = 1.0
    for j in range(1, T + 1):
        all_terms[j] = all_terms[j - 1] * Z / j
    suffix = np.cumsum(all_terms[::-1], axis=0)[::-1]
    return suffix[:n]                           # (n, N, M, M)


# ─────────────────────────────────────────────────────────────────────────────
# Main class
# ─────────────────────────────────────────────────────────────────────────────

class RBFGAStencil:
    """
    1D RBF-GA stencil computer for Gaussian RBFs  φ(r) = exp(−(ε·r)²).

    Instantiate with a fixed shape parameter ε, then call .weights() or
    .predict() for any local stencil.

    Parameters
    ----------
    eps           : float > 0
        Gaussian shape parameter ε.
    eps_threshold : float
        Above this value, weights_auto() falls back to the faster RBF-Direct.
    """

    def __init__(self, eps, eps_threshold=0.3):
        self.eps           = float(eps)
        self.eps_threshold = float(eps_threshold)

    # ── RBF-GA (stable, always) ───────────────────────────────────────────────

    def weights(self, x_nodes, x_c, deriv=0):
        """
        RBF-GA FD weights in 1D. Stable for all ε > 0, including ε → 0.

        Algorithm (1D specialisation of Fornberg et al. 2013, §5)
        ---------------------------------------------------------
        Stage k (k = 0 … n−1) produces one well-conditioned basis function:

            ψ_{k+1}(x) = exp(−ε²x²) · (1/ε^{2k}) · B_k · [G_k(2ε²x·x₁) … G_k(2ε²x·x_{k+1})]

        B_k = null vector of the (k × k+1) monomial Vandermonde of the first
        k+1 centred nodes.  The null space analytically removes the ill-conditioned
        Taylor terms of order 0 … k−1, leaving only the stable G_k remainder.

        Nodes are centred at x_c internally (FD weights are translation-invariant).

        Parameters
        ----------
        x_nodes : (n,) array_like  — stencil node positions
        x_c     : float            — evaluation point
        deriv   : 0 → interpolation L[f] = f(x_c)
                  1 → first derivative L[f] = f'(x_c)

        Returns
        -------
        w : (n,) ndarray  s.t.  w @ f_nodes ≈ L[f](x_c)
        """
        eps     = self.eps
        x_nodes = np.asarray(x_nodes, dtype=float)
        n       = len(x_nodes)
        xn      = x_nodes - x_c          # centred nodes; x_c maps to 0

        B_pad      = _build_Bpad(xn)                               # (n, n)
        Z          = 2.0 * eps**2 * np.outer(xn, xn)              # (n, n)
        GG         = _build_Gk_tensor(Z, n)                       # (n, n, n)
        scales     = eps**(-2.0 * np.arange(n, dtype=float))      # (n,)
        scales[0]  = 1.0
        gauss      = np.exp(-eps**2 * xn**2)                      # (n,)

        # contracted[k, m] = Σ_i GG[k, m, i] * B_pad[k, i]
        contracted = np.einsum('kmi,ki->km', GG, B_pad)           # (n, n)
        Psi        = scales[:, None] * gauss[None, :] * contracted # (n, n)

        # b vector — sparse in centred coords (x_c maps to 0)
        # G_k(0) = 0^k/k! tail = 1 for k=0, 0 for k≥1.
        # deriv=0: b[k] = scale * B_pad[k,:k+1] @ G_k(0)^{k+1} → b[0]=1, rest 0
        # deriv=1: b[k] = scale * B_pad[k,:k+1] @ (2ε²xn[:k+1] * G_{k-1}(0))
        #          G_{k-1}(0)=1 only for k-1=0 → k=1; also G_{max(0,-1)}(0)=1 → k=0
        b = np.zeros(n)
        if deriv == 0:
            b[0] = 1.0
        elif deriv == 1:
            # At z=0: G_k(0)=1 iff k=0, else 0.  G_{k-1}(0)=1 iff k<=1.
            Gkm1_0    = np.zeros(n)
            Gkm1_0[0] = 1.0          # k=0: G_{-1} → G_0(0) = 1
            if n > 1:
                Gkm1_0[1] = 1.0      # k=1: G_0(0) = 1
            b = scales * (B_pad @ (2.0 * eps**2 * xn)) * Gkm1_0
        elif deriv == 2:
            # d²/dx² ψ_k|_{x=0} = s_k * [g''(0)*h_k(0) + h_k''(0)]
            # g''(0) = -2ε²,  h_k(0) = 1 only for k=0
            # h_k''(0) = G_{k-2}(0) * (2ε²)² * (B_pad @ xn²)[k]
            # G_{k-2}(0) = 1 for k=0,1,2 (same negative-index convention as deriv=1)
            Gkm2_0    = np.zeros(n)
            for i in range(min(3, n)):
                Gkm2_0[i] = 1.0
            b  = scales * (4.0 * eps**4 * (B_pad @ xn**2) * Gkm2_0)
            b[0] -= 2.0 * eps**2     # g''(0)*h_0(0) contribution (scales[0]=1)

        return np.linalg.solve(Psi, b)

    # ── RBF-GA batch (N stencils at once) ────────────────────────────────────

    def weights_batch(self, x_stencil, x_centers, deriv=0, nugget=0.0):
        """
        Compute RBF-GA FD weights for N stencils simultaneously.

        All stencils must have the same padded width M.  Invalid (padding) slots
        should be pre-replaced by safe dummy positions near each x_center before
        calling this method (see gp_interp.py); u values for those slots must be
        zero so they do not contribute to the prediction.

        Parameters
        ----------
        x_stencil : (N, M) — stencil node positions (pre-sanitised, no sentinels)
        x_centers : (N,)   — evaluation points
        deriv     : 0 or 1
        nugget    : float  — small diagonal regularisation for near-singular Ψ

        Returns
        -------
        W : (N, M)  weights s.t.  (W * u_stencil).sum(axis=1) ≈ L[u](x_centers)
        """
        eps       = self.eps
        x_stencil = np.asarray(x_stencil, dtype=float)
        x_centers = np.asarray(x_centers, dtype=float)
        N, M      = x_stencil.shape

        xn = x_stencil - x_centers[:, None]                      # (N, M)

        B_pad = _build_Bpad_batch(xn, M)                         # (N, M, M)
        Z     = 2.0 * eps**2 * xn[:, :, None] * xn[:, None, :]  # (N, M, M)
        GG    = _build_Gk_tensor_batch(Z, M)                     # (M, N, M, M)

        scales    = eps**(-2.0 * np.arange(M, dtype=float))      # (M,)
        scales[0] = 1.0
        gauss     = np.exp(-eps**2 * xn**2)                      # (N, M)

        # contracted[k, i, m] = Σ_j GG[k, i, m, j] * B_pad[i, k, j]
        contracted = np.einsum('knmj,nkj->knm', GG, B_pad)       # (M, N, M)

        # Psi[i, k, m] = scales[k] * gauss[i, m] * contracted[k, i, m]
        Psi  = (scales[:, None, None] * gauss[None, :, :] * contracted
                ).transpose(1, 0, 2)                              # (N, M, M)
        Psi += nugget * np.eye(M)

        b = np.zeros((N, M))
        if deriv == 0:
            b[:, 0] = 1.0
        elif deriv == 1:
            Gkm1_0    = np.zeros(M)
            Gkm1_0[0] = 1.0
            if M > 1:
                Gkm1_0[1] = 1.0
            b = (scales[None, :]
                 * np.einsum('nkj,nj->nk', B_pad, 2.0 * eps**2 * xn)
                 * Gkm1_0[None, :])
        elif deriv == 2:
            Gkm2_0    = np.zeros(M)
            for i in range(min(3, M)):
                Gkm2_0[i] = 1.0
            b  = (scales[None, :]
                  * 4.0 * eps**4 * np.einsum('nkj,nj->nk', B_pad, xn**2)
                  * Gkm2_0[None, :])
            b[:, 0] -= 2.0 * eps**2  # g''(0)*h_0(0) contribution

        return np.linalg.solve(Psi, b[..., None]).squeeze(-1)    # (N, M)

    # ── RBF-Direct (fast but unstable for small ε) ────────────────────────────

    def weights_direct(self, x_nodes, x_c, deriv=0):
        """RBF-Direct: solve A w = b directly. Fast but fails for small ε."""
        eps     = self.eps
        x_nodes = np.asarray(x_nodes, dtype=float)
        r2      = (x_nodes[:, None] - x_nodes[None, :])**2
        A       = np.exp(-eps**2 * r2)
        rc2     = (x_c - x_nodes)**2
        if deriv == 0:
            b_rhs = np.exp(-eps**2 * rc2)
        else:
            b_rhs = -2.0 * eps**2 * (x_c - x_nodes) * np.exp(-eps**2 * rc2)
        return np.linalg.solve(A, b_rhs)

    def weights_direct_batch(self, x_stencil, x_centers, deriv=0):
        """RBF-Direct batch for N stencils. Stable for large ε."""
        eps       = self.eps
        N, M      = x_stencil.shape
        diff2     = (x_stencil[:, :, None] - x_stencil[:, None, :]) ** 2  # (N, M, M)
        A         = np.exp(-eps**2 * diff2)                                # (N, M, M)
        rc2       = (x_centers[:, None] - x_stencil) ** 2                 # (N, M)
        if deriv == 0:
            b = np.exp(-eps**2 * rc2)
        else:
            b = -2.0 * eps**2 * (x_centers[:, None] - x_stencil) * np.exp(-eps**2 * rc2)
        return np.linalg.solve(A, b[..., None]).squeeze(-1)               # (N, M)

    # ── Auto: GA for small ε, Direct for large ε ─────────────────────────────

    def weights_auto(self, x_nodes, x_c, deriv=0):
        """Use RBF-GA if ε ≤ eps_threshold, else RBF-Direct."""
        if self.eps <= self.eps_threshold:
            return self.weights(x_nodes, x_c, deriv)
        return self.weights_direct(x_nodes, x_c, deriv)

    def weights_auto_batch(self, x_stencil, x_centers, deriv=0):
        """Batch auto-dispatch: RBF-GA for small ε, RBF-Direct for large ε."""
        if self.eps <= self.eps_threshold:
            return self.weights_batch(x_stencil, x_centers, deriv)
        return self.weights_direct_batch(x_stencil, x_centers, deriv)

    # ── Convenience: predict directly ────────────────────────────────────────

    def predict(self, x_nodes, f_nodes, x_c, deriv=0):
        """Return  weights(x_nodes, x_c, deriv) @ f_nodes."""
        return self.weights(x_nodes, x_c, deriv) @ np.asarray(f_nodes, dtype=float)


# ─────────────────────────────────────────────────────────────────────────────
# Uniform-grid stencil (weights computed once, reused everywhere)
# ─────────────────────────────────────────────────────────────────────────────

class RBFGAUniformGrid:
    """
    Precomputed RBF-GA FD weights for a uniform 1D grid.

    On a uniform grid with spacing dx every interior stencil of width M is a
    translate of one prototype, so the weight vector is identical for every
    interior node and needs to be computed only once.

    Parameters
    ----------
    dx       : float — uniform grid spacing
    M        : int   — stencil width (odd recommended, e.g. 2*r + 1)
    eps      : float — Gaussian shape parameter ε
    deriv    : 0, 1, or 2 — operator to precompute (default 0)
    x_offset : float — where to evaluate inside the stencil, relative to the
                       centre node (default 0.0 = on-node).

    Attributes
    ----------
    w : (M,) ndarray — precomputed weights for the chosen deriv

    Quick-start
    -----------
    >>> g2 = RBFGAUniformGrid(dx=0.01, M=7, eps=1.0, deriv=2)
    >>> d2f = g2.apply(f)                    # second derivative on uniform grid
    >>> A2  = g2.build_matrix(N)             # full N×N matrix (build once)
    >>> d2f = A2 @ f                         # fast repeated application
    """

    def __init__(self, dx, M, eps, deriv=0, x_offset=0.0):
        self.dx       = float(dx)
        self.M        = int(M)
        self.eps      = float(eps)
        self.deriv    = int(deriv)
        self.x_offset = float(x_offset)

        r       = (M - 1) // 2
        x_proto = np.arange(-r, M - r, dtype=float) * dx   # (M,) prototype nodes
        self.w  = RBFGAStencil(eps).weights(x_proto, x_offset, deriv=self.deriv)

    # ── Point-wise application ────────────────────────────────────────────────

    def apply(self, f, bc='clamp'):
        """
        Apply the precomputed weights to every point of a 1D array.

        Parameters
        ----------
        f  : (N,) array_like
        bc : 'clamp' or 'periodic'

        Returns
        -------
        out : (N,) ndarray  s.t.  out[i] ≈ L[f](x_i + x_offset)
        """
        f    = np.asarray(f, dtype=float)
        N    = len(f)
        r    = (self.M - 1) // 2
        out  = np.empty(N)
        base = np.arange(-r, self.M - r)   # (M,) relative offsets

        for i in range(N):
            if bc == 'periodic':
                idx = (i + base) % N
            else:
                idx = np.clip(i + base, 0, N - 1)
            out[i] = self.w @ f[idx]
        return out

    # ── Full matrix build (amortise index work for repeated application) ───────

    def build_matrix(self, N, bc='clamp'):
        """
        Build a dense (N, N) application matrix for a grid of N points.

        Parameters
        ----------
        N  : int
        bc : 'clamp' or 'periodic'

        Returns
        -------
        A : (N, N) ndarray
        """
        r    = (self.M - 1) // 2
        base = np.arange(-r, r + 1)        # (M,)
        A    = np.zeros((N, N))

        for i in range(N):
            if bc == 'periodic':
                idx = (i + base) % N
            else:
                idx = np.clip(i + base, 0, N - 1)
            np.add.at(A[i], idx, self.w)

        return A


# ─────────────────────────────────────────────────────────────────────────────
# Validation / demo
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import matplotlib.pyplot as plt

    rng = np.random.default_rng(42)
    f   = lambda x: np.sin(x) + 0.5 * np.cos(2 * x)
    df  = lambda x: np.cos(x) - np.sin(2 * x)

    x_c     = 0.3
    n       = 14
    xmin = -1.#x_c - 2 / spacing
    xmax = 1.#x_c + 2 / spacing
    x_nodes = x_c + np.linspace(xmin, xmax, n) + 0.05 * rng.standard_normal(n)
    f_nodes = f(x_nodes)

    epsilons = np.logspace(-4, 1, 100)
    results  = {lbl: {"f0": [], "d1": []} for lbl in ["RBF-GA", "RBF-Direct"]}

    for eps in epsilons:
        stencil = RBFGAStencil(eps)
        for lbl, method in [("RBF-GA", stencil.weights), ("RBF-Direct", stencil.weights_direct)]:
            for key, deriv, exact in [("f0", 0, f(x_c)), ("d1", 1, df(x_c))]:
                try:
                    w = method(x_nodes, x_c, deriv)
                    results[lbl][key].append(abs(w @ f_nodes - exact))
                except Exception:
                    results[lbl][key].append(np.nan)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, key, title in zip(axes, ["f0", "d1"],
                               ["Interpolation  f(x_c)", "Derivative  f'(x_c)"]):
        for lbl in ["RBF-GA", "RBF-Direct"]:
            ax.semilogy(epsilons, results[lbl][key], label=lbl)
        ax.set_xscale("log"); ax.set_xlabel("ε"); ax.set_ylabel("absolute error")
        ax.set_title(title); ax.legend(); ax.grid(True, which="both", alpha=0.4)
    fig.suptitle(f"n = {n} nodes,  x_c = {x_c}")
    plt.tight_layout(); plt.show()

    print(f"\nn={n}, x_c={x_c}")
    print(f"{'eps':>8}  {'GA f0':>10}  {'GA d1':>10}  {'Dir f0':>10}  {'Dir d1':>10}")
    for eps in [1e-4, 1e-3, 1e-2, 0.1, 1.0]:
        s = RBFGAStencil(eps)
        w0 = s.weights(x_nodes, x_c, 0);  w1 = s.weights(x_nodes, x_c, 1)
        try:
            wd0 = s.weights_direct(x_nodes, x_c, 0)
            wd1 = s.weights_direct(x_nodes, x_c, 1)
            ed0, ed1 = abs(wd0 @ f_nodes - f(x_c)), abs(wd1 @ f_nodes - df(x_c))
        except Exception:
            ed0 = ed1 = float("nan")
        print(f"{eps:>8.0e}  {abs(w0@f_nodes-f(x_c)):>10.2e}  {abs(w1@f_nodes-df(x_c)):>10.2e}"
              f"  {ed0:>10.2e}  {ed1:>10.2e}")
        
# %%

# Does the RBF-GA stencil reproduce degree-k polynomials exactly?
from rbf_ga_weights_1d import RBFGAStencil
import numpy as np

dx   = 0.01
r    = 2
M    = 2*r + 1
x_src = np.linspace(-2, 2, 40) + 0.4 * dx  # shifted source grid
x_tgt = np.linspace(-2, 2, 40)             # uniform target

stencil = RBFGAStencil(1e-3)
# build stencil as gp_interp does: M nearest per target
D = np.abs(x_tgt[:, None] - x_src[None, :])
order = np.argpartition(D, M-1, axis=1)[:, :M]
xs = x_src[order]   # (N, M)

for deg in range(M):
    f_src = xs**deg
    W = stencil.weights_batch(xs, x_tgt, deriv=0)
    result = (W * f_src).sum(axis=1)
    error  = np.max(np.abs(result - x_tgt**deg))
   
    print(f"degree {deg}: max reproduction error = {error:.2e}")
# %%

# ──────────────────────────────────────────────────────────────────────────────
# Spatial convergence: freeze deg=3, vary h and r
# For a polynomial of degree d, a stencil of size M=2r+1 reproduces it exactly
# when d <= 2r (error ~ machine eps). When d > 2r the error decays as O(h^(2r+1)).
# Here we freeze d=3 and sweep r=1,2,3 to see the transition.
# ──────────────────────────────────────────────────────────────────────────────


deg    = 3
f_test = lambda x: x**deg
x_c    = 0.3                         # fixed target
delta  = 0.8                         # fractional offset so x_c is never a node
h_list = [0.20, 0.10, 0.05, 0.025, 0.0125, 0.00625]

print(f"\n── Spatial convergence  (f = x^{deg},  x_c = {x_c},  node offset δ={delta}) ──")
print(f"{'r':>3}  {'M':>3}  {'repro deg≤':>10}  {'h':>8}  {'error':>10}  {'order':>6}")
print("-" * 55)

for r in [1, 2, 3]:
    M       = 2 * r + 1
    stencil = RBFGAStencil(1e-2)
    errors  = []

    for h in h_list:
        # nodes at x_c + (i + delta)*h  →  x_c is never a stencil node
        x_nodes = x_c + (np.arange(-r, r + 1) + delta) * h
        w       = stencil.weights(x_nodes, x_c, deriv=0)
        errors.append(abs(w @ f_test(x_nodes) - f_test(x_c)))

    repro = 2 * r
    exact = deg <= repro
    for i, h in enumerate(h_list):
        if i == 0:
            order_str = "   -"
        elif errors[i] < 1e-12:
            order_str = "  ~0"
        else:
            order_str = f"{np.log(errors[i-1]/errors[i]) / np.log(h_list[i-1]/h_list[i]):6.2f}"
        suffix = f"  ← exact (deg {deg} ≤ {repro})" if (exact and i == 0) else ""
        print(f"{r:>3}  {M:>3}  {repro:>10}  {h:>8.5f}  {errors[i]:>10.2e}  {order_str}{suffix}")
    print()

# %%


