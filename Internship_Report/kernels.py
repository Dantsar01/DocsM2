import numpy as np


def get_kernel(name, **kwargs):
    """
    Return a kernel function k(diff, sigma) for the requested kernel.

    The returned callable has signature
        k(diff, sigma) -> ndarray, same shape as diff
    where diff = x - x' is the signed difference (scalar or array).

    Parameters
    ----------
    name : str
        'se'         Squared Exponential (Gaussian)
        'mq'         Multiquadric                      (cond. neg. definite)
        'imq'        Inverse Multiquadric
        'iq'         Inverse Quadratic
        'rq'         Rational Quadratic                (needs alpha=, default 1)
        'matern12'   Matérn ν=1/2  (exponential)
        'matern32'   Matérn ν=3/2
        'matern52'   Matérn ν=5/2
        'wendland2'  Wendland C²   (compact support on |diff| < sigma)
        'wendland4'  Wendland C⁴   (compact support on |diff| < sigma)
        'bump'       C∞ bump       (compact support on |diff| < sigma)
        'cubic'      Cubic |r|³    (alias for phs3)
        'tps'        Thin-plate spline r² log r  (alias for phs2)
        'periodic'   Periodic SE                       (needs period=, default 1)
        'phs'        Polyharmonic spline, general      (needs k=, default 3)
                       odd  k : φ(r) = r^k
                       even k : φ(r) = r^k log(r),  φ(0) = 0
        'phs1'       |r|           (linear, k=1)
        'phs3'       |r|³          (cubic,  k=3)
        'phs5'       |r|⁵          (k=5)
        'phs7'       |r|⁷          (k=7)
        'phs2'       r² log r      (thin-plate, k=2)
        'phs4'       r⁴ log r      (k=4)
        'phs6'       r⁶ log r      (k=6)

    **kwargs
        alpha  : float – shape parameter for 'rq'       (default 1.0)
        period : float – domain period for 'periodic'   (default 1.0)
        k      : int   – degree for 'phs'               (default 3)
    """
    name = name.lower().strip()

    if name == 'se':
        def k(diff, sigma):
            return np.exp(-diff**2 / (2.0 * sigma**2))

    elif name == 'mq':
        def k(diff, sigma):
            return np.sqrt(1.0 + diff**2 / sigma**2)

    elif name == 'imq':
        def k(diff, sigma):
            return 1.0 / np.sqrt(1.0 + diff**2 / sigma**2)

    elif name == 'iq':
        def k(diff, sigma):
            return 1.0 / (1.0 + diff**2 / sigma**2)

    elif name == 'rq':
        alpha = float(kwargs.get('alpha', 1.0))
        def k(diff, sigma):
            return (1.0 + diff**2 / (2.0 * alpha * sigma**2)) ** (-alpha)

    elif name == 'matern12':
        def k(diff, sigma):
            return np.exp(-np.abs(diff) / sigma)

    elif name == 'matern32':
        def k(diff, sigma):
            s = np.sqrt(3.0) * np.abs(diff) / sigma
            return (1.0 + s) * np.exp(-s)

    elif name == 'matern52':
        def k(diff, sigma):
            s = np.sqrt(5.0) * np.abs(diff) / sigma
            return (1.0 + s + s**2 / 3.0) * np.exp(-s)

    elif name == 'wendland2':
        # (1 - r)_+^4 (4r + 1),  C² on R,  compact support |diff| < sigma
        def k(diff, sigma):
            r = np.abs(diff) / sigma
            t = np.maximum(1.0 - r, 0.0)
            return t**4 * (4.0 * r + 1.0)

    elif name == 'wendland4':
        # (1 - r)_+^6 (35r² + 18r + 3) / 3,  C⁴ on R,  compact support
        def k(diff, sigma):
            r = np.abs(diff) / sigma
            t = np.maximum(1.0 - r, 0.0)
            return t**6 * (35.0 * r**2 + 18.0 * r + 3.0) / 3.0

    elif name == 'bump':
        # exp(-1 / (1 - r²)) for |diff| < sigma, else 0  — C∞, compact support
        def k(diff, sigma):
            r = np.abs(diff) / sigma
            inside = r < 1.0
            out = np.zeros(np.shape(diff), dtype=float)
            out[inside] = np.exp(-1.0 / (1.0 - r[inside]**2))
            return out

    elif name == 'cubic' or name == 'phs3':
        def k(diff, sigma):
            return (np.abs(diff) / sigma) ** 3

    elif name == 'tps' or name == 'phs2':
        def k(diff, sigma):
            r = np.abs(diff) / sigma
            return np.where(r == 0.0, 0.0, r**2 * np.log(r))

    elif name == 'phs1':
        def k(diff, sigma):
            return np.abs(diff) / sigma

    elif name == 'phs5':
        def k(diff, sigma):
            return (np.abs(diff) / sigma) ** 5

    elif name == 'phs7':
        def k(diff, sigma):
            return (np.abs(diff) / sigma) ** 7

    elif name == 'phs4':
        def k(diff, sigma):
            r = np.abs(diff) / sigma
            return np.where(r == 0.0, 0.0, r**4 * np.log(r))

    elif name == 'phs6':
        def k(diff, sigma):
            r = np.abs(diff) / sigma
            return np.where(r == 0.0, 0.0, r**6 * np.log(r))

    elif name == 'phs':
        deg = int(kwargs.get('k', 3))
        if deg % 2 == 1:
            def k(diff, sigma):
                return (np.abs(diff) / sigma) ** deg
        else:
            def k(diff, sigma):
                r = np.abs(diff) / sigma
                return np.where(r == 0.0, 0.0, r**deg * np.log(r))

    elif name == 'periodic':
        # exp(-2 sin²(π diff / period) / sigma²) — SE on the circle
        period = float(kwargs.get('period', 1.0))
        def k(diff, sigma):
            return np.exp(-2.0 * np.sin(np.pi * diff / period)**2 / sigma**2)

    else:
        supported = ('se, mq, imq, iq, rq, matern12, matern32, matern52, '
                     'wendland2, wendland4, bump, periodic, '
                     'phs1, phs2, phs3, phs4, phs5, phs6, phs7, phs(k=), '
                     'cubic (=phs3), tps (=phs2)')
        raise ValueError(f"Unknown kernel '{name}'. Supported: {supported}")

    return k
