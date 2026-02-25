import numpy as np



def grad_potential(x, t, well_separation=2.5, speed=7):
    """
    Gradient of a time-dependent potential with complex branching dynamics.
    Version corrigée avec transitions continues.
    
    Dynamics:
    - t ∈ [0, 0.5]: Bifurcation from one well into two wells (left and right)
    - t ∈ [0.5, 1.0]: 
        * Left well: moves along x_0 axis (drifts to the left)
        * Right well: splits again into two sub-wells (up and down)
    
    Parameters
    ----------
    x : array-like, shape (d,)
        Particle position (d >= 2).
    t : float
        Time (should be in [0, 1]).
    well_separation : float
        Distance between wells.
    speed : float
        Speed of drift along x_1 axis.
    
    Returns
    -------
    grad_V : ndarray, shape (d,)
        Gradient of the potential at (x, t).
    """
    x = np.asarray(x)
    dim = x.size
    grad_V = np.zeros(dim)
    
    # Paramètres de transition
    t_transition = 0.5
    scale_depth = 4.0
    
    # Phase 1: Initial bifurcation (t < 0.5)
    if t < t_transition:
        # Profondeur croissante pour la bifurcation initiale
        # On veut que depth augmente de 0 à une valeur élevée
        depth = (t / t_transition)  # De 0 à 1
    
        
        # Deux puits fixes
        x0_left = -well_separation
        x0_right = well_separation
        
        grad_V[0] = scale_depth * depth/2 * ((x[0] - x0_left) * (x[0] - x0_right)**2 + 
                                            (x[0] - x0_left)**2 * (x[0] - x0_right))
    
    # Phase 2: Complex branching (t >= 0.5)
    else:
        # Paramètre de transition (0 à 0.5 quand t va de 0.5 à 1.0)
        tau = (t - t_transition) / (1.0 - t_transition)  # tau ∈ [0, 1]
        
        # Potentiel à 2 puits (transition douce)
        x0_left = -well_separation
        x0_right = well_separation
        grad_2wells_x0 = x[0] * (x[0] - x0_left) * (x[0] - x0_right)

        x_left_up = 2 * x0_left
        x_left_down = 0 * x0_left
        x_right_up = 2 * x0_right
        x_right_down = 0 * x0_right
        
        # Au début de la phase 2 (t=0.5), on a encore 2 puits
        # À la fin (t=1.0), on a 3 puits
        # Poids de transition: au début on garde le potentiel à 2 puits,
        # à la fin on a un potentiel à 3 puits
        weight_2wells = (1 - tau)**2  # Décroît de 1 à 0
        weight_3wells = (tau)**(1/4)        # Croît de 0 à 1
        
        # Potentiel à 3 puits
        # Distances aux trois puits (dans l'espace 2D)
        d_left_up = np.sqrt((x[0] - x_left_up)**2)
        d_left_down = np.sqrt((x[0] - x_left_down)**2)
        d_right_up = np.sqrt((x[0] - x_right_up)**2)
        d_right_down = np.sqrt((x[0] - x_right_down)**2)
        
        # Poids softmax pour chaque puits
        temperature = .5
        w_left_down = np.exp(-d_left_down**2 / temperature)
        w_left_up = np.exp(-d_left_up**2 / temperature)
        w_right_up = np.exp(-d_right_up**2 / temperature)
        w_right_down = np.exp(-d_right_down**2 / temperature)
        total_w = w_left_up + w_right_up + w_right_down + w_left_down
        
        w_left_up /= total_w
        w_left_down /= total_w
        w_right_up /= total_w
        w_right_down /= total_w
        
        # Forces attractives (gradient négatif = attraction)
        grad_3wells_x0 = (
            w_left_down * (x[0] - x_left_down) +
            w_left_up * (x[0] - x_left_up) +
            w_right_up * (x[0] - x_right_up) +
            w_right_down * (x[0] - x_right_down))
        
        # Combinaison pondérée
        grad_V[0] = scale_depth * (weight_2wells/2 * grad_2wells_x0 + weight_3wells * grad_3wells_x0)
        
    if dim >= 2:
        # Pour x_1, on combine la dérive globale avec l'attraction des puits
        grad_V[1] = speed * (1 + t)
    
    # Higher dimensions: add damping
    for i in range(2, dim):
        grad_V[i] = speed * x[i]
    
    return grad_V


def simulate_sde(
    n_particles=100*np.ones(6),
    dim=2,
    t0=0.0,
    t1=1.0,
    dt=1e-3,
    sigma=0.5,
    snapshot_times=None,
    seed=0,
):
    """
    Simulate a stochastic differential equation using Euler-Maruyama.

    dX_t = -∇V(X_t, t) dt + sqrt(2 * sigma) dB_t

    Parameters
    ----------
    n_particles : int
        Number of particles.
    dim : int
        Dimension of the state space.
    t0, t1 : float
        Initial and final times.
    dt : float
        Time step.
    sigma : float
        Diffusion coefficient.
    snapshot_times : list or array of floats
        Times at which distributions are recorded.
    seed : int
        Random seed.

    Returns
    -------
    snapshots : dict
        Dictionary {t: array of shape (n_particles, dim)}
    """
    rng = np.random.default_rng(seed)

    if snapshot_times is None:
        snapshot_times = np.linspace(t0, t1, len(n_particles))

    snapshot_times = np.array(snapshot_times)
    snapshots = {}

    # Initial condition: centered Gaussian

    for snapshot_idx, time in enumerate(snapshot_times):
        t = t0
        X = rng.normal(0.0, 0.5, size=(n_particles[snapshot_idx], dim))
        while t <= time:
            if np.isclose(t, snapshot_times[snapshot_idx], atol=dt):
                snapshots[snapshot_times[snapshot_idx]] = X.copy()
                snapshot_idx += 1

            # Euler-Maruyama step
            drift = np.array([grad_potential(x, t) for x in X])

            noise = rng.normal(0.0, 1.0, size=X.shape)

            X += -drift * dt + np.sqrt(2 * sigma * dt) * noise   

            t += dt   

    return snapshots
