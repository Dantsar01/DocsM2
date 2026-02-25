import numpy as np
from scipy.spatial.distance import cdist
import ot

def compute_ot_coupling_manual(X_source, X_target, epsilon, max_iter=1000, tol=1e-6):
    """
    Implémenter l'algorithme de Sinkhorn à la main.
    
    Parameters
    ----------
    X_source : ndarray, shape (n_source, d)
        Particules sources
    X_target : ndarray, shape (n_target, d)
        Particules cibles
    epsilon : float
        Paramètre de régularisation entropique
        
    Returns
    -------
    gamma : ndarray, shape (n_source, n_target)
        Plan de transport optimal (matrice de couplage)
    """

    n_source = len(X_source)
    n_target = len(X_target)
    
    # Distributions uniformes
    a = np.ones(n_source) / n_source
    b = np.ones(n_target) / n_target
    
    # Matrice de coût et noyau de Gibbs
    C = cdist(X_source, X_target, metric='sqeuclidean')
    K = np.exp(-C / epsilon)
    
    # Initialisation
    v = np.ones(n_target)
    
    # Itérations de Sinkhorn
    for iteration in range(max_iter):
        u = a / (K @ v)
        v = b / (K.T @ u)
        
        if iteration > 0 and iteration % 10 == 0:
            # Vérification de convergence tous les 10 itérations
            marginal_error = np.max(np.abs((K @ v) * u - a))
            if marginal_error < tol:
                break
    
    # Plan de transport: gamma = diag(u) K diag(v)
    gamma = u[:, None] * K * v[None, :]
    
    return gamma


def compute_ot_coupling(X_source, X_target, epsilon):
    """
    Implémenter Sinkhorn avec ot.sinkhorn.
    
    Parameters
    ----------
    X_source : ndarray, shape (n_source, d)
        Particules sources
    X_target : ndarray, shape (n_target, d)
        Particules cibles
    epsilon : float
        Paramètre de régularisation entropique
    
    Returns
    -------
    gamma : ndarray, shape (n_source, n_target)
        Plan de transport optimal
    """
    n_source = len(X_source)
    n_target = len(X_target)
    
    # Distributions uniformes
    a = np.ones(n_source) / n_source
    b = np.ones(n_target) / n_target
    
    # Matrice de coût
    C = cdist(X_source, X_target, metric='sqeuclidean')
    
    # Sinkhorn
    gamma = ot.sinkhorn(a, b, C, reg=epsilon, numItermax=1000, stopThr=1e-9)
    
    return gamma


def build_trajectories(snapshots_dict, couplings, n_trajectories=100):
    """
    Construit des trajectoires en chaînant les couplages OT.
    
    Parameters
    ----------
    snapshots_dict : dict
        Dictionnaire {temps: array de particules}
    couplings : dict
        Dictionnaire {(t_start, t_end): matrice de couplage}
    n_trajectories : int
        Nombre de trajectoires à construire
    
    Returns
    -------
    trajectories : list of list of ndarray
        Liste de trajectoires, chaque trajectoire est une liste de positions
    """
    times = sorted(snapshots_dict.keys())
    
    # Indices de départ aléatoires
    n_particles_start = len(snapshots_dict[times[0]])
    start_indices = np.random.choice(n_particles_start, size=n_trajectories, replace=True)
    
    trajectories = []
    
    for start_idx in start_indices:
        trajectory = []
        current_idx = start_idx
        
        for k in range(len(times)):
            t_current = times[k]
            trajectory.append(snapshots_dict[t_current][current_idx])
            
            # Si ce n'est pas le dernier temps, transporter vers le suivant
            if k < len(times) - 1:
                t_next = times[k + 1]
                gamma = couplings[(t_current, t_next)]
                
                # Distribution conditionnelle: gamma[current_idx, :] (renormalisée)
                probs = gamma[current_idx, :]
                probs = probs / probs.sum()  # Normalisation
                
                # Échantillonner l'indice suivant
                n_next = len(snapshots_dict[t_next])
                current_idx = np.random.choice(n_next, p=probs)
        
        trajectories.append(trajectory)
    
    return trajectories, times


def mccann_interpolation(X_source, X_target, gamma, t, n_samples=1000):
    """
    Calcule l'interpolation de McCann au temps t entre deux distributions.
    
    Parameters
    ----------
    X_source : ndarray, shape (n_source, d)
        Particules sources (au temps 0)
    X_target : ndarray, shape (n_target, d)
        Particules cibles (au temps 1)
    gamma : ndarray, shape (n_source, n_target)
        Plan de transport optimal
    t : float
        Temps d'interpolation (entre 0 et 1)
    n_samples : int
        Nombre d'échantillons à générer
    
    Returns
    -------
    X_interp : ndarray, shape (n_samples, d)
        Particules interpolées au temps t
    """
    n_source, n_target = gamma.shape
    
    # Aplatir gamma pour avoir une distribution sur les paires
    gamma_flat = gamma.flatten()
    gamma_flat = gamma_flat / gamma_flat.sum()  # Normalisation
    
    # Échantillonner n_samples paires (i, j) selon gamma
    pair_indices = np.random.choice(n_source * n_target, size=n_samples, p=gamma_flat)
    
    # Convertir les indices plats en indices (i, j)
    source_indices = pair_indices // n_target
    target_indices = pair_indices % n_target
    
    # Interpolation: (1-t) * x + t * y
    X_interp = (1 - t) * X_source[source_indices] + t * X_target[target_indices]
    
    return X_interp


def distribution_distance(X, Y, epsilon=0):
    """
    Calcule la distance de Wasserstein (si epsilon = 0) ou la divergence entropique (si epsilon > 0) entre deux distributions empiriques.
    
    Utilise l'OT exact ou la sinkhorn divergence.
    """

    n, m = len(X), len(Y)
    a = np.ones(n) / n
    b = np.ones(m) / m
    
    C = cdist(X, Y, metric='sqeuclidean')
    if not epsilon:
        w2_squared = ot.emd2(a, b, C, numItermax=100000)
        return np.sqrt(w2_squared)
    
    def sinkhorn_divergence(X, Y):
        Cxy = cdist(X, Y, metric="sqeuclidean")
        Cxx = cdist(X, X, metric="sqeuclidean")
        Cyy = cdist(Y, Y, metric="sqeuclidean")
        return (
            ot.sinkhorn2(a, b, Cxy, reg=epsilon)
            - 0.5 * ot.sinkhorn2(a, a, Cxx, reg=epsilon)
            - 0.5 * ot.sinkhorn2(b, b, Cyy, reg=epsilon)
        )
    
    if epsilon:
        w2_squared = sinkhorn_divergence(X, Y)
        return np.sqrt(w2_squared)