# run_gwot_on_simulation.py
"""
Example script:
- simulate a branching SDE with sparse sampling (100, 10, ..., 10, 100),
- simulate in parallel a "full" SDE with 100 cells at all timepoints,
- run global Waddington-OT (gWOT) on the sparse data to infer 100 particles per time,
- compare inferred marginals with the full simulation (100 everywhere),
- plot:
    - left: gWOT loss vs iterations,
    - right: distances per time:
        * inferred vs full,
        * sparse vs full.
"""

import numpy as np
import matplotlib.pyplot as plt
import torch

import sys
from pathlib import Path

base_dir = Path(__file__).resolve().parent
sys.path.append(str(base_dir / "src"))

from simulation import simulate_sde
from functions import distribution_distance
from models_gwot import TrajLoss, optimize_model

from tqdm import tqdm


# %%



plt.rcParams["figure.figsize"] = (10, 4)
plt.rcParams["font.size"] = 12


def resample_points(X, n_samples):
    """
    Resample n_samples points from X (with replacement).

    Parameters
    ----------
    X : np.ndarray, shape (n, d)
    n_samples : int
    seed : int or None

    Returns
    -------
    Y : np.ndarray, shape (n_samples, d)
    """
    n, d = X.shape
    n_runs = int(n_samples / n) + 1
    idx = np.ones(n_samples, dtype=int)
    cnt = 0
    for _ in range(n_runs):
        add = min(n_samples - cnt, n)
        idx[cnt:cnt+add] = np.random.choice(n, add, replace=False).astype(int)
        cnt += add
    return X[idx]

# %%



def main():
    # ======================
    # 1. SDE simulations (sparse vs full)
    # ======================
    dim = 2
    t0, t1 = 0.0, 1.0
    dt = 1e-3
    sigma = 3.0
    n_times = 10
    snapshot_times = np.linspace(t0, t1, n_times)
    n_cells = 100

    # Sparse: 100 at endpoints, 10 at intermediate times
    n_sparse = int(n_cells) * np.ones(len(snapshot_times), dtype=int)
    n_sparse[1:-1] = 10

    # Full: 100 at all timepoints
    n_full = n_cells * np.ones(len(snapshot_times), dtype=int)

    # Simulate sparse data (used for inference)
    snapshots_sparse = simulate_sde(
        n_particles=n_sparse,
        dim=dim,
        t0=t0,
        t1=t1,
        dt=dt,
        sigma=sigma,
        snapshot_times=snapshot_times,
        seed=42,
    )

    # Simulate full data (reference for evaluation)
    snapshots_full = simulate_sde(
        n_particles=n_full,
        dim=dim,
        t0=t0,
        t1=t1,
        dt=dt,
        sigma=sigma,
        snapshot_times=snapshot_times,
    )

    times = sorted(snapshots_sparse.keys())
    #print(f"Times: {times}")
    #print(f"Shapes sparse: {[snapshots_sparse[t].shape for t in times]}")
    #print(f"Shapes full:   {[snapshots_full[t].shape for t in times]}")

    
    # ======================
    # 2. Define and optimize TrajLoss on sparse data
    # ======================
    lam_reg = 10
    lam_fit = 1
    eps = sigma / n_times       
    sigma_fit = 0.2 
    n_epochs = 300
    lr = 1e-2

    
    # ======================
    # 3. Prepare data for gWOT
    # ======================
    device = "cuda" if torch.cuda.is_available() else "cpu"
    M = 100  # number of model particles per time

    obs_list_sparse = []   # data used for inference (built from sparse simulation)
    x0_list = []           # initial model particles (100 per time)

    for i, t in enumerate(times):
        X_sparse = snapshots_sparse[t]  # (n_sparse_t, dim)
        # We want 100 obs points per time for the data-fit term, based on sparse samples
        X_obs_sparse_resampled = resample_points(X_sparse, n_samples=M)
        obs_list_sparse.append(torch.from_numpy(X_sparse).float())

        # Initialize model particles close to these resampled sparse observations
        x0 = X_obs_sparse_resampled + np.random.randn(*X_obs_sparse_resampled.shape) * (t > 0)
        x0_list.append(torch.from_numpy(x0).float())
        

    
    traj_model = TrajLoss(
        x0_list=x0_list,
        obs_list=obs_list_sparse,
        lam_reg=lam_reg,
        lam_fit=lam_fit,
        eps=eps,
        sigma_fit=sigma_fit,
        device=device,
    ).to(device)
    
    history, best_positions = optimize_model(
        traj_model,
        n_epochs=n_epochs,
        lr=lr,
        print_every=50,
    )

    print([b.shape for b in best_positions])
    
    # ======================
    # 4. Compare marginals (inferred vs full, sparse vs full)
    # ======================
    distances_inferred_vs_full = []
    distances_sparse_vs_full = []

    for i, t in enumerate(times):
        # Reference: full simulation, 100 points at each time
        X_full = snapshots_full[t]          # (100, dim)

        # Inferred: best_positions[i] from gWOT (100, dim)
        X_inferred = best_positions[i]

        # Sparse: original sparse snapshot (n_sparse_t, dim)
        X_sparse = snapshots_sparse[t]

        # Distances:
        d_inf_full = distribution_distance(X_inferred, X_full, epsilon=0.0)
        distances_inferred_vs_full.append(d_inf_full)

        d_sparse_full = distribution_distance(X_sparse, X_full, epsilon=0.0)
        distances_sparse_vs_full.append(d_sparse_full)

    # ======================
    # 5. Figure: loss (left) + distances (right)
    # ======================
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    # Left: loss vs iterations
    ax_loss = axes[0]
    ax_loss.plot(history)
    ax_loss.set_xlabel("Iteration")
    ax_loss.set_ylabel("gWOT loss")
    ax_loss.set_title("Convergence of gWOT (Adam)")

    # Right: distances vs time
    ax_dist = axes[1]
    ax_dist.plot(times, distances_inferred_vs_full, marker="o", label="Inferred vs full")
    ax_dist.plot(times, distances_sparse_vs_full, marker="s", label="Sparse vs full")
    ax_dist.set_xlabel("Time")
    ax_dist.set_ylabel("W2 distance")
    ax_dist.set_title("Reconstruction quality per timepoint")
    ax_dist.legend()

    fig.tight_layout()
    plt.show()

    # ======================
    # 6. Figure: sparse vs inferred vs full in phase space
    # ======================
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))

    # Helper to build stacked positions and times
    def stack_snapshots(snapshot_dict):
        all_X = []
        all_t = []
        for t in times:
            X = snapshot_dict[t]  # np.ndarray (n_t, dim)
            all_X.append(X)
            all_t.append(t * np.ones(len(X)))
        all_X = np.vstack(all_X)
        all_t = np.concatenate(all_t)
        return all_X, all_t

    # Left: sparse data (variable n_t)
    all_X_sparse, all_t_sparse = stack_snapshots(snapshots_sparse)
    sc0 = axes[0].scatter(
        all_X_sparse[:, 0],
        all_X_sparse[:, 1],
        c=all_t_sparse,
        cmap="viridis",
        s=4,
        alpha=0.5,
    )
    axes[0].set_title("Sparse data")
    axes[0].set_xlabel("x₀")
    axes[0].set_ylabel("x₁")
    axes[0].axis("equal")

    # Middle: inferred data (100 per time, best_positions)
    all_X_inf = []
    all_t_inf = []
    for i, t in enumerate(times):
        X_inf = best_positions[i]  # (100, dim)
        all_X_inf.append(X_inf)
        all_t_inf.append(t * np.ones(len(X_inf)))
    all_X_inf = np.vstack(all_X_inf)
    all_t_inf = np.concatenate(all_t_inf)

    sc1 = axes[1].scatter(
        all_X_inf[:, 0],
        all_X_inf[:, 1],
        c=all_t_inf,
        cmap="viridis",
        s=4,
        alpha=0.5,
    )
    axes[1].set_title("Inferred data (gWOT)")
    axes[1].set_xlabel("x₀")
    axes[1].set_ylabel("x₁")
    axes[1].axis("equal")

    # Right: full data (100 per time)
    all_X_full, all_t_full = stack_snapshots(snapshots_full)
    sc2 = axes[2].scatter(
        all_X_full[:, 0],
        all_X_full[:, 1],
        c=all_t_full,
        cmap="viridis",
        s=4,
        alpha=0.5,
    )
    axes[2].set_title("Full data")
    axes[2].set_xlabel("x₀")
    axes[2].set_ylabel("x₁")
    axes[2].axis("equal")

    fig.tight_layout()
    plt.show()
    

if __name__ == "__main__":
    main()
    
# %%

def run_experiment(
    dim=2,
    t0=0.0,
    t1=1.0,
    dt=1e-3,
    sigma=3.0,
    n_times=10,
    n_cells=100,
    lam_reg=0.05,
    lam_fit=1,
    sigma_fit=0.2,
    n_epochs=100,
    lr=1e-2,
    seed=None,
):
    # ======================
    # 1. Simulate data
    # ======================
    snapshot_times = np.linspace(t0, t1, n_times)

    n_sparse = n_cells * np.ones(len(snapshot_times), dtype=int)
    n_sparse[1:-1] = 10  # sparse in the middle

    n_full = n_cells * np.ones(len(snapshot_times), dtype=int)

    snapshots_sparse = simulate_sde(
        n_particles=n_sparse,
        dim=dim,
        t0=t0,
        t1=t1,
        dt=dt,
        sigma=sigma,
        snapshot_times=snapshot_times,
        seed=seed,
    )

    snapshots_full = simulate_sde(
        n_particles=n_full,
        dim=dim,
        t0=t0,
        t1=t1,
        dt=dt,
        sigma=sigma,
        snapshot_times=snapshot_times,
    )

    times = sorted(snapshots_sparse.keys())

    # ======================
    # 2. Prepare gWOT model
    # ======================
    device = "cuda" if torch.cuda.is_available() else "cpu"
    M = n_cells
    eps = sigma / n_times

    obs_list_sparse = []
    x0_list = []

    for t in times:
        X_sparse = snapshots_sparse[t]
        X_resampled = resample_points(X_sparse, n_samples=M)

        obs_list_sparse.append(torch.from_numpy(X_sparse).float())

        x0 = X_resampled + np.random.randn(*X_resampled.shape) * (t > 0)
        x0_list.append(torch.from_numpy(x0).float())

    traj_model = TrajLoss(
        x0_list=x0_list,
        obs_list=obs_list_sparse,
        lam_reg=lam_reg,
        lam_fit=lam_fit,
        eps=eps,
        sigma_fit=sigma_fit,
        device=device,
    ).to(device)

    _, best_positions = optimize_model(
        traj_model,
        n_epochs=n_epochs,
        lr=lr,
        print_every=0,
    )

    # ======================
    # 3. Compute summed W2 distance
    # ======================
    total_w2 = []

    for i, t in enumerate(times):
        X_full = snapshots_full[t]
        X_inferred = best_positions[i]

        d = distribution_distance(X_inferred, X_full, epsilon=0.0)
        total_w2 += [d]
    
    return np.mean(total_w2)


# ======================
# 4. Run 10 repetitions
# ======================
if __name__ == "__main__":

    lam = np.array([1e-2,1e-1,1,1e1])
    results = []

    for l in lam:
        w2_runs = []
        
        for k in tqdm(range(5)):
            w2 = run_experiment(lam_reg=l, seed=42+k)
            w2_runs.append(w2)
        print(l, np.mean(w2_runs), np.std(w2_runs))
# %%


