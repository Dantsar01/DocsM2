#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 15:18:10 2026

@author: dantsar
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2D semi-Lagrangian (particle) method for linear advection
    u_t + a(x,y) · ∇u = 0    on [0,1]^2, periodic BCs

Particles are translated along characteristics, then snapped back to the
background Cartesian grid (grid-based collision / merging).  Empty cells
are replenished by nearest-neighbour interpolation.

Two Eulerian reference methods run in parallel: upwind (Godunov) and
Rusanov (local Lax-Friedrichs).  Both support spatially varying a(x,y).

Velocity modes (set VEL_MODE below)
------------------------------------
  "x"        : uniform rightward    a = (1, 0)
  "y"        : uniform upward       a = (0, 1)
  "diagonal" : uniform diagonal     a = (1, 1) / sqrt(2)
  "rotation" : solid-body rotation  a = (-(y-0.5), x-0.5)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata


# =============================================================================
# Velocity field   a : (x, y) → (ax, ay)
# =============================================================================

def velocity_field(x, y, mode="x"):
    """
    Evaluate the advection velocity a(x,y) at positions (x, y).
    Both x and y may be arrays of any shape; returned arrays are same shape.

    To add a custom field, add an elif branch here.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if mode == "x":
        return np.ones_like(x), np.zeros_like(y)

    elif mode == "y":
        return np.zeros_like(x), np.ones_like(y)

    elif mode == "diagonal":
        c = 1.0 / np.sqrt(2.0)
        return np.full_like(x, c), np.full_like(y, c)

    elif mode == "rotation":
        # Solid-body rotation around the domain centre (0.5, 0.5), |ω| = 1
        cx, cy = 0.5, 0.5
        return -(y - cy), (x - cx)

    else:
        raise ValueError(f"Unknown velocity mode: {mode!r}")


# =============================================================================
# Domain & background Cartesian grid
# =============================================================================

xmin, xmax = 0.0, 1.0
ymin, ymax = 0.0, 1.0
Nx, Ny     = 64, 64

dx = (xmax - xmin) / Nx
dy = (ymax - ymin) / Ny

x_ref = np.linspace(xmin + dx / 2, xmax - dx / 2, Nx)   # cell centres (Nx,)
y_ref = np.linspace(ymin + dy / 2, ymax - dy / 2, Ny)   # cell centres (Ny,)
X_ref, Y_ref = np.meshgrid(x_ref, y_ref, indexing="ij")  # (Nx, Ny)


# =============================================================================
# Initial condition
# =============================================================================

def u0_2d(x, y):
    return np.sin(2 * np.pi * x) * np.sin(2 * np.pi * y)


# =============================================================================
# Periodic boundary conditions (particles)
# =============================================================================

def wrap_periodic_2d(x, y):
    """Map positions back into [xmin, xmax) × [ymin, ymax) periodically."""
    Lx = xmax - xmin
    Ly = ymax - ymin
    return xmin + np.mod(x - xmin, Lx), ymin + np.mod(y - ymin, Ly)


# =============================================================================
# Grid-based collision: snap particles to cell centres & merge duplicates
#
# Replaces the 1D remove_collisions + argsort logic.
#
# After translating, each particle (x_p, y_p) is assigned to a grid cell
# (ix, iy).  Multiple particles landing in the same cell are merged into
# one: their u values are averaged, and a single particle is placed at the
# cell centre.
# =============================================================================

def snap_and_merge(x_p, y_p, u_p):
    """
    Map each particle to its grid cell; average u for colliding particles.
    Returns one particle per occupied cell, sitting at the cell centre.
    """
    ix = np.clip(np.floor((x_p - xmin) / dx).astype(int), 0, Nx - 1)
    iy = np.clip(np.floor((y_p - ymin) / dy).astype(int), 0, Ny - 1)

    cell_id = ix * Ny + iy                            # flat (unique) cell index

    unique_ids, inverse = np.unique(cell_id, return_inverse=True)

    # Sum u within each occupied cell, then divide by the count
    u_merged = np.zeros(len(unique_ids))
    np.add.at(u_merged, inverse, u_p)
    u_merged /= np.bincount(inverse)

    # Place surviving particle at the cell centre
    ix_u = unique_ids // Ny
    iy_u = unique_ids  % Ny
    x_out = x_ref[ix_u]
    y_out = y_ref[iy_u]

    return x_out, y_out, u_merged


# =============================================================================
# Replenishment: inject a new particle into every empty cell
#
# Replaces the 1D mesh_interpolation.
#
# For non-uniform velocity, particles bunch into some cells and leave
# others empty.  We fill empty cells by nearest-neighbour interpolation
# from currently occupied cells so the particle count stays ≈ Nx × Ny.
# =============================================================================

def replenish_empty_cells(x_p, y_p, u_p):
    """
    Find cells that hold no particle and fill them with the nearest
    occupied neighbour's u value (nearest-neighbour interpolation).
    """
    ix = np.clip(np.floor((x_p - xmin) / dx).astype(int), 0, Nx - 1)
    iy = np.clip(np.floor((y_p - ymin) / dy).astype(int), 0, Ny - 1)

    occupied = np.zeros((Nx, Ny), dtype=bool)
    occupied[ix, iy] = True

    if occupied.all():
        return x_p, y_p, u_p                # nothing to do

    empty_ix, empty_iy = np.where(~occupied)
    x_empty = x_ref[empty_ix]
    y_empty = y_ref[empty_iy]

    # Nearest-neighbour interpolation (exact for uniform constant flow,
    # first-order approximation otherwise)
    u_fill = griddata(
        (x_p, y_p), u_p,
        (x_empty, y_empty),
        method="nearest"
    )

    return (np.concatenate([x_p, x_empty]),
            np.concatenate([y_p, y_empty]),
            np.concatenate([u_p, u_fill]))


# =============================================================================
# Particle translation step   (2D analogue of mesh_translation_step)
# =============================================================================

def particle_step_2d(x_p, y_p, u_p, dt, vel_mode="x"):
    """
    One time step of the particle (semi-Lagrangian) method:
      1. Translate each particle along its characteristic: x += a(x,y)*dt
      2. Wrap periodically back into the domain
      3. Snap to grid & merge colliding particles
      4. Replenish empty cells
    """
    # 1 — translate
    
    ax, ay = velocity_field(x_p, y_p, mode=vel_mode)
    
    x_p = x_p + ax * dt
    y_p = y_p + ay * dt
    '''
    
    # RK2 (midpoint method)
    ax1, ay1 = velocity_field(x_p, y_p, mode=vel_mode)

    x_mid = x_p + 0.5 * dt * ax1
    y_mid = y_p + 0.5 * dt * ay1

    ax2, ay2 = velocity_field(x_mid, y_mid, mode=vel_mode)
    
    x_p = x_p + dt * ax2
    y_p = y_p + dt * ay2
    '''
    # 2 — periodic wrap
    x_p, y_p = wrap_periodic_2d(x_p, y_p)

    # 3 — grid-based collision handling (replaces 1D remove_collisions)
    x_p, y_p, u_p = snap_and_merge(x_p, y_p, u_p)

    # 4 — replenish empty cells (replaces 1D mesh_interpolation)
    #x_p, y_p, u_p = replenish_empty_cells(x_p, y_p, u_p)

    return x_p, y_p, u_p


# =============================================================================
# CFL time step   (2D analogue of compute_dt_particle)
# =============================================================================

def compute_dt_2d(x_p, y_p, vel_mode, cfl=1):
    """
    dt ≤ CFL / (|ax|_max/dx + |ay|_max/dy)
    This is the multi-dimensional CFL stability condition.
    """
    ax, ay = velocity_field(x_p, y_p, mode=vel_mode)
    speed_x = np.max(np.abs(ax)) / dx
    speed_y = np.max(np.abs(ay)) / dy
    return cfl / (speed_x + speed_y + 1e-15)


# =============================================================================
# 2D Upwind (Godunov) for   u_t + ax(x,y)*u_x + ay(x,y)*u_y = 0
#
# Dimension-by-dimension upwind; supports spatially varying velocity.
# Periodic BCs via np.roll.
# =============================================================================

def upwind_2d_step(u, Ax, Ay, dt):
    """
    Ax, Ay : (Nx, Ny) arrays of velocity components at cell centres.
    Uses backward / forward differences depending on the sign of each
    component, chosen cell-by-cell.
    """
    # --- x direction upwind
    flux_x_pos = Ax * (u - np.roll(u,  1, axis=0)) / dx   # ax > 0: backward diff
    flux_x_neg = Ax * (np.roll(u, -1, axis=0) - u) / dx   # ax < 0: forward  diff
    flux_x = np.where(Ax >= 0, flux_x_pos, flux_x_neg)

    # --- y direction upwind
    flux_y_pos = Ay * (u - np.roll(u,  1, axis=1)) / dy
    flux_y_neg = Ay * (np.roll(u, -1, axis=1) - u) / dy
    flux_y = np.where(Ay >= 0, flux_y_pos, flux_y_neg)

    return u - dt * (flux_x + flux_y)


# =============================================================================
# 2D Rusanov (local Lax-Friedrichs) for linear advection
#
# Interface flux:
#   F_{i+½} = ½ a_int (uL + uR) − ½ |a_int| (uR − uL)
# where a_int is the average of the two adjacent cell velocities.
# =============================================================================

def rusanov_2d_step(u, Ax, Ay, dt):
    """
    Ax, Ay : (Nx, Ny) velocity arrays.
    """
    # --- x interfaces (i+½)
    uL_x   = u
    uR_x   = np.roll(u, -1, axis=0)
    ax_int = 0.5 * (Ax + np.roll(Ax, -1, axis=0))
    Fx     = 0.5 * ax_int * (uL_x + uR_x) - 0.5 * np.abs(ax_int) * (uR_x - uL_x)
    div_x  = (Fx - np.roll(Fx, 1, axis=0)) / dx

    # --- y interfaces (j+½)
    uL_y   = u
    uR_y   = np.roll(u, -1, axis=1)
    ay_int = 0.5 * (Ay + np.roll(Ay, -1, axis=1))
    Fy     = 0.5 * ay_int * (uL_y + uR_y) - 0.5 * np.abs(ay_int) * (uR_y - uL_y)
    div_y  = (Fy - np.roll(Fy, 1, axis=1)) / dy

    return u - dt * (div_x + div_y)


# =============================================================================
# Exact solution
#
# For constant-velocity modes: u(x,y,t) = u0(x − ax*t, y − ay*t)  (periodic).
# For rotation: back-rotate each point by angle t around (0.5, 0.5).
# =============================================================================

def u_exact_2d(t, vel_mode):
    """Evaluate the exact solution on the background grid at time t."""
    if vel_mode == "rotation":
        cx, cy = 0.5, 0.5
        cos_t, sin_t = np.cos(t), np.sin(t)
        # Characteristic back-tracking: rotate by −t
        X0 = cx + (X_ref - cx) * cos_t + (Y_ref - cy) * sin_t
        Y0 = cy - (X_ref - cx) * sin_t + (Y_ref - cy) * cos_t
        X0 = xmin + np.mod(X0 - xmin, xmax - xmin)
        Y0 = ymin + np.mod(Y0 - ymin, ymax - ymin)
        return u0_2d(X0, Y0)
    else:
        ax, ay = velocity_field(X_ref, Y_ref, mode=vel_mode)
        X0 = xmin + np.mod(X_ref - ax * t - xmin, xmax - xmin)
        Y0 = ymin + np.mod(Y_ref - ay * t - ymin, ymax - ymin)
        return u0_2d(X0, Y0)


# =============================================================================
# Plotting utilities
# =============================================================================

def reconstruct_grid(x_p, y_p, u_p):
    """Place particle u values onto the background grid (NaN where empty)."""
    ix = np.clip(np.floor((x_p - xmin) / dx).astype(int), 0, Nx - 1)
    iy = np.clip(np.floor((y_p - ymin) / dy).astype(int), 0, Ny - 1)
    u_grid = np.full((Nx, Ny), np.nan)
    u_grid[ix, iy] = u_p
    return u_grid


def plot_2d(x_p, y_p, u_p, u_up, u_rus, t, vel_mode, k=0, freq=20):
    if k % freq != 0:
        return

    u_part = reconstruct_grid(x_p, y_p, u_p)
    u_ex   = u_exact_2d(t, vel_mode)

    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    titles = ["Particle method", "Upwind (Godunov)", "Rusanov", "Exact solution"]
    data   = [u_part, u_up, u_rus, u_ex]
    vmin, vmax = -1.0, 1.0

    for ax_plot, title, dat in zip(axes.ravel(), titles, data):
        im = ax_plot.imshow(
            dat.T, origin="lower",
            extent=[xmin, xmax, ymin, ymax],
            vmin=vmin, vmax=vmax, cmap="RdBu_r", aspect="equal"
        )
        ax_plot.set_title(title)
        ax_plot.set_xlabel("x")
        ax_plot.set_ylabel("y")
        plt.colorbar(im, ax=ax_plot)

    fig.suptitle(f"t = {t:.4f}  |  mode = {vel_mode!r}")
    plt.tight_layout()
    plt.show()


def plot_l2_error(t_hist, err_part, err_up, err_rus, vel_mode):
    plt.figure(figsize=(7, 4))
    plt.semilogy(t_hist, err_part, 'r-',  label="Particle method")
    plt.semilogy(t_hist, err_up,   'b--', label="Upwind (Godunov)")
    plt.semilogy(t_hist, err_rus,  'g:',  label="Rusanov")
    plt.xlabel("t")
    plt.ylabel("L² error")
    plt.title(f"L² error over time  |  mode = {vel_mode!r}")
    plt.legend()
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.tight_layout()
    plt.show()


# =============================================================================
# Main loop
# =============================================================================

VEL_MODE  = "rotation"    # ← change here:  "x" | "y" | "diagonal" | "rotation"
T         = 1.0    # end time
PLOT_FREQ = 10    # plot every PLOT_FREQ time steps

# --- Velocity grids for Eulerian methods (time-independent, evaluated once)
ax_max, ay_max = velocity_field(X_ref, Y_ref, mode=VEL_MODE)

# --- Particle initialisation: one particle per cell centre
x_p = X_ref.ravel().copy()
y_p = Y_ref.ravel().copy()
u_p = u0_2d(x_p, y_p)

# --- Eulerian initialisations
u_up  = u0_2d(X_ref, Y_ref).copy()
u_rus = u0_2d(X_ref, Y_ref).copy()

t = 0.0
k = 0

t_hist, err_part_hist, err_up_hist, err_rus_hist = [], [], [], []


while t <= T:

    # --- Shared time step (2D CFL)
    dt = compute_dt_2d(ax_max, ay_max, VEL_MODE)
   
    # --- Particle step
    x_p, y_p, u_p = particle_step_2d(x_p, y_p, u_p, dt, vel_mode=VEL_MODE)

    # --- Eulerian steps
    u_up  = upwind_2d_step(u_up,  ax_max, ay_max, dt)
    u_rus = rusanov_2d_step(u_rus, ax_max, ay_max, dt)

    t += dt
    k += 1

    # --- L² error tracking
    u_ex = u_exact_2d(t, VEL_MODE)
    u_part_grid = reconstruct_grid(x_p, y_p, u_p)
    # Use 0 for NaN cells when computing error (shouldn't happen after replenish)
    u_part_grid = np.nan_to_num(u_part_grid)

    dA = dx * dy
    err_part = np.sqrt(np.sum((u_part_grid - u_ex)**2) * dA)
    err_up   = np.sqrt(np.sum((u_up         - u_ex)**2) * dA)
    err_rus  = np.sqrt(np.sum((u_rus         - u_ex)**2) * dA)

    t_hist.append(t)
    err_part_hist.append(err_part)
    err_up_hist.append(err_up)
    err_rus_hist.append(err_rus)

    # --- Plot
    plot_2d(x_p, y_p, u_p, u_up, u_rus, t, VEL_MODE, k=k, freq=PLOT_FREQ)

# --- Final L² error plot
plot_l2_error(t_hist, err_part_hist, err_up_hist, err_rus_hist, VEL_MODE)

# %%
