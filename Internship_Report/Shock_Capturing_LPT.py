#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  7 15:44:44 2026

@author: dantsar
"""

import numpy as np
import matplotlib.pyplot as plt

def f_burgers(u):
    return u**2/2

def f_p_burgers(u):
    return u

def u0(ul, ur, x):
    x = np.asarray(x)
    y = np.ones_like(x, dtype=float)

    mask_left = x < 1/4
    mask_right = x > 3/4
    mask_middle = (~mask_left) & (~mask_right)

    # Left and right states
    y[mask_left] = ul
    y[mask_right] = ur

    # Affine (linear) part between 1/4 and 3/4
    y[mask_middle] = ul + (ur - ul) * (x[mask_middle] - 1/4) / (3/4 - 1/4)

    return y

def u0_bis(ul,uc, ur, x):
    x = np.asarray(x)
    y = np.ones_like(x, dtype=float)

    mask_left = x < 1/4
    mask_right = x > 3/4
    mask_middle = (~mask_left) & (~mask_right)
    
    y[mask_left] = ul
    y[mask_right] = ur
    y[mask_middle] = uc
    
    return y

# %%

def apply_bc(u, bc_type):
    if bc_type == "periodic":
        u_ext = np.concatenate(([u[-1]], u, [u[0]]))

    elif bc_type == "neumann":  # zero gradient
        u_ext = np.concatenate(([u[0]], u, [u[-1]]))

    elif bc_type == "dirichlet":
        # example: fixed values (modify if needed)
        uL_val = u[0]
        uR_val = u[-1]
        u_ext = np.concatenate(([uL_val], u, [uR_val]))

    else:
        raise ValueError("Unknown BC type")

    return u_ext

def godunov_flux(uL, uR):
    # Exact Riemann solver flux for Burgers (Upwind/Godunov)
    flux = np.zeros_like(uL)

    mask = uL <= uR  # rarefaction
    flux[mask] = np.where(
        uL[mask] >= 0, f_burgers(uL[mask]),
        np.where(uR[mask] <= 0, f_burgers(uR[mask]), 0)
    )

    mask = uL > uR  # shock
    s = (uL[mask] + uR[mask]) / 2
    flux[mask] = np.where(s >= 0, f_burgers(uL[mask]), f_burgers(uR[mask]))

    return flux

def upwind_step(u, dt, dx, bc_type="periodic"):
    u_ext = apply_bc(u, bc_type)

    uL = u_ext[:-1]
    uR = u_ext[1:]

    F = godunov_flux(uL, uR)

    return u - dt/dx * (F[1:] - F[:-1])


def rusanov_step(u, dt, dx, bc_type="periodic"):
    u_ext = apply_bc(u, bc_type)

    uL = u_ext[:-1]
    uR = u_ext[1:]

    alpha = np.maximum(np.abs(uL), np.abs(uR))

    F = 0.5 * (f_burgers(uL) + f_burgers(uR)) \
        - 0.5 * alpha * (uR - uL)

    return u - dt/dx * (F[1:] - F[:-1])


# %%


def mesh_out(x_next):
    right = np.sum(x_next >= xmax)
    left =  np.sum(x_next <= xmin)
    return left, right


#def u_ex(ul,ur,t, x):
    #return u0(ul, ur, x - (ur + ul)/2 * t)
    

def u_ex(ul, ur, t, x):
    x = np.asarray(x, dtype=float)
    u = np.zeros_like(x)

    a = 2 * (ur - ul)        # slope
    T = -1 / a              # shock time (a < 0 required)

    if t < T:
        mask_left = x < 1/4 + ul * t
        mask_right = x > 3/4 + ur * t

        u[mask_left] = ul
        u[mask_right] = ur

        mask_middle = (~mask_left) & (~mask_right)

        x0 = (x[mask_middle] - (ul - a/4) * t) / (1 + a * t)
        u[mask_middle] = ul + a * (x0 - 1/4)

    else:
        # shock solution
        s = (ul + ur)/2

        x_s_T = 1/4 - ul / a     # <<< THIS WAS MISSING
        x_s = x_s_T + s * (t - T)

        u = np.where(x < x_s, ul, ur)

    return u


# %%

def plot_(x, u, t, k, sol = False, freq = 1):
    if k%freq == 0 :
        plt.figure()
        plt.grid()
        plt.plot(x,u, 'r+', label = f"numerical at t = {t}")
        if sol == True:
            u_sol = u_ex(ul,ur,t,x)
            plt.plot(x, u_sol, 'b-', label = 'ground truth')
        plt.xlim([xmin,xmax])
        plt.ylim([-1,1])
        plt.legend()
        plt.show()

# %%
'''
def remove_collisions(
    x, u, f, remove_right_if, tol=1e-6
):
    x = np.asarray(x)
    u = np.asarray(u)

    n = len(x)

    dx = x[1:] - x[:-1]
    du = u[:-1] - u[1:]

    # collision after translation
    collide = (dx == 0)

    # Rankine–Hugoniot speed for all pairs
    uL = u[:-1]
    uR = u[1:]
    
    
    s  = (f(uL) - f(uR)) / (uL - uR + 1e-15)

    # decide which one to remove (only where collide)
    remove_right = np.zeros(n-1, dtype=bool)
    remove_left  = np.zeros(n-1, dtype=bool)

    decision = remove_right_if(uL, uR)

    remove_right[collide] = decision[collide]
    remove_left[collide]  = ~decision[collide]

    # update surviving particle speeds (RH)
    u_new = u.copy()

    # if right removed → update left
    u_new[:-1][remove_right] = s[remove_right]

    # if left removed → update right
    u_new[1:][remove_left] = s[remove_left]

    # build keep mask
    keep = np.ones(n, dtype=bool)
    keep[1:]  &= ~remove_right
    keep[:-1] &= ~remove_left
    print(keep)
    

    return x[keep], u_new[keep]
'''
# %%
'''
def remove_collisions(x, u, f, remove_right_if, tol=1e-6):
    x = np.asarray(x)
    u = np.asarray(u)
    n = len(x)

    dx = x[1:] - x[:-1]
    collide = np.abs(dx) < tol           # (n-1)

    keep = np.ones(n, dtype=bool)
    u_new = u.copy()

    # ------------------------------------------------------------------
    # Identify collision blocks in dx-space
    # ------------------------------------------------------------------
    block_start = collide & np.concatenate(([True], ~collide[:-1]))
    block_end   = collide & np.concatenate((~collide[1:], [True]))

    block_id = np.cumsum(block_start) - 1
    block_id[~collide] = -1

    valid = block_id >= 0
    block_sizes = np.bincount(block_id[valid])

    size_per_dx = np.zeros_like(block_id)
    size_per_dx[valid] = block_sizes[block_id[valid]]

    # ------------------------------------------------------------------
    # Remove middle particles for blocks with ≥ 3 particles
    # (i.e. block has ≥ 2 consecutive dx == 0)
    # ------------------------------------------------------------------
    multi = size_per_dx >= 2                     # dx-space mask

    # middle particles are (i+1) for dx i that is NOT the last dx of block
    middle_dx = multi & ~block_end               # dx indices
    middle_particles = np.where(middle_dx)[0] + 1
    keep[middle_particles] = False

    # ------------------------------------------------------------------
    # Endpoints of each collision block
    # ------------------------------------------------------------------
    left_dx  = np.where(block_start)[0]
    right_dx = np.where(block_end)[0]

    left_particles  = left_dx
    right_particles = right_dx + 1

    # ------------------------------------------------------------------
    # Rankine–Hugoniot update on surviving endpoints
    # ------------------------------------------------------------------
    uL = u[left_particles]
    uR = u[right_particles]
    s  = 0.5 * (uL + uR) #(f(uR) - f(uL)) / (uR - uL + 1e-15)

    decision = remove_right_if(uL, uR)

    # remove chosen endpoint
    keep[right_particles[decision]] = False
    keep[left_particles[~decision]] = False

    # update speed of surviving endpoint
    u_new[left_particles[decision]]  = s[decision]
    u_new[right_particles[~decision]] = s[~decision]

    return x[keep], u_new[keep]
'''

def remove_collisions(x, u_l, u_r, u_speed, f, tol=1e-10):
    """
    Each particle carries:
      u_l     – left state  (characteristic value, or left side of shock)
      u_r     – right state (characteristic value, or right side of shock)
      u_speed – propagation speed (= u_l = u_r for characteristics,
                                   = RH speed for shocks)

    When a collision block [iL … iR] is detected we keep the leftmost
    particle and update:
      new_u_r[iL]     = u_r[iR]            (right state of rightmost particle)
      new_u_speed[iL] = RH(u_l[iL], u_r[iR])

    This is correct whether the leftmost is a fresh characteristic or an
    existing shock that is absorbing a new particle from the left.
    """
    x       = np.asarray(x)
    u_l     = np.asarray(u_l)
    u_r     = np.asarray(u_r)
    u_speed = np.asarray(u_speed)

    n = len(x)
    if n <= 1:
        return x, u_l, u_r, u_speed

    dx = x[1:] - x[:-1]
    collide = dx < tol   # (n-1)

    # --- Identify block starts and ends
    block_start = np.zeros(n, dtype=bool)
    block_end   = np.zeros(n, dtype=bool)

    block_start[0]  = True
    block_start[1:] = ~collide

    block_end[:-1] = ~collide
    block_end[-1]  = True

    starts = np.where(block_start)[0]
    ends   = np.where(block_end)[0]

    is_collision = ends > starts

    # --- Build output (keep leftmost of each block)
    new_x       = x[starts].copy()
    new_u_l     = u_l[starts].copy()
    new_u_r     = u_r[starts].copy()
    new_u_speed = u_speed[starts].copy()

    # --- Resolve collision blocks
    print(np.any(is_collision))
    if np.any(is_collision):
        iL = starts[is_collision]
        iR = ends[is_collision]

        uL = u_l[iL]           # actual left state of the leftmost particle
        uR = u_r[iR]           # actual right state of the rightmost particle

        s = (f(uR) - f(uL)) / (uR - uL + 1e-15)

        new_u_r[is_collision]     = uR
        new_u_speed[is_collision] = s
        # new_u_l is already correct: u_l[iL] (left state never changes)

    return new_x, new_u_l, new_u_r, new_u_speed

# %%

def mesh_interpolation(x, u):
    x = np.asarray(x)
    u = np.asarray(u)

    gap_left = x[0] - xmin
    gap_right = xmax - x[-1]

    if gap_left <= dx and gap_right <= dx:
        return x, u

    n_add_left = int(gap_left // dx)
    n_add_right = int(gap_right // dx)

    if n_add_left > 0:
        x_new = x[0] - dx * np.arange(n_add_left, 0, -1)
        u_new = np.full(n_add_left, u[0])
        x = np.concatenate((x_new, x))
        u = np.concatenate((u_new, u))

    if n_add_right > 0:
        x_new = x[-1] + dx * np.arange(1, n_add_right + 1)
        u_new = np.full(n_add_right, u[-1])
        x = np.concatenate((x, x_new))
        u = np.concatenate((u, u_new))

    return x, u
    
# %%
'''
def mesh_translation(x, u, f = lambda u: u**2/2 , f_p = lambda u: u):
    
    
    dx = x[1:] - x[:-1]
    du = u[1:] - u[:-1]
    mask_for_one_collision = (du < 0)
    dt_to_shock = dx[mask_for_one_collision] / du[mask_for_one_collision]
    if dt_to_shock.size > 0:
        dt_to_shock = np.min(abs(dt_to_shock))
    if any(mask_for_one_collision) == True and dt_to_shock < dt_ref:
        dt = dt_to_shock
        print('dt_shock = ',dt)
    else:
        dt = dt_ref
        print('dt_ref = ', dt)
    
    x = x + f_p(u) * dt

    left, right = mesh_out(x)

    L = xmax - xmin
    x = xmin + np.mod(x - xmin, L)
    
     
    idx = np.argsort(x)
    x = x[idx]
    u = u[idx]
   
    x, u = remove_collisions(
        x, u, f, 
        #lambda ul, ur: np.logical_and(ur < ul, ul < 0) ,   # remove lower velocity
        lambda ul, ur: ul > ur,
        dt
    )
   

    return x, u, dt
'''

def compute_dt_particle(x, u_speed, CFL = 1):
    dx = x[1:] - x[:-1]
    du = u_speed[1:] - u_speed[:-1]

    mask = (du < 0)
    if np.any(mask):
        dt_shock = np.min(np.abs(dx[mask] / du[mask]))  * CFL
        return min(dt_shock, dt_ref)
    else:
        return dt_ref

def mesh_translation_step(x, u_l, u_r, u_speed, dt, f=lambda u: u**2/2, tol = 1e-6):

    x = x + u_speed * dt

    L = xmax - xmin
    x = xmin + np.mod(x - xmin, L)

    idx = np.argsort(x)
    x       = x[idx]
    u_l     = u_l[idx]
    u_r     = u_r[idx]
    u_speed = u_speed[idx]

    x, u_l, u_r, u_speed = remove_collisions(x, u_l, u_r, u_speed, f, tol)

    return x, u_l, u_r, u_speed
# %%

xmin,xmax = 0., 1.
N = 256
dx = (xmax - xmin)/N

x_ref = np.linspace(xmin + dx/2, xmax - dx/2,N)
x = x_ref.copy()

t = 0 

ul, ur = 1, 0
uc = 0.9


dt_ref = dx

#u = u0_bis(ul, uc, ur, x)
u = np.sin(2 * np.pi * x)+ 1
#u = u0(ul, ur, x)
#u = u0_single_hump_compact(x)

u_up, u_rus = u.copy(), u.copy()
bc_type = "periodic"

# Each particle carries its left state, right state, and propagation speed.
# For characteristics all three start equal; shocks diverge over time.
u_l     = u.copy()
u_r     = u.copy()
u_speed = u.copy()

T = 1
l2_err = []
freq = 10
k = 0

%matplotlib inline
while t <= T:

    # --- Compute BOTH candidate time steps FIRST
    dt_particle = compute_dt_particle(x, u_speed, CFL = 0.99)

    umax = max(np.max(np.abs(u_up)), np.max(np.abs(u_rus)), 1e-8)
    dt_euler = 0.9 * dx / umax

    # --- SINGLE shared dt
    dt = min(dt_particle, dt_euler)
    #dt = dt_particle
    # --- Advance particle method
    x, u_l, u_r, u_speed = mesh_translation_step(x, u_l, u_r, u_speed, dt, tol= 1e-6)

    # --- Advance Eulerian methods
    u_up = upwind_step(u_up, dt, dx, bc_type)
    u_rus = rusanov_step(u_rus, dt, dx, bc_type)

    # --- Plot
    if k % freq == 0:
        plt.figure()
        plt.grid()

        plt.plot(x, u_speed, 'r+', label="Particle method")
        plt.plot(x_ref, u_up, 'b-', label="Upwind (Godunov)")
        plt.plot(x_ref, u_rus, 'g--', label="Rusanov")
        #plt.plot(t, np.sqrt(2) / (10 + t) + 0.5, marker = 'x')

        plt.xlim([0, 1])
        plt.xlim([xmin, xmax])
        plt.legend()
        plt.title(f"t = {t:.4f}")
        plt.show()

    t += dt
    k += 1

# %%
'''
x_next = x + u * dt # translation

sortie = np.sum(x_next >= xmax) # on compte le nombre de point sortant du domaine
#x_next = np.roll(x_next - (x_next // xmax),sortie) # on fait rentrer par periodicite

# y-coordinates for the two discretizations
y = np.vstack((np.zeros_like(x), np.ones_like(x)))

# x-coordinates (your stack)
X = np.vstack((x, x_next))

# Plot all segments at once
plt.plot(X, y, 'k-')

# Optional: show points
plt.plot(x, np.zeros_like(x), 'ko')
plt.plot(x_next, np.ones_like(x_next), 'ko')

plt.xlabel("x")
plt.ylabel("step")
plt.yticks([0, 1], ["original", "translated"])
plt.axis('equal')
plt.show()

print(x_next)

'''
# %%