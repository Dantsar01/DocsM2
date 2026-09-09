# LPT Development Log

---

## 2026-05-08 (Thursday)

### Session 1 — GP_CIP tuning + exact solution infrastructure

**Parameter tuning (`1D_Euler_LPTpy.py`)**
- GP_CIP: changed `r = 5 → 4`, `nugget = 1e-7 → 1e-10`
- Initial conditions updated: `rho_0 = 1.5 + exp(-x²)`, `u_0 = -0.2*exp(-x²)`
- Grid bumped from `Nx=128` to `Nx=256`, `T_max=0.5 → 1.0`

**Accuracy / error estimation — discussion**
Options considered for measuring LPT scheme accuracy across interpolation methods:
1. Spatial convergence study (L2 vs HLLC reference, vary Nx)
2. Exact solution comparison via simple wave (avoids HLLC reference error)
3. Interpolation unit test (freeze grid, test interp in isolation)
4. Conservation diagnostics (track Σρ·dx drift each step)
5. Temporal vs spatial separation (fix Nx, vary CFL)

**New file: `exact_simple_wave.py`**
Exact solution for 1D isentropic Euler via right-going simple wave (R⁻ = const).
- `make_simple_wave_ic(x, rho_bg, K, gamma, perturbation)` — build consistent (ρ₀, u₀)
- `SimpleWaveExact` — traces λ⁺ characteristics forward, inverts x(t)=x₀+λ⁺t for each query point
  - `.shock_time()` — estimates t* = 1/|min dλ⁺/dx₀|, raises error if t ≥ t*
  - `.solve(x_eval, t)` — returns (ρ, u) exact
  - `.l2_error(x, rho_num, u_num, t)` — L2 vs exact
- `convergence_study(...)` + `plot_convergence()` — automated log-log convergence table

---

## 2026-05-07 (Wednesday)

- Implemented **spatially variable σ** in `gp_interpolate` (Fornberg & Zuev 2007)
- Merged into main via worktree `romantic-heyrovsky-bb13be` / `pedantic-euler-2c967d`

---

## 2026-05-05 (Monday)

- Refactored `1D_Euler_LPTpy.py` into classes: `HLLCSolver`, `LPTSolver`, `EulerSimulation`
  (worktrees: `nifty-fermat`, `quizzical-williams`, `stupefied-mccarthy`)
- Added `pdfs/` folder with analysis documents (worktree `suspicious-dewdney`)

---

## 2026-05-03 (Saturday)

- Diagnosed and fixed structural issue: plateau regions causing wrong shocks
  (worktree `nifty-beaver`)

---

## TODO

- [ ] **Add Riemann Problem (RP) test case**
  Standard shock-tube IC (e.g. Sod) to test LPT in the presence of discontinuities.
  Plug into `exact_simple_wave.py` or a separate `exact_riemann.py` solver for reference.

- [ ] **2×2 → 3×3 system (full isentropic → full Euler)**
  Add energy equation: state becomes (ρ, ρu, E).
  Third Riemann invariant (entropy) + entropy wave (λ⁰ = u) alongside λ±.
  Requires a third characteristic grid tracking the contact wave.

- [ ] **RK4 time integration for GP-CIP**
  Currently: forward-Euler advection of characteristics + single interp step.
  Replace with classical 4-stage RK4 on the characteristic ODE (dx/dt = λ±)
  to reduce temporal error and let the GP/CIP spatial accuracy dominate.
