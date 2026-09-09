The python code for the report is located in `Internship_Report`, the other folders' content are irrelevant.

## Important files

- **`1D_Euler_RBFGA_ImplicitSL.py`** — Implicit semi-Lagrangian solver for the 1D isentropic Euler equations in Riemann-invariant form, using RBF-GA interpolation to evaluate R± at departure points. Supports implicit Euler (`ck=0`) and implicit trapezoidal (`ck=1`) departure schemes.
- **`MGFM_2fluid.py`** — Modified Ghost Fluid Method (MGFM) solver for the 1D isentropic Euler equations with two fluids (Liu, Khoo & Yeo 2003 formulation), coupling a level set, interface Riemann solver, and ghost cell construction.

## Other files

- `Full_Euler_RBFGA_ImplicitSL.py` — draft extension of the RBF-GA implicit semi-Lagrangian solver to the full Euler system (adds the entropy field).
- `MGFM_1D.py` — earlier/base version of the MGFM two-fluid solver.
- `NOTES.md` — development log tracking parameter tuning, design decisions, and TODOs.
- Remaining `.py` files are earlier solver variants (LPT, CIP, RBF-GA) and supporting utilities (kernels, exact solutions, convergence studies, benchmarks).
