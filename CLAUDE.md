# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a computational physics project that computes **transport cross-sections from intermolecular potentials** with rotational (angular) dependence. The workflow lives in plain `.py` modules and a set of input data files (`PC.in.XXX`):

- `transport_from_potential.py` — single-angle pipeline: 1D `Potential` class (short/long extrapolation + cubic spline in r), turning-point solve, deflection angle, and cross-section integral.
- `potential_surface.py` *(added for angle-resolved work)* — 2D `PotentialSurface(V(r, phi))` class over the full angle grid.
- `v_tilde.py` *(added for angle-resolved work)* — the effective 1D potential V̂tilde(r) that collapses the angle-resolved deflection integral back into the atom–atom form so `CrossSectionSolver` can be reused unchanged.
- `prototypes/` — earlier 2D-interpolation prototype and its test harness, kept for reference (not part of the production pipeline).

The physical pipeline is:
1. Read angle-dependent pair potentials V(r, phi) from input files
2. Extrapolate potentials to short range (power-law fit) and long range
3. Determine orbiting conditions (where classical trajectories spiral indefinitely)
4. Compute deflection angles chi(epsilon, b) via numerical integration
5. Compute transport cross-sections Q^(ell)(epsilon) from deflection angles

## Input Data Format

Files are named `PC.in.XXX` where `XXX` is the angle in degrees (000 through 180, in steps of 10). Each file:
- Line 1: number of data points
- Subsequent lines: tab-separated `r V(r)` pairs (distance and potential energy)
- Distances are in atomic units, not necessarily equally spaced

## Running

Requires Python 3 with: `numpy`, `scipy`, `matplotlib`. Run the modules and tests directly:

```
python transport_from_potential.py          # single-angle demo driver
pytest test_cross_sections.py               # single-angle regression tests
pytest test_potential_surface.py            # 2D V(r, phi) tests (angle-resolved)
pytest test_v_tilde.py                      # V̂tilde construction tests (angle-resolved)
python prototypes/test_interpolation_2d.py  # prototype angular-interp comparison
```

## Key Physics / Code Concepts

- **Short-range extrapolation**: Fits first two data points to `V(r) = C * r^k` (power law)
- **Long-range behavior**: Controlled by `long_range_pow` parameter (e.g., r^{-6})
- **Cubic spline interpolation**: `scipy.interpolate.CubicSpline` used between data points; derivatives computed analytically from the spline
- **Orbiting detection**: Checks where effective potential has a maximum (centrifugal barrier), indicating classical orbiting. Uses conditions on V, V', V'' at separation r
- **Deflection angle integral**: Uses substitution w = b/r to transform the integral to finite limits, then numerical quadrature (trapezoidal rule on uniform grid)
- **Cross-section integral**: Split into [0,1] and [1, infinity) pieces; the high-b tail uses an approximation
- **Angular symmetry**: `PotentialSurface` takes a required `symmetry` argument. `"heteronuclear"` expects data on [0°, 180°]; `"homonuclear"` expects data on [0°, 90°]. The class handles folding φ back into the canonical range internally.
- **V̂tilde construction**: For angle-resolved collisions, `v_tilde.build_collision_frame(...)` sets up a gauge-fixed frame at closest approach (**r_m** along ẑ, **v_m** along x̂; ℓ̂_m and **L** parametrized by physical angles ψ and α_L) and `make_v_tilde(surface, frame)` returns a 1D closure that plugs into the existing single-angle `CrossSectionSolver`.

## Numerical Considerations

- Root-finding uses `scipy.optimize.root_scalar` with bracketing methods
- The deflection integrand has a square-root singularity at w = b/r_m (the turning point); the current implementation uses uniform grid quadrature which can lose accuracy near this singularity
- Orbiting parameters are found by scanning over separations, then refined with `minimize_scalar` and `root_scalar`
