# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a computational physics project that computes **transport cross-sections from intermolecular potentials** with rotational (angular) dependence. The main workflow lives in a single Jupyter notebook (`transport_from_potental.ipynb`) and a set of input data files (`PC.in.XXX`).

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

Requires Python 3 with: `numpy`, `scipy`, `matplotlib`. Run with Jupyter:

```
jupyter notebook transport_from_potental.ipynb
```

## Key Physics / Code Concepts

- **Short-range extrapolation**: Fits first two data points to `V(r) = C * r^k` (power law)
- **Long-range behavior**: Controlled by `long_range_pow` parameter (e.g., r^{-6})
- **Cubic spline interpolation**: `scipy.interpolate.CubicSpline` used between data points; derivatives computed analytically from the spline
- **Orbiting detection**: Checks where effective potential has a maximum (centrifugal barrier), indicating classical orbiting. Uses conditions on V, V', V'' at separation r
- **Deflection angle integral**: Uses substitution w = b/r to transform the integral to finite limits, then numerical quadrature (trapezoidal rule on uniform grid)
- **Cross-section integral**: Split into [0,1] and [1, infinity) pieces; the high-b tail uses an approximation

## Numerical Considerations

- Root-finding uses `scipy.optimize.root_scalar` with bracketing methods
- The deflection integrand has a square-root singularity at w = b/r_m (the turning point); the current implementation uses uniform grid quadrature which can lose accuracy near this singularity
- Orbiting parameters are found by scanning over separations, then refined with `minimize_scalar` and `root_scalar`
