# TODO List

## Input Validation
- [ ] Check that R values increase monotonically in each input file
- [ ] Verify potential is well-approximated by R^(long_range_pow) for large R
- [ ] Check that the first three potential energies (small R) are all positive
- [ ] Check that the last three potential energies (large R) are all negative
- [ ] Check that the largest given potential energy equals or exceeds e_max

## Code Improvements
- [ ] Refactor `create_extrapolated_potential`, `create_derivative_potential`, and `create_second_derivative_potential` to share a single cubic spline object instead of each independently rebuilding the same spline
- [ ] Create a function that calls find_orbiting, find_energy_decrease, find_energy_decrease_precisely, find_energy_value_precisely and returns an interpolated mapping of each energy to a separation
- [ ] Add documentation for `C_ell` function

## Physics
- [ ] Support multiple orbiting (potentials with more than one centrifugal barrier maximum)

## Numerical Methods
- [ ] Investigate using Clenshaw-Curtis integration (available via `scipy.integrate.nquad`)
- [ ] Investigate higher precision arithmetic (e.g., `heyoka.py` for arbitrary precision)
