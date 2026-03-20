import numpy as np
from scipy.interpolate import CubicSpline, RectBivariateSpline
from scipy.special import legendre


def _build_radial_spline(distances, energies, long_range_pow):
    """Build a single cubic spline for a given angle, with clamped boundary
    conditions matching the short-range and long-range extrapolations.

    Returns (spline, short_coef, short_pow, long_coef) so that the caller
    can evaluate the potential, its first derivative, or its second derivative
    in any region without rebuilding the spline.
    """
    short_pow = (np.log(energies[0]) - np.log(energies[1])) / \
                (np.log(distances[0]) - np.log(distances[1]))
    short_coef = energies[0] / distances[0] ** short_pow
    short_deriv = short_coef * short_pow * distances[0] ** (short_pow - 1)

    long_coef = energies[-1] / distances[-1] ** long_range_pow
    long_deriv = long_coef * long_range_pow * distances[-1] ** (long_range_pow - 1)

    spline = CubicSpline(
        distances, energies, bc_type=((1, short_deriv), (1, long_deriv)))

    return spline, short_coef, short_pow, long_coef


def _eval_radial(r, spline, short_coef, short_pow, long_coef, long_range_pow,
                 r_min, r_max, deriv_order=0):
    """Evaluate the potential (or its 1st/2nd derivative) at separation r,
    using the three-region piecewise definition."""
    if deriv_order == 0:
        if r < r_min:
            return short_coef * r ** short_pow
        elif r > r_max:
            return long_coef * r ** long_range_pow
        else:
            return float(spline(r))
    elif deriv_order == 1:
        if r < r_min:
            return short_coef * short_pow * r ** (short_pow - 1)
        elif r > r_max:
            return long_coef * long_range_pow * r ** (long_range_pow - 1)
        else:
            return float(spline.derivative(1)(r))
    elif deriv_order == 2:
        if r < r_min:
            return short_coef * short_pow * (short_pow - 1) * r ** (short_pow - 2)
        elif r > r_max:
            return long_coef * long_range_pow * (long_range_pow - 1) * r ** (long_range_pow - 2)
        else:
            return float(spline.derivative(2)(r))


def _build_all_radial_splines(angles_given, distances_given, energies_given,
                              long_range_pow):
    """Build a radial spline for every angle. Returns a dict keyed by angle."""
    radial = {}
    for angle in angles_given:
        d = distances_given[angle]
        e = energies_given[angle]
        spline, sc, sp, lc = _build_radial_spline(d, e, long_range_pow)
        radial[angle] = {
            'spline': spline,
            'short_coef': sc,
            'short_pow': sp,
            'long_coef': lc,
            'r_min': min(d),
            'r_max': max(d),
        }
    return radial


# ---------------------------------------------------------------------------
# Method 1: Legendre polynomial expansion in angle
# ---------------------------------------------------------------------------

def _build_legendre(angles_given, distances_given, energies_given,
                    long_range_pow, num_legendre_terms, num_r_points):
    """Fit a Legendre expansion V(r, phi) = sum_l V_l(r) P_l(cos phi).

    At each point on a common radial grid, the 19 angular values are
    projected onto Legendre polynomials. Each radial coefficient V_l(r)
    is then itself stored as a 1D cubic spline (with short/long-range
    extrapolation).
    """
    radial = _build_all_radial_splines(
        angles_given, distances_given, energies_given, long_range_pow)
    angles_sorted = sorted(radial.keys())

    # Common radial grid spanning the union of all data ranges
    r_min_global = min(info['r_min'] for info in radial.values())
    r_max_global = max(info['r_max'] for info in radial.values())
    r_common = np.linspace(r_min_global, r_max_global, num_r_points)

    # Evaluate V(r_i, phi_j) on the common grid
    cos_angles = np.cos(np.radians(angles_sorted))
    V_grid = np.zeros((len(r_common), len(angles_sorted)))
    for j, angle in enumerate(angles_sorted):
        info = radial[angle]
        for i, r in enumerate(r_common):
            V_grid[i, j] = _eval_radial(
                r, info['spline'], info['short_coef'], info['short_pow'],
                info['long_coef'], long_range_pow, info['r_min'], info['r_max'])

    # Build the Legendre projection matrix: P[j, l] = P_l(cos(phi_j))
    n_terms = min(num_legendre_terms, len(angles_sorted))
    P_matrix = np.zeros((len(angles_sorted), n_terms))
    for l in range(n_terms):
        Pl = legendre(l)
        P_matrix[:, l] = Pl(cos_angles)

    # Least-squares fit at each radial point: V_grid[i,:] ≈ P_matrix @ coeffs[i,:]
    # coeffs[i, l] = V_l(r_i)
    coeffs, _, _, _ = np.linalg.lstsq(P_matrix, V_grid.T, rcond=None)
    # coeffs has shape (n_terms, num_r_points), transpose to (num_r_points, n_terms)
    coeffs = coeffs.T

    # Build a 1D cubic spline for each Legendre coefficient V_l(r)
    coeff_splines = []
    for l in range(n_terms):
        cs = CubicSpline(r_common, coeffs[:, l])
        coeff_splines.append(cs)

    return coeff_splines, n_terms, r_min_global, r_max_global, long_range_pow, radial


def _make_legendre_potential(build_data):
    coeff_splines, n_terms, r_min, r_max, long_range_pow, radial = build_data

    def potential(r, phi, deriv_order=0):
        cos_phi = np.cos(np.radians(phi))
        result = 0.0
        for l in range(n_terms):
            Pl = legendre(l)
            Pl_val = float(Pl(cos_phi))
            if deriv_order == 0:
                Vl = float(coeff_splines[l](r))
            else:
                Vl = float(coeff_splines[l].derivative(deriv_order)(r))
            result += Vl * Pl_val
        return result

    return potential


# ---------------------------------------------------------------------------
# Method 2: Radial-first interpolation (evaluate all radial splines, then
#            interpolate across angles)
# ---------------------------------------------------------------------------

def _make_radial_first_potential(angles_given, distances_given, energies_given,
                                 long_range_pow):
    radial = _build_all_radial_splines(
        angles_given, distances_given, energies_given, long_range_pow)
    angles_sorted = sorted(radial.keys())
    angles_array = np.array(angles_sorted, dtype=float)

    def potential(r, phi, deriv_order=0):
        # Evaluate the radial potential/derivative at every known angle
        vals = np.zeros(len(angles_sorted))
        for j, angle in enumerate(angles_sorted):
            info = radial[angle]
            vals[j] = _eval_radial(
                r, info['spline'], info['short_coef'], info['short_pow'],
                info['long_coef'], long_range_pow, info['r_min'], info['r_max'],
                deriv_order=deriv_order)

        # Interpolate in angle with a cubic spline (periodic-safe: data goes
        # 0 to 180, so we use a clamped spline rather than periodic)
        angle_spline = CubicSpline(angles_array, vals)
        return float(angle_spline(phi))

    return potential


# ---------------------------------------------------------------------------
# Method 3: 2D tensor-product spline on a regular grid
# ---------------------------------------------------------------------------

def _make_tensor_product_potential(angles_given, distances_given, energies_given,
                                   long_range_pow, num_r_points):
    radial = _build_all_radial_splines(
        angles_given, distances_given, energies_given, long_range_pow)
    angles_sorted = sorted(radial.keys())
    angles_array = np.array(angles_sorted, dtype=float)

    r_min_global = min(info['r_min'] for info in radial.values())
    r_max_global = max(info['r_max'] for info in radial.values())
    r_common = np.linspace(r_min_global, r_max_global, num_r_points)

    # Build V, dV/dr, d2V/dr2 grids: shape (num_r_points, num_angles)
    grids = {}
    for deriv_order in [0, 1, 2]:
        G = np.zeros((len(r_common), len(angles_sorted)))
        for j, angle in enumerate(angles_sorted):
            info = radial[angle]
            for i, r in enumerate(r_common):
                G[i, j] = _eval_radial(
                    r, info['spline'], info['short_coef'], info['short_pow'],
                    info['long_coef'], long_range_pow,
                    info['r_min'], info['r_max'], deriv_order=deriv_order)
        # RectBivariateSpline wants (x, y, z) where x and y are strictly increasing
        grids[deriv_order] = RectBivariateSpline(r_common, angles_array, G)

    def potential(r, phi, deriv_order=0):
        r_clamped = np.clip(r, r_min_global, r_max_global)
        phi_clamped = np.clip(phi, angles_array[0], angles_array[-1])
        return float(grids[deriv_order](r_clamped, phi_clamped))

    return potential


# ---------------------------------------------------------------------------
# Public interface
# ---------------------------------------------------------------------------

def create_2d_potential(angles_given, distances_given, energies_given,
                        long_range_pow,
                        method='radial_first',
                        num_legendre_terms=12,
                        num_r_points=200):
    """Create a single function V(r, phi, deriv_order=0) that returns the
    interpolated potential (or its 1st/2nd radial derivative) at any
    separation r and angle phi (in degrees).

    Parameters
    ----------
    angles_given : list of int/float
        The angles (in degrees) at which data is provided.
    distances_given : dict
        distances_given[angle] = list of r values.
    energies_given : dict
        energies_given[angle] = list of V(r) values.
    long_range_pow : float
        Power law exponent for long-range extrapolation (e.g. -4).
    method : str
        'legendre'       — Legendre polynomial expansion in angle (default)
        'radial_first'   — evaluate radial splines then interpolate in angle
        'tensor_product'  — 2D tensor-product spline on a regular grid
    num_legendre_terms : int
        Number of Legendre terms (only used by 'legendre' method).
    num_r_points : int
        Number of points on the common radial grid (used by 'legendre'
        and 'tensor_product' methods).

    Returns
    -------
    potential : callable
        potential(r, phi, deriv_order=0) -> float
        where deriv_order is 0, 1, or 2 (radial derivative order).
    """
    if method == 'legendre':
        build_data = _build_legendre(
            angles_given, distances_given, energies_given,
            long_range_pow, num_legendre_terms, num_r_points)
        return _make_legendre_potential(build_data)

    elif method == 'radial_first':
        return _make_radial_first_potential(
            angles_given, distances_given, energies_given, long_range_pow)

    elif method == 'tensor_product':
        return _make_tensor_product_potential(
            angles_given, distances_given, energies_given,
            long_range_pow, num_r_points)

    else:
        raise ValueError(
            f"Unknown method '{method}'. "
            "Choose from 'legendre', 'radial_first', or 'tensor_product'.")
