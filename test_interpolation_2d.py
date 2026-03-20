import numpy as np
from scipy.interpolate import CubicSpline
from interpolation_2d import (
    create_2d_potential, _build_all_radial_splines, _eval_radial
)

# ---------------------------------------------------------------------------
# Load data (same logic as notebook)
# ---------------------------------------------------------------------------
potential_files_and_angles = [
    ("PC.in.000", 0),   ("PC.in.010", 10),  ("PC.in.020", 20),  ("PC.in.030", 30),
    ("PC.in.040", 40),  ("PC.in.050", 50),  ("PC.in.060", 60),  ("PC.in.070", 70),
    ("PC.in.080", 80),  ("PC.in.090", 90),  ("PC.in.100", 100), ("PC.in.110", 110),
    ("PC.in.120", 120), ("PC.in.130", 130), ("PC.in.140", 140), ("PC.in.150", 150),
    ("PC.in.160", 160), ("PC.in.170", 170), ("PC.in.180", 180)]
long_range_pow = -4

given_potentials = {}
for filename, angle in potential_files_and_angles:
    with open(filename) as f:
        lines = [line.strip() for line in f.readlines()][1:]
        given_potentials[angle] = [
            (float(x.split()[0]), float(x.split()[1])) for x in lines]

angles_given = []
distances_given = {}
energies_given = {}
for angle, pv in given_potentials.items():
    angles_given.append(angle)
    distances_given[angle] = [d for d, _ in pv]
    energies_given[angle] = [e for _, e in pv]


def build_original_1d(angle):
    """Reproduce the notebook's original 1D spline for a given angle."""
    d = distances_given[angle]
    e = energies_given[angle]
    sp = (np.log(e[0]) - np.log(e[1])) / (np.log(d[0]) - np.log(d[1]))
    sc = e[0] / d[0] ** sp
    sd = sc * sp * d[0] ** (sp - 1)
    lc = e[-1] / d[-1] ** long_range_pow
    ld = lc * long_range_pow * d[-1] ** (long_range_pow - 1)
    spline = CubicSpline(d, e, bc_type=((1, sd), (1, ld)))
    return spline, sc, sp, lc, min(d), max(d)


def eval_original_1d(r, angle, deriv_order=0):
    """Evaluate the original notebook potential at (r, angle)."""
    spline, sc, sp, lc, r_min, r_max = build_original_1d(angle)
    if deriv_order == 0:
        if r < r_min:
            return sc * r ** sp
        elif r > r_max:
            return lc * r ** long_range_pow
        else:
            return float(spline(r))
    elif deriv_order == 1:
        if r < r_min:
            return sc * sp * r ** (sp - 1)
        elif r > r_max:
            return lc * long_range_pow * r ** (long_range_pow - 1)
        else:
            return float(spline.derivative(1)(r))
    elif deriv_order == 2:
        if r < r_min:
            return sc * sp * (sp - 1) * r ** (sp - 2)
        elif r > r_max:
            return lc * long_range_pow * (long_range_pow - 1) * r ** (long_range_pow - 2)
        else:
            return float(spline.derivative(2)(r))


# ---------------------------------------------------------------------------
# Test 1: radial_first reproduces original 1D at all known angles
# ---------------------------------------------------------------------------
def test_radial_first_matches_original():
    """At every known angle, radial_first should exactly match the original
    1D spline for V, dV/dr, and d2V/dr2 across the mid-range."""
    V = create_2d_potential(angles_given, distances_given, energies_given,
                            long_range_pow, method='radial_first')

    test_rs = [4.0, 5.0, 6.0, 8.0, 10.0, 15.0, 25.0]
    max_err = {0: 0.0, 1: 0.0, 2: 0.0}

    for angle in angles_given:
        for r in test_rs:
            for d_order in [0, 1, 2]:
                val_new = V(r, angle, deriv_order=d_order)
                val_orig = eval_original_1d(r, angle, deriv_order=d_order)
                err = abs(val_new - val_orig)
                max_err[d_order] = max(max_err[d_order], err)
                if err > 1e-12:
                    print(f"  MISMATCH: angle={angle}, r={r}, deriv={d_order}: "
                          f"new={val_new:.15e}, orig={val_orig:.15e}, err={err:.2e}")

    for d_order in [0, 1, 2]:
        status = "PASS" if max_err[d_order] < 1e-12 else "FAIL"
        print(f"  [{status}] deriv_order={d_order}: max error = {max_err[d_order]:.2e}")


# ---------------------------------------------------------------------------
# Test 2: short-range extrapolation (r < r_min)
# ---------------------------------------------------------------------------
def test_short_range_extrapolation():
    """All methods should produce reasonable values for r below the data range.
    Compare against original 1D at known angles."""
    test_rs = [1.0, 2.0, 3.0, 3.4]

    for method in ['legendre', 'radial_first', 'tensor_product']:
        V = create_2d_potential(angles_given, distances_given, energies_given,
                                long_range_pow, method=method, num_r_points=500)
        max_rel_err = 0.0
        worst_case = None

        for angle in angles_given:
            for r in test_rs:
                val_new = V(r, angle, deriv_order=0)
                val_orig = eval_original_1d(r, angle, deriv_order=0)
                if abs(val_orig) > 1e-15:
                    rel_err = abs(val_new - val_orig) / abs(val_orig)
                    if rel_err > max_rel_err:
                        max_rel_err = rel_err
                        worst_case = (angle, r, val_new, val_orig)

        status = "PASS" if max_rel_err < 0.05 else "FAIL"
        detail = ""
        if worst_case:
            a, r, vn, vo = worst_case
            detail = f" (worst: angle={a}, r={r}, new={vn:.6e}, orig={vo:.6e})"
        print(f"  [{status}] {method:16s}: max relative error = {max_rel_err:.4e}{detail}")


# ---------------------------------------------------------------------------
# Test 3: long-range extrapolation (r > r_max)
# ---------------------------------------------------------------------------
def test_long_range_extrapolation():
    """All methods should produce reasonable values for r above the data range.
    Compare against original 1D at known angles."""
    test_rs = [35.0, 50.0, 100.0]

    for method in ['legendre', 'radial_first', 'tensor_product']:
        V = create_2d_potential(angles_given, distances_given, energies_given,
                                long_range_pow, method=method, num_r_points=500)
        max_rel_err = 0.0
        worst_case = None

        for angle in angles_given:
            for r in test_rs:
                val_new = V(r, angle, deriv_order=0)
                val_orig = eval_original_1d(r, angle, deriv_order=0)
                if abs(val_orig) > 1e-15:
                    rel_err = abs(val_new - val_orig) / abs(val_orig)
                    if rel_err > max_rel_err:
                        max_rel_err = rel_err
                        worst_case = (angle, r, val_new, val_orig)

        status = "PASS" if max_rel_err < 0.05 else "FAIL"
        detail = ""
        if worst_case:
            a, r, vn, vo = worst_case
            detail = f" (worst: angle={a}, r={r}, new={vn:.6e}, orig={vo:.6e})"
        print(f"  [{status}] {method:16s}: max relative error = {max_rel_err:.4e}{detail}")


# ---------------------------------------------------------------------------
# Test 4: all methods agree at known angles across full r range
# ---------------------------------------------------------------------------
def test_methods_agree_at_known_angles():
    """At known angles all three methods should produce similar values."""
    methods = {}
    for m in ['legendre', 'radial_first', 'tensor_product']:
        methods[m] = create_2d_potential(
            angles_given, distances_given, energies_given,
            long_range_pow, method=m, num_r_points=500)

    test_rs = [4.0, 5.0, 6.5, 8.0, 10.0, 15.0, 25.0]
    max_rel_diff = 0.0
    worst_case = None

    for angle in angles_given:
        for r in test_rs:
            vals = {m: methods[m](r, angle) for m in methods}
            ref = vals['radial_first']
            if abs(ref) > 1e-15:
                for m in ['legendre', 'tensor_product']:
                    rel_diff = abs(vals[m] - ref) / abs(ref)
                    if rel_diff > max_rel_diff:
                        max_rel_diff = rel_diff
                        worst_case = (m, angle, r, vals[m], ref)

    status = "PASS" if max_rel_diff < 0.01 else "FAIL"
    detail = ""
    if worst_case:
        m, a, r, vn, vr = worst_case
        detail = f" (worst: {m}, angle={a}, r={r}, val={vn:.6e}, ref={vr:.6e})"
    print(f"  [{status}] max relative difference = {max_rel_diff:.4e}{detail}")


# ---------------------------------------------------------------------------
# Test 5: angular interpolation is monotone/smooth between known angles
# ---------------------------------------------------------------------------
def test_angular_smoothness():
    """At a fixed r, V(r, phi) should vary smoothly with phi. Check that
    there are no wild oscillations by verifying the interpolated values at
    5-degree steps stay between (or near) neighboring known-angle values."""
    r_test = 6.0
    n_violations = 0

    for method in ['legendre', 'radial_first', 'tensor_product']:
        V = create_2d_potential(angles_given, distances_given, energies_given,
                                long_range_pow, method=method, num_r_points=500)

        # Evaluate at known angles
        known_vals = {a: V(r_test, a) for a in angles_given}
        angles_sorted = sorted(angles_given)

        for i in range(len(angles_sorted) - 1):
            a1 = angles_sorted[i]
            a2 = angles_sorted[i + 1]
            v1 = known_vals[a1]
            v2 = known_vals[a2]
            lo = min(v1, v2)
            hi = max(v1, v2)
            margin = abs(hi - lo) * 0.5  # allow 50% overshoot

            mid_angle = (a1 + a2) / 2.0
            v_mid = V(r_test, mid_angle)

            if v_mid < lo - margin or v_mid > hi + margin:
                n_violations += 1
                print(f"  OSCILLATION: {method}, phi={mid_angle}, "
                      f"V={v_mid:.6e}, bracket=[{lo:.6e}, {hi:.6e}]")

    status = "PASS" if n_violations == 0 else "FAIL"
    print(f"  [{status}] angular oscillation violations: {n_violations}")


# ---------------------------------------------------------------------------
# Test 6: Legendre convergence — error should decrease with more terms
# ---------------------------------------------------------------------------
def test_legendre_convergence():
    """Increasing num_legendre_terms should reduce error at known angles."""
    test_rs = [5.0, 8.0, 15.0]
    V_ref = create_2d_potential(angles_given, distances_given, energies_given,
                                long_range_pow, method='radial_first')

    prev_err = float('inf')
    monotonic = True

    for n_terms in [4, 8, 12, 16]:
        V = create_2d_potential(angles_given, distances_given, energies_given,
                                long_range_pow, method='legendre',
                                num_legendre_terms=n_terms, num_r_points=500)
        max_rel_err = 0.0
        for angle in angles_given:
            for r in test_rs:
                ref = V_ref(r, angle)
                val = V(r, angle)
                if abs(ref) > 1e-15:
                    max_rel_err = max(max_rel_err, abs(val - ref) / abs(ref))

        improved = max_rel_err <= prev_err
        if not improved:
            monotonic = False
        print(f"  n_terms={n_terms:2d}: max relative error = {max_rel_err:.4e}"
              f"  {'<=' if improved else '> '} previous")
        prev_err = max_rel_err

    status = "PASS" if monotonic else "FAIL"
    print(f"  [{status}] error decreases monotonically with more terms")


# ---------------------------------------------------------------------------
# Test 7: derivatives are consistent with finite differences
# ---------------------------------------------------------------------------
def test_derivatives_finite_difference():
    """Check that the returned derivatives are consistent with finite
    differences of the returned potential."""
    h = 1e-5

    for method in ['legendre', 'radial_first', 'tensor_product']:
        V = create_2d_potential(angles_given, distances_given, energies_given,
                                long_range_pow, method=method, num_r_points=500)
        max_rel_err_d1 = 0.0
        max_rel_err_d2 = 0.0

        for angle in [0, 45, 90, 135, 180]:
            for r in [5.0, 8.0, 12.0, 20.0]:
                # First derivative: central difference
                fd1 = (V(r + h, angle) - V(r - h, angle)) / (2 * h)
                d1 = V(r, angle, deriv_order=1)
                if abs(fd1) > 1e-15:
                    max_rel_err_d1 = max(max_rel_err_d1,
                                         abs(d1 - fd1) / abs(fd1))

                # Second derivative: central difference
                fd2 = (V(r + h, angle) - 2 * V(r, angle) + V(r - h, angle)) / h**2
                d2 = V(r, angle, deriv_order=2)
                if abs(fd2) > 1e-15:
                    max_rel_err_d2 = max(max_rel_err_d2,
                                         abs(d2 - fd2) / abs(fd2))

        status_d1 = "PASS" if max_rel_err_d1 < 1e-4 else "FAIL"
        status_d2 = "PASS" if max_rel_err_d2 < 1e-3 else "FAIL"
        print(f"  [{status_d1}] {method:16s} d1: max relative error vs finite diff = {max_rel_err_d1:.4e}")
        print(f"  [{status_d2}] {method:16s} d2: max relative error vs finite diff = {max_rel_err_d2:.4e}")


# ---------------------------------------------------------------------------
# Run all tests
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    print("Test 1: radial_first matches original 1D at all known angles")
    test_radial_first_matches_original()
    print()

    print("Test 2: short-range extrapolation (r < r_min)")
    test_short_range_extrapolation()
    print()

    print("Test 3: long-range extrapolation (r > r_max)")
    test_long_range_extrapolation()
    print()

    print("Test 4: all methods agree at known angles")
    test_methods_agree_at_known_angles()
    print()

    print("Test 5: angular smoothness (no wild oscillations)")
    test_angular_smoothness()
    print()

    print("Test 6: Legendre convergence")
    test_legendre_convergence()
    print()

    print("Test 7: derivatives consistent with finite differences")
    test_derivatives_finite_difference()
    print()
