import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import minimize_scalar, root_scalar


# ---------------------------------------------------------------------------
# Input data
# ---------------------------------------------------------------------------

potential_files_and_angles = [
    ("PC.in.000", 0),   ("PC.in.010", 10),  ("PC.in.020", 20),  ("PC.in.030", 30),
    ("PC.in.040", 40),  ("PC.in.050", 50),  ("PC.in.060", 60),  ("PC.in.070", 70),
    ("PC.in.080", 80),  ("PC.in.090", 90),  ("PC.in.100", 100), ("PC.in.110", 110),
    ("PC.in.120", 120), ("PC.in.130", 130), ("PC.in.140", 140), ("PC.in.150", 150),
    ("PC.in.160", 160), ("PC.in.170", 170), ("PC.in.180", 180)]

long_range_pow = -4
minimum_energy = 1e-9
maximum_energy = 1.0


# ---------------------------------------------------------------------------
# Read potential data from files
# ---------------------------------------------------------------------------

def read_potentials(files_and_angles):
    given_potentials = {}
    for filename, angle in files_and_angles:
        with open(filename, "r") as f:
            lines = [line.strip() for line in f.readlines()]
            lines = lines[1:]
            given_potentials[angle] = []
            for line in lines:
                distance, energy = [float(x) for x in line.split()]
                given_potentials[angle].append((distance, energy))

    angles_given = []
    distances_given = {}
    energies_given = {}
    for angle, potential_values in given_potentials.items():
        angles_given.append(angle)
        distances_given[angle] = [distance for distance, _ in potential_values]
        energies_given[angle] = [energy for _, energy in potential_values]

    return angles_given, distances_given, energies_given


def read_single_potential(filename):
    """Read a single potential file in PC.in format (Fortran PC.f95 style).

    Format:
      Line 1: comment
      Line 2: accuracy, emin, emax
      Line 3: NV (number of data points)
      Line 4: NLONG
      Lines 5..5+NV-1: R  V(R)
    """
    with open(filename, "r") as f:
        comment = f.readline().strip()
        accuracy, emin, emax = [float(x) for x in f.readline().split()]
        nv = int(f.readline().strip())
        nlong = int(f.readline().strip())
        distances = []
        energies = []
        for _ in range(nv):
            parts = f.readline().split()
            distances.append(float(parts[0]))
            energies.append(float(parts[1]))
    return comment, accuracy, emin, emax, nlong, distances, energies


# ---------------------------------------------------------------------------
# Short-range extrapolation: fit first two points to V(r) = C * r^k
# ---------------------------------------------------------------------------

def get_short_range_power(distances, energies):
    delta_log_potential = np.log(energies[0]) - np.log(energies[1])
    delta_log_distance = np.log(distances[0]) - np.log(distances[1])
    return delta_log_potential / delta_log_distance


def get_short_range_coef(distances, energies):
    power = get_short_range_power(distances, energies)
    return energies[0] / distances[0]**power


# ---------------------------------------------------------------------------
# Extrapolated potential, first derivative, and second derivative
# ---------------------------------------------------------------------------

def _build_spline(distances, energies, long_range_pow):
    """Build the clamped cubic spline used by all three potential functions."""
    short_coef = get_short_range_coef(distances, energies)
    short_pow = get_short_range_power(distances, energies)
    short_deriv = short_coef * short_pow * distances[0]**(short_pow - 1)
    long_coef = energies[-1] / distances[-1]**long_range_pow
    long_deriv = long_coef * long_range_pow * distances[-1]**(long_range_pow - 1)
    spline = CubicSpline(
        distances, energies, bc_type=((1, short_deriv), (1, long_deriv)))
    return spline, short_coef, short_pow, long_coef


def create_extrapolated_potential(distances, energies, long_range_pow):
    spline, short_coef, short_pow, long_coef = _build_spline(
        distances, energies, long_range_pow)
    rmin = min(distances)
    rmax = max(distances)

    def vf(r):
        if r < rmin:
            return short_coef * r ** short_pow
        elif r > rmax:
            return long_coef * r ** long_range_pow
        else:
            return float(spline(r))
    return vf


def create_derivative_potential(distances, energies, long_range_pow):
    spline, short_coef, short_pow, long_coef = _build_spline(
        distances, energies, long_range_pow)
    deriv_spline = spline.derivative(1)
    rmin = min(distances)
    rmax = max(distances)

    def vd(r):
        if r < rmin:
            return short_coef * short_pow * r ** (short_pow - 1)
        elif r > rmax:
            return long_coef * long_range_pow * r ** (long_range_pow - 1)
        else:
            return float(deriv_spline(r))
    return vd


def create_second_derivative_potential(distances, energies, long_range_pow):
    spline, short_coef, short_pow, long_coef = _build_spline(
        distances, energies, long_range_pow)
    deriv2_spline = spline.derivative(2)
    rmin = min(distances)
    rmax = max(distances)

    def vd2(r):
        if r < rmin:
            return short_coef * short_pow * (short_pow - 1) * r ** (short_pow - 2)
        elif r > rmax:
            return long_coef * long_range_pow * (long_range_pow - 1) * r ** (long_range_pow - 2)
        else:
            return float(deriv2_spline(r))
    return vd2


# ---------------------------------------------------------------------------
# Orbiting detection (basic helpers)
# ---------------------------------------------------------------------------

def is_orbiting(r, vf, vd, vd2):
    E = vf(r) + r * vd(r) / 2.0
    if E <= 0.0:
        return False
    if vd(r) <= 0.0:
        return False
    if vd2(r) + 3 * vd(r) / r >= 0.0:
        return False
    return True


def orbiting_params(r, vf, vd):
    E = vf(r) + r * vd(r) / 2.0
    b = np.sqrt(r**3 * vd(r) / (E * 2))
    return (E, b, r)


# ---------------------------------------------------------------------------
# Step 1: Orbiting scan
# ---------------------------------------------------------------------------

def orbiting_scan(vf, vd, vd2, distances, emin, nlong, clong, n_samples=300):
    """Find orbiting parameters by root-finding the orbiting region boundaries
    and sampling on a fixed log-spaced grid.

    Returns a list of (E, b, R) tuples ordered by increasing energy, and ED.
    """
    # The orbiting conditions are:
    #   C1: V'(r) > 0
    #   C2: V''(r) + 3*V'(r)/r < 0
    #   C3: Veff = V(r) + r*V'(r)/2 > 0
    # The orbiting energy at separation r is Eorb(r) = V(r) + r*V'(r)/2.

    # Step 1: Find the inner boundary r_min where C2 = 0.
    # This is where orbiting first becomes possible (going outward).
    # Eorb(r_min) = EC (maximum orbiting energy).
    def c2_func(r):
        return vd2(r) + 3 * vd(r) / r

    r_scan = np.logspace(np.log10(distances[0]), 3, 5000)
    r_min = None
    for i in range(len(r_scan) - 1):
        if c2_func(r_scan[i]) >= 0 and c2_func(r_scan[i + 1]) < 0:
            r_min = root_scalar(c2_func, bracket=(r_scan[i], r_scan[i + 1]),
                                method='brentq', xtol=1e-14).root
            break

    if r_min is None:
        raise RuntimeError("Could not find inner orbiting boundary")

    EC = vf(r_min) + r_min * vd(r_min) / 2

    # Step 2: Find the outer boundary r_max where Eorb(r) = emin.
    def eorb_minus_emin(r):
        return vf(r) + r * vd(r) / 2 - emin

    # Eorb decreases with r (for standard potentials), so scan outward
    r_max = None
    for i in range(len(r_scan) - 1):
        if r_scan[i] > r_min:
            ea = eorb_minus_emin(r_scan[i])
            eb = eorb_minus_emin(r_scan[i + 1])
            if ea > 0 and eb <= 0:
                r_max = root_scalar(eorb_minus_emin,
                                    bracket=(r_scan[i], r_scan[i + 1]),
                                    method='brentq', xtol=1e-10).root
                break

    if r_max is None:
        # Eorb may not reach emin within scan range; use the largest r
        # where orbiting conditions still hold
        r_max = r_scan[-1]

    # Step 3: Sample on a fixed log-spaced grid within [r_min, r_max].
    r_grid = np.logspace(np.log10(r_min + 1e-10), np.log10(r_max), n_samples)

    orbit_list = []
    for r in r_grid:
        c1 = vd(r)
        c2 = vd2(r) + 3 * c1 / r
        v = vf(r)
        veff = v + r * c1 / 2

        if c1 > 0 and c2 < 0 and veff > 0:
            b = np.sqrt(r**3 * c1 / (2 * veff))
            orbit_list.append((veff, b, r))

    # Reverse so energy increases (r increases => E decreases, so reverse)
    orbit_list.reverse()

    # Handle NLONG=3 with CLONG<0: find ED (minimum orbiting energy)
    if nlong == 3 and clong < 0:
        ed = orbit_list[0][0] if orbit_list else 0.0
    else:
        ed = 0.0

    return orbit_list, ed


# ---------------------------------------------------------------------------
# Step 2: Orbiting regions
# ---------------------------------------------------------------------------

def orbiting_regions(orbit_list):
    """Partition orbiting parameters into monotonically increasing energy regions.

    Returns:
        regions: list of (E_low, E_high, idx_low, idx_high) tuples
        EC: maximum orbiting energy
        ED: minimum orbiting energy (0 if not applicable)
    """
    if not orbit_list:
        return [], 0.0

    regions = []
    increasing = True
    start_idx = 0
    start_E = orbit_list[0][0]

    for j in range(1, len(orbit_list)):
        if increasing and orbit_list[j][0] < orbit_list[j - 1][0]:
            regions.append((start_E, orbit_list[j - 1][0], start_idx, j - 1))
            increasing = False
        elif not increasing and orbit_list[j][0] > orbit_list[j - 1][0]:
            start_idx = j - 1
            start_E = orbit_list[j - 1][0]
            increasing = True

    # Final region
    regions.append((start_E, orbit_list[-1][0], start_idx, len(orbit_list) - 1))

    EC = max(r[1] for r in regions)  # max of all E_high values
    return regions


# ---------------------------------------------------------------------------
# Step 3: Turning point finder
# ---------------------------------------------------------------------------

def find_turning_point(E, b, vf, vd, distances, nlong, clong, EC=0.0, acc1=8e-4,
                       r_guess=None):
    """Find the turning point R where E - V(R) - E*b^2/R^2 = 0.

    Uses scipy brentq with robust bracket finding. r_guess provides an
    initial estimate for the search region.
    """
    def Y(r):
        return E - vf(r) - E * b**2 / r**2

    rmin_tab = distances[0]
    rmax_tab = distances[-1]

    # Analytical shortcut for NLONG=4 with large impact parameter
    if nlong == 4 and b > 0:
        if E * b**4 >= 4 * clong:
            ra = np.sqrt(b * b / 2 * (1 + np.sqrt(1 - 4 * clong / E / b**4)))
            if ra >= rmax_tab:
                return ra, b

    # Strategy: scan to find a bracket [r_lo, r_hi] where Y changes sign.
    # Y(r) = E - V(r) - E*b^2/r^2
    # At the turning point, Y = 0. For small r, V(r) is large positive
    # so Y < 0. For large r, V(r) -> 0 and E*b^2/r^2 -> 0, so Y -> E > 0.
    # We want the innermost root (smallest r where Y = 0).

    # Build a scan grid centered on the guess
    r0 = r_guess if r_guess is not None else rmax_tab
    r0 = max(r0, rmin_tab)

    # Scan inward from r0 to find where Y goes negative
    r_neg = None
    r_pos = None
    r = r0
    scale = 0.95
    if EC > 0 and 0.99 * EC < E < 1.01 * EC:
        scale = 0.999

    # First, find a point where Y > 0 (should be easy at large r)
    for _ in range(500):
        if Y(r) > 0:
            r_pos = r
            break
        r /= scale
        if r > 1e8:
            break
    if r_pos is None:
        # Try large r
        r_pos = 1e4
        if Y(r_pos) <= 0:
            r_pos = None

    # Then scan inward from r_pos to find where Y < 0
    if r_pos is not None:
        r = r_pos
        for _ in range(2000):
            r *= scale
            if r < 0.3 * rmin_tab:
                break
            if Y(r) < 0:
                r_neg = r
                break

    # If that didn't work, scan from the small-r end outward
    if r_neg is None or r_pos is None:
        r = 0.5 * rmin_tab
        if Y(r) < 0:
            r_neg = r
            r = rmin_tab
            for _ in range(2000):
                r *= 1.0 / scale
                if Y(r) > 0:
                    r_pos = r
                    break
                if r > 1e6:
                    break

    if r_neg is None or r_pos is None:
        raise RuntimeError(
            f"Could not bracket turning point: E={E}, b={b}, "
            f"r_neg={r_neg}, r_pos={r_pos}")

    # Ensure r_neg < r_pos and they bracket a root
    if r_neg > r_pos:
        r_neg, r_pos = r_pos, r_neg
    if Y(r_neg) > 0 and Y(r_pos) < 0:
        r_neg, r_pos = r_pos, r_neg

    # Use brentq for precise root finding
    r = root_scalar(Y, bracket=[min(r_neg, r_pos), max(r_neg, r_pos)],
                    method='brentq', xtol=1e-14, rtol=1e-14).root

    # Adjust b to compensate for small errors in R (critical for accuracy)
    bc = 1 - vf(r) / E
    if bc >= 0:
        b_adj = r * np.sqrt(bc)
    else:
        b_adj = r * np.sqrt(abs(bc))
    return r, b_adj


# ---------------------------------------------------------------------------
# Step 4: Deflection angle
# ---------------------------------------------------------------------------

def deflection_angle(E, b, vf, vd, distances, nlong, clong, EC, acc1,
                     naitk=0, BO=None, RO=None, ROMAX=0):
    """Compute the scattering angle theta at energy E, impact parameter b.

    Uses the Fortran ANGLE algorithm with GA/GB integrands.
    """
    HPI = np.pi / 2

    if b <= 0:
        raise ValueError(f"Impact parameter b={b} must be positive")

    # Determine which case applies
    ED = 0.0  # Will be passed in context
    case1 = True
    if naitk >= 2 and BO is not None:
        for i in range(naitk - 1):
            if BO[i] <= b < BO[i + 1]:
                case1 = False
                break

    if case1:
        # Section 6.3.1: non-orbiting deflection angle
        r_guess = ROMAX if ROMAX > 0 else distances[-1]
        rm, b_adj = find_turning_point(E, b, vf, vd, distances, nlong, clong, EC, acc1,
                                       r_guess=r_guess)

        # Endpoint values
        denom = 1 - rm**3 * vd(rm) / (2 * b_adj**2 * E)
        if denom <= 0:
            denom = abs(denom) + 1e-30
        ea = 1 - 1 / np.sqrt(denom)
        eb = 1 - b_adj / rm

        # GA integrand
        def ga(y):
            r = rm / np.cos(np.pi * (y + 1) / 4)
            val = 1 - (b_adj / r)**2 - vf(r) / E
            if val <= 0:
                # Taylor series near turning point
                val = (r - rm) * (2 * b_adj**2 / rm**3 - vf(rm) / E)
            if val <= 0:
                return 0.0
            return 1 - b_adj / rm * np.sin(np.pi * (y + 1) / 4) / np.sqrt(val)

        # Use adaptive quadrature
        result, _ = quad(ga, -1, 1, limit=200, epsrel=acc1)

        # Check for small-angle case
        theta = HPI * result
        return theta

    else:
        # Section 6.3.2: orbiting deflection angle
        rbar = max(RO[:naitk])
        r_guess = 0.5 * rbar
        rm, b_adj = find_turning_point(E, b, vf, vd, distances, nlong, clong, EC, acc1)
        if rm >= rbar:
            rm, rbar = rbar, rm
            rm, b_adj = find_turning_point(E, b_adj, vf, vd, distances, nlong, clong, EC, acc1)

        # Endpoint values
        ea = 1.0
        denom = 1 - rm**3 * vd(rm) / (2 * E * b_adj**2)
        if denom <= 0:
            # retry
            rm, b_adj = find_turning_point(E, b_adj, vf, vd, distances, nlong, clong, EC, acc1)
            denom = 1 - rm**3 * vd(rm) / (2 * E * b_adj**2)
            if denom <= 0:
                denom = abs(denom) + 1e-30
        eb = 1 - b_adj / rbar - np.arccos(rm / rbar) / np.sqrt(denom)

        def gb(y):
            # Part 1: R from rbar outward
            r8 = rbar / np.cos(np.pi * (y + 1) / 4)
            z3 = 1 - (b_adj / r8)**2 - vf(r8) / E
            if z3 != 0:
                z3_abs = abs(z3)
                func = 1 - b_adj / rbar * np.sin(np.pi * (y + 1) / 4) / np.sqrt(z3_abs)
            else:
                return ea

            # Part 2: R from rm to rbar
            z1 = np.arccos(rm / rbar)
            z2 = z1 * np.cos(np.pi * (y + 1) / 4)
            r7 = rm / np.cos(z2)
            z3b = 1 - (b_adj / r7)**2 - vf(r7) / E
            if z3b > 0:
                func -= b_adj / rm * z1 * np.sin(z2) * np.sin(np.pi * (1 + y) / 4) / np.sqrt(z3b)
            else:
                if abs(y) < 1e-15:
                    func = ea
                else:
                    func = eb

            return func

        result, _ = quad(gb, -1, 1, limit=200, epsrel=acc1)
        theta = HPI * result
        return theta


# ---------------------------------------------------------------------------
# Step 5: Cross-section integrand (qint) and EofB
# ---------------------------------------------------------------------------

def eofb(b, temp, E, ell_max, vf, vd, distances, nlong, clong, EC, acc1,
         naitk, BO, RO, ROMAX):
    """Compute temp * (1 - cos(theta)^L) for L=1..ell_max."""
    theta = deflection_angle(E, b, vf, vd, distances, nlong, clong, EC, acc1,
                             naitk, BO, RO, ROMAX)
    cos_theta = np.cos(theta)
    return np.array([temp * (1 - cos_theta**L) for L in range(1, ell_max + 1)])


def qint(y, E, ell_max, vf, vd, distances, nlong, clong,
         orbit_list, regions, EC, ED, acc1,
         naitk, BO, RO, ROMAX, BOMAX):
    """Cross-section integrand at integration point y in [-1, 1]."""
    HPI = np.pi / 2
    NO = len(orbit_list)
    EC2 = 10 * EC

    fun = np.zeros(ell_max)

    # Regime 1: ED <= E <= EC (orbiting present)
    if E >= ED and E <= EC and len(orbit_list) > 0:
        # Piece 1: smallest orbiting impact parameter
        if y <= -1 or y >= 1:
            return fun
        b2 = BO[0] * np.cos(np.pi * (y + 1) / 4)
        temp = (HPI * BO[0])**2 * np.sin(HPI * (y + 1))
        fun = eofb(b2, temp, E, ell_max, vf, vd, distances, nlong, clong, EC, acc1,
                    naitk, BO, RO, ROMAX)

        # Piece 2: intermediate orbiting impact parameters
        if y < 1:
            for i in range(naitk - 1):
                rbar = (RO[i] + RO[i + 1]) / 2
                x11 = 4 / np.pi * np.arcsin(min(RO[i], RO[i + 1]) / rbar) - 1
                x1 = ((1 - x11) * y + x11 + 1) / 2
                r1 = rbar * np.sin(np.pi * (x1 + 1) / 4)
                bc = 1 - vf(r1) / E
                if bc < 0:
                    continue
                b2 = r1 * np.sqrt(bc)
                temp = (HPI * rbar)**2 / 2 * (1 - x11) * np.sin(HPI * (x1 + 1)) * \
                    (1 - vf(r1) / E - r1 / 2 * vd(r1) / E)
                ans = eofb(b2, temp, E, ell_max, vf, vd, distances, nlong, clong, EC, acc1,
                           naitk, BO, RO, ROMAX)
                if RO[i] < RO[i + 1]:
                    fun += ans
                else:
                    fun -= ans

        if y > -1:
            for i in range(naitk - 1):
                rbar = (RO[i] + RO[i + 1]) / 2
                ro_max_i = max(RO[i], RO[i + 1])
                x12 = 4 / np.pi * np.arccos(rbar / ro_max_i) - 1
                x2 = ((1 + x12) * y + x12 - 1) / 2
                r2 = ro_max_i * np.cos(np.pi * (x2 + 1) / 4)
                bc = 1 - vf(r2) / E
                if bc < 0:
                    continue
                b2 = r2 * np.sqrt(bc)
                temp = (HPI * ro_max_i)**2 / 2 * (1 + x12) * np.sin(HPI * (x2 + 1)) * \
                    (1 - vf(r2) / E - r2 / 2 * vd(r2) / E)
                ans = eofb(b2, temp, E, ell_max, vf, vd, distances, nlong, clong, EC, acc1,
                           naitk, BO, RO, ROMAX)
                if RO[i] < RO[i + 1]:
                    fun += ans
                else:
                    fun -= ans

        # Piece 3: large orbiting impact parameters
        if -1 < y < 1:
            x3 = (1 + y) / 2
            sin_val = np.sin(HPI * x3)
            if sin_val > 0:
                r3 = RO[naitk - 1] / sin_val
                bc = 1 - vf(r3) / E
                if bc > 0:
                    b2 = r3 * np.sqrt(bc)
                    temp = (HPI * RO[naitk - 1])**2 * 2 * (r3 / RO[naitk - 1])**3 * \
                        np.cos(HPI * x3) * (1 - vf(r3) / E - r3 / 2 * vd(r3) / E)
                    ans = eofb(b2, temp, E, ell_max, vf, vd, distances, nlong, clong, EC, acc1,
                               naitk, BO, RO, ROMAX)
                    fun += ans

    # Regime 2: EC < E <= EC2 (above orbiting)
    elif E > ED and E <= EC2:
        # Find R1: root of E - V(R) - E*ROMAX^2/R^2 = 0, initial guess = ROMAX
        r1_guess = ROMAX if ROMAX > 0 else distances[-1]
        r1, _ = find_turning_point(E, ROMAX, vf, vd, distances, nlong, clong, EC, acc1,
                                   r_guess=r1_guess)
        # Find R2: root of E - V(R) = 0 (b=0), initial guess = R1
        r2, _ = find_turning_point(E, 0.0, vf, vd, distances, nlong, clong, EC, acc1,
                                   r_guess=r1)

        if y <= -1:
            fun = np.zeros(ell_max)
            for L in range(2, ell_max + 1, 2):
                fun[L - 1] = -np.pi * r2**2 * (r1 - r2) * vd(r2) / E
        elif y >= 1:
            temp = np.pi * r1 * (r1 - r2) * (1 - vf(r1) / E - r1 / 2 * vd(r1) / E)
            fun = eofb(ROMAX, temp, E, ell_max, vf, vd, distances, nlong, clong, EC, acc1,
                       naitk, BO, RO, ROMAX)
        else:
            # Piece 1: linear mapping
            r4 = ((r1 - r2) * y + r1 + r2) / 2
            bc = 1 - vf(r4) / E
            if bc < 0:
                bc = abs(bc)
            b2 = r4 * np.sqrt(bc)
            temp = np.pi * (r1 - r2) * r4 * (1 - vf(r4) / E - r4 / 2 * vd(r4) / E)
            fun = eofb(b2, temp, E, ell_max, vf, vd, distances, nlong, clong, EC, acc1,
                       naitk, BO, RO, ROMAX)

            # Piece 2: arcsine mapping
            x3 = (1 + y) / 2
            r5 = r1 / np.sin(HPI * x3) if np.sin(HPI * x3) > 0 else 1e10
            bc = 1 - vf(r5) / E
            if bc > 0:
                b2 = r5 * np.sqrt(bc)
                temp = (HPI * r1)**2 * 2 * (r5 / r1)**3 * np.cos(HPI * x3) * \
                    (1 - vf(r5) / E - r5 / 2 * vd(r5) / E)
                ans = eofb(b2, temp, E, ell_max, vf, vd, distances, nlong, clong, EC, acc1,
                           naitk, BO, RO, ROMAX)
                fun += ans

    # Regime 3: E < ED or E > EC2 (far from orbiting)
    else:
        if NO == 0:
            return fun
        if E < ED:
            b_ref = orbit_list[0][1]
        else:
            b_ref = orbit_list[-1][1]

        if y <= -1:
            fun = np.zeros(ell_max)
        elif y >= 1:
            temp = 2 * np.pi * b_ref**2
            fun = eofb(b_ref, temp, E, ell_max, vf, vd, distances, nlong, clong, EC, acc1,
                       naitk, BO, RO, ROMAX)
        else:
            # Piece 1: linear mapping
            b2 = b_ref / 2 * (y + 1)
            temp = HPI * b_ref**2 * (y + 1)
            fun = eofb(b2, temp, E, ell_max, vf, vd, distances, nlong, clong, EC, acc1,
                       naitk, BO, RO, ROMAX)

            # Piece 2: inverse mapping
            b2 = b_ref * 2 / (y + 1)
            temp = np.pi * b_ref**2 * (2 / (y + 1))**3
            ans = eofb(b2, temp, E, ell_max, vf, vd, distances, nlong, clong, EC, acc1,
                       naitk, BO, RO, ROMAX)
            fun += ans

    return fun


# ---------------------------------------------------------------------------
# Step 6: Orbiting interpolation at given energy
# ---------------------------------------------------------------------------

def interpolate_orbiting(E, orbit_list, regions, vf):
    """Interpolate orbiting parameters at energy E.

    Returns BO, RO arrays (sorted by increasing BO) and naitk.
    """
    BO_list = []
    RO_list = []

    for e_low, e_high, idx_low, idx_high in regions:
        if E < e_low or E > e_high:
            BO_list.append(0.0)
            RO_list.append(0.0)
        elif E <= e_low:
            BO_list.append(orbit_list[idx_low][1])
            RO_list.append(orbit_list[idx_low][2])
        elif E >= e_high:
            BO_list.append(orbit_list[idx_high][1])
            RO_list.append(orbit_list[idx_high][2])
        else:
            # Interpolation within this region
            n_pts = idx_high - idx_low + 1
            exx = np.array([orbit_list[idx_low + k][0] for k in range(n_pts)])
            bxx = np.array([orbit_list[idx_low + k][1] for k in range(n_pts)])
            rxx = np.array([orbit_list[idx_low + k][2] for k in range(n_pts)])

            if n_pts >= 4:
                b_spline = CubicSpline(exx, bxx, bc_type='natural')
                r_spline = CubicSpline(exx, rxx, bc_type='natural')
                BO_list.append(float(b_spline(E)))
                RO_list.append(float(r_spline(E)))
            elif n_pts >= 2:
                BO_list.append(float(np.interp(E, exx, bxx)))
                RO_list.append(float(np.interp(E, exx, rxx)))
            else:
                BO_list.append(bxx[0])
                RO_list.append(rxx[0])
                continue

    # Eliminate invalid entries, adjust BO from RO
    BO_valid = []
    RO_valid = []
    for i in range(len(RO_list)):
        if RO_list[i] > 0:
            bc = 1 - vf(RO_list[i]) / E
            if bc >= 0:
                BO_valid.append(RO_list[i] * np.sqrt(bc))
            else:
                BO_valid.append(BO_list[i])
            RO_valid.append(RO_list[i])

    # Sort by increasing BO
    if BO_valid:
        pairs = sorted(zip(BO_valid, RO_valid), key=lambda x: x[0])
        BO_sorted = [p[0] for p in pairs]
        RO_sorted = [p[1] for p in pairs]
    else:
        BO_sorted = []
        RO_sorted = []

    return BO_sorted, RO_sorted, len(BO_sorted)


def find_orbiting_at_energy(E, vf, vd, vd2, distances):
    """Directly find the orbiting separation and impact parameter at energy E
    by solving V(r) + r*V'(r)/2 = E for r.

    This avoids interpolation errors from a pre-computed table.
    Returns (BO, RO, naitk) or empty if no orbiting at this energy.
    """
    # The orbiting energy as a function of r is: Eorb(r) = V(r) + r*V'(r)/2
    # We need to find r where Eorb(r) = E, subject to the orbiting conditions.

    # Scan to find all r values where orbiting conditions hold and Eorb crosses E
    rmin_tab = distances[0]
    rmax_tab = distances[-1]

    r = rmin_tab
    crossings = []
    prev_diff = None

    while r < 1000:
        c1 = vd(r)
        c2 = vd2(r) + 3 * c1 / r
        v = vf(r)
        veff = v + r * c1 / 2.0

        if c1 > 0 and c2 < 0 and veff > 0:
            diff = veff - E
            if prev_diff is not None and prev_diff * diff < 0:
                # Sign change — there's a root between prev_r and r
                crossings.append((prev_r, r))
            prev_diff = diff
            prev_r = r
        else:
            prev_diff = None

        if r > rmax_tab:
            r *= 1.05
        else:
            r *= 1.01

    # Refine each crossing
    BO_list = []
    RO_list = []
    for r_lo, r_hi in crossings:
        def eorb_minus_E(r):
            return vf(r) + r * vd(r) / 2.0 - E
        try:
            r_root = root_scalar(eorb_minus_E, bracket=[r_lo, r_hi],
                                 method='brentq', xtol=1e-12).root
            # Verify orbiting conditions at the root
            c1 = vd(r_root)
            c2 = vd2(r_root) + 3 * c1 / r_root
            veff = vf(r_root) + r_root * c1 / 2.0
            if c1 > 0 and c2 < 0 and veff > 0:
                bc = 1 - vf(r_root) / E
                if bc >= 0:
                    BO_list.append(r_root * np.sqrt(bc))
                else:
                    BO_list.append(np.sqrt(r_root**3 * c1 / (2 * veff)))
                RO_list.append(r_root)
        except (ValueError, RuntimeError):
            continue

    if BO_list:
        pairs = sorted(zip(BO_list, RO_list), key=lambda x: x[0])
        return [p[0] for p in pairs], [p[1] for p in pairs], len(pairs)
    return [], [], 0


# ---------------------------------------------------------------------------
# Step 7: Cross-section at single energy
# ---------------------------------------------------------------------------

def compute_cross_sections(E, ell_max, vf, vd, vd2, distances, nlong, clong,
                           orbit_list, regions, EC, ED, accuracy):
    """Compute transport cross-sections Q^(1)..Q^(ell_max) at energy E."""
    acc1 = 0.8 * accuracy
    EC2 = 10 * EC

    # Set up orbiting data for this energy
    naitk = 0
    BO = []
    RO = []
    ROMAX = 0.0
    BOMAX = 0.0

    if E >= ED and E <= EC and len(orbit_list) > 0:
        # Use direct solver for orbiting parameters at this energy
        BO, RO, naitk = find_orbiting_at_energy(E, vf, vd, vd2, distances)
        if not BO:
            # Fall back to interpolation
            BO, RO, naitk = interpolate_orbiting(E, orbit_list, regions, vf)
        if RO:
            ROMAX = max(RO)
            BOMAX = max(BO)
    else:
        if orbit_list:
            ROMAX = orbit_list[-1][2]
            BOMAX = orbit_list[-1][1]

    # Define the integrand as a function of y only
    def integrand_component(y, L_index):
        result = qint(y, E, ell_max, vf, vd, distances, nlong, clong,
                      orbit_list, regions, EC, ED, acc1,
                      naitk, BO, RO, ROMAX, BOMAX)
        return result[L_index]

    # Integrate for each L
    cross_sections = np.zeros(ell_max)
    # More efficient: integrate all L at once by computing the full vector at each y
    # Use a vectorized approach with quad for each component
    for L_idx in range(ell_max):
        val, err = quad(integrand_component, -1, 1, args=(L_idx,),
                        limit=200, epsrel=acc1)
        cross_sections[L_idx] = val

    return cross_sections


# ---------------------------------------------------------------------------
# Step 8: Energy grid and full computation
# ---------------------------------------------------------------------------

def compute_all_cross_sections(vf, vd, vd2, distances, emin, emax, accuracy,
                               nlong, clong, ell_max=30):
    """Compute transport cross-sections across the full energy range."""
    # Step 1: Find orbiting parameters
    orbit_list, ED = orbiting_scan(vf, vd, vd2, distances, emin, nlong, clong)
    print(f"Found {len(orbit_list)} orbiting parameter sets")

    # Step 2: Determine regions
    regions = orbiting_regions(orbit_list)
    EC = max(r[1] for r in regions) if regions else 0.0
    EC2 = 10 * EC
    print(f"ED = {ED:.6e}, EC = {EC:.6e}")
    print(f"Number of orbiting regions: {len(regions)}")

    # Adjust emax
    emax = min(emax, vf(distances[0]))

    # Build energy grid using log-Chebyshev nodes in 3 regions
    all_results = []

    region_bounds = []
    if EC > emin:
        region_bounds.append((emin, min(EC, emax)))
    if EC < emax and EC2 > emin:
        region_bounds.append((max(EC, emin), min(EC2, emax)))
    if EC2 < emax:
        region_bounds.append((max(EC2, emin), emax))

    for e1, e2 in region_bounds:
        if e1 >= e2:
            continue
        log_e1 = np.log10(e1)
        log_e2 = np.log10(e2)

        # Start with a modest number of points
        nm = 5
        for idx in range(nm):
            log_e = (log_e2 + log_e1) / 2 - (log_e2 - log_e1) / 2 * \
                np.cos(idx * np.pi / (nm - 1))
            energy = 10**log_e
            print(f"Computing cross-sections at E = {energy:.6e}")
            cs = compute_cross_sections(energy, ell_max, vf, vd, vd2,
                                        distances, nlong, clong,
                                        orbit_list, regions, EC, ED, accuracy)
            all_results.append((energy, cs))

    # Sort by energy
    all_results.sort(key=lambda x: x[0])
    return all_results, orbit_list, regions, EC, ED


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import sys

    # Read input
    comment, accuracy, emin, emax, nlong, distances, energies = \
        read_single_potential("PC.in")

    long_range_pow = -nlong

    # Build potential functions
    vf = create_extrapolated_potential(distances, energies, long_range_pow)
    vd = create_derivative_potential(distances, energies, long_range_pow)
    vd2 = create_second_derivative_potential(distances, energies, long_range_pow)

    # Compute CLONG (same as Fortran: -V(R_last)*R_last^NLONG)
    clong = -energies[-1] * distances[-1]**nlong

    print(f"Comment: {comment}")
    print(f"Accuracy: {accuracy}, EMIN: {emin}, EMAX: {emax}, NLONG: {nlong}")
    print(f"CLONG: {clong}")
    print(f"Number of potential points: {len(distances)}")

    # Step 1: Orbiting scan
    orbit_list, ED = orbiting_scan(vf, vd, vd2, distances, emin, nlong, clong)
    print(f"\nFound {len(orbit_list)} sets of orbiting parameters")
    print(f"{'E':>24s} {'B':>24s} {'RM':>24s}")
    for E_orb, b_orb, r_orb in orbit_list:
        print(f"{E_orb:24.14e} {b_orb:24.14e} {r_orb:24.14e}")

    # Step 2: Regions
    regions = orbiting_regions(orbit_list)
    EC = max(r[1] for r in regions) if regions else 0.0
    print(f"\nED = {ED:.14e}")
    print(f"EC = {EC:.14e}")
    print("Orbiting energy regions:")
    for i, (e_low, e_high, idx_low, idx_high) in enumerate(regions):
        print(f"  Region {i + 1}: {e_low:.14e} - {e_high:.14e}")

    # Step 7-8: Compute cross-sections at a few test energies from PC.out
    test_energies = [1e-9, 6.48e-9, 5.9e-7, 5.37e-5, 3.48e-4]
    print("\nCross-section results:")
    for E_test in test_energies:
        if E_test > emax:
            continue
        print(f"\n  E = {E_test:.6e}")
        cs = compute_cross_sections(E_test, 30, vf, vd, vd2, distances,
                                    nlong, clong, orbit_list, regions,
                                    EC, ED, accuracy)
        # Print first 6 values (Q1 through Q6)
        for i in range(0, min(30, len(cs)), 2):
            print(f"    Q({i+1:2d})={cs[i]:12.4f}  Q({i+2:2d})={cs[i+1]:12.4f}")
