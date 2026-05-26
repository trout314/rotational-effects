"""Compute transport cross-sections from intermolecular potentials.

Reproduces the results of the Fortran program PC.f95 by L. A. Viehland.
"""

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import quad
from scipy.optimize import root_scalar


# Cross-cutting constants used by multiple integrand / quadrature call sites.
# Floor added to denominators near turning-point singularities to avoid
# division by zero in the deflection-angle integrands.
DENOM_FLOOR = 1.0e-30
# Maximum number of subintervals scipy.quad is allowed to use for any of
# the deflection / cross-section integrals.
QUAD_SUBINTERVAL_LIMIT = 200


# ---------------------------------------------------------------------------
# Potential with extrapolation
# ---------------------------------------------------------------------------

class Potential:
    """Intermolecular potential V(r) with short- and long-range extrapolation.

    Short range: V(r) = C * r^k  (power law fit to first two data points)
    Mid range:   Clamped cubic spline through tabulated data
    Long range:  V(r) = long_coef * r^long_range_pow

    Provides V(r), V'(r), and V''(r) from a single spline.
    """

    def __init__(self, distances, energies, long_range_pow):
        self.distances = np.asarray(distances, dtype=float)
        self.energies = np.asarray(energies, dtype=float)
        self.long_range_pow = long_range_pow
        self.rmin = self.distances[0]
        self.rmax = self.distances[-1]

        # Short-range power law: V = short_coef * r^short_pow
        self.short_pow = (np.log(energies[0] / energies[1])
                          / np.log(distances[0] / distances[1]))
        self.short_coef = energies[0] / distances[0]**self.short_pow

        # Long-range: V = long_coef * r^long_range_pow
        self.long_coef = energies[-1] / distances[-1]**long_range_pow

        # Clamped cubic spline with matching derivatives at boundaries
        short_deriv = self.short_coef * self.short_pow * distances[0]**(self.short_pow - 1)
        long_deriv = self.long_coef * long_range_pow * distances[-1]**(long_range_pow - 1)
        self._spline = CubicSpline(
            distances, energies, bc_type=((1, short_deriv), (1, long_deriv)))
        self._spline_d1 = self._spline.derivative(1)
        self._spline_d2 = self._spline.derivative(2)

    def __call__(self, r):
        """Evaluate V(r)."""
        if r < self.rmin:
            return self.short_coef * r ** self.short_pow
        elif r > self.rmax:
            return self.long_coef * r ** self.long_range_pow
        return float(self._spline(r))

    def deriv(self, r):
        """Evaluate V'(r)."""
        if r < self.rmin:
            return self.short_coef * self.short_pow * r ** (self.short_pow - 1)
        elif r > self.rmax:
            return self.long_coef * self.long_range_pow * r ** (self.long_range_pow - 1)
        return float(self._spline_d1(r))

    def deriv2(self, r):
        """Evaluate V''(r)."""
        if r < self.rmin:
            return (self.short_coef * self.short_pow * (self.short_pow - 1)
                    * r ** (self.short_pow - 2))
        elif r > self.rmax:
            n = self.long_range_pow
            return self.long_coef * n * (n - 1) * r ** (n - 2)
        return float(self._spline_d2(r))


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def read_single_potential(filename):
    """Read a potential file in PC.in format.

    Returns (comment, accuracy, emin, emax, nlong, distances, energies).
    """
    with open(filename) as f:
        comment = f.readline().strip()
        accuracy, emin, emax = map(float, f.readline().split())
        nv = int(f.readline())
        nlong = int(f.readline())
        data = np.loadtxt(f, max_rows=nv)
    return comment, accuracy, emin, emax, nlong, data[:, 0].tolist(), data[:, 1].tolist()


def read_potentials(files_and_angles):
    """Read angle-dependent potentials from PC.in.XXX files."""
    angles, distances, energies = [], {}, {}
    for filename, angle in files_and_angles:
        data = np.loadtxt(filename, skiprows=1)
        angles.append(angle)
        distances[angle] = data[:, 0].tolist()
        energies[angle] = data[:, 1].tolist()
    return angles, distances, energies


# ---------------------------------------------------------------------------
# Orbiting scan
# ---------------------------------------------------------------------------

def orbiting_scan(pot, emin, nlong, clong):
    """Find orbiting parameters by outward walk + refinement + densification.

    Port of SUBROUTINE ORBITS in PC.f95 by L. A. Viehland.  Returns
    (orbit_list, ED) where orbit_list is [(E, b, R), ...] in (mostly)
    ascending energy order, with localized E dips marking multi-orbit
    regions that orbiting_regions then splits on.
    """
    # Hard ceiling on r during the outward search (matches PC.f95 line 401).
    WALK_R_LIMIT = 1.0e3
    # Multiplicative step taken while r is still within the tabulated data.
    WALK_STEP_NEAR = 1.01
    # Larger step once r has exited the tabulated data range.
    WALK_STEP_FAR = 1.05
    # Base of the geometric nudge factor (1 - REFINE_FACTOR_BASE**level) used
    # to polish each orbiting separation; level runs 1..REFINE_LEVELS.
    REFINE_FACTOR_BASE = 0.001
    REFINE_LEVELS = 7
    # Number of intermediate r samples inserted between two consecutive orbits
    # when the orbit-list energy goes "the wrong way" (multi-orbit signature).
    DENSIFY_STEPS = 19

    def eval_at(r):
        c1 = pot.deriv(r)
        c2 = pot.deriv2(r) + 3 * c1 / r
        veff = pot(r) + r * c1 / 2
        return c1, c2, veff

    def make_entry(r, c1, veff):
        return (veff, (r**3 * c1 / (2 * veff))**0.5, r)

    # Phase 1: outward walk from pot.rmin with inward refinement at each orbit.
    # orbits is filled in increasing-R order (innermost first, outermost last).
    orbits = []
    r = pot.rmin
    while True:
        c1, c2, veff = eval_at(r)
        if c1 > 0 and c2 < 0 and veff > 0:
            orbits.append(make_entry(r, c1, veff))
            is_first = (len(orbits) == 1)
            for level in range(1, REFINE_LEVELS + 1):
                while True:
                    r_try = (1 - REFINE_FACTOR_BASE**level) * orbits[-1][2]
                    c1_t, c2_t, veff_t = eval_at(r_try)
                    if (c1_t > 0 and c2_t < 0
                            and veff_t > orbits[-1][0]):
                        orbits[-1] = make_entry(r_try, c1_t, veff_t)
                        if is_first:
                            continue
                    break
        if orbits and 0 < orbits[-1][0] <= emin:
            break
        r *= WALK_STEP_FAR if r > pot.rmax else WALK_STEP_NEAR
        if r > WALK_R_LIMIT:
            if not orbits:
                raise RuntimeError(
                    f"No orbiting found up to r={WALK_R_LIMIT}.")
            break

    if not orbits:
        return [], 0.0

    # Phase 2: special ED handling for NLONG=3 with attractive long range.
    # Refines the outermost orbit by nudging OUTWARD, looking for a local
    # veff MINIMUM (the lower edge of the orbiting energy band).
    ed = 0.0
    nlong3_special = False
    if nlong == 3 and clong < 0:
        outer = orbits[-1]
        for level in range(1, REFINE_LEVELS + 1):
            r_try = (1 + REFINE_FACTOR_BASE**level) * outer[2]
            c1_t, c2_t, veff_t = eval_at(r_try)
            if (c1_t > 0 and c2_t < 0
                    and 0 < veff_t < outer[0]):
                outer = make_entry(r_try, c1_t, veff_t)
        orbits[-1] = outer
        ed = outer[0]
        nlong3_special = True

    # Phase 3: walk outermost->innermost, densifying multi-orbit regions.
    # walk_order[0] is outermost (lowest E), walk_order[-1] is innermost.
    walk_order = list(reversed(orbits))
    new_orbits = []
    start_k = 0
    if nlong3_special:
        new_orbits.append(walk_order[0])
        start_k = 1

    n = len(walk_order)
    for k in range(start_k, n):
        new_orbits.append(walk_order[k])
        if k == n - 1:
            break
        if k < n - 2:
            # Standard case: interpolate between walk_order[k] and walk_order[k+1]
            # whenever E is about to dip (multi-orbit signature).
            if new_orbits[-1][0] > walk_order[k + 1][0]:
                r_K = new_orbits[-1][2]
                r_NO = walk_order[k + 1][2]
                divisions = DENSIFY_STEPS + 1
                for i in range(1, divisions):
                    r_test = ((divisions - i) * r_K + i * r_NO) / divisions
                    c1_t, c2_t, veff_t = eval_at(r_test)
                    if (c1_t > 0 and c2_t < 0
                            and veff_t > new_orbits[-1][0]):
                        new_orbits.append(make_entry(r_test, c1_t, veff_t))
        else:
            # k == n - 2: exactly one entry remains in walk_order.
            # Extrapolate inward from the previous step (PC.f95:484-499).
            if (len(new_orbits) >= 2
                    and new_orbits[-1][0] < new_orbits[-2][0]):
                r_J = new_orbits[-1][2]
                r_Jm1 = new_orbits[-2][2]
                divisions = DENSIFY_STEPS + 1
                for i in range(1, divisions):
                    r_test = r_J + i / divisions * (r_J - r_Jm1)
                    c1_t, c2_t, veff_t = eval_at(r_test)
                    if (c1_t > 0 and c2_t < 0
                            and veff_t > new_orbits[-1][0]):
                        new_orbits.append(make_entry(r_test, c1_t, veff_t))

    # Phase 4: trim trailing entries that haven't started ascending again.
    while len(new_orbits) >= 2 and new_orbits[-1][0] < new_orbits[-2][0]:
        new_orbits.pop()

    return new_orbits, ed


def orbiting_regions(orbit_list):
    """Partition orbiting parameters into monotonically increasing energy regions.

    Returns list of (E_low, E_high, idx_low, idx_high) tuples.
    """
    if not orbit_list:
        return []

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

    regions.append((start_E, orbit_list[-1][0], start_idx, len(orbit_list) - 1))
    return regions


# ---------------------------------------------------------------------------
# Cross-section solver
# ---------------------------------------------------------------------------

class CrossSectionSolver:
    """Computes transport cross-sections Q^(1)..Q^(ell_max) at given energies.

    Encapsulates the potential, orbiting data, and numerical parameters.
    """

    def __init__(self, pot, orbit_list, regions, nlong, clong, accuracy):
        # Integrand-level accuracy is held this much tighter than the requested
        # overall cross-section accuracy.
        INTEGRAND_ACCURACY_FRACTION = 0.8
        # Cross-section energy bands split at EMIN, EC, and EC2 = EC2_FACTOR * EC.
        EC2_FACTOR = 10

        self.pot = pot
        self.orbit_list = orbit_list
        self.regions = regions
        self.nlong = nlong
        self.clong = clong
        self.accuracy = accuracy
        self.acc1 = INTEGRAND_ACCURACY_FRACTION * accuracy
        self.EC = max(r[1] for r in regions) if regions else 0.0
        self.ED = 0.0  # Set externally for nlong=3 case
        self.EC2 = EC2_FACTOR * self.EC

    # -- Turning point finder ------------------------------------------------

    def find_turning_point(self, E, b, r_guess=None):
        """Find R where E - V(R) - E*b²/R² = 0."""
        # Half-width of the "near EC" window (1% on either side of EC) where
        # the bracket-shrink factor must be made very small.
        EC_PROXIMITY_FRACTION = 0.01
        # Bracket-shrink factor inside the near-EC window (slow & careful).
        BRACKET_SHRINK_NEAR_EC = 0.999
        # Default bracket-shrink factor away from EC.
        BRACKET_SHRINK_DEFAULT = 0.95
        # Iteration cap on every bracket-search loop below.
        BRACKET_SEARCH_ITER_MAX = 2000
        # Upper-r fallback used when the outward bracket couldn't find Y > 0.
        BRACKET_UPPER_FALLBACK_R = 1.0e4
        # Abort the search if r ever exceeds this during outward bracketing.
        BRACKET_ABORT_R = 1.0e6
        # Stop the inward bracket search once r drops below this fraction of rmin.
        BRACKET_R_NEG_RMIN_FRACTION = 0.3
        # Initial r used in the bracket-finding fallback branch (fraction of rmin).
        BRACKET_FALLBACK_R_FRACTION = 0.5
        # brentq xtol/rtol for the final turning-point root.
        TURNING_POINT_TOL = 1.0e-14

        pot = self.pot

        def Y(r):
            return E - pot(r) - E * b**2 / r**2

        # Analytical shortcut for NLONG=4
        if self.nlong == 4 and b > 0 and E * b**4 >= 4 * self.clong:
            ra = np.sqrt(b**2 / 2 * (1 + np.sqrt(1 - 4 * self.clong / E / b**4)))
            if ra >= pot.rmax:
                return ra, b

        # Bracket finding: scan from initial guess
        r0 = r_guess if r_guess is not None else pot.rmax
        r0 = max(r0, pot.rmin)
        near_ec = (self.EC > 0
                   and (1 - EC_PROXIMITY_FRACTION) * self.EC < E
                                                            < (1 + EC_PROXIMITY_FRACTION) * self.EC)
        scale = BRACKET_SHRINK_NEAR_EC if near_ec else BRACKET_SHRINK_DEFAULT

        # Find r_pos (Y > 0) and r_neg (Y < 0)
        r_pos = r_neg = None
        r = r0
        for _ in range(BRACKET_SEARCH_ITER_MAX):
            if Y(r) > 0:
                r_pos = r
                break
            r /= scale
        if r_pos is None:
            if Y(BRACKET_UPPER_FALLBACK_R) > 0:
                r_pos = BRACKET_UPPER_FALLBACK_R

        if r_pos is not None:
            r = r_pos
            for _ in range(BRACKET_SEARCH_ITER_MAX):
                r *= scale
                if r < BRACKET_R_NEG_RMIN_FRACTION * pot.rmin:
                    break
                if Y(r) < 0:
                    r_neg = r
                    break

        if r_neg is None or r_pos is None:
            r = BRACKET_FALLBACK_R_FRACTION * pot.rmin
            if Y(r) < 0:
                r_neg = r
                r = pot.rmin
                for _ in range(BRACKET_SEARCH_ITER_MAX):
                    r /= scale
                    if Y(r) > 0:
                        r_pos = r
                        break
                    if r > BRACKET_ABORT_R:
                        break

        if r_neg is None or r_pos is None:
            raise RuntimeError(f"Could not bracket turning point: E={E}, b={b}")

        lo, hi = min(r_neg, r_pos), max(r_neg, r_pos)
        r = root_scalar(Y, bracket=[lo, hi], method='brentq',
                        xtol=TURNING_POINT_TOL, rtol=TURNING_POINT_TOL).root

        # Adjust b for consistency
        bc = 1 - pot(r) / E
        b_adj = r * np.sqrt(abs(bc))
        return r, b_adj

    # -- Orbiting at specific energy -----------------------------------------

    def find_orbiting_at_energy(self, E):
        """Directly solve for orbiting (b, R) at energy E."""
        # Number of log-spaced samples used to bracket orbiting roots at fixed E.
        ORBIT_AT_ENERGY_SCAN_POINTS = 3000
        # brentq xtol for the orbit-at-energy root.
        ORBIT_AT_ENERGY_TOL = 1.0e-12

        pot = self.pot
        r_scan = np.logspace(np.log10(pot.rmin), 3, ORBIT_AT_ENERGY_SCAN_POINTS)

        # Evaluate Eorb - E across scan, tracking sign changes within orbiting region
        prev_diff = prev_r = None
        crossings = []
        for r in r_scan:
            c1 = pot.deriv(r)
            c2 = pot.deriv2(r) + 3 * c1 / r
            veff = pot(r) + r * c1 / 2
            if c1 > 0 and c2 < 0 and veff > 0:
                diff = veff - E
                if prev_diff is not None and prev_diff * diff < 0:
                    crossings.append((prev_r, r))
                prev_diff, prev_r = diff, r
            else:
                prev_diff = None

        BO, RO = [], []
        for r_lo, r_hi in crossings:
            try:
                r_root = root_scalar(
                    lambda r: pot(r) + r * pot.deriv(r) / 2 - E,
                    bracket=[r_lo, r_hi], method='brentq',
                    xtol=ORBIT_AT_ENERGY_TOL).root
                c1 = pot.deriv(r_root)
                c2 = pot.deriv2(r_root) + 3 * c1 / r_root
                veff = pot(r_root) + r_root * c1 / 2
                if c1 > 0 and c2 < 0 and veff > 0:
                    bc = 1 - pot(r_root) / E
                    BO.append(r_root * np.sqrt(abs(bc)) if bc >= 0
                              else np.sqrt(r_root**3 * c1 / (2 * veff)))
                    RO.append(r_root)
            except (ValueError, RuntimeError):
                continue

        if BO:
            pairs = sorted(zip(BO, RO))
            return [p[0] for p in pairs], [p[1] for p in pairs], len(pairs)
        return [], [], 0

    def interpolate_orbiting(self, E):
        """Interpolate orbiting parameters from pre-computed table."""
        BO_list, RO_list = [], []
        for e_low, e_high, idx_low, idx_high in self.regions:
            if E < e_low or E > e_high:
                BO_list.append(0.0); RO_list.append(0.0)
            elif E <= e_low:
                BO_list.append(self.orbit_list[idx_low][1])
                RO_list.append(self.orbit_list[idx_low][2])
            elif E >= e_high:
                BO_list.append(self.orbit_list[idx_high][1])
                RO_list.append(self.orbit_list[idx_high][2])
            else:
                n = idx_high - idx_low + 1
                exx = np.array([self.orbit_list[idx_low + k][0] for k in range(n)])
                bxx = np.array([self.orbit_list[idx_low + k][1] for k in range(n)])
                rxx = np.array([self.orbit_list[idx_low + k][2] for k in range(n)])
                if n >= 4:
                    BO_list.append(float(CubicSpline(exx, bxx, bc_type='natural')(E)))
                    RO_list.append(float(CubicSpline(exx, rxx, bc_type='natural')(E)))
                else:
                    BO_list.append(float(np.interp(E, exx, bxx)))
                    RO_list.append(float(np.interp(E, exx, rxx)))

        # Filter, adjust, sort
        BO_valid, RO_valid = [], []
        for bo, ro in zip(BO_list, RO_list):
            if ro > 0:
                bc = 1 - self.pot(ro) / E
                BO_valid.append(ro * np.sqrt(bc) if bc >= 0 else bo)
                RO_valid.append(ro)

        if BO_valid:
            pairs = sorted(zip(BO_valid, RO_valid))
            return [p[0] for p in pairs], [p[1] for p in pairs], len(pairs)
        return [], [], 0

    # -- Deflection angle ----------------------------------------------------

    def deflection_angle(self, E, b, naitk, BO, RO, ROMAX):
        """Compute scattering angle theta at energy E, impact parameter b."""
        HPI = np.pi / 2
        pot = self.pot

        # Determine case
        case1 = True
        if naitk >= 2 and BO is not None:
            for i in range(naitk - 1):
                if BO[i] <= b < BO[i + 1]:
                    case1 = False
                    break

        if case1:
            r_guess = ROMAX if ROMAX > 0 else pot.rmax
            rm, b_adj = self.find_turning_point(E, b, r_guess=r_guess)

            denom = 1 - rm**3 * pot.deriv(rm) / (2 * b_adj**2 * E)
            if denom <= 0:
                denom = abs(denom) + DENOM_FLOOR

            def ga(y):
                r = rm / np.cos(np.pi * (y + 1) / 4)
                val = 1 - (b_adj / r)**2 - pot(r) / E
                if val <= 0:
                    val = (r - rm) * (2 * b_adj**2 / rm**3 - pot(rm) / E)
                if val <= 0:
                    return 0.0
                return 1 - b_adj / rm * np.sin(np.pi * (y + 1) / 4) / np.sqrt(val)

            result, _ = quad(ga, -1, 1,
                             limit=QUAD_SUBINTERVAL_LIMIT, epsrel=self.acc1)
            return HPI * result

        else:
            # Treat |y| below this as exactly 0 when evaluating gb (machine-eps).
            Y_EPSILON = 1.0e-15

            rbar = max(RO[:naitk])
            rm, b_adj = self.find_turning_point(E, b)
            if rm >= rbar:
                rm, rbar = rbar, rm
                rm, b_adj = self.find_turning_point(E, b_adj)

            denom = 1 - rm**3 * pot.deriv(rm) / (2 * E * b_adj**2)
            if denom <= 0:
                rm, b_adj = self.find_turning_point(E, b_adj)
                denom = 1 - rm**3 * pot.deriv(rm) / (2 * E * b_adj**2)
                if denom <= 0:
                    denom = abs(denom) + DENOM_FLOOR

            ea = 1.0
            eb = 1 - b_adj / rbar - np.arccos(rm / rbar) / np.sqrt(denom)

            def gb(y):
                r8 = rbar / np.cos(np.pi * (y + 1) / 4)
                z3 = 1 - (b_adj / r8)**2 - pot(r8) / E
                if z3 == 0:
                    return ea
                func = 1 - b_adj / rbar * np.sin(np.pi * (y + 1) / 4) / np.sqrt(abs(z3))

                z1 = np.arccos(rm / rbar)
                z2 = z1 * np.cos(np.pi * (y + 1) / 4)
                r7 = rm / np.cos(z2)
                z3b = 1 - (b_adj / r7)**2 - pot(r7) / E
                if z3b > 0:
                    func -= (b_adj / rm * z1 * np.sin(z2)
                             * np.sin(np.pi * (1 + y) / 4) / np.sqrt(z3b))
                else:
                    func = ea if abs(y) < Y_EPSILON else eb
                return func

            result, _ = quad(gb, -1, 1,
                             limit=QUAD_SUBINTERVAL_LIMIT, epsrel=self.acc1)
            return HPI * result

    # -- EofB: deflection -> cross-section integrand -------------------------

    def eofb(self, b, temp, E, ell_max, naitk, BO, RO, ROMAX):
        """Return temp * (1 - cos(theta)^L) for L=1..ell_max."""
        theta = self.deflection_angle(E, b, naitk, BO, RO, ROMAX)
        ct = np.cos(theta)
        # Vectorized power computation
        ct_powers = ct ** np.arange(1, ell_max + 1)
        return temp * (1 - ct_powers)

    # -- QINT: cross-section integrand ---------------------------------------

    def qint(self, y, E, ell_max, naitk, BO, RO, ROMAX, BOMAX):
        """Cross-section integrand at point y in [-1, 1].

        Returns array of length ell_max.
        """
        HPI = np.pi / 2
        pot = self.pot
        NO = len(self.orbit_list)
        fun = np.zeros(ell_max)

        # Regime 1: orbiting (ED <= E <= EC)
        if E >= self.ED and E <= self.EC and NO > 0:
            if y <= -1 or y >= 1:
                return fun

            # Piece 1: small b (0 to BO[0])
            b2 = BO[0] * np.cos(np.pi * (y + 1) / 4)
            temp = (HPI * BO[0])**2 * np.sin(HPI * (y + 1))
            fun = self.eofb(b2, temp, E, ell_max, naitk, BO, RO, ROMAX)

            # Piece 2a: intermediate (lower half of each RO pair)
            if y < 1:
                for i in range(naitk - 1):
                    rbar = (RO[i] + RO[i + 1]) / 2
                    x11 = 4 / np.pi * np.arcsin(min(RO[i], RO[i + 1]) / rbar) - 1
                    x1 = ((1 - x11) * y + x11 + 1) / 2
                    r1 = rbar * np.sin(np.pi * (x1 + 1) / 4)
                    bc = 1 - pot(r1) / E
                    if bc < 0:
                        continue
                    b2 = r1 * np.sqrt(bc)
                    temp = ((HPI * rbar)**2 / 2 * (1 - x11)
                            * np.sin(HPI * (x1 + 1))
                            * (1 - pot(r1) / E - r1 / 2 * pot.deriv(r1) / E))
                    ans = self.eofb(b2, temp, E, ell_max, naitk, BO, RO, ROMAX)
                    fun += ans if RO[i] < RO[i + 1] else -ans

            # Piece 2b: intermediate (upper half)
            if y > -1:
                for i in range(naitk - 1):
                    rbar = (RO[i] + RO[i + 1]) / 2
                    ro_max_i = max(RO[i], RO[i + 1])
                    x12 = 4 / np.pi * np.arccos(rbar / ro_max_i) - 1
                    x2 = ((1 + x12) * y + x12 - 1) / 2
                    r2 = ro_max_i * np.cos(np.pi * (x2 + 1) / 4)
                    bc = 1 - pot(r2) / E
                    if bc < 0:
                        continue
                    b2 = r2 * np.sqrt(bc)
                    temp = ((HPI * ro_max_i)**2 / 2 * (1 + x12)
                            * np.sin(HPI * (x2 + 1))
                            * (1 - pot(r2) / E - r2 / 2 * pot.deriv(r2) / E))
                    ans = self.eofb(b2, temp, E, ell_max, naitk, BO, RO, ROMAX)
                    fun += ans if RO[i] < RO[i + 1] else -ans

            # Piece 3: large b (RO outward)
            if -1 < y < 1:
                x3 = (1 + y) / 2
                sin_val = np.sin(HPI * x3)
                if sin_val > 0:
                    r3 = RO[naitk - 1] / sin_val
                    bc = 1 - pot(r3) / E
                    if bc > 0:
                        b2 = r3 * np.sqrt(bc)
                        temp = ((HPI * RO[naitk - 1])**2 * 2
                                * (r3 / RO[naitk - 1])**3
                                * np.cos(HPI * x3)
                                * (1 - pot(r3) / E - r3 / 2 * pot.deriv(r3) / E))
                        fun += self.eofb(b2, temp, E, ell_max, naitk, BO, RO, ROMAX)

        # Regime 2: above orbiting (EC < E <= 10*EC)
        elif E > self.ED and E <= self.EC2:
            r1_guess = ROMAX if ROMAX > 0 else pot.rmax
            r1, _ = self.find_turning_point(E, ROMAX, r_guess=r1_guess)
            r2, _ = self.find_turning_point(E, 0.0, r_guess=r1)

            if y <= -1:
                for L in range(2, ell_max + 1, 2):
                    fun[L - 1] = -np.pi * r2**2 * (r1 - r2) * pot.deriv(r2) / E
            elif y >= 1:
                temp = np.pi * r1 * (r1 - r2) * (1 - pot(r1) / E - r1 / 2 * pot.deriv(r1) / E)
                fun = self.eofb(ROMAX, temp, E, ell_max, naitk, BO, RO, ROMAX)
            else:
                r4 = ((r1 - r2) * y + r1 + r2) / 2
                bc = 1 - pot(r4) / E
                b2 = r4 * np.sqrt(abs(bc))
                temp = np.pi * (r1 - r2) * r4 * (1 - pot(r4) / E - r4 / 2 * pot.deriv(r4) / E)
                fun = self.eofb(b2, temp, E, ell_max, naitk, BO, RO, ROMAX)

                x3 = (1 + y) / 2
                sin_val = np.sin(HPI * x3)
                if sin_val > 0:
                    r5 = r1 / sin_val
                    bc = 1 - pot(r5) / E
                    if bc > 0:
                        b2 = r5 * np.sqrt(bc)
                        temp = ((HPI * r1)**2 * 2 * (r5 / r1)**3
                                * np.cos(HPI * x3)
                                * (1 - pot(r5) / E - r5 / 2 * pot.deriv(r5) / E))
                        fun += self.eofb(b2, temp, E, ell_max, naitk, BO, RO, ROMAX)

        # Regime 3: far from orbiting
        else:
            if NO == 0:
                return fun
            b_ref = self.orbit_list[0][1] if E < self.ED else self.orbit_list[-1][1]

            if y <= -1:
                pass
            elif y >= 1:
                fun = self.eofb(b_ref, 2 * np.pi * b_ref**2, E, ell_max,
                                naitk, BO, RO, ROMAX)
            else:
                b2 = b_ref / 2 * (y + 1)
                temp = HPI * b_ref**2 * (y + 1)
                fun = self.eofb(b2, temp, E, ell_max, naitk, BO, RO, ROMAX)

                b2 = b_ref * 2 / (y + 1)
                temp = np.pi * b_ref**2 * (2 / (y + 1))**3
                fun += self.eofb(b2, temp, E, ell_max, naitk, BO, RO, ROMAX)

        return fun

    # -- Main entry point ----------------------------------------------------

    def compute(self, E, ell_max=30):
        """Compute transport cross-sections Q^(1)..Q^(ell_max) at energy E."""
        # Set up orbiting data
        naitk, BO, RO, ROMAX, BOMAX = 0, [], [], 0.0, 0.0

        if E >= self.ED and E <= self.EC and self.orbit_list:
            BO, RO, naitk = self.find_orbiting_at_energy(E)
            if not BO:
                BO, RO, naitk = self.interpolate_orbiting(E)
            if RO:
                ROMAX, BOMAX = max(RO), max(BO)
        elif self.orbit_list:
            ROMAX = self.orbit_list[-1][2]
            BOMAX = self.orbit_list[-1][1]

        # Integrate all ell values simultaneously
        def integrand(y):
            return self.qint(y, E, ell_max, naitk, BO, RO, ROMAX, BOMAX)

        # quad only handles scalar functions, so we integrate component by
        # component but cache the integrand evaluation to avoid redundant
        # deflection angle computations at the same y.
        cache = {}

        def cached_integrand(y, L_idx):
            if y not in cache:
                cache[y] = integrand(y)
            return cache[y][L_idx]

        cross_sections = np.zeros(ell_max)
        for L_idx in range(ell_max):
            cache.clear()
            val, _ = quad(cached_integrand, -1, 1, args=(L_idx,),
                          limit=QUAD_SUBINTERVAL_LIMIT, epsrel=self.acc1)
            cross_sections[L_idx] = val

        return cross_sections


# ---------------------------------------------------------------------------
# Full computation driver
# ---------------------------------------------------------------------------

def compute_all_cross_sections(pot, emin, emax, accuracy, nlong, clong, ell_max=30):
    """Compute transport cross-sections across the full energy range."""
    orbit_list, ED = orbiting_scan(pot, emin, nlong, clong)
    regions = orbiting_regions(orbit_list)

    solver = CrossSectionSolver(pot, orbit_list, regions, nlong, clong, accuracy)
    solver.ED = ED

    EC = solver.EC
    EC2 = solver.EC2
    emax = min(emax, pot(pot.rmin))

    # Log-Chebyshev energy grids across 3 regions
    region_bounds = []
    if EC > emin:
        region_bounds.append((emin, min(EC, emax)))
    if EC < emax and EC2 > emin:
        region_bounds.append((max(EC, emin), min(EC2, emax)))
    if EC2 < emax:
        region_bounds.append((max(EC2, emin), emax))

    all_results = []
    for e1, e2 in region_bounds:
        if e1 >= e2:
            continue
        log_e1, log_e2 = np.log10(e1), np.log10(e2)
        nm = 5
        for idx in range(nm):
            energy = 10**((log_e2 + log_e1) / 2
                          - (log_e2 - log_e1) / 2 * np.cos(idx * np.pi / (nm - 1)))
            cs = solver.compute(energy, ell_max)
            all_results.append((energy, cs))

    all_results.sort(key=lambda x: x[0])
    return all_results, solver


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    comment, accuracy, emin, emax, nlong, distances, energies = \
        read_single_potential("data/PC.in")

    pot = Potential(distances, energies, -nlong)
    clong = -energies[-1] * distances[-1]**nlong

    print(f"Comment: {comment}")
    print(f"Accuracy: {accuracy}, EMIN: {emin}, EMAX: {emax}, NLONG: {nlong}")
    print(f"CLONG: {clong}")

    orbit_list, ED = orbiting_scan(pot, emin, nlong, clong)
    regions = orbiting_regions(orbit_list)
    EC = max(r[1] for r in regions) if regions else 0.0

    print(f"\n{len(orbit_list)} orbiting sets, EC = {EC:.14e}")
    print(f"{'E':>24s} {'B':>24s} {'RM':>24s}")
    for E_orb, b_orb, r_orb in orbit_list[:5]:
        print(f"{E_orb:24.14e} {b_orb:24.14e} {r_orb:24.14e}")
    print("...")
    for E_orb, b_orb, r_orb in orbit_list[-3:]:
        print(f"{E_orb:24.14e} {b_orb:24.14e} {r_orb:24.14e}")

    solver = CrossSectionSolver(pot, orbit_list, regions, nlong, clong, accuracy)
    solver.ED = ED

    # Test at energies from Fortran output
    test_energies = [1e-9, 6.48e-9, 5.9e-7, 5.37e-5, 3.48e-4, 5e-4]
    print("\nCross-section results:")
    for E_test in test_energies:
        cs = solver.compute(E_test, 30)
        print(f"  E={E_test:.2e}  Q(1)={cs[0]:10.1f}  Q(2)={cs[1]:10.1f}")
