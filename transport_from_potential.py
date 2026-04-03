"""Compute transport cross-sections from intermolecular potentials.

Reproduces the results of the Fortran program PC.f95 by L. A. Viehland.
"""

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import quad
from scipy.optimize import root_scalar


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

def orbiting_scan(pot, emin, nlong, clong, n_samples=300):
    """Find orbiting parameters by root-finding the region boundaries.

    Returns (orbit_list, ED) where orbit_list is [(E, b, R), ...] sorted
    by increasing energy.
    """
    def c2_func(r):
        return pot.deriv2(r) + 3 * pot.deriv(r) / r

    def eorb(r):
        return pot(r) + r * pot.deriv(r) / 2

    # Find inner boundary where c2 = 0 (this gives EC)
    r_scan = np.logspace(np.log10(pot.rmin), 3, 5000)
    c2_vals = np.array([c2_func(r) for r in r_scan])
    sign_changes = np.where((c2_vals[:-1] >= 0) & (c2_vals[1:] < 0))[0]

    if len(sign_changes) == 0:
        raise RuntimeError("Could not find inner orbiting boundary")

    r_min = root_scalar(c2_func, bracket=(r_scan[sign_changes[0]], r_scan[sign_changes[0] + 1]),
                        method='brentq', xtol=1e-14).root
    EC = eorb(r_min)

    # Find outer boundary where Eorb(r) = emin
    eorb_vals = np.array([eorb(r) for r in r_scan])
    mask = (r_scan > r_min) & (eorb_vals[:-1] if len(eorb_vals) > len(r_scan) else np.ones(len(r_scan), dtype=bool))
    r_outer = r_scan[r_scan > r_min]
    eorb_outer = np.array([eorb(r) for r in r_outer])
    crossings = np.where((eorb_outer[:-1] - emin > 0) & (eorb_outer[1:] - emin <= 0))[0]

    if len(crossings) > 0:
        r_max = root_scalar(lambda r: eorb(r) - emin,
                            bracket=(r_outer[crossings[0]], r_outer[crossings[0] + 1]),
                            method='brentq', xtol=1e-10).root
    else:
        r_max = r_scan[-1]

    # Sample on fixed log-spaced grid
    r_grid = np.logspace(np.log10(r_min + 1e-10), np.log10(r_max), n_samples)
    orbit_list = []
    for r in r_grid:
        c1 = pot.deriv(r)
        c2 = pot.deriv2(r) + 3 * c1 / r
        veff = pot(r) + r * c1 / 2
        if c1 > 0 and c2 < 0 and veff > 0:
            orbit_list.append((veff, np.sqrt(r**3 * c1 / (2 * veff)), r))

    orbit_list.reverse()  # increasing energy

    ed = orbit_list[0][0] if (nlong == 3 and clong < 0 and orbit_list) else 0.0
    return orbit_list, ed


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
        self.pot = pot
        self.orbit_list = orbit_list
        self.regions = regions
        self.nlong = nlong
        self.clong = clong
        self.accuracy = accuracy
        self.acc1 = 0.8 * accuracy
        self.EC = max(r[1] for r in regions) if regions else 0.0
        self.ED = 0.0  # Set externally for nlong=3 case
        self.EC2 = 10 * self.EC

    # -- Turning point finder ------------------------------------------------

    def find_turning_point(self, E, b, r_guess=None):
        """Find R where E - V(R) - E*b²/R² = 0."""
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
        scale = 0.999 if (self.EC > 0 and 0.99 * self.EC < E < 1.01 * self.EC) else 0.95

        # Find r_pos (Y > 0) and r_neg (Y < 0)
        r_pos = r_neg = None
        r = r0
        for _ in range(500):
            if Y(r) > 0:
                r_pos = r
                break
            r /= scale
        if r_pos is None:
            if Y(1e4) > 0:
                r_pos = 1e4

        if r_pos is not None:
            r = r_pos
            for _ in range(2000):
                r *= scale
                if r < 0.3 * pot.rmin:
                    break
                if Y(r) < 0:
                    r_neg = r
                    break

        if r_neg is None or r_pos is None:
            r = 0.5 * pot.rmin
            if Y(r) < 0:
                r_neg = r
                r = pot.rmin
                for _ in range(2000):
                    r /= scale
                    if Y(r) > 0:
                        r_pos = r
                        break
                    if r > 1e6:
                        break

        if r_neg is None or r_pos is None:
            raise RuntimeError(f"Could not bracket turning point: E={E}, b={b}")

        lo, hi = min(r_neg, r_pos), max(r_neg, r_pos)
        r = root_scalar(Y, bracket=[lo, hi], method='brentq',
                        xtol=1e-14, rtol=1e-14).root

        # Adjust b for consistency
        bc = 1 - pot(r) / E
        b_adj = r * np.sqrt(abs(bc))
        return r, b_adj

    # -- Orbiting at specific energy -----------------------------------------

    def find_orbiting_at_energy(self, E):
        """Directly solve for orbiting (b, R) at energy E."""
        pot = self.pot
        r_scan = np.logspace(np.log10(pot.rmin), 3, 3000)

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
                    bracket=[r_lo, r_hi], method='brentq', xtol=1e-12).root
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
                denom = abs(denom) + 1e-30

            def ga(y):
                r = rm / np.cos(np.pi * (y + 1) / 4)
                val = 1 - (b_adj / r)**2 - pot(r) / E
                if val <= 0:
                    val = (r - rm) * (2 * b_adj**2 / rm**3 - pot(rm) / E)
                if val <= 0:
                    return 0.0
                return 1 - b_adj / rm * np.sin(np.pi * (y + 1) / 4) / np.sqrt(val)

            result, _ = quad(ga, -1, 1, limit=200, epsrel=self.acc1)
            return HPI * result

        else:
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
                    denom = abs(denom) + 1e-30

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
                    func = ea if abs(y) < 1e-15 else eb
                return func

            result, _ = quad(gb, -1, 1, limit=200, epsrel=self.acc1)
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
                          limit=200, epsrel=self.acc1)
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
