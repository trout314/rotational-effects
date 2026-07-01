"""2D potential surface V(r, φ) for angle-resolved transport cross sections.

Wraps a set of per-angle 1D `Potential` objects (from `transport_from_potential`)
with a `symmetry`-aware fold and angular interpolation. Reuses the per-angle
short/long-range extrapolation and analytic derivatives from `Potential` —
does not duplicate that code.
"""

from __future__ import annotations

import glob
import os
import re

import numpy as np
from scipy.interpolate import CubicSpline

from transport_from_potential import Potential


_ANGLE_FROM_FILENAME = re.compile(r"PC\.in\.(\d{3})$")


def _read_pc_in_file(filename):
    """Read one PC.in.XXX file.

    Format: line 1 = integer count of data points; following lines =
    whitespace-separated (r, V) pairs. Returns (distances, energies) as
    plain Python lists (matching what `Potential.__init__` expects).
    """
    with open(filename) as f:
        n = int(f.readline())
        data = np.loadtxt(f, max_rows=n)
    return data[:, 0].tolist(), data[:, 1].tolist()


class PotentialSurface:
    """V(r, φ) built from tabulated per-angle 1D potentials.

    Parameters
    ----------
    per_angle_data : dict[float, tuple[list, list]]
        Mapping angle-in-degrees → (distances, energies) for that angle.
    long_range_pow : float
        Long-range extrapolation exponent, passed through to each per-angle
        `Potential`.
    symmetry : {"heteronuclear", "homonuclear"}
        Selects the canonical angle range and fold rule. `"heteronuclear"`
        expects data on [0°, 180°] and only applies the reflection
        φ ↔ −φ (implemented as 360°−φ for φ > 180°). `"homonuclear"`
        expects data on [0°, 90°] and additionally applies φ ↔ 180°−φ.
    method : {"radial_first"}
        Angular interpolation method. `"radial_first"` is the only method
        supported in Phase 1 (see the plan for Phase 2's `"legendre"`).

    Notes
    -----
    Angular interpolation is a cubic spline in φ (degrees) with clamped
    boundary conditions dV/dφ = 0 at each endpoint, matching the reflection
    symmetries V(r, −φ) = V(r, φ) at the low end and (for heteronuclear)
    V(r, 180°+δ) = V(r, 180°−δ) at the high end; homonuclear adds the same
    at 90°.
    """

    _CANONICAL_MAX = {"heteronuclear": 180.0, "homonuclear": 90.0}
    _SUPPORTED_METHODS = ("radial_first",)

    def __init__(self, per_angle_data, long_range_pow, symmetry,
                 method="radial_first"):
        if symmetry not in self._CANONICAL_MAX:
            raise ValueError(
                f"symmetry must be one of {list(self._CANONICAL_MAX)}, "
                f"got {symmetry!r}")
        if method not in self._SUPPORTED_METHODS:
            raise NotImplementedError(
                f"method={method!r}: only {self._SUPPORTED_METHODS} "
                f"is supported in Phase 1.")

        canonical_max = self._CANONICAL_MAX[symmetry]
        angles = sorted(float(a) for a in per_angle_data.keys())
        if len(angles) < 4:
            raise ValueError(
                "need at least 4 angles for a cubic spline; "
                f"got {len(angles)}.")
        if angles[0] < -1e-9 or angles[-1] > canonical_max + 1e-9:
            raise ValueError(
                f"symmetry={symmetry!r} requires all data angles in "
                f"[0°, {canonical_max}°]; got range "
                f"[{angles[0]}°, {angles[-1]}°]")

        self.symmetry = symmetry
        self.method = method
        self.long_range_pow = long_range_pow
        self.angles_deg = np.array(angles, dtype=float)
        self.potentials = [
            Potential(per_angle_data[a][0], per_angle_data[a][1],
                      long_range_pow)
            for a in angles
        ]

        # 1-slot cache: at a given r, V̂tilde queries three φ samples in a
        # row, and χ quadrature revisits each r for both integrand branches.
        # We cache the angular spline(s) built at the last r seen.
        self._cache_r = None
        self._cache_V_spline = None
        self._cache_dV_dr_spline = None

    def _fold(self, phi):
        """Fold an arbitrary φ (degrees) into the canonical range."""
        phi = phi % 360.0
        if phi > 180.0:
            phi = 360.0 - phi
        if self.symmetry == "homonuclear" and phi > 90.0:
            phi = 180.0 - phi
        return phi

    def _fold_array(self, phi_array):
        """Vectorized version of `_fold` for `V_at_r`."""
        folded = np.mod(phi_array, 360.0)
        folded = np.where(folded > 180.0, 360.0 - folded, folded)
        if self.symmetry == "homonuclear":
            folded = np.where(folded > 90.0, 180.0 - folded, folded)
        return folded

    def _angular_spline(self, r, kind):
        """Build (or return the cached) angular cubic spline at radius r.

        `kind ∈ {"V", "dV_dr"}` selects value or first radial derivative.
        The last-built spline of each kind is cached and reused as long as r
        is unchanged.
        """
        if self._cache_r != r:
            self._cache_r = r
            self._cache_V_spline = None
            self._cache_dV_dr_spline = None

        attr = f"_cache_{kind}_spline"
        cached = getattr(self, attr)
        if cached is not None:
            return cached

        if kind == "V":
            vals = np.array([p(r) for p in self.potentials])
        elif kind == "dV_dr":
            vals = np.array([p.deriv(r) for p in self.potentials])
        else:
            raise ValueError(f"unknown spline kind {kind!r}")

        spline = CubicSpline(
            self.angles_deg, vals,
            bc_type=((1, 0.0), (1, 0.0)))
        setattr(self, attr, spline)
        return spline

    def V(self, r, phi):
        """V(r, φ) with φ in degrees."""
        return float(self._angular_spline(r, "V")(self._fold(phi)))

    def dV_dr(self, r, phi):
        """∂V/∂r at (r, φ) with φ in degrees."""
        return float(self._angular_spline(r, "dV_dr")(self._fold(phi)))

    def V_at_r(self, r, phi_array):
        """Vectorized V(r, φ) over an array of φ (degrees) at fixed r.

        Fast path for V̂tilde's ⅓ Σᵢ V(r, φₘ + sᵢΔφ) average — builds the
        angular spline once and evaluates at all three φ samples.
        """
        folded = self._fold_array(np.asarray(phi_array, dtype=float))
        return np.asarray(self._angular_spline(r, "V")(folded), dtype=float)


def read_potential_surface(directory, symmetry, long_range_pow,
                           pattern="PC.in.[0-9][0-9][0-9]"):
    """Load a `PotentialSurface` from a directory of PC.in.XXX files.

    Parameters
    ----------
    directory : str
        Path to the directory containing PC.in.XXX files.
    symmetry : {"heteronuclear", "homonuclear"}
        Passed straight through to `PotentialSurface`.
    long_range_pow : float
        Long-range extrapolation exponent (e.g. -4).
    pattern : str
        Glob pattern applied inside `directory`. The integer angle is
        parsed from the last three digits of each filename via
        `PC.in.<NNN>`.
    """
    files = sorted(glob.glob(os.path.join(directory, pattern)))
    if not files:
        raise FileNotFoundError(
            f"no files matching {pattern!r} in {directory!r}")

    per_angle_data = {}
    for filename in files:
        base = os.path.basename(filename)
        match = _ANGLE_FROM_FILENAME.search(base)
        if not match:
            raise ValueError(
                f"cannot parse angle from filename {base!r} (expected "
                f"suffix .NNN with three digits)")
        angle_deg = float(match.group(1))
        per_angle_data[angle_deg] = _read_pc_in_file(filename)
    return PotentialSurface(per_angle_data, long_range_pow, symmetry)
