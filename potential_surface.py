"""2D potential surface V(r, φ) for angle-resolved transport cross sections.

Wraps a set of per-angle 1D `Potential` objects (from `transport_from_potential`)
with a `symmetry`-aware fold and angular interpolation. Reuses the per-angle
short/long-range extrapolation and analytic derivatives from `Potential` —
does not duplicate that code.

Two angular-interpolation methods are supported:

- **radial_first**: at each r, evaluate per-angle 1D potentials, then build
  an angular cubic spline (in degrees, clamped-BC to match reflection
  symmetries). Reproduces the per-angle `Potential` exactly at every grid
  angle. Extrapolation in r is per-angle and always correct.
- **legendre**: at each r, evaluate per-angle 1D potentials, project onto a
  Legendre-in-cos(φ) basis using precomputed pseudoinverse weights, and
  reconstruct V(r, φ). Smoother in φ (analytic derivatives with respect to
  cos φ) which makes ∂ᵩV cheap. Extrapolation in r is per-angle and always
  correct because the projection uses per-angle `Potential.__call__` /
  `Potential.deriv` directly, not a spline in r.

Both methods expose the same API: `V(r, φ)`, `dV_dr(r, φ)`, `dV_dphi(r, φ)`,
`dV_dphi2(r, φ)`, `V_at_r(r, phi_array)`. Angles φ are in **degrees**;
`dV_dphi` returns ∂V/∂φ in **radians⁻¹** and `dV_dphi2` returns ∂²V/∂φ² in
**radians⁻²**, matching the physics convention in `v_tilde`.
"""

from __future__ import annotations

import glob
import os
import re

import numpy as np
from numpy.polynomial.legendre import legder, legval
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
    method : {"radial_first", "legendre"}
        Angular interpolation method (see module docstring).
    n_legendre_terms : int, optional
        Number of Legendre polynomial terms used in `legendre` projection.
        Defaults to the number of grid angles (exact interpolation at
        grid points). Ignored for `radial_first`.

    Notes
    -----
    For `radial_first`, the angular cubic spline uses clamped boundary
    conditions dV/dφ = 0 at each endpoint, matching the reflection
    symmetries V(r, −φ) = V(r, φ) at 0° and (for heteronuclear)
    V(r, 180°+δ) = V(r, 180°−δ) at 180°; homonuclear adds the same at 90°.
    For `legendre` with `symmetry="homonuclear"`, the expansion is
    restricted to even-ℓ Legendre polynomials only, which enforces
    V(r, π−φ) = V(r, φ) exactly (P_ℓ(−x) = (−1)^ℓ P_ℓ(x)).
    """

    _CANONICAL_MAX = {"heteronuclear": 180.0, "homonuclear": 90.0}
    _SUPPORTED_METHODS = ("radial_first", "legendre")

    def __init__(self, per_angle_data, long_range_pow, symmetry,
                 method="radial_first", n_legendre_terms=None):
        if symmetry not in self._CANONICAL_MAX:
            raise ValueError(
                f"symmetry must be one of {list(self._CANONICAL_MAX)}, "
                f"got {symmetry!r}")
        if method not in self._SUPPORTED_METHODS:
            raise NotImplementedError(
                f"method={method!r}: only {self._SUPPORTED_METHODS} "
                f"are supported.")

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

        if method == "legendre":
            self._setup_legendre(n_legendre_terms)
        else:
            self.n_legendre_terms = None

        # 1-slot cache: at a given r, V̂tilde queries three φ samples in a
        # row, and χ quadrature revisits each r for both integrand branches.
        # Both methods evaluate the 19 per-angle 1D potentials at r as their
        # first step; those grids are cached and shared. `radial_first`
        # additionally caches the angular cubic spline it builds from each
        # grid.
        self._cache_r = None
        self._cache_V_grid = None
        self._cache_dV_dr_grid = None
        self._cache_V_spline = None
        self._cache_dV_dr_spline = None

    # ------------------------------------------------------------------
    # Legendre-method setup
    # ------------------------------------------------------------------

    def _setup_legendre(self, n_legendre_terms):
        n_angles = len(self.angles_deg)
        n_terms = n_legendre_terms if n_legendre_terms is not None else n_angles
        if not (1 <= n_terms <= n_angles):
            raise ValueError(
                f"n_legendre_terms must be in [1, {n_angles}]; "
                f"got {n_terms}.")
        self.n_legendre_terms = n_terms

        # Active ℓ indices. For homonuclear V(r, π−φ) = V(r, φ), only
        # even-ℓ Legendre polynomials contribute (P_ℓ(−x) = (−1)^ℓ P_ℓ(x)),
        # so we restrict the basis accordingly. For heteronuclear we
        # use every ℓ up to n_terms.
        if self.symmetry == "homonuclear":
            self._legendre_ells = 2 * np.arange(n_terms, dtype=int)
        else:
            self._legendre_ells = np.arange(n_terms, dtype=int)
        # Padding length: highest ℓ used + 1. `legval` / `legder` need
        # a coefficient array indexed by ℓ, so we store the projected
        # coefficients padded with zeros at the unused (odd) positions.
        self._legendre_pad_length = int(self._legendre_ells[-1]) + 1

        # Projection matrix P[j, k] = P_{ell[k]}(cos φ_j).
        cos_phis = np.cos(np.radians(self.angles_deg))
        P_matrix = np.zeros((n_angles, n_terms))
        for k, ell in enumerate(self._legendre_ells):
            coef = np.zeros(ell + 1)
            coef[ell] = 1.0
            P_matrix[:, k] = legval(cos_phis, coef)

        # Pseudoinverse gives the least-squares projection weights.
        # For n_terms == n_angles this is the plain inverse (exact
        # interpolation at grid angles).
        self._legendre_weights = np.linalg.pinv(P_matrix)  # (n_terms, n_angles)

    # ------------------------------------------------------------------
    # Fold rule
    # ------------------------------------------------------------------

    def _fold(self, phi):
        """Fold an arbitrary φ (degrees) into the canonical range."""
        phi = phi % 360.0
        if phi > 180.0:
            phi = 360.0 - phi
        if self.symmetry == "homonuclear" and phi > 90.0:
            phi = 180.0 - phi
        return phi

    def _fold_array(self, phi_array):
        """Vectorized version of `_fold`."""
        folded = np.mod(phi_array, 360.0)
        folded = np.where(folded > 180.0, 360.0 - folded, folded)
        if self.symmetry == "homonuclear":
            folded = np.where(folded > 90.0, 180.0 - folded, folded)
        return folded

    # ------------------------------------------------------------------
    # Per-angle grid cache (shared by both methods)
    # ------------------------------------------------------------------

    def _get_grid(self, r, kind):
        """Return the array of per-angle V or dV/dr values at radius r.

        Uses the 1-slot r cache. On an r change, invalidates the grid
        arrays and (for `radial_first`) the angular splines built from
        them; on a hit for the same kind, returns the cached array.
        """
        if self._cache_r != r:
            self._cache_r = r
            self._cache_V_grid = None
            self._cache_dV_dr_grid = None
            self._cache_V_spline = None
            self._cache_dV_dr_spline = None

        if kind == "V":
            if self._cache_V_grid is None:
                self._cache_V_grid = np.array([p(r) for p in self.potentials])
            return self._cache_V_grid
        if kind == "dV_dr":
            if self._cache_dV_dr_grid is None:
                self._cache_dV_dr_grid = np.array(
                    [p.deriv(r) for p in self.potentials])
            return self._cache_dV_dr_grid
        raise ValueError(f"unknown grid kind {kind!r}")

    def _get_angular_spline(self, r, kind):
        """`radial_first`: build (or return cached) angular cubic spline."""
        grid = self._get_grid(r, kind)
        cache_attr = f"_cache_{kind}_spline"
        cached = getattr(self, cache_attr)
        if cached is not None:
            return cached
        spline = CubicSpline(
            self.angles_deg, grid,
            bc_type=((1, 0.0), (1, 0.0)))
        setattr(self, cache_attr, spline)
        return spline

    def _legendre_project(self, grid):
        """`legendre`: project a per-angle value array onto the Legendre
        basis. Returns a coefficient array padded to length
        `_legendre_pad_length` with zeros at unused (odd, for homonuclear)
        positions — the form `legval` / `legder` expects.
        """
        V_active = self._legendre_weights @ grid  # shape (n_terms,)
        V_padded = np.zeros(self._legendre_pad_length)
        V_padded[self._legendre_ells] = V_active
        return V_padded

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def V(self, r, phi):
        """V(r, φ) with φ in degrees."""
        phi_folded = self._fold(phi)
        if self.method == "radial_first":
            return float(self._get_angular_spline(r, "V")(phi_folded))
        # legendre
        V_ell = self._legendre_project(self._get_grid(r, "V"))
        return float(legval(np.cos(np.radians(phi_folded)), V_ell))

    def dV_dr(self, r, phi):
        """∂V/∂r at (r, φ) with φ in degrees."""
        phi_folded = self._fold(phi)
        if self.method == "radial_first":
            return float(self._get_angular_spline(r, "dV_dr")(phi_folded))
        # legendre
        V_ell = self._legendre_project(self._get_grid(r, "dV_dr"))
        return float(legval(np.cos(np.radians(phi_folded)), V_ell))

    def dV_dphi(self, r, phi):
        """∂V/∂φ at (r, φ), in radians⁻¹. Input φ in degrees.

        For `radial_first`, the angular cubic spline is a function of φ
        in degrees; its scipy-provided derivative is dV/dφ_deg, which is
        rescaled by 180/π to convert to dV/dφ_rad.

        For `legendre`, uses the closed-form
            ∂V/∂φ = -sin(φ) · Σ_ℓ V_ℓ(r) · P_ℓ'(cos φ),
        with V_ℓ(r) obtained from the same fresh projection used by
        `V(r, φ)`.
        """
        phi_folded = self._fold(phi)
        if self.method == "radial_first":
            spline = self._get_angular_spline(r, "V")
            deriv_per_deg = float(spline(phi_folded, 1))
            return deriv_per_deg * (180.0 / np.pi)
        # legendre
        V_ell = self._legendre_project(self._get_grid(r, "V"))
        V_ell_deriv = legder(V_ell)
        phi_rad = np.radians(phi_folded)
        return float(-np.sin(phi_rad) * legval(np.cos(phi_rad), V_ell_deriv))

    def dV_dphi2(self, r, phi):
        """∂²V/∂φ² at (r, φ), in radians⁻². Input φ in degrees.

        For `radial_first`, the second derivative of the angular cubic
        spline in φ_deg, rescaled by (180/π)² to convert to per-radian².

        For `legendre`, uses the closed-form second derivative of
            V(r, φ) = Σ_ℓ V_ℓ(r) · P_ℓ(cos φ)
        with respect to φ:
            ∂²V/∂φ² = -cos(φ) · Σ_ℓ V_ℓ · P_ℓ'(cos φ)
                      + sin²(φ) · Σ_ℓ V_ℓ · P_ℓ''(cos φ),
        obtained by differentiating -sin(φ)·Σ V_ℓ P_ℓ'(cos φ) once more
        with respect to φ.
        """
        phi_folded = self._fold(phi)
        if self.method == "radial_first":
            spline = self._get_angular_spline(r, "V")
            deriv2_per_deg2 = float(spline(phi_folded, 2))
            return deriv2_per_deg2 * (180.0 / np.pi) ** 2
        # legendre
        V_ell = self._legendre_project(self._get_grid(r, "V"))
        V_ell_d1 = legder(V_ell)
        V_ell_d2 = legder(V_ell_d1)
        phi_rad = np.radians(phi_folded)
        cos_phi = np.cos(phi_rad)
        sin_phi = np.sin(phi_rad)
        term1 = -cos_phi * legval(cos_phi, V_ell_d1)
        term2 = sin_phi * sin_phi * legval(cos_phi, V_ell_d2)
        return float(term1 + term2)

    def V_at_r(self, r, phi_array):
        """Vectorized V(r, φ) over an array of φ (degrees) at fixed r."""
        folded = self._fold_array(np.asarray(phi_array, dtype=float))
        if self.method == "radial_first":
            return np.asarray(self._get_angular_spline(r, "V")(folded),
                              dtype=float)
        # legendre
        V_ell = self._legendre_project(self._get_grid(r, "V"))
        return np.asarray(legval(np.cos(np.radians(folded)), V_ell),
                          dtype=float)


def read_potential_surface(directory, symmetry, long_range_pow,
                           pattern="PC.in.[0-9][0-9][0-9]",
                           method="radial_first", n_legendre_terms=None):
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
    method : str
        Angular interpolation method for `PotentialSurface`.
    n_legendre_terms : int, optional
        Passed straight through when method="legendre".
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
    return PotentialSurface(per_angle_data, long_range_pow, symmetry,
                            method=method,
                            n_legendre_terms=n_legendre_terms)
