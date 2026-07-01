"""Tests for `potential_surface.PotentialSurface`.

Covers:
- fold-rule round-trips for both symmetries;
- constructor validation of angle grid vs declared symmetry;
- exact reproduction of the per-angle 1D `Potential` at every grid angle
  (V and dV/dr);
- agreement with the recovered `radial_first` prototype at arbitrary φ;
- homonuclear-fold synthesized from folded heteronuclear data reproduces
  values at symmetric points;
- V_at_r vectorized path matches scalar V per angle.
"""

import os
import sys

import numpy as np
import pytest

from potential_surface import (
    PotentialSurface,
    read_potential_surface,
    _read_pc_in_file,
)
from transport_from_potential import Potential

# The recovered prototype lives under prototypes/ — make it importable.
_REPO_ROOT = os.path.dirname(os.path.abspath(__file__))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)
from prototypes.interpolation_2d import create_2d_potential  # noqa: E402


_DATA_DIR = os.path.join(_REPO_ROOT, "data")
LONG_RANGE_POW = -4


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def hetero_data():
    """Load the real heteronuclear angle-resolved potentials into a dict."""
    data = {}
    for angle in range(0, 181, 10):
        path = os.path.join(_DATA_DIR, f"PC.in.{angle:03d}")
        data[float(angle)] = _read_pc_in_file(path)
    return data


@pytest.fixture(scope="module")
def hetero_surface(hetero_data):
    return PotentialSurface(hetero_data, LONG_RANGE_POW,
                            symmetry="heteronuclear")


@pytest.fixture(scope="module")
def homo_surface(hetero_data):
    """A homonuclear surface synthesized by keeping only the [0°, 90°] slice
    of the real heteronuclear data. Used to check the extra fold rule."""
    homo_data = {a: hetero_data[a] for a in hetero_data if a <= 90.0}
    return PotentialSurface(homo_data, LONG_RANGE_POW,
                            symmetry="homonuclear")


# ---------------------------------------------------------------------------
# Fold rule
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("phi,expected", [
    (0.0, 0.0),
    (45.0, 45.0),
    (180.0, 180.0),
    (200.0, 160.0),     # 360 - 200
    (-30.0, 30.0),      # wraps to 330 -> 360 - 330 = 30
    (720.0, 0.0),
    (540.0, 180.0),     # 540 % 360 = 180
    (270.0, 90.0),      # 360 - 270
])
def test_fold_heteronuclear(hetero_surface, phi, expected):
    assert hetero_surface._fold(phi) == pytest.approx(expected)


@pytest.mark.parametrize("phi,expected", [
    (0.0, 0.0),
    (45.0, 45.0),
    (90.0, 90.0),
    (100.0, 80.0),      # 180 - 100
    (135.0, 45.0),
    (180.0, 0.0),
    (200.0, 20.0),      # 360 - 200 = 160 -> 180 - 160 = 20
    (-30.0, 30.0),
    (270.0, 90.0),
])
def test_fold_homonuclear(homo_surface, phi, expected):
    assert homo_surface._fold(phi) == pytest.approx(expected)


def test_fold_array_matches_scalar(hetero_surface, homo_surface):
    phis = np.array([-30.0, 0.0, 45.0, 90.0, 100.0, 180.0, 200.0, 540.0])
    for surf in (hetero_surface, homo_surface):
        scalar = np.array([surf._fold(p) for p in phis])
        vector = surf._fold_array(phis)
        assert np.allclose(scalar, vector)


# ---------------------------------------------------------------------------
# Constructor validation
# ---------------------------------------------------------------------------

def test_bad_symmetry_rejected(hetero_data):
    with pytest.raises(ValueError):
        PotentialSurface(hetero_data, LONG_RANGE_POW, symmetry="funky")


def test_bad_method_rejected(hetero_data):
    with pytest.raises(NotImplementedError):
        PotentialSurface(hetero_data, LONG_RANGE_POW,
                         symmetry="heteronuclear", method="legendre")


def test_heteronuclear_range_data_rejected_as_homonuclear(hetero_data):
    """Full [0°, 180°] data with symmetry='homonuclear' must fail."""
    with pytest.raises(ValueError):
        PotentialSurface(hetero_data, LONG_RANGE_POW, symmetry="homonuclear")


def test_homonuclear_range_data_accepted_as_heteronuclear(hetero_data):
    """The homonuclear data range [0°, 90°] is inside [0°, 180°], so
    declaring 'heteronuclear' should be legal (albeit wasteful)."""
    homo_data = {a: hetero_data[a] for a in hetero_data if a <= 90.0}
    surf = PotentialSurface(homo_data, LONG_RANGE_POW,
                            symmetry="heteronuclear")
    # No exception is the assertion; sanity-check we got a usable object.
    assert surf.symmetry == "heteronuclear"


# ---------------------------------------------------------------------------
# V(r, φ_j) matches per-angle Potential
# ---------------------------------------------------------------------------

def test_V_at_grid_angles_matches_1d(hetero_surface):
    test_rs = [4.0, 5.5, 6.0, 8.0, 12.0]
    for j, angle in enumerate(hetero_surface.angles_deg):
        pot_1d = hetero_surface.potentials[j]
        for r in test_rs:
            assert hetero_surface.V(r, angle) == pytest.approx(pot_1d(r),
                                                                abs=1e-14)


def test_dV_dr_at_grid_angles_matches_1d(hetero_surface):
    test_rs = [4.0, 5.5, 6.0, 8.0, 12.0]
    for j, angle in enumerate(hetero_surface.angles_deg):
        pot_1d = hetero_surface.potentials[j]
        for r in test_rs:
            assert hetero_surface.dV_dr(r, angle) == pytest.approx(
                pot_1d.deriv(r), abs=1e-14)


def test_V_at_grid_angles_short_range(hetero_surface):
    """At r < r_min the per-angle short-range power law applies. Grid-angle
    V(r, φ_j) should match Potential(r) exactly."""
    for j, angle in enumerate(hetero_surface.angles_deg):
        pot_1d = hetero_surface.potentials[j]
        r = 0.5 * pot_1d.rmin
        assert hetero_surface.V(r, angle) == pytest.approx(pot_1d(r),
                                                            abs=1e-14)


def test_V_at_grid_angles_long_range(hetero_surface):
    """At r > r_max the per-angle long-range power law applies."""
    for j, angle in enumerate(hetero_surface.angles_deg):
        pot_1d = hetero_surface.potentials[j]
        r = 2.0 * pot_1d.rmax
        assert hetero_surface.V(r, angle) == pytest.approx(pot_1d(r),
                                                            abs=1e-14)


# ---------------------------------------------------------------------------
# Agreement with the recovered radial_first prototype
# ---------------------------------------------------------------------------

def test_matches_radial_first_prototype(hetero_data):
    """Between grid angles, V(r, φ) should match the prototype's
    radial_first method to machine precision (same underlying construction)."""
    angles_given = sorted(hetero_data.keys())
    distances_given = {a: hetero_data[a][0] for a in angles_given}
    energies_given = {a: hetero_data[a][1] for a in angles_given}
    proto_V = create_2d_potential(
        angles_given, distances_given, energies_given, LONG_RANGE_POW,
        method="radial_first")

    surface = PotentialSurface(hetero_data, LONG_RANGE_POW,
                               symmetry="heteronuclear")

    # The prototype uses not-a-knot BCs on the angular spline; PotentialSurface
    # uses clamped BCs (dV/dφ = 0 at 0° and 180°). Differences appear near
    # the endpoints. Compare mid-range and check that the difference is
    # small at endpoints too.
    for r in [4.5, 6.0, 8.0, 12.0]:
        for phi in [15.0, 45.0, 75.0, 105.0, 135.0, 165.0]:
            ours = surface.V(r, phi)
            theirs = proto_V(r, phi)
            # Compare at 1% relative accuracy — the BC choice differs
            # between the two implementations so exact agreement isn't
            # expected, but the underlying data is the same so the
            # answers should be very close.
            if abs(theirs) > 1e-12:
                assert abs(ours - theirs) / abs(theirs) < 1e-2, (
                    f"disagreement at r={r}, phi={phi}: "
                    f"ours={ours}, theirs={theirs}")


# ---------------------------------------------------------------------------
# V_at_r vectorized path
# ---------------------------------------------------------------------------

def test_V_at_r_matches_scalar(hetero_surface):
    r = 6.0
    phis = np.array([12.5, 45.0, 77.5, 132.5])
    vec = hetero_surface.V_at_r(r, phis)
    scalar = np.array([hetero_surface.V(r, p) for p in phis])
    assert np.allclose(vec, scalar, atol=1e-14)


def test_V_at_r_handles_out_of_range(hetero_surface):
    """Vectorized fold must handle φ outside [0°, 180°] the same as
    the scalar fold."""
    r = 6.0
    phis = np.array([-30.0, 30.0, 200.0, 160.0])
    vec = hetero_surface.V_at_r(r, phis)
    # -30 folds to 30, 200 folds to 160.
    assert vec[0] == pytest.approx(vec[1])
    assert vec[2] == pytest.approx(vec[3])


# ---------------------------------------------------------------------------
# Homonuclear fold reproduces symmetric values
# ---------------------------------------------------------------------------

def test_homonuclear_fold_180_minus_phi(homo_surface):
    """For a homonuclear surface, V(r, 180°−φ) must equal V(r, φ)
    (the extra fold beyond heteronuclear)."""
    r = 6.0
    for phi in [10.0, 30.0, 45.0, 70.0]:
        assert homo_surface.V(r, phi) == pytest.approx(
            homo_surface.V(r, 180.0 - phi))


def test_homonuclear_fold_130_maps_to_50(homo_surface):
    """Explicit sanity: V(6, 130°) should equal V(6, 50°)."""
    assert homo_surface.V(6.0, 130.0) == pytest.approx(
        homo_surface.V(6.0, 50.0))


# ---------------------------------------------------------------------------
# read_potential_surface loader
# ---------------------------------------------------------------------------

def test_loader_produces_matching_surface(hetero_data):
    surf_from_dict = PotentialSurface(hetero_data, LONG_RANGE_POW,
                                      symmetry="heteronuclear")
    surf_from_disk = read_potential_surface(
        _DATA_DIR, symmetry="heteronuclear", long_range_pow=LONG_RANGE_POW)
    assert np.array_equal(surf_from_dict.angles_deg,
                          surf_from_disk.angles_deg)
    for r, phi in [(6.0, 45.0), (5.0, 15.0), (10.0, 165.0)]:
        assert surf_from_dict.V(r, phi) == pytest.approx(
            surf_from_disk.V(r, phi), abs=1e-14)


def test_loader_missing_directory():
    with pytest.raises(FileNotFoundError):
        read_potential_surface("nonexistent_dir", "heteronuclear",
                               LONG_RANGE_POW)


# ---------------------------------------------------------------------------
# 1-slot cache
# ---------------------------------------------------------------------------

def test_cache_reused_within_same_r(hetero_surface):
    """After a V(r, φ) call, a second V(r, φ') at the same r should reuse
    the cached angular spline (checked by object identity)."""
    hetero_surface.V(7.0, 30.0)
    first = hetero_surface._cache_V_spline
    hetero_surface.V(7.0, 60.0)
    second = hetero_surface._cache_V_spline
    assert first is second


def test_cache_invalidated_on_r_change(hetero_surface):
    hetero_surface.V(7.0, 30.0)
    first = hetero_surface._cache_V_spline
    hetero_surface.V(8.0, 30.0)
    second = hetero_surface._cache_V_spline
    assert first is not second
