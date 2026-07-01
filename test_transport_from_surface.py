"""Tests for `transport_from_surface.compute_Q_ell`.

Uses small quadrature grids so the whole file runs in a few seconds.
The main correctness check is the **isotropic limit**: if V(r, φ) doesn't
depend on φ, then V̂tilde(r) = V(r) regardless of orientation, so
`compute_Q_ell` on the isotropic surface must equal a direct
single-angle Gauss-Legendre b-integration on the underlying 1D
`Potential` — up to machine precision when the same b nodes are used.
"""

import os

import numpy as np
import pytest

from potential_surface import PotentialSurface, _read_pc_in_file
from transport_from_potential import CrossSectionSolver, Potential
from transport_from_surface import compute_Q_ell, compute_Q_ell_grid


_DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
LONG_RANGE_POW = -4
MU_IJ = 1000.0
MU_I = 1000.0
ELL_LEN = 2.0


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def isotropic_surface():
    """Surface built from the φ=90° slice replicated across every grid
    angle — so V(r, φ) is φ-independent."""
    ref = _read_pc_in_file(os.path.join(_DATA_DIR, "PC.in.090"))
    per_angle = {float(a): (list(ref[0]), list(ref[1]))
                 for a in range(0, 181, 10)}
    return PotentialSurface(per_angle, LONG_RANGE_POW,
                             symmetry="heteronuclear")


@pytest.fixture(scope="module")
def isotropic_1d_potential():
    """The underlying 1D `Potential` for the isotropic-surface fixture."""
    ref = _read_pc_in_file(os.path.join(_DATA_DIR, "PC.in.090"))
    return Potential(ref[0], ref[1], LONG_RANGE_POW)


@pytest.fixture(scope="module")
def anisotropic_surface():
    """Real heteronuclear surface — for basic runtime / positivity checks."""
    per_angle = {}
    for a in range(0, 181, 10):
        per_angle[float(a)] = _read_pc_in_file(
            os.path.join(_DATA_DIR, f"PC.in.{a:03d}"))
    return PotentialSurface(per_angle, LONG_RANGE_POW,
                             symmetry="heteronuclear")


def _q_ell_1d_reference(pot_1d, epsilon, ell_max, n_b, b_max, accuracy):
    """Direct single-angle Q^(ℓ) with the same Gauss-Legendre b-quadrature
    that `compute_Q_ell` uses internally. Reference for the isotropic
    consistency check.
    """
    x_b, w_b = np.polynomial.legendre.leggauss(n_b)
    b_nodes = 0.5 * (x_b + 1.0) * b_max
    solver = CrossSectionSolver(pot_1d, orbit_list=[], regions=[],
                                 nlong=4, clong=0.0, accuracy=accuracy)
    ell_range = np.arange(1, ell_max + 1)
    Q = np.zeros(ell_max)
    for i, b in enumerate(b_nodes):
        if b == 0.0:
            continue
        chi = solver.deflection_angle(epsilon, b, naitk=0, BO=None,
                                       RO=None, ROMAX=0)
        Q += b * w_b[i] * (1.0 - np.cos(chi) ** ell_range)
    return Q * np.pi * b_max


# ---------------------------------------------------------------------------
# Shape / type
# ---------------------------------------------------------------------------

def test_return_shape(anisotropic_surface):
    Q = compute_Q_ell(
        anisotropic_surface, epsilon=1e-6, L_mag=0.0,
        mu_ij=MU_IJ, mu_i=MU_I, ell_len=ELL_LEN,
        ell_max=5,
        n_b=4, n_cos_phi=2, n_psi=2, n_alpha=2)
    assert Q.shape == (5,)
    assert Q.dtype == np.float64


def test_all_positive(anisotropic_surface):
    Q = compute_Q_ell(
        anisotropic_surface, epsilon=1e-6, L_mag=0.0,
        mu_ij=MU_IJ, mu_i=MU_I, ell_len=ELL_LEN,
        ell_max=4,
        n_b=6, n_cos_phi=2, n_psi=2, n_alpha=2)
    assert np.all(Q > 0)


def test_grid_wrapper_shape(anisotropic_surface):
    """The batched wrapper returns (n_energies, ell_max)."""
    Qs = compute_Q_ell_grid(
        anisotropic_surface, epsilon_array=[1e-6, 5e-6],
        L_mag=0.0, mu_ij=MU_IJ, mu_i=MU_I, ell_len=ELL_LEN,
        ell_max=3,
        n_b=4, n_cos_phi=2, n_psi=2, n_alpha=2)
    assert Qs.shape == (2, 3)


# ---------------------------------------------------------------------------
# Isotropic-limit consistency
# ---------------------------------------------------------------------------

def test_isotropic_matches_1d_reference(isotropic_surface,
                                          isotropic_1d_potential):
    """With an isotropic surface, V̂tilde(r) = V(r) at every orientation,
    so `compute_Q_ell` with any orientation grid must reproduce a direct
    single-angle Gauss-Legendre b-integration on the underlying 1D
    `Potential`, up to floating-point round-off."""
    epsilon = 1e-6
    ell_max = 5
    n_b = 6
    b_max = 40.0     # covers most of the tail contribution at low ε
    accuracy = 1e-4

    Q_2d = compute_Q_ell(
        isotropic_surface, epsilon=epsilon, L_mag=0.0,
        mu_ij=MU_IJ, mu_i=MU_I, ell_len=ELL_LEN,
        ell_max=ell_max,
        n_b=n_b, n_cos_phi=2, n_psi=2, n_alpha=2,
        b_max=b_max, accuracy=accuracy)
    Q_1d = _q_ell_1d_reference(
        isotropic_1d_potential, epsilon, ell_max, n_b, b_max, accuracy)

    for ell in range(1, ell_max + 1):
        assert Q_2d[ell - 1] == pytest.approx(Q_1d[ell - 1],
                                                rel=1e-10, abs=1e-14), \
            f"ell={ell}: Q_2d={Q_2d[ell-1]}, Q_1d={Q_1d[ell-1]}"


def test_isotropic_invariant_under_orientation_grid(isotropic_surface):
    """For an isotropic surface, changing orientation-grid resolution
    must not change Q^(ℓ) — every orientation gives the same χ(b), so
    the average is exact regardless of the number of nodes used."""
    epsilon = 1e-6
    ell_max = 3
    kwargs = dict(
        epsilon=epsilon, L_mag=0.0,
        mu_ij=MU_IJ, mu_i=MU_I, ell_len=ELL_LEN,
        ell_max=ell_max, n_b=6, b_max=20.0, accuracy=1e-4)
    Q_coarse = compute_Q_ell(isotropic_surface,
                               n_cos_phi=2, n_psi=2, n_alpha=2, **kwargs)
    Q_fine = compute_Q_ell(isotropic_surface,
                             n_cos_phi=4, n_psi=6, n_alpha=6, **kwargs)
    for ell in range(1, ell_max + 1):
        assert Q_coarse[ell - 1] == pytest.approx(Q_fine[ell - 1],
                                                    rel=1e-10, abs=1e-14)


# ---------------------------------------------------------------------------
# Deriv-mode independence in the isotropic limit
# ---------------------------------------------------------------------------

def test_isotropic_fd_and_analytic_agree(isotropic_surface):
    """In the isotropic limit the analytic and FD derivative modes of
    Vtilde reduce to the same value (∂ᵩV = 0 and Δφ multiplies zero),
    so Q^(ℓ) must agree exactly across the two modes."""
    common = dict(
        surface=isotropic_surface, epsilon=1e-6, L_mag=0.0,
        mu_ij=MU_IJ, mu_i=MU_I, ell_len=ELL_LEN,
        ell_max=3, n_b=6, n_cos_phi=2, n_psi=2, n_alpha=2,
        b_max=20.0, accuracy=1e-4)
    Q_fd = compute_Q_ell(derivative="fd", **common)
    Q_an = compute_Q_ell(derivative="analytic", **common)
    for ell in range(1, 4):
        assert Q_fd[ell - 1] == pytest.approx(Q_an[ell - 1],
                                                rel=1e-10, abs=1e-14)
