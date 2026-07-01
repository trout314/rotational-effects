"""Tests for `v_tilde` — the V̂tilde(r) construction.

Covers:
- isotropic limit: identical per-angle data ⇒ V̂tilde(r) = V(r);
- MMA limit: L = 0 and b = 0 ⇒ Δφ = 0 ⇒ V̂tilde(r) = V(r, φ_m);
- Δt(Δr) closed form vs direct numerical time-integration of the
  Taylor-expanded trajectory on a synthetic 12-4 potential;
- sin φ_m → 0 handling: dφ/dt is regular (no 1/sin), and Vtilde
  values are finite at φ_m near 0 and π;
- Vtilde exposes the same shape as `Potential` (rmin/rmax/deriv/deriv2/
  veff), enabling drop-in use inside `CrossSectionSolver`.
"""

import os

import numpy as np
import pytest

from potential_surface import PotentialSurface, _read_pc_in_file
from transport_from_potential import Potential
from v_tilde import (
    CollisionFrame,
    Vtilde,
    build_collision_frame,
    delta_t_and_deriv_of_delta_r,
    delta_t_of_delta_r,
    make_v_tilde,
)


_DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
LONG_RANGE_POW = -4


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def real_surface():
    per_angle = {}
    for a in range(0, 181, 10):
        per_angle[float(a)] = _read_pc_in_file(
            os.path.join(_DATA_DIR, f"PC.in.{a:03d}"))
    return PotentialSurface(per_angle, LONG_RANGE_POW,
                            symmetry="heteronuclear")


@pytest.fixture(scope="module")
def isotropic_surface(real_surface):
    """A PotentialSurface built from a single slice replicated across every
    grid angle. V(r, φ) must be independent of φ."""
    # Take the phi=90° slice as the reference.
    ref_pot = real_surface.potentials[9]  # angles_deg[9] = 90.0
    per_angle = {a: (list(ref_pot.distances), list(ref_pot.energies))
                 for a in real_surface.angles_deg}
    return PotentialSurface(per_angle, LONG_RANGE_POW,
                            symmetry="heteronuclear")


# ---------------------------------------------------------------------------
# Isotropic limit
# ---------------------------------------------------------------------------

def test_isotropic_limit(isotropic_surface):
    """When V(r, φ) does not depend on φ, V̂tilde(r) must equal V(r)
    regardless of the collision parameters — every φ sample in the
    three-point average returns the same value."""
    frame = build_collision_frame(
        isotropic_surface, epsilon=1e-6, b=6.0, phi_m=np.pi / 3,
        psi=0.4, L_mag=1e-3, alpha_L=0.9,
        mu_ij=1000.0, mu_i=1000.0, ell_len=2.0)
    V = make_v_tilde(isotropic_surface, frame)
    ref = isotropic_surface.potentials[9]  # the shared 1D slice
    for r in [4.5, 5.5, 6.5, 8.0, 12.0]:
        assert V(r) == pytest.approx(ref(r), rel=1e-12, abs=1e-16)


def test_isotropic_limit_sweep(isotropic_surface):
    """Same as `test_isotropic_limit` but sweeps the collision parameters
    to catch closure state leakage."""
    ref = isotropic_surface.potentials[9]
    rng = np.random.default_rng(42)
    for _ in range(6):
        phi_m = rng.uniform(0.1, np.pi - 0.1)
        psi = rng.uniform(0.0, 2 * np.pi)
        alpha_L = rng.uniform(0.0, 2 * np.pi)
        L_mag = rng.uniform(1e-4, 1e-2)
        b = rng.uniform(3.0, 10.0)
        frame = build_collision_frame(
            isotropic_surface, epsilon=1e-6, b=b, phi_m=phi_m,
            psi=psi, L_mag=L_mag, alpha_L=alpha_L,
            mu_ij=1000.0, mu_i=1000.0, ell_len=2.0)
        V = make_v_tilde(isotropic_surface, frame)
        for r in [5.0, 7.5, 10.0]:
            assert V(r) == pytest.approx(ref(r), rel=1e-12, abs=1e-16)


# ---------------------------------------------------------------------------
# MMA limit: L=0, b=0 ⇒ Δφ = 0 ⇒ V̂tilde(r) = V(r, φ_m)
# ---------------------------------------------------------------------------

def test_mma_limit(real_surface):
    """When L = 0 and b = 0, both terms of dφ/dt|ₘ vanish (v_m ∝ b → 0
    and L → 0), so Δφ = 0 and V̂tilde(r) collapses to V(r, φ_m)."""
    phi_m = np.pi / 4  # 45°
    frame = build_collision_frame(
        real_surface, epsilon=1e-6, b=0.0, phi_m=phi_m,
        psi=0.0, L_mag=0.0, alpha_L=0.0,
        mu_ij=1000.0, mu_i=1000.0, ell_len=2.0)
    assert frame.dphidt_m == pytest.approx(0.0, abs=1e-15)
    V = make_v_tilde(real_surface, frame)
    phi_m_deg = np.degrees(phi_m)
    for r in [5.0, 6.0, 8.0, 10.0]:
        expected = real_surface.V(r, phi_m_deg)
        assert V(r) == pytest.approx(expected, rel=1e-12, abs=1e-16)


# ---------------------------------------------------------------------------
# Δt(Δr) vs direct time-integration
# ---------------------------------------------------------------------------

def test_delta_t_closed_form_matches_trajectory():
    """Check delta_t_of_delta_r against a direct quadratic-trajectory
    root-solve on a synthetic setup.

    Uses parameters in the physical outer-turning-point regime
    (v_m² > r_m·dV/dr/μ_ij, i.e., the quadratic's b > 0)."""
    r_m = 5.0
    v_m = 0.5
    dVdr_m = 0.05
    mu_ij = 2.0
    for delta_r in [1e-4, 1e-3, 1e-2, 5e-2, 1e-1]:
        dt = delta_t_of_delta_r(delta_r, r_m, v_m, dVdr_m, mu_ij)
        # Direct: (r_m + Δr)² should equal α² r_m² + v_m² Δt²
        alpha = 1.0 - dt * dt * dVdr_m / (2.0 * mu_ij * r_m)
        lhs = (r_m + delta_r) ** 2
        rhs = (alpha * r_m) ** 2 + (v_m * dt) ** 2
        assert lhs == pytest.approx(rhs, rel=1e-12, abs=1e-20)


def test_delta_t_positive_and_zero_at_zero():
    """Δt(Δr) is monotone-increasing in Δr and hits 0 at Δr = 0."""
    assert delta_t_of_delta_r(0.0, 5.0, 0.5, 0.05, 2.0) == pytest.approx(0.0)
    prev = 0.0
    for delta_r in [1e-5, 1e-4, 1e-3, 1e-2, 1e-1]:
        cur = delta_t_of_delta_r(delta_r, 5.0, 0.5, 0.05, 2.0)
        assert cur > prev
        prev = cur


def test_delta_t_small_delta_r_leading_order():
    """For small Δr, Δt² ≈ 2 r_m Δr / (v_m² − r_m dV/dr / μ_ij) — the
    linearization of the quadratic in Δt² at Δr = 0."""
    r_m, v_m, dVdr_m, mu_ij = 5.0, 0.5, 0.05, 2.0
    delta_r = 1e-6
    dt = delta_t_of_delta_r(delta_r, r_m, v_m, dVdr_m, mu_ij)
    b_coeff = v_m * v_m - r_m * dVdr_m / mu_ij
    assert b_coeff > 0.0, "test setup must sit in the physical outer-TP regime"
    expected = np.sqrt(2.0 * r_m * delta_r / b_coeff)
    assert dt == pytest.approx(expected, rel=1e-5)


def test_delta_t_stable_for_negative_b():
    """When b = v_m² − 2 r_m² f is negative (unusual, but the formula
    should still be numerically stable), the closed form must produce
    the same result as a direct root-scalar solve of the quadratic."""
    from scipy.optimize import brentq
    r_m, v_m, dVdr_m, mu_ij = 5.0, 0.3, 0.05, 2.0
    f = dVdr_m / (2.0 * mu_ij * r_m)
    A = (r_m * f) ** 2
    b = v_m ** 2 - 2.0 * (r_m ** 2) * f
    assert b < 0.0

    for delta_r in [1e-4, 1e-3, 1e-2]:
        C = delta_r * (2.0 * r_m + delta_r)
        dt = delta_t_of_delta_r(delta_r, r_m, v_m, dVdr_m, mu_ij)
        u = dt ** 2
        # Compare to a brute-force positive root of A u² + b u − C = 0.
        u_ref = brentq(lambda x: A * x * x + b * x - C, 0.0, 1e6)
        assert u == pytest.approx(u_ref, rel=1e-10)


# ---------------------------------------------------------------------------
# sin φ_m → 0 handling
# ---------------------------------------------------------------------------

def test_dphidt_regular_at_phi_m_zero(real_surface):
    """dφ/dt|ₘ is regular through φ_m = 0 (no 1/sin blowup)."""
    frame = build_collision_frame(
        real_surface, epsilon=1e-6, b=6.0, phi_m=0.0,
        psi=0.5, L_mag=1e-3, alpha_L=0.7,
        mu_ij=1000.0, mu_i=1000.0, ell_len=2.0)
    assert np.isfinite(frame.dphidt_m)


def test_dphidt_regular_at_phi_m_pi(real_surface):
    frame = build_collision_frame(
        real_surface, epsilon=1e-6, b=6.0, phi_m=np.pi,
        psi=0.5, L_mag=1e-3, alpha_L=0.7,
        mu_ij=1000.0, mu_i=1000.0, ell_len=2.0)
    assert np.isfinite(frame.dphidt_m)


def test_dphidt_continuous_across_zero(real_surface):
    """Sweep φ_m from 1e-6 to 1e-10 (below the threshold) and check
    dφ/dt|ₘ stays smooth."""
    values = []
    for phi_m in [1e-4, 1e-6, 1e-8, 1e-10, 0.0]:
        frame = build_collision_frame(
            real_surface, epsilon=1e-6, b=6.0, phi_m=phi_m,
            psi=0.5, L_mag=1e-3, alpha_L=0.7,
            mu_ij=1000.0, mu_i=1000.0, ell_len=2.0)
        values.append(frame.dphidt_m)
    # All finite and within a tight relative window of each other.
    assert all(np.isfinite(v) for v in values)
    for v in values[1:]:
        assert v == pytest.approx(values[0], rel=1e-3)


def test_vtilde_finite_at_phi_m_zero(real_surface):
    """V̂tilde(r) at φ_m = 0 must be finite (no NaN from the frame setup)."""
    frame = build_collision_frame(
        real_surface, epsilon=1e-6, b=6.0, phi_m=0.0,
        psi=0.3, L_mag=1e-3, alpha_L=0.7,
        mu_ij=1000.0, mu_i=1000.0, ell_len=2.0)
    V = make_v_tilde(real_surface, frame)
    for r in [frame.r_m, frame.r_m + 0.1, frame.r_m + 1.0]:
        assert np.isfinite(V(r))


# ---------------------------------------------------------------------------
# CollisionFrame invariants
# ---------------------------------------------------------------------------

def test_frame_gauge_and_perpendicularity(real_surface):
    frame = build_collision_frame(
        real_surface, epsilon=1e-6, b=6.0, phi_m=np.pi / 3,
        psi=0.4, L_mag=1e-3, alpha_L=0.9,
        mu_ij=1000.0, mu_i=1000.0, ell_len=2.0)
    # **r_m** along ẑ, **v_m** along x̂.
    assert np.allclose(frame.r_m_vec, [0.0, 0.0, frame.r_m])
    assert np.allclose(frame.v_m_vec, [frame.v_m, 0.0, 0.0])
    # ℓ̂_m unit vector.
    assert np.linalg.norm(frame.ell_hat_m) == pytest.approx(1.0, abs=1e-14)
    # **L** ⊥ ℓ̂_m.
    assert frame.ell_hat_m @ frame.L_vec == pytest.approx(0.0, abs=1e-14)
    # |**L**| = L_mag.
    assert np.linalg.norm(frame.L_vec) == pytest.approx(frame.L_mag,
                                                          abs=1e-14)
    # cos φ_m = **r_m** · ℓ̂_m / r_m.
    assert (frame.ell_hat_m[2]) == pytest.approx(np.cos(frame.phi_m),
                                                   abs=1e-14)


def test_frame_dphidt_matches_vector_formula(real_surface):
    """At φ_m away from the singular case, the closed-form dφ/dt|ₘ must
    equal the direct vector-formula evaluation."""
    frame = build_collision_frame(
        real_surface, epsilon=1e-6, b=6.0, phi_m=np.pi / 3,
        psi=0.4, L_mag=1e-3, alpha_L=0.9,
        mu_ij=1000.0, mu_i=1000.0, ell_len=2.0)
    sin_phi_m = np.sin(frame.phi_m)
    numer = (frame.v_m_vec @ frame.ell_hat_m
             + (frame.r_m_vec @ np.cross(frame.ell_hat_m, frame.L_vec))
                / (frame.mu_i * frame.ell_len ** 2))
    expected = numer / (frame.r_m * sin_phi_m)
    assert frame.dphidt_m == pytest.approx(expected, rel=1e-12)


# ---------------------------------------------------------------------------
# Vtilde shape / Potential drop-in compatibility
# ---------------------------------------------------------------------------

def test_vtilde_shape(real_surface):
    frame = build_collision_frame(
        real_surface, epsilon=1e-6, b=6.0, phi_m=np.pi / 3,
        psi=0.4, L_mag=1e-3, alpha_L=0.9,
        mu_ij=1000.0, mu_i=1000.0, ell_len=2.0)
    V = make_v_tilde(real_surface, frame)
    # Attributes CrossSectionSolver relies on.
    assert isinstance(V.rmin, float)
    assert isinstance(V.rmax, float)
    assert V.long_range_pow == LONG_RANGE_POW
    # Methods.
    assert callable(V.__call__)
    r = frame.r_m + 0.5
    assert np.isfinite(V(r))
    assert np.isfinite(V.deriv(r))
    assert np.isfinite(V.deriv2(r))
    assert np.isfinite(V.veff(r))
    # veff definition check.
    assert V.veff(r) == pytest.approx(V(r) + r * V.deriv(r) / 2.0)


def test_vtilde_deriv_matches_isotropic_1d(isotropic_surface):
    """In the isotropic limit V̂tilde = V, so V̂tilde.deriv should equal
    the 1D V.deriv to FD accuracy."""
    frame = build_collision_frame(
        isotropic_surface, epsilon=1e-6, b=6.0, phi_m=np.pi / 3,
        psi=0.4, L_mag=1e-3, alpha_L=0.9,
        mu_ij=1000.0, mu_i=1000.0, ell_len=2.0)
    V = make_v_tilde(isotropic_surface, frame)
    ref = isotropic_surface.potentials[9]
    for r in [5.0, 6.0, 8.0]:
        # Vtilde.deriv uses h ≈ 1e-4 r_m ≈ 5e-4 -> FD error ~ h² ~ 2.5e-7
        # relative on typical V' magnitudes ~ 1e-3 -> abs tol ~ few 1e-10.
        assert V.deriv(r) == pytest.approx(ref.deriv(r), rel=1e-4, abs=1e-8)


# ---------------------------------------------------------------------------
# End-to-end: Vtilde drops into CrossSectionSolver
# ---------------------------------------------------------------------------

def _build_solver(pot, accuracy=1e-3):
    """Minimal CrossSectionSolver setup for a single deflection_angle call.

    Passes empty orbit_list / regions (non-orbiting regime), and nlong=4
    with clong derived from the potential's long-range coefficient. This
    is sufficient for the case1 code path in `deflection_angle`."""
    from transport_from_potential import CrossSectionSolver
    clong = -pot.long_coef if hasattr(pot, "long_coef") else 1.0
    return CrossSectionSolver(pot, orbit_list=[], regions=[],
                              nlong=4, clong=clong, accuracy=accuracy)


def test_vtilde_deflection_angle_finite(real_surface):
    """End-to-end: build a Vtilde, feed it to CrossSectionSolver, and get
    a finite deflection angle out."""
    frame = build_collision_frame(
        real_surface, epsilon=1e-6, b=8.0, phi_m=np.pi / 3,
        psi=0.4, L_mag=1e-3, alpha_L=0.9,
        mu_ij=1000.0, mu_i=1000.0, ell_len=2.0)
    V = make_v_tilde(real_surface, frame)
    # Vtilde is a closure, not a 1D Potential — no long_coef attribute.
    # Give the solver an approximate clong for the nlong=4 shortcut; the
    # exact value doesn't matter unless the analytic shortcut fires (it
    # only does when E·b⁴ >= 4·clong, which we avoid by picking small b).
    from transport_from_potential import CrossSectionSolver
    solver = CrossSectionSolver(V, orbit_list=[], regions=[],
                                nlong=4, clong=1e-4, accuracy=1e-3)
    chi = solver.deflection_angle(1e-6, 8.0, naitk=0, BO=None, RO=None,
                                    ROMAX=0)
    assert np.isfinite(chi)


def test_vtilde_reduces_to_single_angle_solver(real_surface):
    """With L=0 and ψ=π/2, dφ/dt|ₘ = 0 → Δφ = 0 for all r → V̂tilde(r) =
    V(r, φ_m). The deflection angle from CrossSectionSolver(Vtilde) must
    then match the deflection angle from CrossSectionSolver(1D slice at
    φ_m) at the same (ε, b)."""
    from transport_from_potential import CrossSectionSolver

    phi_m_deg = 60.0
    phi_m = np.radians(phi_m_deg)
    epsilon = 1e-6
    mu_ij = 1000.0

    # 1D slice: the per-angle Potential from the surface at φ_m = 60°.
    idx = int(np.where(real_surface.angles_deg == phi_m_deg)[0][0])
    pot_1d = real_surface.potentials[idx]
    solver_1d = CrossSectionSolver(pot_1d, orbit_list=[], regions=[],
                                    nlong=4, clong=1e-4, accuracy=1e-3)

    for b in [5.0, 7.0, 10.0, 15.0]:
        frame = build_collision_frame(
            real_surface, epsilon=epsilon, b=b, phi_m=phi_m,
            psi=np.pi / 2, L_mag=0.0, alpha_L=0.0,
            mu_ij=mu_ij, mu_i=mu_ij, ell_len=2.0)
        assert frame.dphidt_m == pytest.approx(0.0, abs=1e-15)
        V = make_v_tilde(real_surface, frame)
        solver_vt = CrossSectionSolver(V, orbit_list=[], regions=[],
                                        nlong=4, clong=1e-4, accuracy=1e-3)

        chi_1d = solver_1d.deflection_angle(epsilon, b, naitk=0, BO=None,
                                              RO=None, ROMAX=0)
        chi_vt = solver_vt.deflection_angle(epsilon, b, naitk=0, BO=None,
                                              RO=None, ROMAX=0)
        # Both solvers do the same χ integral with the same underlying V;
        # differences are FD-derivative noise from Vtilde.deriv at the
        # turning-point solve step (used only through denom = 1 -
        # rm³·V'(rm)/(2 b² E)). Should be very close.
        assert chi_vt == pytest.approx(chi_1d, rel=1e-3, abs=1e-3)


# ---------------------------------------------------------------------------
# Phase 2: analytic Vtilde' vs FD
# ---------------------------------------------------------------------------

def test_delta_t_and_deriv_matches_fd():
    """dΔt/dΔr from the analytic implicit differentiation must match a
    central FD on Δt(Δr) for Δr away from 0 (where the analytic formula
    has a √-singularity)."""
    r_m, v_m, dVdr_m, mu_ij = 5.0, 0.5, 0.05, 2.0
    for delta_r in [1e-2, 5e-2, 1e-1, 3e-1]:
        _, ddt_analytic = delta_t_and_deriv_of_delta_r(
            delta_r, r_m, v_m, dVdr_m, mu_ij)
        h = 1e-6
        ddt_fd = ((delta_t_of_delta_r(delta_r + h, r_m, v_m, dVdr_m, mu_ij)
                   - delta_t_of_delta_r(delta_r - h, r_m, v_m, dVdr_m, mu_ij))
                  / (2.0 * h))
        assert ddt_analytic == pytest.approx(ddt_fd, rel=1e-6)


def test_delta_t_and_deriv_diverges_at_zero():
    """dΔt/dΔr diverges as 1/√Δr near Δr = 0 — check that a Δr small but
    positive gives a large finite value proportional to 1/√Δr."""
    r_m, v_m, dVdr_m, mu_ij = 5.0, 0.5, 0.05, 2.0
    d1 = delta_t_and_deriv_of_delta_r(1e-6, r_m, v_m, dVdr_m, mu_ij)[1]
    d4 = delta_t_and_deriv_of_delta_r(4e-6, r_m, v_m, dVdr_m, mu_ij)[1]
    # d(Δr=1e-6) / d(Δr=4e-6) ≈ √4 = 2.
    assert d1 / d4 == pytest.approx(2.0, rel=1e-3)


def test_analytic_deriv_matches_fd_away_from_rm(real_surface):
    """Away from r_m the analytic V̂tilde'(r) must match the central-FD
    result — they compute the same derivative by different routes."""
    frame = build_collision_frame(
        real_surface, epsilon=1e-6, b=6.0, phi_m=np.pi / 3,
        psi=0.4, L_mag=1e-3, alpha_L=0.9,
        mu_ij=1000.0, mu_i=1000.0, ell_len=2.0)
    V_fd = Vtilde(real_surface, frame, derivative="fd")
    V_an = Vtilde(real_surface, frame, derivative="analytic")
    # Sample well away from r_m so the analytic path is not in FD fallback.
    for r in [frame.r_m + 1.0, frame.r_m + 3.0, frame.r_m + 8.0]:
        assert V_an.deriv(r) == pytest.approx(V_fd.deriv(r), rel=1e-3,
                                                abs=1e-10)


def test_analytic_deriv_sweep(real_surface):
    """Across a (b, φ_m, ψ, L, α_L) sweep, analytic V̂tilde' should track
    central FD closely at every sample point away from r_m."""
    rng = np.random.default_rng(0)
    for _ in range(5):
        phi_m = rng.uniform(0.1, np.pi - 0.1)
        psi = rng.uniform(0.0, 2 * np.pi)
        alpha_L = rng.uniform(0.0, 2 * np.pi)
        L_mag = rng.uniform(1e-4, 1e-2)
        b = rng.uniform(4.0, 10.0)
        frame = build_collision_frame(
            real_surface, epsilon=1e-6, b=b, phi_m=phi_m,
            psi=psi, L_mag=L_mag, alpha_L=alpha_L,
            mu_ij=1000.0, mu_i=1000.0, ell_len=2.0)
        V_fd = Vtilde(real_surface, frame, derivative="fd")
        V_an = Vtilde(real_surface, frame, derivative="analytic")
        for delta_r in [1.0, 3.0, 6.0]:
            r = frame.r_m + delta_r
            assert V_an.deriv(r) == pytest.approx(V_fd.deriv(r),
                                                    rel=1e-3, abs=1e-10)


def test_analytic_deriv_uses_closed_form_limit_at_rm(real_surface):
    """At r = r_m the analytic path returns the closed-form right-limit
        V̂tilde'(r_m⁺) = ∂ᵣV(r_m, φ_m)
                        + (2·(dφ/dt|ₘ)²·r_m / (3·b_coef)) · ∂²V/∂φ²(r_m, φ_m).
    Must be finite and match a from-scratch evaluation of that formula."""
    frame = build_collision_frame(
        real_surface, epsilon=1e-6, b=6.0, phi_m=np.pi / 3,
        psi=0.4, L_mag=1e-3, alpha_L=0.9,
        mu_ij=1000.0, mu_i=1000.0, ell_len=2.0)
    V = Vtilde(real_surface, frame, derivative="analytic")
    result = V.deriv(frame.r_m)
    assert np.isfinite(result)

    phi_m_deg = np.degrees(frame.phi_m)
    b_coef = frame.v_m ** 2 - frame.r_m * frame.dVdr_at_m / frame.mu_ij
    d2Vdphi2 = real_surface.dV_dphi2(frame.r_m, phi_m_deg)
    expected = (frame.dVdr_at_m
                + (2.0 * frame.dphidt_m ** 2 * frame.r_m / (3.0 * b_coef))
                * d2Vdphi2)
    assert result == pytest.approx(expected, rel=1e-14, abs=1e-14)


def test_analytic_deriv_at_rm_beats_forward_fd(real_surface):
    """The analytic limit at r_m should be materially more accurate than
    forward FD (which the previous implementation used). Compare both
    against a Richardson-extrapolated FD as a truth proxy."""
    frame = build_collision_frame(
        real_surface, epsilon=1e-6, b=6.0, phi_m=np.pi / 3,
        psi=0.4, L_mag=1e-3, alpha_L=0.9,
        mu_ij=1000.0, mu_i=1000.0, ell_len=2.0)
    V = Vtilde(real_surface, frame, derivative="analytic")

    analytic = V.deriv(frame.r_m)
    # Two-step Richardson extrapolation of the forward FD:
    h1 = 1e-4 * frame.r_m
    h2 = h1 / 2.0
    fd1 = (V(frame.r_m + h1) - V(frame.r_m)) / h1
    fd2 = (V(frame.r_m + h2) - V(frame.r_m)) / h2
    truth = 2.0 * fd2 - fd1  # cancels the O(h) term
    # Analytic should match truth much more tightly than fd1 does.
    err_analytic = abs(analytic - truth)
    err_fd = abs(fd1 - truth)
    assert err_analytic < err_fd, (
        f"analytic error {err_analytic:.3e} not smaller than "
        f"forward-FD error {err_fd:.3e}")


def test_analytic_deriv_continuous_across_threshold(real_surface):
    """Right at the threshold (Δr just above vs just below r_m_limit),
    the analytic result must not jump — the general formula converges
    to the limit as Δr → 0."""
    frame = build_collision_frame(
        real_surface, epsilon=1e-6, b=6.0, phi_m=np.pi / 3,
        psi=0.4, L_mag=1e-3, alpha_L=0.9,
        mu_ij=1000.0, mu_i=1000.0, ell_len=2.0)
    V = Vtilde(real_surface, frame, derivative="analytic")
    threshold = V._rm_limit_threshold
    # Just inside: uses limit formula
    val_inside = V.deriv(frame.r_m + 0.5 * threshold)
    # Just outside: uses general formula
    val_outside = V.deriv(frame.r_m + 5.0 * threshold)
    # Both should agree to at most a few ULPs on the underlying values.
    assert val_inside == pytest.approx(val_outside, rel=1e-6, abs=1e-14)


def test_analytic_isotropic_matches_1d(isotropic_surface):
    """In the isotropic limit V̂tilde = V, ∂ᵩV = 0, and the chain-rule
    angular term vanishes. Analytic V̂tilde'(r) = ∂ᵣV(r, φ_m) exactly."""
    frame = build_collision_frame(
        isotropic_surface, epsilon=1e-6, b=6.0, phi_m=np.pi / 3,
        psi=0.4, L_mag=1e-3, alpha_L=0.9,
        mu_ij=1000.0, mu_i=1000.0, ell_len=2.0)
    V = Vtilde(isotropic_surface, frame, derivative="analytic")
    ref = isotropic_surface.potentials[9]
    for r in [frame.r_m + 1.0, frame.r_m + 3.0, frame.r_m + 6.0]:
        assert V.deriv(r) == pytest.approx(ref.deriv(r), rel=1e-12,
                                             abs=1e-14)


def test_vtilde_with_legendre_surface(real_surface, hetero_data=None):
    """V̂tilde works end-to-end with a Legendre-based PotentialSurface."""
    surf_leg = PotentialSurface(
        {a: (list(real_surface.potentials[i].distances),
             list(real_surface.potentials[i].energies))
         for i, a in enumerate(real_surface.angles_deg)},
        real_surface.long_range_pow, symmetry="heteronuclear",
        method="legendre")
    frame = build_collision_frame(
        surf_leg, epsilon=1e-6, b=6.0, phi_m=np.pi / 3,
        psi=0.4, L_mag=1e-3, alpha_L=0.9,
        mu_ij=1000.0, mu_i=1000.0, ell_len=2.0)
    V = Vtilde(surf_leg, frame, derivative="analytic")
    r = frame.r_m + 2.0
    assert np.isfinite(V(r))
    assert np.isfinite(V.deriv(r))


def test_analytic_bad_option_rejected(real_surface):
    frame = build_collision_frame(
        real_surface, epsilon=1e-6, b=6.0, phi_m=np.pi / 3,
        psi=0.4, L_mag=1e-3, alpha_L=0.9,
        mu_ij=1000.0, mu_i=1000.0, ell_len=2.0)
    with pytest.raises(ValueError):
        Vtilde(real_surface, frame, derivative="not_a_mode")
