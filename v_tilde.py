"""V̂tilde(r) construction for angle-resolved transport cross sections.

Given collision parameters at closest approach, `build_collision_frame`
returns a `CollisionFrame`; `make_v_tilde(surface, frame)` returns a
`Potential`-shaped closure V̂tilde(r) that packs the three-point angular
average from Eq. `v0_tilde` of `paper/main.tex`. That closure plugs into
the existing `transport_from_potential.CrossSectionSolver` in place of
the 1D `Potential`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import root_scalar

from transport_from_potential import DENOM_FLOOR
from potential_surface import PotentialSurface


@dataclass
class CollisionFrame:
    """Gauge-fixed kinematic state at closest approach.

    Conventions
    -----------
    - Coordinate axes: **r_m** = r_m ẑ, **v_m** = v_m x̂ (gauge-fixed;
      the physical problem is rotation-invariant about the lab origin).
    - Diatom orientation: ℓ̂_m at polar angle φ_m from ẑ, azimuth ψ measured
      from x̂ = v̂_m about ẑ. Formally
      ℓ̂_m = (sin φ_m cos ψ, sin φ_m sin ψ, cos φ_m).
    - Angular momentum: **L** ⊥ ℓ̂_m, magnitude L_mag, direction α_L in the
      tangent basis
      ê₁ = (ẑ - cos φ_m ℓ̂_m)/sin φ_m,  ê₂ = ℓ̂_m × ê₁,
      i.e., **L** = L_mag (cos α_L ê₁ + sin α_L ê₂). The (ẑ - cos φ_m ℓ̂_m)
      formula is 0/0 at sin φ_m = 0 but has a finite limit; the limit is
      applied when |sin φ_m| < `_SIN_PHI_THRESHOLD`.
    - Angles φ_m, ψ, α_L are all in **radians**.

    Fields
    ------
    Scalars consumed on every V̂tilde(r) call: `r_m`, `v_m`, `phi_m`,
    `mu_ij`, `dVdr_at_m`, `dphidt_m`. The rest are kept for testing,
    debugging, and the future analytic V̂tilde'(r) work in Phase 2.
    """

    # Hot-path scalars.
    r_m: float
    v_m: float
    phi_m: float
    mu_ij: float
    dVdr_at_m: float
    dphidt_m: float

    # Frame vectors (Phase-2-useful; also handy for tests).
    r_m_vec: np.ndarray
    v_m_vec: np.ndarray
    ell_hat_m: np.ndarray
    L_vec: np.ndarray

    # Physical parameters echoed for downstream use.
    L_mag: float
    psi: float
    alpha_L: float
    mu_i: float
    ell_len: float


_SIN_PHI_THRESHOLD = 1.0e-8


def delta_t_of_delta_r(delta_r, r_m, v_m, dVdr_m, mu_ij):
    """Closed-form Δt(Δr) at closest approach.

    Derivation. The trajectory near t_m is
    r(t_m+Δt) ≈ **r_m** + **v_m**·Δt + ½·**a_m**·Δt² with **a_m** =
    -(1/μ_ij)(dV/dr) r̂_m. Since **v_m** ⊥ **r_m** at the turning point,
    |r(t_m+Δt)|² = α² r_m² + v_m²·Δt² where
    α = 1 - Δt²·(dV/dr)/(2 μ_ij r_m). Setting |r(t_m+Δt)| = r_m + Δr and
    squaring gives a quadratic in Δt²:

        A · Δt⁴ + b · Δt² - Δr·(2·r_m + Δr) = 0,

    with A = r_m²·f², b = v_m² - 2·r_m²·f, f = (dV/dr)/(2·μ_ij·r_m). The
    physical root (Δt² > 0, → 0 as Δr → 0) is returned in its
    rationalized form

        Δt² = 2·Δr(2·r_m+Δr) / (b + √(b² + 4·A·Δr(2·r_m+Δr))),

    which avoids the near-cancellation of -b + √(…) that appears in the
    naive quadratic-formula answer when b > 0 (the physical outer-
    turning-point regime).

    Note. `paper/main.tex` Eq. `delta_t_from_delta_r` writes
    C = 4·Δr·(4·r_m+Δr)/r_m². Direct numerical time-integration of the
    quadratic trajectory reproduces this implementation's formula
    (equivalent to C = 4·Δr·(2·r_m+Δr)/r_m²) instead — the paper's `4·r_m`
    looks like a typo for `2·r_m`.
    """
    f = dVdr_m / (2.0 * mu_ij * r_m)
    A = (r_m * f) ** 2
    b = v_m ** 2 - 2.0 * (r_m ** 2) * f
    C_star = delta_r * (2.0 * r_m + delta_r)

    # Degenerate quadratic (dV/dr = 0 ⇒ A = 0): b·u = C_star.
    if A < DENOM_FLOOR:
        if abs(b) < DENOM_FLOOR:
            return 0.0
        u = C_star / b
        return np.sqrt(u) if u > 0.0 else 0.0

    disc = b * b + 4.0 * A * C_star
    sqrt_disc = np.sqrt(max(disc, 0.0))
    # Two algebraically-equivalent forms of the positive root of
    # A·u² + b·u − C_star = 0. Pick the one that avoids subtractive
    # cancellation based on sign(b):
    #   b > 0: 2·C_star / (b + √disc)  (denominator is a sum of positives)
    #   b < 0: (−b + √disc) / (2·A)    (numerator is a sum of positives)
    if b >= 0.0:
        u = 2.0 * C_star / (b + sqrt_disc)
    else:
        u = (-b + sqrt_disc) / (2.0 * A)
    return np.sqrt(u) if u > 0.0 else 0.0


def delta_t_and_deriv_of_delta_r(delta_r, r_m, v_m, dVdr_m, mu_ij):
    """Return (Δt, dΔt/dΔr) at the given Δr.

    dΔt/dΔr is obtained by implicit differentiation of the quadratic
    A·u² + b·u − C(Δr) = 0 in u = Δt². With C(Δr) = Δr·(2·r_m+Δr) and
    the identity 2·A·u + b = √disc for the positive root, we get

        dΔt/dΔr = C'(Δr) / (2·Δt·√disc),

    where C'(Δr) = 2·(r_m + Δr). This diverges as Δt → 0 (i.e., Δr → 0);
    callers must guard against evaluation at Δr = 0.
    """
    f = dVdr_m / (2.0 * mu_ij * r_m)
    A = (r_m * f) ** 2
    b = v_m ** 2 - 2.0 * (r_m ** 2) * f
    C_star = delta_r * (2.0 * r_m + delta_r)
    dCdΔr = 2.0 * (r_m + delta_r)

    if A < DENOM_FLOOR:
        if abs(b) < DENOM_FLOOR:
            return 0.0, 0.0
        u = C_star / b
        if u <= 0.0:
            return 0.0, 0.0
        dt = np.sqrt(u)
        return dt, (dCdΔr / b) / (2.0 * dt)

    disc = b * b + 4.0 * A * C_star
    sqrt_disc = np.sqrt(max(disc, 0.0))
    if b >= 0.0:
        u = 2.0 * C_star / (b + sqrt_disc)
    else:
        u = (-b + sqrt_disc) / (2.0 * A)
    if u <= 0.0:
        return 0.0, 0.0
    dt = np.sqrt(u)
    if sqrt_disc < DENOM_FLOOR:
        return dt, 0.0
    return dt, (dCdΔr / sqrt_disc) / (2.0 * dt)


def _find_turning_point_slice(surface, epsilon, b, phi_m_deg):
    """Solve V(r, φ_m) + ε·b²/r² = ε for the outermost r.

    A minimal, dependency-free bracketer: start at a large r where
    Y(r) = ε − V − ε·b²/r² > 0 (well outside the well) and walk inward
    until Y flips sign. Suits Phase-1 needs; the heavier machinery in
    `CrossSectionSolver.find_turning_point` (which handles multi-orbit
    regions via `EC`/`clong`) is not needed here because V̂tilde does not
    itself compute cross sections — it is *fed to* the solver.
    """
    def Y(r):
        return epsilon - surface.V(r, phi_m_deg) - epsilon * b * b / (r * r)

    rmax_data = max(p.rmax for p in surface.potentials)
    rmin_data = min(p.rmin for p in surface.potentials)

    r_hi = 5.0 * rmax_data
    if Y(r_hi) <= 0.0:
        raise RuntimeError(
            f"outer boundary r={r_hi:.3g} has Y <= 0 at (ε={epsilon}, b={b}, "
            f"φ_m={phi_m_deg}°); trajectory does not reach a turning point.")

    r_lo = r_hi
    for _ in range(400):
        r_lo *= 0.9
        if Y(r_lo) < 0.0:
            break
        if r_lo < 0.1 * rmin_data:
            raise RuntimeError(
                f"could not bracket turning point above r={r_lo:.3g}")
    else:
        raise RuntimeError("bracket search did not converge")

    result = root_scalar(Y, bracket=[r_lo, r_hi], method='brentq',
                          xtol=1e-14, rtol=1e-14)
    return float(result.root)


def build_collision_frame(potential_surface, epsilon, b, phi_m,
                          psi, L_mag, alpha_L,
                          mu_ij, mu_i, ell_len):
    """Assemble a `CollisionFrame` at closest approach.

    Parameters
    ----------
    potential_surface : PotentialSurface
    epsilon : float
        Collision energy ε (positive).
    b : float
        Impact parameter (non-negative).
    phi_m : float
        Angle between **r_m** and ℓ̂_m, in radians. Range [0, π].
    psi : float
        Azimuth of ℓ̂_m about ẑ, measured from x̂ = v̂_m, in radians.
    L_mag : float
        Angular-momentum magnitude |**L**|.
    alpha_L : float
        Direction of **L** in the plane ⊥ ℓ̂_m, in radians, measured from
        ê₁ toward ê₂.
    mu_ij : float
        Collision reduced mass μ_ij (used in the radial equation and in
        Δt(Δr)).
    mu_i : float
        Diatom reduced mass μ_i (moment of inertia = μ_i ℓ²).
    ell_len : float
        Diatom bond length ℓ (moment of inertia parameter).
    """
    # 1. Solve for r_m using V(·, φ_m) as a 1D slice.
    phi_m_deg = float(np.degrees(phi_m))
    r_m = _find_turning_point_slice(potential_surface, epsilon, b, phi_m_deg)

    # 2. Kinematic scalars at closest approach.
    g = np.sqrt(2.0 * epsilon / mu_ij)
    v_m = b * g / r_m
    dVdr_at_m = potential_surface.dV_dr(r_m, phi_m_deg)

    # 3. Frame vectors: **r_m** = r_m ẑ, **v_m** = v_m x̂.
    r_m_vec = np.array([0.0, 0.0, r_m])
    v_m_vec = np.array([v_m, 0.0, 0.0])

    # 4. Diatom orientation ℓ̂_m.
    sin_phi_m = np.sin(phi_m)
    cos_phi_m = np.cos(phi_m)
    ell_hat_m = np.array([
        sin_phi_m * np.cos(psi),
        sin_phi_m * np.sin(psi),
        cos_phi_m,
    ])

    # 5. Tangent basis {ê₁, ê₂} at ℓ̂_m.
    if abs(sin_phi_m) > _SIN_PHI_THRESHOLD:
        e1 = (np.array([0.0, 0.0, 1.0]) - cos_phi_m * ell_hat_m) / sin_phi_m
    else:
        # Limit as sin φ_m → 0: ê₁ = -sign(cos φ_m)·(cos ψ, sin ψ, 0).
        # Near φ_m = 0 (cos φ_m = +1): -(cos ψ, sin ψ, 0).
        # Near φ_m = π (cos φ_m = -1): +(cos ψ, sin ψ, 0).
        sign = -1.0 if cos_phi_m >= 0 else 1.0
        e1 = sign * np.array([np.cos(psi), np.sin(psi), 0.0])
    e2 = np.cross(ell_hat_m, e1)
    L_vec = L_mag * (np.cos(alpha_L) * e1 + np.sin(alpha_L) * e2)

    # 6. dφ/dt|ₘ. General formula:
    #    dφ/dt|ₘ = (1/(r_m sin φ_m))·[**v_m**·ℓ̂_m
    #              + (1/(μ_i ℓ²))·**r_m**·(ℓ̂_m × **L**)]
    #    but analytic cancellation gives the regular form
    #    dφ/dt|ₘ = (v_m cos ψ)/r_m - L_mag sin α_L / (μ_i ℓ²)
    #    which is valid for all φ_m (including sin φ_m = 0).
    dphidt_m = (v_m * np.cos(psi) / r_m
                - L_mag * np.sin(alpha_L) / (mu_i * ell_len ** 2))

    return CollisionFrame(
        r_m=r_m,
        v_m=v_m,
        phi_m=phi_m,
        mu_ij=mu_ij,
        dVdr_at_m=dVdr_at_m,
        dphidt_m=dphidt_m,
        r_m_vec=r_m_vec,
        v_m_vec=v_m_vec,
        ell_hat_m=ell_hat_m,
        L_vec=L_vec,
        L_mag=L_mag,
        psi=psi,
        alpha_L=alpha_L,
        mu_i=mu_i,
        ell_len=ell_len,
    )


class Vtilde:
    """`Potential`-shaped closure implementing Eq. `v0_tilde`.

    Exposes __call__, `deriv`, `deriv2`, `veff`, `rmin`, `rmax`,
    `long_range_pow` so it drops into `CrossSectionSolver` in place of
    the 1D `Potential`.

    Parameters
    ----------
    surface : PotentialSurface
    frame : CollisionFrame
    fd_step_fraction : float
        h = fd_step_fraction · r_m is the finite-difference step used
        by the FD derivative paths.
    derivative : {"fd", "analytic"}
        `"fd"` (default): centered finite difference on V̂tilde itself —
        Phase-1 behaviour, always applicable. `"analytic"`: chain-rule
        formula that consumes surface.dV_dr and surface.dV_dphi, with a
        forward-FD fallback in a small window around r_m where the
        analytic dΔt/dΔr diverges as 1/√Δr.
    """

    _SUPPORTED_DERIVATIVES = ("fd", "analytic")

    def __init__(self, surface, frame, fd_step_fraction=1.0e-4,
                 derivative="fd"):
        if derivative not in self._SUPPORTED_DERIVATIVES:
            raise ValueError(
                f"derivative must be one of {self._SUPPORTED_DERIVATIVES}, "
                f"got {derivative!r}")
        self.surface = surface
        self.frame = frame
        self.h = fd_step_fraction * frame.r_m
        self.derivative = derivative
        # Analytic dΔt/dΔr diverges as Δr → 0 (square-root singularity of
        # the trajectory in radial displacement). Use forward FD when
        # Δr < this threshold; forward FD gives the physical right-limit
        # derivative at r_m, which is what CrossSectionSolver.deflection_angle
        # uses in the deflection-integrand normalization.
        self._analytic_fd_threshold = 5.0 * self.h
        # Expose the same shape as `Potential` for CrossSectionSolver.
        self.rmin = min(p.rmin for p in surface.potentials)
        self.rmax = max(p.rmax for p in surface.potentials)
        self.long_range_pow = surface.long_range_pow

    def __call__(self, r):
        f = self.frame
        dt = delta_t_of_delta_r(r - f.r_m, f.r_m, f.v_m, f.dVdr_at_m, f.mu_ij)
        dphi = f.dphidt_m * dt
        phi_m_deg = np.degrees(f.phi_m)
        dphi_deg = np.degrees(dphi)
        phis = np.array([
            phi_m_deg - dphi_deg,
            phi_m_deg,
            phi_m_deg + dphi_deg,
        ])
        vals = self.surface.V_at_r(r, phis)
        return float(np.mean(vals))

    def deriv(self, r):
        if self.derivative == "analytic":
            return self._deriv_analytic(r)
        # "fd": centered finite difference.
        h = self.h
        return (self(r + h) - self(r - h)) / (2.0 * h)

    def _deriv_analytic(self, r):
        """Analytic chain-rule V̂tilde'(r).

        V̂tilde(r) = ⅓ Σ_s V(r, φ_m + s·Δφ(r)) with s ∈ {−1, 0, +1}.
        d/dr [V(r, φ_m + s·Δφ)] = ∂ᵣV + s·(dΔφ/dr)·∂ᵩV. Summing:

            V̂tilde'(r) = ⅓ Σ_s ∂ᵣV(r, φ_m + s·Δφ)
                       + (dΔφ/dr / 3) · [∂ᵩV(r, φ_m+Δφ) - ∂ᵩV(r, φ_m−Δφ)],

        with dΔφ/dr = (dφ/dt|ₘ) · (dΔt/dΔr).
        """
        f = self.frame
        delta_r = r - f.r_m

        # Forward-FD fallback near r_m (where dΔt/dΔr → ∞).
        if delta_r < self._analytic_fd_threshold:
            h = self.h
            return (self(r + h) - self(r)) / h

        dt, ddt_ddr = delta_t_and_deriv_of_delta_r(
            delta_r, f.r_m, f.v_m, f.dVdr_at_m, f.mu_ij)
        dphi = f.dphidt_m * dt
        dphi_ddr = f.dphidt_m * ddt_ddr

        phi_m_deg = np.degrees(f.phi_m)
        dphi_deg = np.degrees(dphi)
        phi_minus = phi_m_deg - dphi_deg
        phi_center = phi_m_deg
        phi_plus = phi_m_deg + dphi_deg

        dVdr_avg = (self.surface.dV_dr(r, phi_minus)
                    + self.surface.dV_dr(r, phi_center)
                    + self.surface.dV_dr(r, phi_plus)) / 3.0

        dVdphi_plus = self.surface.dV_dphi(r, phi_plus)
        dVdphi_minus = self.surface.dV_dphi(r, phi_minus)
        angular_term = (dphi_ddr / 3.0) * (dVdphi_plus - dVdphi_minus)

        return dVdr_avg + angular_term

    def deriv2(self, r):
        # FD-on-FD (independent of derivative mode).
        h = self.h
        return (self(r + h) - 2.0 * self(r) + self(r - h)) / (h * h)

    def veff(self, r):
        return self(r) + r * self.deriv(r) / 2.0


def make_v_tilde(surface, frame, fd_step_fraction=1.0e-4, derivative="fd"):
    """Return a `Vtilde` closure over `surface` and `frame`.

    Convenience factory; equivalent to
    `Vtilde(surface, frame, fd_step_fraction=..., derivative=...)`.
    """
    return Vtilde(surface, frame, fd_step_fraction=fd_step_fraction,
                  derivative=derivative)
