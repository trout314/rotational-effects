"""Angle-resolved transport cross sections Q^(ℓ)(ε, L) via V̂tilde.

Extends the single-angle machinery in `transport_from_potential.py` to
angle-dependent potentials by:

1. Wrapping each collision configuration (ε, b, φ_m, ψ, L_mag, α_L)
   in a V̂tilde closure (`v_tilde.make_v_tilde`).
2. Instantiating a fresh `CrossSectionSolver` per (b, orientation) sample
   (Option A of the design decision — reuses the tested `deflection_angle`
   integrator, at the cost of one solver construction per sample).
3. Averaging (1 - cos^ℓ χ) over the classical isotropic ensemble:
   ℓ̂_m uniform on S² and **L** direction uniform on the circle ⊥ ℓ̂_m,
   which is equivalent to sampling cos φ_m ∼ Uniform[−1, 1],
   ψ ∼ Uniform[0, 2π], α_L ∼ Uniform[0, 2π] (Q1 resolution).
4. Integrating over b with the standard 2π·b weight.

Quadrature (Q3):
- cos φ_m: Gauss-Legendre on [−1, 1].
- ψ, α_L: periodic trapezoid on [0, 2π] (spectrally accurate for smooth
  periodic integrands).
- b: Gauss-Legendre on [0, b_max] via linear change of variables.
- Inner χ integrand: scipy.integrate.quad, delegated to
  `CrossSectionSolver.deflection_angle`.

The current version handles the non-orbiting regime only (empty
`orbit_list` and `regions` passed to `CrossSectionSolver`). Orbiting
support in the multi-angle case would require a per-orientation orbit
scan and is deferred.
"""

from __future__ import annotations

import numpy as np

from potential_surface import PotentialSurface
from transport_from_potential import CrossSectionSolver
from v_tilde import build_collision_frame, make_v_tilde


def compute_Q_ell(
    surface: PotentialSurface,
    epsilon: float,
    L_mag: float,
    mu_ij: float,
    mu_i: float,
    ell_len: float,
    ell_max: int = 30,
    n_b: int = 20,
    n_cos_phi: int = 20,
    n_psi: int = 40,
    n_alpha: int = 40,
    b_max: float | None = None,
    accuracy: float = 1.0e-3,
    derivative: str = "analytic",
) -> np.ndarray:
    """Compute orientation-averaged transport cross sections at (ε, L_mag).

    Parameters
    ----------
    surface : PotentialSurface
        2D V(r, φ) surface, built with either method="radial_first" or
        method="legendre".
    epsilon : float
        Collision energy ε (positive).
    L_mag : float
        Angular-momentum magnitude of the diatom, held fixed (single-J
        calculation; rotational ground state ⇒ L_mag = 0).
    mu_ij, mu_i, ell_len : float
        Reduced masses and diatom bond length; see `build_collision_frame`.
    ell_max : int
        Maximum ℓ index. Returns Q^(1), …, Q^(ell_max) in one shot; the
        cost per (b, orientation) sample is dominated by χ, so extracting
        all ell values from the same χ is essentially free.
    n_b, n_cos_phi : int
        Number of Gauss-Legendre nodes for the b-integral and the
        cos(φ_m)-integral respectively.
    n_psi, n_alpha : int
        Number of periodic-trapezoid nodes on [0, 2π] for ψ and α_L.
    b_max : float, optional
        Upper limit of the b-integral. Defaults to 5 · max(potentials.rmax).
    accuracy : float
        Passed to CrossSectionSolver for its inner χ quadrature (epsrel).
    derivative : {"fd", "analytic"}
        Passed to `make_v_tilde`.

    Returns
    -------
    Q : ndarray of shape (ell_max,)
        Q^(1), Q^(2), …, Q^(ell_max) in the same units as
        (length)² · (energy)⁰ (whatever the surface uses).

    Notes
    -----
    The ensemble average is written in the (φ_m, ψ, α_L) parametrization,
    where sampling cos φ_m ∼ Uniform[−1, 1] and ψ ∼ Uniform[0, 2π] gives
    ℓ̂_m uniform on S², and α_L ∼ Uniform[0, 2π] gives **L** uniform on
    the circle ⊥ ℓ̂_m. This is equivalent to sampling **L** uniform on S²
    with ℓ̂_m uniform on the circle ⊥ **L** — the two descriptions are
    related by an SO(3) automorphism.
    """
    rmax_data = max(p.rmax for p in surface.potentials)
    if b_max is None:
        b_max = 5.0 * rmax_data

    # ---- Quadrature nodes and weights --------------------------------
    # cos φ_m ∈ [-1, 1] (Gauss-Legendre)
    x_cos, w_cos = np.polynomial.legendre.leggauss(n_cos_phi)
    # b ∈ [0, b_max] via b = (b_max/2)(x + 1) (Gauss-Legendre on [-1, 1])
    x_b, w_b = np.polynomial.legendre.leggauss(n_b)
    b_nodes = 0.5 * (x_b + 1.0) * b_max
    # ψ, α_L ∈ [0, 2π] (periodic trapezoid)
    psi_nodes = 2.0 * np.pi * np.arange(n_psi) / n_psi
    alpha_nodes = 2.0 * np.pi * np.arange(n_alpha) / n_alpha

    # ---- Overall prefactor -------------------------------------------
    # Q^(ℓ) = 2π ∫ b db · (1/(8π²)) ∫∫∫ d(cos φ_m) dψ dα_L (1 - cos^ℓ χ)
    #       = π · b_max / (2 · N_ψ · N_α) · Σ b_i·w_b_i · w_c_j · Σ_k Σ_l …
    prefactor = np.pi * b_max / (2.0 * n_psi * n_alpha)

    # ---- Accumulator -------------------------------------------------
    Q = np.zeros(ell_max)
    ell_range = np.arange(1, ell_max + 1)

    # ---- Main 4-D loop ----------------------------------------------
    for i_b in range(n_b):
        b = b_nodes[i_b]
        w_b_i = w_b[i_b]
        if b == 0.0:
            # Head-on collision; χ integrand degenerates. Contribution
            # vanishes because of the b·db weight.
            continue
        for i_c in range(n_cos_phi):
            phi_m = float(np.arccos(x_cos[i_c]))
            w_c_j = w_cos[i_c]
            for psi in psi_nodes:
                for alpha_L in alpha_nodes:
                    frame = build_collision_frame(
                        surface, epsilon, b, phi_m,
                        float(psi), L_mag, float(alpha_L),
                        mu_ij, mu_i, ell_len)
                    Vt = make_v_tilde(surface, frame,
                                      derivative=derivative)
                    solver = CrossSectionSolver(
                        Vt, orbit_list=[], regions=[],
                        nlong=4, clong=0.0, accuracy=accuracy)
                    chi = solver.deflection_angle(
                        epsilon, b, naitk=0, BO=None, RO=None, ROMAX=0)

                    cos_chi = np.cos(chi)
                    integrand = 1.0 - cos_chi ** ell_range
                    Q += b * w_b_i * w_c_j * integrand

    Q *= prefactor
    return Q


def compute_Q_ell_grid(
    surface: PotentialSurface,
    epsilon_array,
    L_mag: float,
    mu_ij: float,
    mu_i: float,
    ell_len: float,
    ell_max: int = 30,
    **kwargs,
) -> np.ndarray:
    """Compute Q^(ℓ)(ε, L_mag) on an array of energies.

    Thin loop over `compute_Q_ell`. Kept as a separate function so that
    a future parallel version (multiprocessing / threading) can be added
    without changing either API.

    Returns
    -------
    Q : ndarray of shape (n_energies, ell_max)
        Q[k, ℓ - 1] = Q^(ℓ)(epsilon_array[k], L_mag).
    """
    epsilon_array = np.asarray(epsilon_array, dtype=float)
    out = np.empty((epsilon_array.size, ell_max))
    for k, eps in enumerate(epsilon_array):
        out[k] = compute_Q_ell(
            surface, float(eps), L_mag, mu_ij, mu_i, ell_len,
            ell_max=ell_max, **kwargs)
    return out
