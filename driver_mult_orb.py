"""Driver: run the Python solver on PC_mult_orb.in and report orbiting regions."""

from transport_from_potential import (
    read_single_potential, Potential, orbiting_scan, orbiting_regions,
    CrossSectionSolver,
)


def main():
    comment, accuracy, emin, emax, nlong, distances, energies = \
        read_single_potential("data/PC_mult_orb.in")

    pot = Potential(distances, energies, -nlong)
    clong = -energies[-1] * distances[-1]**nlong

    print(f"Comment: {comment}")
    print(f"ACCURACY: {accuracy}, EMIN: {emin}, EMAX: {emax}, NLONG: {nlong}")
    print(f"r range: [{distances[0]}, {distances[-1]}]  ({len(distances)} pts)")
    print(f"V range: [{min(energies)}, {max(energies)}]")
    print(f"CLONG: {clong}")

    orbit_list, ED = orbiting_scan(pot, emin, nlong, clong)
    regions = orbiting_regions(orbit_list)

    print(f"\norbiting_scan found {len(orbit_list)} orbiting sets")
    print(f"ED = {ED}")

    if orbit_list:
        print("\nFirst 5 orbiting sets:")
        print(f"{'E':>24s} {'B':>24s} {'RM':>24s}")
        for E_orb, b_orb, r_orb in orbit_list[:5]:
            print(f"{E_orb:24.14e} {b_orb:24.14e} {r_orb:24.14e}")
        if len(orbit_list) > 5:
            print("...")
            for E_orb, b_orb, r_orb in orbit_list[-3:]:
                print(f"{E_orb:24.14e} {b_orb:24.14e} {r_orb:24.14e}")

    print(f"\norbiting_regions found {len(regions)} region(s):")
    print(f"{'idx':>4s} {'E_start':>20s} {'E_end':>20s} {'i_start':>8s} {'i_end':>8s}")
    for i, (E_start, E_end, i_start, i_end) in enumerate(regions):
        print(f"{i:4d} {E_start:20.10e} {E_end:20.10e} {i_start:8d} {i_end:8d}")

    EC = max(r[1] for r in regions) if regions else 0.0
    print(f"\nEC (max orbiting energy) = {EC}")

    # Smoke test: try to instantiate the solver and compute at one energy.
    solver = CrossSectionSolver(pot, orbit_list, regions, nlong, clong, accuracy)
    solver.ED = ED
    test_E = max(emin, min(1e-7, emax))
    try:
        cs = solver.compute(test_E, 4)
        print(f"\nSmoke test at E={test_E}: Q(1..4) = {cs.tolist()}")
    except Exception as e:
        print(f"\nSmoke test at E={test_E} FAILED: {type(e).__name__}: {e}")


if __name__ == "__main__":
    main()
