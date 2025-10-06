"""
Crystal Structure Builder - Comprehensive Showcase

Demonstrates all 9+ crystal structure types available in Kanad.
"""

from kanad.io import build_crystal, build_binary_crystal, get_kpath, get_lattice_info


def showcase_elemental_crystals():
    """Showcase elemental crystal structures."""
    print("\n" + "="*70)
    print("ELEMENTAL CRYSTAL STRUCTURES")
    print("="*70)

    structures = [
        ('H', 'sc', 3.0, 'Simple Cubic (Polonium-like)'),
        ('Fe', 'bcc', 2.87, 'Body-Centered Cubic'),
        ('Cu', 'fcc', 3.61, 'Face-Centered Cubic'),
        ('Mg', 'hcp', 3.21, 'Hexagonal Close-Packed'),
        ('Si', 'diamond', 5.43, 'Diamond Structure'),
    ]

    for element, lattice_type, a, description in structures:
        crystal = build_crystal(element, lattice_type, lattice_constant=a)
        info = get_lattice_info(lattice_type)

        print(f"\n{description}:")
        print(f"  Element: {element}")
        print(f"  Lattice type: {lattice_type.upper()}")
        print(f"  Lattice constant: {a} Å")
        print(f"  Atoms in unit cell: {crystal.n_atoms}")
        print(f"  Unit cell volume: {crystal.lattice.volume:.3f} Ų")
        print(f"  Coordination number: {info.get('coordination', 'N/A')}")
        print(f"  Common examples: {', '.join(info.get('examples', ['None']))}")


def showcase_binary_crystals():
    """Showcase binary compound structures."""
    print("\n" + "="*70)
    print("BINARY COMPOUND STRUCTURES")
    print("="*70)

    compounds = [
        ('Na', 'Cl', 'rocksalt', 5.64, 'Rocksalt (NaCl)'),
        ('Ga', 'As', 'zincblende', 5.65, 'Zinc Blende (GaAs)'),
        ('Zn', 'O', 'wurtzite', 3.25, 'Wurtzite (ZnO)'),
        ('Cs', 'Cl', 'cscl', 4.11, 'CsCl Structure'),
    ]

    for elem_a, elem_b, lattice_type, a, description in compounds:
        crystal = build_binary_crystal(elem_a, elem_b, lattice_type,
                                      lattice_constant=a)
        info = get_lattice_info(lattice_type)

        print(f"\n{description}:")
        print(f"  Composition: {elem_a}{elem_b}")
        print(f"  Lattice type: {lattice_type.upper()}")
        print(f"  Lattice constant: {a} Å")
        print(f"  Atoms in unit cell: {crystal.n_atoms}")
        print(f"  Unit cell volume: {crystal.lattice.volume:.3f} Ų")
        print(f"  Coordination number: {info.get('coordination', 'N/A')}")
        print(f"  Common examples: {', '.join(info.get('examples', ['None']))}")


def showcase_supercells():
    """Demonstrate supercell expansion."""
    print("\n" + "="*70)
    print("SUPERCELL EXPANSION")
    print("="*70)

    # Start with silicon unit cell
    si_unit = build_crystal('Si', 'diamond', lattice_constant=5.43)

    print(f"\nSilicon Diamond Unit Cell:")
    print(f"  Atoms: {si_unit.n_atoms}")
    print(f"  Volume: {si_unit.lattice.volume:.3f} Ų")

    # Create different supercells
    supercells = [
        ((2, 2, 2), '2×2×2'),
        ((3, 1, 1), '3×1×1'),
        ((2, 2, 1), '2×2×1 (slab)'),
    ]

    for size, label in supercells:
        supercell = si_unit.make_supercell(size)
        nx, ny, nz = size
        expected_atoms = si_unit.n_atoms * nx * ny * nz
        expected_volume = si_unit.lattice.volume * nx * ny * nz

        print(f"\n{label} Supercell:")
        print(f"  Atoms: {supercell.n_atoms} (expected: {expected_atoms})")
        print(f"  Volume: {supercell.lattice.volume:.3f} Ų (expected: {expected_volume:.3f})")
        print(f"  Lattice vectors scaled correctly: ", end='')

        # Check lattice vectors
        ratio_a = supercell.lattice.a / si_unit.lattice.a
        ratio_b = supercell.lattice.b / si_unit.lattice.b
        ratio_c = supercell.lattice.c / si_unit.lattice.c
        print(f"a×{ratio_a:.1f}, b×{ratio_b:.1f}, c×{ratio_c:.1f}")


def showcase_kpaths():
    """Show high-symmetry k-paths for different lattices."""
    print("\n" + "="*70)
    print("HIGH-SYMMETRY k-PATHS FOR BAND STRUCTURE")
    print("="*70)

    lattice_types = ['fcc', 'bcc', 'sc', 'hexagonal']

    for lattice_type in lattice_types:
        try:
            k_path, labels, positions = get_kpath(lattice_type, n_points=50)

            print(f"\n{lattice_type.upper()} Lattice:")
            print(f"  k-path: {' → '.join(labels)}")
            print(f"  Total k-points: {len(k_path)}")
            print(f"  Label positions: {positions}")
        except Exception as e:
            print(f"\n{lattice_type.upper()} Lattice: {e}")


def complete_workflow_example():
    """Show complete workflow from structure to calculation."""
    print("\n" + "="*70)
    print("COMPLETE WORKFLOW EXAMPLE: Gallium Arsenide")
    print("="*70)

    print("\nStep 1: Build crystal structure")
    gaas = build_binary_crystal('Ga', 'As', 'zincblende',
                                lattice_constant=5.65,
                                k_points=(4, 4, 4))

    print(f"  ✓ Created GaAs (zinc blende)")
    print(f"  ✓ Atoms in unit cell: {gaas.n_atoms}")
    print(f"  ✓ K-points: 4×4×4")
    print(f"  ✓ Is periodic: {gaas.is_periodic}")

    print("\nStep 2: Check lattice properties")
    print(f"  Lattice vectors:")
    for i, vec in enumerate(gaas.lattice.lattice_vectors):
        print(f"    a{i+1} = [{vec[0]:.3f}, {vec[1]:.3f}, {vec[2]:.3f}] Å")
    print(f"  Volume: {gaas.lattice.volume:.3f} Ų")
    print(f"  PBC: {gaas.lattice.pbc}")

    print("\nStep 3: Get k-path for band structure")
    k_path, labels, positions = get_kpath('fcc', n_points=50)
    print(f"  ✓ k-path: {' → '.join(labels)}")
    print(f"  ✓ Number of k-points: {len(k_path)}")

    print("\nStep 4: What you would do next...")
    print("""
    # Solve SCF
    result = gaas.solve_scf_pbc()

    # Compute band structure
    bands = gaas.compute_band_structure(k_path, n_bands=8)

    # Compute DOS
    from kanad.analysis import DOSCalculator
    dos_calc = DOSCalculator(gaas.hamiltonian)
    dos = dos_calc.compute_dos(energy_range=(-15, 15))

    # Get band gap
    gap = gaas.get_band_gap()
    print(f"Band gap: {gap['gap']:.2f} eV")
    """)


def main():
    """Run all showcases."""
    print("\n" + "="*70)
    print("KANAD CRYSTAL STRUCTURE BUILDER - COMPREHENSIVE SHOWCASE")
    print("="*70)

    showcase_elemental_crystals()
    showcase_binary_crystals()
    showcase_supercells()
    showcase_kpaths()
    complete_workflow_example()

    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)

    print("""
Elemental Structures (5):
  ✓ Simple Cubic (sc)
  ✓ Body-Centered Cubic (bcc)
  ✓ Face-Centered Cubic (fcc)
  ✓ Hexagonal Close-Packed (hcp)
  ✓ Diamond

Binary Compound Structures (5):
  ✓ Rocksalt (NaCl)
  ✓ Zinc Blende (GaAs)
  ✓ Wurtzite (ZnO)
  ✓ CsCl
  ✓ Fluorite (CaF₂)

Additional Features:
  ✓ Supercell expansion (any size)
  ✓ High-symmetry k-paths (FCC, BCC, SC, Hexagonal)
  ✓ Lattice information database
  ✓ Automatic k-point setup
  ✓ Seamless integration with periodic calculations

Total Structures: 10+
Total Lines of Code: 580 lines (kanad/io/crystal_builder.py)

For actual calculations, see:
  • examples/periodic_systems_demo.py
  • examples/dos_demo.py
    """)


if __name__ == '__main__':
    main()
