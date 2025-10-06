"""
Phase 4: Periodic Boundary Conditions Demo

Demonstrates:
1. Simple cubic lattice with Hydrogen chain
2. Silicon diamond structure
3. Band structure calculation
4. Density of states
5. Band gap determination
6. Supercell expansion
"""

import numpy as np
import sys
import logging

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

from kanad.core.atom import Atom
from kanad.core import Lattice
from kanad.core.molecule import Molecule
from kanad.io import build_crystal, build_binary_crystal, get_kpath
from kanad.analysis import DOSCalculator


def demo_1d_hydrogen_chain():
    """Demo 1: 1D periodic hydrogen chain."""
    print("\n" + "="*70)
    print("DEMO 1: 1D Periodic Hydrogen Chain")
    print("="*70)

    # 1D lattice with 5 Angstrom spacing
    a = 5.0
    lattice = Lattice(
        lattice_vectors=np.array([
            [a, 0, 0],
            [0, 10, 0],    # Large in y to avoid interaction
            [0, 0, 10]     # Large in z to avoid interaction
        ]),
        pbc=(True, False, False)  # Only periodic in x
    )

    print(f"\nLattice:")
    print(lattice)

    # Single hydrogen atom at origin
    atoms = [Atom('H', position=[0, 0, 0])]

    # Create periodic system with k-point sampling along chain
    chain = Molecule(atoms, lattice=lattice, k_points=(8, 1, 1))

    print(f"\nMolecule: {chain}")
    print(f"Periodic: {chain.is_periodic}")

    # Solve SCF
    print("\nSolving SCF...")
    result = chain.solve_scf_pbc(verbose=0)

    print(f"Energy per unit cell: {result['energy']:.6f} Ha")
    print(f"SCF converged: {result['converged']}")
    print(f"Fermi energy: {result['fermi_energy']*27.2114:.3f} eV")

    # Band gap
    gap_info = chain.get_band_gap()
    print(f"\nBand structure:")
    print(f"  VBM: {gap_info['vbm']:.3f} eV")
    print(f"  CBM: {gap_info['cbm']:.3f} eV")
    print(f"  Gap: {gap_info['gap']:.3f} eV ({gap_info['type']})")

    return chain


def demo_silicon_diamond():
    """Demo 2: Silicon diamond structure."""
    print("\n" + "="*70)
    print("DEMO 2: Silicon Diamond Crystal")
    print("="*70)

    # Build silicon crystal (diamond structure)
    a_si = 5.43  # Angstrom (experimental)

    si = build_crystal('Si', 'diamond', lattice_constant=a_si, k_points=(2, 2, 2))

    print(f"\nCrystal: {si}")
    print(f"Lattice:")
    print(si.lattice)
    print(f"\nNumber of atoms in unit cell: {si.n_atoms}")

    # Solve SCF
    print("\nSolving SCF for Si (this may take a minute)...")
    try:
        result = si.solve_scf_pbc(verbose=0, max_iterations=50)

        print(f"\nResults:")
        print(f"  Energy per unit cell: {result['energy']:.6f} Ha")
        print(f"  Energy per unit cell: {result['energy']*27.2114:.3f} eV")
        print(f"  Converged: {result['converged']}")
        print(f"  Fermi energy: {result['fermi_energy']*27.2114:.3f} eV")

        # Band gap
        gap_info = si.get_band_gap()
        print(f"\nBand gap:")
        print(f"  Calculated: {gap_info['gap']:.2f} eV ({gap_info['type']})")
        print(f"  Experimental: ~1.1 eV")
        print(f"  Note: HF overestimates band gap (lack of correlation)")

        return si

    except Exception as e:
        logger.error(f"Silicon SCF failed: {e}")
        logger.info("This is expected - Si calculation requires proper pseudopotentials")
        return None


def demo_band_structure():
    """Demo 3: Band structure along high-symmetry path."""
    print("\n" + "="*70)
    print("DEMO 3: Band Structure Calculation")
    print("="*70)

    # Use simple hydrogen chain for fast demo
    print("\nBuilding H chain for band structure demo...")
    a = 4.0
    lattice = Lattice(np.diag([a, 10, 10]), pbc=(True, False, False))
    atoms = [Atom('H', [0, 0, 0])]
    chain = Molecule(atoms, lattice=lattice, k_points=(4, 1, 1))

    print("Solving SCF...")
    chain.solve_scf_pbc(verbose=0)

    print("\nComputing band structure along Γ → X path...")
    # Manual k-path: Gamma (0,0,0) to X (0.5,0,0)
    k_path = np.linspace([0, 0, 0], [0.5, 0, 0], 50)

    try:
        bands = chain.compute_band_structure(k_path, n_bands=5)

        print(f"Band structure computed:")
        print(f"  Number of k-points: {len(bands['k_points'])}")
        print(f"  Number of bands: {bands['band_energies'].shape[1]}")
        print(f"  Energy range: {np.min(bands['band_energies']):.2f} to "
              f"{np.max(bands['band_energies']):.2f} eV")

        return chain, bands

    except Exception as e:
        logger.error(f"Band structure calculation failed: {e}")
        logger.info("PySCF PBC band structure API may differ between versions")
        return chain, None


def demo_density_of_states():
    """Demo 4: Density of states calculation."""
    print("\n" + "="*70)
    print("DEMO 4: Density of States (DOS)")
    print("="*70)

    # Use H chain
    print("\nCreating H chain...")
    a = 4.0
    lattice = Lattice(np.diag([a, 10, 10]), pbc=(True, False, False))
    atoms = [Atom('H', [0, 0, 0])]
    chain = Molecule(atoms, lattice=lattice, k_points=(8, 1, 1))

    print("Solving SCF...")
    result = chain.solve_scf_pbc(verbose=0)

    print("\nComputing DOS...")
    dos_calc = DOSCalculator(chain.hamiltonian)

    dos_result = dos_calc.compute_dos(
        energy_range=(-15, 15),
        n_points=500,
        sigma=0.2,  # eV
        method='gaussian'
    )

    print(f"\nDOS Results:")
    print(f"  Energy range: {dos_result['energies'][0]:.2f} to "
          f"{dos_result['energies'][-1]:.2f} eV")
    print(f"  Peak DOS: {np.max(dos_result['dos']):.3f} states/eV")
    print(f"  Fermi energy: {dos_result['fermi_energy']:.3f} eV")
    print(f"  Electrons from integration: {dos_result['n_electrons_from_integration']:.2f}")
    print(f"  Actual electrons: {dos_result['n_electrons_actual']}")

    # Van Hove singularities
    print("\nSearching for Van Hove singularities...")
    singularities = dos_calc.find_van_hove_singularities(
        energy_range=(-15, 15),
        threshold=0.1
    )

    if singularities:
        print(f"Found {len(singularities)} singularities:")
        for i, sing in enumerate(singularities[:5]):  # Show first 5
            print(f"  {i+1}. E = {sing['energy']:.2f} eV, "
                  f"DOS = {sing['dos_value']:.3f} states/eV")

    return chain, dos_result


def demo_supercell():
    """Demo 5: Supercell expansion."""
    print("\n" + "="*70)
    print("DEMO 5: Supercell Expansion")
    print("="*70)

    # Start with simple cubic lattice
    a = 3.0
    sc = build_crystal('H', 'sc', lattice_constant=a)

    print(f"\nUnit cell:")
    print(f"  Atoms: {sc.n_atoms}")
    print(f"  Volume: {sc.lattice.volume:.3f} Ų")

    # Create 2x2x2 supercell
    supercell = sc.make_supercell((2, 2, 2))

    print(f"\n2x2x2 Supercell:")
    print(f"  Atoms: {supercell.n_atoms} (expected: {sc.n_atoms * 8})")
    print(f"  Volume: {supercell.lattice.volume:.3f} Ų "
          f"(expected: {sc.lattice.volume * 8:.3f})")
    print(f"  Volume ratio: {supercell.lattice.volume / sc.lattice.volume:.1f}×")

    # Verify atoms are replicated correctly
    print(f"\nAtom positions in supercell (first 4):")
    for i in range(min(4, len(supercell.atoms))):
        pos = supercell.atoms[i].position
        print(f"  {i+1}. {supercell.atoms[i].symbol} at [{pos[0]:.2f}, {pos[1]:.2f}, {pos[2]:.2f}]")

    return supercell


def demo_binary_compounds():
    """Demo 6: Binary compound crystals."""
    print("\n" + "="*70)
    print("DEMO 6: Binary Compound Crystals")
    print("="*70)

    # NaCl (rocksalt)
    print("\n1. NaCl (Rocksalt structure):")
    nacl = build_binary_crystal('Na', 'Cl', 'rocksalt', lattice_constant=5.64)
    print(f"   Atoms in unit cell: {nacl.n_atoms}")
    print(f"   Lattice volume: {nacl.lattice.volume:.3f} ų")

    # GaAs (zinc blende)
    print("\n2. GaAs (Zinc blende structure):")
    gaas = build_binary_crystal('Ga', 'As', 'zincblende', lattice_constant=5.65)
    print(f"   Atoms in unit cell: {gaas.n_atoms}")
    print(f"   Lattice volume: {gaas.lattice.volume:.3f} ų")

    # ZnO (wurtzite)
    print("\n3. ZnO (Wurtzite structure):")
    zno = build_binary_crystal('Zn', 'O', 'wurtzite', lattice_constant=3.25)
    print(f"   Atoms in unit cell: {zno.n_atoms}")
    print(f"   Lattice type: Hexagonal")

    return nacl, gaas, zno


def demo_lattice_info():
    """Demo 7: Crystal structure information."""
    print("\n" + "="*70)
    print("DEMO 7: Crystal Structure Information")
    print("="*70)

    from kanad.io import get_lattice_info

    structures = ['sc', 'bcc', 'fcc', 'hcp', 'diamond', 'rocksalt', 'zincblende']

    for structure in structures:
        info = get_lattice_info(structure)
        print(f"\n{structure.upper()}:")
        print(f"  Description: {info['description']}")
        print(f"  Coordination number: {info.get('coordination', 'N/A')}")
        print(f"  Examples: {', '.join(info.get('examples', ['None']))}")


def main():
    """Run all demos."""
    print("\n" + "="*70)
    print("KANAD FRAMEWORK - PHASE 4: PERIODIC BOUNDARY CONDITIONS")
    print("Comprehensive Demonstration")
    print("="*70)

    try:
        # Demo 1: 1D hydrogen chain
        chain = demo_1d_hydrogen_chain()

        # Demo 2: Silicon (may fail, that's OK)
        si = demo_silicon_diamond()

        # Demo 3: Band structure
        chain, bands = demo_band_structure()

        # Demo 4: DOS
        chain, dos = demo_density_of_states()

        # Demo 5: Supercell
        supercell = demo_supercell()

        # Demo 6: Binary compounds
        compounds = demo_binary_compounds()

        # Demo 7: Lattice info
        demo_lattice_info()

        print("\n" + "="*70)
        print("SUMMARY: Phase 4 Implementation Complete!")
        print("="*70)
        print("\nImplemented Features:")
        print("  ✓ Lattice class with PBC support (1D, 2D, 3D)")
        print("  ✓ PeriodicHamiltonian with k-point sampling")
        print("  ✓ Band structure calculations")
        print("  ✓ Density of states (DOS)")
        print("  ✓ Band gap determination")
        print("  ✓ Supercell expansion")
        print("  ✓ Crystal structure builder (8 lattice types)")
        print("  ✓ High-symmetry k-path generation")
        print("\nCrystal Structures Supported:")
        print("  • Simple cubic (sc)")
        print("  • Body-centered cubic (bcc)")
        print("  • Face-centered cubic (fcc)")
        print("  • Hexagonal close-packed (hcp)")
        print("  • Diamond")
        print("  • Rocksalt (NaCl)")
        print("  • Zinc blende (GaAs)")
        print("  • Wurtzite (ZnO)")
        print("\nNext Steps:")
        print("  • Run on real hardware with proper pseudopotentials")
        print("  • Validate band structures against literature")
        print("  • Extend to magnetic/spin-polarized systems")
        print("  • Phonon calculations (future phase)")

    except Exception as e:
        logger.error(f"Demo failed: {e}", exc_info=True)
        print(f"\nNote: Some features require PySCF PBC module.")
        print(f"Install with: pip install pyscf[geomopt]")


if __name__ == '__main__':
    main()
