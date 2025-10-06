"""
Density of States (DOS) - Quick API Demonstration

Shows the DOS calculator API without running expensive SCF calculations.
"""

import numpy as np
from kanad.core.atom import Atom
from kanad.core import Lattice
from kanad.core.molecule import Molecule
from kanad.analysis import DOSCalculator


def show_dos_api():
    """Demonstrate DOS calculator API."""
    print("\n" + "="*70)
    print("DENSITY OF STATES (DOS) - API DEMONSTRATION")
    print("="*70)

    print("\nDOS Calculator is already fully implemented in:")
    print("  kanad/analysis/dos_calculator.py (430 lines)")
    print("\nFeatures available:")
    print("  ✓ Total DOS with Gaussian/Lorentzian broadening")
    print("  ✓ Integrated DOS (IDOS)")
    print("  ✓ Projected DOS (PDOS) - framework ready")
    print("  ✓ Fermi energy calculation")
    print("  ✓ Van Hove singularity detection")
    print("  ✓ Band gap from DOS")
    print("  ✓ DOS plotting utilities")

    print("\n" + "="*70)
    print("USAGE EXAMPLES")
    print("="*70)

    print("\n1. Basic Usage:")
    print("""
    from kanad.io import build_crystal
    from kanad.analysis import DOSCalculator

    # Build a crystal
    si = build_crystal('Si', 'diamond', lattice_constant=5.43, k_points=(4,4,4))

    # Solve SCF
    si.solve_scf_pbc()

    # Create DOS calculator
    dos_calc = DOSCalculator(si.hamiltonian)

    # Compute DOS
    dos = dos_calc.compute_dos(
        energy_range=(-15, 15),  # eV
        n_points=1000,           # Energy grid points
        sigma=0.1,               # Broadening in eV
        method='gaussian'        # or 'lorentzian'
    )
    """)

    print("\n2. DOS Results:")
    print("""
    # Result is a dictionary with:
    dos['energies']              # Energy grid (eV)
    dos['dos']                   # DOS values (states/eV)
    dos['idos']                  # Integrated DOS
    dos['fermi_energy']          # Fermi level (eV)
    dos['n_electrons_from_integration']
    dos['n_electrons_actual']
    """)

    print("\n3. Advanced Features:")
    print("""
    # Find Van Hove singularities
    singularities = dos_calc.find_van_hove_singularities(
        energy_range=(-10, 10),
        threshold=0.5  # Minimum DOS value
    )

    # Band gap from DOS
    gap_info = dos_calc.find_band_gap()
    print(f"Band gap: {gap_info['gap']:.2f} eV")
    print(f"Type: {gap_info['type']}")  # 'direct' or 'indirect'

    # Projected DOS (PDOS)
    pdos = dos_calc.compute_pdos(
        atom_indices=[0, 1, 2]  # Specific atoms
    )
    """)

    print("\n4. Plotting:")
    print("""
    # Simple DOS plot
    dos_calc.plot_dos(dos, show_fermi=True, show_gap=True)

    # Combined band structure + DOS
    bands = si.compute_band_structure(k_path, n_bands=8)
    dos_calc.plot_band_structure_with_dos(bands, dos)
    """)

    print("\n5. Different Broadening Methods:")
    print("""
    # Gaussian (default)
    dos_gauss = dos_calc.compute_dos(sigma=0.1, method='gaussian')

    # Lorentzian (heavier tails)
    dos_lorentz = dos_calc.compute_dos(sigma=0.1, method='lorentzian')
    """)

    print("\n" + "="*70)
    print("THEORY")
    print("="*70)

    print("""
Total DOS:
    DOS(E) = Σ_{nk} w_k δ(E - E_nk)

With Gaussian broadening:
    DOS(E) ≈ Σ_{nk} w_k / (σ√2π) exp(-(E - E_nk)² / 2σ²)

With Lorentzian broadening:
    DOS(E) ≈ Σ_{nk} w_k × (γ/π) / ((E - E_nk)² + γ²)

Fermi Energy:
    ∫_{-∞}^{E_F} DOS(E) dE = N_electrons

Van Hove Singularities:
    Points where ∇_k E(k) = 0 → peaks in DOS
    """)

    print("\n" + "="*70)
    print("COMPUTATIONAL NOTES")
    print("="*70)

    print("""
K-point Convergence:
    - Denser k-grids → smoother DOS
    - Metals: Need dense grids (8×8×8 or more)
    - Semiconductors: 4×4×4 often sufficient

Broadening Parameter:
    - σ too small: Noisy, many sharp peaks
    - σ too large: Over-smoothed, lose features
    - Typical: 0.05-0.3 eV depending on system

Energy Range:
    - Include several eV above/below Fermi level
    - For optical properties: -10 to +10 eV typical
    - For total energy: May need wider range

Grid Points:
    - n_points=1000 usually sufficient
    - For high resolution: 2000-5000 points
    """)

    print("\n" + "="*70)
    print("EXAMPLE OUTPUT")
    print("="*70)

    # Create a minimal example (no SCF needed, just show structure)
    print("\nCreating example periodic system...")

    a = 4.0
    lattice = Lattice(
        lattice_vectors=np.array([
            [a, 0, 0],
            [0, 10, 0],
            [0, 0, 10]
        ]),
        pbc=(True, False, False)
    )

    atoms = [Atom('H', position=[0, 0, 0])]
    chain = Molecule(atoms, lattice=lattice, k_points=(8, 1, 1))

    print(f"\nSystem created:")
    print(f"  Type: 1D periodic H chain")
    print(f"  Lattice constant: {a} Å")
    print(f"  K-points: 8")
    print(f"  Periodic: {chain.is_periodic}")

    print(f"\nTo compute DOS, you would:")
    print(f"  1. chain.solve_scf_pbc()")
    print(f"  2. dos_calc = DOSCalculator(chain.hamiltonian)")
    print(f"  3. dos = dos_calc.compute_dos(...)")

    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)

    print("""
The DOSCalculator class is fully implemented and ready to use!

Key Methods:
  • compute_dos()                  - Total DOS
  • compute_pdos()                 - Projected DOS
  • find_band_gap()                - Extract band gap
  • find_van_hove_singularities()  - Peak detection
  • plot_dos()                     - Visualization
  • plot_band_structure_with_dos() - Combined plots

For working examples that run SCF and compute DOS, see:
  • examples/dos_demo.py           - Full DOS demonstrations
  • examples/periodic_systems_demo.py - Complete periodic demo

To run the full DOS demo (requires SCF calculations):
  $ python examples/dos_demo.py
    """)


if __name__ == '__main__':
    show_dos_api()
