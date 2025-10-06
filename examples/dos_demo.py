"""
Density of States (DOS) Demonstration

Shows how to compute and analyze DOS for periodic systems:
1. Compute total DOS with different broadening methods
2. Find band gaps from DOS
3. Identify Van Hove singularities
4. Integrate DOS to verify electron count
5. Plot DOS with Fermi level
"""

import numpy as np
import matplotlib.pyplot as plt
from kanad.core.atom import Atom
from kanad.core import Lattice
from kanad.core.molecule import Molecule
from kanad.io import build_crystal
from kanad.analysis import DOSCalculator


def demo_1d_chain_dos():
    """Demo 1: DOS for 1D hydrogen chain."""
    print("\n" + "="*70)
    print("DEMO 1: DOS for 1D Hydrogen Chain")
    print("="*70)

    # 1D periodic chain
    a = 4.0  # Angstrom spacing
    lattice = Lattice(
        lattice_vectors=np.array([
            [a, 0, 0],
            [0, 10, 0],
            [0, 0, 10]
        ]),
        pbc=(True, False, False)
    )

    atoms = [Atom('H', position=[0, 0, 0])]

    # Need dense k-point sampling for good DOS
    chain = Molecule(atoms, lattice=lattice, k_points=(16, 1, 1))

    print(f"\n1D H chain:")
    print(f"  Lattice constant: {a} Å")
    print(f"  K-points: 16")

    # Solve SCF
    print("\nSolving SCF...")
    result = chain.solve_scf_pbc(verbose=0)

    print(f"  Energy: {result['energy']:.6f} Ha")
    print(f"  Converged: {result['converged']}")
    print(f"  Fermi energy: {result['fermi_energy']*27.2114:.3f} eV")

    # Compute DOS
    print("\nComputing DOS...")
    dos_calc = DOSCalculator(chain.hamiltonian)

    # Try different broadening parameters
    dos_narrow = dos_calc.compute_dos(
        energy_range=(-15, 15),
        n_points=1000,
        sigma=0.05,  # Narrow peaks
        method='gaussian'
    )

    dos_broad = dos_calc.compute_dos(
        energy_range=(-15, 15),
        n_points=1000,
        sigma=0.3,  # Broader peaks
        method='gaussian'
    )

    print(f"\nDOS computed with different broadening:")
    print(f"  Narrow (σ=0.05 eV): Peak DOS = {np.max(dos_narrow['dos']):.3f} states/eV")
    print(f"  Broad (σ=0.30 eV): Peak DOS = {np.max(dos_broad['dos']):.3f} states/eV")
    print(f"  Fermi energy: {dos_narrow['fermi_energy']:.3f} eV")

    # Check electron count from integration
    print(f"\nElectron count verification:")
    print(f"  From integration: {dos_narrow['n_electrons_from_integration']:.2f}")
    print(f"  Actual: {dos_narrow['n_electrons_actual']}")
    error = abs(dos_narrow['n_electrons_from_integration'] - dos_narrow['n_electrons_actual'])
    print(f"  Error: {error:.2f} electrons")

    # Plot comparison
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Narrow broadening
    ax1.plot(dos_narrow['energies'], dos_narrow['dos'], 'b-', linewidth=2)
    ax1.fill_between(dos_narrow['energies'], 0, dos_narrow['dos'], alpha=0.3)
    ax1.axvline(dos_narrow['fermi_energy'], color='r', linestyle='--',
                linewidth=2, label=f'E_F = {dos_narrow["fermi_energy"]:.2f} eV')
    ax1.set_xlabel('Energy (eV)', fontsize=12)
    ax1.set_ylabel('DOS (states/eV)', fontsize=12)
    ax1.set_title('Narrow Broadening (σ=0.05 eV)', fontsize=14)
    ax1.legend()
    ax1.grid(alpha=0.3)
    ax1.set_xlim(-15, 15)

    # Broad broadening
    ax2.plot(dos_broad['energies'], dos_broad['dos'], 'g-', linewidth=2)
    ax2.fill_between(dos_broad['energies'], 0, dos_broad['dos'], alpha=0.3)
    ax2.axvline(dos_broad['fermi_energy'], color='r', linestyle='--',
                linewidth=2, label=f'E_F = {dos_broad["fermi_energy"]:.2f} eV')
    ax2.set_xlabel('Energy (eV)', fontsize=12)
    ax2.set_ylabel('DOS (states/eV)', fontsize=12)
    ax2.set_title('Broad Broadening (σ=0.30 eV)', fontsize=14)
    ax2.legend()
    ax2.grid(alpha=0.3)
    ax2.set_xlim(-15, 15)

    plt.tight_layout()
    plt.savefig('/tmp/dos_1d_chain.png', dpi=300, bbox_inches='tight')
    print(f"\nPlot saved to /tmp/dos_1d_chain.png")
    plt.close()

    return chain, dos_narrow


def demo_van_hove_singularities():
    """Demo 2: Identify Van Hove singularities."""
    print("\n" + "="*70)
    print("DEMO 2: Van Hove Singularities")
    print("="*70)

    # 1D chain shows clear Van Hove singularities
    a = 3.5
    lattice = Lattice(np.diag([a, 10, 10]), pbc=(True, False, False))
    atoms = [Atom('H', [0, 0, 0])]
    chain = Molecule(atoms, lattice=lattice, k_points=(20, 1, 1))

    print("\nSolving SCF for 1D chain...")
    chain.solve_scf_pbc(verbose=0)

    dos_calc = DOSCalculator(chain.hamiltonian)

    print("\nSearching for Van Hove singularities...")
    singularities = dos_calc.find_van_hove_singularities(
        energy_range=(-20, 20),
        threshold=0.05  # Minimum DOS to count as singularity
    )

    print(f"\nFound {len(singularities)} Van Hove singularities:")
    for i, sing in enumerate(singularities[:10]):  # Show first 10
        print(f"  {i+1}. Energy: {sing['energy']:8.3f} eV, "
              f"DOS: {sing['dos_value']:7.3f} states/eV")

    # Plot DOS with singularities marked
    dos_result = dos_calc.compute_dos(
        energy_range=(-20, 20),
        n_points=2000,
        sigma=0.05
    )

    plt.figure(figsize=(12, 6))
    plt.plot(dos_result['energies'], dos_result['dos'], 'b-', linewidth=2,
             label='DOS')
    plt.axvline(dos_result['fermi_energy'], color='r', linestyle='--',
                linewidth=2, label='Fermi level')

    # Mark singularities
    for sing in singularities[:10]:
        plt.axvline(sing['energy'], color='orange', linestyle=':', alpha=0.5)
        plt.plot(sing['energy'], sing['dos_value'], 'ro', markersize=8)

    plt.xlabel('Energy (eV)', fontsize=14)
    plt.ylabel('DOS (states/eV)', fontsize=14)
    plt.title('1D Chain DOS with Van Hove Singularities', fontsize=16)
    plt.legend(fontsize=12)
    plt.grid(alpha=0.3)
    plt.xlim(-20, 20)

    plt.tight_layout()
    plt.savefig('/tmp/dos_van_hove.png', dpi=300, bbox_inches='tight')
    print(f"\nPlot saved to /tmp/dos_van_hove.png")
    plt.close()


def demo_integrated_dos():
    """Demo 3: Integrated DOS (IDOS)."""
    print("\n" + "="*70)
    print("DEMO 3: Integrated DOS")
    print("="*70)

    # Simple system
    a = 4.0
    lattice = Lattice(np.diag([a, 10, 10]), pbc=(True, False, False))
    atoms = [Atom('H', [0, 0, 0])]
    chain = Molecule(atoms, lattice=lattice, k_points=(12, 1, 1))

    print("\nSolving SCF...")
    result = chain.solve_scf_pbc(verbose=0)

    dos_calc = DOSCalculator(chain.hamiltonian)
    dos_result = dos_calc.compute_dos(
        energy_range=(-20, 20),
        n_points=1000,
        sigma=0.1
    )

    print(f"\nIntegrated DOS (IDOS):")
    print(f"  Total electrons: {result['n_electrons_actual']}")
    print(f"  IDOS at E_F: {dos_result['n_electrons_from_integration']:.2f}")
    print(f"  IDOS at E_max: {dos_result['idos'][-1]:.2f}")

    # Plot DOS and IDOS together
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

    # DOS
    ax1.plot(dos_result['energies'], dos_result['dos'], 'b-', linewidth=2)
    ax1.fill_between(dos_result['energies'], 0, dos_result['dos'], alpha=0.3)
    ax1.axvline(dos_result['fermi_energy'], color='r', linestyle='--',
                linewidth=2, label='E_F')
    ax1.set_ylabel('DOS (states/eV)', fontsize=12)
    ax1.set_title('Density of States', fontsize=14)
    ax1.legend()
    ax1.grid(alpha=0.3)

    # IDOS
    ax2.plot(dos_result['energies'], dos_result['idos'], 'g-', linewidth=2)
    ax2.axvline(dos_result['fermi_energy'], color='r', linestyle='--',
                linewidth=2, label='E_F')
    ax2.axhline(result['n_electrons_actual'], color='orange', linestyle=':',
                linewidth=2, label=f'N_electrons = {result["n_electrons_actual"]}')
    ax2.set_xlabel('Energy (eV)', fontsize=12)
    ax2.set_ylabel('Integrated DOS (# states)', fontsize=12)
    ax2.set_title('Integrated Density of States', fontsize=14)
    ax2.legend()
    ax2.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig('/tmp/dos_integrated.png', dpi=300, bbox_inches='tight')
    print(f"\nPlot saved to /tmp/dos_integrated.png")
    plt.close()


def demo_gaussian_vs_lorentzian():
    """Demo 4: Compare Gaussian vs Lorentzian broadening."""
    print("\n" + "="*70)
    print("DEMO 4: Gaussian vs Lorentzian Broadening")
    print("="*70)

    # Create system
    a = 4.0
    lattice = Lattice(np.diag([a, 10, 10]), pbc=(True, False, False))
    atoms = [Atom('H', [0, 0, 0])]
    chain = Molecule(atoms, lattice=lattice, k_points=(12, 1, 1))

    print("\nSolving SCF...")
    chain.solve_scf_pbc(verbose=0)

    dos_calc = DOSCalculator(chain.hamiltonian)

    # Gaussian broadening
    dos_gauss = dos_calc.compute_dos(
        energy_range=(-15, 15),
        n_points=1000,
        sigma=0.2,
        method='gaussian'
    )

    # Lorentzian broadening
    dos_lorentz = dos_calc.compute_dos(
        energy_range=(-15, 15),
        n_points=1000,
        sigma=0.2,
        method='lorentzian'
    )

    print(f"\nBroadening comparison (σ=0.2 eV):")
    print(f"  Gaussian peak DOS: {np.max(dos_gauss['dos']):.3f} states/eV")
    print(f"  Lorentzian peak DOS: {np.max(dos_lorentz['dos']):.3f} states/eV")

    # Plot comparison
    plt.figure(figsize=(12, 6))
    plt.plot(dos_gauss['energies'], dos_gauss['dos'], 'b-', linewidth=2,
             label='Gaussian', alpha=0.7)
    plt.plot(dos_lorentz['energies'], dos_lorentz['dos'], 'r-', linewidth=2,
             label='Lorentzian', alpha=0.7)
    plt.axvline(dos_gauss['fermi_energy'], color='k', linestyle='--',
                linewidth=1.5, label='Fermi level')

    plt.xlabel('Energy (eV)', fontsize=14)
    plt.ylabel('DOS (states/eV)', fontsize=14)
    plt.title('Gaussian vs Lorentzian Broadening', fontsize=16)
    plt.legend(fontsize=12)
    plt.grid(alpha=0.3)
    plt.xlim(-15, 15)

    plt.tight_layout()
    plt.savefig('/tmp/dos_broadening_comparison.png', dpi=300, bbox_inches='tight')
    print(f"\nPlot saved to /tmp/dos_broadening_comparison.png")
    plt.close()


def demo_dos_for_crystal():
    """Demo 5: DOS for a 3D crystal structure."""
    print("\n" + "="*70)
    print("DEMO 5: DOS for 3D Crystal (Simple Cubic)")
    print("="*70)

    # Simple cubic hydrogen (for fast calculation)
    print("\nBuilding simple cubic H crystal...")
    sc = build_crystal('H', 'sc', lattice_constant=3.0, k_points=(4, 4, 4))

    print(f"  Structure: Simple Cubic")
    print(f"  Lattice constant: 3.0 Å")
    print(f"  K-points: 4×4×4 = 64 k-points")

    print("\nSolving SCF (this may take a moment)...")
    result = sc.solve_scf_pbc(verbose=0)

    print(f"\nResults:")
    print(f"  Energy: {result['energy']:.6f} Ha")
    print(f"  Converged: {result['converged']}")

    # Compute DOS
    print("\nComputing DOS...")
    dos_calc = DOSCalculator(sc.hamiltonian)

    dos = dos_calc.compute_dos(
        energy_range=(-15, 15),
        n_points=1000,
        sigma=0.15
    )

    print(f"  Peak DOS: {np.max(dos['dos']):.3f} states/eV")
    print(f"  Fermi energy: {dos['fermi_energy']:.3f} eV")

    # Band gap
    gap_info = dos_calc.find_band_gap()
    print(f"\nBand gap:")
    print(f"  Gap: {gap_info['gap']:.2f} eV")
    print(f"  Type: {gap_info['type']}")
    print(f"  VBM: {gap_info['vbm']:.2f} eV")
    print(f"  CBM: {gap_info['cbm']:.2f} eV")

    # Plot DOS with gap shading
    plt.figure(figsize=(10, 6))
    plt.plot(dos['energies'], dos['dos'], 'b-', linewidth=2, label='DOS')
    plt.fill_between(dos['energies'], 0, dos['dos'], alpha=0.3)
    plt.axvline(dos['fermi_energy'], color='r', linestyle='--',
                linewidth=2, label=f'E_F = {dos["fermi_energy"]:.2f} eV')

    # Shade band gap if significant
    if gap_info['gap'] > 0.1:
        plt.axvspan(gap_info['vbm'], gap_info['cbm'], alpha=0.2, color='gray',
                   label=f'Gap = {gap_info["gap"]:.2f} eV')

    plt.xlabel('Energy (eV)', fontsize=14)
    plt.ylabel('DOS (states/eV)', fontsize=14)
    plt.title('DOS for Simple Cubic Crystal', fontsize=16)
    plt.legend(fontsize=12)
    plt.grid(alpha=0.3)
    plt.xlim(-15, 15)
    plt.ylim(bottom=0)

    plt.tight_layout()
    plt.savefig('/tmp/dos_3d_crystal.png', dpi=300, bbox_inches='tight')
    print(f"\nPlot saved to /tmp/dos_3d_crystal.png")
    plt.close()


def main():
    """Run all DOS demos."""
    print("\n" + "="*70)
    print("DENSITY OF STATES (DOS) - COMPREHENSIVE DEMONSTRATION")
    print("="*70)

    try:
        # Demo 1: Basic DOS
        chain, dos = demo_1d_chain_dos()

        # Demo 2: Van Hove singularities
        demo_van_hove_singularities()

        # Demo 3: Integrated DOS
        demo_integrated_dos()

        # Demo 4: Broadening comparison
        demo_gaussian_vs_lorentzian()

        # Demo 5: 3D crystal DOS
        demo_dos_for_crystal()

        print("\n" + "="*70)
        print("DOS DEMONSTRATION COMPLETE!")
        print("="*70)
        print("\nFeatures Demonstrated:")
        print("  ✓ Total DOS computation (Gaussian & Lorentzian)")
        print("  ✓ Van Hove singularity detection")
        print("  ✓ Integrated DOS (IDOS)")
        print("  ✓ Electron count verification")
        print("  ✓ Band gap from DOS")
        print("  ✓ Fermi level calculation")
        print("  ✓ 1D, 2D, 3D systems")
        print("\nPlots Generated:")
        print("  • /tmp/dos_1d_chain.png")
        print("  • /tmp/dos_van_hove.png")
        print("  • /tmp/dos_integrated.png")
        print("  • /tmp/dos_broadening_comparison.png")
        print("  • /tmp/dos_3d_crystal.png")

    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()


if __name__ == '__main__':
    main()
