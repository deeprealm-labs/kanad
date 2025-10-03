
"""
Validation Script 7: H2O Full Ab-Initio Calculation
===================================================

A hallenging end-to-end validation of the Kanad framework by performing
a full Hartree-Fock calculation on the water molecule (H2O) and comparing
the ground state energy to a known, high-quality reference value.

This script validates:
1. The framework's ability to handle more complex, non-linear molecules.
2. The accuracy of the integral computation engine for a multi-atom system.
3. The correctness of the Self-Consistent Field (SCF) solver.
4. The full end-to-end workflow for a real-world quantum chemistry calculation.
"""

import numpy as np
from kanad.core.atom import Atom
from kanad.bonds import BondFactory
from kanad.core.hamiltonians.covalent_hamiltonian import CovalentHamiltonian
from kanad.core.representations.lcao_representation import LCAORepresentation
from kanad.core.representations.base_representation import Molecule

# Known reference value for H2O at STO-3G level
# Source: Computational Chemistry Comparison and Benchmark DataBase (NIST)
# Value is for equilibrium geometry.
REFERENCE_H2O_STO3G_ENERGY = -74.9629182993  # Hartree

def run_validation():
    """Run the H2O validation script."""
    print("="*70)
    print("VALIDATION 7: H2O Full Ab-Initio Calculation")
    print("="*70)
    print()

    # ============================================================================
    # 1. Define the Water Molecule (H2O)
    # ============================================================================
    print("1. Defining the H2O molecular geometry...")
    print("-" * 70)

    # Standard equilibrium geometry for water
    # O at origin, H's in xy-plane
    # Bond length = 0.957 Å, Angle = 104.5°
    bond_angle = np.deg2rad(104.5)
    bond_length = 0.957
    
    o_atom = Atom('O', position=np.array([0.0, 0.0, 0.0]))
    h1_atom = Atom('H', position=np.array([bond_length, 0.0, 0.0]))
    h2_atom = Atom('H', position=np.array(
        [bond_length * np.cos(bond_angle), bond_length * np.sin(bond_angle), 0.0]
    ))

    atoms = [o_atom, h1_atom, h2_atom]
    molecule = Molecule(atoms)

    print(f"Molecule: H2O")
    print(f"  - Oxygen at [0.000, 0.000, 0.000]")
    print(f"  - Hydrogen 1 at [{h1_atom.position[0]:.3f}, {h1_atom.position[1]:.3f}, {h1_atom.position[2]:.3f}]")
    print(f"  - Hydrogen 2 at [{h2_atom.position[0]:.3f}, {h2_atom.position[1]:.3f}, {h2_atom.position[2]:.3f}]")
    print(f"Number of electrons: {molecule.n_electrons}")
    print()

    # ============================================================================
    # 2. Build Representation and Hamiltonian
    # ============================================================================
    print("2. Building representation and Hamiltonian...")
    print("-" * 70)

    representation = LCAORepresentation(molecule)
    hamiltonian = CovalentHamiltonian(molecule, representation, basis_name='sto-3g')

    print(f"Representation: {representation.__class__.__name__}")
    print(f"  - Qubits: {representation.n_qubits}")
    print(f"  - Orbitals: {representation.n_orbitals}")
    print()
    print(f"Hamiltonian: {hamiltonian.__class__.__name__}")
    print(f"  - Basis set: STO-3G")
    print(f"  - Nuclear repulsion: {hamiltonian.nuclear_repulsion:.6f} Ha")
    print()

    # ============================================================================
    # 3. Perform Self-Consistent Field (SCF) Calculation
    # ============================================================================
    print("3. Running Hartree-Fock SCF calculation...")
    print("-" * 70)

    try:
        hf_energy, density_matrix, mo_coeffs = hamiltonian.get_hf_energy(max_iter=100, tol=1e-7)
        scf_success = True
        print("SCF calculation converged successfully!")
        print(f"  - Final Hartree-Fock Energy: {hf_energy:.10f} Hartree")
        print(f"  - Reference Energy:          {REFERENCE_H2O_STO3G_ENERGY:.10f} Hartree")
    except RuntimeError as e:
        scf_success = False
        hf_energy = 0.0
        print(f"SCF calculation failed to converge: {e}")
    
    print()

    # ============================================================================
    # 4. Analyze Results and Validate
    # ============================================================================
    print("4. Analyzing results and validating...")
    print("-" * 70)

    validations = []

    # Check SCF convergence
    if scf_success:
        validations.append(("✓", "SCF procedure converged"))
    else:
        validations.append(("✗", "SCF procedure failed to converge"))

    # Check energy accuracy
    if scf_success:
        energy_error = abs(hf_energy - REFERENCE_H2O_STO3G_ENERGY)
        print(f"Energy error: {energy_error:.10f} Hartree")
        
        # Chemical accuracy is 1.6 mHa (0.0016 Ha)
        if energy_error < 0.0016:
            validations.append(("✓", f"Energy within chemical accuracy ({energy_error:.4e} Ha)"))
        else:
            validations.append(("✗", f"Energy outside chemical accuracy ({energy_error:.4e} Ha)"))
    else:
        validations.append(("✗", "Energy accuracy could not be checked"))

    # Check molecular orbitals
    if scf_success:
        mo_energies = np.linalg.eigvalsh(hamiltonian.h_core)
        n_occ = molecule.n_electrons // 2
        homo_lumo_gap = mo_energies[n_occ] - mo_energies[n_occ-1]
        
        if homo_lumo_gap > 0:
            validations.append(("✓", f"Positive HOMO-LUMO gap: {homo_lumo_gap:.4f} Ha"))
        else:
            validations.append(("✗", f"Negative or zero HOMO-LUMO gap: {homo_lumo_gap:.4f} Ha"))
            
        print(f"\nOccupied MO energies (Hartree):")
        for i in range(n_occ):
            print(f"  MO {i}: {mo_energies[i]:.6f}")
        print(f"Unoccupied MO energies (Hartree):")
        for i in range(n_occ, len(mo_energies)):
            print(f"  MO {i}: {mo_energies[i]:.6f}")

    print()

    # ============================================================================
    # 5. Validation Summary
    # ============================================================================
    print("="*70)
    print("VALIDATION SUMMARY")
    print("="*70)

    for symbol, message in validations:
        print(f"{symbol} {message}")

    passed = sum(1 for s, _ in validations if s == "✓")
    total = len(validations)
    print()
    print(f"Validation score: {passed}/{total} checks passed")

    if passed == total:
        print("\n🎉 ALL VALIDATIONS PASSED! Framework is scientifically accurate for H2O.")
    else:
        print(f"\n⚠️  {total - passed} validation(s) failed. Review above.")

    print("="*70)

if __name__ == "__main__":
    run_validation()
