"""
Governance Validation: IONIC BONDING (LiH)

Tests ionic governance protocol validation on LiH molecule.
Validates: weak transfer, large EN difference, large Hubbard U, charge transfer.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../..'))

import numpy as np
from kanad.core.molecule import Molecule
from kanad.core.atom import Atom
from kanad.core.representations.second_quantization import SecondQuantizationRepresentation
from kanad.core.hamiltonians.ionic_hamiltonian import IonicHamiltonian


def test_ionic_lih():
    """Test ionic bonding governance on LiH."""
    print("=" * 80)
    print("GOVERNANCE VALIDATION: IONIC BONDING (LiH)")
    print("=" * 80)

    # Create LiH molecule
    print("\n1. Creating LiH molecule...")
    li = Atom('Li', position=np.array([0.0, 0.0, 0.0]))
    h = Atom('H', position=np.array([1.6, 0.0, 0.0]))

    molecule = Molecule(atoms=[li, h], charge=0, spin=0)

    print(f"   ✓ LiH created (bond length: 1.6 Å)")
    print(f"   - Li electronegativity: {li.electronegativity:.2f}")
    print(f"   - H electronegativity: {h.electronegativity:.2f}")
    print(f"   - ΔEN = {abs(h.electronegativity - li.electronegativity):.2f}")

    # Build Hamiltonian with governance
    print("\n2. Building Ionic Hamiltonian with Governance...")
    representation = SecondQuantizationRepresentation(molecule=molecule)

    hamiltonian = IonicHamiltonian(
        molecule=molecule,
        representation=representation,
        use_governance=True
    )

    print(f"   ✓ Hamiltonian built")
    print(f"   - Orbitals: {hamiltonian.n_orbitals}")
    print(f"   - Electrons: {hamiltonian.n_electrons}")
    print(f"   - Nuclear repulsion: {hamiltonian.nuclear_repulsion:.6f} Ha")

    # Validate with governance
    print("\n3. Running Governance Validation...")
    validation = hamiltonian.validate_with_governance()

    print("\n" + "=" * 80)
    print("GOVERNANCE VALIDATION RESULTS")
    print("=" * 80)

    for check in validation['checks']:
        status = "✓" if check['passed'] else "✗"
        print(f"[{status}] {check['name']}")
        print(f"    {check['message']}")

    print("=" * 80)

    # Overall result
    if validation['all_checks_passed']:
        print("\n🎉 ALL GOVERNANCE CHECKS PASSED - IONIC CHARACTER CONFIRMED")
        print("   → LiH exhibits correct ionic bonding physics")
    else:
        print("\n⚠ Some governance checks failed")
        print("   → Review bonding character")

    # Detailed metrics
    print("\nDetailed Metrics:")
    print(f"  Transfer integral: {validation['max_transfer_integral']:.6f} Ha "
          f"({validation['max_transfer_integral'] * 27.211:.2f} eV)")
    print(f"  Energy spread: {validation['energy_spread']:.6f} Ha "
          f"({validation['energy_spread'] * 27.211:.2f} eV)")
    print(f"  Average Hubbard U: {validation['average_hubbard_u']:.6f} Ha "
          f"({validation['average_hubbard_u'] * 27.211:.1f} eV)")

    # SCF analysis
    print("\n4. Classical SCF Analysis...")
    density, energy = hamiltonian.solve_scf(max_iterations=50)

    print(f"   SCF Energy: {energy:.6f} Ha ({energy * 27.211:.4f} eV)")

    # Charge transfer
    charges = hamiltonian.compute_charge_transfer(density)
    print(f"\nCharge Transfer:")
    print(f"  Li: {charges[0]:+.4f} e (expected: +0.8 to +0.9)")
    print(f"  H:  {charges[1]:+.4f} e (expected: -0.8 to -0.9)")

    charge_transfer = abs(charges[0])
    if charge_transfer > 0.5:
        print(f"  ✓ IONIC character confirmed (q = {charge_transfer:.2f} e)")
    else:
        print(f"  ✗ Weak charge transfer (q = {charge_transfer:.2f} e)")

    print("\n" + "=" * 80)
    print("IONIC BONDING VALIDATION COMPLETE")
    print("=" * 80)

    return validation['all_checks_passed']


if __name__ == "__main__":
    success = test_ionic_lih()
    sys.exit(0 if success else 1)
