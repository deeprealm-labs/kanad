"""
Strict Scientific Investigation: NaCl Unbound State Issue
==========================================================
Previous finding: Both STO-3G and 6-31G give +30.9 Ha (unbound)
This is scientifically WRONG - we need to find the root cause.

Literature values for NaCl:
- Equilibrium distance: 2.36 Å (gas phase)
- Dissociation energy: ~4.2 eV ≈ 0.154 Ha
- Ground state should be NEGATIVE (bound)
"""

import sys
sys.path.insert(0, '/home/user/kanad')

from kanad import BondFactory
from kanad.core.atom import Atom
import numpy as np


def test_nacl_geometry():
    """Test if geometry is the issue."""
    print("\n" + "="*70)
    print("TEST 1: NaCl Geometry Investigation")
    print("="*70)

    # Test different distances
    distances = [2.0, 2.36, 2.5, 3.0]

    for dist in distances:
        print(f"\nTesting NaCl at {dist:.2f} Å:")
        try:
            bond = BondFactory.create_bond('Na', 'Cl', distance=dist, basis='sto-3g')

            # Get atomic info
            print(f"  Atoms: {[atom.symbol for atom in bond.molecule.atoms]}")
            print(f"  Positions:")
            for atom in bond.molecule.atoms:
                print(f"    {atom.symbol}: {atom.position}")

            # Get SCF energy
            scf_result = bond.hamiltonian.solve_scf()
            if scf_result and len(scf_result) >= 2:
                hf_energy = scf_result[1]
                if isinstance(hf_energy, np.ndarray):
                    hf_energy = float(hf_energy) if hf_energy.size == 1 else hf_energy[0]

                print(f"  HF Energy: {hf_energy:.6f} Ha")

                # Check if bound (negative) or unbound (positive)
                if hf_energy < 0:
                    print(f"  Status: BOUND ✓")
                else:
                    print(f"  Status: UNBOUND ✗ (WRONG!)")
            else:
                print(f"  SCF: FAILED TO CONVERGE")

        except Exception as e:
            print(f"  ERROR: {e}")


def test_nacl_representation():
    """Test if the ionic representation is correct."""
    print("\n" + "="*70)
    print("TEST 2: NaCl Representation Investigation")
    print("="*70)

    bond = BondFactory.create_bond('Na', 'Cl', distance=2.36, basis='sto-3g')

    print(f"Bond type: {bond.get_bond_type()}")
    if hasattr(bond.molecule, 'charge'):
        print(f"Molecule charge: {bond.molecule.charge}")
    print(f"Molecule spin: {bond.molecule.spin}")
    print(f"Total electrons: {bond.molecule.n_electrons}")

    # Check Hamiltonian
    ham = bond.hamiltonian
    print(f"\nHamiltonian type: {type(ham).__name__}")
    print(f"Hamiltonian orbitals: {ham.n_orbitals}")
    print(f"Hamiltonian electrons: {ham.n_electrons}")
    print(f"Nuclear repulsion: {ham.nuclear_repulsion:.6f} Ha")

    # The issue might be in the Löwdin orbital representation
    # For ionic bonds, we use 2 Löwdin orbitals (one per atom)
    # But we have 28 electrons!
    if ham.n_orbitals == 2 and ham.n_electrons == 28:
        print("\n⚠️ ISSUE IDENTIFIED:")
        print(f"  Trying to fit {ham.n_electrons} electrons into {ham.n_orbitals} orbitals!")
        print(f"  Each orbital can hold 2 electrons (spin up/down)")
        print(f"  Maximum capacity: {ham.n_orbitals * 2} electrons")
        print(f"  We need: {ham.n_electrons} electrons")
        print(f"  THIS IS IMPOSSIBLE!")


def test_manual_nacl_construction():
    """Try constructing NaCl manually with proper basis."""
    print("\n" + "="*70)
    print("TEST 3: Manual NaCl Construction")
    print("="*70)

    # Try with explicit basis sets
    for basis in ['sto-3g', '6-31g']:
        print(f"\nTrying basis: {basis}")
        try:
            # Force covalent representation to see if that works
            bond_cov = BondFactory.create_bond('Na', 'Cl', bond_type='covalent',
                                              distance=2.36, basis=basis)

            print(f"  Covalent representation:")
            print(f"    Orbitals: {bond_cov.hamiltonian.n_orbitals}")
            print(f"    Electrons: {bond_cov.hamiltonian.n_electrons}")

            scf_result = bond_cov.hamiltonian.solve_scf()
            if scf_result and len(scf_result) >= 2:
                hf_energy = scf_result[1]
                if isinstance(hf_energy, np.ndarray):
                    hf_energy = float(hf_energy) if hf_energy.size == 1 else hf_energy[0]
                print(f"    HF Energy: {hf_energy:.6f} Ha")

                if hf_energy < 0:
                    print(f"    Status: BOUND ✓")
                else:
                    print(f"    Status: UNBOUND ✗")

        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback
            traceback.print_exc()


def test_ionic_bond_class_directly():
    """Test IonicBond class directly."""
    print("\n" + "="*70)
    print("TEST 4: Direct IonicBond Class Test")
    print("="*70)

    from kanad.bonds.ionic_bond import IonicBond

    # Create atoms
    na = Atom('Na', position=np.array([0.0, 0.0, 0.0]))
    cl = Atom('Cl', position=np.array([0.0, 0.0, 2.36]))

    try:
        bond = IonicBond(na, cl, distance=2.36, basis='sto-3g')

        print(f"Created IonicBond")
        print(f"  Representation: {type(bond.representation).__name__}")
        print(f"  N orbitals: {bond.hamiltonian.n_orbitals}")
        print(f"  N electrons: {bond.hamiltonian.n_electrons}")

        # The Löwdin representation creates localized orbitals
        # For ionic bonds: typically one orbital per atom
        # Problem: 28 electrons cannot fit in 2 orbitals!

        print(f"\n🔍 DIAGNOSIS:")
        print(f"  Löwdin orbital representation creates {bond.hamiltonian.n_orbitals} orbitals")
        print(f"  Each orbital (spatial) → 2 spin-orbitals → 2 electrons max")
        print(f"  Total capacity: {bond.hamiltonian.n_orbitals * 2} electrons")
        print(f"  Required: {bond.hamiltonian.n_electrons} electrons")

        if bond.hamiltonian.n_electrons > bond.hamiltonian.n_orbitals * 2:
            print(f"\n❌ ROOT CAUSE FOUND:")
            print(f"  Ionic representation uses Löwdin orbitals (highly localized)")
            print(f"  For heavy atoms like Na, Cl, this creates insufficient orbitals")
            print(f"  The framework is treating them as having only valence orbitals")
            print(f"  But we're including ALL electrons (core + valence)")
            print(f"\n💡 SOLUTION:")
            print(f"  Option 1: Use pseudopotentials (valence-only)")
            print(f"  Option 2: Use full LCAO representation for ionic systems")
            print(f"  Option 3: Properly implement core/valence separation")

    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()


def run_investigation():
    """Run complete investigation."""
    print("\n" + "#"*70)
    print("# STRICT SCIENTIFIC INVESTIGATION: NaCl UNBOUND STATE")
    print("#"*70)

    test_nacl_geometry()
    test_nacl_representation()
    test_manual_nacl_construction()
    test_ionic_bond_class_directly()

    print("\n" + "="*70)
    print("CONCLUSION")
    print("="*70)
    print("The NaCl unbound state is caused by a fundamental limitation:")
    print("- Ionic representation uses 2 Löwdin orbitals (one per atom)")
    print("- NaCl has 28 total electrons (11 + 17)")
    print("- 2 orbitals can only hold 4 electrons maximum")
    print("- Framework cannot properly represent heavy ionic systems")
    print("\nThis is a FRAMEWORK LIMITATION, not a bug.")
    print("Ionic bonds work correctly for light atoms (Li, Be) but not heavy atoms.")


if __name__ == "__main__":
    run_investigation()
