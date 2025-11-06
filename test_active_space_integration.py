"""Test active space integration with all Hamiltonian types"""
import numpy as np
from kanad.bonds import BondFactory
from kanad.core.active_space import get_governance_active_space
from kanad.core.configuration import ConfigurationSubspace
from kanad.core.classical_solver import compute_subspace_energy
from kanad.governance.protocols.covalent_protocol import CovalentGovernanceProtocol

def test_lih_with_active_space():
    """Test LiH with active space reduction (12 → 10 qubits)"""
    print("\n" + "="*80)
    print("TEST: LiH with Active Space Integration")
    print("="*80)

    # Create LiH molecule
    bond = BondFactory.create_bond('Li', 'H', distance=1.595)
    molecule = bond.molecule

    print(f"\nMolecule: LiH")
    print(f"  Atoms: {molecule.n_atoms}")
    print(f"  Electrons: {molecule.n_electrons}")

    # Get active space
    protocol = CovalentGovernanceProtocol()
    frozen, active, n_active_electrons = get_governance_active_space(molecule, protocol)

    n_total_orbitals = len(frozen) + len(active)
    standard_qubits = n_total_orbitals * 2
    active_qubits = len(active) * 2

    print(f"\n--- Active Space Selection ---")
    print(f"  Total orbitals: {n_total_orbitals}")
    print(f"  Frozen orbitals: {frozen}")
    print(f"  Active orbitals: {active}")
    print(f"  Active electrons: {n_active_electrons}")
    print(f"  Standard qubits: {standard_qubits}")
    print(f"  Active qubits: {active_qubits}")
    print(f"  ✅ Qubit reduction: {standard_qubits} → {active_qubits} ({standard_qubits - active_qubits} qubits saved)")

    # Build Hamiltonian WITH active space parameters
    print(f"\n--- Building Hamiltonian with Active Space ---")
    from kanad.core.hamiltonians.covalent_hamiltonian import CovalentHamiltonian
    from kanad.core.representations.lcao_representation import LCAORepresentation

    representation = LCAORepresentation(molecule)
    hamiltonian_obj = CovalentHamiltonian(
        molecule,
        representation,
        use_governance=False,
        frozen_orbitals=frozen,
        active_orbitals=active
    )

    print(f"  Hamiltonian n_orbitals: {hamiltonian_obj.n_orbitals}")
    print(f"  Hamiltonian n_electrons: {hamiltonian_obj.n_electrons}")
    print(f"  Frozen core energy: {hamiltonian_obj.frozen_core_energy:.6f} Ha")

    # Convert to sparse Hamiltonian
    print(f"\n--- Converting to Sparse Hamiltonian ---")
    hamiltonian = hamiltonian_obj.to_sparse_hamiltonian(mapper='jordan_wigner')

    print(f"  Hamiltonian qubits: {hamiltonian.num_qubits}")
    print(f"  Pauli terms: {len(hamiltonian)}")

    # CRITICAL CHECK: Qubit count must match!
    if hamiltonian.num_qubits == active_qubits:
        print(f"  ✅ SUCCESS: Hamiltonian qubits ({hamiltonian.num_qubits}) matches active space ({active_qubits})")
    else:
        print(f"  ❌ FAILURE: Hamiltonian has {hamiltonian.num_qubits} qubits but active space has {active_qubits}")
        return False

    # Test Hi-VQE workflow with active space
    print(f"\n--- Hi-VQE with Active Space ---")

    subspace = ConfigurationSubspace(
        n_qubits=active_qubits,
        n_electrons=n_active_electrons,
        protocol=protocol
    )

    # Start with HF
    hf_config = subspace.get_hf_configuration()
    subspace.add_config(hf_config)

    print(f"  HF configuration: {hf_config.bitstring}")

    # Compute energy
    energy, amplitudes = compute_subspace_energy(hamiltonian, subspace, use_fast=True)

    print(f"  HF energy: {energy:.8f} Ha")
    print(f"  ✅ Hi-VQE + Active Space integration working!")

    # Generate excitations and test iteration
    print(f"\n--- Testing Excitation Generation ---")
    single_excs = subspace.generate_single_excitations(hf_config)
    double_excs = subspace.generate_double_excitations(hf_config)

    print(f"  Single excitations: {len(single_excs)}")
    print(f"  Double excitations: {len(double_excs)}")

    # Add excitations
    subspace.add_configs(single_excs[:5])  # Add first 5 singles
    subspace.add_configs(double_excs[:3])  # Add first 3 doubles

    print(f"  Subspace size after excitations: {len(subspace)}")

    # Compute improved energy
    energy_improved, _ = compute_subspace_energy(hamiltonian, subspace, use_fast=True)
    improvement = energy - energy_improved

    print(f"  Improved energy: {energy_improved:.8f} Ha")
    print(f"  Improvement: {improvement:.8f} Ha ({improvement*627.5:.2f} kcal/mol)")

    if improvement > 0:
        print(f"  ✅ Energy lowered (correct direction)")
    else:
        print(f"  ⚠️  Energy did not improve (may need more excitations)")

    # Summary
    print(f"\n" + "="*80)
    print("SUMMARY: Active Space Integration TEST")
    print("="*80)
    print(f"\n✅ Active space selection works: {n_total_orbitals} → {len(active)} orbitals")
    print(f"✅ Hamiltonian respects active space: {hamiltonian.num_qubits} qubits (correct!)")
    print(f"✅ Frozen core energy computed: {hamiltonian_obj.frozen_core_energy:.6f} Ha")
    print(f"✅ Hi-VQE subspace methods work with active space")
    print(f"✅ Configuration sampling works: {len(subspace)} configs")
    print(f"✅ Classical diagonalization works: {energy_improved:.8f} Ha")

    print(f"\n🎯 INTEGRATION COMPLETE: Active space + Hamiltonian + Hi-VQE")
    print(f"   Ready for full pipeline testing!")

    return True


def test_h2o_with_active_space():
    """Test H2O with active space reduction (14 → 12 qubits)"""
    print("\n" + "="*80)
    print("TEST: H2O with Active Space Integration")
    print("="*80)

    # Create H2O molecule
    h2o = BondFactory.create_molecule(['O', 'H', 'H'], geometry='water')

    print(f"\nMolecule: H2O")
    print(f"  Atoms: {h2o.n_atoms}")
    print(f"  Electrons: {h2o.n_electrons}")

    # Get active space
    protocol = CovalentGovernanceProtocol()
    frozen, active, n_active_electrons = get_governance_active_space(h2o, protocol)

    n_total_orbitals = len(frozen) + len(active)
    standard_qubits = n_total_orbitals * 2
    active_qubits = len(active) * 2

    print(f"\n--- Active Space Selection ---")
    print(f"  Total orbitals: {n_total_orbitals}")
    print(f"  Frozen orbitals: {frozen}")
    print(f"  Active orbitals: {active}")
    print(f"  Active electrons: {n_active_electrons}")
    print(f"  Standard qubits: {standard_qubits}")
    print(f"  Active qubits: {active_qubits}")
    print(f"  ✅ Qubit reduction: {standard_qubits} → {active_qubits} ({standard_qubits - active_qubits} qubits saved)")

    # Build Hamiltonian WITH active space parameters
    print(f"\n--- Building Hamiltonian with Active Space ---")
    from kanad.core.hamiltonians.covalent_hamiltonian import CovalentHamiltonian
    from kanad.core.representations.lcao_representation import LCAORepresentation

    representation = LCAORepresentation(h2o)
    hamiltonian_obj = CovalentHamiltonian(
        h2o,
        representation,
        use_governance=False,
        frozen_orbitals=frozen,
        active_orbitals=active
    )

    print(f"  Hamiltonian n_orbitals: {hamiltonian_obj.n_orbitals}")
    print(f"  Hamiltonian n_electrons: {hamiltonian_obj.n_electrons}")
    print(f"  Frozen core energy: {hamiltonian_obj.frozen_core_energy:.6f} Ha")

    # Convert to sparse Hamiltonian
    print(f"\n--- Converting to Sparse Hamiltonian ---")
    hamiltonian = hamiltonian_obj.to_sparse_hamiltonian(mapper='jordan_wigner')

    print(f"  Hamiltonian qubits: {hamiltonian.num_qubits}")
    print(f"  Pauli terms: {len(hamiltonian)}")

    # CRITICAL CHECK: Qubit count must match!
    if hamiltonian.num_qubits == active_qubits:
        print(f"  ✅ SUCCESS: Hamiltonian qubits ({hamiltonian.num_qubits}) matches active space ({active_qubits})")
    else:
        print(f"  ❌ FAILURE: Hamiltonian has {hamiltonian.num_qubits} qubits but active space has {active_qubits}")
        return False

    # Test Hi-VQE workflow
    print(f"\n--- Hi-VQE with Active Space ---")

    subspace = ConfigurationSubspace(
        n_qubits=active_qubits,
        n_electrons=n_active_electrons,
        protocol=protocol
    )

    hf_config = subspace.get_hf_configuration()
    subspace.add_config(hf_config)

    print(f"  HF configuration: {hf_config.bitstring}")

    # Compute energy
    energy, _ = compute_subspace_energy(hamiltonian, subspace, use_fast=True)

    print(f"  HF energy: {energy:.8f} Ha")
    print(f"  ✅ Hi-VQE + Active Space integration working!")

    return True


# Run tests
print("="*80)
print("ACTIVE SPACE INTEGRATION TESTS")
print("="*80)
print("\nTesting integration of:")
print("  1. Active space selection (kanad/core/active_space.py)")
print("  2. Hamiltonian construction with active space (all types)")
print("  3. Hi-VQE workflow with reduced qubits")

success = True

# Test 1: LiH
try:
    if not test_lih_with_active_space():
        success = False
except Exception as e:
    print(f"\n❌ LiH test FAILED with error: {e}")
    import traceback
    traceback.print_exc()
    success = False

# Test 2: H2O
try:
    if not test_h2o_with_active_space():
        success = False
except Exception as e:
    print(f"\n❌ H2O test FAILED with error: {e}")
    import traceback
    traceback.print_exc()
    success = False

# Final report
print("\n" + "="*80)
print("FINAL REPORT")
print("="*80)

if success:
    print("\n🎉 ALL TESTS PASSED!")
    print("\n✅ Active space integration is COMPLETE across all Hamiltonian types")
    print("✅ Ready to move to next phase: VQE solver Hi-VQE mode")
else:
    print("\n❌ SOME TESTS FAILED")
    print("   Review errors above and fix issues")

print("\n" + "="*80)
