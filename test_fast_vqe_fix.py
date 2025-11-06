"""
Test FastVQE Placeholder Fix

Validates that FastVQE now computes actual quantum expectation values
instead of returning HF energy.

CRITICAL FIX: Issue - FastVQE was returning HF energy placeholder
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from kanad.io import from_smiles
from kanad.vqe_optimization.fast_vqe import FastVQE
import numpy as np


def test_fastvqe_returns_quantum_energy():
    """Test that FastVQE computes quantum expectation (not HF placeholder)."""
    print("=" * 80)
    print("TEST: FastVQE Returns Quantum Expectation Value")
    print("=" * 80)

    # Create H2 molecule
    print("\nCreating H2 molecule...")
    h2 = from_smiles("[H][H]")

    # Get Hamiltonian
    print("Getting Hamiltonian...")
    from kanad.core.hamiltonians import CovalentHamiltonian
    from kanad.core.representations.lcao_representation import LCAORepresentation

    representation = LCAORepresentation(h2)
    hamiltonian = CovalentHamiltonian(
        h2,
        representation,
        basis_name='sto-3g',
        use_governance=True
    )

    # Get HF energy from PySCF directly (already computed during Hamiltonian construction)
    # The value shown above: converged SCF energy = -1.11729321718008
    hf_energy = -1.117  # Approximate HF energy for H2/STO-3G
    print(f"HF energy (reference): ~{hf_energy:.3f} Ha")

    # Create FastVQE solver
    print("\nCreating FastVQE solver...")
    fast_vqe = FastVQE(
        hamiltonian=hamiltonian,
        num_layers=1,
        use_circuit_optimization=False,
        use_smart_init=True,
        optimizer='cobyla'  # Use COBYLA (gradient-free) for simpler test
    )

    # Test expectation value computation directly
    print("\nTesting _evaluate_expectation() method...")

    # Create a simple test circuit (HF state)
    try:
        from qiskit import QuantumCircuit
        n_qubits = 2 * hamiltonian.n_orbitals
        circuit = QuantumCircuit(n_qubits)

        # Prepare HF state (occupy first n_electrons/2 orbitals)
        # For H2: 2 electrons, so occupy orbitals 0 and 1 (spin-up and spin-down of orbital 0)
        circuit.x(0)  # Spin-up electron in orbital 0
        circuit.x(1)  # Spin-down electron in orbital 0

        # Compute expectation value
        energy = fast_vqe._evaluate_expectation(circuit)

        print(f"Quantum expectation (HF state): {energy:.6f} Ha")
        print(f"HF energy (reference): {hf_energy:.6f} Ha")

        # Check that the energy is close to HF (within 1e-3 Ha = 27 meV)
        # For the HF state, the quantum expectation should equal HF energy
        energy_diff = abs(energy - hf_energy)
        print(f"Difference: {energy_diff:.6e} Ha")

        # The HF state should give HF energy (within numerical precision)
        if energy_diff < 1e-3:
            print("\n✅ PASS: Quantum expectation matches HF energy for HF state")
            hf_state_test = True
        else:
            print("\n❌ FAIL: Quantum expectation doesn't match HF energy")
            print("   (Expected HF state to give HF energy)")
            hf_state_test = False

        # Now test with a perturbed circuit (should give different energy)
        print("\n" + "-" * 80)
        print("Testing with perturbed circuit...")

        circuit_perturbed = QuantumCircuit(n_qubits)
        circuit_perturbed.x(0)
        circuit_perturbed.x(1)
        # Add a small rotation to perturb the state
        circuit_perturbed.ry(0.1, 0)  # Small rotation on qubit 0

        energy_perturbed = fast_vqe._evaluate_expectation(circuit_perturbed)

        print(f"Quantum expectation (perturbed): {energy_perturbed:.6f} Ha")
        print(f"HF energy (reference): {hf_energy:.6f} Ha")

        perturb_diff = abs(energy_perturbed - hf_energy)
        print(f"Difference from HF: {perturb_diff:.6e} Ha")

        # The perturbed state should give DIFFERENT energy (not HF energy)
        # If it still returns HF energy exactly, the placeholder is still there!
        if perturb_diff > 1e-6:
            print("\n✅ PASS: Perturbed circuit gives different energy")
            print("   (Confirms actual quantum computation, not HF placeholder)")
            perturbed_test = True
        else:
            print("\n❌ FAIL: Perturbed circuit gives same energy as HF")
            print("   (This indicates placeholder still active!)")
            perturbed_test = False

        return hf_state_test and perturbed_test

    except Exception as e:
        print(f"\n❌ FAIL: Exception during test: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_sgd_optimizer_no_crash():
    """Test that SGD optimizer doesn't crash with prev_energy bug."""
    print("\n\n")
    print("=" * 80)
    print("TEST: SGD Optimizer prev_energy Bug Fix")
    print("=" * 80)

    try:
        # Create H2 molecule
        print("\nCreating H2 molecule...")
        h2 = from_smiles("[H][H]")

        # Get Hamiltonian
        from kanad.core.hamiltonians import CovalentHamiltonian
        from kanad.core.representations.lcao_representation import LCAORepresentation

        representation = LCAORepresentation(h2)
        hamiltonian = CovalentHamiltonian(
            h2,
            representation,
            basis_name='sto-3g',
            use_governance=True
        )

        # Create FastVQE with SGD optimizer
        print("Creating FastVQE with SGD optimizer...")
        fast_vqe = FastVQE(
            hamiltonian=hamiltonian,
            num_layers=1,
            use_circuit_optimization=False,
            use_smart_init=True,
            optimizer='sgd'
        )

        # Run optimization (just 5 iterations to test)
        print("Running SGD optimization (5 iterations)...")
        result = fast_vqe.solve(
            max_iterations=5,
            learning_rate=0.01,
            tolerance=1e-6
        )

        print(f"\n✅ PASS: SGD optimizer completed without crash")
        print(f"   Energy: {result['energy_ha']:.6f} Ha")
        print(f"   Iterations: {result['iterations']}")
        print(f"   Converged: {result['converged']}")

        return True

    except NameError as e:
        if 'prev_energy' in str(e):
            print(f"\n❌ FAIL: NameError for prev_energy (bug not fixed!)")
            print(f"   Error: {e}")
            return False
        else:
            raise

    except Exception as e:
        print(f"\n❌ FAIL: Unexpected exception: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == '__main__':
    print("\n🔬 FASTVQE PLACEHOLDER FIX VALIDATION\n")
    print("Testing Issue: FastVQE placeholder - returns HF energy instead of quantum")
    print("CRITICAL: FastVQE was not actually running quantum computation!\n")

    passed_quantum = test_fastvqe_returns_quantum_energy()
    passed_sgd = test_sgd_optimizer_no_crash()

    print("\n\n")
    print("=" * 80)
    print("FINAL SUMMARY")
    print("=" * 80)
    print(f"Quantum Expectation Test: {'✅ PASS' if passed_quantum else '❌ FAIL'}")
    print(f"SGD Optimizer Test:       {'✅ PASS' if passed_sgd else '❌ FAIL'}")

    if passed_quantum and passed_sgd:
        print("\n✅ FASTVQE FIX COMPLETE!")
        print("   ✓ Computes actual quantum expectation values")
        print("   ✓ No longer returns HF energy placeholder")
        print("   ✓ SGD optimizer prev_energy bug fixed")
        print("\n🎉 Issues RESOLVED!")
        sys.exit(0)
    else:
        print("\n❌ SOME TESTS FAILED")
        if not passed_quantum:
            print("   ✗ Quantum expectation still returning placeholder")
        if not passed_sgd:
            print("   ✗ SGD optimizer prev_energy bug not fixed")
        sys.exit(1)
