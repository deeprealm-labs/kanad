#!/usr/bin/env python3
"""
Local VQE Diagnostic Suite
Tests Hamiltonian and Ansatz implementations with simple molecules
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from api.services.experiment_service import create_molecule_from_config
from kanad.ansatze.ucc_ansatz import UCCAnsatz
from qiskit.quantum_info import Statevector
import numpy as np

def test_molecule(smiles, basis, expected_hf, name):
    """Test a single molecule"""
    print(f"\n{'='*80}")
    print(f"Testing {name} ({smiles})")
    print(f"{'='*80}")

    # Create molecule
    molecule_config = {
        "smiles": smiles,
        "basis": basis,
        "charge": 0,
        "multiplicity": 1
    }

    molecule = create_molecule_from_config(molecule_config)
    ham = molecule.hamiltonian

    print(f"Molecule: {molecule.n_electrons} electrons, {ham.n_orbitals} orbitals")
    print(f"Nuclear repulsion: {ham.nuclear_repulsion:.8f} Ha")

    # Get HF energy from SCF
    _, hf_energy = ham.solve_scf(max_iterations=50, conv_tol=1e-6)
    print(f"\nHF SCF energy: {hf_energy:.8f} Ha")

    if expected_hf:
        error = abs(hf_energy - expected_hf)
        print(f"Expected HF:   {expected_hf:.8f} Ha")
        print(f"Error:         {error:.8f} Ha ({error*627.5:.2f} kcal/mol)")

    # Build Pauli Hamiltonian
    pauli_ham = ham.to_sparse_hamiltonian(mapper="jordan_wigner")
    print(f"\nPauli Hamiltonian: {len(pauli_ham)} terms")

    # Test 1: UCC ansatz with zero parameters (should give HF state)
    print("\n" + "-"*80)
    print("TEST 1: UCC Ansatz HF State (zero parameters)")
    print("-"*80)

    n_qubits = 2 * ham.n_orbitals
    n_electrons = molecule.n_electrons

    ansatz = UCCAnsatz(n_qubits=n_qubits, n_electrons=n_electrons)
    ansatz.build_circuit()

    print(f"Ansatz: {ansatz.n_parameters} parameters")

    # Get HF state (zero parameters)
    zero_params = np.zeros(ansatz.n_parameters)
    qiskit_circuit = ansatz.circuit.to_qiskit()

    if qiskit_circuit.num_parameters > 0:
        param_dict = {qiskit_circuit.parameters[i]: zero_params[i] for i in range(len(zero_params))}
        bound_circuit = qiskit_circuit.assign_parameters(param_dict)
    else:
        bound_circuit = qiskit_circuit

    print(f"Circuit: {bound_circuit.num_qubits} qubits, depth {bound_circuit.depth()}")

    # Get statevector
    statevector = Statevector.from_instruction(bound_circuit)

    # Find dominant components
    nonzero = [(i, abs(statevector.data[i])) for i in range(len(statevector.data)) if abs(statevector.data[i]) > 1e-10]
    print(f"\nNon-zero components: {len(nonzero)}")
    for i, amp in nonzero[:3]:
        binary = format(i, f"0{n_qubits}b")
        print(f"  |{binary}> : {amp:.6f}")

    # Compute energy
    ucc_hf_energy = statevector.expectation_value(pauli_ham).real
    print(f"\nUCC HF state energy: {ucc_hf_energy:.8f} Ha")
    print(f"SCF HF energy:       {hf_energy:.8f} Ha")

    hf_diff = abs(ucc_hf_energy - hf_energy)
    print(f"Difference:          {hf_diff:.8f} Ha ({hf_diff*627.5:.2f} kcal/mol)")

    if hf_diff < 0.001:
        print("✅ PASS: UCC ansatz creates correct HF state")
        test1_pass = True
    else:
        print("❌ FAIL: UCC ansatz HF state is wrong!")
        test1_pass = False

    # Test 2: Run simple VQE optimization
    print("\n" + "-"*80)
    print("TEST 2: Simple VQE Optimization (10 iterations)")
    print("-"*80)

    from kanad.utils.vqe_solver import VQESolver

    vqe = VQESolver(
        hamiltonian=ham,
        ansatz=ansatz,
        optimizer_method='COBYLA',
        max_iterations=10
    )

    try:
        result = vqe.solve()
        vqe_energy = result['energy']
        iterations = result.get('iteration_count', 0)

        print(f"\nVQE energy:      {vqe_energy:.8f} Ha")
        print(f"HF energy:       {hf_energy:.8f} Ha")
        print(f"Iterations:      {iterations}")

        correlation_energy = vqe_energy - hf_energy
        print(f"Correlation:     {correlation_energy:.8f} Ha ({correlation_energy*627.5:.2f} kcal/mol)")

        if vqe_energy < hf_energy - 0.0001:
            print("✅ PASS: VQE energy below HF (correlation energy recovered)")
            test2_pass = True
        else:
            print("⚠️  WARNING: VQE energy not below HF")
            test2_pass = False

    except Exception as e:
        print(f"❌ FAIL: VQE optimization failed: {e}")
        test2_pass = False

    # Summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print(f"Test 1 (HF state):     {'PASS' if test1_pass else 'FAIL'}")
    print(f"Test 2 (VQE optim):    {'PASS' if test2_pass else 'FAIL'}")

    return test1_pass and test2_pass


def main():
    """Run diagnostic tests on benchmark molecules"""
    print("\n" + "="*80)
    print("KANAD VQE DIAGNOSTIC SUITE - LOCAL TESTING")
    print("="*80)

    # Test molecules with known reference energies (STO-3G basis)
    test_cases = [
        # (SMILES, basis, expected_HF, name)
        ("[H][H]", "sto-3g", -1.1167, "H2 (Hydrogen)"),
        ("O", "sto-3g", -74.9659, "H2O (Water)"),
        # ("C", "sto-3g", -39.7266, "CH4 (Methane)"),  # Skip CH4 for now (too big)
    ]

    results = []
    for smiles, basis, expected_hf, name in test_cases:
        try:
            passed = test_molecule(smiles, basis, expected_hf, name)
            results.append((name, passed))
        except Exception as e:
            print(f"\n❌ ERROR testing {name}: {e}")
            import traceback
            traceback.print_exc()
            results.append((name, False))

    # Final summary
    print("\n" + "="*80)
    print("FINAL RESULTS")
    print("="*80)
    for name, passed in results:
        status = "✅ PASS" if passed else "❌ FAIL"
        print(f"{status}: {name}")

    all_pass = all(passed for _, passed in results)
    if all_pass:
        print("\n🎉 All tests passed!")
        return 0
    else:
        print("\n⚠️  Some tests failed - needs investigation")
        return 1


if __name__ == "__main__":
    sys.exit(main())
