#!/usr/bin/env python3
"""
Debug why UCC ansatz fails to find correlation energy
"""

import numpy as np
from kanad.bonds import BondFactory
from kanad.ansatze.ucc_ansatz import UCCAnsatz
from qiskit.quantum_info import Statevector

print("="*80)
print("DEBUG UCC ANSATZ FAILURE")
print("="*80)

# Create H2
h2 = BondFactory.create_bond('H', 'H', distance=0.74, basis='sto-3g')
ham = h2.hamiltonian

# Get reference
dm, hf_energy = ham.solve_scf()
ham_matrix = ham.to_matrix(n_qubits=4)
exact_energy = np.linalg.eigvalsh(ham_matrix)[0]

print(f"\nReference:")
print(f"  HF:    {hf_energy:.8f} Ha")
print(f"  Exact: {exact_energy:.8f} Ha")
print(f"  Correlation needed: {exact_energy - hf_energy:.8f} Ha")

# Create UCC ansatz
ansatz = UCCAnsatz(n_qubits=4, n_electrons=2, include_singles=True, include_doubles=True)

print(f"\nUCC Ansatz:")
print(f"  n_qubits: {ansatz.n_qubits}")
print(f"  n_electrons: {ansatz.n_electrons}")
print(f"  Excitations: {ansatz.excitations}")
print(f"  n_parameters: {ansatz.n_parameters}")

# Build circuit with zero parameters
print(f"\n📊 Test 1: Zero parameters (should give HF)")
params_zero = np.zeros(ansatz.n_parameters)
circuit = ansatz.build_circuit()

qiskit_circuit = circuit.to_qiskit()

if qiskit_circuit.num_parameters > 0:
    param_dict = {qiskit_circuit.parameters[i]: params_zero[i] for i in range(len(params_zero))}
    bound_circuit = qiskit_circuit.assign_parameters(param_dict)
else:
    bound_circuit = qiskit_circuit

statevector = Statevector.from_instruction(bound_circuit)
psi = statevector.data

energy_zero = np.real(np.conj(psi) @ ham_matrix @ psi)
print(f"  Energy: {energy_zero:.8f} Ha")
print(f"  Error from HF: {abs(energy_zero - hf_energy):.8f} Ha")

# Test with parameters that should give correlation
print(f"\n📊 Test 2: Grid search for best single parameter")

# Test each parameter individually
best_energy = hf_energy
best_params = params_zero.copy()
best_param_idx = -1

for param_idx in range(ansatz.n_parameters):
    # Test range of values for this parameter
    test_values = np.linspace(-np.pi, np.pi, 20)

    for val in test_values:
        params_test = params_zero.copy()
        params_test[param_idx] = val

        # Build and evaluate
        circuit_test = ansatz.build_circuit()
        qiskit_circuit_test = circuit_test.to_qiskit()

        if qiskit_circuit_test.num_parameters > 0:
            param_dict_test = {qiskit_circuit_test.parameters[i]: params_test[i] for i in range(len(params_test))}
            bound_circuit_test = qiskit_circuit_test.assign_parameters(param_dict_test)
        else:
            bound_circuit_test = qiskit_circuit_test

        statevector_test = Statevector.from_instruction(bound_circuit_test)
        psi_test = statevector_test.data

        energy_test = np.real(np.conj(psi_test) @ ham_matrix @ psi_test)

        if energy_test < best_energy:
            best_energy = energy_test
            best_params = params_test.copy()
            best_param_idx = param_idx

print(f"\n  Best single-parameter result:")
print(f"    Parameter index: {best_param_idx}")
if best_param_idx >= 0:
    print(f"    Parameter value: {best_params[best_param_idx]:.6f}")
    print(f"    Excitation: {ansatz.excitations[best_param_idx] if best_param_idx < len(ansatz.excitations) else 'N/A'}")
print(f"    Energy: {best_energy:.8f} Ha")
print(f"    Correlation: {best_energy - hf_energy:.8f} Ha")
print(f"    Error from exact: {abs(best_energy - exact_energy):.8f} Ha")

if best_energy < hf_energy:
    print(f"  ✓ Found correlation with single parameter!")
else:
    print(f"  ✗ Could NOT find correlation even with grid search!")
    print(f"\n  This suggests UCC ansatz is NOT expressive enough for H2")
    print(f"  OR there's a bug in the UCC implementation")

# Manual check: what is the HF state in this basis?
print(f"\n📊 Test 3: Analyze HF state")
print(f"  HF state amplitudes (from zero-param circuit):")
for i, amp in enumerate(psi):
    if np.abs(amp) > 0.01:
        print(f"    |{i:04b}⟩: {amp:.6f}")

# Check what the exact ground state looks like
eigenvalues, eigenvectors = np.linalg.eigh(ham_matrix)
exact_gs = eigenvectors[:, 0]
print(f"\n  Exact ground state amplitudes:")
for i, amp in enumerate(exact_gs):
    if np.abs(amp) > 0.01:
        print(f"    |{i:04b}⟩: {amp:.6f}")

print(f"\n{'='*80}")
