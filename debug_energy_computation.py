#!/usr/bin/env python3
"""
Debug VQE Energy Computation - Find WHERE the energy calculation goes wrong
"""

import numpy as np
from kanad.bonds import BondFactory
from qiskit.quantum_info import Statevector

print("="*80)
print("DEBUG: VQE ENERGY COMPUTATION")
print("="*80)

# Create H2
h2 = BondFactory.create_bond('H', 'H', distance=0.74, basis='sto-3g')
ham = h2.hamiltonian

# Get exact values
dm, hf_energy = ham.solve_scf()
ham_matrix = ham.to_matrix(n_qubits=4)
eigenvalues = np.linalg.eigvalsh(ham_matrix)

print(f"\n📊 Reference Values:")
print(f"  HF Energy:      {hf_energy:.6f} Ha")
print(f"  Ground state:   {eigenvalues[0]:.6f} Ha")
print(f"  Correlation:    {eigenvalues[0] - hf_energy:.6f} Ha")

# Create governance ansatz
from kanad.ansatze.governance_aware_ansatz import CovalentGovernanceAnsatz

ansatz = CovalentGovernanceAnsatz(
    n_qubits=4,
    n_electrons=2,
    n_layers=2
)

print(f"\n📊 Ansatz:")
print(f"  n_qubits: {ansatz.n_qubits}")
print(f"  n_parameters: {ansatz.n_parameters}")

# Test 1: Zero parameters (should give HF-like state)
print(f"\n" + "="*80)
print("TEST 1: Zero Parameters")
print("="*80)

params_zero = np.zeros(ansatz.n_parameters)
circuit = ansatz.build_circuit(params_zero)

print(f"\nCircuit:")
print(f"  n_qubits: {circuit.n_qubits}")
print(f"  depth: {circuit.depth}")
print(f"  gates: {len(circuit.gates)}")

# Convert to Qiskit and get statevector
qiskit_circuit = circuit.to_qiskit()
print(f"\nQiskit circuit:")
print(f"  n_qubits: {qiskit_circuit.num_qubits}")
print(f"  num_parameters: {qiskit_circuit.num_parameters}")

# Bind parameters if needed
if qiskit_circuit.num_parameters > 0:
    param_dict = {qiskit_circuit.parameters[i]: params_zero[i] for i in range(len(params_zero))}
    bound_circuit = qiskit_circuit.assign_parameters(param_dict)
else:
    bound_circuit = qiskit_circuit

# Get statevector
statevector = Statevector.from_instruction(bound_circuit)
psi = statevector.data

print(f"\nStatevector:")
print(f"  Dimension: {len(psi)}")
print(f"  Norm: {np.linalg.norm(psi):.6f}")
print(f"  Non-zero elements: {np.count_nonzero(np.abs(psi) > 1e-10)}")

# Show which basis states have amplitude
print(f"\n  Basis states with amplitude > 0.01:")
for i, amp in enumerate(psi):
    if np.abs(amp) > 0.01:
        print(f"    |{i:04b}⟩: {amp:.4f} (prob: {np.abs(amp)**2:.4f})")

# Compute energy manually
energy_manual = np.real(np.conj(psi) @ ham_matrix @ psi)
print(f"\n📊 Energy (zero params):")
print(f"  Manual calculation: {energy_manual:.6f} Ha")
print(f"  HF Energy:          {hf_energy:.6f} Ha")
print(f"  Difference:         {abs(energy_manual - hf_energy):.6f} Ha")

if abs(energy_manual - hf_energy) < 0.01:
    print(f"  ✓ Close to HF (as expected for zero params)")
else:
    print(f"  ⚠️  NOT close to HF!")

# Test 2: Use VQE solver's energy function
print(f"\n" + "="*80)
print("TEST 2: VQE Solver Energy Function")
print("="*80)

from kanad.solvers.vqe_solver import VQESolver

vqe = VQESolver(
    hamiltonian=ham,
    ansatz_type='governance',
    mapper_type='jordan_wigner',
    max_iterations=1
)

# Call the energy function directly
energy_vqe = vqe._compute_energy_statevector(params_zero)

print(f"\n📊 VQE Energy Function:")
print(f"  VQE result:         {energy_vqe:.6f} Ha")
print(f"  Manual result:      {energy_manual:.6f} Ha")
print(f"  Difference:         {abs(energy_vqe - energy_manual):.6f} Ha")

if abs(energy_vqe - energy_manual) < 1e-6:
    print(f"  ✓ VQE matches manual calculation")
else:
    print(f"  ✗ VQE DOES NOT MATCH manual calculation!")
    print(f"\n  Debugging VQE internals:")
    print(f"    _hamiltonian_matrix cached: {vqe._hamiltonian_matrix is not None}")
    if vqe._hamiltonian_matrix is not None:
        print(f"    _hamiltonian_matrix shape: {vqe._hamiltonian_matrix.shape}")
        print(f"    _needs_padding: {getattr(vqe, '_needs_padding', 'N/A')}")
        print(f"    _ansatz_qubits: {getattr(vqe, '_ansatz_qubits', 'N/A')}")
        print(f"    _full_qubits: {getattr(vqe, '_full_qubits', 'N/A')}")

# Test 3: Check Hamiltonian matrix used by VQE
print(f"\n" + "="*80)
print("TEST 3: Compare Hamiltonian Matrices")
print("="*80)

ham_matrix_vqe = vqe._hamiltonian_matrix
ham_matrix_direct = ham.to_matrix(n_qubits=4, use_mo_basis=True)

print(f"\nHamiltonian Matrix Comparison:")
print(f"  Direct shape: {ham_matrix_direct.shape}")
print(f"  VQE shape:    {ham_matrix_vqe.shape}")
print(f"  Matrices equal: {np.allclose(ham_matrix_direct, ham_matrix_vqe)}")

if not np.allclose(ham_matrix_direct, ham_matrix_vqe):
    print(f"  ✗ MATRICES ARE DIFFERENT!")
    print(f"  Max difference: {np.max(np.abs(ham_matrix_direct - ham_matrix_vqe)):.6e}")

    # Check diagonal elements
    diag_direct = np.diag(ham_matrix_direct)
    diag_vqe = np.diag(ham_matrix_vqe)

    print(f"\n  Diagonal comparison (first 8 elements):")
    for i in range(min(8, len(diag_direct))):
        print(f"    [{i:02d}] Direct: {diag_direct[i]:.6f}, VQE: {diag_vqe[i]:.6f}, Diff: {abs(diag_direct[i] - diag_vqe[i]):.6e}")

# Test 4: Check if ansatz actually changes with parameters
print(f"\n" + "="*80)
print("TEST 4: Parameter Sensitivity")
print("="*80)

# Small random parameters
params_small = np.random.random(ansatz.n_parameters) * 0.01

circuit_small = ansatz.build_circuit(params_small)
qiskit_circuit_small = circuit_small.to_qiskit()

if qiskit_circuit_small.num_parameters > 0:
    param_dict = {qiskit_circuit_small.parameters[i]: params_small[i] for i in range(len(params_small))}
    bound_circuit_small = qiskit_circuit_small.assign_parameters(param_dict)
else:
    bound_circuit_small = qiskit_circuit_small

statevector_small = Statevector.from_instruction(bound_circuit_small)
psi_small = statevector_small.data

energy_small = np.real(np.conj(psi_small) @ ham_matrix @ psi_small)

print(f"\n📊 Energy with small random parameters:")
print(f"  Parameters norm: {np.linalg.norm(params_small):.6f}")
print(f"  Energy: {energy_small:.6f} Ha")
print(f"  Zero params energy: {energy_manual:.6f} Ha")
print(f"  Change: {energy_small - energy_manual:.6f} Ha")

if abs(energy_small - energy_manual) > 1e-6:
    print(f"  ✓ Ansatz is sensitive to parameters")
else:
    print(f"  ✗ Ansatz NOT sensitive to parameters!")

print(f"\n{'='*80}")
print("DEBUG COMPLETE")
print("="*80)
