#!/usr/bin/env python3
"""
Diagnose VQE Issues - Inspect actual values step by step
"""

import numpy as np
from kanad.bonds import BondFactory

print("="*80)
print("VQE DIAGNOSTIC - STEP BY STEP VALUE INSPECTION")
print("="*80)

# =============================================================================
# STEP 1: Create H2 and verify HF energy
# =============================================================================
print("\nSTEP 1: Create H2 Bond and Verify Hartree-Fock")
print("-"*80)

h2 = BondFactory.create_bond('H', 'H', distance=0.74, basis='sto-3g')
ham = h2.hamiltonian

print(f"Hamiltonian created:")
print(f"  n_orbitals: {ham.n_orbitals}")
print(f"  n_electrons: {ham.n_electrons}")
print(f"  n_qubits (2*n_orbitals): {2 * ham.n_orbitals}")

# Get integrals
print(f"\n📊 Integral Values:")
print(f"  Nuclear repulsion: {ham.nuclear_repulsion:.6f} Ha")

h_one = ham.get_one_body_tensor()
h_two = ham.get_two_body_tensor()
print(f"  One-electron matrix shape: {h_one.shape}")
print(f"  One-electron matrix:\n{h_one}")
print(f"  Two-electron matrix shape: {h_two.shape}")

# HF energy
dm, hf_energy = ham.solve_scf()
print(f"\n📊 Hartree-Fock Results:")
print(f"  HF Energy: {hf_energy:.6f} Ha")
print(f"  Expected (H2 @ 0.74Å): ~-1.117 Ha")
print(f"  Error: {abs(hf_energy - (-1.117)):.6f} Ha")

# =============================================================================
# STEP 2: Build Hamiltonian operator and check eigenvalues
# =============================================================================
print("\nSTEP 2: Build Qubit Hamiltonian and Check Spectrum")
print("-"*80)

# Get full Hamiltonian matrix directly
try:
    n_qubits = 2 * ham.n_orbitals
    ham_matrix = ham.to_matrix(n_qubits=n_qubits)
    print(f"Hamiltonian matrix:")
    print(f"  Shape: {ham_matrix.shape}")
    print(f"  Hermitian: {np.allclose(ham_matrix, ham_matrix.conj().T)}")

    # Get eigenvalues
    eigenvalues = np.linalg.eigvalsh(ham_matrix)
    print(f"\n📊 Exact Spectrum (from diagonalization):")
    print(f"  Ground state:    {eigenvalues[0]:.6f} Ha")
    print(f"  1st excited:     {eigenvalues[1]:.6f} Ha")
    print(f"  2nd excited:     {eigenvalues[2]:.6f} Ha")
    print(f"  3rd excited:     {eigenvalues[3]:.6f} Ha")

    print(f"\n  HF Energy:       {hf_energy:.6f} Ha")
    print(f"  Ground state:    {eigenvalues[0]:.6f} Ha")
    print(f"  Correlation:     {eigenvalues[0] - hf_energy:.6f} Ha")

    if eigenvalues[0] < hf_energy:
        print(f"  ✓ Ground state below HF (correlation energy is negative)")
    else:
        print(f"  ✗ Ground state ABOVE HF (something is wrong!)")

except Exception as e:
    print(f"\n✗ Failed to build matrix: {e}")
    import traceback
    traceback.print_exc()

# =============================================================================
# STEP 3: Test ansatz circuit evaluation
# =============================================================================
print("\nSTEP 3: Test Ansatz Circuit Energy Evaluation")
print("-"*80)

from kanad.ansatze.governance_aware_ansatz import CovalentGovernanceAnsatz

# Create governance ansatz
ansatz = CovalentGovernanceAnsatz(
    n_qubits=2 * ham.n_orbitals,
    n_electrons=ham.n_electrons,
    n_layers=2
)

print(f"Ansatz created:")
print(f"  Type: {type(ansatz).__name__}")
print(f"  n_qubits: {ansatz.n_qubits}")
print(f"  n_parameters: {ansatz.n_parameters}")

# Test with HF initial state (should give HF energy)
print(f"\n📊 Test 1: Zero parameters (should give HF-like energy)")
params_zero = np.zeros(ansatz.n_parameters)
circuit = ansatz.build_circuit(params_zero)
print(f"  Circuit built: {circuit.n_qubits} qubits, depth {circuit.depth}")

# Evaluate energy
try:
    from kanad.backends.statevector_simulator import StatevectorSimulator
    backend = StatevectorSimulator()

    # Get statevector
    statevector = backend.run(circuit)
    print(f"  Statevector shape: {statevector.shape}")
    print(f"  Statevector norm: {np.linalg.norm(statevector):.6f}")

    # Compute expectation value <ψ|H|ψ>
    energy = np.real(np.conj(statevector) @ ham_matrix @ statevector)
    print(f"  Energy (zero params): {energy:.6f} Ha")
    print(f"  HF Energy: {hf_energy:.6f} Ha")
    print(f"  Difference: {abs(energy - hf_energy):.6f} Ha")

except Exception as e:
    print(f"  ✗ Energy evaluation failed: {e}")
    import traceback
    traceback.print_exc()

# Test with random small parameters
print(f"\n📊 Test 2: Small random parameters")
params_random = np.random.random(ansatz.n_parameters) * 0.1
circuit = ansatz.build_circuit(params_random)

try:
    statevector = backend.run(circuit)
    energy = np.real(np.conj(statevector) @ ham_matrix @ statevector)
    print(f"  Energy (random params): {energy:.6f} Ha")
    print(f"  HF Energy: {hf_energy:.6f} Ha")
    print(f"  Ground state: {eigenvalues[0]:.6f} Ha")

    if energy < hf_energy:
        print(f"  ✓ Energy below HF (correlation)")
    else:
        print(f"  ⚠️  Energy above HF (not optimal)")

    if energy > eigenvalues[0]:
        print(f"  ✓ Energy above ground state (variational principle)")
    else:
        print(f"  ✗ Energy BELOW ground state (VIOLATION - something wrong!)")

except Exception as e:
    print(f"  ✗ Energy evaluation failed: {e}")
    import traceback
    traceback.print_exc()

# =============================================================================
# STEP 4: Test VQE energy function
# =============================================================================
print("\nSTEP 4: Test VQE Energy Function")
print("-"*80)

from kanad.solvers.vqe_solver import VQESolver

vqe = VQESolver(
    hamiltonian=ham,
    ansatz_type='governance',
    mapper_type='jordan_wigner',
    max_iterations=10
)

print(f"VQE Solver created:")
print(f"  Ansatz: {type(vqe.ansatz).__name__}")
print(f"  Mapper: {type(vqe.mapper).__name__}")
print(f"  Parameters: {vqe.ansatz.n_parameters}")

# Test energy function directly
print(f"\n📊 Test VQE energy function:")
try:
    # Get the energy function
    if hasattr(vqe, '_energy_evaluation'):
        energy_func = vqe._energy_evaluation
    else:
        print(f"  Available methods: {[m for m in dir(vqe) if not m.startswith('_')]}")

    # Test with zero parameters
    params_test = np.zeros(vqe.ansatz.n_parameters)
    # Need to find the actual energy evaluation method

    print(f"  VQE object attributes:")
    for attr in dir(vqe):
        if 'energy' in attr.lower() or 'evaluate' in attr.lower():
            print(f"    - {attr}")

except Exception as e:
    print(f"  ✗ Failed: {e}")
    import traceback
    traceback.print_exc()

print(f"\n{'='*80}")
print("DIAGNOSTIC COMPLETE")
print("="*80)
