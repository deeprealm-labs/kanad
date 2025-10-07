#!/usr/bin/env python3
"""
Test VQE energy function with multiple calls (different parameters)
"""

import numpy as np
from kanad.bonds import BondFactory
from kanad.solvers.vqe_solver import VQESolver

print("="*80)
print("TEST: VQE Energy Function - Multiple Calls")
print("="*80)

# Create H2
h2 = BondFactory.create_bond('H', 'H', distance=0.74, basis='sto-3g')
ham = h2.hamiltonian

# Get reference energy
dm, hf_energy = ham.solve_scf()
ham_matrix = ham.to_matrix(n_qubits=4)
eigenvalues = np.linalg.eigvalsh(ham_matrix)

print(f"\n📊 Reference:")
print(f"  HF Energy:    {hf_energy:.6f} Ha")
print(f"  Ground state: {eigenvalues[0]:.6f} Ha")

# Create VQE
vqe = VQESolver(
    hamiltonian=ham,
    ansatz_type='governance',
    mapper_type='jordan_wigner',
    max_iterations=5
)

print(f"\n📊 VQE Setup:")
print(f"  Ansatz: {type(vqe.ansatz).__name__}")
print(f"  n_parameters: {vqe.ansatz.n_parameters}")

# Call energy function multiple times with different parameters
print(f"\n📊 Energy Evaluations:")

params_test = [
    np.zeros(vqe.ansatz.n_parameters),
    np.ones(vqe.ansatz.n_parameters) * 0.1,
    np.random.random(vqe.ansatz.n_parameters) * 0.1,
    np.random.random(vqe.ansatz.n_parameters) * 0.5,
]

for i, params in enumerate(params_test):
    energy = vqe._compute_energy_statevector(params)
    print(f"\n  Call {i+1}:")
    print(f"    params[:3]: {params[:3]}")
    print(f"    Energy: {energy:.6f} Ha")

    # Check if energy is physical
    if energy < eigenvalues[0]:
        print(f"    ✗ BELOW GROUND STATE (violation!)")
    elif energy > hf_energy + 10:
        print(f"    ✗ WAY TOO HIGH (unphysical)")
    elif energy > hf_energy:
        print(f"    ⚠️  Above HF (not finding correlation)")
    else:
        print(f"    ✓ Between ground state and HF")

print(f"\n{'='*80}")
