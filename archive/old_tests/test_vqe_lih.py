#!/usr/bin/env python3
"""Test VQE with LiH - has significant correlation energy"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from api.services.experiment_service import create_molecule_from_config
from kanad.ansatze.ucc_ansatz import UCCAnsatz
from kanad.utils.vqe_solver import VQESolver

# Test LiH - known to have significant correlation energy
print("="*80)
print("Testing LiH VQE")
print("Expected HF: ~-7.86 Ha, Exact FCI: ~-7.88 Ha (correlation ~0.02 Ha)")
print("="*80)

molecule_config = {
    "smiles": "[Li]",  # Will form LiH
    "basis": "sto-3g",
    "charge": 0,
    "multiplicity": 1
}

molecule = create_molecule_from_config(molecule_config)
ham = molecule.hamiltonian

print(f"\nMolecule: {molecule.n_electrons} electrons, {ham.n_orbitals} orbitals")

# Get HF energy
_, hf_energy = ham.solve_scf(max_iterations=50, conv_tol=1e-6)
print(f"HF energy: {hf_energy:.8f} Ha")

# Create ansatz
n_qubits = 2 * ham.n_orbitals
ansatz = UCCAnsatz(n_qubits=n_qubits, n_electrons=molecule.n_electrons)

print(f"Ansatz parameters: {ansatz.n_parameters}")

# Run VQE with COBYLA
print(f"\nRunning VQE with COBYLA optimizer (50 iterations)...")

vqe = VQESolver(
    hamiltonian=ham,
    ansatz=ansatz,
    optimizer_method='COBYLA',
    max_iterations=50
)

result = vqe.solve()

print(f"\n" + "="*80)
print("RESULTS")
print("="*80)
print(f"HF energy:          {hf_energy:.8f} Ha")
print(f"VQE energy:         {result['energy']:.8f} Ha")
print(f"Correlation energy: {result['energy'] - hf_energy:.8f} Ha ({(result['energy'] - hf_energy)*627.5:.2f} kcal/mol)")
print(f"Function evals:     {result.get('function_evaluations', 'N/A')}")
print(f"Iterations:         {result.get('iterations', 'N/A')}")
print(f"Converged:          {result.get('converged', False)}")

# Check if we recovered correlation energy
if result['energy'] < hf_energy - 0.001:
    print(f"\n✅ SUCCESS: VQE recovered correlation energy!")
    print(f"   Energy is {hf_energy - result['energy']:.6f} Ha below HF")
else:
    print(f"\n❌ FAIL: VQE did not recover correlation energy")
    print(f"   Energy difference from HF: {result['energy'] - hf_energy:.6f} Ha")

# Show final energy comparison
print(f"\n" + "-"*80)
print(f"Energy lowering: {(hf_energy - result['energy'])*627.5:.2f} kcal/mol")
