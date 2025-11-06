#!/usr/bin/env python3
"""
Diagnose API VQE Issue
Test VQE using the same path as the API/frontend
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

print("="*80)
print("API VQE DIAGNOSTIC")
print("Testing VQE using high-level API (same as frontend)")
print("="*80)

# Import API service
from api.services.experiment_service import create_molecule_from_config, run_vqe_experiment

# Create H2 molecule
molecule_config = {
    "smiles": "[H][H]",
    "basis": "sto-3g",
    "charge": 0,
    "multiplicity": 1
}

vqe_config = {
    "ansatz": "covalent_governance",
    "optimizer": "COBYLA",
    "max_iterations": 50,  # Short test
    "mapper": "jordan_wigner"
}

print("\n📋 Configuration:")
print(f"  Molecule: {molecule_config['smiles']}")
print(f"  Basis: {molecule_config['basis']}")
print(f"  Ansatz: {vqe_config['ansatz']}")
print(f"  Optimizer: {vqe_config['optimizer']}")
print(f"  Max iterations: {vqe_config['max_iterations']}")

# Create molecule
print("\n🔬 Creating molecule...")
molecule = create_molecule_from_config(molecule_config)
print(f"✓ Molecule created: {molecule.n_atoms} atoms, {molecule.n_electrons} electrons")

# Track energy history
energy_history = []

def progress_callback(iteration, energy, parameters):
    energy_history.append({'iter': iteration, 'energy': energy})
    if iteration % 5 == 0:
        print(f"  Iter {iteration}: E = {energy:.8f} Ha")

# Run VQE
print("\n🚀 Running VQE...")
try:
    result = run_vqe_experiment(
        molecule=molecule,
        ansatz_type=vqe_config['ansatz'],
        optimizer=vqe_config['optimizer'],
        max_iterations=vqe_config['max_iterations'],
        mapper_type=vqe_config['mapper'],
        backend='classical',
        callback=progress_callback
    )

    print("\n✅ VQE Completed")
    print(f"  Energy: {result['energy']:.8f} Ha")
    print(f"  HF Energy: {result.get('hf_energy', 0):.8f} Ha")
    print(f"  Correlation: {result.get('correlation_energy', 0):.8f} Ha")
    print(f"  Converged: {result.get('converged', False)}")
    print(f"  Iterations: {result.get('iterations', 0)}")

    # Analyze energy history
    if len(energy_history) > 0:
        energies = [e['energy'] for e in energy_history]
        print(f"\n📊 Energy History Analysis:")
        print(f"  Total evaluations: {len(energies)}")
        print(f"  Initial energy: {energies[0]:.8f} Ha")
        print(f"  Final energy: {energies[-1]:.8f} Ha")
        print(f"  Min energy: {min(energies):.8f} Ha")
        print(f"  Max energy: {max(energies):.8f} Ha")
        print(f"  Energy range: {max(energies) - min(energies):.8f} Ha")

        # Check if energy is changing
        energy_std = np.std(energies)
        print(f"  Energy std dev: {energy_std:.8f} Ha")

        if energy_std < 1e-6:
            print(f"\n  ❌ PROBLEM: Energy barely changing!")
            print(f"     All energies essentially identical")
            print(f"     Optimizer not exploring parameter space")
        elif energy_std < 0.001:
            print(f"\n  ⚠️  WARNING: Small energy variation")
            print(f"     Optimizer may be stuck")
        else:
            print(f"\n  ✓ Energy varying significantly")

        # Check for improvement
        improvement = energies[0] - min(energies)
        print(f"  Energy improvement: {improvement:.8f} Ha ({improvement*627.5:.4f} kcal/mol)")

        if improvement < 1e-6:
            print(f"     ❌ No improvement! VQE stuck at initial state")
        elif improvement < 0.001:
            print(f"     ⚠️  Minimal improvement")
        else:
            print(f"     ✓ Significant improvement")

except Exception as e:
    print(f"\n❌ VQE Failed: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "="*80)

import numpy as np
