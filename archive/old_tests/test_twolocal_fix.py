#!/usr/bin/env python3
"""Quick test to verify TwoLocal bug fix."""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from api.services.experiment_service import create_molecule_from_config
from kanad.ansatze.two_local_ansatz import TwoLocalAnsatz
from kanad.utils.vqe_solver import VQESolver
import numpy as np

print("=" * 80)
print("TESTING TwoLocal ANSATZ BUG FIX")
print("=" * 80)

# Create H2 molecule
print("\n1. Creating H2 molecule (STO-3G basis)...")
molecule_config = {
    'smiles': '[H][H]',
    'basis': 'sto-3g',
    'charge': 0,
    'multiplicity': 1
}
h2 = create_molecule_from_config(molecule_config)

ham = h2.hamiltonian
n_qubits = 2 * ham.n_orbitals
print(f"   Molecule: {h2.n_electrons} electrons, {ham.n_orbitals} orbitals, {n_qubits} qubits")

# Solve SCF
print("\n2. Running SCF calculation...")
scf_results, hf_energy = ham.solve_scf(max_iterations=50, conv_tol=1e-6)
print(f"   HF Energy: {hf_energy:.8f} Ha")

# Test TwoLocal ansatz
print("\n3. Testing TwoLocalAnsatz...")
try:
    ansatz = TwoLocalAnsatz(
        n_qubits=n_qubits,
        n_electrons=h2.n_electrons,
        n_layers=2,
        rotation_gates='ry',
        entanglement='linear'
    )

    # Test n_parameters property (this was the bug!)
    print(f"   ✅ n_parameters property exists: {ansatz.n_parameters} parameters")

    # Build circuit
    circuit = ansatz.build_circuit()
    print(f"   ✅ Circuit built successfully")

    # Create VQE solver
    print("\n4. Running VQE with TwoLocal ansatz...")
    vqe = VQESolver(
        hamiltonian=ham,
        ansatz=ansatz,
        optimizer='COBYLA',
        max_iterations=50
    )

    result = vqe.solve()

    print(f"\n{'=' * 80}")
    print("RESULTS:")
    print(f"{'=' * 80}")
    print(f"HF Energy:          {hf_energy:.8f} Ha")
    print(f"VQE Energy:         {result['energy']:.8f} Ha")
    print(f"Correlation:        {(result['energy'] - hf_energy) * 1000:.3f} mHa")
    print(f"Function Evals:     {result['function_evaluations']}")
    print(f"Iterations:         {result['iterations']}")
    print(f"Status:             {result['status']}")
    print(f"Circuit Depth:      {result.get('circuit_depth', 'N/A')}")
    print(f"Parameters:         {result.get('n_parameters', ansatz.n_parameters)}")

    correlation_mha = (result['energy'] - hf_energy) * 1000

    if correlation_mha < -5:
        print(f"\n✅ TwoLocal WORKS! Recovered {abs(correlation_mha):.3f} mHa correlation")
        print("✅ Bug fix SUCCESSFUL!")
    elif abs(correlation_mha) < 0.1:
        print(f"\n⚠️  WARNING: No correlation recovered (stuck at HF)")
    else:
        print(f"\n✅ TwoLocal works (some correlation recovered)")

except AttributeError as e:
    print(f"\n❌ AttributeError: {e}")
    print("❌ Bug fix FAILED!")
    sys.exit(1)
except Exception as e:
    print(f"\n❌ Error: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

print("\n" + "=" * 80)
print("TEST COMPLETE")
print("=" * 80)
