#!/usr/bin/env python3
"""
Test ALL ansatze in the framework
"""

import numpy as np
from kanad.bonds import BondFactory
from kanad.solvers.vqe_solver import VQESolver

print("="*80)
print("TEST ALL ANSATZE")
print("="*80)

# Create H2
h2 = BondFactory.create_bond('H', 'H', distance=0.74, basis='sto-3g')
ham = h2.hamiltonian

# Get reference
dm, hf_energy = ham.solve_scf()
ham_matrix = ham.to_matrix(n_qubits=4)
exact_energy = np.linalg.eigvalsh(ham_matrix)[0]

print(f"\n📊 Reference:")
print(f"  HF:    {hf_energy:.8f} Ha")
print(f"  Exact: {exact_energy:.8f} Ha")

# Test all ansatz types
ansatz_types = ['ucc', 'governance', 'hardware_efficient']

results = []

for ansatz_type in ansatz_types:
    print(f"\n{'='*80}")
    print(f"Testing Ansatz: {ansatz_type}")
    print("="*80)

    try:
        # Test with best optimizer (SLSQP)
        vqe = VQESolver(
            hamiltonian=ham,
            ansatz_type=ansatz_type,
            mapper_type='jordan_wigner',
            optimizer='SLSQP',
            max_iterations=500
        )

        print(f"  Created successfully")
        print(f"  Parameters: {vqe.ansatz.n_parameters}")

        # Test circuit builds
        print(f"  Testing circuit build...")
        params_test = np.zeros(vqe.ansatz.n_parameters)
        circuit = vqe.ansatz.build_circuit(params_test) if ansatz_type != 'ucc' else vqe.ansatz.build_circuit()

        print(f"  Circuit: {circuit.n_qubits} qubits, depth {circuit.depth}")

        # Test energy evaluation
        print(f"  Testing energy evaluation...")
        energy_test = vqe._compute_energy_statevector(params_test)
        print(f"  Energy (zero params): {energy_test:.8f} Ha")

        # Run VQE
        print(f"  Running VQE optimization...")
        result = vqe.solve()

        energy = result['energy']
        error = abs(energy - exact_energy)
        correlation = energy - hf_energy

        print(f"\n  Results:")
        print(f"    Final Energy:  {energy:.8f} Ha")
        print(f"    Error:         {error:.8f} Ha ({error*1000:.3f} mHa)")
        print(f"    Correlation:   {correlation:.8f} Ha")
        print(f"    Converged:     {result.get('converged', False)}")
        print(f"    Iterations:    {result.get('iterations', 'N/A')}")

        if error < 0.001:
            print(f"    ✓✓ EXCELLENT")
        elif error < 0.01:
            print(f"    ✓ GOOD")
        else:
            print(f"    ✗ POOR")

        results.append({
            'ansatz': ansatz_type,
            'success': True,
            'energy': energy,
            'error': error,
            'correlation': correlation,
            'n_params': vqe.ansatz.n_parameters
        })

    except Exception as e:
        print(f"  ✗ FAILED: {e}")
        import traceback
        traceback.print_exc()

        results.append({
            'ansatz': ansatz_type,
            'success': False,
            'error_msg': str(e)
        })

# Summary
print(f"\n{'='*80}")
print("SUMMARY")
print("="*80)

successful = [r for r in results if r['success']]
failed = [r for r in results if not r['success']]

print(f"\n✓ Successful: {len(successful)}/{len(results)}")
for r in successful:
    print(f"  - {r['ansatz']}: {r['error']*1000:.3f} mHa error, {r['n_params']} params")

if failed:
    print(f"\n✗ Failed: {len(failed)}/{len(results)}")
    for r in failed:
        print(f"  - {r['ansatz']}: {r.get('error_msg', 'Unknown error')}")

print(f"\n{'='*80}")
