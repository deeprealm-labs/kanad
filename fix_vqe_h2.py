#!/usr/bin/env python3
"""
Fix VQE to get correct H2 ground state energy: -1.137284 Ha

Test different optimizers systematically
"""

import numpy as np
from kanad.bonds import BondFactory
from kanad.solvers.vqe_solver import VQESolver

print("="*80)
print("FIX VQE - SYSTEMATIC OPTIMIZER TESTING")
print("="*80)

# Create H2
h2 = BondFactory.create_bond('H', 'H', distance=0.74, basis='sto-3g')
ham = h2.hamiltonian

# Get reference energies
dm, hf_energy = ham.solve_scf()
ham_matrix = ham.to_matrix(n_qubits=4)
eigenvalues = np.linalg.eigvalsh(ham_matrix)

print(f"\n📊 Reference Energies:")
print(f"  HF Energy:       {hf_energy:.8f} Ha")
print(f"  Exact Ground:    {eigenvalues[0]:.8f} Ha")
print(f"  Correlation:     {eigenvalues[0] - hf_energy:.8f} Ha ({abs(eigenvalues[0] - hf_energy)/abs(hf_energy)*100:.3f}%)")
print(f"  Target Error:    < 1 mHa (0.001 Ha)")

# Test different optimizers
optimizers = [
    ('COBYLA', {}),
    ('SLSQP', {}),
    ('L-BFGS-B', {}),
    ('Nelder-Mead', {}),
    ('Powell', {}),
]

results = []

for opt_name, opt_kwargs in optimizers:
    print(f"\n{'='*80}")
    print(f"Testing Optimizer: {opt_name}")
    print("="*80)

    for ansatz_type in ['ucc', 'governance']:
        print(f"\n  Ansatz: {ansatz_type}")

        try:
            vqe = VQESolver(
                hamiltonian=ham,
                ansatz_type=ansatz_type,
                mapper_type='jordan_wigner',
                optimizer=opt_name,
                max_iterations=500,  # Increase iterations
                conv_threshold=1e-6
            )

            print(f"    Parameters: {vqe.ansatz.n_parameters}")

            # Run optimization
            result = vqe.solve()

            energy = result['energy']
            error = abs(energy - eigenvalues[0])
            iterations = result.get('iterations', len(vqe.energy_history))
            converged = result.get('converged', False)

            print(f"    Final Energy:  {energy:.8f} Ha")
            print(f"    Error:         {error:.8f} Ha ({error*1000:.3f} mHa)")
            print(f"    Iterations:    {iterations}")
            print(f"    Converged:     {converged}")

            # Check if we found correlation
            correlation = energy - hf_energy
            if correlation < -0.001:  # At least 1 mHa correlation
                print(f"    ✓ Found correlation: {correlation:.8f} Ha")
            else:
                print(f"    ✗ No correlation found: {correlation:.8f} Ha")

            # Check accuracy
            if error < 0.001:  # 1 mHa accuracy
                print(f"    ✓✓ EXCELLENT ACCURACY!")
            elif error < 0.01:  # 10 mHa
                print(f"    ✓ Good accuracy")
            else:
                print(f"    ✗ Poor accuracy")

            results.append({
                'optimizer': opt_name,
                'ansatz': ansatz_type,
                'energy': energy,
                'error': error,
                'correlation': correlation,
                'iterations': iterations,
                'converged': converged
            })

        except Exception as e:
            print(f"    ✗ FAILED: {e}")
            results.append({
                'optimizer': opt_name,
                'ansatz': ansatz_type,
                'energy': None,
                'error': None,
                'correlation': None,
                'iterations': 0,
                'converged': False,
                'error_msg': str(e)
            })

# Summary
print(f"\n{'='*80}")
print("SUMMARY - Best Results")
print("="*80)

# Sort by error
valid_results = [r for r in results if r['error'] is not None]
if valid_results:
    valid_results.sort(key=lambda x: x['error'])

    print(f"\nTop 5 Best Configurations:")
    for i, r in enumerate(valid_results[:5]):
        print(f"\n{i+1}. {r['optimizer']} + {r['ansatz']}")
        print(f"   Energy:      {r['energy']:.8f} Ha")
        print(f"   Error:       {r['error']:.8f} Ha ({r['error']*1000:.3f} mHa)")
        print(f"   Correlation: {r['correlation']:.8f} Ha")
        print(f"   Iterations:  {r['iterations']}")

    # Check if ANY configuration achieved target
    best = valid_results[0]
    if best['error'] < 0.001:
        print(f"\n✓✓ SUCCESS! Achieved < 1 mHa accuracy with {best['optimizer']} + {best['ansatz']}")
    elif best['error'] < 0.01:
        print(f"\n✓ GOOD! Achieved < 10 mHa accuracy with {best['optimizer']} + {best['ansatz']}")
    else:
        print(f"\n✗ FAILED to achieve 1 mHa accuracy. Best error: {best['error']*1000:.3f} mHa")
        print(f"\nPossible issues:")
        print(f"  1. Ansatz not expressive enough")
        print(f"  2. Optimizer stuck in local minimum")
        print(f"  3. Need better initialization")
        print(f"  4. Need gradient-based optimization")
else:
    print(f"\n✗ ALL configurations failed!")

print(f"\n{'='*80}")
