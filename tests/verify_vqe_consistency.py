#!/usr/bin/env python3
"""
Verify VQE Consistency After Fix
Run VQE 10 times and check that all runs achieve good correlation recovery
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from kanad.bonds.bond_factory import BondFactory
from kanad.core.atom import Atom
from kanad.utils.vqe_solver import VQESolver

print("="*80)
print("VQE CONSISTENCY TEST - 10 CONSECUTIVE RUNS")
print("="*80)
print("\nThis test verifies the initialization fix improves consistency")
print("Expected: >90% of runs should achieve >80% correlation recovery")
print()

# Create H2 molecule
atom1 = Atom('H', [0.0, 0.0, 0.0])
atom2 = Atom('H', [0.0, 0.0, 0.74])
bond = BondFactory.create_bond(atom1, atom2, distance=0.74, basis='sto-3g')

# Get reference energies
hf_energy = bond.hamiltonian.solve_scf(max_iterations=50)[1]
fci_energy = -1.137284  # Known FCI energy for H2/STO-3G @ 0.74 Å
expected_correlation = fci_energy - hf_energy

print(f"Reference energies:")
print(f"  HF:  {hf_energy:.8f} Ha")
print(f"  FCI: {fci_energy:.8f} Ha")
print(f"  Expected correlation: {expected_correlation:.8f} Ha")
print()

# Run VQE 10 times
n_runs = 10
results = []

print(f"Running VQE {n_runs} times with COBYLA optimizer...")
print("-"*80)

for run in range(1, n_runs + 1):
    # Create fresh solver for each run
    solver = VQESolver(
        bond=bond,
        ansatz_type='governance',
        mapper_type='jordan_wigner',
        optimizer='COBYLA',
        max_iterations=200,
        backend='statevector',
        shots=None
    )

    # Run VQE
    result = solver.solve()

    vqe_energy = result['energy']
    correlation = vqe_energy - hf_energy
    recovery = abs(correlation / expected_correlation) * 100 if abs(expected_correlation) > 1e-6 else 0

    results.append({
        'run': run,
        'energy': vqe_energy,
        'correlation': correlation,
        'recovery': recovery,
        'iterations': result.get('iterations', 0)
    })

    # Print result
    status = "✅" if recovery > 80 else "⚠️" if recovery > 50 else "❌"
    print(f"Run {run:2d}: E = {vqe_energy:.8f} Ha | Corr = {correlation:.8f} Ha | Recovery = {recovery:5.1f}% {status}")

print("-"*80)

# Analysis
print("\nSTATISTICAL ANALYSIS:")
print("="*80)

energies = [r['energy'] for r in results]
correlations = [r['correlation'] for r in results]
recoveries = [r['recovery'] for r in results]

print(f"\nEnergy Statistics:")
print(f"  Mean:   {np.mean(energies):.8f} Ha")
print(f"  Std:    {np.std(energies):.8f} Ha")
print(f"  Min:    {np.min(energies):.8f} Ha")
print(f"  Max:    {np.max(energies):.8f} Ha")

print(f"\nCorrelation Recovery Statistics:")
print(f"  Mean:   {np.mean(recoveries):.2f}%")
print(f"  Std:    {np.std(recoveries):.2f}%")
print(f"  Min:    {np.min(recoveries):.2f}%")
print(f"  Max:    {np.max(recoveries):.2f}%")

# Success criteria
excellent_runs = sum(1 for r in recoveries if r > 80)
good_runs = sum(1 for r in recoveries if 50 < r <= 80)
poor_runs = sum(1 for r in recoveries if r <= 50)

print(f"\nSuccess Rate:")
print(f"  Excellent (>80%): {excellent_runs}/{n_runs} ({excellent_runs/n_runs*100:.1f}%)")
print(f"  Good (50-80%):    {good_runs}/{n_runs} ({good_runs/n_runs*100:.1f}%)")
print(f"  Poor (<50%):      {poor_runs}/{n_runs} ({poor_runs/n_runs*100:.1f}%)")

print("\n" + "="*80)
print("VERDICT")
print("="*80)

if excellent_runs >= 9:
    print("\n✅ EXCELLENT: Fix is working perfectly!")
    print(f"   {excellent_runs}/{n_runs} runs achieved >80% recovery")
    print("   VQE is now highly consistent")
elif excellent_runs >= 7:
    print("\n✅ GOOD: Fix significantly improved consistency")
    print(f"   {excellent_runs}/{n_runs} runs achieved >80% recovery")
    print("   Much better than before (was ~1-2/10)")
elif excellent_runs >= 5:
    print("\n⚠️  IMPROVED: Some improvement but not ideal")
    print(f"   {excellent_runs}/{n_runs} runs achieved >80% recovery")
    print("   May need multi-start approach for better reliability")
else:
    print("\n❌ INSUFFICIENT: Fix did not improve consistency enough")
    print(f"   Only {excellent_runs}/{n_runs} runs achieved >80% recovery")
    print("   Need to investigate further or implement multi-start")

print("\n" + "="*80)
