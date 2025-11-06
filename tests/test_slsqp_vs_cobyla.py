#!/usr/bin/env python3
"""
Test SLSQP vs COBYLA
Does gradient-based optimization avoid HF trap better?
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from kanad.bonds.bond_factory import BondFactory
from kanad.core.atom import Atom
from kanad.utils.vqe_solver import VQESolver

print("="*80)
print("SLSQP vs COBYLA Comparison")
print("="*80)

# Create H2
atom1 = Atom('H', [0.0, 0.0, 0.0])
atom2 = Atom('H', [0.0, 0.0, 0.74])
bond = BondFactory.create_bond(atom1, atom2, distance=0.74, basis='sto-3g')

hf_energy = bond.hamiltonian.solve_scf(max_iterations=50)[1]
fci_energy = -1.137284

print(f"\nReference:")
print(f"  HF:  {hf_energy:.8f} Ha")
print(f"  FCI: {fci_energy:.8f} Ha")

optimizers = ['COBYLA', 'SLSQP', 'Powell']
n_runs = 5

for opt in optimizers:
    print(f"\n{'='*80}")
    print(f"Testing {opt} ({n_runs} runs)")
    print(f"{'='*80}")

    energies = []
    recoveries = []

    for i in range(1, n_runs + 1):
        solver = VQESolver(
            bond=bond,
            ansatz_type='governance',
            mapper_type='jordan_wigner',
            optimizer=opt,
            max_iterations=200,
            backend='statevector',
            shots=None
        )

        result = solver.solve()
        energy = result['energy']
        correlation = energy - hf_energy
        expected_corr = fci_energy - hf_energy
        recovery = abs(correlation / expected_corr) * 100 if abs(expected_corr) > 1e-6 else 0

        energies.append(energy)
        recoveries.append(recovery)

        status = "✅" if recovery > 50 else "❌"
        print(f"  Run {i}: {energy:.8f} Ha | {recovery:5.1f}% recovery {status}")

    print(f"\nStatistics:")
    print(f"  Mean energy:    {np.mean(energies):.8f} Ha")
    print(f"  Mean recovery:  {np.mean(recoveries):.1f}%")
    print(f"  Success rate:   {sum(1 for r in recoveries if r > 50)}/{n_runs} (>50% recovery)")

print("\n" + "="*80)
