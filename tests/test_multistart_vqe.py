#!/usr/bin/env python3
"""
Test Multi-Start VQE
Demonstrate that running VQE 3 times and keeping best result gives reliable convergence
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from kanad.bonds.bond_factory import BondFactory
from kanad.core.atom import Atom
from kanad.utils.vqe_solver import VQESolver

print("="*80)
print("MULTI-START VQE TEST")
print("="*80)
print("\nThis demonstrates that multi-start VQE reliably finds good solutions")
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

# Create solver
solver = VQESolver(
    bond=bond,
    ansatz_type='governance',
    mapper_type='jordan_wigner',
    optimizer='COBYLA',
    max_iterations=200,
    backend='statevector',
    shots=None
)

print("Running multi-start VQE (3 restarts)...")
print("-"*80)

# Run with multi-start
result = solver.solve_with_restarts(n_restarts=3)

print("\n" + "="*80)
print("RESULTS")
print("="*80)

vqe_energy = result['energy']
correlation = vqe_energy - hf_energy
recovery = abs(correlation / expected_correlation) * 100 if abs(expected_correlation) > 1e-6 else 0

print(f"\nBest VQE Result:")
print(f"  Energy:            {vqe_energy:.8f} Ha")
print(f"  HF energy:         {hf_energy:.8f} Ha")
print(f"  Correlation:       {correlation:.8f} Ha")
print(f"  Recovery:          {recovery:.2f}%")

# Multi-start statistics
ms_info = result.get('multi_start', {})
if ms_info:
    print(f"\nMulti-Start Statistics:")
    print(f"  Restarts:          {ms_info['n_restarts']}")
    print(f"  Best attempt:      {ms_info['best_attempt']}")
    print(f"  Energy range:      {ms_info['energy_range']:.8f} Ha")
    print(f"  Energy std:        {ms_info['energy_std']:.8f} Ha")
    print(f"\n  All energies:")
    for i, e in enumerate(ms_info['all_energies'], 1):
        marker = " ⭐ BEST" if i == ms_info['best_attempt'] else ""
        print(f"    Attempt {i}: {e:.8f} Ha{marker}")

print("\n" + "="*80)
print("VERDICT")
print("="*80)

if recovery > 80:
    print(f"\n✅ EXCELLENT: {recovery:.1f}% correlation recovery!")
    print("   Multi-start VQE successfully found good solution")
elif recovery > 50:
    print(f"\n⚠️  GOOD: {recovery:.1f}% correlation recovery")
    print("   Multi-start helped but could be better")
else:
    print(f"\n❌ POOR: Only {recovery:.1f}% correlation recovery")
    print("   May need more restarts or different optimizer")

print("\n" + "="*80)
