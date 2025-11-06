#!/usr/bin/env python3
"""
Test Bond API VQE - Replicate Dashboard Behavior
This uses the exact same code path as the dashboard
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from kanad.bonds.bond_factory import BondFactory
from kanad.core.atom import Atom
from kanad.utils.vqe_solver import VQESolver

print("="*80)
print("BOND API VQE TEST - Replicating Dashboard")
print("="*80)

# Create H2 using Bond API (same as dashboard for diatomic molecules)
print("\n1. Creating H2 bond...")
atom1 = Atom('H', [0.0, 0.0, 0.0])
atom2 = Atom('H', [0.0, 0.0, 0.74])  # 0.74 Angstrom

bond = BondFactory.create_bond(
    atom1,
    atom2,
    distance=0.74,
    basis='sto-3g'
)

print(f"✓ Bond created: {bond.bond_type} bond")
print(f"  Atoms: {bond.atom_1.symbol} - {bond.atom_2.symbol}")
print(f"  Distance: {bond.distance:.4f} Å")
print(f"  Basis: sto-3g")

# Create VQE solver using bond API (same as execute_vqe for diatomic)
print("\n2. Creating VQE solver with Bond API...")

# This replicates lines 367-378 in experiment_service.py
solver = VQESolver(
    bond=bond,
    ansatz_type='governance',  # Maps to covalent for H2
    mapper_type='jordan_wigner',
    optimizer='COBYLA',
    max_iterations=100,
    backend='statevector',
    shots=None
)

print(f"✓ VQE solver created")
print(f"  Ansatz type: governance (auto-detects covalent)")
print(f"  Mapper: jordan_wigner")
print(f"  Optimizer: COBYLA")
print(f"  Max iterations: 100")

# Get HF reference for comparison
print("\n3. Getting HF reference energy...")
hf_energy_ref = bond.hamiltonian.solve_scf(max_iterations=50)[1]
print(f"✓ HF energy: {hf_energy_ref:.8f} Ha")

# Track iterations
iteration_log = []

def callback(iteration, energy, parameters):
    iteration_log.append({'iter': iteration, 'energy': energy})
    if iteration % 10 == 0:
        print(f"  Iter {iteration}: E = {energy:.8f} Ha")

# Run VQE
print("\n4. Running VQE...")
print("🚀 Starting optimization...")

result = solver.solve(callback=callback)

print(f"\n✅ VQE Completed")
print(f"  Final energy:      {result['energy']:.8f} Ha")
print(f"  HF energy:         {result.get('hf_energy', hf_energy_ref):.8f} Ha")
print(f"  Correlation:       {result.get('correlation_energy', 0):.8f} Ha")
print(f"  Converged:         {result.get('converged', False)}")
print(f"  Iterations:        {result.get('iterations', 0)}")
print(f"  Function evals:    {len(iteration_log)}")

# Analyze results
print("\n5. Analysis:")

correlation = result['energy'] - hf_energy_ref
corr_kcal = correlation * 627.5

print(f"  Correlation energy: {correlation:.8f} Ha ({corr_kcal:.4f} kcal/mol)")

# Expected FCI energy for H2/STO-3G at 0.74 Å
fci_energy_expected = -1.137284  # From our previous test
expected_correlation = fci_energy_expected - hf_energy_ref

print(f"  Expected correlation (FCI): {expected_correlation:.8f} Ha")

if abs(correlation) > 1e-6:
    recovery_percent = abs(correlation / expected_correlation) * 100
    print(f"  Correlation recovery: {recovery_percent:.2f}%")

    if recovery_percent > 90:
        print(f"  ✅ EXCELLENT: >90% recovery")
    elif recovery_percent > 50:
        print(f"  ⚠️  GOOD: 50-90% recovery")
    else:
        print(f"  ❌ POOR: <50% recovery")
else:
    print(f"  ❌ CRITICAL: No correlation energy recovered!")
    print(f"     VQE stuck at HF state")

# Check energy history
if len(iteration_log) > 0:
    energies = [e['energy'] for e in iteration_log]
    energy_range = max(energies) - min(energies)

    print(f"\n  Energy history:")
    print(f"    Initial: {energies[0]:.8f} Ha")
    print(f"    Final:   {energies[-1]:.8f} Ha")
    print(f"    Min:     {min(energies):.8f} Ha")
    print(f"    Range:   {energy_range:.8f} Ha")

    if energy_range < 1e-6:
        print(f"    ❌ Energy flat! Optimizer not moving")
    elif energy_range < 0.001:
        print(f"    ⚠️  Small variation")
    else:
        print(f"    ✓ Energy varying")

print("\n" + "="*80)
print("VERDICT")
print("="*80)

if abs(correlation) < 1e-6:
    print("\n❌ BUG CONFIRMED: VQE using Bond API gives HF energy!")
    print("   This matches your dashboard behavior.")
    print("\n   Possible causes:")
    print("   1. Bond API creates different Hamiltonian than direct API")
    print("   2. Ansatz not being applied properly with Bond API")
    print("   3. Initial state preparation issue")
    print("\n   Need to investigate Bond API code path...")
elif abs(correlation / expected_correlation) > 0.9:
    print("\n✅ VQE working correctly with Bond API!")
    print(f"   Recovering {abs(correlation / expected_correlation) * 100:.1f}% of correlation")
    print("   Dashboard issue must be elsewhere (frontend? database?)")
else:
    print("\n⚠️  VQE partially working")
    print(f"   Recovering {abs(correlation / expected_correlation) * 100:.1f}% of correlation")
    print("   Could be improved but not broken")

print("\n" + "="*80)
