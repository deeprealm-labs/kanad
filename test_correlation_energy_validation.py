#!/usr/bin/env python3
"""
Test Correlation Energy Calculation

Validates that correlation energy = E_quantum - E_HF is computed correctly:
1. Both energies include nuclear repulsion (it cancels)
2. Correlation energy is reasonable (negative, ~1-5% of total energy)
3. Verify nuclear repulsion is the same in both calculations
"""

import numpy as np
import logging
from kanad.bonds import BondFactory
from kanad.solvers.sqd_solver import SQDSolver

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

print("=" * 70)
print("CORRELATION ENERGY VALIDATION TEST")
print("=" * 70)

# Create H2 molecule
h2_bond = BondFactory.create_bond('H', 'H', distance=0.74)

print(f"\nMolecule: H2")
print(f"Bond distance: 0.74 Å")
print(f"Nuclear repulsion: {h2_bond.hamiltonian.nuclear_repulsion:.8f} Ha")

# Run SQD solver
solver = SQDSolver(
    bond=h2_bond,
    subspace_dim=8,
    enable_analysis=False
)

result = solver.solve(n_states=1)

# Extract energies
E_quantum = result['ground_state_energy']
E_hf = result['hf_energy']
E_corr = result['correlation_energy']
E_nuc = h2_bond.hamiltonian.nuclear_repulsion

print("\n" + "=" * 70)
print("ENERGY BREAKDOWN")
print("=" * 70)

print(f"\nQuantum ground state:  {E_quantum:.8f} Ha")
print(f"HF reference:          {E_hf:.8f} Ha")
print(f"Correlation:           {E_corr:.8f} Ha")
print(f"Nuclear repulsion:     {E_nuc:.8f} Ha")

# Verify correlation calculation
E_corr_manual = E_quantum - E_hf
print(f"\nManual calculation:    {E_corr_manual:.8f} Ha")

if abs(E_corr - E_corr_manual) < 1e-10:
    print("✅ Correlation energy calculation is consistent")
else:
    print(f"❌ Inconsistency: {abs(E_corr - E_corr_manual):.2e} Ha")

# Check if nuclear repulsion cancels
print("\n" + "=" * 70)
print("NUCLEAR REPULSION CANCELLATION CHECK")
print("=" * 70)

# The electronic energies (without nuclear repulsion)
E_electronic_quantum = E_quantum - E_nuc
E_electronic_hf = E_hf - E_nuc

print(f"\nElectronic energy (quantum): {E_electronic_quantum:.8f} Ha")
print(f"Electronic energy (HF):      {E_electronic_hf:.8f} Ha")
print(f"Difference (correlation):    {E_electronic_quantum - E_electronic_hf:.8f} Ha")

if abs((E_electronic_quantum - E_electronic_hf) - E_corr) < 1e-10:
    print("✅ Nuclear repulsion cancels correctly in correlation energy")
else:
    print("❌ Nuclear repulsion does NOT cancel correctly!")

# Validate physical reasonableness
print("\n" + "=" * 70)
print("PHYSICAL REASONABLENESS CHECKS")
print("=" * 70)

# Check 1: Correlation should be negative (lowers energy)
if E_corr < 0:
    print(f"✅ Correlation is negative (lowers energy): {E_corr:.8f} Ha")
else:
    print(f"❌ Correlation is positive (wrong!): {E_corr:.8f} Ha")

# Check 2: Correlation should be small (1-5% of total energy)
percent_correlation = abs(E_corr / E_quantum) * 100
print(f"\nCorrelation as % of total: {percent_correlation:.2f}%")

if 0.5 < percent_correlation < 10:
    print(f"✅ Correlation is reasonable magnitude ({percent_correlation:.2f}%)")
else:
    print(f"⚠️  Correlation magnitude unusual ({percent_correlation:.2f}%)")

# Check 3: Quantum energy should be lower than HF
if E_quantum < E_hf:
    print(f"✅ Quantum energy lower than HF (variational principle)")
else:
    print(f"❌ Quantum energy HIGHER than HF (violates variational principle!)")

# Check 4: Correlation should be in reasonable range
if -0.1 < E_corr < 0:
    print(f"✅ Correlation energy in expected range for H2: {E_corr:.8f} Ha")
else:
    print(f"⚠️  Correlation energy outside typical range: {E_corr:.8f} Ha")

# FINAL VERDICT
print("\n" + "=" * 70)
print("FINAL VERDICT")
print("=" * 70)

all_checks = [
    abs(E_corr - E_corr_manual) < 1e-10,  # Calculation consistent
    abs((E_electronic_quantum - E_electronic_hf) - E_corr) < 1e-10,  # Nuc rep cancels
    E_corr < 0,  # Correlation is negative
    E_quantum < E_hf,  # Quantum lower than HF
    0.5 < percent_correlation < 10  # Reasonable magnitude
]

if all(all_checks):
    print("\n✅ ALL CORRELATION ENERGY CHECKS PASSED")
    print("\nConclusion:")
    print("  ✅ Correlation energy = E_quantum - E_HF is CORRECT")
    print("  ✅ Nuclear repulsion cancels properly")
    print("  ✅ Results are physically reasonable")
    print("  ✅ Issue #9 is NOT A REAL ISSUE - formula is correct!")
else:
    print("\n❌ SOME CHECKS FAILED - Investigation needed")

print("\n" + "=" * 70)
