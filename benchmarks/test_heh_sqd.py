#!/usr/bin/env python3
"""Test HeH+ with SQD (Subspace Quantum Diagonalization) solver"""
import time
import numpy as np
from kanad.bonds.bond_factory import BondFactory
from kanad.solvers.sqd_solver import SQDSolver

print("=" * 80)
print("SQD Solver Test: HeH Molecule")
print("=" * 80)

# Create HeH bond using BondFactory (neutral, not ionic)
# Note: BondFactory doesn't directly support charged species, so we test neutral HeH
bond = BondFactory.create_bond('He', 'H', distance=0.772, bond_type='covalent')
mol = bond.molecule

print(f"\nMolecule: HeH (neutral)")
print(f"Geometry: He-H = 0.772 Å")
print(f"Electrons: {mol.n_electrons}")

# Run SQD solver
print("\nRunning SQD solver...")
print("Configuration: subspace_dim=10, circuit_depth=3")

t0 = time.time()
solver = SQDSolver(
    bond=bond,
    subspace_dim=10,
    circuit_depth=3,
    backend='statevector',
    enable_analysis=False,
    enable_optimization=False,
    random_seed=42
)
result = solver.solve(n_states=3)
dt = time.time() - t0

E_sqd_ground = result['ground_state_energy']
E_hf = result.get('hf_energy', None)

print(f"\n{'=' * 80}")
print("RESULTS")
print("=" * 80)
if E_hf is not None:
    print(f"HF Energy:     {E_hf:.8f} Ha")
    print(f"SQD Energy:    {E_sqd_ground:.8f} Ha")
    print(f"Correlation:   {(E_sqd_ground - E_hf)*1000:+.3f} mHa")
else:
    print(f"SQD Energy:    {E_sqd_ground:.8f} Ha")
print(f"Time:          {dt:.2f} s")
print(f"Converged:     {result['converged']}")

if 'excited_state_energies' in result and len(result['excited_state_energies']) > 0:
    print(f"\nExcited States:")
    for i, E_ex in enumerate(result['excited_state_energies'], 1):
        excitation_ev = (E_ex - E_sqd_ground) * 27.2114
        print(f"  State {i}: {E_ex:.8f} Ha  (ΔE = {excitation_ev:.4f} eV)")

# Validate
print(f"\n{'=' * 80}")
if result['converged'] and E_sqd_ground < -2.0:  # Reasonable ground state energy for HeH+
    print("✓ PASS - SQD calculation successful")
    success = True
else:
    print("✗ FAIL - SQD calculation failed")
    success = False

print("=" * 80)

assert success, "SQD solver failed"
