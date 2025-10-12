#!/usr/bin/env python3
"""Test LiH with Excited States solver (CIS method)"""
import time
import numpy as np
from kanad.core.atom import Atom
from kanad.core.molecule import Molecule
from kanad.solvers.excited_states_solver import ExcitedStatesSolver
from kanad.bonds.bond_factory import BondFactory

print("=" * 80)
print("Excited States Solver Test: LiH Molecule")
print("=" * 80)

# Create LiH bond using BondFactory
bond = BondFactory.create_bond('Li', 'H', distance=1.595, bond_type='covalent')
mol = bond.molecule

print(f"\nMolecule: LiH")
print(f"Geometry: Li-H = 1.595 Å")
print(f"Electrons: {mol.n_electrons}")

# Get HF reference
print("\nComputing HF reference...")
density_matrix, hf_energy = bond.hamiltonian.solve_scf(
    max_iterations=100,
    conv_tol=1e-8
)
print(f"HF Energy: {hf_energy:.8f} Ha")

# Run Excited States solver
print("\nRunning CIS excited states solver...")
print("Configuration: method='cis', n_states=5")

t0 = time.time()
solver = ExcitedStatesSolver(
    bond=bond,
    method='cis',
    n_states=5,
    enable_analysis=False,
    enable_optimization=False
)
result = solver.solve()
dt = time.time() - t0

E_ground = result['ground_state_energy']
excitations_ev = result.get('excitation_energies_ev', [])

print(f"\n{'=' * 80}")
print("RESULTS")
print("=" * 80)
print(f"Ground State:  {E_ground:.8f} Ha")
print(f"Method:        {result.get('method', 'CIS')}")
print(f"Time:          {dt:.2f} s")
print(f"Converged:     {result.get('converged', False)}")

if len(excitations_ev) > 0:
    print(f"\nExcited States ({len(excitations_ev)} found):")
    oscillators = result.get('oscillator_strengths', [0]*len(excitations_ev))
    transitions = result.get('dominant_transitions', ['?']*len(excitations_ev))

    for i, (E_ev, f, trans) in enumerate(zip(excitations_ev, oscillators, transitions), 1):
        wavelength_nm = 1240 / E_ev if E_ev > 0 else 0
        print(f"  State {i}: {E_ev:8.4f} eV ({wavelength_nm:6.1f} nm)  f={f:.4f}  {trans}")
else:
    print("\nNo excited states found")

# Validate
print(f"\n{'=' * 80}")
if result.get('converged', False):
    print("✓ PASS - CIS calculation converged")
    success = True
else:
    print("✗ FAIL - CIS calculation did not converge")
    success = False

print("=" * 80)

assert success, "Excited states solver failed"
