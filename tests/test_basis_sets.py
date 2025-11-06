#!/usr/bin/env python3
"""Test different basis sets to see which ones PySCF supports."""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from api.services.experiment_service import create_molecule_from_config

# Test basis sets
BASIS_SETS = [
    'sto-3g',      # Minimal - KNOWN TO WORK
    '3-21g',       # Split-valence
    '6-31g',       # Popular split-valence
    '6-31g*',      # With polarization (d)
    '6-31g**',     # With polarization (d,p)
    '6-311g',      # Triple-zeta
    'cc-pvdz',     # Correlation-consistent double-zeta
    'cc-pvtz',     # Correlation-consistent triple-zeta
    'cc-pvqz',     # Correlation-consistent quadruple-zeta
    'aug-cc-pvdz', # Augmented cc-pVDZ
]

print("=" * 80)
print("TESTING BASIS SET AVAILABILITY IN PYSCF")
print("=" * 80)
print("\nTest molecule: H2")
print()

results = []

for basis in BASIS_SETS:
    print(f"Testing {basis:15s}... ", end='', flush=True)
    try:
        molecule_config = {
            'smiles': '[H][H]',
            'basis': basis,
            'charge': 0,
            'multiplicity': 1
        }
        mol = create_molecule_from_config(molecule_config)
        ham = mol.hamiltonian
        scf_results, hf_energy = ham.solve_scf(max_iterations=50, conv_tol=1e-6)

        print(f"✅ WORKS   HF Energy: {hf_energy:.8f} Ha, {ham.n_orbitals} orbitals")
        results.append((basis, True, hf_energy, ham.n_orbitals))
    except Exception as e:
        error_msg = str(e)
        if len(error_msg) > 50:
            error_msg = error_msg[:50] + "..."
        print(f"❌ FAILED  {error_msg}")
        results.append((basis, False, None, None))

print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)
print()
print(f"{'Basis Set':<15s} {'Status':<10s} {'HF Energy (Ha)':<20s} {'N Orbitals'}")
print("-" * 80)

working = []
failed = []

for basis, works, energy, n_orb in results:
    status = "✅ WORKS" if works else "❌ FAILED"
    energy_str = f"{energy:.8f}" if works else "N/A"
    orb_str = str(n_orb) if works else "N/A"
    print(f"{basis:<15s} {status:<10s} {energy_str:<20s} {orb_str}")

    if works:
        working.append(basis)
    else:
        failed.append(basis)

print("\n" + "=" * 80)
print(f"WORKING: {len(working)}/{len(BASIS_SETS)} basis sets")
print("=" * 80)
for basis in working:
    print(f"  ✅ {basis}")

if failed:
    print("\n" + "=" * 80)
    print(f"FAILED: {len(failed)}/{len(BASIS_SETS)} basis sets")
    print("=" * 80)
    for basis in failed:
        print(f"  ❌ {basis}")

print("\n" + "=" * 80)
print("RECOMMENDATION:")
print("=" * 80)
if '6-31g' in working:
    print("✅ 6-31g is AVAILABLE - ready for comprehensive testing!")
if 'cc-pvdz' in working:
    print("✅ cc-pVDZ is AVAILABLE - ready for high-accuracy benchmarks!")
if len(working) >= 5:
    print(f"✅ {len(working)} basis sets available - excellent coverage!")
else:
    print(f"⚠️  Only {len(working)} basis sets available - may need to install more basis set data")

print()
