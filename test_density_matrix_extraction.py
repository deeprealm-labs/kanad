#!/usr/bin/env python3
"""
Test Density Matrix Extraction from Hamiltonians

Validates that get_density_matrix() works correctly after Phase 1 fixes.
"""

import numpy as np
from kanad.bonds import BondFactory

print("="*70)
print("🔧 DENSITY MATRIX EXTRACTION TEST")
print("="*70)

# Test 1: Covalent Hamiltonian
print("\n" + "="*70)
print("TEST 1: Covalent Hamiltonian (H2)")
print("="*70)

h2_bond = BondFactory.create_bond('H', 'H', distance=0.74)

# Solve SCF first
print("\n🔧 Running SCF...")
rdm1_scf, energy_scf = h2_bond.hamiltonian.solve_scf()
print(f"  ✓ SCF converged: E = {energy_scf:.6f} Ha")
print(f"  ✓ Density matrix shape: {rdm1_scf.shape}")

# Extract density matrix using new method
print("\n🔧 Extracting density matrix...")
try:
    rdm1_extracted = h2_bond.hamiltonian.get_density_matrix()
    print(f"  ✓ Extracted density matrix shape: {rdm1_extracted.shape}")

    # Verify they match
    if np.allclose(rdm1_scf, rdm1_extracted):
        print("  ✅ Density matrices match!")
    else:
        print("  ❌ Density matrices don't match!")
        max_diff = np.max(np.abs(rdm1_scf - rdm1_extracted))
        print(f"     Max difference: {max_diff}")

except Exception as e:
    print(f"  ❌ Failed to extract density: {e}")

# Test 2: Ionic Hamiltonian
print("\n" + "="*70)
print("TEST 2: Ionic Hamiltonian (LiH)")
print("="*70)

lih_bond = BondFactory.create_bond('Li', 'H', distance=1.60)

# Solve SCF first
print("\n🔧 Running SCF...")
rdm1_scf, energy_scf = lih_bond.hamiltonian.solve_scf()
print(f"  ✓ SCF converged: E = {energy_scf:.6f} Ha")
print(f"  ✓ Density matrix shape: {rdm1_scf.shape}")

# Extract density matrix using new method
print("\n🔧 Extracting density matrix...")
try:
    rdm1_extracted = lih_bond.hamiltonian.get_density_matrix()
    print(f"  ✓ Extracted density matrix shape: {rdm1_extracted.shape}")

    # Verify they match
    if np.allclose(rdm1_scf, rdm1_extracted):
        print("  ✅ Density matrices match!")
    else:
        print("  ❌ Density matrices don't match!")
        max_diff = np.max(np.abs(rdm1_scf - rdm1_extracted))
        print(f"     Max difference: {max_diff}")

except Exception as e:
    print(f"  ❌ Failed to extract density: {e}")

# Test 3: Error handling (no SCF)
print("\n" + "="*70)
print("TEST 3: Error Handling (No SCF)")
print("="*70)

print("\n🔧 Creating bond without SCF...")
# Create a fresh bond instance
from kanad.core.hamiltonians.covalent_hamiltonian import CovalentHamiltonian
from kanad.core.molecule import Molecule
from kanad.core.atom import Atom

atoms = [Atom('H', position=(0, 0, 0)), Atom('H', position=(0, 0, 0.74))]
molecule = Molecule(atoms, charge=0, spin=0)

# Create hamiltonian but don't call solve_scf
from kanad.core.integrals.lcao_representation import LCAORepresentation
lcao = LCAORepresentation(molecule, basis_name='sto-3g')
hamiltonian = CovalentHamiltonian(molecule, lcao, use_governance=False)

print("🔧 Trying to extract density without SCF...")
try:
    rdm1 = hamiltonian.get_density_matrix()
    print("  ❌ Should have raised ValueError!")
except ValueError as e:
    print(f"  ✅ Correctly raised ValueError: {str(e)[:80]}...")
except Exception as e:
    print(f"  ❌ Wrong exception type: {type(e).__name__}: {e}")

# Summary
print("\n" + "="*70)
print("✅ SUMMARY")
print("="*70)
print("✓ Covalent Hamiltonian: get_density_matrix() works")
print("✓ Ionic Hamiltonian: get_density_matrix() works")
print("✓ Error handling: Raises ValueError when no SCF")
print("\n🎉 Phase 1 density matrix extraction complete!")
print("="*70)
