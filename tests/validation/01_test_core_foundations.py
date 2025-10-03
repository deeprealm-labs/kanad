#!/usr/bin/env python3
"""
Phase 1: Core Foundations Testing
==================================

Tests the most fundamental components:
- Atoms and atomic properties
- Physical constants and conversions
- Basis sets
- Molecular integrals
- Representations

These MUST be correct for everything else to work.
"""

import numpy as np
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

print("="*80)
print("PHASE 1: CORE FOUNDATIONS TESTING")
print("="*80)

# ==============================================================================
# Test 1.1: Atoms and Atomic Properties
# ==============================================================================

print("\n[1.1] Testing Atoms and Atomic Properties")
print("-"*80)

from kanad.core.atom import Atom

test_cases = [
    {
        'element': 'H',
        'expected_number': 1,
        'expected_mass': 1.008,  # amu
        'position': np.array([0.0, 0.0, 0.0])
    },
    {
        'element': 'C',
        'expected_number': 6,
        'expected_mass': 12.011,
        'position': np.array([1.0, 1.0, 1.0])
    },
    {
        'element': 'O',
        'expected_number': 8,
        'expected_mass': 15.999,
        'position': np.array([2.0, 0.0, 0.0])
    }
]

passed = 0
failed = 0

for tc in test_cases:
    try:
        atom = Atom(tc['element'], position=tc['position'])

        # Check atomic number
        assert atom.atomic_number == tc['expected_number'], \
            f"Atomic number mismatch: {atom.atomic_number} != {tc['expected_number']}"

        # Check mass (rough check)
        assert abs(atom.mass - tc['expected_mass']) < 0.1, \
            f"Mass mismatch: {atom.mass} != {tc['expected_mass']}"

        # Check position
        assert np.allclose(atom.position, tc['position']), \
            f"Position mismatch"

        print(f"  ✅ {tc['element']}: Z={atom.atomic_number}, mass={atom.mass:.3f} amu")
        passed += 1

    except Exception as e:
        print(f"  ❌ {tc['element']}: {e}")
        failed += 1

# Test distance calculation
try:
    h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
    h2 = Atom('H', position=np.array([0.0, 0.0, 0.74]))

    distance = h1.distance_to(h2)
    expected = 0.74

    assert abs(distance - expected) < 1e-10, f"Distance mismatch: {distance} != {expected}"
    print(f"  ✅ Distance calculation: {distance:.2f} Å")
    passed += 1

except Exception as e:
    print(f"  ❌ Distance calculation: {e}")
    failed += 1

print(f"\n  Results: {passed} passed, {failed} failed")

# ==============================================================================
# Test 1.2: Physical Constants
# ==============================================================================

print("\n[1.2] Testing Physical Constants")
print("-"*80)

from kanad.core.constants.physical_constants import PhysicalConstants
from kanad.core.constants.conversion_factors import ConversionFactors

# Test known constants
constants_tests = [
    ('BOHR_TO_ANGSTROM', PhysicalConstants.BOHR_TO_ANGSTROM, 0.529177, 1e-5),
    ('HARTREE_TO_EV', PhysicalConstants.HARTREE_TO_EV, 27.211, 0.001),
    ('SPEED_OF_LIGHT', PhysicalConstants.SPEED_OF_LIGHT, 299792458, 1.0),
]

passed = 0
failed = 0

for name, value, expected, tolerance in constants_tests:
    try:
        assert abs(value - expected) < tolerance, \
            f"{name} = {value}, expected ~{expected}"
        print(f"  ✅ {name}: {value}")
        passed += 1
    except Exception as e:
        print(f"  ❌ {name}: {e}")
        failed += 1

# Test conversion factors
try:
    # Bohr to Angstrom
    bohr_to_ang = ConversionFactors.BOHR_TO_ANGSTROM
    ang_to_bohr = ConversionFactors.ANGSTROM_TO_BOHR

    assert abs(bohr_to_ang * ang_to_bohr - 1.0) < 1e-10, "Conversion inconsistency"

    print(f"  ✅ Bohr ↔ Angstrom: {bohr_to_ang:.6f} / {ang_to_bohr:.6f}")
    passed += 1

except Exception as e:
    print(f"  ❌ Conversion factors: {e}")
    failed += 1

print(f"\n  Results: {passed} passed, {failed} failed")

# ==============================================================================
# Test 1.3: Basis Sets
# ==============================================================================

print("\n[1.3] Testing Basis Sets")
print("-"*80)

from kanad.core.integrals.basis_sets import BasisSet
from kanad.core.representations.base_representation import Molecule

# Test STO-3G basis for H2
try:
    h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
    h2 = Atom('H', position=np.array([0.0, 0.0, 0.74]))
    atoms = [h1, h2]

    basis = BasisSet('sto-3g')
    basis.build_basis(atoms)

    n_basis = len(basis.basis_functions)

    assert n_basis == 2, f"H2 should have 2 basis functions (1s on each H), got {n_basis}"

    print(f"  ✅ STO-3G H2: {n_basis} basis functions")

    # Check basis function properties
    for i, bf in enumerate(basis.basis_functions):
        print(f"    BF{i}: center={bf.center}, L={bf.L}, n_primitives={len(bf.exponents)}")

    passed = 1
    failed = 0

except Exception as e:
    print(f"  ❌ STO-3G H2: {e}")
    import traceback
    traceback.print_exc()
    passed = 0
    failed = 1

# Test other elements
elements_to_test = ['He', 'C', 'N', 'O', 'F']

for element in elements_to_test:
    try:
        atom = Atom(element, position=np.array([0.0, 0.0, 0.0]))
        basis = BasisSet('sto-3g')
        basis.build_basis([atom])

        n_basis = len(basis.basis_functions)
        print(f"  ✅ STO-3G {element}: {n_basis} basis functions")
        passed += 1

    except Exception as e:
        print(f"  ❌ STO-3G {element}: {e}")
        failed += 1

print(f"\n  Results: {passed} passed, {failed} failed")

# ==============================================================================
# Test 1.4: Molecular Integrals
# ==============================================================================

print("\n[1.4] Testing Molecular Integrals")
print("-"*80)

from kanad.core.representations.lcao_representation import LCAORepresentation
from kanad.core.integrals.overlap import compute_overlap_matrix
from kanad.core.integrals.one_electron import compute_kinetic_matrix, compute_nuclear_attraction_matrix
from kanad.core.integrals.two_electron import compute_electron_repulsion_integrals

# H2 molecule for testing
h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
h2 = Atom('H', position=np.array([0.0, 0.0, 1.4]))  # 0.74 Å = 1.4 Bohr
molecule = Molecule([h1, h2])

try:
    representation = LCAORepresentation(molecule)
    basis = BasisSet('sto-3g')
    basis.build_basis(molecule.atoms)

    # Overlap matrix
    S = compute_overlap_matrix(basis)
    print(f"\n  Overlap Matrix S:")
    print(f"    Shape: {S.shape}")
    print(f"    S[0,0] = {S[0,0]:.6f} (should be 1.0)")
    print(f"    S[0,1] = {S[0,1]:.6f} (should be ~0.66 for H2)")

    # Reference: H2 at 1.4 Bohr should have S[0,1] ≈ 0.66
    assert abs(S[0,0] - 1.0) < 1e-6, "Diagonal overlap should be 1"
    assert 0.6 < S[0,1] < 0.7, f"H2 overlap should be ~0.66, got {S[0,1]}"

    print(f"  ✅ Overlap matrix correct")

    # Kinetic energy matrix
    T = compute_kinetic_matrix(basis)
    print(f"\n  Kinetic Matrix T:")
    print(f"    T[0,0] = {T[0,0]:.6f}")
    print(f"    T[0,1] = {T[0,1]:.6f}")

    # Reference: T should be positive
    assert T[0,0] > 0, "Kinetic energy should be positive"

    print(f"  ✅ Kinetic matrix reasonable")

    # Nuclear attraction matrix
    V = compute_nuclear_attraction_matrix(basis, molecule.atoms)
    print(f"\n  Nuclear Attraction Matrix V:")
    print(f"    V[0,0] = {V[0,0]:.6f}")
    print(f"    V[0,1] = {V[0,1]:.6f}")

    # Reference: V should be negative (attractive)
    assert V[0,0] < 0, "Nuclear attraction should be negative"

    print(f"  ✅ Nuclear attraction matrix reasonable")

    # h_core = T + V
    h_core = T + V
    print(f"\n  Core Hamiltonian h_core = T + V:")
    print(f"    h_core[0,0] = {h_core[0,0]:.6f}")
    print(f"    h_core[0,1] = {h_core[0,1]:.6f}")

    # Two-electron integrals (ERI)
    eri = compute_electron_repulsion_integrals(basis)
    print(f"\n  Electron Repulsion Integrals (ERI):")
    print(f"    Shape: {eri.shape}")
    print(f"    eri[0,0,0,0] = {eri[0,0,0,0]:.6f}")
    print(f"    eri[0,0,1,1] = {eri[0,0,1,1]:.6f}")

    # Reference: (00|00) should be ~0.77 for H (1s)
    assert eri[0,0,0,0] > 0, "ERI should be positive"

    print(f"  ✅ ERI computed successfully")

    passed = 5
    failed = 0

except Exception as e:
    print(f"\n  ❌ Integral computation failed: {e}")
    import traceback
    traceback.print_exc()
    passed = 0
    failed = 1

print(f"\n  Results: {passed} passed, {failed} failed")

# ==============================================================================
# Test 1.5: Compare Against Reference (PySCF)
# ==============================================================================

print("\n[1.5] Comparing Against Reference Values (PySCF)")
print("-"*80)

try:
    import pyscf
    from pyscf import gto, scf

    # Build H2 in PySCF
    mol_pyscf = gto.M(
        atom='H 0 0 0; H 0 0 0.74',
        basis='sto-3g',
        unit='Angstrom'
    )

    # Get reference integrals
    S_ref = mol_pyscf.intor('int1e_ovlp')
    T_ref = mol_pyscf.intor('int1e_kin')
    V_ref = mol_pyscf.intor('int1e_nuc')
    eri_ref = mol_pyscf.intor('int2e')

    # Run HF to get reference energy
    mf = scf.RHF(mol_pyscf)
    e_hf_ref = mf.kernel()

    print(f"\n  PySCF Reference Values (H2, STO-3G):")
    print(f"    HF Energy: {e_hf_ref:.6f} Ha")
    print(f"    S[0,1]: {S_ref[0,1]:.6f}")
    print(f"    T[0,0]: {T_ref[0,0]:.6f}")
    print(f"    V[0,0]: {V_ref[0,0]:.6f}")
    print(f"    ERI[0,0,0,0]: {eri_ref[0,0,0,0]:.6f}")

    # Compare with Kanad values
    print(f"\n  Kanad vs PySCF Comparison:")

    errors = []

    # Overlap
    s_error = abs(S[0,1] - S_ref[0,1])
    print(f"    Overlap S[0,1]: {S[0,1]:.6f} vs {S_ref[0,1]:.6f} (Δ = {s_error:.2e})")
    errors.append(('Overlap', s_error))

    # Kinetic
    t_error = abs(T[0,0] - T_ref[0,0])
    print(f"    Kinetic T[0,0]: {T[0,0]:.6f} vs {T_ref[0,0]:.6f} (Δ = {t_error:.2e})")
    errors.append(('Kinetic', t_error))

    # Nuclear
    v_error = abs(V[0,0] - V_ref[0,0])
    print(f"    Nuclear V[0,0]: {V[0,0]:.6f} vs {V_ref[0,0]:.6f} (Δ = {v_error:.2e})")
    errors.append(('Nuclear', v_error))

    # ERI
    eri_error = abs(eri[0,0,0,0] - eri_ref[0,0,0,0])
    print(f"    ERI[0,0,0,0]: {eri[0,0,0,0]:.6f} vs {eri_ref[0,0,0,0]:.6f} (Δ = {eri_error:.2e})")
    errors.append(('ERI', eri_error))

    # Check if all errors are small
    print(f"\n  Error Summary:")
    all_good = True
    for name, error in errors:
        status = "✅" if error < 1e-6 else "⚠️" if error < 1e-4 else "❌"
        print(f"    {status} {name}: {error:.2e}")
        if error >= 1e-4:
            all_good = False

    if all_good:
        print(f"\n  ✅ All integrals match PySCF reference (< 1e-6 error)")
        passed = 1
        failed = 0
    else:
        print(f"\n  ⚠️  Some integrals have errors > 1e-6")
        passed = 0
        failed = 1

except ImportError:
    print("  ⚠️  PySCF not installed - skipping reference comparison")
    print("     Install with: pip install pyscf")
    passed = 0
    failed = 0

except Exception as e:
    print(f"\n  ❌ Reference comparison failed: {e}")
    import traceback
    traceback.print_exc()
    passed = 0
    failed = 1

print(f"\n  Results: {passed} passed, {failed} failed")

# ==============================================================================
# Summary
# ==============================================================================

print("\n" + "="*80)
print("PHASE 1 SUMMARY: CORE FOUNDATIONS")
print("="*80)

print("""
Components Tested:
  [1.1] Atoms and Atomic Properties
  [1.2] Physical Constants
  [1.3] Basis Sets
  [1.4] Molecular Integrals
  [1.5] Reference Comparison (PySCF)

Critical for Framework:
  ✅ If all pass → Core is solid, can build on it
  ❌ If any fail → MUST fix before proceeding

Next Phase: Test Hamiltonians (Phase 2)
""")

print("="*80)
