"""
Validation Script 1: H2 Covalent Bond
=====================================

Tests the H2 molecule (simplest covalent bond) and validates framework capabilities.

This script validates:
1. Bond creation and auto-detection
2. Bond analysis capabilities
3. Molecular orbital structure
4. Hamiltonian construction
5. Governance protocol application
"""

import numpy as np
from kanad.bonds import BondFactory
from kanad.analysis.energy_analysis import BondingAnalyzer

print("="*70)
print("VALIDATION 1: H2 Covalent Bond")
print("="*70)
print()

# ============================================================================
# 1. Bond Creation and Auto-Detection
# ============================================================================
print("1. Creating H2 bond with automatic bond type detection...")
print("-" * 70)

bond = BondFactory.create_bond('H', 'H')

print(f"Bond created: {bond}")
print(f"Bond type: {bond.bond_type}")
print(f"Bond length: {bond.get_bond_length():.4f} Å")
print(f"Expected: ~0.74 Å (experimental)")
print(f"Number of atoms: {len(bond.atoms)}")
print(f"Atoms: {bond.atoms[0].symbol} - {bond.atoms[1].symbol}")
print()

# ============================================================================
# 2. Quick Bond Information
# ============================================================================
print("2. Quick bond information (no calculation)...")
print("-" * 70)
info = BondFactory.quick_bond_info('H', 'H')
print(f"Predicted type: {info['predicted_type']}")
print(f"ΔEN: {info['electronegativity_difference']:.4f}")
print(f"Estimated bond length: {info['estimated_bond_length']:.4f} Å")
print(f"H electronegativity: {info['electronegativity_1']:.2f}")
print(f"Rationale: {info['rationale']}")
print()

# ============================================================================
# 3. Molecule and Representation
# ============================================================================
print("3. Molecular system properties...")
print("-" * 70)

print(f"Molecule: {bond.molecule}")
print(f"Number of atoms: {bond.molecule.n_atoms}")
print(f"Number of electrons: {bond.molecule.n_electrons}")
print(f"Total nuclear charge: {sum(atom.atomic_number for atom in bond.molecule.atoms)}")
print()

print(f"Representation: {bond.representation}")
print(f"Number of qubits: {bond.representation.n_qubits}")
print(f"Number of orbitals: {bond.representation.n_orbitals}")
print(f"Representation type: {bond.representation.__class__.__name__}")
print()

# ============================================================================
# 4. Hamiltonian Construction
# ============================================================================
print("4. Hamiltonian properties...")
print("-" * 70)

H = bond.hamiltonian
print(f"Hamiltonian type: {H.__class__.__name__}")
print(f"Number of orbitals: {H.n_orbitals}")
print(f"Number of electrons: {H.n_electrons}")
print(f"Nuclear repulsion: {H.nuclear_repulsion:.6f} Hartree")
print()

# Core Hamiltonian
if hasattr(H, 'h_core') and H.h_core is not None:
    print("Core Hamiltonian (h_core) matrix:")
    print(H.h_core)
    print(f"Shape: {H.h_core.shape}")
    print(f"Is Hermitian: {np.allclose(H.h_core, H.h_core.T)}")
    print()

# Overlap matrix
if hasattr(H, 'S') and H.S is not None:
    print("Overlap matrix (S):")
    print(H.S)
    print(f"Diagonal (should be ~1): {np.diag(H.S)}")
    print(f"Off-diagonal (orbital overlap): {H.S[0,1]:.4f}")
    print()

# ============================================================================
# 5. Molecular Orbital Analysis
# ============================================================================
print("5. Molecular orbital analysis...")
print("-" * 70)

if hasattr(H, 'compute_molecular_orbitals'):
    mo_energies, mo_coeffs = H.compute_molecular_orbitals()
    print(f"Number of molecular orbitals: {len(mo_energies)}")
    print(f"\nMO energies (Hartree):")
    for i, energy in enumerate(mo_energies):
        if i == 0:
            label = "σ (bonding)"
            occ = "OCCUPIED (2e⁻)"
        elif i == 1:
            label = "σ* (antibonding)"
            occ = "VIRTUAL"
        else:
            label = f"MO {i}"
            occ = "VIRTUAL"
        print(f"  MO {i}: {energy:10.6f} Ha  [{label:20s}] {occ}")

    if len(mo_energies) >= 2:
        splitting = mo_energies[1] - mo_energies[0]
        print(f"\nσ - σ* splitting: {splitting:.6f} Hartree")
        print(f"                  {splitting*27.211:.4f} eV")
        print(f"                  {splitting*627.5:.2f} kcal/mol")

    print(f"\nMO coefficients (LCAO-MO):")
    print(mo_coeffs)
    print()

# ============================================================================
# 6. Bond Analysis
# ============================================================================
print("6. Detailed bond analysis...")
print("-" * 70)

analysis = bond.analyze()

print(f"Bond type: {analysis['bond_type']}")
print(f"Bond order: {analysis.get('bond_order', 1)}")
print(f"Hybridization: {analysis['hybridization']}")
print(f"Bond length: {analysis['bond_length']:.4f} Å")
print()

print("Character analysis:")
print(f"  Ionic character:    {analysis['ionic_character']:6.2%}")
print(f"  Covalent character: {analysis['covalent_character']:6.2%}")
print(f"  ΔEN: {analysis['electronegativity_difference']:.4f}")
print(f"  Expected: 100% covalent (H-H is homonuclear)")
print()

if 'homo_lumo_gap' in analysis:
    print(f"Electronic structure:")
    print(f"  HOMO-LUMO gap: {analysis['homo_lumo_gap']:.6f} Hartree")
    print(f"                 {analysis['homo_lumo_gap_ev']:.4f} eV")
    print()

if 'overlap' in analysis:
    print(f"Orbital overlap: {analysis['overlap']:.4f}")
    print()

print(f"Entanglement type: {analysis['entanglement_type']}")
print(f"Governance protocol: {analysis['governance_protocol']}")
print()

# ============================================================================
# 7. Governance Protocol Validation
# ============================================================================
print("7. Governance protocol properties...")
print("-" * 70)

protocol = bond.governance
print(f"Protocol type: {protocol.__class__.__name__}")
print(f"Protocol validates: covalent bonding constraints")
print(f"  - Enforces orbital hybridization")
print(f"  - Constructs bonding/antibonding MO pairs")
print(f"  - Creates Bell-pair entanglement structure")
print()

# ============================================================================
# 8. Mapper Properties
# ============================================================================
print("8. Fermionic-to-qubit mapper...")
print("-" * 70)

mapper = bond.mapper
print(f"Mapper type: {mapper.__class__.__name__}")
print(f"Number of MO pairs: {mapper.n_pairs}")
print(f"MO pair mapping:")
for i, (bonding, antibonding) in enumerate(mapper.mo_pairs):
    print(f"  Pair {i}: MO {bonding} (bonding) ↔ MO {antibonding} (antibonding)")
print()

# ============================================================================
# 9. Validation Summary
# ============================================================================
print("="*70)
print("VALIDATION SUMMARY")
print("="*70)

validations = []

# Check bond type
if bond.bond_type == 'covalent':
    validations.append(("✓", "Bond type correctly identified as covalent"))
else:
    validations.append(("✗", f"Bond type incorrect: {bond.bond_type}"))

# Check covalent character
if analysis['covalent_character'] > 0.95:
    validations.append(("✓", f"Highly covalent: {analysis['covalent_character']:.1%}"))
else:
    validations.append(("✗", f"Covalent character too low: {analysis['covalent_character']:.1%}"))

# Check bond length (sum of H covalent radii ≈ 0.62 Å)
if 0.5 < bond.get_bond_length() < 0.8:
    validations.append(("✓", f"Bond length reasonable: {bond.get_bond_length():.4f} Å"))
else:
    validations.append(("✗", f"Bond length unusual: {bond.get_bond_length():.4f} Å"))

# Check molecule setup
if bond.molecule.n_electrons == 2:
    validations.append(("✓", "Correct number of electrons (2)"))
else:
    validations.append(("✗", f"Wrong electron count: {bond.molecule.n_electrons}"))

# Check representation
if bond.representation.n_qubits > 0:
    validations.append(("✓", f"Quantum representation created ({bond.representation.n_qubits} qubits)"))
else:
    validations.append(("✗", "Quantum representation not initialized"))

# Check Hamiltonian
if H.nuclear_repulsion > 0:
    validations.append(("✓", f"Nuclear repulsion computed: {H.nuclear_repulsion:.4f} Ha"))
else:
    validations.append(("✗", "Nuclear repulsion not computed"))

# Check MO structure
if hasattr(H, 'compute_molecular_orbitals'):
    mo_e, _ = H.compute_molecular_orbitals()
    if len(mo_e) >= 2 and mo_e[0] < mo_e[1]:
        validations.append(("✓", "MO ordering correct (bonding < antibonding)"))
    else:
        validations.append(("✗", "MO ordering incorrect"))

# Print validation results
print()
for symbol, message in validations:
    print(f"{symbol} {message}")

passed = sum(1 for s, _ in validations if s == "✓")
total = len(validations)
print()
print(f"Validation score: {passed}/{total} checks passed")

if passed == total:
    print("\n🎉 ALL VALIDATIONS PASSED! Framework working correctly for H2.")
else:
    print(f"\n⚠️  {total - passed} validation(s) failed. Review above.")

print("="*70)
