"""
Validation Script 2: NaCl Ionic Bond
====================================

Tests the NaCl molecule (classic ionic bond) and validates framework capabilities.

Expected properties for NaCl:
- Large electronegativity difference (2.23)
- Significant charge transfer (~0.8 e)
- Long bond length (~2.36 Å)
- Predominantly ionic character (>70%)

This script validates:
1. Ionic bond auto-detection
2. Charge transfer analysis
3. Electrostatic interactions
4. Localized orbital representation
"""

import numpy as np
from kanad.bonds import BondFactory

print("="*70)
print("VALIDATION 2: NaCl Ionic Bond")
print("="*70)
print()

# ============================================================================
# 1. Bond Creation and Auto-Detection
# ============================================================================
print("1. Creating NaCl bond with automatic bond type detection...")
print("-" * 70)

bond = BondFactory.create_bond('Na', 'Cl')

print(f"Bond created: {bond}")
print(f"Bond type: {bond.bond_type}")
print(f"Bond length: {bond.get_bond_length():.4f} Å")
print(f"Expected: ~2.36 Å (experimental)")
print()

# ============================================================================
# 2. Quick Bond Information
# ============================================================================
print("2. Quick bond information analysis...")
print("-" * 70)
info = BondFactory.quick_bond_info('Na', 'Cl')
print(f"Predicted type: {info['predicted_type']}")
print(f"ΔEN: {info['electronegativity_difference']:.4f}")
print(f"  Na electronegativity: {info['electronegativity_1']:.2f}")
print(f"  Cl electronegativity: {info['electronegativity_2']:.2f}")
print(f"Expected: ΔEN > 1.7 for ionic bonding")
print(f"Estimated bond length: {info['estimated_bond_length']:.4f} Å")
print(f"Rationale: {info['rationale']}")
print()

# ============================================================================
# 3. Ionic Bond-Specific Properties
# ============================================================================
print("3. Ionic bond specific properties...")
print("-" * 70)

print(f"Donor atom: {bond.donor.symbol}")
print(f"  Atomic number: {bond.donor.atomic_number}")
print(f"  Valence electrons: {bond.donor.n_valence}")
print(f"  Electronegativity: {bond.donor.electronegativity:.2f}")
print(f"  Is metal: {bond.donor.is_metal}")
print()

print(f"Acceptor atom: {bond.acceptor.symbol}")
print(f"  Atomic number: {bond.acceptor.atomic_number}")
print(f"  Valence electrons: {bond.acceptor.n_valence}")
print(f"  Electronegativity: {bond.acceptor.electronegativity:.2f}")
print(f"  Is metal: {bond.acceptor.is_metal}")
print()

# ============================================================================
# 4. Molecule and Representation
# ============================================================================
print("4. Molecular system properties...")
print("-" * 70)

print(f"Molecule: {bond.molecule}")
print(f"Number of atoms: {bond.molecule.n_atoms}")
print(f"Number of electrons: {bond.molecule.n_electrons}")
print(f"  Na: 11 electrons")
print(f"  Cl: 17 electrons")
print(f"  Total: {bond.molecule.n_electrons} electrons")
print()

print(f"Representation: {bond.representation}")
print(f"Representation type: {bond.representation.__class__.__name__}")
print(f"Number of qubits: {bond.representation.n_qubits}")
print(f"Expected: Second quantization (localized orbitals)")
print()

# ============================================================================
# 5. Hamiltonian Analysis
# ============================================================================
print("5. Ionic Hamiltonian properties...")
print("-" * 70)

H = bond.hamiltonian
print(f"Hamiltonian type: {H.__class__.__name__}")
print(f"Number of orbitals: {H.n_orbitals}")
print(f"Number of electrons: {H.n_electrons}")
print()

# Check for ionic-specific properties
if hasattr(H, 'on_site_energies'):
    print("On-site energies (atomic orbital energies):")
    for i, energy in enumerate(H.on_site_energies):
        atom = bond.atoms[i]
        print(f"  {atom.symbol}: {energy:.4f} Hartree")
    print()

if hasattr(H, 'transfer_integrals'):
    print("Transfer integrals (hopping parameters):")
    print(H.transfer_integrals)
    print()

if hasattr(H, 'U_parameters'):
    print("Hubbard U parameters (on-site repulsion):")
    for i, U in enumerate(H.U_parameters):
        atom = bond.atoms[i]
        print(f"  U_{atom.symbol}: {U:.4f} Hartree")
    print()

# ============================================================================
# 6. Charge Transfer Analysis
# ============================================================================
print("6. Charge transfer analysis...")
print("-" * 70)

charges = bond.get_charge_distribution()

print("Charge distribution:")
for symbol, charge in charges.items():
    sign = "+" if charge > 0 else ""
    print(f"  {symbol}: {sign}{charge:.4f} e")
print()

print("Interpretation:")
print(f"  Na loses ~{charges[bond.donor.symbol]:.2f} electrons → Na⁺")
print(f"  Cl gains ~{abs(charges[bond.acceptor.symbol]):.2f} electrons → Cl⁻")
print(f"  Expected: ~0.8 e transfer for NaCl")
print()

# ============================================================================
# 7. Bond Analysis
# ============================================================================
print("7. Detailed bond analysis...")
print("-" * 70)

analysis = bond.analyze()

print(f"Bond type: {analysis['bond_type']}")
print(f"Donor: {analysis['donor']}")
print(f"Acceptor: {analysis['acceptor']}")
print(f"Bond length: {analysis['bond_length']:.4f} Å")
print()

print("Character analysis:")
print(f"  ΔEN: {analysis['electronegativity_difference']:.4f}")
print(f"  Charge transfer: {analysis['charge_transfer']:.4f} e")
print(f"  Ionic character: {analysis['ionic_character']:.2%}")
print(f"  Expected: >70% ionic for NaCl")
print()

if 'coulomb_energy' in analysis:
    print(f"Electrostatic (Coulomb) energy: {analysis['coulomb_energy']:.4f} Hartree")
    print(f"                                {analysis['coulomb_energy']*27.211:.2f} eV")
    print(f"                                {analysis['coulomb_energy']*627.5:.2f} kcal/mol")
    print()

print(f"Entanglement type: {analysis['entanglement_type']}")
print(f"Governance protocol: {analysis['governance_protocol']}")
print()

# ============================================================================
# 8. Governance Protocol Validation
# ============================================================================
print("8. Governance protocol properties...")
print("-" * 70)

protocol = bond.governance
print(f"Protocol type: {protocol.__class__.__name__}")
print(f"Protocol enforces:")
print(f"  - Localized atomic orbitals")
print(f"  - Minimal entanglement (charge transfer only)")
print(f"  - Electrostatic interactions dominant")
print()

# ============================================================================
# 9. Mapper Properties
# ============================================================================
print("9. Fermionic-to-qubit mapper...")
print("-" * 70)

mapper = bond.mapper
print(f"Mapper type: {mapper.__class__.__name__}")
print(f"Expected: Jordan-Wigner (for localized states)")
print(f"Property: Sequential mapping, local operators")
print()

# ============================================================================
# 10. Validation Summary
# ============================================================================
print("="*70)
print("VALIDATION SUMMARY")
print("="*70)

validations = []

# Check bond type
if bond.bond_type == 'ionic':
    validations.append(("✓", "Bond type correctly identified as ionic"))
else:
    validations.append(("✗", f"Bond type incorrect: {bond.bond_type}"))

# Check ΔEN
if info['electronegativity_difference'] > 1.7:
    validations.append(("✓", f"ΔEN > 1.7: {info['electronegativity_difference']:.4f}"))
else:
    validations.append(("✗", f"ΔEN too small: {info['electronegativity_difference']:.4f}"))

# Check donor/acceptor assignment
if bond.donor.symbol == 'Na' and bond.acceptor.symbol == 'Cl':
    validations.append(("✓", "Donor/acceptor correctly assigned (Na→Cl)"))
else:
    validations.append(("✗", "Donor/acceptor assignment incorrect"))

# Check charge transfer
if analysis['charge_transfer'] > 0.6:
    validations.append(("✓", f"Significant charge transfer: {analysis['charge_transfer']:.2f} e"))
else:
    validations.append(("✗", f"Insufficient charge transfer: {analysis['charge_transfer']:.2f} e"))

# Check ionic character
if analysis['ionic_character'] > 0.6:
    validations.append(("✓", f"Predominantly ionic: {analysis['ionic_character']:.1%}"))
else:
    validations.append(("✗", f"Ionic character too low: {analysis['ionic_character']:.1%}"))

# Check representation type
if 'SecondQuantization' in bond.representation.__class__.__name__:
    validations.append(("✓", "Second quantization representation (localized)"))
else:
    validations.append(("~", f"Representation: {bond.representation.__class__.__name__}"))

# Check mapper type
if 'JordanWigner' in mapper.__class__.__name__:
    validations.append(("✓", "Jordan-Wigner mapper (for ionic bonding)"))
else:
    validations.append(("~", f"Mapper: {mapper.__class__.__name__}"))

# Check governance
if 'Ionic' in protocol.__class__.__name__:
    validations.append(("✓", "Ionic governance protocol applied"))
else:
    validations.append(("✗", f"Wrong protocol: {protocol.__class__.__name__}"))

# Print validation results
print()
for symbol, message in validations:
    print(f"{symbol} {message}")

passed = sum(1 for s, _ in validations if s == "✓")
total = len(validations)
print()
print(f"Validation score: {passed}/{total} checks passed")

if passed == total:
    print("\n🎉 ALL VALIDATIONS PASSED! Framework working correctly for NaCl.")
elif passed >= total - 2:
    print("\n✅ MOSTLY PASSED! Framework working well for ionic bonds.")
else:
    print(f"\n⚠️  {total - passed} validation(s) failed. Review above.")

print("="*70)
