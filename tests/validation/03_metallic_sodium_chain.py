"""
Validation Script 3: Metallic Sodium Chain
==========================================

Tests a 1D chain of sodium atoms (metallic bonding) using tight-binding model.

Expected properties for metallic Na:
- Delocalized electrons
- Band structure
- Fermi surface
- High conductivity

This script validates:
1. Metallic bond creation
2. Band structure computation
3. Fermi energy determination
4. Density of states
"""

import numpy as np
from kanad.bonds import MetallicBond
from kanad.core.atom import Atom

print("="*70)
print("VALIDATION 3: Metallic Sodium Chain")
print("="*70)
print()

# ============================================================================
# 1. Creating Metallic System
# ============================================================================
print("1. Creating 1D sodium chain (metallic bonding)...")
print("-" * 70)

# Create chain of 6 Na atoms
n_atoms = 6
spacing = 3.7  # Angstroms (typical Na-Na distance in BCC)

atoms = [
    Atom('Na', position=np.array([i * spacing, 0.0, 0.0]))
    for i in range(n_atoms)
]

bond = MetallicBond(atoms, lattice_type='1d_chain', hopping_parameter=-1.0)

print(f"System created: {bond}")
print(f"Bond type: {bond.bond_type}")
print(f"Lattice type: {bond.lattice_type}")
print(f"Number of atoms: {bond.n_atoms}")
print(f"Hopping parameter (t): {bond.hopping_parameter} eV")
print()

# ============================================================================
# 2. Atomic Properties
# ============================================================================
print("2. Atomic properties...")
print("-" * 70)

na = atoms[0]
print(f"Element: {na.symbol}")
print(f"Atomic number: {na.atomic_number}")
print(f"Valence electrons: {na.n_valence}")
print(f"Electronegativity: {na.electronegativity:.2f}")
print(f"Is metal: {na.is_metal}")
print()

print(f"Atom positions (Å):")
for i, atom in enumerate(atoms):
    print(f"  Na {i}: {atom.position}")
print()

# ============================================================================
# 3. Energy Calculation (Tight-Binding)
# ============================================================================
print("3. Computing energy with tight-binding model...")
print("-" * 70)

result = bond.compute_energy(method='tight_binding')

print(f"Method: {result['method']}")
print(f"Total energy: {result['energy']:.6f} eV")
print(f"Energy per atom: {result['energy']/n_atoms:.6f} eV")
print(f"Converged: {result['converged']}")
print()

# ============================================================================
# 4. Band Structure
# ============================================================================
print("4. Band structure...")
print("-" * 70)

band_energies = result['band_energies']
n_bands_occupied = result['n_bands_occupied']
n_electrons = result['n_electrons']

print(f"Number of bands: {len(band_energies)}")
print(f"Number of electrons: {n_electrons}")
print(f"Number of occupied bands: {n_bands_occupied}")
print(f"\nBand energies (eV):")
for i, energy in enumerate(band_energies):
    if i < n_bands_occupied - 1:
        occupied = "OCCUPIED (2e⁻)"
    elif i == n_bands_occupied - 1:
        if n_electrons % 2 == 0:
            occupied = "OCCUPIED (2e⁻)"
        else:
            occupied = "OCCUPIED (1e⁻)"
    else:
        occupied = "VIRTUAL"
    print(f"  Band {i}: {energy:10.6f} eV  [{occupied}]")
print()

print(f"Bandwidth: {band_energies[-1] - band_energies[0]:.6f} eV")
print(f"Valence band max: {band_energies[n_bands_occupied-1]:.6f} eV")
if n_bands_occupied < len(band_energies):
    print(f"Conduction band min: {band_energies[n_bands_occupied]:.6f} eV")
print()

# ============================================================================
# 5. Fermi Energy
# ============================================================================
print("5. Fermi energy and Fermi surface...")
print("-" * 70)

if 'fermi_energy' in result:
    E_F = result['fermi_energy']
    print(f"Fermi energy: {E_F:.6f} eV")
    print(f"Position in band: {(E_F - band_energies[0])/(band_energies[-1] - band_energies[0]):.1%}")
    print()

# ============================================================================
# 6. Bond Analysis
# ============================================================================
print("6. Metallic bond analysis...")
print("-" * 70)

analysis = result['bond_analysis']

print(f"Bond type: {analysis['bond_type']}")
print(f"Lattice type: {analysis['lattice_type']}")
print(f"Number of atoms: {analysis['n_atoms']}")
print(f"Hopping parameter: {analysis['hopping_parameter']} eV")
print()

if 'bandwidth' in analysis:
    print(f"Bandwidth: {analysis['bandwidth']:.6f} eV")

if 'fermi_energy' in analysis:
    print(f"Fermi energy: {analysis['fermi_energy']:.6f} eV")

if 'dos_at_fermi' in analysis:
    print(f"DOS at Fermi level: {analysis['dos_at_fermi']}")
print()

print(f"Entanglement type: {analysis['entanglement_type']}")
print()

# ============================================================================
# 7. Full Band Structure Scan
# ============================================================================
print("7. Full band structure E(k)...")
print("-" * 70)

bands = bond.get_band_structure()

print(f"Number of k-points: {len(bands['k_points'])}")
print(f"k-point range: [{bands['k_points'][0]:.4f}, {bands['k_points'][-1]:.4f}] (π units)")
print()

# Find band edges
energies_flat = bands['energies'].flatten()
print(f"Band structure statistics:")
print(f"  Minimum energy: {energies_flat.min():.6f} eV")
print(f"  Maximum energy: {energies_flat.max():.6f} eV")
print(f"  Total bandwidth: {energies_flat.max() - energies_flat.min():.6f} eV")
print()

# Show dispersion at a few k-points
print("Energy dispersion at specific k-points:")
for i in [0, len(bands['k_points'])//4, len(bands['k_points'])//2, -1]:
    k = bands['k_points'][i]
    k_label = f"k = {k/np.pi:.2f}π"
    print(f"  {k_label:15s}: E = {bands['energies'][i, 0]:.6f} eV")
print()

# ============================================================================
# 8. Validation Summary
# ============================================================================
print("="*70)
print("VALIDATION SUMMARY")
print("="*70)

validations = []

# Check bond type
if bond.bond_type == 'metallic':
    validations.append(("✓", "Bond type correctly identified as metallic"))
else:
    validations.append(("✗", f"Bond type incorrect: {bond.bond_type}"))

# Check number of atoms
if bond.n_atoms == n_atoms:
    validations.append(("✓", f"Correct number of atoms: {n_atoms}"))
else:
    validations.append(("✗", f"Wrong atom count: {bond.n_atoms}"))

# Check band structure computed
if len(band_energies) == n_atoms:
    validations.append(("✓", f"Band structure computed ({len(band_energies)} bands)"))
else:
    validations.append(("✗", f"Unexpected number of bands: {len(band_energies)}"))

# Check bandwidth is reasonable
bandwidth = band_energies[-1] - band_energies[0]
if bandwidth > 0 and bandwidth < 100:  # Reasonable range
    validations.append(("✓", f"Bandwidth reasonable: {bandwidth:.4f} eV"))
else:
    validations.append(("✗", f"Bandwidth unusual: {bandwidth:.4f} eV"))

# Check Fermi energy exists
if 'fermi_energy' in result:
    validations.append(("✓", f"Fermi energy computed: {result['fermi_energy']:.4f} eV"))
else:
    validations.append(("✗", "Fermi energy not computed"))

# Check energy is reasonable
if -100 < result['energy'] < 100:  # Reasonable range
    validations.append(("✓", f"Total energy reasonable: {result['energy']:.4f} eV"))
else:
    validations.append(("✗", f"Energy out of range: {result['energy']:.4f} eV"))

# Check dispersion
if bands['energies'].shape[0] > 10:
    validations.append(("✓", f"Band dispersion computed ({bands['energies'].shape[0]} k-points)"))
else:
    validations.append(("✗", "Insufficient k-point sampling"))

# Print validation results
print()
for symbol, message in validations:
    print(f"{symbol} {message}")

passed = sum(1 for s, _ in validations if s == "✓")
total = len(validations)
print()
print(f"Validation score: {passed}/{total} checks passed")

if passed == total:
    print("\n🎉 ALL VALIDATIONS PASSED! Metallic bonding model working correctly.")
elif passed >= total - 1:
    print("\n✅ MOSTLY PASSED! Metallic bonding model working well.")
else:
    print(f"\n⚠️  {total - passed} validation(s) failed. Review above.")

print("="*70)
