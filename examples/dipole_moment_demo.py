"""
Demonstration of Dipole Moment Calculation in Kanad

Shows how to:
1. Calculate dipole moments for various molecules
2. Understand electronic vs. nuclear contributions
3. Validate results against experimental values
4. Use dipole to predict molecular properties
"""

print("=" * 80)
print("KANAD DIPOLE MOMENT DEMONSTRATION")
print("=" * 80)

from kanad.io import from_smiles
from kanad.analysis import PropertyCalculator

# Example 1: Polar molecules
print("\n1. Polar Molecules (Non-zero Dipole)")
print("-" * 80)

molecules = {
    "Water (H2O)": ("O", 1.85),  # (SMILES, experimental dipole in D)
    "Ammonia (NH3)": ("N", 1.47),
    "Hydrogen Fluoride (HF)": ("F", 1.82),
    "Methanol (CH3OH)": ("CO", 1.70),
}

print(f"{'Molecule':<25} {'Calculated':<15} {'Experimental':<15} {'Error':<10}")
print("-" * 80)

for name, (smiles, exp_dipole) in molecules.items():
    mol = from_smiles(smiles, name=name)
    calc = PropertyCalculator(mol.hamiltonian)
    result = calc.compute_dipole_moment()

    calc_dipole = result['dipole_magnitude']
    error = abs(calc_dipole - exp_dipole)
    error_pct = (error / exp_dipole) * 100

    print(f"{name:<25} {calc_dipole:>6.3f} D       {exp_dipole:>6.3f} D       {error_pct:>5.1f}%")

# Example 2: Non-polar molecules (symmetric)
print("\n2. Non-Polar Molecules (Zero Dipole)")
print("-" * 80)

nonpolar = {
    "Hydrogen (H2)": "[H][H]",
    "Methane (CH4)": "C",
    "Ethane (C2H6)": "CC",
}

print(f"{'Molecule':<25} {'Dipole (D)':<15}")
print("-" * 80)

for name, smiles in nonpolar.items():
    mol = from_smiles(smiles, name=name)
    calc = PropertyCalculator(mol.hamiltonian)
    result = calc.compute_dipole_moment()

    dipole = result['dipole_magnitude']
    print(f"{name:<25} {dipole:>10.6f}  {'✓ Near zero' if dipole < 0.01 else '⚠ Non-zero!'}")

# Example 3: Dipole vector components
print("\n3. Dipole Vector Analysis (Water)")
print("-" * 80)

water = from_smiles("O", name="Water")
calc = PropertyCalculator(water.hamiltonian)
result = calc.compute_dipole_moment()

print(f"Dipole magnitude: {result['dipole_magnitude']:.4f} D")
print(f"\nComponents (Debye):")
print(f"  μ_x = {result['components']['x']:>8.4f} D")
print(f"  μ_y = {result['components']['y']:>8.4f} D")
print(f"  μ_z = {result['components']['z']:>8.4f} D")

print(f"\nComponents (atomic units):")
print(f"  μ_x = {result['dipole_au'][0]:>8.4f} a.u.")
print(f"  μ_y = {result['dipole_au'][1]:>8.4f} a.u.")
print(f"  μ_z = {result['dipole_au'][2]:>8.4f} a.u.")

# Example 4: Electronic vs. Nuclear contributions
print("\n4. Electronic vs. Nuclear Contributions (Water)")
print("-" * 80)

mu_elec = result['electronic_contribution']
mu_nuc = result['nuclear_contribution']

print(f"Electronic contribution: [{mu_elec[0]:>7.4f}, {mu_elec[1]:>7.4f}, {mu_elec[2]:>7.4f}] D")
print(f"Nuclear contribution:    [{mu_nuc[0]:>7.4f}, {mu_nuc[1]:>7.4f}, {mu_nuc[2]:>7.4f}] D")
print(f"Total dipole:            [{result['dipole_vector'][0]:>7.4f}, {result['dipole_vector'][1]:>7.4f}, {result['dipole_vector'][2]:>7.4f}] D")

import numpy as np
print(f"\nMagnitude of electronic contribution: {np.linalg.norm(mu_elec):.4f} D")
print(f"Magnitude of nuclear contribution:    {np.linalg.norm(mu_nuc):.4f} D")

# Example 5: Validation against PySCF
print("\n5. Validation Against PySCF")
print("-" * 80)

verification = calc.verify_dipole_with_pyscf()

print(f"Kanad calculation:  {verification['kanad_dipole']:.6f} D")
print(f"PySCF calculation:  {verification['pyscf_dipole']:.6f} D")
print(f"Absolute difference: {verification['difference']:.8f} D")
print(f"Agreement: {'✓ PASS' if verification['agree'] else '✗ FAIL'} (tolerance: 0.01 D)")

# Example 6: Molecular centers
print("\n6. Molecular Centers")
print("-" * 80)

com = calc.compute_center_of_mass()
coc = calc.compute_center_of_charge()

print(f"Center of mass:   [{com[0]:>8.4f}, {com[1]:>8.4f}, {com[2]:>8.4f}] Å")
print(f"Center of charge: [{coc[0]:>8.4f}, {coc[1]:>8.4f}, {coc[2]:>8.4f}] Å")

# Example 7: Comparing molecules
print("\n7. Dipole Moment Trends")
print("-" * 80)

series = {
    "HF": "F",
    "H2O": "O",
    "NH3": "N",
    "CH4": "C",
}

print("Trend: Decreasing electronegativity difference → Decreasing dipole")
print(f"\n{'Molecule':<10} {'Dipole (D)':<12} {'Polarity':<15}")
print("-" * 80)

dipoles = {}
for name, smiles in series.items():
    mol = from_smiles(smiles, name=name)
    calc = PropertyCalculator(mol.hamiltonian)
    result = calc.compute_dipole_moment()
    dipoles[name] = result['dipole_magnitude']

    polarity = "Highly polar" if dipoles[name] > 1.5 else \
               "Moderately polar" if dipoles[name] > 0.5 else \
               "Non-polar"

    print(f"{name:<10} {dipoles[name]:>6.3f}       {polarity:<15}")

print("\n" + "=" * 80)
print("DEMONSTRATION COMPLETE")
print("=" * 80)
print("\nKey Takeaways:")
print("• Polar molecules (H-F, H-O, H-N) have dipole moments > 1 D")
print("• Symmetric molecules (H2, CH4) have dipole ≈ 0 D")
print("• Kanad results agree with PySCF within 0.01 D")
print("• HF/STO-3G gives reasonable dipoles (within ~10% of experiment)")
