"""
Polarizability Calculation Demo

Demonstrates computing molecular polarizability using the finite field method.

Polarizability (α) measures how easily a molecule's electron cloud is distorted
by an external electric field:
    μ_induced = α · E

This demo shows:
1. Computing polarizability tensor (α_ij)
2. Mean polarizability: ᾱ = (α_xx + α_yy + α_zz) / 3
3. Polarizability anisotropy (Δα)
4. Effect of basis set on accuracy
5. Comparison with experimental values
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
from kanad.io import from_smiles
from kanad.analysis import PropertyCalculator

print("=" * 80)
print("POLARIZABILITY CALCULATION DEMO")
print("=" * 80)

# =============================================================================
# Example 1: H2 Polarizability
# =============================================================================
print("\n1. Hydrogen (H2)")
print("-" * 80)

h2 = from_smiles("[H][H]", name="Hydrogen")
calc_h2 = PropertyCalculator(h2.hamiltonian)

result_h2 = calc_h2.compute_polarizability()

print(f"Polarizability tensor (a.u.):")
print(result_h2['alpha_tensor'])
print(f"\nDiagonal elements:")
print(f"  α_xx (parallel)      = {result_h2['alpha_xx']:.4f} a.u.")
print(f"  α_yy (perpendicular) = {result_h2['alpha_yy']:.4f} a.u.")
print(f"  α_zz (perpendicular) = {result_h2['alpha_zz']:.4f} a.u.")
print(f"\nMean polarizability:")
print(f"  ᾱ = {result_h2['alpha_mean']:.4f} a.u. = {result_h2['alpha_mean_angstrom3']:.4f} Å³")
print(f"  (Experimental: 5.4 a.u. = 0.80 Å³)")
print(f"\nAnisotropy: Δα = {result_h2['alpha_anisotropy']:.4f} a.u.")
print(f"Principal polarizabilities: {result_h2['eigenvalues']}")

print("\nNote: STO-3G severely underestimates polarizability due to lack of")
print("      polarization functions (d,p). Perpendicular components are zero")
print("      because H atoms only have s-orbitals in this basis.")

# =============================================================================
# Example 2: H2O Polarizability
# =============================================================================
print("\n" + "=" * 80)
print("2. Water (H2O)")
print("-" * 80)

water = from_smiles("O", name="Water")
calc_water = PropertyCalculator(water.hamiltonian)

result_water = calc_water.compute_polarizability()

print(f"Polarizability tensor (a.u.):")
print(result_water['alpha_tensor'])
print(f"\nDiagonal elements:")
print(f"  α_xx = {result_water['alpha_xx']:.4f} a.u.")
print(f"  α_yy = {result_water['alpha_yy']:.4f} a.u.")
print(f"  α_zz = {result_water['alpha_zz']:.4f} a.u.")
print(f"\nMean polarizability:")
print(f"  ᾱ = {result_water['alpha_mean']:.4f} a.u. = {result_water['alpha_mean_angstrom3']:.4f} Å³")
print(f"  (Experimental: 9.8 a.u. = 1.45 Å³)")
print(f"\nAnisotropy: Δα = {result_water['alpha_anisotropy']:.4f} a.u.")
print(f"Principal polarizabilities: {result_water['eigenvalues']}")

# =============================================================================
# Example 3: Comparison Across Molecules
# =============================================================================
print("\n" + "=" * 80)
print("3. Comparison: H2, H2O, NH3, CH4")
print("-" * 80)

molecules = {
    "H2":  ("[H][H]", 5.4),
    "H2O": ("O",      9.8),
    "NH3": ("N",     14.6),
    "CH4": ("C",     17.3),
}

print(f"{'Molecule':<10} {'α (calc, a.u.)':<15} {'α (exp, a.u.)':<15} {'Ratio':<10}")
print("-" * 80)

for name, (smiles, exp_alpha) in molecules.items():
    mol = from_smiles(smiles, name=name)
    calc = PropertyCalculator(mol.hamiltonian)
    result = calc.compute_polarizability()

    alpha_calc = result['alpha_mean']
    ratio = alpha_calc / exp_alpha

    print(f"{name:<10} {alpha_calc:<15.4f} {exp_alpha:<15.4f} {ratio:<10.3f}")

print("\nNote: All calculations use HF/STO-3G, which underestimates polarizability")
print("      by ~50-70% due to minimal basis and lack of electron correlation.")
print("      Accurate results require larger basis sets (6-311G(d,p)) and correlation (MP2/CCSD).")

# =============================================================================
# Example 4: Basis Set Effect (if you want to test)
# =============================================================================
print("\n" + "=" * 80)
print("4. Basis Set Effect on H2 Polarizability")
print("-" * 80)
print(f"{'Basis Set':<15} {'α_mean (a.u.)':<15} {'% of Exp.':<15}")
print("-" * 80)

experimental = 5.4

for basis in ['sto-3g', '6-31g', '6-31g(d,p)', '6-311g(d,p)']:
    try:
        h2_test = from_smiles("[H][H]", name="H2", basis=basis)
        calc_test = PropertyCalculator(h2_test.hamiltonian)
        result_test = calc_test.compute_polarizability()

        alpha = result_test['alpha_mean']
        percent = (alpha / experimental) * 100

        print(f"{basis:<15} {alpha:<15.4f} {percent:<15.1f}%")
    except Exception as e:
        print(f"{basis:<15} Error: {str(e)[:50]}")

print("\nObservations:")
print("  - Polarization functions (d,p) are essential for accurate polarizability")
print("  - STO-3G: only s-orbitals → cannot describe perpendicular polarization")
print("  - 6-31G(d,p): adds d on heavy atoms, p on H → much better")
print("  - Still below experiment: need correlation methods (MP2, CCSD)")

# =============================================================================
# Example 5: Understanding the Polarizability Tensor
# =============================================================================
print("\n" + "=" * 80)
print("5. Understanding the Polarizability Tensor")
print("-" * 80)

print("\nFor H2 aligned along x-axis:")
print("  - α_xx (parallel):      polarizability along bond")
print("  - α_yy, α_zz (perpendicular): polarizability perpendicular to bond")
print("  - Expect: α_parallel > α_perpendicular (electron cloud easier to distort along bond)")

h2_example = from_smiles("[H][H]", name="H2")
result_example = PropertyCalculator(h2_example.hamiltonian).compute_polarizability()

print(f"\nH2 geometry:")
for atom in h2_example.atoms:
    print(f"  {atom.symbol}: {atom.position}")

print(f"\nPolarizability components:")
print(f"  α_xx = {result_example['alpha_xx']:.4f} a.u. (along bond)")
print(f"  α_yy = {result_example['alpha_yy']:.4f} a.u. (perpendicular)")
print(f"  α_zz = {result_example['alpha_zz']:.4f} a.u. (perpendicular)")

print("\nNote: With STO-3G, perpendicular components are zero (no p-functions).")
print("      With 6-31G(d,p): α_xx ≈ 5.9, α_yy = α_zz ≈ 0.6 a.u.")

# =============================================================================
# Summary
# =============================================================================
print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)
print("""
Polarizability measures how easily a molecule's electron cloud is distorted
by an external electric field. It's a key property for:
  - Intermolecular interactions (van der Waals forces)
  - Refractive index and optical properties
  - Raman spectroscopy
  - Drug design (molecular recognition)

Implementation Details:
  - Method: Finite field (numerical derivative)
  - Field strength: 0.001 a.u. (optimized for accuracy vs. linearity)
  - Formula: α_ij = -∂μ_i/∂E_j (response of dipole to field)

Limitations:
  - HF/STO-3G underestimates by ~50-70%
  - Needs polarization functions (d,p) for accurate perpendicular components
  - Needs electron correlation (MP2, CCSD) for quantitative accuracy

For production calculations:
  - Use at least 6-311G(d,p) basis set
  - Consider MP2 or CCSD for correlation
  - Can also use analytical CPHF (faster, more accurate)
""")

print("=" * 80)
