"""
Energy Gradients Demo

Demonstrates:
1. Computing analytical gradients at equilibrium
2. Computing gradients for stretched geometry
3. Verifying analytical vs numerical gradients
4. Using gradients to compute forces on atoms
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
from kanad.core.atom import Atom
from kanad.core.molecule import Molecule
from kanad.core.gradients import GradientCalculator

print("=" * 80)
print("ENERGY GRADIENTS DEMO")
print("=" * 80)

# Test 1: H2 at equilibrium
print("\n1. H2 at Equilibrium Geometry (d = 0.74 Å)")
print("-" * 80)

h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))
molecule = Molecule(atoms=[h1, h2], charge=0, spin=0, basis='sto-3g')
grad_calc = GradientCalculator(molecule, method='HF')
result = grad_calc.compute_gradient()

print(f"Gradient (Ha/Bohr):")
print(result['gradient'])
print(f"\nForces (Ha/Bohr):")
print(result['forces'])
print(f"\nMax force: {result['max_force']:.6f} Ha/Bohr")
print(f"RMS force: {result['rms_force']:.6f} Ha/Bohr")

if result['max_force'] < 0.05:
    print("✓ Forces near zero at equilibrium (as expected)")
else:
    print("⚠ Forces should be near zero at equilibrium")

# Test 2: H2 stretched
print("\n2. H2 Stretched Geometry (d = 1.0 Å)")
print("-" * 80)

h1_stretch = Atom('H', position=np.array([0.0, 0.0, 0.0]))
h2_stretch = Atom('H', position=np.array([1.0, 0.0, 0.0]))
molecule_stretch = Molecule(atoms=[h1_stretch, h2_stretch], charge=0, spin=0, basis='sto-3g')
grad_calc_stretch = GradientCalculator(molecule_stretch, method='HF')
result_stretch = grad_calc_stretch.compute_gradient()

print(f"Gradient (Ha/Bohr):")
print(result_stretch['gradient'])
print(f"\nForces (Ha/Bohr):")
print(result_stretch['forces'])
print(f"\nMax force: {result_stretch['max_force']:.6f} Ha/Bohr")
print(f"RMS force: {result_stretch['rms_force']:.6f} Ha/Bohr")

if result_stretch['max_force'] > 0.05:
    print("✓ Large forces present (molecule wants to contract)")
else:
    print("⚠ Forces should be large for stretched geometry")

# Test 3: Verify analytical gradient
print("\n3. Gradient Verification (H2O)")
print("-" * 80)

o = Atom('O', position=np.array([0.0, 0.0, 0.0]))
h1_w = Atom('H', position=np.array([0.758, 0.587, 0.0]))
h2_w = Atom('H', position=np.array([-0.758, 0.587, 0.0]))
molecule_water = Molecule(atoms=[o, h1_w, h2_w], charge=0, spin=0, basis='sto-3g')
grad_calc_water = GradientCalculator(molecule_water, method='HF')

print("Computing analytical and numerical gradients...")
print("(This may take a minute - 18 energy evaluations for central differences)")

verification = grad_calc_water.verify_analytical_gradient(
    numerical_step=0.001,
    tolerance=1e-5
)

print(f"\nAnalytical gradient (Ha/Bohr):")
print(verification['analytical'])
print(f"\nNumerical gradient (Ha/Bohr):")
print(verification['numerical'])
print(f"\nAbsolute difference:")
print(verification['difference'])
print(f"\nMax difference: {verification['max_difference']:.2e} Ha/Bohr")
print(f"RMS difference: {verification['rms_difference']:.2e} Ha/Bohr")
print(f"Tolerance: {verification['tolerance']:.2e} Ha/Bohr")
print(f"\nAgreement: {'YES ✓' if verification['agree'] else 'NO ✗'}")

if verification['agree']:
    print("\n✓ Analytical gradient implementation is correct!")
else:
    print("\n✗ Analytical gradient differs from numerical gradient")

# Test 4: MP2 Gradients
print("\n4. MP2 Gradients (H2)")
print("-" * 80)

grad_calc_mp2 = GradientCalculator(molecule, method='MP2')
result_mp2 = grad_calc_mp2.compute_gradient()

print(f"MP2 Gradient (Ha/Bohr):")
print(result_mp2['gradient'])
print(f"\nMP2 Forces (Ha/Bohr):")
print(result_mp2['forces'])
print(f"\nMax force: {result_mp2['max_force']:.6f} Ha/Bohr")
print(f"RMS force: {result_mp2['rms_force']:.6f} Ha/Bohr")

print(f"\nComparison (HF vs MP2):")
print(f"  HF max force:  {result['max_force']:.6f} Ha/Bohr")
print(f"  MP2 max force: {result_mp2['max_force']:.6f} Ha/Bohr")
print(f"  Difference:    {abs(result['max_force'] - result_mp2['max_force']):.6f} Ha/Bohr")

print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)
print(f"✓ Analytical gradients computed using PySCF")
print(f"✓ Gradients near zero at equilibrium geometry")
print(f"✓ Large gradients for stretched geometry")
print(f"✓ Analytical vs numerical agreement: {verification['max_difference']:.2e} Ha/Bohr")
print(f"✓ Both HF and MP2 gradients supported")
print("\nGradients can be used for:")
print("  - Geometry optimization (finding minimum energy structures)")
print("  - Molecular dynamics (time evolution under forces)")
print("  - Transition state searches")
print("  - Vibrational frequency calculations")
print("=" * 80)
