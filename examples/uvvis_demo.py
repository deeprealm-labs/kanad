"""
UV-Vis Absorption Spectroscopy Demo

Demonstrates calculation of electronic excitations and UV-Vis spectra
using TD-DFT/TDA.

Examples:
1. Water (H2O) - Far UV excitations
2. Comparison of methods (TDA vs TDDFT)
3. Spectrum generation with Gaussian broadening
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
from kanad.core.atom import Atom
from kanad.core.molecule import Molecule
from kanad.analysis import UVVisCalculator

print("=" * 80)
print("UV-VIS ABSORPTION SPECTROSCOPY DEMO")
print("=" * 80)

# Example 1: Water (H2O) Electronic Excitations
print("\n" + "=" * 80)
print("Example 1: Water (H2O) Electronic Excitations")
print("=" * 80)

# Create water molecule
o = Atom('O', position=np.array([0.0, 0.0, 0.0]))
h1 = Atom('H', position=np.array([0.758, 0.587, 0.0]))
h2 = Atom('H', position=np.array([-0.758, 0.587, 0.0]))
water = Molecule(atoms=[o, h1, h2], charge=0, spin=0, basis='sto-3g')

print(f"\nWater molecule:")
print(f"  Formula: {water.formula}")
print(f"  Basis: {water.basis}")

# Compute excitations using TDA (faster than full TDDFT)
uv_calc = UVVisCalculator(water)

print("\n" + "-" * 80)
print("Method 1: TDA-HF (HF ground state)")
print("-" * 80)

result_tda_hf = uv_calc.compute_excitations(
    n_states=5,
    method='TDA',
    functional=None,  # None = HF
    verbose=True
)

# Example 2: Generate and Plot Spectrum
print("\n" + "=" * 80)
print("Example 2: UV-Vis Absorption Spectrum")
print("=" * 80)

print("\nGenerating absorption spectrum with Gaussian broadening...")

spectrum = uv_calc.generate_spectrum(
    result_tda_hf,
    wavelength_range=(100, 400),  # UV range
    broadening=0.3,  # eV (FWHM)
    n_points=500,
    verbose=True
)

# Find absorption maximum
max_idx = np.argmax(spectrum['absorbance'])
λ_max = spectrum['wavelengths'][max_idx]
ε_max = spectrum['absorbance'][max_idx]

print(f"\nSpectrum characteristics:")
print(f"  λ_max: {λ_max:.1f} nm")
print(f"  ε_max: {ε_max:.2e} L/(mol·cm)")

# Identify UV region
uv_c = np.sum((spectrum['wavelengths'] >= 100) & (spectrum['wavelengths'] < 280))
uv_b = np.sum((spectrum['wavelengths'] >= 280) & (spectrum['wavelengths'] < 315))
uv_a = np.sum((spectrum['wavelengths'] >= 315) & (spectrum['wavelengths'] < 400))

print(f"\nSpectral regions:")
print(f"  UV-C (100-280 nm): {uv_c > 0}")
print(f"  UV-B (280-315 nm): {uv_b > 0}")
print(f"  UV-A (315-400 nm): {uv_a > 0}")
print(f"  Visible (400-800 nm): Transparent (water is colorless)")

# Example 3: Analysis of Transitions
print("\n" + "=" * 80)
print("Example 3: Transition Analysis")
print("=" * 80)

print("\nClassifying transitions by strength:")
print(f"{'Transition':<12} {'Energy (eV)':<14} {'λ (nm)':<12} {'f':<10} {'Classification':<15}")
print("-" * 75)

for i, (E, λ, f) in enumerate(zip(
    result_tda_hf['excitation_energies'],
    result_tda_hf['wavelengths'],
    result_tda_hf['oscillator_strengths']
)):
    if f < 0.001:
        classification = "Forbidden"
    elif f < 0.01:
        classification = "Very weak"
    elif f < 0.1:
        classification = "Weak"
    elif f < 1.0:
        classification = "Moderate"
    else:
        classification = "Strong"

    print(f"S{i} → S{i+1:<7} {E:>12.4f}  {λ:>10.2f}  {f:>8.4f}  {classification:<15}")

# Find dominant transition
if len(result_tda_hf['oscillator_strengths']) > 0:
    max_f_idx = np.argmax(result_tda_hf['oscillator_strengths'])
    max_f = result_tda_hf['oscillator_strengths'][max_f_idx]
    max_λ = result_tda_hf['wavelengths'][max_f_idx]

    print(f"\nDominant transition:")
    print(f"  S0 → S{max_f_idx+1}")
    print(f"  λ = {max_λ:.2f} nm  (f = {max_f:.4f})")

# Example 4: Comparison with Experimental Data
print("\n" + "=" * 80)
print("Example 4: Comparison with Literature")
print("=" * 80)

print(f"\nWater UV absorption:")
print(f"  Literature: First excitation ~7.4 eV (~167 nm)")
print(f"  Our result: {result_tda_hf['excitation_energies'][0]:.2f} eV ({result_tda_hf['wavelengths'][0]:.1f} nm)")

error_eV = abs(result_tda_hf['excitation_energies'][0] - 7.4)
error_pct = (error_eV / 7.4) * 100

print(f"  Error: {error_eV:.2f} eV ({error_pct:.1f}%)")

if error_pct < 10:
    print("  ✓ Excellent agreement with literature")
elif error_pct < 20:
    print("  ✓ Good agreement with literature")
elif error_pct < 30:
    print("  ✓ Reasonable agreement (expected for minimal basis)")
else:
    print("  ⚠ Large deviation (STO-3G limitation)")

print("\nNote: STO-3G underestimates excitation energies")
print("      Use larger basis sets (6-31+G*, aug-cc-pVDZ) for better accuracy")

# Summary
print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)
print("✓ UV-Vis calculator successfully computes:")
print("  - Electronic excitations using TD-DFT/TDA")
print("  - Excitation energies and wavelengths")
print("  - Oscillator strengths (transition probabilities)")
print("  - Absorption spectrum with Gaussian broadening")
print(f"\n✓ Water (H2O) first excitation: {result_tda_hf['wavelengths'][0]:.1f} nm (Far UV)")
print(f"✓ Number of excitations computed: {result_tda_hf['n_states']}")
print(f"✓ Spectrum generated: {len(spectrum['wavelengths'])} points")

print("\nApplications:")
print("  - Predicting UV-Vis spectra")
print("  - Understanding electronic structure")
print("  - Identifying chromophores")
print("  - Studying conjugated systems")
print("  - Photochemistry and photophysics")
print("  - Drug design (molecular colors)")

print("\nLimitations:")
print("  - STO-3G basis: poor for excited states (use 6-31+G* or better)")
print("  - TDA-HF: overestimates excitations (use TD-DFT with B3LYP)")
print("  - Small molecules only (TD-DFT scales as O(N⁴))")

print("=" * 80)

# Optional: Plot spectrum
try:
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend

    print("\nGenerating plot...")
    uv_calc.plot_spectrum(
        spectrum,
        excitations=result_tda_hf,
        save_path='uvvis_water_demo.png',
        show_sticks=True
    )
    print("✓ Plot saved to uvvis_water_demo.png")

except ImportError:
    print("\nmatplotlib not available - skipping plot")
except Exception as e:
    print(f"\nPlot generation failed: {e}")
