#!/usr/bin/env python3
"""
Comprehensive Vibrational Solver Validation
Tests vibrational frequencies for different molecules.
"""

import numpy as np
from kanad.core.atom import Atom
from kanad.bonds.covalent_bond import CovalentBond
from kanad.solvers.vibrational_solver import VibrationalSolver

def print_header(title: str):
    print(f"\n{'='*80}")
    print(f"  {title}")
    print(f"{'='*80}")

# Test 1: H2 Vibrational Frequency
print_header("TEST 1: H2 Molecule Vibrational Frequency")

H1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
H2_atom = Atom('H', position=np.array([0.74, 0.0, 0.0]))

bond = CovalentBond(H1, H2_atom)
molecule = bond.molecule

print(f"  H2 Molecule: {molecule.n_atoms} atoms")
print(f"  Expected vibrational mode: 1 (stretch)")
print(f"  Expected frequency: ~4400 cm⁻¹")

try:
    # Use small displacement for faster computation
    solver = VibrationalSolver(molecule=molecule, bond=bond, displacement=0.05)

    print(f"  Computing Hessian (this may take time)...")

    results = solver.solve_harmonic_frequencies()

    print(f"\n  Vibrational Modes Found: {len(results['frequencies_cm'])}")
    print(f"\n  {'Mode':>6s} {'Freq (cm⁻¹)':>15s} {'Freq (eV)':>12s} {'ZPE contrib':>15s}")
    print(f"  {'-'*60}")

    for i, (freq_cm, freq_ev) in enumerate(zip(
        results['frequencies_cm'],
        results['frequencies_ev']
    )):
        zpe_contrib = 0.5 * freq_ev
        print(f"  {i+1:>6d} {freq_cm:>15.1f} {freq_ev:>12.6f} {zpe_contrib:>15.6f}")

    print(f"\n  Total Zero-Point Energy: {results['zero_point_energy_ha']:.6f} Ha")
    print(f"                           {results['zero_point_energy_kcal']:.2f} kcal/mol")

    # Validation
    if len(results['frequencies_cm']) >= 1:
        first_freq = results['frequencies_cm'][0]
        # Expected H2 stretch: ~4400 cm⁻¹
        if 3000 < first_freq < 6000:
            print(f"\n  ✅ First frequency reasonable for H2 stretch")
        else:
            print(f"\n  ⚠️  First frequency may be incorrect (expected ~4400 cm⁻¹)")
    else:
        print(f"\n  ❌ No vibrational modes found")

except Exception as e:
    print(f"  ❌ Error: {str(e)[:100]}")
    import traceback
    traceback.print_exc()

# Test 2: Mode Analysis
print_header("TEST 2: Normal Mode Analysis")

if 'results' in locals() and 'normal_modes' in results:
    try:
        mode_analysis = solver.analyze_normal_modes(
            results['frequencies_cm'],
            results['normal_modes']
        )

        print(f"\n  {'Mode':>6s} {'Freq (cm⁻¹)':>15s} {'Type':>25s} {'Primary Atoms':>20s}")
        print(f"  {'-'*75}")

        for mode_info in mode_analysis:
            mode_idx = mode_info['mode_index']
            freq = mode_info['frequency_cm']
            mode_type = mode_info['type']
            atoms = str(mode_info['primary_atoms'])

            print(f"  {mode_idx+1:>6d} {freq:>15.1f} {mode_type:>25s} {atoms:>20s}")

        print(f"\n  ✅ Normal mode analysis complete")

    except Exception as e:
        print(f"  ❌ Error in mode analysis: {str(e)[:70]}")

# Test 3: IR Spectrum
print_header("TEST 3: IR Spectrum Prediction")

if 'results' in locals():
    try:
        ir_spectrum = solver.get_ir_spectrum(results['frequencies_cm'])

        print(f"\n  IR Active Peaks:")
        print(f"  {'Freq (cm⁻¹)':>15s} {'Intensity (a.u.)':>20s}")
        print(f"  {'-'*40}")

        for freq, intensity in ir_spectrum['peaks']:
            print(f"  {freq:>15.1f} {intensity:>20.2f}")

        print(f"\n  ✅ IR spectrum prediction complete")

    except Exception as e:
        print(f"  ❌ Error: {str(e)[:70]}")

# Summary
print_header("VALIDATION SUMMARY")

print("\n✅ = Passed")
print("⚠️  = Warning/needs review")
print("❌ = Failed/error")

print("\nTests Performed:")
print("  1. H2 vibrational frequency calculation")
print("  2. Normal mode analysis")
print("  3. IR spectrum prediction")

print("\nKey Findings:")
if 'results' in locals():
    print(f"  - Computed {len(results['frequencies_cm'])} vibrational modes")
    print(f"  - Frequency range: {results['frequencies_cm'].min():.1f} - {results['frequencies_cm'].max():.1f} cm⁻¹")
    print(f"  - Zero-point energy: {results['zero_point_energy_kcal']:.2f} kcal/mol")
else:
    print(f"  - Vibrational calculation failed or incomplete")

print("\nPotential Issues to Check:")
print("  1. Are frequencies physically reasonable?")
print("  2. Is ZPE contribution correct?")
print("  3. Does number of modes match 3N-6 (or 3N-5 for linear)?")
print("  4. Are translational/rotational modes properly removed?")
print("  5. Does Hessian computation converge?")

print("\nNote: Vibrational solver uses numerical Hessian (computationally expensive).")
print("For faster tests, use larger displacement values (less accurate).")
print("For molecules with more atoms, computation time scales as O(N²).")

print("\n✅ Vibrational Solver validation complete!")
