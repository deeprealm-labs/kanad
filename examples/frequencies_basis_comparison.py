"""
Vibrational Frequencies - Basis Set Comparison

Test how different basis sets affect frequency accuracy.

Basis sets to test:
- STO-3G (minimal)
- 6-31G (split-valence)
- 6-311G (triple-zeta)
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
from kanad.core.atom import Atom
from kanad.core.molecule import Molecule
from kanad.analysis import FrequencyCalculator

print("=" * 80)
print("VIBRATIONAL FREQUENCIES - BASIS SET COMPARISON")
print("=" * 80)

# Experimental values for H2
FREQ_EXP_H2 = 4401.2  # cm⁻¹
ZPE_EXP_H2 = 6.29     # kcal/mol

# Basis sets to test
basis_sets = ['sto-3g', '6-31g', '6-311g']

print(f"\nTesting H2 vibrational frequency with different basis sets")
print(f"Experimental: {FREQ_EXP_H2:.1f} cm⁻¹")
print(f"Expected ZPE: {ZPE_EXP_H2:.2f} kcal/mol")
print("\n" + "=" * 80)

results = {}

for basis in basis_sets:
    print(f"\n{'='*80}")
    print(f"BASIS SET: {basis.upper()}")
    print(f"{'='*80}")

    try:
        # Create H2 molecule with current basis
        h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        h2 = Atom('H', position=np.array([0.741, 0.0, 0.0]))
        h2_molecule = Molecule(atoms=[h1, h2], charge=0, spin=0, basis=basis)

        print(f"Molecule: {h2_molecule.formula}")
        print(f"Basis: {basis}")
        print(f"Bond length: 0.741 Å")

        # Compute frequencies
        freq_calc = FrequencyCalculator(h2_molecule)
        result = freq_calc.compute_frequencies(method='HF', step_size=0.01, verbose=False)

        freq = result['frequencies'][0] if len(result['frequencies']) > 0 else 0.0
        zpe_kcal = result['zpe'] * 627.509
        force_const = result['force_constants'][0] if len(result['force_constants']) > 0 else 0.0

        # Calculate errors
        freq_error = abs(freq - FREQ_EXP_H2)
        freq_error_pct = (freq_error / FREQ_EXP_H2) * 100
        zpe_error = abs(zpe_kcal - ZPE_EXP_H2)

        results[basis] = {
            'frequency': freq,
            'zpe': zpe_kcal,
            'force_constant': force_const,
            'freq_error': freq_error,
            'freq_error_pct': freq_error_pct,
            'zpe_error': zpe_error
        }

        print(f"\nResults:")
        print(f"  Frequency:       {freq:10.2f} cm⁻¹")
        print(f"  Error:           {freq_error:10.2f} cm⁻¹  ({freq_error_pct:5.2f}%)")
        print(f"  ZPE:             {zpe_kcal:10.2f} kcal/mol")
        print(f"  ZPE error:       {zpe_error:10.2f} kcal/mol")
        print(f"  Force constant:  {force_const:10.4f} mdyn/Å")

        if freq_error_pct < 5:
            print(f"  ✓ Excellent agreement (<5% error)")
        elif freq_error_pct < 10:
            print(f"  ✓ Very good agreement (<10% error)")
        elif freq_error_pct < 20:
            print(f"  ✓ Good agreement (<20% error)")
        else:
            print(f"  ⚠ Moderate agreement (>20% error)")

    except Exception as e:
        print(f"\n✗ Failed with {basis}: {e}")
        results[basis] = None

# Summary comparison
print("\n" + "=" * 80)
print("SUMMARY COMPARISON")
print("=" * 80)

print(f"\n{'Basis Set':<12} {'Frequency':<12} {'Error':<10} {'% Error':<10} {'ZPE':<12} {'Force k':<10}")
print("-" * 80)

for basis in basis_sets:
    if results.get(basis):
        r = results[basis]
        print(f"{basis:<12} {r['frequency']:>10.1f}  {r['freq_error']:>8.1f}  {r['freq_error_pct']:>8.2f}%  {r['zpe']:>10.2f}  {r['force_constant']:>8.4f}")
    else:
        print(f"{basis:<12} {'FAILED':<50}")

print(f"\n{'Experimental':<12} {FREQ_EXP_H2:>10.1f}  {'---':>8}  {'---':>8}  {ZPE_EXP_H2:>10.2f}  {'5.75':>8}")

# Analysis
print("\n" + "=" * 80)
print("ANALYSIS")
print("=" * 80)

successful_results = {k: v for k, v in results.items() if v is not None}

if len(successful_results) > 1:
    best_basis = min(successful_results.keys(),
                     key=lambda k: successful_results[k]['freq_error_pct'])
    worst_basis = max(successful_results.keys(),
                      key=lambda k: successful_results[k]['freq_error_pct'])

    print(f"\nBest basis set: {best_basis.upper()}")
    print(f"  Error: {successful_results[best_basis]['freq_error_pct']:.2f}%")

    print(f"\nWorst basis set: {worst_basis.upper()}")
    print(f"  Error: {successful_results[worst_basis]['freq_error_pct']:.2f}%")

    # Improvement
    if worst_basis == 'sto-3g' and best_basis != 'sto-3g':
        improvement = successful_results[worst_basis]['freq_error_pct'] - successful_results[best_basis]['freq_error_pct']
        print(f"\nImprovement from STO-3G to {best_basis.upper()}: {improvement:.2f}% reduction in error")

print("\nObservations:")
print("  - Larger basis sets generally give better frequencies")
print("  - HF frequencies are typically 10-15% too high (systematic error)")
print("  - Scaling factors (e.g., 0.89 for HF/6-31G) can correct this")
print("  - For production, use 6-31G* or larger with MP2 or DFT")

print("\n" + "=" * 80)
print("RECOMMENDATION")
print("=" * 80)

if 'sto-3g' in successful_results:
    sto3g_error = successful_results['sto-3g']['freq_error_pct']
    print(f"\nSTO-3G error: {sto3g_error:.2f}%")

    better_bases = [b for b in successful_results
                    if successful_results[b]['freq_error_pct'] < sto3g_error]

    if better_bases:
        print(f"\nBetter alternatives: {', '.join(b.upper() for b in better_bases)}")
        print("\nFor production calculations, consider:")
        print("  - 6-31G* or 6-311G* for moderate accuracy")
        print("  - cc-pVDZ or cc-pVTZ for high accuracy")
        print("  - Use frequency scaling factor (0.89-0.95 for HF)")
    else:
        print("\nSTO-3G provides reasonable estimates for minimal cost")

print("=" * 80)
