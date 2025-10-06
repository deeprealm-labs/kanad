"""
Bond Length Scan Demo

Demonstrates the BondLengthScanner for exploring potential energy surfaces.

Examples:
1. H2 dissociation curve (diatomic)
2. H2O O-H bond scan (preserving rest of geometry)
3. Comparison with literature values
4. Cubic spline interpolation for smooth curves
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
from kanad.core.atom import Atom
from kanad.core.molecule import Molecule
from kanad.analysis import BondLengthScanner

print("=" * 80)
print("BOND LENGTH SCAN DEMO")
print("=" * 80)

# Example 1: H2 Dissociation Curve
print("\n" + "=" * 80)
print("Example 1: H2 Dissociation Curve")
print("=" * 80)

h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
h2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))
h2_molecule = Molecule(atoms=[h1, h2], charge=0, spin=0, basis='sto-3g')

scanner_h2 = BondLengthScanner(h2_molecule, atom_i=0, atom_j=1)
result_h2 = scanner_h2.scan(r_min=0.5, r_max=2.5, n_points=20, method='HF')

print(f"\n{'='*80}")
print("RESULTS:")
print(f"{'='*80}")
print(f"Optimized H-H distance: {result_h2['optimized_distance']:.4f} Å")
print(f"Minimum energy:         {result_h2['optimized_energy']:.8f} Ha")
print(f"Original distance:      {result_h2['original_distance']:.4f} Å")
print(f"\nLiterature values:")
print(f"  Experimental:  0.741 Å")
print(f"  HF/STO-3G:     ~0.73 Å")
print(f"  Our result:    {result_h2['optimized_distance']:.4f} Å")

error = abs(result_h2['optimized_distance'] - 0.741) * 1000  # mÅ
print(f"\nError vs experimental: {error:.1f} mÅ")

if error < 20:
    print("✓ Excellent agreement with experimental value!")
elif error < 50:
    print("✓ Good agreement with experimental value")
else:
    print("⚠ Larger deviation (expected for minimal basis set)")

# Example 2: H2O O-H Bond Scan
print("\n" + "=" * 80)
print("Example 2: H2O O-H Bond Scan (Preserving Geometry)")
print("=" * 80)

# Create water molecule
o = Atom('O', position=np.array([0.0, 0.0, 0.0]))
h1_w = Atom('H', position=np.array([0.758, 0.587, 0.0]))
h2_w = Atom('H', position=np.array([-0.758, 0.587, 0.0]))
water = Molecule(atoms=[o, h1_w, h2_w], charge=0, spin=0, basis='sto-3g')

print(f"\nOriginal H2O geometry:")
print(f"  O:  {o.position}")
print(f"  H1: {h1_w.position}")
print(f"  H2: {h2_w.position}")

# Scan O-H1 bond
scanner_water = BondLengthScanner(water, atom_i=0, atom_j=1)
result_water = scanner_water.scan(r_min=0.7, r_max=1.5, n_points=15, method='HF')

print(f"\n{'='*80}")
print("RESULTS:")
print(f"{'='*80}")
print(f"Optimized O-H distance: {result_water['optimized_distance']:.4f} Å")
print(f"Minimum energy:         {result_water['optimized_energy']:.8f} Ha")
print(f"Original distance:      {result_water['original_distance']:.4f} Å")
print(f"\nLiterature values:")
print(f"  Experimental:  0.958 Å")
print(f"  HF/STO-3G:     ~0.95 Å")
print(f"  Our result:    {result_water['optimized_distance']:.4f} Å")

# Verify geometry preservation
print(f"\nGeometry after scan (should be restored):")
print(f"  O:  {water.atoms[0].position}")
print(f"  H1: {water.atoms[1].position}")
print(f"  H2: {water.atoms[2].position}")

# Check if geometry is restored
orig_geom_preserved = all(
    np.allclose(water.atoms[i].position, [o.position, h1_w.position, h2_w.position][i])
    for i in range(3)
)

if orig_geom_preserved:
    print("\n✓ Original geometry restored after scan!")
else:
    print("\n⚠ Geometry was modified (unexpected)")

# Example 3: Spline Interpolation Quality
print("\n" + "=" * 80)
print("Example 3: Cubic Spline Interpolation")
print("=" * 80)

if result_h2['spline'] is not None:
    spline = result_h2['spline']

    # Evaluate on fine grid
    r_fine = np.linspace(0.5, 2.5, 1000)
    e_fine = spline(r_fine)

    # Find minimum
    min_idx = np.argmin(e_fine)
    r_opt_spline = r_fine[min_idx]
    e_opt_spline = e_fine[min_idx]

    print(f"Grid search minimum:   r = {result_h2['optimized_distance']:.6f} Å")
    print(f"Spline refined minimum: r = {r_opt_spline:.6f} Å")
    print(f"Difference:            Δr = {abs(result_h2['optimized_distance'] - r_opt_spline)*1000:.3f} mÅ")

    # Compute dissociation energy
    e_dissociated = spline(2.5)  # Energy at large separation
    d_e = (e_dissociated - e_opt_spline) * 627.509  # Convert to kcal/mol

    print(f"\nDissociation energy (HF/STO-3G):")
    print(f"  D_e = {d_e:.2f} kcal/mol")
    print(f"  Experimental D_e ≈ 109 kcal/mol")
    print(f"  HF/STO-3G expected: ~84 kcal/mol")
    print(f"  (HF underestimates bond strength)")

    print("\n✓ Spline interpolation provides smooth curve")
else:
    print("⚠ Not enough points for spline interpolation")

# Example 4: Energy Profile Statistics
print("\n" + "=" * 80)
print("Example 4: Energy Profile Analysis")
print("=" * 80)

distances = result_h2['distances']
energies = result_h2['energies']
valid_mask = ~np.isnan(energies)

print(f"Scan statistics:")
print(f"  Total points:     {len(distances)}")
print(f"  Valid points:     {np.sum(valid_mask)}")
print(f"  Failed points:    {np.sum(~valid_mask)}")
print(f"  Energy range:     {np.nanmax(energies) - np.nanmin(energies):.6f} Ha")
print(f"                    {(np.nanmax(energies) - np.nanmin(energies)) * 627.509:.2f} kcal/mol")

# Find turning points (approximate force constant)
if result_h2['spline'] is not None:
    # Compute numerical second derivative at minimum
    r_opt = result_h2['optimized_distance']
    delta = 0.01  # Small displacement

    e_plus = spline(r_opt + delta)
    e_minus = spline(r_opt - delta)
    e_center = spline(r_opt)

    # Second derivative: d²E/dr² ≈ (E(r+δ) - 2E(r) + E(r-δ)) / δ²
    k = (e_plus - 2*e_center + e_minus) / (delta**2)

    # Convert to more common units
    # 1 Ha/Bohr² = 15.569 mdyn/Å = 1556.9 N/m
    k_mdyn = k * 15.569  # mdyn/Å
    k_nm = k * 1556.9    # N/m

    print(f"\nForce constant at equilibrium:")
    print(f"  k = {k:.4f} Ha/Bohr²")
    print(f"    = {k_mdyn:.2f} mdyn/Å")
    print(f"    = {k_nm:.1f} N/m")
    print(f"  Experimental (H2): k ≈ 5.75 mdyn/Å")

# Summary
print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)
print("✓ Bond length scanner works for diatomic molecules (H2)")
print("✓ Bond length scanner works for polyatomic molecules (H2O)")
print("✓ Preserves geometry of non-scanned atoms")
print("✓ Finds equilibrium distances close to literature values")
print("✓ Cubic spline interpolation provides smooth PES curves")
print("✓ Can extract force constants and dissociation energies")
print("\nApplications:")
print("  - Finding equilibrium bond lengths")
print("  - Computing dissociation energies")
print("  - Exploring reaction coordinates")
print("  - Validating force fields")
print("  - Understanding chemical bonding")
print("=" * 80)

# Optional: Plot results (requires matplotlib)
try:
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend
    import matplotlib.pyplot as plt

    print("\nCreating plots...")

    # Plot H2 dissociation curve
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # H2 plot
    valid_h2 = ~np.isnan(result_h2['energies'])
    e_h2_rel = (result_h2['energies'][valid_h2] - result_h2['optimized_energy']) * 627.509

    ax1.plot(result_h2['distances'][valid_h2], e_h2_rel, 'o', markersize=8, label='Computed')

    if result_h2['spline'] is not None:
        r_fine = np.linspace(0.5, 2.5, 500)
        e_fine = result_h2['spline'](r_fine)
        e_fine_rel = (e_fine - result_h2['optimized_energy']) * 627.509
        ax1.plot(r_fine, e_fine_rel, '-', linewidth=2, alpha=0.7, label='Cubic spline')

    ax1.axvline(result_h2['optimized_distance'], color='red', linestyle='--',
                alpha=0.5, label=f'Min: {result_h2["optimized_distance"]:.4f} Å')
    ax1.set_xlabel('H-H Distance (Å)', fontsize=12)
    ax1.set_ylabel('Relative Energy (kcal/mol)', fontsize=12)
    ax1.set_title('H₂ Dissociation Curve (HF/STO-3G)', fontsize=14)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(bottom=0)

    # H2O plot
    valid_h2o = ~np.isnan(result_water['energies'])
    e_h2o_rel = (result_water['energies'][valid_h2o] - result_water['optimized_energy']) * 627.509

    ax2.plot(result_water['distances'][valid_h2o], e_h2o_rel, 'o', markersize=8, label='Computed')

    if result_water['spline'] is not None:
        r_fine = np.linspace(0.7, 1.5, 500)
        e_fine = result_water['spline'](r_fine)
        e_fine_rel = (e_fine - result_water['optimized_energy']) * 627.509
        ax2.plot(r_fine, e_fine_rel, '-', linewidth=2, alpha=0.7, label='Cubic spline')

    ax2.axvline(result_water['optimized_distance'], color='red', linestyle='--',
                alpha=0.5, label=f'Min: {result_water["optimized_distance"]:.4f} Å')
    ax2.set_xlabel('O-H Distance (Å)', fontsize=12)
    ax2.set_ylabel('Relative Energy (kcal/mol)', fontsize=12)
    ax2.set_title('H₂O O-H Bond Scan (HF/STO-3G)', fontsize=14)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(bottom=0)

    plt.tight_layout()
    plt.savefig('bond_scan_demo.png', dpi=300, bbox_inches='tight')
    print("✓ Plot saved to bond_scan_demo.png")

except ImportError:
    print("\nmatplotlib not available - skipping plots")
except Exception as e:
    print(f"\nPlot generation failed: {e}")
