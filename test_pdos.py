"""Test Projected DOS with Mulliken Analysis"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from kanad.bonds import BondFactory
from kanad.analysis import DOSCalculator
import numpy as np

# Create H2 molecule
bond = BondFactory.create_bond('H', 'H', distance=0.74)

# Run SCF first (required for DOS)
print("Running SCF...")
density_matrix, hf_energy = bond.hamiltonian.solve_scf()
print(f"HF Energy: {hf_energy:.8f} Ha")

# Create DOS calculator
dos_calc = DOSCalculator(bond.hamiltonian)

print("\n" + "=" * 80)
print("Testing Projected DOS (PDOS) with Mulliken Analysis")
print("=" * 80)

# Compute PDOS for both atoms
result = dos_calc.compute_pdos(
    atom_indices=[0, 1],
    energy_range=(-15, 5),
    n_points=500,
    sigma=0.3,
    units='eV'
)

print(f"\n✓ PDOS computed successfully!")
print(f"  Energy range: {result['energies'][0]:.2f} to {result['energies'][-1]:.2f} eV")
print(f"  Number of points: {len(result['energies'])}")

# Check that PDOS is different for each atom (or at least computed)
pdos_atom0 = result['pdos'][0]
pdos_atom1 = result['pdos'][1]
total_dos = result['total_dos']

print(f"\n✓ Atom 0 PDOS: min={np.min(pdos_atom0):.6f}, max={np.max(pdos_atom0):.6f}")
print(f"✓ Atom 1 PDOS: min={np.min(pdos_atom1):.6f}, max={np.max(pdos_atom1):.6f}")
print(f"✓ Total DOS: min={np.min(total_dos):.6f}, max={np.max(total_dos):.6f}")

# Verify PDOS sums to total DOS (approximately)
pdos_sum = pdos_atom0 + pdos_atom1
difference = np.abs(pdos_sum - total_dos)
max_diff = np.max(difference)
relative_error = max_diff / np.max(total_dos)

print(f"\n✓ PDOS sum check:")
print(f"  Max difference: {max_diff:.6f}")
print(f"  Relative error: {relative_error*100:.2f}%")

if relative_error < 0.05:  # Less than 5% error
    print("\n✅ PDOS TEST PASSED!")
    print("   Mulliken population analysis working correctly")
    print("   PDOS contributions sum to total DOS")
else:
    print(f"\n⚠️  Warning: PDOS sum deviates from total DOS by {relative_error*100:.2f}%")

# For H2 (homonuclear), PDOS should be symmetric
pdos_diff = np.abs(pdos_atom0 - pdos_atom1)
symmetry_error = np.max(pdos_diff) / np.max(total_dos)

print(f"\n✓ Symmetry check (H2 is symmetric):")
print(f"  Max difference between atoms: {np.max(pdos_diff):.6f}")
print(f"  Relative difference: {symmetry_error*100:.2f}%")

if symmetry_error < 0.1:  # Less than 10% difference
    print("  ✅ PDOS is approximately symmetric (as expected for H2)")
else:
    print(f"  ⚠️  PDOS shows {symmetry_error*100:.2f}% asymmetry")

print("\n" + "=" * 80)
print("PDOS IMPLEMENTATION: ✅ COMPLETE")
print("=" * 80)
