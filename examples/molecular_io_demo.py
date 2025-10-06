"""
Demonstration of Kanad Molecular I/O Module

Shows how to:
1. Load molecules from SMILES strings
2. Load molecules from XYZ files
3. Save molecules to XYZ files
4. Run calculations on loaded molecules
"""

print("=" * 80)
print("KANAD MOLECULAR I/O DEMONSTRATION")
print("=" * 80)

# Example 1: Load from SMILES
print("\n1. Loading molecules from SMILES strings")
print("-" * 80)

from kanad.io import from_smiles

# Simple molecules
water = from_smiles("O", name="Water")
print(f"✓ Water: {water.formula}, {water.n_atoms} atoms, {water.n_electrons} electrons")

ethanol = from_smiles("CCO", name="Ethanol")
print(f"✓ Ethanol: {ethanol.formula}, {ethanol.n_atoms} atoms")

methane = from_smiles("C", name="Methane")
print(f"✓ Methane: {methane.formula}, {methane.n_atoms} atoms")

# Aromatic molecule
benzene = from_smiles("c1ccccc1", name="Benzene")
print(f"✓ Benzene: {benzene.formula}, {benzene.n_atoms} atoms")

# Charged species
ammonium = from_smiles("[NH4+]", name="Ammonium")
print(f"✓ Ammonium: {ammonium.formula}, charge={ammonium.charge}, spin={ammonium.spin}")

chloride = from_smiles("[Cl-]", name="Chloride")
print(f"✓ Chloride: {chloride.formula}, charge={chloride.charge}")

# Example 2: Save to XYZ file
print("\n2. Saving molecules to XYZ files")
print("-" * 80)

from kanad.io import to_xyz
import tempfile
import os

temp_dir = tempfile.mkdtemp()

# Save water molecule
water_file = os.path.join(temp_dir, "water.xyz")
to_xyz(water, water_file, comment="Water molecule from SMILES")
print(f"✓ Saved water to: {water_file}")

# Save benzene
benzene_file = os.path.join(temp_dir, "benzene.xyz")
to_xyz(benzene, benzene_file, comment="Benzene (aromatic)")
print(f"✓ Saved benzene to: {benzene_file}")

# Example 3: Load from XYZ file
print("\n3. Loading molecules from XYZ files")
print("-" * 80)

from kanad.io import from_xyz

# Load water back
loaded_water = from_xyz(water_file)
print(f"✓ Loaded: {loaded_water.formula}, {loaded_water.n_atoms} atoms")

# Load benzene back
loaded_benzene = from_xyz(benzene_file)
print(f"✓ Loaded: {loaded_benzene.formula}, {loaded_benzene.n_atoms} atoms")

# Example 4: Run calculations on loaded molecules
print("\n4. Running calculations on loaded molecules")
print("-" * 80)

# Hartree-Fock calculation on water
print(f"\nCalculating HF energy for {water.formula}...")
result = water.compute_energy(method='HF')
print(f"✓ HF Energy: {result['energy_ha']:.6f} Ha ({result['energy']:.4f} eV)")
print(f"  Converged: {result['converged']}")

# Hartree-Fock on methane
print(f"\nCalculating HF energy for {methane.formula}...")
result_ch4 = methane.compute_energy(method='HF')
print(f"✓ HF Energy: {result_ch4['energy_ha']:.6f} Ha ({result_ch4['energy']:.4f} eV)")

# Example 5: Round-trip test
print("\n5. Round-trip test (SMILES → Molecule → XYZ → Molecule)")
print("-" * 80)

# Use ethanol for quick demo
ethanol_saved = os.path.join(temp_dir, "ethanol.xyz")
to_xyz(ethanol, ethanol_saved, comment="Ethanol from SMILES round-trip")
print(f"✓ Saved ethanol to: {ethanol_saved}")

# Load back
ethanol_reloaded = from_xyz(ethanol_saved)
print(f"✓ Loaded from XYZ: {ethanol_reloaded.formula}")
print(f"  Atoms match: {ethanol_reloaded.n_atoms == ethanol.n_atoms}")

# Example 6: Molecule library
print("\n6. Building a molecule library")
print("-" * 80)

molecules = {
    "Water": from_smiles("O"),
    "Methane": from_smiles("C"),
    "Ethanol": from_smiles("CCO"),
    "Ammonia": from_smiles("N"),
    "Benzene": from_smiles("c1ccccc1"),
}

print(f"✓ Created library of {len(molecules)} molecules:")
for name, mol in molecules.items():
    print(f"  • {name:15s}: {mol.formula:15s} ({mol.n_atoms:2d} atoms, {mol.n_electrons:3d} electrons)")

# Cleanup
print("\n" + "=" * 80)
print("DEMONSTRATION COMPLETE")
print("=" * 80)
print(f"\nTemporary files saved to: {temp_dir}")
print("Clean up with: rm -rf " + temp_dir)
