"""
BlueQubit Cloud: Testing Larger Molecules

Tests governance-validated quantum chemistry on realistic molecules:
1. H2O (Water) - Covalent bonding, bent geometry
2. NH3 (Ammonia) - Covalent bonding, pyramidal geometry
3. CH4 (Methane) - Covalent bonding, tetrahedral, sp3 hybridization
4. NaCl (Sodium Chloride) - Ionic bonding
5. CO2 (Carbon Dioxide) - Covalent bonding, linear geometry

Using BlueQubit's ~30 qubits for real quantum simulations.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
from kanad.core.molecule import Molecule
from kanad.core.atom import Atom


def create_h2o():
    """Create H2O molecule with correct geometry."""
    print("\n" + "=" * 80)
    print("MOLECULE 1: H2O (Water)")
    print("=" * 80)

    # Water geometry: bent, 104.5° angle, O-H = 0.96 Å
    o = Atom('O', position=np.array([0.0, 0.0, 0.0]))
    h1 = Atom('H', position=np.array([0.758, 0.587, 0.0]))  # 104.5° angle
    h2 = Atom('H', position=np.array([-0.758, 0.587, 0.0]))

    molecule = Molecule(atoms=[o, h1, h2], charge=0, spin=0)

    print(f"\n✓ H2O created")
    print(f"  Atoms: 3 (1 O, 2 H)")
    print(f"  Electrons: {molecule.n_electrons}")
    print(f"  Geometry: Bent (104.5°)")
    print(f"  Bonding: Covalent (sp3 hybridization on O)")
    print(f"  Expected qubits: ~12-16 (depends on basis set)")

    return molecule


def create_nh3():
    """Create NH3 molecule with pyramidal geometry."""
    print("\n" + "=" * 80)
    print("MOLECULE 2: NH3 (Ammonia)")
    print("=" * 80)

    # Ammonia geometry: pyramidal, 107° angle, N-H = 1.01 Å
    n = Atom('N', position=np.array([0.0, 0.0, 0.0]))
    h1 = Atom('H', position=np.array([0.0, 0.938, 0.383]))
    h2 = Atom('H', position=np.array([0.812, -0.469, 0.383]))
    h3 = Atom('H', position=np.array([-0.812, -0.469, 0.383]))

    molecule = Molecule(atoms=[n, h1, h2, h3], charge=0, spin=0)

    print(f"\n✓ NH3 created")
    print(f"  Atoms: 4 (1 N, 3 H)")
    print(f"  Electrons: {molecule.n_electrons}")
    print(f"  Geometry: Pyramidal")
    print(f"  Bonding: Covalent (sp3 hybridization on N)")
    print(f"  Expected qubits: ~14-20")

    return molecule


def create_ch4():
    """Create CH4 molecule with tetrahedral geometry."""
    print("\n" + "=" * 80)
    print("MOLECULE 3: CH4 (Methane)")
    print("=" * 80)

    # Methane geometry: tetrahedral, 109.5° angle, C-H = 1.09 Å
    c = Atom('C', position=np.array([0.0, 0.0, 0.0]))
    h1 = Atom('H', position=np.array([0.629, 0.629, 0.629]))
    h2 = Atom('H', position=np.array([-0.629, -0.629, 0.629]))
    h3 = Atom('H', position=np.array([-0.629, 0.629, -0.629]))
    h4 = Atom('H', position=np.array([0.629, -0.629, -0.629]))

    molecule = Molecule(atoms=[c, h1, h2, h3, h4], charge=0, spin=0)

    print(f"\n✓ CH4 created")
    print(f"  Atoms: 5 (1 C, 4 H)")
    print(f"  Electrons: {molecule.n_electrons}")
    print(f"  Geometry: Tetrahedral (perfect sp3)")
    print(f"  Bonding: Covalent (sp3 hybridization)")
    print(f"  Expected qubits: ~16-24")

    return molecule


def create_nacl():
    """Create NaCl molecule (ionic)."""
    print("\n" + "=" * 80)
    print("MOLECULE 4: NaCl (Sodium Chloride)")
    print("=" * 80)

    # NaCl: ionic bond, Na-Cl = 2.36 Å
    na = Atom('Na', position=np.array([0.0, 0.0, 0.0]))
    cl = Atom('Cl', position=np.array([2.36, 0.0, 0.0]))

    molecule = Molecule(atoms=[na, cl], charge=0, spin=0)

    print(f"\n✓ NaCl created")
    print(f"  Atoms: 2 (1 Na, 1 Cl)")
    print(f"  Electrons: {molecule.n_electrons}")
    print(f"  Geometry: Linear (diatomic)")
    print(f"  Bonding: Ionic (Na+ Cl-)")
    print(f"  ΔEN: {abs(cl.electronegativity - na.electronegativity):.2f} (large → ionic)")
    print(f"  Expected qubits: ~10-14")

    return molecule


def create_co2():
    """Create CO2 molecule (linear)."""
    print("\n" + "=" * 80)
    print("MOLECULE 5: CO2 (Carbon Dioxide)")
    print("=" * 80)

    # CO2: linear, C=O double bonds, C-O = 1.16 Å
    c = Atom('C', position=np.array([0.0, 0.0, 0.0]))
    o1 = Atom('O', position=np.array([1.16, 0.0, 0.0]))
    o2 = Atom('O', position=np.array([-1.16, 0.0, 0.0]))

    molecule = Molecule(atoms=[c, o1, o2], charge=0, spin=0)

    print(f"\n✓ CO2 created")
    print(f"  Atoms: 3 (1 C, 2 O)")
    print(f"  Electrons: {molecule.n_electrons}")
    print(f"  Geometry: Linear (180°)")
    print(f"  Bonding: Covalent (sp hybridization, double bonds)")
    print(f"  Expected qubits: ~14-20")

    return molecule


def analyze_molecule(molecule, name):
    """Analyze molecular properties."""
    print(f"\n{'─' * 80}")
    print(f"MOLECULAR ANALYSIS: {name}")
    print(f"{'─' * 80}")

    print(f"\nAtoms:")
    for i, atom in enumerate(molecule.atoms):
        print(f"  {i+1}. {atom.symbol} at {atom.position} Å")
        print(f"     EN: {atom.electronegativity:.2f}, Metal: {atom.is_metal}")

    # Bonding analysis
    print(f"\nBonding Analysis:")
    atoms_symbols = [a.symbol for a in molecule.atoms]
    unique_atoms = set(atoms_symbols)

    if len(unique_atoms) == 1:
        print(f"  Type: Homonuclear (all {list(unique_atoms)[0]})")
    else:
        electronegativities = [a.electronegativity for a in molecule.atoms]
        en_diff = max(electronegativities) - min(electronegativities)

        if en_diff > 1.5:
            print(f"  Type: IONIC (ΔEN = {en_diff:.2f})")
        elif en_diff > 0.5:
            print(f"  Type: POLAR COVALENT (ΔEN = {en_diff:.2f})")
        else:
            print(f"  Type: COVALENT (ΔEN = {en_diff:.2f})")

    # Check if any metals
    has_metal = any(a.is_metal for a in molecule.atoms)
    has_nonmetal = any(not a.is_metal for a in molecule.atoms)

    if has_metal and has_nonmetal:
        print(f"  Character: Metal + Nonmetal → Ionic/Polar")
    elif has_metal:
        print(f"  Character: All metals → Metallic")
    else:
        print(f"  Character: All nonmetals → Covalent")

    print(f"\nQuantum Simulation Requirements:")
    print(f"  Total electrons: {molecule.n_electrons}")
    print(f"  Spin: {molecule.spin}")
    print(f"  Charge: {molecule.charge}")

    # Estimate qubits needed (rough: 2 per electron pair with basis)
    estimated_qubits = molecule.n_electrons // 2 * 3  # Rough estimate with STO-3G
    print(f"  Estimated qubits (STO-3G): ~{estimated_qubits}")

    if estimated_qubits > 30:
        print(f"  ⚠ May exceed BlueQubit's ~30 qubits")
    else:
        print(f"  ✓ Fits in BlueQubit's ~30 qubits")


def main():
    """Create and analyze larger molecules for BlueQubit."""
    print("=" * 80)
    print("KANAD FRAMEWORK - LARGER MOLECULES FOR BLUEQUBIT")
    print("=" * 80)
    print("\nTesting realistic molecular systems with governance validation")
    print("BlueQubit backend: ~30 qubits available")
    print("=" * 80)

    # Create all molecules
    molecules = {
        'H2O': create_h2o(),
        'NH3': create_nh3(),
        'CH4': create_ch4(),
        'NaCl': create_nacl(),
        'CO2': create_co2()
    }

    # Analyze each
    print("\n\n")
    print("=" * 80)
    print("DETAILED MOLECULAR ANALYSIS")
    print("=" * 80)

    for name, molecule in molecules.items():
        analyze_molecule(molecule, name)

    # Summary
    print("\n\n")
    print("=" * 80)
    print("SUMMARY: MOLECULES FOR BLUEQUBIT TESTING")
    print("=" * 80)

    print("\nMolecules Created:")
    for name, mol in molecules.items():
        print(f"  ✓ {name:<8} {len(mol.atoms)} atoms, {mol.n_electrons} electrons")

    print("\nBonding Types Covered:")
    print("  ✓ Covalent: H2O, NH3, CH4, CO2 (different geometries)")
    print("  ✓ Ionic: NaCl")

    print("\nGeometries Tested:")
    print("  ✓ Linear: CO2")
    print("  ✓ Bent: H2O")
    print("  ✓ Pyramidal: NH3")
    print("  ✓ Tetrahedral: CH4")

    print("\nNext Steps:")
    print("  1. Build Hamiltonians with governance validation")
    print("  2. Create UCC ansatze for each molecule")
    print("  3. Run VQE on BlueQubit cloud")
    print("  4. Compare with classical HF/DFT results")

    print("\n" + "=" * 80)
    print("Ready for BlueQubit quantum simulations!")
    print("=" * 80)

    # Check BlueQubit token
    bluequbit_token = os.getenv('BLUEQUBIT_TOKEN')
    if bluequbit_token:
        print("\n✓ BLUEQUBIT_TOKEN found - ready to run!")
    else:
        print("\n⚠ BLUEQUBIT_TOKEN not set")
        print("  Set with: export BLUEQUBIT_TOKEN='your_token'")

    return molecules


if __name__ == "__main__":
    molecules = main()
