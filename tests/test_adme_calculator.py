#!/usr/bin/env python3
"""
Test ADME Calculator
Tests the ADME property prediction for drug discovery
"""
import sys
sys.path.insert(0, '/home/mk/deeprealm/kanad')

from kanad.analysis.adme_calculator import ADMECalculator
import json


def test_aspirin():
    """Test ADME calculation for Aspirin"""
    print("=" * 80)
    print("Testing ADME Calculator with Aspirin (Acetylsalicylic Acid)")
    print("=" * 80)

    # Aspirin geometry (simplified)
    geometry = [
        ('C', (0.0, 0.0, 0.0)),
        ('C', (1.4, 0.0, 0.0)),
        ('C', (2.1, 1.2, 0.0)),
        ('C', (1.4, 2.4, 0.0)),
        ('C', (0.0, 2.4, 0.0)),
        ('C', (-0.7, 1.2, 0.0)),
        ('O', (1.5, 2.3, 0.0)),
        ('C', (3.4, 1.2, 0.0)),
        ('O', (1.5, -1.2, 0.0)),
    ]

    # Create calculator
    calc = ADMECalculator(
        geometry=geometry,
        smiles="CC(=O)Oc1ccccc1C(=O)O",  # Aspirin SMILES
        charge=0
    )

    print("\n📊 Calculating molecular descriptors...")
    # Calculate descriptors (without quantum data for now)
    descriptors = calc.calculate_descriptors()

    print("\n✅ Molecular Descriptors:")
    print(f"   MW: {descriptors.molecular_weight:.2f} g/mol")
    print(f"   Heavy atoms: {descriptors.heavy_atom_count}")
    print(f"   H-bond donors: {descriptors.h_bond_donors}")
    print(f"   H-bond acceptors: {descriptors.h_bond_acceptors}")
    print(f"   Rotatable bonds: {descriptors.rotatable_bonds}")
    print(f"   Aromatic rings: {descriptors.aromatic_rings}")
    print(f"   TPSA: {descriptors.polar_surface_area:.2f} Ų" if descriptors.polar_surface_area else "   TPSA: N/A")

    print("\n📊 Predicting ADME properties...")
    # Predict ADME properties
    adme = calc.predict_adme(descriptors)

    print("\n✅ ADME Properties:")
    print(f"\n   Absorption:")
    print(f"   - logP (lipophilicity): {adme.logP:.2f}")
    print(f"   - logS (solubility): {adme.logS:.2f}")
    print(f"   - Caco-2 permeability: {adme.caco2_permeability:.2e} cm/s")
    print(f"   - PAMPA permeability: {adme.pampa_permeability:.2e} cm/s")
    print(f"   - Classification: {adme.absorption_class}")

    print(f"\n   Distribution:")
    print(f"   - logBB (BBB penetration): {adme.log_bb:.2f}")
    print(f"   - BBB penetration: {adme.bbb_penetration}")
    print(f"   - Plasma protein binding: {adme.plasma_protein_binding:.1f}%")
    print(f"   - P-gp substrate: {adme.pgp_substrate}")

    print(f"\n   Drug-Likeness:")
    print(f"   - Lipinski violations: {adme.lipinski_violations}/4")
    print(f"   - Veber violations: {adme.veber_violations}/2")
    print(f"   - Ghose violations: {adme.ghose_violations}/4")

    print("\n📋 Generating full report...")
    # Generate comprehensive report
    report = calc.generate_report(descriptors, adme)

    print("\n✅ Overall Assessment:")
    print(f"   {report['overall_assessment']}")

    # Save report
    with open('/tmp/adme_test_aspirin.json', 'w') as f:
        json.dump(report, f, indent=2)
    print("\n💾 Report saved to: /tmp/adme_test_aspirin.json")

    return report


def test_caffeine():
    """Test ADME calculation for Caffeine"""
    print("\n\n" + "=" * 80)
    print("Testing ADME Calculator with Caffeine")
    print("=" * 80)

    # Caffeine geometry (simplified)
    geometry = [
        ('C', (0.0, 0.0, 0.0)),
        ('N', (1.4, 0.0, 0.0)),
        ('C', (2.1, 1.2, 0.0)),
        ('N', (1.4, 2.4, 0.0)),
        ('C', (0.0, 2.4, 0.0)),
        ('N', (-0.7, 1.2, 0.0)),
        ('C', (3.5, 1.2, 0.0)),
        ('O', (4.2, 2.3, 0.0)),
        ('N', (4.2, 0.0, 0.0)),
    ]

    calc = ADMECalculator(
        geometry=geometry,
        smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # Caffeine SMILES
        charge=0
    )

    descriptors = calc.calculate_descriptors()
    adme = calc.predict_adme(descriptors)

    print(f"\n✅ Caffeine ADME Summary:")
    print(f"   logP: {adme.logP:.2f}, BBB: {adme.bbb_penetration}")
    print(f"   Lipinski violations: {adme.lipinski_violations}/4")
    print(f"   Absorption: {adme.absorption_class}")

    return adme


def test_with_quantum_data():
    """Test ADME with quantum-enhanced descriptors"""
    print("\n\n" + "=" * 80)
    print("Testing ADME Calculator with Quantum Data")
    print("=" * 80)

    # Water molecule (simple test)
    geometry = [
        ('O', (0.0, 0.0, 0.0)),
        ('H', (0.96, 0.0, 0.0)),
        ('H', (-0.24, 0.93, 0.0))
    ]

    calc = ADMECalculator(geometry=geometry, smiles="O", charge=0)

    # Simulate quantum data (from VQE results)
    import numpy as np
    rdm1 = np.eye(7) * 0.5  # Simplified density matrix
    orbital_energies = [-20.24, -1.26, -0.61, -0.45, -0.39, 0.59, 0.73]
    dipole = [-0.03, -1.71, 0.0]
    polarizability = 1.5

    print("\n📊 Calculating with quantum descriptors...")
    descriptors = calc.calculate_descriptors(
        rdm1=rdm1,
        orbital_energies=orbital_energies,
        dipole=dipole,
        polarizability=polarizability
    )

    print(f"\n✅ Quantum-enhanced descriptors:")
    print(f"   HOMO: {descriptors.homo_energy:.2f} eV")
    print(f"   LUMO: {descriptors.lumo_energy:.2f} eV")
    print(f"   HOMO-LUMO gap: {descriptors.homo_lumo_gap:.2f} eV")
    print(f"   Dipole moment: {descriptors.dipole_moment:.2f} D")
    print(f"   Polarizability: {descriptors.polarizability:.2f} ų")

    return descriptors


if __name__ == '__main__':
    print("\n🧪 ADME Calculator Test Suite\n")

    try:
        # Test 1: Aspirin
        aspirin_report = test_aspirin()

        # Test 2: Caffeine
        caffeine_adme = test_caffeine()

        # Test 3: Quantum data
        quantum_desc = test_with_quantum_data()

        print("\n\n" + "=" * 80)
        print("✅ All ADME tests passed!")
        print("=" * 80)

    except Exception as e:
        print(f"\n\n❌ Test failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
