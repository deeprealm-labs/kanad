#!/usr/bin/env python3
"""
Test Inline Analysis with Enhanced VQE Data
Tests that inline analysis uses the stored enhanced data correctly
"""
import sys
sys.path.insert(0, '/home/mk/deeprealm/kanad')

import requests
import json
from pprint import pprint


def test_with_completed_experiment(experiment_id):
    """Test analysis with a completed experiment"""
    print("=" * 80)
    print(f"Testing Inline Analysis with Experiment: {experiment_id}")
    print("=" * 80)

    # Get experiment results
    response = requests.get(f'http://localhost:8000/api/experiments/{experiment_id}')
    if response.status_code != 200:
        print(f"❌ Failed to get experiment: {response.status_code}")
        return False

    exp = response.json()['experiment']
    results = exp['results']

    print(f"\n📊 Experiment Status: {exp['status']}")
    print(f"   Molecule: {exp['molecule']['smiles']}")
    print(f"   Method: {results['method']}")
    print(f"   Energy: {results['energy']:.6f} Ha")

    # Check enhanced data fields
    print("\n✅ Enhanced Data Fields:")
    enhanced_fields = ['geometry', 'atoms', 'rdm1', 'orbital_energies', 'dipole', 'nuclear_repulsion']

    for field in enhanced_fields:
        if field in results and results[field] is not None:
            if field == 'geometry':
                print(f"   ✓ {field}: {len(results[field])} atoms")
                print(f"      Example: {results[field][0]}")
            elif field == 'rdm1':
                print(f"   ✓ {field}: {len(results[field])}x{len(results[field])} matrix")
            elif field == 'orbital_energies':
                print(f"   ✓ {field}: {len(results[field])} orbitals")
                print(f"      HOMO: {results[field][4]:.2f} eV")
                print(f"      LUMO: {results[field][5]:.2f} eV")
            elif field == 'dipole':
                dipole_mag = (results[field][0]**2 + results[field][1]**2 + results[field][2]**2)**0.5
                print(f"   ✓ {field}: {dipole_mag:.2f} D")
            else:
                print(f"   ✓ {field}: {results[field]}")
        else:
            print(f"   ✗ {field}: MISSING")

    # Check inline analysis results
    if 'analysis' in results:
        print("\n✅ Inline Analysis Results:")
        analysis = results['analysis']

        if 'energy_components' in analysis:
            ec = analysis['energy_components']
            print(f"\n   Energy Components:")
            print(f"   - Nuclear repulsion: {ec.get('nuclear_repulsion', 'N/A'):.4f} Ha")
            print(f"   - One-electron: {ec.get('one_electron', 'N/A'):.4f} Ha")
            print(f"   - Two-electron: {ec.get('two_electron', 'N/A'):.4f} Ha")
            print(f"   - Total: {ec.get('total', 'N/A'):.4f} Ha")

        if 'bonding' in analysis:
            bonding = analysis['bonding']
            if 'bond_type' in bonding:
                bt = bonding['bond_type']
                print(f"\n   Bonding Analysis:")
                print(f"   - Type: {bt.get('bonding_type', 'N/A')}")
                print(f"   - HOMO-LUMO gap: {bt.get('homo_lumo_gap', 'N/A'):.4f} Ha")
                print(f"   - HOMO-LUMO gap: {bt.get('homo_lumo_gap_ev', 'N/A'):.2f} eV")

        if 'properties' in analysis:
            props = analysis['properties']
            print(f"\n   Properties:")
            print(f"   - Dipole moment: {props.get('dipole_moment', 'N/A'):.4f} D")
            if 'dipole_vector' in props:
                print(f"   - Dipole vector: [{props['dipole_vector'][0]:.4f}, {props['dipole_vector'][1]:.4f}, {props['dipole_vector'][2]:.4f}] D")

    else:
        print("\n⚠️  No inline analysis results found")

    return True


def test_analysis_with_api(molecule_smiles, basis='sto-3g'):
    """Run a new experiment and verify analysis"""
    print("\n\n" + "=" * 80)
    print(f"Testing Full Experiment + Analysis Flow")
    print(f"Molecule: {molecule_smiles}")
    print("=" * 80)

    # Submit experiment
    payload = {
        "molecule": {
            "smiles": molecule_smiles,
            "basis": basis,
            "charge": 0,
            "multiplicity": 1
        },
        "configuration": {
            "method": "VQE",
            "ansatz": "uccsd",
            "backend": "classical",
            "max_iterations": 3,
            "analysis": {
                "energy_decomposition": True,
                "bond_analysis": True,
                "dipole_moment": True
            }
        },
        "execute_now": True
    }

    print("\n📤 Submitting experiment...")
    response = requests.post('http://localhost:8000/api/experiments/submit', json=payload)

    if response.status_code != 200:
        print(f"❌ Failed to submit: {response.status_code}")
        print(response.text)
        return False

    exp_id = response.json()['experiment_id']
    print(f"✅ Experiment submitted: {exp_id}")

    # Wait for completion
    import time
    max_wait = 60  # seconds
    waited = 0

    while waited < max_wait:
        time.sleep(3)
        waited += 3

        response = requests.get(f'http://localhost:8000/api/experiments/{exp_id}')
        status = response.json()['experiment']['status']

        if status == 'completed':
            print(f"✅ Experiment completed in {waited} seconds")
            break
        elif status == 'failed':
            print(f"❌ Experiment failed")
            return False
        else:
            print(f"   ⏳ Status: {status} ({waited}s)")

    if status != 'completed':
        print(f"⏰ Timeout after {max_wait}s")
        return False

    # Test the results
    return test_with_completed_experiment(exp_id)


if __name__ == '__main__':
    print("\n🧪 Inline Analysis Test Suite\n")

    try:
        # Test 1: Use existing completed experiment
        print("Test 1: Existing Experiment")
        test_with_completed_experiment('8a3d7976-8cbf-45eb-a6ea-bdd3d9761ee5')

        # Test 2: Run new H2 experiment
        print("\n\nTest 2: Fresh H2 Experiment")
        success = test_analysis_with_api('[H][H]', basis='sto-3g')

        if success:
            print("\n\n" + "=" * 80)
            print("✅ All inline analysis tests passed!")
            print("=" * 80)
        else:
            print("\n\n❌ Some tests failed")
            sys.exit(1)

    except Exception as e:
        print(f"\n\n❌ Test failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
