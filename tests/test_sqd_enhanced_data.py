#!/usr/bin/env python3
"""
Test SQD Solver Enhanced Data Storage
Tests that SQD solver stores geometry, rdm1, orbital_energies, and dipole
"""
import sys
sys.path.insert(0, '/home/mk/deeprealm/kanad')

import requests
import json
import time
from pprint import pprint


def submit_sqd_experiment(molecule_smiles='[H][H]', backend='classical'):
    """Submit SQD experiment to test enhanced data storage"""

    payload = {
        "molecule": {
            "smiles": molecule_smiles,
            "basis": "sto-3g",
            "charge": 0,
            "multiplicity": 1
        },
        "configuration": {
            "method": "SQD",
            "backend": backend,
            "subspace_dim": 3,  # Ground state + 2 excited states
            "circuit_depth": 3,
            "max_iterations": 5,  # Keep it short
            "analysis": {
                "energy_decomposition": True,
                "bond_analysis": True,
                "dipole_moment": True
            }
        },
        "execute_now": True
    }

    print(f"\n📤 Submitting SQD experiment...")
    print(f"   Molecule: {molecule_smiles}")
    print(f"   Backend: {backend}")
    print(f"   Subspace dimension: 3 (ground + 2 excited states)")

    response = requests.post('http://localhost:8000/api/experiments/submit', json=payload)

    if response.status_code != 200:
        print(f"❌ Failed to submit: {response.status_code}")
        print(response.text)
        return None

    result = response.json()
    exp_id = result['experiment_id']
    job_id = result['job_id']

    print(f"✅ Experiment submitted:")
    print(f"   Experiment ID: {exp_id}")
    print(f"   Job ID: {job_id}")

    return exp_id


def wait_for_completion(exp_id, timeout=120):
    """Wait for experiment to complete"""
    print(f"\n⏳ Waiting for experiment {exp_id} to complete...")

    start_time = time.time()
    last_status = None

    while time.time() - start_time < timeout:
        response = requests.get(f'http://localhost:8000/api/experiments/{exp_id}')
        if response.status_code != 200:
            print(f"❌ Failed to get status: {response.status_code}")
            return False

        exp = response.json()['experiment']
        status = exp['status']

        # Print status changes
        if status != last_status:
            elapsed = int(time.time() - start_time)
            if status == 'running' and 'job' in exp:
                job = exp['job']
                progress = job.get('progress', 0)
                iteration = job.get('current_iteration', 0)
                energy = job.get('current_energy')
                print(f"   [{elapsed}s] Status: {status} | Progress: {progress:.0f}% | Iter: {iteration} | E: {energy:.4f} Ha" if energy else f"   [{elapsed}s] Status: {status} | Progress: {progress:.0f}%")
            else:
                print(f"   [{elapsed}s] Status: {status}")
            last_status = status

        if status == 'completed':
            elapsed = int(time.time() - start_time)
            print(f"✅ Experiment completed in {elapsed} seconds")
            return True
        elif status == 'failed':
            print(f"❌ Experiment failed")
            if 'error_message' in exp:
                print(f"   Error: {exp['error_message']}")
            return False

        time.sleep(2)

    print(f"⏰ Timeout after {timeout} seconds")
    return False


def verify_enhanced_data(exp_id):
    """Verify enhanced data fields are present and valid"""
    print(f"\n🔍 Verifying SQD enhanced data for experiment {exp_id}...")

    response = requests.get(f'http://localhost:8000/api/experiments/{exp_id}')
    if response.status_code != 200:
        print(f"❌ Failed to get experiment: {response.status_code}")
        return False

    exp = response.json()['experiment']
    results = exp['results']

    if not results:
        print("❌ No results found")
        return False

    print(f"\n✅ SQD Experiment Results:")
    print(f"   Energy: {results.get('energy', 'N/A'):.6f} Ha")
    print(f"   HF Energy: {results.get('hf_energy', 'N/A'):.6f} Ha")
    print(f"   Converged: {results.get('converged', 'N/A')}")
    print(f"   Subspace Dim: {results.get('subspace_dim', 'N/A')}")

    # Show excited states if present
    if 'excited_state_energies' in results and results['excited_state_energies']:
        print(f"   Excited States:")
        for i, energy in enumerate(results['excited_state_energies']):
            print(f"      State {i+1}: {energy:.6f} Ha")

    # Check enhanced data fields
    enhanced_fields = {
        'geometry': 'Molecular geometry',
        'atoms': 'Atom symbols',
        'rdm1': 'Density matrix',
        'orbital_energies': 'Orbital energies',
        'dipole': 'Dipole moment',
        'nuclear_repulsion': 'Nuclear repulsion'
    }

    print(f"\n📊 Enhanced Data Fields:")
    all_present = True

    for field, description in enhanced_fields.items():
        if field in results and results[field] is not None:
            if field == 'geometry':
                print(f"   ✅ {field}: {len(results[field])} atoms")
                # Show geometry
                for atom_symbol, coords in results[field]:
                    print(f"      {atom_symbol}: ({coords[0]:.4f}, {coords[1]:.4f}, {coords[2]:.4f})")
            elif field == 'rdm1':
                print(f"   ✅ {field}: {len(results[field])}x{len(results[field])} matrix")
            elif field == 'orbital_energies':
                orbs = results[field]
                print(f"   ✅ {field}: {len(orbs)} orbitals")
                if len(orbs) >= 2:
                    # Find HOMO and LUMO
                    print(f"      HOMO: {orbs[-1]:.3f} eV")
                    print(f"      LUMO: {orbs[0]:.3f} eV" if len(orbs) > 1 else "")
                    # Show gap
                    if len(orbs) > 1:
                        gap = abs(orbs[0] - orbs[-1])
                        print(f"      HOMO-LUMO gap: {gap:.3f} eV")
            elif field == 'dipole':
                dipole_mag = (results[field][0]**2 + results[field][1]**2 + results[field][2]**2)**0.5
                print(f"   ✅ {field}: {dipole_mag:.3f} D")
                print(f"      Vector: ({results[field][0]:.3f}, {results[field][1]:.3f}, {results[field][2]:.3f})")
            elif field == 'atoms':
                print(f"   ✅ {field}: {results[field]}")
            else:
                print(f"   ✅ {field}: {results[field]}")
        else:
            print(f"   ❌ {field}: MISSING")
            all_present = False

    # Check inline analysis
    if 'analysis' in results and results['analysis']:
        print(f"\n📈 Inline Analysis:")
        analysis = results['analysis']

        if 'energy_components' in analysis:
            print(f"   ✅ Energy decomposition present")
            components = analysis['energy_components']
            if 'nuclear_repulsion' in components:
                print(f"      Nuclear repulsion: {components['nuclear_repulsion']:.4f} Ha")
            if 'one_electron' in components:
                print(f"      One-electron: {components['one_electron']:.4f} Ha")
            if 'two_electron' in components:
                print(f"      Two-electron: {components['two_electron']:.4f} Ha")

        if 'bonding' in analysis:
            print(f"   ✅ Bonding analysis present")
            if 'bond_type' in analysis['bonding']:
                bond_type = analysis['bonding']['bond_type']
                if 'homo_lumo_gap_ev' in bond_type:
                    print(f"      HOMO-LUMO gap: {bond_type['homo_lumo_gap_ev']:.2f} eV")
                if 'type' in bond_type:
                    print(f"      Bond type: {bond_type['type']}")

        if 'properties' in analysis:
            print(f"   ✅ Properties analysis present")
            props = analysis['properties']
            if 'dipole_moment' in props:
                print(f"      Dipole moment: {props['dipole_moment']:.3f} D")
    else:
        print(f"\n⚠️  No inline analysis results")

    return all_present


def test_sqd_h2():
    """Test SQD with H2 molecule"""
    print("\n" + "=" * 80)
    print("Testing SQD Enhanced Data Storage - H2 Molecule")
    print("=" * 80)

    exp_id = submit_sqd_experiment(
        molecule_smiles="[H][H]",
        backend="classical"
    )

    if not exp_id:
        return False

    # Wait for completion
    success = wait_for_completion(exp_id, timeout=120)

    if not success:
        return False

    # Verify enhanced data
    return verify_enhanced_data(exp_id)


def test_sqd_h2o():
    """Test SQD with H2O molecule"""
    print("\n\n" + "=" * 80)
    print("Testing SQD Enhanced Data Storage - H2O Molecule")
    print("=" * 80)

    exp_id = submit_sqd_experiment(
        molecule_smiles="O",
        backend="classical"
    )

    if not exp_id:
        return False

    # Wait for completion (H2O may take longer)
    success = wait_for_completion(exp_id, timeout=180)

    if not success:
        return False

    # Verify enhanced data
    return verify_enhanced_data(exp_id)


if __name__ == '__main__':
    print("\n🧪 SQD Enhanced Data Storage Test Suite\n")
    print("This will test that SQD solver stores enhanced data correctly")
    print("Tests will run with classical backend for speed\n")

    results = {}

    try:
        # Test H2 (small molecule)
        print("\n🚀 Starting H2 test...")
        results['h2'] = test_sqd_h2()

        # Test H2O (slightly larger)
        print("\n🚀 Starting H2O test...")
        results['h2o'] = test_sqd_h2o()

        # Summary
        print("\n\n" + "=" * 80)
        print("Test Summary")
        print("=" * 80)

        for molecule, passed in results.items():
            status = "✅ PASSED" if passed else "❌ FAILED"
            print(f"{molecule.upper()}: {status}")

        all_passed = all(results.values())

        if all_passed:
            print("\n✅ All SQD enhanced data tests passed!")
            print("SQD solver is correctly storing geometry, rdm1, orbital_energies, and dipole!")
            sys.exit(0)
        else:
            print("\n⚠️  Some tests failed - review the output above")
            sys.exit(1)

    except KeyboardInterrupt:
        print("\n\n⚠️  Tests interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n\n❌ Test failed with exception: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
