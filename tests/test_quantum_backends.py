#!/usr/bin/env python3
"""
Test Quantum Backends with Enhanced Data Storage
Tests IBM Quantum and BlueQubit backends to verify enhanced data is stored correctly
"""
import sys
sys.path.insert(0, '/home/mk/deeprealm/kanad')

import requests
import json
import time
from pprint import pprint


def submit_experiment(backend, molecule_smiles='[H][H]', backend_name=None, device=None):
    """Submit experiment to specified backend"""

    payload = {
        "molecule": {
            "smiles": molecule_smiles,
            "basis": "sto-3g",
            "charge": 0,
            "multiplicity": 1
        },
        "configuration": {
            "method": "VQE",
            "ansatz": "ucc",
            "backend": backend,
            "max_iterations": 3,  # Keep it short for quantum hardware
            "analysis": {
                "energy_decomposition": True,
                "bond_analysis": True,
                "dipole_moment": True
            }
        },
        "execute_now": True
    }

    # Add backend-specific config
    if backend == "ibm_quantum":
        payload["configuration"]["backend_name"] = backend_name or "ibm_torino"
    elif backend == "bluequbit":
        payload["configuration"]["bluequbit_device"] = device or "cpu"

    print(f"\n📤 Submitting experiment to {backend}...")
    print(f"   Molecule: {molecule_smiles}")
    print(f"   Backend config: {backend_name or device or 'default'}")

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


def wait_for_completion(exp_id, timeout=300):
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

        time.sleep(3)

    print(f"⏰ Timeout after {timeout} seconds")
    return False


def verify_enhanced_data(exp_id):
    """Verify enhanced data fields are present and valid"""
    print(f"\n🔍 Verifying enhanced data for experiment {exp_id}...")

    response = requests.get(f'http://localhost:8000/api/experiments/{exp_id}')
    if response.status_code != 200:
        print(f"❌ Failed to get experiment: {response.status_code}")
        return False

    exp = response.json()['experiment']
    results = exp['results']

    if not results:
        print("❌ No results found")
        return False

    print(f"\n✅ Experiment Results:")
    print(f"   Energy: {results.get('energy', 'N/A'):.6f} Ha")
    print(f"   HF Energy: {results.get('hf_energy', 'N/A'):.6f} Ha")
    print(f"   Converged: {results.get('converged', 'N/A')}")
    print(f"   Iterations: {results.get('iterations', 'N/A')}")

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
            elif field == 'rdm1':
                print(f"   ✅ {field}: {len(results[field])}x{len(results[field])} matrix")
            elif field == 'orbital_energies':
                orbs = results[field]
                print(f"   ✅ {field}: {len(orbs)} orbitals")
                if len(orbs) >= 2:
                    # Find HOMO and LUMO
                    print(f"      HOMO: {orbs[-1]:.3f} eV")
                    print(f"      LUMO: {orbs[0]:.3f} eV" if len(orbs) > 1 else "")
            elif field == 'dipole':
                dipole_mag = (results[field][0]**2 + results[field][1]**2 + results[field][2]**2)**0.5
                print(f"   ✅ {field}: {dipole_mag:.3f} D")
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
        if 'bonding' in analysis:
            print(f"   ✅ Bonding analysis present")
            if 'bond_type' in analysis['bonding']:
                gap = analysis['bonding']['bond_type'].get('homo_lumo_gap_ev')
                if gap:
                    print(f"      HOMO-LUMO gap: {gap:.2f} eV")
        if 'properties' in analysis:
            print(f"   ✅ Properties analysis present")
            dipole = analysis['properties'].get('dipole_moment')
            if dipole:
                print(f"      Dipole moment: {dipole:.3f} D")
    else:
        print(f"\n⚠️  No inline analysis results")

    return all_present


def test_ibm_quantum():
    """Test with IBM Quantum backend"""
    print("\n" + "=" * 80)
    print("Testing IBM Quantum Backend")
    print("=" * 80)

    exp_id = submit_experiment(
        backend="ibm_quantum",
        molecule_smiles="[H][H]",  # H2 - small molecule
        backend_name="ibm_torino"
    )

    if not exp_id:
        return False

    # Wait for completion (IBM Quantum can take several minutes)
    success = wait_for_completion(exp_id, timeout=600)  # 10 minutes max

    if not success:
        return False

    # Verify enhanced data
    return verify_enhanced_data(exp_id)


def test_bluequbit():
    """Test with BlueQubit backend"""
    print("\n\n" + "=" * 80)
    print("Testing BlueQubit Backend")
    print("=" * 80)

    exp_id = submit_experiment(
        backend="bluequbit",
        molecule_smiles="[H][H]",  # H2 - small molecule
        device="cpu"
    )

    if not exp_id:
        return False

    # Wait for completion (BlueQubit CPU is usually fast)
    success = wait_for_completion(exp_id, timeout=300)  # 5 minutes max

    if not success:
        return False

    # Verify enhanced data
    return verify_enhanced_data(exp_id)


if __name__ == '__main__':
    print("\n🧪 Quantum Backend Testing Suite\n")
    print("This will test enhanced data storage with real quantum backends")
    print("Note: These tests may take several minutes to complete\n")

    results = {}

    try:
        # Test IBM Quantum
        print("\n🚀 Starting IBM Quantum test...")
        results['ibm_quantum'] = test_ibm_quantum()

        # Test BlueQubit
        print("\n🚀 Starting BlueQubit test...")
        results['bluequbit'] = test_bluequbit()

        # Summary
        print("\n\n" + "=" * 80)
        print("Test Summary")
        print("=" * 80)

        for backend, passed in results.items():
            status = "✅ PASSED" if passed else "❌ FAILED"
            print(f"{backend}: {status}")

        all_passed = all(results.values())

        if all_passed:
            print("\n✅ All quantum backend tests passed!")
            print("Enhanced data storage is working correctly with quantum hardware!")
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
