#!/usr/bin/env python3
"""
Comprehensive API test script.

Tests all endpoints to ensure complete implementation.
"""

import requests
import json
import sys

BASE_URL = "http://localhost:8000"
API_V1 = f"{BASE_URL}/api/v1"


def print_test(name, passed, details=""):
    """Print test result."""
    status = "✓" if passed else "✗"
    print(f"{status} {name}")
    if details and not passed:
        print(f"  Error: {details}")


def test_info_endpoint():
    """Test /api/v1/info endpoint."""
    print("\n=== Testing Info Endpoint ===")

    try:
        response = requests.get(f"{API_V1}/info")
        data = response.json()

        # Check all required fields
        tests = [
            ("Info endpoint accessible", response.status_code == 200),
            ("Has solvers list", "solvers" in data),
            ("Has ansatze list", "ansatze" in data),
            ("Has mappers list", "mappers" in data),
            ("Has hamiltonians list", "hamiltonians" in data),
            ("Has optimizers list", "optimizers" in data),
            ("Has backends list", "backends" in data),
            ("VQE in solvers", "VQE" in data.get("solvers", [])),
            ("SQD in solvers", "SQD" in data.get("solvers", [])),
            ("ExcitedStates in solvers", "ExcitedStates" in data.get("solvers", [])),
            ("ucc in ansatze", "ucc" in data.get("ansatze", [])),
            ("uccsd in ansatze", "uccsd" in data.get("ansatze", [])),
            ("hardware_efficient in ansatze", "hardware_efficient" in data.get("ansatze", [])),
            ("governance_aware in ansatze", "governance_aware" in data.get("ansatze", [])),
            ("two_local in ansatze", "two_local" in data.get("ansatze", [])),
            ("jordan_wigner in mappers", "jordan_wigner" in data.get("mappers", [])),
            ("bravyi_kitaev in mappers", "bravyi_kitaev" in data.get("mappers", [])),
            ("hybrid_orbital in mappers", "hybrid_orbital" in data.get("mappers", [])),
            ("parity in mappers", "parity" in data.get("mappers", [])),
            ("covalent in hamiltonians", "covalent" in data.get("hamiltonians", [])),
            ("ionic in hamiltonians", "ionic" in data.get("hamiltonians", [])),
            ("metallic in hamiltonians", "metallic" in data.get("hamiltonians", [])),
            ("periodic in hamiltonians", "periodic" in data.get("hamiltonians", [])),
            ("molecular in hamiltonians", "molecular" in data.get("hamiltonians", [])),
            ("classical in backends", "classical" in data.get("backends", [])),
            ("ibm_quantum in backends", "ibm_quantum" in data.get("backends", [])),
            ("bluequbit_gpu in backends", "bluequbit_gpu" in data.get("backends", [])),
        ]

        for name, result in tests:
            print_test(name, result)

        return all(result for _, result in tests)

    except Exception as e:
        print_test("Info endpoint", False, str(e))
        return False


def test_cloud_credentials_endpoints():
    """Test cloud credentials endpoints."""
    print("\n=== Testing Cloud Credentials Endpoints ===")

    try:
        # Test status endpoint
        response = requests.get(f"{API_V1}/cloud-credentials/status")
        print_test("Credentials status endpoint", response.status_code == 200)

        # Test IBM endpoints
        response = requests.get(f"{API_V1}/cloud-credentials/ibm")
        print_test("IBM credentials GET endpoint", response.status_code == 200)

        # Test BlueQubit endpoints
        response = requests.get(f"{API_V1}/cloud-credentials/bluequbit")
        print_test("BlueQubit credentials GET endpoint", response.status_code == 200)

        # Test IBM POST (with dummy data)
        ibm_data = {
            "crn": "test_crn",
            "api_key": "test_key"
        }
        response = requests.post(f"{API_V1}/cloud-credentials/ibm", json=ibm_data)
        print_test("IBM credentials POST endpoint", response.status_code == 200)

        # Test BlueQubit POST (with dummy data)
        bq_data = {
            "api_token": "test_token"
        }
        response = requests.post(f"{API_V1}/cloud-credentials/bluequbit", json=bq_data)
        print_test("BlueQubit credentials POST endpoint", response.status_code == 200)

        # Verify credentials were stored
        response = requests.get(f"{API_V1}/cloud-credentials/ibm")
        data = response.json()
        print_test("IBM credentials stored", data.get("configured") == True)

        response = requests.get(f"{API_V1}/cloud-credentials/bluequbit")
        data = response.json()
        print_test("BlueQubit credentials stored", data.get("configured") == True)

        # Clean up - delete test credentials
        requests.delete(f"{API_V1}/cloud-credentials/ibm")
        requests.delete(f"{API_V1}/cloud-credentials/bluequbit")

        return True

    except Exception as e:
        print_test("Cloud credentials endpoints", False, str(e))
        return False


def test_experiment_creation():
    """Test experiment creation with all solvers."""
    print("\n=== Testing Experiment Creation ===")

    experiments = [
        {
            "name": "Test VQE",
            "solver": "VQE",
            "ansatz": "ucc"
        },
        {
            "name": "Test SQD",
            "solver": "SQD",
            "ansatz": "hardware_efficient"
        },
    ]

    all_passed = True

    for exp_config in experiments:
        try:
            data = {
                "name": exp_config["name"],
                "molecule": {
                    "smiles": "H2",
                    "basis": "sto-3g"
                },
                "configuration": {
                    "method": exp_config["solver"],
                    "ansatz": exp_config["ansatz"],
                    "backend": "classical"
                },
                "execute_immediately": False
            }

            response = requests.post(f"{API_V1}/experiments/", json=data)
            passed = response.status_code == 201
            print_test(f"Create {exp_config['name']}", passed,
                      response.text if not passed else "")

            all_passed = all_passed and passed

        except Exception as e:
            print_test(f"Create {exp_config['name']}", False, str(e))
            all_passed = False

    return all_passed


def test_validators():
    """Test validators accept all valid options."""
    print("\n=== Testing Validators ===")

    all_passed = True

    # Test all ansatze
    ansatze = ["ucc", "uccsd", "hardware_efficient", "governance_aware", "two_local"]
    for ansatz in ansatze:
        try:
            data = {
                "name": f"Test {ansatz}",
                "molecule": {"smiles": "H2", "basis": "sto-3g"},
                "configuration": {
                    "method": "VQE",
                    "ansatz": ansatz,
                    "backend": "classical"
                },
                "execute_immediately": False
            }
            response = requests.post(f"{API_V1}/experiments/", json=data)
            passed = response.status_code == 201
            print_test(f"Ansatz: {ansatz}", passed, response.text if not passed else "")
            all_passed = all_passed and passed
        except Exception as e:
            print_test(f"Ansatz: {ansatz}", False, str(e))
            all_passed = False

    # Test all mappers
    mappers = ["jordan_wigner", "bravyi_kitaev", "hybrid_orbital", "parity"]
    for mapper in mappers:
        try:
            data = {
                "name": f"Test {mapper}",
                "molecule": {"smiles": "H2", "basis": "sto-3g"},
                "configuration": {
                    "method": "VQE",
                    "mapper": mapper,
                    "backend": "classical"
                },
                "execute_immediately": False
            }
            response = requests.post(f"{API_V1}/experiments/", json=data)
            passed = response.status_code == 201
            print_test(f"Mapper: {mapper}", passed, response.text if not passed else "")
            all_passed = all_passed and passed
        except Exception as e:
            print_test(f"Mapper: {mapper}", False, str(e))
            all_passed = False

    return all_passed


def main():
    """Run all tests."""
    print("=" * 60)
    print("Kanad API Complete Implementation Test")
    print("=" * 60)

    # Check if server is running
    try:
        response = requests.get(f"{BASE_URL}/health")
        if response.status_code != 200:
            print("✗ Server is not healthy!")
            sys.exit(1)
        print("✓ Server is running and healthy")
    except requests.exceptions.ConnectionError:
        print("✗ Cannot connect to server at", BASE_URL)
        print("  Please start the server with: uvicorn main:app --reload")
        sys.exit(1)

    # Run tests
    results = []
    results.append(("Info Endpoint", test_info_endpoint()))
    results.append(("Cloud Credentials", test_cloud_credentials_endpoints()))
    results.append(("Experiment Creation", test_experiment_creation()))
    results.append(("Validators", test_validators()))

    # Summary
    print("\n" + "=" * 60)
    print("Test Summary")
    print("=" * 60)

    total = len(results)
    passed = sum(1 for _, result in results if result)

    for name, result in results:
        status = "PASS" if result else "FAIL"
        print(f"{name:30} {status}")

    print(f"\nTotal: {passed}/{total} test suites passed")

    if passed == total:
        print("\n✓ All tests passed!")
        sys.exit(0)
    else:
        print(f"\n✗ {total - passed} test suite(s) failed")
        sys.exit(1)


if __name__ == "__main__":
    main()
