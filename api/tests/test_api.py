"""
Test script for Kanad API endpoints.

Run this after starting the API server to verify all endpoints work correctly.
"""

import requests
import json
import time
from typing import Dict, Any

# API base URL
BASE_URL = "http://localhost:8000"
API_V1 = f"{BASE_URL}/api/v1"


def print_response(title: str, response: requests.Response):
    """Print formatted response."""
    print(f"\n{'='*80}")
    print(f"{title}")
    print(f"{'='*80}")
    print(f"Status Code: {response.status_code}")
    print(f"Response:")
    print(json.dumps(response.json(), indent=2))


def test_health_check():
    """Test health check endpoint."""
    response = requests.get(f"{BASE_URL}/health")
    print_response("Health Check", response)
    assert response.status_code == 200


def test_api_info():
    """Test API info endpoint."""
    response = requests.get(f"{API_V1}/info")
    print_response("API Info", response)
    assert response.status_code == 200


def test_validate_smiles():
    """Test SMILES validation."""
    # Valid SMILES
    response = requests.post(
        f"{API_V1}/molecules/validate",
        json={"smiles": "CCO"}
    )
    print_response("Validate SMILES - Valid", response)
    assert response.status_code == 200
    assert response.json()["valid"] is True

    # Invalid SMILES
    response = requests.post(
        f"{API_V1}/molecules/validate",
        json={"smiles": "invalid_smiles"}
    )
    print_response("Validate SMILES - Invalid", response)
    assert response.status_code == 200
    assert response.json()["valid"] is False


def test_molecule_library():
    """Test molecule library endpoints."""
    # Get full library
    response = requests.get(f"{API_V1}/molecules/library")
    print_response("Molecule Library", response)
    assert response.status_code == 200

    # Get categories
    response = requests.get(f"{API_V1}/molecules/library/categories")
    print_response("Molecule Categories", response)
    assert response.status_code == 200

    # Get specific molecule
    response = requests.get(f"{API_V1}/molecules/library/h2")
    print_response("Get H2 from Library", response)
    assert response.status_code == 200


def test_settings():
    """Test settings endpoints."""
    # Get settings
    response = requests.get(f"{API_V1}/settings")
    print_response("Get Settings", response)
    assert response.status_code == 200

    # Update settings
    response = requests.put(
        f"{API_V1}/settings",
        json={
            "method": "VQE",
            "ansatz": "hardware_efficient",
            "optimizer": "COBYLA",
            "circuit_optimization": True
        }
    )
    print_response("Update Settings", response)
    assert response.status_code == 200


def test_create_experiment():
    """Test experiment creation."""
    experiment_data = {
        "name": "H2 VQE Test",
        "molecule": {
            "smiles": "[H][H]",
            "basis": "sto-3g",
            "charge": 0,
            "multiplicity": 1
        },
        "configuration": {
            "method": "VQE",
            "ansatz": "ucc",
            "mapper": "jordan_wigner",
            "optimizer": "SLSQP",
            "backend": "classical",
            "max_iterations": 100,
            "conv_threshold": 1e-6
        },
        "execute_immediately": True
    }

    response = requests.post(f"{API_V1}/experiments/", json=experiment_data)
    print_response("Create Experiment", response)
    assert response.status_code == 201

    return response.json()["id"]


def test_list_experiments():
    """Test experiment listing."""
    response = requests.get(f"{API_V1}/experiments/")
    print_response("List Experiments", response)
    assert response.status_code == 200


def test_get_experiment(experiment_id: int):
    """Test getting experiment details."""
    # Poll until experiment completes (with timeout)
    max_wait = 120  # 2 minutes
    start_time = time.time()

    while time.time() - start_time < max_wait:
        response = requests.get(f"{API_V1}/experiments/{experiment_id}")
        data = response.json()

        print(f"\nExperiment Status: {data['status']}")

        if data['status'] in ['completed', 'failed']:
            print_response(f"Experiment {experiment_id} - Final", response)
            break

        time.sleep(2)  # Wait 2 seconds before next poll

    assert response.status_code == 200


def test_convergence_data(experiment_id: int):
    """Test getting convergence data."""
    response = requests.get(f"{API_V1}/experiments/{experiment_id}/convergence")
    print_response(f"Convergence Data - Experiment {experiment_id}", response)
    assert response.status_code == 200


def test_queue_operations():
    """Test queue operations."""
    # Create experiment without immediate execution
    experiment_data = {
        "name": "Queued Experiment",
        "molecule": {
            "smiles": "O",
            "basis": "sto-3g",
            "charge": 0,
            "multiplicity": 1
        },
        "configuration": {
            "method": "VQE",
            "ansatz": "ucc",
            "mapper": "jordan_wigner",
            "optimizer": "SLSQP",
            "backend": "classical"
        },
        "execute_immediately": False
    }

    response = requests.post(f"{API_V1}/experiments/", json=experiment_data)
    experiment_id = response.json()["id"]

    # Add to queue
    queue_data = {
        "experiment_id": experiment_id,
        "priority": 5
    }
    response = requests.post(f"{API_V1}/queue/", json=queue_data)
    print_response("Add to Queue", response)
    assert response.status_code == 201
    queue_id = response.json()["id"]

    # List queue
    response = requests.get(f"{API_V1}/queue/")
    print_response("List Queue", response)
    assert response.status_code == 200

    # Update queue item
    response = requests.put(
        f"{API_V1}/queue/{queue_id}",
        json={"priority": 10}
    )
    print_response("Update Queue Priority", response)
    assert response.status_code == 200

    # Delete queue item
    response = requests.delete(f"{API_V1}/queue/{queue_id}")
    print_response("Delete Queue Item", response)
    assert response.status_code == 200


def run_all_tests():
    """Run all API tests."""
    print("\n" + "="*80)
    print("KANAD API TEST SUITE")
    print("="*80)

    try:
        # Basic endpoints
        test_health_check()
        test_api_info()

        # Molecules
        test_validate_smiles()
        test_molecule_library()

        # Settings
        test_settings()

        # Experiments
        test_list_experiments()
        experiment_id = test_create_experiment()

        print(f"\n\nWaiting for experiment {experiment_id} to complete...")
        test_get_experiment(experiment_id)
        test_convergence_data(experiment_id)

        # Queue
        test_queue_operations()

        print("\n" + "="*80)
        print("ALL TESTS PASSED!")
        print("="*80)

    except AssertionError as e:
        print(f"\n\nTEST FAILED: {e}")
    except Exception as e:
        print(f"\n\nERROR: {e}")


if __name__ == "__main__":
    # Check if server is running
    try:
        response = requests.get(f"{BASE_URL}/health", timeout=5)
        if response.status_code == 200:
            print("Server is running. Starting tests...")
            run_all_tests()
        else:
            print("Server returned unexpected status. Please check if it's running correctly.")
    except requests.exceptions.ConnectionError:
        print("\nERROR: Cannot connect to API server.")
        print(f"Please start the server first:")
        print(f"  cd /home/mk/deeprealm/kanad/api")
        print(f"  python main.py")
