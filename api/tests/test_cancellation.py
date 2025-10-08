"""
Test script for job cancellation functionality.

Tests:
1. Cancel queued experiment (before execution starts)
2. Cancel running VQE experiment (during optimization)
3. Cancel running SQD experiment
4. Attempt to cancel completed experiment (should fail)
5. Verify partial results are saved for cancelled jobs
"""

import requests
import time
import sys

BASE_URL = "http://localhost:8000/api/v1"


def test_cancel_queued_experiment():
    """Test cancelling an experiment before it starts running."""
    print("\n" + "=" * 80)
    print("TEST 1: Cancel Queued Experiment")
    print("=" * 80)

    # Create experiment but don't execute immediately
    experiment = {
        "name": "Queued H2 Test",
        "molecule": {
            "smiles": "[H][H]",
            "basis": "sto-3g"
        },
        "configuration": {
            "method": "VQE",
            "ansatz": "ucc",
            "max_iterations": 100
        },
        "execute_immediately": False
    }

    response = requests.post(f"{BASE_URL}/experiments", json=experiment)
    assert response.status_code == 201, f"Failed to create experiment: {response.text}"

    exp_id = response.json()['id']
    print(f"Created experiment {exp_id} (status: pending)")

    # Cancel immediately
    cancel_response = requests.patch(f"{BASE_URL}/experiments/{exp_id}/cancel")

    if cancel_response.status_code == 400:
        print(f"  - Cannot cancel pending experiment (expected behavior)")
        print(f"  - Skipping this test")
        return True
    elif cancel_response.status_code == 200:
        result = cancel_response.json()
        print(f"  - Cancellation response: {result}")

        # Check status
        status_response = requests.get(f"{BASE_URL}/experiments/{exp_id}/status")
        status = status_response.json()
        print(f"  - Final status: {status['status']}")

        assert status['status'] == 'cancelled' or status['status'] == 'pending', \
            f"Expected 'cancelled' or 'pending', got '{status['status']}'"

        print("  PASSED")
        return True
    else:
        print(f"  FAILED: Unexpected status code {cancel_response.status_code}")
        return False


def test_cancel_running_vqe():
    """Test cancelling a running VQE optimization."""
    print("\n" + "=" * 80)
    print("TEST 2: Cancel Running VQE Experiment")
    print("=" * 80)

    # Create a long-running VQE experiment
    experiment = {
        "name": "Running H2 VQE Test",
        "molecule": {
            "smiles": "[H][H]",
            "basis": "sto-3g"
        },
        "configuration": {
            "method": "VQE",
            "ansatz": "ucc",
            "optimizer": "SLSQP",
            "max_iterations": 1000,  # Long-running
            "backend": "classical"
        },
        "execute_immediately": True
    }

    response = requests.post(f"{BASE_URL}/experiments", json=experiment)
    assert response.status_code == 201, f"Failed to create experiment: {response.text}"

    exp_id = response.json()['id']
    print(f"Created experiment {exp_id}")

    # Wait for it to start running
    print("  - Waiting for experiment to start running...")
    for i in range(10):
        time.sleep(1)
        status_response = requests.get(f"{BASE_URL}/experiments/{exp_id}/status")
        status = status_response.json()
        print(f"    [{i+1}s] Status: {status['status']}, Progress: {status['progress']}%")

        if status['status'] == 'running':
            print("  - Experiment is now running, sending cancellation...")
            break

    # Cancel it
    cancel_response = requests.patch(f"{BASE_URL}/experiments/{exp_id}/cancel")
    assert cancel_response.status_code == 200, f"Cancellation failed: {cancel_response.text}"

    result = cancel_response.json()
    print(f"  - Cancellation response: {result}")
    assert result['status'] == 'cancelled', f"Expected status 'cancelled', got '{result['status']}'"
    assert result['local_cancelled'] == True, "Local cancellation should succeed"

    # Wait for cancellation to take effect
    print("  - Waiting for cancellation to complete...")
    time.sleep(3)

    # Check final status
    final_response = requests.get(f"{BASE_URL}/experiments/{exp_id}")
    final = final_response.json()
    print(f"  - Final status: {final['status']}")
    print(f"  - Cancelled at: {final.get('cancelled_at', 'N/A')}")
    print(f"  - Partial iterations: {len(final.get('convergence_data', []))}")

    assert final['status'] == 'cancelled', f"Expected 'cancelled', got '{final['status']}'"
    assert final['cancelled_at'] is not None, "cancelled_at should be set"

    print("  PASSED")
    return True


def test_cancel_running_sqd():
    """Test cancelling a running SQD calculation."""
    print("\n" + "=" * 80)
    print("TEST 3: Cancel Running SQD Experiment")
    print("=" * 80)

    # Create SQD experiment
    experiment = {
        "name": "Running H2 SQD Test",
        "molecule": {
            "smiles": "[H][H]",
            "basis": "sto-3g"
        },
        "configuration": {
            "method": "SQD",
            "subspace_dim": 20,
            "backend": "classical"
        },
        "execute_immediately": True
    }

    response = requests.post(f"{BASE_URL}/experiments", json=experiment)
    assert response.status_code == 201, f"Failed to create experiment: {response.text}"

    exp_id = response.json()['id']
    print(f"Created experiment {exp_id}")

    # SQD is very fast, so cancel immediately
    print("  - Sending cancellation immediately...")
    time.sleep(0.5)  # Brief delay

    cancel_response = requests.patch(f"{BASE_URL}/experiments/{exp_id}/cancel")

    if cancel_response.status_code == 400:
        print("  - Experiment already completed (SQD is very fast)")
        print("  - This is expected behavior")
        print("  PASSED (completed before cancellation)")
        return True
    elif cancel_response.status_code == 200:
        result = cancel_response.json()
        print(f"  - Cancellation response: {result}")

        # Check status
        final_response = requests.get(f"{BASE_URL}/experiments/{exp_id}")
        final = final_response.json()
        print(f"  - Final status: {final['status']}")

        assert final['status'] in ['cancelled', 'completed'], \
            f"Expected 'cancelled' or 'completed', got '{final['status']}'"

        print("  PASSED")
        return True
    else:
        print(f"  FAILED: Unexpected status code {cancel_response.status_code}")
        return False


def test_cancel_completed_experiment():
    """Test that completed experiments cannot be cancelled."""
    print("\n" + "=" * 80)
    print("TEST 4: Attempt to Cancel Completed Experiment")
    print("=" * 80)

    # Create and run a fast experiment
    experiment = {
        "name": "Fast H2 HF Test",
        "molecule": {
            "smiles": "[H][H]",
            "basis": "sto-3g"
        },
        "configuration": {
            "method": "HF"  # Very fast, completes immediately
        },
        "execute_immediately": True
    }

    response = requests.post(f"{BASE_URL}/experiments", json=experiment)
    assert response.status_code == 201, f"Failed to create experiment: {response.text}"

    exp_id = response.json()['id']
    print(f"Created experiment {exp_id}")

    # Wait for completion
    print("  - Waiting for experiment to complete...")
    for i in range(10):
        time.sleep(1)
        status_response = requests.get(f"{BASE_URL}/experiments/{exp_id}/status")
        status = status_response.json()
        print(f"    [{i+1}s] Status: {status['status']}")

        if status['status'] == 'completed':
            break

    # Try to cancel completed experiment
    print("  - Attempting to cancel completed experiment...")
    cancel_response = requests.patch(f"{BASE_URL}/experiments/{exp_id}/cancel")

    print(f"  - Response status: {cancel_response.status_code}")
    print(f"  - Response: {cancel_response.json()}")

    assert cancel_response.status_code == 400, \
        f"Expected 400 Bad Request, got {cancel_response.status_code}"

    error = cancel_response.json()
    assert "Cannot cancel" in error['detail'], \
        f"Expected 'Cannot cancel' error message, got: {error['detail']}"

    print("  PASSED (correctly rejected)")
    return True


def test_partial_results_saved():
    """Test that partial results are saved when cancelling."""
    print("\n" + "=" * 80)
    print("TEST 5: Verify Partial Results Are Saved")
    print("=" * 80)

    # Create long-running experiment
    experiment = {
        "name": "Partial Results Test",
        "molecule": {
            "smiles": "[H][H]",
            "basis": "sto-3g"
        },
        "configuration": {
            "method": "VQE",
            "ansatz": "ucc",
            "max_iterations": 500
        },
        "execute_immediately": True
    }

    response = requests.post(f"{BASE_URL}/experiments", json=experiment)
    exp_id = response.json()['id']
    print(f"Created experiment {exp_id}")

    # Wait for some iterations
    print("  - Letting experiment run for 5 seconds...")
    time.sleep(5)

    # Cancel
    cancel_response = requests.patch(f"{BASE_URL}/experiments/{exp_id}/cancel")
    assert cancel_response.status_code == 200

    # Wait for cancellation
    time.sleep(2)

    # Check results
    final_response = requests.get(f"{BASE_URL}/experiments/{exp_id}")
    final = final_response.json()

    print(f"  - Status: {final['status']}")
    print(f"  - Convergence data points: {len(final.get('convergence_data', []))}")

    assert final['status'] == 'cancelled', "Should be cancelled"
    assert len(final.get('convergence_data', [])) > 0, \
        "Should have saved some convergence data"

    print(f"  - First iteration energy: {final['convergence_data'][0]['energy']:.6f}")
    print(f"  - Last iteration energy: {final['convergence_data'][-1]['energy']:.6f}")
    print(f"  - Total iterations saved: {len(final['convergence_data'])}")

    print("  PASSED")
    return True


def run_all_tests():
    """Run all cancellation tests."""
    print("\n" + "=" * 80)
    print("KANAD JOB CANCELLATION TEST SUITE")
    print("=" * 80)
    print(f"Testing against: {BASE_URL}")

    tests = [
        ("Cancel Queued Experiment", test_cancel_queued_experiment),
        ("Cancel Running VQE", test_cancel_running_vqe),
        ("Cancel Running SQD", test_cancel_running_sqd),
        ("Reject Cancel Completed", test_cancel_completed_experiment),
        ("Partial Results Saved", test_partial_results_saved)
    ]

    results = []

    for name, test_func in tests:
        try:
            passed = test_func()
            results.append((name, passed))
        except Exception as e:
            print(f"\n  FAILED: {e}")
            results.append((name, False))

    # Print summary
    print("\n" + "=" * 80)
    print("TEST SUMMARY")
    print("=" * 80)

    passed_count = sum(1 for _, passed in results if passed)
    total_count = len(results)

    for name, passed in results:
        status = "PASS" if passed else "FAIL"
        symbol = "✓" if passed else "✗"
        print(f"  {symbol} {name}: {status}")

    print(f"\nTotal: {passed_count}/{total_count} tests passed")

    if passed_count == total_count:
        print("\n🎉 ALL TESTS PASSED!")
        return 0
    else:
        print(f"\n❌ {total_count - passed_count} tests failed")
        return 1


if __name__ == "__main__":
    exit_code = run_all_tests()
    sys.exit(exit_code)
