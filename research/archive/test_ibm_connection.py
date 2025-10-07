"""
Simple IBM Quantum Connection Test
===================================

Quick test to verify IBM Quantum connection and list available backends.
"""

import os
import sys
from pathlib import Path
from dotenv import load_dotenv

# Load credentials
env_path = Path(__file__).parent.parent / '.env'
load_dotenv(env_path)

print("=" * 80)
print("IBM QUANTUM CONNECTION TEST")
print("=" * 80)

# Check credentials
ibm_api = os.getenv('IBM_API')
if not ibm_api:
    print("✗ IBM_API not found in .env")
    sys.exit(1)

print(f"✓ IBM_API loaded: {ibm_api[:15]}...{ibm_api[-10:]}\n")

# Try to connect
print("Connecting to IBM Quantum Platform...")
print("(This may take 30-60 seconds on first connection...)\n")

try:
    from qiskit_ibm_runtime import QiskitRuntimeService

    # Initialize service
    service = QiskitRuntimeService(
        channel='ibm_quantum_platform',
        token=ibm_api
    )

    print("✓ Successfully connected to IBM Quantum!\n")

    # List available backends
    print("=" * 80)
    print("AVAILABLE QUANTUM BACKENDS")
    print("=" * 80)

    backends = service.backends()

    print(f"\nFound {len(backends)} backends:\n")

    for backend in backends:
        status = backend.status()
        config = backend.configuration()

        print(f"Backend: {backend.name}")
        print(f"  Qubits: {config.n_qubits}")
        print(f"  Status: {'✓ Operational' if status.operational else '✗ Down'}")
        print(f"  Queue: {status.pending_jobs} jobs")
        print(f"  Simulator: {config.simulator}")
        print()

    # Try to get a specific backend
    print("=" * 80)
    print("TESTING SPECIFIC BACKEND: ibm_brisbane")
    print("=" * 80)

    try:
        backend = service.backend('ibm_brisbane')
        print(f"\n✓ Successfully accessed: {backend.name}")
        print(f"  Qubits: {backend.num_qubits}")
        print(f"  Version: {backend.version}")

        status = backend.status()
        print(f"  Operational: {status.operational}")
        print(f"  Pending jobs: {status.pending_jobs}")

    except Exception as e:
        print(f"\n✗ Could not access ibm_brisbane: {e}")
        print("\nTry one of the available backends listed above")

    print("\n" + "=" * 80)
    print("✓✓✓ CONNECTION TEST SUCCESSFUL ✓✓✓")
    print("=" * 80)
    print("\nYou can now run: python research/ibm_large_molecules.py")

except Exception as e:
    print(f"✗ Connection failed: {e}")
    import traceback
    traceback.print_exc()

    print("\nTroubleshooting:")
    print("1. Check your IBM_API token is valid at https://quantum.ibm.com")
    print("2. Verify network connectivity")
    print("3. Try: pip install --upgrade qiskit-ibm-runtime")
    sys.exit(1)
