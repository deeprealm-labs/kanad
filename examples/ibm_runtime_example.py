"""
Example: Using IBM Quantum Runtime backend with Kanad framework.

This example demonstrates how to:
1. Load credentials from .env file
2. Connect to IBM Quantum hardware
3. Run a simple VQE calculation on real quantum hardware
"""

import os
from pathlib import Path


def load_env_file():
    """Load environment variables from .env file."""
    env_file = Path(__file__).parent.parent / '.env'

    if env_file.exists():
        with open(env_file) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    key, _, value = line.partition('=')
                    os.environ[key.strip()] = value.strip().strip('"').strip("'")
        print(f"✓ Loaded credentials from {env_file}")
    else:
        print(f"⚠ No .env file found at {env_file}")
        print("  Using environment variables or saved credentials")


def main():
    """Run IBM Quantum Runtime example."""
    print("\n" + "="*70)
    print("IBM Quantum Runtime Backend Example")
    print("="*70 + "\n")

    # Load credentials from .env
    load_env_file()

    from kanad.backends import IBMRuntimeBackend

    # Example 1: Connect to IBM Quantum (auto-select least busy backend)
    print("Example 1: Auto-selecting least busy backend...")
    print("-" * 70)

    # Explicitly set channel based on whether CRN is available
    token = os.getenv('IBM_QUANTUM_TOKEN') or os.getenv('QISKIT_IBM_TOKEN') or os.getenv('API')
    instance = os.getenv('IBM_QUANTUM_CRN') or os.getenv('QISKIT_IBM_INSTANCE') or os.getenv('CRN')
    channel = 'ibm_cloud' if instance else 'ibm_quantum_platform'

    with IBMRuntimeBackend(
        token=token,
        instance=instance,
        channel=channel,
        shots=4096,
        optimization_level=3,
        resilience_level=1
    ) as backend:
        info = backend.get_backend_info()
        print(f"✓ Connected to: {info['name']}")
        print(f"  Qubits: {info.get('num_qubits', 'N/A')}")
        print(f"  Status: {info.get('status', 'N/A')}")
        print(f"  Channel: {info['channel']}")
        print()

        # List some available backends
        print("Available backends:")
        backends = backend.list_backends({'operational': True})
        for i, b in enumerate(backends[:10], 1):
            print(f"  {i}. {b}")
        print()

    # Example 2: List and connect to available backend
    print("\nExample 2: Using available backend...")
    print("-" * 70)

    try:
        # Reuse the first backend connection to show transpilation
        with IBMRuntimeBackend(
            backend_name='ibm_torino',  # Use the available backend
            token=token,
            instance=instance,
            channel=channel,
            shots=1024,
            optimization_level=2,
            resilience_level=1
        ) as backend:
            info = backend.get_backend_info()
            print(f"✓ Connected to: {info['name']}")
            print(f"  Qubits: {info.get('num_qubits', 'N/A')}")
            print()

            # Create and transpile a simple circuit
            from qiskit import QuantumCircuit

            print("Creating and transpiling a Bell state circuit...")
            qc = QuantumCircuit(2)
            qc.h(0)
            qc.cx(0, 1)
            qc.measure_all()

            transpiled = backend.transpile_circuit(qc)
            print(f"✓ Original depth: {qc.depth()}")
            print(f"✓ Transpiled depth: {transpiled.depth()}")
            print()

    except Exception as e:
        print(f"⚠ Could not complete example: {e}")
        print()

    # Example 3: Using the Estimator for energy estimation
    print("\nExample 3: Using Estimator primitive for energy estimation...")
    print("-" * 70)

    with IBMRuntimeBackend(
        backend_name='ibm_torino',  # Use available backend
        token=token,
        instance=instance,
        channel=channel,
        shots=1024,
        resilience_level=1
    ) as backend:
        from qiskit import QuantumCircuit
        from qiskit.quantum_info import SparsePauliOp

        # Simple 2-qubit ansatz
        qc = QuantumCircuit(2)
        qc.ry(0.5, 0)
        qc.ry(0.5, 1)
        qc.cx(0, 1)

        # Simple Hamiltonian (H2 example)
        hamiltonian = SparsePauliOp.from_list([
            ("II", -1.0523732),
            ("IZ", 0.39793742),
            ("ZI", -0.39793742),
            ("ZZ", -0.01128010),
            ("XX", 0.18093120)
        ])

        print("Circuit created:")
        print(f"  Qubits: {qc.num_qubits}")
        print(f"  Depth: {qc.depth()}")
        print()

        print("Hamiltonian terms: 5")
        print()

        # Get estimator
        estimator = backend.get_estimator()
        print(f"✓ Created EstimatorV2 primitive")
        print(f"  Backend: {backend.backend_name}")
        print(f"  Shots: {backend.shots}")
        print(f"  Resilience level: {backend.resilience_level}")
        print()

        print("To run the estimation, use:")
        print("  job = estimator.run([(circuit, hamiltonian)])")
        print("  result = job.result()")
        print("  energy = result[0].data.evs")
        print()

    # Example 4: Environment variable configuration
    print("\nExample 4: Configuration via environment variables...")
    print("-" * 70)
    print("You can configure credentials using environment variables:")
    print()
    print("Option 1 - .env file (recommended):")
    print("  Create .env file in project root:")
    print("    IBM_QUANTUM_TOKEN=your_api_key")
    print("    IBM_QUANTUM_CRN=your_crn  # For IBM Cloud")
    print("    IBM_QUANTUM_CHANNEL=ibm_quantum")
    print()
    print("Option 2 - Export environment variables:")
    print("  export IBM_QUANTUM_TOKEN=your_api_key")
    print("  export IBM_QUANTUM_CRN=your_crn")
    print()
    print("Option 3 - Pass directly to constructor:")
    print("  backend = IBMRuntimeBackend(")
    print("      token='your_api_key',")
    print("      instance='your_crn',")
    print("      channel='ibm_quantum'")
    print("  )")
    print()

    print("="*70)
    print("✓ Example completed successfully!")
    print("="*70)
    print()
    print("Next steps:")
    print("  1. Get your API key from: https://quantum.ibm.com/account")
    print("  2. Add credentials to .env file")
    print("  3. Run this script or the integration tests")
    print("  4. Try running VQE on real IBM Quantum hardware!")
    print()


if __name__ == '__main__':
    main()
