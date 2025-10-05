"""
Integration test for IBM Quantum Runtime backend.

Requires IBM Quantum credentials in .env file or environment variables.
"""

import pytest
import os
from pathlib import Path


def load_env_file():
    """Load environment variables from .env file if it exists."""
    env_file = Path(__file__).parent.parent.parent / '.env'

    if env_file.exists():
        with open(env_file) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    key, _, value = line.partition('=')
                    os.environ[key.strip()] = value.strip().strip('"').strip("'")
        print(f"Loaded environment from {env_file}")
    else:
        print(f"No .env file found at {env_file}")


# Load .env before importing backend
load_env_file()


from kanad.backends import IBMRuntimeBackend


@pytest.fixture
def ibm_credentials():
    """Get IBM Quantum credentials from environment."""
    token = (
        os.getenv('IBM_QUANTUM_TOKEN') or
        os.getenv('QISKIT_IBM_TOKEN') or
        os.getenv('API')  # Legacy/custom naming
    )
    instance = (
        os.getenv('IBM_QUANTUM_CRN') or
        os.getenv('QISKIT_IBM_INSTANCE') or
        os.getenv('CRN')  # Legacy/custom naming
    )
    channel = (
        os.getenv('IBM_QUANTUM_CHANNEL') or
        os.getenv('QISKIT_IBM_CHANNEL') or
        'ibm_quantum_platform'
    )

    if not token:
        pytest.skip("IBM Quantum token not found in environment")

    return {
        'token': token,
        'instance': instance,
        'channel': channel
    }


@pytest.fixture
def ibm_backend(ibm_credentials):
    """Create IBM Runtime backend instance."""
    backend = IBMRuntimeBackend(
        token=ibm_credentials['token'],
        instance=ibm_credentials['instance'],
        channel=ibm_credentials['channel'],
        shots=1024,
        optimization_level=3,
        resilience_level=1
    )
    yield backend
    backend.close()


def test_backend_initialization(ibm_backend):
    """Test that backend initializes correctly."""
    assert ibm_backend is not None
    assert ibm_backend.backend is not None
    assert ibm_backend._service is not None
    print(f"✓ Connected to backend: {ibm_backend.backend_name}")


def test_list_backends(ibm_backend):
    """Test listing available backends."""
    backends = ibm_backend.list_backends()
    assert len(backends) > 0
    print(f"✓ Found {len(backends)} available backends")
    print(f"  Examples: {backends[:5]}")


def test_backend_info(ibm_backend):
    """Test getting backend information."""
    info = ibm_backend.get_backend_info()

    assert 'name' in info
    assert 'channel' in info
    assert 'num_qubits' in info

    print(f"✓ Backend info:")
    print(f"  Name: {info['name']}")
    print(f"  Channel: {info['channel']}")
    print(f"  Qubits: {info.get('num_qubits', 'N/A')}")
    print(f"  Status: {info.get('status', 'N/A')}")


def test_get_estimator(ibm_backend):
    """Test creating an Estimator primitive."""
    estimator = ibm_backend.get_estimator()
    assert estimator is not None
    print(f"✓ Created EstimatorV2 primitive")


def test_get_sampler(ibm_backend):
    """Test creating a Sampler primitive."""
    sampler = ibm_backend.get_sampler()
    assert sampler is not None
    print(f"✓ Created SamplerV2 primitive")


def test_transpile_simple_circuit(ibm_backend):
    """Test circuit transpilation."""
    from qiskit import QuantumCircuit

    # Create simple circuit
    qc = QuantumCircuit(2)
    qc.h(0)
    qc.cx(0, 1)
    qc.measure_all()

    # Transpile for backend
    transpiled = ibm_backend.transpile_circuit(qc)

    assert transpiled is not None
    assert transpiled.num_qubits >= qc.num_qubits

    print(f"✓ Transpiled circuit:")
    print(f"  Original depth: {qc.depth()}")
    print(f"  Transpiled depth: {transpiled.depth()}")


def test_context_manager(ibm_credentials):
    """Test using backend as context manager."""
    with IBMRuntimeBackend(
        token=ibm_credentials['token'],
        instance=ibm_credentials['instance'],
        channel=ibm_credentials['channel']
    ) as backend:
        assert backend is not None
        info = backend.get_backend_info()
        assert info['name'] is not None
        print(f"✓ Context manager works correctly")


@pytest.mark.slow
def test_simple_vqe_energy_estimation(ibm_backend):
    """
    Test simple energy estimation using Estimator.

    This test runs a simple H2 VQE circuit (without full optimization).
    Marked as slow since it submits a job to IBM Quantum.
    """
    from qiskit import QuantumCircuit
    from qiskit.quantum_info import SparsePauliOp

    # Create simple H2 ansatz circuit (2 qubits)
    qc = QuantumCircuit(2)
    qc.ry(0.5, 0)
    qc.ry(0.5, 1)
    qc.cx(0, 1)

    # Simple H2 Hamiltonian (example)
    hamiltonian = SparsePauliOp.from_list([
        ("II", -1.0523),
        ("IZ", 0.3979),
        ("ZI", -0.3979),
        ("ZZ", -0.0112),
        ("XX", 0.1809)
    ])

    # Get estimator and run
    estimator = ibm_backend.get_estimator()

    # Note: EstimatorV2 has different API
    # Submit job (this will queue on IBM Quantum)
    job = estimator.run([(qc, hamiltonian)])

    print(f"✓ Submitted job to IBM Quantum")
    print(f"  Job ID: {job.job_id()}")
    print(f"  Status: {job.status()}")

    # Don't wait for result in test (can take a long time)
    # In real usage, you would call: result = job.result()


if __name__ == '__main__':
    """Run tests directly for quick validation."""
    print("\n" + "="*60)
    print("IBM Quantum Runtime Backend Integration Test")
    print("="*60 + "\n")

    # Load credentials
    load_env_file()

    token = (
        os.getenv('IBM_QUANTUM_TOKEN') or
        os.getenv('QISKIT_IBM_TOKEN') or
        os.getenv('API')
    )
    instance = (
        os.getenv('IBM_QUANTUM_CRN') or
        os.getenv('QISKIT_IBM_INSTANCE') or
        os.getenv('CRN')
    )
    channel = (
        os.getenv('IBM_QUANTUM_CHANNEL') or
        os.getenv('QISKIT_IBM_CHANNEL') or
        'ibm_quantum_platform'
    )

    if not token:
        print("❌ ERROR: IBM Quantum token not found!")
        print("\nPlease create a .env file in the project root with:")
        print("  IBM_QUANTUM_TOKEN=your_api_key_here")
        print("  IBM_QUANTUM_CRN=your_crn_here  # For IBM Cloud channel")
        print("  IBM_QUANTUM_CHANNEL=ibm_quantum  # or ibm_cloud")
        print("\nGet credentials from: https://quantum.ibm.com/account")
        exit(1)

    print(f"✓ Found credentials")
    print(f"  Channel: {channel}")
    print(f"  Token: {token[:20]}..." if token else "  Token: None")
    print(f"  CRN: {instance[:30]}..." if instance else "  CRN: None")
    print()

    # Run tests
    try:
        print("Test 1: Backend initialization...")
        with IBMRuntimeBackend(token=token, instance=instance, channel=channel) as backend:
            test_backend_initialization(backend)
            print()

            print("Test 2: List backends...")
            test_list_backends(backend)
            print()

            print("Test 3: Backend info...")
            test_backend_info(backend)
            print()

            print("Test 4: Get estimator...")
            test_get_estimator(backend)
            print()

            print("Test 5: Get sampler...")
            test_get_sampler(backend)
            print()

            print("Test 6: Transpile circuit...")
            test_transpile_simple_circuit(backend)
            print()

        print("\n" + "="*60)
        print("✓ All tests passed!")
        print("="*60)

    except Exception as e:
        print(f"\n❌ Test failed: {e}")
        import traceback
        traceback.print_exc()
        exit(1)
