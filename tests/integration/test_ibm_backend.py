"""
Integration tests for IBM Quantum backend.

These tests require valid IBM Quantum credentials.
Set environment variables:
  export IBM_API="your_api_token"
  export IBM_CRN="your_cloud_resource_name"

Run with:
  python -m pytest tests/integration/test_ibm_backend.py -v -s
"""

import os
import pytest
import numpy as np
from kanad.core.molecule import Molecule
from kanad.core.atom import Atom
from kanad.core.hamiltonians.molecular_hamiltonian import MolecularHamiltonian
from kanad.solvers.vqe_solver import VQESolver
from kanad.backends.ibm.backend import IBMBackend


# Skip all tests if IBM credentials not available
pytestmark = pytest.mark.skipif(
    not os.getenv('IBM_API'),
    reason="IBM_API environment variable not set"
)


class TestIBMBackendBasic:
    """Test IBM backend initialization and basic operations."""

    def test_ibm_backend_init(self):
        """Test IBM backend initialization."""
        api_token = os.getenv('IBM_API')
        crn = os.getenv('IBM_CRN')

        backend = IBMBackend(
            backend_name='ibm_brisbane',  # Use a specific backend
            api_token=api_token,
            instance=crn
        )

        assert backend is not None
        assert backend.backend_name == 'ibm_brisbane'
        assert backend.service is not None
        assert backend.backend is not None

        print(f"\n✓ Backend initialized: {backend.backend.name}")
        print(f"  Qubits: {backend.backend.num_qubits}")

    def test_ibm_backend_info(self):
        """Test getting backend information."""
        api_token = os.getenv('IBM_API')
        crn = os.getenv('IBM_CRN')

        backend = IBMBackend(
            backend_name='ibm_brisbane',
            api_token=api_token,
            instance=crn
        )

        info = backend.get_backend_info()

        assert 'name' in info
        assert 'num_qubits' in info
        assert 'is_simulator' in info
        assert 'is_operational' in info

        print(f"\n✓ Backend info:")
        print(f"  Name: {info['name']}")
        print(f"  Qubits: {info['num_qubits']}")
        print(f"  Simulator: {info['is_simulator']}")
        print(f"  Operational: {info['is_operational']}")
        print(f"  Pending jobs: {info['pending_jobs']}")


class TestIBMBackendCircuitExecution:
    """Test circuit execution on IBM backend."""

    def test_simple_circuit_submission(self):
        """Test submitting a simple circuit to IBM."""
        from qiskit import QuantumCircuit
        from qiskit.quantum_info import SparsePauliOp

        api_token = os.getenv('IBM_API')
        crn = os.getenv('IBM_CRN')

        backend = IBMBackend(
            backend_name='ibm_brisbane',
            api_token=api_token,
            instance=crn
        )

        # Create simple 2-qubit circuit
        circuit = QuantumCircuit(2)
        circuit.h(0)
        circuit.cx(0, 1)

        # Create simple observable (Z₀)
        observable = SparsePauliOp.from_list([("ZI", 1.0)])

        # Submit job
        print("\n✓ Submitting circuit to IBM Quantum...")
        result = backend.run_batch(
            circuits=[circuit],
            observables=[observable],
            shots=1024
        )

        assert 'job_id' in result
        assert 'status' in result
        assert 'backend' in result

        job_id = result['job_id']
        print(f"  Job ID: {job_id}")
        print(f"  Status: {result['status']}")
        print(f"  Backend: {result['backend']}")

        # Check job status
        status = backend.get_job_status(job_id)
        print(f"  Current status: {status}")

        # Note: We don't wait for completion in this test
        # as IBM jobs can take several minutes


class TestIBMBackendWithVQE:
    """Test VQE solver with IBM backend."""

    def test_h2_vqe_with_ibm(self):
        """
        Test H2 molecule VQE with IBM backend.

        WARNING: This test submits real jobs to IBM Quantum!
        It will use IBM compute credits and may take 10-30 minutes.
        """
        api_token = os.getenv('IBM_API')
        crn = os.getenv('IBM_CRN')

        # Create IBM backend
        ibm_backend = IBMBackend(
            backend_name='ibm_brisbane',
            api_token=api_token,
            instance=crn
        )

        # Create H2 molecule
        h1 = Atom('H', position=(0.0, 0.0, 0.0))
        h2 = Atom('H', position=(0.74, 0.0, 0.0))  # Bond length ~0.74 Angstrom

        molecule = Molecule(
            atoms=[h1, h2],
            charge=0,
            spin=0,
            basis='sto-3g'
        )

        # Create Hamiltonian
        hamiltonian = MolecularHamiltonian(molecule)

        # Create VQE solver with IBM backend
        print("\n✓ Creating VQE solver with IBM backend...")
        solver = VQESolver(
            hamiltonian=hamiltonian,
            ansatz_type='ucc',
            mapper_type='jordan_wigner',
            molecule=molecule,
            optimizer='SLSQP',
            max_iterations=3,  # Use very few iterations for testing!
            backend=ibm_backend,
            shots=1024,
            enable_analysis=True
        )

        print(f"  Solver initialized: {solver.ansatz_type} ansatz")
        print(f"  Backend: {solver.backend}")
        print(f"  Max iterations: {solver.max_iterations}")

        # Run VQE (WARNING: This will submit jobs to IBM!)
        print("\n⚠️  Starting VQE optimization on IBM Quantum...")
        print("   This will submit jobs to IBM and may take 10-30 minutes!")
        print("   Press Ctrl+C to cancel if you don't want to use IBM credits.")

        try:
            result = solver.solve()

            print(f"\n✓ VQE completed!")
            print(f"  Energy: {result['energy']:.8f} Hartree")
            print(f"  HF Energy: {result['hf_energy']:.8f} Hartree")
            print(f"  Correlation: {result['correlation_energy']:.8f} Hartree")
            print(f"  Iterations: {result['iterations']}")
            print(f"  Converged: {result['converged']}")

            # Check that energy is reasonable for H2
            assert result['energy'] < 0  # Should be negative
            assert result['energy'] > -2.0  # Shouldn't be too negative
            assert result['iterations'] <= 3  # Should respect max_iterations

        except KeyboardInterrupt:
            print("\n⚠️  Test cancelled by user")
            pytest.skip("Test cancelled - user interrupted")


class TestIBMBackendJobManagement:
    """Test IBM job management features."""

    def test_list_backends(self):
        """Test listing available IBM backends."""
        api_token = os.getenv('IBM_API')
        crn = os.getenv('IBM_CRN')

        backend = IBMBackend(
            backend_name='ibm_brisbane',
            api_token=api_token,
            instance=crn
        )

        # List real hardware backends
        hardware_backends = backend.list_backends(simulator=False, operational=True)
        print(f"\n✓ Available hardware backends: {len(hardware_backends)}")
        for name in hardware_backends[:5]:  # Show first 5
            print(f"  - {name}")

        assert len(hardware_backends) > 0

        # List simulators
        simulator_backends = backend.list_backends(simulator=True, operational=True)
        print(f"\n✓ Available simulators: {len(simulator_backends)}")
        for name in simulator_backends[:3]:
            print(f"  - {name}")


if __name__ == '__main__':
    # Run tests with verbose output
    pytest.main([__file__, '-v', '-s'])
