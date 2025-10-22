#!/usr/bin/env python3
"""
Integration test for IBM Quantum backend with VQE solver.
Tests the full stack: Molecule -> Hamiltonian -> VQE -> IBM Backend
"""

import os
import pytest
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class TestIBMIntegration:
    """Integration tests for IBM Quantum backend."""

    @pytest.fixture(autouse=True)
    def setup(self):
        """Setup test environment."""
        # Check if credentials are available
        self.has_ibm_credentials = (
            os.getenv('IBM_API') is not None and
            os.getenv('IBM_CRN') is not None
        )

        if not self.has_ibm_credentials:
            pytest.skip("IBM credentials not available")

    def test_ibm_backend_initialization(self):
        """Test IBM backend initializes correctly."""
        from kanad.backends.ibm import IBMBackend

        backend = IBMBackend(
            backend_name='ibm_brisbane',
            api_token=os.environ['IBM_API'],
            instance=os.environ['IBM_CRN']
        )

        assert backend.backend.num_qubits == 127
        assert backend.backend.name == 'ibm_brisbane'
        logger.info(f"✓ IBM backend initialized: {backend.backend.name}")

    def test_ibm_backend_info(self):
        """Test getting backend information."""
        from kanad.backends.ibm import IBMBackend

        backend = IBMBackend(
            backend_name='ibm_brisbane',
            api_token=os.environ['IBM_API'],
            instance=os.environ['IBM_CRN']
        )

        info = backend.get_backend_info()

        assert info['name'] == 'ibm_brisbane'
        assert info['num_qubits'] == 127
        assert isinstance(info['is_operational'], bool)
        assert isinstance(info['pending_jobs'], int)
        logger.info(f"✓ Backend operational: {info['is_operational']}, pending jobs: {info['pending_jobs']}")

    def test_ibm_circuit_submission(self):
        """Test submitting a simple circuit to IBM."""
        from kanad.backends.ibm import IBMBackend
        from qiskit import QuantumCircuit
        from qiskit.quantum_info import SparsePauliOp

        backend = IBMBackend(
            backend_name='ibm_brisbane',
            api_token=os.environ['IBM_API'],
            instance=os.environ['IBM_CRN']
        )

        # Create Bell state circuit
        circuit = QuantumCircuit(2)
        circuit.h(0)
        circuit.cx(0, 1)

        # Observable
        observable = SparsePauliOp.from_list([("ZI", 1.0)])

        # Submit job
        result = backend.run_batch(
            circuits=[circuit],
            observables=[observable],
            shots=1024
        )

        assert 'job_id' in result
        assert 'status' in result
        assert 'backend' in result
        assert result['backend'] == 'ibm_brisbane'

        logger.info(f"✓ Circuit submitted, Job ID: {result['job_id']}")

        # Test job status retrieval
        job_id = result['job_id']
        status = backend.get_job_status(job_id)
        assert isinstance(status, str)
        assert status in ['QUEUED', 'RUNNING', 'DONE', 'ERROR', 'CANCELLED']
        logger.info(f"✓ Job status: {status}")

    def test_vqe_ibm_integration(self):
        """Test VQE solver with IBM backend (minimal test)."""
        from kanad.core.molecule import Molecule
        from kanad.core.atom import Atom
        from kanad.solvers.vqe_solver import VQESolver

        # Create H2 molecule
        h1 = Atom('H', position=(0.0, 0.0, 0.0))
        h2 = Atom('H', position=(0.74, 0.0, 0.0))

        molecule = Molecule(
            atoms=[h1, h2],
            charge=0,
            spin=0,
            basis='sto-3g'
        )

        # Get Hamiltonian from molecule
        hamiltonian = molecule.hamiltonian

        # Create VQE solver with IBM backend
        solver = VQESolver(
            hamiltonian=hamiltonian,
            ansatz_type='ucc',
            mapper_type='jordan_wigner',
            molecule=molecule,
            optimizer='SLSQP',
            max_iterations=1,  # Only 1 iteration for testing
            backend='ibm',
            shots=1024,
            enable_analysis=False,
            # IBM credentials
            backend_name='ibm_brisbane',
            api_token=os.environ['IBM_API'],
            instance=os.environ['IBM_CRN']
        )

        # Verify backend is properly configured
        assert solver.backend == 'ibm'
        assert solver._use_statevector is False
        assert hasattr(solver, '_ibm_backend')
        assert solver._ibm_backend is not None

        logger.info(f"✓ VQE solver configured with IBM backend")
        logger.info(f"  Backend: {solver._ibm_backend}")
        logger.info(f"  Use statevector: {solver._use_statevector}")

        # NOTE: We don't actually run solve() here to avoid IBM usage
        # The full solve would submit multiple jobs to IBM
        # This test just verifies the setup is correct


if __name__ == '__main__':
    # Run tests
    pytest.main([__file__, '-v', '-s'])
