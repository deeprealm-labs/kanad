#!/usr/bin/env python3
"""
Integration test for BlueQubit backend with VQE solver.
Tests the full stack: Molecule -> Hamiltonian -> VQE -> BlueQubit Backend
"""

import os
import pytest
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class TestBlueQubitIntegration:
    """Integration tests for BlueQubit backend."""

    @pytest.fixture(autouse=True)
    def setup(self):
        """Setup test environment."""
        # Check if credentials are available
        self.has_bluequbit_credentials = os.getenv('BLUE_TOKEN') is not None

        if not self.has_bluequbit_credentials:
            pytest.skip("BlueQubit credentials not available")

    def test_bluequbit_backend_initialization(self):
        """Test BlueQubit backend initializes correctly."""
        from kanad.backends.bluequbit import BlueQubitBackend

        backend = BlueQubitBackend(
            api_token=os.environ['BLUE_TOKEN']
        )

        assert backend.api_token is not None
        logger.info(f"✓ BlueQubit backend initialized")

    def test_vqe_bluequbit_integration(self):
        """Test VQE solver with BlueQubit backend (minimal test)."""
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

        # Create VQE solver with BlueQubit backend
        solver = VQESolver(
            hamiltonian=hamiltonian,
            ansatz_type='ucc',
            mapper_type='jordan_wigner',
            molecule=molecule,
            optimizer='SLSQP',
            max_iterations=1,  # Only 1 iteration for testing
            backend='bluequbit',
            shots=1024,
            enable_analysis=False,
            # BlueQubit credentials
            api_token=os.environ['BLUE_TOKEN']
        )

        # Verify backend is properly configured
        assert solver.backend == 'bluequbit'
        assert solver._use_statevector is False
        assert hasattr(solver, '_bluequbit_backend')
        assert solver._bluequbit_backend is not None

        logger.info(f"✓ VQE solver configured with BlueQubit backend")
        logger.info(f"  Backend: {solver._bluequbit_backend}")
        logger.info(f"  Use statevector: {solver._use_statevector}")

        # NOTE: We don't actually run solve() here to avoid BlueQubit credits usage
        # The full solve would consume credits
        # This test just verifies the setup is correct


if __name__ == '__main__':
    # Run tests
    pytest.main([__file__, '-v', '-s'])
