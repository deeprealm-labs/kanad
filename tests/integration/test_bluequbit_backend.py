#!/usr/bin/env python3
"""
Integration Tests for BlueQubit Cloud Backend

Tests the BlueQubit cloud quantum computing integration with real API calls.
Requires valid BLUEQUBIT_API_TOKEN in environment or .env file.
"""

import pytest
import numpy as np
import sys
import os

# Add project root to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from kanad.backends.bluequbit_backend import (
    BlueQubitBackend,
    BlueQubitVQESolver,
    get_bluequbit_backend
)
from kanad.core.atom import Atom
from kanad.bonds.covalent_bond import CovalentBond


class TestBlueQubitBackend:
    """Test BlueQubit backend initialization and configuration."""

    def test_backend_initialization(self):
        """Test backend can be initialized with token from environment."""
        try:
            backend = BlueQubitBackend(device='cpu')
            assert backend.device == 'cpu'
            assert backend.execution_mode == 'cloud'
            assert backend.shots == 1024
        except ValueError as e:
            pytest.skip(f"No BlueQubit token available: {e}")

    def test_get_device_info(self):
        """Test device information retrieval."""
        try:
            backend = BlueQubitBackend(device='cpu')
            info = backend.get_device_info()

            assert 'current_device' in info
            assert 'available_devices' in info
            assert 'device_info' in info

            assert info['current_device'] == 'cpu'
            assert info['device_info']['qubits'] == 34
            assert info['device_info']['cost'] == 'free'

        except ValueError as e:
            pytest.skip(f"No BlueQubit token available: {e}")

    def test_convenience_function(self):
        """Test get_bluequbit_backend convenience function."""
        try:
            backend = get_bluequbit_backend(device='cpu', shots=2048)
            assert backend.device == 'cpu'
            assert backend.shots == 2048
        except ValueError as e:
            pytest.skip(f"No BlueQubit token available: {e}")


class TestBlueQubitCircuitExecution:
    """Test quantum circuit execution on BlueQubit cloud."""

    @pytest.mark.slow
    def test_simple_circuit(self):
        """Test running a simple quantum circuit on BlueQubit."""
        try:
            from qiskit import QuantumCircuit

            backend = BlueQubitBackend(device='cpu')

            # Create simple Bell state
            qc = QuantumCircuit(2)
            qc.h(0)
            qc.cx(0, 1)
            qc.measure_all()

            # Run on BlueQubit
            result = backend.run_circuit(qc, shots=1024)

            # Check results
            assert 'counts' in result
            assert 'job_id' in result
            assert 'device' in result
            assert result['device'] == 'cpu'
            assert result['shots'] == 1024
            assert result['success'] is True

            # Check measurement outcomes
            counts = result['counts']
            assert len(counts) > 0

        except ValueError as e:
            pytest.skip(f"No BlueQubit token available: {e}")
        except Exception as e:
            pytest.fail(f"Circuit execution failed: {e}")


class TestBlueQubitMolecules:
    """Test molecular calculations on BlueQubit."""

    def test_h2_molecule_setup(self):
        """Test H2 molecule can be prepared for BlueQubit."""
        atoms = [
            Atom('H', position=np.array([0.0, 0.0, 0.0])),
            Atom('H', position=np.array([0.74, 0.0, 0.0]))
        ]

        bond = CovalentBond(atoms[0], atoms[1])
        hamiltonian = bond.hamiltonian

        assert hamiltonian.n_orbitals == 2
        assert hamiltonian.n_electrons == 2

        # Should need 4 qubits (2 orbitals * 2)
        n_qubits = 2 * hamiltonian.n_orbitals
        assert n_qubits == 4

        # This fits on BlueQubit CPU (34 qubits)
        assert n_qubits <= 34

    @pytest.mark.slow
    def test_h2_circuit_on_cloud(self):
        """Test H2 molecular circuit execution on BlueQubit."""
        try:
            from qiskit import QuantumCircuit

            backend = BlueQubitBackend(device='cpu')

            # Create simple H2 circuit
            atoms = [
                Atom('H', position=np.array([0.0, 0.0, 0.0])),
                Atom('H', position=np.array([0.74, 0.0, 0.0]))
            ]
            bond = CovalentBond(atoms[0], atoms[1])
            hamiltonian = bond.hamiltonian

            n_qubits = 2 * hamiltonian.n_orbitals

            # Create circuit with HF initial state
            qc = QuantumCircuit(n_qubits)
            qc.x(0)  # First electron
            qc.x(1)  # Second electron
            qc.measure_all()

            # Run on BlueQubit
            result = backend.run_circuit(qc, shots=1024)

            assert result['success'] is True
            assert 'counts' in result

        except ValueError as e:
            pytest.skip(f"No BlueQubit token available: {e}")


class TestBlueQubitDevices:
    """Test different BlueQubit devices."""

    def test_cpu_device(self):
        """Test CPU device (free tier)."""
        try:
            backend = BlueQubitBackend(device='cpu')
            info = backend.get_device_info()

            assert info['device_info']['cost'] == 'free'
            assert info['device_info']['qubits'] == 34

        except ValueError as e:
            pytest.skip(f"No BlueQubit token available: {e}")

    def test_invalid_device(self):
        """Test handling of invalid device names."""
        try:
            backend = BlueQubitBackend(device='invalid_device')
            # Should still initialize, but execution may fail
            assert backend.device == 'invalid_device'

        except ValueError as e:
            pytest.skip(f"No BlueQubit token available: {e}")


class TestBlueQubitErrorHandling:
    """Test error handling and edge cases."""

    def test_missing_token(self):
        """Test error when token is missing."""
        # Temporarily clear token
        import os
        old_token = os.environ.get('BLUEQUBIT_API_TOKEN')
        old_token2 = os.environ.get('TOKEN')

        try:
            if old_token:
                del os.environ['BLUEQUBIT_API_TOKEN']
            if old_token2:
                del os.environ['TOKEN']

            with pytest.raises(ValueError, match="BlueQubit API token required"):
                BlueQubitBackend(api_token=None)

        finally:
            # Restore tokens
            if old_token:
                os.environ['BLUEQUBIT_API_TOKEN'] = old_token
            if old_token2:
                os.environ['TOKEN'] = old_token2

    def test_sdk_not_installed(self, monkeypatch):
        """Test error when BlueQubit SDK is not installed."""
        # Skip this test - mock import causes recursion issues in Python 3.13
        pytest.skip("SDK import mocking causes recursion issues in Python 3.13")


# Mark tests that require valid API token
pytestmark = pytest.mark.integration


if __name__ == '__main__':
    # Run tests
    pytest.main([__file__, '-v', '-s'])
