"""
Integration tests for BlueQubit backend.

These tests require a valid BlueQubit API token.
Set environment variable:
  export BLUE_TOKEN="your_api_token"

Run with:
  python -m pytest tests/integration/test_bluequbit_backend.py -v -s
"""

import os
import pytest
import numpy as np
from kanad.core.molecule import Molecule
from kanad.core.atom import Atom
from kanad.core.hamiltonians.molecular_hamiltonian import MolecularHamiltonian
from kanad.solvers.vqe_solver import VQESolver
from kanad.backends.bluequbit.backend import BlueQubitBackend


# Skip all tests if BlueQubit token not available
pytestmark = pytest.mark.skipif(
    not os.getenv('BLUE_TOKEN'),
    reason="BLUE_TOKEN environment variable not set"
)


class TestBlueQubitBackendBasic:
    """Test BlueQubit backend initialization and basic operations."""

    def test_bluequbit_backend_init(self):
        """Test BlueQubit backend initialization."""
        api_token = os.getenv('BLUE_TOKEN')

        backend = BlueQubitBackend(
            device='gpu',
            api_token=api_token
        )

        assert backend is not None
        assert backend.device == 'gpu'
        assert backend.bq is not None

        print(f"\n✓ BlueQubit backend initialized")
        print(f"  Device: {backend.device}")

    def test_bluequbit_device_info(self):
        """Test getting device information."""
        api_token = os.getenv('BLUE_TOKEN')

        backend = BlueQubitBackend(
            device='gpu',
            api_token=api_token
        )

        info = backend.get_device_info()

        assert 'device' in info
        assert 'max_qubits' in info
        assert 'supports_statevector' in info

        print(f"\n✓ Device info:")
        print(f"  Device: {info['device']}")
        print(f"  Max qubits: {info['max_qubits']}")
        print(f"  GPU accelerated: {info['is_gpu_accelerated']}")
        print(f"  Supports statevector: {info['supports_statevector']}")


class TestBlueQubitCircuitExecution:
    """Test circuit execution on BlueQubit."""

    def test_simple_circuit_execution(self):
        """Test running a simple circuit on BlueQubit."""
        from qiskit import QuantumCircuit

        api_token = os.getenv('BLUE_TOKEN')

        backend = BlueQubitBackend(
            device='gpu',
            api_token=api_token
        )

        # Create simple 2-qubit circuit
        circuit = QuantumCircuit(2)
        circuit.h(0)
        circuit.cx(0, 1)

        # Run circuit (statevector mode)
        print("\n✓ Running circuit on BlueQubit...")
        result = backend.run_circuit(
            circuit=circuit,
            shots=None,  # Statevector mode
            asynchronous=False
        )

        assert 'statevector' in result
        statevector = np.array(result['statevector'])

        print(f"  Statevector shape: {statevector.shape}")
        print(f"  Statevector norm: {np.linalg.norm(statevector):.6f}")

        # Check statevector is normalized
        assert np.abs(np.linalg.norm(statevector) - 1.0) < 1e-6

        # Check Bell state: |00⟩ + |11⟩ (with normalization)
        expected = np.zeros(4, dtype=complex)
        expected[0] = 1 / np.sqrt(2)  # |00⟩
        expected[3] = 1 / np.sqrt(2)  # |11⟩

        # Allow some numerical error
        assert np.allclose(statevector, expected, atol=1e-6)

        print(f"  ✓ Bell state verified!")


class TestBlueQubitWithVQE:
    """Test VQE solver with BlueQubit backend."""

    def test_h2_vqe_with_bluequbit(self):
        """
        Test H2 molecule VQE with BlueQubit backend.

        This test uses BlueQubit GPU simulator for fast execution.
        """
        api_token = os.getenv('BLUE_TOKEN')

        # Create BlueQubit backend
        bq_backend = BlueQubitBackend(
            device='gpu',
            api_token=api_token
        )

        # Create H2 molecule
        h1 = Atom('H', position=(0.0, 0.0, 0.0))
        h2 = Atom('H', position=(0.74, 0.0, 0.0))

        molecule = Molecule(
            atoms=[h1, h2],
            charge=0,
            spin=0,
            basis='sto-3g'
        )

        # Create Hamiltonian
        hamiltonian = MolecularHamiltonian(molecule)

        # Create VQE solver with BlueQubit backend
        print("\n✓ Creating VQE solver with BlueQubit backend...")
        solver = VQESolver(
            hamiltonian=hamiltonian,
            ansatz_type='ucc',
            mapper_type='jordan_wigner',
            molecule=molecule,
            optimizer='SLSQP',
            max_iterations=10,
            backend=bq_backend,
            shots=None,  # Statevector mode for faster execution
            enable_analysis=True
        )

        print(f"  Solver initialized: {solver.ansatz_type} ansatz")
        print(f"  Backend: {solver.backend}")
        print(f"  Max iterations: {solver.max_iterations}")

        # Run VQE
        print("\n✓ Starting VQE optimization on BlueQubit GPU...")
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
        assert result['energy'] < result['hf_energy']  # VQE should improve on HF

        # Known H2 exact energy at this geometry: ~-1.137 Ha
        # VQE with UCC should get close
        assert abs(result['energy'] - (-1.137)) < 0.1  # Within 0.1 Ha


class TestBlueQubitPerformance:
    """Test BlueQubit performance characteristics."""

    def test_circuit_depth_scaling(self):
        """Test performance with different circuit depths."""
        from qiskit import QuantumCircuit
        import time

        api_token = os.getenv('BLUE_TOKEN')

        backend = BlueQubitBackend(
            device='gpu',
            api_token=api_token
        )

        print(f"\n✓ Testing circuit depth scaling...")

        for n_qubits in [2, 4, 6]:
            circuit = QuantumCircuit(n_qubits)

            # Create layered circuit
            for layer in range(5):
                for i in range(n_qubits):
                    circuit.h(i)
                for i in range(n_qubits - 1):
                    circuit.cx(i, i + 1)

            start_time = time.time()
            result = backend.run_circuit(circuit, shots=None, asynchronous=False)
            elapsed = time.time() - start_time

            print(f"  {n_qubits} qubits, depth {circuit.depth()}: {elapsed:.3f}s")

            assert 'statevector' in result


if __name__ == '__main__':
    # Run tests with verbose output
    pytest.main([__file__, '-v', '-s'])
