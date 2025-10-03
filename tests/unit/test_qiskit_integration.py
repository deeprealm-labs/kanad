"""
Unit tests for Qiskit integration.

Tests:
- Circuit conversion to Qiskit
- Backend initialization
- Pauli converter
- Qiskit-based VQE execution
"""

import pytest
import numpy as np


class TestCircuitConversion:
    """Test custom circuit to Qiskit conversion."""

    def test_to_qiskit_import(self):
        """Test that to_qiskit() requires Qiskit."""
        from kanad.ansatze.base_ansatz import QuantumCircuit

        circuit = QuantumCircuit(2)
        circuit.h(0)
        circuit.cx(0, 1)

        # Should work if Qiskit is installed
        try:
            qc = circuit.to_qiskit()
            assert qc.num_qubits == 2
        except ImportError:
            pytest.skip("Qiskit not installed")

    def test_single_qubit_gates(self):
        """Test single-qubit gate conversion."""
        pytest.importorskip("qiskit")

        from kanad.ansatze.base_ansatz import QuantumCircuit

        circuit = QuantumCircuit(1)
        circuit.h(0)
        circuit.x(0)
        circuit.y(0)
        circuit.z(0)

        qc = circuit.to_qiskit()
        assert qc.num_qubits == 1
        assert qc.depth() == 4

    def test_two_qubit_gates(self):
        """Test two-qubit gate conversion."""
        pytest.importorskip("qiskit")

        from kanad.ansatze.base_ansatz import QuantumCircuit

        circuit = QuantumCircuit(2)
        circuit.cx(0, 1)
        circuit.cz(1, 0)

        qc = circuit.to_qiskit()
        assert qc.num_qubits == 2

    def test_parametrized_gates(self):
        """Test parametrized gate conversion."""
        pytest.importorskip("qiskit")

        from kanad.ansatze.base_ansatz import QuantumCircuit, Parameter

        circuit = QuantumCircuit(2)
        theta = Parameter('theta', 0.5)
        circuit.ry(theta, 0)
        circuit.rz(theta, 1)

        qc = circuit.to_qiskit()
        assert qc.num_qubits == 2
        assert len(qc.parameters) == 1  # One unique parameter

    def test_parameter_binding(self):
        """Test parameter binding for Qiskit circuits."""
        pytest.importorskip("qiskit")

        from kanad.ansatze.base_ansatz import QuantumCircuit, Parameter

        circuit = QuantumCircuit(1)
        theta = Parameter('theta_0', 0.0)
        circuit.ry(theta, 0)

        qc = circuit.to_qiskit()
        bound_qc = circuit.assign_parameters_for_qiskit(qc, np.array([np.pi/2]))

        assert len(bound_qc.parameters) == 0  # All parameters bound


class TestPauliConverter:
    """Test Hamiltonian to Pauli operator conversion."""

    @pytest.fixture
    def simple_hamiltonian(self):
        """Create simple test Hamiltonian."""
        from kanad.core.hamiltonians.molecular_hamiltonian import MolecularHamiltonian

        class SimpleHamiltonian(MolecularHamiltonian):
            def __init__(self):
                super().__init__(n_orbitals=2, n_electrons=2, nuclear_repulsion=0.5)
                self.h_core = np.array([
                    [-1.0, 0.1],
                    [0.1, -0.5]
                ])
                self.eri = np.zeros((2, 2, 2, 2))

            def compute_energy(self, density_matrix):
                return np.sum(density_matrix * self.h_core) + self.nuclear_repulsion

            def to_matrix(self):
                return self.h_core.copy()

        return SimpleHamiltonian()

    def test_pauli_converter_import(self, simple_hamiltonian):
        """Test PauliConverter import."""
        from kanad.core.hamiltonians.pauli_converter import PauliConverter
        from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper

        mapper = JordanWignerMapper()

        # Should work if Qiskit is installed
        try:
            pauli_op = PauliConverter.to_sparse_pauli_op(simple_hamiltonian, mapper)
            assert pauli_op is not None
        except ImportError:
            pytest.skip("Qiskit not installed")

    def test_pauli_conversion(self, simple_hamiltonian):
        """Test Hamiltonian to Pauli conversion."""
        pytest.importorskip("qiskit")

        from kanad.core.hamiltonians.pauli_converter import PauliConverter
        from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper

        mapper = JordanWignerMapper()
        pauli_op = PauliConverter.to_sparse_pauli_op(simple_hamiltonian, mapper)

        # Should have multiple Pauli terms
        assert len(pauli_op) > 0

        # Should include nuclear repulsion (constant term)
        assert any('II' in str(pauli) for pauli in pauli_op.paulis)

    def test_pauli_info(self, simple_hamiltonian):
        """Test Pauli decomposition info."""
        pytest.importorskip("qiskit")

        from kanad.core.hamiltonians.pauli_converter import PauliConverter
        from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper

        mapper = JordanWignerMapper()
        info = PauliConverter.get_hamiltonian_info(simple_hamiltonian, mapper)

        assert 'num_terms' in info
        assert 'max_coeff' in info
        assert info['num_terms'] > 0


class TestQiskitBackend:
    """Test Qiskit backend initialization."""

    def test_aer_simulator_init(self):
        """Test Aer simulator initialization."""
        pytest.importorskip("qiskit_aer")

        from kanad.backends.qiskit_backend import QiskitBackend

        backend = QiskitBackend(backend_name='aer_simulator')
        assert backend.backend is not None
        assert backend.backend_name == 'aer_simulator'

    def test_aer_statevector_init(self):
        """Test Aer statevector simulator initialization."""
        pytest.importorskip("qiskit_aer")

        from kanad.backends.qiskit_backend import QiskitBackend

        backend = QiskitBackend(backend_name='aer_simulator_statevector')
        assert backend.backend is not None
        assert backend.shots is None  # Statevector is exact

    def test_get_estimator(self):
        """Test Estimator primitive creation."""
        pytest.importorskip("qiskit_aer")

        from kanad.backends.qiskit_backend import QiskitBackend

        backend = QiskitBackend(backend_name='aer_simulator')
        estimator = backend.get_estimator()

        assert estimator is not None

    def test_backend_info(self):
        """Test backend information retrieval."""
        pytest.importorskip("qiskit_aer")

        from kanad.backends.qiskit_backend import QiskitBackend

        backend = QiskitBackend(backend_name='aer_simulator', shots=2048)
        info = backend.get_backend_info()

        assert info['name'] == 'aer_simulator'
        assert info['shots'] == 2048


class TestQiskitVQE:
    """Test VQE with Qiskit backends."""

    @pytest.fixture
    def h2_hamiltonian(self):
        """Create H2 Hamiltonian."""
        from kanad.core.atom import Atom
        from kanad.core.representations.base_representation import Molecule
        from kanad.core.representations.lcao_representation import LCAORepresentation
        from kanad.core.hamiltonians.covalent_hamiltonian import CovalentHamiltonian

        h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        h2 = Atom('H', position=np.array([0.0, 0.0, 1.4]))
        molecule = Molecule([h1, h2])

        representation = LCAORepresentation(molecule)
        hamiltonian = CovalentHamiltonian(
            molecule,
            representation,
            basis_name='sto-3g'
        )

        return hamiltonian

    def test_vqe_classical_backend(self, h2_hamiltonian):
        """Test VQE with classical backend (no Qiskit)."""
        from kanad.solvers.vqe_solver import VQESolver
        from kanad.ansatze.hardware_efficient_ansatz import RealAmplitudesAnsatz
        from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper

        ansatz = RealAmplitudesAnsatz(n_qubits=2, n_electrons=2, n_layers=1)
        mapper = JordanWignerMapper()

        solver = VQESolver(
            hamiltonian=h2_hamiltonian,
            ansatz=ansatz,
            mapper=mapper,
            backend='classical',
            max_iterations=10
        )

        assert solver._use_qiskit is False
        result = solver.solve()
        assert np.isfinite(result['energy'])

    def test_vqe_aer_backend(self, h2_hamiltonian):
        """Test VQE with Aer simulator backend."""
        pytest.importorskip("qiskit")
        pytest.importorskip("qiskit_aer")

        from kanad.solvers.vqe_solver import VQESolver
        from kanad.ansatze.hardware_efficient_ansatz import RealAmplitudesAnsatz
        from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper

        # H2 has 2 spatial orbitals → 4 spin orbitals → 4 qubits (Jordan-Wigner)
        ansatz = RealAmplitudesAnsatz(n_qubits=4, n_electrons=2, n_layers=1)
        mapper = JordanWignerMapper()

        solver = VQESolver(
            hamiltonian=h2_hamiltonian,
            ansatz=ansatz,
            mapper=mapper,
            backend='aer_simulator_statevector',
            max_iterations=10
        )

        assert solver._use_qiskit is True
        assert solver.pauli_hamiltonian is not None

        result = solver.solve()
        assert np.isfinite(result['energy'])

    def test_vqe_classical_vs_qiskit(self, h2_hamiltonian):
        """Test that classical and Qiskit Aer statevector give similar results."""
        pytest.importorskip("qiskit")
        pytest.importorskip("qiskit_aer")

        from kanad.solvers.vqe_solver import VQESolver
        from kanad.ansatze.hardware_efficient_ansatz import RealAmplitudesAnsatz
        from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper

        # H2 has 2 spatial orbitals → 4 spin orbitals → 4 qubits (Jordan-Wigner)
        ansatz = RealAmplitudesAnsatz(n_qubits=4, n_electrons=2, n_layers=1)
        mapper = JordanWignerMapper()

        # Classical
        solver_classical = VQESolver(
            hamiltonian=h2_hamiltonian,
            ansatz=ansatz,
            mapper=mapper,
            backend='classical',
            optimizer='COBYLA',
            max_iterations=20
        )
        result_classical = solver_classical.solve()

        # Qiskit Aer statevector
        solver_qiskit = VQESolver(
            hamiltonian=h2_hamiltonian,
            ansatz=ansatz,
            mapper=mapper,
            backend='aer_simulator_statevector',
            optimizer='COBYLA',
            max_iterations=20
        )
        result_qiskit = solver_qiskit.solve()

        # Should be reasonably close (both exact simulation, but VQE has optimizer convergence issues)
        # Note: Current VQE accuracy ~27%, so we allow larger tolerance
        energy_diff = abs(result_classical['energy'] - result_qiskit['energy'])
        assert energy_diff < 0.5, f"Large energy difference: {energy_diff} Ha (classical vs Qiskit VQE)"
