"""
Unit tests for VQE solver and energy estimation.

Tests:
- VQE initialization
- Energy computation
- Optimization
- Circuit simulation
- Integration with ansätze and Hamiltonians
"""

import pytest
import numpy as np

from kanad.solvers.vqe_solver import VQESolver
from kanad.ansatze.hardware_efficient_ansatz import RealAmplitudesAnsatz
from kanad.ansatze.ucc_ansatz import UCC_S_Ansatz
from kanad.core.hamiltonians.covalent_hamiltonian import CovalentHamiltonian
from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper
from kanad.core.atom import Atom
from kanad.core.representations.base_representation import Molecule
from kanad.core.representations.lcao_representation import LCAORepresentation


class TestVQESolver:
    """Test VQE solver functionality."""

    @pytest.fixture
    def h2_molecule(self):
        """Create H2 molecule for testing."""
        h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
        h2 = Atom('H', position=np.array([0.0, 0.0, 1.4]))  # 1.4 bohr ≈ 0.74 Å
        return Molecule([h1, h2])

    @pytest.fixture
    def h2_hamiltonian(self, h2_molecule):
        """Create H2 Hamiltonian."""
        representation = LCAORepresentation(h2_molecule)
        hamiltonian = CovalentHamiltonian(
            h2_molecule,
            representation,
            basis_name='sto-3g'
        )
        return hamiltonian

    def test_vqe_solver_creation(self, h2_hamiltonian):
        """Test VQE solver initialization."""
        ansatz = RealAmplitudesAnsatz(n_qubits=2, n_electrons=2, n_layers=1)
        mapper = JordanWignerMapper()

        solver = VQESolver(
            hamiltonian=h2_hamiltonian,
            ansatz=ansatz,
            mapper=mapper
        )

        assert solver.n_parameters > 0
        assert len(solver.energy_history) == 0

    def test_vqe_energy_computation(self, h2_hamiltonian):
        """Test energy computation for given parameters."""
        ansatz = RealAmplitudesAnsatz(n_qubits=2, n_electrons=2, n_layers=1)
        mapper = JordanWignerMapper()

        solver = VQESolver(
            hamiltonian=h2_hamiltonian,
            ansatz=ansatz,
            mapper=mapper
        )

        # Compute energy with zero parameters
        params = np.zeros(solver.n_parameters)
        energy = solver.compute_energy(params)

        # Energy should be real and finite
        assert np.isfinite(energy)
        assert isinstance(energy, (float, np.floating))

    def test_vqe_optimization(self, h2_hamiltonian):
        """Test VQE optimization."""
        ansatz = RealAmplitudesAnsatz(n_qubits=2, n_electrons=2, n_layers=1)
        mapper = JordanWignerMapper()

        solver = VQESolver(
            hamiltonian=h2_hamiltonian,
            ansatz=ansatz,
            mapper=mapper,
            optimizer='COBYLA',  # Fast for testing
            max_iterations=50
        )

        result = solver.solve()

        # Check result structure
        assert 'energy' in result
        assert 'parameters' in result
        assert 'converged' in result
        assert 'iterations' in result

        # Energy should be finite
        assert np.isfinite(result['energy'])

        # Should have run some iterations
        assert result['iterations'] > 0

    def test_vqe_with_ucc_ansatz(self, h2_hamiltonian):
        """Test VQE with UCC ansatz."""
        ansatz = UCC_S_Ansatz(n_qubits=2, n_electrons=2)
        mapper = JordanWignerMapper()

        solver = VQESolver(
            hamiltonian=h2_hamiltonian,
            ansatz=ansatz,
            mapper=mapper,
            max_iterations=30
        )

        result = solver.solve()
        assert np.isfinite(result['energy'])

    def test_vqe_callback(self, h2_hamiltonian):
        """Test VQE callback functionality."""
        ansatz = RealAmplitudesAnsatz(n_qubits=2, n_electrons=2, n_layers=1)
        mapper = JordanWignerMapper()

        solver = VQESolver(
            hamiltonian=h2_hamiltonian,
            ansatz=ansatz,
            mapper=mapper,
            max_iterations=10
        )

        callback_data = []

        def callback(iteration, energy, params):
            callback_data.append({
                'iteration': iteration,
                'energy': energy,
                'params': params.copy()
            })

        result = solver.solve(callback=callback)

        # Callback should have been called
        assert len(callback_data) > 0
        assert callback_data[0]['iteration'] == 1

    def test_vqe_energy_history(self, h2_hamiltonian):
        """Test energy history tracking."""
        ansatz = RealAmplitudesAnsatz(n_qubits=2, n_electrons=2, n_layers=1)
        mapper = JordanWignerMapper()

        solver = VQESolver(
            hamiltonian=h2_hamiltonian,
            ansatz=ansatz,
            mapper=mapper,
            max_iterations=20
        )

        result = solver.solve()

        # Energy history should exist
        assert len(result['energy_history']) > 0
        assert len(result['energy_history']) == result['iterations']

        # Final energy should match
        assert np.isclose(result['energy_history'][-1], result['energy'])

    def test_vqe_variance_computation(self, h2_hamiltonian):
        """Test energy variance computation."""
        ansatz = RealAmplitudesAnsatz(n_qubits=2, n_electrons=2, n_layers=1)
        mapper = JordanWignerMapper()

        solver = VQESolver(
            hamiltonian=h2_hamiltonian,
            ansatz=ansatz,
            mapper=mapper
        )

        params = np.zeros(solver.n_parameters)
        variance = solver.get_energy_variance(params)

        # Variance should be non-negative
        assert variance >= 0


class TestCircuitSimulation:
    """Test quantum circuit simulation."""

    def test_single_qubit_gates(self):
        """Test single-qubit gate simulation."""
        from kanad.ansatze.base_ansatz import QuantumCircuit, Parameter

        circuit = QuantumCircuit(1)
        circuit.x(0)

        # Simple test - just ensure no errors
        # Full simulation tested through VQE

    def test_two_qubit_gates(self):
        """Test two-qubit gate simulation."""
        from kanad.ansatze.base_ansatz import QuantumCircuit

        circuit = QuantumCircuit(2)
        circuit.h(0)
        circuit.cx(0, 1)  # Create Bell state

        # Full simulation tested through VQE

    def test_parametrized_gates(self):
        """Test parametrized gate simulation."""
        from kanad.ansatze.base_ansatz import QuantumCircuit, Parameter

        circuit = QuantumCircuit(2)
        theta = Parameter('theta', np.pi/4)
        circuit.ry(theta, 0)

        # Full simulation tested through VQE


class TestVQEIntegration:
    """Integration tests for VQE with different components."""

    @pytest.fixture
    def simple_hamiltonian(self):
        """Create simple test Hamiltonian."""
        from kanad.core.hamiltonians.molecular_hamiltonian import MolecularHamiltonian

        class SimpleHamiltonian(MolecularHamiltonian):
            def __init__(self):
                super().__init__(n_orbitals=2, n_electrons=2, nuclear_repulsion=0.5)
                # Simple Hamiltonian
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

    def test_vqe_with_different_optimizers(self, simple_hamiltonian):
        """Test VQE with different optimizers."""
        ansatz = RealAmplitudesAnsatz(n_qubits=2, n_electrons=2, n_layers=1)
        mapper = JordanWignerMapper()

        for optimizer in ['COBYLA', 'SLSQP']:
            solver = VQESolver(
                hamiltonian=simple_hamiltonian,
                ansatz=ansatz,
                mapper=mapper,
                optimizer=optimizer,
                max_iterations=10
            )

            result = solver.solve()
            assert np.isfinite(result['energy'])

    def test_vqe_parameter_initialization(self, simple_hamiltonian):
        """Test VQE with different parameter initializations."""
        ansatz = RealAmplitudesAnsatz(n_qubits=2, n_electrons=2, n_layers=1)
        mapper = JordanWignerMapper()

        solver = VQESolver(
            hamiltonian=simple_hamiltonian,
            ansatz=ansatz,
            mapper=mapper,
            max_iterations=10
        )

        # Random initialization
        result1 = solver.solve()
        assert np.isfinite(result1['energy'])

        # Zero initialization
        params_zero = np.zeros(solver.n_parameters)
        result2 = solver.solve(initial_parameters=params_zero)
        assert np.isfinite(result2['energy'])

    def test_vqe_reproducibility(self, simple_hamiltonian):
        """Test VQE gives reproducible results with same initialization."""
        ansatz = RealAmplitudesAnsatz(n_qubits=2, n_electrons=2, n_layers=1)
        mapper = JordanWignerMapper()

        solver = VQESolver(
            hamiltonian=simple_hamiltonian,
            ansatz=ansatz,
            mapper=mapper,
            max_iterations=20
        )

        # Same initial parameters
        n_params = solver.n_parameters
        params_init = np.random.uniform(-0.1, 0.1, n_params)

        result1 = solver.solve(initial_parameters=params_init.copy())
        result2 = solver.solve(initial_parameters=params_init.copy())

        # Should give same results (deterministic optimization)
        # Allow some tolerance for numerical differences
        assert np.allclose(result1['energy'], result2['energy'], rtol=1e-3)

    def test_vqe_converges_below_initial(self, simple_hamiltonian):
        """Test VQE converges to energy lower than initial."""
        ansatz = RealAmplitudesAnsatz(n_qubits=2, n_electrons=2, n_layers=1)
        mapper = JordanWignerMapper()

        solver = VQESolver(
            hamiltonian=simple_hamiltonian,
            ansatz=ansatz,
            mapper=mapper,
            max_iterations=50
        )

        # Start from random parameters
        params_init = ansatz.initialize_parameters('random')

        result = solver.solve(initial_parameters=params_init)

        # Final energy should be lower than or equal to initial
        # (VQE minimizes energy)
        initial_energy = result['energy_history'][0]
        final_energy = result['energy_history'][-1]

        # Allow for numerical tolerance and potential local minima
        assert final_energy <= initial_energy + 1e-6
