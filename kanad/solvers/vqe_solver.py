"""
Variational Quantum Eigensolver (VQE) for molecular ground state energy.

VQE is a hybrid quantum-classical algorithm that combines:
1. Quantum circuit (ansatz) for state preparation
2. Classical optimizer for parameter optimization
3. Energy expectation value measurement
"""

from typing import Optional, Callable, Dict, List, Tuple
import numpy as np
from scipy.optimize import minimize

from kanad.ansatze.base_ansatz import BaseAnsatz
from kanad.core.hamiltonians.molecular_hamiltonian import MolecularHamiltonian
from kanad.core.mappers.base_mapper import BaseMapper


class VQESolver:
    """
    Variational Quantum Eigensolver.

    Algorithm:
    1. Prepare parametrized quantum state |ψ(θ)⟩ using ansatz
    2. Measure energy E(θ) = ⟨ψ(θ)|H|ψ(θ)⟩
    3. Optimize θ to minimize E(θ)
    4. Return optimized energy and parameters
    """

    def __init__(
        self,
        hamiltonian: MolecularHamiltonian,
        ansatz: BaseAnsatz,
        mapper: BaseMapper,
        optimizer: str = 'SLSQP',
        max_iterations: int = 1000,
        convergence_threshold: float = 1e-6
    ):
        """
        Initialize VQE solver.

        Args:
            hamiltonian: Molecular Hamiltonian
            ansatz: Variational ansatz
            mapper: Fermionic-to-qubit mapper
            optimizer: Classical optimizer ('SLSQP', 'COBYLA', 'L-BFGS-B')
            max_iterations: Maximum optimization iterations
            convergence_threshold: Energy convergence threshold (Hartree)
        """
        self.hamiltonian = hamiltonian
        self.ansatz = ansatz
        self.mapper = mapper
        self.optimizer = optimizer
        self.max_iterations = max_iterations
        self.convergence_threshold = convergence_threshold

        # Build circuit
        self.circuit = ansatz.build_circuit()
        self.n_parameters = self.circuit.get_num_parameters()

        # Optimization history
        self.energy_history: List[float] = []
        self.parameter_history: List[np.ndarray] = []
        self.iteration_count = 0

    def solve(
        self,
        initial_parameters: Optional[np.ndarray] = None,
        callback: Optional[Callable] = None
    ) -> Dict:
        """
        Solve for ground state energy using VQE.

        Args:
            initial_parameters: Initial parameter values
            callback: Optional callback function(iteration, energy, parameters)

        Returns:
            Dictionary with results:
                - energy: Ground state energy (Hartree)
                - parameters: Optimized parameters
                - converged: Whether optimization converged
                - iterations: Number of iterations
                - energy_history: Energy at each iteration
        """
        # Initialize parameters
        if initial_parameters is None:
            initial_parameters = self.ansatz.initialize_parameters('small_random')

        # Reset history
        self.energy_history = []
        self.parameter_history = []
        self.iteration_count = 0

        # Define objective function
        def objective(params):
            energy = self.compute_energy(params)
            self.energy_history.append(energy)
            self.parameter_history.append(params.copy())
            self.iteration_count += 1

            if callback is not None:
                callback(self.iteration_count, energy, params)

            return energy

        # Optimize
        result = minimize(
            objective,
            initial_parameters,
            method=self.optimizer,
            options={'maxiter': self.max_iterations}
        )

        # Check convergence
        converged = result.success and result.fun < self.energy_history[0]

        return {
            'energy': result.fun,
            'parameters': result.x,
            'converged': converged,
            'iterations': self.iteration_count,
            'energy_history': np.array(self.energy_history),
            'parameter_history': self.parameter_history,
            'optimizer_message': result.message
        }

    def compute_energy(self, parameters: np.ndarray) -> float:
        """
        Compute energy expectation value for given parameters.

        E(θ) = ⟨ψ(θ)|H|ψ(θ)⟩

        Args:
            parameters: Ansatz parameters

        Returns:
            Energy expectation value (Hartree)
        """
        # Bind parameters to circuit
        self.circuit.bind_parameters(parameters)

        # Get state vector (classical simulation)
        state_vector = self._simulate_circuit()

        # Compute energy expectation
        energy = self._compute_expectation_value(state_vector)

        return energy

    def _simulate_circuit(self) -> np.ndarray:
        """
        Simulate quantum circuit classically.

        Returns:
            State vector (complex amplitudes)
        """
        n_qubits = self.circuit.n_qubits
        state = np.zeros(2**n_qubits, dtype=complex)
        state[0] = 1.0  # |0...0⟩

        # Apply gates sequentially
        for gate in self.circuit.gates:
            gate_type = gate['type']
            qubits = gate['qubits']
            params = gate.get('params', [])

            if gate_type == 'barrier':
                continue

            # Get gate matrix
            matrix = self._get_gate_matrix(gate_type, qubits, params, n_qubits)

            # Apply to state
            state = matrix @ state

        return state

    def _get_gate_matrix(
        self,
        gate_type: str,
        qubits: List[int],
        params: List,
        n_qubits: int
    ) -> np.ndarray:
        """
        Get full gate matrix for n-qubit system.

        Args:
            gate_type: Gate type
            qubits: Qubit indices
            params: Gate parameters
            n_qubits: Total number of qubits

        Returns:
            2^n × 2^n gate matrix
        """
        # Single-qubit gates
        if gate_type in ['h', 'x', 'y', 'z', 'rx', 'ry', 'rz']:
            base_matrix = self._get_single_qubit_gate(gate_type, params)
            return self._expand_gate_to_full_system(base_matrix, qubits[0], n_qubits)

        # Two-qubit gates
        elif gate_type in ['cx', 'cz', 'rxx', 'ryy', 'rzz', 'swap']:
            base_matrix = self._get_two_qubit_gate(gate_type, params)
            return self._expand_two_qubit_gate(base_matrix, qubits[0], qubits[1], n_qubits)

        else:
            raise ValueError(f"Unknown gate type: {gate_type}")

    def _get_single_qubit_gate(self, gate_type: str, params: List) -> np.ndarray:
        """Get 2×2 single-qubit gate matrix."""
        if gate_type == 'h':
            return np.array([[1, 1], [1, -1]]) / np.sqrt(2)
        elif gate_type == 'x':
            return np.array([[0, 1], [1, 0]])
        elif gate_type == 'y':
            return np.array([[0, -1j], [1j, 0]])
        elif gate_type == 'z':
            return np.array([[1, 0], [0, -1]])
        elif gate_type == 'rx':
            theta = params[0].value if hasattr(params[0], 'value') else params[0]
            return np.array([
                [np.cos(theta/2), -1j*np.sin(theta/2)],
                [-1j*np.sin(theta/2), np.cos(theta/2)]
            ])
        elif gate_type == 'ry':
            theta = params[0].value if hasattr(params[0], 'value') else params[0]
            return np.array([
                [np.cos(theta/2), -np.sin(theta/2)],
                [np.sin(theta/2), np.cos(theta/2)]
            ])
        elif gate_type == 'rz':
            theta = params[0].value if hasattr(params[0], 'value') else params[0]
            return np.array([
                [np.exp(-1j*theta/2), 0],
                [0, np.exp(1j*theta/2)]
            ])

    def _get_two_qubit_gate(self, gate_type: str, params: List) -> np.ndarray:
        """Get 4×4 two-qubit gate matrix."""
        if gate_type == 'cx':
            return np.array([
                [1, 0, 0, 0],
                [0, 1, 0, 0],
                [0, 0, 0, 1],
                [0, 0, 1, 0]
            ])
        elif gate_type == 'cz':
            return np.array([
                [1, 0, 0, 0],
                [0, 1, 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, -1]
            ])
        elif gate_type == 'swap':
            return np.array([
                [1, 0, 0, 0],
                [0, 0, 1, 0],
                [0, 1, 0, 0],
                [0, 0, 0, 1]
            ])
        elif gate_type == 'rxx':
            theta = params[0].value if hasattr(params[0], 'value') else params[0]
            c = np.cos(theta/2)
            s = np.sin(theta/2)
            return np.array([
                [c, 0, 0, -1j*s],
                [0, c, -1j*s, 0],
                [0, -1j*s, c, 0],
                [-1j*s, 0, 0, c]
            ])
        elif gate_type == 'ryy':
            theta = params[0].value if hasattr(params[0], 'value') else params[0]
            c = np.cos(theta/2)
            s = np.sin(theta/2)
            return np.array([
                [c, 0, 0, 1j*s],
                [0, c, -1j*s, 0],
                [0, -1j*s, c, 0],
                [1j*s, 0, 0, c]
            ])
        elif gate_type == 'rzz':
            theta = params[0].value if hasattr(params[0], 'value') else params[0]
            return np.array([
                [np.exp(-1j*theta/2), 0, 0, 0],
                [0, np.exp(1j*theta/2), 0, 0],
                [0, 0, np.exp(1j*theta/2), 0],
                [0, 0, 0, np.exp(-1j*theta/2)]
            ])

    def _expand_gate_to_full_system(
        self,
        gate_matrix: np.ndarray,
        target_qubit: int,
        n_qubits: int
    ) -> np.ndarray:
        """Expand single-qubit gate to full n-qubit system."""
        # Build tensor product: I ⊗ ... ⊗ gate ⊗ ... ⊗ I
        I = np.eye(2)
        result = np.array([1.0])

        for qubit in range(n_qubits):
            if qubit == target_qubit:
                result = np.kron(result, gate_matrix)
            else:
                result = np.kron(result, I)

        return result

    def _expand_two_qubit_gate(
        self,
        gate_matrix: np.ndarray,
        control: int,
        target: int,
        n_qubits: int
    ) -> np.ndarray:
        """Expand two-qubit gate to full n-qubit system."""
        # Simplified: only handle adjacent qubits
        # Full implementation would use SWAP networks
        if abs(control - target) != 1:
            raise NotImplementedError("Non-adjacent two-qubit gates require SWAP network")

        I = np.eye(2)
        result = np.array([1.0])
        min_qubit = min(control, target)

        for qubit in range(n_qubits):
            if qubit == min_qubit:
                result = np.kron(result, gate_matrix)
                qubit += 1  # Skip next qubit (already included)
            elif qubit != min_qubit + 1:
                result = np.kron(result, I)

        return result

    def _compute_expectation_value(self, state_vector: np.ndarray) -> float:
        """
        Compute Hamiltonian expectation value.

        E = ⟨ψ|H|ψ⟩

        Args:
            state_vector: Quantum state

        Returns:
            Energy (Hartree)
        """
        # For simplicity, use matrix form of Hamiltonian
        # Full implementation would use Pauli decomposition from mapper

        # Get Hamiltonian matrix
        H_matrix = self._build_hamiltonian_matrix()

        # Compute expectation
        energy = np.real(np.conjugate(state_vector) @ H_matrix @ state_vector)

        return energy

    def _build_hamiltonian_matrix(self) -> np.ndarray:
        """
        Build full Hamiltonian matrix in qubit basis.

        Returns:
            2^n × 2^n Hamiltonian matrix
        """
        n_qubits = self.circuit.n_qubits
        dim = 2**n_qubits

        # Simplified: use second quantization approach
        # Map fermionic Hamiltonian to qubit Hamiltonian using mapper

        # For now, return identity (placeholder)
        # Full implementation would use mapper.map_hamiltonian()
        H = np.zeros((dim, dim), dtype=complex)

        # Add one-body terms
        h_core = self.hamiltonian.h_core
        for i in range(len(h_core)):
            for j in range(len(h_core)):
                if abs(h_core[i, j]) > 1e-10:
                    # Map a†_i a_j to qubits
                    H += h_core[i, j] * self._build_excitation_operator(i, j, n_qubits)

        # Add nuclear repulsion (constant term)
        H += self.hamiltonian.nuclear_repulsion * np.eye(dim)

        return H

    def _build_excitation_operator(
        self,
        i: int,
        j: int,
        n_qubits: int
    ) -> np.ndarray:
        """
        Build excitation operator a†_i a_j in qubit basis.

        Uses Jordan-Wigner transformation (simplified).

        Args:
            i: Creation orbital
            j: Annihilation orbital
            n_qubits: Number of qubits

        Returns:
            Operator matrix
        """
        dim = 2**n_qubits

        if i == j:
            # Number operator: a†_i a_i → (I - Z_i)/2
            Z = np.array([[1, 0], [0, -1]])
            I = np.eye(2)

            op = np.array([1.0])
            for q in range(n_qubits):
                if q == i:
                    op = np.kron(op, (I - Z) / 2)
                else:
                    op = np.kron(op, I)

            return op
        else:
            # Simplified: return zeros
            # Full implementation would use proper Jordan-Wigner transformation
            return np.zeros((dim, dim), dtype=complex)

    def get_energy_variance(self, parameters: np.ndarray) -> float:
        """
        Compute energy variance for given parameters.

        Var(E) = ⟨ψ|H²|ψ⟩ - ⟨ψ|H|ψ⟩²

        Args:
            parameters: Ansatz parameters

        Returns:
            Energy variance
        """
        self.circuit.bind_parameters(parameters)
        state = self._simulate_circuit()
        H = self._build_hamiltonian_matrix()

        E = np.real(np.conjugate(state) @ H @ state)
        E2 = np.real(np.conjugate(state) @ H @ H @ state)

        return E2 - E**2

    def __repr__(self) -> str:
        return (f"VQESolver(ansatz={self.ansatz.__class__.__name__}, "
                f"optimizer={self.optimizer}, n_params={self.n_parameters})")
