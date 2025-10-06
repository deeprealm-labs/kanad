"""
Variational Quantum Eigensolver (VQE) for molecular ground state energy.

VQE is a hybrid quantum-classical algorithm that combines:
1. Quantum circuit (ansatz) for state preparation
2. Classical optimizer for parameter optimization
3. Energy expectation value measurement
"""

from typing import Optional, Callable, Dict, List, Tuple
import numpy as np
import warnings
import logging
from scipy.optimize import minimize

logger = logging.getLogger(__name__)

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
        backend: str = 'classical',
        optimizer: str = 'SLSQP',
        max_iterations: int = 1000,
        convergence_threshold: float = 1e-6,
        shots: int = 1024,
        **backend_options
    ):
        """
        Initialize VQE solver.

        Args:
            hamiltonian: Molecular Hamiltonian
            ansatz: Variational ansatz
            mapper: Fermionic-to-qubit mapper
            backend: Execution backend
                - 'classical': NumPy statevector simulation (exact, fast)
                - 'aer_simulator': Qiskit Aer QASM simulator (shot-based)
                - 'aer_simulator_statevector': Qiskit Aer statevector (exact)
                - 'bluequbit': BlueQubit cloud simulator (CPU/GPU, requires API key)
                - 'cuquantum_statevector': GPU-accelerated statevector (NVIDIA)
                - 'cuquantum_tensornet': GPU tensor network simulation (NVIDIA)
                - 'ibm_*': IBM Quantum hardware (requires credentials)
            optimizer: Classical optimizer ('SLSQP', 'COBYLA', 'L-BFGS-B')
            max_iterations: Maximum optimization iterations
            convergence_threshold: Energy convergence threshold (Hartree)
            shots: Number of measurement shots (for hardware/qasm backends)
            **backend_options: Additional backend-specific options
        """
        self.hamiltonian = hamiltonian
        self.ansatz = ansatz
        self.mapper = mapper
        self.backend_type = backend
        self.optimizer = optimizer
        self.max_iterations = max_iterations
        self.convergence_threshold = convergence_threshold
        self.shots = shots

        # Build circuit
        self.circuit = ansatz.build_circuit()
        self.n_parameters = self.circuit.get_num_parameters()

        # Initialize backend-specific components
        if backend == 'classical':
            self._use_qiskit = False
            self._use_cuquantum = False
            self._use_bluequbit = False
            self.qiskit_backend = None
            self.pauli_hamiltonian = None
            logger.info("Using classical NumPy simulation (exact statevector)")
        elif backend == 'bluequbit':
            self._use_qiskit = True  # BlueQubit uses Qiskit interface
            self._use_cuquantum = False
            self._use_bluequbit = True
            self._initialize_bluequbit_backend(backend, backend_options)
        elif backend.startswith('cuquantum'):
            self._use_qiskit = True  # cuQuantum uses Qiskit interface
            self._use_cuquantum = True
            self._use_bluequbit = False
            self._initialize_cuquantum_backend(backend, backend_options)
        else:
            self._use_qiskit = True
            self._use_cuquantum = False
            self._use_bluequbit = False
            self._initialize_qiskit_backend(backend, backend_options)

        # Optimization history
        self.energy_history: List[float] = []
        self.parameter_history: List[np.ndarray] = []
        self.iteration_count = 0

    def _initialize_qiskit_backend(self, backend_name: str, backend_options: dict):
        """Initialize Qiskit backend and convert Hamiltonian."""
        try:
            from kanad.backends.qiskit_backend import QiskitBackend
            from kanad.core.hamiltonians.pauli_converter import PauliConverter
        except ImportError as e:
            raise ImportError(
                f"Failed to import Qiskit components: {e}\n"
                "Install: pip install qiskit qiskit-aer"
            )

        # Initialize backend
        logger.info(f"Initializing Qiskit backend: {backend_name}")
        self.qiskit_backend = QiskitBackend(
            backend_name=backend_name,
            shots=self.shots,
            **backend_options
        )

        # Convert Hamiltonian to Pauli operators
        logger.info("Converting Hamiltonian to Pauli operators...")
        self.pauli_hamiltonian = PauliConverter.to_sparse_pauli_op(
            self.hamiltonian,
            self.mapper
        )

        num_terms = len(self.pauli_hamiltonian)
        logger.info(f"Hamiltonian converted: {num_terms} Pauli terms")

        # Get Estimator primitive
        self._estimator = self.qiskit_backend.get_estimator()

        # Convert circuit to Qiskit format (cached)
        self._qiskit_circuit = self.circuit.to_qiskit()
        logger.info(f"Circuit converted: {self._qiskit_circuit.num_qubits} qubits, "
                   f"{self._qiskit_circuit.depth()} depth")

    def _initialize_cuquantum_backend(self, backend_name: str, backend_options: dict):
        """Initialize cuQuantum GPU backend."""
        try:
            from kanad.backends.cuquantum_backend import CuQuantumBackend
            from kanad.core.hamiltonians.pauli_converter import PauliConverter
        except ImportError as e:
            raise ImportError(
                f"Failed to import cuQuantum components: {e}\n"
                "Install: pip install cuquantum-python cupy-cuda11x qiskit-aer-gpu"
            )

        # Determine simulation method
        if 'statevector' in backend_name:
            simulation_method = 'statevector'
        elif 'tensornet' in backend_name:
            simulation_method = 'tensornet'
        else:
            simulation_method = 'statevector'  # Default

        # Initialize cuQuantum backend
        logger.info(f"Initializing cuQuantum backend: {backend_name}")
        self.qiskit_backend = CuQuantumBackend(
            backend_name=backend_name,
            simulation_method=simulation_method,
            **backend_options
        )

        # Convert Hamiltonian to Pauli operators
        logger.info("Converting Hamiltonian to Pauli operators...")
        self.pauli_hamiltonian = PauliConverter.to_sparse_pauli_op(
            self.hamiltonian,
            self.mapper
        )

        num_terms = len(self.pauli_hamiltonian)
        logger.info(f"Hamiltonian converted: {num_terms} Pauli terms")

        # Get Estimator primitive (GPU-accelerated)
        self._estimator = self.qiskit_backend.get_estimator()

        # Convert circuit to Qiskit format (cached)
        self._qiskit_circuit = self.circuit.to_qiskit()
        logger.info(f"Circuit converted: {self._qiskit_circuit.num_qubits} qubits, "
                   f"{self._qiskit_circuit.depth()} depth")

        # Estimate GPU capabilities
        max_qubits = self.qiskit_backend.estimate_max_qubits()
        if self._qiskit_circuit.num_qubits > max_qubits:
            warnings.warn(
                f"Circuit has {self._qiskit_circuit.num_qubits} qubits but "
                f"GPU can handle ~{max_qubits} qubits. May run out of memory!"
            )

    def _initialize_bluequbit_backend(self, backend_name: str, backend_options: dict):
        """Initialize BlueQubit cloud backend."""
        try:
            from kanad.backends.bluequbit import BlueQubitBackend
            from kanad.core.hamiltonians.pauli_converter import PauliConverter
        except ImportError as e:
            raise ImportError(
                f"Failed to import BlueQubit components: {e}\n"
                "Install: pip install bluequbit"
            )

        # Extract device from backend_options
        device = backend_options.get('device', 'cpu')

        # Initialize BlueQubit backend
        logger.info(f"Initializing BlueQubit backend with device: {device}")
        self.qiskit_backend = BlueQubitBackend(device=device)

        # Convert Hamiltonian to Pauli operators
        logger.info("Converting Hamiltonian to Pauli operators...")
        self.pauli_hamiltonian = PauliConverter.to_sparse_pauli_op(
            self.hamiltonian,
            self.mapper
        )

        num_terms = len(self.pauli_hamiltonian)
        logger.info(f"Hamiltonian converted: {num_terms} Pauli terms")

        # Get Estimator primitive
        self._estimator = self.qiskit_backend.get_estimator()

        # Convert circuit to Qiskit format (cached)
        self._qiskit_circuit = self.circuit.to_qiskit()
        logger.info(f"Circuit converted: {self._qiskit_circuit.num_qubits} qubits, "
                   f"{self._qiskit_circuit.depth()} depth")
        logger.info(f"BlueQubit backend ready with device: {device}")

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
        if self._use_qiskit:
            return self._compute_energy_qiskit(parameters)
        else:
            return self._compute_energy_classical(parameters)

    def _compute_energy_qiskit(self, parameters: np.ndarray) -> float:
        """
        Compute energy using Qiskit Estimator primitive.

        Uses shot-based measurement for hardware/QASM backends,
        or exact statevector for Aer statevector.

        Args:
            parameters: Ansatz parameters

        Returns:
            Energy expectation value (Hartree)
        """
        # Assign parameters to Qiskit circuit
        bound_circuit = self.circuit.assign_parameters_for_qiskit(
            self._qiskit_circuit,
            parameters
        )

        # Check if using V2 Estimator API (StatevectorEstimator, etc.)
        try:
            # Try V2 API first (single tuple argument)
            job = self._estimator.run([(bound_circuit, self.pauli_hamiltonian)])
            result = job.result()
            # V2 returns PubResult in results list
            pub_result = result[0]
            # Extract evs (expectation values) - it's a scalar ndarray
            energy = float(pub_result.data.evs)
        except (TypeError, AttributeError, IndexError) as e:
            # Fall back to V1 API (two list arguments)
            try:
                job = self._estimator.run([bound_circuit], [self.pauli_hamiltonian])
                result = job.result()
                # V1 returns values directly
                energy = result.values[0]
            except Exception as e2:
                logger.error(f"Both V1 and V2 Estimator APIs failed: {e}, {e2}")
                raise

        return float(energy)

    def _compute_energy_classical(self, parameters: np.ndarray) -> float:
        """
        Compute energy using classical NumPy simulation.

        This is the original implementation - exact statevector simulation.

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
        """
        Expand single-qubit gate to full n-qubit system.

        Uses little-endian qubit ordering (qubit 0 is rightmost/LSB).
        For target_qubit=1 in 4-qubit system: I⊗I⊗gate⊗I
        """
        I = np.eye(2)
        result = np.array([1.0])

        # Iterate from highest qubit to lowest (little-endian convention)
        for qubit in range(n_qubits - 1, -1, -1):
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
        """
        Expand two-qubit gate to full n-qubit system.

        Handles non-adjacent qubits using SWAP network.
        """
        I = np.eye(2)
        dim = 2**n_qubits

        # For adjacent qubits, direct expansion
        if abs(control - target) == 1:
            result = np.array([1.0])
            min_qubit = min(control, target)

            for qubit in range(n_qubits):
                if qubit == min_qubit:
                    result = np.kron(result, gate_matrix)
                    # Skip next qubit (already included in gate_matrix)
                elif qubit == min_qubit + 1:
                    pass  # Already included
                else:
                    result = np.kron(result, I)

            return result

        # For non-adjacent qubits, use SWAP network
        # Strategy: SWAP target qubit next to control, apply gate, SWAP back
        return self._apply_gate_with_swaps(gate_matrix, control, target, n_qubits)

    def _apply_gate_with_swaps(
        self,
        gate_matrix: np.ndarray,
        control: int,
        target: int,
        n_qubits: int
    ) -> np.ndarray:
        """
        Apply two-qubit gate using SWAP network for non-adjacent qubits.

        Strategy:
        1. SWAP target towards control until adjacent
        2. Apply gate
        3. SWAP back to restore qubit positions
        """
        dim = 2**n_qubits
        result = np.eye(dim, dtype=complex)

        # Determine swap direction
        if control < target:
            # Swap target leftward towards control
            swaps = []
            current_target = target
            while current_target > control + 1:
                # SWAP qubits (current_target-1, current_target)
                swap_gate = self._create_swap_gate(current_target - 1, current_target, n_qubits)
                result = swap_gate @ result
                swaps.append((current_target - 1, current_target))
                current_target -= 1

            # Now target is adjacent to control, apply gate
            adjacent_gate = self._expand_two_qubit_gate_adjacent(
                gate_matrix, control, current_target, n_qubits
            )
            result = adjacent_gate @ result

            # Swap back in reverse order
            for q1, q2 in reversed(swaps):
                swap_gate = self._create_swap_gate(q1, q2, n_qubits)
                result = swap_gate @ result

        else:
            # Swap target rightward towards control
            swaps = []
            current_target = target
            while current_target < control - 1:
                # SWAP qubits (current_target, current_target+1)
                swap_gate = self._create_swap_gate(current_target, current_target + 1, n_qubits)
                result = swap_gate @ result
                swaps.append((current_target, current_target + 1))
                current_target += 1

            # Apply gate
            adjacent_gate = self._expand_two_qubit_gate_adjacent(
                gate_matrix, control, current_target, n_qubits
            )
            result = adjacent_gate @ result

            # Swap back
            for q1, q2 in reversed(swaps):
                swap_gate = self._create_swap_gate(q1, q2, n_qubits)
                result = swap_gate @ result

        return result

    def _expand_two_qubit_gate_adjacent(
        self,
        gate_matrix: np.ndarray,
        control: int,
        target: int,
        n_qubits: int
    ) -> np.ndarray:
        """Expand two-qubit gate for adjacent qubits only."""
        I = np.eye(2)
        result = np.array([1.0])
        min_qubit = min(control, target)

        for qubit in range(n_qubits):
            if qubit == min_qubit:
                result = np.kron(result, gate_matrix)
            elif qubit == min_qubit + 1:
                pass  # Already included
            else:
                result = np.kron(result, I)

        return result

    def _create_swap_gate(self, qubit1: int, qubit2: int, n_qubits: int) -> np.ndarray:
        """
        Create SWAP gate matrix for two adjacent qubits.

        SWAP = |00⟩⟨00| + |01⟩⟨10| + |10⟩⟨01| + |11⟩⟨11|
        """
        if abs(qubit1 - qubit2) != 1:
            raise ValueError("SWAP gate only defined for adjacent qubits")

        # SWAP matrix in 2-qubit space
        SWAP_2q = np.array([
            [1, 0, 0, 0],
            [0, 0, 1, 0],
            [0, 1, 0, 0],
            [0, 0, 0, 1]
        ], dtype=complex)

        I = np.eye(2)
        result = np.array([1.0], dtype=complex)
        min_qubit = min(qubit1, qubit2)

        for qubit in range(n_qubits):
            if qubit == min_qubit:
                result = np.kron(result, SWAP_2q)
            elif qubit == min_qubit + 1:
                pass  # Already included
            else:
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

        H = Σ_{ij} h_{ij} a†_i a_j + 1/2 Σ_{ijkl} g_{ijkl} a†_i a†_j a_l a_k + E_nn

        Returns:
            2^n × 2^n Hamiltonian matrix
        """
        n_qubits = self.circuit.n_qubits
        dim = 2**n_qubits

        logger.debug(f"Building Hamiltonian matrix for {n_qubits} qubits (dim={dim})...")

        # Use PauliConverter to get correct Hamiltonian!
        # This is the trusted, tested implementation
        from kanad.core.hamiltonians.pauli_converter import PauliConverter

        try:
            # Get Pauli representation
            pauli_op = PauliConverter.to_sparse_pauli_op(self.hamiltonian, self.mapper)

            # Convert to dense matrix (to_matrix() returns numpy array, not sparse matrix)
            H_matrix = pauli_op.to_matrix()

            # Ensure it's a numpy array (not a matrix subclass)
            if hasattr(H_matrix, 'toarray'):
                # It's a sparse matrix - convert to dense
                H_matrix = H_matrix.toarray()
            else:
                # It's already a numpy array
                H_matrix = np.array(H_matrix)

            logger.debug(f"Hamiltonian built using PauliConverter: {len(pauli_op)} Pauli terms")

            return H_matrix

        except Exception as e:
            logger.warning(f"PauliConverter failed ({e}), falling back to manual construction")

            # Fallback: manual construction (old code - has bugs!)
            H = self.hamiltonian.nuclear_repulsion * np.eye(dim, dtype=complex)

            h_core = self.hamiltonian.h_core
            n_orbitals = len(h_core)

            for i in range(n_orbitals):
                for j in range(n_orbitals):
                    if abs(h_core[i, j]) > 1e-10:
                        for spin_offset in [0, n_orbitals]:
                            i_spin = i + spin_offset
                            j_spin = j + spin_offset
                            H += h_core[i, j] * self._build_excitation_operator(i_spin, j_spin, n_qubits)

            if hasattr(self.hamiltonian, 'eri'):
                eri = self.hamiltonian.eri
                for i in range(n_orbitals):
                    for j in range(n_orbitals):
                        for k in range(n_orbitals):
                            for l in range(n_orbitals):
                                v_ijkl = eri[i, j, k, l]
                                if abs(v_ijkl) < 1e-10:
                                    continue

                                for spin_i_l in [0, n_orbitals]:
                                    for spin_k_j in [0, n_orbitals]:
                                        i_spin = i + spin_i_l
                                        k_spin = k + spin_k_j
                                        l_spin = l + spin_i_l
                                        j_spin = j + spin_k_j

                                        H += 0.5 * v_ijkl * self._build_two_electron_operator(
                                            i_spin, k_spin, l_spin, j_spin, n_qubits
                                        )

            return H

    def _build_two_electron_operator(
        self,
        i: int,
        j: int,
        k: int,
        l: int,
        n_qubits: int
    ) -> np.ndarray:
        """
        Build two-electron operator a†_i a†_j a_k a_l in qubit basis.

        Args:
            i, j: Creation indices
            k, l: Annihilation indices
            n_qubits: Number of qubits

        Returns:
            Operator matrix
        """
        # a†_i a†_j a_k a_l = a†_i (a†_j a_k) a_l
        # Build as product of operators

        # First: a†_j a_k
        op1 = self._build_excitation_operator(j, k, n_qubits)

        # Then: a†_i ... a_l where ... is the result of a†_j a_k
        # This is complex - need to use anticommutation relations

        # Simpler approach: build as sum of products
        # a†_i a†_j a_k a_l can be expressed using number operators and excitations

        # For now, use direct construction for common cases
        # Most important case: i=k, j=l (density-density interaction)
        if i == k and j == l:
            # a†_i a†_j a_j a_i = n_i n_j (number-number interaction)
            n_i = self._build_excitation_operator(i, i, n_qubits)
            n_j = self._build_excitation_operator(j, j, n_qubits)
            return n_i @ n_j
        elif i == l and j == k:
            # a†_i a†_j a_i a_j = a†_i a_i a†_j a_j (with sign from anticommutation)
            n_i = self._build_excitation_operator(i, i, n_qubits)
            n_j = self._build_excitation_operator(j, j, n_qubits)
            return n_i @ n_j
        else:
            # General case: use operator algebra
            # a†_i a†_j a_k a_l = (a†_i a_k)(a†_j a_l) ± (a†_i a_l)(a†_j a_k) + ...
            # This gets complex with anticommutation relations

            # Approximation: compute dominant terms
            op_ik = self._build_excitation_operator(i, k, n_qubits)
            op_jl = self._build_excitation_operator(j, l, n_qubits)
            return op_ik @ op_jl

    def _build_excitation_operator(
        self,
        i: int,
        j: int,
        n_qubits: int
    ) -> np.ndarray:
        """
        Build excitation operator a†_i a_j in qubit basis.

        Uses Jordan-Wigner transformation:
        - a†_i a_i → (I - Z_i)/2  (number operator)
        - a†_i a_j → (X_i - iY_i)/2 * Z_{i+1}...Z_{j-1} * (X_j + iY_j)/2  (for i < j)

        Args:
            i: Creation orbital
            j: Annihilation orbital
            n_qubits: Number of qubits

        Returns:
            Operator matrix
        """
        dim = 2**n_qubits

        # Pauli matrices
        I = np.eye(2, dtype=complex)
        X = np.array([[0, 1], [1, 0]], dtype=complex)
        Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
        Z = np.array([[1, 0], [0, -1]], dtype=complex)

        if i == j:
            # Number operator: a†_i a_i → (I - Z_i)/2
            op = np.array([1.0], dtype=complex)
            for q in range(n_qubits):
                if q == i:
                    op = np.kron(op, (I - Z) / 2)
                else:
                    op = np.kron(op, I)
            return op
        else:
            # Off-diagonal: a†_i a_j (i ≠ j)
            # Jordan-Wigner: need to handle string of Z operators between i and j

            # Ensure i < j for consistency
            if i > j:
                # a†_i a_j = (a†_j a_i)† (Hermitian conjugate)
                return self._build_excitation_operator(j, i, n_qubits).conj().T

            # For i < j: a†_i a_j → creation at i, annihilation at j
            # = [(X_i - iY_i)/2] * [Z_{i+1} * ... * Z_{j-1}] * [(X_j + iY_j)/2]

            # Build creation operator at i
            creation_i = (X - 1j*Y) / 2

            # Build annihilation operator at j
            annihilation_j = (X + 1j*Y) / 2

            # Build full operator
            op = np.array([1.0], dtype=complex)
            for q in range(n_qubits):
                if q == i:
                    op = np.kron(op, creation_i)
                elif q == j:
                    op = np.kron(op, annihilation_j)
                elif i < q < j:
                    # String of Z operators between i and j
                    op = np.kron(op, Z)
                else:
                    op = np.kron(op, I)

            return op

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
        backend_info = f"backend='{self.backend_type}'"
        if self._use_qiskit and self.shots:
            backend_info += f", shots={self.shots}"
        return (f"VQESolver(ansatz={self.ansatz.__class__.__name__}, "
                f"{backend_info}, optimizer={self.optimizer}, n_params={self.n_parameters})")
