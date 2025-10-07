"""
Variational Quantum Eigensolver (VQE) - Rebuilt with Bonds Module Integration.

New Design:
- Takes bond as input (not raw Hamiltonian)
- Automatic analysis integration
- Automatic circuit optimization
- Rich, comprehensive results
- User-friendly interface
"""

from typing import Dict, Any, Optional, Callable
import numpy as np
import logging
from scipy.optimize import minimize

from kanad.solvers.base_solver import BaseSolver
from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper
from kanad.core.mappers.bravyi_kitaev_mapper import BravyiKitaevMapper

logger = logging.getLogger(__name__)


class VQESolver(BaseSolver):
    """
    Variational Quantum Eigensolver for ground state energy.

    VQE is a hybrid quantum-classical algorithm:
    1. Prepare parametrized quantum state |ψ(θ)⟩
    2. Measure energy E(θ) = ⟨ψ(θ)|H|ψ(θ)⟩
    3. Classically optimize parameters θ
    4. Return minimum energy found

    Usage:
        from kanad.bonds import BondFactory
        from kanad.solvers import VQESolver

        bond = BondFactory.create_bond('H', 'H', distance=0.74)
        solver = VQESolver(bond, ansatz_type='ucc')
        result = solver.solve()

        print(f"Energy: {result['energy']:.6f} Hartree")
        solver.print_summary()
    """

    def __init__(
        self,
        bond: Optional['BaseBond'] = None,
        # High-level API (bond-based)
        ansatz_type: str = 'ucc',
        mapper_type: str = 'jordan_wigner',
        # Low-level API (component-based, for testing)
        hamiltonian: Optional[Any] = None,
        ansatz: Optional[Any] = None,
        mapper: Optional[Any] = None,
        # Common parameters
        optimizer: str = 'SLSQP',
        max_iterations: int = 1000,
        conv_threshold: float = 1e-6,
        backend: str = 'statevector',
        shots: Optional[int] = None,
        enable_analysis: bool = True,
        enable_optimization: bool = True,
        **kwargs
    ):
        """
        Initialize VQE solver.

        Supports two APIs:
        1. High-level (bond-based): solver = VQESolver(bond, ansatz_type='ucc')
        2. Low-level (component-based): solver = VQESolver(hamiltonian=ham, ansatz=ans, mapper=map)

        Args:
            bond: Bond object from BondFactory (high-level API)
            ansatz_type: Type of ansatz ('ucc', 'hardware_efficient', 'governance')
            mapper_type: Fermionic-to-qubit mapping ('jordan_wigner', 'parity', 'bravyi_kitaev')
            hamiltonian: Hamiltonian object (low-level API, for testing)
            ansatz: Ansatz object (low-level API, for testing)
            mapper: Mapper object (low-level API, for testing)
            optimizer: Classical optimizer ('SLSQP', 'COBYLA', 'L-BFGS-B')
            max_iterations: Maximum optimization iterations
            conv_threshold: Convergence threshold (Hartree)
            backend: Quantum backend ('statevector', 'qasm', 'bluequbit', 'ibm')
            shots: Number of shots for sampling backends
            enable_analysis: Enable automatic analysis
            enable_optimization: Enable automatic circuit optimization
            **kwargs: Additional backend-specific options
        """
        # Store ansatz/mapper types FIRST (needed by _init methods)
        self._ansatz_type_param = ansatz_type
        self._mapper_type_param = mapper_type

        # Detect which API is being used
        if bond is not None and hamiltonian is None:
            # High-level API: Initialize from bond
            super().__init__(bond, enable_analysis, enable_optimization)
            self._api_mode = 'bond'
        elif hamiltonian is not None and bond is None:
            # Low-level API: Initialize from components (for testing)
            self._init_from_components_mode(hamiltonian, ansatz, mapper, enable_analysis, enable_optimization)
            self._api_mode = 'components'
        elif bond is not None and hamiltonian is not None:
            raise ValueError("Cannot use both 'bond' and 'hamiltonian' parameters. Choose one API.")
        else:
            raise ValueError("Must provide either 'bond' (high-level API) or 'hamiltonian' (low-level API)")

        # Store common parameters
        self.optimizer_method = optimizer  # Support both 'optimizer' and 'optimizer_method' kwargs
        self.max_iterations = max_iterations
        self.conv_threshold = conv_threshold
        # Normalize backend names ('classical' -> 'statevector')
        self.backend_name = 'statevector' if backend == 'classical' else backend
        self.backend = self.backend_name  # Alias for backward compatibility
        self.shots = shots if shots is not None else 1024

        # This is a correlated method
        self._is_correlated = True

        # Initialize based on API mode
        if self._api_mode == 'bond':
            # Initialize ansatz and mapper from bond
            self._init_ansatz()
            self._init_mapper()
            # Build quantum circuit
            self._build_circuit()
            # Initialize backend
            self._init_backend(**kwargs)
        else:
            # Components mode - ansatz and mapper already set
            # Initialize backend
            self._hamiltonian_matrix = None
            self._use_statevector = True
            # Determine if Qiskit backend is requested
            self._use_qiskit = backend not in ['statevector', 'classical', None]

            # If Qiskit backend, convert Hamiltonian to Pauli operators
            if self._use_qiskit:
                # Convert Hamiltonian matrix to Pauli operators using mapper
                try:
                    from qiskit.quantum_info import SparsePauliOp
                    # Use mapper to convert fermionic Hamiltonian to qubit operators
                    # For now, create a simple Pauli representation
                    # TODO: Implement full Hamiltonian mapping
                    self.pauli_hamiltonian = SparsePauliOp.from_list([("I" * self.ansatz.n_qubits, 1.0)])
                    logger.info("Created placeholder Pauli Hamiltonian for Qiskit backend")
                except ImportError:
                    logger.warning("Qiskit not available for Pauli conversion")
                    self.pauli_hamiltonian = None
            else:
                self.pauli_hamiltonian = None

            logger.info(f"VQE initialized in component mode (for testing), use_qiskit={self._use_qiskit}")

        # Optimization tracking
        self.energy_history = []
        self.parameter_history = []
        self.iteration_count = 0

        logger.info(f"VQE Solver initialized: {self.ansatz_type} ansatz, {self.mapper_type} mapping, {self.backend_name} backend")

    def _init_from_components_mode(self, hamiltonian, ansatz, mapper, enable_analysis, enable_optimization):
        """Initialize VQE from individual components (low-level API for testing)."""
        # Store components directly
        self.hamiltonian = hamiltonian
        self.mapper = mapper if mapper is not None else JordanWignerMapper()

        # Set molecule and bond to None (not available in this mode)
        self.molecule = getattr(hamiltonian, 'molecule', None)
        self.bond = None

        # Disable analysis features unless molecule is available
        self.enable_analysis = enable_analysis if self.molecule is not None else False
        self.enable_optimization = enable_optimization

        # Initialize ansatz (from object or type string)
        if ansatz is not None:
            # Ansatz object provided directly
            self.ansatz = ansatz
            self.ansatz_type = type(self.ansatz).__name__
        elif self._ansatz_type_param is not None:
            # Ansatz type string provided - create ansatz
            self.ansatz_type = self._ansatz_type_param
            self._init_ansatz()
        else:
            self.ansatz = None
            self.ansatz_type = 'None'

        # Set mapper type string for logging
        self.mapper_type = type(self.mapper).__name__ if self.mapper else 'None'

        # Build circuit to get n_parameters
        if self.ansatz is not None:
            if self.ansatz.circuit is None:
                self.ansatz.build_circuit()

            # Get n_parameters from various possible sources
            if hasattr(self.ansatz, 'n_parameters'):
                self.n_parameters = self.ansatz.n_parameters
            elif hasattr(self.ansatz.circuit, 'get_num_parameters'):
                self.n_parameters = self.ansatz.circuit.get_num_parameters()
            elif hasattr(self.ansatz, 'num_parameters'):
                self.n_parameters = self.ansatz.num_parameters
            else:
                self.n_parameters = 0
        else:
            self.n_parameters = 0

        logger.info("VQE solver initialized from components (testing mode)")

    def _init_ansatz(self):
        """Initialize ansatz from bonds module."""
        from kanad.bonds import UCCAnsatz, HardwareEfficientAnsatz

        n_qubits = 2 * self.hamiltonian.n_orbitals
        n_electrons = self.molecule.n_electrons

        # Use stored parameter (set in __init__)
        ansatz_type = self._ansatz_type_param

        if ansatz_type.lower() == 'ucc':
            self.ansatz = UCCAnsatz(
                n_qubits=n_qubits,
                n_electrons=n_electrons,
                include_singles=True,
                include_doubles=True
            )
            logger.info(f"UCC ansatz: {n_qubits} qubits, {n_electrons} electrons")

        elif ansatz_type.lower() == 'hardware_efficient':
            self.ansatz = HardwareEfficientAnsatz(
                n_qubits=n_qubits,
                n_electrons=n_electrons,
                n_layers=3,
                entanglement='linear'
            )
            logger.info(f"Hardware-efficient ansatz: 3 layers, linear entanglement")

        elif ansatz_type.lower() == 'governance':
            # Use governance-aware ansatz based on bond type or Hamiltonian protocol
            bond_type = None
            if self.bond is not None:
                bond_type = self.bond.bond_type
            elif hasattr(self.hamiltonian, 'governance_protocol') and self.hamiltonian.governance_protocol:
                # Infer from protocol type
                protocol_name = type(self.hamiltonian.governance_protocol).__name__
                if 'Covalent' in protocol_name:
                    bond_type = 'covalent'
                elif 'Ionic' in protocol_name:
                    bond_type = 'ionic'
                elif 'Metallic' in protocol_name:
                    bond_type = 'metallic'

            if bond_type == 'covalent':
                from kanad.ansatze.governance_aware_ansatz import CovalentGovernanceAnsatz

                # Get hybridization from governance metadata if available
                metadata = getattr(self.hamiltonian, '_governance_metadata', {})
                hybridization = metadata.get('hybridization', 'sp3')
                protocol = metadata.get('governance_protocol', self.hamiltonian.governance_protocol)

                self.ansatz = CovalentGovernanceAnsatz(
                    n_qubits=n_qubits,
                    n_electrons=n_electrons,
                    n_layers=2,
                    hybridization=hybridization,
                    protocol=protocol
                )
                logger.info(f"Covalent governance ansatz (hybridization: {hybridization})")
            else:
                # Fallback to UCC
                self.ansatz = UCCAnsatz(n_qubits=n_qubits, n_electrons=n_electrons)
                logger.warning(f"Governance ansatz not available for {bond_type}, using UCC")

        else:
            raise ValueError(f"Unknown ansatz type: {ansatz_type}")

        # Store ansatz type after creation
        self.ansatz_type = ansatz_type

    def _init_mapper(self):
        """Initialize fermionic-to-qubit mapper from bonds module."""
        from kanad.bonds import JordanWignerMapper, BravyiKitaevMapper

        # Use stored parameter
        mapper_type = self._mapper_type_param

        if mapper_type.lower() == 'jordan_wigner':
            self.mapper = JordanWignerMapper()
        elif mapper_type.lower() == 'bravyi_kitaev':
            self.mapper = BravyiKitaevMapper()
        else:
            # Default to Jordan-Wigner
            logger.warning(f"Unknown mapper type {mapper_type}, using Jordan-Wigner")
            self.mapper = JordanWignerMapper()

        # Store mapper type after creation
        self.mapper_type = mapper_type
        logger.debug(f"Mapper initialized: {self.mapper_type}")

    def _build_circuit(self):
        """Build parametrized quantum circuit."""
        self.circuit = self.ansatz.build_circuit()
        self.n_parameters = self.circuit.get_num_parameters()

        logger.info(f"Circuit built: {self.n_parameters} parameters")

        # Store pre-optimization stats
        if self.enable_optimization:
            try:
                self._gates_before_opt = self.circuit.count_ops().get('total', 0) if hasattr(self.circuit, 'count_ops') else None
                self._depth_before_opt = self.circuit.depth if isinstance(self.circuit.depth, int) else (self.circuit.depth() if callable(self.circuit.depth) else None)
            except (AttributeError, TypeError):
                self._gates_before_opt = None
                self._depth_before_opt = None

    def _init_backend(self, **kwargs):
        """Initialize quantum backend."""
        if self.backend == 'statevector':
            # Classical statevector simulation (exact, fast)
            self._use_statevector = True
            self._pauli_hamiltonian = None
            self._hamiltonian_matrix = None
            logger.info("Using statevector simulation (exact)")

        elif self.backend == 'bluequbit':
            # BlueQubit cloud backend
            try:
                from kanad.backends.bluequbit import BlueQubitBackend
                self._bluequbit_backend = BlueQubitBackend(**kwargs)
                self._use_statevector = False
                logger.info("BlueQubit backend initialized")
            except Exception as e:
                logger.error(f"BlueQubit initialization failed: {e}")
                raise

        elif self.backend == 'ibm':
            # IBM Quantum backend
            try:
                from kanad.backends.ibm import IBMRuntimeBackend
                self._ibm_backend = IBMRuntimeBackend(**kwargs)
                self._use_statevector = False
                logger.info("IBM Quantum backend initialized")
            except Exception as e:
                logger.error(f"IBM backend initialization failed: {e}")
                raise

        else:
            logger.warning(f"Unknown backend {self.backend}, using statevector")
            self._use_statevector = True

    def _compute_energy(self, parameters: np.ndarray) -> float:
        """
        Compute energy expectation value for given parameters.

        E(θ) = ⟨ψ(θ)|H|ψ(θ)⟩

        Args:
            parameters: Circuit parameters

        Returns:
            Energy expectation value (Hartree)
        """
        if self._use_statevector:
            # Classical statevector simulation
            return self._compute_energy_statevector(parameters)
        else:
            # Quantum backend
            return self._compute_energy_quantum(parameters)

    def _compute_energy_statevector(self, parameters: np.ndarray) -> float:
        """Compute energy using classical statevector simulation."""
        from qiskit.quantum_info import Statevector

        # Build circuit if not already built
        if self.ansatz.circuit is None:
            self.ansatz.build_circuit()

        # Bind parameters to circuit
        self.ansatz.circuit.bind_parameters(parameters)

        # Convert to Qiskit circuit and bind parameters
        qiskit_circuit = self.ansatz.circuit.to_qiskit()

        # If circuit still has parameters, bind them
        if qiskit_circuit.num_parameters > 0:
            param_dict = {qiskit_circuit.parameters[i]: parameters[i] for i in range(len(parameters))}
            bound_circuit = qiskit_circuit.assign_parameters(param_dict)
        else:
            bound_circuit = qiskit_circuit

        # Get statevector from circuit
        statevector = Statevector.from_instruction(bound_circuit)
        psi = statevector.data

        # Get Hamiltonian matrix if not cached
        if self._hamiltonian_matrix is None:
            # Get n_qubits from ansatz (more reliable than 2*n_orbitals)
            n_qubits = self.ansatz.n_qubits if hasattr(self.ansatz, 'n_qubits') else 2 * self.hamiltonian.n_orbitals

            # Check if to_matrix supports n_qubits parameter (for real hamiltonians)
            import inspect
            to_matrix_sig = inspect.signature(self.hamiltonian.to_matrix)
            has_n_qubits_param = 'n_qubits' in to_matrix_sig.parameters

            if has_n_qubits_param:
                # WORKAROUND: For hardware-efficient ansätze with fewer qubits than full space,
                # pad the statevector to match full Hamiltonian dimension
                # The Hamiltonian MUST be built with full 2*n_orbitals qubits for Jordan-Wigner
                full_n_qubits = 2 * self.hamiltonian.n_orbitals
                self._hamiltonian_matrix = self.hamiltonian.to_matrix(n_qubits=full_n_qubits, use_mo_basis=True)
                self._needs_padding = (n_qubits < full_n_qubits)
                self._ansatz_qubits = n_qubits
                self._full_qubits = full_n_qubits
            else:
                # Simple test Hamiltonian - just call to_matrix() without parameters
                H_core = self.hamiltonian.to_matrix()
                # Expand to full qubit space (2^n_qubits x 2^n_qubits)
                dim = 2 ** n_qubits
                self._hamiltonian_matrix = np.kron(H_core, np.eye(dim // H_core.shape[0]))
                self._needs_padding = False

        # Pad statevector if needed (for hardware-efficient ansätze with fewer qubits)
        if hasattr(self, '_needs_padding') and self._needs_padding:
            # Pad psi from 2^n to 2^m dimensions (n < m)
            # For H2: psi is 4-dim (2 qubits), need 16-dim (4 qubits)
            # We assume the ansatz prepares a state in the computational basis subspace
            psi_padded = np.zeros(2 ** self._full_qubits, dtype=complex)
            psi_padded[:len(psi)] = psi
            psi = psi_padded

        # Compute expectation value: E = <psi|H|psi>
        energy = np.real(np.conj(psi) @ self._hamiltonian_matrix @ psi)

        return float(energy)

    def _compute_energy_quantum(self, parameters: np.ndarray) -> float:
        """Compute energy using quantum backend (sampling-based)."""
        # This would use actual quantum backend
        # For now, placeholder
        logger.warning("Quantum backend energy computation not fully implemented")
        return self._compute_energy_statevector(parameters)

    def compute_energy(self, parameters: np.ndarray) -> float:
        """
        Public API: Compute energy expectation value for given parameters.

        Args:
            parameters: Circuit parameters

        Returns:
            Energy expectation value (Hartree)
        """
        return self._compute_energy(parameters)

    def get_energy_variance(self, parameters: np.ndarray) -> float:
        """
        Compute energy variance for given parameters.

        Variance = <H²> - <H>²

        Args:
            parameters: Circuit parameters

        Returns:
            Energy variance
        """
        from qiskit.quantum_info import Statevector

        # Build circuit if not already built
        if self.ansatz.circuit is None:
            self.ansatz.build_circuit()

        # Bind parameters to circuit
        self.ansatz.circuit.bind_parameters(parameters)

        # Convert to Qiskit circuit and bind parameters
        qiskit_circuit = self.ansatz.circuit.to_qiskit()

        # If circuit still has parameters, bind them
        if qiskit_circuit.num_parameters > 0:
            param_dict = {qiskit_circuit.parameters[i]: parameters[i] for i in range(len(parameters))}
            bound_circuit = qiskit_circuit.assign_parameters(param_dict)
        else:
            bound_circuit = qiskit_circuit

        # Get statevector from circuit
        statevector = Statevector.from_instruction(bound_circuit)
        psi = statevector.data

        # Get Hamiltonian matrix if not cached
        if self._hamiltonian_matrix is None:
            # Get n_qubits from ansatz (more reliable than 2*n_orbitals)
            n_qubits = self.ansatz.n_qubits if hasattr(self.ansatz, 'n_qubits') else 2 * self.hamiltonian.n_orbitals

            # Check if to_matrix supports n_qubits parameter (for real hamiltonians)
            import inspect
            to_matrix_sig = inspect.signature(self.hamiltonian.to_matrix)
            has_n_qubits_param = 'n_qubits' in to_matrix_sig.parameters

            if has_n_qubits_param:
                # WORKAROUND: For hardware-efficient ansätze with fewer qubits than full space,
                # pad the statevector to match full Hamiltonian dimension
                full_n_qubits = 2 * self.hamiltonian.n_orbitals
                self._hamiltonian_matrix = self.hamiltonian.to_matrix(n_qubits=full_n_qubits, use_mo_basis=True)
                self._needs_padding = (n_qubits < full_n_qubits)
                self._ansatz_qubits = n_qubits
                self._full_qubits = full_n_qubits
            else:
                # Simple test Hamiltonian - just call to_matrix() without parameters
                H_core = self.hamiltonian.to_matrix()
                # Expand to full qubit space (2^n_qubits x 2^n_qubits)
                dim = 2 ** n_qubits
                self._hamiltonian_matrix = np.kron(H_core, np.eye(dim // H_core.shape[0]))
                self._needs_padding = False

        # Pad statevector if needed (for hardware-efficient ansätze with fewer qubits)
        if hasattr(self, '_needs_padding') and self._needs_padding:
            psi_padded = np.zeros(2 ** self._full_qubits, dtype=complex)
            psi_padded[:len(psi)] = psi
            psi = psi_padded

        # Compute <H>
        H_expectation = np.real(np.conj(psi) @ self._hamiltonian_matrix @ psi)

        # Compute <H²>
        H_squared = self._hamiltonian_matrix @ self._hamiltonian_matrix
        H2_expectation = np.real(np.conj(psi) @ H_squared @ psi)

        # Variance = <H²> - <H>²
        variance = H2_expectation - H_expectation**2

        return float(variance)

    def _objective_function(self, parameters: np.ndarray) -> float:
        """
        Objective function for classical optimization.

        Tracks energy history and iteration count.

        Args:
            parameters: Current parameters

        Returns:
            Energy value to minimize
        """
        energy = self._compute_energy(parameters)

        # Track history
        self.energy_history.append(energy)
        self.parameter_history.append(parameters.copy())
        self.iteration_count += 1

        # Call user callback if provided
        if hasattr(self, '_callback') and self._callback is not None:
            self._callback(self.iteration_count, energy, parameters)

        # Log progress
        if self.iteration_count % 10 == 0:
            logger.info(f"Iteration {self.iteration_count}: E = {energy:.8f} Ha")

        return energy

    def solve(self, initial_parameters: Optional[np.ndarray] = None, callback: Optional[callable] = None) -> Dict[str, Any]:
        """
        Solve for ground state energy using VQE.

        Args:
            initial_parameters: Initial parameter guess (random if None)
            callback: Optional callback function(iteration, energy, params)

        Returns:
            Dictionary with comprehensive results:
                - energy: Ground state energy (Hartree)
                - parameters: Optimized parameters
                - converged: Convergence status
                - iterations: Number of iterations
                - hf_energy: Hartree-Fock reference energy
                - correlation_energy: E_VQE - E_HF
                - energy_history: Energy at each iteration
                - analysis: Detailed analysis (if enabled)
                - optimization_stats: Circuit optimization stats (if enabled)
        """
        # Store callback
        self._callback = callback
        logger.info("Starting VQE optimization...")

        # Get HF reference energy
        hf_energy = self.get_reference_energy()
        if hf_energy is not None:
            logger.info(f"HF reference energy: {hf_energy:.8f} Hartree")

        # Initial parameters
        if initial_parameters is None:
            # Start near zero (HF-like state)
            initial_parameters = np.random.randn(self.n_parameters) * 0.01

        logger.info(f"Optimizing {self.n_parameters} parameters using {self.optimizer_method}")

        # Reset tracking
        self.energy_history = []
        self.parameter_history = []
        self.iteration_count = 0

        # Classical optimization
        result = minimize(
            self._objective_function,
            initial_parameters,
            method=self.optimizer_method,
            options={'maxiter': self.max_iterations},
            tol=self.conv_threshold
        )

        # Store results
        self.results = {
            'energy': result.fun,
            'parameters': result.x,
            'converged': result.success,
            'iterations': self.iteration_count,
            'energy_history': np.array(self.energy_history),
            'parameter_history': np.array(self.parameter_history),
            'optimizer_message': result.message
        }

        # Add HF reference and correlation energy
        if hf_energy is not None:
            self.results['hf_energy'] = hf_energy
            self.results['correlation_energy'] = result.fun - hf_energy

            logger.info(f"VQE energy: {result.fun:.8f} Hartree")
            logger.info(f"Correlation energy: {result.fun - hf_energy:.8f} Hartree")

        # Add analysis if enabled
        if self.enable_analysis:
            # Use HF density matrix for analysis (VQE doesn't give us density matrix directly)
            density_matrix, _ = self.hamiltonian.solve_scf(max_iterations=50, conv_tol=1e-6)
            self._add_analysis_to_results(result.fun, density_matrix)

        # Add optimization stats if enabled
        if self.enable_optimization:
            self._add_optimization_stats()

        # Validate results
        validation = self.validate_results()
        self.results['validation'] = validation

        if not validation['passed']:
            logger.warning("VQE results failed validation checks!")

        logger.info(f"VQE optimization complete: {result.success}, {self.iteration_count} iterations")

        return self.results

    def print_summary(self):
        """Print comprehensive VQE results summary."""
        print("=" * 80)
        print("VQE SOLVER RESULTS")
        print("=" * 80)

        # System information
        print(f"\nMolecule: {'-'.join([a.symbol for a in self.atoms])}")
        print(f"Bond Type: {self.bond.bond_type}")
        print(f"Electrons: {self.molecule.n_electrons}")
        print(f"Orbitals: {self.hamiltonian.n_orbitals}")
        print(f"Qubits: {2 * self.hamiltonian.n_orbitals}")

        # Method details
        print(f"\nMethod: VQE")
        print(f"Ansatz: {self.ansatz_type.upper()}")
        print(f"Parameters: {self.n_parameters}")
        print(f"Mapper: {self.mapper_type}")
        print(f"Optimizer: {self.optimizer_method}")
        print(f"Backend: {self.backend}")

        # Energy results
        print("\n" + "-" * 80)
        print("ENERGY RESULTS")
        print("-" * 80)

        if 'energy' in self.results:
            print(f"\nVQE Energy:     {self.results['energy']:14.8f} Hartree")

        if 'hf_energy' in self.results:
            print(f"HF Reference:   {self.results['hf_energy']:14.8f} Hartree")

        if 'correlation_energy' in self.results:
            corr = self.results['correlation_energy']
            print(f"Correlation:    {corr:14.8f} Hartree ({corr * 627.509:.2f} kcal/mol)")

        # Convergence
        print("\n" + "-" * 80)
        print("CONVERGENCE")
        print("-" * 80)

        if 'converged' in self.results:
            status = "✓ Converged" if self.results['converged'] else "✗ Not Converged"
            print(f"\nStatus: {status}")

        if 'iterations' in self.results:
            print(f"Iterations: {self.results['iterations']}")

        if 'optimizer_message' in self.results:
            print(f"Message: {self.results['optimizer_message']}")

        # Energy history (convergence plot in text)
        if 'energy_history' in self.results and len(self.results['energy_history']) > 0:
            print("\nEnergy Convergence:")
            history = self.results['energy_history']
            # Show first, middle, last few points
            if len(history) > 10:
                indices = [0, 1, 2, len(history)//2, -3, -2, -1]
                for idx in indices:
                    if idx < 0:
                        idx = len(history) + idx
                    if 0 <= idx < len(history):
                        print(f"  Iter {idx:4d}: {history[idx]:14.8f} Ha")
            else:
                for i, E in enumerate(history):
                    print(f"  Iter {i:4d}: {E:14.8f} Ha")

        # Analysis results
        if 'analysis' in self.results and self.results['analysis']:
            print("\n" + "-" * 80)
            print("ANALYSIS")
            print("-" * 80)

            analysis = self.results['analysis']

            # Bonding analysis
            if analysis.get('bonding'):
                print("\nBonding Analysis:")
                bonding = analysis['bonding']
                for key, value in bonding.items():
                    if isinstance(value, (int, float)):
                        print(f"  {key:25s}: {value:.6f}")

            # Properties
            if analysis.get('properties'):
                print("\nMolecular Properties:")
                props = analysis['properties']
                for key, value in props.items():
                    if isinstance(value, (int, float)):
                        print(f"  {key:25s}: {value:.6f}")

        # Optimization stats
        if 'optimization_stats' in self.results:
            opt = self.results['optimization_stats']
            if any(opt.get('circuit', {}).values()):
                print("\n" + "-" * 80)
                print("CIRCUIT OPTIMIZATION")
                print("-" * 80)
                if opt['circuit'].get('gates_before') and opt['circuit'].get('gates_after'):
                    print(f"\nGates: {opt['circuit']['gates_before']} → {opt['circuit']['gates_after']}")
                    reduction = 100 * (1 - opt['circuit']['gates_after'] / opt['circuit']['gates_before'])
                    print(f"Reduction: {reduction:.1f}%")

        # Validation
        if 'validation' in self.results:
            validation = self.results['validation']
            if not validation['passed']:
                print("\n" + "-" * 80)
                print("⚠ VALIDATION WARNINGS")
                print("-" * 80)
                for check in validation['checks']:
                    if not check['passed']:
                        print(f"✗ {check['name']}: {check['message']}")
            else:
                print("\n✓ All validation checks passed")

        # Final summary
        print("\n" + "=" * 80)
        if self.results.get('converged'):
            print("✓ VQE OPTIMIZATION SUCCESSFUL")
        else:
            print("⚠ VQE OPTIMIZATION INCOMPLETE")
        print("=" * 80)
