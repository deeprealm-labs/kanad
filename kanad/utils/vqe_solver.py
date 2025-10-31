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
        molecule: Optional[Any] = None,  # Molecule for hamiltonian-based API
        # Common parameters
        optimizer: str = 'SLSQP',
        max_iterations: int = 100,  # Default, but should be set from frontend
        conv_threshold: float = 1e-6,
        backend: str = 'statevector',
        shots: Optional[int] = None,
        enable_analysis: bool = True,
        enable_optimization: bool = True,
        experiment_id: Optional[str] = None,  # For WebSocket broadcasting and cancellation
        job_id: Optional[str] = None,  # For cancellation checking
        callback: Optional[Callable] = None,  # Progress callback
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
            if ansatz is None and mapper is None:
                # Hamiltonian-based API with types (for polyatomic molecules from API)
                # Use provided molecule or extract from hamiltonian
                if molecule is not None:
                    molecule_obj = molecule
                elif hasattr(hamiltonian, 'molecule'):
                    molecule_obj = hamiltonian.molecule
                else:
                    # Hamiltonian might not have molecule reference, try to get it
                    molecule_obj = getattr(hamiltonian, '_molecule', None)

                if molecule_obj is None:
                    raise ValueError("Must provide 'molecule' parameter or hamiltonian with molecule reference for type-based initialization")

                # Initialize with molecule and hamiltonian
                self.molecule = molecule_obj
                self.hamiltonian = hamiltonian
                self.bond = None
                self._enable_analysis = enable_analysis
                self._enable_optimization = enable_optimization
                self.enable_analysis = enable_analysis  # Public attribute for _build_circuit
                self.enable_optimization = enable_optimization  # Public attribute for _build_circuit
                self._api_mode = 'hamiltonian_types'
            else:
                # Low-level API: Initialize from components (for testing)
                self._init_from_components_mode(hamiltonian, ansatz, mapper, molecule, enable_analysis, enable_optimization)
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

        # Store experiment_id and job_id for WebSocket broadcasting and cancellation
        self.experiment_id = experiment_id
        self.job_id = job_id

        # Track cloud provider job information
        self.cloud_job_ids = []  # List of cloud job IDs from IBM/BlueQubit
        self.cloud_provider = None  # 'ibm' or 'bluequbit'
        self.execution_mode = None  # 'batch', 'session', 'instance' for IBM

        # Store callback (don't pass to backend via kwargs)
        self._callback = callback

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
        elif self._api_mode == 'hamiltonian_types':
            # Hamiltonian-based API with types - initialize like bond mode
            self._init_ansatz()
            self._init_mapper()
            # Build quantum circuit
            self._build_circuit()
            # Initialize backend
            self._init_backend(**kwargs)
        else:
            # Components mode - ansatz and mapper already set
            # Initialize backend (this will set _use_statevector correctly)
            self._hamiltonian_matrix = None
            self._init_backend(**kwargs)
            logger.info(f"VQE initialized in components mode, backend={self.backend}, use_statevector={self._use_statevector}")

        # Optimization tracking
        self.energy_history = []
        self.parameter_history = []
        self.iteration_count = 0

        # Performance optimization: Cache sparse Pauli operator
        # CRITICAL: Cache must be invalidated if mapper changes!
        self._sparse_pauli_op = None
        self._cached_mapper = None  # Track which mapper was used for cache
        self._use_sparse = False

        logger.info(f"VQE Solver initialized: {self.ansatz_type} ansatz, {self.mapper_type} mapping, {self.backend_name} backend")

    def _init_from_components_mode(self, hamiltonian, ansatz, mapper, molecule, enable_analysis, enable_optimization):
        """Initialize VQE from individual components (low-level API for testing)."""
        # Store components directly
        self.hamiltonian = hamiltonian
        self.mapper = mapper if mapper is not None else JordanWignerMapper()

        # Try to get molecule from parameter, then hamiltonian, then None
        if molecule is not None:
            self.molecule = molecule
        else:
            self.molecule = getattr(hamiltonian, 'molecule', None)

        self.bond = None

        # Enable analysis features if molecule is available
        self.enable_analysis = enable_analysis if self.molecule is not None else False
        self.enable_optimization = enable_optimization

        # Initialize analysis tools if enabled and molecule available
        if self.enable_analysis and self.molecule is not None:
            try:
                from kanad.analysis import EnergyAnalyzer, BondingAnalyzer, PropertyCalculator
                self.energy_analyzer = EnergyAnalyzer(self.hamiltonian)
                self.bonding_analyzer = BondingAnalyzer(self.molecule.hamiltonian)
                self.property_calculator = PropertyCalculator(self.hamiltonian)
                self.atoms = self.molecule.atoms
                logger.info("✅ Analysis tools initialized in components mode")
            except Exception as e:
                logger.warning(f"❌ Failed to initialize analysis tools: {e}")
                import traceback
                traceback.print_exc()
                self.enable_analysis = False

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
                entanglement='linear',
                mapper=self._mapper_type_param
            )
            logger.info(f"Hardware-efficient ansatz: 3 layers, linear entanglement, {self._mapper_type_param} mapper")

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
                    protocol=protocol,
                    mapper=self._mapper_type_param
                )
                logger.info(f"Covalent governance ansatz (hybridization: {hybridization}, {self._mapper_type_param} mapper)")
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
        logger.info(f"Initializing backend: {self.backend}")
        print(f"🔧 Initializing backend: {self.backend}")

        if self.backend == 'statevector':
            # Classical statevector simulation (exact, fast)
            self._use_statevector = True
            self._pauli_hamiltonian = None
            self._hamiltonian_matrix = None
            logger.info("Using statevector simulation (exact)")
            print(f"📍 Using statevector simulation")

        elif self.backend == 'bluequbit':
            # BlueQubit cloud backend
            try:
                print(f"🌐 Initializing BlueQubit backend with kwargs: {list(kwargs.keys())}")
                from kanad.backends.bluequbit import BlueQubitBackend
                self._bluequbit_backend = BlueQubitBackend(**kwargs)
                self._use_statevector = False
                device = kwargs.get('device', 'gpu')
                logger.info(f"BlueQubit backend initialized: device={device}")
                print(f"✅ Connected to BlueQubit cloud: device={device}")
                print(f"🔗 Track your jobs at: https://app.bluequbit.io/jobs")
            except Exception as e:
                logger.error(f"BlueQubit initialization failed: {e}")
                print(f"❌ BlueQubit initialization failed: {e}")
                raise

        elif self.backend == 'ibm':
            # IBM Quantum backend
            try:
                print(f"🌐 Initializing IBM backend with kwargs: {list(kwargs.keys())}")
                from kanad.backends.ibm import IBMBackend
                self._ibm_backend = IBMBackend(**kwargs)
                self._use_statevector = False
                backend_name = kwargs.get('backend_name', 'ibm_torino')
                logger.info(f"IBM Quantum backend initialized: {backend_name}")
                print(f"✅ Connected to IBM Quantum: {backend_name}")
                print(f"🔗 Track your jobs at: https://quantum.ibm.com/jobs")
            except Exception as e:
                logger.error(f"IBM backend initialization failed: {e}")
                print(f"❌ IBM backend initialization failed: {e}")
                raise

        else:
            logger.warning(f"Unknown backend {self.backend}, using statevector")
            print(f"⚠️ Unknown backend {self.backend}, falling back to statevector")
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
        """
        Compute energy using classical statevector simulation.

        PERFORMANCE: Uses sparse Pauli operators for efficiency (100-1000x faster).
        MEMORY SAFE: No dense matrix construction for large molecules.
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

        # Get n_qubits from ansatz
        n_qubits = self.ansatz.n_qubits if hasattr(self.ansatz, 'n_qubits') else 2 * self.hamiltonian.n_orbitals

        # CRITICAL PERFORMANCE IMPROVEMENT: Use sparse Pauli operators instead of dense matrices
        # Check if Hamiltonian has sparse method (covalent, ionic, molecular hamiltonians)
        # CRITICAL BUG FIX: Must check if cached Hamiltonian matches current mapper!

        # Determine mapper type to pass to Hamiltonian
        mapper_name = getattr(self, 'mapper_type', 'jordan_wigner')
        if mapper_name.lower() in ['bravyikitaevmapper', 'bravyi_kitaev']:
            mapper_arg = 'bravyi_kitaev'
        else:
            mapper_arg = 'jordan_wigner'

        # Check if we need to rebuild Hamiltonian (cache miss or mapper changed)
        cache_invalid = (self._sparse_pauli_op is None or
                        self._cached_mapper != mapper_arg)

        if hasattr(self.hamiltonian, 'to_sparse_hamiltonian') and cache_invalid:
            # Build sparse Pauli operator (FAST, memory-efficient)
            logger.info(f"Building sparse Pauli Hamiltonian with {mapper_arg} mapper")

            # Build Hamiltonian with correct mapper
            self._sparse_pauli_op = self.hamiltonian.to_sparse_hamiltonian(mapper=mapper_arg)
            self._cached_mapper = mapper_arg  # Update cache tracker
            self._use_sparse = True

            # Check qubit count consistency
            if self._sparse_pauli_op.num_qubits != n_qubits:
                logger.warning(f"Qubit count mismatch: Hamiltonian={self._sparse_pauli_op.num_qubits}, Ansatz={n_qubits}")
                # Pad/truncate if needed
                if self._sparse_pauli_op.num_qubits > n_qubits:
                    self._needs_padding = True
                    self._ansatz_qubits = n_qubits
                    self._full_qubits = self._sparse_pauli_op.num_qubits
                else:
                    self._needs_padding = False
            else:
                self._needs_padding = False

        # Compute energy using sparse or dense method
        if hasattr(self, '_use_sparse') and self._use_sparse:
            # FAST PATH: Sparse Pauli operator
            # Pad statevector if needed
            if self._needs_padding:
                psi = statevector.data
                psi_padded = np.zeros(2 ** self._full_qubits, dtype=complex)
                psi_padded[:len(psi)] = psi
                statevector_padded = Statevector(psi_padded)
                energy = statevector_padded.expectation_value(self._sparse_pauli_op).real
            else:
                # Direct expectation value computation (FAST!)
                energy = statevector.expectation_value(self._sparse_pauli_op).real
        else:
            # FALLBACK: Dense matrix path (only for small test systems)
            logger.warning("Using SLOW dense matrix Hamiltonian - consider using sparse method")

            # Get Hamiltonian matrix if not cached
            if self._hamiltonian_matrix is None:
                # MEMORY SAFETY CHECK
                full_n_qubits = 2 * self.hamiltonian.n_orbitals
                required_memory_gb = (2 ** full_n_qubits) ** 2 * 16 / 1e9

                if required_memory_gb > 16:  # 16 GB limit
                    raise MemoryError(
                        f"Dense Hamiltonian matrix requires {required_memory_gb:.1f} GB RAM!\n"
                        f"System: {self.hamiltonian.n_orbitals} orbitals → {full_n_qubits} qubits\n"
                        f"Matrix size: {2**full_n_qubits} × {2**full_n_qubits} = {(2**full_n_qubits)**2:,} elements\n"
                        f"SOLUTION: Your Hamiltonian should implement to_sparse_hamiltonian() method"
                    )

                logger.info(f"Building dense Hamiltonian matrix ({required_memory_gb:.2f} GB)")

                # Check if to_matrix supports n_qubits parameter
                import inspect
                to_matrix_sig = inspect.signature(self.hamiltonian.to_matrix)
                has_n_qubits_param = 'n_qubits' in to_matrix_sig.parameters

                if has_n_qubits_param:
                    self._hamiltonian_matrix = self.hamiltonian.to_matrix(n_qubits=full_n_qubits, use_mo_basis=True)
                    self._needs_padding = (n_qubits < full_n_qubits)
                    self._ansatz_qubits = n_qubits
                    self._full_qubits = full_n_qubits
                else:
                    # Simple test Hamiltonian
                    H_core = self.hamiltonian.to_matrix()
                    dim = 2 ** n_qubits
                    self._hamiltonian_matrix = np.kron(H_core, np.eye(dim // H_core.shape[0]))
                    self._needs_padding = False

            # Pad statevector if needed
            psi = statevector.data
            if hasattr(self, '_needs_padding') and self._needs_padding:
                psi_padded = np.zeros(2 ** self._full_qubits, dtype=complex)
                psi_padded[:len(psi)] = psi
                psi = psi_padded

            # Compute expectation value: E = <psi|H|psi>
            energy = np.real(np.conj(psi) @ self._hamiltonian_matrix @ psi)

        return float(energy)

    def _compute_energy_from_counts(self, counts: dict, pauli_hamiltonian) -> float:
        """
        Compute energy expectation value from measurement counts.

        IMPORTANT: This is a simplified implementation that only accurately measures
        Z-basis Pauli terms. For X and Y terms, proper basis rotation is needed.

        For molecular Hamiltonians with Jordan-Wigner mapping, the majority of terms
        are Z-basis, so this provides a reasonable approximation. For production use,
        consider using statevector mode (shots=None) or implementing full Pauli measurements.

        For a Hamiltonian H = Σ c_i P_i where P_i are Pauli strings,
        we compute <H> = Σ c_i <P_i> where <P_i> is estimated from counts.

        Args:
            counts: Measurement counts dict (bitstring -> count)
            pauli_hamiltonian: SparsePauliOp representing the Hamiltonian

        Returns:
            Estimated energy expectation value (approximate for non-Z terms)
        """
        # Get total number of shots
        total_shots = sum(counts.values())

        # Convert counts to probabilities
        probabilities = {bitstring: count / total_shots for bitstring, count in counts.items()}

        # Energy is weighted sum of Pauli term expectations
        energy = 0.0

        # Separate Z-only terms from terms with X/Y (for diagnostic purposes)
        z_terms_weight = 0.0
        xy_terms_weight = 0.0

        # Iterate over Pauli terms in Hamiltonian
        for pauli_str, coeff in zip(pauli_hamiltonian.paulis, pauli_hamiltonian.coeffs):
            # Compute expectation value of this Pauli term from counts
            expectation = 0.0

            pauli_string = str(pauli_str)
            has_xy = any(c in pauli_string for c in ['X', 'Y'])

            if has_xy:
                xy_terms_weight += abs(float(np.real(coeff)))
            else:
                z_terms_weight += abs(float(np.real(coeff)))

            for bitstring, prob in probabilities.items():
                # Compute eigenvalue of Pauli string for this bitstring
                # Pauli operators have eigenvalues ±1
                eigenvalue = 1.0

                # Convert bitstring to list of bits (reverse order for Qiskit convention)
                bits = [int(b) for b in bitstring[::-1]]

                # For each qubit, check if Pauli operator is Z
                # Z|0> = +|0>, Z|1> = -|1>
                for i, pauli_char in enumerate(pauli_string[::-1]):  # Reverse to match qubit ordering
                    if pauli_char == 'Z' and i < len(bits):
                        if bits[i] == 1:
                            eigenvalue *= -1
                    elif pauli_char in ['X', 'Y'] and i < len(bits):
                        # X and Y terms need basis rotation for accurate measurement
                        # For now, we treat them as averaging to 0 contribution
                        # This is inaccurate but prevents complete nonsense results
                        eigenvalue = 0.0
                        break

                expectation += prob * eigenvalue

            # Add weighted contribution to energy
            energy += float(np.real(coeff)) * expectation

        # Log diagnostic info
        if xy_terms_weight > 0.01:
            logger.warning(
                f"Hamiltonian has significant X/Y terms (weight: {xy_terms_weight:.4f} vs Z: {z_terms_weight:.4f}). "
                f"Counts-based estimation may be inaccurate. Consider using shots=None for statevector mode."
            )

        return energy

    def _compute_energy_quantum(self, parameters: np.ndarray) -> float:
        """
        Compute energy using quantum backend (sampling-based).

        Supports IBM Quantum and BlueQubit backends.
        """
        from qiskit.quantum_info import SparsePauliOp

        # Build circuit if not already built
        if self.ansatz.circuit is None:
            self.ansatz.build_circuit()

        # Bind parameters to circuit
        self.ansatz.circuit.bind_parameters(parameters)

        # Convert to Qiskit circuit
        qiskit_circuit = self.ansatz.circuit.to_qiskit()

        # Bind parameters if circuit has them
        if qiskit_circuit.num_parameters > 0:
            param_dict = {qiskit_circuit.parameters[i]: parameters[i]
                         for i in range(len(parameters))}
            bound_circuit = qiskit_circuit.assign_parameters(param_dict)
        else:
            bound_circuit = qiskit_circuit

        # Get Pauli representation of Hamiltonian
        try:
            from kanad.core.hamiltonians.pauli_converter import PauliConverter
            pauli_hamiltonian = PauliConverter.to_sparse_pauli_op(
                self.hamiltonian,
                self.mapper,
                use_qiskit_nature=True
            )
        except Exception as e:
            logger.error(f"Failed to convert Hamiltonian to Pauli operators: {e}")
            logger.warning("Falling back to statevector simulation")
            import traceback
            traceback.print_exc()
            return self._compute_energy_statevector(parameters)

        # Use IBM backend if available
        if hasattr(self, '_ibm_backend') and self._ibm_backend is not None:
            logger.info(f"Submitting to IBM Quantum (iteration {self.iteration_count})")

            # Broadcast to frontend if experiment_id is available
            if self.experiment_id:
                try:
                    from api.utils import broadcast_log_sync
                    broadcast_log_sync(self.experiment_id, f"🚀 Submitting job to IBM Quantum (function eval {self.iteration_count})")
                except Exception:
                    print(f"🚀 Submitting job to IBM Quantum (function eval {self.iteration_count})")
            else:
                print(f"🚀 Submitting job to IBM Quantum (function eval {self.iteration_count})")

            try:
                # Submit job to IBM
                result = self._ibm_backend.run_batch(
                    circuits=[bound_circuit],
                    observables=[pauli_hamiltonian],
                    shots=self.shots
                )

                job_id = result['job_id']
                logger.info(f"IBM job submitted: {job_id}")

                # Track cloud job information
                self.cloud_job_ids.append(job_id)
                self.cloud_provider = 'ibm'
                self.execution_mode = 'batch'  # IBM uses batch mode

                # Broadcast job ID to frontend
                if self.experiment_id:
                    try:
                        from api.utils import broadcast_log_sync
                        broadcast_log_sync(self.experiment_id, f"✅ IBM job submitted: {job_id}")
                        broadcast_log_sync(self.experiment_id, f"🔗 Track job at: https://quantum.ibm.com/jobs/{job_id}")
                    except Exception:
                        print(f"✅ IBM job submitted: {job_id}")
                        print(f"🔗 Track job at: https://quantum.ibm.com/jobs/{job_id}")
                else:
                    print(f"✅ IBM job submitted: {job_id}")
                    print(f"🔗 Track job at: https://quantum.ibm.com/jobs/{job_id}")

                # Wait for job to complete with cancellation checking
                job = self._ibm_backend.service.job(job_id)
                logger.info(f"Waiting for IBM job {job_id}...")

                if self.experiment_id:
                    try:
                        from api.utils import broadcast_log_sync
                        broadcast_log_sync(self.experiment_id, f"⏳ Waiting for IBM job to complete...")
                    except Exception:
                        print(f"⏳ Waiting for IBM job to complete...")
                else:
                    print(f"⏳ Waiting for IBM job to complete...")

                # Poll job status with cancellation checks (instead of blocking job.result())
                import time
                poll_interval = 5  # Check every 5 seconds
                job_result = None

                while True:
                    # Check for cancellation FIRST before checking job status
                    if self.experiment_id and self.job_id:
                        try:
                            from api.services.experiment_service import check_cancellation
                            check_cancellation(self.experiment_id, self.job_id)
                        except Exception as cancel_error:
                            # Cancellation detected - cancel IBM job and raise
                            logger.warning(f"Cancellation detected, cancelling IBM job {job_id}...")
                            try:
                                from api.utils import broadcast_log_sync
                                broadcast_log_sync(self.experiment_id, f"🚫 Cancelling IBM job {job_id}...")
                            except:
                                pass

                            try:
                                self._ibm_backend.cancel_job(job_id)
                                logger.info(f"IBM job {job_id} cancelled successfully")
                            except Exception as e:
                                logger.error(f"Failed to cancel IBM job {job_id}: {e}")

                            # Re-raise the cancellation exception
                            raise cancel_error

                    # Check job status
                    status = self._ibm_backend.get_job_status(job_id)
                    logger.debug(f"IBM job {job_id} status: {status}")

                    # Check if job is done (handles different status formats)
                    status_upper = str(status).upper()
                    if 'DONE' in status_upper or 'COMPLETED' in status_upper:
                        logger.info(f"IBM job {job_id} completed, retrieving result...")
                        job_result = job.result()
                        break
                    elif 'ERROR' in status_upper or 'CANCELLED' in status_upper or 'FAILED' in status_upper:
                        logger.error(f"IBM job {job_id} failed with status: {status}")
                        raise RuntimeError(f"IBM job failed: {status}")

                    # Job still running, wait before next poll
                    time.sleep(poll_interval)

                # Extract energy from Estimator result
                # Handle both V1 and V2 primitive formats
                if hasattr(job_result, 'values'):
                    # V1 Estimator result format
                    energy = float(job_result.values[0])
                elif hasattr(job_result, '__getitem__'):
                    # V2 Estimator result format (EstimatorV2)
                    # Result is a PrimitiveResult object with PubResult items
                    # Each PubResult has .data.evs (expectation values)
                    energy = float(job_result[0].data.evs[0])
                else:
                    raise AttributeError(f"Unknown Estimator result format: {type(job_result)}")

                logger.info(f"IBM job {job_id} completed: E = {energy:.8f} Ha")
                return energy

            except Exception as e:
                # Check if this is a cancellation exception - if so, re-raise it
                if 'ExperimentCancelledException' in type(e).__name__ or 'cancelled' in str(e).lower():
                    logger.warning(f"🚫 Cancellation detected during IBM backend execution - stopping")
                    raise  # Re-raise cancellation exception to stop optimizer

                # For other errors, fall back to statevector
                logger.error(f"IBM backend execution failed: {e}")
                logger.warning("Falling back to statevector simulation")
                return self._compute_energy_statevector(parameters)

        # Use BlueQubit backend if available
        elif hasattr(self, '_bluequbit_backend') and self._bluequbit_backend is not None:
            logger.info(f"Submitting to BlueQubit (iteration {self.iteration_count})")

            # Broadcast to frontend if experiment_id is available
            if self.experiment_id:
                try:
                    from api.utils import broadcast_log_sync
                    broadcast_log_sync(self.experiment_id, f"🚀 Submitting job to BlueQubit (function eval {self.iteration_count})")
                except Exception:
                    print(f"🚀 Submitting job to BlueQubit (function eval {self.iteration_count})")
            else:
                print(f"🚀 Submitting job to BlueQubit (function eval {self.iteration_count})")

            try:
                # For BlueQubit CPU/MPS devices, use statevector mode for accuracy
                # (counts-based measurement requires basis rotation which is complex)
                use_shots = self.shots
                if hasattr(self._bluequbit_backend, 'device'):
                    device = self._bluequbit_backend.device
                    if device in ['cpu', 'mps.cpu', 'mps.gpu']:
                        use_shots = None  # Force statevector mode
                        logger.info(f"Using statevector mode for {device} (more accurate than sampling)")

                # Submit to BlueQubit (asynchronous to allow cancellation)
                job_info = self._bluequbit_backend.run_circuit(
                    circuit=bound_circuit,
                    shots=use_shots,  # None = statevector mode
                    asynchronous=True  # Changed to async for cancellation support
                )

                job_id = job_info['job_id']
                logger.info(f"BlueQubit job submitted: {job_id}")

                # Track cloud job information
                self.cloud_job_ids.append(job_id)
                self.cloud_provider = 'bluequbit'

                # Broadcast job submission
                if self.experiment_id:
                    try:
                        from api.utils import broadcast_log_sync
                        broadcast_log_sync(self.experiment_id, f"✅ BlueQubit job submitted: {job_id}")
                        broadcast_log_sync(self.experiment_id, f"🔗 Track job at: https://app.bluequbit.io/jobs/{job_id}")
                        broadcast_log_sync(self.experiment_id, f"⏳ Waiting for BlueQubit job to complete...")
                    except Exception:
                        print(f"✅ BlueQubit job submitted: {job_id}")
                        print(f"🔗 Track job at: https://app.bluequbit.io/jobs/{job_id}")
                else:
                    print(f"✅ BlueQubit job submitted: {job_id}")
                    print(f"🔗 Track job at: https://app.bluequbit.io/jobs/{job_id}")

                # Check for cancellation before waiting for job
                # (Note: BlueQubit's wait_for_job is blocking, so we check before calling it)
                if self.experiment_id and self.job_id:
                    try:
                        from api.services.experiment_service import check_cancellation
                        check_cancellation(self.experiment_id, self.job_id)
                    except Exception as cancel_error:
                        # Cancellation detected - cancel BlueQubit job and raise
                        logger.warning(f"Cancellation detected, cancelling BlueQubit job {job_id}...")
                        try:
                            from api.utils import broadcast_log_sync
                            broadcast_log_sync(self.experiment_id, f"🚫 Cancelling BlueQubit job {job_id}...")
                        except:
                            pass

                        try:
                            self._bluequbit_backend.cancel_job(job_id)
                            logger.info(f"BlueQubit job {job_id} cancelled successfully")
                        except Exception as e:
                            logger.error(f"Failed to cancel BlueQubit job {job_id}: {e}")

                        # Re-raise the cancellation exception
                        raise cancel_error

                # Wait for job to complete (blocking call)
                result = self._bluequbit_backend.wait_for_job(job_id)
                logger.info(f"BlueQubit job {job_id} completed successfully")

                # Extract statevector or counts
                if 'statevector' in result:
                    statevector = np.array(result['statevector'])

                    # Convert Pauli operator to matrix
                    pauli_matrix = pauli_hamiltonian.to_matrix()

                    # Compute expectation value: <ψ|H|ψ>
                    energy = float(np.real(np.conj(statevector) @ pauli_matrix @ statevector))

                elif 'counts' in result:
                    # Sampling-based energy estimation
                    counts = result['counts']
                    logger.info(f"Computing energy from counts: {sum(counts.values())} total shots")

                    # Compute expectation value from Pauli measurements
                    energy = self._compute_energy_from_counts(
                        counts=counts,
                        pauli_hamiltonian=pauli_hamiltonian
                    )
                    logger.info(f"Counts-based energy: {energy:.8f} Ha")
                else:
                    raise ValueError("BlueQubit result missing statevector or counts")

                logger.info(f"BlueQubit energy: {energy:.8f} Ha")
                return energy

            except Exception as e:
                # Check if this is a cancellation exception - if so, re-raise it
                if 'ExperimentCancelledException' in type(e).__name__ or 'cancelled' in str(e).lower():
                    logger.warning(f"🚫 Cancellation detected during BlueQubit backend execution - stopping")
                    raise  # Re-raise cancellation exception to stop optimizer

                # For other errors, fall back to statevector
                logger.error(f"BlueQubit backend execution failed: {e}")
                logger.warning("Falling back to statevector simulation")
                return self._compute_energy_statevector(parameters)

        else:
            logger.warning("No quantum backend available, using statevector")
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
        electronic_energy = self._compute_energy(parameters)

        # Get nuclear repulsion for total energy (for broadcasting and display)
        nuclear_repulsion = getattr(self.hamiltonian, 'nuclear_repulsion', 0.0)
        total_energy = electronic_energy + nuclear_repulsion

        # Track history (store total energy for consistency)
        self.energy_history.append(total_energy)
        self.parameter_history.append(parameters.copy())
        self.iteration_count += 1

        # Call user callback if provided (used for progress broadcasting from API layer)
        # IMPORTANT: Pass total energy (with nuclear repulsion) to callback for correct display
        if hasattr(self, '_callback') and self._callback is not None:
            try:
                self._callback(self.iteration_count, total_energy, parameters)
            except Exception as e:
                # Check if this is a cancellation exception - if so, re-raise it to stop optimizer
                if 'ExperimentCancelledException' in type(e).__name__ or 'cancelled' in str(e).lower():
                    logger.warning(f"🚫 Experiment cancelled - stopping optimizer")
                    raise  # Re-raise to stop optimizer
                else:
                    # For other callback errors, log but don't stop the optimizer
                    print(f"⚠️ Callback failed: {e}")
                    import traceback
                    traceback.print_exc()
        elif self.iteration_count == 1:
            # Debug: print on first iteration if callback is missing
            print(f"🔍 DEBUG: _callback exists={hasattr(self, '_callback')}, is None={getattr(self, '_callback', 'MISSING') is None}")

        # Log progress (every 10 function evals)
        if self.iteration_count % 10 == 0:
            # Estimate optimizer iteration (rough approximation)
            est_iter = self.iteration_count // 40 if self.optimizer_method in ['SLSQP', 'L-BFGS-B'] else self.iteration_count
            logger.info(f"Function eval {self.iteration_count} (~iter {est_iter}): E = {total_energy:.8f} Ha")
            print(f"📊 Progress: {self.iteration_count} function evals (~{est_iter} iterations), E = {total_energy:.8f} Ha")

        # Return electronic energy (no nuclear repulsion) for optimizer
        # The optimizer should minimize ONLY the electronic part
        return electronic_energy

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
        # Store callback (only if explicitly provided, don't overwrite __init__ callback)
        if callback is not None:
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

        # Warn if using cloud backend with gradient-based optimizer
        if self.backend in ['ibm', 'bluequbit'] and self.optimizer_method in ['SLSQP', 'L-BFGS-B']:
            estimated_jobs = self.max_iterations * 40
            warning_msg = (
                f"⚠️  OPTIMIZER WARNING ⚠️\n"
                f"   Optimizer: {self.optimizer_method} (gradient-based)\n"
                f"   Max iterations: {self.max_iterations}\n"
                f"   Estimated quantum jobs: ~{estimated_jobs} function evaluations\n"
                f"   \n"
                f"   SLSQP uses ~40 function evaluations per iteration.\n"
                f"   Consider using COBYLA or POWELL for fewer job submissions (~1-3x iterations)."
            )
            logger.warning(warning_msg)
            print(warning_msg)

        # Classical optimization
        # Prepare optimizer options
        opt_options = {'maxiter': self.max_iterations}

        # Set maxfun (max function evaluations) for ALL optimizers
        # This is critical - scipy defaults are HUGE (SLSQP=15000, POWELL=unlimited!)
        # User sets max_iterations, we calculate reasonable maxfun based on optimizer

        # Normalize optimizer name to handle case variations
        optimizer_upper = self.optimizer_method.upper()

        # Derivative-free optimizers (low cost per iteration)
        if optimizer_upper == 'COBYLA':
            # COBYLA: ~1-3 function evals per iteration
            # Uses 'maxfun' option
            min_maxfun = self.n_parameters + 2
            maxfun = max(self.max_iterations * 3, min_maxfun)
            opt_options['maxfun'] = maxfun
            logger.info(f"COBYLA: maxiter={self.max_iterations}, maxfun={maxfun} (~3x)")

        elif optimizer_upper == 'POWELL':
            # Powell's method: ~2-5 function evals per iteration (uses line searches)
            # Uses 'maxfev' option (NOT maxfun!)
            maxfev = self.max_iterations * 5
            opt_options['maxfev'] = maxfev
            logger.info(f"Powell: maxiter={self.max_iterations}, maxfev={maxfev} (~5x)")

        elif optimizer_upper == 'NELDER-MEAD':
            # Nelder-Mead: simplex method, ~5-10 evals per iteration
            # Uses 'maxfev' option
            maxfev = self.max_iterations * 10
            opt_options['maxfev'] = maxfev
            logger.info(f"Nelder-Mead: maxiter={self.max_iterations}, maxfev={maxfev} (~10x)")

        # Gradient-based optimizers (HIGH cost per iteration - need finite differences)
        elif optimizer_upper == 'SLSQP':
            # SLSQP: ~(2*n_params) to 100 function evals per iteration
            # Uses 'maxiter' only - no separate function eval limit in scipy
            # Estimate: 2*n_params for gradient + some for line search
            logger.info(f"⚠️ SLSQP: maxiter={self.max_iterations}, ~{self.n_parameters * 2} function evals per iter - HIGH COST!")

        elif optimizer_upper == 'L-BFGS-B':
            # L-BFGS-B: ~(2*n_params) to 100 function evals per iteration
            # Uses 'maxfun' option
            maxfun = self.max_iterations * (self.n_parameters * 2 + 10)
            opt_options['maxfun'] = maxfun
            logger.info(f"⚠️ L-BFGS-B: maxiter={self.max_iterations}, maxfun={maxfun} (~{self.n_parameters * 2}x per iter) - HIGH COST!")

        elif optimizer_upper == 'BFGS':
            # BFGS: ~(2*n_params) function evals per iteration
            # Uses 'maxiter' only - no maxfun option
            logger.info(f"⚠️ BFGS: maxiter={self.max_iterations}, ~{self.n_parameters * 2} function evals per iter - HIGH COST!")

        elif optimizer_upper == 'CG':
            # Conjugate Gradient: ~(2*n_params) function evals per iteration
            # Uses 'maxiter' only - no maxfun option
            logger.info(f"⚠️ CG: maxiter={self.max_iterations}, ~{self.n_parameters * 2} function evals per iter - HIGH COST!")

        elif optimizer_upper == 'TNC':
            # TNC (Truncated Newton): ~(2*n_params) to 50 function evals per iteration
            # TNC uses 'maxfun' but doesn't accept 'maxiter' - remove it
            del opt_options['maxiter']
            maxfun = self.max_iterations * (self.n_parameters * 2 + 10)
            opt_options['maxfun'] = maxfun
            logger.info(f"⚠️ TNC: maxfun={maxfun} (~{self.n_parameters * 2}x per iter) - HIGH COST!")

        else:
            # Other optimizers: try setting maxfun as fallback
            logger.warning(f"⚠️ Optimizer '{self.optimizer_method}' may not respect iteration limits")

        result = minimize(
            self._objective_function,
            initial_parameters,
            method=self.optimizer_method,
            options=opt_options,
            tol=self.conv_threshold
        )

        # Get nuclear repulsion energy (constant energy shift)
        nuclear_repulsion = getattr(self.hamiltonian, 'nuclear_repulsion', 0.0)

        # Total energy = electronic energy + nuclear repulsion
        total_energy = result.fun + nuclear_repulsion

        # Store results
        self.results = {
            'energy': total_energy,  # FIXED: Now includes nuclear repulsion!
            'electronic_energy': result.fun,  # Raw electronic energy from VQE
            'parameters': result.x,
            'converged': result.success,
            'iterations': result.nit if hasattr(result, 'nit') else self.iteration_count,  # Use optimizer iterations, not function evals
            'function_evaluations': self.iteration_count,  # Track function evaluations separately
            'energy_history': np.array(self.energy_history),
            'parameter_history': np.array(self.parameter_history),
            'optimizer_message': result.message
        }

        # Add HF reference and correlation energy
        if hf_energy is not None:
            self.results['hf_energy'] = hf_energy
            self.results['correlation_energy'] = total_energy - hf_energy

            logger.info(f"VQE energy (total): {total_energy:.8f} Hartree")
            logger.info(f"Electronic energy: {result.fun:.8f} Hartree")
            logger.info(f"Nuclear repulsion: {nuclear_repulsion:.8f} Hartree")
            logger.info(f"Correlation energy: {total_energy - hf_energy:.8f} Hartree")

        # Add analysis if enabled
        if self.enable_analysis:
            logger.info(f"✅ Analysis enabled - generating analysis data...")
            try:
                # Use HF density matrix for analysis (VQE doesn't give us density matrix directly)
                density_matrix, _ = self.hamiltonian.solve_scf(max_iterations=50, conv_tol=1e-6)
                self._add_analysis_to_results(result.fun, density_matrix)
                logger.info(f"✅ Analysis data added to results")
            except Exception as e:
                logger.error(f"❌ Analysis generation failed: {e}")
                import traceback
                traceback.print_exc()
        else:
            logger.warning(f"⚠️  Analysis disabled - enable_analysis={self.enable_analysis}, molecule={self.molecule is not None}")

        # Add optimization stats if enabled
        if self.enable_optimization:
            self._add_optimization_stats()

        # Validate results
        validation = self.validate_results()
        self.results['validation'] = validation

        if not validation['passed']:
            logger.warning("VQE results failed validation checks!")

        # ADD ENHANCED DATA FOR ANALYSIS SERVICE
        try:
            # Store molecule geometry for ADME and other analyses
            if self.molecule is not None:
                self.results['geometry'] = [
                    (atom.symbol, tuple(atom.position))
                    for atom in self.molecule.atoms
                ]
                self.results['atoms'] = [atom.symbol for atom in self.molecule.atoms]
                self.results['n_atoms'] = self.molecule.n_atoms
                self.results['n_electrons'] = self.molecule.n_electrons
                self.results['charge'] = getattr(self.molecule, 'charge', 0)
                self.results['multiplicity'] = getattr(self.molecule, 'multiplicity', 1)
                logger.info(f"✅ Stored molecule geometry for analysis")

            # Store nuclear repulsion energy
            if hasattr(self.hamiltonian, 'nuclear_repulsion'):
                self.results['nuclear_repulsion'] = float(self.hamiltonian.nuclear_repulsion)

            # Try to get density matrix from HF reference
            try:
                if hasattr(self.hamiltonian, 'mf'):
                    # PySCF mean-field object has density matrix
                    if hasattr(self.hamiltonian.mf, 'make_rdm1'):
                        rdm1 = self.hamiltonian.mf.make_rdm1()
                        self.results['rdm1'] = rdm1.tolist()
                        logger.info(f"✅ Stored RDM1 for bonding analysis")
            except Exception as e:
                logger.warning(f"Could not extract RDM1: {e}")

            # Try to get orbital energies
            try:
                logger.info(f"🔍 Checking hamiltonian for orbital energies: hasattr(mf)={hasattr(self.hamiltonian, 'mf')}")
                if hasattr(self.hamiltonian, 'mf'):
                    logger.info(f"🔍 Hamiltonian has mf attribute, checking mo_energy: hasattr(mo_energy)={hasattr(self.hamiltonian.mf, 'mo_energy')}")
                    if hasattr(self.hamiltonian.mf, 'mo_energy'):
                        orb_energies = self.hamiltonian.mf.mo_energy
                        logger.info(f"🔍 Found orbital energies: shape={orb_energies.shape}, dtype={orb_energies.dtype}")
                        self.results['orbital_energies'] = orb_energies.tolist()
                        logger.info(f"✅ Stored orbital energies for DOS analysis")
                    else:
                        logger.warning(f"⚠️  Hamiltonian.mf does not have mo_energy attribute")
                else:
                    logger.warning(f"⚠️  Hamiltonian does not have mf attribute (type: {type(self.hamiltonian).__name__})")
            except Exception as e:
                logger.warning(f"Could not extract orbital energies: {e}")

            # Try to get dipole moment
            try:
                if hasattr(self.hamiltonian, 'mf'):
                    from pyscf import scf
                    if hasattr(scf, 'hf') and hasattr(scf.hf, 'dip_moment'):
                        dipole = scf.hf.dip_moment(self.hamiltonian.mf.mol, self.hamiltonian.mf.make_rdm1())
                        self.results['dipole'] = dipole.tolist()
                        logger.info(f"✅ Stored dipole moment")
            except Exception as e:
                logger.warning(f"Could not calculate dipole: {e}")

        except Exception as e:
            logger.error(f"Error storing enhanced data: {e}")

        # Add cloud provider job information to results
        if self.cloud_provider:
            self.results['cloud_provider'] = self.cloud_provider
            if self.cloud_job_ids:
                self.results['cloud_job_ids'] = self.cloud_job_ids
                # Generate job URLs
                if self.cloud_provider == 'bluequbit':
                    # BlueQubit job URLs
                    self.results['cloud_job_urls'] = [
                        f"https://app.bluequbit.io/jobs/{jid}" for jid in self.cloud_job_ids
                    ]
                elif self.cloud_provider == 'ibm':
                    # IBM Quantum job URLs
                    self.results['cloud_job_urls'] = [
                        f"https://quantum.ibm.com/jobs/{jid}" for jid in self.cloud_job_ids
                    ]
            if self.execution_mode:
                self.results['execution_mode'] = self.execution_mode
            logger.info(f"✅ Added cloud job info: provider={self.cloud_provider}, jobs={len(self.cloud_job_ids)}")

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
