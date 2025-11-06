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
from kanad.utils.hivqe_solver_mixin import HiVQESolverMixin

logger = logging.getLogger(__name__)


class VQESolver(HiVQESolverMixin, BaseSolver):
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
        # Hi-VQE parameters
        mode: str = 'standard',  # 'standard' or 'hivqe'
        hivqe_max_iterations: int = 10,  # Hi-VQE subspace expansion iterations
        hivqe_subspace_threshold: float = 0.05,  # Amplitude threshold for important configs
        use_active_space: bool = False,  # Enable governance-aware active space reduction
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

        # Hi-VQE parameters
        self.mode = mode.lower()  # 'standard' or 'hivqe'
        self.hivqe_max_iterations = hivqe_max_iterations
        self.hivqe_subspace_threshold = hivqe_subspace_threshold
        self.use_active_space = use_active_space

        if self.mode not in ['standard', 'hivqe']:
            raise ValueError(f"Invalid mode '{mode}'. Must be 'standard' or 'hivqe'")

        if self.mode == 'hivqe':
            logger.info(f"🔥 Hi-VQE mode enabled: {hivqe_max_iterations} iterations, threshold={hivqe_subspace_threshold}")
            if self.use_active_space:
                logger.info(f"   ✅ Active space reduction enabled")

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
        self.iteration_count = 0  # Function evaluation counter
        self.optimizer_iteration = 0  # Real optimizer iteration counter
        self._last_energy = None  # Track energy changes to detect optimizer iterations

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

        elif ansatz_type.lower() in ['governance', 'adaptive_governance']:
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
                # Get hybridization from governance metadata if available
                metadata = getattr(self.hamiltonian, '_governance_metadata', {})
                hybridization = metadata.get('hybridization', 'sp3')
                protocol = metadata.get('governance_protocol', self.hamiltonian.governance_protocol)

                # Choose between regular and adaptive governance
                if ansatz_type.lower() == 'adaptive_governance':
                    from kanad.ansatze import AdaptiveGovernanceOptimized

                    self.ansatz = AdaptiveGovernanceOptimized(
                        n_qubits=n_qubits,
                        n_electrons=n_electrons,
                        initial_layers=2,
                        max_layers=3,
                        hybridization=hybridization,
                        protocol=protocol,
                        mapper=self._mapper_type_param,
                        use_mp2_init=True  # Auto-enable MP2 init
                    )
                    logger.info(f"✨ Adaptive Governance ansatz (MP2 init enabled, {self._mapper_type_param} mapper)")
                else:
                    from kanad.ansatze.governance_aware_ansatz import CovalentGovernanceAnsatz

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

        # CRITICAL FIX: Filter ansatz excitations using governance protocol validation
        # This ensures VQE uses same excitations as SQD (consistency!)
        if hasattr(self.hamiltonian, 'governance_protocol') and self.hamiltonian.governance_protocol:
            protocol = self.hamiltonian.governance_protocol

            # Only filter if ansatz has excitations attribute (UCC, not hardware-efficient)
            if hasattr(self.ansatz, 'excitations') and hasattr(self.ansatz, 'get_excitation_list'):
                from kanad.core.configuration import Configuration

                original_excitations = self.ansatz.get_excitation_list()
                n_original = len(original_excitations)

                # HF reference bitstring
                hf_bitstring = '1' * n_electrons + '0' * (n_qubits - n_electrons)

                valid_excitations = []
                for exc in original_excitations:
                    # Convert excitation to configuration
                    occ, virt = exc

                    # Build bitstring by applying excitation to HF reference
                    bitlist = list(hf_bitstring)
                    for i in occ:
                        bitlist[i] = '0'  # Remove electron
                    for a in virt:
                        bitlist[a] = '1'  # Add electron

                    excited_bitstring = ''.join(bitlist)

                    # Check if valid according to governance
                    if protocol.is_valid_configuration(excited_bitstring):
                        valid_excitations.append(exc)

                # Update ansatz excitations
                self.ansatz.excitations = valid_excitations

                logger.info(f"✅ Governance filtering: {n_original} → {len(valid_excitations)} valid excitations")
                logger.info(f"   Protocol: {type(protocol).__name__}")
            else:
                logger.debug(f"Ansatz {ansatz_type} does not have excitations - skipping governance filtering")

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
        """Initialize quantum backend (uses base class implementation)."""
        super()._init_backend(self.backend, **kwargs)

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

        # DEBUG: Check statevector first time
        if not hasattr(self, '_statevector_checked'):
            print(f"🔍 Statevector: {len(statevector.data)} components")
            print(f"🔍 Bound circuit: {bound_circuit.num_qubits} qubits, depth {bound_circuit.depth()}")
            self._statevector_checked = True

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

        # Check if Hamiltonian is already a SparsePauliOp
        from qiskit.quantum_info import SparsePauliOp
        if isinstance(self.hamiltonian, SparsePauliOp):
            # Already sparse - just use it directly
            if self._sparse_pauli_op is None:
                self._sparse_pauli_op = self.hamiltonian
                self._cached_mapper = mapper_arg
                self._use_sparse = True
                print(f"📊 Using provided SparsePauliOp: {len(self._sparse_pauli_op)} Pauli terms")
        elif hasattr(self.hamiltonian, 'to_sparse_hamiltonian') and cache_invalid:
            # Build sparse Pauli operator (FAST, memory-efficient)
            logger.info(f"Building sparse Pauli Hamiltonian with {mapper_arg} mapper")
            print(f"🔧 Building sparse Pauli Hamiltonian with {mapper_arg} mapper")

            # Build Hamiltonian with correct mapper
            self._sparse_pauli_op = self.hamiltonian.to_sparse_hamiltonian(mapper=mapper_arg)
            self._cached_mapper = mapper_arg  # Update cache tracker
            self._use_sparse = True

            print(f"📊 Sparse Hamiltonian built: {len(self._sparse_pauli_op)} Pauli terms")

            # DEBUG: Check identity coefficient
            identity_coeff = 0.0
            for i, pauli_str in enumerate(self._sparse_pauli_op.paulis):
                if all(c == 'I' for c in str(pauli_str)):
                    identity_coeff += self._sparse_pauli_op.coeffs[i]
            print(f"🔍 VQE Hamiltonian identity coefficient: {identity_coeff:.8f} Ha")
            print(f"🔍 Nuclear repulsion from molecule: {self.hamiltonian.nuclear_repulsion:.8f} Ha")

            # Check qubit count consistency
            if self._sparse_pauli_op.num_qubits != n_qubits:
                print(f"⚠️  QUBIT MISMATCH: Hamiltonian={self._sparse_pauli_op.num_qubits}, Ansatz={n_qubits}")
                logger.warning(f"Qubit count mismatch: Hamiltonian={self._sparse_pauli_op.num_qubits}, Ansatz={n_qubits}")
                # Pad/truncate if needed
                if self._sparse_pauli_op.num_qubits > n_qubits:
                    self._needs_padding = True
                    self._ansatz_qubits = n_qubits
                    self._full_qubits = self._sparse_pauli_op.num_qubits
                    print(f"⚠️  PADDING ENABLED: Will pad {n_qubits}-qubit statevector to {self._full_qubits} qubits")
                    print(f"⚠️  THIS IS LIKELY THE BUG!")
                else:
                    self._needs_padding = False
            else:
                print(f"✅ Qubit counts match: {n_qubits} qubits")
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

                # DEBUG: Print first few energy evaluations
                if not hasattr(self, '_debug_counter'):
                    self._debug_counter = 0
                if self._debug_counter < 3:
                    print(f"🔍 Energy evaluation {self._debug_counter + 1}: {energy:.8f} Ha")
                    self._debug_counter += 1
        else:
            # FALLBACK: Dense matrix path (only for small test systems)
            logger.warning("Using SLOW dense matrix Hamiltonian - consider using sparse method")

            # Get Hamiltonian matrix if not cached
            if self._hamiltonian_matrix is None:
                # MEMORY SAFETY CHECK
                # Handle different Hamiltonian types
                from qiskit.quantum_info import SparsePauliOp
                if isinstance(self.hamiltonian, SparsePauliOp):
                    # SparsePauliOp passed directly (from tests)
                    # Get n_qubits from ansatz or Hamiltonian
                    if hasattr(self, 'circuit') and self.circuit:
                        full_n_qubits = self.circuit.num_qubits
                    elif hasattr(self.ansatz, 'n_qubits'):
                        full_n_qubits = self.ansatz.n_qubits
                    else:
                        # Get from Hamiltonian
                        full_n_qubits = self.hamiltonian.num_qubits
                elif hasattr(self.hamiltonian, 'n_orbitals'):
                    # MolecularHamiltonian
                    full_n_qubits = 2 * self.hamiltonian.n_orbitals
                else:
                    raise TypeError(
                        f"Unsupported Hamiltonian type: {type(self.hamiltonian)}. "
                        f"Expected MolecularHamiltonian or SparsePauliOp"
                    )

                required_memory_gb = (2 ** full_n_qubits) ** 2 * 16 / 1e9

                if required_memory_gb > 16:  # 16 GB limit
                    raise MemoryError(
                        f"Dense Hamiltonian matrix requires {required_memory_gb:.1f} GB RAM!\n"
                        f"System: {self.hamiltonian.n_orbitals} orbitals → {full_n_qubits} qubits\n"
                        f"Matrix size: {2**full_n_qubits} × {2**full_n_qubits} = {(2**full_n_qubits)**2:,} elements\n"
                        f"SOLUTION: Your Hamiltonian should implement to_sparse_hamiltonian() method"
                    )

                logger.info(f"Building dense Hamiltonian matrix ({required_memory_gb:.2f} GB)")

                # Handle different Hamiltonian types for matrix conversion
                from qiskit.quantum_info import SparsePauliOp
                if isinstance(self.hamiltonian, SparsePauliOp):
                    # SparsePauliOp can be converted to matrix directly
                    self._hamiltonian_matrix = self.hamiltonian.to_matrix()
                    self._needs_padding = False
                else:
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
        energy = self._compute_energy(parameters)

        # Track history
        self.energy_history.append(energy)
        self.parameter_history.append(parameters.copy())
        self.iteration_count += 1

        # Detect optimizer iterations: energy changes significantly or it's the first eval
        # This heuristic detects when optimizer completes an iteration and starts a new one
        energy_changed = (self._last_energy is None or
                         abs(energy - self._last_energy) > 1e-12)

        if energy_changed and self.iteration_count > 1:
            self.optimizer_iteration += 1
            self._last_energy = energy
        elif self.iteration_count == 1:
            self.optimizer_iteration = 1
            self._last_energy = energy

        # Call user callback if provided (used for progress broadcasting from API layer)
        # Pass: (optimizer_iteration, energy, parameters, function_eval_count)
        if hasattr(self, '_callback') and self._callback is not None:
            try:
                # Check if callback accepts 4 arguments (new signature)
                import inspect
                sig = inspect.signature(self._callback)
                if len(sig.parameters) >= 4:
                    # New signature: (iteration, energy, parameters, function_evals)
                    self._callback(self.optimizer_iteration, energy, parameters, self.iteration_count)
                else:
                    # Old signature: (iteration, energy, parameters) - pass optimizer iteration
                    self._callback(self.optimizer_iteration, energy, parameters)
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
            logger.info(f"Function eval {self.iteration_count} (~iter {est_iter}): E = {energy:.8f} Ha")
            print(f"📊 Progress: {self.iteration_count} function evals (~{est_iter} iterations), E = {energy:.8f} Ha")

        return energy

    def _spsa_minimize(self, initial_parameters: np.ndarray) -> tuple:
        """
        SPSA (Simultaneous Perturbation Stochastic Approximation) optimizer.

        Key advantage: Only 2 function evaluations per iteration regardless of parameter count!
        Perfect for cloud quantum backends where function evaluations are expensive.

        Args:
            initial_parameters: Initial parameter values

        Returns:
            (optimized_parameters, final_energy)
        """
        params = initial_parameters.copy()
        best_energy = float('inf')
        best_params = params.copy()

        # SPSA hyperparameters (standard values from Spall 1998)
        a = 0.16  # Step size coefficient
        c = 0.1   # Perturbation size coefficient
        A = 0.1 * self.max_iterations  # Stability constant
        alpha = 0.602  # Step size decay
        gamma = 0.101  # Perturbation decay

        prev_energy = None

        for k in range(self.max_iterations):
            # Compute decay schedules
            a_k = a / (k + 1 + A)**alpha
            c_k = c / (k + 1)**gamma

            # Random perturbation direction (Bernoulli ±1)
            delta = 2 * np.random.randint(0, 2, size=len(params)) - 1

            # Evaluate at perturbed points (only 2 evaluations!)
            params_plus = params + c_k * delta
            params_minus = params - c_k * delta

            energy_plus = self._objective_function(params_plus)
            energy_minus = self._objective_function(params_minus)

            # Gradient approximation
            gradient = (energy_plus - energy_minus) / (2 * c_k) * delta

            # Update parameters
            params = params - a_k * gradient

            # Track best result
            current_energy = min(energy_plus, energy_minus)
            if current_energy < best_energy:
                best_energy = current_energy
                best_params = params_plus if energy_plus < energy_minus else params_minus

            # Check convergence
            if prev_energy is not None:
                energy_change = abs(current_energy - prev_energy)
                if energy_change < self.conv_threshold:
                    print(f"✅ SPSA converged at iteration {k+1}")
                    break

            prev_energy = current_energy

            if (k + 1) % 5 == 0:
                print(f"  SPSA iter {k+1}/{self.max_iterations}: E = {best_energy:.8f} Ha")

        return best_params, best_energy

    def _compute_quantum_rdm1_from_statevector(self, statevector: np.ndarray) -> np.ndarray:
        """
        Compute 1-particle reduced density matrix (1-RDM) from VQE statevector.

        The statevector is |ψ⟩ = Σ_i c_i |i⟩ where |i⟩ are computational basis states.
        Each computational basis state is a Slater determinant.

        For a wavefunction |ψ⟩, the 1-RDM is:
            ρ_pq = ⟨ψ| a†_p a_q |ψ⟩ = Σ_ij c*_i c_j ⟨i| a†_p a_q |j⟩

        Args:
            statevector: Quantum state amplitudes (2^n_qubits,) in computational basis

        Returns:
            Quantum 1-RDM in spatial orbital basis (n_orbitals, n_orbitals)
        """
        n_orbitals = self.hamiltonian.n_orbitals
        n_qubits = 2 * n_orbitals
        hilbert_dim = 2 ** n_qubits

        # Ensure statevector is the right size
        if len(statevector) != hilbert_dim:
            logger.warning(f"Statevector size {len(statevector)} != expected {hilbert_dim}")
            # Try to handle this gracefully
            if len(statevector) < hilbert_dim:
                # Pad with zeros
                padded = np.zeros(hilbert_dim, dtype=complex)
                padded[:len(statevector)] = statevector
                statevector = padded
            else:
                # Truncate
                statevector = statevector[:hilbert_dim]

        # Initialize 1-RDM (spin-summed in spatial orbital basis)
        rdm1 = np.zeros((n_orbitals, n_orbitals), dtype=complex)

        # For each pair of spatial orbitals p, q
        for p in range(n_orbitals):
            for q in range(n_orbitals):
                rdm_element = 0.0 + 0.0j

                # Alpha spin contribution: a†_{p,α} a_{q,α}
                p_alpha = p
                q_alpha = q

                # Beta spin contribution: a†_{p,β} a_{q,β}
                p_beta = p + n_orbitals
                q_beta = q + n_orbitals

                # Sum over all pairs of computational basis states
                for i in range(hilbert_dim):
                    c_i = statevector[i]
                    if np.abs(c_i) < 1e-12:
                        continue  # Skip negligible amplitudes

                    for j in range(hilbert_dim):
                        c_j = statevector[j]
                        if np.abs(c_j) < 1e-12:
                            continue

                        # Alpha contribution
                        matrix_element_alpha = self._slater_condon_1body(i, j, p_alpha, q_alpha)
                        rdm_element += c_i.conj() * c_j * matrix_element_alpha

                        # Beta contribution
                        matrix_element_beta = self._slater_condon_1body(i, j, p_beta, q_beta)
                        rdm_element += c_i.conj() * c_j * matrix_element_beta

                rdm1[p, q] = rdm_element

        return rdm1.real

    def _slater_condon_1body(self, occ_I: int, occ_J: int, p: int, q: int) -> float:
        """
        Compute ⟨φ_I| a†_p a_q |φ_J⟩ using Slater-Condon rules.

        (Same implementation as SQD solver)

        Args:
            occ_I: Occupation bitstring for |φ_I⟩
            occ_J: Occupation bitstring for |φ_J⟩
            p: Creation orbital (spin-orbital index)
            q: Annihilation orbital (spin-orbital index)

        Returns:
            Matrix element (0, +1, or -1)
        """
        # Check if q is occupied in J and p is unoccupied in J
        q_occ_J = (occ_J >> q) & 1
        p_occ_J = (occ_J >> p) & 1

        if q == p:
            # Number operator: a†_p a_p = n_p
            if occ_I == occ_J:
                return float(q_occ_J)
            else:
                return 0.0

        # For p != q
        if q_occ_J == 0:  # Can't annihilate from unoccupied
            return 0.0

        # Apply a_q: remove electron from q
        occ_after_aq = occ_J ^ (1 << q)

        # Check if p is occupied after removing q
        p_occ_after_aq = (occ_after_aq >> p) & 1
        if p_occ_after_aq == 1:  # Can't create in occupied
            return 0.0

        # Apply a†_p: add electron to p
        occ_final = occ_after_aq ^ (1 << p)

        # Check if we got |φ_I⟩
        if occ_final != occ_I:
            return 0.0

        # Compute fermion sign
        sign = self._fermion_sign_1body(occ_J, p, q)
        return float(sign)

    def _fermion_sign_1body(self, occ: int, p: int, q: int) -> int:
        """
        Compute fermion sign for a†_p a_q acting on occupation occ.

        Args:
            occ: Initial occupation bitstring
            p: Creation index
            q: Annihilation index

        Returns:
            +1 or -1
        """
        if p == q:
            return 1

        # Count occupied orbitals between p and q (exclusive)
        min_idx = min(p, q)
        max_idx = max(p, q)

        count = 0
        for i in range(min_idx + 1, max_idx):
            if (occ >> i) & 1:
                count += 1

        return 1 if count % 2 == 0 else -1

    def solve(self, initial_parameters: Optional[np.ndarray] = None, callback: Optional[callable] = None) -> Dict[str, Any]:
        """
        Solve for ground state energy using VQE or Hi-VQE.

        Args:
            initial_parameters: Initial parameter guess (random if None)
            callback: Optional callback function(iteration, energy, params)

        Returns:
            Dictionary with comprehensive results:
                - energy: Ground state energy (Hartree)
                - parameters: Optimized parameters (None for Hi-VQE)
                - converged: Convergence status
                - iterations: Number of iterations
                - hf_energy: Hartree-Fock reference energy
                - correlation_energy: E_VQE - E_HF
                - energy_history: Energy at each iteration
                - analysis: Detailed analysis (if enabled)
                - optimization_stats: Circuit optimization stats (if enabled)
                - mode: 'standard' or 'hivqe'
                - hivqe_stats: Hi-VQE specific stats (if mode='hivqe')
        """
        # Store callback (only if explicitly provided, don't overwrite __init__ callback)
        if callback is not None:
            self._callback = callback

        # Route to appropriate solver based on mode
        if self.mode == 'hivqe':
            logger.info("🔥 Starting Hi-VQE optimization...")
            return self._solve_hivqe()
        else:
            logger.info("Starting standard VQE optimization...")
            return self._solve_standard_vqe(initial_parameters)

    def _solve_standard_vqe(self, initial_parameters: Optional[np.ndarray] = None) -> Dict[str, Any]:
        """Standard VQE solve with variational optimization."""
        # Get HF reference energy
        hf_energy = self.get_reference_energy()
        if hf_energy is not None:
            logger.info(f"HF reference energy: {hf_energy:.8f} Hartree")

        # Initial parameters
        if initial_parameters is None:
            # Check if ansatz supports smart initialization (Adaptive Governance)
            if hasattr(self.ansatz, 'get_smart_initial_params'):
                try:
                    logger.info(f"🎯 Using smart MP2-based initialization...")
                    initial_parameters = self.ansatz.get_smart_initial_params(hamiltonian=self.hamiltonian)
                    logger.info(f"✅ MP2 initialization successful (range: [{np.min(initial_parameters):.4f}, {np.max(initial_parameters):.4f}])")
                except Exception as e:
                    logger.warning(f"MP2 initialization failed: {e}, falling back to random")
                    initial_parameters = np.random.uniform(-0.1, 0.1, size=self.n_parameters)
            else:
                # Use uniform random initialization - provides good exploration
                # Range [-0.1, 0.1] is standard VQE practice
                initial_parameters = np.random.uniform(-0.1, 0.1, size=self.n_parameters)
                logger.info(f"Generated random initial parameters in range [-0.1, 0.1]")

        # Auto-select SPSA for cloud backends (huge cost savings)
        original_optimizer = self.optimizer_method
        if self.backend in ['ibm', 'bluequbit'] and self.optimizer_method not in ['SPSA', 'COBYLA', 'POWELL']:
            switch_msg = (
                f"☁️  CLOUD BACKEND OPTIMIZATION ☁️\n"
                f"   Backend: {self.backend}\n"
                f"   Original optimizer: {self.optimizer_method} (gradient-based)\n"
                f"   Auto-switching to: SPSA\n"
                f"   \n"
                f"   📉 Efficiency gain:\n"
                f"      {self.optimizer_method}: ~40 function evaluations/iteration\n"
                f"      SPSA: 2 function evaluations/iteration\n"
                f"      Expected speedup: 20x fewer quantum jobs\n"
                f"   \n"
                f"   💡 SPSA (Simultaneous Perturbation Stochastic Approximation) is\n"
                f"      specifically designed for noisy quantum hardware and cloud execution.\n"
                f"      It uses finite differences with random perturbations instead of gradients."
            )
            logger.warning(switch_msg)
            print(switch_msg)
            self.optimizer_method = 'SPSA'

            # Adjust max_iterations for SPSA (it typically needs fewer iterations)
            if self.max_iterations > 100:
                old_max_iter = self.max_iterations
                self.max_iterations = min(100, self.max_iterations)
                logger.info(f"📊 Adjusted max_iterations: {old_max_iter} → {self.max_iterations} (SPSA converges faster)")
                print(f"📊 Adjusted max_iterations: {old_max_iter} → {self.max_iterations}")

        logger.info(f"Optimizing {self.n_parameters} parameters using {self.optimizer_method}")

        # Reset tracking
        self.energy_history = []
        self.parameter_history = []
        self.iteration_count = 0  # Function evaluations
        self.optimizer_iteration = 0  # Real optimizer iterations
        self._last_energy = None  # For iteration detection

        # Classical optimization - Simple user-controlled approach
        # User sets max_iterations, optimizer uses natural convergence behavior
        # This is more predictable than trying to control function evaluations

        opt_options = {
            'maxiter': self.max_iterations,
            'disp': False  # Suppress optimizer output
        }

        logger.info(f"Optimizer: {self.optimizer_method}, max_iterations: {self.max_iterations}")
        print(f"🔧 Optimizer: {self.optimizer_method} with {self.max_iterations} iterations")

        # Use SPSA for cloud backends or if explicitly requested
        if self.optimizer_method == 'SPSA':
            print(f"📊 Using SPSA: 2 function evaluations per iteration (efficient for cloud)")
            final_params, final_energy = self._spsa_minimize(initial_parameters)

            # Create scipy-like result object for consistency
            class SPSAResult:
                def __init__(self, x, fun, nit, success=True, message="SPSA converged"):
                    self.x = x
                    self.fun = fun
                    self.nit = nit
                    self.success = success
                    self.message = message

            result = SPSAResult(
                x=final_params,
                fun=final_energy,
                nit=self.optimizer_iteration,
                success=True,
                message=f"SPSA completed {self.optimizer_iteration} iterations"
            )
        else:
            # Standard scipy optimizers
            result = minimize(
                self._objective_function,
                initial_parameters,
                method=self.optimizer_method,
                options=opt_options,
                tol=self.conv_threshold
            )

        # Store results
        self.results = {
            'energy': result.fun,
            'parameters': result.x,
            'converged': result.success,
            'iterations': result.nit if hasattr(result, 'nit') else self.iteration_count,  # Use optimizer iterations, not function evals
            'function_evaluations': self.iteration_count,  # Track function evaluations separately
            'energy_history': np.array(self.energy_history),
            'parameter_history': np.array(self.parameter_history),
            'optimizer_message': result.message,
            'mode': 'standard'  # Indicate this was standard VQE
        }

        # Add HF reference and correlation energy
        if hf_energy is not None:
            self.results['hf_energy'] = hf_energy
            self.results['correlation_energy'] = result.fun - hf_energy

            logger.info(f"VQE energy: {result.fun:.8f} Hartree")
            logger.info(f"Correlation energy: {result.fun - hf_energy:.8f} Hartree")

        # Add analysis if enabled
        if self.enable_analysis:
            logger.info(f"✅ Analysis enabled - generating analysis data...")
            try:
                # CRITICAL FIX: Extract quantum density from VQE statevector
                # Build circuit with optimal parameters
                if self._use_statevector:
                    from qiskit.quantum_info import Statevector

                    # Build circuit with optimal parameters (proper UCCAnsatz handling)
                    # Build circuit if not already built
                    if self.ansatz.circuit is None:
                        self.ansatz.build_circuit()

                    # Bind parameters to internal circuit
                    self.ansatz.circuit.bind_parameters(result.x)

                    # Convert to Qiskit circuit
                    qiskit_circuit = self.ansatz.circuit.to_qiskit()

                    # Bind parameters to Qiskit circuit
                    if qiskit_circuit.num_parameters > 0:
                        param_dict = {qiskit_circuit.parameters[i]: result.x[i] for i in range(len(result.x))}
                        bound_circuit = qiskit_circuit.assign_parameters(param_dict)
                    else:
                        bound_circuit = qiskit_circuit

                    # Get statevector
                    statevector = Statevector.from_instruction(bound_circuit)
                    statevector_array = statevector.data

                    logger.info("Computing quantum 1-RDM from VQE statevector...")
                    quantum_density = self._compute_quantum_rdm1_from_statevector(statevector_array)
                    logger.info(f"✅ Quantum density computed (includes correlation effects)")

                    # Store quantum density in results
                    self.results['quantum_rdm1'] = quantum_density

                    # CRITICAL FIX: Store quantum density in hamiltonian
                    if hasattr(self.hamiltonian, 'set_quantum_density_matrix'):
                        self.hamiltonian.set_quantum_density_matrix(quantum_density)

                    # Use quantum density for analysis
                    self._add_analysis_to_results(result.fun, quantum_density)
                    logger.info(f"✅ Analysis data added to results (using quantum density)")
                else:
                    # For non-statevector backends, fall back to HF for now
                    # TODO: Extract density from sampling results
                    logger.warning("Non-statevector backend: using HF density for analysis")
                    density_matrix, _ = self.hamiltonian.solve_scf(max_iterations=50, conv_tol=1e-6)
                    self._add_analysis_to_results(result.fun, density_matrix)
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

            # Try to get density matrix - prefer quantum, fallback to HF
            try:
                if 'quantum_rdm1' in self.results:
                    # Use quantum density if available (includes correlation)
                    self.results['rdm1'] = self.results['quantum_rdm1'].tolist()
                    logger.info(f"✅ Stored QUANTUM RDM1 for bonding analysis (correlated)")
                elif hasattr(self.hamiltonian, 'mf'):
                    # Fallback to HF density if quantum not available
                    if hasattr(self.hamiltonian.mf, 'make_rdm1'):
                        rdm1 = self.hamiltonian.mf.make_rdm1()
                        self.results['rdm1'] = rdm1.tolist()
                        logger.info(f"⚠️  Stored HF RDM1 (quantum density not computed)")
            except Exception as e:
                logger.warning(f"Could not extract RDM1: {e}")

            # Try to get orbital energies
            try:
                logger.debug(f"🔍 Checking hamiltonian for orbital energies: hasattr(mf)={hasattr(self.hamiltonian, 'mf')}")
                if hasattr(self.hamiltonian, 'mf'):
                    logger.debug(f"🔍 Hamiltonian has mf attribute, checking mo_energy: hasattr(mo_energy)={hasattr(self.hamiltonian.mf, 'mo_energy')}")
                    if hasattr(self.hamiltonian.mf, 'mo_energy'):
                        orb_energies = self.hamiltonian.mf.mo_energy
                        logger.debug(f"🔍 Found orbital energies: shape={orb_energies.shape}, dtype={orb_energies.dtype}")
                        self.results['orbital_energies'] = orb_energies.tolist()
                        logger.info(f"✅ Stored orbital energies for DOS analysis")
                    else:
                        logger.debug(f"Hamiltonian.mf does not have mo_energy attribute")
                else:
                    logger.debug(f"Hamiltonian does not have mf attribute (type: {type(self.hamiltonian).__name__}) - orbital energies not available")
            except Exception as e:
                logger.debug(f"Could not extract orbital energies: {e}")

            # Try to get dipole moment - prefer quantum density
            try:
                if hasattr(self.hamiltonian, 'mf'):
                    from pyscf import scf
                    if hasattr(scf, 'hf') and hasattr(scf.hf, 'dip_moment'):
                        # Use quantum density if available
                        if 'quantum_rdm1' in self.results:
                            dipole = scf.hf.dip_moment(self.hamiltonian.mf.mol, self.results['quantum_rdm1'])
                            self.results['dipole'] = dipole.tolist()
                            logger.info(f"✅ Stored dipole moment (using quantum density)")
                        else:
                            dipole = scf.hf.dip_moment(self.hamiltonian.mf.mol, self.hamiltonian.mf.make_rdm1())
                            self.results['dipole'] = dipole.tolist()
                            logger.info(f"⚠️  Stored dipole moment (using HF density)")
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

    def solve_with_restarts(self, n_restarts=3, callback=None):
        """
        Run VQE multiple times with different random initializations.
        Returns the best result among all attempts.

        This significantly improves reliability for stochastic optimization.

        Args:
            n_restarts: Number of VQE attempts (default: 3)
            callback: Optional callback function for progress tracking

        Returns:
            dict: Best VQE result with lowest energy
        """
        logger.info(f"🔄 Starting multi-start VQE with {n_restarts} restarts")

        best_energy = float('inf')
        best_result = None
        all_energies = []

        for attempt in range(1, n_restarts + 1):
            logger.info(f"🎯 VQE attempt {attempt}/{n_restarts}")

            # Run VQE with new random initialization
            result = self.solve(callback=callback)

            energy = result['energy']
            all_energies.append(energy)

            # Track best result
            if energy < best_energy:
                best_energy = energy
                best_result = result
                logger.info(f"   ✅ New best: {best_energy:.8f} Ha")
            else:
                logger.info(f"   Energy: {energy:.8f} Ha (not better than {best_energy:.8f})")

        # Add multi-start metadata
        best_result['multi_start'] = {
            'n_restarts': n_restarts,
            'all_energies': all_energies,
            'best_attempt': all_energies.index(best_energy) + 1,
            'energy_std': float(np.std(all_energies)),
            'energy_range': float(max(all_energies) - min(all_energies))
        }

        logger.info(f"🏆 Multi-start VQE complete:")
        logger.info(f"   Best energy: {best_energy:.8f} Ha (attempt {all_energies.index(best_energy) + 1})")
        logger.info(f"   Energy range: {max(all_energies) - min(all_energies):.8f} Ha")
        logger.info(f"   Energy std: {np.std(all_energies):.8f} Ha")

        return best_result

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
