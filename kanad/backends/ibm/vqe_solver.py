"""
VQE Solver for IBM Quantum Platform.

Implements Variational Quantum Eigensolver using IBM Runtime primitives
with optimized transpilation, error mitigation, and resilience settings.
"""

import numpy as np
from typing import Optional, Dict, List
import logging

from kanad.backends.ibm.backend import IBMRuntimeBackend
from kanad.core.hamiltonians.molecular_hamiltonian import MolecularHamiltonian
from kanad.ansatze.base_ansatz import BaseAnsatz

logger = logging.getLogger(__name__)


class IBMVQESolver:
    """
    VQE solver optimized for IBM Quantum hardware.

    Features:
    - Uses IBM Runtime Estimator primitive for efficient execution
    - Automatic transpilation with optimization levels
    - Error mitigation (resilience levels)
    - Session management for batch jobs
    - Supports both simulator and real hardware
    """

    def __init__(
        self,
        hamiltonian: MolecularHamiltonian,
        ansatz: BaseAnsatz,
        backend: Optional[IBMRuntimeBackend] = None,
        backend_name: Optional[str] = None,
        token: Optional[str] = None,
        instance: Optional[str] = None,
        optimizer: str = 'COBYLA',
        max_iterations: int = 100,
        shots: int = 4096,
        optimization_level: int = 3,
        resilience_level: int = 1,
        # Advanced error mitigation options
        enable_dynamical_decoupling: bool = True,
        dynamical_decoupling_sequence: str = 'XX',
        enable_twirling: bool = True,
        twirling_strategy: str = 'all',
        zne_extrapolator: Optional[str] = None,
        pec_mitigation: bool = False,
        **backend_options
    ):
        """
        Initialize IBM VQE solver with advanced error mitigation.

        Args:
            hamiltonian: Molecular Hamiltonian to minimize
            ansatz: Variational ansatz circuit
            backend: Pre-initialized IBM backend (optional)
            backend_name: IBM backend name (if backend not provided)
            token: IBM Quantum token
            instance: IBM Quantum instance (CRN)
            optimizer: Classical optimizer ('COBYLA', 'SLSQP', 'SPSA')
            max_iterations: Maximum optimization iterations
            shots: Number of measurement shots
            optimization_level: Transpiler optimization (0-3)
            resilience_level: Base error mitigation level (0-2)
                - 0: No mitigation
                - 1: Measurement error mitigation
                - 2: ZNE (Zero-Noise Extrapolation)

            Advanced Error Mitigation:
            enable_dynamical_decoupling: Insert DD sequences during idle times
            dynamical_decoupling_sequence: 'XX', 'XY4', or 'CPMG'
            enable_twirling: Randomized gate twirling for noise averaging
            twirling_strategy: 'active', 'passive', or 'all'
            zne_extrapolator: ZNE extrapolation method ('linear', 'exponential', etc.)
            pec_mitigation: Probabilistic Error Cancellation (experimental)
            **backend_options: Additional backend options
        """
        self.hamiltonian = hamiltonian
        self.ansatz = ansatz
        self.optimizer_name = optimizer
        self.max_iterations = max_iterations

        # Store error mitigation settings
        self.error_mitigation = {
            'resilience_level': resilience_level,
            'dynamical_decoupling': enable_dynamical_decoupling,
            'dd_sequence': dynamical_decoupling_sequence,
            'twirling': enable_twirling,
            'twirling_strategy': twirling_strategy,
            'zne_extrapolator': zne_extrapolator,
            'pec': pec_mitigation
        }

        # Initialize or use provided backend
        if backend is None:
            if backend_name is None:
                raise ValueError("Either 'backend' or 'backend_name' must be provided")

            logger.info(f"Initializing IBM backend: {backend_name}")
            self.backend = IBMRuntimeBackend(
                backend_name=backend_name,
                token=token,
                instance=instance,
                shots=shots,
                optimization_level=optimization_level,
                resilience_level=resilience_level,
                **backend_options
            )
        else:
            self.backend = backend
            logger.info(f"Using provided IBM backend")

        # Build circuit
        self.circuit = ansatz.build_circuit()
        self.n_parameters = self.circuit.get_num_parameters()

        # Convert Hamiltonian to sparse Pauli operators using PauliConverter
        logger.info("Converting Hamiltonian to sparse Pauli representation...")
        from kanad.core.hamiltonians.pauli_converter import PauliConverter
        from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper

        mapper = JordanWignerMapper()
        self.pauli_hamiltonian = PauliConverter.to_sparse_pauli_op(hamiltonian, mapper)
        logger.info(f"Hamiltonian: {len(self.pauli_hamiltonian)} Pauli terms")

        # Convert circuit to Qiskit
        self.qiskit_circuit = self.circuit.to_qiskit()
        logger.info(f"Circuit: {self.qiskit_circuit.num_qubits} qubits, "
                   f"{self.qiskit_circuit.depth()} depth, "
                   f"{self.n_parameters} parameters")

        # Optimization history
        self.energy_history: List[float] = []
        self.parameter_history: List[np.ndarray] = []
        self.iteration_count = 0

    def solve(
        self,
        initial_parameters: Optional[np.ndarray] = None,
        callback: Optional[callable] = None
    ) -> Dict:
        """
        Run VQE optimization on IBM Quantum.

        Args:
            initial_parameters: Initial parameter values (random if None)
            callback: Optional callback function(iteration, energy, parameters)

        Returns:
            Dictionary with:
                - energy: Optimized ground state energy (Hartree)
                - parameters: Optimized circuit parameters
                - iterations: Number of iterations
                - history: Energy convergence history
                - backend_info: IBM backend information
        """
        from scipy.optimize import minimize as scipy_minimize

        # Initialize parameters
        if initial_parameters is None:
            initial_parameters = np.random.uniform(
                -np.pi, np.pi, self.n_parameters
            )

        logger.info("=" * 80)
        logger.info("IBM QUANTUM VQE OPTIMIZATION")
        logger.info("=" * 80)

        backend_info = self.backend.get_backend_info()
        logger.info(f"Backend: {backend_info['name']}")
        logger.info(f"  Qubits: {backend_info.get('num_qubits', 'N/A')}")
        logger.info(f"  Type: {backend_info.get('backend_type', 'N/A')}")
        logger.info(f"  Shots: {self.backend.shots}")
        logger.info(f"  Optimization Level: {self.backend.optimization_level}")
        logger.info(f"  Resilience Level: {self.backend.resilience_level}")
        logger.info(f"Optimizer: {self.optimizer_name}")
        logger.info(f"Max Iterations: {self.max_iterations}")
        logger.info("=" * 80)

        # Get Estimator primitive with error mitigation options
        estimator_options = {}

        # Add twirling if enabled (Qiskit Runtime V2 feature)
        if self.error_mitigation['twirling']:
            estimator_options['twirling'] = {
                'enable_gates': self.error_mitigation['twirling_strategy'] in ['active', 'all'],
                'enable_measure': self.error_mitigation['twirling_strategy'] in ['passive', 'all'],
                'num_randomizations': 32  # Number of random gates per twirl
            }

        # Add ZNE extrapolator if specified (resilience_level 2 enables ZNE)
        if self.error_mitigation['zne_extrapolator'] and self.error_mitigation['resilience_level'] >= 2:
            estimator_options['zne'] = {
                'extrapolator': self.error_mitigation['zne_extrapolator'],
                'noise_factors': [1, 1.5, 2, 2.5, 3]  # Noise amplification factors
            }

        # Use batch mode (no session) for free/open plans
        # When premium subscription is available, set session=True for better batching
        estimator_options['session'] = False
        estimator = self.backend.get_estimator(**estimator_options)

        # Transpile circuit to ISA with error mitigation passes
        from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
        from qiskit.transpiler import PassManager

        logger.info("\nTranspiling circuit to ISA with error mitigation...")
        logger.info(f"Error Mitigation Settings:")
        logger.info(f"  Dynamical Decoupling: {self.error_mitigation['dynamical_decoupling']}")
        if self.error_mitigation['dynamical_decoupling']:
            logger.info(f"    Sequence: {self.error_mitigation['dd_sequence']}")
        logger.info(f"  Gate Twirling: {self.error_mitigation['twirling']}")
        if self.error_mitigation['twirling']:
            logger.info(f"    Strategy: {self.error_mitigation['twirling_strategy']}")

        # Generate base ISA transpilation
        pm = generate_preset_pass_manager(
            optimization_level=self.backend.optimization_level,
            backend=self.backend.backend
        )
        isa_circuit = pm.run(self.qiskit_circuit)

        # Add Dynamical Decoupling if enabled
        if self.error_mitigation['dynamical_decoupling']:
            from qiskit.transpiler.passes import PadDynamicalDecoupling
            from qiskit.circuit.library import XGate, YGate

            # Choose DD sequence
            dd_sequence = self.error_mitigation['dd_sequence']
            if dd_sequence == 'XX':
                dd_gates = [XGate(), XGate()]
            elif dd_sequence == 'XY4':
                dd_gates = [XGate(), YGate(), XGate(), YGate()]
            elif dd_sequence == 'CPMG':
                dd_gates = [XGate()] * 4  # Simplified CPMG
            else:
                dd_gates = [XGate(), XGate()]  # Default to XX

            try:
                dd_pm = PassManager([
                    PadDynamicalDecoupling(
                        durations=None,  # Use backend timing
                        dd_sequence=dd_gates
                    )
                ])
                isa_circuit = dd_pm.run(isa_circuit)
                logger.info(f"  ✓ Dynamical decoupling applied ({dd_sequence})")
            except Exception as e:
                logger.warning(f"  ⚠️  DD failed: {e}, continuing without DD")

        logger.info(f"ISA Circuit: {isa_circuit.depth()} depth, {isa_circuit.num_qubits} qubits")

        # Transform Hamiltonian to match ISA circuit qubits
        from qiskit.quantum_info import SparsePauliOp

        # Pad Hamiltonian with identity operators to match ISA circuit size
        n_extra_qubits = isa_circuit.num_qubits - len(self.pauli_hamiltonian.paulis[0])
        if n_extra_qubits > 0:
            logger.info(f"Padding Hamiltonian from {len(self.pauli_hamiltonian.paulis[0])} to {isa_circuit.num_qubits} qubits")
            padded_paulis = []
            for pauli, coeff in zip(self.pauli_hamiltonian.paulis, self.pauli_hamiltonian.coeffs):
                # Pad with identity on the right
                padded_pauli = str(pauli) + 'I' * n_extra_qubits
                padded_paulis.append((padded_pauli, coeff))
            isa_hamiltonian = SparsePauliOp.from_list(padded_paulis)
        else:
            isa_hamiltonian = self.pauli_hamiltonian

        # Define cost function
        def cost_function(parameters):
            """Compute energy expectation value."""
            self.iteration_count += 1

            # Bind parameters to ISA circuit
            bound_circuit = isa_circuit.assign_parameters(parameters)

            # Run estimation
            job = estimator.run([(bound_circuit, isa_hamiltonian)])
            result = job.result()

            # Extract energy (evs is now a scalar, not array)
            energy = result[0].data.evs

            # Store history
            self.energy_history.append(energy)
            self.parameter_history.append(parameters.copy())

            # Log progress
            logger.info(f"Iteration {self.iteration_count}: E = {energy:.6f} Ha")

            # Callback
            if callback:
                callback(self.iteration_count, energy, parameters)

            return energy

        # Run optimization
        logger.info("\nStarting optimization...")
        result = scipy_minimize(
            cost_function,
            initial_parameters,
            method=self.optimizer_name,
            options={'maxiter': self.max_iterations}
        )

        final_energy = result.fun
        final_parameters = result.x

        logger.info("=" * 80)
        logger.info("OPTIMIZATION COMPLETE")
        logger.info("=" * 80)
        logger.info(f"Final Energy: {final_energy:.6f} Ha")
        logger.info(f"             {final_energy * 27.211386:.4f} eV")
        logger.info(f"Iterations: {self.iteration_count}")
        logger.info(f"Converged: {result.success}")
        logger.info("=" * 80)

        return {
            'energy': final_energy,
            'energy_hartree': final_energy,
            'energy_ev': final_energy * 27.211386,
            'parameters': final_parameters,
            'iterations': self.iteration_count,
            'converged': result.success,
            'history': self.energy_history,
            'parameter_history': self.parameter_history,
            'backend_info': backend_info,
            'optimizer_result': result
        }
