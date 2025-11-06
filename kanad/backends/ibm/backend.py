"""
IBM Quantum Backend Implementation

Provides interface to IBM Quantum hardware and cloud simulators.
"""

import os
import logging
from typing import Dict, Any, Optional, List, Union
import numpy as np

from kanad.backends.ibm.error_mitigation import ErrorMitigationStrategy

logger = logging.getLogger(__name__)


class IBMBackend:
    """
    IBM Quantum backend for quantum hardware and simulators.

    Supports:
    - Real quantum hardware (127+ qubits)
    - Cloud simulators
    - Batch mode (parallel independent jobs)
    - Session mode (reserved hardware for iterative algorithms)
    - Qiskit Runtime primitives (SamplerV2, EstimatorV2)

    Usage:
        # Batch mode (for independent jobs)
        backend = IBMBackend(backend_name='ibm_brisbane')
        results = backend.run_batch(circuits, observables)

        # Session mode (for Hi-VQE and iterative algorithms)
        results = backend.run_session(circuits, observables, max_time='1h')
    """

    def __init__(
        self,
        backend_name: Optional[str] = None,
        api_token: Optional[str] = None,
        channel: str = 'ibm_quantum_platform',
        crn: Optional[str] = None,
        instance: Optional[str] = None
    ):
        """
        Initialize IBM Quantum backend.

        Args:
            backend_name: Backend name (e.g., 'ibm_brisbane', 'ibmq_qasm_simulator')
                         If None, uses least busy backend
            api_token: IBM Quantum API token (or set IBM_API env var)
            channel: 'ibm_quantum_platform' or 'ibm_cloud'
            crn: Cloud Resource Name (required for ibm_cloud)
            instance: IBM Quantum instance (optional)
        """
        self.backend_name = backend_name
        self.channel = channel
        self.crn = crn or os.getenv('IBM_CRN')
        self.instance = instance

        # Get API token
        self.api_token = api_token or os.getenv('IBM_API')

        if not self.api_token:
            raise ValueError(
                "IBM Quantum API token required. Set IBM_API environment variable "
                "or pass api_token parameter. Get token from https://quantum.ibm.com"
            )

        if channel == 'ibm_cloud' and not self.crn:
            raise ValueError(
                "IBM_CRN (Cloud Resource Name) required for ibm_cloud channel"
            )

        # Initialize service and backend
        self._init_service()
        self._init_backend()

        logger.info(f"IBM backend initialized: {self.backend.name}")

    def _init_service(self):
        """Initialize Qiskit Runtime service."""
        try:
            from qiskit_ibm_runtime import QiskitRuntimeService

            # Save account if not already saved
            try:
                if self.channel == 'ibm_cloud':
                    QiskitRuntimeService.save_account(
                        channel=self.channel,
                        token=self.api_token,
                        instance=self.crn,
                        overwrite=True
                    )
                else:
                    QiskitRuntimeService.save_account(
                        channel=self.channel,
                        token=self.api_token,
                        instance=self.instance,
                        overwrite=True
                    )
            except:
                pass  # Account already exists

            # Initialize service
            self.service = QiskitRuntimeService(channel=self.channel)

            logger.info(f"IBM Quantum service initialized ({self.channel})")

        except ImportError:
            raise ImportError(
                "qiskit-ibm-runtime required. Install with: pip install qiskit-ibm-runtime"
            )
        except Exception as e:
            raise RuntimeError(f"Failed to initialize IBM Quantum service: {e}")

    def _init_backend(self):
        """Initialize quantum backend."""
        if self.backend_name:
            self.backend = self.service.backend(self.backend_name)
        else:
            # Get least busy backend
            self.backend = self.service.least_busy(operational=True, simulator=False)
            self.backend_name = self.backend.name

        logger.info(f"Using backend: {self.backend.name}")
        logger.info(f"  Qubits: {self.backend.num_qubits}")
        logger.info(f"  Quantum: {not self.backend.simulator}")

    def run_session(
        self,
        circuits: Union[List, 'QuantumCircuit'],
        observables: Optional[List] = None,
        shots: int = 1024,
        optimization_level: int = 1,
        resilience_level: int = 1,
        max_time: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        Run circuits in session mode (reserved hardware access).

        Session mode is ideal for Hi-VQE and other iterative algorithms that
        require sequential job execution with priority queue access.

        Args:
            circuits: Quantum circuit(s) to run
            observables: Observable(s) for Estimator (optional)
            shots: Number of measurement shots
            optimization_level: Transpilation optimization (0-3)
            resilience_level: Error mitigation level (0-2)
            max_time: Maximum session time (e.g., '1h', '30m', '2h')
                     If None, uses default session timeout

        Returns:
            Results dictionary with job ID and session context
        """
        from qiskit_ibm_runtime import Session, SamplerV2 as Sampler, EstimatorV2 as Estimator
        from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

        # Ensure circuits is a list
        if not isinstance(circuits, list):
            circuits = [circuits]

        logger.info(f"Running session on {self.backend.name}")
        logger.info(f"  Circuits: {len(circuits)}")
        logger.info(f"  Shots: {shots}")
        logger.info(f"  Max time: {max_time or 'default'}")

        # Transpile circuits for target hardware
        logger.info(f"  Transpiling circuits (optimization_level={optimization_level})...")
        pm = generate_preset_pass_manager(
            backend=self.backend,
            optimization_level=optimization_level
        )
        transpiled_circuits = pm.run(circuits)
        logger.info(f"  Transpilation complete")

        # Expand observables to match transpiled circuit size if needed
        if observables is not None:
            transpiled_observables = []
            for i, (circuit, observable) in enumerate(zip(transpiled_circuits, observables)):
                original_qubits = circuits[i].num_qubits
                transpiled_qubits = circuit.num_qubits

                if transpiled_qubits > original_qubits:
                    # Pad observable with identity operators
                    padding = "I" * (transpiled_qubits - original_qubits)
                    padded_terms = []
                    for pauli_str, coeff in observable.to_list():
                        padded_pauli = pauli_str + padding
                        padded_terms.append((padded_pauli, coeff))

                    from qiskit.quantum_info import SparsePauliOp
                    padded_observable = SparsePauliOp.from_list(padded_terms)
                    transpiled_observables.append(padded_observable)
                else:
                    transpiled_observables.append(observable)

            observables = transpiled_observables

        try:
            # Build session parameters
            session_params = {'backend': self.backend}
            if max_time:
                session_params['max_time'] = max_time

            with Session(**session_params) as session:
                if observables is not None:
                    # Use Estimator for energy expectation values
                    estimator = Estimator(session=session)

                    # Set options (V2 primitives)
                    estimator.options.default_shots = shots

                    # CRITICAL FIX (Issue #8): Use auto_configure() for error mitigation
                    error_mitigation = ErrorMitigationStrategy.auto_configure(self.backend.name)

                    # Apply auto-configured error mitigation strategy
                    estimator.options.resilience_level = error_mitigation.resilience_level
                    resilience_options = error_mitigation.get_resilience_options()

                    # Set resilience options if any
                    if resilience_options:
                        for key, value in resilience_options.items():
                            setattr(estimator.options.resilience, key, value)
                        logger.info(f"  ✅ Auto-configured error mitigation applied")
                        logger.info(f"     Resilience level: {error_mitigation.resilience_level}")
                        if error_mitigation.zne_extrapolation:
                            logger.info(f"     ZNE: {error_mitigation.zne_extrapolation}")
                        if error_mitigation.dynamical_decoupling:
                            logger.info(f"     DD: {error_mitigation.dynamical_decoupling}")

                    logger.info("Using Estimator primitive (session mode)")

                    # Build pub (circuit, observable) tuples with transpiled circuits
                    pubs = [(transpiled_circuits[i], observables[i]) for i in range(len(transpiled_circuits))]

                    job = estimator.run(pubs)

                    # Return job and session info
                    return {
                        'job_id': job.job_id(),
                        'status': str(job.status()),
                        'backend': self.backend.name,
                        'session_id': session.session_id,
                        'mode': 'session'
                    }

                else:
                    # Use Sampler for measurement counts
                    sampler = Sampler(session=session)

                    # Set options
                    sampler.options.default_shots = shots

                    logger.info("Using Sampler primitive (session mode)")

                    job = sampler.run(transpiled_circuits)

                    # Return job and session info
                    return {
                        'job_id': job.job_id(),
                        'status': str(job.status()),
                        'backend': self.backend.name,
                        'session_id': session.session_id,
                        'mode': 'session'
                    }

        except Exception as e:
            logger.error(f"IBM session execution failed: {e}")
            raise

    def run_batch(
        self,
        circuits: Union[List, 'QuantumCircuit'],
        observables: Optional[List] = None,
        shots: int = 1024,
        optimization_level: int = 1,
        resilience_level: int = 1
    ) -> Dict[str, Any]:
        """
        Run circuits in batch mode (for non-premium users).

        Args:
            circuits: Quantum circuit(s) to run
            observables: Observable(s) for Estimator (optional)
            shots: Number of measurement shots
            optimization_level: Transpilation optimization (0-3)
            resilience_level: Error mitigation level (0-2)

        Returns:
            Results dictionary
        """
        from qiskit_ibm_runtime import Batch, SamplerV2 as Sampler, EstimatorV2 as Estimator
        from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

        # Ensure circuits is a list
        if not isinstance(circuits, list):
            circuits = [circuits]

        logger.info(f"Running batch on {self.backend.name}")
        logger.info(f"  Circuits: {len(circuits)}")
        logger.info(f"  Shots: {shots}")

        # Transpile circuits for target hardware
        logger.info(f"  Transpiling circuits (optimization_level={optimization_level})...")
        pm = generate_preset_pass_manager(
            backend=self.backend,
            optimization_level=optimization_level
        )
        transpiled_circuits = pm.run(circuits)
        logger.info(f"  Transpilation complete")

        # Expand observables to match transpiled circuit size if needed
        if observables is not None:
            transpiled_observables = []
            for i, (circuit, observable) in enumerate(zip(transpiled_circuits, observables)):
                original_qubits = circuits[i].num_qubits
                transpiled_qubits = circuit.num_qubits

                if transpiled_qubits > original_qubits:
                    # Pad observable with identity operators
                    padding = "I" * (transpiled_qubits - original_qubits)
                    padded_terms = []
                    for pauli_str, coeff in observable.to_list():
                        padded_pauli = pauli_str + padding
                        padded_terms.append((padded_pauli, coeff))

                    from qiskit.quantum_info import SparsePauliOp
                    padded_observable = SparsePauliOp.from_list(padded_terms)
                    transpiled_observables.append(padded_observable)
                else:
                    transpiled_observables.append(observable)

            observables = transpiled_observables

        try:
            with Batch(backend=self.backend) as batch:
                if observables is not None:
                    # Use Estimator for energy expectation values
                    estimator = Estimator(mode=batch)

                    # Set options (V2 primitives)
                    estimator.options.default_shots = shots

                    # CRITICAL FIX (Issue #8): Use auto_configure() for error mitigation
                    error_mitigation = ErrorMitigationStrategy.auto_configure(self.backend.name)

                    # Apply auto-configured error mitigation strategy
                    estimator.options.resilience_level = error_mitigation.resilience_level
                    resilience_options = error_mitigation.get_resilience_options()

                    # Set resilience options if any
                    if resilience_options:
                        for key, value in resilience_options.items():
                            setattr(estimator.options.resilience, key, value)
                        logger.info(f"  ✅ Auto-configured error mitigation applied")
                        logger.info(f"     Resilience level: {error_mitigation.resilience_level}")
                        if error_mitigation.zne_extrapolation:
                            logger.info(f"     ZNE: {error_mitigation.zne_extrapolation}")
                        if error_mitigation.dynamical_decoupling:
                            logger.info(f"     DD: {error_mitigation.dynamical_decoupling}")

                    logger.info("Using Estimator primitive (batch mode)")

                    # Build pub (circuit, observable) tuples with transpiled circuits
                    pubs = [(transpiled_circuits[i], observables[i]) for i in range(len(transpiled_circuits))]

                    job = estimator.run(pubs)

                    # Return job immediately (non-blocking)
                    return {
                        'job_id': job.job_id(),
                        'status': str(job.status()),
                        'backend': self.backend.name,
                        'mode': 'batch'
                    }

                else:
                    # Use Sampler for measurement counts
                    sampler = Sampler(mode=batch)

                    # Set options
                    sampler.options.default_shots = shots

                    logger.info("Using Sampler primitive")

                    job = sampler.run(transpiled_circuits)

                    # Return job immediately (non-blocking)
                    return {
                        'job_id': job.job_id(),
                        'status': str(job.status()),
                        'backend': self.backend.name,
                        'mode': 'batch'
                    }

        except Exception as e:
            logger.error(f"IBM batch execution failed: {e}")
            raise

    def run_single(
        self,
        circuit: 'QuantumCircuit',
        observable: Optional[Any] = None,
        shots: int = 1024
    ) -> Dict[str, Any]:
        """Run a single circuit (convenience method)."""
        return self.run_batch([circuit], [observable] if observable else None, shots=shots)

    def get_backend_info(self) -> Dict[str, Any]:
        """Get information about current backend."""
        config = self.backend.configuration()

        info = {
            'name': self.backend.name,
            'num_qubits': self.backend.num_qubits,
            'is_simulator': self.backend.simulator,
            'is_operational': self.backend.status().operational,
            'pending_jobs': self.backend.status().pending_jobs,
            'basis_gates': config.basis_gates if hasattr(config, 'basis_gates') else None,
            'coupling_map': self.backend.coupling_map,
            'max_shots': config.max_shots if hasattr(config, 'max_shots') else 'unlimited'
        }

        return info

    def get_job_status(self, job_id: str) -> str:
        """Get status of a submitted job."""
        job = self.service.job(job_id)
        status = job.status()
        # Handle different Qiskit versions - status might be string or enum
        if isinstance(status, str):
            return status
        elif hasattr(status, 'name'):
            return status.name
        else:
            return str(status)

    def get_job_result(self, job_id: str) -> Any:
        """Retrieve results for a completed job."""
        job = self.service.job(job_id)
        return job.result()

    def list_backends(self, simulator: bool = False, operational: bool = True) -> List[str]:
        """List available backends."""
        backends = self.service.backends(simulator=simulator, operational=operational)
        return [b.name for b in backends]

    def cancel_job(self, job_id: str) -> Dict[str, Any]:
        """
        Cancel a running or queued job.

        Args:
            job_id: IBM job ID to cancel

        Returns:
            Dictionary with cancellation status
        """
        try:
            job = self.service.job(job_id)
            job.cancel()
            logger.info(f"IBM job {job_id} cancelled successfully")

            return {
                'status': 'cancelled',
                'job_id': job_id,
                'message': 'Job cancellation requested'
            }
        except Exception as e:
            logger.error(f"Failed to cancel IBM job {job_id}: {e}")
            raise RuntimeError(f"Failed to cancel job: {e}")

    def __repr__(self):
        return f"IBMBackend(backend='{self.backend.name}', qubits={self.backend.num_qubits})"
