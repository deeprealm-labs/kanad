"""
Qiskit backend integration for hardware execution.

Supports Qiskit 2.x primitives (Estimator, Sampler) for both
local simulation and IBM Quantum hardware.
"""

from typing import Optional, Union, Dict, Any
import warnings
import logging

logger = logging.getLogger(__name__)


class QiskitBackend:
    """
    Manages Qiskit backends for quantum circuit execution.

    Supports:
    - Local Aer simulators (statevector, qasm)
    - IBM Quantum hardware
    - Third-party providers

    Uses Qiskit 2.x Estimator primitive for VQE energy estimation.
    """

    def __init__(
        self,
        backend_name: str = 'aer_simulator',
        shots: int = 1024,
        optimization_level: int = 1,
        **backend_options
    ):
        """
        Initialize Qiskit backend.

        Args:
            backend_name: Backend identifier
                - 'aer_simulator': Local Aer simulator (default)
                - 'aer_simulator_statevector': Exact statevector simulation
                - 'ibm_*': IBM Quantum hardware (requires credentials)
            shots: Number of measurement shots (for hardware/qasm simulation)
            optimization_level: Transpiler optimization level (0-3)
            **backend_options: Additional backend-specific options
        """
        self.backend_name = backend_name
        self.shots = shots
        self.optimization_level = optimization_level
        self.backend_options = backend_options

        # Initialize backend
        self.backend = None
        self._service = None
        self._estimator = None

        self._initialize_backend()

    def _initialize_backend(self):
        """Initialize the appropriate Qiskit backend."""
        try:
            if self.backend_name == 'aer_simulator':
                self._initialize_aer_simulator()
            elif self.backend_name == 'aer_simulator_statevector':
                self._initialize_aer_statevector()
            elif self.backend_name.startswith('ibm_'):
                self._initialize_ibm_backend()
            else:
                raise ValueError(
                    f"Unknown backend: {self.backend_name}. "
                    f"Supported: 'aer_simulator', 'aer_simulator_statevector', 'ibm_*'"
                )
        except ImportError as e:
            raise ImportError(
                f"Failed to initialize backend '{self.backend_name}': {e}\n"
                f"Install required packages: pip install qiskit qiskit-aer"
            )

    def _initialize_aer_simulator(self):
        """Initialize Aer QASM simulator."""
        try:
            from qiskit_aer import AerSimulator
        except ImportError:
            raise ImportError(
                "qiskit-aer not installed. Install with: pip install qiskit-aer"
            )

        self.backend = AerSimulator(**self.backend_options)
        logger.info(f"Initialized Aer simulator: {self.backend}")

    def _initialize_aer_statevector(self):
        """Initialize Aer statevector simulator (exact, no shots)."""
        try:
            from qiskit_aer import AerSimulator
        except ImportError:
            raise ImportError(
                "qiskit-aer not installed. Install with: pip install qiskit-aer"
            )

        self.backend = AerSimulator(method='statevector', **self.backend_options)
        self.shots = None  # Statevector is exact
        logger.info("Initialized Aer statevector simulator (exact computation)")

    def _initialize_ibm_backend(self):
        """Initialize IBM Quantum hardware backend."""
        try:
            from qiskit_ibm_runtime import QiskitRuntimeService
        except ImportError:
            raise ImportError(
                "qiskit-ibm-runtime not installed. "
                "Install with: pip install qiskit-ibm-runtime"
            )

        # Load IBM Quantum credentials
        try:
            self._service = QiskitRuntimeService()
            self.backend = self._service.backend(self.backend_name)
            logger.info(f"Connected to IBM backend: {self.backend_name}")
            logger.info(f"Status: {self.backend.status().status_msg}")
            logger.info(f"Qubits: {self.backend.num_qubits}")
        except Exception as e:
            raise RuntimeError(
                f"Failed to connect to IBM backend '{self.backend_name}': {e}\n"
                f"Configure credentials: QiskitRuntimeService.save_account(token='YOUR_TOKEN')"
            )

    def get_estimator(self, **options):
        """
        Get Estimator primitive for energy expectation value computation.

        Args:
            **options: Estimator-specific options (resilience_level, etc.)

        Returns:
            Estimator primitive (Qiskit 2.x)
        """
        if self._estimator is not None:
            return self._estimator

        # Choose appropriate Estimator based on backend
        if self.backend_name.startswith('aer_'):
            self._estimator = self._get_aer_estimator(**options)
        elif self.backend_name.startswith('ibm_'):
            self._estimator = self._get_runtime_estimator(**options)
        else:
            self._estimator = self._get_generic_estimator(**options)

        return self._estimator

    def _get_aer_estimator(self, **options):
        """Get Aer-based Estimator."""
        from qiskit_aer.primitives import Estimator

        # Use Aer's native Estimator primitive
        estimator = Estimator()

        return estimator

    def _get_runtime_estimator(self, **options):
        """Get IBM Runtime Estimator (for real hardware)."""
        try:
            from qiskit_ibm_runtime import EstimatorV2 as Estimator
            from qiskit_ibm_runtime import Session
        except ImportError:
            raise ImportError(
                "qiskit-ibm-runtime not installed. "
                "Install with: pip install qiskit-ibm-runtime"
            )

        # Create session for batching jobs
        session = Session(backend=self.backend)

        # Default runtime options
        estimator_options = {
            'resilience_level': options.get('resilience_level', 1),
            'optimization_level': self.optimization_level,
        }
        estimator_options.update(options)

        estimator = Estimator(session=session, options=estimator_options)

        logger.info(f"Created Runtime Estimator with resilience_level={estimator_options['resilience_level']}")

        return estimator

    def _get_generic_estimator(self, **options):
        """Get generic Estimator for custom backends."""
        from qiskit_aer.primitives import Estimator

        return Estimator()

    def get_sampler(self, **options):
        """
        Get Sampler primitive for circuit sampling.

        Args:
            **options: Sampler-specific options

        Returns:
            Sampler primitive (Qiskit 2.x)
        """
        if self.backend_name.startswith('aer_'):
            from qiskit_aer.primitives import Sampler
            return Sampler()
        elif self.backend_name.startswith('ibm_'):
            from qiskit_ibm_runtime import SamplerV2 as Sampler
            from qiskit_ibm_runtime import Session
            session = Session(backend=self.backend)
            return Sampler(session=session, options=options)
        else:
            from qiskit_aer.primitives import Sampler
            return Sampler()

    def transpile_circuit(self, circuit, **transpile_options):
        """
        Transpile circuit for target backend.

        Args:
            circuit: Qiskit QuantumCircuit
            **transpile_options: Options for qiskit.transpile()

        Returns:
            Transpiled circuit
        """
        from qiskit import transpile

        options = {
            'backend': self.backend,
            'optimization_level': self.optimization_level,
        }
        options.update(transpile_options)

        return transpile(circuit, **options)

    def get_backend_info(self) -> Dict[str, Any]:
        """
        Get backend information and status.

        Returns:
            Dictionary with backend properties
        """
        info = {
            'name': self.backend_name,
            'shots': self.shots,
            'optimization_level': self.optimization_level,
        }

        if hasattr(self.backend, 'num_qubits'):
            info['num_qubits'] = self.backend.num_qubits

        if hasattr(self.backend, 'status'):
            status = self.backend.status()
            info['status'] = status.status_msg
            info['pending_jobs'] = status.pending_jobs

        return info

    def __repr__(self) -> str:
        """String representation."""
        return f"QiskitBackend(name='{self.backend_name}', shots={self.shots})"
