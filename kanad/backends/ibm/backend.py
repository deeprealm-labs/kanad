"""
IBM Quantum Runtime backend for cloud quantum computing.

Provides access to IBM Quantum hardware and cloud simulators via
IBM Quantum Runtime service with proper authentication using
API keys and Cloud Resource Names (CRN).
"""

from typing import Optional, Dict, Any, List
import os
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


class IBMRuntimeBackend:
    """
    IBM Quantum Runtime backend for executing circuits on IBM Quantum hardware.

    Supports:
    - IBM Quantum Platform (cloud-based)
    - IBM Cloud instances (with CRN)
    - Runtime primitives (EstimatorV2, SamplerV2)
    - Session management for job batching
    - Error mitigation and resilience options

    Authentication methods:
    1. Environment variables (recommended for production):
       - IBM_QUANTUM_TOKEN or QISKIT_IBM_TOKEN: API key
       - IBM_QUANTUM_CRN or QISKIT_IBM_INSTANCE: Cloud Resource Name (for IBM Cloud)
       - IBM_QUANTUM_CHANNEL or QISKIT_IBM_CHANNEL: 'ibm_cloud' or 'ibm_quantum_platform'

    2. Direct parameters (for testing/development):
       - Pass token and instance to __init__()

    3. Saved credentials:
       - Previously saved via QiskitRuntimeService.save_account()
    """

    def __init__(
        self,
        backend_name: Optional[str] = None,
        token: Optional[str] = None,
        instance: Optional[str] = None,
        channel: Optional[str] = None,  # Will auto-detect based on CRN presence
        shots: int = 4096,
        optimization_level: int = 3,
        resilience_level: int = 1,
        **backend_options
    ):
        """
        Initialize IBM Runtime backend.

        Args:
            backend_name: Name of IBM backend (e.g., 'ibm_kyoto', 'ibm_osaka')
                         If None, uses least busy backend
            token: IBM Quantum API key (or set IBM_QUANTUM_TOKEN env var)
            instance: IBM Cloud CRN (or set IBM_QUANTUM_CRN env var)
                     Required for IBM Cloud channel
            channel: 'ibm_quantum_platform' (IBM Quantum Platform) or 'ibm_cloud' (IBM Cloud)
            shots: Number of measurement shots
            optimization_level: Transpiler optimization (0-3, default 3 for hardware)
            resilience_level: Error mitigation level (0-2, default 1)
            **backend_options: Additional backend options
        """
        self.backend_name = backend_name
        self.shots = shots
        self.optimization_level = optimization_level
        self.resilience_level = resilience_level
        self.backend_options = backend_options

        # Get credentials from environment or parameters
        # Support multiple environment variable naming conventions
        self.token = (
            token or
            os.getenv('IBM_QUANTUM_TOKEN') or
            os.getenv('QISKIT_IBM_TOKEN') or
            os.getenv('API')  # Legacy/custom naming
        )
        self.instance = (
            instance or
            os.getenv('IBM_QUANTUM_CRN') or
            os.getenv('QISKIT_IBM_INSTANCE') or
            os.getenv('CRN')  # Legacy/custom naming
        )
        # Determine channel (ibm_cloud and ibm_quantum_platform are interchangeable)
        # Default to ibm_cloud if CRN is provided, otherwise ibm_quantum_platform
        default_channel = 'ibm_cloud' if self.instance else 'ibm_quantum_platform'
        channel_env = (
            os.getenv('IBM_QUANTUM_CHANNEL') or
            os.getenv('QISKIT_IBM_CHANNEL')
        )

        # Map legacy channel names to correct ones for both env and parameter
        if channel == 'ibm_quantum':
            channel = 'ibm_quantum_platform'
        if channel_env == 'ibm_quantum':
            channel_env = 'ibm_quantum_platform'

        self.channel = channel or channel_env or default_channel

        # Initialize service and backend
        self._service = None
        self.backend = None
        self._session = None
        self._estimator = None
        self._sampler = None

        self._initialize_service()
        self._initialize_backend()

    def _initialize_service(self):
        """Initialize QiskitRuntimeService with credentials."""
        try:
            from qiskit_ibm_runtime import QiskitRuntimeService
        except ImportError:
            raise ImportError(
                "qiskit-ibm-runtime not installed. "
                "Install with: pip install qiskit-ibm-runtime"
            )

        try:
            # Try to initialize with explicit credentials first
            if self.token:
                kwargs = {'channel': self.channel, 'token': self.token}
                if self.instance:
                    kwargs['instance'] = self.instance

                self._service = QiskitRuntimeService(**kwargs)
                logger.info(f"Initialized IBM Runtime service with {self.channel} channel")
            else:
                # Fall back to saved credentials
                self._service = QiskitRuntimeService()
                logger.info("Initialized IBM Runtime service with saved credentials")

        except Exception as e:
            error_msg = (
                f"Failed to initialize IBM Runtime service: {e}\n\n"
                f"Please provide credentials via:\n"
                f"1. Environment variables:\n"
                f"   - IBM_QUANTUM_TOKEN or QISKIT_IBM_TOKEN (API key)\n"
                f"   - IBM_QUANTUM_CRN or QISKIT_IBM_INSTANCE (CRN, for IBM Cloud)\n"
                f"   - IBM_QUANTUM_CHANNEL or QISKIT_IBM_CHANNEL (ibm_quantum_platform or ibm_cloud)\n"
                f"2. Constructor parameters: token, instance, channel\n"
                f"3. Save account: QiskitRuntimeService.save_account(channel='{self.channel}', "
                f"token='YOUR_TOKEN', instance='YOUR_CRN')\n\n"
                f"Get credentials from: https://quantum.ibm.com/account"
            )
            raise RuntimeError(error_msg)

    def _initialize_backend(self):
        """Initialize the IBM Quantum backend."""
        try:
            if self.backend_name:
                # Use specified backend
                self.backend = self._service.backend(self.backend_name)
                logger.info(f"Connected to backend: {self.backend_name}")
            else:
                # Get least busy backend
                self.backend = self._service.least_busy(
                    operational=True,
                    simulator=False,
                    min_num_qubits=5
                )
                self.backend_name = self.backend.name
                logger.info(f"Selected least busy backend: {self.backend_name}")

            # Log backend info
            if hasattr(self.backend, 'status'):
                status = self.backend.status()
                logger.info(f"Backend status: {status.status_msg}")
                logger.info(f"Pending jobs: {status.pending_jobs}")

            if hasattr(self.backend, 'num_qubits'):
                logger.info(f"Available qubits: {self.backend.num_qubits}")

        except Exception as e:
            raise RuntimeError(
                f"Failed to connect to backend '{self.backend_name}': {e}\n"
                f"Available backends: {self.list_backends()}"
            )

    def list_backends(self, filters: Optional[Dict[str, Any]] = None) -> List[str]:
        """
        List available IBM Quantum backends.

        Args:
            filters: Optional filters (e.g., {'simulator': False, 'operational': True})

        Returns:
            List of backend names
        """
        if self._service is None:
            return []

        try:
            backends = self._service.backends(**(filters or {}))
            return [b.name for b in backends]
        except Exception as e:
            logger.warning(f"Failed to list backends: {e}")
            return []

    def get_estimator(self, **options):
        """
        Get IBM Runtime EstimatorV2 primitive for energy estimation.

        Args:
            **options: Estimator options (override defaults)
                - resilience_level: Error mitigation level (0-2)
                - optimization_level: Transpiler optimization (0-3)
                - default_shots: Number of shots per circuit

        Returns:
            EstimatorV2 primitive configured for this backend
        """
        if self._estimator is not None:
            return self._estimator

        try:
            from qiskit_ibm_runtime import EstimatorV2 as Estimator
            from qiskit_ibm_runtime import Session
        except ImportError:
            raise ImportError(
                "qiskit-ibm-runtime not installed. "
                "Install with: pip install qiskit-ibm-runtime"
            )

        # Try to create session (not available on open/free plans)
        session_mode = options.pop('session', True)  # Default to trying session
        if session_mode and self._session is None:
            try:
                self._session = Session(backend=self.backend)
                logger.info(f"Created runtime session: {self._session.session_id}")
            except Exception as e:
                logger.warning(f"Could not create session (likely using open plan): {e}")
                logger.info("Running without session mode")
                self._session = None

        # Configure estimator options (EstimatorV2 options)
        # Note: optimization_level is not an Estimator option - transpilation happens separately
        estimator_options = {
            'resilience_level': self.resilience_level,
            'default_shots': self.shots,
        }
        estimator_options.update(options)

        # Create estimator (with or without session)
        if self._session:
            self._estimator = Estimator(mode=self._session, options=estimator_options)
        else:
            # Use backend directly without session (for open plan)
            self._estimator = Estimator(mode=self.backend, options=estimator_options)

        logger.info(
            f"Created EstimatorV2 with resilience_level={estimator_options['resilience_level']}, "
            f"shots={estimator_options['default_shots']}"
        )

        return self._estimator

    def get_sampler(self, **options):
        """
        Get IBM Runtime SamplerV2 primitive for circuit sampling.

        Args:
            **options: Sampler options
                - default_shots: Number of shots per circuit

        Returns:
            SamplerV2 primitive configured for this backend
        """
        if self._sampler is not None:
            return self._sampler

        try:
            from qiskit_ibm_runtime import SamplerV2 as Sampler
            from qiskit_ibm_runtime import Session
        except ImportError:
            raise ImportError(
                "qiskit-ibm-runtime not installed. "
                "Install with: pip install qiskit-ibm-runtime"
            )

        # Try to create session (not available on open/free plans)
        session_mode = options.pop('session', True)
        if session_mode and self._session is None:
            try:
                self._session = Session(backend=self.backend)
                logger.info(f"Created runtime session: {self._session.session_id}")
            except Exception as e:
                logger.warning(f"Could not create session (likely using open plan): {e}")
                logger.info("Running without session mode")
                self._session = None

        # Configure sampler options
        sampler_options = {
            'default_shots': self.shots,
        }
        sampler_options.update(options)

        # Create sampler (with or without session)
        if self._session:
            self._sampler = Sampler(mode=self._session, options=sampler_options)
        else:
            self._sampler = Sampler(mode=self.backend, options=sampler_options)

        logger.info(f"Created SamplerV2 with shots={sampler_options['default_shots']}")

        return self._sampler

    def transpile_circuit(self, circuit, **transpile_options):
        """
        Transpile circuit for IBM backend.

        Args:
            circuit: Qiskit QuantumCircuit
            **transpile_options: Options for qiskit.transpile()

        Returns:
            Transpiled circuit optimized for target backend
        """
        from qiskit import transpile

        options = {
            'backend': self.backend,
            'optimization_level': self.optimization_level,
        }
        options.update(transpile_options)

        transpiled = transpile(circuit, **options)

        logger.info(
            f"Transpiled circuit: {circuit.num_qubits} qubits, "
            f"{circuit.depth()} depth → {transpiled.depth()} depth"
        )

        return transpiled

    def get_backend_info(self) -> Dict[str, Any]:
        """
        Get comprehensive backend information.

        Returns:
            Dictionary with backend properties and status
        """
        info = {
            'name': self.backend_name,
            'channel': self.channel,
            'shots': self.shots,
            'optimization_level': self.optimization_level,
            'resilience_level': self.resilience_level,
        }

        if self.backend:
            if hasattr(self.backend, 'num_qubits'):
                info['num_qubits'] = self.backend.num_qubits

            if hasattr(self.backend, 'version'):
                info['version'] = self.backend.version

            if hasattr(self.backend, 'status'):
                status = self.backend.status()
                info['status'] = status.status_msg
                info['pending_jobs'] = status.pending_jobs
                info['operational'] = status.operational

        if self._session:
            info['session_id'] = self._session.session_id

        return info

    def close(self):
        """Close the runtime session and release resources."""
        if self._session:
            try:
                self._session.close()
                logger.info(f"Closed runtime session: {self._session.session_id}")
            except Exception as e:
                logger.warning(f"Error closing session: {e}")
            finally:
                self._session = None
                self._estimator = None
                self._sampler = None

    def __enter__(self):
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit - close session."""
        self.close()

    def __repr__(self) -> str:
        """String representation."""
        return (
            f"IBMRuntimeBackend(name='{self.backend_name}', "
            f"channel='{self.channel}', shots={self.shots})"
        )

    def __del__(self):
        """Destructor - ensure session is closed."""
        self.close()
