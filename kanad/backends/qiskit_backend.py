"""
Qiskit Backend wrapper for Kanad Framework.

Provides a simple interface to Qiskit Aer simulators.
"""

import logging
from typing import Optional, Dict, Any

logger = logging.getLogger(__name__)


class QiskitBackend:
    """
    Wrapper for Qiskit Aer backends.

    Supports:
    - aer_simulator: General purpose simulator
    - aer_simulator_statevector: Exact statevector simulator
    """

    def __init__(
        self,
        backend_name: str = 'aer_simulator',
        shots: Optional[int] = 1024,
        **kwargs
    ):
        """
        Initialize Qiskit backend.

        Args:
            backend_name: Name of Qiskit Aer backend
            shots: Number of shots (None for statevector)
            **kwargs: Additional backend options
        """
        try:
            from qiskit_aer import Aer
        except ImportError:
            raise ImportError(
                "Qiskit Aer is required for QiskitBackend. "
                "Install with: pip install qiskit-aer"
            )

        self.backend_name = backend_name
        self.shots = None if 'statevector' in backend_name else shots
        self.backend_options = kwargs

        # Get backend from Aer
        try:
            self.backend = Aer.get_backend(backend_name)
        except Exception as e:
            logger.error(f"Failed to get Aer backend '{backend_name}': {e}")
            raise

        logger.info(f"Initialized Qiskit backend: {backend_name}")

    def get_estimator(self):
        """
        Get Qiskit Estimator primitive.

        Returns:
            Estimator instance
        """
        try:
            # Try new Qiskit primitives API
            from qiskit.primitives import BackendEstimator
            estimator = BackendEstimator(backend=self.backend)
        except (ImportError, TypeError):
            try:
                # Try Aer primitives
                from qiskit_aer.primitives import Estimator
                estimator = Estimator()
            except ImportError:
                # Fallback to older API
                from qiskit.primitives import Estimator
                estimator = Estimator()

        return estimator

    def get_sampler(self):
        """
        Get Qiskit Sampler primitive.

        Returns:
            Sampler instance
        """
        try:
            from qiskit.primitives import BackendSampler
        except ImportError:
            from qiskit_aer.primitives import Sampler as BackendSampler

        sampler = BackendSampler(backend=self.backend)
        return sampler

    def get_backend_info(self) -> Dict[str, Any]:
        """
        Get backend information.

        Returns:
            Dictionary with backend info
        """
        return {
            'name': self.backend_name,
            'shots': self.shots,
            'backend': str(self.backend),
        }

    def run(self, circuit, **kwargs):
        """
        Run circuit on backend.

        Args:
            circuit: Qiskit QuantumCircuit
            **kwargs: Additional run options

        Returns:
            Job result
        """
        run_kwargs = {'shots': self.shots} if self.shots else {}
        run_kwargs.update(kwargs)

        job = self.backend.run(circuit, **run_kwargs)
        return job.result()
