"""
BlueQubit Backend Implementation

Provides interface to BlueQubit cloud quantum simulators.
"""

import os
import logging
from typing import Dict, Any, Optional, List
import numpy as np

logger = logging.getLogger(__name__)


class BlueQubitBackend:
    """
    BlueQubit cloud backend for quantum simulations.

    Supports:
    - GPU simulators (free, 36 qubits, fast)
    - CPU simulators (34 qubits)
    - MPS tensor network (40+ qubits, requires balance)

    Usage:
        backend = BlueQubitBackend(device='gpu', api_token='your_token')
        result = backend.run_circuit(circuit, shots=1024)
    """

    SUPPORTED_DEVICES = ['cpu', 'gpu', 'mps.cpu', 'mps.gpu', 'pauli-path']

    def __init__(
        self,
        device: str = 'gpu',
        api_token: Optional[str] = None,
        method: Optional[str] = None,  # Legacy parameter, maps to device
        **options
    ):
        """
        Initialize BlueQubit backend.

        Args:
            device: Device type ('cpu', 'gpu', 'mps.cpu', 'mps.gpu', 'pauli-path')
            api_token: BlueQubit API token (or set BLUEQUBIT_API_TOKEN env var)
            method: Legacy alias for device parameter
            **options: Additional device options:
                - mps_bond_dimension: For MPS devices (default: 100)
                - pauli_path_truncation_threshold: For pauli-path (min: 1e-5)
                - pauli_path_circuit_transpilation_level: Transpilation level
        """
        # Support legacy 'method' parameter
        if method is not None:
            device = method

        self.device = device
        self.options = options

        # Get API token
        self.api_token = api_token or os.getenv('BLUE_TOKEN')

        if not self.api_token:
            raise ValueError(
                "BlueQubit API token required. Set BLUE_TOKEN environment variable "
                "or pass api_token parameter. Get token from https://app.bluequbit.io"
            )

        if device not in self.SUPPORTED_DEVICES:
            raise ValueError(f"Device must be one of {self.SUPPORTED_DEVICES}")

        # Initialize client
        self._init_client()

        logger.info(f"BlueQubit backend initialized: device={device}")

    def _init_client(self):
        """Initialize BlueQubit client."""
        try:
            import bluequbit
            self.bq = bluequbit.init(self.api_token)
            logger.info("BlueQubit client initialized successfully")
        except ImportError:
            raise ImportError(
                "bluequbit package required. Install with: pip install bluequbit"
            )
        except Exception as e:
            raise RuntimeError(f"Failed to initialize BlueQubit client: {e}")

    def run_circuit(
        self,
        circuit,
        shots: Optional[int] = None,
        job_name: Optional[str] = None,
        asynchronous: bool = False
    ) -> Dict[str, Any]:
        """
        Run a quantum circuit on BlueQubit.

        Args:
            circuit: Qiskit QuantumCircuit
            shots: Number of measurement shots (None for statevector)
            job_name: Optional job name for tracking
            asynchronous: If True, return job handle without waiting

        Returns:
            Result dictionary with counts/statevector
        """
        job_name = job_name or f"kanad_{self.device}"

        logger.info(f"Submitting circuit to BlueQubit ({self.device})")
        print(f"🔧 BlueQubit device selected: {self.device}")
        print(f"   Circuit: {circuit.num_qubits} qubits, depth {circuit.depth()}")
        logger.info(f"  Qubits: {circuit.num_qubits}, Depth: {circuit.depth()}")

        try:
            # Run circuit
            if asynchronous:
                job = self.bq.run(
                    circuit,
                    device=self.device,
                    job_name=job_name,
                    options=self.options,
                    asynchronous=True
                )
                return {'job': job, 'job_id': job.job_id}
            else:
                result = self.bq.run(
                    circuit,
                    device=self.device,
                    job_name=job_name,
                    options=self.options
                )

                # Extract results
                output = {}

                if shots is None:
                    # Get statevector
                    output['statevector'] = result.get_statevector()
                    logger.info("Retrieved statevector")
                else:
                    # Get counts
                    output['counts'] = result.get_counts()
                    logger.info(f"Retrieved counts ({shots} shots)")

                return output

        except Exception as e:
            logger.error(f"BlueQubit execution failed: {e}")
            raise

    def wait_for_job(self, job_id: str) -> Dict[str, Any]:
        """Wait for asynchronous job to complete."""
        logger.info(f"Waiting for job {job_id}")
        result = self.bq.wait(job_id)

        return {
            'statevector': result.get_statevector(),
            'counts': result.get_counts() if hasattr(result, 'get_counts') else None
        }

    def cancel_job(self, job_id: str):
        """Cancel a pending job."""
        logger.info(f"Cancelling job {job_id}")
        self.bq.cancel(job_id)

    def get_device_info(self) -> Dict[str, Any]:
        """Get information about current device."""
        info = {
            'device': self.device,
            'max_qubits': self._get_max_qubits(),
            'supports_statevector': True,
            'supports_shots': self.device not in ['mps.cpu', 'mps.gpu'],
            'is_gpu_accelerated': 'gpu' in self.device
        }

        return info

    def _get_max_qubits(self) -> int:
        """Get maximum qubits for device."""
        limits = {
            'cpu': 34,
            'gpu': 36,
            'mps.cpu': 40,
            'mps.gpu': 40,
            'pauli-path': 50  # Depends on truncation threshold
        }
        return limits.get(self.device, 30)

    def __repr__(self):
        return f"BlueQubitBackend(device='{self.device}')"
