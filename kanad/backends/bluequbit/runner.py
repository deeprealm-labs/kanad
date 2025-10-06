"""
BlueQubit job runner for Kanad solvers.

Integrates BlueQubit backend with Kanad VQE, SQD, and other solvers.
"""

import logging
from typing import Dict, Any, Optional
import numpy as np

logger = logging.getLogger(__name__)


class BlueQubitRunner:
    """
    High-level runner for Kanad solvers on BlueQubit.

    Automatically handles:
    - Circuit preparation
    - Hamiltonian mapping
    - Job submission
    - Result processing
    """

    def __init__(self, backend: 'BlueQubitBackend'):
        """
        Initialize runner with BlueQubit backend.

        Args:
            backend: BlueQubitBackend instance
        """
        self.backend = backend
        logger.info(f"BlueQubitRunner initialized with {backend.device}")

    def run_vqe(
        self,
        bond: 'BaseBond',
        ansatz_type: str = 'ucc',
        optimizer: str = 'COBYLA',
        max_iterations: int = 100,
        **kwargs
    ) -> Dict[str, Any]:
        """
        Run VQE on BlueQubit.

        Args:
            bond: Molecular bond
            ansatz_type: Ansatz type ('ucc', 'hardware_efficient')
            optimizer: Classical optimizer
            max_iterations: Max optimization iterations
            **kwargs: Additional VQE options

        Returns:
            VQE results dictionary
        """
        from kanad.solvers import VQESolver

        logger.info("Running VQE on BlueQubit")
        logger.info(f"  Device: {self.backend.device}")
        logger.info(f"  Ansatz: {ansatz_type}")
        logger.info(f"  Optimizer: {optimizer}")

        # Create VQE solver with BlueQubit backend
        solver = VQESolver(
            bond=bond,
            ansatz_type=ansatz_type,
            optimizer_method=optimizer,
            max_iterations=max_iterations,
            backend='bluequbit',
            backend_options={'bluequbit_backend': self.backend},
            **kwargs
        )

        # Run solver
        results = solver.solve()

        logger.info(f"VQE complete: E = {results['energy']:.8f} Ha")

        return results

    def run_sqd(
        self,
        bond: 'BaseBond',
        subspace_dim: int = 10,
        n_states: int = 3,
        **kwargs
    ) -> Dict[str, Any]:
        """
        Run SQD on BlueQubit.

        Args:
            bond: Molecular bond
            subspace_dim: Subspace dimension
            n_states: Number of states to compute
            **kwargs: Additional SQD options

        Returns:
            SQD results dictionary
        """
        from kanad.solvers import SQDSolver

        logger.info("Running SQD on BlueQubit")
        logger.info(f"  Device: {self.backend.device}")
        logger.info(f"  Subspace dim: {subspace_dim}")

        # Create SQD solver with BlueQubit backend
        solver = SQDSolver(
            bond=bond,
            subspace_dim=subspace_dim,
            backend='bluequbit',
            backend_options={'bluequbit_backend': self.backend},
            **kwargs
        )

        # Run solver
        results = solver.solve(n_states=n_states)

        logger.info(f"SQD complete: E = {results['energy']:.8f} Ha")

        return results

    def estimate_cost(
        self,
        n_qubits: int,
        circuit_depth: int,
        n_iterations: int = 100
    ) -> Dict[str, Any]:
        """
        Estimate computational cost for BlueQubit job.

        Args:
            n_qubits: Number of qubits
            circuit_depth: Circuit depth
            n_iterations: Number of VQE iterations

        Returns:
            Cost estimate dictionary
        """
        device_info = self.backend.get_device_info()

        # Check limits
        if n_qubits > device_info['max_qubits']:
            logger.warning(
                f"Circuit requires {n_qubits} qubits but device max is {device_info['max_qubits']}"
            )

        # Estimate runtime
        # GPU is ~10x faster than CPU for large circuits
        base_time_per_iter = 0.1  # seconds for small circuits
        scaling = 2 ** (n_qubits / 10)  # Exponential scaling
        gpu_speedup = 10 if device_info['is_gpu_accelerated'] else 1

        estimated_time = (base_time_per_iter * scaling * n_iterations) / gpu_speedup

        cost_estimate = {
            'device': self.backend.device,
            'n_qubits': n_qubits,
            'circuit_depth': circuit_depth,
            'n_iterations': n_iterations,
            'estimated_time_seconds': estimated_time,
            'estimated_time_minutes': estimated_time / 60,
            'within_limits': n_qubits <= device_info['max_qubits'],
            'gpu_accelerated': device_info['is_gpu_accelerated'],
            'recommendation': self._get_device_recommendation(n_qubits)
        }

        return cost_estimate

    def _get_device_recommendation(self, n_qubits: int) -> str:
        """Get device recommendation based on problem size."""
        if n_qubits <= 20:
            return "cpu (sufficient, free)"
        elif n_qubits <= 36:
            return "gpu (recommended, free, faster)"
        else:
            return "mps.gpu (required for >36 qubits, requires balance)"

    def __repr__(self):
        return f"BlueQubitRunner(device='{self.backend.device}')"
