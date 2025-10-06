"""
IBM Quantum job runner for Kanad solvers.
"""

import logging
from typing import Dict, Any, Optional
import numpy as np

logger = logging.getLogger(__name__)


class IBMRunner:
    """
    High-level runner for Kanad solvers on IBM Quantum.

    Automatically handles:
    - Job preparation
    - Circuit transpilation
    - Batch submission
    - Result processing
    """

    def __init__(self, backend: 'IBMBackend'):
        """
        Initialize runner with IBM backend.

        Args:
            backend: IBMBackend instance
        """
        self.backend = backend
        logger.info(f"IBMRunner initialized with {backend.backend.name}")

    def run_vqe(
        self,
        bond: 'BaseBond',
        ansatz_type: str = 'ucc',
        shots: int = 2048,
        optimization_level: int = 1,
        **kwargs
    ) -> Dict[str, Any]:
        """
        Run VQE on IBM Quantum.

        Args:
            bond: Molecular bond
            ansatz_type: Ansatz type
            shots: Number of measurement shots
            optimization_level: Transpilation optimization (0-3)
            **kwargs: Additional options

        Returns:
            VQE results
        """
        from kanad.backends.ibm import IBMPreparation

        logger.info("Running VQE on IBM Quantum")
        logger.info(f"  Backend: {self.backend.backend.name}")
        logger.info(f"  Ansatz: {ansatz_type}")

        # Prepare job
        prep = IBMPreparation(bond, ansatz_type=ansatz_type)

        # Get circuits and observables
        circuits, observables = prep.prepare_vqe_circuits()

        # Transpile for hardware
        circuits = prep.transpile_for_hardware(
            circuits,
            self.backend,
            optimization_level=optimization_level
        )

        # Run on IBM
        results = self.backend.run_batch(
            circuits,
            observables,
            shots=shots,
            optimization_level=optimization_level
        )

        # Process results
        energy = results['values'][0] if results['values'] else None

        output = {
            'energy': float(energy) if energy is not None else None,
            'hf_energy': float(prep.hf_energy),
            'job_id': results['job_id'],
            'metadata': results['metadata'],
            'backend': self.backend.backend.name,
            'shots': shots,
            'preparation_summary': prep.get_preparation_summary()
        }

        if energy is not None:
            logger.info(f"VQE complete: E = {energy:.8f} Ha")

        return output

    def run_sqd(
        self,
        bond: 'BaseBond',
        subspace_dim: int = 10,
        shots: int = 4096,
        **kwargs
    ) -> Dict[str, Any]:
        """
        Run SQD on IBM Quantum.

        Args:
            bond: Molecular bond
            subspace_dim: Subspace dimension
            shots: Number of shots per circuit
            **kwargs: Additional options

        Returns:
            SQD results
        """
        from kanad.backends.ibm import IBMPreparation

        logger.info("Running SQD on IBM Quantum")
        logger.info(f"  Backend: {self.backend.backend.name}")
        logger.info(f"  Subspace dim: {subspace_dim}")

        # Prepare job
        prep = IBMPreparation(bond, ansatz_type='ucc')

        # Get circuits
        circuits = prep.prepare_sqd_circuits(subspace_dim=subspace_dim)

        # Transpile
        circuits = prep.transpile_for_hardware(circuits, self.backend)

        # Run batch
        results = self.backend.run_batch(circuits, shots=shots)

        # Process counts to build subspace
        counts_list = results['counts']

        # Classical post-processing would happen here
        # For now, return raw results

        output = {
            'counts': counts_list,
            'job_id': results['job_id'],
            'metadata': results['metadata'],
            'backend': self.backend.backend.name,
            'subspace_dim': subspace_dim,
            'preparation_summary': prep.get_preparation_summary()
        }

        logger.info("SQD complete")

        return output

    def __repr__(self):
        return f"IBMRunner(backend='{self.backend.backend.name}')"
