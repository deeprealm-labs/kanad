"""
Cloud Service - IBM Quantum and BlueQubit integration.

Handles submission of jobs to cloud quantum backends.
"""

from typing import Dict, Any, Optional
import logging
import asyncio

logger = logging.getLogger(__name__)


class CloudService:
    """Service for cloud quantum backend integration."""

    def __init__(self):
        """Initialize cloud service."""
        self.ibm_runners = {}  # Cache IBM runners per user
        self.bluequbit_runners = {}
        logger.info("CloudService initialized")

    async def submit_to_ibm(
        self,
        molecule: Any,
        config: Dict[str, Any],
        credentials: Dict[str, str]
    ) -> Dict[str, Any]:
        """
        Submit job to IBM Quantum.

        Args:
            molecule: Kanad Molecule object
            config: Simulation configuration
            credentials: IBM credentials (api_token, crn)

        Returns:
            Job submission result with cloud_job_id
        """
        logger.info("Submitting job to IBM Quantum...")

        try:
            from kanad.backends.ibm import IBMBackend, IBMRunner

            # Create IBM backend
            backend = IBMBackend(
                backend_name=config['backend'].get('backend_name', 'ibm_torino'),
                api_token=credentials['api_token'],
                crn=credentials.get('crn'),
                channel=credentials.get('channel', 'ibm_quantum')
            )

            # Create runner
            runner = IBMRunner(backend)

            # Submit VQE job
            result = runner.run_vqe(
                bond=molecule,  # Can accept Molecule or Bond
                ansatz_type=config.get('ansatz', 'ucc'),
                shots=config.get('shots', 2048)
            )

            logger.info(f"Job submitted to IBM: {result['job_id']}")

            return {
                'cloud_job_id': result['job_id'],
                'backend': result['backend'],
                'status': result['status']
            }

        except Exception as e:
            logger.error(f"IBM submission failed: {e}")
            raise

    async def submit_to_bluequbit(
        self,
        molecule: Any,
        config: Dict[str, Any],
        credentials: Dict[str, str]
    ) -> Dict[str, Any]:
        """
        Submit job to BlueQubit.

        Args:
            molecule: Kanad Molecule object
            config: Simulation configuration
            credentials: BlueQubit credentials

        Returns:
            Job submission result
        """
        logger.info("Submitting job to BlueQubit...")

        try:
            from kanad.backends.bluequbit import BlueQubitBackend, BlueQubitRunner

            # Create BlueQubit backend
            backend = BlueQubitBackend(
                api_token=credentials['api_token']
            )

            # Create runner
            runner = BlueQubitRunner(backend)

            # Submit job
            result = runner.run_vqe(
                bond=molecule,
                ansatz_type=config.get('ansatz', 'ucc'),
                shots=config.get('shots', 2048)
            )

            logger.info(f"Job submitted to BlueQubit: {result['job_id']}")

            return {
                'cloud_job_id': result['job_id'],
                'backend': 'bluequbit',
                'status': result['status']
            }

        except Exception as e:
            logger.error(f"BlueQubit submission failed: {e}")
            raise

    async def get_available_backends(
        self,
        provider: str,
        credentials: Dict[str, str]
    ) -> list:
        """
        Get list of available cloud backends.

        Args:
            provider: 'ibm_quantum' or 'bluequbit'
            credentials: Provider credentials

        Returns:
            List of backend information dicts
        """
        if provider == 'ibm_quantum':
            try:
                from qiskit_ibm_runtime import QiskitRuntimeService

                service = QiskitRuntimeService(
                    channel=credentials.get('channel', 'ibm_quantum'),
                    token=credentials['api_token'],
                    instance=credentials.get('crn')
                )

                backends = service.backends()

                return [
                    {
                        'name': backend.name,
                        'qubits': backend.num_qubits,
                        'online': backend.status().operational,
                        'queue_depth': backend.status().pending_jobs
                    }
                    for backend in backends
                ]

            except Exception as e:
                logger.error(f"Failed to get IBM backends: {e}")
                return []

        elif provider == 'bluequbit':
            # BlueQubit has one main backend
            return [
                {
                    'name': 'bluequbit_gpu',
                    'qubits': 36,
                    'speedup': '10x',
                    'online': True
                }
            ]

        return []

    async def check_job_status(
        self,
        cloud_job_id: str,
        provider: str,
        credentials: Dict[str, str]
    ) -> Dict[str, Any]:
        """
        Check status of cloud job.

        Args:
            cloud_job_id: Cloud provider's job ID
            provider: 'ibm_quantum' or 'bluequbit'
            credentials: Provider credentials

        Returns:
            Job status information
        """
        if provider == 'ibm_quantum':
            try:
                from qiskit_ibm_runtime import QiskitRuntimeService

                service = QiskitRuntimeService(
                    channel=credentials.get('channel', 'ibm_quantum'),
                    token=credentials['api_token'],
                    instance=credentials.get('crn')
                )

                job = service.job(cloud_job_id)

                return {
                    'status': job.status().value,
                    'queue_position': getattr(job.status(), 'queue_position', None)
                }

            except Exception as e:
                logger.error(f"Failed to check IBM job status: {e}")
                return {'status': 'unknown'}

        return {'status': 'unknown'}
