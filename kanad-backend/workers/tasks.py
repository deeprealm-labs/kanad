"""
Celery background tasks for quantum chemistry computations.

These tasks run in worker processes and handle long-running calculations.
"""

from workers.celery_app import celery_app
from services.computation_service import ComputationService
from core.database import SessionLocal
from db.models import Job, Result, JobLog, JobStatusEnum
from datetime import datetime
import logging
import asyncio
import traceback
import json

logger = logging.getLogger(__name__)


def log_to_db(job_id: str, message: str, level: str = 'INFO'):
    """Log message to database for WebSocket streaming."""
    db = SessionLocal()
    try:
        log = JobLog(
            job_id=job_id,
            message=message,
            level=level,
            timestamp=datetime.utcnow()
        )
        db.add(log)
        db.commit()
    except Exception as e:
        logger.error(f"Failed to log to DB: {e}")
    finally:
        db.close()


@celery_app.task(bind=True, name='workers.tasks.run_vqe_computation')
def run_vqe_computation(self, job_id: str, molecule_data: dict, config: dict):
    """
    Run VQE computation as background task.

    Args:
        self: Celery task instance (for state updates)
        job_id: Database job ID
        molecule_data: Serialized molecule data
        config: Computation configuration

    Returns:
        Results dictionary
    """
    db = SessionLocal()

    try:
        # Update job status to RUNNING
        job = db.query(Job).filter(Job.job_id == job_id).first()
        job.status = JobStatusEnum.RUNNING
        job.started_at = datetime.utcnow()
        db.commit()

        log_to_db(job_id, "Job started - Initializing computation service", "INFO")

        # Initialize computation service
        comp_service = ComputationService()

        # Create progress callback
        async def progress_callback(progress: int, message: str):
            """Update task state and log progress."""
            self.update_state(
                state='PROGRESS',
                meta={'progress': progress, 'message': message}
            )
            job.progress = progress
            db.commit()
            log_to_db(job_id, message, "INFO")

        # Create molecule
        log_to_db(job_id, f"Creating molecule with {len(molecule_data['atoms'])} atoms", "INFO")

        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)

        molecule_result = loop.run_until_complete(
            comp_service.create_molecule_from_atoms(
                atoms=molecule_data['atoms'],
                basis=config.get('basis', 'sto-3g'),
                charge=molecule_data.get('charge', 0),
                multiplicity=molecule_data.get('multiplicity', 1)
            )
        )

        molecule = molecule_result['molecule']
        log_to_db(job_id, f"Molecule created: {molecule.formula} ({molecule.n_electrons} electrons)", "INFO")

        # Run computation
        log_to_db(job_id, f"Starting {config['method']} computation...", "INFO")

        computation_result = loop.run_until_complete(
            comp_service.run_computation(
                molecule=molecule,
                method=config['method'],
                config=config,
                progress_callback=progress_callback
            )
        )

        log_to_db(job_id, f"Computation complete: E = {computation_result['energy']:.6f} Ha", "INFO")

        # Run analysis if requested
        analysis_result = None
        if config.get('analysis'):
            log_to_db(job_id, "Running analyses...", "INFO")
            analysis_result = loop.run_until_complete(
                comp_service.run_analysis(
                    molecule=molecule,
                    results=computation_result,
                    analysis_requests=config['analysis']
                )
            )
            log_to_db(job_id, "Analysis complete", "INFO")

        loop.close()

        # Save results to database
        result = Result(
            job_id=job_id,
            energy=computation_result['energy'],
            hf_energy=computation_result.get('hf_energy'),
            correlation_energy=computation_result.get('correlation_energy'),
            n_iterations=computation_result.get('n_iterations'),
            converged=computation_result.get('converged'),
            convergence_history=computation_result.get('convergence_history'),
            analysis=analysis_result
        )
        db.add(result)

        # Update job status
        job.status = JobStatusEnum.COMPLETED
        job.completed_at = datetime.utcnow()
        job.progress = 100
        db.commit()

        log_to_db(job_id, "Job completed successfully!", "INFO")

        return {
            'status': 'completed',
            'energy': computation_result['energy'],
            'job_id': job_id
        }

    except Exception as e:
        # Log error
        error_msg = f"Computation failed: {str(e)}\n{traceback.format_exc()}"
        logger.error(error_msg)
        log_to_db(job_id, error_msg, "ERROR")

        # Update job status
        try:
            job = db.query(Job).filter(Job.job_id == job_id).first()
            job.status = JobStatusEnum.FAILED
            job.error_message = str(e)
            db.commit()
        except Exception as db_error:
            logger.error(f"Failed to update job status: {db_error}")

        # Re-raise for Celery to mark as failed
        raise

    finally:
        db.close()


@celery_app.task(bind=True, name='workers.tasks.run_cloud_computation')
def run_cloud_computation(self, job_id: str, molecule_data: dict, config: dict, cloud_credentials: dict):
    """
    Run computation on cloud backend (IBM Quantum or BlueQubit).

    Args:
        self: Celery task instance
        job_id: Database job ID
        molecule_data: Serialized molecule
        config: Configuration including backend details
        cloud_credentials: Encrypted cloud credentials

    Returns:
        Cloud job information
    """
    db = SessionLocal()

    try:
        job = db.query(Job).filter(Job.job_id == job_id).first()
        job.status = JobStatusEnum.RUNNING
        job.started_at = datetime.utcnow()
        db.commit()

        log_to_db(job_id, "Preparing cloud computation...", "INFO")

        # Initialize services
        comp_service = ComputationService()

        from services.cloud_service import CloudService
        cloud_service = CloudService()

        # Create molecule
        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)

        molecule_result = loop.run_until_complete(
            comp_service.create_molecule_from_atoms(
                atoms=molecule_data['atoms'],
                basis=config.get('basis', 'sto-3g')
            )
        )

        molecule = molecule_result['molecule']

        # Submit to cloud
        backend_type = config['backend']['type']
        log_to_db(job_id, f"Submitting to {backend_type} backend...", "INFO")

        if backend_type == 'ibm_quantum':
            cloud_result = loop.run_until_complete(
                cloud_service.submit_to_ibm(
                    molecule=molecule,
                    config=config,
                    credentials=cloud_credentials
                )
            )
        elif backend_type == 'bluequbit':
            cloud_result = loop.run_until_complete(
                cloud_service.submit_to_bluequbit(
                    molecule=molecule,
                    config=config,
                    credentials=cloud_credentials
                )
            )
        else:
            raise ValueError(f"Unknown backend: {backend_type}")

        loop.close()

        # Store cloud job ID
        job.cloud_job_id = cloud_result['cloud_job_id']
        db.commit()

        log_to_db(job_id, f"Job submitted to cloud: {cloud_result['cloud_job_id']}", "INFO")
        log_to_db(job_id, "Waiting for cloud job to complete...", "INFO")

        return {
            'status': 'submitted',
            'cloud_job_id': cloud_result['cloud_job_id'],
            'backend': backend_type
        }

    except Exception as e:
        error_msg = f"Cloud submission failed: {str(e)}\n{traceback.format_exc()}"
        logger.error(error_msg)
        log_to_db(job_id, error_msg, "ERROR")

        job = db.query(Job).filter(Job.job_id == job_id).first()
        job.status = JobStatusEnum.FAILED
        job.error_message = str(e)
        db.commit()

        raise

    finally:
        db.close()
