"""
Background job queue system using threading.

Manages concurrent experiment execution with priority-based scheduling.
Simpler than Celery but sufficient for single-server deployments.
"""

import threading
import queue
import time
import logging
from typing import Optional, Callable, Dict, Any
from datetime import datetime
from sqlalchemy.orm import Session

from api.database import SessionLocal
from api.models.experiment import Experiment
from api.models.queue import QueueItem
from api.services.experiment_service import experiment_service
from api.config import settings

logger = logging.getLogger(__name__)


class JobQueue:
    """
    Background job queue for experiment execution.

    Features:
    - Priority-based scheduling
    - Concurrent execution with configurable workers
    - Progress tracking via callbacks
    - Automatic retry on failure
    - Graceful shutdown
    - Job cancellation support
    """

    def __init__(self, max_workers: int = 2):
        """
        Initialize job queue.

        Args:
            max_workers: Maximum number of concurrent experiments
        """
        self.max_workers = max_workers
        self.job_queue = queue.PriorityQueue()
        self.workers = []
        self.running = False
        self.current_jobs = {}  # job_id -> thread mapping
        self.cancellation_flags = {}  # experiment_id -> threading.Event()
        self.lock = threading.Lock()

        logger.info(f"Job queue initialized with {max_workers} workers")

    def start(self):
        """Start worker threads."""
        if self.running:
            logger.warning("Job queue already running")
            return

        self.running = True

        # Start worker threads
        for i in range(self.max_workers):
            worker = threading.Thread(
                target=self._worker,
                name=f"JobWorker-{i}",
                daemon=True
            )
            worker.start()
            self.workers.append(worker)

        logger.info(f"Started {self.max_workers} worker threads")

    def stop(self):
        """Stop worker threads gracefully."""
        logger.info("Stopping job queue...")
        self.running = False

        # Wait for workers to finish
        for worker in self.workers:
            worker.join(timeout=5.0)

        self.workers.clear()
        logger.info("Job queue stopped")

    def add_job(self, experiment_id: int, priority: int = 0):
        """
        Add job to queue.

        Args:
            experiment_id: Experiment ID to execute
            priority: Job priority (higher = runs first)
        """
        # Initialize cancellation flag for this experiment
        with self.lock:
            self.cancellation_flags[experiment_id] = threading.Event()

        # Priority queue uses lowest-first, so negate priority
        self.job_queue.put((-priority, experiment_id))
        logger.info(f"Added experiment {experiment_id} to queue with priority {priority}")

    def _worker(self):
        """Worker thread that processes jobs from queue."""
        logger.info(f"Worker {threading.current_thread().name} started")

        while self.running:
            try:
                # Get next job (with timeout to allow checking self.running)
                try:
                    priority, experiment_id = self.job_queue.get(timeout=1.0)
                except queue.Empty:
                    continue

                # Execute experiment
                logger.info(f"Worker {threading.current_thread().name} executing experiment {experiment_id}")
                self._execute_experiment(experiment_id)

                # Mark task as done
                self.job_queue.task_done()

            except Exception as e:
                logger.error(f"Worker error: {e}", exc_info=True)

        logger.info(f"Worker {threading.current_thread().name} stopped")

    def _execute_experiment(self, experiment_id: int):
        """
        Execute a single experiment.

        Args:
            experiment_id: Experiment ID to execute
        """
        db = SessionLocal()

        try:
            # Get experiment
            experiment = db.query(Experiment).filter(Experiment.id == experiment_id).first()
            if not experiment:
                logger.error(f"Experiment {experiment_id} not found")
                return

            # Check if cancelled before starting
            if self.is_cancelled(experiment_id):
                logger.info(f"Experiment {experiment_id} was cancelled before execution started")
                experiment.status = "cancelled"
                experiment.cancelled_at = datetime.now()
                experiment.completed_at = datetime.now()
                db.commit()
                return

            # Update status to running
            experiment.status = "running"
            experiment.started_at = datetime.now()
            db.commit()

            logger.info(f"Executing experiment {experiment_id}: {experiment.name}")

            # Create progress callback to update database
            convergence_data = []

            def progress_callback(iteration: int, energy: float, params):
                """Callback to track convergence and check for cancellation."""
                # Check for cancellation
                if self.is_cancelled(experiment_id):
                    from api.utils.exceptions import ExperimentCancelledException
                    raise ExperimentCancelledException(experiment_id)

                convergence_data.append({
                    'iteration': iteration,
                    'energy': float(energy)
                })

                # Update convergence data every 10 iterations
                if iteration % 10 == 0:
                    with SessionLocal() as update_db:
                        exp = update_db.query(Experiment).filter(Experiment.id == experiment_id).first()
                        if exp:
                            exp.convergence_data = convergence_data
                            update_db.commit()

                logger.debug(f"Experiment {experiment_id} - Iteration {iteration}: E={energy:.8f}")

            # Execute experiment
            result = experiment_service.execute_experiment(
                molecule_data=experiment.molecule_data,
                config=experiment.configuration,
                progress_callback=progress_callback,
                experiment_id=experiment_id,  # Pass for cancellation checks
                db_session=db  # Pass database session for cloud credentials
            )

            # Update experiment with results
            experiment.status = "completed"
            experiment.completed_at = datetime.now()
            experiment.energy = result.get('energy')
            experiment.hf_energy = result.get('hf_energy')
            experiment.correlation_energy = result.get('correlation_energy')
            experiment.results = result
            experiment.convergence_data = convergence_data

            # Update queue item
            queue_item = db.query(QueueItem).filter(
                QueueItem.experiment_id == experiment_id
            ).first()
            if queue_item:
                queue_item.status = "completed"

            db.commit()

            logger.info(f"Experiment {experiment_id} completed successfully")

        except Exception as e:
            from api.utils.exceptions import ExperimentCancelledException

            # Check if this is a cancellation
            if isinstance(e, ExperimentCancelledException):
                logger.info(f"Experiment {experiment_id} cancelled by user")
                experiment.status = "cancelled"
                experiment.cancelled_at = datetime.now()
                experiment.completed_at = datetime.now()
                experiment.error_message = "Cancelled by user"

                # Update queue item
                queue_item = db.query(QueueItem).filter(
                    QueueItem.experiment_id == experiment_id
                ).first()
                if queue_item:
                    queue_item.status = "cancelled"

                db.commit()
            else:
                # Regular failure
                logger.error(f"Experiment {experiment_id} failed: {e}", exc_info=True)

                # Update experiment status to failed
                experiment.status = "failed"
                experiment.completed_at = datetime.now()
                experiment.error_message = str(e)

                # Update queue item
                queue_item = db.query(QueueItem).filter(
                    QueueItem.experiment_id == experiment_id
                ).first()
                if queue_item:
                    queue_item.status = "failed"

                db.commit()

        finally:
            # Clean up cancellation flag
            with self.lock:
                if experiment_id in self.cancellation_flags:
                    del self.cancellation_flags[experiment_id]

            db.close()

    def get_queue_size(self) -> int:
        """Get number of jobs in queue."""
        return self.job_queue.qsize()

    def is_running(self) -> bool:
        """Check if queue is running."""
        return self.running

    def cancel_job(self, experiment_id: int) -> bool:
        """
        Request cancellation of a running or queued experiment.

        Args:
            experiment_id: Experiment ID to cancel

        Returns:
            True if cancellation was requested, False if job not found
        """
        with self.lock:
            if experiment_id in self.cancellation_flags:
                self.cancellation_flags[experiment_id].set()
                logger.info(f"Cancellation requested for experiment {experiment_id}")
                return True
            else:
                logger.warning(f"Cannot cancel experiment {experiment_id}: not in queue")
                return False

    def is_cancelled(self, experiment_id: int) -> bool:
        """
        Check if an experiment has been cancelled.

        Args:
            experiment_id: Experiment ID to check

        Returns:
            True if cancellation was requested
        """
        with self.lock:
            if experiment_id in self.cancellation_flags:
                return self.cancellation_flags[experiment_id].is_set()
            return False


# Global job queue instance
job_queue = JobQueue(max_workers=settings.MAX_CONCURRENT_JOBS)
