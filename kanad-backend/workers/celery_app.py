"""
Celery application for background task processing.

Handles computationally intensive quantum chemistry calculations
asynchronously.
"""

from celery import Celery
import os
import logging

logger = logging.getLogger(__name__)

# Get Redis URL from environment
REDIS_URL = os.getenv("REDIS_URL", "redis://localhost:6379/0")

# Create Celery app
celery_app = Celery(
    "kanad",
    broker=REDIS_URL,
    backend=REDIS_URL,
    include=['workers.tasks']
)

# Configure Celery
celery_app.conf.update(
    task_serializer='json',
    accept_content=['json'],
    result_serializer='json',
    timezone='UTC',
    enable_utc=True,
    task_track_started=True,
    task_time_limit=7200,  # 2 hours max per task
    task_soft_time_limit=6600,  # Soft limit 110 minutes
    worker_prefetch_multiplier=1,  # Don't prefetch tasks (important for long tasks)
    worker_max_tasks_per_child=50,  # Restart worker after 50 tasks (prevent memory leaks)
)

logger.info(f"Celery app configured with broker: {REDIS_URL}")
