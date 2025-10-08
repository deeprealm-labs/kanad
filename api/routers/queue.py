"""
Job queue API router.

Handles job queue management, scheduling, and execution.
"""

from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy.orm import Session
from datetime import datetime
from typing import Optional

from api.database import get_db
from api.models.experiment import Experiment
from api.models.queue import QueueItem
from api.utils.validators import QueueItemCreate, QueueItemUpdate
from api.utils.exceptions import QueueItemNotFoundError, ExperimentNotFoundError
from api.services.job_queue import job_queue

router = APIRouter(prefix="/queue", tags=["queue"])


@router.post("/", status_code=201)
def add_to_queue(
    queue_data: QueueItemCreate,
    db: Session = Depends(get_db)
):
    """
    Add experiment to job queue.

    Optionally schedule for later execution.
    """
    # Verify experiment exists
    experiment = db.query(Experiment).filter(Experiment.id == queue_data.experiment_id).first()
    if not experiment:
        raise ExperimentNotFoundError(queue_data.experiment_id)

    # Parse scheduled time if provided
    scheduled_time = None
    if queue_data.scheduled_time:
        try:
            scheduled_time = datetime.fromisoformat(queue_data.scheduled_time)
        except ValueError:
            raise HTTPException(status_code=400, detail="Invalid datetime format. Use ISO 8601.")

    # Create queue item
    queue_item = QueueItem(
        experiment_id=queue_data.experiment_id,
        status="scheduled" if scheduled_time else "queued",
        priority=queue_data.priority,
        scheduled_time=scheduled_time
    )

    db.add(queue_item)

    # Update experiment status
    experiment.status = "scheduled" if scheduled_time else "queued"

    db.commit()
    db.refresh(queue_item)

    # Add to background queue if not scheduled
    if not scheduled_time:
        job_queue.add_job(queue_data.experiment_id, priority=queue_data.priority)

    return queue_item.to_dict()


@router.get("/")
def list_queue(
    status: Optional[str] = None,
    db: Session = Depends(get_db)
):
    """
    List all queue items.

    Optionally filter by status.
    """
    query = db.query(QueueItem)

    if status:
        query = query.filter(QueueItem.status == status)

    # Order by priority (descending) then creation time
    queue_items = query.order_by(
        QueueItem.priority.desc(),
        QueueItem.created_at.asc()
    ).all()

    # Join with experiments to get full info
    result = []
    for item in queue_items:
        experiment = db.query(Experiment).filter(Experiment.id == item.experiment_id).first()
        result.append({
            **item.to_dict(),
            "experiment": experiment.to_dict() if experiment else None
        })

    return {
        "total": len(result),
        "queue": result,
        "queue_size": job_queue.get_queue_size()
    }


@router.get("/{queue_id}")
def get_queue_item(
    queue_id: int,
    db: Session = Depends(get_db)
):
    """
    Get queue item details.
    """
    queue_item = db.query(QueueItem).filter(QueueItem.id == queue_id).first()

    if not queue_item:
        raise QueueItemNotFoundError(queue_id)

    # Get associated experiment
    experiment = db.query(Experiment).filter(Experiment.id == queue_item.experiment_id).first()

    return {
        **queue_item.to_dict(),
        "experiment": experiment.to_dict() if experiment else None
    }


@router.put("/{queue_id}")
def update_queue_item(
    queue_id: int,
    update_data: QueueItemUpdate,
    db: Session = Depends(get_db)
):
    """
    Update queue item.

    Supports updating status, priority, and scheduled time.
    """
    queue_item = db.query(QueueItem).filter(QueueItem.id == queue_id).first()

    if not queue_item:
        raise QueueItemNotFoundError(queue_id)

    # Update fields
    if update_data.status is not None:
        queue_item.status = update_data.status

        # Update associated experiment status
        experiment = db.query(Experiment).filter(Experiment.id == queue_item.experiment_id).first()
        if experiment:
            experiment.status = update_data.status

    if update_data.priority is not None:
        queue_item.priority = update_data.priority

    if update_data.scheduled_time is not None:
        try:
            queue_item.scheduled_time = datetime.fromisoformat(update_data.scheduled_time)
            queue_item.status = "scheduled"
        except ValueError:
            raise HTTPException(status_code=400, detail="Invalid datetime format. Use ISO 8601.")

    db.commit()
    db.refresh(queue_item)

    return queue_item.to_dict()


@router.delete("/{queue_id}")
def delete_queue_item(
    queue_id: int,
    db: Session = Depends(get_db)
):
    """
    Remove item from queue.

    Does not delete the experiment itself.
    """
    queue_item = db.query(QueueItem).filter(QueueItem.id == queue_id).first()

    if not queue_item:
        raise QueueItemNotFoundError(queue_id)

    # Update experiment status
    experiment = db.query(Experiment).filter(Experiment.id == queue_item.experiment_id).first()
    if experiment and experiment.status in ["queued", "scheduled"]:
        experiment.status = "cancelled"

    db.delete(queue_item)
    db.commit()

    return {
        "message": f"Queue item {queue_id} deleted successfully"
    }


@router.post("/{queue_id}/execute")
def execute_queue_item(
    queue_id: int,
    db: Session = Depends(get_db)
):
    """
    Manually trigger execution of a queued item.
    """
    queue_item = db.query(QueueItem).filter(QueueItem.id == queue_id).first()

    if not queue_item:
        raise QueueItemNotFoundError(queue_id)

    # Update status
    queue_item.status = "queued"
    db.commit()

    # Add to job queue
    job_queue.add_job(queue_item.experiment_id, priority=queue_item.priority)

    return {
        "message": f"Queue item {queue_id} added to execution queue"
    }


@router.get("/stats")
def get_queue_statistics(db: Session = Depends(get_db)):
    """
    Get comprehensive queue statistics.

    Returns counts of jobs by status, queue size, and performance metrics.
    """
    # Count by status
    total = db.query(QueueItem).count()
    queued = db.query(QueueItem).filter(QueueItem.status == "queued").count()
    running = db.query(QueueItem).filter(QueueItem.status == "running").count()
    completed = db.query(QueueItem).filter(QueueItem.status == "completed").count()
    failed = db.query(QueueItem).filter(QueueItem.status == "failed").count()
    scheduled = db.query(QueueItem).filter(QueueItem.status == "scheduled").count()

    # Get experiments statistics
    total_experiments = db.query(Experiment).count()
    completed_experiments = db.query(Experiment).filter(Experiment.status == "completed").count()
    failed_experiments = db.query(Experiment).filter(Experiment.status == "failed").count()

    # Calculate success rate
    finished_count = completed_experiments + failed_experiments
    success_rate = (completed_experiments / finished_count * 100) if finished_count > 0 else 0

    return {
        "queue": {
            "total": total,
            "queued": queued,
            "running": running,
            "completed": completed,
            "failed": failed,
            "scheduled": scheduled,
            "active_size": job_queue.get_queue_size()
        },
        "experiments": {
            "total": total_experiments,
            "completed": completed_experiments,
            "failed": failed_experiments,
            "success_rate": round(success_rate, 2)
        }
    }
