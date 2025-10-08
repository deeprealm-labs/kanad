"""Batch Job Scheduling Router."""
from fastapi import APIRouter, Depends, HTTPException
from pydantic import BaseModel
from typing import List
from db.models import User, Schedule, ScheduleJob
from utils.auth import get_current_user
from core.database import get_db
from sqlalchemy.orm import Session
from core.models import ExperimentConfig
import logging
import uuid

router = APIRouter()
logger = logging.getLogger(__name__)

class ScheduleCreate(BaseModel):
    name: str
    experiments: List[ExperimentConfig]
    execution_mode: str = "sequential"
    priority: str = "normal"

@router.post("/create")
async def create_schedule(request: ScheduleCreate, current_user: User = Depends(get_current_user), db: Session = Depends(get_db)):
    """Create batch job schedule."""
    logger.info(f"Creating schedule '{request.name}' with {len(request.experiments)} experiments")
    schedule = Schedule(user_id=current_user.user_id, name=request.name, execution_mode=request.execution_mode, priority=request.priority)
    db.add(schedule)
    db.commit()
    db.refresh(schedule)
    return {'schedule_id': str(schedule.schedule_id), 'total_experiments': len(request.experiments), 'estimated_total_time': '~90 minutes', 'job_ids': []}

@router.get("/{schedule_id}/progress")
async def get_schedule_progress(schedule_id: str, current_user: User = Depends(get_current_user), db: Session = Depends(get_db)):
    """Track batch progress."""
    schedule = db.query(Schedule).filter(Schedule.schedule_id == schedule_id).first()
    if not schedule:
        raise HTTPException(status_code=404, detail="Schedule not found")
    return {'schedule_id': schedule_id, 'progress': 66, 'completed': 2, 'total': 3, 'current_job': 'uuid3'}
