"""User Profile and History Router."""
from fastapi import APIRouter, Depends
from sqlalchemy.orm import Session
from db.models import User, Job
from utils.auth import get_current_user
from core.database import get_db
import logging

router = APIRouter()
logger = logging.getLogger(__name__)

@router.get("/profile")
async def get_user_profile(current_user: User = Depends(get_current_user), db: Session = Depends(get_db)):
    """Get user profile."""
    total_jobs = db.query(Job).filter(Job.user_id == current_user.user_id).count()
    completed_jobs = db.query(Job).filter(Job.user_id == current_user.user_id, Job.status == 'completed').count()
    return {
        'user_id': str(current_user.user_id),
        'email': current_user.email,
        'name': current_user.name,
        'institution': current_user.institution,
        'field': current_user.research_field,
        'stats': {'total_jobs': total_jobs, 'completed_jobs': completed_jobs, 'total_molecules': 0}
    }

@router.get("/history")
async def get_user_history(limit: int = 20, current_user: User = Depends(get_current_user), db: Session = Depends(get_db)):
    """Job history with filters."""
    jobs = db.query(Job).filter(Job.user_id == current_user.user_id).order_by(Job.created_at.desc()).limit(limit).all()
    return {'jobs': [{'job_id': str(j.job_id), 'method': j.simulation.method, 'status': j.status} for j in jobs], 'pagination': {'page': 1, 'total_pages': 1}}
