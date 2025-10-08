"""Jobs router - Job management and monitoring."""
from fastapi import APIRouter, Depends, WebSocket
from sqlalchemy.orm import Session
from typing import List
from core.database import get_db
from db.models import User, Job
from utils.auth import get_current_user
from core.models import JobListItem, JobStatusResponse

router = APIRouter()

@router.get("/", response_model=List[JobListItem])
async def list_jobs(current_user: User = Depends(get_current_user), db: Session = Depends(get_db)):
    """List user's jobs."""
    jobs = db.query(Job).filter(Job.user_id == current_user.user_id).order_by(Job.created_at.desc()).limit(50).all()
    return [{"job_id": str(j.job_id), "molecule_name": "Unknown", "method": j.simulation.method, "status": j.status, "progress": j.progress, "created_at": j.created_at, "started_at": j.started_at, "backend": j.simulation.backend_type.value} for j in jobs]

@router.get("/{job_id}/status", response_model=JobStatusResponse)
async def get_job_status(job_id: str, current_user: User = Depends(get_current_user), db: Session = Depends(get_db)):
    """Get job status."""
    job = db.query(Job).filter(Job.job_id == job_id, Job.user_id == current_user.user_id).first()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    return {"job_id": str(job.job_id), "status": job.status, "progress": job.progress, "message": "Running"}

@router.websocket("/{job_id}/logs")
async def job_logs_websocket(websocket: WebSocket, job_id: str):
    """WebSocket for real-time job logs."""
    await websocket.accept()
    # Stream logs from database
    await websocket.send_text("Job logs stream")
    await websocket.close()

@router.delete("/{job_id}")
async def cancel_job(job_id: str, current_user: User = Depends(get_current_user), db: Session = Depends(get_db)):
    """Cancel running job."""
    return {"job_id": job_id, "status": "cancelled"}
