"""
Job queue management endpoints
"""

from fastapi import APIRouter, HTTPException
from typing import Optional

from api.core.database import JobDB, ExperimentDB

router = APIRouter()


@router.get("/list")
async def list_jobs(status: Optional[str] = None):
    """List all jobs in the queue."""
    try:
        jobs = JobDB.list(status=status)

        # Enrich with experiment data
        enriched_jobs = []
        for job in jobs:
            experiment = ExperimentDB.get(job['experiment_id'])
            if experiment:
                enriched_jobs.append({
                    **job,
                    'molecule': experiment['molecule'],
                    'method': experiment['method'],
                    'backend': experiment['backend'],
                })

        return {
            "jobs": enriched_jobs,
            "total": len(enriched_jobs)
        }

    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/{job_id}")
async def get_job(job_id: str):
    """Get job details."""
    job = JobDB.get(job_id)

    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    # Get associated experiment
    experiment = ExperimentDB.get(job['experiment_id'])

    return {
        "job": job,
        "experiment": experiment
    }


@router.get("/{job_id}/status")
async def get_job_status(job_id: str):
    """Get job status and progress."""
    job = JobDB.get(job_id)

    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    return {
        "job_id": job_id,
        "status": job['status'],
        "progress": job['progress'],
        "current_iteration": job['current_iteration'],
        "max_iterations": job['max_iterations'],
        "current_energy": job['current_energy'],
        "best_energy": job['best_energy'],
    }


@router.post("/{job_id}/cancel")
async def cancel_job(job_id: str):
    """Cancel a job."""
    job = JobDB.get(job_id)

    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    if job['status'] in ['completed', 'failed', 'cancelled']:
        raise HTTPException(
            status_code=400,
            detail=f"Cannot cancel job with status: {job['status']}"
        )

    JobDB.update_status(job_id, 'cancelled')
    ExperimentDB.update_status(job['experiment_id'], 'cancelled')

    return {"message": "Job cancelled successfully"}


@router.delete("/{job_id}")
async def delete_job(job_id: str):
    """Delete a job from queue."""
    job = JobDB.get(job_id)

    if not job:
        raise HTTPException(status_code=404, detail="Job not found")

    JobDB.delete(job_id)

    return {"message": "Job deleted successfully"}
