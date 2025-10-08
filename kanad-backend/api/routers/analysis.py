"""Analysis router - Results and analysis endpoints."""
from fastapi import APIRouter
router = APIRouter()

@router.get("/{job_id}/results")
async def get_results(job_id: str):
    """Get job results with analysis."""
    return {"job_id": job_id, "status": "completed"}
