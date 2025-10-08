"""Simulations router - Computation configuration and job submission."""
from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy.orm import Session
from core.database import get_db
from db.models import User
from utils.auth import get_current_user
from core.models import SimulationConfig, SimulationPreview, SimulationSubmit
import uuid

router = APIRouter()

@router.post("/configure", response_model=SimulationPreview)
async def configure_simulation(config: SimulationConfig, current_user: User = Depends(get_current_user), db: Session = Depends(get_db)):
    """Configure simulation and get preview."""
    # Implementation here - returns circuit preview and T&C
    return {
        "simulation_id": str(uuid.uuid4()),
        "status": "configured",
        "preview": {"n_qubits": 14, "n_parameters": 48, "circuit_depth": 24, "estimated_time": "~30s", "estimated_cost": "$0.00"},
        "terms_and_conditions": {"data_usage": "Results stored securely"}
    }

@router.post("/{simulation_id}/accept-and-run")
async def accept_and_run(simulation_id: str, submission: SimulationSubmit, current_user: User = Depends(get_current_user), db: Session = Depends(get_db)):
    """Accept T&C and submit job."""
    # Submit job to Celery
    return {"job_id": str(uuid.uuid4()), "status": "queued"}
