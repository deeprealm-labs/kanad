"""
Experiment submission and management endpoints
"""

from fastapi import APIRouter, HTTPException, BackgroundTasks
from pydantic import BaseModel
from typing import List, Optional, Dict, Any
import uuid

from api.core.database import ExperimentDB, JobDB
from api.services.experiment_service import execute_experiment

router = APIRouter()


class AtomInput(BaseModel):
    symbol: str
    x: float
    y: float
    z: float


class MoleculeConfig(BaseModel):
    smiles: Optional[str] = None
    atoms: Optional[List[AtomInput]] = None
    basis: str = "sto-3g"
    charge: int = 0
    multiplicity: int = 1


class BackendConfig(BaseModel):
    method: str = "VQE"  # HF, VQE, SQD, etc.
    ansatz: Optional[str] = "hardware_efficient"
    mapper: Optional[str] = "jordan_wigner"
    optimizer: Optional[str] = "SLSQP"
    max_iterations: Optional[int] = 100  # Changed default from 1000 to 100
    backend: str = "classical"  # classical, ibm_quantum, bluequbit
    backend_name: Optional[str] = None
    bluequbit_device: Optional[str] = "gpu"  # BlueQubit device: cpu, gpu, mps.cpu, mps.gpu, pauli-path

    class Config:
        # Allow both camelCase (frontend) and snake_case (backend) field names
        populate_by_name = True
        alias_generator = lambda field_name: ''.join(
            word.capitalize() if i > 0 else word
            for i, word in enumerate(field_name.split('_'))
        )


class AnalysisConfig(BaseModel):
    energy_decomposition: bool = True
    bond_analysis: bool = True
    dipole_moment: bool = True
    polarizability: bool = False
    thermochemistry: bool = False
    spectroscopy: bool = False
    vibrational: bool = False


class OptimizationConfig(BaseModel):
    geometry: bool = False
    orbitals: bool = False
    circuit: bool = True
    adaptive: bool = False


class ExperimentSubmitRequest(BaseModel):
    molecule: MoleculeConfig
    configuration: BackendConfig
    analysis: Optional[AnalysisConfig] = AnalysisConfig()
    optimization: Optional[OptimizationConfig] = OptimizationConfig()
    execute_now: bool = True  # If False, add to queue


@router.post("/submit")
async def submit_experiment(
    request: ExperimentSubmitRequest,
    background_tasks: BackgroundTasks
):
    """
    Submit a new experiment for execution.

    Returns experiment ID and status.
    """
    try:
        # Generate unique ID
        experiment_id = str(uuid.uuid4())

        # Create experiment record
        experiment_data = {
            'id': experiment_id,
            'molecule': request.molecule.model_dump(),
            'configuration': {
                **request.configuration.model_dump(),
                'analysis': request.analysis.model_dump() if request.analysis else {},
                'optimization': request.optimization.model_dump() if request.optimization else {},
            },
            'status': 'queued',
            'method': request.configuration.method,
            'backend': request.configuration.backend,
        }

        ExperimentDB.create(experiment_data)

        # Create job
        job_id = str(uuid.uuid4())
        job_data = {
            'id': job_id,
            'experiment_id': experiment_id,
            'status': 'queued',
            'priority': 0,
            'max_iterations': request.configuration.max_iterations,
        }
        JobDB.create(job_data)

        # Execute immediately or add to queue
        if request.execute_now:
            # Execute in background
            background_tasks.add_task(
                execute_experiment,
                experiment_id,
                job_id
            )
            status = "running"
        else:
            status = "queued"

        return {
            "experiment_id": experiment_id,
            "job_id": job_id,
            "status": status,
            "message": f"Experiment {'started' if request.execute_now else 'queued'} successfully"
        }

    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/list")
async def list_experiments(
    status: Optional[str] = None,
    limit: int = 50,
    offset: int = 0
):
    """List experiments with optional filtering."""
    try:
        experiments = ExperimentDB.list(status=status, limit=limit, offset=offset)

        return {
            "experiments": experiments,
            "total": len(experiments),
            "limit": limit,
            "offset": offset
        }

    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/{experiment_id}")
async def get_experiment(experiment_id: str):
    """Get experiment details."""
    experiment = ExperimentDB.get(experiment_id)

    if not experiment:
        raise HTTPException(status_code=404, detail="Experiment not found")

    return {"experiment": experiment}


@router.get("/{experiment_id}/status")
async def get_experiment_status(experiment_id: str):
    """Get experiment status."""
    experiment = ExperimentDB.get(experiment_id)

    if not experiment:
        raise HTTPException(status_code=404, detail="Experiment not found")

    return {
        "experiment_id": experiment_id,
        "status": experiment['status'],
        "created_at": experiment['created_at'],
        "started_at": experiment['started_at'],
        "completed_at": experiment['completed_at'],
    }


@router.get("/{experiment_id}/results")
async def get_experiment_results(experiment_id: str):
    """Get experiment results."""
    experiment = ExperimentDB.get(experiment_id)

    if not experiment:
        raise HTTPException(status_code=404, detail="Experiment not found")

    if experiment['status'] != 'completed':
        raise HTTPException(
            status_code=400,
            detail=f"Experiment not completed (status: {experiment['status']})"
        )

    return {
        "experiment_id": experiment_id,
        "results": experiment['results'],
        "molecule": experiment['molecule'],
        "configuration": experiment['configuration'],
    }


@router.delete("/{experiment_id}")
async def delete_experiment(experiment_id: str):
    """Delete an experiment."""
    experiment = ExperimentDB.get(experiment_id)

    if not experiment:
        raise HTTPException(status_code=404, detail="Experiment not found")

    ExperimentDB.delete(experiment_id)

    return {"message": "Experiment deleted successfully"}


@router.post("/{experiment_id}/cancel")
async def cancel_experiment(experiment_id: str):
    """Cancel a running or queued experiment."""
    experiment = ExperimentDB.get(experiment_id)

    if not experiment:
        raise HTTPException(status_code=404, detail="Experiment not found")

    # If already cancelled, just return success
    if experiment['status'] == 'cancelled':
        return {"message": "Experiment already cancelled"}

    # Allow cancelling failed experiments (for cleanup)
    # Only reject if completed successfully
    if experiment['status'] == 'completed':
        raise HTTPException(
            status_code=400,
            detail="Cannot cancel a successfully completed experiment"
        )

    ExperimentDB.update_status(experiment_id, 'cancelled')

    # Also cancel associated job
    jobs = JobDB.list()
    for job in jobs:
        if job['experiment_id'] == experiment_id:
            JobDB.update_status(job['id'], 'cancelled')

    return {"message": "Experiment cancelled successfully"}
