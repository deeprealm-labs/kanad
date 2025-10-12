"""
Analysis endpoints
"""

from fastapi import APIRouter, HTTPException

router = APIRouter()


@router.get("/{experiment_id}/bond-analysis")
async def get_bond_analysis(experiment_id: str):
    """Get bond analysis for completed experiment."""
    # TODO: Implement bond analysis extraction from results
    return {
        "experiment_id": experiment_id,
        "message": "Bond analysis endpoint - to be implemented"
    }


@router.get("/{experiment_id}/energy-decomposition")
async def get_energy_decomposition(experiment_id: str):
    """Get energy decomposition."""
    # TODO: Implement energy decomposition
    return {
        "experiment_id": experiment_id,
        "message": "Energy decomposition endpoint - to be implemented"
    }
