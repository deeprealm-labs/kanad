"""
Analysis endpoints
"""

from fastapi import APIRouter, HTTPException, BackgroundTasks
from typing import List, Dict, Any, Optional
from pydantic import BaseModel
import logging

from api.core.database import ExperimentDB, AnalysisDB
from kanad.services.analysis_service import AnalysisService

logger = logging.getLogger(__name__)

router = APIRouter()


class CompareExperimentsRequest(BaseModel):
    experiment_ids: List[str]


class RunAnalysisRequest(BaseModel):
    profile: Optional[str] = None
    analyses: Optional[List[str]] = None
    parameters: Optional[Dict[str, Any]] = None
    async_execution: bool = False


@router.get("/{experiment_id}/bond-analysis")
async def get_bond_analysis(experiment_id: str):
    """Get bond analysis for completed experiment."""
    experiment = ExperimentDB.get(experiment_id)

    if not experiment:
        raise HTTPException(status_code=404, detail="Experiment not found")

    if experiment['status'] != 'completed':
        raise HTTPException(
            status_code=400,
            detail=f"Experiment not completed (status: {experiment['status']})"
        )

    results = experiment.get('results', {})
    analysis = results.get('analysis', {})

    # Extract bond-specific analysis
    bond_data = {
        "experiment_id": experiment_id,
        "bond_order": analysis.get('bond_order'),
        "bond_length": analysis.get('bond_length'),
        "bond_type": analysis.get('bond_type'),
        "bond_energy": analysis.get('bond_energy'),
        "dissociation_energy": analysis.get('dissociation_energy'),
        "vibrational_frequency": analysis.get('vibrational_frequency'),
        "force_constant": analysis.get('force_constant'),
        # Include all bond-related keys from analysis
        **{k: v for k, v in analysis.items() if 'bond' in k.lower()}
    }

    # Remove None values
    bond_data = {k: v for k, v in bond_data.items() if v is not None}

    if len(bond_data) <= 1:  # Only experiment_id present
        return {
            "experiment_id": experiment_id,
            "message": "No bond analysis data available for this experiment"
        }

    return bond_data


@router.get("/{experiment_id}/energy-decomposition")
async def get_energy_decomposition(experiment_id: str):
    """Get energy decomposition."""
    experiment = ExperimentDB.get(experiment_id)

    if not experiment:
        raise HTTPException(status_code=404, detail="Experiment not found")

    if experiment['status'] != 'completed':
        raise HTTPException(
            status_code=400,
            detail=f"Experiment not completed (status: {experiment['status']})"
        )

    results = experiment.get('results', {})
    analysis = results.get('analysis', {})

    # Extract energy decomposition
    energy_data = {
        "experiment_id": experiment_id,
        "total_energy": results.get('energy'),
        "hf_energy": results.get('hf_energy'),
        "correlation_energy": results.get('correlation_energy'),
        "nuclear_repulsion": analysis.get('nuclear_repulsion_energy'),
        "kinetic_energy": analysis.get('kinetic_energy'),
        "potential_energy": analysis.get('potential_energy'),
        "coulomb_energy": analysis.get('coulomb_energy'),
        "exchange_energy": analysis.get('exchange_energy'),
        "one_electron_energy": analysis.get('one_electron_energy'),
        "two_electron_energy": analysis.get('two_electron_energy'),
        # Include all energy-related keys
        **{k: v for k, v in analysis.items() if 'energy' in k.lower()}
    }

    # Remove None values and duplicates
    seen = set()
    energy_data = {k: v for k, v in energy_data.items()
                   if v is not None and not (k in seen or seen.add(k))}

    if len(energy_data) <= 1:  # Only experiment_id present
        return {
            "experiment_id": experiment_id,
            "message": "No energy decomposition data available"
        }

    return energy_data


@router.get("/{experiment_id}/full-analysis")
async def get_full_analysis(experiment_id: str):
    """Get complete analysis data for an experiment."""
    experiment = ExperimentDB.get(experiment_id)

    if not experiment:
        raise HTTPException(status_code=404, detail="Experiment not found")

    if experiment['status'] != 'completed':
        raise HTTPException(
            status_code=400,
            detail=f"Experiment not completed (status: {experiment['status']})"
        )

    results = experiment.get('results', {})
    analysis = results.get('analysis', {})

    if not analysis:
        return {
            "experiment_id": experiment_id,
            "message": "No analysis data available for this experiment"
        }

    return {
        "experiment_id": experiment_id,
        "analysis": analysis,
        "energy": results.get('energy'),
        "method": experiment.get('method'),
        "molecule": experiment.get('molecule')
    }


@router.post("/compare")
async def compare_experiments(request: CompareExperimentsRequest):
    """Compare multiple experiments side-by-side."""
    if len(request.experiment_ids) < 2:
        raise HTTPException(
            status_code=400,
            detail="At least 2 experiments required for comparison"
        )

    if len(request.experiment_ids) > 10:
        raise HTTPException(
            status_code=400,
            detail="Maximum 10 experiments can be compared at once"
        )

    experiments_data = []

    for exp_id in request.experiment_ids:
        experiment = ExperimentDB.get(exp_id)

        if not experiment:
            raise HTTPException(
                status_code=404,
                detail=f"Experiment {exp_id} not found"
            )

        results = experiment.get('results', {})

        # Extract key comparison data
        exp_data = {
            "id": exp_id,
            "status": experiment.get('status'),
            "method": experiment.get('method'),
            "backend": experiment.get('backend'),
            "molecule": experiment.get('molecule'),
            "energy": results.get('energy'),
            "hf_energy": results.get('hf_energy'),
            "correlation_energy": results.get('correlation_energy'),
            "iterations": results.get('iterations'),
            "converged": results.get('converged'),
            "convergence_history": results.get('convergence_history', []),
            "energy_history": results.get('energy_history', []),
            "analysis": results.get('analysis'),
            "created_at": experiment.get('created_at'),
            "completed_at": experiment.get('completed_at'),
            "configuration": experiment.get('configuration')
        }

        experiments_data.append(exp_data)

    # Calculate comparison statistics
    completed_experiments = [e for e in experiments_data if e['status'] == 'completed']

    comparison_stats = {
        "total_compared": len(request.experiment_ids),
        "completed": len(completed_experiments),
        "failed": len([e for e in experiments_data if e['status'] == 'failed'])
    }

    if completed_experiments:
        energies = [e['energy'] for e in completed_experiments if e['energy'] is not None]
        if energies:
            comparison_stats["lowest_energy"] = min(energies)
            comparison_stats["highest_energy"] = max(energies)
            comparison_stats["energy_range"] = max(energies) - min(energies)
            comparison_stats["best_experiment_id"] = completed_experiments[
                energies.index(min(energies))
            ]["id"]

    return {
        "experiments": experiments_data,
        "comparison_stats": comparison_stats
    }


# ============================================================================
# New Analysis Service Endpoints
# ============================================================================

@router.get("/profiles")
async def get_analysis_profiles():
    """Get available analysis profiles and their descriptions."""
    try:
        # Create a dummy service to access profile registry
        dummy_service = AnalysisService(experiment_results={'final_energy': 0.0})
        profiles = dummy_service.get_available_profiles()

        return {
            "profiles": profiles,
            "total": len(profiles)
        }
    except Exception as e:
        logger.error(f"Failed to get analysis profiles: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/modules")
async def get_analysis_modules():
    """Get available analysis modules and their requirements."""
    try:
        # Create a dummy service to access module registry
        dummy_service = AnalysisService(experiment_results={'final_energy': 0.0})
        modules = dummy_service.get_available_analyses()

        return {
            "modules": modules,
            "total": len(modules)
        }
    except Exception as e:
        logger.error(f"Failed to get analysis modules: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/{experiment_id}/run-analysis")
async def run_analysis_on_experiment(
    experiment_id: str,
    request: RunAnalysisRequest,
    background_tasks: BackgroundTasks
):
    """
    Run on-demand analysis on a completed experiment.

    Can run either:
    - A domain-specific profile (Chemistry, Spectroscopy, Materials, etc.)
    - Custom selection of individual analyses

    Supports both synchronous and asynchronous execution.
    """
    # Get experiment
    experiment = ExperimentDB.get(experiment_id)
    if not experiment:
        raise HTTPException(status_code=404, detail="Experiment not found")

    if experiment['status'] != 'completed':
        raise HTTPException(
            status_code=400,
            detail=f"Can only analyze completed experiments (status: {experiment['status']})"
        )

    # Validate request
    if not request.profile and not request.analyses:
        raise HTTPException(
            status_code=400,
            detail="Must specify either 'profile' or 'analyses'"
        )

    if request.profile and request.analyses:
        raise HTTPException(
            status_code=400,
            detail="Cannot specify both 'profile' and 'analyses' - choose one"
        )

    try:
        # Create analysis service
        service = AnalysisService.from_experiment_id(experiment_id)

        if request.async_execution:
            # Run in background
            # TODO: Implement async execution with WebSocket progress updates
            raise HTTPException(
                status_code=501,
                detail="Async execution not yet implemented"
            )
        else:
            # Run synchronously
            if request.profile:
                logger.info(f"Running {request.profile} profile on experiment {experiment_id}")
                results = service.run_analysis_profile(
                    request.profile,
                    request.parameters
                )
            else:
                # Run custom analyses
                logger.info(f"Running custom analyses {request.analyses} on experiment {experiment_id}")
                results = {
                    'profile': 'custom',
                    'analyses': {}
                }
                for analysis_name in request.analyses:
                    params = request.parameters.get(analysis_name, {}) if request.parameters else {}
                    result = service.run_analysis(analysis_name, params)
                    results['analyses'][analysis_name] = result

            # Save results to database
            analysis_data = {
                'experiment_id': experiment_id,
                'profile': request.profile or 'custom',
                'analyses': request.analyses or list(results['analyses'].keys()),
                'parameters': request.parameters or {},
                'results': results,
                'status': 'completed'
            }

            # Calculate total computation time
            computation_times = [
                r.get('computation_time', 0)
                for r in results['analyses'].values()
                if isinstance(r, dict)
            ]
            analysis_data['computation_time'] = sum(computation_times)

            analysis_id = AnalysisDB.create(analysis_data)

            return {
                "analysis_id": analysis_id,
                "experiment_id": experiment_id,
                "profile": request.profile or 'custom',
                "status": "completed",
                "computation_time": analysis_data['computation_time'],
                "results": results
            }

    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.error(f"Analysis failed for experiment {experiment_id}: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/results/{analysis_id}")
async def get_analysis_result(analysis_id: str):
    """Get analysis result by ID."""
    result = AnalysisDB.get(analysis_id)

    if not result:
        raise HTTPException(status_code=404, detail="Analysis result not found")

    return result


@router.get("/{experiment_id}/analyses")
async def get_experiment_analyses(
    experiment_id: str,
    profile: Optional[str] = None
):
    """Get all analysis results for an experiment, optionally filtered by profile."""
    # Verify experiment exists
    experiment = ExperimentDB.get(experiment_id)
    if not experiment:
        raise HTTPException(status_code=404, detail="Experiment not found")

    # Get analyses
    analyses = AnalysisDB.get_by_experiment(experiment_id, profile)

    return {
        "experiment_id": experiment_id,
        "analyses": analyses,
        "total": len(analyses)
    }


@router.delete("/results/{analysis_id}")
async def delete_analysis_result(analysis_id: str):
    """Delete an analysis result."""
    result = AnalysisDB.get(analysis_id)
    if not result:
        raise HTTPException(status_code=404, detail="Analysis result not found")

    AnalysisDB.delete(analysis_id)

    return {"message": "Analysis result deleted successfully"}
