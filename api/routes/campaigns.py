"""
Campaign/batch execution endpoints
"""

from fastapi import APIRouter, HTTPException, BackgroundTasks, Depends
from pydantic import BaseModel
from typing import List, Optional
import threading
import uuid

from api.core.database import CampaignDB, ExperimentDB, JobDB
from api.services.experiment_service import execute_experiment
from api.dependencies.auth import get_optional_user, get_current_verified_user
from api.core.database_postgres import User

router = APIRouter()

# Semaphore for parallel campaign execution (max 3 campaigns at once)
MAX_CONCURRENT_CAMPAIGNS = 3
campaign_semaphore = threading.Semaphore(MAX_CONCURRENT_CAMPAIGNS)

# Lock for individual campaign coordination
campaign_locks = {}
campaign_locks_lock = threading.Lock()


class CreateCampaignRequest(BaseModel):
    name: str
    description: Optional[str] = ""
    experiment_ids: List[str]


class ReorderQueueRequest(BaseModel):
    experiment_id: str
    new_order: int


@router.post("/create")
async def create_campaign(
    request: CreateCampaignRequest,
    current_user: User = Depends(get_current_verified_user)
):
    """
    Create a campaign from queued experiments.
    Requires authentication.
    """
    try:
        # Create campaign
        campaign_id = CampaignDB.create({
            'name': request.name,
            'description': request.description,
            'total_experiments': len(request.experiment_ids),
            'user_id': current_user.id  # Associate campaign with authenticated user
        })

        # Update experiments to be part of this campaign
        for idx, exp_id in enumerate(request.experiment_ids):
            experiment = ExperimentDB.get(exp_id)
            if not experiment:
                continue

            # Update experiment with campaign_id and sequence order
            ExperimentDB.update_campaign_info(exp_id, campaign_id, idx)

        return {
            "campaign_id": campaign_id,
            "message": f"Campaign created with {len(request.experiment_ids)} experiments"
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/{campaign_id}/execute")
async def execute_campaign(campaign_id: str, background_tasks: BackgroundTasks):
    """Execute all experiments in a campaign sequentially."""
    try:
        campaign = CampaignDB.get(campaign_id)
        if not campaign:
            raise HTTPException(status_code=404, detail="Campaign not found")

        if campaign['status'] == 'running':
            raise HTTPException(status_code=400, detail="Campaign already running")

        # Start campaign execution in background
        background_tasks.add_task(run_campaign, campaign_id)

        return {
            "campaign_id": campaign_id,
            "status": "running",
            "message": "Campaign execution started"
        }
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/{campaign_id}")
async def get_campaign(campaign_id: str):
    """Get campaign details with progress."""
    try:
        campaign = CampaignDB.get(campaign_id)
        if not campaign:
            raise HTTPException(status_code=404, detail="Campaign not found")

        # Get all experiments in campaign
        experiments = CampaignDB.get_experiments(campaign_id)

        # Calculate progress
        total = len(experiments)
        completed = sum(1 for e in experiments if e['status'] == 'completed')
        failed = sum(1 for e in experiments if e['status'] == 'failed')
        running = sum(1 for e in experiments if e['status'] == 'running')

        return {
            "campaign": campaign,
            "experiments": experiments,
            "progress": {
                "total": total,
                "completed": completed,
                "failed": failed,
                "running": running,
                "percentage": (completed + failed) / total * 100 if total > 0 else 0
            }
        }
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/list")
async def list_campaigns(
    limit: int = 100,
    offset: int = 0,
    current_user: Optional[User] = Depends(get_optional_user)
):
    """
    List campaigns.
    Returns only user's own campaigns if authenticated.
    """
    try:
        campaigns = CampaignDB.list(limit, offset)

        # Filter by user if authenticated (only show user's own campaigns)
        if current_user:
            campaigns = [
                camp for camp in campaigns
                if camp.get('user_id') == current_user.id
            ]

        return {"campaigns": campaigns}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/{campaign_id}/results")
async def get_campaign_results(campaign_id: str):
    """Get aggregated results and analysis for a completed campaign."""
    try:
        campaign = CampaignDB.get(campaign_id)
        if not campaign:
            raise HTTPException(status_code=404, detail="Campaign not found")

        # Get all experiments in campaign
        experiments = CampaignDB.get_experiments(campaign_id)

        # Filter completed experiments
        completed_experiments = [e for e in experiments if e['status'] == 'completed']

        if not completed_experiments:
            return {
                "campaign_id": campaign_id,
                "message": "No completed experiments in campaign yet",
                "campaign": campaign,
                "experiments_count": len(experiments),
                "completed_count": 0
            }

        # Collect all energies and results
        energies = []
        hf_energies = []
        correlation_energies = []
        iterations_list = []
        convergence_data = []

        for exp in completed_experiments:
            exp_full = ExperimentDB.get(exp['id'])
            if exp_full and exp_full.get('results'):
                results = exp_full['results']

                energy = results.get('energy')
                if energy is not None:
                    energies.append({
                        'experiment_id': exp['id'],
                        'energy': energy,
                        'molecule': exp_full.get('molecule', {}).get('smiles', 'Unknown')
                    })

                if results.get('hf_energy') is not None:
                    hf_energies.append(results['hf_energy'])

                if results.get('correlation_energy') is not None:
                    correlation_energies.append(results['correlation_energy'])

                if results.get('iterations') is not None:
                    iterations_list.append(results['iterations'])

                # Collect convergence data
                if results.get('convergence_history'):
                    convergence_data.append({
                        'experiment_id': exp['id'],
                        'convergence': results['convergence_history']
                    })

        # Calculate aggregate statistics
        aggregate_stats = {
            "total_experiments": len(experiments),
            "completed_experiments": len(completed_experiments),
            "failed_experiments": len([e for e in experiments if e['status'] == 'failed']),
            "cancelled_experiments": len([e for e in experiments if e['status'] == 'cancelled']),
        }

        if energies:
            energy_values = [e['energy'] for e in energies]
            aggregate_stats["lowest_energy"] = min(energy_values)
            aggregate_stats["highest_energy"] = max(energy_values)
            aggregate_stats["average_energy"] = sum(energy_values) / len(energy_values)
            aggregate_stats["energy_range"] = max(energy_values) - min(energy_values)

            # Find best experiment (lowest energy)
            best_idx = energy_values.index(min(energy_values))
            aggregate_stats["best_experiment"] = energies[best_idx]

        if hf_energies:
            aggregate_stats["average_hf_energy"] = sum(hf_energies) / len(hf_energies)

        if correlation_energies:
            aggregate_stats["average_correlation_energy"] = sum(correlation_energies) / len(correlation_energies)

        if iterations_list:
            aggregate_stats["average_iterations"] = sum(iterations_list) / len(iterations_list)
            aggregate_stats["min_iterations"] = min(iterations_list)
            aggregate_stats["max_iterations"] = max(iterations_list)

        return {
            "campaign_id": campaign_id,
            "campaign": campaign,
            "aggregate_stats": aggregate_stats,
            "energies": energies,
            "convergence_data": convergence_data[:5]  # Limit to first 5 for visualization
        }

    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/{campaign_id}/cancel")
async def cancel_campaign(campaign_id: str):
    """Cancel a running campaign."""
    try:
        campaign = CampaignDB.get(campaign_id)
        if not campaign:
            raise HTTPException(status_code=404, detail="Campaign not found")

        if campaign['status'] != 'running':
            raise HTTPException(status_code=400, detail="Campaign is not running")

        # Update campaign status
        CampaignDB.update_status(campaign_id, 'cancelled')

        # Cancel any running experiments
        experiments = CampaignDB.get_experiments(campaign_id)
        for exp in experiments:
            if exp['status'] == 'running':
                ExperimentDB.update_status(exp['id'], 'cancelled')

        return {
            "campaign_id": campaign_id,
            "status": "cancelled",
            "message": "Campaign cancelled"
        }
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


def run_campaign(campaign_id: str):
    """Execute campaign experiments sequentially (runs in background)."""
    # Get or create a lock for this specific campaign
    with campaign_locks_lock:
        if campaign_id not in campaign_locks:
            campaign_locks[campaign_id] = threading.Lock()
        current_campaign_lock = campaign_locks[campaign_id]

    # Acquire semaphore slot (allows MAX_CONCURRENT_CAMPAIGNS to run in parallel)
    with campaign_semaphore:
        # Lock this specific campaign to prevent duplicate execution
        with current_campaign_lock:
            try:
                print(f"🚀 Starting campaign {campaign_id}")
                CampaignDB.update_status(campaign_id, 'running')

                # Get all experiments in sequence order
                experiments = CampaignDB.get_experiments(campaign_id)

                completed_count = 0
                failed_count = 0
                cancelled_count = 0

                for exp in experiments:
                    # Skip if already completed/failed/cancelled
                    if exp['status'] in ['completed', 'failed', 'cancelled']:
                        if exp['status'] == 'completed':
                            completed_count += 1
                        elif exp['status'] == 'failed':
                            failed_count += 1
                        else:
                            cancelled_count += 1
                        continue

                    exp_id = exp['id']
                    print(f"📊 Executing experiment {exp_id} ({completed_count + failed_count + cancelled_count + 1}/{len(experiments)})")

                    # Get or create job for this experiment
                    jobs = JobDB.list()
                    job_id = None
                    for job in jobs:
                        if job['experiment_id'] == exp_id:
                            job_id = job['id']
                            break

                    if not job_id:
                        # Create new job
                        job_id = str(uuid.uuid4())
                        JobDB.create({
                            'id': job_id,
                            'experiment_id': exp_id,
                            'status': 'queued',
                            'priority': 0
                        })

                    # Execute experiment synchronously
                    try:
                        execute_experiment(exp_id, job_id)

                        # Check final status (experiment_service catches cancellation internally)
                        final_exp = ExperimentDB.get(exp_id)
                        if final_exp['status'] == 'completed':
                            completed_count += 1
                        elif final_exp['status'] == 'cancelled':
                            cancelled_count += 1
                            print(f"⚠️  Experiment {exp_id} was cancelled")
                        else:
                            failed_count += 1
                    except Exception as e:
                        print(f"❌ Experiment {exp_id} failed: {e}")
                        failed_count += 1

                    # Update campaign progress
                    CampaignDB.update_progress(campaign_id, completed_count, failed_count)

                # Mark campaign as completed
                if cancelled_count > 0 and completed_count == 0:
                    final_status = 'cancelled'
                elif failed_count == 0 and cancelled_count == 0:
                    final_status = 'completed'
                else:
                    final_status = 'partially_completed'
                CampaignDB.update_status(campaign_id, final_status)

                print(f"✅ Campaign {campaign_id} completed: {completed_count} succeeded, {failed_count} failed, {cancelled_count} cancelled")

            except Exception as e:
                print(f"❌ Campaign {campaign_id} failed: {e}")
                CampaignDB.update_status(campaign_id, 'failed')

            finally:
                # Clean up campaign lock
                with campaign_locks_lock:
                    if campaign_id in campaign_locks:
                        del campaign_locks[campaign_id]
