"""
Campaign/batch execution endpoints
"""

from fastapi import APIRouter, HTTPException, BackgroundTasks
from pydantic import BaseModel
from typing import List, Optional
import threading
import uuid

from api.core.database import CampaignDB, ExperimentDB, JobDB
from api.services.experiment_service import execute_experiment

router = APIRouter()

# Global lock for sequential campaign execution
campaign_lock = threading.Lock()


class CreateCampaignRequest(BaseModel):
    name: str
    description: Optional[str] = ""
    experiment_ids: List[str]


class ReorderQueueRequest(BaseModel):
    experiment_id: str
    new_order: int


@router.post("/create")
async def create_campaign(request: CreateCampaignRequest):
    """Create a campaign from queued experiments."""
    try:
        # Create campaign
        campaign_id = CampaignDB.create({
            'name': request.name,
            'description': request.description,
            'total_experiments': len(request.experiment_ids)
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
async def list_campaigns(limit: int = 100, offset: int = 0):
    """List all campaigns."""
    try:
        campaigns = CampaignDB.list(limit, offset)
        return {"campaigns": campaigns}
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
    with campaign_lock:  # Ensure only one campaign runs at a time
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
