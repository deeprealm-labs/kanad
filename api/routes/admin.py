"""
Admin Routes for Kanad Platform
Priority: 5 (Keys) → 4 (Users) → 3 (Live Monitoring) → 2 (Usage) → 1 (Stats)
"""

from fastapi import APIRouter, Depends, HTTPException, status, Query
from sqlalchemy.orm import Session
from sqlalchemy import func, and_
from pydantic import BaseModel, Field, EmailStr
from typing import List, Optional
from datetime import datetime, timedelta
import secrets

from api.core.database_postgres import (
    get_db,
    User,
    UserRole,
    AccessKey,
    Experiment,
    Campaign,
    Job,
    ExperimentStatus,
)
from api.dependencies.auth import require_admin
from api.auth.password import generate_temporary_password, hash_password

router = APIRouter(prefix="/admin", tags=["Admin"], dependencies=[Depends(require_admin)])


# ============================================================================
# PRIORITY 5: ACCESS KEY MANAGEMENT
# ============================================================================

class CreateAccessKeyRequest(BaseModel):
    description: Optional[str] = None
    max_uses: int = Field(default=1, ge=1, le=1000)
    expires_in_days: Optional[int] = Field(default=None, ge=1, le=365)


class AccessKeyResponse(BaseModel):
    id: int
    key: str
    description: Optional[str]
    max_uses: int
    used_count: int
    is_active: bool
    expires_at: Optional[datetime]
    created_by: Optional[int]
    created_at: datetime


@router.post("/keys", response_model=AccessKeyResponse, status_code=status.HTTP_201_CREATED)
async def create_access_key(
    request: CreateAccessKeyRequest,
    current_admin: User = Depends(require_admin),
    db: Session = Depends(get_db),
):
    """
    [PRIORITY 5] Generate a new early access key

    Admin can specify:
    - max_uses: How many users can register with this key
    - expires_in_days: Optional expiration
    - description: Optional note about this key
    """
    # Generate unique key
    key = f"KANAD-{secrets.token_urlsafe(16).upper()}"

    # Calculate expiration
    expires_at = None
    if request.expires_in_days:
        expires_at = datetime.utcnow() + timedelta(days=request.expires_in_days)

    # Create access key
    access_key = AccessKey(
        key=key,
        description=request.description,
        max_uses=request.max_uses,
        used_count=0,
        is_active=True,
        expires_at=expires_at,
        created_by=current_admin.id,
    )

    db.add(access_key)
    db.commit()
    db.refresh(access_key)

    return access_key


@router.get("/keys", response_model=List[AccessKeyResponse])
async def list_access_keys(
    active_only: bool = Query(False, description="Show only active keys"),
    db: Session = Depends(get_db),
):
    """
    [PRIORITY 5] List all access keys
    """
    query = db.query(AccessKey)

    if active_only:
        query = query.filter(AccessKey.is_active == True)

    keys = query.order_by(AccessKey.created_at.desc()).all()
    return keys


@router.get("/keys/{key_id}", response_model=AccessKeyResponse)
async def get_access_key(key_id: int, db: Session = Depends(get_db)):
    """
    [PRIORITY 5] Get access key details
    """
    key = db.query(AccessKey).filter(AccessKey.id == key_id).first()
    if not key:
        raise HTTPException(status_code=404, detail="Access key not found")
    return key


@router.put("/keys/{key_id}/deactivate")
async def deactivate_access_key(key_id: int, db: Session = Depends(get_db)):
    """
    [PRIORITY 5] Deactivate an access key
    """
    key = db.query(AccessKey).filter(AccessKey.id == key_id).first()
    if not key:
        raise HTTPException(status_code=404, detail="Access key not found")

    key.is_active = False
    db.commit()

    return {"message": "Access key deactivated", "key": key.key}


@router.delete("/keys/{key_id}")
async def delete_access_key(key_id: int, db: Session = Depends(get_db)):
    """
    [PRIORITY 5] Delete an access key

    Warning: This will not affect users who already registered with this key.
    """
    key = db.query(AccessKey).filter(AccessKey.id == key_id).first()
    if not key:
        raise HTTPException(status_code=404, detail="Access key not found")

    # Check if key has been used
    if key.used_count > 0:
        raise HTTPException(
            status_code=400,
            detail=f"Cannot delete key that has been used {key.used_count} times. Deactivate instead.",
        )

    db.delete(key)
    db.commit()

    return {"message": "Access key deleted"}


# ============================================================================
# PRIORITY 4: USER MANAGEMENT
# ============================================================================

class UserResponse(BaseModel):
    id: int
    email: str
    full_name: str
    role: str
    is_verified: bool
    is_active: bool
    has_google_auth: bool
    access_key_id: Optional[int]
    created_at: datetime
    last_login: Optional[datetime]
    experiments_count: int
    campaigns_count: int


class UpdateUserRequest(BaseModel):
    role: Optional[UserRole] = None
    is_active: Optional[bool] = None
    is_verified: Optional[bool] = None


class CreateUserRequest(BaseModel):
    email: EmailStr
    full_name: str
    role: UserRole = UserRole.USER
    send_credentials: bool = Field(default=True, description="Send password via email")


@router.get("/users", response_model=List[UserResponse])
async def list_users(
    role: Optional[UserRole] = None,
    verified_only: bool = False,
    active_only: bool = True,
    skip: int = 0,
    limit: int = 100,
    db: Session = Depends(get_db),
):
    """
    [PRIORITY 4] List all users with filters
    """
    query = db.query(
        User,
        func.count(Experiment.id).label("experiments_count"),
        func.count(Campaign.id).label("campaigns_count"),
    ).outerjoin(Experiment).outerjoin(Campaign).group_by(User.id)

    if role:
        query = query.filter(User.role == role)
    if verified_only:
        query = query.filter(User.is_verified == True)
    if active_only:
        query = query.filter(User.is_active == True)

    users = query.offset(skip).limit(limit).all()

    return [
        {
            "id": user.id,
            "email": user.email,
            "full_name": user.full_name,
            "role": user.role.value,
            "is_verified": user.is_verified,
            "is_active": user.is_active,
            "has_google_auth": bool(user.google_id),
            "access_key_id": user.access_key_id,
            "created_at": user.created_at,
            "last_login": user.last_login,
            "experiments_count": exp_count,
            "campaigns_count": camp_count,
        }
        for user, exp_count, camp_count in users
    ]


@router.get("/users/{user_id}", response_model=UserResponse)
async def get_user(user_id: int, db: Session = Depends(get_db)):
    """
    [PRIORITY 4] Get user details
    """
    user = db.query(User).filter(User.id == user_id).first()
    if not user:
        raise HTTPException(status_code=404, detail="User not found")

    experiments_count = db.query(Experiment).filter(Experiment.user_id == user_id).count()
    campaigns_count = db.query(Campaign).filter(Campaign.user_id == user_id).count()

    return {
        "id": user.id,
        "email": user.email,
        "full_name": user.full_name,
        "role": user.role.value,
        "is_verified": user.is_verified,
        "is_active": user.is_active,
        "has_google_auth": bool(user.google_id),
        "access_key_id": user.access_key_id,
        "created_at": user.created_at,
        "last_login": user.last_login,
        "experiments_count": experiments_count,
        "campaigns_count": campaigns_count,
    }


@router.put("/users/{user_id}")
async def update_user(
    user_id: int,
    request: UpdateUserRequest,
    current_admin: User = Depends(require_admin),
    db: Session = Depends(get_db),
):
    """
    [PRIORITY 4] Update user details

    Admin can change:
    - role: Change user role (admin/user/viewer)
    - is_active: Enable/disable account
    - is_verified: Force verify email
    """
    user = db.query(User).filter(User.id == user_id).first()
    if not user:
        raise HTTPException(status_code=404, detail="User not found")

    # Prevent self-demotion
    if user_id == current_admin.id and request.role and request.role != UserRole.ADMIN:
        raise HTTPException(status_code=400, detail="Cannot change your own admin role")

    # Update fields
    if request.role is not None:
        user.role = request.role
    if request.is_active is not None:
        user.is_active = request.is_active
    if request.is_verified is not None:
        user.is_verified = request.is_verified

    db.commit()
    db.refresh(user)

    return {"message": "User updated", "user_id": user_id}


@router.post("/users", response_model=UserResponse, status_code=status.HTTP_201_CREATED)
async def create_user(
    request: CreateUserRequest,
    db: Session = Depends(get_db),
):
    """
    [PRIORITY 4] Create user manually (admin only)

    Generates temporary password and optionally emails it.
    """
    # Check if email exists
    existing = db.query(User).filter(User.email == request.email).first()
    if existing:
        raise HTTPException(status_code=400, detail="Email already exists")

    # Generate temporary password
    temp_password = generate_temporary_password()
    password_hash = hash_password(temp_password)

    # Create user
    user = User(
        email=request.email,
        full_name=request.full_name,
        password_hash=password_hash,
        role=request.role,
        is_verified=True,  # Admin-created users are pre-verified
        is_active=True,
    )

    db.add(user)
    db.commit()
    db.refresh(user)

    # TODO: Send credentials email if requested
    if request.send_credentials:
        print(f"TODO: Send credentials to {user.email}: Password: {temp_password}")

    return {
        "id": user.id,
        "email": user.email,
        "full_name": user.full_name,
        "role": user.role.value,
        "is_verified": user.is_verified,
        "is_active": user.is_active,
        "has_google_auth": False,
        "access_key_id": None,
        "created_at": user.created_at,
        "last_login": None,
        "experiments_count": 0,
        "campaigns_count": 0,
    }


@router.delete("/users/{user_id}")
async def delete_user(
    user_id: int,
    current_admin: User = Depends(require_admin),
    db: Session = Depends(get_db),
):
    """
    [PRIORITY 4] Delete user

    Warning: This will also delete all user's experiments and campaigns.
    """
    if user_id == current_admin.id:
        raise HTTPException(status_code=400, detail="Cannot delete yourself")

    user = db.query(User).filter(User.id == user_id).first()
    if not user:
        raise HTTPException(status_code=404, detail="User not found")

    # Count related data
    exp_count = db.query(Experiment).filter(Experiment.user_id == user_id).count()
    camp_count = db.query(Campaign).filter(Campaign.user_id == user_id).count()

    # Delete user (cascades will handle related data)
    db.delete(user)
    db.commit()

    return {
        "message": "User deleted",
        "deleted_experiments": exp_count,
        "deleted_campaigns": camp_count,
    }


# ============================================================================
# PRIORITY 3: LIVE EXPERIMENT MONITORING
# ============================================================================

class LiveExperimentResponse(BaseModel):
    id: int
    name: str
    user_email: str
    status: str
    progress: float
    backend: Optional[str]
    method: Optional[str]
    started_at: Optional[datetime]
    running_time_seconds: Optional[int]


@router.get("/experiments/live", response_model=List[LiveExperimentResponse])
async def get_live_experiments(db: Session = Depends(get_db)):
    """
    [PRIORITY 3] Get all currently running experiments

    Real-time monitoring of active quantum computations.
    """
    running_experiments = (
        db.query(Experiment, User)
        .join(User, Experiment.user_id == User.id, isouter=True)
        .filter(Experiment.status == ExperimentStatus.RUNNING)
        .order_by(Experiment.started_at.desc())
        .all()
    )

    results = []
    for exp, user in running_experiments:
        running_time = None
        if exp.started_at:
            running_time = int((datetime.utcnow() - exp.started_at).total_seconds())

        results.append({
            "id": exp.id,
            "name": exp.name,
            "user_email": user.email if user else "System",
            "status": exp.status.value,
            "progress": exp.progress,
            "backend": exp.backend,
            "method": exp.method,
            "started_at": exp.started_at,
            "running_time_seconds": running_time,
        })

    return results


@router.get("/experiments/recent")
async def get_recent_experiments(
    limit: int = Query(50, ge=1, le=500),
    db: Session = Depends(get_db),
):
    """
    [PRIORITY 3] Get recently completed/failed experiments
    """
    experiments = (
        db.query(Experiment, User)
        .join(User, Experiment.user_id == User.id, isouter=True)
        .filter(
            Experiment.status.in_([
                ExperimentStatus.COMPLETED,
                ExperimentStatus.FAILED,
                ExperimentStatus.CANCELLED,
            ])
        )
        .order_by(Experiment.completed_at.desc())
        .limit(limit)
        .all()
    )

    return [
        {
            "id": exp.id,
            "name": exp.name,
            "user_email": user.email if user else "System",
            "status": exp.status.value,
            "progress": exp.progress,
            "completed_at": exp.completed_at,
            "error_message": exp.error_message,
        }
        for exp, user in experiments
    ]


# ============================================================================
# PRIORITY 2: BACKEND USAGE TRACKING
# ============================================================================

@router.get("/usage/by-backend")
async def get_usage_by_backend(
    days: int = Query(30, ge=1, le=365),
    db: Session = Depends(get_db),
):
    """
    [PRIORITY 2] Get compute usage statistics by quantum backend

    Shows which backends are being used most.
    """
    since = datetime.utcnow() - timedelta(days=days)

    usage = (
        db.query(
            Experiment.backend,
            func.count(Experiment.id).label("count"),
            func.avg(Experiment.progress).label("avg_progress"),
        )
        .filter(Experiment.created_at >= since)
        .group_by(Experiment.backend)
        .order_by(func.count(Experiment.id).desc())
        .all()
    )

    return [
        {
            "backend": backend or "Unknown",
            "experiment_count": count,
            "average_progress": round(avg_progress, 2) if avg_progress else 0,
        }
        for backend, count, avg_progress in usage
    ]


@router.get("/usage/by-user")
async def get_usage_by_user(
    days: int = Query(30, ge=1, le=365),
    limit: int = Query(10, ge=1, le=100),
    db: Session = Depends(get_db),
):
    """
    [PRIORITY 2] Get compute usage by user (top users)
    """
    since = datetime.utcnow() - timedelta(days=days)

    usage = (
        db.query(
            User.email,
            User.full_name,
            func.count(Experiment.id).label("experiment_count"),
            func.count(
                and_(Experiment.status == ExperimentStatus.COMPLETED, True)
            ).label("completed_count"),
        )
        .join(Experiment, User.id == Experiment.user_id)
        .filter(Experiment.created_at >= since)
        .group_by(User.id, User.email, User.full_name)
        .order_by(func.count(Experiment.id).desc())
        .limit(limit)
        .all()
    )

    return [
        {
            "user_email": email,
            "user_name": name,
            "total_experiments": exp_count,
            "completed_experiments": comp_count,
            "success_rate": round((comp_count / exp_count * 100), 2) if exp_count > 0 else 0,
        }
        for email, name, exp_count, comp_count in usage
    ]


@router.get("/usage/timeline")
async def get_usage_timeline(
    days: int = Query(30, ge=1, le=365),
    db: Session = Depends(get_db),
):
    """
    [PRIORITY 2] Get daily experiment counts over time
    """
    since = datetime.utcnow() - timedelta(days=days)

    # This is a simplified version - for production, you'd want proper date binning
    daily_usage = (
        db.query(
            func.date(Experiment.created_at).label("date"),
            func.count(Experiment.id).label("count"),
        )
        .filter(Experiment.created_at >= since)
        .group_by(func.date(Experiment.created_at))
        .order_by(func.date(Experiment.created_at).desc())
        .all()
    )

    return [
        {"date": str(date), "experiment_count": count}
        for date, count in daily_usage
    ]


# ============================================================================
# PRIORITY 1: SYSTEM STATISTICS
# ============================================================================

@router.get("/stats/overview")
async def get_system_stats(db: Session = Depends(get_db)):
    """
    [PRIORITY 1] Get overall system statistics

    High-level dashboard metrics.
    """
    total_users = db.query(User).count()
    active_users = db.query(User).filter(User.is_active == True).count()
    verified_users = db.query(User).filter(User.is_verified == True).count()

    total_experiments = db.query(Experiment).count()
    running_experiments = db.query(Experiment).filter(
        Experiment.status == ExperimentStatus.RUNNING
    ).count()
    completed_experiments = db.query(Experiment).filter(
        Experiment.status == ExperimentStatus.COMPLETED
    ).count()

    total_campaigns = db.query(Campaign).count()

    total_access_keys = db.query(AccessKey).count()
    active_access_keys = db.query(AccessKey).filter(AccessKey.is_active == True).count()

    # Recent activity (last 24 hours)
    last_24h = datetime.utcnow() - timedelta(hours=24)
    experiments_today = db.query(Experiment).filter(
        Experiment.created_at >= last_24h
    ).count()
    new_users_today = db.query(User).filter(User.created_at >= last_24h).count()

    return {
        "users": {
            "total": total_users,
            "active": active_users,
            "verified": verified_users,
            "new_today": new_users_today,
        },
        "experiments": {
            "total": total_experiments,
            "running": running_experiments,
            "completed": completed_experiments,
            "started_today": experiments_today,
            "success_rate": round(
                (completed_experiments / total_experiments * 100) if total_experiments > 0 else 0,
                2,
            ),
        },
        "campaigns": {
            "total": total_campaigns,
        },
        "access_keys": {
            "total": total_access_keys,
            "active": active_access_keys,
        },
    }


@router.get("/stats/growth")
async def get_growth_stats(
    days: int = Query(30, ge=1, le=365),
    db: Session = Depends(get_db),
):
    """
    [PRIORITY 1] Get user and experiment growth over time
    """
    since = datetime.utcnow() - timedelta(days=days)

    # User growth
    user_growth = (
        db.query(
            func.date(User.created_at).label("date"),
            func.count(User.id).label("count"),
        )
        .filter(User.created_at >= since)
        .group_by(func.date(User.created_at))
        .order_by(func.date(User.created_at))
        .all()
    )

    # Experiment growth
    exp_growth = (
        db.query(
            func.date(Experiment.created_at).label("date"),
            func.count(Experiment.id).label("count"),
        )
        .filter(Experiment.created_at >= since)
        .group_by(func.date(Experiment.created_at))
        .order_by(func.date(Experiment.created_at))
        .all()
    )

    return {
        "user_growth": [{"date": str(date), "new_users": count} for date, count in user_growth],
        "experiment_growth": [
            {"date": str(date), "new_experiments": count} for date, count in exp_growth
        ],
    }
