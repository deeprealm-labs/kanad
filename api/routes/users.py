"""
User Profile and Account Management Routes
Handles user profile, settings, statistics, and session management
"""

from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy.orm import Session
from pydantic import BaseModel, EmailStr, Field
from typing import Optional, List
from datetime import datetime

from api.core.database_postgres import get_db, User, Session as DBSession
from api.dependencies.auth import get_current_user, get_current_verified_user
from api.auth.password import validate_and_hash_password, verify_password

router = APIRouter(prefix="/users", tags=["Users"])


# Pydantic Models
class UserProfileResponse(BaseModel):
    id: int
    email: str
    full_name: str
    role: str
    is_verified: bool
    is_active: bool
    google_id: Optional[str]
    avatar_url: Optional[str]
    created_at: datetime
    last_login: Optional[datetime]
    access_key_used: bool

    class Config:
        from_attributes = True


class UpdateProfileRequest(BaseModel):
    full_name: Optional[str] = Field(None, min_length=1, max_length=255)
    avatar_url: Optional[str] = None


class ChangePasswordRequest(BaseModel):
    current_password: str
    new_password: str = Field(..., min_length=8, max_length=128)


class UserStatisticsResponse(BaseModel):
    total_experiments: int
    completed_experiments: int
    failed_experiments: int
    total_campaigns: int
    account_age_days: int
    last_experiment_date: Optional[datetime]


class SessionInfo(BaseModel):
    id: int
    created_at: datetime
    expires_at: datetime
    last_used: Optional[datetime]
    is_current: bool


# Routes
@router.get("/me", response_model=UserProfileResponse)
async def get_current_user_profile(
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """Get current user's profile information"""
    return {
        "id": current_user.id,
        "email": current_user.email,
        "full_name": current_user.full_name,
        "role": current_user.role.value,
        "is_verified": current_user.is_verified,
        "is_active": current_user.is_active,
        "google_id": current_user.google_id,
        "avatar_url": current_user.avatar_url,
        "created_at": current_user.created_at,
        "last_login": current_user.last_login,
        "access_key_used": current_user.access_key_id is not None,
    }


@router.put("/me", response_model=UserProfileResponse)
async def update_user_profile(
    profile_data: UpdateProfileRequest,
    current_user: User = Depends(get_current_verified_user),
    db: Session = Depends(get_db)
):
    """Update current user's profile information"""

    # Update fields if provided
    if profile_data.full_name is not None:
        current_user.full_name = profile_data.full_name

    if profile_data.avatar_url is not None:
        current_user.avatar_url = profile_data.avatar_url

    db.commit()
    db.refresh(current_user)

    return {
        "id": current_user.id,
        "email": current_user.email,
        "full_name": current_user.full_name,
        "role": current_user.role.value,
        "is_verified": current_user.is_verified,
        "is_active": current_user.is_active,
        "google_id": current_user.google_id,
        "avatar_url": current_user.avatar_url,
        "created_at": current_user.created_at,
        "last_login": current_user.last_login,
        "access_key_used": current_user.access_key_id is not None,
    }


@router.post("/change-password")
async def change_password(
    password_data: ChangePasswordRequest,
    current_user: User = Depends(get_current_verified_user),
    db: Session = Depends(get_db)
):
    """Change user's password"""

    # Check if user uses OAuth (no password)
    if current_user.google_id and not current_user.password_hash:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Cannot change password for Google OAuth accounts",
        )

    # Verify current password
    if not verify_password(password_data.current_password, current_user.password_hash):
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Current password is incorrect",
        )

    # Validate and hash new password
    try:
        new_password_hash = validate_and_hash_password(password_data.new_password)
    except ValueError as e:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=str(e),
        )

    # Update password
    current_user.password_hash = new_password_hash
    db.commit()

    return {"message": "Password changed successfully"}


@router.get("/me/statistics", response_model=UserStatisticsResponse)
async def get_user_statistics(
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """Get current user's usage statistics"""
    from api.core.database import ExperimentDB, CampaignDB

    # Get experiments statistics
    all_experiments = ExperimentDB.list(limit=10000)
    user_experiments = [exp for exp in all_experiments if exp.get('user_id') == current_user.id]

    total_experiments = len(user_experiments)
    completed_experiments = len([exp for exp in user_experiments if exp.get('status') == 'completed'])
    failed_experiments = len([exp for exp in user_experiments if exp.get('status') == 'failed'])

    # Get campaigns statistics
    all_campaigns = CampaignDB.list(limit=10000)
    user_campaigns = [camp for camp in all_campaigns if camp.get('user_id') == current_user.id]
    total_campaigns = len(user_campaigns)

    # Calculate account age
    account_age = datetime.utcnow() - current_user.created_at
    account_age_days = account_age.days

    # Get last experiment date
    last_experiment_date = None
    if user_experiments:
        sorted_experiments = sorted(
            user_experiments,
            key=lambda x: x.get('created_at', ''),
            reverse=True
        )
        if sorted_experiments:
            last_exp_date_str = sorted_experiments[0].get('created_at')
            if last_exp_date_str:
                try:
                    last_experiment_date = datetime.fromisoformat(last_exp_date_str)
                except:
                    pass

    return {
        "total_experiments": total_experiments,
        "completed_experiments": completed_experiments,
        "failed_experiments": failed_experiments,
        "total_campaigns": total_campaigns,
        "account_age_days": account_age_days,
        "last_experiment_date": last_experiment_date,
    }


@router.get("/me/sessions", response_model=List[SessionInfo])
async def get_user_sessions(
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """Get all active sessions for current user"""

    # Get all sessions for user
    sessions = db.query(DBSession).filter(
        DBSession.user_id == current_user.id,
        DBSession.is_revoked == False
    ).order_by(DBSession.created_at.desc()).all()

    # Get current session token from context (would need to be passed)
    # For now, mark the most recent as current
    current_session_id = sessions[0].id if sessions else None

    return [
        {
            "id": session.id,
            "created_at": session.created_at,
            "expires_at": session.expires_at,
            "last_used": session.last_used,
            "is_current": session.id == current_session_id,
        }
        for session in sessions
    ]


@router.delete("/me/sessions/{session_id}")
async def revoke_session(
    session_id: int,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """Revoke a specific session"""

    # Get session
    session = db.query(DBSession).filter(
        DBSession.id == session_id,
        DBSession.user_id == current_user.id
    ).first()

    if not session:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Session not found",
        )

    # Revoke session
    session.is_revoked = True
    db.commit()

    return {"message": "Session revoked successfully"}


@router.delete("/me")
async def delete_account(
    current_user: User = Depends(get_current_verified_user),
    db: Session = Depends(get_db)
):
    """Deactivate user account (soft delete)"""

    # Soft delete - just deactivate the account
    current_user.is_active = False

    # Revoke all sessions
    db.query(DBSession).filter(
        DBSession.user_id == current_user.id,
        DBSession.is_revoked == False
    ).update({
        "is_revoked": True
    })

    db.commit()

    return {"message": "Account deactivated successfully"}
