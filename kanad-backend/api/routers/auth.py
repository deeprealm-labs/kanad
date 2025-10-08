"""
Authentication router - User registration, login, and token management.
"""

from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy.orm import Session
from datetime import timedelta
import uuid

from core.database import get_db
from core.models import UserRegister, UserLogin, Token, UserProfile
from db.models import User, UserSettings
from utils.auth import hash_password, verify_password, create_access_token
from api.config import settings

router = APIRouter()


@router.post("/register", response_model=Token, status_code=status.HTTP_201_CREATED)
async def register(user_data: UserRegister, db: Session = Depends(get_db)):
    """
    Register a new user.

    Creates user account with hashed password and default settings.
    """
    # Check if user already exists
    existing_user = db.query(User).filter(User.email == user_data.email).first()
    if existing_user:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Email already registered"
        )

    # Create user
    user = User(
        user_id=uuid.uuid4(),
        email=user_data.email,
        hashed_password=hash_password(user_data.password),
        name=user_data.name,
        institution=user_data.institution,
        research_field=user_data.field
    )

    db.add(user)
    db.flush()

    # Create default settings
    user_settings = UserSettings(user_id=user.user_id)
    db.add(user_settings)

    db.commit()

    # Generate access token
    access_token = create_access_token(
        data={"sub": str(user.user_id)},
        expires_delta=timedelta(minutes=settings.JWT_ACCESS_TOKEN_EXPIRE_MINUTES)
    )

    return {"access_token": access_token, "token_type": "bearer"}


@router.post("/login", response_model=Token)
async def login(credentials: UserLogin, db: Session = Depends(get_db)):
    """
    Login with email and password.

    Returns JWT access token.
    """
    # Find user
    user = db.query(User).filter(User.email == credentials.email).first()

    if not user or not verify_password(credentials.password, user.hashed_password):
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Incorrect email or password",
            headers={"WWW-Authenticate": "Bearer"},
        )

    # Generate access token
    access_token = create_access_token(
        data={"sub": str(user.user_id)},
        expires_delta=timedelta(minutes=settings.JWT_ACCESS_TOKEN_EXPIRE_MINUTES)
    )

    return {"access_token": access_token, "token_type": "bearer"}


@router.get("/profile", response_model=UserProfile)
async def get_profile(current_user: User = Depends(lambda: __import__('utils.auth', fromlist=['get_current_user']).get_current_user), db: Session = Depends(get_db)):
    """
    Get current user's profile.

    Requires authentication.
    """
    # Get user stats
    job_count = db.query(__import__('db.models', fromlist=['Job']).Job).filter_by(user_id=current_user.user_id).count()
    completed_jobs = db.query(__import__('db.models', fromlist=['Job']).Job).filter_by(
        user_id=current_user.user_id,
        status='completed'
    ).count()
    molecule_count = db.query(__import__('db.models', fromlist=['Molecule']).Molecule).filter_by(user_id=current_user.user_id).count()

    return {
        "user_id": str(current_user.user_id),
        "email": current_user.email,
        "name": current_user.name,
        "institution": current_user.institution,
        "field": current_user.research_field,
        "stats": {
            "total_jobs": job_count,
            "completed_jobs": completed_jobs,
            "total_molecules": molecule_count,
            "quantum_credits_used": 0  # Placeholder
        }
    }
