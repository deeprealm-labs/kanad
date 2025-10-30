"""
FastAPI Authentication Dependencies
Provides dependency injection for protected routes and role-based access control
"""

from typing import Optional
from fastapi import Depends, HTTPException, status, Header
from fastapi.security import HTTPBearer, HTTPAuthorizationCredentials
from sqlalchemy.orm import Session

from api.core.database_postgres import User, UserRole, get_db
from api.auth.jwt_handler import validate_token

# Security scheme for Swagger UI
security = HTTPBearer(
    scheme_name="JWT Bearer Token",
    description="Enter your JWT token (without 'Bearer' prefix)",
)


async def get_token_from_header(
    authorization: Optional[str] = Header(None),
) -> Optional[str]:
    """
    Extract token from Authorization header

    Args:
        authorization: Authorization header value

    Returns:
        Token string or None
    """
    if not authorization:
        return None

    # Handle "Bearer <token>" format
    if authorization.startswith("Bearer "):
        return authorization[7:]

    # Return as-is if not Bearer format
    return authorization


async def get_current_user(
    credentials: HTTPAuthorizationCredentials = Depends(security),
    db: Session = Depends(get_db),
) -> User:
    """
    Get current authenticated user from JWT token

    Args:
        credentials: HTTP authorization credentials
        db: Database session

    Returns:
        User object

    Raises:
        HTTPException: If authentication fails
    """
    token = credentials.credentials

    # Validate token
    result = validate_token(token, "access")

    if not result.valid:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail=result.error or "Invalid authentication token",
            headers={"WWW-Authenticate": "Bearer"},
        )

    # Get user from database
    user_id = result.payload.get("user_id")
    user = db.query(User).filter(User.id == user_id).first()

    if not user:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="User not found",
            headers={"WWW-Authenticate": "Bearer"},
        )

    # Check if user is active
    if not user.is_active:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="User account is disabled",
        )

    # Update last login
    from datetime import datetime

    user.last_login = datetime.utcnow()
    db.commit()

    return user


async def get_current_active_user(
    current_user: User = Depends(get_current_user),
) -> User:
    """
    Get current active user (alias for get_current_user)

    Args:
        current_user: Current user from dependency

    Returns:
        User object
    """
    return current_user


async def get_current_verified_user(
    current_user: User = Depends(get_current_user),
) -> User:
    """
    Get current verified user (email verified)

    Args:
        current_user: Current user from dependency

    Returns:
        User object

    Raises:
        HTTPException: If email not verified
    """
    if not current_user.is_verified:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Email verification required. Please verify your email address.",
        )

    return current_user


# Role-based access control dependencies
async def require_admin(
    current_user: User = Depends(get_current_user),
) -> User:
    """
    Require admin role

    Args:
        current_user: Current user from dependency

    Returns:
        User object

    Raises:
        HTTPException: If user is not admin
    """
    if current_user.role != UserRole.ADMIN:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Admin access required",
        )

    return current_user


async def require_user_or_admin(
    current_user: User = Depends(get_current_verified_user),
) -> User:
    """
    Require user or admin role (verified)

    Args:
        current_user: Current user from dependency

    Returns:
        User object

    Raises:
        HTTPException: If user is viewer
    """
    if current_user.role == UserRole.VIEWER:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="User or Admin access required. Your account has viewer-only access.",
        )

    return current_user


async def require_viewer_or_above(
    current_user: User = Depends(get_current_verified_user),
) -> User:
    """
    Require any authenticated and verified user (viewer, user, or admin)

    Args:
        current_user: Current user from dependency

    Returns:
        User object
    """
    # All roles allowed, just need to be verified
    return current_user


# Optional authentication (doesn't fail if no token)
async def get_optional_user(
    authorization: Optional[str] = Header(None),
    db: Session = Depends(get_db),
) -> Optional[User]:
    """
    Get current user if authenticated, otherwise None

    Args:
        authorization: Authorization header
        db: Database session

    Returns:
        User object or None
    """
    if not authorization:
        return None

    try:
        # Extract token
        token = authorization[7:] if authorization.startswith("Bearer ") else authorization

        # Validate token
        result = validate_token(token, "access")
        if not result.valid:
            return None

        # Get user
        user_id = result.payload.get("user_id")
        user = db.query(User).filter(User.id == user_id, User.is_active == True).first()

        return user

    except Exception:
        return None


# Helper to check if user has access key
async def require_access_key(
    current_user: User = Depends(get_current_user),
) -> User:
    """
    Require user to have used an access key during registration

    Args:
        current_user: Current user from dependency

    Returns:
        User object

    Raises:
        HTTPException: If no access key associated
    """
    if not current_user.access_key_id:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Access key required. Please use an early access key to unlock features.",
        )

    return current_user


# Check specific permissions
def check_experiment_access(user: User, experiment_user_id: Optional[int]) -> bool:
    """
    Check if user has access to an experiment

    Args:
        user: Current user
        experiment_user_id: Owner user ID of experiment (None = system experiment)

    Returns:
        True if user has access
    """
    # Admin can access everything
    if user.role == UserRole.ADMIN:
        return True

    # System experiments (user_id = None) are accessible to all
    if experiment_user_id is None:
        return True

    # Users can access their own experiments
    if user.id == experiment_user_id:
        return True

    return False


def check_campaign_access(user: User, campaign_user_id: Optional[int]) -> bool:
    """
    Check if user has access to a campaign

    Args:
        user: Current user
        campaign_user_id: Owner user ID of campaign (None = system campaign)

    Returns:
        True if user has access
    """
    # Admin can access everything
    if user.role == UserRole.ADMIN:
        return True

    # System campaigns (user_id = None) are accessible to all
    if campaign_user_id is None:
        return True

    # Users can access their own campaigns
    if user.id == campaign_user_id:
        return True

    return False


# Dependency to check experiment ownership
async def verify_experiment_access(
    experiment_id: int,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db),
) -> bool:
    """
    Verify user has access to experiment

    Args:
        experiment_id: Experiment ID
        current_user: Current user
        db: Database session

    Returns:
        True if access granted

    Raises:
        HTTPException: If access denied
    """
    from api.core.database_postgres import Experiment

    experiment = db.query(Experiment).filter(Experiment.id == experiment_id).first()

    if not experiment:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Experiment not found",
        )

    if not check_experiment_access(current_user, experiment.user_id):
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="You don't have access to this experiment",
        )

    return True
