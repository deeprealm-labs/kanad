"""
Authentication Routes for Kanad Platform
Handles registration, login, email verification, OAuth, and token management
"""

from fastapi import APIRouter, Depends, HTTPException, status, Request
from sqlalchemy.orm import Session
from pydantic import BaseModel, EmailStr, Field
from typing import Optional
from datetime import datetime, timedelta

from api.core.database_postgres import (
    get_db,
    User,
    UserRole,
    AccessKey,
    Session as DBSession,
)
from api.auth.password import validate_and_hash_password, verify_password
from api.auth.jwt_handler import (
    create_token_pair,
    verify_refresh_token,
    get_token_expiration,
)
from api.auth.email_otp import (
    create_otp,
    verify_otp,
    send_otp_email,
    resend_otp,
    send_welcome_email,
)
from api.auth.google_oauth import (
    GoogleOAuthHandler,
    create_or_update_google_user,
    is_google_oauth_configured,
)
from api.dependencies.auth import get_current_user, get_current_active_user
from api.middleware.rate_limit import login_limiter, register_limiter, rate_limit_dependency

router = APIRouter(prefix="/auth", tags=["Authentication"])


# Pydantic models
class RegisterRequest(BaseModel):
    email: EmailStr
    password: str = Field(..., min_length=8, max_length=128)
    full_name: str = Field(..., min_length=1, max_length=255)
    access_key: str = Field(..., min_length=1, description="Early access key required")


class RegisterResponse(BaseModel):
    message: str
    email: str
    requires_verification: bool


class LoginRequest(BaseModel):
    email: EmailStr
    password: str


class LoginResponse(BaseModel):
    access_token: str
    refresh_token: str
    token_type: str = "bearer"
    expires_in: int
    user: dict


class VerifyEmailRequest(BaseModel):
    email: EmailStr
    otp: str = Field(..., min_length=6, max_length=6)


class VerifyEmailResponse(BaseModel):
    message: str
    access_token: str
    refresh_token: str
    token_type: str = "bearer"


class ResendOTPRequest(BaseModel):
    email: EmailStr


class RefreshTokenRequest(BaseModel):
    refresh_token: str


class GoogleAuthRequest(BaseModel):
    id_token: Optional[str] = None
    code: Optional[str] = None
    access_key: Optional[str] = None  # Required for new users


class ChangePasswordRequest(BaseModel):
    old_password: str
    new_password: str = Field(..., min_length=8, max_length=128)


# Rate limit dependencies (must be defined before routes that use them)
async def _login_rate_limit(req: Request):
    await rate_limit_dependency(
        req, login_limiter, error_message="Too many login attempts. Please try again in 15 minutes."
    )

async def _register_rate_limit(req: Request):
    await rate_limit_dependency(
        req, register_limiter, error_message="Too many registration attempts. Please try again in 1 hour."
    )


# Routes
@router.post("/register", response_model=RegisterResponse, status_code=status.HTTP_201_CREATED)
async def register(
    request: RegisterRequest,
    req: Request,
    db: Session = Depends(get_db),
    _rate_limit: None = Depends(_register_rate_limit)
):
    """
    Register a new user with email and password

    Requires early access key for registration.
    After registration, user must verify email before full access.
    """
    # Check if email already exists
    existing_user = db.query(User).filter(User.email == request.email).first()
    if existing_user:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Email already registered",
        )

    # Validate and find access key
    access_key = (
        db.query(AccessKey)
        .filter(AccessKey.key == request.access_key, AccessKey.is_active == True)
        .first()
    )

    if not access_key:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Invalid access key",
        )

    # Check if access key is expired
    if access_key.expires_at and datetime.utcnow() > access_key.expires_at:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Access key has expired",
        )

    # Check if access key has remaining uses
    if access_key.used_count >= access_key.max_uses:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Access key has reached maximum uses",
        )

    # Validate and hash password
    password_hash, error = validate_and_hash_password(request.password)
    if error:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=error,
        )

    # Create user
    user = User(
        email=request.email,
        password_hash=password_hash,
        full_name=request.full_name,
        role=UserRole.USER,
        is_verified=False,
        is_active=True,
        access_key_id=access_key.id,
    )

    db.add(user)

    # Increment access key usage
    access_key.used_count += 1

    db.commit()
    db.refresh(user)

    # Send verification email
    otp, expires_at = create_otp(db, request.email)
    send_otp_email(request.email, otp, request.full_name)

    return RegisterResponse(
        message="Registration successful. Please check your email for verification code.",
        email=request.email,
        requires_verification=True,
    )


@router.post("/verify-email", response_model=VerifyEmailResponse)
async def verify_email(request: VerifyEmailRequest, db: Session = Depends(get_db)):
    """
    Verify email address using OTP code

    After verification, user receives access and refresh tokens.
    """
    # Find user
    user = db.query(User).filter(User.email == request.email).first()
    if not user:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="User not found",
        )

    # Check if already verified
    if user.is_verified:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Email already verified",
        )

    # Verify OTP
    is_valid, error = verify_otp(db, request.email, request.otp)
    if not is_valid:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=error,
        )

    # Mark user as verified
    user.is_verified = True
    user.last_login = datetime.utcnow()
    db.commit()
    db.refresh(user)

    # Create tokens
    tokens = create_token_pair(user.id, user.email, user.role.value)

    # Create session
    session = DBSession(
        user_id=user.id,
        access_token=tokens["access_token"],
        refresh_token=tokens["refresh_token"],
        expires_at=datetime.utcnow() + timedelta(hours=1),
        refresh_expires_at=datetime.utcnow() + timedelta(days=30),
    )
    db.add(session)
    db.commit()

    # Send welcome email
    send_welcome_email(user.email, user.full_name)

    return VerifyEmailResponse(
        message="Email verified successfully",
        access_token=tokens["access_token"],
        refresh_token=tokens["refresh_token"],
        token_type="bearer",
    )


@router.post("/resend-otp")
async def resend_otp_route(request: ResendOTPRequest, db: Session = Depends(get_db)):
    """
    Resend OTP verification code to email
    """
    # Find user
    user = db.query(User).filter(User.email == request.email).first()
    if not user:
        # Don't reveal if email exists
        return {"message": "If the email is registered, a new code has been sent"}

    if user.is_verified:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Email already verified",
        )

    # Resend OTP
    success, error = resend_otp(db, request.email)
    if not success:
        raise HTTPException(
            status_code=status.HTTP_429_TOO_MANY_REQUESTS,
            detail=error,
        )

    return {"message": "Verification code sent"}


@router.post("/login", response_model=LoginResponse)
async def login(
    request: LoginRequest,
    req: Request,
    db: Session = Depends(get_db),
    _rate_limit: None = Depends(_login_rate_limit)
):
    """
    Login with email and password

    Returns access and refresh tokens on success.
    """
    # Find user
    user = db.query(User).filter(User.email == request.email).first()
    if not user:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Invalid email or password",
        )

    # Verify password
    if not verify_password(request.password, user.password_hash):
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Invalid email or password",
        )

    # Check if user is active
    if not user.is_active:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Account is disabled",
        )

    # Check if email is verified
    if not user.is_verified:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Email not verified. Please check your email for verification code.",
        )

    # Update last login
    user.last_login = datetime.utcnow()
    db.commit()

    # Create tokens
    tokens = create_token_pair(user.id, user.email, user.role.value)

    # Create session
    session = DBSession(
        user_id=user.id,
        access_token=tokens["access_token"],
        refresh_token=tokens["refresh_token"],
        expires_at=datetime.utcnow() + timedelta(hours=1),
        refresh_expires_at=datetime.utcnow() + timedelta(days=30),
        user_agent=req.headers.get("user-agent"),
        ip_address=req.client.host if req.client else None,
    )
    db.add(session)
    db.commit()

    return LoginResponse(
        access_token=tokens["access_token"],
        refresh_token=tokens["refresh_token"],
        token_type="bearer",
        expires_in=3600,  # 1 hour in seconds
        user={
            "id": user.id,
            "email": user.email,
            "full_name": user.full_name,
            "role": user.role.value,
            "is_verified": user.is_verified,
            "avatar_url": user.avatar_url,
        },
    )


@router.post("/refresh", response_model=LoginResponse)
async def refresh_token(request: RefreshTokenRequest, db: Session = Depends(get_db)):
    """
    Refresh access token using refresh token

    Returns new access and refresh tokens.
    """
    # Verify refresh token
    payload = verify_refresh_token(request.refresh_token)
    if not payload:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Invalid or expired refresh token",
        )

    # Check if session exists and is not revoked
    session = (
        db.query(DBSession)
        .filter(
            DBSession.refresh_token == request.refresh_token,
            DBSession.is_revoked == False,
        )
        .first()
    )

    if not session:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Session not found or revoked",
        )

    # Get user
    user = db.query(User).filter(User.id == session.user_id).first()
    if not user or not user.is_active:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="User not found or inactive",
        )

    # Create new tokens
    tokens = create_token_pair(user.id, user.email, user.role.value)

    # Update session
    session.access_token = tokens["access_token"]
    session.refresh_token = tokens["refresh_token"]
    session.expires_at = datetime.utcnow() + timedelta(hours=1)
    session.refresh_expires_at = datetime.utcnow() + timedelta(days=30)
    session.last_used = datetime.utcnow()
    db.commit()

    return LoginResponse(
        access_token=tokens["access_token"],
        refresh_token=tokens["refresh_token"],
        token_type="bearer",
        expires_in=3600,
        user={
            "id": user.id,
            "email": user.email,
            "full_name": user.full_name,
            "role": user.role.value,
            "is_verified": user.is_verified,
            "avatar_url": user.avatar_url,
        },
    )


@router.post("/logout")
async def logout(
    current_user: User = Depends(get_current_active_user),
    authorization: str = Depends(lambda: None),
    db: Session = Depends(get_db),
):
    """
    Logout current user (revoke session)
    """
    # Extract token from authorization header
    # Note: In production, you'd get this from the dependency
    # For now, we'll revoke all active sessions for the user

    # Revoke all user's sessions
    db.query(DBSession).filter(
        DBSession.user_id == current_user.id, DBSession.is_revoked == False
    ).update({"is_revoked": True})

    db.commit()

    return {"message": "Logged out successfully"}


@router.get("/me")
async def get_current_user_info(current_user: User = Depends(get_current_active_user)):
    """
    Get current authenticated user information
    """
    return {
        "id": current_user.id,
        "email": current_user.email,
        "full_name": current_user.full_name,
        "role": current_user.role.value,
        "is_verified": current_user.is_verified,
        "is_active": current_user.is_active,
        "avatar_url": current_user.avatar_url,
        "created_at": current_user.created_at.isoformat(),
        "last_login": current_user.last_login.isoformat() if current_user.last_login else None,
        "has_google_auth": bool(current_user.google_id),
    }


@router.post("/change-password")
async def change_password(
    request: ChangePasswordRequest,
    current_user: User = Depends(get_current_active_user),
    db: Session = Depends(get_db),
):
    """
    Change user password
    """
    # Verify old password
    if not verify_password(request.old_password, current_user.password_hash):
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Incorrect current password",
        )

    # Validate and hash new password
    password_hash, error = validate_and_hash_password(request.new_password)
    if error:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=error,
        )

    # Update password
    current_user.password_hash = password_hash
    db.commit()

    # Revoke all sessions (force re-login)
    db.query(DBSession).filter(
        DBSession.user_id == current_user.id, DBSession.is_revoked == False
    ).update({"is_revoked": True})
    db.commit()

    return {"message": "Password changed successfully. Please login again."}


# Google OAuth routes
@router.get("/google/url")
async def get_google_auth_url():
    """
    Get Google OAuth authorization URL
    """
    if not is_google_oauth_configured():
        raise HTTPException(
            status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
            detail="Google OAuth not configured",
        )

    handler = GoogleOAuthHandler()
    auth_url = handler.get_authorization_url()

    return {"authorization_url": auth_url}


@router.post("/google", response_model=LoginResponse)
async def google_auth(
    request: GoogleAuthRequest,
    req: Request,
    db: Session = Depends(get_db),
):
    """
    Authenticate with Google OAuth

    Accepts either:
    - id_token: ID token from client-side Google Sign-In
    - code: Authorization code from server-side OAuth flow
    """
    if not is_google_oauth_configured():
        raise HTTPException(
            status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
            detail="Google OAuth not configured",
        )

    handler = GoogleOAuthHandler()
    google_user_info = None

    # Try ID token first
    if request.id_token:
        google_user_info = handler.verify_id_token(request.id_token)

    # Try authorization code
    elif request.code:
        google_user_info = handler.authenticate(request.code)

    if not google_user_info:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Google authentication failed",
        )

    # Create or update user (pass access_key for validation)
    user = create_or_update_google_user(db, google_user_info, request.access_key)
    if not user:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Valid access key required for new Google account registration",
        )

    # Create tokens
    tokens = create_token_pair(user.id, user.email, user.role.value)

    # Create session
    session = DBSession(
        user_id=user.id,
        access_token=tokens["access_token"],
        refresh_token=tokens["refresh_token"],
        expires_at=datetime.utcnow() + timedelta(hours=1),
        refresh_expires_at=datetime.utcnow() + timedelta(days=30),
        user_agent=req.headers.get("user-agent"),
        ip_address=req.client.host if req.client else None,
    )
    db.add(session)
    db.commit()

    return LoginResponse(
        access_token=tokens["access_token"],
        refresh_token=tokens["refresh_token"],
        token_type="bearer",
        expires_in=3600,
        user={
            "id": user.id,
            "email": user.email,
            "full_name": user.full_name,
            "role": user.role.value,
            "is_verified": user.is_verified,
            "avatar_url": user.avatar_url,
        },
    )


@router.get("/status")
async def get_auth_status():
    """
    Get authentication system status
    """
    return {
        "google_oauth_enabled": is_google_oauth_configured(),
        "email_verification_enabled": True,
        "access_key_required": True,
    }
