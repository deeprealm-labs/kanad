"""
JWT Token Handler for Kanad Platform
Handles generation, validation, and refresh of JSON Web Tokens
"""

import os
from datetime import datetime, timedelta
from typing import Optional, Dict, Any
import jwt
from jwt.exceptions import InvalidTokenError, ExpiredSignatureError
import secrets

# JWT Configuration
SECRET_KEY = os.getenv("JWT_SECRET_KEY", secrets.token_urlsafe(32))
ALGORITHM = "HS256"
ACCESS_TOKEN_EXPIRE_MINUTES = 60  # 1 hour
REFRESH_TOKEN_EXPIRE_DAYS = 30  # 30 days


class TokenPayload:
    """JWT Token Payload Structure"""

    def __init__(
        self,
        user_id: int,
        email: str,
        role: str,
        token_type: str = "access",
        exp: Optional[datetime] = None,
    ):
        self.user_id = user_id
        self.email = email
        self.role = role
        self.token_type = token_type
        self.exp = exp or (
            datetime.utcnow()
            + (
                timedelta(minutes=ACCESS_TOKEN_EXPIRE_MINUTES)
                if token_type == "access"
                else timedelta(days=REFRESH_TOKEN_EXPIRE_DAYS)
            )
        )
        self.iat = datetime.utcnow()

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JWT encoding"""
        return {
            "user_id": self.user_id,
            "email": self.email,
            "role": self.role,
            "token_type": self.token_type,
            "exp": self.exp,
            "iat": self.iat,
        }


def create_access_token(user_id: int, email: str, role: str) -> str:
    """
    Create a new access token (short-lived)

    Args:
        user_id: User's database ID
        email: User's email address
        role: User's role (admin, user, viewer)

    Returns:
        Encoded JWT access token
    """
    payload = TokenPayload(
        user_id=user_id, email=email, role=role, token_type="access"
    )

    token = jwt.encode(payload.to_dict(), SECRET_KEY, algorithm=ALGORITHM)
    return token


def create_refresh_token(user_id: int, email: str, role: str) -> str:
    """
    Create a new refresh token (long-lived)

    Args:
        user_id: User's database ID
        email: User's email address
        role: User's role (admin, user, viewer)

    Returns:
        Encoded JWT refresh token
    """
    payload = TokenPayload(
        user_id=user_id, email=email, role=role, token_type="refresh"
    )

    token = jwt.encode(payload.to_dict(), SECRET_KEY, algorithm=ALGORITHM)
    return token


def create_token_pair(user_id: int, email: str, role: str) -> Dict[str, str]:
    """
    Create both access and refresh tokens

    Args:
        user_id: User's database ID
        email: User's email address
        role: User's role (admin, user, viewer)

    Returns:
        Dictionary with access_token and refresh_token
    """
    return {
        "access_token": create_access_token(user_id, email, role),
        "refresh_token": create_refresh_token(user_id, email, role),
    }


def decode_token(token: str) -> Optional[Dict[str, Any]]:
    """
    Decode and validate a JWT token

    Args:
        token: JWT token string

    Returns:
        Decoded payload dictionary or None if invalid

    Raises:
        ExpiredSignatureError: If token has expired
        InvalidTokenError: If token is invalid
    """
    try:
        payload = jwt.decode(token, SECRET_KEY, algorithms=[ALGORITHM])
        return payload
    except ExpiredSignatureError:
        raise ExpiredSignatureError("Token has expired")
    except InvalidTokenError as e:
        raise InvalidTokenError(f"Invalid token: {str(e)}")


def verify_access_token(token: str) -> Optional[Dict[str, Any]]:
    """
    Verify an access token

    Args:
        token: JWT access token

    Returns:
        Decoded payload if valid, None otherwise
    """
    try:
        payload = decode_token(token)
        if payload.get("token_type") != "access":
            raise InvalidTokenError("Not an access token")
        return payload
    except (ExpiredSignatureError, InvalidTokenError):
        return None


def verify_refresh_token(token: str) -> Optional[Dict[str, Any]]:
    """
    Verify a refresh token

    Args:
        token: JWT refresh token

    Returns:
        Decoded payload if valid, None otherwise
    """
    try:
        payload = decode_token(token)
        if payload.get("token_type") != "refresh":
            raise InvalidTokenError("Not a refresh token")
        return payload
    except (ExpiredSignatureError, InvalidTokenError):
        return None


def get_token_expiration(token: str) -> Optional[datetime]:
    """
    Get expiration time of a token

    Args:
        token: JWT token

    Returns:
        Expiration datetime or None if invalid
    """
    try:
        payload = decode_token(token)
        exp_timestamp = payload.get("exp")
        if exp_timestamp:
            return datetime.fromtimestamp(exp_timestamp)
        return None
    except (ExpiredSignatureError, InvalidTokenError):
        return None


def is_token_expired(token: str) -> bool:
    """
    Check if a token is expired

    Args:
        token: JWT token

    Returns:
        True if expired, False otherwise
    """
    try:
        decode_token(token)
        return False
    except ExpiredSignatureError:
        return True
    except InvalidTokenError:
        return True


def extract_user_from_token(token: str) -> Optional[Dict[str, Any]]:
    """
    Extract user information from token

    Args:
        token: JWT token (access or refresh)

    Returns:
        Dictionary with user_id, email, role or None if invalid
    """
    try:
        payload = decode_token(token)
        return {
            "user_id": payload.get("user_id"),
            "email": payload.get("email"),
            "role": payload.get("role"),
        }
    except (ExpiredSignatureError, InvalidTokenError):
        return None


def generate_session_id() -> str:
    """
    Generate a unique session ID

    Returns:
        Random URL-safe session ID
    """
    return secrets.token_urlsafe(32)


# Token validation results
class TokenValidationResult:
    """Result of token validation"""

    def __init__(self, valid: bool, payload: Optional[Dict] = None, error: Optional[str] = None):
        self.valid = valid
        self.payload = payload
        self.error = error

    def __bool__(self):
        return self.valid


def validate_token(token: str, token_type: str = "access") -> TokenValidationResult:
    """
    Comprehensive token validation

    Args:
        token: JWT token to validate
        token_type: Expected token type ("access" or "refresh")

    Returns:
        TokenValidationResult object
    """
    try:
        payload = decode_token(token)

        # Check token type
        if payload.get("token_type") != token_type:
            return TokenValidationResult(
                valid=False, error=f"Expected {token_type} token, got {payload.get('token_type')}"
            )

        # Check required fields
        required_fields = ["user_id", "email", "role", "exp", "iat"]
        missing_fields = [field for field in required_fields if field not in payload]
        if missing_fields:
            return TokenValidationResult(
                valid=False, error=f"Missing required fields: {', '.join(missing_fields)}"
            )

        return TokenValidationResult(valid=True, payload=payload)

    except ExpiredSignatureError:
        return TokenValidationResult(valid=False, error="Token has expired")
    except InvalidTokenError as e:
        return TokenValidationResult(valid=False, error=f"Invalid token: {str(e)}")
    except Exception as e:
        return TokenValidationResult(valid=False, error=f"Validation error: {str(e)}")


# Helper for FastAPI dependency injection
def get_current_user_id(token: str) -> Optional[int]:
    """
    Extract user ID from token (for FastAPI dependencies)

    Args:
        token: JWT access token

    Returns:
        User ID or None if invalid
    """
    result = validate_token(token, "access")
    if result.valid:
        return result.payload.get("user_id")
    return None


def get_current_user_role(token: str) -> Optional[str]:
    """
    Extract user role from token (for FastAPI dependencies)

    Args:
        token: JWT access token

    Returns:
        User role or None if invalid
    """
    result = validate_token(token, "access")
    if result.valid:
        return result.payload.get("role")
    return None
