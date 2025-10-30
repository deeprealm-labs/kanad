"""
Google OAuth 2.0 Integration for Kanad Platform
Handles Google Sign-In authentication flow
"""

import os
from typing import Optional, Dict, Any
import requests
from datetime import datetime

# Google OAuth configuration
GOOGLE_CLIENT_ID = os.getenv("GOOGLE_CLIENT_ID", "")
GOOGLE_CLIENT_SECRET = os.getenv("GOOGLE_CLIENT_SECRET", "")
GOOGLE_REDIRECT_URI = os.getenv(
    "GOOGLE_REDIRECT_URI", "http://localhost:3000/auth/google/callback"
)

# Google OAuth endpoints
GOOGLE_AUTH_URL = "https://accounts.google.com/o/oauth2/v2/auth"
GOOGLE_TOKEN_URL = "https://oauth2.googleapis.com/token"
GOOGLE_USERINFO_URL = "https://www.googleapis.com/oauth2/v2/userinfo"

# OAuth scopes
GOOGLE_SCOPES = [
    "https://www.googleapis.com/auth/userinfo.email",
    "https://www.googleapis.com/auth/userinfo.profile",
]


def get_google_auth_url(state: Optional[str] = None) -> str:
    """
    Generate Google OAuth authorization URL

    Args:
        state: Optional state parameter for CSRF protection

    Returns:
        Google OAuth authorization URL
    """
    if not GOOGLE_CLIENT_ID:
        raise ValueError("GOOGLE_CLIENT_ID not configured")

    params = {
        "client_id": GOOGLE_CLIENT_ID,
        "redirect_uri": GOOGLE_REDIRECT_URI,
        "response_type": "code",
        "scope": " ".join(GOOGLE_SCOPES),
        "access_type": "offline",
        "prompt": "select_account",
    }

    if state:
        params["state"] = state

    # Build URL
    query_string = "&".join(f"{key}={value}" for key, value in params.items())
    return f"{GOOGLE_AUTH_URL}?{query_string}"


def exchange_code_for_token(code: str) -> Optional[Dict[str, Any]]:
    """
    Exchange authorization code for access token

    Args:
        code: Authorization code from Google

    Returns:
        Token response dictionary or None if failed
    """
    if not GOOGLE_CLIENT_ID or not GOOGLE_CLIENT_SECRET:
        raise ValueError("Google OAuth credentials not configured")

    try:
        response = requests.post(
            GOOGLE_TOKEN_URL,
            data={
                "code": code,
                "client_id": GOOGLE_CLIENT_ID,
                "client_secret": GOOGLE_CLIENT_SECRET,
                "redirect_uri": GOOGLE_REDIRECT_URI,
                "grant_type": "authorization_code",
            },
            timeout=10,
        )

        if response.status_code == 200:
            return response.json()
        else:
            print(f"✗ Google token exchange failed: {response.status_code}")
            print(f"   Response: {response.text}")
            return None

    except Exception as e:
        print(f"✗ Error exchanging Google code for token: {str(e)}")
        return None


def get_google_user_info(access_token: str) -> Optional[Dict[str, Any]]:
    """
    Get user information from Google using access token

    Args:
        access_token: Google OAuth access token

    Returns:
        User info dictionary or None if failed
    """
    try:
        response = requests.get(
            GOOGLE_USERINFO_URL,
            headers={"Authorization": f"Bearer {access_token}"},
            timeout=10,
        )

        if response.status_code == 200:
            return response.json()
        else:
            print(f"✗ Google userinfo request failed: {response.status_code}")
            return None

    except Exception as e:
        print(f"✗ Error fetching Google user info: {str(e)}")
        return None


def verify_google_token(id_token: str) -> Optional[Dict[str, Any]]:
    """
    Verify Google ID token (alternative method)

    Args:
        id_token: Google ID token from client-side

    Returns:
        Decoded token payload or None if invalid
    """
    try:
        # Verify token with Google's tokeninfo endpoint
        response = requests.get(
            f"https://oauth2.googleapis.com/tokeninfo?id_token={id_token}",
            timeout=10,
        )

        if response.status_code == 200:
            payload = response.json()

            # Verify audience (client ID)
            if payload.get("aud") != GOOGLE_CLIENT_ID:
                print("✗ Invalid token audience")
                return None

            # Check expiration
            exp = int(payload.get("exp", 0))
            if exp < datetime.utcnow().timestamp():
                print("✗ Token expired")
                return None

            return payload
        else:
            print(f"✗ Token verification failed: {response.status_code}")
            return None

    except Exception as e:
        print(f"✗ Error verifying Google token: {str(e)}")
        return None


class GoogleOAuthHandler:
    """
    Google OAuth handler class for managing authentication flow
    """

    def __init__(self, client_id: str = None, client_secret: str = None, redirect_uri: str = None):
        self.client_id = client_id or GOOGLE_CLIENT_ID
        self.client_secret = client_secret or GOOGLE_CLIENT_SECRET
        self.redirect_uri = redirect_uri or GOOGLE_REDIRECT_URI

        if not self.client_id or not self.client_secret:
            raise ValueError("Google OAuth credentials not configured")

    def get_authorization_url(self, state: Optional[str] = None) -> str:
        """Generate authorization URL"""
        return get_google_auth_url(state)

    def authenticate(self, code: str) -> Optional[Dict[str, Any]]:
        """
        Complete authentication flow

        Args:
            code: Authorization code from Google

        Returns:
            User information dictionary or None if failed
        """
        # Exchange code for token
        token_response = exchange_code_for_token(code)
        if not token_response:
            return None

        access_token = token_response.get("access_token")
        if not access_token:
            print("✗ No access token in response")
            return None

        # Get user info
        user_info = get_google_user_info(access_token)
        if not user_info:
            return None

        # Extract relevant information
        return {
            "google_id": user_info.get("id"),
            "email": user_info.get("email"),
            "email_verified": user_info.get("verified_email", False),
            "full_name": user_info.get("name"),
            "given_name": user_info.get("given_name"),
            "family_name": user_info.get("family_name"),
            "avatar_url": user_info.get("picture"),
            "locale": user_info.get("locale"),
        }

    def verify_id_token(self, id_token: str) -> Optional[Dict[str, Any]]:
        """
        Verify ID token from client-side Google Sign-In

        Args:
            id_token: Google ID token

        Returns:
            User information or None if invalid
        """
        payload = verify_google_token(id_token)
        if not payload:
            return None

        return {
            "google_id": payload.get("sub"),
            "email": payload.get("email"),
            "email_verified": payload.get("email_verified", False),
            "full_name": payload.get("name"),
            "given_name": payload.get("given_name"),
            "family_name": payload.get("family_name"),
            "avatar_url": payload.get("picture"),
            "locale": payload.get("locale"),
        }


def create_or_update_google_user(db, google_user_info: Dict[str, Any], access_key: Optional[str] = None) -> Optional[Any]:
    """
    Create or update user from Google OAuth data

    Args:
        db: Database session
        google_user_info: User information from Google
        access_key: Access key for new user registration (required for new users)

    Returns:
        User object or None if failed
    """
    from api.core.database_postgres import User, UserRole, AccessKey

    google_id = google_user_info.get("google_id")
    email = google_user_info.get("email")

    if not google_id or not email:
        print("✗ Missing required Google user info")
        return None

    # Check if user exists by Google ID
    user = db.query(User).filter(User.google_id == google_id).first()

    if user:
        # Update existing user
        user.email = email
        user.full_name = google_user_info.get("full_name") or user.full_name
        user.avatar_url = google_user_info.get("avatar_url") or user.avatar_url
        user.is_verified = True  # Google emails are pre-verified
        user.last_login = datetime.utcnow()
        db.commit()
        db.refresh(user)
        print(f"✓ Updated existing Google user: {email}")
        return user

    # Check if user exists by email
    user = db.query(User).filter(User.email == email).first()

    if user:
        # Link Google account to existing user
        user.google_id = google_id
        user.avatar_url = google_user_info.get("avatar_url") or user.avatar_url
        user.is_verified = True
        user.last_login = datetime.utcnow()
        db.commit()
        db.refresh(user)
        print(f"✓ Linked Google account to existing user: {email}")
        return user

    # Create new user - REQUIRES ACCESS KEY
    if not access_key:
        print(f"✗ Access key required for new Google user: {email}")
        return None

    # Validate access key
    key = db.query(AccessKey).filter(
        AccessKey.key == access_key,
        AccessKey.is_active == True
    ).first()

    if not key:
        print(f"✗ Invalid access key for Google user: {email}")
        return None

    # Check if key has remaining uses
    if key.max_uses is not None and key.used_count >= key.max_uses:
        print(f"✗ Access key exhausted for Google user: {email}")
        return None

    # Create user with access key
    user = User(
        email=email,
        google_id=google_id,
        full_name=google_user_info.get("full_name"),
        avatar_url=google_user_info.get("avatar_url"),
        password_hash="",  # No password for OAuth users
        role=UserRole.USER,
        is_verified=True,  # Google emails are verified
        is_active=True,
        access_key_id=key.id,  # Link to access key
        last_login=datetime.utcnow(),
    )

    # Increment key usage
    key.used_count += 1

    db.add(user)
    db.commit()
    db.refresh(user)

    print(f"✓ Created new Google user with access key: {email}")
    return user


def is_google_oauth_configured() -> bool:
    """
    Check if Google OAuth is properly configured

    Returns:
        True if configured, False otherwise
    """
    return bool(GOOGLE_CLIENT_ID and GOOGLE_CLIENT_SECRET)


def get_oauth_config_status() -> Dict[str, Any]:
    """
    Get OAuth configuration status

    Returns:
        Dictionary with configuration status
    """
    return {
        "google_oauth_enabled": is_google_oauth_configured(),
        "client_id_configured": bool(GOOGLE_CLIENT_ID),
        "client_secret_configured": bool(GOOGLE_CLIENT_SECRET),
        "redirect_uri": GOOGLE_REDIRECT_URI,
    }
