"""
Cloud backend management endpoints
"""

from fastapi import APIRouter, HTTPException, Depends
from pydantic import BaseModel
from typing import Optional
import json
from datetime import datetime

from api.core.database import get_db
from api.dependencies.auth import get_optional_user, get_current_verified_user
from api.core.database_postgres import User

router = APIRouter()


class CredentialsUpdate(BaseModel):
    provider: str  # 'ibm' or 'bluequbit'
    credentials: dict


@router.get("/backends")
async def get_available_backends():
    """Get list of available quantum backends."""
    return {
        "backends": [
            {
                "provider": "classical",
                "name": "Local Classical Simulator",
                "type": "simulator",
                "available": True
            },
            {
                "provider": "ibm_quantum",
                "name": "IBM Quantum",
                "type": "hardware",
                "backends": [
                    {"name": "ibm_torino", "qubits": 133, "available": True},
                    {"name": "ibm_brisbane", "qubits": 127, "available": True},
                    {"name": "ibm_kyoto", "qubits": 127, "available": True},
                ]
            },
            # DISABLED - BlueQubit backend has freezing issues
            # See IMMEDIATE_FIXES_REQUIRED.md for details
            # {
            #     "provider": "bluequbit",
            #     "name": "BlueQubit",
            #     "type": "gpu_simulator",
            #     "available": True
            # }
        ]
    }


@router.post("/credentials")
async def store_credentials(
    request: CredentialsUpdate,
    current_user: User = Depends(get_current_verified_user)
):
    """
    Store cloud provider credentials.
    Credentials are scoped to the authenticated user.
    """
    try:
        with get_db() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                INSERT INTO cloud_credentials (provider, user_id, credentials, updated_at)
                VALUES (?, ?, ?, ?)
                ON CONFLICT(provider, user_id) DO UPDATE SET
                    credentials = excluded.credentials,
                    updated_at = excluded.updated_at
            """, (
                request.provider,
                current_user.id,
                json.dumps(request.credentials),
                datetime.utcnow().isoformat()
            ))
            conn.commit()

        return {"message": f"Credentials for {request.provider} stored successfully"}

    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/credentials/status")
async def get_credentials_status(
    current_user: Optional[User] = Depends(get_optional_user)
):
    """
    Get configuration status for all cloud providers.
    Returns user's credentials if authenticated, otherwise global credentials.
    """
    with get_db() as conn:
        cursor = conn.cursor()

        user_id = current_user.id if current_user else None

        # Check IBM - try user-specific first, then fall back to global
        cursor.execute("""
            SELECT updated_at FROM cloud_credentials
            WHERE provider = 'ibm' AND (user_id = ? OR user_id IS NULL)
            ORDER BY user_id DESC NULLS LAST
            LIMIT 1
        """, (user_id,))
        ibm_row = cursor.fetchone()

        # Check BlueQubit - try user-specific first, then fall back to global
        cursor.execute("""
            SELECT updated_at FROM cloud_credentials
            WHERE provider = 'bluequbit' AND (user_id = ? OR user_id IS NULL)
            ORDER BY user_id DESC NULLS LAST
            LIMIT 1
        """, (user_id,))
        bluequbit_row = cursor.fetchone()

        return {
            "ibm": {
                "configured": ibm_row is not None,
                "updated_at": ibm_row['updated_at'] if ibm_row else None
            },
            "bluequbit": {
                "configured": bluequbit_row is not None,
                "updated_at": bluequbit_row['updated_at'] if bluequbit_row else None
            }
        }


@router.get("/credentials/{provider}")
async def get_credentials(provider: str):
    """Get stored credentials (returns only metadata, not actual credentials)."""
    with get_db() as conn:
        cursor = conn.cursor()
        cursor.execute(
            "SELECT updated_at FROM cloud_credentials WHERE provider = ?",
            (provider,)
        )
        row = cursor.fetchone()

        if row:
            return {
                "provider": provider,
                "has_credentials": True,
                "updated_at": row['updated_at']
            }

    return {
        "provider": provider,
        "has_credentials": False
    }


@router.delete("/credentials/{provider}")
async def delete_credentials(provider: str):
    """Delete stored credentials for a provider."""
    try:
        with get_db() as conn:
            cursor = conn.cursor()
            cursor.execute(
                "DELETE FROM cloud_credentials WHERE provider = ?",
                (provider,)
            )
            conn.commit()

            if cursor.rowcount == 0:
                raise HTTPException(status_code=404, detail=f"No credentials found for {provider}")

        return {"message": f"Credentials for {provider} deleted successfully"}

    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
