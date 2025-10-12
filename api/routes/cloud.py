"""
Cloud backend management endpoints
"""

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
import json
from datetime import datetime

from api.core.database import get_db

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
            {
                "provider": "bluequbit",
                "name": "BlueQubit",
                "type": "gpu_simulator",
                "available": True
            }
        ]
    }


@router.post("/credentials")
async def store_credentials(request: CredentialsUpdate):
    """Store cloud provider credentials."""
    try:
        with get_db() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                INSERT INTO cloud_credentials (provider, credentials, updated_at)
                VALUES (?, ?, ?)
                ON CONFLICT(provider) DO UPDATE SET
                    credentials = excluded.credentials,
                    updated_at = excluded.updated_at
            """, (
                request.provider,
                json.dumps(request.credentials),
                datetime.utcnow().isoformat()
            ))
            conn.commit()

        return {"message": f"Credentials for {request.provider} stored successfully"}

    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


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
