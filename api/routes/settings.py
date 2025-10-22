"""
User settings endpoints
"""

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import Dict, Any
import json
from datetime import datetime

from api.core.database import get_db

router = APIRouter()


class SettingsUpdate(BaseModel):
    settings: Dict[str, Any]


@router.get("/defaults")
async def get_defaults():
    """Get default settings."""
    with get_db() as conn:
        cursor = conn.cursor()
        cursor.execute("SELECT settings FROM user_settings WHERE id = 1")
        row = cursor.fetchone()

        if row:
            return {"settings": json.loads(row['settings'])}

    # Return default settings if none saved
    return {
        "settings": {
            "method": "VQE",
            "ansatz": "hardware_efficient",
            "mapper": "jordan_wigner",
            "optimizer": "COBYLA",  # Default to COBYLA - better for cloud backends
            "backend": "classical",
            "backend_name": "ibm_torino",
            "bluequbitDevice": "gpu",
            "maxIterations": 100,
            "optimization": {
                "geometry": False,
                "orbitals": False,
                "circuit": True,
                "adaptive": False
            }
        }
    }


@router.put("/defaults")
async def update_defaults(request: SettingsUpdate):
    """Update default settings."""
    try:
        with get_db() as conn:
            cursor = conn.cursor()

            # Upsert settings
            cursor.execute("""
                INSERT INTO user_settings (id, settings, updated_at)
                VALUES (1, ?, ?)
                ON CONFLICT(id) DO UPDATE SET
                    settings = excluded.settings,
                    updated_at = excluded.updated_at
            """, (json.dumps(request.settings), datetime.utcnow().isoformat()))

            conn.commit()

        return {"message": "Settings updated successfully"}

    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
