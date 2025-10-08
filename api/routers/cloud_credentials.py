"""
Cloud credentials API endpoints.

Manage cloud provider credentials for IBM Quantum and BlueQubit.
"""

from fastapi import APIRouter, HTTPException, Depends
from sqlalchemy.orm import Session
from pydantic import BaseModel
import json

from api.database import get_db
from api.models.cloud_credentials import CloudCredentials


router = APIRouter(prefix="/api/v1/cloud-credentials", tags=["Cloud Credentials"])


# ===== Pydantic Models =====

class IBMCredentialsRequest(BaseModel):
    """Request model for IBM Quantum credentials."""
    crn: str
    api_key: str


class BlueQubitCredentialsRequest(BaseModel):
    """Request model for BlueQubit credentials."""
    api_token: str


class CredentialsResponse(BaseModel):
    """Response model for credentials status."""
    provider: str
    configured: bool
    message: str


# ===== Endpoints =====

@router.post("/ibm", response_model=CredentialsResponse)
async def save_ibm_credentials(
    request: IBMCredentialsRequest,
    db: Session = Depends(get_db)
):
    """
    Save IBM Quantum credentials.

    Stores CRN and API key for IBM Quantum backend access.
    """
    try:
        # Create credentials JSON
        credentials_data = {
            "crn": request.crn,
            "api_key": request.api_key
        }

        # Check if credentials already exist
        existing = db.query(CloudCredentials).filter(
            CloudCredentials.provider == "ibm",
            CloudCredentials.user_id == 1
        ).first()

        if existing:
            # Update existing credentials
            existing.credentials = json.dumps(credentials_data)
            db.commit()
            message = "IBM Quantum credentials updated successfully"
        else:
            # Create new credentials
            credentials = CloudCredentials(
                provider="ibm",
                user_id=1,
                credentials=json.dumps(credentials_data)
            )
            db.add(credentials)
            db.commit()
            message = "IBM Quantum credentials saved successfully"

        return CredentialsResponse(
            provider="ibm",
            configured=True,
            message=message
        )

    except Exception as e:
        db.rollback()
        raise HTTPException(status_code=500, detail=f"Failed to save credentials: {str(e)}")


@router.post("/bluequbit", response_model=CredentialsResponse)
async def save_bluequbit_credentials(
    request: BlueQubitCredentialsRequest,
    db: Session = Depends(get_db)
):
    """
    Save BlueQubit credentials.

    Stores API token for BlueQubit backend access.
    """
    try:
        # Create credentials JSON
        credentials_data = {
            "api_token": request.api_token
        }

        # Check if credentials already exist
        existing = db.query(CloudCredentials).filter(
            CloudCredentials.provider == "bluequbit",
            CloudCredentials.user_id == 1
        ).first()

        if existing:
            # Update existing credentials
            existing.credentials = json.dumps(credentials_data)
            db.commit()
            message = "BlueQubit credentials updated successfully"
        else:
            # Create new credentials
            credentials = CloudCredentials(
                provider="bluequbit",
                user_id=1,
                credentials=json.dumps(credentials_data)
            )
            db.add(credentials)
            db.commit()
            message = "BlueQubit credentials saved successfully"

        return CredentialsResponse(
            provider="bluequbit",
            configured=True,
            message=message
        )

    except Exception as e:
        db.rollback()
        raise HTTPException(status_code=500, detail=f"Failed to save credentials: {str(e)}")


@router.get("/ibm")
async def get_ibm_credentials_status(db: Session = Depends(get_db)):
    """
    Get IBM Quantum credentials status.

    Returns whether credentials are configured (not the actual credentials).
    """
    try:
        credentials = db.query(CloudCredentials).filter(
            CloudCredentials.provider == "ibm",
            CloudCredentials.user_id == 1
        ).first()

        if not credentials:
            return {
                "provider": "ibm",
                "configured": False,
                "message": "No IBM Quantum credentials configured"
            }

        creds_data = json.loads(credentials.credentials)

        return {
            "provider": "ibm",
            "configured": True,
            "has_api_key": bool(creds_data.get("api_key")),
            "has_crn": bool(creds_data.get("crn")),
            "created_at": credentials.created_at.isoformat() if credentials.created_at else None,
            "updated_at": credentials.updated_at.isoformat() if credentials.updated_at else None
        }

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to get credentials status: {str(e)}")


@router.get("/bluequbit")
async def get_bluequbit_credentials_status(db: Session = Depends(get_db)):
    """
    Get BlueQubit credentials status.

    Returns whether credentials are configured (not the actual credentials).
    """
    try:
        credentials = db.query(CloudCredentials).filter(
            CloudCredentials.provider == "bluequbit",
            CloudCredentials.user_id == 1
        ).first()

        if not credentials:
            return {
                "provider": "bluequbit",
                "configured": False,
                "message": "No BlueQubit credentials configured"
            }

        creds_data = json.loads(credentials.credentials)

        return {
            "provider": "bluequbit",
            "configured": True,
            "has_api_token": bool(creds_data.get("api_token")),
            "created_at": credentials.created_at.isoformat() if credentials.created_at else None,
            "updated_at": credentials.updated_at.isoformat() if credentials.updated_at else None
        }

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to get credentials status: {str(e)}")


@router.get("/status")
async def get_credentials_status(db: Session = Depends(get_db)):
    """
    Get status of all cloud credentials.

    Returns whether IBM and BlueQubit credentials are configured.
    """
    try:
        ibm_configured = db.query(CloudCredentials).filter(
            CloudCredentials.provider == "ibm",
            CloudCredentials.user_id == 1
        ).first() is not None

        bluequbit_configured = db.query(CloudCredentials).filter(
            CloudCredentials.provider == "bluequbit",
            CloudCredentials.user_id == 1
        ).first() is not None

        return {
            "ibm": {
                "configured": ibm_configured
            },
            "bluequbit": {
                "configured": bluequbit_configured
            }
        }

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to get credentials status: {str(e)}")


@router.delete("/ibm")
async def delete_ibm_credentials(db: Session = Depends(get_db)):
    """Delete IBM Quantum credentials."""
    try:
        credentials = db.query(CloudCredentials).filter(
            CloudCredentials.provider == "ibm",
            CloudCredentials.user_id == 1
        ).first()

        if not credentials:
            raise HTTPException(status_code=404, detail="IBM credentials not found")

        db.delete(credentials)
        db.commit()

        return {"message": "IBM Quantum credentials deleted successfully"}

    except HTTPException:
        raise
    except Exception as e:
        db.rollback()
        raise HTTPException(status_code=500, detail=f"Failed to delete credentials: {str(e)}")


@router.delete("/bluequbit")
async def delete_bluequbit_credentials(db: Session = Depends(get_db)):
    """Delete BlueQubit credentials."""
    try:
        credentials = db.query(CloudCredentials).filter(
            CloudCredentials.provider == "bluequbit",
            CloudCredentials.user_id == 1
        ).first()

        if not credentials:
            raise HTTPException(status_code=404, detail="BlueQubit credentials not found")

        db.delete(credentials)
        db.commit()

        return {"message": "BlueQubit credentials deleted successfully"}

    except HTTPException:
        raise
    except Exception as e:
        db.rollback()
        raise HTTPException(status_code=500, detail=f"Failed to delete credentials: {str(e)}")
