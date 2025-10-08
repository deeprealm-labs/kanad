"""
Credential retrieval utilities.

Helper functions to retrieve cloud credentials from database or environment variables.
"""

import json
import logging
from typing import Optional, Dict, Any
from sqlalchemy.orm import Session

from api.models.cloud_credentials import CloudCredentials
from api.config import settings

logger = logging.getLogger(__name__)


def get_ibm_credentials(db: Session) -> Optional[Dict[str, Any]]:
    """
    Retrieve IBM Quantum credentials.

    Checks database first, falls back to environment variables.

    Args:
        db: Database session

    Returns:
        Dictionary with api_key and crn, or None if not configured
    """
    # Try database first
    creds = db.query(CloudCredentials).filter(
        CloudCredentials.provider == "ibm",
        CloudCredentials.user_id == 1
    ).first()

    if creds:
        try:
            creds_data = json.loads(creds.credentials)
            logger.debug("Retrieved IBM credentials from database")
            return {
                "api_key": creds_data.get("api_key"),
                "crn": creds_data.get("crn"),
                "channel": creds_data.get("channel", "ibm_quantum")
            }
        except Exception as e:
            logger.error(f"Failed to parse IBM credentials from database: {e}")

    # Fallback to environment variables
    if settings.IBM_QUANTUM_TOKEN:
        logger.debug("Using IBM credentials from environment variables")
        return {
            "api_key": settings.IBM_QUANTUM_TOKEN,
            "crn": settings.IBM_QUANTUM_CRN,
            "channel": "ibm_quantum"
        }

    logger.warning("No IBM credentials configured")
    return None


def get_bluequbit_credentials(db: Session) -> Optional[Dict[str, Any]]:
    """
    Retrieve BlueQubit credentials.

    Checks database first, falls back to environment variables.

    Args:
        db: Database session

    Returns:
        Dictionary with api_token, or None if not configured
    """
    # Try database first
    creds = db.query(CloudCredentials).filter(
        CloudCredentials.provider == "bluequbit",
        CloudCredentials.user_id == 1
    ).first()

    if creds:
        try:
            creds_data = json.loads(creds.credentials)
            logger.debug("Retrieved BlueQubit credentials from database")
            return {
                "api_token": creds_data.get("api_token")
            }
        except Exception as e:
            logger.error(f"Failed to parse BlueQubit credentials from database: {e}")

    # Fallback to environment variables
    if settings.BLUEQUBIT_API_KEY:
        logger.debug("Using BlueQubit credentials from environment variables")
        return {
            "api_token": settings.BLUEQUBIT_API_KEY
        }

    logger.warning("No BlueQubit credentials configured")
    return None
