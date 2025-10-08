"""
Settings API router.

Handles user settings and preferences.
"""

from fastapi import APIRouter, Depends
from sqlalchemy.orm import Session

from api.database import get_db
from api.models.settings import UserSettings
from api.utils.validators import SettingsUpdate

router = APIRouter(prefix="/settings", tags=["settings"])


@router.get("/")
def get_settings(db: Session = Depends(get_db)):
    """
    Get user settings.

    Returns default settings if none exist.
    """
    # For now, single-user mode (user_id=1)
    settings = db.query(UserSettings).filter(UserSettings.user_id == 1).first()

    if not settings:
        # Create default settings
        settings = UserSettings(
            user_id=1,
            method="VQE",
            ansatz="ucc",
            mapper="jordan_wigner",
            optimizer="SLSQP",
            backend="classical",
            circuit_optimization=True
        )
        db.add(settings)
        db.commit()
        db.refresh(settings)

    return settings.to_dict()


@router.put("/")
def update_settings(
    update_data: SettingsUpdate,
    db: Session = Depends(get_db)
):
    """
    Update user settings.

    Only provided fields will be updated.
    """
    # Get or create settings
    settings = db.query(UserSettings).filter(UserSettings.user_id == 1).first()

    if not settings:
        settings = UserSettings(user_id=1)
        db.add(settings)

    # Update provided fields
    update_dict = update_data.dict(exclude_unset=True)
    for field, value in update_dict.items():
        if hasattr(settings, field):
            setattr(settings, field, value)

    db.commit()
    db.refresh(settings)

    return settings.to_dict()


@router.delete("/")
def reset_settings(db: Session = Depends(get_db)):
    """
    Reset settings to defaults.
    """
    settings = db.query(UserSettings).filter(UserSettings.user_id == 1).first()

    if settings:
        # Reset to defaults
        settings.method = "VQE"
        settings.ansatz = "ucc"
        settings.mapper = "jordan_wigner"
        settings.optimizer = "SLSQP"
        settings.backend = "classical"
        settings.backend_name = None
        settings.geometry_optimization = False
        settings.orbital_optimization = False
        settings.circuit_optimization = True
        settings.adaptive_vqe = False
        settings.advanced_settings = None

        db.commit()
        db.refresh(settings)

    return {
        "message": "Settings reset to defaults"
    }
