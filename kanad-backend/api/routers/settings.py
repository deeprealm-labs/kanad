"""Settings router - User settings."""
from fastapi import APIRouter
router = APIRouter()

@router.get("/defaults")
async def get_defaults():
    """Get default settings."""
    return {"computation": {}, "analysis": {}}
