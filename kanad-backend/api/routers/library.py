"""Library router - Molecule library."""
from fastapi import APIRouter
router = APIRouter()

@router.get("/molecules")
async def get_library():
    """Get molecule library."""
    return {"categories": []}
