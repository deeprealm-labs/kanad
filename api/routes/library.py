"""
Molecule library endpoints
"""

from fastapi import APIRouter

router = APIRouter()


MOLECULE_LIBRARY = [
    {
        "id": "h2",
        "name": "Hydrogen (H₂)",
        "formula": "H2",
        "smiles": "[H][H]",
        "category": "diatomic",
        "description": "Simplest molecule, commonly used for testing"
    },
    {
        "id": "h2o",
        "name": "Water (H₂O)",
        "formula": "H2O",
        "smiles": "O",
        "category": "polyatomic",
        "description": "Water molecule"
    },
    {
        "id": "nh3",
        "name": "Ammonia (NH₃)",
        "formula": "NH3",
        "smiles": "N",
        "category": "polyatomic",
        "description": "Ammonia molecule"
    },
    {
        "id": "ch4",
        "name": "Methane (CH₄)",
        "formula": "CH4",
        "smiles": "C",
        "category": "organic",
        "description": "Simplest hydrocarbon"
    },
    {
        "id": "lih",
        "name": "Lithium Hydride (LiH)",
        "formula": "LiH",
        "smiles": "[Li]H",
        "category": "ionic",
        "description": "Simple ionic compound"
    },
    {
        "id": "nacl",
        "name": "Sodium Chloride (NaCl)",
        "formula": "NaCl",
        "smiles": "[Na+].[Cl-]",
        "category": "ionic",
        "description": "Table salt"
    },
]


@router.get("/molecules")
async def get_molecule_library():
    """Get pre-defined molecule library."""
    return {
        "molecules": MOLECULE_LIBRARY,
        "total": len(MOLECULE_LIBRARY)
    }


@router.get("/molecules/{molecule_id}")
async def get_molecule(molecule_id: str):
    """Get specific molecule from library."""
    for mol in MOLECULE_LIBRARY:
        if mol['id'] == molecule_id:
            return {"molecule": mol}

    return {"error": "Molecule not found"}, 404
