"""
Molecules API router.

Handles SMILES validation and molecule library.
"""

from fastapi import APIRouter, HTTPException
from typing import List, Dict, Any

from api.utils.validators import SMILESValidation
from kanad.io.smiles_parser import validate_smiles, smiles_to_formula

router = APIRouter(prefix="/molecules", tags=["molecules"])


# Predefined molecule library
MOLECULE_LIBRARY = [
    {
        "id": "h2",
        "name": "Hydrogen",
        "smiles": "[H][H]",
        "formula": "H2",
        "category": "diatomic",
        "description": "Simplest molecule, commonly used for testing"
    },
    {
        "id": "h2o",
        "name": "Water",
        "smiles": "O",
        "formula": "H2O",
        "category": "small_molecules",
        "description": "Water molecule"
    },
    {
        "id": "ch4",
        "name": "Methane",
        "smiles": "C",
        "formula": "CH4",
        "category": "hydrocarbons",
        "description": "Simplest hydrocarbon"
    },
    {
        "id": "nh3",
        "name": "Ammonia",
        "smiles": "N",
        "formula": "NH3",
        "category": "small_molecules",
        "description": "Ammonia molecule"
    },
    {
        "id": "co2",
        "name": "Carbon Dioxide",
        "smiles": "O=C=O",
        "formula": "CO2",
        "category": "small_molecules",
        "description": "Carbon dioxide molecule"
    },
    {
        "id": "ethanol",
        "name": "Ethanol",
        "smiles": "CCO",
        "formula": "C2H6O",
        "category": "organic",
        "description": "Ethanol (drinking alcohol)"
    },
    {
        "id": "ethane",
        "name": "Ethane",
        "smiles": "CC",
        "formula": "C2H6",
        "category": "hydrocarbons",
        "description": "Two-carbon alkane"
    },
    {
        "id": "propane",
        "name": "Propane",
        "smiles": "CCC",
        "formula": "C3H8",
        "category": "hydrocarbons",
        "description": "Three-carbon alkane"
    },
    {
        "id": "benzene",
        "name": "Benzene",
        "smiles": "c1ccccc1",
        "formula": "C6H6",
        "category": "aromatic",
        "description": "Aromatic hydrocarbon with resonance"
    },
    {
        "id": "acetylene",
        "name": "Acetylene",
        "smiles": "C#C",
        "formula": "C2H2",
        "category": "hydrocarbons",
        "description": "Triple bond hydrocarbon"
    },
    {
        "id": "formaldehyde",
        "name": "Formaldehyde",
        "smiles": "C=O",
        "formula": "CH2O",
        "category": "organic",
        "description": "Simplest aldehyde"
    },
    {
        "id": "acetic_acid",
        "name": "Acetic Acid",
        "smiles": "CC(=O)O",
        "formula": "C2H4O2",
        "category": "organic",
        "description": "Vinegar component"
    },
    {
        "id": "lih",
        "name": "Lithium Hydride",
        "smiles": "[LiH]",
        "formula": "LiH",
        "category": "inorganic",
        "description": "Ionic bond example"
    },
    {
        "id": "beh2",
        "name": "Beryllium Hydride",
        "smiles": "[BeH2]",
        "formula": "BeH2",
        "category": "inorganic",
        "description": "Linear molecule"
    },
]


@router.post("/validate")
def validate_smiles_string(validation: SMILESValidation):
    """
    Validate a SMILES string.

    Returns validation status, formula, and molecule info if valid.
    """
    is_valid, error = validate_smiles(validation.smiles)

    if not is_valid:
        return {
            "valid": False,
            "smiles": validation.smiles,
            "error": error
        }

    # Get formula
    try:
        formula = smiles_to_formula(validation.smiles)
    except Exception as e:
        formula = None

    return {
        "valid": True,
        "smiles": validation.smiles,
        "formula": formula,
        "error": None
    }


@router.get("/library")
def get_molecule_library(
    category: str = None
) -> List[Dict[str, Any]]:
    """
    Get predefined molecule library.

    Optionally filter by category.
    """
    if category:
        filtered = [mol for mol in MOLECULE_LIBRARY if mol["category"] == category]
        return filtered

    return MOLECULE_LIBRARY


@router.get("/library/categories")
def get_categories():
    """
    Get list of available molecule categories.
    """
    categories = list(set(mol["category"] for mol in MOLECULE_LIBRARY))
    return {
        "categories": sorted(categories)
    }


@router.get("/library/{molecule_id}")
def get_molecule_from_library(molecule_id: str):
    """
    Get specific molecule from library by ID.
    """
    for mol in MOLECULE_LIBRARY:
        if mol["id"] == molecule_id:
            return mol

    raise HTTPException(status_code=404, detail=f"Molecule {molecule_id} not found in library")
