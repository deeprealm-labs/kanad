"""
Molecule creation and validation endpoints
"""

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import List, Optional
import numpy as np

from kanad.core.atom import Atom
from kanad.core.molecule import Molecule
from kanad.bonds import BondFactory

router = APIRouter()


class AtomInput(BaseModel):
    symbol: str
    x: float
    y: float
    z: float


class MoleculeCreateRequest(BaseModel):
    smiles: Optional[str] = None
    atoms: Optional[List[AtomInput]] = None
    basis: str = "sto-3g"
    charge: int = 0
    multiplicity: int = 1


class ValidateSmilesRequest(BaseModel):
    smiles: str


@router.post("/create")
async def create_molecule(request: MoleculeCreateRequest):
    """
    Create a molecule from SMILES or atom coordinates.

    Returns basic molecule properties.
    """
    try:
        if request.smiles:
            # Create from SMILES
            molecule = Molecule.from_smiles(
                request.smiles,
                basis=request.basis,
                charge=request.charge,
                spin=request.multiplicity - 1
            )
        elif request.atoms:
            # Create from atoms
            atoms = [
                Atom(
                    atom.symbol,
                    position=np.array([atom.x, atom.y, atom.z])
                )
                for atom in request.atoms
            ]
            molecule = Molecule(
                atoms,
                basis=request.basis,
                charge=request.charge,
                spin=request.multiplicity - 1
            )
        else:
            raise HTTPException(
                status_code=400,
                detail="Must provide either SMILES or atoms"
            )

        return {
            "formula": molecule.formula,
            "n_atoms": molecule.n_atoms,
            "n_electrons": molecule.n_electrons,
            "n_orbitals": molecule.n_orbitals,
            "basis": request.basis,
            "charge": request.charge,
            "multiplicity": request.multiplicity,
            "success": True
        }

    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.post("/validate-smiles")
async def validate_smiles(request: ValidateSmilesRequest):
    """Validate a SMILES string."""
    try:
        # Try to create molecule
        molecule = Molecule.from_smiles(request.smiles, basis='sto-3g')

        return {
            "valid": True,
            "molecular_formula": molecule.formula,
            "n_atoms": molecule.n_atoms,
            "n_electrons": molecule.n_electrons,
            "message": "Valid SMILES string"
        }

    except Exception as e:
        return {
            "valid": False,
            "message": str(e)
        }


@router.post("/bond-info")
async def get_bond_info(atom_1: str, atom_2: str):
    """Get quick bond information without expensive calculations."""
    try:
        info = BondFactory.quick_bond_info(atom_1, atom_2)
        return info

    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))
