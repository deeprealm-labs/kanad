"""
Molecules router - Molecule creation and management.
"""

from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy.orm import Session
from typing import List
import uuid
import asyncio

from core.database import get_db
from core.models import MoleculeCreate, MoleculeResponse, SMILESCreate
from db.models import User, Molecule
from utils.auth import get_current_user
from services.computation_service import ComputationService

router = APIRouter()
comp_service = ComputationService()


@router.post("/create", response_model=MoleculeResponse, status_code=status.HTTP_201_CREATED)
async def create_molecule(
    mol_data: MoleculeCreate,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Create a molecule from atoms.

    Supports multiple creation methods: atoms, SMILES, library, XYZ.
    """
    try:
        if mol_data.method.value == "atoms":
            # Create from atom list
            atoms = mol_data.data['atoms']
            basis = mol_data.data.get('basis', 'sto-3g')
            charge = mol_data.data.get('charge', 0)
            multiplicity = mol_data.data.get('multiplicity', 1)

            result = await comp_service.create_molecule_from_atoms(
                atoms=atoms,
                basis=basis,
                charge=charge,
                multiplicity=multiplicity
            )

        elif mol_data.method.value == "smiles":
            # Create from SMILES
            smiles = mol_data.data['smiles']
            basis = mol_data.data.get('basis', 'sto-3g')

            result = await comp_service.create_molecule_from_smiles(
                smiles=smiles,
                basis=basis
            )

        else:
            raise HTTPException(
                status_code=status.HTTP_501_NOT_IMPLEMENTED,
                detail=f"Method {mol_data.method} not yet implemented"
            )

        # Save to database
        molecule_id = uuid.uuid4()
        molecule = Molecule(
            molecule_id=molecule_id,
            user_id=current_user.user_id,
            name=mol_data.data.get('name'),
            formula=result['formula'],
            smiles=result.get('smiles'),
            geometry=result['geometry'],
            basis=basis,
            charge=mol_data.data.get('charge', 0),
            multiplicity=mol_data.data.get('multiplicity', 1),
            n_electrons=result['n_electrons'],
            n_orbitals=result['n_orbitals'],
            n_qubits=result['n_qubits'],
            lewis_structure=result.get('lewis_structure')
        )

        db.add(molecule)
        db.commit()
        db.refresh(molecule)

        return {
            "molecule_id": str(molecule.molecule_id),
            "name": molecule.name,
            "formula": molecule.formula,
            "geometry": molecule.geometry,
            "lewis_structure": molecule.lewis_structure,
            "n_electrons": molecule.n_electrons,
            "n_orbitals": molecule.n_orbitals,
            "n_qubits": molecule.n_qubits,
            "created_at": molecule.created_at
        }

    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to create molecule: {str(e)}"
        )


@router.post("/from-smiles", response_model=MoleculeResponse)
async def create_from_smiles(
    smiles_data: SMILESCreate,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Create molecule from SMILES string.

    Optionally optimizes geometry.
    """
    try:
        result = await comp_service.create_molecule_from_smiles(
            smiles=smiles_data.smiles,
            basis=smiles_data.basis.value,
            optimize_geometry=smiles_data.optimize_geometry
        )

        # Save to database
        molecule_id = uuid.uuid4()
        molecule = Molecule(
            molecule_id=molecule_id,
            user_id=current_user.user_id,
            formula=result['formula'],
            smiles=smiles_data.smiles,
            geometry=result['geometry'],
            basis=smiles_data.basis.value,
            n_electrons=result['n_electrons'],
            n_orbitals=result['n_orbitals'],
            n_qubits=result['n_qubits']
        )

        db.add(molecule)
        db.commit()

        return {
            "molecule_id": str(molecule.molecule_id),
            "formula": molecule.formula,
            "smiles": molecule.smiles,
            "geometry": molecule.geometry,
            "n_electrons": molecule.n_electrons,
            "n_orbitals": molecule.n_orbitals,
            "n_qubits": molecule.n_qubits,
            "created_at": molecule.created_at
        }

    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to create molecule from SMILES: {str(e)}"
        )


@router.get("/{molecule_id}", response_model=MoleculeResponse)
async def get_molecule(
    molecule_id: str,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """Get molecule by ID."""
    molecule = db.query(Molecule).filter(
        Molecule.molecule_id == molecule_id,
        Molecule.user_id == current_user.user_id
    ).first()

    if not molecule:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Molecule not found"
        )

    return {
        "molecule_id": str(molecule.molecule_id),
        "name": molecule.name,
        "formula": molecule.formula,
        "geometry": molecule.geometry,
        "lewis_structure": molecule.lewis_structure,
        "n_electrons": molecule.n_electrons,
        "n_orbitals": molecule.n_orbitals,
        "n_qubits": molecule.n_qubits,
        "created_at": molecule.created_at
    }


@router.get("/", response_model=List[MoleculeResponse])
async def list_molecules(
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db),
    limit: int = 50
):
    """List user's molecules."""
    molecules = db.query(Molecule).filter(
        Molecule.user_id == current_user.user_id
    ).order_by(Molecule.created_at.desc()).limit(limit).all()

    return [
        {
            "molecule_id": str(mol.molecule_id),
            "name": mol.name,
            "formula": mol.formula,
            "geometry": mol.geometry,
            "n_electrons": mol.n_electrons,
            "n_orbitals": mol.n_orbitals,
            "n_qubits": mol.n_qubits,
            "created_at": mol.created_at
        }
        for mol in molecules
    ]
