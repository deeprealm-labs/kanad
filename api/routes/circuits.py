"""
Circuit preview and visualization endpoints
"""

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import Dict, Any, List, Optional

router = APIRouter()


class AtomInput(BaseModel):
    symbol: str
    x: float
    y: float
    z: float


class MoleculeConfig(BaseModel):
    atoms: Optional[List[AtomInput]] = None
    smiles: Optional[str] = None
    basis: str = "sto-3g"
    charge: int = 0
    multiplicity: int = 1


class BackendConfig(BaseModel):
    method: str = "VQE"
    ansatz: Optional[str] = "hardware_efficient"
    mapper: Optional[str] = "jordan_wigner"


class CircuitPreviewRequest(BaseModel):
    molecule: MoleculeConfig
    configuration: BackendConfig


@router.post("/preview")
async def get_circuit_preview(request: CircuitPreviewRequest):
    """
    Generate circuit preview for given molecule and configuration.

    Returns ASCII diagram and circuit statistics.
    """
    try:
        from kanad.visualization import create_preview_circuit

        # Convert request to dictionaries
        molecule_dict = {
            'basis': request.molecule.basis,
            'charge': request.molecule.charge,
            'multiplicity': request.molecule.multiplicity
        }

        # Add either atoms or SMILES
        if request.molecule.smiles:
            molecule_dict['smiles'] = request.molecule.smiles
        elif request.molecule.atoms:
            molecule_dict['atoms'] = [atom.model_dump() for atom in request.molecule.atoms]

        backend_dict = request.configuration.model_dump()

        # Generate preview
        preview_data = create_preview_circuit(molecule_dict, backend_dict)

        # Normalize response format for frontend
        if preview_data.get('has_circuit'):
            preview = {
                'circuit_diagram': preview_data.get('diagram', ''),
                'statistics': preview_data.get('stats', {}),
                'gates': preview_data.get('gates', []),  # Include structured gates
                'circuit_image': preview_data.get('circuit_image'),  # Include Matplotlib image
                'n_qubits': preview_data.get('n_qubits', 0),
                'n_electrons': preview_data.get('n_electrons', 0),
                'n_parameters': preview_data.get('n_parameters', 0),
                'ansatz_type': preview_data.get('ansatz_type', ''),
                'mapper_type': preview_data.get('mapper_type', ''),
                'method': preview_data.get('method', 'VQE'),
                'note': preview_data.get('note', '')
            }
        else:
            preview = {
                'has_circuit': False,
                'message': preview_data.get('message', 'Circuit preview not available')
            }

        return {
            "success": True,
            "preview": preview
        }

    except Exception as e:
        import traceback
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=f"Failed to generate circuit preview: {str(e)}")
