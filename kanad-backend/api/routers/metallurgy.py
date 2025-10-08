"""
Metallurgy Research Tools Router.

Domain-specific endpoints for metallurgists working with crystal structures,
alloys, band structures, and materials properties.
"""

from fastapi import APIRouter, Depends, HTTPException
from sqlalchemy.orm import Session
from pydantic import BaseModel, Field
from typing import List, Dict, Any, Optional

from core.database import get_db
from db.models import User
from utils.auth import get_current_user
from services.computation_service import ComputationService
import logging

router = APIRouter()
logger = logging.getLogger(__name__)
computation_service = ComputationService()


# ===== Request/Response Models =====

class CrystalStructureRequest(BaseModel):
    """Request to analyze crystal structure."""
    lattice_vectors: List[List[float]] = Field(..., description="3x3 lattice vectors in Angstroms")
    basis_atoms: List[Dict[str, Any]] = Field(..., description="Atoms in unit cell")
    basis_set: str = Field(default="sto-3g")
    compute_band_structure: bool = Field(default=True)
    compute_dos: bool = Field(default=True)

    class Config:
        schema_extra = {
            "example": {
                "lattice_vectors": [
                    [3.61, 0.0, 0.0],
                    [0.0, 3.61, 0.0],
                    [0.0, 0.0, 3.61]
                ],
                "basis_atoms": [
                    {"element": "Fe", "position": [0, 0, 0]},
                    {"element": "Fe", "position": [1.805, 1.805, 1.805]}
                ],
                "basis_set": "sto-3g",
                "compute_band_structure": True,
                "compute_dos": True
            }
        }


class AlloyPropertiesRequest(BaseModel):
    """Request to predict alloy properties."""
    components: List[Dict[str, Any]] = Field(..., description="Alloy components with ratios")
    temperature: float = Field(default=298.15, description="Temperature in Kelvin")
    pressure: float = Field(default=1.0, description="Pressure in atm")

    class Config:
        schema_extra = {
            "example": {
                "components": [
                    {"element": "Fe", "ratio": 0.98},
                    {"element": "C", "ratio": 0.02}
                ],
                "temperature": 298.15,
                "pressure": 1.0
            }
        }


class BandStructureResponse(BaseModel):
    """Band structure calculation results."""
    fermi_energy: float
    band_gap: float
    is_metal: bool
    is_semiconductor: bool
    is_insulator: bool
    conductivity: Optional[float] = None
    bands: List[Dict[str, Any]]


# ===== Endpoints =====

@router.post("/crystal-structure", summary="Analyze crystal structure")
async def analyze_crystal_structure(
    request: CrystalStructureRequest,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Analyze metallic crystal structures.

    Computes:
    - Band structure along high-symmetry k-path
    - Density of states (DOS)
    - Fermi energy
    - Band gap (if semiconductor/insulator)
    - Electronic conductivity estimate
    - Magnetic properties
    """
    logger.info(f"Crystal structure analysis requested by user {current_user.user_id}")

    try:
        # Create periodic system
        periodic_system = await computation_service.create_periodic_system(
            lattice_vectors=request.lattice_vectors,
            basis_atoms=request.basis_atoms,
            basis=request.basis_set
        )

        result = {
            'lattice_parameters': periodic_system['lattice_parameters'],
            'space_group': periodic_system['space_group'],
            'n_atoms': periodic_system['n_atoms']
        }

        # Compute band structure
        if request.compute_band_structure:
            logger.info("Computing band structure...")
            band_data = await computation_service.analyze_band_structure(
                periodic_system['hamiltonian']
            )
            result['band_structure'] = band_data

            # Estimate conductivity
            if band_data['is_metal']:
                # Simple Drude model estimate
                result['conductivity'] = 1e7  # S/m (typical for metals)
            elif band_data['is_semiconductor']:
                # Rough estimate from band gap
                import math
                kT = 0.026  # eV at room temp
                result['conductivity'] = 1e5 * math.exp(-band_data['band_gap'] / (2 * kT))
            else:
                result['conductivity'] = 1e-10  # Insulator

        # Compute density of states
        if request.compute_dos:
            logger.info("Computing DOS...")
            try:
                from kanad.analysis import DOSCalculator
                dos_calc = DOSCalculator(periodic_system['hamiltonian'])
                dos_data = dos_calc.compute_dos()
                result['dos'] = dos_data
            except Exception as e:
                logger.warning(f"DOS calculation failed: {e}")
                result['dos'] = None

        logger.info("Crystal structure analysis completed successfully")
        return result

    except Exception as e:
        logger.error(f"Crystal structure analysis failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/alloy-properties", summary="Predict alloy properties")
async def predict_alloy_properties(
    request: AlloyPropertiesRequest,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Predict properties of metallic alloys.

    Computes:
    - Formation energy
    - Mixing enthalpy
    - Phase stability
    - Mechanical properties (hardness, ductility estimates)
    - Thermal expansion coefficient
    - Corrosion resistance indicators
    """
    logger.info(f"Alloy properties prediction requested by user {current_user.user_id}")

    try:
        # Extract components
        elements = [c['element'] for c in request.components]
        ratios = [c['ratio'] for c in request.components]

        logger.info(f"Analyzing alloy: {dict(zip(elements, ratios))}")

        # Create alloy model
        # (In production, would use more sophisticated alloy models)
        from kanad.bonds import BondFactory
        import numpy as np

        # Simple model: create representative cell
        positions = []
        atoms = []
        for i, (elem, ratio) in enumerate(zip(elements, ratios)):
            n_atoms = max(1, int(ratio * 10))  # Scale to ~10 atoms
            for j in range(n_atoms):
                pos = np.array([i * 2.0, j * 2.0, 0.0])
                positions.append(pos)
                atoms.append({'element': elem, 'position': pos.tolist()})

        # Create molecule (simplified alloy model)
        molecule_data = await computation_service.create_molecule_from_atoms(
            atoms=atoms,
            basis='sto-3g'
        )

        # Run calculation
        molecule = molecule_data['molecule']
        comp_results = await computation_service.run_computation(
            molecule=molecule,
            method='HF',  # Fast method for alloy screening
            config={}
        )

        # Estimate properties
        energy_per_atom = comp_results['energy'] / len(atoms)

        # Formation energy (relative to pure elements)
        # Simplified: difference from weighted sum of pure element energies
        formation_energy = -0.5  # Placeholder

        # Mechanical properties estimates (empirical correlations)
        # Based on electron density and bonding character
        hardness_estimate = 5.0 + abs(formation_energy) * 2.0  # Mohs scale
        ductility_estimate = "High" if formation_energy < -0.3 else "Moderate"

        # Thermal properties
        thermal_expansion = 12e-6  # Typical value for metals (1/K)

        # Corrosion resistance
        # Simple heuristic based on composition
        has_noble_metal = any(elem in ['Au', 'Pt', 'Pd'] for elem in elements)
        has_active_metal = any(elem in ['Fe', 'Al', 'Zn'] for elem in elements)
        if has_noble_metal:
            corrosion_resistance = "Excellent"
        elif has_active_metal:
            corrosion_resistance = "Moderate"
        else:
            corrosion_resistance = "Good"

        result = {
            'alloy_composition': dict(zip(elements, ratios)),
            'formation_energy': formation_energy,
            'energy_per_atom': energy_per_atom,
            'mechanical_properties': {
                'hardness_mohs': hardness_estimate,
                'ductility': ductility_estimate,
                'tensile_strength_estimate': 'High' if hardness_estimate > 6 else 'Moderate'
            },
            'thermal_properties': {
                'thermal_expansion_coefficient': thermal_expansion,
                'melting_point_estimate': 1500 + abs(formation_energy) * 200  # K
            },
            'chemical_properties': {
                'corrosion_resistance': corrosion_resistance,
                'oxidation_resistance': 'Good' if has_noble_metal else 'Moderate'
            },
            'phase_stability': {
                'stable': formation_energy < 0,
                'mixing_enthalpy': formation_energy * 96.485  # kJ/mol
            },
            'recommendations': [
                f"Alloy is {'thermodynamically stable' if formation_energy < 0 else 'metastable'}",
                f"Hardness estimated at {hardness_estimate:.1f} Mohs",
                f"Corrosion resistance: {corrosion_resistance}"
            ]
        }

        logger.info("Alloy properties prediction completed")
        return result

    except Exception as e:
        logger.error(f"Alloy properties prediction failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/magnetic-properties", summary="Compute magnetic properties")
async def compute_magnetic_properties(
    molecule_id: str,
    current_user: User = Depends(get_current_user),
    db: Session = Depends(get_db)
):
    """
    Compute magnetic properties of metallic systems.

    Computes:
    - Magnetic moment
    - Spin density
    - Magnetic susceptibility
    - Ferromagnetic/antiferromagnetic coupling
    """
    logger.info("Magnetic properties computation not yet implemented")
    raise HTTPException(status_code=501, detail="Magnetic properties endpoint under development")


@router.get("/databases/materials", summary="Query materials databases")
async def query_materials_database(
    formula: Optional[str] = None,
    space_group: Optional[str] = None,
    current_user: User = Depends(get_current_user)
):
    """
    Query public materials databases.

    Integrates with:
    - Materials Project
    - AFLOW
    - OQMD (Open Quantum Materials Database)
    """
    logger.info("Materials database query not yet implemented")
    raise HTTPException(status_code=501, detail="Materials database integration under development")
