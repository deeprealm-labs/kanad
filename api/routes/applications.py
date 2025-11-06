"""
Application Domain Platforms API

Exposes 4 domain-specific platforms:
1. DrugDiscoveryPlatform - Pharmaceutical research
2. AlloyDesigner - Materials science
3. CatalystOptimizer - Chemical catalysis
4. MaterialsScout - Optoelectronic materials
"""

from fastapi import APIRouter, HTTPException, Depends
from pydantic import BaseModel, Field
from typing import List, Dict, Any, Optional, Tuple
import logging
import numpy as np

from api.core.database import ExperimentDB
from api.dependencies.auth import get_current_verified_user
from api.core.database_postgres import User

# Import application platforms
# Import service layer
from api.services.application_service import (
    DrugDiscoveryService,
    AlloyDesignService,
    CatalystOptimizationService,
    MaterialsAnalysisService
)

logger = logging.getLogger(__name__)

router = APIRouter()

# ================================
# REQUEST/RESPONSE MODELS
# ================================

# -------- Drug Discovery --------
class DrugAnalysisRequest(BaseModel):
    experiment_id: str
    target_protein: Optional[str] = None
    ph: float = Field(default=7.4, description="pH for binding calculation")
    temperature: float = Field(default=310.15, description="Temperature in Kelvin (default: body temp)")
    calculate_adme: bool = Field(default=True, description="Calculate ADME properties")
    calculate_metabolites: bool = Field(default=False, description="Predict metabolites (slow)")
    ionic_strength: float = Field(default=0.15, description="Ionic strength in mol/L")

class DrugBindingRequest(BaseModel):
    ligand_smiles: str
    target_protein: Optional[str] = None
    ph: float = 7.4
    temperature: float = 310.15
    basis: str = "sto-3g"

class DrugScreeningRequest(BaseModel):
    ligand_smiles_list: List[str]
    target_protein: Optional[str] = None
    ph: float = 7.4
    temperature: float = 310.15
    top_n: int = Field(default=5, description="Return top N candidates")

# -------- Alloy Designer --------
class AlloyDesignRequest(BaseModel):
    composition: Dict[str, float]  # {'Ti': 0.9, 'Al': 0.1}
    temperature: float = Field(default=298.15, description="Temperature in Kelvin")
    pressure: float = Field(default=1.0, description="Pressure in atm")
    calculate_phase_diagram: bool = False
    t_range: Optional[Tuple[float, float]] = None
    p_range: Optional[Tuple[float, float]] = None

class AlloyCompositionPredictionRequest(BaseModel):
    elements: List[str]  # ['Ti', 'Al', 'V']
    target_properties: Dict[str, Any]  # {'strength': '>1000 MPa', 'density': '<4.5'}
    n_candidates: int = Field(default=10, ge=1, le=50)
    temperature: float = 298.15

# -------- Catalyst Optimizer --------
class CatalystOptimizationRequest(BaseModel):
    catalyst_smiles: str
    reactant_smiles: str
    product_smiles: Optional[str] = None
    temperature: float = Field(default=298.15, description="Reaction temperature")
    pressure: float = Field(default=1.0, description="Reaction pressure")
    find_transition_state: bool = Field(default=True, description="Find transition state (slow)")
    calculate_selectivity: bool = Field(default=False, description="Calculate selectivity")

class CatalystScreeningRequest(BaseModel):
    catalyst_smiles_list: List[str]
    reactant_smiles: str
    product_smiles: Optional[str] = None
    temperature: float = 298.15
    top_n: int = Field(default=5, description="Return top N catalysts")

# -------- Materials Scout --------
class MaterialAnalysisRequest(BaseModel):
    experiment_id: str
    calculate_bandgap: bool = True
    calculate_optical: bool = True
    calculate_led_color: bool = False
    wavelength_range: Optional[Tuple[float, float]] = Field(default=(400, 700), description="nm")

class MaterialDesignRequest(BaseModel):
    material_smiles: str
    target_bandgap: Optional[float] = None  # eV
    temperature: float = 298.15
    pressure: float = 1.0
    doping_elements: Optional[List[str]] = None
    doping_concentration: float = 0.01

class MaterialScreeningRequest(BaseModel):
    material_smiles_list: List[str]
    target_bandgap: Optional[float] = None
    target_property: Optional[str] = None  # 'bandgap', 'optical', 'led'
    top_n: int = 5


# ================================
# DRUG DISCOVERY ENDPOINTS
# ================================

@router.post("/drug-discovery/analyze")
async def analyze_drug_candidate(
    request: DrugAnalysisRequest,
    current_user: User = Depends(get_current_verified_user)
):
    """
    Analyze completed experiment as drug candidate.

    Returns druglikeness, ADME properties, binding affinity, etc.
    Competitive with SwissADME, cheaper than Schrödinger.
    """
    try:
        # Get experiment
        experiment = ExperimentDB.get(request.experiment_id)
        if not experiment:
            raise HTTPException(status_code=404, detail="Experiment not found")

        if experiment.get('user_id') != current_user.id:
            raise HTTPException(status_code=403, detail="Not authorized")

        if experiment.get('status') != 'completed':
            raise HTTPException(status_code=400, detail="Experiment not completed")

        # Initialize platform
        # Use service layer

        # Get molecule SMILES
        molecule_config = experiment.get('molecule', {})
        smiles = molecule_config.get('smiles')
        if not smiles:
            raise HTTPException(status_code=400, detail="Experiment has no SMILES string")

        # Perform drug discovery analysis
        results = experiment.get('results', {})

        # Calculate druglikeness and ADME
        analysis_result = DrugDiscoveryService.analyze_drug_candidate(
            smiles=smiles,
            quantum_energy=results.get('energy'),
            homo_lumo_gap=results.get('homo_lumo_gap'),
            dipole_moment=results.get('dipole'),
            molecular_properties=results.get('analysis', {}),
            ph=request.ph,
            temperature=request.temperature,
            calculate_adme=request.calculate_adme,
            predict_metabolites=request.calculate_metabolites
        )

        return {
            "experiment_id": request.experiment_id,
            "application": "drug_discovery",
            "results": analysis_result
        }

    except Exception as e:
        logger.error(f"Drug discovery analysis failed: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"Analysis failed: {str(e)}")


@router.post("/drug-discovery/binding-energy")
async def calculate_binding_energy(
    request: DrugBindingRequest,
    current_user: User = Depends(get_current_verified_user)
):
    """
    Calculate protein-ligand binding affinity.

    Uses quantum calculations for higher accuracy than force fields.
    pH-dependent binding included (unique to Kanad).
    """
    try:
        # Use service layer

        binding_result = DrugDiscoveryService.calculate_binding_affinity(
            ligand_smiles=request.ligand_smiles,
            target_protein=request.target_protein,
            ph=request.ph,
            temperature=request.temperature,
            basis=request.basis
        )

        return {
            "application": "drug_discovery",
            "analysis_type": "binding_energy",
            "results": binding_result
        }

    except Exception as e:
        logger.error(f"Binding calculation failed: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"Calculation failed: {str(e)}")


@router.post("/drug-discovery/screen-library")
async def screen_drug_library(
    request: DrugScreeningRequest,
    current_user: User = Depends(get_current_verified_user)
):
    """
    Screen multiple compounds for drug candidacy.

    Returns ranked list by druglikeness and predicted binding.
    """
    try:
        # Use service layer

        results = DrugDiscoveryService.screen_compound_library(
            smiles_list=request.ligand_smiles_list,
            target_protein=request.target_protein,
            ph=request.ph,
            temperature=request.temperature,
            top_n=request.top_n
        )

        return {
            "application": "drug_discovery",
            "analysis_type": "library_screening",
            "candidates_screened": len(request.ligand_smiles_list),
            "top_candidates": results
        }

    except Exception as e:
        logger.error(f"Library screening failed: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"Screening failed: {str(e)}")


@router.get("/drug-discovery/adme/{smiles}")
async def calculate_adme_properties(
    smiles: str,
    current_user: User = Depends(get_current_verified_user)
):
    """
    Calculate ADME (Absorption, Distribution, Metabolism, Excretion) properties.

    Includes Lipinski Rule of 5, solubility, permeability, etc.
    """
    try:
        # Use service layer

        adme_props = DrugDiscoveryService.calculate_adme_properties(smiles)

        return {
            "application": "drug_discovery",
            "smiles": smiles,
            "adme_properties": adme_props
        }

    except Exception as e:
        logger.error(f"ADME calculation failed: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"ADME calculation failed: {str(e)}")


# ================================
# ALLOY DESIGNER ENDPOINTS
# ================================

@router.post("/alloy/design")
async def design_alloy(
    request: AlloyDesignRequest,
    current_user: User = Depends(get_current_verified_user)
):
    """
    Design and analyze alloy composition.

    Calculates formation energy, mechanical properties, phase stability.
    Competes with CALPHAD/Thermo-Calc at fraction of the cost.
    """
    try:
        # Use service layer

        # Validate composition sums to ~1.0
        total = sum(request.composition.values())
        if abs(total - 1.0) > 0.01:
            raise HTTPException(
                status_code=400,
                detail=f"Composition fractions must sum to 1.0, got {total}"
            )

        # Design alloy
        alloy_result = AlloyDesignService.design_alloy(
            composition=request.composition,
            temperature=request.temperature,
            pressure=request.pressure
        )

        # Optionally calculate phase diagram
        phase_diagram = None
        if request.calculate_phase_diagram and request.t_range and request.p_range:
            phase_diagram = AlloyDesignService.compute_phase_diagram(
                composition=request.composition,
                t_range=request.t_range,
                p_range=request.p_range
            )

        return {
            "application": "alloy_design",
            "composition": request.composition,
            "properties": alloy_result,
            "phase_diagram": phase_diagram
        }

    except Exception as e:
        logger.error(f"Alloy design failed: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"Alloy design failed: {str(e)}")


@router.post("/alloy/predict-composition")
async def predict_optimal_composition(
    request: AlloyCompositionPredictionRequest,
    current_user: User = Depends(get_current_verified_user)
):
    """
    Predict optimal alloy compositions for target properties.

    Explores composition space using quantum calculations.
    Unique capability - CALPHAD can't predict NEW alloys.
    """
    try:
        # Use service layer

        # Validate elements
        if len(request.elements) < 2:
            raise HTTPException(status_code=400, detail="Need at least 2 elements")
        if len(request.elements) > 5:
            raise HTTPException(status_code=400, detail="Maximum 5 elements supported")

        # Predict compositions
        candidates = AlloyDesignService.predict_optimal_compositions(
            elements=request.elements,
            target_properties=request.target_properties,
            n_candidates=request.n_candidates,
            temperature=request.temperature
        )

        return {
            "application": "alloy_design",
            "analysis_type": "composition_prediction",
            "elements": request.elements,
            "target_properties": request.target_properties,
            "candidates": candidates
        }

    except Exception as e:
        logger.error(f"Composition prediction failed: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"Prediction failed: {str(e)}")


@router.post("/alloy/phase-diagram")
async def compute_phase_diagram(
    composition: Dict[str, float],
    t_range: Tuple[float, float] = (300, 1800),
    p_range: Tuple[float, float] = (1, 100),
    resolution: int = 20,
    current_user: User = Depends(get_current_verified_user)
):
    """
    Compute phase diagram for alloy composition.

    Includes pressure dependence (rare in commercial tools).
    """
    try:
        # Use service layer

        phase_diagram = AlloyDesignService.compute_phase_diagram(
            composition=composition,
            t_range=t_range,
            p_range=p_range,
            resolution=resolution
        )

        return {
            "application": "alloy_design",
            "composition": composition,
            "phase_diagram": phase_diagram
        }

    except Exception as e:
        logger.error(f"Phase diagram computation failed: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"Computation failed: {str(e)}")


# ================================
# CATALYST OPTIMIZER ENDPOINTS
# ================================

@router.post("/catalyst/optimize")
async def optimize_catalyst(
    request: CatalystOptimizationRequest,
    current_user: User = Depends(get_current_verified_user)
):
    """
    Optimize catalyst for reaction.

    Calculates activation energy, selectivity, TOF.
    800M-1.9B x speedup vs DFT validated.
    """
    try:
        # Use service layer

        # Optimize catalyst
        result = CatalystOptimizationService.optimize_catalyst(
            catalyst_smiles=request.catalyst_smiles,
            reactant_smiles=request.reactant_smiles,
            product_smiles=request.product_smiles,
            temperature=request.temperature,
            pressure=request.pressure,
            find_transition_state=request.find_transition_state,
            calculate_selectivity=request.calculate_selectivity
        )

        return {
            "application": "catalyst_optimization",
            "catalyst": request.catalyst_smiles,
            "reaction": {
                "reactant": request.reactant_smiles,
                "product": request.product_smiles
            },
            "results": result
        }

    except Exception as e:
        logger.error(f"Catalyst optimization failed: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"Optimization failed: {str(e)}")


@router.post("/catalyst/transition-state")
async def find_transition_state(
    catalyst_smiles: str,
    reactant_smiles: str,
    product_smiles: Optional[str] = None,
    temperature: float = 298.15,
    current_user: User = Depends(get_current_verified_user)
):
    """
    Find transition state for catalytic reaction.

    Uses quantum TS search with governance acceleration.
    """
    try:
        # Use service layer

        ts_result = CatalystOptimizationService.find_transition_state(
            catalyst_smiles=catalyst_smiles,
            reactant_smiles=reactant_smiles,
            product_smiles=product_smiles,
            temperature=temperature
        )

        return {
            "application": "catalyst_optimization",
            "analysis_type": "transition_state",
            "results": ts_result
        }

    except Exception as e:
        logger.error(f"TS search failed: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"TS search failed: {str(e)}")


@router.post("/catalyst/screen")
async def screen_catalysts(
    request: CatalystScreeningRequest,
    current_user: User = Depends(get_current_verified_user)
):
    """
    Screen multiple catalysts for reaction.

    Returns ranked list by predicted activity.
    """
    try:
        # Use service layer

        results = CatalystOptimizationService.screen_catalysts(
            catalyst_smiles_list=request.catalyst_smiles_list,
            reactant_smiles=request.reactant_smiles,
            product_smiles=request.product_smiles,
            temperature=request.temperature,
            top_n=request.top_n
        )

        return {
            "application": "catalyst_optimization",
            "analysis_type": "screening",
            "catalysts_screened": len(request.catalyst_smiles_list),
            "top_catalysts": results
        }

    except Exception as e:
        logger.error(f"Catalyst screening failed: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"Screening failed: {str(e)}")


# ================================
# MATERIALS SCOUT ENDPOINTS
# ================================

@router.post("/materials/analyze")
async def analyze_material(
    request: MaterialAnalysisRequest,
    current_user: User = Depends(get_current_verified_user)
):
    """
    Analyze material for optoelectronic properties.

    Calculates band gap, absorption, LED color, etc.
    """
    try:
        # Get experiment
        experiment = ExperimentDB.get(request.experiment_id)
        if not experiment:
            raise HTTPException(status_code=404, detail="Experiment not found")

        if experiment.get('user_id') != current_user.id:
            raise HTTPException(status_code=403, detail="Not authorized")

        if experiment.get('status') != 'completed':
            raise HTTPException(status_code=400, detail="Experiment not completed")

        # Use service layer

        # Get results
        results = experiment.get('results', {})

        # Analyze material
        analysis = MaterialsAnalysisService.analyze_material(
            quantum_energy=results.get('energy'),
            homo_lumo_gap=results.get('homo_lumo_gap'),
            orbital_energies=results.get('orbital_energies', []),
            molecular_properties=results.get('analysis', {}),
            calculate_bandgap=request.calculate_bandgap,
            calculate_optical=request.calculate_optical,
            calculate_led_color=request.calculate_led_color,
            wavelength_range=request.wavelength_range
        )

        return {
            "experiment_id": request.experiment_id,
            "application": "materials_science",
            "results": analysis
        }

    except Exception as e:
        logger.error(f"Materials analysis failed: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"Analysis failed: {str(e)}")


@router.post("/materials/bandgap")
async def calculate_bandgap(
    material_smiles: str,
    temperature: float = 298.15,
    current_user: User = Depends(get_current_verified_user)
):
    """
    Calculate band gap for material.

    Includes temperature dependence.
    """
    try:
        # Use service layer

        bandgap_result = MaterialsAnalysisService.calculate_bandgap(
            material_smiles=material_smiles,
            temperature=temperature
        )

        return {
            "application": "materials_science",
            "analysis_type": "bandgap",
            "smiles": material_smiles,
            "results": bandgap_result
        }

    except Exception as e:
        logger.error(f"Bandgap calculation failed: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"Calculation failed: {str(e)}")


@router.post("/materials/optical-properties")
async def calculate_optical_properties(
    material_smiles: str,
    wavelength_range: Tuple[float, float] = (400, 700),
    current_user: User = Depends(get_current_verified_user)
):
    """
    Calculate optical absorption spectrum.

    Returns wavelength-dependent absorption.
    """
    try:
        # Use service layer

        optical_result = MaterialsAnalysisService.calculate_optical_properties(
            material_smiles=material_smiles,
            wavelength_range=wavelength_range
        )

        return {
            "application": "materials_science",
            "analysis_type": "optical_properties",
            "smiles": material_smiles,
            "results": optical_result
        }

    except Exception as e:
        logger.error(f"Optical calculation failed: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"Calculation failed: {str(e)}")


@router.post("/materials/led-color")
async def predict_led_color(
    material_smiles: str,
    current_user: User = Depends(get_current_verified_user)
):
    """
    Predict LED emission color.

    Returns RGB color and wavelength.
    """
    try:
        # Use service layer

        led_result = MaterialsAnalysisService.predict_led_color(material_smiles=material_smiles)

        return {
            "application": "materials_science",
            "analysis_type": "led_color",
            "smiles": material_smiles,
            "results": led_result
        }

    except Exception as e:
        logger.error(f"LED prediction failed: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"Prediction failed: {str(e)}")


@router.post("/materials/screen")
async def screen_materials(
    request: MaterialScreeningRequest,
    current_user: User = Depends(get_current_verified_user)
):
    """
    Screen multiple materials for target properties.

    Returns ranked list by target property.
    """
    try:
        # Use service layer

        results = MaterialsAnalysisService.screen_materials(
            material_smiles_list=request.material_smiles_list,
            target_bandgap=request.target_bandgap,
            target_property=request.target_property,
            top_n=request.top_n
        )

        return {
            "application": "materials_science",
            "analysis_type": "screening",
            "materials_screened": len(request.material_smiles_list),
            "top_materials": results
        }

    except Exception as e:
        logger.error(f"Materials screening failed: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"Screening failed: {str(e)}")


# ================================
# UTILITY ENDPOINTS
# ================================

@router.get("/domains")
async def get_available_domains():
    """
    Get list of available application domains.
    """
    return {
        "domains": [
            {
                "id": "drug_discovery",
                "name": "Drug Discovery",
                "description": "Pharmaceutical research - binding, ADME, druglikeness",
                "competitive_advantage": "800M-1.9B x speedup, <1 kcal/mol binding accuracy",
                "use_cases": ["lead optimization", "virtual screening", "ADME prediction"]
            },
            {
                "id": "alloy_design",
                "name": "Alloy Designer",
                "description": "Materials science - alloy composition, phase diagrams, mechanical properties",
                "competitive_advantage": "NEW alloy prediction, pressure-dependent phases",
                "use_cases": ["aerospace alloys", "automotive materials", "armor design"]
            },
            {
                "id": "catalyst_optimization",
                "name": "Catalyst Optimizer",
                "description": "Chemical catalysis - activation energy, selectivity, TOF",
                "competitive_advantage": "800M-1.9B x speedup vs DFT, quantum TS",
                "use_cases": ["industrial catalysis", "green chemistry", "reaction optimization"]
            },
            {
                "id": "materials_science",
                "name": "Materials Scout",
                "description": "Optoelectronic materials - band gap, optical properties, LED color",
                "competitive_advantage": "Quantum accuracy for semiconductors, fast screening",
                "use_cases": ["LED design", "solar cells", "photodetectors"]
            }
        ]
    }


@router.get("/capabilities/{domain}")
async def get_domain_capabilities(domain: str):
    """
    Get detailed capabilities for specific domain.
    """
    capabilities = {
        "drug_discovery": {
            "analyses": [
                "binding_affinity",
                "adme_properties",
                "druglikeness_score",
                "lipinski_violations",
                "metabolite_prediction",
                "ph_dependent_binding"
            ],
            "competitive_tools": ["SwissADME (free)", "Schrödinger ($$$)", "DataWarrior (free)"],
            "our_advantages": [
                "Quantum binding accuracy",
                "pH-dependent analysis",
                "Free tier available"
            ]
        },
        "alloy_design": {
            "analyses": [
                "formation_energy",
                "mechanical_properties",
                "phase_stability",
                "phase_diagrams",
                "composition_prediction"
            ],
            "competitive_tools": ["CALPHAD", "Thermo-Calc ($30k-100k/yr)", "Materials Project (free)"],
            "our_advantages": [
                "NEW alloy discovery",
                "Pressure-dependent phases",
                "Quantum first-principles"
            ]
        },
        "catalyst_optimization": {
            "analyses": [
                "activation_energy",
                "transition_state",
                "selectivity",
                "turnover_frequency",
                "reaction_mechanism"
            ],
            "competitive_tools": ["Gaussian ($$)", "VASP ($$)", "DFT codes (slow)"],
            "our_advantages": [
                "800M-1.9B x speedup",
                "Governance acceleration",
                "Real-time optimization"
            ]
        },
        "materials_science": {
            "analyses": [
                "band_gap",
                "optical_absorption",
                "led_color",
                "doping_effects",
                "carrier_concentration"
            ],
            "competitive_tools": ["Materials Project (free)", "VASP ($$)", "Quantum Espresso (free)"],
            "our_advantages": [
                "Fast screening",
                "LED color prediction",
                "Quantum accuracy"
            ]
        }
    }

    if domain not in capabilities:
        raise HTTPException(status_code=404, detail=f"Domain '{domain}' not found")

    return {
        "domain": domain,
        **capabilities[domain]
    }
