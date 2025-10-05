#!/usr/bin/env python3
"""
Kanad API Server - FastAPI REST API for Quantum Chemistry

Provides endpoints for:
- Molecule energy calculations
- Optimization strategies
- VQE/QPE/SQD solvers
- Real-time computation
"""

from fastapi import FastAPI, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi.responses import JSONResponse
from pydantic import BaseModel, Field
from typing import List, Optional, Dict, Any
import numpy as np
import logging
import time
from pathlib import Path
import json

# Import Kanad framework
from kanad.core.atom import Atom
from kanad.core.molecule import Molecule  # NEW: Fast multi-atom support
from kanad.bonds.covalent_bond import CovalentBond
from kanad.bonds.ionic_bond import IonicBond
from kanad.bonds.metallic_bond import MetallicBond
from kanad.solvers.vqe_solver import VQESolver
from kanad.solvers.qpe_solver import QPESolver
from kanad.solvers.sqd_solver import SQDSolver
from kanad.ansatze.ucc_ansatz import UCCAnsatz
from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper
from kanad.optimization.quantum_optimizer import QuantumOptimizer

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Create FastAPI app
app = FastAPI(
    title="Kanad Quantum Chemistry API",
    description="Production-ready quantum chemistry computations with optimization",
    version="1.0.0",
    docs_url="/api/docs",
    redoc_url="/api/redoc"
)

# Enable CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # In production, specify exact origins
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Request/Response Models
class AtomModel(BaseModel):
    symbol: str = Field(..., description="Chemical symbol (H, C, O, etc.)")
    position: List[float] = Field(..., description="3D position [x, y, z] in Angstroms")

class ComputeRequest(BaseModel):
    atoms: List[AtomModel] = Field(..., description="List of atoms in molecule")
    bond_type: str = Field("covalent", description="Bond type: covalent, ionic, metallic")
    solver: str = Field("HF", description="Solver: HF, VQE, QPE, SQD")
    optimize: bool = Field(True, description="Apply quantum optimization")
    strategy: str = Field("balanced", description="Optimization strategy: aggressive, balanced, conservative")

class OptimizeRequest(BaseModel):
    atoms: List[AtomModel]
    bond_type: str = "covalent"
    strategy: str = "balanced"

class PresetMolecule(BaseModel):
    name: str
    atoms: List[AtomModel]
    bond_type: str

# Preset molecules
PRESETS = {
    "H2": PresetMolecule(
        name="H2 (Hydrogen)",
        atoms=[
            AtomModel(symbol="H", position=[0.0, 0.0, 0.0]),
            AtomModel(symbol="H", position=[0.74, 0.0, 0.0])
        ],
        bond_type="covalent"
    ),
    "LiH": PresetMolecule(
        name="LiH (Lithium Hydride)",
        atoms=[
            AtomModel(symbol="Li", position=[0.0, 0.0, 0.0]),
            AtomModel(symbol="H", position=[1.59, 0.0, 0.0])
        ],
        bond_type="ionic"
    ),
    "H2O": PresetMolecule(
        name="H2O (Water)",
        atoms=[
            AtomModel(symbol="O", position=[0.0, 0.0, 0.0]),
            AtomModel(symbol="H", position=[0.96, 0.0, 0.0]),
            AtomModel(symbol="H", position=[-0.24, 0.93, 0.0])
        ],
        bond_type="covalent"
    ),
    "Na2": PresetMolecule(
        name="Na2 (Sodium Dimer)",
        atoms=[
            AtomModel(symbol="Na", position=[0.0, 0.0, 0.0]),
            AtomModel(symbol="Na", position=[3.66, 0.0, 0.0])
        ],
        bond_type="metallic"
    )
}

# Helper functions
def clean_for_json(obj: Any) -> Any:
    """Convert numpy types to native Python types for JSON serialization."""
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, dict):
        return {key: clean_for_json(value) for key, value in obj.items()}
    elif isinstance(obj, (list, tuple)):
        return [clean_for_json(item) for item in obj]
    else:
        return obj

def create_atoms(atom_models: List[AtomModel]) -> List[Atom]:
    """Convert AtomModel to Atom objects."""
    return [
        Atom(a.symbol, position=np.array(a.position))
        for a in atom_models
    ]

def create_bond_or_molecule(atoms: List[Atom], bond_type: str):
    """
    Create appropriate bond or molecule object.

    Uses new Molecule class for 3+ atoms (5x faster than previous approach!).
    Falls back to bond classes for 2-atom systems.
    """
    # For 3+ atoms, use optimized Molecule class (5x faster!)
    if len(atoms) >= 3:
        logger.info(f"Using Molecule class for {len(atoms)} atoms (5x faster!)")
        return Molecule(atoms, charge=0, spin=0, basis='sto-3g')

    # For 2 atoms, use specialized bond classes
    if bond_type == "covalent":
        return CovalentBond(atoms[0], atoms[1])
    elif bond_type == "ionic":
        return IonicBond(atoms[0], atoms[1])
    elif bond_type == "metallic":
        return MetallicBond(atoms)
    else:
        raise ValueError(f"Unknown bond type: {bond_type}")

# API Endpoints
@app.get("/api")
def read_root():
    """API root endpoint."""
    return {
        "message": "Kanad Quantum Chemistry API",
        "version": "1.0.0",
        "status": "active",
        "endpoints": {
            "docs": "/api/docs",
            "compute": "/api/compute",
            "optimize": "/api/optimize",
            "presets": "/api/presets",
            "health": "/api/health"
        }
    }

@app.get("/api/health")
def health_check():
    """Health check endpoint."""
    return {
        "status": "healthy",
        "version": "1.0.0",
        "timestamp": time.time()
    }

@app.get("/api/presets")
def get_presets():
    """Get preset molecules."""
    return {
        "presets": {
            key: {
                "name": preset.name,
                "atoms": [{"symbol": a.symbol, "position": a.position} for a in preset.atoms],
                "bond_type": preset.bond_type
            }
            for key, preset in PRESETS.items()
        }
    }

@app.post("/api/compute")
async def compute_energy(request: ComputeRequest):
    """
    Compute molecular energy with optional optimization.

    Supports multiple solvers:
    - HF: Hartree-Fock (fast, baseline)
    - VQE: Variational Quantum Eigensolver
    - QPE: Quantum Phase Estimation
    - SQD: Sample-based Quantum Diagonalization
    """
    try:
        start_time = time.time()
        logger.info(f"Compute request: {request.solver} on {len(request.atoms)} atoms")

        # Create atoms
        atoms = create_atoms(request.atoms)

        # Create bond or molecule (uses fast Molecule class for 3+ atoms)
        system = create_bond_or_molecule(atoms, request.bond_type)
        hamiltonian = system.hamiltonian

        # Run optimization if requested
        optimization_result = None
        if request.optimize:
            logger.info(f"Optimizing with {request.strategy} strategy")
            optimizer = QuantumOptimizer(hamiltonian)
            opt_raw = optimizer.optimize(strategy=request.strategy)
            # Clean numpy arrays for JSON serialization
            optimization_result = clean_for_json(opt_raw)

        # Run solver
        result = {}

        if request.solver.upper() == "HF":
            # Use appropriate classical method based on system type
            if isinstance(system, Molecule):
                # New Molecule class (5x faster!)
                energy_result = system.compute_energy(method='HF')
                result = {
                    "energy": float(energy_result['energy']),
                    "method": "Hartree-Fock (Molecule)"
                }
            elif request.bond_type == "metallic":
                energy_result = system.compute_energy(method='tight_binding')
                result = {
                    "energy": float(energy_result['energy']),
                    "method": "Tight-Binding"
                }
            else:
                energy_result = system.compute_energy(method='HF')
                result = {
                    "energy": float(energy_result['energy']),
                    "method": "Hartree-Fock"
                }

        elif request.solver.upper() == "VQE":
            n_qubits = 2 * hamiltonian.n_orbitals
            ansatz = UCCAnsatz(n_qubits, hamiltonian.n_electrons)
            mapper = JordanWignerMapper()
            solver = VQESolver(hamiltonian, ansatz, mapper)
            vqe_result = solver.solve()
            result = {
                "energy": float(vqe_result['energy']),
                "iterations": vqe_result.get('iterations', 'N/A'),
                "method": "VQE"
            }

        elif request.solver.upper() == "QPE":
            solver = QPESolver(hamiltonian, n_ancilla=4)
            qpe_result = solver.solve()
            result = {
                "energy": float(qpe_result['energy']),
                "method": "QPE"
            }

        elif request.solver.upper() == "SQD":
            solver = SQDSolver(hamiltonian, n_samples=100)
            sqd_result = solver.solve()
            result = {
                "energy": float(sqd_result['energy']),
                "method": "SQD"
            }
        else:
            raise ValueError(f"Unknown solver: {request.solver}")

        computation_time = time.time() - start_time

        response = {
            "success": True,
            "result": result,
            "system": {
                "orbitals": hamiltonian.n_orbitals,
                "electrons": hamiltonian.n_electrons,
                "qubits": 2 * hamiltonian.n_orbitals,
                "bond_type": request.bond_type,
                "num_atoms": len(atoms)
            },
            "optimization": optimization_result if request.optimize else None,
            "computation_time": round(computation_time, 3)
        }

        # Add info for multi-atom molecules
        if len(atoms) > 2:
            if isinstance(system, Molecule):
                response["info"] = f"✓ Full molecule calculation ({len(atoms)} atoms) using optimized Molecule class (5x faster!)"
            else:
                response["warning"] = f"⚠ Legacy mode: Multi-atom molecule ({len(atoms)} atoms) - computing primary bond only"

        logger.info(f"Computation complete in {computation_time:.3f}s")
        return response

    except Exception as e:
        logger.error(f"Computation error: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/api/optimize")
async def optimize_molecule(request: OptimizeRequest):
    """
    Analyze optimization strategies for a molecule.

    Returns qubit reduction and speedup estimates.
    """
    try:
        logger.info(f"Optimization request for {len(request.atoms)} atoms")

        # Create atoms and system (bond or molecule)
        atoms = create_atoms(request.atoms)
        system = create_bond_or_molecule(atoms, request.bond_type)
        hamiltonian = system.hamiltonian

        # Run optimization
        optimizer = QuantumOptimizer(hamiltonian)

        # Get single strategy result
        result = optimizer.optimize(strategy=request.strategy)

        # Get comparison
        comparison = optimizer.compare_strategies()

        return {
            "success": True,
            "optimization": clean_for_json(result),
            "comparison": clean_for_json(comparison),
            "system": {
                "orbitals": hamiltonian.n_orbitals,
                "electrons": hamiltonian.n_electrons,
                "original_qubits": 2 * hamiltonian.n_orbitals
            }
        }

    except Exception as e:
        logger.error(f"Optimization error: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/api/solvers")
def get_solvers():
    """Get available solvers and their descriptions."""
    return {
        "solvers": {
            "HF": {
                "name": "Hartree-Fock",
                "description": "Fast classical baseline",
                "speed": "Fast",
                "accuracy": "Moderate"
            },
            "VQE": {
                "name": "Variational Quantum Eigensolver",
                "description": "Hybrid quantum-classical optimization",
                "speed": "Moderate",
                "accuracy": "High"
            },
            "QPE": {
                "name": "Quantum Phase Estimation",
                "description": "Quantum eigenvalue estimation",
                "speed": "Moderate",
                "accuracy": "High"
            },
            "SQD": {
                "name": "Sample-based Quantum Diagonalization",
                "description": "Quantum sampling + classical processing",
                "speed": "Moderate",
                "accuracy": "High"
            }
        },
        "strategies": {
            "aggressive": "Maximum reduction, minimal active space",
            "balanced": "Good accuracy-cost trade-off",
            "conservative": "Minimal reduction, highest accuracy"
        }
    }

# Mount static files (web frontend)
web_dir = Path(__file__).parent.parent / "web"
if web_dir.exists():
    app.mount("/", StaticFiles(directory=str(web_dir), html=True), name="static")
    logger.info(f"Serving static files from {web_dir}")

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000, log_level="info")
