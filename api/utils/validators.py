"""
Validation utilities for API requests.
"""

from typing import Optional
from pydantic import BaseModel, Field, validator


class MoleculeData(BaseModel):
    """Molecule data validation model."""

    smiles: Optional[str] = None
    atoms: Optional[list] = None
    basis: str = Field(default="sto-3g", description="Basis set")
    charge: int = Field(default=0, description="Molecular charge")
    multiplicity: int = Field(default=1, description="Spin multiplicity (2S+1)")

    @validator('basis')
    def validate_basis(cls, v):
        """Validate basis set."""
        valid_bases = ['sto-3g', '6-31g', '6-31g*', '6-31g**', 'cc-pvdz', 'cc-pvtz']
        if v.lower() not in valid_bases:
            raise ValueError(f"Invalid basis set. Valid options: {', '.join(valid_bases)}")
        return v.lower()

    @validator('multiplicity')
    def validate_multiplicity(cls, v):
        """Validate multiplicity."""
        if v < 1:
            raise ValueError("Multiplicity must be >= 1")
        return v

    class Config:
        schema_extra = {
            "example": {
                "smiles": "CCO",
                "basis": "sto-3g",
                "charge": 0,
                "multiplicity": 1
            }
        }


class ExperimentConfig(BaseModel):
    """Experiment configuration validation model."""

    method: str = Field(default="VQE", description="Computation method")
    ansatz: str = Field(default="ucc", description="Ansatz type")
    mapper: str = Field(default="jordan_wigner", description="Qubit mapper")
    optimizer: str = Field(default="SLSQP", description="Classical optimizer")
    backend: str = Field(default="classical", description="Backend type")
    backend_name: Optional[str] = Field(None, description="Specific backend name")
    max_iterations: int = Field(default=1000, description="Max iterations")
    conv_threshold: float = Field(default=1e-6, description="Convergence threshold")

    @validator('method')
    def validate_method(cls, v):
        """Validate computation method."""
        valid_methods = ['VQE', 'HF', 'SQD', 'QPE']
        if v.upper() not in valid_methods:
            raise ValueError(f"Invalid method. Valid options: {', '.join(valid_methods)}")
        return v.upper()

    @validator('ansatz')
    def validate_ansatz(cls, v):
        """Validate ansatz type."""
        valid_ansatze = ['ucc', 'hardware_efficient', 'governance']
        if v.lower() not in valid_ansatze:
            raise ValueError(f"Invalid ansatz. Valid options: {', '.join(valid_ansatze)}")
        return v.lower()

    @validator('mapper')
    def validate_mapper(cls, v):
        """Validate mapper type."""
        valid_mappers = ['jordan_wigner', 'bravyi_kitaev', 'hybrid_orbital']
        if v.lower() not in valid_mappers:
            raise ValueError(f"Invalid mapper. Valid options: {', '.join(valid_mappers)}")
        return v.lower()

    @validator('optimizer')
    def validate_optimizer(cls, v):
        """Validate optimizer."""
        valid_optimizers = ['SLSQP', 'COBYLA', 'L-BFGS-B', 'ADAM', 'POWELL']
        if v.upper() not in valid_optimizers:
            raise ValueError(f"Invalid optimizer. Valid options: {', '.join(valid_optimizers)}")
        return v.upper()

    @validator('backend')
    def validate_backend(cls, v):
        """Validate backend type."""
        valid_backends = ['classical', 'ibm_quantum', 'bluequbit']
        if v.lower() not in valid_backends:
            raise ValueError(f"Invalid backend. Valid options: {', '.join(valid_backends)}")
        return v.lower()

    class Config:
        schema_extra = {
            "example": {
                "method": "VQE",
                "ansatz": "ucc",
                "mapper": "jordan_wigner",
                "optimizer": "SLSQP",
                "backend": "classical",
                "max_iterations": 1000,
                "conv_threshold": 1e-6
            }
        }


class ExperimentCreate(BaseModel):
    """Request model for creating an experiment."""

    name: Optional[str] = Field(None, description="Experiment name")
    molecule: MoleculeData
    configuration: ExperimentConfig
    execute_immediately: bool = Field(default=False, description="Execute now or add to queue")

    class Config:
        schema_extra = {
            "example": {
                "name": "Ethanol VQE Calculation",
                "molecule": {
                    "smiles": "CCO",
                    "basis": "sto-3g",
                    "charge": 0,
                    "multiplicity": 1
                },
                "configuration": {
                    "method": "VQE",
                    "ansatz": "ucc",
                    "mapper": "jordan_wigner",
                    "optimizer": "SLSQP",
                    "backend": "classical"
                },
                "execute_immediately": False
            }
        }


class QueueItemCreate(BaseModel):
    """Request model for adding item to queue."""

    experiment_id: int
    priority: int = Field(default=0, description="Queue priority")
    scheduled_time: Optional[str] = Field(None, description="ISO 8601 datetime")


class QueueItemUpdate(BaseModel):
    """Request model for updating queue item."""

    status: Optional[str] = None
    priority: Optional[int] = None
    scheduled_time: Optional[str] = None


class SettingsUpdate(BaseModel):
    """Request model for updating settings."""

    method: Optional[str] = None
    ansatz: Optional[str] = None
    mapper: Optional[str] = None
    optimizer: Optional[str] = None
    backend: Optional[str] = None
    backend_name: Optional[str] = None
    geometry_optimization: Optional[bool] = None
    orbital_optimization: Optional[bool] = None
    circuit_optimization: Optional[bool] = None
    adaptive_vqe: Optional[bool] = None


class SMILESValidation(BaseModel):
    """Request model for SMILES validation."""

    smiles: str = Field(..., description="SMILES string to validate")

    class Config:
        schema_extra = {
            "example": {
                "smiles": "CCO"
            }
        }
