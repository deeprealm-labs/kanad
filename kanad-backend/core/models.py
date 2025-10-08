"""
Pydantic models for API requests and responses.

These models define the structure of data exchanged between the frontend
and backend, with automatic validation.
"""

from pydantic import BaseModel, Field, validator
from typing import List, Optional, Dict, Any, Union
from datetime import datetime
from enum import Enum
import re


# ===== Enums =====

class BasisSet(str, Enum):
    """Available basis sets."""
    STO_3G = "sto-3g"
    STO_6G = "sto-6g"
    E_631G = "6-31g"
    E_631Gd = "6-31g*"
    E_631Gdp = "6-31g**"
    CC_PVDZ = "cc-pvdz"
    CC_PVTZ = "cc-pvtz"


class ComputationMethod(str, Enum):
    """Supported quantum chemistry methods."""
    HF = "HF"
    VQE = "VQE"
    MP2 = "MP2"
    FCI = "FCI"
    SQD = "SQD"  # Subspace Quantum Diagonalization
    EXCITED_STATES = "EXCITED_STATES"  # Excited states calculation


class AnsatzType(str, Enum):
    """VQE ansatz types."""
    UCC = "ucc"
    HARDWARE_EFFICIENT = "hardware_efficient"
    GOVERNANCE = "governance"
    TWO_LOCAL = "two_local"  # Two-local ansatz
    UCC_CORRECT_DOUBLE = "ucc_correct_double"  # UCC with corrected doubles


class MapperType(str, Enum):
    """Fermion-to-qubit mapping strategies."""
    JORDAN_WIGNER = "jordan_wigner"
    BRAVYI_KITAEV = "bravyi_kitaev"
    PARITY = "parity"
    HYBRID_ORBITAL = "hybrid_orbital"  # Hybrid orbital mapper for metallurgy


class BackendType(str, Enum):
    """Quantum backend types."""
    CLASSICAL = "classical"
    IBM_QUANTUM = "ibm_quantum"
    BLUEQUBIT = "bluequbit"


class JobStatus(str, Enum):
    """Job execution status."""
    QUEUED = "queued"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


class MoleculeCreationMethod(str, Enum):
    """Methods for creating molecules."""
    ATOMS = "atoms"
    SMILES = "smiles"
    LIBRARY = "library"
    XYZ = "xyz"


# ===== Molecule Models =====

class Atom(BaseModel):
    """Single atom representation."""
    element: str = Field(..., description="Element symbol (e.g., 'H', 'O', 'Fe')")
    position: List[float] = Field(..., description="3D position [x, y, z] in Angstroms")

    @validator('position')
    def validate_position(cls, v):
        if len(v) != 3:
            raise ValueError('Position must have exactly 3 coordinates')
        return v


class MoleculeCreate(BaseModel):
    """Request to create a molecule."""
    method: MoleculeCreationMethod
    data: Dict[str, Any] = Field(..., description="Method-specific data")

    class Config:
        schema_extra = {
            "example": {
                "method": "atoms",
                "data": {
                    "atoms": [
                        {"element": "H", "position": [0, 0, 0]},
                        {"element": "O", "position": [0, 0, 0.96]},
                        {"element": "H", "position": [0, 0.76, 0.59]}
                    ],
                    "basis": "sto-3g",
                    "charge": 0,
                    "multiplicity": 1
                }
            }
        }


class LewisStructure(BaseModel):
    """Lewis structure representation."""
    bonds: List[List[int]] = Field(..., description="Pairs of atom indices forming bonds")
    lone_pairs: Dict[str, int] = Field(default_factory=dict, description="Lone pairs per atom")
    formal_charges: Dict[str, int] = Field(default_factory=dict, description="Formal charges")


class MoleculeResponse(BaseModel):
    """Response after creating a molecule."""
    molecule_id: str
    name: Optional[str] = None
    formula: str
    geometry: Dict[str, Any]
    lewis_structure: Optional[LewisStructure] = None
    n_electrons: int
    n_orbitals: int
    n_qubits: int
    preview: Optional[str] = Field(None, description="Base64-encoded SVG image")
    created_at: datetime = Field(default_factory=datetime.utcnow)


class SMILESCreate(BaseModel):
    """Create molecule from SMILES."""
    smiles: str = Field(..., description="SMILES string")
    basis: BasisSet = Field(default=BasisSet.STO_3G)
    optimize_geometry: bool = Field(default=False)


# ===== Simulation Configuration =====

class BackendConfig(BaseModel):
    """Backend configuration."""
    type: BackendType = Field(default=BackendType.CLASSICAL)
    use_user_credentials: bool = Field(default=False)
    backend_name: Optional[str] = Field(None, description="Specific backend (e.g., 'ibm_torino')")


class AnalysisConfig(BaseModel):
    """Analysis options."""
    energy_decomposition: bool = Field(default=True)
    bond_analysis: bool = Field(default=True)
    dipole_moment: bool = Field(default=True)
    polarizability: bool = Field(default=False)
    thermochemistry: bool = Field(default=False)
    spectroscopy: bool = Field(default=False)
    vibrational: bool = Field(default=False)


class OptimizationConfig(BaseModel):
    """Optimization settings."""
    geometry: bool = Field(default=False)
    orbitals: bool = Field(default=False)
    circuit: bool = Field(default=True)
    adaptive: bool = Field(default=False)


class AdvancedSettings(BaseModel):
    """Advanced computation settings (Pro mode)."""
    active_space: Optional[Dict[str, int]] = Field(None, description="{'n_electrons': 4, 'n_orbitals': 6}")
    frozen_core: bool = Field(default=True)
    symmetry: str = Field(default="auto")
    initial_state: str = Field(default="hf")


class SimulationConfig(BaseModel):
    """Simulation configuration request."""
    molecule_id: str
    method: ComputationMethod = Field(default=ComputationMethod.VQE)

    # VQE-specific
    ansatz: Optional[AnsatzType] = Field(default=AnsatzType.HARDWARE_EFFICIENT)
    mapper: Optional[MapperType] = Field(default=MapperType.JORDAN_WIGNER)
    optimizer: str = Field(default="SLSQP")
    max_iterations: int = Field(default=1000, ge=1, le=10000)
    convergence_threshold: float = Field(default=1e-6, gt=0)

    # Backend
    backend: BackendConfig = Field(default_factory=BackendConfig)

    # Analysis and optimization
    analysis: AnalysisConfig = Field(default_factory=AnalysisConfig)
    optimization: OptimizationConfig = Field(default_factory=OptimizationConfig)
    advanced: AdvancedSettings = Field(default_factory=AdvancedSettings)

    @validator('molecule_id')
    def validate_molecule_id(cls, v):
        # UUID format validation
        uuid_pattern = r'^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$'
        if not re.match(uuid_pattern, v):
            raise ValueError('Invalid molecule_id format')
        return v


class CircuitPreview(BaseModel):
    """Circuit preview information."""
    n_qubits: int
    n_parameters: int
    circuit_depth: int
    estimated_time: str
    estimated_cost: str


class SimulationPreview(BaseModel):
    """Preview of simulation before execution."""
    simulation_id: str
    status: str = "configured"
    preview: CircuitPreview
    circuit_preview: Optional[str] = Field(None, description="Base64-encoded circuit diagram")
    terms_and_conditions: Dict[str, str]


class UserCredentials(BaseModel):
    """Encrypted user credentials for cloud providers."""
    ibm_api: Optional[str] = None
    ibm_crn: Optional[str] = None
    blue_token: Optional[str] = None


class SimulationSubmit(BaseModel):
    """Submit simulation for execution."""
    accepted_terms: bool = Field(..., description="Must accept T&C")
    user_credentials: Optional[UserCredentials] = None


# ===== Job Models =====

class JobListItem(BaseModel):
    """Single job in job list."""
    job_id: str
    molecule_name: str
    method: ComputationMethod
    status: JobStatus
    progress: int = Field(..., ge=0, le=100)
    created_at: datetime
    started_at: Optional[datetime] = None
    backend: str


class JobStatusResponse(BaseModel):
    """Real-time job status."""
    job_id: str
    status: JobStatus
    progress: int = Field(..., ge=0, le=100)
    current_iteration: Optional[int] = None
    max_iterations: Optional[int] = None
    current_energy: Optional[float] = None
    best_energy: Optional[float] = None
    message: str


# ===== Results Models =====

class EnergyDecomposition(BaseModel):
    """Energy component breakdown."""
    kinetic: float
    nuclear_attraction: float
    electron_repulsion: float


class BondInfo(BaseModel):
    """Single bond information."""
    atoms: List[int]
    order: float
    length: float


class BondAnalysis(BaseModel):
    """Bonding analysis results."""
    bonds: List[BondInfo]
    homo_lumo_gap: float


class DipoleMoment(BaseModel):
    """Dipole moment."""
    magnitude: float = Field(..., description="Dipole magnitude in Debye")
    direction: List[float]


class Thermochemistry(BaseModel):
    """Thermochemistry results."""
    enthalpy: float
    entropy: float
    gibbs_free_energy: float


class AnalysisResults(BaseModel):
    """Complete analysis results."""
    energy_decomposition: Optional[EnergyDecomposition] = None
    bond_analysis: Optional[BondAnalysis] = None
    dipole_moment: Optional[DipoleMoment] = None
    thermochemistry: Optional[Thermochemistry] = None


class ConvergencePoint(BaseModel):
    """Single point in convergence history."""
    iteration: int
    energy: float


class ComputationResults(BaseModel):
    """Core computation results."""
    method: ComputationMethod
    energy: float
    hf_energy: float
    correlation_energy: float
    n_iterations: int
    converged: bool
    convergence_history: List[ConvergencePoint]


class LLMReport(BaseModel):
    """AI-generated report."""
    summary: str
    key_findings: List[str]
    interpretation: str
    recommendations: List[str]


class JobResults(BaseModel):
    """Complete job results."""
    job_id: str
    status: JobStatus
    molecule: Dict[str, str]
    results: ComputationResults
    analysis: Optional[AnalysisResults] = None
    llm_report: Optional[LLMReport] = None


# ===== Cloud Provider Models =====

class CloudCredentials(BaseModel):
    """Cloud provider credentials."""
    provider: str = Field(..., description="'ibm_quantum' or 'bluequbit'")
    credentials: Dict[str, str]


class BackendInfo(BaseModel):
    """Cloud backend information."""
    name: str
    qubits: int
    queue_depth: Optional[int] = None
    online: bool
    estimated_wait: Optional[str] = None


class CloudProviderInfo(BaseModel):
    """Cloud provider status."""
    name: str
    status: str
    backends: List[BackendInfo]


# ===== User Models =====

class UserRegister(BaseModel):
    """User registration."""
    email: str
    password: str
    name: Optional[str] = None
    institution: Optional[str] = None
    field: Optional[str] = None


class UserLogin(BaseModel):
    """User login."""
    email: str
    password: str


class Token(BaseModel):
    """JWT token response."""
    access_token: str
    token_type: str = "bearer"


class UserProfile(BaseModel):
    """User profile information."""
    user_id: str
    email: str
    name: Optional[str]
    institution: Optional[str]
    field: Optional[str]
    stats: Dict[str, int]


# ===== Batch/Schedule Models =====

class ExperimentConfig(BaseModel):
    """Single experiment configuration."""
    molecule_id: str
    method: ComputationMethod
    basis: BasisSet
    ansatz: Optional[AnsatzType] = None
    mapper: Optional[MapperType] = None


class ScheduleCreate(BaseModel):
    """Create batch schedule."""
    name: str
    experiments: List[ExperimentConfig]
    execution_mode: str = Field(default="sequential", description="'sequential' or 'parallel'")
    priority: str = Field(default="normal", description="'low', 'normal', or 'high'")


class ScheduleResponse(BaseModel):
    """Schedule creation response."""
    schedule_id: str
    total_experiments: int
    estimated_total_time: str
    job_ids: List[str]


# ===== Settings Models =====

class DefaultSettings(BaseModel):
    """User default settings."""
    computation: Dict[str, str]
    optimization: Dict[str, bool]
    analysis: Dict[str, Any]
    cloud: Dict[str, str]


# ===== Library Models =====

class LibraryMolecule(BaseModel):
    """Molecule in library."""
    id: str
    name: str
    formula: str
    smiles: Optional[str] = None


class LibraryCategory(BaseModel):
    """Category in molecule library."""
    name: str
    molecules: List[LibraryMolecule]


class MoleculeLibrary(BaseModel):
    """Complete molecule library."""
    categories: List[LibraryCategory]
