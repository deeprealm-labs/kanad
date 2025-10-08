"""
Experiment database model.

Stores quantum chemistry experiments with configuration, results, and convergence data.
"""

from sqlalchemy import Column, Integer, String, Text, DateTime, JSON, Float
from sqlalchemy.sql import func
from datetime import datetime

from api.database import Base


class Experiment(Base):
    """Experiment database model."""

    __tablename__ = "experiments"

    # Primary key
    id = Column(Integer, primary_key=True, index=True)

    # Experiment metadata
    name = Column(String(255), nullable=True)
    status = Column(String(50), nullable=False, default="queued", index=True)
    # Status values: "queued", "running", "completed", "failed", "cancelled"

    # Molecule data
    smiles = Column(String(500), nullable=True)
    molecule_data = Column(JSON, nullable=False)
    # JSON structure: {atoms: [...], smiles: "...", basis: "...", charge: 0, multiplicity: 1}

    # Computation configuration
    configuration = Column(JSON, nullable=False)
    # JSON structure: {method: "VQE", ansatz: "ucc", mapper: "jordan_wigner",
    #                  optimizer: "SLSQP", backend: "classical", ...}

    # Results
    energy = Column(Float, nullable=True)
    hf_energy = Column(Float, nullable=True)
    correlation_energy = Column(Float, nullable=True)
    results = Column(JSON, nullable=True)
    # JSON structure: {energy: -1.137, properties: {...}, analysis: {...}}

    # Convergence tracking
    convergence_data = Column(JSON, nullable=True)
    # JSON structure: [{iteration: 1, energy: -1.1}, {iteration: 2, energy: -1.13}, ...]

    # Error handling
    error_message = Column(Text, nullable=True)

    # Cloud job tracking for cancellation
    cloud_job_id = Column(String(255), nullable=True)  # IBM/BlueQubit job ID
    cloud_backend = Column(String(50), nullable=True)  # 'ibm' or 'bluequbit'

    # Timestamps
    created_at = Column(DateTime, default=func.now(), nullable=False)
    started_at = Column(DateTime, nullable=True)
    completed_at = Column(DateTime, nullable=True)
    cancelled_at = Column(DateTime, nullable=True)
    updated_at = Column(DateTime, default=func.now(), onupdate=func.now(), nullable=False)

    def __repr__(self):
        return f"<Experiment(id={self.id}, name={self.name}, status={self.status})>"

    def to_dict(self):
        """Convert to dictionary for API response."""
        return {
            "id": self.id,
            "name": self.name,
            "status": self.status,
            "smiles": self.smiles,
            "molecule_data": self.molecule_data,
            "configuration": self.configuration,
            "energy": self.energy,
            "hf_energy": self.hf_energy,
            "correlation_energy": self.correlation_energy,
            "results": self.results,
            "convergence_data": self.convergence_data,
            "error_message": self.error_message,
            "cloud_job_id": self.cloud_job_id,
            "cloud_backend": self.cloud_backend,
            "created_at": self.created_at.isoformat() if self.created_at else None,
            "started_at": self.started_at.isoformat() if self.started_at else None,
            "completed_at": self.completed_at.isoformat() if self.completed_at else None,
            "cancelled_at": self.cancelled_at.isoformat() if self.cancelled_at else None,
            "updated_at": self.updated_at.isoformat() if self.updated_at else None,
        }
