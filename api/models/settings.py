"""
User settings database model.

Stores default computation settings and preferences.
"""

from sqlalchemy import Column, Integer, String, Boolean, JSON, DateTime
from sqlalchemy.sql import func

from api.database import Base


class UserSettings(Base):
    """User settings database model."""

    __tablename__ = "settings"

    # Primary key
    id = Column(Integer, primary_key=True, index=True)

    # For future multi-user support, default to single user for now
    user_id = Column(Integer, nullable=True, default=1, index=True)

    # Computation defaults
    method = Column(String(50), default="VQE", nullable=False)
    ansatz = Column(String(50), default="ucc", nullable=False)
    mapper = Column(String(50), default="jordan_wigner", nullable=False)
    optimizer = Column(String(50), default="SLSQP", nullable=False)
    backend = Column(String(50), default="classical", nullable=False)
    backend_name = Column(String(100), nullable=True)  # e.g., "ibm_torino"

    # Optimization toggles
    geometry_optimization = Column(Boolean, default=False)
    orbital_optimization = Column(Boolean, default=False)
    circuit_optimization = Column(Boolean, default=True)
    adaptive_vqe = Column(Boolean, default=False)

    # Advanced settings (JSON for flexibility)
    advanced_settings = Column(JSON, nullable=True)
    # JSON structure: {max_iterations: 1000, conv_threshold: 1e-6, ...}

    # Timestamps
    created_at = Column(DateTime, default=func.now(), nullable=False)
    updated_at = Column(DateTime, default=func.now(), onupdate=func.now(), nullable=False)

    def __repr__(self):
        return f"<UserSettings(id={self.id}, method={self.method}, backend={self.backend})>"

    def to_dict(self):
        """Convert to dictionary for API response."""
        return {
            "id": self.id,
            "user_id": self.user_id,
            "method": self.method,
            "ansatz": self.ansatz,
            "mapper": self.mapper,
            "optimizer": self.optimizer,
            "backend": self.backend,
            "backend_name": self.backend_name,
            "geometry_optimization": self.geometry_optimization,
            "orbital_optimization": self.orbital_optimization,
            "circuit_optimization": self.circuit_optimization,
            "adaptive_vqe": self.adaptive_vqe,
            "advanced_settings": self.advanced_settings,
            "created_at": self.created_at.isoformat() if self.created_at else None,
            "updated_at": self.updated_at.isoformat() if self.updated_at else None,
        }
