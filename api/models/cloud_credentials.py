"""
Cloud credentials database model.

Stores encrypted cloud provider credentials for IBM Quantum and BlueQubit.
"""

from sqlalchemy import Column, Integer, String, DateTime, Text
from sqlalchemy.sql import func

from api.database import Base


class CloudCredentials(Base):
    """Cloud credentials database model."""

    __tablename__ = "cloud_credentials"

    # Primary key
    id = Column(Integer, primary_key=True, index=True)

    # Provider type
    provider = Column(String(50), nullable=False, index=True)  # "ibm" or "bluequbit"

    # For future multi-user support
    user_id = Column(Integer, nullable=True, default=1, index=True)

    # Credentials (stored as encrypted text in production)
    # For IBM: stores CRN and API key
    # For BlueQubit: stores API token
    credentials = Column(Text, nullable=False)

    # Timestamps
    created_at = Column(DateTime, default=func.now(), nullable=False)
    updated_at = Column(DateTime, default=func.now(), onupdate=func.now(), nullable=False)

    def __repr__(self):
        return f"<CloudCredentials(id={self.id}, provider={self.provider})>"

    def to_dict(self):
        """Convert to dictionary for API response."""
        return {
            "id": self.id,
            "provider": self.provider,
            "user_id": self.user_id,
            "has_credentials": bool(self.credentials),
            "created_at": self.created_at.isoformat() if self.created_at else None,
            "updated_at": self.updated_at.isoformat() if self.updated_at else None,
        }
