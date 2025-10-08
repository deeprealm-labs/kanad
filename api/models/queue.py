"""
Job queue database model.

Manages queued, scheduled, and running experiments.
"""

from sqlalchemy import Column, Integer, String, DateTime, ForeignKey
from sqlalchemy.sql import func
from sqlalchemy.orm import relationship

from api.database import Base


class QueueItem(Base):
    """Job queue database model."""

    __tablename__ = "queue"

    # Primary key
    id = Column(Integer, primary_key=True, index=True)

    # Reference to experiment
    experiment_id = Column(Integer, ForeignKey("experiments.id"), nullable=False, index=True)

    # Queue management
    status = Column(String(50), nullable=False, default="queued", index=True)
    # Status values: "queued", "scheduled", "running", "paused", "completed", "failed", "cancelled"

    priority = Column(Integer, nullable=False, default=0, index=True)
    # Higher priority = runs first (0 is lowest priority)

    # Scheduling
    scheduled_time = Column(DateTime, nullable=True)
    # If set, experiment will run at this time (null = run ASAP)

    # Timestamps
    created_at = Column(DateTime, default=func.now(), nullable=False)
    updated_at = Column(DateTime, default=func.now(), onupdate=func.now(), nullable=False)

    def __repr__(self):
        return f"<QueueItem(id={self.id}, experiment_id={self.experiment_id}, status={self.status}, priority={self.priority})>"

    def to_dict(self):
        """Convert to dictionary for API response."""
        return {
            "id": self.id,
            "experiment_id": self.experiment_id,
            "status": self.status,
            "priority": self.priority,
            "scheduled_time": self.scheduled_time.isoformat() if self.scheduled_time else None,
            "created_at": self.created_at.isoformat() if self.created_at else None,
            "updated_at": self.updated_at.isoformat() if self.updated_at else None,
        }
