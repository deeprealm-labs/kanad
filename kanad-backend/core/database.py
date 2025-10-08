"""
Database configuration and session management.

Handles PostgreSQL connections, session lifecycle, and dependency injection.
"""

from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker, Session
from typing import Generator
import os
import logging

logger = logging.getLogger(__name__)

# Database URL from environment
DATABASE_URL = os.getenv(
    "DATABASE_URL",
    "postgresql://kanad_user:kanad_password@localhost:5432/kanad_experiments"
)

# Create SQLAlchemy engine
engine = create_engine(
    DATABASE_URL,
    pool_pre_ping=True,  # Verify connections before using
    pool_size=10,  # Connection pool size
    max_overflow=20,  # Max overflow connections
    echo=False  # Set to True for SQL query logging
)

# Create session factory
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

# Base class for models (imported in db/models.py)
Base = declarative_base()


def get_db() -> Generator[Session, None, None]:
    """
    FastAPI dependency for database sessions.

    Usage:
        @app.get("/users/{user_id}")
        def get_user(user_id: str, db: Session = Depends(get_db)):
            user = db.query(User).filter(User.user_id == user_id).first()
            return user
    """
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()


def init_db():
    """
    Initialize database (create tables).

    Call this at application startup.
    """
    from db.models import Base
    Base.metadata.create_all(bind=engine)
    logger.info("Database tables created successfully")


def drop_db():
    """Drop all tables (use with caution!)."""
    from db.models import Base
    Base.metadata.drop_all(bind=engine)
    logger.warning("All database tables dropped")
