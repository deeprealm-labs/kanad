"""
Database connection and session management.

Uses SQLAlchemy with SQLite for development (easy setup, no external dependencies).
Can be switched to PostgreSQL for production by changing DATABASE_URL.
"""

from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker, Session
from typing import Generator

from api.config import settings

# Create SQLAlchemy engine
engine = create_engine(
    settings.DATABASE_URL,
    connect_args={"check_same_thread": False} if "sqlite" in settings.DATABASE_URL else {},
    echo=settings.DATABASE_ECHO
)

# Session factory
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

# Base class for models
Base = declarative_base()


def get_db() -> Generator[Session, None, None]:
    """
    Dependency for getting database sessions.

    Usage in FastAPI endpoints:
        @app.get("/experiments")
        def get_experiments(db: Session = Depends(get_db)):
            ...
    """
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()


def init_db():
    """Initialize database tables."""
    # Import all models to register them with Base
    import api.models.experiment
    import api.models.queue
    import api.models.settings

    # Create all tables
    Base.metadata.create_all(bind=engine)
