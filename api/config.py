"""
Configuration management for Kanad API.

Handles environment variables, database connections, and application settings.
"""

import os
from pathlib import Path
from typing import Optional
from pydantic_settings import BaseSettings
from pydantic import ConfigDict


class Settings(BaseSettings):
    """Application settings loaded from environment variables."""

    # Application
    APP_NAME: str = "Kanad Quantum Chemistry API"
    APP_VERSION: str = "1.0.0"
    DEBUG: bool = True

    # API Configuration
    API_V1_PREFIX: str = "/api/v1"

    # CORS Origins
    CORS_ORIGINS: list[str] = [
        "http://localhost:3000",
        "http://localhost:3001",
        "https://kanad.deeprealm.in"
    ]

    # Database
    DATABASE_URL: str = "sqlite:///./kanad.db"
    DATABASE_ECHO: bool = False  # Set to True for SQL logging

    # Job Queue
    MAX_CONCURRENT_JOBS: int = 2
    JOB_TIMEOUT_SECONDS: int = 3600  # 1 hour default timeout

    # Quantum Backends - Environment variable fallbacks
    IBM_QUANTUM_TOKEN: Optional[str] = None  # Fallback if not in DB
    IBM_QUANTUM_CRN: Optional[str] = None  # Fallback if not in DB
    BLUEQUBIT_API_KEY: Optional[str] = None  # Fallback if not in DB

    # Default Computation Settings
    DEFAULT_METHOD: str = "VQE"
    DEFAULT_ANSATZ: str = "ucc"
    DEFAULT_MAPPER: str = "jordan_wigner"
    DEFAULT_OPTIMIZER: str = "SLSQP"
    DEFAULT_BACKEND: str = "classical"
    DEFAULT_BASIS: str = "sto-3g"

    # File Storage
    RESULTS_DIR: Path = Path("./results")
    UPLOAD_DIR: Path = Path("./uploads")

    model_config = ConfigDict(
        env_file=".env",
        case_sensitive=True,
        extra="ignore"  # Ignore extra environment variables
    )


# Global settings instance
settings = Settings()

# Create necessary directories
settings.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
settings.UPLOAD_DIR.mkdir(parents=True, exist_ok=True)
