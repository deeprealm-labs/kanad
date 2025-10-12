"""
Configuration settings for the API server
"""

import os
from functools import lru_cache
from typing import List


class Settings:
    """Application settings."""

    VERSION: str = "0.1.0"

    # Server configuration
    HOST: str = os.getenv("API_HOST", "0.0.0.0")
    PORT: int = int(os.getenv("API_PORT", "8000"))
    DEBUG: bool = os.getenv("DEBUG", "true").lower() == "true"

    # CORS
    CORS_ORIGINS: List[str] = [
        "http://localhost:3000",
        "http://localhost:3001",
        "http://127.0.0.1:3000",
    ]

    # Database
    DATABASE_PATH: str = os.getenv("DATABASE_PATH", "kanad_experiments.db")

    # Cloud credentials
    IBM_API_TOKEN: str = os.getenv("IBM_API", "")
    IBM_CRN: str = os.getenv("IBM_CRN", "")
    BLUEQUBIT_TOKEN: str = os.getenv("BLUE_TOKEN", "")

    # Job processing
    MAX_CONCURRENT_JOBS: int = int(os.getenv("MAX_CONCURRENT_JOBS", "2"))
    JOB_TIMEOUT: int = int(os.getenv("JOB_TIMEOUT", "3600"))  # 1 hour

    # Default computation settings
    DEFAULT_BASIS: str = "sto-3g"
    DEFAULT_OPTIMIZER: str = "SLSQP"
    DEFAULT_MAX_ITERATIONS: int = 1000


@lru_cache()
def get_settings() -> Settings:
    """Get cached settings instance."""
    return Settings()
