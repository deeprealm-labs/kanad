"""
Application configuration management.

Loads environment variables and provides configuration constants.
"""

from pydantic_settings import BaseSettings
from typing import Optional
import os


class Settings(BaseSettings):
    """Application settings loaded from environment variables."""

    # Application
    APP_NAME: str = "Kanad Quantum Chemistry API"
    APP_VERSION: str = "1.0.0"
    DEBUG: bool = False
    ENVIRONMENT: str = "production"  # 'development', 'staging', 'production'

    # API
    API_PREFIX: str = "/api"
    ALLOWED_ORIGINS: str = "https://your-app.vercel.app,http://localhost:3000"

    # Database
    DATABASE_URL: str = "postgresql://kanad_user:kanad_password@localhost:5432/kanad_experiments"

    # Redis & Celery
    REDIS_URL: str = "redis://localhost:6379/0"
    CELERY_BROKER_URL: str = "redis://localhost:6379/0"
    CELERY_RESULT_BACKEND: str = "redis://localhost:6379/0"

    # Security
    JWT_SECRET_KEY: str = os.getenv("JWT_SECRET_KEY", "CHANGE_THIS_IN_PRODUCTION")
    JWT_ALGORITHM: str = "HS256"
    JWT_ACCESS_TOKEN_EXPIRE_MINUTES: int = 1440  # 24 hours
    ENCRYPTION_KEY: str = os.getenv("ENCRYPTION_KEY", "CHANGE_THIS_IN_PRODUCTION")

    # Cloud Providers (default credentials)
    IBM_API_TOKEN: Optional[str] = None
    IBM_CRN: Optional[str] = None
    BLUEQUBIT_TOKEN: Optional[str] = None

    # AI Services
    ANTHROPIC_API_KEY: Optional[str] = None

    # Storage
    UPLOAD_DIR: str = "/tmp/kanad_uploads"
    MAX_UPLOAD_SIZE: int = 10 * 1024 * 1024  # 10 MB

    # Rate Limiting
    RATE_LIMIT_ENABLED: bool = True
    RATE_LIMIT_PER_MINUTE: int = 60

    # Logging
    LOG_LEVEL: str = "INFO"
    LOG_FORMAT: str = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

    class Config:
        env_file = ".env"
        case_sensitive = True


# Global settings instance
settings = Settings()


def get_allowed_origins() -> list:
    """Parse ALLOWED_ORIGINS into list."""
    return [origin.strip() for origin in settings.ALLOWED_ORIGINS.split(",")]
