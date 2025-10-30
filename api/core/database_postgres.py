"""
PostgreSQL Database Configuration for Production
Supports user authentication, access keys, and multi-user experiments
"""

from sqlalchemy import (
    create_engine,
    Column,
    Integer,
    String,
    Float,
    Text,
    DateTime,
    Boolean,
    ForeignKey,
    JSON,
    Enum as SQLEnum,
)
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker, relationship
from datetime import datetime
import enum
import os

# Database configuration
DATABASE_URL = os.getenv(
    "DATABASE_URL",
    "postgresql://kanad_user:your_password@localhost:5432/kanad_db"
)

engine = create_engine(DATABASE_URL, pool_pre_ping=True, pool_size=10, max_overflow=20)
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)
Base = declarative_base()


# Enums
class UserRole(str, enum.Enum):
    ADMIN = "admin"
    USER = "user"
    VIEWER = "viewer"


class ExperimentStatus(str, enum.Enum):
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


# Authentication Tables
class User(Base):
    __tablename__ = "users"

    id = Column(Integer, primary_key=True, index=True)
    email = Column(String(255), unique=True, index=True, nullable=False)
    password_hash = Column(String(255), nullable=False)
    full_name = Column(String(255), nullable=True)
    role = Column(SQLEnum(UserRole, values_callable=lambda x: [e.value for e in x]), default=UserRole.USER, nullable=False)
    is_verified = Column(Boolean, default=False, nullable=False)
    is_active = Column(Boolean, default=True, nullable=False)
    access_key_id = Column(Integer, ForeignKey("access_keys.id"), nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)
    last_login = Column(DateTime, nullable=True)

    # OAuth fields
    google_id = Column(String(255), unique=True, nullable=True, index=True)
    avatar_url = Column(Text, nullable=True)

    # Relationships
    access_key = relationship("AccessKey", foreign_keys=[access_key_id], back_populates="users")
    sessions = relationship("Session", back_populates="user", cascade="all, delete-orphan")
    experiments = relationship("Experiment", back_populates="user")
    campaigns = relationship("Campaign", back_populates="user")


class AccessKey(Base):
    __tablename__ = "access_keys"

    id = Column(Integer, primary_key=True, index=True)
    key = Column(String(64), unique=True, index=True, nullable=False)
    description = Column(Text, nullable=True)
    max_uses = Column(Integer, default=1, nullable=False)
    used_count = Column(Integer, default=0, nullable=False)
    expires_at = Column(DateTime, nullable=True)
    is_active = Column(Boolean, default=True, nullable=False)
    created_by = Column(Integer, ForeignKey("users.id"), nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)

    # Relationships
    users = relationship("User", foreign_keys="[User.access_key_id]", back_populates="access_key")


class Session(Base):
    __tablename__ = "sessions"

    id = Column(Integer, primary_key=True, index=True)
    user_id = Column(Integer, ForeignKey("users.id"), nullable=False)
    access_token = Column(String(512), unique=True, index=True, nullable=False)
    refresh_token = Column(String(512), unique=True, index=True, nullable=False)
    expires_at = Column(DateTime, nullable=False)
    refresh_expires_at = Column(DateTime, nullable=False)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    last_used = Column(DateTime, default=datetime.utcnow, nullable=False)
    user_agent = Column(Text, nullable=True)
    ip_address = Column(String(45), nullable=True)
    is_revoked = Column(Boolean, default=False, nullable=False)

    # Relationships
    user = relationship("User", back_populates="sessions")


class EmailVerification(Base):
    __tablename__ = "email_verifications"

    id = Column(Integer, primary_key=True, index=True)
    email = Column(String(255), index=True, nullable=False)
    otp = Column(String(6), nullable=False)
    expires_at = Column(DateTime, nullable=False)
    is_used = Column(Boolean, default=False, nullable=False)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)


# Existing Tables with User Association
class Campaign(Base):
    __tablename__ = "campaigns"

    id = Column(Integer, primary_key=True, index=True)
    user_id = Column(Integer, ForeignKey("users.id"), nullable=True)  # Nullable for backward compatibility
    name = Column(String(255), nullable=False)
    description = Column(Text, nullable=True)
    parameters = Column(JSON, nullable=True)
    status = Column(String(50), default="active", nullable=False)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)

    # Relationships
    user = relationship("User", back_populates="campaigns")
    experiments = relationship("Experiment", back_populates="campaign", cascade="all, delete-orphan")


class Experiment(Base):
    __tablename__ = "experiments"

    id = Column(Integer, primary_key=True, index=True)
    user_id = Column(Integer, ForeignKey("users.id"), nullable=True)  # Nullable for backward compatibility
    campaign_id = Column(Integer, ForeignKey("campaigns.id"), nullable=True)
    name = Column(String(255), nullable=False)
    molecule_name = Column(String(255), nullable=True)
    smiles = Column(Text, nullable=True)
    xyz_data = Column(Text, nullable=True)
    atoms = Column(JSON, nullable=True)
    basis_set = Column(String(50), nullable=True)
    charge = Column(Integer, default=0)
    multiplicity = Column(Integer, default=1)
    backend = Column(String(50), nullable=True)
    method = Column(String(50), nullable=True)
    status = Column(SQLEnum(ExperimentStatus, values_callable=lambda x: [e.value for e in x]), default=ExperimentStatus.PENDING, nullable=False)
    progress = Column(Float, default=0.0)
    error_message = Column(Text, nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    started_at = Column(DateTime, nullable=True)
    completed_at = Column(DateTime, nullable=True)

    # Relationships
    user = relationship("User", back_populates="experiments")
    campaign = relationship("Campaign", back_populates="experiments")
    jobs = relationship("Job", back_populates="experiment", cascade="all, delete-orphan")
    analysis_results = relationship("AnalysisResult", back_populates="experiment", cascade="all, delete-orphan")


class Job(Base):
    __tablename__ = "jobs"

    id = Column(Integer, primary_key=True, index=True)
    experiment_id = Column(Integer, ForeignKey("experiments.id"), nullable=False)
    job_id = Column(String(255), unique=True, index=True, nullable=False)
    backend = Column(String(50), nullable=True)
    status = Column(String(50), default="pending", nullable=False)
    progress = Column(Float, default=0.0)
    result = Column(JSON, nullable=True)
    error = Column(Text, nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    started_at = Column(DateTime, nullable=True)
    completed_at = Column(DateTime, nullable=True)

    # Relationships
    experiment = relationship("Experiment", back_populates="jobs")


class AnalysisResult(Base):
    __tablename__ = "analysis_results"

    id = Column(Integer, primary_key=True, index=True)
    experiment_id = Column(Integer, ForeignKey("experiments.id"), nullable=False)
    analysis_type = Column(String(100), nullable=False)
    result_data = Column(JSON, nullable=True)
    file_path = Column(Text, nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)

    # Relationships
    experiment = relationship("Experiment", back_populates="analysis_results")


class UserSettings(Base):
    __tablename__ = "user_settings"

    id = Column(Integer, primary_key=True, index=True)
    user_id = Column(Integer, ForeignKey("users.id"), nullable=True)  # NULL = global settings
    settings = Column(JSON, nullable=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)


class CloudCredentials(Base):
    __tablename__ = "cloud_credentials"

    id = Column(Integer, primary_key=True, index=True)
    user_id = Column(Integer, ForeignKey("users.id"), nullable=True)  # NULL = global credentials
    provider = Column(String(50), nullable=False)
    credentials = Column(JSON, nullable=False)
    is_active = Column(Boolean, default=True, nullable=False)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)


# Database initialization
def init_db():
    """Create all tables in the database"""
    Base.metadata.create_all(bind=engine)
    print("✓ Database tables created successfully")


def get_db():
    """FastAPI dependency for database sessions"""
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()


# Utility function to create admin user
def create_admin_user(db, email: str, password: str, full_name: str = "Admin"):
    """Create the first admin user (for initialization only)"""
    from api.auth.password import hash_password

    existing_user = db.query(User).filter(User.email == email).first()
    if existing_user:
        print(f"✗ Admin user {email} already exists")
        return existing_user

    admin = User(
        email=email,
        password_hash=hash_password(password),
        full_name=full_name,
        role=UserRole.ADMIN,
        is_verified=True,
        is_active=True
    )

    db.add(admin)
    db.commit()
    db.refresh(admin)

    print(f"✓ Admin user created: {email}")
    return admin
