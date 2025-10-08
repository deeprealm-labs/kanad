"""
SQLAlchemy database models for PostgreSQL.

These models define the database schema for persisting molecules, jobs,
results, and user data.
"""

from sqlalchemy import Column, String, Integer, Float, Boolean, Text, DateTime, JSON, ForeignKey, Enum as SQLEnum
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy.dialects.postgresql import UUID
from datetime import datetime
import uuid
import enum

Base = declarative_base()


# ===== Enums =====

class JobStatusEnum(str, enum.Enum):
    """Job status enumeration."""
    QUEUED = "queued"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


class ComputationMethodEnum(str, enum.Enum):
    """Computation method enumeration."""
    HF = "HF"
    VQE = "VQE"
    MP2 = "MP2"
    FCI = "FCI"
    SQD = "SQD"


class BackendTypeEnum(str, enum.Enum):
    """Backend type enumeration."""
    CLASSICAL = "classical"
    IBM_QUANTUM = "ibm_quantum"
    BLUEQUBIT = "bluequbit"


# ===== User Models =====

class User(Base):
    """User account."""
    __tablename__ = "users"

    user_id = Column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    email = Column(String(255), unique=True, nullable=False, index=True)
    hashed_password = Column(String(255), nullable=False)
    name = Column(String(255))
    institution = Column(String(255))
    research_field = Column(String(100))
    created_at = Column(DateTime, default=datetime.utcnow)
    last_login = Column(DateTime)

    # Relationships
    molecules = relationship("Molecule", back_populates="user", cascade="all, delete-orphan")
    simulations = relationship("Simulation", back_populates="user", cascade="all, delete-orphan")
    jobs = relationship("Job", back_populates="user", cascade="all, delete-orphan")
    credentials = relationship("CloudCredential", back_populates="user", cascade="all, delete-orphan")
    settings = relationship("UserSettings", back_populates="user", uselist=False, cascade="all, delete-orphan")


class UserSettings(Base):
    """User preferences and default settings."""
    __tablename__ = "user_settings"

    user_id = Column(UUID(as_uuid=True), ForeignKey("users.user_id", ondelete="CASCADE"), primary_key=True)
    default_basis = Column(String(50), default='sto-3g')
    default_method = Column(String(50), default='VQE')
    default_ansatz = Column(String(50), default='hardware_efficient')
    default_mapper = Column(String(50), default='jordan_wigner')
    default_backend = Column(String(50), default='classical')
    auto_analyze = Column(Boolean, default=True)
    theme = Column(String(50), default='light')
    preferences = Column(JSON)

    # Relationships
    user = relationship("User", back_populates="settings")


# ===== Molecule Models =====

class Molecule(Base):
    """Molecular structure."""
    __tablename__ = "molecules"

    molecule_id = Column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    user_id = Column(UUID(as_uuid=True), ForeignKey("users.user_id", ondelete="CASCADE"), nullable=False)
    name = Column(String(255))
    formula = Column(String(100), nullable=False)
    smiles = Column(String(500))
    inchi = Column(Text)
    geometry = Column(JSON, nullable=False)  # List of atoms with positions
    basis = Column(String(50), nullable=False)
    charge = Column(Integer, default=0)
    multiplicity = Column(Integer, default=1)
    created_at = Column(DateTime, default=datetime.utcnow)
    n_electrons = Column(Integer)
    n_orbitals = Column(Integer)
    n_qubits = Column(Integer)
    lewis_structure = Column(JSON)  # Bonds, lone pairs, charges
    preview_svg = Column(Text)  # Base64-encoded SVG

    # Relationships
    user = relationship("User", back_populates="molecules")
    simulations = relationship("Simulation", back_populates="molecule", cascade="all, delete-orphan")


# ===== Simulation & Job Models =====

class Simulation(Base):
    """Simulation configuration."""
    __tablename__ = "simulations"

    simulation_id = Column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    user_id = Column(UUID(as_uuid=True), ForeignKey("users.user_id", ondelete="CASCADE"), nullable=False)
    molecule_id = Column(UUID(as_uuid=True), ForeignKey("molecules.molecule_id", ondelete="CASCADE"), nullable=False)
    method = Column(SQLEnum(ComputationMethodEnum), nullable=False)
    ansatz = Column(String(50))
    mapper = Column(String(50))
    optimizer = Column(String(50))
    backend_type = Column(SQLEnum(BackendTypeEnum), nullable=False)
    backend_name = Column(String(100))
    configuration = Column(JSON, nullable=False)  # Full config dict
    created_at = Column(DateTime, default=datetime.utcnow)

    # Relationships
    user = relationship("User", back_populates="simulations")
    molecule = relationship("Molecule", back_populates="simulations")
    jobs = relationship("Job", back_populates="simulation", cascade="all, delete-orphan")


class Job(Base):
    """Computation job."""
    __tablename__ = "jobs"

    job_id = Column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    user_id = Column(UUID(as_uuid=True), ForeignKey("users.user_id", ondelete="CASCADE"), nullable=False)
    simulation_id = Column(UUID(as_uuid=True), ForeignKey("simulations.simulation_id", ondelete="CASCADE"), nullable=False)
    molecule_id = Column(UUID(as_uuid=True), ForeignKey("molecules.molecule_id", ondelete="CASCADE"), nullable=False)
    status = Column(SQLEnum(JobStatusEnum), default=JobStatusEnum.QUEUED, nullable=False, index=True)
    progress = Column(Integer, default=0)  # 0-100
    created_at = Column(DateTime, default=datetime.utcnow, index=True)
    started_at = Column(DateTime)
    completed_at = Column(DateTime)
    cloud_job_id = Column(String(255))  # External job ID (IBM/BlueQubit)
    error_message = Column(Text)
    celery_task_id = Column(String(255))  # Celery task ID for tracking

    # Relationships
    user = relationship("User", back_populates="jobs")
    simulation = relationship("Simulation", back_populates="jobs")
    molecule = relationship("Molecule")
    result = relationship("Result", back_populates="job", uselist=False, cascade="all, delete-orphan")
    logs = relationship("JobLog", back_populates="job", cascade="all, delete-orphan")


class Result(Base):
    """Computation results."""
    __tablename__ = "results"

    result_id = Column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    job_id = Column(UUID(as_uuid=True), ForeignKey("jobs.job_id", ondelete="CASCADE"), unique=True, nullable=False)
    energy = Column(Float)
    hf_energy = Column(Float)
    correlation_energy = Column(Float)
    n_iterations = Column(Integer)
    converged = Column(Boolean)
    convergence_history = Column(JSON)  # List of {iteration, energy}
    analysis = Column(JSON)  # All analysis results
    llm_report = Column(JSON)  # AI-generated report
    created_at = Column(DateTime, default=datetime.utcnow)

    # Relationships
    job = relationship("Job", back_populates="result")


class JobLog(Base):
    """Real-time job logs for WebSocket streaming."""
    __tablename__ = "job_logs"

    log_id = Column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    job_id = Column(UUID(as_uuid=True), ForeignKey("jobs.job_id", ondelete="CASCADE"), nullable=False, index=True)
    timestamp = Column(DateTime, default=datetime.utcnow)
    level = Column(String(20), default='INFO')  # INFO, WARNING, ERROR
    message = Column(Text, nullable=False)

    # Relationships
    job = relationship("Job", back_populates="logs")


# ===== Cloud Provider Models =====

class CloudCredential(Base):
    """Encrypted cloud provider credentials."""
    __tablename__ = "cloud_credentials"

    credential_id = Column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    user_id = Column(UUID(as_uuid=True), ForeignKey("users.user_id", ondelete="CASCADE"), nullable=False)
    provider = Column(String(50), nullable=False)  # 'ibm_quantum', 'bluequbit'
    encrypted_token = Column(Text, nullable=False)  # AES-256 encrypted
    encrypted_crn = Column(Text)  # For IBM Quantum
    verified = Column(Boolean, default=False)
    last_verified = Column(DateTime)
    created_at = Column(DateTime, default=datetime.utcnow)

    # Relationships
    user = relationship("User", back_populates="credentials")


# ===== Molecule Library =====

class MoleculeLibrary(Base):
    """Pre-built molecule library (public and user-specific)."""
    __tablename__ = "molecule_library"

    library_id = Column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    name = Column(String(255), nullable=False)
    category = Column(String(100), nullable=False, index=True)  # 'Small', 'Metallurgy', 'Bioscience', etc.
    formula = Column(String(100), nullable=False)
    smiles = Column(String(500))
    geometry = Column(JSON, nullable=False)
    is_public = Column(Boolean, default=True, index=True)
    user_id = Column(UUID(as_uuid=True), ForeignKey("users.user_id", ondelete="CASCADE"))  # NULL for public
    description = Column(Text)
    tags = Column(JSON)  # List of tags


# ===== Batch Scheduling =====

class Schedule(Base):
    """Batch job schedule."""
    __tablename__ = "schedules"

    schedule_id = Column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    user_id = Column(UUID(as_uuid=True), ForeignKey("users.user_id", ondelete="CASCADE"), nullable=False)
    name = Column(String(255), nullable=False)
    execution_mode = Column(String(50), default='sequential')  # 'sequential', 'parallel'
    priority = Column(String(50), default='normal')  # 'low', 'normal', 'high'
    status = Column(SQLEnum(JobStatusEnum), default=JobStatusEnum.QUEUED)
    created_at = Column(DateTime, default=datetime.utcnow)
    completed_at = Column(DateTime)

    # Relationships
    jobs = relationship("ScheduleJob", back_populates="schedule", cascade="all, delete-orphan")


class ScheduleJob(Base):
    """Association between schedules and jobs."""
    __tablename__ = "schedule_jobs"

    schedule_id = Column(UUID(as_uuid=True), ForeignKey("schedules.schedule_id", ondelete="CASCADE"), primary_key=True)
    job_id = Column(UUID(as_uuid=True), ForeignKey("jobs.job_id", ondelete="CASCADE"), primary_key=True)
    sequence_order = Column(Integer, nullable=False)

    # Relationships
    schedule = relationship("Schedule", back_populates="jobs")
    job = relationship("Job")
