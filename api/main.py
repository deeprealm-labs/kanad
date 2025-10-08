"""
Kanad Quantum Chemistry API - FastAPI Application

Main entry point for the backend server.
"""

import logging
from contextlib import asynccontextmanager

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.exceptions import RequestValidationError

from api.config import settings
from api.database import init_db
from api.services.job_queue import job_queue
from api.utils.exceptions import validation_exception_handler, general_exception_handler

# Import routers
from api.routers import experiments, queue, molecules, settings as settings_router

# Configure logging
logging.basicConfig(
    level=logging.INFO if not settings.DEBUG else logging.DEBUG,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@asynccontextmanager
async def lifespan(app: FastAPI):
    """
    Lifespan context manager for startup and shutdown events.
    """
    # Startup
    logger.info("Starting Kanad API server...")

    # Initialize database
    init_db()
    logger.info("Database initialized")

    # Start job queue
    job_queue.start()
    logger.info("Job queue started")

    yield

    # Shutdown
    logger.info("Shutting down Kanad API server...")
    job_queue.stop()
    logger.info("Job queue stopped")


# Create FastAPI application
app = FastAPI(
    title=settings.APP_NAME,
    version=settings.APP_VERSION,
    description="FastAPI backend for Kanad quantum chemistry calculations",
    lifespan=lifespan,
    docs_url="/docs",
    redoc_url="/redoc"
)

# Configure CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=settings.CORS_ORIGINS,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Register exception handlers
app.add_exception_handler(RequestValidationError, validation_exception_handler)
app.add_exception_handler(Exception, general_exception_handler)

# Include routers
app.include_router(experiments.router, prefix=settings.API_V1_PREFIX)
app.include_router(queue.router, prefix=settings.API_V1_PREFIX)
app.include_router(molecules.router, prefix=settings.API_V1_PREFIX)
app.include_router(settings_router.router, prefix=settings.API_V1_PREFIX)


@app.get("/")
def root():
    """Root endpoint - API information."""
    return {
        "name": settings.APP_NAME,
        "version": settings.APP_VERSION,
        "status": "running",
        "docs": "/docs",
        "api_prefix": settings.API_V1_PREFIX
    }


@app.get("/health")
def health_check():
    """Health check endpoint."""
    return {
        "status": "healthy",
        "job_queue_running": job_queue.is_running(),
        "queue_size": job_queue.get_queue_size()
    }


@app.get(f"{settings.API_V1_PREFIX}/info")
def api_info():
    """API information and available capabilities."""
    return {
        "version": settings.APP_VERSION,
        "endpoints": {
            "experiments": f"{settings.API_V1_PREFIX}/experiments",
            "queue": f"{settings.API_V1_PREFIX}/queue",
            "queue_stats": f"{settings.API_V1_PREFIX}/queue/stats",
            "molecules": f"{settings.API_V1_PREFIX}/molecules",
            "settings": f"{settings.API_V1_PREFIX}/settings",
            "circuit_visualization": f"{settings.API_V1_PREFIX}/experiments/{{id}}/circuit",
            "experiment_report": f"{settings.API_V1_PREFIX}/experiments/{{id}}/report",
            "convergence_data": f"{settings.API_V1_PREFIX}/experiments/{{id}}/convergence"
        },
        "capabilities": {
            "methods": {
                "VQE": "Variational Quantum Eigensolver (ground state)",
                "SQD": "Subspace Quantum Diagonalization (ground + excited states)",
                "EXCITED_STATES": "Excited states solver (CIS, TDDFT)",
                "HF": "Hartree-Fock (classical reference)"
            },
            "ansatze": {
                "ucc": "Unitary Coupled Cluster (UCCSD)",
                "hardware_efficient": "Hardware-Efficient Ansatz (layered)",
                "governance": "Governance-Aware Ansatz (bond-type specific)",
                "two_local": "Two-Local Ansatz (customizable)"
            },
            "mappers": {
                "jordan_wigner": "Jordan-Wigner transformation",
                "bravyi_kitaev": "Bravyi-Kitaev transformation",
                "hybrid_orbital": "Hybrid Orbital mapping"
            },
            "excited_methods": {
                "cis": "Configuration Interaction Singles",
                "tddft": "Time-Dependent DFT"
            },
            "optimizers": ["SLSQP", "COBYLA", "L-BFGS-B", "ADAM", "POWELL"],
            "backends": {
                "classical": "Classical statevector simulation (exact)",
                "ibm_quantum": "IBM Quantum hardware/simulators",
                "bluequbit": "BlueQubit cloud platform"
            },
            "basis_sets": ["sto-3g", "6-31g", "6-31g*", "6-31g**", "cc-pvdz", "cc-pvtz"]
        },
        "features": {
            "circuit_visualization": ["json", "ascii", "qasm"],
            "report_formats": ["json", "markdown"],
            "real_time_convergence": True,
            "job_queue": True,
            "queue_statistics": True
        }
    }


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(
        "main:app",
        host="0.0.0.0",
        port=8000,
        reload=settings.DEBUG,
        log_level="info"
    )
