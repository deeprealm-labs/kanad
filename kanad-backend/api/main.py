"""
Kanad FastAPI Application Entry Point.

Main application with routing, middleware, and configuration.
"""

from fastapi import FastAPI, Request, status
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from fastapi.exceptions import RequestValidationError
from slowapi import Limiter, _rate_limit_exceeded_handler
from slowapi.util import get_remote_address
from slowapi.errors import RateLimitExceeded
import logging
import time
from contextlib import asynccontextmanager

from api.config import settings, get_allowed_origins
from core.database import init_db

# Import routers
from api.routers import (
    auth,
    molecules,
    simulations,
    jobs,
    analysis,
    cloud,
    library,
    settings as settings_router,
    metallurgy,
    bioscience,
    chemical_engineering,
    schedules,
    user
)

# Configure logging
logging.basicConfig(
    level=getattr(logging, settings.LOG_LEVEL),
    format=settings.LOG_FORMAT
)
logger = logging.getLogger(__name__)


# Rate limiter
limiter = Limiter(key_func=get_remote_address)


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Lifecycle management (startup/shutdown)."""
    # Startup
    logger.info("Starting Kanad API server...")
    init_db()  # Initialize database tables
    logger.info(f"Database initialized: {settings.DATABASE_URL}")
    logger.info("Kanad API is ready!")

    yield

    # Shutdown
    logger.info("Shutting down Kanad API...")


# Create FastAPI app
app = FastAPI(
    title=settings.APP_NAME,
    version=settings.APP_VERSION,
    description="FastAPI backend for Kanad quantum chemistry framework",
    docs_url=f"{settings.API_PREFIX}/docs",
    redoc_url=f"{settings.API_PREFIX}/redoc",
    openapi_url=f"{settings.API_PREFIX}/openapi.json",
    lifespan=lifespan
)

# Add rate limiter state
app.state.limiter = limiter
app.add_exception_handler(RateLimitExceeded, _rate_limit_exceeded_handler)


# ===== Middleware =====

# CORS Middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=get_allowed_origins(),
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
    expose_headers=["*"]
)


# Request timing middleware
@app.middleware("http")
async def add_process_time_header(request: Request, call_next):
    """Add X-Process-Time header with request duration."""
    start_time = time.time()
    response = await call_next(request)
    process_time = time.time() - start_time
    response.headers["X-Process-Time"] = str(process_time)
    return response


# ===== Exception Handlers =====

@app.exception_handler(RequestValidationError)
async def validation_exception_handler(request: Request, exc: RequestValidationError):
    """Handle Pydantic validation errors with detailed messages."""
    return JSONResponse(
        status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
        content={
            "detail": exc.errors(),
            "message": "Request validation failed"
        }
    )


@app.exception_handler(Exception)
async def general_exception_handler(request: Request, exc: Exception):
    """Catch-all exception handler."""
    logger.error(f"Unhandled exception: {exc}", exc_info=True)
    return JSONResponse(
        status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
        content={
            "detail": "Internal server error",
            "message": str(exc) if settings.DEBUG else "An error occurred"
        }
    )


# ===== Root Routes =====

@app.get("/")
async def root():
    """API root endpoint."""
    return {
        "name": settings.APP_NAME,
        "version": settings.APP_VERSION,
        "status": "online",
        "docs": f"{settings.API_PREFIX}/docs"
    }


@app.get("/health")
async def health_check():
    """Health check endpoint for monitoring."""
    return {
        "status": "healthy",
        "environment": settings.ENVIRONMENT,
        "version": settings.APP_VERSION
    }


@app.get("/ready")
async def readiness_check():
    """Readiness check for Kubernetes."""
    # Check database connection
    try:
        from core.database import SessionLocal
        db = SessionLocal()
        db.execute("SELECT 1")
        db.close()
        db_status = "connected"
    except Exception as e:
        logger.error(f"Database connection failed: {e}")
        return JSONResponse(
            status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
            content={"status": "not ready", "database": "disconnected"}
        )

    return {"status": "ready", "database": db_status}


# ===== Include Routers =====

app.include_router(auth.router, prefix=f"{settings.API_PREFIX}/auth", tags=["Authentication"])
app.include_router(molecules.router, prefix=f"{settings.API_PREFIX}/molecules", tags=["Molecules"])
app.include_router(simulations.router, prefix=f"{settings.API_PREFIX}/simulations", tags=["Simulations"])
app.include_router(jobs.router, prefix=f"{settings.API_PREFIX}/jobs", tags=["Jobs"])
app.include_router(analysis.router, prefix=f"{settings.API_PREFIX}/analysis", tags=["Analysis"])
app.include_router(cloud.router, prefix=f"{settings.API_PREFIX}/cloud", tags=["Cloud"])
app.include_router(library.router, prefix=f"{settings.API_PREFIX}/library", tags=["Library"])
app.include_router(settings_router.router, prefix=f"{settings.API_PREFIX}/settings", tags=["Settings"])

# Domain-specific routers
app.include_router(metallurgy.router, prefix=f"{settings.API_PREFIX}/metallurgy", tags=["Metallurgy"])
app.include_router(bioscience.router, prefix=f"{settings.API_PREFIX}/bioscience", tags=["Bioscience"])
app.include_router(chemical_engineering.router, prefix=f"{settings.API_PREFIX}/chemical-engineering", tags=["Chemical Engineering"])

# Additional routers
app.include_router(schedules.router, prefix=f"{settings.API_PREFIX}/schedules", tags=["Schedules"])
app.include_router(user.router, prefix=f"{settings.API_PREFIX}/user", tags=["User"])


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(
        "api.main:app",
        host="0.0.0.0",
        port=8000,
        reload=settings.DEBUG,
        log_level=settings.LOG_LEVEL.lower()
    )
