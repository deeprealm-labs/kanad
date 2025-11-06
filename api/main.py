"""
FastAPI Backend Server for Kanad Quantum Chemistry Framework

This server provides REST API endpoints for the Kanad web frontend.
"""

import os
import sys
from typing import Optional
from contextlib import asynccontextmanager
from pathlib import Path
from dotenv import load_dotenv

# Load environment variables from .env file in project root
env_path = Path(__file__).parent.parent / '.env'
load_dotenv(dotenv_path=env_path, override=True)

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
import uvicorn

# Add parent directory to path to import kanad
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from api.routes import (
    molecules,
    experiments,
    jobs,
    analysis,
    settings as settings_router,
    library,
    cloud,
    health,
    configuration,
    campaigns,
    circuits,
    websockets,
    auth,  # Authentication routes
    admin,  # Admin dashboard routes
    users,  # User profile and account management
    applications,  # Application domain platforms
)
from api.core.config import get_settings
from api.core.database import init_db, cleanup_stuck_experiments
from api.core.database_postgres import init_db as init_postgres_db


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Initialize resources on startup, cleanup on shutdown."""
    # Startup
    import asyncio
    from api.utils import set_main_event_loop

    config = get_settings()
    print(f"🚀 Starting Kanad API Server v{config.VERSION}")
    print(f"📁 Database: {config.DATABASE_PATH}")

    # Store reference to main event loop for background tasks
    loop = asyncio.get_running_loop()
    set_main_event_loop(loop)

    # Initialize SQLite database (experiments, jobs, campaigns)
    init_db()
    print("✅ SQLite Database initialized")

    # Initialize PostgreSQL database (users, auth, admin)
    init_postgres_db()
    print("✅ PostgreSQL Database initialized")

    # Clean up stuck experiments from previous session
    cleanup_stuck_experiments()

    # Start rate limit cleanup task
    from api.middleware.rate_limit import start_cleanup_task, stop_cleanup_task
    start_cleanup_task()

    yield

    # Shutdown
    print("👋 Shutting down Kanad API Server")
    stop_cleanup_task()


# Create FastAPI app
app = FastAPI(
    title="Kanad Quantum Chemistry API",
    description="REST API for Kanad quantum chemistry calculations",
    version="0.1.0",
    lifespan=lifespan,
)

# Security and compute limits middleware
from api.middleware.security import SecurityHeadersMiddleware, RequestLoggingMiddleware
from api.middleware.compute_limits import ComputeLimitsMiddleware

app.add_middleware(ComputeLimitsMiddleware)  # Check compute limits first
app.add_middleware(SecurityHeadersMiddleware)
app.add_middleware(RequestLoggingMiddleware)

# CORS configuration
config = get_settings()
app.add_middleware(
    CORSMiddleware,
    allow_origins=config.CORS_ORIGINS,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
    expose_headers=["X-RateLimit-Limit", "X-RateLimit-Remaining", "X-RateLimit-Reset"],
)

# Include routers
app.include_router(health.router, tags=["Health"])
app.include_router(websockets.router, prefix="/api", tags=["WebSockets"])  # WebSocket endpoint

# Authentication & Admin routes (public)
app.include_router(auth.router, prefix="/api", tags=["Authentication"])
app.include_router(admin.router, prefix="/api", tags=["Admin"])
app.include_router(users.router, prefix="/api", tags=["Users"])

# Application routes
app.include_router(molecules.router, prefix="/api/molecules", tags=["Molecules"])
app.include_router(experiments.router, prefix="/api/experiments", tags=["Experiments"])
app.include_router(jobs.router, prefix="/api/jobs", tags=["Jobs"])
app.include_router(campaigns.router, prefix="/api/campaigns", tags=["Campaigns"])
app.include_router(analysis.router, prefix="/api/analysis", tags=["Analysis"])
app.include_router(settings_router.router, prefix="/api/settings", tags=["Settings"])
app.include_router(library.router, prefix="/api/library", tags=["Library"])
app.include_router(cloud.router, prefix="/api/cloud", tags=["Cloud"])
app.include_router(configuration.router, prefix="/api/configuration", tags=["Configuration"])
app.include_router(circuits.router, prefix="/api/circuits", tags=["Circuits"])
app.include_router(applications.router, prefix="/api/applications", tags=["Applications"])


@app.exception_handler(HTTPException)
async def http_exception_handler(request, exc):
    """Custom HTTP exception handler."""
    return JSONResponse(
        status_code=exc.status_code,
        content={
            "error": exc.detail,
            "status_code": exc.status_code,
        },
    )


@app.exception_handler(Exception)
async def general_exception_handler(request, exc):
    """Handle unexpected exceptions."""
    import traceback
    print(f"❌ Unexpected error: {exc}")
    traceback.print_exc()

    return JSONResponse(
        status_code=500,
        content={
            "error": "Internal server error",
            "message": str(exc),
            "status_code": 500,
        },
    )


if __name__ == "__main__":
    config = get_settings()
    uvicorn.run(
        "main:app",
        host=config.HOST,
        port=config.PORT,
        reload=config.DEBUG,
        log_level="info",
    )
