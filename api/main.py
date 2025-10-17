"""
FastAPI Backend Server for Kanad Quantum Chemistry Framework

This server provides REST API endpoints for the Kanad web frontend.
"""

import os
import sys
from typing import Optional
from contextlib import asynccontextmanager

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
    circuits
)
from api.core.config import get_settings
from api.core.database import init_db, cleanup_stuck_experiments


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Initialize resources on startup, cleanup on shutdown."""
    # Startup
    config = get_settings()
    print(f"🚀 Starting Kanad API Server v{config.VERSION}")
    print(f"📁 Database: {config.DATABASE_PATH}")

    # Initialize database
    init_db()
    print("✅ Database initialized")

    # Clean up stuck experiments from previous session
    cleanup_stuck_experiments()

    yield

    # Shutdown
    print("👋 Shutting down Kanad API Server")


# Create FastAPI app
app = FastAPI(
    title="Kanad Quantum Chemistry API",
    description="REST API for Kanad quantum chemistry calculations",
    version="0.1.0",
    lifespan=lifespan,
)

# CORS configuration
config = get_settings()
app.add_middleware(
    CORSMiddleware,
    allow_origins=config.CORS_ORIGINS,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Include routers
app.include_router(health.router, tags=["Health"])
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
