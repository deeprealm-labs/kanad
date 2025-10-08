"""
Custom exceptions and error handlers for the API.
"""

from fastapi import HTTPException, Request, status
from fastapi.responses import JSONResponse
from fastapi.exceptions import RequestValidationError


class ExperimentNotFoundError(HTTPException):
    """Raised when experiment is not found."""

    def __init__(self, experiment_id: int):
        super().__init__(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Experiment with ID {experiment_id} not found"
        )


class QueueItemNotFoundError(HTTPException):
    """Raised when queue item is not found."""

    def __init__(self, queue_id: int):
        super().__init__(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Queue item with ID {queue_id} not found"
        )


class ExperimentExecutionError(HTTPException):
    """Raised when experiment execution fails."""

    def __init__(self, message: str):
        super().__init__(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Experiment execution failed: {message}"
        )


class InvalidMoleculeError(HTTPException):
    """Raised when molecule data is invalid."""

    def __init__(self, message: str):
        super().__init__(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Invalid molecule data: {message}"
        )


class InvalidConfigurationError(HTTPException):
    """Raised when configuration is invalid."""

    def __init__(self, message: str):
        super().__init__(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Invalid configuration: {message}"
        )


class ExperimentCancelledException(Exception):
    """Raised when an experiment is cancelled by user request."""

    def __init__(self, experiment_id: int, message: str = "Experiment cancelled by user"):
        self.experiment_id = experiment_id
        self.message = message
        super().__init__(self.message)


# Error handlers

async def validation_exception_handler(request: Request, exc: RequestValidationError):
    """Handle validation errors with detailed messages."""
    errors = []
    for error in exc.errors():
        errors.append({
            "field": " -> ".join(str(x) for x in error["loc"]),
            "message": error["msg"],
            "type": error["type"]
        })

    return JSONResponse(
        status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
        content={
            "detail": "Validation error",
            "errors": errors
        }
    )


async def general_exception_handler(request: Request, exc: Exception):
    """Handle general exceptions."""
    return JSONResponse(
        status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
        content={
            "detail": "Internal server error",
            "message": str(exc)
        }
    )
