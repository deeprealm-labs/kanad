"""
WebSocket utility functions for broadcasting logs to frontend.
"""

import asyncio
from typing import Optional

# Global reference to main event loop (set by main.py on startup)
_main_event_loop: Optional[asyncio.AbstractEventLoop] = None


def set_main_event_loop(loop: asyncio.AbstractEventLoop):
    """Store reference to main event loop for cross-thread WebSocket broadcasting."""
    global _main_event_loop
    _main_event_loop = loop


def broadcast_log_sync(experiment_id: str, message: str, level: str = "info"):
    """
    Broadcast log message to frontend from sync context.

    This function safely broadcasts log messages via WebSocket updates from synchronous code
    by submitting the coroutine to the main event loop.

    Args:
        experiment_id: Experiment ID to broadcast to
        message: Log message to send
        level: Log level ("info", "warning", "error")
    """
    global _main_event_loop

    # Always print to console for backend logs
    print(message)

    try:
        # Try to get the stored main loop first
        loop = _main_event_loop

        # Fallback: try to get running loop (works in async context)
        if loop is None:
            try:
                loop = asyncio.get_running_loop()
            except RuntimeError:
                pass

        # If still no loop, we can't broadcast (but we already printed)
        if loop is None:
            return

        # Import here to avoid circular dependency
        from api.routes.websockets import manager as ws_manager

        # Create the coroutine
        coro = ws_manager.broadcast_log(
            experiment_id,
            level=level,
            message=message
        )

        # Schedule it on the main loop (thread-safe)
        future = asyncio.run_coroutine_threadsafe(coro, loop)

        # Wait briefly for completion (non-blocking)
        future.result(timeout=0.5)

    except Exception as e:
        # Silently fail - we already printed to console
        pass
