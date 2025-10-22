"""
WebSocket endpoints for real-time experiment updates
"""

from fastapi import APIRouter, WebSocket, WebSocketDisconnect
from typing import Dict, Set
import asyncio
import json
import logging

logger = logging.getLogger(__name__)

router = APIRouter()

# Store active WebSocket connections per experiment
active_connections: Dict[str, Set[WebSocket]] = {}


class ConnectionManager:
    """Manage WebSocket connections for experiment updates."""

    def __init__(self):
        self.active_connections: Dict[str, Set[WebSocket]] = {}

    async def connect(self, websocket: WebSocket, experiment_id: str):
        """Accept and register a new WebSocket connection."""
        await websocket.accept()

        if experiment_id not in self.active_connections:
            self.active_connections[experiment_id] = set()

        self.active_connections[experiment_id].add(websocket)
        logger.info(f"✅ WebSocket connected for experiment {experiment_id}")
        print(f"✅ WebSocket connected for experiment {experiment_id}")
        print(f"   Total connections: {len(self.active_connections[experiment_id])}")

    def disconnect(self, websocket: WebSocket, experiment_id: str):
        """Remove a WebSocket connection."""
        if experiment_id in self.active_connections:
            self.active_connections[experiment_id].discard(websocket)

            # Clean up empty sets
            if not self.active_connections[experiment_id]:
                del self.active_connections[experiment_id]

        logger.info(f"WebSocket disconnected for experiment {experiment_id}")

    async def send_update(self, experiment_id: str, data: dict):
        """Send update to all connected clients for an experiment."""
        if experiment_id not in self.active_connections:
            return

        # Create copy to avoid modification during iteration
        connections = list(self.active_connections[experiment_id])

        for connection in connections:
            try:
                await connection.send_json(data)
            except Exception as e:
                logger.error(f"Error sending to WebSocket: {e}")
                # Remove dead connection
                self.disconnect(connection, experiment_id)

    async def broadcast_log(self, experiment_id: str, level: str, message: str):
        """Broadcast a log message to all clients."""
        await self.send_update(experiment_id, {
            "type": "log",
            "level": level,
            "message": message,
            "timestamp": asyncio.get_event_loop().time()
        })

    async def broadcast_convergence(self, experiment_id: str, iteration: int, energy: float,
                                   parameters: list = None, is_optimizer_iteration: bool = False):
        """Broadcast convergence data point."""
        data = {
            "type": "convergence",
            "iteration": iteration,
            "energy": energy,
            "is_optimizer_iteration": is_optimizer_iteration
        }

        if parameters is not None:
            data["parameters"] = parameters

        print(f"📊 Broadcasting convergence: iter={iteration}, E={energy:.8f}")
        await self.send_update(experiment_id, data)

    async def broadcast_status(self, experiment_id: str, status: str, progress: float = None):
        """Broadcast experiment status update."""
        data = {
            "type": "status",
            "status": status
        }

        if progress is not None:
            data["progress"] = progress

        await self.send_update(experiment_id, data)


# Global connection manager instance
manager = ConnectionManager()


@router.websocket("/ws/experiments/{experiment_id}")
async def experiment_websocket(websocket: WebSocket, experiment_id: str):
    """
    WebSocket endpoint for real-time experiment updates.

    Sends:
    - Convergence data (iteration, energy, parameters)
    - Log messages (info, warning, error)
    - Status updates (running, paused, completed, failed)
    """
    await manager.connect(websocket, experiment_id)

    try:
        # Keep connection alive and handle incoming messages if needed
        while True:
            # Wait for messages (currently we only send, but could receive commands)
            data = await websocket.receive_text()

            # Echo back for testing
            logger.info(f"Received from client: {data}")

    except WebSocketDisconnect:
        manager.disconnect(websocket, experiment_id)
        logger.info(f"Client disconnected from experiment {experiment_id}")
    except Exception as e:
        logger.error(f"WebSocket error for experiment {experiment_id}: {e}")
        manager.disconnect(websocket, experiment_id)


# Export manager for use in other modules
__all__ = ['router', 'manager']
