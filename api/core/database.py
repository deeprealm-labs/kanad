"""
SQLite database for storing experiments and jobs
"""

import sqlite3
import json
from datetime import datetime
from typing import Optional, List, Dict, Any
from contextlib import contextmanager

from api.core.config import get_settings


def init_db():
    """Initialize database schema."""
    settings = get_settings()

    with get_db() as conn:
        cursor = conn.cursor()

        # Experiments table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS experiments (
                id TEXT PRIMARY KEY,
                molecule_data TEXT NOT NULL,
                configuration TEXT NOT NULL,
                status TEXT NOT NULL,
                method TEXT NOT NULL,
                backend TEXT NOT NULL,
                results TEXT,
                error_message TEXT,
                created_at TEXT NOT NULL,
                started_at TEXT,
                completed_at TEXT
            )
        """)

        # Jobs queue table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS jobs (
                id TEXT PRIMARY KEY,
                experiment_id TEXT NOT NULL,
                status TEXT NOT NULL,
                priority INTEGER DEFAULT 0,
                progress REAL DEFAULT 0.0,
                current_iteration INTEGER,
                max_iterations INTEGER,
                current_energy REAL,
                best_energy REAL,
                scheduled_time TEXT,
                created_at TEXT NOT NULL,
                started_at TEXT,
                completed_at TEXT,
                FOREIGN KEY (experiment_id) REFERENCES experiments(id)
            )
        """)

        # Settings table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS user_settings (
                id INTEGER PRIMARY KEY CHECK (id = 1),
                settings TEXT NOT NULL,
                updated_at TEXT NOT NULL
            )
        """)

        # Cloud credentials table (encrypted in production)
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS cloud_credentials (
                provider TEXT PRIMARY KEY,
                credentials TEXT NOT NULL,
                updated_at TEXT NOT NULL
            )
        """)

        conn.commit()


@contextmanager
def get_db():
    """Get database connection context manager."""
    settings = get_settings()
    conn = sqlite3.connect(settings.DATABASE_PATH)
    conn.row_factory = sqlite3.Row
    try:
        yield conn
    finally:
        conn.close()


class ExperimentDB:
    """Database operations for experiments."""

    @staticmethod
    def create(experiment_data: Dict[str, Any]) -> str:
        """Create a new experiment."""
        with get_db() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                INSERT INTO experiments (
                    id, molecule_data, configuration, status, method, backend, created_at
                ) VALUES (?, ?, ?, ?, ?, ?, ?)
            """, (
                experiment_data['id'],
                json.dumps(experiment_data['molecule']),
                json.dumps(experiment_data['configuration']),
                experiment_data['status'],
                experiment_data['method'],
                experiment_data['backend'],
                datetime.utcnow().isoformat(),
            ))
            conn.commit()
            return experiment_data['id']

    @staticmethod
    def get(experiment_id: str) -> Optional[Dict[str, Any]]:
        """Get experiment by ID."""
        with get_db() as conn:
            cursor = conn.cursor()
            cursor.execute("SELECT * FROM experiments WHERE id = ?", (experiment_id,))
            row = cursor.fetchone()

            if not row:
                return None

            return {
                'id': row['id'],
                'molecule': json.loads(row['molecule_data']),
                'configuration': json.loads(row['configuration']),
                'status': row['status'],
                'method': row['method'],
                'backend': row['backend'],
                'results': json.loads(row['results']) if row['results'] else None,
                'error_message': row['error_message'],
                'created_at': row['created_at'],
                'started_at': row['started_at'],
                'completed_at': row['completed_at'],
            }

    @staticmethod
    def list(
        status: Optional[str] = None,
        limit: int = 50,
        offset: int = 0
    ) -> List[Dict[str, Any]]:
        """List experiments with optional filtering."""
        with get_db() as conn:
            cursor = conn.cursor()

            if status:
                cursor.execute("""
                    SELECT * FROM experiments
                    WHERE status = ?
                    ORDER BY created_at DESC
                    LIMIT ? OFFSET ?
                """, (status, limit, offset))
            else:
                cursor.execute("""
                    SELECT * FROM experiments
                    ORDER BY created_at DESC
                    LIMIT ? OFFSET ?
                """, (limit, offset))

            rows = cursor.fetchall()
            return [
                {
                    'id': row['id'],
                    'molecule': json.loads(row['molecule_data']),
                    'configuration': json.loads(row['configuration']),
                    'status': row['status'],
                    'method': row['method'],
                    'backend': row['backend'],
                    'results': json.loads(row['results']) if row['results'] else None,
                    'error_message': row['error_message'],
                    'created_at': row['created_at'],
                    'started_at': row['started_at'],
                    'completed_at': row['completed_at'],
                }
                for row in rows
            ]

    @staticmethod
    def update_status(
        experiment_id: str,
        status: str,
        results: Optional[Dict[str, Any]] = None,
        error_message: Optional[str] = None
    ):
        """Update experiment status and results."""
        with get_db() as conn:
            cursor = conn.cursor()

            if status == 'running' and not results:
                cursor.execute("""
                    UPDATE experiments
                    SET status = ?, started_at = ?
                    WHERE id = ?
                """, (status, datetime.utcnow().isoformat(), experiment_id))

            elif status == 'completed' and results:
                cursor.execute("""
                    UPDATE experiments
                    SET status = ?, results = ?, completed_at = ?
                    WHERE id = ?
                """, (
                    status,
                    json.dumps(results),
                    datetime.utcnow().isoformat(),
                    experiment_id
                ))

            elif status == 'failed':
                cursor.execute("""
                    UPDATE experiments
                    SET status = ?, error_message = ?, completed_at = ?
                    WHERE id = ?
                """, (
                    status,
                    error_message,
                    datetime.utcnow().isoformat(),
                    experiment_id
                ))

            conn.commit()

    @staticmethod
    def delete(experiment_id: str):
        """Delete experiment."""
        with get_db() as conn:
            cursor = conn.cursor()
            cursor.execute("DELETE FROM experiments WHERE id = ?", (experiment_id,))
            conn.commit()


class JobDB:
    """Database operations for jobs."""

    @staticmethod
    def create(job_data: Dict[str, Any]) -> str:
        """Create a new job."""
        with get_db() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                INSERT INTO jobs (
                    id, experiment_id, status, priority, max_iterations, created_at
                ) VALUES (?, ?, ?, ?, ?, ?)
            """, (
                job_data['id'],
                job_data['experiment_id'],
                job_data['status'],
                job_data.get('priority', 0),
                job_data.get('max_iterations'),
                datetime.utcnow().isoformat(),
            ))
            conn.commit()
            return job_data['id']

    @staticmethod
    def get(job_id: str) -> Optional[Dict[str, Any]]:
        """Get job by ID."""
        with get_db() as conn:
            cursor = conn.cursor()
            cursor.execute("SELECT * FROM jobs WHERE id = ?", (job_id,))
            row = cursor.fetchone()

            if not row:
                return None

            return {
                'id': row['id'],
                'experiment_id': row['experiment_id'],
                'status': row['status'],
                'priority': row['priority'],
                'progress': row['progress'],
                'current_iteration': row['current_iteration'],
                'max_iterations': row['max_iterations'],
                'current_energy': row['current_energy'],
                'best_energy': row['best_energy'],
                'scheduled_time': row['scheduled_time'],
                'created_at': row['created_at'],
                'started_at': row['started_at'],
                'completed_at': row['completed_at'],
            }

    @staticmethod
    def list(status: Optional[str] = None) -> List[Dict[str, Any]]:
        """List jobs."""
        with get_db() as conn:
            cursor = conn.cursor()

            if status:
                cursor.execute("""
                    SELECT * FROM jobs
                    WHERE status = ?
                    ORDER BY priority DESC, created_at ASC
                """, (status,))
            else:
                cursor.execute("""
                    SELECT * FROM jobs
                    ORDER BY priority DESC, created_at ASC
                """)

            rows = cursor.fetchall()
            return [
                {
                    'id': row['id'],
                    'experiment_id': row['experiment_id'],
                    'status': row['status'],
                    'priority': row['priority'],
                    'progress': row['progress'],
                    'current_iteration': row['current_iteration'],
                    'max_iterations': row['max_iterations'],
                    'current_energy': row['current_energy'],
                    'best_energy': row['best_energy'],
                    'scheduled_time': row['scheduled_time'],
                    'created_at': row['created_at'],
                    'started_at': row['started_at'],
                    'completed_at': row['completed_at'],
                }
                for row in rows
            ]

    @staticmethod
    def update_progress(
        job_id: str,
        progress: float,
        current_iteration: Optional[int] = None,
        current_energy: Optional[float] = None,
        best_energy: Optional[float] = None
    ):
        """Update job progress."""
        with get_db() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                UPDATE jobs
                SET progress = ?, current_iteration = ?, current_energy = ?, best_energy = ?
                WHERE id = ?
            """, (progress, current_iteration, current_energy, best_energy, job_id))
            conn.commit()

    @staticmethod
    def update_status(job_id: str, status: str):
        """Update job status."""
        with get_db() as conn:
            cursor = conn.cursor()

            if status == 'running':
                cursor.execute("""
                    UPDATE jobs
                    SET status = ?, started_at = ?
                    WHERE id = ?
                """, (status, datetime.utcnow().isoformat(), job_id))
            elif status in ['completed', 'failed', 'cancelled']:
                cursor.execute("""
                    UPDATE jobs
                    SET status = ?, completed_at = ?
                    WHERE id = ?
                """, (status, datetime.utcnow().isoformat(), job_id))
            else:
                cursor.execute("""
                    UPDATE jobs SET status = ? WHERE id = ?
                """, (status, job_id))

            conn.commit()

    @staticmethod
    def delete(job_id: str):
        """Delete job."""
        with get_db() as conn:
            cursor = conn.cursor()
            cursor.execute("DELETE FROM jobs WHERE id = ?", (job_id,))
            conn.commit()
