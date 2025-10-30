"""
SQLite database for storing experiments and jobs
"""

import sqlite3
import json
from datetime import datetime
from typing import Optional, List, Dict, Any
from contextlib import contextmanager

from api.core.config import get_settings


def cleanup_stuck_experiments():
    """
    Mark experiments stuck in 'running' or 'queued' status as 'failed'.
    This is called on server startup to clean up zombie experiments from previous crashes.
    """
    from datetime import datetime, timedelta

    with get_db() as conn:
        cursor = conn.cursor()

        # Mark old running/queued experiments as failed (older than 10 minutes)
        ten_min_ago = (datetime.now() - timedelta(minutes=10)).isoformat()

        cursor.execute('''
            UPDATE experiments
            SET status = 'failed',
                error_message = 'Experiment was interrupted by server restart'
            WHERE status IN ('running', 'queued')
            AND created_at < ?
        ''', (ten_min_ago,))

        stuck_count = cursor.rowcount

        # Also mark corresponding jobs as failed
        cursor.execute('''
            UPDATE jobs
            SET status = 'failed'
            WHERE status IN ('running', 'queued')
            AND created_at < ?
        ''', (ten_min_ago,))

        conn.commit()

        if stuck_count > 0:
            print(f"🧹 Cleaned up {stuck_count} stuck experiments from previous session")


def init_db():
    """Initialize database schema."""
    settings = get_settings()

    with get_db() as conn:
        cursor = conn.cursor()

        # Campaigns table (for batch/queue executions)
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS campaigns (
                id TEXT PRIMARY KEY,
                name TEXT NOT NULL,
                description TEXT,
                status TEXT NOT NULL,
                total_experiments INTEGER DEFAULT 0,
                completed_experiments INTEGER DEFAULT 0,
                failed_experiments INTEGER DEFAULT 0,
                created_at TEXT NOT NULL,
                started_at TEXT,
                completed_at TEXT
            )
        """)

        # Experiments table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS experiments (
                id TEXT PRIMARY KEY,
                campaign_id TEXT,
                molecule_data TEXT NOT NULL,
                configuration TEXT NOT NULL,
                status TEXT NOT NULL,
                method TEXT NOT NULL,
                backend TEXT NOT NULL,
                results TEXT,
                error_message TEXT,
                sequence_order INTEGER DEFAULT 0,
                created_at TEXT NOT NULL,
                started_at TEXT,
                completed_at TEXT,
                FOREIGN KEY (campaign_id) REFERENCES campaigns(id)
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

        # Analysis results table (for on-demand analysis)
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS analysis_results (
                id TEXT PRIMARY KEY,
                experiment_id TEXT NOT NULL,
                profile TEXT NOT NULL,
                analyses TEXT NOT NULL,
                parameters TEXT,
                results TEXT NOT NULL,
                status TEXT NOT NULL,
                computation_time REAL,
                created_at TEXT NOT NULL,
                FOREIGN KEY (experiment_id) REFERENCES experiments(id)
            )
        """)

        # Index for faster lookups
        cursor.execute("""
            CREATE INDEX IF NOT EXISTS idx_analysis_experiment
            ON analysis_results(experiment_id)
        """)
        cursor.execute("""
            CREATE INDEX IF NOT EXISTS idx_analysis_profile
            ON analysis_results(profile)
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
                    id, molecule_data, configuration, status, method, backend, user_id, created_at
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                experiment_data['id'],
                json.dumps(experiment_data['molecule']),
                json.dumps(experiment_data['configuration']),
                experiment_data['status'],
                experiment_data['method'],
                experiment_data['backend'],
                experiment_data.get('user_id'),  # Add user_id
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
                'campaign_id': row['campaign_id'],
                'sequence_order': row['sequence_order'],
                'results': json.loads(row['results']) if row['results'] else None,
                'error_message': row['error_message'],
                'user_id': row['user_id'],  # Add user_id
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
                    'campaign_id': row['campaign_id'],
                    'sequence_order': row['sequence_order'],
                    'results': json.loads(row['results']) if row['results'] else None,
                    'error_message': row['error_message'],
                    'user_id': row['user_id'],  # Add user_id
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

            elif status == 'cancelled':
                # Handle cancellation with optional results (partial results if available)
                if results:
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
                else:
                    cursor.execute("""
                        UPDATE experiments
                        SET status = ?, completed_at = ?
                        WHERE id = ?
                    """, (
                        status,
                        datetime.utcnow().isoformat(),
                        experiment_id
                    ))

            conn.commit()

    @staticmethod
    def update_campaign_info(experiment_id: str, campaign_id: str, sequence_order: int):
        """Update experiment campaign info."""
        with get_db() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                UPDATE experiments
                SET campaign_id = ?, sequence_order = ?
                WHERE id = ?
            """, (campaign_id, sequence_order, experiment_id))
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
    def get_by_experiment_id(experiment_id: str) -> Optional[Dict[str, Any]]:
        """Get job by experiment ID."""
        with get_db() as conn:
            cursor = conn.cursor()
            cursor.execute("SELECT * FROM jobs WHERE experiment_id = ? ORDER BY created_at DESC LIMIT 1", (experiment_id,))
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


class CampaignDB:
    """Database operations for experiment campaigns."""

    @staticmethod
    def create(campaign_data: Dict[str, Any]) -> str:
        """Create new campaign."""
        import uuid
        campaign_id = str(uuid.uuid4())

        with get_db() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                INSERT INTO campaigns (
                    id, name, description, status, total_experiments, user_id, created_at
                ) VALUES (?, ?, ?, ?, ?, ?, ?)
            """, (
                campaign_id,
                campaign_data.get('name', 'Unnamed Campaign'),
                campaign_data.get('description', ''),
                'pending',
                campaign_data.get('total_experiments', 0),
                campaign_data.get('user_id'),  # Add user_id
                datetime.utcnow().isoformat()
            ))
            conn.commit()

        return campaign_id

    @staticmethod
    def get(campaign_id: str) -> Optional[Dict[str, Any]]:
        """Get campaign by ID."""
        with get_db() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                SELECT * FROM campaigns WHERE id = ?
            """, (campaign_id,))
            row = cursor.fetchone()

            if not row:
                return None

            return {
                'id': row['id'],
                'name': row['name'],
                'description': row['description'],
                'status': row['status'],
                'total_experiments': row['total_experiments'],
                'completed_experiments': row['completed_experiments'],
                'failed_experiments': row['failed_experiments'],
                'user_id': row['user_id'],  # Add user_id
                'created_at': row['created_at'],
                'started_at': row['started_at'],
                'completed_at': row['completed_at']
            }

    @staticmethod
    def list(limit: int = 100, offset: int = 0) -> List[Dict[str, Any]]:
        """List all campaigns."""
        with get_db() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                SELECT * FROM campaigns
                ORDER BY created_at DESC
                LIMIT ? OFFSET ?
            """, (limit, offset))

            return [
                {
                    'id': row['id'],
                    'name': row['name'],
                    'description': row['description'],
                    'status': row['status'],
                    'total_experiments': row['total_experiments'],
                    'completed_experiments': row['completed_experiments'],
                    'failed_experiments': row['failed_experiments'],
                    'user_id': row['user_id'],  # Add user_id
                    'created_at': row['created_at'],
                    'started_at': row['started_at'],
                    'completed_at': row['completed_at']
                }
                for row in cursor.fetchall()
            ]

    @staticmethod
    def update_status(campaign_id: str, status: str):
        """Update campaign status."""
        with get_db() as conn:
            cursor = conn.cursor()

            if status == 'running':
                cursor.execute("""
                    UPDATE campaigns
                    SET status = ?, started_at = ?
                    WHERE id = ?
                """, (status, datetime.utcnow().isoformat(), campaign_id))
            elif status in ['completed', 'failed', 'cancelled']:
                cursor.execute("""
                    UPDATE campaigns
                    SET status = ?, completed_at = ?
                    WHERE id = ?
                """, (status, datetime.utcnow().isoformat(), campaign_id))
            else:
                cursor.execute("""
                    UPDATE campaigns SET status = ? WHERE id = ?
                """, (status, campaign_id))

            conn.commit()

    @staticmethod
    def update_progress(campaign_id: str, completed: int, failed: int):
        """Update campaign progress."""
        with get_db() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                UPDATE campaigns
                SET completed_experiments = ?, failed_experiments = ?
                WHERE id = ?
            """, (completed, failed, campaign_id))
            conn.commit()

    @staticmethod
    def get_experiments(campaign_id: str) -> List[Dict[str, Any]]:
        """Get all experiments in a campaign."""
        with get_db() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                SELECT * FROM experiments
                WHERE campaign_id = ?
                ORDER BY sequence_order ASC
            """, (campaign_id,))

            return [
                {
                    'id': row['id'],
                    'campaign_id': row['campaign_id'],
                    'status': row['status'],
                    'method': row['method'],
                    'backend': row['backend'],
                    'sequence_order': row['sequence_order'],
                    'created_at': row['created_at'],
                    'started_at': row['started_at'],
                    'completed_at': row['completed_at']
                }
                for row in cursor.fetchall()
            ]


class AnalysisDB:
    """Database operations for analysis results."""

    @staticmethod
    def create(analysis_data: Dict[str, Any]) -> str:
        """Save analysis results."""
        import uuid
        analysis_id = str(uuid.uuid4())

        with get_db() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                INSERT INTO analysis_results (
                    id, experiment_id, profile, analyses, parameters,
                    results, status, computation_time, created_at
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                analysis_id,
                analysis_data['experiment_id'],
                analysis_data['profile'],
                json.dumps(analysis_data.get('analyses', [])),
                json.dumps(analysis_data.get('parameters', {})),
                json.dumps(analysis_data['results']),
                analysis_data.get('status', 'completed'),
                analysis_data.get('computation_time', 0.0),
                datetime.utcnow().isoformat()
            ))
            conn.commit()

        return analysis_id

    @staticmethod
    def get(analysis_id: str) -> Optional[Dict[str, Any]]:
        """Get analysis result by ID."""
        with get_db() as conn:
            cursor = conn.cursor()
            cursor.execute("""
                SELECT * FROM analysis_results WHERE id = ?
            """, (analysis_id,))
            row = cursor.fetchone()

            if not row:
                return None

            return {
                'id': row['id'],
                'experiment_id': row['experiment_id'],
                'profile': row['profile'],
                'analyses': json.loads(row['analyses']),
                'parameters': json.loads(row['parameters']) if row['parameters'] else {},
                'results': json.loads(row['results']),
                'status': row['status'],
                'computation_time': row['computation_time'],
                'created_at': row['created_at']
            }

    @staticmethod
    def get_by_experiment(experiment_id: str, profile: Optional[str] = None) -> List[Dict[str, Any]]:
        """Get all analysis results for an experiment, optionally filtered by profile."""
        with get_db() as conn:
            cursor = conn.cursor()

            if profile:
                cursor.execute("""
                    SELECT * FROM analysis_results
                    WHERE experiment_id = ? AND profile = ?
                    ORDER BY created_at DESC
                """, (experiment_id, profile))
            else:
                cursor.execute("""
                    SELECT * FROM analysis_results
                    WHERE experiment_id = ?
                    ORDER BY created_at DESC
                """, (experiment_id,))

            return [
                {
                    'id': row['id'],
                    'experiment_id': row['experiment_id'],
                    'profile': row['profile'],
                    'analyses': json.loads(row['analyses']),
                    'parameters': json.loads(row['parameters']) if row['parameters'] else {},
                    'results': json.loads(row['results']),
                    'status': row['status'],
                    'computation_time': row['computation_time'],
                    'created_at': row['created_at']
                }
                for row in cursor.fetchall()
            ]

    @staticmethod
    def list(limit: int = 100, offset: int = 0, profile: Optional[str] = None) -> List[Dict[str, Any]]:
        """List analysis results with optional profile filter."""
        with get_db() as conn:
            cursor = conn.cursor()

            if profile:
                cursor.execute("""
                    SELECT * FROM analysis_results
                    WHERE profile = ?
                    ORDER BY created_at DESC
                    LIMIT ? OFFSET ?
                """, (profile, limit, offset))
            else:
                cursor.execute("""
                    SELECT * FROM analysis_results
                    ORDER BY created_at DESC
                    LIMIT ? OFFSET ?
                """, (limit, offset))

            return [
                {
                    'id': row['id'],
                    'experiment_id': row['experiment_id'],
                    'profile': row['profile'],
                    'analyses': json.loads(row['analyses']),
                    'parameters': json.loads(row['parameters']) if row['parameters'] else {},
                    'results': json.loads(row['results']),
                    'status': row['status'],
                    'computation_time': row['computation_time'],
                    'created_at': row['created_at']
                }
                for row in cursor.fetchall()
            ]

    @staticmethod
    def delete(analysis_id: str):
        """Delete analysis result."""
        with get_db() as conn:
            cursor = conn.cursor()
            cursor.execute("DELETE FROM analysis_results WHERE id = ?", (analysis_id,))
            conn.commit()

    @staticmethod
    def delete_by_experiment(experiment_id: str):
        """Delete all analysis results for an experiment."""
        with get_db() as conn:
            cursor = conn.cursor()
            cursor.execute("DELETE FROM analysis_results WHERE experiment_id = ?", (experiment_id,))
            conn.commit()
