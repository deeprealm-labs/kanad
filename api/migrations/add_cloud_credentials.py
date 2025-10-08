"""
Migration: Add cloud_credentials table

Adds table for storing encrypted cloud provider credentials.
"""

import sqlite3
from pathlib import Path


def run_migration():
    """Add cloud_credentials table."""

    # Connect to database
    db_path = Path(__file__).parent.parent.parent / "kanad.db"
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()

    # Check if table exists
    cursor.execute("""
        SELECT name FROM sqlite_master
        WHERE type='table' AND name='cloud_credentials'
    """)

    if cursor.fetchone():
        print("cloud_credentials table already exists, skipping migration")
        conn.close()
        return

    # Create cloud_credentials table
    cursor.execute("""
        CREATE TABLE cloud_credentials (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            provider VARCHAR(50) NOT NULL,
            user_id INTEGER DEFAULT 1,
            credentials TEXT NOT NULL,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP NOT NULL,
            updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP NOT NULL
        )
    """)

    # Create indexes
    cursor.execute("""
        CREATE INDEX ix_cloud_credentials_provider
        ON cloud_credentials(provider)
    """)

    cursor.execute("""
        CREATE INDEX ix_cloud_credentials_user_id
        ON cloud_credentials(user_id)
    """)

    # Create unique constraint on provider + user_id
    cursor.execute("""
        CREATE UNIQUE INDEX ix_cloud_credentials_provider_user
        ON cloud_credentials(provider, user_id)
    """)

    conn.commit()
    conn.close()

    print("✓ Added cloud_credentials table")


if __name__ == '__main__':
    run_migration()
