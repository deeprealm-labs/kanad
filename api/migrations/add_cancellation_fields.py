"""
Database migration to add job cancellation support.

Adds the following columns to the experiments table:
- cloud_job_id: Store cloud provider job ID (IBM/BlueQubit)
- cloud_backend: Store cloud backend type ('ibm' or 'bluequbit')
- cancelled_at: Timestamp when experiment was cancelled

Run this migration after updating the Experiment model.
"""

from sqlalchemy import create_engine, text
from api.config import settings

def run_migration():
    """Run the migration to add cancellation fields."""
    engine = create_engine(settings.DATABASE_URL)

    with engine.connect() as conn:
        # Check if columns already exist using SQLite's PRAGMA
        result = conn.execute(text("PRAGMA table_info(experiments)"))
        existing_columns = {row[1] for row in result}  # row[1] is the column name

        # Add cloud_job_id if not exists
        if 'cloud_job_id' not in existing_columns:
            print("Adding cloud_job_id column...")
            conn.execute(text("ALTER TABLE experiments ADD COLUMN cloud_job_id TEXT"))
            conn.commit()
            print("  - cloud_job_id added successfully")
        else:
            print("  - cloud_job_id already exists, skipping")

        # Add cloud_backend if not exists
        if 'cloud_backend' not in existing_columns:
            print("Adding cloud_backend column...")
            conn.execute(text("ALTER TABLE experiments ADD COLUMN cloud_backend TEXT"))
            conn.commit()
            print("  - cloud_backend added successfully")
        else:
            print("  - cloud_backend already exists, skipping")

        # Add cancelled_at if not exists
        if 'cancelled_at' not in existing_columns:
            print("Adding cancelled_at column...")
            conn.execute(text("ALTER TABLE experiments ADD COLUMN cancelled_at TEXT"))
            conn.commit()
            print("  - cancelled_at added successfully")
        else:
            print("  - cancelled_at already exists, skipping")

        print("\nMigration completed successfully!")

if __name__ == "__main__":
    print("Running database migration: add_cancellation_fields")
    print("=" * 60)
    run_migration()
