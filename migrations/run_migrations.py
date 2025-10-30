#!/usr/bin/env python3
"""
Database Migration Runner for Kanad Platform
Executes SQL migrations in order and tracks migration history
"""

import os
import sys
import psycopg2
from psycopg2.extensions import ISOLATION_LEVEL_AUTOCOMMIT
from datetime import datetime
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))


def get_db_connection(database_url: str = None):
    """Create PostgreSQL connection"""
    if database_url is None:
        database_url = os.getenv("DATABASE_URL")

    if not database_url:
        raise ValueError("DATABASE_URL environment variable not set")

    # Parse connection string
    # Format: postgresql://user:password@host:port/database
    import re
    match = re.match(r'postgresql://([^:]+):([^@]+)@([^:]+):(\d+)/(.+)', database_url)
    if not match:
        raise ValueError(f"Invalid DATABASE_URL format: {database_url}")

    user, password, host, port, database = match.groups()

    return psycopg2.connect(
        host=host,
        port=int(port),
        database=database,
        user=user,
        password=password
    )


def create_migration_table(conn):
    """Create migrations tracking table if it doesn't exist"""
    cursor = conn.cursor()
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS schema_migrations (
            id SERIAL PRIMARY KEY,
            migration_name VARCHAR(255) UNIQUE NOT NULL,
            executed_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
            success BOOLEAN NOT NULL DEFAULT TRUE,
            error_message TEXT
        );
    """)
    conn.commit()
    cursor.close()
    print("✓ Migration tracking table ready")


def get_executed_migrations(conn):
    """Get list of already executed migrations"""
    cursor = conn.cursor()
    cursor.execute("SELECT migration_name FROM schema_migrations WHERE success = TRUE ORDER BY id")
    executed = [row[0] for row in cursor.fetchall()]
    cursor.close()
    return executed


def execute_migration(conn, migration_file: Path):
    """Execute a single migration file"""
    migration_name = migration_file.name

    print(f"\n{'='*60}")
    print(f"Executing migration: {migration_name}")
    print(f"{'='*60}")

    cursor = conn.cursor()

    try:
        # Read migration SQL
        with open(migration_file, 'r') as f:
            sql = f.read()

        # Execute migration
        cursor.execute(sql)
        conn.commit()

        # Record success
        cursor.execute(
            "INSERT INTO schema_migrations (migration_name, success) VALUES (%s, %s)",
            (migration_name, True)
        )
        conn.commit()

        print(f"✓ Migration {migration_name} completed successfully")
        return True

    except Exception as e:
        conn.rollback()

        # Record failure
        cursor.execute(
            "INSERT INTO schema_migrations (migration_name, success, error_message) VALUES (%s, %s, %s)",
            (migration_name, False, str(e))
        )
        conn.commit()

        print(f"✗ Migration {migration_name} failed:")
        print(f"  Error: {str(e)}")
        return False

    finally:
        cursor.close()


def run_migrations(database_url: str = None):
    """Run all pending migrations"""
    print("\n" + "="*60)
    print("Kanad Platform - Database Migration Runner")
    print("="*60)

    # Connect to database
    try:
        conn = get_db_connection(database_url)
        print(f"✓ Connected to PostgreSQL database")
    except Exception as e:
        print(f"✗ Failed to connect to database:")
        print(f"  Error: {str(e)}")
        print(f"\nPlease ensure:")
        print(f"  1. PostgreSQL is running")
        print(f"  2. DATABASE_URL environment variable is set")
        print(f"  3. Database exists and credentials are correct")
        sys.exit(1)

    try:
        # Create migration tracking table
        create_migration_table(conn)

        # Get executed migrations
        executed = get_executed_migrations(conn)
        print(f"✓ Found {len(executed)} previously executed migrations")

        # Find all migration files
        migrations_dir = Path(__file__).parent
        migration_files = sorted(migrations_dir.glob("*.sql"))

        if not migration_files:
            print("✓ No migration files found")
            return

        print(f"✓ Found {len(migration_files)} migration files")

        # Execute pending migrations
        pending = [f for f in migration_files if f.name not in executed]

        if not pending:
            print("\n✓ All migrations are up to date!")
            return

        print(f"\n→ Executing {len(pending)} pending migrations...")

        success_count = 0
        for migration_file in pending:
            if execute_migration(conn, migration_file):
                success_count += 1
            else:
                print(f"\n✗ Migration failed, stopping execution")
                break

        print(f"\n{'='*60}")
        print(f"Migration Summary:")
        print(f"  Total: {len(pending)}")
        print(f"  Success: {success_count}")
        print(f"  Failed: {len(pending) - success_count}")
        print(f"{'='*60}\n")

        if success_count == len(pending):
            print("✓ All migrations completed successfully!")
        else:
            print("✗ Some migrations failed. Please check errors above.")
            sys.exit(1)

    finally:
        conn.close()


def rollback_migration(migration_name: str, database_url: str = None):
    """Rollback a specific migration (manual intervention required)"""
    print(f"\n⚠️  Rolling back migration: {migration_name}")
    print(f"⚠️  Note: Automatic rollback is not implemented.")
    print(f"⚠️  You need to manually write and execute the rollback SQL.")

    conn = get_db_connection(database_url)
    cursor = conn.cursor()

    try:
        # Mark as not executed
        cursor.execute(
            "DELETE FROM schema_migrations WHERE migration_name = %s",
            (migration_name,)
        )
        conn.commit()
        print(f"✓ Migration {migration_name} marked as not executed")

    except Exception as e:
        conn.rollback()
        print(f"✗ Rollback failed: {str(e)}")

    finally:
        cursor.close()
        conn.close()


def check_migration_status(database_url: str = None):
    """Check current migration status"""
    print("\n" + "="*60)
    print("Migration Status")
    print("="*60)

    try:
        conn = get_db_connection(database_url)
        cursor = conn.cursor()

        cursor.execute("""
            SELECT migration_name, executed_at, success, error_message
            FROM schema_migrations
            ORDER BY id
        """)

        migrations = cursor.fetchall()

        if not migrations:
            print("No migrations executed yet")
        else:
            print(f"\nExecuted migrations: {len(migrations)}\n")
            for name, executed_at, success, error in migrations:
                status = "✓" if success else "✗"
                print(f"{status} {name}")
                print(f"   Executed: {executed_at}")
                if not success:
                    print(f"   Error: {error}")
                print()

        cursor.close()
        conn.close()

    except Exception as e:
        print(f"✗ Failed to check status: {str(e)}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Kanad Database Migration Tool")
    parser.add_argument("command", choices=["migrate", "status", "rollback"], help="Command to execute")
    parser.add_argument("--migration", help="Migration name (for rollback)")
    parser.add_argument("--database-url", help="PostgreSQL connection string (overrides DATABASE_URL env)")

    args = parser.parse_args()

    if args.command == "migrate":
        run_migrations(args.database_url)
    elif args.command == "status":
        check_migration_status(args.database_url)
    elif args.command == "rollback":
        if not args.migration:
            print("✗ --migration argument required for rollback")
            sys.exit(1)
        rollback_migration(args.migration, args.database_url)
