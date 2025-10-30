-- Migration 001: Create Authentication Tables
-- Description: Creates users, access_keys, sessions, and email_verifications tables
-- Date: 2025-10-30

-- User roles enum
CREATE TYPE user_role AS ENUM ('admin', 'user', 'viewer');

-- Experiment status enum
CREATE TYPE experiment_status AS ENUM ('pending', 'running', 'completed', 'failed', 'cancelled');

-- Access Keys Table
CREATE TABLE IF NOT EXISTS access_keys (
    id SERIAL PRIMARY KEY,
    key VARCHAR(64) UNIQUE NOT NULL,
    description TEXT,
    max_uses INTEGER NOT NULL DEFAULT 1,
    used_count INTEGER NOT NULL DEFAULT 0,
    expires_at TIMESTAMP,
    is_active BOOLEAN NOT NULL DEFAULT TRUE,
    created_by INTEGER,
    created_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,

    -- Indexes
    CONSTRAINT access_keys_key_idx UNIQUE (key)
);

CREATE INDEX idx_access_keys_active ON access_keys(is_active) WHERE is_active = TRUE;
CREATE INDEX idx_access_keys_expires ON access_keys(expires_at) WHERE expires_at IS NOT NULL;

-- Users Table
CREATE TABLE IF NOT EXISTS users (
    id SERIAL PRIMARY KEY,
    email VARCHAR(255) UNIQUE NOT NULL,
    password_hash VARCHAR(255) NOT NULL,
    full_name VARCHAR(255),
    role user_role NOT NULL DEFAULT 'user',
    is_verified BOOLEAN NOT NULL DEFAULT FALSE,
    is_active BOOLEAN NOT NULL DEFAULT TRUE,
    access_key_id INTEGER REFERENCES access_keys(id) ON DELETE SET NULL,
    created_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    last_login TIMESTAMP,

    -- OAuth fields
    google_id VARCHAR(255) UNIQUE,
    avatar_url TEXT,

    -- Indexes
    CONSTRAINT users_email_idx UNIQUE (email)
);

CREATE INDEX idx_users_email ON users(email);
CREATE INDEX idx_users_google_id ON users(google_id) WHERE google_id IS NOT NULL;
CREATE INDEX idx_users_role ON users(role);
CREATE INDEX idx_users_active ON users(is_active) WHERE is_active = TRUE;

-- Add foreign key constraint for created_by in access_keys
ALTER TABLE access_keys ADD CONSTRAINT fk_access_keys_created_by
    FOREIGN KEY (created_by) REFERENCES users(id) ON DELETE SET NULL;

-- Sessions Table
CREATE TABLE IF NOT EXISTS sessions (
    id SERIAL PRIMARY KEY,
    user_id INTEGER NOT NULL REFERENCES users(id) ON DELETE CASCADE,
    access_token VARCHAR(512) UNIQUE NOT NULL,
    refresh_token VARCHAR(512) UNIQUE NOT NULL,
    expires_at TIMESTAMP NOT NULL,
    refresh_expires_at TIMESTAMP NOT NULL,
    created_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    last_used TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    user_agent TEXT,
    ip_address VARCHAR(45),
    is_revoked BOOLEAN NOT NULL DEFAULT FALSE,

    -- Indexes
    CONSTRAINT sessions_access_token_idx UNIQUE (access_token),
    CONSTRAINT sessions_refresh_token_idx UNIQUE (refresh_token)
);

CREATE INDEX idx_sessions_user_id ON sessions(user_id);
CREATE INDEX idx_sessions_access_token ON sessions(access_token);
CREATE INDEX idx_sessions_refresh_token ON sessions(refresh_token);
CREATE INDEX idx_sessions_expires ON sessions(expires_at);
CREATE INDEX idx_sessions_active ON sessions(is_revoked) WHERE is_revoked = FALSE;

-- Email Verifications Table
CREATE TABLE IF NOT EXISTS email_verifications (
    id SERIAL PRIMARY KEY,
    email VARCHAR(255) NOT NULL,
    otp VARCHAR(6) NOT NULL,
    expires_at TIMESTAMP NOT NULL,
    is_used BOOLEAN NOT NULL DEFAULT FALSE,
    created_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX idx_email_verifications_email ON email_verifications(email);
CREATE INDEX idx_email_verifications_otp ON email_verifications(otp);
CREATE INDEX idx_email_verifications_expires ON email_verifications(expires_at);
CREATE INDEX idx_email_verifications_unused ON email_verifications(is_used) WHERE is_used = FALSE;

-- Trigger to update updated_at timestamp
CREATE OR REPLACE FUNCTION update_updated_at_column()
RETURNS TRIGGER AS $$
BEGIN
    NEW.updated_at = CURRENT_TIMESTAMP;
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;

CREATE TRIGGER update_users_updated_at
    BEFORE UPDATE ON users
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

-- Comments for documentation
COMMENT ON TABLE users IS 'User accounts with authentication and role-based access';
COMMENT ON TABLE access_keys IS 'Early access keys for user registration';
COMMENT ON TABLE sessions IS 'User sessions with JWT tokens';
COMMENT ON TABLE email_verifications IS 'Email verification OTP codes';

COMMENT ON COLUMN users.role IS 'User role: admin (full access), user (experiments), viewer (read-only)';
COMMENT ON COLUMN users.is_verified IS 'Email verification status';
COMMENT ON COLUMN access_keys.max_uses IS 'Maximum number of users allowed to use this key';
COMMENT ON COLUMN access_keys.used_count IS 'Number of times this key has been used';

-- Success message
DO $$
BEGIN
    RAISE NOTICE '✓ Migration 001 completed: Authentication tables created';
END $$;
