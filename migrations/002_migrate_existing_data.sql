-- Migration 002: Migrate Existing Data and Add User Foreign Keys
-- Description: Adds user_id columns to existing tables and preserves current data
-- Date: 2025-10-30

-- Create campaigns table with user_id (if migrating from SQLite)
CREATE TABLE IF NOT EXISTS campaigns (
    id SERIAL PRIMARY KEY,
    user_id INTEGER REFERENCES users(id) ON DELETE SET NULL,
    name VARCHAR(255) NOT NULL,
    description TEXT,
    parameters JSONB,
    status VARCHAR(50) NOT NULL DEFAULT 'active',
    created_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX idx_campaigns_user_id ON campaigns(user_id);
CREATE INDEX idx_campaigns_status ON campaigns(status);
CREATE INDEX idx_campaigns_created_at ON campaigns(created_at);

-- Create experiments table with user_id
CREATE TABLE IF NOT EXISTS experiments (
    id SERIAL PRIMARY KEY,
    user_id INTEGER REFERENCES users(id) ON DELETE SET NULL,
    campaign_id INTEGER REFERENCES campaigns(id) ON DELETE SET NULL,
    name VARCHAR(255) NOT NULL,
    molecule_name VARCHAR(255),
    smiles TEXT,
    xyz_data TEXT,
    atoms JSONB,
    basis_set VARCHAR(50),
    charge INTEGER DEFAULT 0,
    multiplicity INTEGER DEFAULT 1,
    backend VARCHAR(50),
    method VARCHAR(50),
    status experiment_status NOT NULL DEFAULT 'pending',
    progress FLOAT DEFAULT 0.0,
    error_message TEXT,
    created_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    started_at TIMESTAMP,
    completed_at TIMESTAMP
);

CREATE INDEX idx_experiments_user_id ON experiments(user_id);
CREATE INDEX idx_experiments_campaign_id ON experiments(campaign_id);
CREATE INDEX idx_experiments_status ON experiments(status);
CREATE INDEX idx_experiments_created_at ON experiments(created_at);

-- Create jobs table
CREATE TABLE IF NOT EXISTS jobs (
    id SERIAL PRIMARY KEY,
    experiment_id INTEGER NOT NULL REFERENCES experiments(id) ON DELETE CASCADE,
    job_id VARCHAR(255) UNIQUE NOT NULL,
    backend VARCHAR(50),
    status VARCHAR(50) NOT NULL DEFAULT 'pending',
    progress FLOAT DEFAULT 0.0,
    result JSONB,
    error TEXT,
    created_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    started_at TIMESTAMP,
    completed_at TIMESTAMP
);

CREATE INDEX idx_jobs_experiment_id ON jobs(experiment_id);
CREATE INDEX idx_jobs_job_id ON jobs(job_id);
CREATE INDEX idx_jobs_status ON jobs(status);

-- Create analysis_results table
CREATE TABLE IF NOT EXISTS analysis_results (
    id SERIAL PRIMARY KEY,
    experiment_id INTEGER NOT NULL REFERENCES experiments(id) ON DELETE CASCADE,
    analysis_type VARCHAR(100) NOT NULL,
    result_data JSONB,
    file_path TEXT,
    created_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX idx_analysis_results_experiment_id ON analysis_results(experiment_id);
CREATE INDEX idx_analysis_results_type ON analysis_results(analysis_type);

-- Create user_settings table
CREATE TABLE IF NOT EXISTS user_settings (
    id SERIAL PRIMARY KEY,
    user_id INTEGER REFERENCES users(id) ON DELETE CASCADE,
    settings JSONB,
    created_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX idx_user_settings_user_id ON user_settings(user_id);

-- Create cloud_credentials table
CREATE TABLE IF NOT EXISTS cloud_credentials (
    id SERIAL PRIMARY KEY,
    user_id INTEGER REFERENCES users(id) ON DELETE CASCADE,
    provider VARCHAR(50) NOT NULL,
    credentials JSONB NOT NULL,
    is_active BOOLEAN NOT NULL DEFAULT TRUE,
    created_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX idx_cloud_credentials_user_id ON cloud_credentials(user_id);
CREATE INDEX idx_cloud_credentials_provider ON cloud_credentials(provider);

-- Triggers for updated_at
CREATE TRIGGER update_campaigns_updated_at
    BEFORE UPDATE ON campaigns
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_user_settings_updated_at
    BEFORE UPDATE ON user_settings
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_cloud_credentials_updated_at
    BEFORE UPDATE ON cloud_credentials
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

-- Data Migration Strategy:
-- Note: Existing experiments and campaigns will have NULL user_id initially
-- This preserves backward compatibility and allows viewing as "system data"
-- Admin can later assign ownership if needed

-- Comments
COMMENT ON TABLE campaigns IS 'Experiment campaigns organized by users';
COMMENT ON TABLE experiments IS 'Quantum chemistry experiments with results';
COMMENT ON TABLE jobs IS 'Quantum computing jobs from various backends';
COMMENT ON TABLE analysis_results IS 'Analysis results from completed experiments';
COMMENT ON TABLE user_settings IS 'Per-user settings and preferences';
COMMENT ON TABLE cloud_credentials IS 'Cloud provider credentials per user';

COMMENT ON COLUMN campaigns.user_id IS 'NULL = system/shared campaign';
COMMENT ON COLUMN experiments.user_id IS 'NULL = system/shared experiment (existing data)';
COMMENT ON COLUMN user_settings.user_id IS 'NULL = global settings';
COMMENT ON COLUMN cloud_credentials.user_id IS 'NULL = global credentials';

-- Success message
DO $$
BEGIN
    RAISE NOTICE '✓ Migration 002 completed: Existing data tables created with user relationships';
END $$;
