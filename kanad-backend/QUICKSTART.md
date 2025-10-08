# Kanad Backend - Quick Start Guide

Get your Kanad backend running in **5 minutes**!

## Prerequisites

- Python 3.11+
- PostgreSQL 15+
- Redis 7+

Or just use **Docker** (easiest).

---

## Option 1: Docker (Recommended)

### Step 1: Navigate to backend directory
```bash
cd /home/mk/deeprealm/kanad/kanad-backend
```

### Step 2: Create .env file
```bash
cp .env.example .env
```

### Step 3: Generate security keys
```bash
# Generate JWT secret
python3 -c "import secrets; print('JWT_SECRET_KEY=' + secrets.token_urlsafe(32))" >> .env

# Generate encryption key
python3 -c "from cryptography.fernet import Fernet; print('ENCRYPTION_KEY=' + Fernet.generate_key().decode())" >> .env
```

### Step 4: Start all services
```bash
docker-compose up --build
```

### Step 5: Access the API
- **API Docs**: http://localhost:8000/api/docs
- **Health Check**: http://localhost:8000/health

**That's it!** The backend is running with PostgreSQL, Redis, FastAPI, and Celery workers.

---

## Option 2: Manual Setup (Development)

### Step 1: Install system dependencies (Ubuntu/Debian)
```bash
sudo apt update
sudo apt install -y python3.11 python3.11-venv postgresql redis-server
sudo systemctl start postgresql redis-server
```

### Step 2: Setup database
```bash
sudo -u postgres psql
```

In PostgreSQL shell:
```sql
CREATE DATABASE kanad_experiments;
CREATE USER kanad_user WITH PASSWORD 'kanad_password';
GRANT ALL PRIVILEGES ON DATABASE kanad_experiments TO kanad_user;
\q
```

### Step 3: Create Python environment
```bash
cd /home/mk/deeprealm/kanad/kanad-backend
python3.11 -m venv venv
source venv/bin/activate
```

### Step 4: Install dependencies
```bash
# Install Kanad framework first
cd ..
pip install -e .

# Install backend dependencies
cd kanad-backend
pip install -r requirements.txt
```

### Step 5: Configure environment
```bash
cp .env.example .env

# Generate keys
python -c "import secrets; print('JWT_SECRET_KEY=' + secrets.token_urlsafe(32))" >> .env
python -c "from cryptography.fernet import Fernet; print('ENCRYPTION_KEY=' + Fernet.generate_key().decode())" >> .env
```

### Step 6: Initialize database
```bash
python -c "from core.database import init_db; init_db()"
```

### Step 7: Start services

**Terminal 1 - FastAPI**:
```bash
uvicorn api.main:app --reload --host 0.0.0.0 --port 8000
```

**Terminal 2 - Celery Worker**:
```bash
celery -A workers.celery_app worker --loglevel=info --concurrency=4
```

### Step 8: Test the API
```bash
curl http://localhost:8000/health
# Should return: {"status": "healthy", ...}
```

---

## Quick API Test

### 1. Register a user
```bash
curl -X POST http://localhost:8000/api/auth/register \
  -H "Content-Type: application/json" \
  -d '{
    "email": "test@example.com",
    "password": "secure_password123",
    "name": "Test User",
    "institution": "University",
    "field": "chemistry"
  }'
```

**Response**:
```json
{
  "access_token": "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9...",
  "token_type": "bearer"
}
```

**Save the token** - you'll need it for authenticated requests.

### 2. Create a molecule (H₂)
```bash
export TOKEN="<your_access_token_from_step_1>"

curl -X POST http://localhost:8000/api/molecules/create \
  -H "Content-Type: application/json" \
  -H "Authorization: Bearer $TOKEN" \
  -d '{
    "method": "atoms",
    "data": {
      "atoms": [
        {"element": "H", "position": [0, 0, 0]},
        {"element": "H", "position": [0, 0, 0.74]}
      ],
      "basis": "sto-3g",
      "charge": 0,
      "multiplicity": 1,
      "name": "Hydrogen Molecule"
    }
  }'
```

**Response**:
```json
{
  "molecule_id": "a1b2c3d4-...",
  "name": "Hydrogen Molecule",
  "formula": "H2",
  "n_electrons": 2,
  "n_orbitals": 2,
  "n_qubits": 4,
  ...
}
```

### 3. Create a molecule from SMILES (Water)
```bash
curl -X POST http://localhost:8000/api/molecules/from-smiles \
  -H "Content-Type: application/json" \
  -H "Authorization: Bearer $TOKEN" \
  -d '{
    "smiles": "O",
    "basis": "sto-3g",
    "optimize_geometry": false
  }'
```

### 4. List your molecules
```bash
curl -X GET http://localhost:8000/api/molecules/ \
  -H "Authorization: Bearer $TOKEN"
```

### 5. Check API documentation
Open browser: http://localhost:8000/api/docs

Interactive Swagger UI with all endpoints!

---

## Testing Computation

### Submit a VQE calculation
```bash
# First, get your molecule_id from step 2 or 3

curl -X POST http://localhost:8000/api/simulations/configure \
  -H "Content-Type: application/json" \
  -H "Authorization: Bearer $TOKEN" \
  -d '{
    "molecule_id": "<your_molecule_id>",
    "method": "VQE",
    "ansatz": "hardware_efficient",
    "mapper": "jordan_wigner",
    "optimizer": "SLSQP",
    "max_iterations": 100,
    "backend": {
      "type": "classical"
    },
    "analysis": {
      "energy_decomposition": true,
      "bond_analysis": true,
      "dipole_moment": true
    }
  }'
```

This returns a simulation preview. Then accept and run:

```bash
curl -X POST http://localhost:8000/api/simulations/<simulation_id>/accept-and-run \
  -H "Content-Type: application/json" \
  -H "Authorization: Bearer $TOKEN" \
  -d '{
    "accepted_terms": true
  }'
```

Returns a `job_id`. Check status:

```bash
curl -X GET http://localhost:8000/api/jobs/<job_id>/status \
  -H "Authorization: Bearer $TOKEN"
```

---

## Troubleshooting

### Database connection error
```bash
# Check PostgreSQL is running
sudo systemctl status postgresql

# Check credentials in .env match database
cat .env | grep DATABASE_URL
```

### Redis connection error
```bash
# Check Redis is running
sudo systemctl status redis-server

# Test Redis
redis-cli ping  # Should return "PONG"
```

### Celery worker not starting
```bash
# Check Redis connection
# Make sure DATABASE_URL and REDIS_URL are set in .env
# Try running with verbose logging
celery -A workers.celery_app worker --loglevel=debug
```

### Import errors
```bash
# Make sure Kanad is installed
pip list | grep kanad

# If not found, install from parent directory
cd /home/mk/deeprealm/kanad
pip install -e .
```

### Port already in use
```bash
# Check what's using port 8000
sudo lsof -i :8000

# Kill the process or use different port
uvicorn api.main:app --port 8001
```

---

## Next Steps

1. ✅ **Explore API docs**: http://localhost:8000/api/docs
2. ✅ **Test all endpoints**: Use Swagger UI or Postman
3. ✅ **Monitor logs**: Check terminal outputs for errors
4. ✅ **Build frontend**: Connect Next.js app to API
5. ✅ **Deploy to cloud**: Follow README.md for Azure deployment

---

## Useful Commands

### Check service status
```bash
# FastAPI
curl http://localhost:8000/health

# Celery workers
celery -A workers.celery_app inspect active

# Database
psql -U kanad_user -d kanad_experiments -c "SELECT COUNT(*) FROM users;"
```

### View logs
```bash
# FastAPI (if using systemd)
sudo journalctl -u kanad-api -f

# Celery
sudo journalctl -u kanad-worker -f

# Docker
docker-compose logs -f api
docker-compose logs -f worker
```

### Reset database
```bash
python -c "from core.database import drop_db, init_db; drop_db(); init_db()"
```

**WARNING**: This deletes ALL data!

---

## Support

- **Documentation**: `/api/docs` endpoint
- **README**: See README.md for detailed info
- **Issues**: Check logs first, then file issue

---

**Ready to compute! 🚀⚛️**
