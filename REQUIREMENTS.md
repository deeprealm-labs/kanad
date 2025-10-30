# Kanad Dependencies Guide

This document explains how to install dependencies for the Kanad quantum chemistry platform.

## Requirements Structure

Kanad uses a modular requirements structure to separate concerns:

```
kanad/
├── requirements.txt        # Core quantum chemistry & scientific computing
└── api/
    └── requirements.txt    # Backend API server & authentication
```

## Installation

### Quick Start - Full Installation

For a complete development setup with all features:

```bash
# Install core dependencies
pip install -r requirements.txt

# Install API server dependencies
pip install -r api/requirements.txt
```

### Minimal Installation Options

#### 1. Core Quantum Chemistry Only

If you only need the quantum chemistry computation engine without the web API:

```bash
pip install -r requirements.txt
```

**Includes:**
- Quantum computing frameworks (Qiskit 2.x)
- Quantum simulators (Aer, BlueQubit)
- Chemistry packages (PySCF, RDKit)
- Scientific computing (NumPy, SciPy)
- Visualization (Matplotlib, Plotly)

#### 2. API Server Only

If you're running just the backend API (assuming core dependencies are already installed):

```bash
pip install -r api/requirements.txt
```

**Includes:**
- FastAPI web framework
- Authentication & security (JWT, bcrypt, OAuth)
- Database support (PostgreSQL, SQLite)
- Email support (SMTP)
- API utilities

## Detailed Dependency Information

### Core Dependencies (`requirements.txt`)

| Category | Packages | Purpose |
|----------|----------|---------|
| **Scientific Computing** | numpy, scipy | Numerical operations |
| **Chemistry** | rdkit, pyscf | Molecular manipulation & quantum chemistry |
| **Quantum Computing** | qiskit, qiskit-aer, qiskit-ibm-runtime, qiskit-addon-sqd | Quantum circuit simulation & execution |
| **Cloud Backends** | bluequbit | Cloud quantum backend integration |
| **Visualization** | matplotlib, plotly | 2D/3D molecular and data visualization |
| **Testing** | pytest, pytest-cov | Unit testing and coverage |
| **Utilities** | dataclasses-json, typing-extensions | Data serialization & type hints |

### API Dependencies (`api/requirements.txt`)

| Category | Packages | Purpose |
|----------|----------|---------|
| **Web Framework** | fastapi, uvicorn, python-multipart | REST API server & WebSocket support |
| **Validation** | pydantic, email-validator | Request/response validation |
| **Database** | psycopg2-binary, sqlalchemy, alembic | PostgreSQL ORM & migrations |
| **Authentication** | pyjwt, bcrypt, python-jose, passlib | JWT tokens & password hashing |
| **OAuth** | requests, httpx | Google OAuth integration |
| **Email** | aiosmtplib | Async email sending |
| **Configuration** | python-dotenv | Environment variable management |

## Optional GPU Acceleration

For NVIDIA GPU support (requires CUDA), uncomment these lines in `requirements.txt`:

```bash
# Uncomment in requirements.txt:
cuquantum-python>=23.10.0
cupy-cuda11x>=12.0.0     # For CUDA 11.x
# OR
cupy-cuda12x>=12.0.0     # For CUDA 12.x
qiskit-aer-gpu>=0.13.0
```

## Database Setup

### SQLite (Automatic)
SQLite is used for experiments/jobs storage and requires no additional setup.

### PostgreSQL (Required for Authentication)

1. Install PostgreSQL:
```bash
# Ubuntu/Debian
sudo apt-get install postgresql postgresql-contrib

# macOS
brew install postgresql
```

2. Create database:
```bash
sudo -u postgres psql
CREATE DATABASE kanad_db;
CREATE USER kanad_user WITH PASSWORD 'your_password';
GRANT ALL PRIVILEGES ON DATABASE kanad_db TO kanad_user;
```

3. Configure environment variables in `.env`:
```bash
DATABASE_URL="postgresql://kanad_user:your_password@localhost:5432/kanad_db"
```

## Verification

After installation, verify your setup:

```bash
# Test core imports
python -c "import qiskit; import pyscf; import rdkit; print('Core: OK')"

# Test API imports
python -c "import fastapi; import sqlalchemy; import jwt; print('API: OK')"

# Run tests
pytest kanad/tests/
```

## Troubleshooting

### RDKit Installation Issues
If RDKit fails to install via pip, use conda:
```bash
conda install -c conda-forge rdkit
```

### PostgreSQL Connection Issues
Ensure PostgreSQL is running:
```bash
sudo systemctl status postgresql    # Linux
brew services list | grep postgresql # macOS
```

### GPU Support
Verify CUDA is installed:
```bash
nvidia-smi
nvcc --version
```

## Development Setup

For active development, install in editable mode:

```bash
# Install core in editable mode
pip install -e .

# Install with all extras
pip install -r requirements.txt -r api/requirements.txt
```

## Production Deployment

For production, pin all versions to ensure reproducibility:

```bash
# Generate pinned requirements
pip freeze > requirements-lock.txt

# Install from pinned file
pip install -r requirements-lock.txt
```

## Version Compatibility

- **Python**: 3.10 or higher required
- **Qiskit**: 2.2.0+ (Qiskit 2.x series)
- **FastAPI**: 0.115.5+
- **PostgreSQL**: 12.0+ recommended

## Support

For dependency issues or questions:
- GitHub Issues: https://github.com/yourusername/kanad/issues
- Documentation: See project README.md
