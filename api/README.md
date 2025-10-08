# Kanad Quantum Chemistry API

FastAPI backend for the Kanad web application, providing quantum chemistry calculations through a RESTful API.

## Features

- **Quantum Methods**: VQE, HF, SQD calculations
- **Multiple Ansatze**: UCC, Hardware-Efficient, Governance-aware
- **Backend Support**: Classical simulation, IBM Quantum, BlueQubit
- **SMILES Parsing**: Automatic molecule creation from SMILES strings
- **Job Queue**: Priority-based background execution
- **Real-time Progress**: Convergence tracking and status updates
- **Persistent Storage**: SQLite database for experiments and results

## Architecture

```
/api/
├── main.py                 # FastAPI application entry point
├── config.py              # Configuration and environment variables
├── database.py            # Database connection and session management
├── models/                # SQLAlchemy database models
│   ├── experiment.py
│   ├── queue.py
│   └── settings.py
├── routers/               # API endpoints
│   ├── experiments.py     # Experiment CRUD operations
│   ├── queue.py           # Job queue management
│   ├── molecules.py       # SMILES validation and library
│   └── settings.py        # User settings
├── services/              # Business logic
│   ├── experiment_service.py  # Kanad framework integration
│   └── job_queue.py       # Background job execution
├── utils/                 # Utilities
│   ├── validators.py      # Pydantic request models
│   └── exceptions.py      # Custom exceptions
└── tests/                 # Test scripts
    └── test_api.py        # API endpoint tests
```

## Installation

### Prerequisites

- Python 3.9+
- Kanad quantum chemistry framework
- RDKit (for SMILES parsing)

### Setup

1. **Install dependencies**:
```bash
cd /home/mk/deeprealm/kanad
pip install -r requirements.txt
pip install fastapi uvicorn sqlalchemy pydantic-settings
```

2. **Create environment file** (optional):
```bash
cd /home/mk/deeprealm/kanad/api
cat > .env << EOF
DATABASE_URL=sqlite:///./kanad.db
MAX_CONCURRENT_JOBS=2
DEBUG=True
IBM_QUANTUM_TOKEN=your_token_here
BLUEQUBIT_API_KEY=your_key_here
EOF
```

3. **Start the server**:
```bash
cd /home/mk/deeprealm/kanad/api
python main.py
```

Or use uvicorn directly:
```bash
uvicorn api.main:app --reload --host 0.0.0.0 --port 8000
```

4. **Access the API**:
- API Documentation: http://localhost:8000/docs
- Alternative Docs: http://localhost:8000/redoc
- Health Check: http://localhost:8000/health

## API Endpoints

### Core Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/` | GET | API information |
| `/health` | GET | Health check |
| `/api/v1/info` | GET | API capabilities |

### Experiments

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/v1/experiments/` | POST | Create experiment |
| `/api/v1/experiments/` | GET | List experiments |
| `/api/v1/experiments/{id}` | GET | Get experiment details |
| `/api/v1/experiments/{id}` | DELETE | Delete experiment |
| `/api/v1/experiments/{id}/status` | GET | Get experiment status |
| `/api/v1/experiments/{id}/convergence` | GET | Get convergence data |

### Queue

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/v1/queue/` | POST | Add to queue |
| `/api/v1/queue/` | GET | List queue |
| `/api/v1/queue/{id}` | GET | Get queue item |
| `/api/v1/queue/{id}` | PUT | Update queue item |
| `/api/v1/queue/{id}` | DELETE | Delete queue item |
| `/api/v1/queue/{id}/execute` | POST | Execute queue item |

### Molecules

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/v1/molecules/validate` | POST | Validate SMILES |
| `/api/v1/molecules/library` | GET | Get molecule library |
| `/api/v1/molecules/library/categories` | GET | Get categories |
| `/api/v1/molecules/library/{id}` | GET | Get specific molecule |

### Settings

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/v1/settings/` | GET | Get settings |
| `/api/v1/settings/` | PUT | Update settings |
| `/api/v1/settings/` | DELETE | Reset settings |

## Usage Examples

### 1. Create and Execute VQE Experiment

```bash
curl -X POST http://localhost:8000/api/v1/experiments/ \
  -H "Content-Type: application/json" \
  -d '{
    "name": "H2 VQE Calculation",
    "molecule": {
      "smiles": "[H][H]",
      "basis": "sto-3g",
      "charge": 0,
      "multiplicity": 1
    },
    "configuration": {
      "method": "VQE",
      "ansatz": "ucc",
      "mapper": "jordan_wigner",
      "optimizer": "SLSQP",
      "backend": "classical",
      "max_iterations": 1000,
      "conv_threshold": 1e-6
    },
    "execute_immediately": true
  }'
```

### 2. Check Experiment Status

```bash
curl http://localhost:8000/api/v1/experiments/1/status
```

### 3. Get Convergence Data

```bash
curl http://localhost:8000/api/v1/experiments/1/convergence
```

### 4. Validate SMILES

```bash
curl -X POST http://localhost:8000/api/v1/molecules/validate \
  -H "Content-Type: application/json" \
  -d '{"smiles": "CCO"}'
```

### 5. Get Molecule Library

```bash
curl http://localhost:8000/api/v1/molecules/library
```

### 6. Update Settings

```bash
curl -X PUT http://localhost:8000/api/v1/settings/ \
  -H "Content-Type: application/json" \
  -d '{
    "method": "VQE",
    "ansatz": "hardware_efficient",
    "optimizer": "COBYLA",
    "backend": "classical"
  }'
```

## Python Client Example

```python
import requests
import time

# Create experiment
response = requests.post(
    "http://localhost:8000/api/v1/experiments/",
    json={
        "name": "Ethanol VQE",
        "molecule": {
            "smiles": "CCO",
            "basis": "sto-3g",
            "charge": 0,
            "multiplicity": 1
        },
        "configuration": {
            "method": "VQE",
            "ansatz": "ucc",
            "mapper": "jordan_wigner",
            "optimizer": "SLSQP",
            "backend": "classical"
        },
        "execute_immediately": True
    }
)

experiment_id = response.json()["id"]
print(f"Created experiment {experiment_id}")

# Poll for completion
while True:
    status_response = requests.get(
        f"http://localhost:8000/api/v1/experiments/{experiment_id}/status"
    )
    status = status_response.json()["status"]
    progress = status_response.json()["progress"]

    print(f"Status: {status}, Progress: {progress}%")

    if status in ["completed", "failed"]:
        break

    time.sleep(2)

# Get results
result = requests.get(
    f"http://localhost:8000/api/v1/experiments/{experiment_id}"
).json()

print(f"Energy: {result['energy']} Hartree")
print(f"Converged: {result['results']['converged']}")
```

## Kanad Framework Capabilities

### Supported Methods

- **VQE**: Variational Quantum Eigensolver
- **HF**: Hartree-Fock
- **SQD**: Subspace Quantum Diagonalization

### Ansatze

- **UCC**: Unitary Coupled Cluster (singles + doubles)
- **Hardware-Efficient**: Parameterized quantum circuits
- **Governance-Aware**: Bond-type specific ansatze

### Mappers

- **Jordan-Wigner**: Standard fermionic-to-qubit mapping
- **Bravyi-Kitaev**: More efficient qubit mapping
- **Hybrid Orbital**: Advanced orbital-based mapping

### Optimizers

- **SLSQP**: Sequential Least Squares Programming
- **COBYLA**: Constrained Optimization BY Linear Approximation
- **L-BFGS-B**: Limited-memory BFGS
- **ADAM**: Adaptive Moment Estimation
- **POWELL**: Powell's method

### Backends

- **Classical**: Local statevector simulation
- **IBM Quantum**: IBM Quantum cloud backends
- **BlueQubit**: GPU-accelerated cloud backend

### Basis Sets

- **STO-3G**: Minimal basis
- **6-31G**: Split-valence basis
- **6-31G***, **6-31G***: Polarized basis
- **cc-pVDZ**, **cc-pVTZ**: Correlation-consistent basis

## Testing

Run the test suite:

```bash
# Start the server first
python /home/mk/deeprealm/kanad/api/main.py

# In another terminal, run tests
python /home/mk/deeprealm/kanad/api/tests/test_api.py
```

## Configuration

### Environment Variables

- `DATABASE_URL`: Database connection string (default: sqlite:///./kanad.db)
- `MAX_CONCURRENT_JOBS`: Max concurrent experiments (default: 2)
- `DEBUG`: Enable debug mode (default: True)
- `IBM_QUANTUM_TOKEN`: IBM Quantum API token
- `BLUEQUBIT_API_KEY`: BlueQubit API key
- `CORS_ORIGINS`: Allowed CORS origins (comma-separated)

### Database

The API uses SQLite by default for simplicity. For production, use PostgreSQL:

```bash
DATABASE_URL=postgresql://user:password@localhost:5432/kanad
```

## Troubleshooting

### Server won't start

1. Check Python version: `python --version` (need 3.9+)
2. Install dependencies: `pip install -r requirements.txt`
3. Check port availability: `lsof -i :8000`

### Experiments fail

1. Check logs in terminal where server is running
2. Verify SMILES string is valid: `/api/v1/molecules/validate`
3. Check basis set is supported
4. Ensure RDKit is installed: `pip install rdkit`

### Queue not processing

1. Check health endpoint: `/health`
2. Verify `job_queue_running: true`
3. Check server logs for worker errors

## Frontend Integration

Update Next.js frontend to use this API:

```typescript
// /web/src/lib/api.ts
const API_BASE = process.env.NEXT_PUBLIC_API_URL || 'http://localhost:8000/api/v1';

export async function submitExperiment(config: any) {
  const response = await fetch(`${API_BASE}/experiments/`, {
    method: 'POST',
    headers: {'Content-Type': 'application/json'},
    body: JSON.stringify(config)
  });
  return response.json();
}
```

## Production Deployment

### Docker

```dockerfile
FROM python:3.11-slim

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY api/ ./api/
COPY kanad/ ./kanad/

CMD ["uvicorn", "api.main:app", "--host", "0.0.0.0", "--port", "8000"]
```

### Kubernetes

See deployment manifests in `/deployments/kubernetes/`

## License

See main Kanad project license.

## Support

For issues or questions:
- GitHub Issues: https://github.com/yourusername/kanad/issues
- Documentation: http://localhost:8000/docs
