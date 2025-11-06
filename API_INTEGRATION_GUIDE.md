# Kanad Framework - API Integration Guide

**Version**: 1.0.0
**Date**: November 7, 2025
**Status**: Production Ready

---

## Overview

This guide provides complete instructions for integrating the Kanad quantum chemistry framework with a backend API (FastAPI) and web application (Next.js).

---

## Architecture

```
┌─────────────────┐
│  Web Frontend   │  Next.js + React
│  (TypeScript)   │
└────────┬────────┘
         │ HTTP/REST
         ▼
┌─────────────────┐
│  Backend API    │  FastAPI + Python
│  (api/main.py)  │
└────────┬────────┘
         │ Direct Import
         ▼
┌─────────────────┐
│ Kanad Framework │  Python Package
│  (kanad/)       │
└─────────────────┘
```

---

## Backend API Implementation

### 1. Core API Structure

The API is already set up in `api/main.py`. Key additions needed:

```python
# api/main.py
from fastapi import FastAPI, BackgroundTasks
from pydantic import BaseModel
from typing import List, Optional
import uuid
from datetime import datetime

from kanad.bonds import BondFactory
from kanad.solvers import VQESolver, SQDSolver
from kanad.dynamics import MDSimulator
from kanad.applications import DrugDiscoveryPlatform, MaterialsScout
from kanad.analysis import PropertyCalculator, DOSCalculator

app = FastAPI(title="Kanad Quantum Chemistry API", version="1.0.0")

# Job storage (use Redis/database in production)
jobs_db = {}
```

### 2. Request/Response Models

```python
# api/models.py
from pydantic import BaseModel
from typing import List, Optional, Dict, Any

class Atom(BaseModel):
    element: str
    position: List[float]  # [x, y, z] in Bohr

class Molecule(BaseModel):
    atoms: List[Atom]
    charge: int = 0
    multiplicity: int = 1
    basis_set: str = "sto-3g"

class VQERequest(BaseModel):
    molecule: Molecule
    backend: str = "statevector"
    use_governance: bool = True
    optimizer: str = "SLSQP"
    max_iterations: int = 100

class VQEResponse(BaseModel):
    job_id: str
    status: str
    energy: Optional[float] = None
    converged: Optional[bool] = None
    iterations: Optional[int] = None
    properties: Optional[Dict] = None
    error: Optional[str] = None

class MDRequest(BaseModel):
    molecule: Molecule
    temperature: float  # Kelvin
    timestep: float  # femtoseconds
    n_steps: int
    force_method: str = "hf"  # hf, mp2, vqe, sqd
    use_governance: bool = True
    backend: str = "statevector"
    thermostat: Optional[str] = None  # berendsen, nose-hoover, langevin
    save_trajectory: bool = True

class MDResponse(BaseModel):
    job_id: str
    status: str
    initial_energy: Optional[float] = None
    final_energy: Optional[float] = None
    energy_drift: Optional[float] = None
    average_temperature: Optional[float] = None
    trajectory_file: Optional[str] = None
    error: Optional[str] = None
```

### 3. API Endpoints

#### A. Quantum Solvers

```python
@app.post("/api/v1/solvers/vqe", response_model=VQEResponse)
async def run_vqe(request: VQERequest, background_tasks: BackgroundTasks):
    """
    Run VQE calculation.

    Supports both synchronous (small molecules) and asynchronous (large molecules) execution.
    """
    job_id = str(uuid.uuid4())

    # For small molecules (< 4 qubits), run synchronously
    if len(request.molecule.atoms) <= 2:
        try:
            # Create bond
            atoms = request.molecule.atoms
            bond = BondFactory.create_bond(
                atoms[0].element,
                atoms[1].element,
                distance=calculate_distance(atoms[0].position, atoms[1].position)
            )

            # Run VQE
            solver = VQESolver(
                bond,
                backend=request.backend,
                use_governance=request.use_governance,
                optimizer=request.optimizer
            )
            result = solver.solve()

            return VQEResponse(
                job_id=job_id,
                status="completed",
                energy=result['energy'],
                converged=result['converged'],
                iterations=result.get('iterations', 0),
                properties=result.get('analysis', {})
            )

        except Exception as e:
            return VQEResponse(
                job_id=job_id,
                status="failed",
                error=str(e)
            )

    # For larger molecules, run in background
    else:
        jobs_db[job_id] = {"status": "queued", "created_at": datetime.now()}
        background_tasks.add_task(run_vqe_background, job_id, request)

        return VQEResponse(
            job_id=job_id,
            status="queued"
        )

def calculate_distance(pos1: List[float], pos2: List[float]) -> float:
    """Calculate distance between two positions."""
    import numpy as np
    return np.linalg.norm(np.array(pos1) - np.array(pos2))

async def run_vqe_background(job_id: str, request: VQERequest):
    """Run VQE in background."""
    try:
        jobs_db[job_id]["status"] = "running"

        # Create molecule and run VQE
        # ... implementation ...

        jobs_db[job_id].update({
            "status": "completed",
            "result": result
        })

    except Exception as e:
        jobs_db[job_id].update({
            "status": "failed",
            "error": str(e)
        })
```

#### B. Molecular Dynamics

```python
@app.post("/api/v1/md/classical", response_model=MDResponse)
async def run_classical_md(request: MDRequest):
    """Run classical MD simulation with HF or MP2 forces."""
    job_id = str(uuid.uuid4())

    try:
        # Create bond
        atoms = request.molecule.atoms
        bond = BondFactory.create_bond(
            atoms[0].element,
            atoms[1].element,
            distance=calculate_distance(atoms[0].position, atoms[1].position)
        )

        # Setup MD
        md = MDSimulator(
            bond,
            temperature=request.temperature,
            timestep=request.timestep,
            force_method=request.force_method,
            thermostat=request.thermostat
        )

        # Run simulation
        result = md.run(n_steps=request.n_steps)

        # Save trajectory
        trajectory_file = None
        if request.save_trajectory:
            trajectory_file = f"data/trajectories/{job_id}.h5"
            result.trajectory.save(trajectory_file)

        return MDResponse(
            job_id=job_id,
            status="completed",
            initial_energy=result.initial_energy,
            final_energy=result.final_energy,
            energy_drift=abs(result.final_energy - result.initial_energy),
            average_temperature=result.average_temperature,
            trajectory_file=trajectory_file
        )

    except Exception as e:
        return MDResponse(
            job_id=job_id,
            status="failed",
            error=str(e)
        )

@app.post("/api/v1/md/quantum", response_model=MDResponse)
async def run_quantum_md(request: MDRequest, background_tasks: BackgroundTasks):
    """
    Run quantum MD simulation with VQE or SQD forces.

    Always runs asynchronously due to computational cost.
    """
    job_id = str(uuid.uuid4())
    jobs_db[job_id] = {"status": "queued", "created_at": datetime.now()}

    background_tasks.add_task(run_quantum_md_background, job_id, request)

    return MDResponse(
        job_id=job_id,
        status="queued"
    )

async def run_quantum_md_background(job_id: str, request: MDRequest):
    """Run quantum MD in background."""
    try:
        jobs_db[job_id]["status"] = "running"

        # Create bond
        atoms = request.molecule.atoms
        bond = BondFactory.create_bond(
            atoms[0].element,
            atoms[1].element,
            distance=calculate_distance(atoms[0].position, atoms[1].position)
        )

        # Setup MD
        md = MDSimulator(
            bond,
            temperature=request.temperature,
            timestep=request.timestep,
            force_method=request.force_method,
            use_governance=request.use_governance,
            backend=request.backend
        )

        # Run simulation
        result = md.run(n_steps=request.n_steps)

        # Save trajectory
        trajectory_file = f"data/trajectories/{job_id}.h5"
        result.trajectory.save(trajectory_file)

        jobs_db[job_id].update({
            "status": "completed",
            "result": {
                "initial_energy": result.initial_energy,
                "final_energy": result.final_energy,
                "energy_drift": abs(result.final_energy - result.initial_energy),
                "average_temperature": result.average_temperature,
                "trajectory_file": trajectory_file
            }
        })

    except Exception as e:
        jobs_db[job_id].update({
            "status": "failed",
            "error": str(e)
        })
```

#### C. Job Status

```python
@app.get("/api/v1/jobs/{job_id}")
async def get_job_status(job_id: str):
    """Get status of a background job."""
    if job_id not in jobs_db:
        raise HTTPException(status_code=404, detail="Job not found")

    return jobs_db[job_id]
```

#### D. Applications

```python
@app.post("/api/v1/drug-discovery/evaluate")
async def evaluate_drug_molecule(smiles: str, use_quantum: bool = False):
    """Evaluate molecule for drug-like properties."""
    try:
        platform = DrugDiscoveryPlatform(backend='statevector')
        result = platform.evaluate_molecule(smiles, use_quantum=use_quantum)

        return {
            "status": "success",
            "smiles": smiles,
            "energy": result['energy'],
            "properties": result['properties'],
            "drug_score": result.get('drug_score', 0.0)
        }

    except Exception as e:
        return {
            "status": "error",
            "error": str(e)
        }

@app.post("/api/v1/materials/alloy-design")
async def design_alloy(elements: List[str], distances: List[float]):
    """Design alloy with specified elements."""
    try:
        scout = MaterialsScout(backend='statevector')
        result = scout.evaluate_material(elements, distances, use_quantum=True)

        return {
            "status": "success",
            "elements": elements,
            "energy": result['energy'],
            "properties": result.get('properties', {}),
            "band_gap": result.get('band_gap', 0.0)
        }

    except Exception as e:
        return {
            "status": "error",
            "error": str(e)
        }
```

---

## Web Frontend Integration

### 1. API Client Setup

```typescript
// web/src/lib/api.ts
import axios from 'axios';

const API_BASE_URL = process.env.NEXT_PUBLIC_API_URL || 'http://localhost:8000';

export const api = axios.create({
  baseURL: API_BASE_URL,
  headers: {
    'Content-Type': 'application/json',
  },
});

// Add authentication interceptor
api.interceptors.request.use((config) => {
  const token = localStorage.getItem('authToken');
  if (token) {
    config.headers.Authorization = `Bearer ${token}`;
  }
  return config;
});
```

### 2. API Service Functions

```typescript
// web/src/services/quantum.ts
import { api } from '@/lib/api';

export interface Atom {
  element: string;
  position: number[];
}

export interface Molecule {
  atoms: Atom[];
  charge?: number;
  multiplicity?: number;
  basis_set?: string;
}

export interface VQERequest {
  molecule: Molecule;
  backend?: string;
  use_governance?: boolean;
  optimizer?: string;
}

export interface MDRequest {
  molecule: Molecule;
  temperature: number;
  timestep: number;
  n_steps: number;
  force_method?: string;
  use_governance?: boolean;
  backend?: string;
}

// VQE calculation
export async function runVQE(request: VQERequest) {
  const response = await api.post('/api/v1/solvers/vqe', request);
  return response.data;
}

// Classical MD
export async function runClassicalMD(request: MDRequest) {
  const response = await api.post('/api/v1/md/classical', request);
  return response.data;
}

// Quantum MD
export async function runQuantumMD(request: MDRequest) {
  const response = await api.post('/api/v1/md/quantum', request);
  return response.data;
}

// Poll job status
export async function getJobStatus(jobId: string) {
  const response = await api.get(`/api/v1/jobs/${jobId}`);
  return response.data;
}

// Drug discovery
export async function evaluateDrugMolecule(smiles: string, useQuantum: boolean = false) {
  const response = await api.post('/api/v1/drug-discovery/evaluate', null, {
    params: { smiles, use_quantum: useQuantum }
  });
  return response.data;
}
```

### 3. React Components

```typescript
// web/src/components/QuantumSimulation.tsx
'use client';

import { useState } from 'react';
import { runVQE, getJobStatus } from '@/services/quantum';

export function QuantumSimulation() {
  const [isRunning, setIsRunning] = useState(false);
  const [result, setResult] = useState(null);

  const handleRunVQE = async () => {
    setIsRunning(true);

    try {
      const request = {
        molecule: {
          atoms: [
            { element: 'H', position: [0, 0, 0] },
            { element: 'H', position: [0, 0, 0.74] }
          ]
        },
        use_governance: true
      };

      const response = await runVQE(request);

      if (response.status === 'queued') {
        // Poll for results
        const jobId = response.job_id;
        const pollInterval = setInterval(async () => {
          const status = await getJobStatus(jobId);

          if (status.status === 'completed') {
            clearInterval(pollInterval);
            setResult(status.result);
            setIsRunning(false);
          } else if (status.status === 'failed') {
            clearInterval(pollInterval);
            alert(`Simulation failed: ${status.error}`);
            setIsRunning(false);
          }
        }, 2000);
      } else {
        // Completed immediately
        setResult(response);
        setIsRunning(false);
      }

    } catch (error) {
      console.error('VQE error:', error);
      setIsRunning(false);
    }
  };

  return (
    <div className="simulation-panel">
      <h2>VQE Simulation</h2>

      <button
        onClick={handleRunVQE}
        disabled={isRunning}
        className="btn-primary"
      >
        {isRunning ? 'Running...' : 'Run VQE'}
      </button>

      {result && (
        <div className="results">
          <h3>Results</h3>
          <p>Energy: {result.energy} Ha</p>
          <p>Converged: {result.converged ? 'Yes' : 'No'}</p>
          <p>Iterations: {result.iterations}</p>
        </div>
      )}
    </div>
  );
}
```

---

## Database Schema (Optional)

For production, use PostgreSQL or MongoDB:

```sql
-- Jobs table
CREATE TABLE jobs (
    id UUID PRIMARY KEY,
    user_id UUID REFERENCES users(id),
    job_type VARCHAR(50) NOT NULL,
    status VARCHAR(20) NOT NULL,
    created_at TIMESTAMP DEFAULT NOW(),
    started_at TIMESTAMP,
    completed_at TIMESTAMP,
    input_data JSONB,
    result_data JSONB,
    error_message TEXT
);

-- Molecules table
CREATE TABLE molecules (
    id UUID PRIMARY KEY,
    user_id UUID REFERENCES users(id),
    name VARCHAR(255),
    smiles TEXT,
    atoms JSONB,
    properties JSONB,
    created_at TIMESTAMP DEFAULT NOW()
);

-- Simulations table
CREATE TABLE simulations (
    id UUID PRIMARY KEY,
    job_id UUID REFERENCES jobs(id),
    molecule_id UUID REFERENCES molecules(id),
    simulation_type VARCHAR(50),
    parameters JSONB,
    trajectory_file TEXT,
    results JSONB,
    created_at TIMESTAMP DEFAULT NOW()
);
```

---

## Deployment

### 1. Backend (FastAPI)

```bash
# Install dependencies
pip install -r requirements.txt

# Run development server
uvicorn api.main:app --reload --host 0.0.0.0 --port 8000

# Production with Gunicorn
gunicorn api.main:app -w 4 -k uvicorn.workers.UvicornWorker --bind 0.0.0.0:8000
```

### 2. Frontend (Next.js)

```bash
# Install dependencies
cd web
npm install

# Development
npm run dev

# Production build
npm run build
npm start
```

### 3. Docker Deployment

```dockerfile
# Dockerfile
FROM python:3.9

WORKDIR /app

# Install dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application
COPY . .

# Expose port
EXPOSE 8000

# Run application
CMD ["gunicorn", "api.main:app", "-w", "4", "-k", "uvicorn.workers.UvicornWorker", "--bind", "0.0.0.0:8000"]
```

```yaml
# docker-compose.yml
version: '3.8'

services:
  backend:
    build: .
    ports:
      - "8000:8000"
    environment:
      - DATABASE_URL=postgresql://user:pass@db:5432/kanad
      - REDIS_URL=redis://redis:6379
    volumes:
      - ./data:/app/data
    depends_on:
      - db
      - redis

  frontend:
    build: ./web
    ports:
      - "3000:3000"
    environment:
      - NEXT_PUBLIC_API_URL=http://backend:8000
    depends_on:
      - backend

  db:
    image: postgres:14
    environment:
      - POSTGRES_USER=user
      - POSTGRES_PASSWORD=pass
      - POSTGRES_DB=kanad
    volumes:
      - postgres_data:/var/lib/postgresql/data

  redis:
    image: redis:7
    volumes:
      - redis_data:/data

volumes:
  postgres_data:
  redis_data:
```

---

## Testing

```python
# tests/test_api.py
from fastapi.testclient import TestClient
from api.main import app

client = TestClient(app)

def test_run_vqe():
    response = client.post("/api/v1/solvers/vqe", json={
        "molecule": {
            "atoms": [
                {"element": "H", "position": [0, 0, 0]},
                {"element": "H", "position": [0, 0, 0.74]}
            ]
        },
        "use_governance": True
    })

    assert response.status_code == 200
    data = response.json()
    assert data["status"] in ["completed", "queued"]
    if data["status"] == "completed":
        assert "energy" in data
        assert data["energy"] < 0  # H2 energy should be negative

def test_run_classical_md():
    response = client.post("/api/v1/md/classical", json={
        "molecule": {
            "atoms": [
                {"element": "H", "position": [0, 0, 0]},
                {"element": "H", "position": [0, 0, 0.74]}
            ]
        },
        "temperature": 300.0,
        "timestep": 0.5,
        "n_steps": 10,
        "force_method": "hf"
    })

    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "completed"
    assert "initial_energy" in data
    assert "final_energy" in data
```

---

## Security Considerations

1. **Authentication**: Use JWT tokens
2. **Rate Limiting**: Prevent abuse
3. **Input Validation**: Sanitize all inputs
4. **CORS**: Configure properly for production
5. **API Keys**: For quantum hardware access

---

## Monitoring

1. **Logging**: Use structured logging (Loguru)
2. **Metrics**: Track API response times, job completion rates
3. **Alerts**: Set up alerts for failures
4. **Dashboard**: Grafana + Prometheus

---

## Next Steps

1. ✅ Create API models
2. ✅ Implement core endpoints
3. ⏳ Add authentication
4. ⏳ Setup database
5. ⏳ Deploy to production
6. ⏳ Add monitoring
7. ⏳ Write API documentation (Swagger/OpenAPI)

---

**Status**: READY FOR IMPLEMENTATION
**Estimated Time**: 2-3 days for full deployment
**Dependencies**: FastAPI, PostgreSQL, Redis, Docker
