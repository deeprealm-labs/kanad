# Backend Integration Plan for Kanad Web Application

## Current Status

### ✅ Completed Frontend Features
1. **Dashboard Home** - Statistics, recent experiments, quick actions
2. **Molecule Creator** - Drag-drop atoms, SMILES input, pre-built molecules
3. **Lewis Structure View** - Visual representation using SmilesDrawer
4. **Preview/Configuration** - Review settings, select analysis properties
5. **Experiment Monitor** - Real-time progress tracking with live convergence charts
6. **Experimentation History** - Browse, search, filter past experiments
7. **Job Queue Management** - Schedule, prioritize, pause/resume experiments
8. **Settings Modal** - Configure methods, backends, optimizers
9. **Local Persistence** - localStorage for experiments and queue

### 🔧 Frontend TODO Items

#### High Priority
1. **Settings Persistence** ⚠️
   - Load settings from localStorage on mount
   - Save settings to localStorage on save
   - Use saved settings as defaults for new experiments
   - File: `/web/src/components/settings/SettingsModal.tsx`

2. **Export/Download Functionality**
   - Implement actual download of experiment results
   - Support formats: JSON, CSV, HDF5
   - Include: molecule data, configuration, results, convergence data
   - Files:
     - `/web/src/app/dashboard/history/page.tsx` (line 202)
     - `/web/src/components/simulation/ExperimentMonitor.tsx` (line 189)

3. **Error Handling & Validation**
   - Validate SMILES input
   - Handle invalid molecule configurations
   - Show error messages for failed experiments
   - Add loading states and spinners

#### Medium Priority
4. **Enhanced Molecule Library**
   - Add more pre-built molecules
   - Categorize by complexity, use case
   - Add 3D structure preview
   - File: `/web/src/data/molecule-library.ts`

5. **User Authentication UI**
   - Login/Register pages
   - User profile management
   - API key management for IBM/BlueQubit

6. **Real-time Updates**
   - WebSocket connection for live experiment updates
   - Multi-tab synchronization
   - Notifications for completed experiments

#### Low Priority
7. **Improved Visualizations**
   - 3D molecular structure viewer (Three.js/Mol*)
   - Interactive convergence charts
   - Energy level diagrams
   - Orbital visualizations

8. **Batch Operations**
   - Select multiple experiments for deletion
   - Bulk export
   - Compare multiple experiments

9. **Dark Mode Polish**
   - Ensure all components look good in dark mode
   - Chart colors adapt to theme

---

## Backend Architecture Plan

### Tech Stack
- **Framework**: FastAPI (already in requirements.txt)
- **Server**: Uvicorn with async/await
- **Database**: PostgreSQL or SQLite for development
- **Queue**: Celery + Redis for job processing
- **WebSocket**: FastAPI WebSocket support
- **Authentication**: JWT tokens
- **File Storage**: Local filesystem or S3

### API Structure

```
/api/v1/
├── /auth/
│   ├── POST /register
│   ├── POST /login
│   ├── POST /refresh
│   └── GET /me
├── /experiments/
│   ├── GET /experiments (list with pagination)
│   ├── POST /experiments (create/submit)
│   ├── GET /experiments/{id}
│   ├── DELETE /experiments/{id}
│   ├── GET /experiments/{id}/results
│   ├── GET /experiments/{id}/download (export)
│   └── GET /experiments/{id}/status
├── /queue/
│   ├── GET /queue (list queued jobs)
│   ├── POST /queue (add to queue)
│   ├── PUT /queue/{id} (update priority/schedule)
│   ├── DELETE /queue/{id}
│   └── POST /queue/{id}/execute
├── /molecules/
│   ├── POST /molecules/validate (validate SMILES)
│   ├── POST /molecules/parse (parse molecule)
│   └── GET /molecules/library
├── /settings/
│   ├── GET /settings
│   └── PUT /settings
└── /ws/
    └── /ws/experiments/{id} (WebSocket for live updates)
```

### Database Schema

#### Users Table
```sql
CREATE TABLE users (
    id SERIAL PRIMARY KEY,
    email VARCHAR(255) UNIQUE NOT NULL,
    password_hash VARCHAR(255) NOT NULL,
    name VARCHAR(255),
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
```

#### User Settings Table
```sql
CREATE TABLE user_settings (
    id SERIAL PRIMARY KEY,
    user_id INTEGER REFERENCES users(id),
    method VARCHAR(50) DEFAULT 'VQE',
    ansatz VARCHAR(50) DEFAULT 'hardware_efficient',
    mapper VARCHAR(50) DEFAULT 'jordan_wigner',
    optimizer VARCHAR(50) DEFAULT 'SLSQP',
    backend VARCHAR(50) DEFAULT 'classical',
    backend_name VARCHAR(100),
    optimization_settings JSONB,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);
```

#### Experiments Table
```sql
CREATE TABLE experiments (
    id SERIAL PRIMARY KEY,
    user_id INTEGER REFERENCES users(id),
    name VARCHAR(255),
    molecule_data JSONB NOT NULL,
    configuration JSONB NOT NULL,
    status VARCHAR(50) NOT NULL,
    results JSONB,
    convergence_data JSONB,
    error_message TEXT,
    started_at TIMESTAMP,
    completed_at TIMESTAMP,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX idx_experiments_user_status ON experiments(user_id, status);
CREATE INDEX idx_experiments_created ON experiments(created_at DESC);
```

#### Job Queue Table
```sql
CREATE TABLE job_queue (
    id SERIAL PRIMARY KEY,
    user_id INTEGER REFERENCES users(id),
    experiment_id INTEGER REFERENCES experiments(id),
    priority INTEGER DEFAULT 0,
    status VARCHAR(50) NOT NULL,
    scheduled_time TIMESTAMP,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX idx_queue_status_priority ON job_queue(status, priority);
CREATE INDEX idx_queue_scheduled ON job_queue(scheduled_time);
```

### Backend Implementation Plan

#### Phase 1: Core API (Week 1-2)
1. **Setup FastAPI Application**
   - Create `/api/main.py`
   - Configure CORS for Next.js frontend
   - Setup middleware (auth, logging, error handling)
   - Database connection pooling

2. **Implement Experiment Endpoints**
   - Create experiment submission endpoint
   - Validate molecule data (SMILES, atoms)
   - Store experiment in database
   - Return experiment ID

3. **Basic Job Queue**
   - Queue system using Celery
   - Worker process to execute experiments
   - Status updates

#### Phase 2: Quantum Integration (Week 3-4)
1. **Integrate Kanad Framework**
   - Wrapper functions for VQE, QPE, SQD solvers
   - Backend selection (Classical, IBM, BlueQubit)
   - Error handling and timeouts

2. **Example Integration Code**:
```python
# /api/services/experiment_service.py
from kanad.core.molecule import Molecule
from kanad.solvers.vqe_solver import VQESolver
from kanad.backends.ibm import IBMRuntimeBackend
from kanad.backends.bluequbit import BlueQubitBackend

async def run_experiment(experiment_config):
    # Parse molecule
    if experiment_config['molecule']['smiles']:
        mol = Molecule.from_smiles(
            experiment_config['molecule']['smiles'],
            basis=experiment_config['molecule']['basis'],
            charge=experiment_config['molecule']['charge'],
            multiplicity=experiment_config['molecule']['multiplicity']
        )
    else:
        # Create from atoms
        mol = Molecule(
            atoms=experiment_config['molecule']['atoms'],
            basis=experiment_config['molecule']['basis'],
            charge=experiment_config['molecule']['charge'],
            multiplicity=experiment_config['molecule']['multiplicity']
        )

    # Select backend
    backend = None
    if experiment_config['backend'] == 'ibm_quantum':
        backend = IBMRuntimeBackend(
            backend_name=experiment_config['backend_name']
        )
    elif experiment_config['backend'] == 'bluequbit':
        backend = BlueQubitBackend()

    # Run VQE
    solver = VQESolver(
        molecule=mol,
        ansatz_type=experiment_config['ansatz'],
        mapper=experiment_config['mapper'],
        optimizer=experiment_config['optimizer'],
        backend=backend
    )

    result = solver.solve()

    return {
        'energy': result.energy,
        'iterations': result.iterations,
        'converged': result.converged,
        'dipole_moment': result.dipole_moment,
        'convergence_data': result.convergence_history
    }
```

#### Phase 3: Real-time Updates (Week 5)
1. **WebSocket Implementation**
   - Live progress updates during execution
   - Convergence data streaming
   - Status changes

2. **Frontend WebSocket Connection**
```typescript
// /web/src/hooks/useExperimentWebSocket.ts
export function useExperimentWebSocket(experimentId: string) {
  const [status, setStatus] = useState<string>('queued');
  const [progress, setProgress] = useState(0);
  const [convergenceData, setConvergenceData] = useState([]);

  useEffect(() => {
    const ws = new WebSocket(
      `ws://localhost:8000/api/v1/ws/experiments/${experimentId}`
    );

    ws.onmessage = (event) => {
      const data = JSON.parse(event.data);
      setStatus(data.status);
      setProgress(data.progress);
      if (data.convergence_point) {
        setConvergenceData(prev => [...prev, data.convergence_point]);
      }
    };

    return () => ws.close();
  }, [experimentId]);

  return { status, progress, convergenceData };
}
```

#### Phase 4: Authentication & Security (Week 6)
1. **JWT Authentication**
   - User registration/login
   - Token refresh mechanism
   - Protected endpoints

2. **API Key Management**
   - Store IBM Quantum API keys securely
   - BlueQubit API keys
   - Encrypt sensitive data

#### Phase 5: Production Features (Week 7-8)
1. **Rate Limiting**
   - Prevent abuse
   - Queue management per user

2. **Monitoring & Logging**
   - Structured logging
   - Performance metrics
   - Error tracking (Sentry)

3. **Deployment**
   - Docker containers
   - Kubernetes deployment
   - CI/CD pipeline

### File Structure

```
/api/
├── main.py                 # FastAPI app entry point
├── config.py              # Configuration management
├── database.py            # Database connection
├── models/
│   ├── user.py
│   ├── experiment.py
│   ├── job_queue.py
│   └── settings.py
├── routers/
│   ├── auth.py
│   ├── experiments.py
│   ├── queue.py
│   ├── molecules.py
│   └── settings.py
├── services/
│   ├── auth_service.py
│   ├── experiment_service.py
│   ├── queue_service.py
│   └── kanad_wrapper.py
├── tasks/
│   └── celery_tasks.py    # Background job processing
├── websockets/
│   └── experiment_ws.py
├── utils/
│   ├── validators.py
│   ├── serializers.py
│   └── exceptions.py
└── tests/
    ├── test_experiments.py
    └── test_queue.py
```

### Environment Variables

```bash
# .env
DATABASE_URL=postgresql://user:password@localhost:5432/kanad
REDIS_URL=redis://localhost:6379/0
JWT_SECRET_KEY=your-secret-key-here
JWT_ALGORITHM=HS256
IBM_QUANTUM_TOKEN=your-ibm-token
BLUEQUBIT_API_KEY=your-bluequbit-key
CORS_ORIGINS=http://localhost:3000,https://kanad.deeprealm.in
```

### Next.js API Integration

Update frontend to use backend API instead of localStorage:

```typescript
// /web/src/lib/api.ts
const API_BASE = process.env.NEXT_PUBLIC_API_URL || 'http://localhost:8000/api/v1';

export async function submitExperiment(config: any) {
  const response = await fetch(`${API_BASE}/experiments`, {
    method: 'POST',
    headers: {
      'Content-Type': 'application/json',
      'Authorization': `Bearer ${getToken()}`
    },
    body: JSON.stringify(config)
  });

  if (!response.ok) throw new Error('Failed to submit experiment');
  return response.json();
}

export async function getExperiments() {
  const response = await fetch(`${API_BASE}/experiments`, {
    headers: {
      'Authorization': `Bearer ${getToken()}`
    }
  });

  if (!response.ok) throw new Error('Failed to fetch experiments');
  return response.json();
}

// ... more API functions
```

---

## Deployment Strategy

### Development Environment
- Next.js: `npm run dev` (port 3000)
- FastAPI: `uvicorn api.main:app --reload` (port 8000)
- PostgreSQL: Docker container
- Redis: Docker container
- Celery Worker: `celery -A api.tasks worker`

### Production Environment
- **Frontend**: Vercel or AWS Amplify
- **Backend**: AWS ECS/EKS or DigitalOcean App Platform
- **Database**: AWS RDS PostgreSQL
- **Queue**: AWS ElastiCache Redis
- **File Storage**: AWS S3

### Docker Compose for Development

```yaml
version: '3.8'
services:
  postgres:
    image: postgres:15
    environment:
      POSTGRES_DB: kanad
      POSTGRES_USER: kanad
      POSTGRES_PASSWORD: kanad123
    ports:
      - "5432:5432"
    volumes:
      - postgres_data:/var/lib/postgresql/data

  redis:
    image: redis:7
    ports:
      - "6379:6379"

  api:
    build: ./api
    command: uvicorn main:app --host 0.0.0.0 --port 8000 --reload
    volumes:
      - ./api:/app
    ports:
      - "8000:8000"
    environment:
      DATABASE_URL: postgresql://kanad:kanad123@postgres:5432/kanad
      REDIS_URL: redis://redis:6379/0
    depends_on:
      - postgres
      - redis

  celery:
    build: ./api
    command: celery -A tasks worker --loglevel=info
    volumes:
      - ./api:/app
    environment:
      DATABASE_URL: postgresql://kanad:kanad123@postgres:5432/kanad
      REDIS_URL: redis://redis:6379/0
    depends_on:
      - postgres
      - redis

volumes:
  postgres_data:
```

---

## Summary & Recommendations

### Immediate Actions (This Week)
1. **Fix Settings Persistence** - 30 minutes
2. **Implement Export Functionality** - 2 hours
3. **Add Basic Error Handling** - 2 hours

### Backend Development (Next 2 Months)
1. **Week 1-2**: Core API with PostgreSQL
2. **Week 3-4**: Integrate Kanad quantum computing
3. **Week 5**: WebSocket real-time updates
4. **Week 6**: Authentication & security
5. **Week 7-8**: Production readiness & deployment

### Success Metrics
- API response time < 200ms (non-compute endpoints)
- Experiment submission to queue < 1 second
- Support 100+ concurrent users
- 99.9% uptime for API
- Real-time updates with < 500ms latency

### Estimated Effort
- **Frontend Polish**: 1 week (40 hours)
- **Backend Development**: 6-8 weeks (240-320 hours)
- **Testing & QA**: 2 weeks (80 hours)
- **Deployment & DevOps**: 1 week (40 hours)

**Total**: ~10-12 weeks for production-ready system

---

## Contact & Next Steps

Once you're ready to begin backend integration:
1. Set up development database
2. Create API skeleton with FastAPI
3. Implement authentication first
4. Build experiments endpoint
5. Integrate Kanad framework
6. Test end-to-end workflow
7. Deploy to staging environment

Let me know which phase you'd like to start with!
