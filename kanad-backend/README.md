# Kanad Backend API

Comprehensive FastAPI backend for the Kanad quantum chemistry framework, enabling researchers to perform quantum chemistry calculations through a modern web interface with cloud quantum backend integration.

## Features

- **Complete REST API** for molecule management, quantum computations, and analysis
- **Real-time job monitoring** via WebSockets
- **Cloud quantum backend integration** (IBM Quantum, BlueQubit)
- **Async task processing** with Celery workers
- **Secure credential management** with encryption
- **PostgreSQL database** for data persistence
- **JWT authentication** for secure access
- **Comprehensive analysis tools** (energy decomposition, bonding, thermochemistry, spectroscopy)
- **Docker deployment** ready

## Architecture

```
┌─────────────────────────────────────────────────────────┐
│                   Frontend (Vercel)                      │
│               Next.js + React + WebSockets               │
└─────────────────────────────────────────────────────────┘
                         │ HTTPS
                         ▼
┌─────────────────────────────────────────────────────────┐
│               FastAPI Backend (Port 8000)                │
│  ┌────────────┐  ┌─────────────┐  ┌──────────────────┐ │
│  │  Routers   │  │  Services   │  │  Kanad Framework │ │
│  │  (REST +   │──│  (Business  │──│  (Quantum Chem)  │ │
│  │  WebSocket)│  │   Logic)    │  │                  │ │
│  └────────────┘  └─────────────┘  └──────────────────┘ │
└─────────────────────────────────────────────────────────┘
         │                    │                    │
         ▼                    ▼                    ▼
┌──────────────────┐  ┌──────────────┐  ┌────────────────┐
│   PostgreSQL     │  │    Redis     │  │ Celery Workers │
│  (Experiments)   │  │ (Queue+Cache)│  │ (Heavy Compute)│
└──────────────────┘  └──────────────┘  └────────────────┘
                                                  │
                                                  ▼
                                        ┌──────────────────┐
                                        │ Cloud Backends   │
                                        │ IBM / BlueQubit  │
                                        └──────────────────┘
```

## Quick Start

### Prerequisites

- Python 3.11+
- PostgreSQL 15+
- Redis 7+
- Docker & Docker Compose (optional)

### Installation

1. **Clone the repository**
```bash
cd /home/mk/deeprealm/kanad/kanad-backend
```

2. **Create virtual environment**
```bash
python3.11 -m venv venv
source venv/bin/activate
```

3. **Install dependencies**
```bash
pip install -r requirements.txt
```

4. **Install Kanad framework**
```bash
cd ..
pip install -e .
cd kanad-backend
```

5. **Configure environment**
```bash
cp .env.example .env
# Edit .env with your settings
```

**IMPORTANT**: Generate secure keys:
```bash
# JWT Secret Key
python -c "import secrets; print('JWT_SECRET_KEY=' + secrets.token_urlsafe(32))"

# Encryption Key (for credentials)
python -c "from cryptography.fernet import Fernet; print('ENCRYPTION_KEY=' + Fernet.generate_key().decode())"
```

6. **Setup database**
```bash
# Start PostgreSQL (if not running)
sudo systemctl start postgresql

# Create database
sudo -u postgres psql
CREATE DATABASE kanad_experiments;
CREATE USER kanad_user WITH PASSWORD 'kanad_password';
GRANT ALL PRIVILEGES ON DATABASE kanad_experiments TO kanad_user;
\q
```

7. **Run database migrations**
```bash
python -c "from core.database import init_db; init_db()"
```

8. **Start Redis**
```bash
sudo systemctl start redis-server
```

9. **Run the application**

Terminal 1 - FastAPI server:
```bash
uvicorn api.main:app --reload --host 0.0.0.0 --port 8000
```

Terminal 2 - Celery worker:
```bash
celery -A workers.celery_app worker --loglevel=info --concurrency=4
```

10. **Access the API**
- API Docs: http://localhost:8000/api/docs
- Health Check: http://localhost:8000/health

## Docker Deployment

### Development
```bash
docker-compose up --build
```

### Production
```bash
docker-compose -f docker-compose.yml up -d
```

## API Endpoints

### Authentication
- `POST /api/auth/register` - Register new user
- `POST /api/auth/login` - Login and get JWT token
- `GET /api/auth/profile` - Get current user profile

### Molecules
- `POST /api/molecules/create` - Create molecule (atoms/SMILES/library)
- `POST /api/molecules/from-smiles` - Create from SMILES
- `GET /api/molecules/{id}` - Get molecule by ID
- `GET /api/molecules/` - List user's molecules

### Simulations
- `POST /api/simulations/configure` - Configure computation
- `POST /api/simulations/{id}/accept-and-run` - Submit job

### Jobs
- `GET /api/jobs` - List jobs
- `GET /api/jobs/{id}/status` - Get job status
- `WS /api/jobs/{id}/logs` - Real-time logs (WebSocket)
- `DELETE /api/jobs/{id}` - Cancel job

### Analysis
- `GET /api/analysis/{job_id}/results` - Get results with analysis
- `GET /api/analysis/{job_id}/report` - Download PDF report
- `GET /api/analysis/{job_id}/export` - Export data (JSON/CSV/XYZ)

### Cloud
- `POST /api/cloud/credentials` - Store encrypted credentials
- `GET /api/cloud/backends` - List available backends

## Project Structure

```
kanad-backend/
├── api/
│   ├── main.py                 # FastAPI app entry point
│   ├── config.py               # Configuration management
│   └── routers/                # API endpoints
│       ├── auth.py             # Authentication
│       ├── molecules.py        # Molecule management
│       ├── simulations.py      # Computation configuration
│       ├── jobs.py             # Job management
│       ├── analysis.py         # Analysis endpoints
│       ├── cloud.py            # Cloud providers
│       └── library.py          # Molecule library
├── core/
│   ├── models.py               # Pydantic data models
│   ├── database.py             # Database connection
│   └── schemas.py              # API schemas
├── db/
│   ├── models.py               # SQLAlchemy database models
│   └── migrations/             # Alembic migrations
├── services/
│   ├── computation_service.py  # Kanad framework integration
│   ├── cloud_service.py        # Cloud backend integration
│   └── analysis_service.py     # Analysis computations
├── workers/
│   ├── celery_app.py           # Celery configuration
│   └── tasks.py                # Background tasks
├── utils/
│   ├── auth.py                 # JWT authentication
│   └── credentials_manager.py  # Credential encryption
├── tests/                      # Unit and integration tests
├── requirements.txt            # Python dependencies
├── Dockerfile                  # Docker image definition
├── docker-compose.yml          # Multi-container setup
└── README.md                   # This file
```

## Configuration

Key environment variables (see `.env.example`):

```bash
# Database
DATABASE_URL=postgresql://kanad_user:password@localhost:5432/kanad_experiments

# Redis & Celery
REDIS_URL=redis://localhost:6379/0

# Security
JWT_SECRET_KEY=<generate_with_secrets>
ENCRYPTION_KEY=<generate_with_fernet>

# Cloud Providers (optional defaults)
IBM_API_TOKEN=your_ibm_token
IBM_CRN=your_ibm_crn
BLUEQUBIT_TOKEN=your_bluequbit_token

# AI Services
ANTHROPIC_API_KEY=your_anthropic_key
```

## Testing

```bash
# Install test dependencies
pip install pytest pytest-asyncio pytest-cov

# Run tests
pytest

# With coverage
pytest --cov=. --cov-report=html
```

## Production Deployment

### Azure VM Setup

1. **Create Azure VM**
   - Size: Standard_D8s_v5 (8 vCPUs, 32 GB RAM)
   - Image: Ubuntu 22.04 LTS
   - Open ports: 22 (SSH), 80 (HTTP), 443 (HTTPS)

2. **Install dependencies**
```bash
sudo apt update && sudo apt upgrade -y
sudo apt install -y python3.11 python3.11-venv postgresql redis-server nginx
```

3. **Deploy application**
```bash
git clone <your-repo>
cd kanad-backend
python3.11 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

4. **Configure systemd services**

Create `/etc/systemd/system/kanad-api.service`:
```ini
[Unit]
Description=Kanad FastAPI Application
After=network.target

[Service]
User=azureuser
WorkingDirectory=/home/azureuser/kanad-backend
Environment="PATH=/home/azureuser/kanad-backend/venv/bin"
EnvironmentFile=/home/azureuser/kanad-backend/.env
ExecStart=/home/azureuser/kanad-backend/venv/bin/uvicorn api.main:app --host 0.0.0.0 --port 8000 --workers 4
Restart=always

[Install]
WantedBy=multi-user.target
```

Create `/etc/systemd/system/kanad-worker.service`:
```ini
[Unit]
Description=Kanad Celery Worker
After=network.target

[Service]
User=azureuser
WorkingDirectory=/home/azureuser/kanad-backend
Environment="PATH=/home/azureuser/kanad-backend/venv/bin"
EnvironmentFile=/home/azureuser/kanad-backend/.env
ExecStart=/home/azureuser/kanad-backend/venv/bin/celery -A workers.celery_app worker --loglevel=info --concurrency=4
Restart=always

[Install]
WantedBy=multi-user.target
```

Enable and start:
```bash
sudo systemctl enable kanad-api kanad-worker
sudo systemctl start kanad-api kanad-worker
```

5. **Configure Nginx**

Create `/etc/nginx/sites-available/kanad-api`:
```nginx
server {
    listen 80;
    server_name your-domain.com;

    location / {
        proxy_pass http://127.0.0.1:8000;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
    }

    location /api/jobs/ {
        proxy_pass http://127.0.0.1:8000;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
    }
}
```

Enable and restart:
```bash
sudo ln -s /etc/nginx/sites-available/kanad-api /etc/nginx/sites-enabled/
sudo nginx -t
sudo systemctl restart nginx
```

6. **HTTPS with Let's Encrypt**
```bash
sudo apt install certbot python3-certbot-nginx
sudo certbot --nginx -d your-domain.com
```

## Monitoring

- **Health check**: `curl http://localhost:8000/health`
- **Service status**: `sudo systemctl status kanad-api kanad-worker`
- **Logs**: `sudo journalctl -u kanad-api -f`

## Security Best Practices

1. ✅ Use strong JWT secrets and encryption keys
2. ✅ Enable HTTPS in production
3. ✅ Restrict database access to localhost
4. ✅ Use environment variables for secrets
5. ✅ Enable rate limiting
6. ✅ Regular security updates
7. ✅ Firewall configuration (allow only 22, 80, 443)

## Contributing

1. Fork the repository
2. Create feature branch (`git checkout -b feature/amazing-feature`)
3. Commit changes (`git commit -m 'Add amazing feature'`)
4. Push to branch (`git push origin feature/amazing-feature`)
5. Open Pull Request

## License

See Kanad framework license.

## Support

- Documentation: `/api/docs`
- Issues: GitHub Issues
- Email: support@kanad.example.com

---

**Built with Kanad Quantum Chemistry Framework**
