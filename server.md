## ☁️ Recommended Cloud Architecture

### **Best Setup for Your Use Case**

```
┌─────────────────────────────────────────────────────┐
│                    USER ACCESS                       │
└─────────────────────────────────────────────────────┘
                        ↓
┌─────────────────────────────────────────────────────┐
│  VERCEL (Frontend)                                   │
│  • Next.js 14 app                                    │
│  • Free tier / $20/month Pro                         │
│  • Auto-scaling, global CDN                          │
│  • HTTPS included                                    │
└─────────────────────────────────────────────────────┘
                        ↓ HTTPS
┌─────────────────────────────────────────────────────┐
│  AZURE VM (Backend API + Workers)                   │
│  • FastAPI application                               │
│  • Celery workers for computations                   │
│  • Redis for task queue                              │
│  • Local PostgreSQL for experiments                  │
└─────────────────────────────────────────────────────┘
                        ↓
┌─────────────────────────────────────────────────────┐
│  SUPABASE / NEON (User Database - Free tier)        │
│  • PostgreSQL for user data                          │
│  • Automatic backups                                 │
│  • Can migrate easily                                │
└─────────────────────────────────────────────────────┘
```

---

## 🎯 Why This Architecture?

### **1. Vercel for Frontend (FREE or $20/month)**

✅ **Zero configuration deployment** ✅ **Automatic HTTPS & CDN** ✅ **Serverless functions for light tasks** ✅ **Git integration (auto-deploy on push)**

### **2. Azure VM for Backend (Use your $1000 credits)**

✅ **Full control over compute resources** ✅ **High CPU/RAM for Hamiltonian computations** ✅ **Can install scientific libraries (NumPy, SciPy, PySCF)** ✅ **Run long-running Celery workers**

### **3. Separate User DB (Supabase/Neon - FREE)**

✅ **User data separate from computation data** ✅ **Easy migration if you upgrade/move** ✅ **Automatic backups** ✅ **PostgreSQL compatible with Prisma**

### **4. Experiment DB on VM**

✅ **Large result datasets stay close to compute** ✅ **Faster data access during computations** ✅ **Can backup/export when needed**

---

## 💰 Cost Breakdown (Maximize $1000 Azure Credits)

### **Recommended Azure VM Specs**

For quantum chemistry computations, you need **high CPU + RAM**:

|VM Type|vCPUs|RAM|Storage|Cost/Month|Best For|
|---|---|---|---|---|---|
|**Standard_D8s_v5**|8|32 GB|256 GB SSD|~$330|**Recommended** - Good balance|
|Standard_D16s_v5|16|64 GB|512 GB SSD|~$660|Heavy computations|
|Standard_D4s_v5|4|16 GB|128 GB SSD|~$165|Development/testing|

**My Recommendation**: Start with **Standard_D8s_v5**

- $330/month = **3 months** of runtime with your credits
- Enough power for H₂O, CO₂, moderate molecules
- Can upgrade to D16s_v5 if needed

### **Additional Azure Costs**

|Service|Cost|Notes|
|---|---|---|
|Bandwidth (outbound)|~$0.05/GB|Minimal for API responses|
|Static IP|~$3/month|Recommended for stable DNS|
|Azure Backup|~$10/month|Important for experiment data|
|**Total Estimate**|**~$343/month**|**~3 months** of credits|

### **Free Tier Services**

|Service|Tier|Notes|
|---|---|---|
|**Vercel**|Free (Hobby)|Unlimited projects, 100 GB bandwidth/month|
|**Supabase**|Free|500 MB database, 1 GB file storage|
|**Neon**|Free|3 GB storage, auto-suspend when idle|

---

## 🏗️ Detailed Setup Guide

### **PHASE 1: Azure VM Setup**

#### **1.1 Create VM on Azure Portal**

```bash
# Via Azure Portal:
1. Go to portal.azure.com
2. Create Resource → Virtual Machine
3. Select:
   - Region: East US (cheapest) or your nearest region
   - Image: Ubuntu 22.04 LTS
   - Size: Standard_D8s_v5 (8 vCPUs, 32 GB RAM)
   - Authentication: SSH key (generate new)
   - Disk: 256 GB Premium SSD

4. Networking:
   - Allow SSH (22)
   - Allow HTTPS (443)
   - Allow HTTP (80) - will redirect to HTTPS
   - Allow custom port 8000 (FastAPI during dev)

5. Create
```

#### **1.2 SSH into VM & Initial Setup**

```bash
# Connect to VM
ssh azureuser@<your-vm-ip>

# Update system
sudo apt update && sudo apt upgrade -y

# Install Python 3.13 (if not available, use 3.11)
sudo apt install -y python3.11 python3.11-venv python3-pip

# Install system dependencies
sudo apt install -y build-essential cmake git curl wget
sudo apt install -y libopenblas-dev liblapack-dev gfortran
sudo apt install -y redis-server postgresql postgresql-contrib nginx

# Install Docker (for easier deployment later)
curl -fsSL https://get.docker.com -o get-docker.sh
sudo sh get-docker.sh
sudo usermod -aG docker $USER
```

#### **1.3 Install Kanad & Backend Dependencies**

```bash
# Clone your Kanad repository
cd /home/azureuser
git clone https://github.com/yourusername/kanad.git
cd kanad

# Create virtual environment
python3.11 -m venv venv
source venv/bin/activate

# Install Kanad
pip install -e .

# Install backend dependencies
pip install fastapi uvicorn[standard] celery redis pydantic sqlalchemy asyncpg
pip install python-jose[cryptography] passlib[bcrypt] python-multipart
pip install anthropic  # For LLM reports
pip install prisma  # If using Prisma ORM
```

#### **1.4 Setup PostgreSQL (Experiment Database)**

```bash
# Start PostgreSQL
sudo systemctl start postgresql
sudo systemctl enable postgresql

# Create database
sudo -u postgres psql
```

```sql
-- In PostgreSQL shell
CREATE DATABASE kanad_experiments;
CREATE USER kanad_user WITH PASSWORD 'your_secure_password';
GRANT ALL PRIVILEGES ON DATABASE kanad_experiments TO kanad_user;
\q
```

#### **1.5 Setup Redis (Task Queue)**

```bash
# Start Redis
sudo systemctl start redis-server
sudo systemctl enable redis-server

# Test Redis
redis-cli ping  # Should return PONG
```

---

### **PHASE 2: Backend API Setup**

#### **2.1 Project Structure**

```bash
cd /home/azureuser
mkdir kanad-backend && cd kanad-backend

# Create directory structure
mkdir -p api/{routers,core,services,workers,utils,db}
mkdir -p tests
```

#### **2.2 FastAPI App (`api/main.py`)**

```python
# api/main.py
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from api.routers import molecules, simulations, jobs, auth
import os

app = FastAPI(
    title="Kanad Quantum Chemistry API",
    version="1.0.0",
    docs_url="/docs",
    redoc_url="/redoc"
)

# CORS - allow Vercel frontend
app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "https://your-app.vercel.app",
        "http://localhost:3000"  # For local development
    ],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Include routers
app.include_router(auth.router, prefix="/api/auth", tags=["auth"])
app.include_router(molecules.router, prefix="/api/molecules", tags=["molecules"])
app.include_router(simulations.router, prefix="/api/simulations", tags=["simulations"])
app.include_router(jobs.router, prefix="/api/jobs", tags=["jobs"])

@app.get("/")
async def root():
    return {"message": "Kanad API is running", "version": "1.0.0"}

@app.get("/health")
async def health():
    return {"status": "healthy"}
```

#### **2.3 Environment Variables (`.env`)**

```bash
# .env
DATABASE_URL=postgresql://kanad_user:your_secure_password@localhost:5432/kanad_experiments
USER_DATABASE_URL=postgresql://user:pass@supabase-or-neon-url:5432/users  # From Supabase/Neon
REDIS_URL=redis://localhost:6379/0

JWT_SECRET_KEY=your_super_secret_jwt_key_generate_with_openssl_rand_hex_32
ENCRYPTION_KEY=your_encryption_key_for_credentials

ANTHROPIC_API_KEY=your_anthropic_api_key

# IBM Quantum (default - users can override)
IBM_API_TOKEN=your_default_ibm_token
IBM_CRN=your_default_ibm_crn

# BlueQubit (default)
BLUEQUBIT_TOKEN=your_default_bluequbit_token

# Environment
ENVIRONMENT=production
LOG_LEVEL=INFO
```

#### **2.4 Run with Uvicorn**

```bash
# Development
uvicorn api.main:app --reload --host 0.0.0.0 --port 8000

# Production with systemd service
sudo nano /etc/systemd/system/kanad-api.service
```

```ini
[Unit]
Description=Kanad FastAPI Application
After=network.target postgresql.service redis.service

[Service]
User=azureuser
Group=azureuser
WorkingDirectory=/home/azureuser/kanad-backend
Environment="PATH=/home/azureuser/kanad-backend/venv/bin"
EnvironmentFile=/home/azureuser/kanad-backend/.env
ExecStart=/home/azureuser/kanad-backend/venv/bin/uvicorn api.main:app --host 0.0.0.0 --port 8000 --workers 4

Restart=always
RestartSec=10

[Install]
WantedBy=multi-user.target
```

```bash
# Enable and start service
sudo systemctl daemon-reload
sudo systemctl enable kanad-api
sudo systemctl start kanad-api
sudo systemctl status kanad-api
```

---

### **PHASE 3: Celery Workers for Heavy Computations**

#### **3.1 Celery App (`workers/celery_app.py`)**

```python
# workers/celery_app.py
from celery import Celery
import os

celery_app = Celery(
    "kanad",
    broker=os.getenv("REDIS_URL", "redis://localhost:6379/0"),
    backend=os.getenv("REDIS_URL", "redis://localhost:6379/0")
)

celery_app.conf.update(
    task_serializer='json',
    accept_content=['json'],
    result_serializer='json',
    timezone='UTC',
    enable_utc=True,
    task_track_started=True,
    task_time_limit=3600,  # 1 hour max per task
)
```

#### **3.2 Computation Tasks (`workers/tasks.py`)**

```python
# workers/tasks.py
from workers.celery_app import celery_app
from kanad.bonds import BondFactory
from kanad.solvers import VQESolver
import asyncio

@celery_app.task(bind=True)
def run_vqe_computation(self, molecule_data, config):
    """
    Heavy VQE computation task
    self.update_state() allows real-time progress updates
    """
    try:
        # Update state: PENDING → STARTED
        self.update_state(state='PROGRESS', meta={'progress': 0, 'message': 'Initializing...'})
        
        # Create molecule
        molecule = BondFactory.create_molecule(
            molecule_data['atoms'],
            basis=config['basis']
        )
        
        self.update_state(state='PROGRESS', meta={'progress': 10, 'message': 'Building Hamiltonian...'})
        
        # Create VQE solver
        solver = VQESolver(
            molecule,
            ansatz_type=config['ansatz'],
            mapper_type=config['mapper'],
            max_iterations=config['max_iterations']
        )
        
        # Run VQE with progress callbacks
        def progress_callback(iteration, energy):
            progress = 10 + int(80 * iteration / config['max_iterations'])
            self.update_state(
                state='PROGRESS',
                meta={
                    'progress': progress,
                    'message': f'Iteration {iteration}: E={energy:.6f} Ha'
                }
            )
        
        result = solver.solve(callback=progress_callback)
        
        self.update_state(state='PROGRESS', meta={'progress': 90, 'message': 'Running analysis...'})
        
        # Run analysis
        # ... analysis code ...
        
        self.update_state(state='PROGRESS', meta={'progress': 100, 'message': 'Complete!'})
        
        return {
            'status': 'completed',
            'energy': result['energy'],
            'results': result
        }
        
    except Exception as e:
        self.update_state(state='FAILURE', meta={'error': str(e)})
        raise
```

#### **3.3 Start Celery Workers**

```bash
# Create systemd service
sudo nano /etc/systemd/system/kanad-worker.service
```

```ini
[Unit]
Description=Kanad Celery Worker
After=network.target redis.service

[Service]
User=azureuser
Group=azureuser
WorkingDirectory=/home/azureuser/kanad-backend
Environment="PATH=/home/azureuser/kanad-backend/venv/bin"
EnvironmentFile=/home/azureuser/kanad-backend/.env
ExecStart=/home/azureuser/kanad-backend/venv/bin/celery -A workers.celery_app worker --loglevel=info --concurrency=4

Restart=always
RestartSec=10

[Install]
WantedBy=multi-user.target
```

```bash
sudo systemctl enable kanad-worker
sudo systemctl start kanad-worker
sudo systemctl status kanad-worker
```

---

### **PHASE 4: Nginx Reverse Proxy & HTTPS**

#### **4.1 Install Certbot (Let's Encrypt)**

```bash
sudo apt install -y certbot python3-certbot-nginx
```

#### **4.2 Setup Domain (Optional but Recommended)**

```bash
# Option 1: Buy domain ($10/year from Namecheap, GoDaddy)
# Point A record to your Azure VM IP

# Option 2: Use free subdomain from services like:
# - DuckDNS.org
# - FreeDNS.afraid.org
# Example: kanad-api.duckdns.org → <your-vm-ip>
```

#### **4.3 Nginx Configuration**

```bash
sudo nano /etc/nginx/sites-available/kanad-api
```

```nginx
server {
    listen 80;
    server_name your-domain.com;  # or your-vm-ip

    location / {
        proxy_pass http://127.0.0.1:8000;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
    }

    # WebSocket support for real-time job updates
    location /api/jobs/ {
        proxy_pass http://127.0.0.1:8000;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
        proxy_set_header Host $host;
    }
}
```

```bash
# Enable site
sudo ln -s /etc/nginx/sites-available/kanad-api /etc/nginx/sites-enabled/
sudo nginx -t
sudo systemctl restart nginx

# Get SSL certificate (if using domain)
sudo certbot --nginx -d your-domain.com
```

---

### **PHASE 5: User Database (Supabase - FREE)**

#### **5.1 Setup Supabase**

```bash
# Go to https://supabase.com
# 1. Sign up (free account)
# 2. Create new project
# 3. Get connection string from Settings → Database
```

**Connection String**:

```
postgresql://postgres:[YOUR-PASSWORD]@db.[PROJECT-REF].supabase.co:5432/postgres
```

#### **5.2 Create User Schema**

```sql
-- In Supabase SQL Editor
CREATE TABLE users (
    user_id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    email VARCHAR(255) UNIQUE NOT NULL,
    name VARCHAR(255),
    hashed_password VARCHAR(255) NOT NULL,
    institution VARCHAR(255),
    research_field VARCHAR(100),
    created_at TIMESTAMP DEFAULT NOW(),
    last_login TIMESTAMP
);

CREATE TABLE cloud_credentials (
    credential_id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID REFERENCES users(user_id) ON DELETE CASCADE,
    provider VARCHAR(50),
    encrypted_token TEXT,
    encrypted_crn TEXT,
    verified BOOLEAN DEFAULT FALSE,
    created_at TIMESTAMP DEFAULT NOW()
);

CREATE TABLE user_settings (
    user_id UUID PRIMARY KEY REFERENCES users(user_id) ON DELETE CASCADE,
    default_basis VARCHAR(50) DEFAULT 'sto-3g',
    default_method VARCHAR(50) DEFAULT 'VQE',
    preferences JSONB
);
```

#### **5.3 Use Prisma with Supabase**

```bash
# In your backend directory
npm install -g prisma
npm install @prisma/client

# Initialize Prisma
npx prisma init
```

**`prisma/schema.prisma`**:

```prisma
datasource db {
  provider = "postgresql"
  url      = env("USER_DATABASE_URL")  // Supabase URL
}

generator client {
  provider = "prisma-client-py"  // For Python
}

model User {
  user_id         String   @id @default(uuid())
  email           String   @unique
  name            String?
  hashed_password String
  institution     String?
  research_field  String?
  created_at      DateTime @default(now())
  last_login      DateTime?
  
  credentials     CloudCredential[]
  settings        UserSettings?
}

model CloudCredential {
  credential_id   String   @id @default(uuid())
  user_id         String
  provider        String
  encrypted_token String
  encrypted_crn   String?
  verified        Boolean  @default(false)
  created_at      DateTime @default(now())
  
  user            User     @relation(fields: [user_id], references: [user_id], onDelete: Cascade)
}

model UserSettings {
  user_id        String  @id
  default_basis  String  @default("sto-3g")
  default_method String  @default("VQE")
  preferences    Json?
  
  user           User    @relation(fields: [user_id], references: [user_id], onDelete: Cascade)
}
```

```bash
# Generate client
npx prisma generate

# Push schema to Supabase
npx prisma db push
```

---

### **PHASE 6: Frontend on Vercel**

#### **6.1 Create Next.js App**

```bash
# On your local machine
npx create-next-app@latest kanad-frontend
cd kanad-frontend

# Install dependencies
npm install @tanstack/react-query zustand axios socket.io-client
npm install shadcn-ui @radix-ui/react-dialog @radix-ui/react-dropdown-menu
npm install recharts three @react-three/fiber  # For visualizations
```

#### **6.2 Environment Variables (`.env.local`)**

```bash
NEXT_PUBLIC_API_URL=https://your-azure-vm-domain.com
NEXT_PUBLIC_WS_URL=wss://your-azure-vm-domain.com
```

#### **6.3 Deploy to Vercel**

```bash
# Install Vercel CLI
npm install -g vercel

# Login
vercel login

# Deploy
vercel

# Production deployment
vercel --prod
```

**Or via GitHub**:

1. Push code to GitHub
2. Go to vercel.com → Import Project
3. Connect GitHub repo
4. Auto-deploys on every push to main branch

---

## 🔒 Security Best Practices

### **1. Firewall Rules (Azure NSG)**

```bash
# Allow only necessary ports
- Port 22 (SSH) - Restrict to your IP only
- Port 80 (HTTP) - Redirect to HTTPS
- Port 443 (HTTPS) - Public
- Block all other ports
```

### **2. Environment Variables**

```bash
# Never commit .env files
echo ".env" >> .gitignore

# Use Azure Key Vault for production secrets (optional)
# https://azure.microsoft.com/en-us/services/key-vault/
```

### **3. Database Backups**

```bash
# Automatic PostgreSQL backups
sudo nano /etc/cron.daily/backup-kanad-db
```

```bash
#!/bin/bash
pg_dump -U kanad_user kanad_experiments > /home/azureuser/backups/kanad_$(date +\%Y\%m\%d).sql
# Upload to Azure Blob Storage or S3
```

---

## 📊 Monitoring & Logging

### **1. Setup Logging**

```python
# api/config.py
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('/var/log/kanad/api.log'),
        logging.StreamHandler()
    ]
)
```

### **2. Monitor Services**

```bash
# Check service status
systemctl status kanad-api
systemctl status kanad-worker
systemctl status redis
systemctl status postgresql

# View logs
journalctl -u kanad-api -f
journalctl -u kanad-worker -f
```

---

## 💡 Cost Optimization Tips

1. **Auto-shutdown during low usage**:
    
    ```bash
    # Azure Portal → VM → Auto-shutdown
    # Schedule shutdown at night if not running 24/7 jobs
    ```
    
2. **Use Spot Instances** (70-90% discount):
    
    - For non-critical workloads
    - Azure can reclaim VM with 30s notice
    - Great for development/testing
3. **Monitor spending**:
    
    - Azure Cost Management + Billing dashboard
    - Set budget alerts

---

## 🚀 Final Deployment Checklist

- [ ]  Azure VM created (Standard_D8s_v5)
- [ ]  Python 3.11+ installed
- [ ]  Kanad framework installed
- [ ]  PostgreSQL setup (experiments DB)
- [ ]  Redis setup (task queue)
- [ ]  FastAPI app running
- [ ]  Celery workers running
- [ ]  Nginx reverse proxy configured
- [ ]  HTTPS certificate installed
- [ ]  Supabase project created (user DB)
- [ ]  Prisma schema pushed
- [ ]  Frontend deployed to Vercel
- [ ]  Environment variables configured
- [ ]  Firewall rules configured
- [ ]  Backup script created

---

This architecture gives you **3 months of runtime** with your $1000 credits, with the flexibility to scale up to D16s_v5 if you need more power. The separation of user data (Supabase) and experiment data (VM) makes future migrations easy. 🎯