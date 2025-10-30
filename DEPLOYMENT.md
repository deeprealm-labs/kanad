# Kanad Deployment Guide

Complete guide for deploying the Kanad quantum chemistry platform to production environments.

## Architecture Overview

```
┌─────────────────┐         ┌──────────────────┐         ┌─────────────────┐
│                 │         │                  │         │                 │
│  Vercel         │────────▶│  Azure Web App   │────────▶│  Azure Database │
│  (Frontend)     │  HTTPS  │  (Backend API)   │         │  for PostgreSQL │
│  Next.js 15.5.4 │         │  FastAPI         │         │                 │
└─────────────────┘         └──────────────────┘         └─────────────────┘
```

**Stack:**
- **Frontend**: Next.js 15.5.4 on Vercel
- **Backend**: FastAPI on Azure Web Apps (Docker container)
- **Database**: Azure Database for PostgreSQL (or Azure SQL)
- **Storage**: SQLite for experiments (Azure Blob Storage optional)

---

## Table of Contents

1. [Prerequisites](#prerequisites)
2. [Frontend Deployment - Vercel](#frontend-deployment---vercel)
3. [Backend Deployment - Azure](#backend-deployment---azure)
4. [Database Setup - Azure PostgreSQL](#database-setup---azure-postgresql)
5. [Environment Variables](#environment-variables)
6. [Docker Deployment](#docker-deployment)
7. [Post-Deployment](#post-deployment)
8. [Monitoring & Troubleshooting](#monitoring--troubleshooting)

---

## Prerequisites

### Required Accounts
- [x] **Vercel Account** (free tier available)
- [x] **Azure Account** (student credits or pay-as-you-go)
- [x] **GitHub Account** (for CI/CD)

### Required Tools
```bash
# Install Azure CLI
curl -sL https://aka.ms/InstallAzureCLIDeb | sudo bash

# Install Vercel CLI
npm install -g vercel

# Install Docker (optional)
sudo apt-get install docker.io docker-compose
```

### Domain & SSL
- Custom domain (optional but recommended)
- SSL certificates (provided by Vercel & Azure)

---

## Frontend Deployment - Vercel

### 1. Prepare the Frontend

Ensure the build works locally:
```bash
cd web
npm install
npm run build
```

### 2. Configure Environment Variables

Create `.env.production` in the `web/` directory:
```env
NEXT_PUBLIC_API_URL=https://kanad-api.azurewebsites.net
```

### 3. Deploy to Vercel

#### Option A: Vercel CLI
```bash
cd web
vercel login
vercel --prod
```

#### Option B: Vercel Dashboard

1. Go to [vercel.com/new](https://vercel.com/new)
2. Import your GitHub repository
3. Set **Root Directory**: `web`
4. **Framework Preset**: Next.js
5. Add environment variables:
   - `NEXT_PUBLIC_API_URL`: `https://kanad-api.azurewebsites.net`
6. Click **Deploy**

### 4. Configure Custom Domain (Optional)

1. Go to **Project Settings** → **Domains**
2. Add your custom domain (e.g., `kanad.app`)
3. Update DNS records as instructed
4. SSL certificate will be auto-provisioned

### 5. Update vercel.json

The [`web/vercel.json`](web/vercel.json) file is pre-configured with:
- Security headers
- API rewrites
- Build settings

Update the `rewrites` section with your actual backend URL:
```json
{
  "rewrites": [
    {
      "source": "/api/:path*",
      "destination": "https://kanad-api.azurewebsites.net/api/:path*"
    }
  ]
}
```

---

## Backend Deployment - Azure

### 1. Create Azure Resources

#### Login to Azure
```bash
az login
```

#### Create Resource Group
```bash
az group create \
  --name kanad-rg \
  --location eastus
```

#### Create Container Registry (ACR)
```bash
az acr create \
  --resource-group kanad-rg \
  --name kanadregistry \
  --sku Basic \
  --admin-enabled true
```

### 2. Build and Push Docker Image

#### Build the Image
```bash
# Navigate to project root
cd /home/mk/deeprealm/kanad

# Build Docker image
docker build -t kanad-api:latest .
```

#### Tag and Push to ACR
```bash
# Login to ACR
az acr login --name kanadregistry

# Tag the image
docker tag kanad-api:latest kanadregistry.azurecr.io/kanad-api:latest

# Push to ACR
docker push kanadregistry.azurecr.io/kanad-api:latest
```

### 3. Create App Service Plan

```bash
az appservice plan create \
  --name kanad-plan \
  --resource-group kanad-rg \
  --is-linux \
  --sku B1  # Change to P1V2 for production
```

### 4. Create Web App

```bash
az webapp create \
  --resource-group kanad-rg \
  --plan kanad-plan \
  --name kanad-api \
  --deployment-container-image-name kanadregistry.azurecr.io/kanad-api:latest
```

### 5. Configure Web App

#### Set Container Registry Credentials
```bash
# Get ACR credentials
ACR_USERNAME=$(az acr credential show --name kanadregistry --query username -o tsv)
ACR_PASSWORD=$(az acr credential show --name kanadregistry --query passwords[0].value -o tsv)

# Configure Web App to use ACR
az webapp config container set \
  --name kanad-api \
  --resource-group kanad-rg \
  --docker-custom-image-name kanadregistry.azurecr.io/kanad-api:latest \
  --docker-registry-server-url https://kanadregistry.azurecr.io \
  --docker-registry-server-user $ACR_USERNAME \
  --docker-registry-server-password $ACR_PASSWORD
```

#### Configure Application Settings
```bash
az webapp config appsettings set \
  --resource-group kanad-rg \
  --name kanad-api \
  --settings \
    DATABASE_URL="postgresql://user:password@server.postgres.database.azure.com:5432/kanad_db" \
    SECRET_KEY="your-secret-key-here" \
    ENVIRONMENT="production" \
    CORS_ORIGINS="https://your-frontend.vercel.app" \
    GOOGLE_CLIENT_ID="your-google-client-id" \
    GOOGLE_CLIENT_SECRET="your-google-client-secret"
```

See [Environment Variables](#environment-variables) section for complete list.

### 6. Enable Continuous Deployment

```bash
az webapp deployment container config \
  --name kanad-api \
  --resource-group kanad-rg \
  --enable-cd true
```

### 7. Configure CORS

```bash
az webapp cors add \
  --resource-group kanad-rg \
  --name kanad-api \
  --allowed-origins https://your-frontend.vercel.app http://localhost:3000
```

---

## Database Setup - Azure PostgreSQL

### 1. Create PostgreSQL Server

```bash
az postgres flexible-server create \
  --resource-group kanad-rg \
  --name kanad-postgres \
  --location eastus \
  --admin-user kanaduser \
  --admin-password "YourSecurePassword123!" \
  --sku-name Standard_B1ms \
  --tier Burstable \
  --version 15 \
  --storage-size 32 \
  --public-access 0.0.0.0-255.255.255.255
```

### 2. Create Database

```bash
az postgres flexible-server db create \
  --resource-group kanad-rg \
  --server-name kanad-postgres \
  --database-name kanad_db
```

### 3. Configure Firewall

```bash
# Allow Azure services
az postgres flexible-server firewall-rule create \
  --resource-group kanad-rg \
  --name kanad-postgres \
  --rule-name AllowAzureServices \
  --start-ip-address 0.0.0.0 \
  --end-ip-address 0.0.0.0

# Allow your IP (for database management)
az postgres flexible-server firewall-rule create \
  --resource-group kanad-rg \
  --name kanad-postgres \
  --rule-name AllowMyIP \
  --start-ip-address YOUR_IP_ADDRESS \
  --end-ip-address YOUR_IP_ADDRESS
```

### 4. Get Connection String

```bash
az postgres flexible-server show-connection-string \
  --server-name kanad-postgres \
  --admin-user kanaduser \
  --database-name kanad_db
```

Update the `DATABASE_URL` in your Web App settings with this connection string.

---

## Environment Variables

### Backend (.env)

Create a `.env` file in the project root (for local development) or set via Azure App Settings (for production):

```env
# ============================================================================
# Kanad Backend Environment Variables
# ============================================================================

# Database
DATABASE_URL=postgresql://user:password@host:5432/kanad_db

# API Configuration
API_VERSION=0.1.0
ENVIRONMENT=production
DEBUG=false

# Security
SECRET_KEY=your-secret-key-minimum-32-characters-long
ACCESS_TOKEN_EXPIRE_MINUTES=60
REFRESH_TOKEN_EXPIRE_DAYS=7

# CORS - Frontend URL
CORS_ORIGINS=https://your-frontend.vercel.app,https://kanad.app

# OAuth - Google
GOOGLE_CLIENT_ID=your-google-client-id.apps.googleusercontent.com
GOOGLE_CLIENT_SECRET=your-google-client-secret
GOOGLE_REDIRECT_URI=https://your-frontend.vercel.app/auth/callback/google

# Email (SMTP)
SMTP_HOST=smtp.gmail.com
SMTP_PORT=587
SMTP_USER=your-email@gmail.com
SMTP_PASSWORD=your-app-password
SMTP_FROM=noreply@kanad.app

# Admin Panel Password
ADMIN_PASSWORD=your-secure-admin-password

# Quantum Backends (Optional)
IBM_QUANTUM_TOKEN=your-ibm-quantum-token
BLUEQUBIT_API_KEY=your-bluequbit-api-key
```

### Frontend (.env.local)

Create `web/.env.local` for local development or set in Vercel dashboard:

```env
# API Endpoint
NEXT_PUBLIC_API_URL=https://kanad-api.azurewebsites.net

# Google OAuth (must match backend)
NEXT_PUBLIC_GOOGLE_CLIENT_ID=your-google-client-id.apps.googleusercontent.com
```

---

## Docker Deployment

### Local Testing

```bash
# Start all services
docker-compose up -d

# Check logs
docker-compose logs -f api

# Stop services
docker-compose down
```

### Manual Docker Deployment

```bash
# Build
docker build -t kanad-api:latest .

# Run with environment file
docker run -d \
  --name kanad-api \
  -p 8000:8000 \
  --env-file .env \
  kanad-api:latest

# Check logs
docker logs -f kanad-api
```

---

## Post-Deployment

### 1. Database Migrations

Run migrations on the Azure Web App:

```bash
# SSH into the container
az webapp ssh --resource-group kanad-rg --name kanad-api

# Run migrations (if using Alembic)
alembic upgrade head
```

### 2. Create Admin User

```bash
# Access the application
curl -X POST https://kanad-api.azurewebsites.net/api/admin/setup \
  -H "Content-Type: application/json" \
  -d '{"password": "your-admin-password"}'
```

### 3. Verify Health

```bash
curl https://kanad-api.azurewebsites.net/api/health
```

### 4. Test Endpoints

```bash
# Test public endpoint
curl https://kanad-api.azurewebsites.net/api/settings/defaults

# Test with authentication
curl https://kanad-api.azurewebsites.net/api/users/me \
  -H "Authorization: Bearer YOUR_TOKEN"
```

### 5. Update Frontend API URL

Ensure Vercel environment variable is set correctly:
```bash
vercel env add NEXT_PUBLIC_API_URL production
# Enter: https://kanad-api.azurewebsites.net
```

---

## Monitoring & Troubleshooting

### Azure Application Insights

Enable monitoring:
```bash
az monitor app-insights component create \
  --app kanad-insights \
  --location eastus \
  --resource-group kanad-rg
```

### View Logs

```bash
# Stream logs
az webapp log tail \
  --resource-group kanad-rg \
  --name kanad-api

# Download logs
az webapp log download \
  --resource-group kanad-rg \
  --name kanad-api
```

### Common Issues

#### 1. Container Failed to Start
```bash
# Check Docker logs
az webapp log tail --resource-group kanad-rg --name kanad-api

# Common fixes:
# - Verify environment variables are set
# - Check DATABASE_URL format
# - Ensure port 8000 is exposed
```

#### 2. Database Connection Errors
```bash
# Test connection from Web App
az webapp ssh --resource-group kanad-rg --name kanad-api
psql $DATABASE_URL

# Verify:
# - Firewall rules allow Azure services
# - Connection string is correct
# - Database exists
```

#### 3. CORS Errors
```bash
# Update CORS settings
az webapp cors add \
  --resource-group kanad-rg \
  --name kanad-api \
  --allowed-origins https://your-new-domain.com
```

#### 4. Frontend Can't Reach Backend
- Verify `NEXT_PUBLIC_API_URL` in Vercel
- Check Azure Web App is running: `az webapp show --name kanad-api --resource-group kanad-rg`
- Test backend health: `curl https://kanad-api.azurewebsites.net/api/health`

### Performance Optimization

1. **Upgrade App Service Plan**: B1 → P1V2 for better performance
2. **Enable Auto-scaling**: Configure in Azure Portal
3. **Use Azure CDN**: For static assets
4. **Enable Redis Cache**: For session management

---

## Cost Estimates

### Minimal Setup (for testing)
- **Vercel**: Free tier
- **Azure B1 App Service**: ~$13/month
- **Azure PostgreSQL Burstable**: ~$12/month
- **Total**: ~$25/month

### Production Setup
- **Vercel Pro**: $20/month
- **Azure P1V2 App Service**: ~$120/month
- **Azure PostgreSQL Standard**: ~$80/month
- **Azure Container Registry**: ~$5/month
- **Total**: ~$225/month

---

## Security Checklist

- [ ] Use HTTPS for all communications
- [ ] Set strong `SECRET_KEY` (32+ characters)
- [ ] Enable Azure Web App authentication
- [ ] Configure PostgreSQL firewall rules
- [ ] Use environment variables for secrets (never commit)
- [ ] Enable Azure DDoS protection
- [ ] Set up rate limiting on API endpoints
- [ ] Regular security updates via Docker image rebuilds
- [ ] Configure backup strategy for databases

---

## CI/CD Setup (GitHub Actions)

See `.github/workflows/` for automated deployment pipelines:
- `deploy-frontend.yml` - Auto-deploy to Vercel on push
- `deploy-backend.yml` - Auto-deploy to Azure on push

---

## Support

For deployment issues:
- GitHub Issues: https://github.com/yourusername/kanad/issues
- Azure Support: https://azure.microsoft.com/support/
- Vercel Support: https://vercel.com/support
