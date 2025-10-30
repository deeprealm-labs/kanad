# Kanad Platform - Deployment Guide
## Authentication & Admin System Setup

---

## Table of Contents
1. [Prerequisites](#prerequisites)
2. [Local Development Setup](#local-development-setup)
3. [Azure PostgreSQL Setup](#azure-postgresql-setup)
4. [Environment Configuration](#environment-configuration)
5. [Database Migration](#database-migration)
6. [Creating First Admin User](#creating-first-admin-user)
7. [Azure Backend Deployment](#azure-backend-deployment)
8. [Vercel Frontend Deployment](#vercel-frontend-deployment)
9. [Testing Authentication](#testing-authentication)
10. [Troubleshooting](#troubleshooting)

---

## Prerequisites

### Required Services
- ✅ Azure Account (for PostgreSQL and App Service)
- ✅ Vercel Account (for frontend hosting)
- ✅ Google Cloud Console (for OAuth)
- ✅ SMTP Server (for email verification)

### Local Development
- Python 3.10+
- Node.js 18+
- PostgreSQL 14+ (for local testing)

---

## Local Development Setup

### 1. Install Backend Dependencies

```bash
cd /home/mk/deeprealm/kanad
pip install -r requirements_auth.txt
```

### 2. Install Frontend Dependencies

```bash
cd web
npm install
```

### 3. Setup PostgreSQL Locally

```bash
# Install PostgreSQL (Ubuntu/Debian)
sudo apt update
sudo apt install postgresql postgresql-contrib

# Create database
sudo -u postgres psql
CREATE DATABASE kanad_db;
CREATE USER kanad_user WITH PASSWORD 'your_secure_password';
GRANT ALL PRIVILEGES ON DATABASE kanad_db TO kanad_user;
\q
```

---

## Azure PostgreSQL Setup

### 1. Create Azure PostgreSQL Flexible Server

```bash
# Login to Azure
az login

# Create resource group
az group create --name kanad-rg --location eastus

# Create PostgreSQL Flexible Server
az postgres flexible-server create \
  --name kanad-db-server \
  --resource-group kanad-rg \
  --location eastus \
  --admin-user kanadmin \
  --admin-password 'YourSecurePassword123!' \
  --sku-name Standard_B2s \
  --tier Burstable \
  --storage-size 32 \
  --version 14

# Create database
az postgres flexible-server db create \
  --resource-group kanad-rg \
  --server-name kanad-db-server \
  --database-name kanad_db

# Configure firewall (allow Azure services)
az postgres flexible-server firewall-rule create \
  --resource-group kanad-rg \
  --name kanad-db-server \
  --rule-name AllowAzureServices \
  --start-ip-address 0.0.0.0 \
  --end-ip-address 0.0.0.0
```

### 2. Get Connection String

```bash
# Format: postgresql://username:password@host:port/database
postgresql://kanadmin:YourSecurePassword123!@kanad-db-server.postgres.database.azure.com:5432/kanad_db
```

---

## Environment Configuration

### 1. Backend Environment Variables

Create `.env` file in project root:

```env
# Database
DATABASE_URL=postgresql://kanadmin:YourSecurePassword123!@kanad-db-server.postgres.database.azure.com:5432/kanad_db

# JWT
JWT_SECRET_KEY=your-super-secret-jwt-key-min-32-chars-random-string-here
JWT_ALGORITHM=HS256

# Google OAuth
GOOGLE_CLIENT_ID=your-google-client-id.apps.googleusercontent.com
GOOGLE_CLIENT_SECRET=your-google-client-secret
GOOGLE_REDIRECT_URI=https://your-domain.com/auth/google/callback

# SMTP Email
SMTP_HOST=smtp.gmail.com
SMTP_PORT=587
SMTP_USER=your-email@gmail.com
SMTP_PASSWORD=your-app-specific-password
SMTP_FROM_EMAIL=noreply@kanad.com
SMTP_FROM_NAME=Kanad Platform

# Application
ENVIRONMENT=production
DEBUG=false
```

### 2. Generate Secure JWT Secret

```python
import secrets
print(secrets.token_urlsafe(32))
# Copy output to JWT_SECRET_KEY
```

### 3. Setup Google OAuth

1. Go to [Google Cloud Console](https://console.cloud.google.com/)
2. Create new project: "Kanad Platform"
3. Enable **Google+ API**
4. Create OAuth 2.0 credentials:
   - Application type: Web application
   - Authorized redirect URIs:
     - `http://localhost:3000/auth/google/callback` (dev)
     - `https://kanad.vercel.app/auth/google/callback` (prod)
5. Copy Client ID and Client Secret to `.env`

### 4. Setup Gmail SMTP

1. Enable 2-factor authentication on Gmail
2. Generate App Password:
   - Google Account → Security → 2-Step Verification → App passwords
   - Select "Mail" and "Other" (Kanad Platform)
3. Copy generated password to `SMTP_PASSWORD`

---

## Database Migration

### 1. Run Migrations

```bash
cd /home/mk/deeprealm/kanad

# Set DATABASE_URL environment variable
export DATABASE_URL="postgresql://kanadmin:password@kanad-db-server.postgres.database.azure.com:5432/kanad_db"

# Run migrations
python migrations/run_migrations.py migrate

# Check status
python migrations/run_migrations.py status
```

Expected output:
```
✓ Connected to PostgreSQL database
✓ Migration tracking table ready
✓ Found 0 previously executed migrations
✓ Found 2 migration files

→ Executing 2 pending migrations...

Executing migration: 001_create_auth_tables.sql
✓ Migration 001_create_auth_tables.sql completed successfully

Executing migration: 002_migrate_existing_data.sql
✓ Migration 002_migrate_existing_data.sql completed successfully

✓ All migrations completed successfully!
```

---

## Creating First Admin User

### Method 1: Using Python Script

Create `create_admin.py`:

```python
import os
from sqlalchemy.orm import Session
from api.core.database_postgres import SessionLocal, create_admin_user

# Set environment
os.environ["DATABASE_URL"] = "your-connection-string"

# Create admin
db = SessionLocal()
try:
    admin = create_admin_user(
        db,
        email="admin@kanad.com",
        password="Admin123!SecurePassword",
        full_name="Admin User"
    )
    print(f"✓ Admin created: {admin.email}")
finally:
    db.close()
```

Run:
```bash
python create_admin.py
```

### Method 2: Using SQL

```sql
-- Insert admin user directly
INSERT INTO users (
    email,
    password_hash,
    full_name,
    role,
    is_verified,
    is_active,
    created_at
) VALUES (
    'admin@kanad.com',
    '$2b$12$hashed_password_here',  -- Generate using bcrypt
    'Admin User',
    'admin',
    true,
    true,
    NOW()
);
```

---

## Azure Backend Deployment

### Option 1: Azure App Service (Recommended for ease)

```bash
# Create App Service plan
az appservice plan create \
  --name kanad-plan \
  --resource-group kanad-rg \
  --sku P2V3 \
  --is-linux

# Create Web App
az webapp create \
  --name kanad-api \
  --resource-group kanad-rg \
  --plan kanad-plan \
  --runtime "PYTHON:3.10"

# Configure environment variables
az webapp config appsettings set \
  --name kanad-api \
  --resource-group kanad-rg \
  --settings \
    DATABASE_URL="postgresql://..." \
    JWT_SECRET_KEY="..." \
    GOOGLE_CLIENT_ID="..." \
    # ... other variables

# Deploy code
cd /home/mk/deeprealm/kanad
zip -r kanad.zip . -x "*.git*" "web/*" "*.pyc"
az webapp deployment source config-zip \
  --name kanad-api \
  --resource-group kanad-rg \
  --src kanad.zip
```

### Option 2: Azure Container Instances (Recommended for high compute)

Create `Dockerfile`:

```dockerfile
FROM python:3.10-slim

WORKDIR /app

COPY requirements.txt requirements_auth.txt ./
RUN pip install --no-cache-dir -r requirements.txt -r requirements_auth.txt

COPY . .

EXPOSE 8000

CMD ["uvicorn", "api.main:app", "--host", "0.0.0.0", "--port", "8000"]
```

Deploy:

```bash
# Build and push to Azure Container Registry
az acr create --name kanadregistry --resource-group kanad-rg --sku Basic
az acr build --registry kanadregistry --image kanad-api:latest .

# Create container instance
az container create \
  --name kanad-api \
  --resource-group kanad-rg \
  --image kanadregistry.azurecr.io/kanad-api:latest \
  --cpu 4 \
  --memory 8 \
  --registry-login-server kanadregistry.azurecr.io \
  --registry-username $(az acr credential show --name kanadregistry --query username -o tsv) \
  --registry-password $(az acr credential show --name kanadregistry --query passwords[0].value -o tsv) \
  --ports 8000 \
  --environment-variables \
    DATABASE_URL="..." \
    JWT_SECRET_KEY="..." \
  --secure-environment-variables \
    SMTP_PASSWORD="..."
```

---

## Vercel Frontend Deployment

### 1. Update Frontend Environment

Create `web/.env.production`:

```env
NEXT_PUBLIC_API_URL=https://kanad-api.azurewebsites.net
NEXT_PUBLIC_GOOGLE_CLIENT_ID=your-google-client-id
```

### 2. Deploy to Vercel

```bash
cd web

# Install Vercel CLI
npm i -g vercel

# Login
vercel login

# Deploy
vercel --prod

# Set environment variables via Vercel dashboard
# Project Settings → Environment Variables → Add:
#   - NEXT_PUBLIC_API_URL
#   - NEXT_PUBLIC_GOOGLE_CLIENT_ID
```

---

## Testing Authentication

### 1. Test Registration

```bash
curl -X POST https://kanad-api.azurewebsites.net/api/auth/register \
  -H "Content-Type: application/json" \
  -d '{
    "email": "test@example.com",
    "password": "Test123!Password",
    "full_name": "Test User",
    "access_key": "KANAD-xxxxxxxxxxxxx"
  }'
```

### 2. Test Login

```bash
curl -X POST https://kanad-api.azurewebsites.net/api/auth/login \
  -H "Content-Type: application/json" \
  -d '{
    "email": "test@example.com",
    "password": "Test123!Password"
  }'
```

### 3. Test Admin Access

```bash
# Login as admin
TOKEN="your-admin-jwt-token"

# Generate access key
curl -X POST https://kanad-api.azurewebsites.net/api/admin/keys \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/json" \
  -d '{
    "description": "Beta tester key",
    "max_uses": 10,
    "expires_in_days": 30
  }'
```

---

## Troubleshooting

### Issue: Database connection fails

**Solution:**
```bash
# Check firewall rules
az postgres flexible-server firewall-rule list \
  --resource-group kanad-rg \
  --name kanad-db-server

# Add your IP
az postgres flexible-server firewall-rule create \
  --resource-group kanad-rg \
  --name kanad-db-server \
  --rule-name AllowMyIP \
  --start-ip-address YOUR_IP \
  --end-ip-address YOUR_IP
```

### Issue: Migrations fail

**Solution:**
```bash
# Check connection
psql "postgresql://kanadmin:password@kanad-db-server.postgres.database.azure.com:5432/kanad_db"

# Rollback migration
python migrations/run_migrations.py rollback --migration 002_migrate_existing_data.sql

# Re-run
python migrations/run_migrations.py migrate
```

### Issue: SMTP emails not sending

**Solution:**
- Verify Gmail app password is correct
- Check firewall allows port 587
- Enable "Less secure app access" (if not using app password)
- Check logs: Emails are printed to console in development

### Issue: Google OAuth fails

**Solution:**
- Verify redirect URI matches exactly
- Check Google Console → APIs & Services → Credentials
- Enable Google+ API
- Clear browser cookies and try again

---

## Production Checklist

- [ ] PostgreSQL database provisioned and migrated
- [ ] Admin user created
- [ ] First access key generated
- [ ] Environment variables configured
- [ ] Google OAuth configured
- [ ] SMTP email working
- [ ] Backend deployed to Azure
- [ ] Frontend deployed to Vercel
- [ ] HTTPS enabled on both
- [ ] CORS configured correctly
- [ ] Firewall rules set
- [ ] Monitoring enabled
- [ ] Backup strategy in place

---

## Next Steps

1. ✅ Generate first access keys via admin panel
2. ✅ Test complete registration → verification → login flow
3. ✅ Invite beta users with access keys
4. ✅ Monitor live experiments in admin dashboard
5. ✅ Review usage statistics weekly

---

## Security Best Practices

### Required for Production:
1. **Use strong passwords**: Minimum 12 characters for admin
2. **Rotate JWT secrets**: Every 90 days
3. **Enable SSL/TLS**: HTTPS only
4. **Rate limiting**: Prevent brute force attacks
5. **Monitor logs**: Check for suspicious activity
6. **Backup database**: Daily automated backups
7. **Update dependencies**: Monthly security patches
8. **Review access keys**: Deactivate unused keys

### Recommended:
- Enable Azure AD authentication for admin
- Set up alerts for failed login attempts
- Implement IP whitelisting for admin routes
- Use Azure Key Vault for secrets
- Enable database encryption at rest

---

## Support

For issues or questions:
- Check logs: `az webapp log tail --name kanad-api --resource-group kanad-rg`
- Review API docs: `https://kanad-api.azurewebsites.net/docs`
- Contact: admin@kanad.com

---

**Deployment completed!** 🚀

Your Kanad Platform is now live with full authentication and admin capabilities.
