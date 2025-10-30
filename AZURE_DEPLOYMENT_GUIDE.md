# Azure Deployment Guide - Step by Step

## Prerequisites Setup

### 1. Install Azure CLI

```bash
# For Arch Linux
sudo pacman -S azure-cli

# For Ubuntu/Debian
curl -sL https://aka.ms/InstallAzureCLIDeb | sudo bash

# For macOS
brew install azure-cli
```

### 2. Login to Azure

```bash
# Login with your Microsoft account
az login

# This will open a browser window for authentication
# Follow the prompts to login with your Azure account
```

### 3. Verify Your Subscription

```bash
# List all subscriptions
az account list --output table

# Set the correct subscription (if you have multiple)
az account set --subscription "Azure for Students" # or your subscription name

# Verify current subscription
az account show
```

### 4. Check Your Credits

```bash
# View subscription details including remaining credits
az consumption usage list --query '[0]'

# Or check in Azure Portal:
# https://portal.azure.com → Cost Management → Credits
```

## Deployment Steps

### Option 1: Automated Deployment (Recommended)

Run the deployment script:

```bash
cd /home/mk/deeprealm/kanad
./deploy-azure.sh
```

This script will:
- ✅ Create all Azure resources
- ✅ Build and push Docker image
- ✅ Setup PostgreSQL database
- ✅ Deploy the API
- ✅ Configure environment variables
- ✅ Enable logging and monitoring

**Estimated time:** 15-20 minutes
**Estimated cost:** ~$395/month (P2V2 App Service + PostgreSQL)

### Option 2: Manual Deployment

If you prefer manual control, follow these steps:

#### Step 1: Create Resource Group

```bash
RESOURCE_GROUP="kanad-rg"
LOCATION="eastus"  # Or choose: westus, westeurope, etc.

az group create --name $RESOURCE_GROUP --location $LOCATION
```

#### Step 2: Create Container Registry

```bash
ACR_NAME="kanadregistry"  # Must be globally unique, lowercase

az acr create \
  --resource-group $RESOURCE_GROUP \
  --name $ACR_NAME \
  --sku Basic \
  --admin-enabled true
```

#### Step 3: Build and Push Docker Image

```bash
# Login to ACR
az acr login --name $ACR_NAME

# Build in Azure (recommended - faster)
cd /home/mk/deeprealm/kanad
az acr build \
  --registry $ACR_NAME \
  --image kanad-api:latest \
  --file Dockerfile \
  .
```

#### Step 4: Create PostgreSQL Database

```bash
POSTGRES_SERVER="kanad-postgres-server"  # Must be globally unique
POSTGRES_ADMIN="kanadmin"
POSTGRES_PASSWORD="YourSecurePassword123!"  # Change this!
POSTGRES_DB="kanad_db"

# Create server (takes 5-10 minutes)
az postgres flexible-server create \
  --resource-group $RESOURCE_GROUP \
  --name $POSTGRES_SERVER \
  --location $LOCATION \
  --admin-user $POSTGRES_ADMIN \
  --admin-password "$POSTGRES_PASSWORD" \
  --sku-name GP_Standard_D2s_v3 \
  --tier GeneralPurpose \
  --version 15 \
  --storage-size 128 \
  --public-access 0.0.0.0-255.255.255.255

# Create database
az postgres flexible-server db create \
  --resource-group $RESOURCE_GROUP \
  --server-name $POSTGRES_SERVER \
  --database-name $POSTGRES_DB

# Allow Azure services
az postgres flexible-server firewall-rule create \
  --resource-group $RESOURCE_GROUP \
  --name $POSTGRES_SERVER \
  --rule-name AllowAzureServices \
  --start-ip-address 0.0.0.0 \
  --end-ip-address 0.0.0.0
```

#### Step 5: Create App Service

```bash
APP_SERVICE_PLAN="kanad-plan"
WEB_APP_NAME="kanad-api"  # Must be globally unique

# Create plan with P2V2 (4 cores, 7GB RAM)
az appservice plan create \
  --name $APP_SERVICE_PLAN \
  --resource-group $RESOURCE_GROUP \
  --is-linux \
  --sku P2V2

# Create web app
az webapp create \
  --resource-group $RESOURCE_GROUP \
  --plan $APP_SERVICE_PLAN \
  --name $WEB_APP_NAME \
  --deployment-container-image-name ${ACR_NAME}.azurecr.io/kanad-api:latest
```

#### Step 6: Configure Container

```bash
# Get ACR credentials
ACR_USERNAME=$(az acr credential show --name $ACR_NAME --query username -o tsv)
ACR_PASSWORD=$(az acr credential show --name $ACR_NAME --query passwords[0].value -o tsv)

# Configure container
az webapp config container set \
  --name $WEB_APP_NAME \
  --resource-group $RESOURCE_GROUP \
  --docker-custom-image-name ${ACR_NAME}.azurecr.io/kanad-api:latest \
  --docker-registry-server-url https://${ACR_NAME}.azurecr.io \
  --docker-registry-server-user $ACR_USERNAME \
  --docker-registry-server-password $ACR_PASSWORD
```

#### Step 7: Configure Environment Variables

```bash
# Generate SECRET_KEY
SECRET_KEY=$(openssl rand -base64 48)

# Build DATABASE_URL
DATABASE_URL="postgresql://${POSTGRES_ADMIN}:${POSTGRES_PASSWORD}@${POSTGRES_SERVER}.postgres.database.azure.com:5432/${POSTGRES_DB}?sslmode=require"

# Set app settings
az webapp config appsettings set \
  --resource-group $RESOURCE_GROUP \
  --name $WEB_APP_NAME \
  --settings \
    DATABASE_URL="$DATABASE_URL" \
    SECRET_KEY="$SECRET_KEY" \
    ENVIRONMENT="production" \
    DEBUG="false" \
    CORS_ORIGINS="https://your-frontend.vercel.app" \
    WEBSITES_PORT="8000" \
    WEBSITES_ENABLE_APP_SERVICE_STORAGE="false"
```

#### Step 8: Add Additional Secrets via Azure Portal

1. Go to https://portal.azure.com
2. Navigate to: App Services → kanad-api → Configuration → Application settings
3. Add these secrets:

```
GOOGLE_CLIENT_ID=<your-google-client-id>
GOOGLE_CLIENT_SECRET=<your-google-client-secret>
GOOGLE_REDIRECT_URI=https://your-frontend.vercel.app/auth/callback/google

SMTP_HOST=smtp.gmail.com
SMTP_PORT=587
SMTP_USER=<your-email@gmail.com>
SMTP_PASSWORD=<your-gmail-app-password>
SMTP_FROM=noreply@kanad.app

ADMIN_PASSWORD=<secure-admin-password>

# Optional: Quantum backends
IBM_QUANTUM_TOKEN=<your-ibm-token>
BLUEQUBIT_API_KEY=<your-bluequbit-key>
```

4. Click "Save"

## Post-Deployment

### 1. Test the API

```bash
WEB_APP_URL="https://${WEB_APP_NAME}.azurewebsites.net"

# Health check
curl $WEB_APP_URL/api/health

# Should return: {"status": "ok", "version": "0.1.0"}
```

### 2. View Logs

```bash
# Stream logs in real-time
az webapp log tail \
  --resource-group $RESOURCE_GROUP \
  --name $WEB_APP_NAME

# Or view in Azure Portal:
# App Services → kanad-api → Monitoring → Log stream
```

### 3. Run Database Migrations

```bash
# SSH into the container
az webapp ssh \
  --resource-group $RESOURCE_GROUP \
  --name $WEB_APP_NAME

# Inside the container, run migrations if needed
# alembic upgrade head
```

### 4. Create Admin User

Use the admin panel at: `https://${WEB_APP_NAME}.azurewebsites.net/mkkisarkar`

### 5. Update Frontend

Update your frontend `.env` or Vercel environment variables:

```bash
NEXT_PUBLIC_API_URL=https://kanad-api.azurewebsites.net
```

## Monitoring & Management

### View Resource Usage

```bash
# CPU and memory metrics
az monitor metrics list \
  --resource /subscriptions/{subscription-id}/resourceGroups/$RESOURCE_GROUP/providers/Microsoft.Web/sites/$WEB_APP_NAME \
  --metric "CpuPercentage,MemoryPercentage" \
  --start-time 2024-01-01T00:00:00Z

# Or use Azure Portal:
# App Services → kanad-api → Monitoring → Metrics
```

### Scale Up/Down

```bash
# Scale to P3V2 (8 cores, 14GB RAM) for heavy workloads
az appservice plan update \
  --resource-group $RESOURCE_GROUP \
  --name $APP_SERVICE_PLAN \
  --sku P3V2

# Scale back to P1V2 to save credits
az appservice plan update \
  --resource-group $RESOURCE_GROUP \
  --name $APP_SERVICE_PLAN \
  --sku P1V2
```

### Enable Auto-Scaling

```bash
# Auto-scale based on CPU usage
az monitor autoscale create \
  --resource-group $RESOURCE_GROUP \
  --resource $APP_SERVICE_PLAN \
  --resource-type Microsoft.Web/serverfarms \
  --name autoscale-kanad \
  --min-count 1 \
  --max-count 5 \
  --count 2

az monitor autoscale rule create \
  --resource-group $RESOURCE_GROUP \
  --autoscale-name autoscale-kanad \
  --condition "CpuPercentage > 70 avg 5m" \
  --scale out 1
```

## Cost Management

### View Current Costs

```bash
# View cost analysis
az consumption usage list \
  --start-date $(date -d "30 days ago" +%Y-%m-%d) \
  --end-date $(date +%Y-%m-%d)
```

### Set Budget Alerts

1. Go to Azure Portal
2. Navigate to: Cost Management + Billing → Budgets
3. Create budget:
   - Name: "Kanad Monthly Budget"
   - Amount: $500 (adjust based on your credits)
   - Alert at: 80%, 90%, 100%

### Estimated Monthly Costs

| Resource | SKU | Cost |
|----------|-----|------|
| App Service | P2V2 | ~$240/month |
| PostgreSQL | GP_Standard_D2s_v3 | ~$150/month |
| Container Registry | Basic | ~$5/month |
| **Total** | | **~$395/month** |

**With $1000 credits:** ~2.5 months
**With $5000 credits:** ~12 months

## Troubleshooting

### Container Won't Start

```bash
# Check logs
az webapp log tail --name $WEB_APP_NAME --resource-group $RESOURCE_GROUP

# Common issues:
# - DATABASE_URL not set correctly
# - Missing environment variables
# - Port 8000 not exposed (check Dockerfile)
```

### Database Connection Errors

```bash
# Test connection from local machine
psql "postgresql://${POSTGRES_ADMIN}:${POSTGRES_PASSWORD}@${POSTGRES_SERVER}.postgres.database.azure.com:5432/${POSTGRES_DB}?sslmode=require"

# Check firewall rules
az postgres flexible-server firewall-rule list \
  --resource-group $RESOURCE_GROUP \
  --name $POSTGRES_SERVER
```

### High CPU Usage

```bash
# Check if computations are running
curl https://${WEB_APP_NAME}.azurewebsites.net/api/admin/experiments/live

# View quota usage
curl https://${WEB_APP_NAME}.azurewebsites.net/api/users/me/quota \
  -H "Authorization: Bearer YOUR_TOKEN"
```

## Cleanup (Delete All Resources)

**WARNING:** This will delete everything!

```bash
az group delete --name $RESOURCE_GROUP --yes --no-wait
```

## Next Steps

1. ✅ Run `./deploy-azure.sh`
2. ✅ Set additional environment variables in Azure Portal
3. ✅ Test API endpoints
4. ✅ Deploy frontend to Vercel
5. ✅ Update frontend with Azure API URL
6. ✅ Test end-to-end workflow
7. ✅ Set up monitoring and alerts
8. ✅ Configure backup strategy

## Support

- Azure Documentation: https://docs.microsoft.com/azure
- Kanad GitHub: https://github.com/yourusername/kanad
- Azure Support: https://portal.azure.com → Help + support
