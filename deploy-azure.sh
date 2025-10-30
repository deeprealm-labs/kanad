#!/bin/bash
# ============================================================================
# Azure Deployment Script for Kanad Backend
# ============================================================================
# This script deploys the Kanad backend API to Azure with compute-optimized settings
#
# Prerequisites:
# - Azure CLI installed (`az`)
# - Docker installed
# - Logged into Azure (`az login`)
#
# Usage: ./deploy-azure.sh
# ============================================================================

set -e  # Exit on error

# ============================================================================
# Configuration
# ============================================================================

# Resource names (customize these)
RESOURCE_GROUP="kanad-rg"
LOCATION="eastus"  # Change if needed
ACR_NAME="kanadregistry"  # Must be globally unique, lowercase only
POSTGRES_SERVER="kanad-postgres-server"  # Must be globally unique
POSTGRES_ADMIN="kanadmin"
POSTGRES_DB="kanad_db"
APP_SERVICE_PLAN="kanad-plan"
WEB_APP_NAME="kanad-api"  # Must be globally unique

# Compute tier - HPC-level for Quantum Chemistry Calculations
# Options for Memory-Optimized (E-series - for Hamiltonian construction, VQE, SQD):
# - E2_v3: 2 vCores, 16GB RAM - ~$150/month (minimum recommended)
# - E4_v3: 4 vCores, 32GB RAM - ~$300/month (good for medium molecules)
# - E8_v3: 8 vCores, 64GB RAM - ~$600/month (large molecules, multiple jobs)
# - E16_v3: 16 vCores, 128GB RAM - ~$1200/month (HPC workloads)
#
# Options for Compute-Optimized (D-series - balanced compute/memory):
# - D4s_v3: 4 vCores, 16GB RAM - ~$200/month
# - D8s_v3: 8 vCores, 32GB RAM - ~$400/month
# - D16s_v3: 16 vCores, 64GB RAM - ~$800/month

# Recommended for Kanad (VQE/SQD Hamiltonian construction):
# E4_v3 provides 32GB RAM for quantum chemistry calculations
# Good balance: handles medium-large molecules, multiple concurrent jobs
APP_SERVICE_SKU="E4_v3"  # 4 vCores, 32GB RAM - HPC-level for quantum simulations

# PostgreSQL tier - Upgraded to match compute tier
# Options:
# - Burstable: $12/month (B_Standard_B1ms) - for testing
# - GeneralPurpose: $80/month (GP_Standard_D2s_v3) - 2 vCores, 8GB RAM
# - GeneralPurpose: $150/month (GP_Standard_D4s_v3) - 4 vCores, 16GB RAM
# - MemoryOptimized: $300/month (MO_Standard_E4s_v3) - 4 vCores, 32GB RAM
POSTGRES_SKU="GP_Standard_D4s_v3"  # 4 vCores, 16GB RAM - handles large molecular data

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# ============================================================================
# Helper Functions
# ============================================================================

print_header() {
    echo -e "\n${BLUE}===================================================================${NC}"
    echo -e "${BLUE}  $1${NC}"
    echo -e "${BLUE}===================================================================${NC}\n"
}

print_success() {
    echo -e "${GREEN}✅ $1${NC}"
}

print_warning() {
    echo -e "${YELLOW}⚠️  $1${NC}"
}

print_error() {
    echo -e "${RED}❌ $1${NC}"
}

# ============================================================================
# Pre-flight Checks
# ============================================================================

print_header "Pre-flight Checks"

# Check if Azure CLI is installed
if ! command -v az &> /dev/null; then
    print_error "Azure CLI not found. Please install it first:"
    echo "  curl -sL https://aka.ms/InstallAzureCLIDeb | sudo bash"
    exit 1
fi
print_success "Azure CLI found"

# Check if Docker is installed
if ! command -v docker &> /dev/null; then
    print_error "Docker not found. Please install it first:"
    echo "  sudo apt-get install docker.io"
    exit 1
fi
print_success "Docker found"

# Check if logged into Azure
if ! az account show &> /dev/null; then
    print_error "Not logged into Azure. Please run: az login"
    exit 1
fi
print_success "Logged into Azure"

# Show current subscription
SUBSCRIPTION=$(az account show --query name -o tsv)
print_success "Using subscription: $SUBSCRIPTION"

# ============================================================================
# Step 1: Create Resource Group
# ============================================================================

print_header "Step 1: Creating Resource Group"

if az group exists --name $RESOURCE_GROUP | grep -q "true"; then
    print_warning "Resource group '$RESOURCE_GROUP' already exists"
else
    az group create \
        --name $RESOURCE_GROUP \
        --location $LOCATION
    print_success "Resource group created: $RESOURCE_GROUP"
fi

# ============================================================================
# Step 2: Create Azure Container Registry (ACR)
# ============================================================================

print_header "Step 2: Creating Azure Container Registry"

if az acr show --name $ACR_NAME --resource-group $RESOURCE_GROUP &> /dev/null; then
    print_warning "ACR '$ACR_NAME' already exists"
else
    az acr create \
        --resource-group $RESOURCE_GROUP \
        --name $ACR_NAME \
        --sku Basic \
        --admin-enabled true
    print_success "ACR created: $ACR_NAME"
fi

# ============================================================================
# Step 3: Build and Push Docker Image
# ============================================================================

print_header "Step 3: Building and Pushing Docker Image"

# Login to ACR
az acr login --name $ACR_NAME
print_success "Logged into ACR"

# Build and push using ACR build (faster, happens in Azure)
print_warning "Building Docker image in Azure (this may take 5-10 minutes)..."
az acr build \
    --registry $ACR_NAME \
    --image kanad-api:latest \
    --file Dockerfile \
    .
print_success "Docker image built and pushed to ACR"

# ============================================================================
# Step 4: Create PostgreSQL Database
# ============================================================================

print_header "Step 4: Creating PostgreSQL Database"

# Generate random password if not set
if [ -z "$POSTGRES_PASSWORD" ]; then
    POSTGRES_PASSWORD=$(openssl rand -base64 32)
    echo -e "${YELLOW}Generated PostgreSQL password (save this!):${NC}"
    echo -e "${GREEN}$POSTGRES_PASSWORD${NC}"
    echo ""
fi

if az postgres flexible-server show --resource-group $RESOURCE_GROUP --name $POSTGRES_SERVER &> /dev/null; then
    print_warning "PostgreSQL server '$POSTGRES_SERVER' already exists"
else
    print_warning "Creating PostgreSQL server (this may take 5-10 minutes)..."
    az postgres flexible-server create \
        --resource-group $RESOURCE_GROUP \
        --name $POSTGRES_SERVER \
        --location $LOCATION \
        --admin-user $POSTGRES_ADMIN \
        --admin-password "$POSTGRES_PASSWORD" \
        --sku-name $POSTGRES_SKU \
        --tier GeneralPurpose \
        --version 15 \
        --storage-size 128 \
        --public-access 0.0.0.0-255.255.255.255
    print_success "PostgreSQL server created"

    # Create database
    az postgres flexible-server db create \
        --resource-group $RESOURCE_GROUP \
        --server-name $POSTGRES_SERVER \
        --database-name $POSTGRES_DB
    print_success "Database created: $POSTGRES_DB"

    # Configure firewall for Azure services
    az postgres flexible-server firewall-rule create \
        --resource-group $RESOURCE_GROUP \
        --name $POSTGRES_SERVER \
        --rule-name AllowAzureServices \
        --start-ip-address 0.0.0.0 \
        --end-ip-address 0.0.0.0
    print_success "Firewall configured"
fi

# Build connection string
DATABASE_URL="postgresql://${POSTGRES_ADMIN}:${POSTGRES_PASSWORD}@${POSTGRES_SERVER}.postgres.database.azure.com:5432/${POSTGRES_DB}?sslmode=require"

# ============================================================================
# Step 5: Create App Service Plan
# ============================================================================

print_header "Step 5: Creating App Service Plan"

if az appservice plan show --name $APP_SERVICE_PLAN --resource-group $RESOURCE_GROUP &> /dev/null; then
    print_warning "App Service Plan '$APP_SERVICE_PLAN' already exists"
else
    print_warning "Creating App Service Plan with SKU: $APP_SERVICE_SKU"
    az appservice plan create \
        --name $APP_SERVICE_PLAN \
        --resource-group $RESOURCE_GROUP \
        --is-linux \
        --sku $APP_SERVICE_SKU
    print_success "App Service Plan created: $APP_SERVICE_PLAN ($APP_SERVICE_SKU)"
fi

# ============================================================================
# Step 6: Create Web App
# ============================================================================

print_header "Step 6: Creating Web App"

# Get ACR credentials
ACR_USERNAME=$(az acr credential show --name $ACR_NAME --query username -o tsv)
ACR_PASSWORD=$(az acr credential show --name $ACR_NAME --query passwords[0].value -o tsv)
ACR_LOGIN_SERVER="${ACR_NAME}.azurecr.io"

if az webapp show --name $WEB_APP_NAME --resource-group $RESOURCE_GROUP &> /dev/null; then
    print_warning "Web App '$WEB_APP_NAME' already exists, updating..."
else
    az webapp create \
        --resource-group $RESOURCE_GROUP \
        --plan $APP_SERVICE_PLAN \
        --name $WEB_APP_NAME \
        --deployment-container-image-name ${ACR_LOGIN_SERVER}/kanad-api:latest
    print_success "Web App created: $WEB_APP_NAME"
fi

# Configure container settings
az webapp config container set \
    --name $WEB_APP_NAME \
    --resource-group $RESOURCE_GROUP \
    --docker-custom-image-name ${ACR_LOGIN_SERVER}/kanad-api:latest \
    --docker-registry-server-url https://${ACR_LOGIN_SERVER} \
    --docker-registry-server-user $ACR_USERNAME \
    --docker-registry-server-password $ACR_PASSWORD
print_success "Container configured"

# Enable continuous deployment
az webapp deployment container config \
    --name $WEB_APP_NAME \
    --resource-group $RESOURCE_GROUP \
    --enable-cd true
print_success "Continuous deployment enabled"

# ============================================================================
# Step 7: Configure Application Settings (Environment Variables)
# ============================================================================

print_header "Step 7: Configuring Environment Variables"

# Generate SECRET_KEY if not set
if [ -z "$SECRET_KEY" ]; then
    SECRET_KEY=$(openssl rand -base64 48)
    echo -e "${YELLOW}Generated SECRET_KEY (save this!):${NC}"
    echo -e "${GREEN}$SECRET_KEY${NC}"
    echo ""
fi

# Configure app settings
az webapp config appsettings set \
    --resource-group $RESOURCE_GROUP \
    --name $WEB_APP_NAME \
    --settings \
        DATABASE_URL="$DATABASE_URL" \
        SECRET_KEY="$SECRET_KEY" \
        ENVIRONMENT="production" \
        DEBUG="false" \
        CORS_ORIGINS="https://your-frontend.vercel.app,http://localhost:3000" \
        WEBSITES_PORT="8000" \
        WEBSITES_ENABLE_APP_SERVICE_STORAGE="false" \
        DOCKER_ENABLE_CI="true"

print_success "Environment variables configured"

print_warning "⚠️  Remember to set these additional variables in Azure Portal:"
echo "  - GOOGLE_CLIENT_ID"
echo "  - GOOGLE_CLIENT_SECRET"
echo "  - SMTP_HOST, SMTP_USER, SMTP_PASSWORD"
echo "  - ADMIN_PASSWORD"
echo "  - IBM_QUANTUM_TOKEN (optional)"
echo "  - BLUEQUBIT_API_KEY (optional)"

# ============================================================================
# Step 8: Configure CORS
# ============================================================================

print_header "Step 8: Configuring CORS"

az webapp cors add \
    --resource-group $RESOURCE_GROUP \
    --name $WEB_APP_NAME \
    --allowed-origins https://your-frontend.vercel.app http://localhost:3000

print_success "CORS configured"

# ============================================================================
# Step 9: Enable Logging
# ============================================================================

print_header "Step 9: Enabling Logging"

az webapp log config \
    --name $WEB_APP_NAME \
    --resource-group $RESOURCE_GROUP \
    --docker-container-logging filesystem

print_success "Logging enabled"

# ============================================================================
# Deployment Complete!
# ============================================================================

print_header "Deployment Complete! 🎉"

WEB_APP_URL="https://${WEB_APP_NAME}.azurewebsites.net"

echo -e "${GREEN}Your Kanad API is deployed!${NC}\n"
echo -e "📍 URL: ${BLUE}${WEB_APP_URL}${NC}"
echo -e "🔍 Health Check: ${BLUE}${WEB_APP_URL}/api/health${NC}"
echo -e "📊 Resource Group: ${BLUE}${RESOURCE_GROUP}${NC}"
echo -e "💻 App Service: ${BLUE}${WEB_APP_NAME}${NC}"
echo -e "🐘 PostgreSQL: ${BLUE}${POSTGRES_SERVER}.postgres.database.azure.com${NC}"
echo -e "🐳 Container Registry: ${BLUE}${ACR_LOGIN_SERVER}${NC}"
echo ""
echo -e "${YELLOW}Important Credentials (save these securely!):${NC}"
echo -e "  PostgreSQL Password: ${GREEN}${POSTGRES_PASSWORD}${NC}"
echo -e "  SECRET_KEY: ${GREEN}${SECRET_KEY}${NC}"
echo ""
echo -e "${YELLOW}Next Steps:${NC}"
echo "  1. Set additional environment variables in Azure Portal"
echo "  2. Update CORS_ORIGINS with your frontend URL"
echo "  3. Test the API: curl ${WEB_APP_URL}/api/health"
echo "  4. View logs: az webapp log tail --name ${WEB_APP_NAME} --resource-group ${RESOURCE_GROUP}"
echo "  5. Update your frontend .env with: NEXT_PUBLIC_API_URL=${WEB_APP_URL}"
echo ""
echo -e "${BLUE}Estimated Monthly Cost (HPC Configuration):${NC}"
echo "  App Service (E4_v3 - 4 vCores, 32GB RAM): ~\$300/month"
echo "  PostgreSQL (GP_Standard_D4s_v3 - 4 vCores, 16GB RAM): ~\$150/month"
echo "  Container Registry (Basic): ~\$5/month"
echo "  Total: ~\$455/month"
echo ""
echo -e "${GREEN}Credit Usage Projection:${NC}"
echo "  - With \$1000 credits: ~2.2 months of runtime"
echo "  - With \$5000 startup credits: ~11 months of runtime"
echo "  - Recommended: Monitor usage and scale down to E2_v3 during idle periods"
echo ""
echo -e "${YELLOW}Note: This HPC-level configuration supports:${NC}"
echo "  - Hamiltonian construction for SQD and VQE"
echo "  - Large molecular simulations (up to 50 atoms)"
echo "  - Multiple concurrent quantum chemistry jobs"
echo "  - 32GB RAM handles memory-intensive PySCF calculations"
echo ""
print_success "Deployment successful! 🚀"
