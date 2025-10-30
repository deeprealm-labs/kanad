#!/bin/bash
# ============================================================================
# Kanad Azure FX-Series VM Deployment Script
# ============================================================================
# Deploys Kanad quantum chemistry platform on Azure FX-series VM
# Architecture:
#   - Single VM hosting PostgreSQL + FastAPI + Nginx
#   - FX16-4mds_v2: 4 vCPUs, 336 GB RAM, 672 GiB local NVMe
#   - PostgreSQL (port 5432) on localhost
#   - FastAPI/Uvicorn (port 8000) on localhost
#   - Nginx (port 80/443) as reverse proxy
# ============================================================================

set -e  # Exit on error

# ============================================================================
# Configuration
# ============================================================================

# Basic settings
RESOURCE_GROUP="kanad-vm-rg"
LOCATION="eastus"  # Change to your preferred region
VM_NAME="kanad-fx16"
ADMIN_USERNAME="kanadmin"

# VM Configuration - FX16-4mds_v2 (4 vCPUs, 336 GB RAM)
# Options:
# - Standard_FX4mds_v2: 4 vCPUs, 84 GB RAM - $451/month (basic)
# - Standard_FX8-4mds_v2: 4 vCPUs, 168 GB RAM - $902/month (good)
# - Standard_FX16-4mds_v2: 4 vCPUs, 336 GB RAM - $1,805/month (RECOMMENDED)
# - Standard_FX32-8mds_v2: 8 vCPUs, 672 GB RAM - $3,609/month (maximum)
VM_SIZE="Standard_FX16-4mds_v2"  # 4 vCPUs, 336 GB RAM

# Network settings
NSG_NAME="kanad-nsg"
VNET_NAME="kanad-vnet"
SUBNET_NAME="kanad-subnet"
PUBLIC_IP_NAME="kanad-public-ip"
NIC_NAME="kanad-nic"

# Disk configuration
OS_DISK_SIZE_GB=128  # OS disk
DATA_DISK_SIZE_GB=512  # Data disk for PostgreSQL and molecular data

# PostgreSQL configuration
POSTGRES_VERSION=15
POSTGRES_PASSWORD=""  # Will be generated if empty

# Application configuration
DOMAIN_NAME=""  # Optional: your domain (e.g., kanad.yourdomain.com)
ADMIN_EMAIL=""  # For Let's Encrypt SSL certificates

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

# Check if logged into Azure
if ! az account show &> /dev/null; then
    print_error "Not logged into Azure. Please run: az login"
    exit 1
fi
print_success "Logged into Azure"

# Show current subscription
SUBSCRIPTION=$(az account show --query name -o tsv)
SUBSCRIPTION_ID=$(az account show --query id -o tsv)
print_success "Using subscription: $SUBSCRIPTION"

# Check if SSH key exists
if [ ! -f ~/.ssh/id_rsa.pub ]; then
    print_warning "SSH key not found. Generating new SSH key..."
    ssh-keygen -t rsa -b 4096 -f ~/.ssh/id_rsa -N "" -C "kanad-azure-vm"
    print_success "SSH key generated"
fi
print_success "SSH key found"

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
# Step 2: Create Virtual Network
# ============================================================================

print_header "Step 2: Creating Virtual Network"

# Create VNet
if az network vnet show --name $VNET_NAME --resource-group $RESOURCE_GROUP &> /dev/null; then
    print_warning "VNet '$VNET_NAME' already exists"
else
    az network vnet create \
        --resource-group $RESOURCE_GROUP \
        --name $VNET_NAME \
        --address-prefix 10.0.0.0/16 \
        --subnet-name $SUBNET_NAME \
        --subnet-prefix 10.0.1.0/24
    print_success "VNet created: $VNET_NAME"
fi

# ============================================================================
# Step 3: Create Network Security Group
# ============================================================================

print_header "Step 3: Creating Network Security Group"

if az network nsg show --name $NSG_NAME --resource-group $RESOURCE_GROUP &> /dev/null; then
    print_warning "NSG '$NSG_NAME' already exists"
else
    az network nsg create \
        --resource-group $RESOURCE_GROUP \
        --name $NSG_NAME
    print_success "NSG created: $NSG_NAME"
fi

# Create NSG rules
print_warning "Configuring NSG rules..."

# Allow SSH (port 22)
az network nsg rule create \
    --resource-group $RESOURCE_GROUP \
    --nsg-name $NSG_NAME \
    --name AllowSSH \
    --priority 1000 \
    --source-address-prefixes '*' \
    --source-port-ranges '*' \
    --destination-address-prefixes '*' \
    --destination-port-ranges 22 \
    --access Allow \
    --protocol Tcp \
    --description "Allow SSH" \
    2>/dev/null || print_warning "SSH rule already exists"

# Allow HTTP (port 80)
az network nsg rule create \
    --resource-group $RESOURCE_GROUP \
    --nsg-name $NSG_NAME \
    --name AllowHTTP \
    --priority 1010 \
    --source-address-prefixes '*' \
    --source-port-ranges '*' \
    --destination-address-prefixes '*' \
    --destination-port-ranges 80 \
    --access Allow \
    --protocol Tcp \
    --description "Allow HTTP" \
    2>/dev/null || print_warning "HTTP rule already exists"

# Allow HTTPS (port 443)
az network nsg rule create \
    --resource-group $RESOURCE_GROUP \
    --nsg-name $NSG_NAME \
    --name AllowHTTPS \
    --priority 1020 \
    --source-address-prefixes '*' \
    --source-port-ranges '*' \
    --destination-address-prefixes '*' \
    --destination-port-ranges 443 \
    --access Allow \
    --protocol Tcp \
    --description "Allow HTTPS" \
    2>/dev/null || print_warning "HTTPS rule already exists"

print_success "NSG rules configured"

# ============================================================================
# Step 4: Create Public IP
# ============================================================================

print_header "Step 4: Creating Public IP"

if az network public-ip show --name $PUBLIC_IP_NAME --resource-group $RESOURCE_GROUP &> /dev/null; then
    print_warning "Public IP '$PUBLIC_IP_NAME' already exists"
    PUBLIC_IP=$(az network public-ip show --name $PUBLIC_IP_NAME --resource-group $RESOURCE_GROUP --query ipAddress -o tsv)
else
    az network public-ip create \
        --resource-group $RESOURCE_GROUP \
        --name $PUBLIC_IP_NAME \
        --sku Standard \
        --allocation-method Static \
        --dns-name $VM_NAME

    PUBLIC_IP=$(az network public-ip show --name $PUBLIC_IP_NAME --resource-group $RESOURCE_GROUP --query ipAddress -o tsv)
    print_success "Public IP created: $PUBLIC_IP"
fi

# ============================================================================
# Step 5: Create Network Interface
# ============================================================================

print_header "Step 5: Creating Network Interface"

if az network nic show --name $NIC_NAME --resource-group $RESOURCE_GROUP &> /dev/null; then
    print_warning "NIC '$NIC_NAME' already exists"
else
    az network nic create \
        --resource-group $RESOURCE_GROUP \
        --name $NIC_NAME \
        --vnet-name $VNET_NAME \
        --subnet $SUBNET_NAME \
        --public-ip-address $PUBLIC_IP_NAME \
        --network-security-group $NSG_NAME
    print_success "NIC created: $NIC_NAME"
fi

# ============================================================================
# Step 6: Create Virtual Machine
# ============================================================================

print_header "Step 6: Creating FX16-4mds_v2 Virtual Machine"
echo "This may take 5-10 minutes..."

if az vm show --name $VM_NAME --resource-group $RESOURCE_GROUP &> /dev/null; then
    print_warning "VM '$VM_NAME' already exists"
else
    print_warning "Creating VM with:"
    echo "  - Size: $VM_SIZE (4 vCPUs, 336 GB RAM)"
    echo "  - OS: Ubuntu 22.04 LTS"
    echo "  - OS Disk: $OS_DISK_SIZE_GB GB"
    echo "  - Admin: $ADMIN_USERNAME"

    az vm create \
        --resource-group $RESOURCE_GROUP \
        --name $VM_NAME \
        --size $VM_SIZE \
        --image Ubuntu2204 \
        --admin-username $ADMIN_USERNAME \
        --ssh-key-values ~/.ssh/id_rsa.pub \
        --nics $NIC_NAME \
        --os-disk-size-gb $OS_DISK_SIZE_GB \
        --storage-sku Premium_LRS \
        --priority Regular

    print_success "VM created: $VM_NAME"
fi

# ============================================================================
# Step 7: Create and Attach Data Disk
# ============================================================================

print_header "Step 7: Creating Data Disk for PostgreSQL"

DISK_NAME="${VM_NAME}-data-disk"

if az disk show --name $DISK_NAME --resource-group $RESOURCE_GROUP &> /dev/null; then
    print_warning "Data disk '$DISK_NAME' already exists"
else
    az disk create \
        --resource-group $RESOURCE_GROUP \
        --name $DISK_NAME \
        --size-gb $DATA_DISK_SIZE_GB \
        --sku Premium_LRS

    print_success "Data disk created: $DATA_DISK_SIZE_GB GB"
fi

# Attach disk to VM
if az vm disk list --resource-group $RESOURCE_GROUP --vm-name $VM_NAME --query "[?name=='$DISK_NAME']" -o tsv | grep -q "$DISK_NAME"; then
    print_warning "Data disk already attached"
else
    az vm disk attach \
        --resource-group $RESOURCE_GROUP \
        --vm-name $VM_NAME \
        --name $DISK_NAME

    print_success "Data disk attached to VM"
fi

# ============================================================================
# Step 8: Generate Configuration Files
# ============================================================================

print_header "Step 8: Generating Configuration Files"

# Generate PostgreSQL password if not set
if [ -z "$POSTGRES_PASSWORD" ]; then
    POSTGRES_PASSWORD=$(openssl rand -base64 32)
    print_success "Generated PostgreSQL password"
fi

# Generate secret key for FastAPI
SECRET_KEY=$(openssl rand -base64 32)

# Save credentials to file
CREDENTIALS_FILE="/tmp/kanad_vm_credentials.txt"
cat > $CREDENTIALS_FILE << EOF
KANAD AZURE VM DEPLOYMENT CREDENTIALS
=====================================
Generated: $(date)

VM Information:
- Resource Group: $RESOURCE_GROUP
- VM Name: $VM_NAME
- VM Size: $VM_SIZE (4 vCPUs, 336 GB RAM)
- Location: $LOCATION
- Public IP: $PUBLIC_IP
- Admin Username: $ADMIN_USERNAME

SSH Access:
ssh $ADMIN_USERNAME@$PUBLIC_IP

PostgreSQL:
- Version: $POSTGRES_VERSION
- Host: localhost
- Port: 5432
- Database: kanad_db
- Username: kanad_user
- Password: $POSTGRES_PASSWORD

Connection String:
postgresql://kanad_user:$POSTGRES_PASSWORD@localhost:5432/kanad_db

Application:
- API URL: http://$PUBLIC_IP (or https://$DOMAIN_NAME)
- Secret Key: $SECRET_KEY

IMPORTANT: Save this file securely and delete it after copying the credentials!
EOF

print_success "Credentials saved to: $CREDENTIALS_FILE"
cat $CREDENTIALS_FILE

# Create setup script for VM
SETUP_SCRIPT="/tmp/kanad_vm_setup.sh"
cat > $SETUP_SCRIPT << 'SETUP_EOF'
#!/bin/bash
# Kanad VM Setup Script - Runs on the VM

set -e

echo "========================================="
echo "  Kanad VM Setup Starting..."
echo "========================================="

# Update system
echo "Updating system packages..."
sudo apt-get update
sudo DEBIAN_FRONTEND=noninteractive apt-get upgrade -y

# Install essential packages
echo "Installing essential packages..."
sudo DEBIAN_FRONTEND=noninteractive apt-get install -y \
    build-essential \
    git \
    curl \
    wget \
    vim \
    htop \
    tmux \
    jq \
    unzip

# Install Python 3.11
echo "Installing Python 3.11..."
sudo add-apt-repository -y ppa:deadsnakes/ppa
sudo apt-get update
sudo DEBIAN_FRONTEND=noninteractive apt-get install -y \
    python3.11 \
    python3.11-venv \
    python3.11-dev \
    python3-pip

# Set Python 3.11 as default
sudo update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 1
sudo update-alternatives --install /usr/bin/python python /usr/bin/python3.11 1

# Install scientific computing libraries
echo "Installing scientific computing libraries..."
sudo DEBIAN_FRONTEND=noninteractive apt-get install -y \
    gfortran \
    libopenblas-dev \
    liblapack-dev \
    libhdf5-dev \
    libatlas-base-dev

# Format and mount data disk
echo "Setting up data disk..."
if ! mount | grep -q "/mnt/data"; then
    sudo parted /dev/sdc --script mklabel gpt mkpart primary ext4 0% 100%
    sudo mkfs.ext4 /dev/sdc1
    sudo mkdir -p /mnt/data
    sudo mount /dev/sdc1 /mnt/data
    echo '/dev/sdc1 /mnt/data ext4 defaults,nofail 0 2' | sudo tee -a /etc/fstab
    echo "✓ Data disk mounted at /mnt/data"
else
    echo "✓ Data disk already mounted"
fi

# Create directories
sudo mkdir -p /mnt/data/postgresql
sudo mkdir -p /mnt/data/kanad
sudo mkdir -p /opt/kanad
sudo chown -R $USER:$USER /opt/kanad

# Install PostgreSQL 15
echo "Installing PostgreSQL 15..."
sudo sh -c 'echo "deb http://apt.postgresql.org/pub/repos/apt $(lsb_release -cs)-pgdg main" > /etc/apt/sources.list.d/pgdg.list'
wget --quiet -O - https://www.postgresql.org/media/keys/ACCC4CF8.asc | sudo apt-key add -
sudo apt-get update
sudo DEBIAN_FRONTEND=noninteractive apt-get install -y postgresql-15 postgresql-contrib-15

# Stop PostgreSQL to move data directory
sudo systemctl stop postgresql

# Move PostgreSQL data to data disk
if [ ! -d "/mnt/data/postgresql/15/main" ]; then
    echo "Moving PostgreSQL data to data disk..."
    sudo mkdir -p /mnt/data/postgresql/15
    sudo cp -R /var/lib/postgresql/15/main /mnt/data/postgresql/15/
    sudo chown -R postgres:postgres /mnt/data/postgresql

    # Update PostgreSQL config
    sudo sed -i "s|data_directory = '/var/lib/postgresql/15/main'|data_directory = '/mnt/data/postgresql/15/main'|" /etc/postgresql/15/main/postgresql.conf
fi

# Configure PostgreSQL for high performance
echo "Configuring PostgreSQL..."
sudo tee -a /etc/postgresql/15/main/postgresql.conf > /dev/null << EOF

# Kanad Performance Tuning
# Total RAM: 336 GB, Allocating 100 GB to PostgreSQL
shared_buffers = 25GB
effective_cache_size = 75GB
maintenance_work_mem = 2GB
checkpoint_completion_target = 0.9
wal_buffers = 16MB
default_statistics_target = 100
random_page_cost = 1.1
effective_io_concurrency = 200
work_mem = 104857kB
min_wal_size = 2GB
max_wal_size = 8GB
max_worker_processes = 4
max_parallel_workers_per_gather = 2
max_parallel_workers = 4
max_parallel_maintenance_workers = 2
EOF

# Allow local connections
sudo sed -i "s/#listen_addresses = 'localhost'/listen_addresses = 'localhost'/g" /etc/postgresql/15/main/postgresql.conf

# Start PostgreSQL
sudo systemctl start postgresql
sudo systemctl enable postgresql

echo "✓ PostgreSQL installed and configured"

# Install Nginx
echo "Installing Nginx..."
sudo DEBIAN_FRONTEND=noninteractive apt-get install -y nginx
sudo systemctl enable nginx
echo "✓ Nginx installed"

# Install Certbot for SSL
echo "Installing Certbot..."
sudo DEBIAN_FRONTEND=noninteractive apt-get install -y certbot python3-certbot-nginx
echo "✓ Certbot installed"

# Install Docker (optional, for future use)
echo "Installing Docker..."
curl -fsSL https://get.docker.com -o get-docker.sh
sudo sh get-docker.sh
sudo usermod -aG docker $USER
rm get-docker.sh
echo "✓ Docker installed"

echo "========================================="
echo "  System setup complete!"
echo "========================================="
echo ""
echo "Next steps:"
echo "1. Create PostgreSQL database and user"
echo "2. Clone Kanad repository"
echo "3. Install Python dependencies"
echo "4. Configure Nginx"
echo "5. Set up systemd service"
echo ""
SETUP_EOF

chmod +x $SETUP_SCRIPT
print_success "Setup script created: $SETUP_SCRIPT"

# ============================================================================
# Step 9: Copy Files and Run Setup on VM
# ============================================================================

print_header "Step 9: Running Setup on VM"
echo "This may take 10-15 minutes..."

# Wait for VM to be ready
print_warning "Waiting for VM to be fully ready..."
sleep 30

# Copy setup script to VM
print_warning "Copying setup script to VM..."
scp -o StrictHostKeyChecking=no $SETUP_SCRIPT $ADMIN_USERNAME@$PUBLIC_IP:/tmp/

# Run setup script on VM
print_warning "Running setup script on VM..."
ssh -o StrictHostKeyChecking=no $ADMIN_USERNAME@$PUBLIC_IP "bash /tmp/kanad_vm_setup.sh"

print_success "VM setup complete!"

# ============================================================================
# Step 10: Create PostgreSQL Database
# ============================================================================

print_header "Step 10: Setting up PostgreSQL Database"

# Create database setup script
DB_SETUP_SCRIPT="/tmp/kanad_db_setup.sql"
cat > $DB_SETUP_SCRIPT << EOF
-- Create database
CREATE DATABASE kanad_db;

-- Create user
CREATE USER kanad_user WITH ENCRYPTED PASSWORD '$POSTGRES_PASSWORD';

-- Grant privileges
GRANT ALL PRIVILEGES ON DATABASE kanad_db TO kanad_user;

-- Connect to database
\\c kanad_db

-- Grant schema privileges
GRANT ALL ON SCHEMA public TO kanad_user;
EOF

# Copy and run DB setup
scp -o StrictHostKeyChecking=no $DB_SETUP_SCRIPT $ADMIN_USERNAME@$PUBLIC_IP:/tmp/
ssh -o StrictHostKeyChecking=no $ADMIN_USERNAME@$PUBLIC_IP "sudo -u postgres psql -f /tmp/kanad_db_setup.sql"

print_success "PostgreSQL database created"

# ============================================================================
# Step 11: Deploy Kanad Application
# ============================================================================

print_header "Step 11: Deploying Kanad Application"

# Create deployment script
DEPLOY_SCRIPT="/tmp/kanad_deploy.sh"
cat > $DEPLOY_SCRIPT << DEPLOY_EOF
#!/bin/bash
set -e

echo "Deploying Kanad application..."

cd /opt/kanad

# Clone repository if not exists
if [ ! -d ".git" ]; then
    echo "Cloning Kanad repository..."
    # Replace with your actual git repository
    # git clone https://github.com/yourusername/kanad.git .
    echo "NOTE: Please clone your repository manually to /opt/kanad"
    echo "For now, you'll need to copy your code to the VM"
fi

# Create virtual environment
if [ ! -d "env" ]; then
    echo "Creating virtual environment..."
    python3.11 -m venv env
fi

# Activate virtual environment and install dependencies
source env/bin/activate

echo "Installing Python dependencies..."
pip install --upgrade pip wheel setuptools

# Install requirements (assuming requirements.txt exists)
if [ -f "requirements.txt" ]; then
    pip install -r requirements.txt
fi

if [ -f "api/requirements.txt" ]; then
    pip install -r api/requirements.txt
fi

# Create .env file
echo "Creating .env file..."
cat > /opt/kanad/.env << EOF
# Database
DATABASE_URL=postgresql://kanad_user:$POSTGRES_PASSWORD@localhost:5432/kanad_db

# Security
SECRET_KEY=$SECRET_KEY
ALGORITHM=HS256
ACCESS_TOKEN_EXPIRE_MINUTES=30

# CORS
CORS_ORIGINS=["http://localhost:3000","http://$PUBLIC_IP","https://$DOMAIN_NAME"]

# Environment
ENVIRONMENT=production
DEBUG=False

# Compute Limits
MAX_ATOMS_PER_MOLECULE=50
MAX_VQE_ITERATIONS=1000
MAX_SQD_ITERATIONS=200
MAX_EXPERIMENTS_PER_HOUR=20
MAX_EXPERIMENTS_PER_DAY=100
EOF

echo "✓ Kanad application deployed"
DEPLOY_EOF

chmod +x $DEPLOY_SCRIPT
scp -o StrictHostKeyChecking=no $DEPLOY_SCRIPT $ADMIN_USERNAME@$PUBLIC_IP:/tmp/
ssh -o StrictHostKeyChecking=no $ADMIN_USERNAME@$PUBLIC_IP "bash /tmp/kanad_deploy.sh"

print_success "Application deployment script executed"

# ============================================================================
# Step 12: Configure Nginx
# ============================================================================

print_header "Step 12: Configuring Nginx"

NGINX_CONFIG="/tmp/kanad_nginx.conf"
cat > $NGINX_CONFIG << 'NGINX_EOF'
# Kanad Nginx Configuration

upstream kanad_api {
    server localhost:8000;
}

# HTTP server - redirect to HTTPS
server {
    listen 80;
    server_name _;

    # Allow Let's Encrypt challenges
    location /.well-known/acme-challenge/ {
        root /var/www/html;
    }

    # Redirect all other traffic to HTTPS
    location / {
        return 301 https://$host$request_uri;
    }
}

# HTTPS server
server {
    listen 443 ssl http2;
    server_name _;

    # SSL certificate paths (will be configured by Certbot)
    ssl_certificate /etc/ssl/certs/ssl-cert-snakeoil.pem;
    ssl_certificate_key /etc/ssl/private/ssl-cert-snakeoil.key;

    # SSL configuration
    ssl_protocols TLSv1.2 TLSv1.3;
    ssl_ciphers HIGH:!aNULL:!MD5;
    ssl_prefer_server_ciphers on;

    # Security headers
    add_header Strict-Transport-Security "max-age=31536000; includeSubDomains" always;
    add_header X-Frame-Options "DENY" always;
    add_header X-Content-Type-Options "nosniff" always;
    add_header X-XSS-Protection "1; mode=block" always;

    # Increase timeouts for long-running quantum computations
    proxy_read_timeout 3600s;
    proxy_connect_timeout 600s;
    proxy_send_timeout 600s;

    # API proxy
    location /api {
        proxy_pass http://kanad_api;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection 'upgrade';
        proxy_set_header Host $host;
        proxy_cache_bypass $http_upgrade;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;

        # CORS headers
        add_header 'Access-Control-Allow-Origin' '*' always;
        add_header 'Access-Control-Allow-Methods' 'GET, POST, PUT, DELETE, OPTIONS' always;
        add_header 'Access-Control-Allow-Headers' 'DNT,User-Agent,X-Requested-With,If-Modified-Since,Cache-Control,Content-Type,Range,Authorization' always;

        if ($request_method = 'OPTIONS') {
            return 204;
        }
    }

    # WebSocket support
    location /ws {
        proxy_pass http://kanad_api;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
    }

    # Health check endpoint
    location /health {
        proxy_pass http://kanad_api/api/health;
        proxy_http_version 1.1;
        proxy_set_header Host $host;
    }

    # Root
    location / {
        return 200 "Kanad Quantum Chemistry Platform API\n";
        add_header Content-Type text/plain;
    }
}
NGINX_EOF

scp -o StrictHostKeyChecking=no $NGINX_CONFIG $ADMIN_USERNAME@$PUBLIC_IP:/tmp/
ssh -o StrictHostKeyChecking=no $ADMIN_USERNAME@$PUBLIC_IP "sudo mv /tmp/kanad_nginx.conf /etc/nginx/sites-available/kanad && sudo ln -sf /etc/nginx/sites-available/kanad /etc/nginx/sites-enabled/ && sudo rm -f /etc/nginx/sites-enabled/default && sudo nginx -t && sudo systemctl reload nginx"

print_success "Nginx configured"

# ============================================================================
# Step 13: Create Systemd Service
# ============================================================================

print_header "Step 13: Creating Systemd Service"

SYSTEMD_SERVICE="/tmp/kanad.service"
cat > $SYSTEMD_SERVICE << SERVICE_EOF
[Unit]
Description=Kanad Quantum Chemistry Platform API
After=network.target postgresql.service

[Service]
Type=simple
User=$ADMIN_USERNAME
WorkingDirectory=/opt/kanad
Environment="PATH=/opt/kanad/env/bin"
ExecStart=/opt/kanad/env/bin/uvicorn api.main:app --host 0.0.0.0 --port 8000 --workers 4
Restart=always
RestartSec=10

# Performance optimizations
LimitNOFILE=65536
LimitNPROC=4096

[Install]
WantedBy=multi-user.target
SERVICE_EOF

scp -o StrictHostKeyChecking=no $SYSTEMD_SERVICE $ADMIN_USERNAME@$PUBLIC_IP:/tmp/
ssh -o StrictHostKeyChecking=no $ADMIN_USERNAME@$PUBLIC_IP "sudo mv /tmp/kanad.service /etc/systemd/system/ && sudo systemctl daemon-reload && sudo systemctl enable kanad.service"

print_success "Systemd service created"

# ============================================================================
# Deployment Complete
# ============================================================================

print_header "Deployment Complete! 🚀"

echo ""
echo -e "${GREEN}Your Kanad VM is ready!${NC}\n"
echo -e "📍 VM Name: ${BLUE}$VM_NAME${NC}"
echo -e "💻 VM Size: ${BLUE}$VM_SIZE${NC} (4 vCPUs, 336 GB RAM)"
echo -e "🌍 Public IP: ${BLUE}$PUBLIC_IP${NC}"
echo -e "🔑 SSH Access: ${BLUE}ssh $ADMIN_USERNAME@$PUBLIC_IP${NC}"
echo ""
echo -e "${BLUE}Estimated Monthly Cost:${NC}"
echo "  VM (FX16-4mds_v2): ~\$1,805/month"
echo "  Data Disk (512 GB): ~\$80/month"
echo "  Public IP: ~\$4/month"
echo "  Total: ~\$1,889/month"
echo ""
echo -e "${GREEN}Credit Usage Projection:${NC}"
echo "  - With \$5000 credits: ~2.6 months of runtime"
echo "  - Perfect for 3-4 month MVP deployment"
echo ""
echo -e "${YELLOW}Next Steps:${NC}"
echo "  1. SSH into VM: ssh $ADMIN_USERNAME@$PUBLIC_IP"
echo "  2. Copy your code: scp -r /path/to/kanad/* $ADMIN_USERNAME@$PUBLIC_IP:/opt/kanad/"
echo "  3. Or clone your repo: cd /opt/kanad && git clone <your-repo> ."
echo "  4. Install dependencies: cd /opt/kanad && source env/bin/activate && pip install -r requirements.txt"
echo "  5. Run database migrations: cd /opt/kanad && source env/bin/activate && alembic upgrade head"
echo "  6. Start the service: sudo systemctl start kanad"
echo "  7. Check status: sudo systemctl status kanad"
echo "  8. Test API: curl http://$PUBLIC_IP/api/health"
echo ""
echo -e "${YELLOW}Optional - Set up SSL:${NC}"
if [ -n "$DOMAIN_NAME" ] && [ -n "$ADMIN_EMAIL" ]; then
    echo "  sudo certbot --nginx -d $DOMAIN_NAME --email $ADMIN_EMAIL --agree-tos --non-interactive"
else
    echo "  sudo certbot --nginx -d your-domain.com --email your-email@example.com"
fi
echo ""
echo -e "${YELLOW}Monitoring:${NC}"
echo "  - View logs: sudo journalctl -u kanad -f"
echo "  - Monitor resources: htop"
echo "  - PostgreSQL logs: sudo tail -f /var/log/postgresql/postgresql-15-main.log"
echo ""
echo -e "${GREEN}Credentials saved to: $CREDENTIALS_FILE${NC}"
echo -e "${RED}IMPORTANT: Save credentials file securely and delete it after copying!${NC}"
echo ""
print_success "Deployment successful! 🎉"
