# Kanad Azure FX-Series VM Deployment Guide

## Overview

This guide walks you through deploying the Kanad quantum chemistry platform on an Azure FX16-4mds_v2 virtual machine with 336 GB RAM - perfect for intensive quantum chemistry workloads!

## Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│              Azure FX16-4mds_v2 VM (336 GB RAM)                 │
│                                                                 │
│  ┌─────────────────────────────────────────────────────────┐  │
│  │  Nginx (Port 80/443)                                     │  │
│  │  - Reverse Proxy                                         │  │
│  │  - SSL Termination (Let's Encrypt)                       │  │
│  │  - Load Balancing                                        │  │
│  └────────────────────┬────────────────────────────────────┘  │
│                       │                                         │
│  ┌────────────────────▼────────────────┐                       │
│  │  FastAPI/Uvicorn (Port 8000)        │                       │
│  │  - Quantum Chemistry API             │                       │
│  │  - VQE & SQD Solvers                 │                       │
│  │  - PySCF, Qiskit Integration         │                       │
│  │  Memory Allocation: ~200 GB          │                       │
│  └──────────────────────────────────────┘                       │
│                                                                 │
│  ┌─────────────────────────────────────┐                       │
│  │  PostgreSQL 15 (Port 5432)          │                       │
│  │  - Molecular Data Storage            │                       │
│  │  - Experiment Results                │                       │
│  │  - User Authentication               │                       │
│  │  Memory Allocation: ~100 GB          │                       │
│  │  Data Location: /mnt/data/postgresql │                       │
│  └─────────────────────────────────────┘                       │
│                                                                 │
│  ┌─────────────────────────────────────┐                       │
│  │  Data Disk (512 GB Premium SSD)     │                       │
│  │  - PostgreSQL Data                   │                       │
│  │  - Molecular Files (.xyz, .mol)      │                       │
│  │  - Computation Results               │                       │
│  │  Mount Point: /mnt/data              │                       │
│  └─────────────────────────────────────┘                       │
│                                                                 │
│  System Resources:                                              │
│  - vCPUs: 4 (Intel Xeon 5th Gen @ 4.0 GHz)                    │
│  - RAM: 336 GB (84 GB per vCPU!)                               │
│  - Local NVMe: 672 GiB (ultra-fast temp storage)               │
│  - OS Disk: 128 GB Premium SSD                                 │
└─────────────────────────────────────────────────────────────────┘
         │
         │ HTTPS (Port 443)
         ▼
    Public Internet
```

## Prerequisites

### 1. Azure Account
- Active Azure subscription
- $5,000 startup credits activated
- Azure CLI installed

### 2. Local Machine Requirements
- Linux, macOS, or WSL2
- SSH client
- Git

### 3. Install Azure CLI

**Ubuntu/Debian:**
```bash
curl -sL https://aka.ms/InstallAzureCLIDeb | sudo bash
```

**macOS:**
```bash
brew install azure-cli
```

**Verify installation:**
```bash
az --version
az login
```

## Deployment Steps

### Step 1: Clone Kanad Repository

```bash
cd ~
git clone <your-kanad-repo-url>
cd kanad
```

### Step 2: Review Configuration

Edit `deploy-azure-vm.sh` to customize:

```bash
# Basic settings
RESOURCE_GROUP="kanad-vm-rg"
LOCATION="eastus"  # Change to your preferred region
VM_NAME="kanad-fx16"
ADMIN_USERNAME="kanadmin"

# VM size - FX16-4mds_v2 recommended
VM_SIZE="Standard_FX16-4mds_v2"  # 4 vCPUs, 336 GB RAM

# Optional: Your domain name
DOMAIN_NAME="kanad.yourdomain.com"
ADMIN_EMAIL="your-email@example.com"
```

### Step 3: Run Deployment Script

```bash
chmod +x deploy-azure-vm.sh
./deploy-azure-vm.sh
```

**Expected Duration**: 20-30 minutes

The script will:
1. ✅ Create resource group
2. ✅ Create virtual network and security groups
3. ✅ Create public IP address
4. ✅ Create and configure FX16-4mds_v2 VM
5. ✅ Attach 512 GB data disk
6. ✅ Install and configure PostgreSQL 15
7. ✅ Install Python 3.11 and dependencies
8. ✅ Install and configure Nginx
9. ✅ Set up systemd service
10. ✅ Generate credentials

### Step 4: Save Credentials

The script will display and save credentials to `/tmp/kanad_vm_credentials.txt`:

```
VM Information:
- Public IP: X.X.X.X
- SSH: ssh kanadmin@X.X.X.X

PostgreSQL:
- Connection: postgresql://kanad_user:PASSWORD@localhost:5432/kanad_db

Application:
- API URL: http://X.X.X.X
- Secret Key: XXXXXXXXX
```

**IMPORTANT**: Copy these credentials to a secure location and delete the temp file:
```bash
cp /tmp/kanad_vm_credentials.txt ~/kanad_credentials_backup.txt
rm /tmp/kanad_vm_credentials.txt
```

### Step 5: Deploy Application Code

#### Option A: Copy from Local Machine

```bash
# From your local machine
cd /path/to/kanad
scp -r * kanadmin@<VM_PUBLIC_IP>:/opt/kanad/
```

#### Option B: Clone from Git

```bash
# SSH into VM
ssh kanadmin@<VM_PUBLIC_IP>

# Clone repository
cd /opt/kanad
git clone <your-repo-url> .
```

### Step 6: Install Python Dependencies

```bash
# SSH into VM
ssh kanadmin@<VM_PUBLIC_IP>

# Navigate to project
cd /opt/kanad

# Activate virtual environment
source env/bin/activate

# Install dependencies
pip install --upgrade pip wheel setuptools
pip install -r requirements.txt
pip install -r api/requirements.txt
```

### Step 7: Run Database Migrations

```bash
cd /opt/kanad
source env/bin/activate

# Run Alembic migrations
alembic upgrade head

# Or manually create tables if using SQLAlchemy
python -c "from api.core.database import Base, engine; Base.metadata.create_all(bind=engine)"
```

### Step 8: Start the Service

```bash
# Start Kanad service
sudo systemctl start kanad

# Check status
sudo systemctl status kanad

# View logs
sudo journalctl -u kanad -f
```

### Step 9: Test the API

```bash
# Health check
curl http://<VM_PUBLIC_IP>/api/health

# Should return:
# {"status": "healthy", "database": "connected"}
```

### Step 10: Set Up SSL (Optional but Recommended)

#### If you have a domain name:

```bash
# Point your domain DNS A record to VM_PUBLIC_IP
# Wait for DNS propagation (5-30 minutes)

# SSH into VM
ssh kanadmin@<VM_PUBLIC_IP>

# Run Certbot
sudo certbot --nginx -d kanad.yourdomain.com --email your-email@example.com --agree-tos --non-interactive

# Test auto-renewal
sudo certbot renew --dry-run
```

#### Test HTTPS:
```bash
curl https://kanad.yourdomain.com/api/health
```

## Post-Deployment Configuration

### 1. Configure Environment Variables

Edit `/opt/kanad/.env`:

```bash
ssh kanadmin@<VM_PUBLIC_IP>
sudo nano /opt/kanad/.env
```

Add/update:
```env
# Google OAuth (optional)
GOOGLE_CLIENT_ID=your-google-client-id
GOOGLE_CLIENT_SECRET=your-google-client-secret

# Email settings (optional)
SMTP_HOST=smtp.gmail.com
SMTP_PORT=587
SMTP_USERNAME=your-email@gmail.com
SMTP_PASSWORD=your-app-password

# Quantum backends (optional)
IBM_QUANTUM_TOKEN=your-ibm-token
BLUEQUBIT_API_KEY=your-bluequbit-key
```

Restart service:
```bash
sudo systemctl restart kanad
```

### 2. Configure Firewall (Optional - Extra Security)

```bash
# Only allow specific IPs to access SSH
az network nsg rule update \
  --resource-group kanad-vm-rg \
  --nsg-name kanad-nsg \
  --name AllowSSH \
  --source-address-prefixes YOUR_IP_ADDRESS/32
```

### 3. Set Up Monitoring

```bash
# Install monitoring tools
sudo apt-get install -y prometheus-node-exporter

# View resource usage
htop

# Monitor PostgreSQL
sudo -u postgres psql -c "SELECT * FROM pg_stat_activity;"
```

### 4. Configure Backups

Create backup script `/opt/kanad/backup.sh`:

```bash
#!/bin/bash
BACKUP_DIR="/mnt/data/backups"
DATE=$(date +%Y%m%d_%H%M%S)

# Backup PostgreSQL
sudo -u postgres pg_dump kanad_db | gzip > $BACKUP_DIR/db_$DATE.sql.gz

# Backup application data
tar -czf $BACKUP_DIR/app_$DATE.tar.gz /opt/kanad

# Keep only last 7 days
find $BACKUP_DIR -name "*.gz" -mtime +7 -delete

echo "Backup completed: $DATE"
```

Set up cron job:
```bash
crontab -e

# Add daily backup at 2 AM
0 2 * * * /opt/kanad/backup.sh >> /var/log/kanad_backup.log 2>&1
```

## Cost Management

### Monthly Costs (FX16-4mds_v2)

| Resource | Specification | Monthly Cost |
|----------|--------------|--------------|
| VM | FX16-4mds_v2 (4 vCPUs, 336 GB RAM) | $1,805 |
| OS Disk | 128 GB Premium SSD | $20 |
| Data Disk | 512 GB Premium SSD | $80 |
| Public IP | Static IP | $4 |
| Bandwidth | ~100 GB/month | $10 |
| **Total** | | **~$1,920/month** |

### Credit Usage

With **$5,000 startup credits**:
- **Runtime**: ~2.6 months
- **Perfect for**: 3-4 month MVP development and testing

### Cost Optimization Tips

#### 1. Scale Down During Idle Periods

```bash
# Stop VM when not in use (nights/weekends)
az vm deallocate --resource-group kanad-vm-rg --name kanad-fx16

# Start VM when needed
az vm start --resource-group kanad-vm-rg --name kanad-fx16
```

**Savings**: Only pay for storage (~$100/month) when stopped

#### 2. Use Auto-Shutdown

```bash
# Set up auto-shutdown at 10 PM
az vm auto-shutdown --resource-group kanad-vm-rg --name kanad-fx16 --time 2200
```

#### 3. Resize VM for Different Workloads

```bash
# Scale down to FX4mds_v2 for light workloads ($451/month)
az vm deallocate --resource-group kanad-vm-rg --name kanad-fx16
az vm resize --resource-group kanad-vm-rg --name kanad-fx16 --size Standard_FX4mds_v2
az vm start --resource-group kanad-vm-rg --name kanad-fx16

# Scale back up when needed
az vm deallocate --resource-group kanad-vm-rg --name kanad-fx16
az vm resize --resource-group kanad-vm-rg --name kanad-fx16 --size Standard_FX16-4mds_v2
az vm start --resource-group kanad-vm-rg --name kanad-fx16
```

#### 4. Use Spot Instances (Advanced)

For non-critical workloads:
```bash
# Create spot instance (80% discount!)
az vm create \
  --resource-group kanad-vm-rg \
  --name kanad-fx16-spot \
  --size Standard_FX16-4mds_v2 \
  --priority Spot \
  --max-price 0.5 \
  --eviction-policy Deallocate
```

## Monitoring & Maintenance

### View Application Logs

```bash
# Real-time logs
sudo journalctl -u kanad -f

# Last 100 lines
sudo journalctl -u kanad -n 100

# Today's logs
sudo journalctl -u kanad --since today
```

### View Nginx Logs

```bash
# Access logs
sudo tail -f /var/log/nginx/access.log

# Error logs
sudo tail -f /var/log/nginx/error.log
```

### View PostgreSQL Logs

```bash
sudo tail -f /var/log/postgresql/postgresql-15-main.log
```

### Monitor System Resources

```bash
# CPU, RAM, processes
htop

# Disk usage
df -h

# PostgreSQL connections
sudo -u postgres psql -c "SELECT count(*) FROM pg_stat_activity;"

# Check service status
sudo systemctl status kanad postgresql nginx
```

### Update Application

```bash
# SSH into VM
ssh kanadmin@<VM_PUBLIC_IP>

# Pull latest code
cd /opt/kanad
git pull

# Update dependencies if needed
source env/bin/activate
pip install -r requirements.txt --upgrade

# Restart service
sudo systemctl restart kanad
```

## Troubleshooting

### Service Won't Start

```bash
# Check logs
sudo journalctl -u kanad -n 50

# Check if port 8000 is in use
sudo netstat -tulpn | grep 8000

# Test manually
cd /opt/kanad
source env/bin/activate
uvicorn api.main:app --host 0.0.0.0 --port 8000
```

### Database Connection Issues

```bash
# Check PostgreSQL status
sudo systemctl status postgresql

# Test connection
psql postgresql://kanad_user:PASSWORD@localhost:5432/kanad_db

# Check PostgreSQL logs
sudo tail -f /var/log/postgresql/postgresql-15-main.log
```

### Nginx Issues

```bash
# Check nginx status
sudo systemctl status nginx

# Test configuration
sudo nginx -t

# Reload configuration
sudo systemctl reload nginx

# Check logs
sudo tail -f /var/log/nginx/error.log
```

### Out of Memory

```bash
# Check memory usage
free -h

# Check which process is using memory
ps aux --sort=-%mem | head -n 10

# Restart service to free memory
sudo systemctl restart kanad
```

### Disk Full

```bash
# Check disk usage
df -h

# Find large files
sudo du -sh /* | sort -hr | head -n 10

# Clean up PostgreSQL logs
sudo find /var/log/postgresql/ -name "*.log" -mtime +7 -delete

# Clean up old backups
sudo find /mnt/data/backups/ -mtime +30 -delete
```

## Performance Optimization

### PostgreSQL Tuning

The deployment script already configures PostgreSQL for 336 GB RAM, but you can further optimize:

```sql
-- Connect to PostgreSQL
sudo -u postgres psql

-- Check current settings
SHOW shared_buffers;
SHOW effective_cache_size;

-- Monitor query performance
SELECT * FROM pg_stat_statements ORDER BY total_time DESC LIMIT 10;
```

### Python Application Tuning

Edit `/etc/systemd/system/kanad.service`:

```ini
[Service]
# Increase workers based on CPU count
ExecStart=/opt/kanad/env/bin/uvicorn api.main:app --host 0.0.0.0 --port 8000 --workers 4

# Increase file descriptors
LimitNOFILE=65536

# Increase process limit
LimitNPROC=4096
```

Reload and restart:
```bash
sudo systemctl daemon-reload
sudo systemctl restart kanad
```

### Nginx Tuning

Edit `/etc/nginx/sites-available/kanad`:

```nginx
# Increase timeouts for long-running computations
proxy_read_timeout 3600s;  # 1 hour
proxy_connect_timeout 600s;
proxy_send_timeout 600s;

# Increase buffer sizes
proxy_buffer_size 16k;
proxy_buffers 32 16k;
proxy_busy_buffers_size 64k;
```

## Security Best Practices

1. **Keep System Updated**
   ```bash
   sudo apt-get update && sudo apt-get upgrade -y
   ```

2. **Enable Firewall**
   ```bash
   sudo ufw enable
   sudo ufw allow 22/tcp
   sudo ufw allow 80/tcp
   sudo ufw allow 443/tcp
   ```

3. **Disable Root Login**
   ```bash
   sudo sed -i 's/PermitRootLogin yes/PermitRootLogin no/' /etc/ssh/sshd_config
   sudo systemctl restart sshd
   ```

4. **Set Up Fail2Ban**
   ```bash
   sudo apt-get install -y fail2ban
   sudo systemctl enable fail2ban
   sudo systemctl start fail2ban
   ```

5. **Regular Backups**
   - Set up automated daily backups
   - Test restore procedures monthly
   - Store backups in Azure Blob Storage

## Scaling Strategy

### When to Scale Up

Scale to **FX32-8mds_v2** (8 vCPUs, 672 GB RAM) if:
- Molecules consistently exceed 40 atoms
- Using advanced basis sets (cc-pvtz, cc-pvqz)
- Running 10+ concurrent experiments
- RAM usage consistently > 80%

### When to Scale Down

Scale to **FX4mds_v2** (4 vCPUs, 84 GB RAM) if:
- Usage drops after MVP phase
- Working with smaller molecules (< 20 atoms)
- Single user or light testing
- RAM usage consistently < 25%

## Support

For issues or questions:
1. Check the logs first
2. Review this documentation
3. Contact your development team
4. Open an issue on GitHub

## Summary

You now have a **production-ready quantum chemistry platform** running on:
- ✅ 336 GB RAM for large molecular calculations
- ✅ High-speed NVMe storage for fast I/O
- ✅ PostgreSQL database on dedicated disk
- ✅ Nginx reverse proxy with SSL support
- ✅ Systemd service for automatic restarts
- ✅ Comprehensive monitoring and logging

**Estimated cost**: $1,920/month (~2.6 months on $5,000 credits)

Perfect for your 3-4 month MVP deployment! 🚀
