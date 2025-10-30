# Azure App Service Tier Selection for Quantum Chemistry Workloads

## Overview

This guide helps you choose the right Azure App Service tier for Kanad based on your computational needs and budget.

## Recommended Tiers for Quantum Chemistry (Memory-Optimized E-series)

### E2_v3 - Entry Level HPC
- **Specs**: 2 vCores, 16GB RAM
- **Cost**: ~$150/month
- **Use case**:
  - Small molecules (< 10 atoms)
  - Testing and development
  - Single user workloads
  - Light VQE runs with small basis sets

### E4_v3 - Recommended for MVP ✅
- **Specs**: 4 vCores, 32GB RAM
- **Cost**: ~$300/month
- **Use case**:
  - Medium molecules (10-30 atoms)
  - Multiple concurrent users (2-3)
  - VQE and SQD calculations
  - Hamiltonian construction for medium-sized systems
  - **Best balance of cost and performance for MVP**

### E8_v3 - High Performance
- **Specs**: 8 vCores, 64GB RAM
- **Cost**: ~$600/month
- **Use case**:
  - Large molecules (30-50 atoms)
  - Multiple concurrent jobs (5+)
  - Heavy Hamiltonian construction
  - Production workloads with high traffic
  - Advanced basis sets (6-31g*, cc-pvdz)

### E16_v3 - Maximum Performance
- **Specs**: 16 vCores, 128GB RAM
- **Cost**: ~$1200/month
- **Use case**:
  - Very large molecules (50+ atoms)
  - Maximum concurrent users
  - Enterprise-level workloads
  - Research-grade calculations
  - **Only use if absolutely necessary**

## Alternative: Compute-Optimized D-series

### D4s_v3
- **Specs**: 4 vCores, 16GB RAM
- **Cost**: ~$200/month
- **Use case**: Good CPU/memory balance, cheaper than E-series
- **Trade-off**: Less RAM than E4_v3, may struggle with large molecules

### D8s_v3
- **Specs**: 8 vCores, 32GB RAM
- **Cost**: ~$400/month
- **Use case**: More CPU power, same RAM as E4_v3
- **Trade-off**: Good for CPU-bound tasks, but E8_v3 is better for memory-bound quantum chemistry

## Cost Optimization Strategies

### 1. Dynamic Scaling
```bash
# Scale up during heavy workloads
az appservice plan update --resource-group kanad-rg --name kanad-plan --sku E8_v3

# Scale down during idle periods (nights/weekends)
az appservice plan update --resource-group kanad-rg --name kanad-plan --sku E2_v3
```

### 2. Credit Projections

| Tier | Monthly Cost | $1000 Credits | $5000 Credits |
|------|-------------|---------------|---------------|
| E2_v3 | $150 | 6.7 months | 33 months |
| E4_v3 | $300 | 3.3 months | 16.7 months |
| E8_v3 | $600 | 1.7 months | 8.3 months |
| E16_v3 | $1200 | 0.8 months | 4.2 months |

**Note**: Costs include PostgreSQL (~$150/month) and ACR (~$5/month)

### 3. Auto-Scaling Rules

Enable auto-scaling to automatically adjust based on CPU/Memory usage:

```bash
# Auto-scale between 2 and 5 instances based on CPU
az monitor autoscale create \
  --resource-group kanad-rg \
  --resource kanad-plan \
  --resource-type Microsoft.Web/serverfarms \
  --name autoscale-kanad \
  --min-count 1 \
  --max-count 3 \
  --count 1

# Scale out when CPU > 70%
az monitor autoscale rule create \
  --resource-group kanad-rg \
  --autoscale-name autoscale-kanad \
  --condition "CpuPercentage > 70 avg 5m" \
  --scale out 1

# Scale in when CPU < 30%
az monitor autoscale rule create \
  --resource-group kanad-rg \
  --autoscale-name autoscale-kanad \
  --condition "CpuPercentage < 30 avg 10m" \
  --scale in 1
```

## Workload-Based Recommendations

### MVP Testing (3 months)
- **Start with**: E4_v3 ($300/month)
- **Total cost**: ~$900 over 3 months (well within $1000 credits)
- **Scale to**: E8_v3 if needed for specific heavy workloads

### MVP with Startup Credits (11 months)
- **Primary tier**: E4_v3 ($300/month)
- **Peak tier**: E8_v3 for heavy computation days
- **Total cost**: ~$3300 over 11 months (within $5000 credits)
- **Buffer**: $1700 remaining for overages

### Production (After MVP)
- **Base tier**: E4_v3
- **Auto-scaling**: Enable to handle traffic spikes
- **Manual scaling**: E8_v3 for large molecule campaigns

## Memory Requirements by Molecule Size

Based on PySCF and Qiskit memory usage:

| Atoms | Basis Set | Estimated RAM | Recommended Tier |
|-------|-----------|---------------|------------------|
| 5-10 | sto-3g | 4-8 GB | E2_v3 (16GB) |
| 10-20 | sto-3g/3-21g | 8-16 GB | E2_v3/E4_v3 |
| 20-30 | 6-31g | 16-32 GB | E4_v3 (32GB) |
| 30-40 | 6-31g* | 32-64 GB | E8_v3 (64GB) |
| 40-50 | cc-pvdz | 64-128 GB | E16_v3 (128GB) |

## Current Configuration

The `deploy-azure.sh` script is configured with:

```bash
APP_SERVICE_SKU="E4_v3"  # 4 vCores, 32GB RAM
POSTGRES_SKU="GP_Standard_D4s_v3"  # 4 vCores, 16GB RAM
```

**This is optimized for**:
- Molecules up to 30 atoms
- Basis sets up to 6-31g*
- 2-3 concurrent users
- VQE and SQD calculations
- ~$455/month total cost

## Changing the Tier

### Before Deployment
Edit `deploy-azure.sh`:
```bash
# Line 46 - Change to desired tier
APP_SERVICE_SKU="E8_v3"  # For more power
```

### After Deployment
```bash
# Via CLI
az appservice plan update \
  --resource-group kanad-rg \
  --name kanad-plan \
  --sku E8_v3

# Via Azure Portal
# App Services → kanad-plan → Scale up → Select tier → Apply
```

## Monitoring and Alerts

Set up alerts to monitor resource usage:

```bash
# Alert when memory > 80%
az monitor metrics alert create \
  --name "High Memory Usage" \
  --resource-group kanad-rg \
  --scopes /subscriptions/{subscription-id}/resourceGroups/kanad-rg/providers/Microsoft.Web/sites/kanad-api \
  --condition "avg MemoryPercentage > 80" \
  --window-size 5m \
  --evaluation-frequency 1m

# Alert when CPU > 80%
az monitor metrics alert create \
  --name "High CPU Usage" \
  --resource-group kanad-rg \
  --scopes /subscriptions/{subscription-id}/resourceGroups/kanad-rg/providers/Microsoft.Web/sites/kanad-api \
  --condition "avg CpuPercentage > 80" \
  --window-size 5m \
  --evaluation-frequency 1m
```

## Summary

**For your MVP with $1000-5000 credits:**

1. **Recommended**: E4_v3 (32GB RAM, $300/month)
   - Handles most quantum chemistry workloads
   - 2.2 months on $1000 credits
   - 11 months on $5000 startup credits

2. **Budget option**: E2_v3 (16GB RAM, $150/month)
   - Good for testing and small molecules
   - 6.7 months on $1000 credits

3. **Performance option**: E8_v3 (64GB RAM, $600/month)
   - Use for specific heavy workloads
   - Scale up temporarily when needed

**The current configuration (E4_v3) is the best balance for your use case!**
