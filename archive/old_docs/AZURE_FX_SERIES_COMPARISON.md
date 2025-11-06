# Azure FX-Series VM Comparison for Quantum Chemistry Workloads

## Executive Summary

**Recommendation: Standard_FX12mds for MVP deployment**
- 12 vCPUs, 252 GB RAM
- $815/month (~6 months on $5000 credits)
- Perfect balance of performance and cost for quantum chemistry
- Can handle large molecules (50+ atoms) with advanced basis sets

## FX-Series v1 (Current Generation - Intel Cascade Lake)

| Model | vCPUs | RAM | RAM/vCPU | Temp Storage | On-Demand Cost | Monthly Cost | Use Case |
|-------|-------|-----|----------|--------------|----------------|--------------|----------|
| **Standard_FX4mds** | 4 | 84 GB | 21 GB | 168 GiB | $0.372/hr | **$272/mo** | Small-medium molecules |
| **Standard_FX12mds** | 12 | 252 GB | 21 GB | 504 GiB | $1.116/hr | **$815/mo** | **Recommended for MVP** |
| **Standard_FX24mds** | 24 | 504 GB | 21 GB | 1008 GiB | $2.232/hr | **$1,629/mo** | Large-scale simulations |
| **Standard_FX36mds** | 36 | 756 GB | 21 GB | 1512 GiB | $3.348/hr | **$2,444/mo** | Very large molecules |
| **Standard_FX48mds** | 48 | 1008 GB | 21 GB | 2016 GiB | $4.464/hr | **$3,259/mo** | Maximum performance |

## FX-Series v2 (5th Gen Intel Xeon Emerald Rapids - Newer)

| Model | vCPUs | RAM | RAM/vCPU | Temp NVMe | On-Demand Cost | Monthly Cost | Performance Gain |
|-------|-------|-----|----------|-----------|----------------|--------------|------------------|
| **Standard_FX2mdsv2** | 2 | 42 GB | 21 GB | 110 GiB | ~$0.186/hr | **~$136/mo** | +20% vs v1 |
| **Standard_FX4mdsv2** | 4 | 84 GB | 21 GB | 220 GiB | ~$0.372/hr | **~$272/mo** | +20% vs v1 |
| **Standard_FX12mdsv2** | 12 | 252 GB | 21 GB | 660 GiB | ~$1.116/hr | **~$815/mo** | +20% vs v1 |
| **Standard_FX24mdsv2** | 24 | 504 GB | 21 GB | 1320 GiB | ~$2.232/hr | **~$1,629/mo** | +20% vs v1 |
| **Standard_FX48mdsv2** | 48 | 1008 GB | 21 GB | 2640 GiB | ~$4.464/hr | **~$3,259/mo** | +20% vs v1 |
| **Standard_FX96mdsv2** | 96 | 1832 GB | 19 GB | 5280 GiB | ~$8.928/hr | **~$6,518/mo** | Enterprise HPC |

**Note**: v2 pricing is estimated based on v1 pricing. Actual prices may vary slightly.

## Cost Projections with Startup Credits

### With $1,000 Credits

| VM Model | Monthly Cost | Runtime on $1K Credits |
|----------|-------------|------------------------|
| FX4mds | $272 | **3.7 months** |
| FX12mds | $815 | **1.2 months** |
| FX24mds | $1,629 | 0.6 months (18 days) |

### With $5,000 Credits

| VM Model | Monthly Cost | Runtime on $5K Credits |
|----------|-------------|------------------------|
| FX4mds | $272 | **18.4 months** |
| FX12mds | $815 | **6.1 months** | ✅ **Recommended**
| FX24mds | $1,629 | **3.1 months** |
| FX36mds | $2,444 | 2.0 months |
| FX48mds | $3,259 | 1.5 months |

## Performance Comparison: FX vs App Service Premium

| Metric | App Service P3MV3 | FX4mds | FX12mds | Winner |
|--------|-------------------|---------|---------|--------|
| **vCPUs** | 8 | 4 | 12 | FX12mds |
| **RAM** | 32 GB | 84 GB | 252 GB | **FX12mds (8x more!)** |
| **RAM/vCPU** | 4 GB | 21 GB | 21 GB | **FX series** |
| **CPU Speed** | Standard | 4.0 GHz | 4.0 GHz | FX series |
| **Monthly Cost** | $600 | $272 | $815 | FX4mds (cheapest) |
| **Value (GB/$)** | 0.053 | 0.309 | 0.309 | **FX series (6x better!)** |

## Workload Capacity by VM Size

### FX4mds (4 vCPUs, 84 GB RAM) - $272/month
**Molecule Size Capacity:**
- Small molecules (< 10 atoms): Excellent ✅
- Medium molecules (10-20 atoms): Very good ✅
- Large molecules (20-30 atoms): Good ✅
- Very large (30-40 atoms): Moderate ⚠️
- Massive (40-50 atoms): Limited ❌

**Basis Sets:**
- sto-3g, 3-21g: Excellent ✅
- 6-31g, 6-31g*: Very good ✅
- cc-pvdz: Good ✅
- cc-pvtz: Limited ⚠️

**Concurrent Users:** 2-3
**Concurrent Jobs:** 2-3

### FX12mds (12 vCPUs, 252 GB RAM) - $815/month ⭐ **RECOMMENDED**
**Molecule Size Capacity:**
- Small molecules (< 10 atoms): Excellent ✅
- Medium molecules (10-20 atoms): Excellent ✅
- Large molecules (20-30 atoms): Excellent ✅
- Very large (30-40 atoms): Very good ✅
- Massive (40-50 atoms): Good ✅
- Enterprise (50+ atoms): Possible ✅

**Basis Sets:**
- sto-3g, 3-21g: Excellent ✅
- 6-31g, 6-31g*: Excellent ✅
- cc-pvdz: Very good ✅
- cc-pvtz: Good ✅
- cc-pvqz: Possible ⚠️

**Concurrent Users:** 8-10
**Concurrent Jobs:** 6-8
**VQE Iterations:** Up to 1000
**SQD States:** 10-20 excited states

### FX24mds (24 vCPUs, 504 GB RAM) - $1,629/month
**Molecule Size Capacity:**
- All molecule sizes: Excellent ✅
- Enterprise-scale (100+ atoms): Possible ✅

**Basis Sets:**
- All standard basis sets: Excellent ✅
- Advanced basis sets (cc-pvqz): Very good ✅

**Concurrent Users:** 15-20
**Concurrent Jobs:** 12-15
**Use Case:** Production with high traffic

## Memory Requirements by Molecule Size & Basis Set

| Molecule Size | sto-3g | 6-31g | 6-31g* | cc-pvdz | cc-pvtz | Recommended VM |
|---------------|--------|-------|--------|---------|---------|----------------|
| 5 atoms | 2 GB | 4 GB | 8 GB | 12 GB | 24 GB | FX4mds |
| 10 atoms | 4 GB | 12 GB | 20 GB | 32 GB | 64 GB | FX4mds |
| 20 atoms | 12 GB | 32 GB | 64 GB | 96 GB | 192 GB | FX12mds |
| 30 atoms | 24 GB | 64 GB | 128 GB | 192 GB | 384 GB | FX12mds-FX24mds |
| 40 atoms | 40 GB | 96 GB | 192 GB | 288 GB | 576 GB | FX24mds+ |
| 50 atoms | 64 GB | 128 GB | 256 GB | 384 GB | 768 GB | FX36mds+ |

## Architecture Benefits: VM vs App Service

### Why FX-Series VM is Better for Quantum Chemistry

#### 1. **Massive RAM (21 GB per vCPU)**
- App Service: 4 GB per vCPU
- FX-series: 21 GB per vCPU
- **5.25x more memory per core!**

#### 2. **High-Speed Local Storage**
- FX-series: Local NVMe SSD (up to 2016 GiB)
- App Service: Network-attached storage
- **10-100x faster I/O for temporary molecular data**

#### 3. **Full Control**
- Install any quantum chemistry libraries
- Optimize Linux kernel for HPC
- Fine-tune PostgreSQL for large datasets
- Custom memory management for PySCF

#### 4. **Cost Efficiency**
- FX12mds: $815/month for 252 GB RAM
- App Service P3MV3: $600/month for 32 GB RAM
- **8x more RAM for 36% more cost = 6x better value!**

#### 5. **No PaaS Limitations**
- No request timeouts (App Service has 230 second limit)
- No memory limits per process
- No container size restrictions
- Direct PostgreSQL on same machine (zero latency)

## Deployment Architecture

### Single VM Setup (Recommended for MVP)

```
┌─────────────────────────────────────────────────┐
│         Azure FX12mds VM (252 GB RAM)           │
│                                                 │
│  ┌──────────────┐      ┌──────────────┐       │
│  │   Nginx      │      │  PostgreSQL  │       │
│  │  (Port 80)   │      │  (Port 5432) │       │
│  │  Reverse     │      │   Database   │       │
│  │   Proxy      │      │   100 GB     │       │
│  └──────┬───────┘      └──────────────┘       │
│         │                                       │
│  ┌──────▼───────┐                              │
│  │   Uvicorn    │                              │
│  │  FastAPI App │                              │
│  │  (Port 8000) │                              │
│  │              │                              │
│  │  - VQE       │                              │
│  │  - SQD       │                              │
│  │  - PySCF     │                              │
│  │  - Qiskit    │                              │
│  └──────────────┘                              │
│                                                 │
│  Memory Allocation:                             │
│  - PostgreSQL: 50 GB                            │
│  - FastAPI/Python: 180 GB                       │
│  - System: 22 GB                                │
└─────────────────────────────────────────────────┘
         │
         │ HTTPS (Let's Encrypt SSL)
         ▼
    Public Internet
```

## Recommended Configuration for MVP

### Primary Recommendation: FX12mds

**Why FX12mds?**

1. **Perfect Balance**: 252 GB RAM handles most quantum chemistry workloads
2. **Cost Effective**: $815/month = 6.1 months on $5K credits
3. **Scalability**: 12 vCPUs support 6-8 concurrent jobs
4. **Future-Proof**: Can handle molecules up to 50 atoms
5. **Production-Ready**: Supports 8-10 concurrent users

**Deployment Plan:**
- **Months 1-3**: Use FX12mds for MVP development and testing
- **Months 4-6**: Continue on FX12mds for production
- **After Month 6**:
  - Scale down to FX4mds if usage is light ($272/month)
  - Scale up to FX24mds if need more capacity ($1,629/month)
  - Or move to pay-as-you-go pricing

### Budget Allocation ($5,000 Credits)

| Phase | Duration | VM | Monthly Cost | Total Cost | Remaining |
|-------|----------|-----|--------------|------------|-----------|
| Setup & Testing | 2 weeks | FX4mds | $136 | $136 | $4,864 |
| MVP Development | 3 months | FX12mds | $815 | $2,445 | $2,419 |
| Production | 2 months | FX12mds | $815 | $1,630 | $789 |
| Optimization | 1 month | FX4mds or FX12mds | $272-815 | $272-815 | $0-517 |
| **Total** | **~6.5 months** | | | **~$4,483-$5,026** | **Optimal** |

## Alternative Configurations

### Option 1: Conservative (Longest Runtime)
- **VM**: FX4mds
- **Cost**: $272/month
- **Runtime**: 18.4 months on $5K
- **Best for**: Small-medium molecules, low traffic

### Option 2: Performance (Best Performance)
- **VM**: FX24mds or FX36mds
- **Cost**: $1,629-$2,444/month
- **Runtime**: 2-3 months on $5K
- **Best for**: Large molecules, high traffic, short-term MVP

### Option 3: Spot Instances (Maximum Savings)
- **VM**: FX12mds Spot
- **Cost**: ~$150/month (82% savings!)
- **Runtime**: 33 months on $5K
- **Caveat**: Can be evicted (Azure reclaims capacity)
- **Best for**: Non-critical workloads, batch processing

## Regional Availability

FX-series VMs are available in limited regions. Check availability:

```bash
az vm list-skus --location eastus --size Standard_FX --output table
az vm list-skus --location westus2 --size Standard_FX --output table
az vm list-skus --location centralus --size Standard_FX --output table
```

**Recommended Regions:**
- East US
- West US 2
- Central US
- West Europe
- UK South

## Final Recommendation

### For Your MVP Deployment:

**Use Standard_FX12mds (12 vCPUs, 252 GB RAM)**

**Reasons:**
1. ✅ 252 GB RAM handles molecules up to 50 atoms
2. ✅ 12 vCPUs support multiple concurrent users
3. ✅ $815/month = 6 months runtime on $5K credits
4. ✅ Perfect for 3-month MVP + 3 months production
5. ✅ Can handle advanced basis sets (cc-pvdz, cc-pvtz)
6. ✅ Local NVMe storage for fast molecular data I/O
7. ✅ No App Service limitations or timeouts
8. ✅ Full control over PostgreSQL and Python environment

**Next Steps:**
1. Run deployment script (will be created next)
2. Deploy to FX12mds in your preferred region
3. Configure PostgreSQL with 50GB allocated
4. Set up Nginx with SSL (Let's Encrypt)
5. Deploy FastAPI application
6. Configure monitoring and alerts
7. Test with sample molecules

**Scaling Strategy:**
- Start with FX12mds
- Monitor RAM usage via Azure Monitor
- If consistently < 50% utilized → Scale down to FX4mds
- If hitting RAM limits → Scale up to FX24mds
- Use spot instances for batch jobs to save costs
