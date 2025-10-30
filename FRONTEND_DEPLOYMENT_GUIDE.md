# Kanad Frontend Deployment Guide

## Overview

This guide covers deploying the Kanad Next.js frontend application to Vercel, connecting it to the Azure backend API.

## Architecture

```
┌─────────────────────────────────────────────┐
│         Vercel (Frontend Hosting)           │
│                                             │
│  ┌─────────────────────────────────────┐  │
│  │   Next.js App (Static + SSR)        │  │
│  │   - React Components                │  │
│  │   - TailwindCSS                     │  │
│  │   - Recharts (3D Visualization)     │  │
│  │   - WebSocket Support               │  │
│  └──────────────┬──────────────────────┘  │
│                 │                           │
└─────────────────┼───────────────────────────┘
                  │ HTTPS API Calls
                  ▼
┌─────────────────────────────────────────────┐
│    Azure FX16-4mds_v2 VM (Backend)          │
│                                             │
│  ┌─────────────────────────────────────┐  │
│  │   Nginx (172.171.222.16:80)         │  │
│  │   ↓                                  │  │
│  │   FastAPI (localhost:8000)          │  │
│  │   - Quantum Chemistry API           │  │
│  │   - PostgreSQL Database             │  │
│  └─────────────────────────────────────┘  │
└─────────────────────────────────────────────┘
```

## Prerequisites

### 1. Vercel Account
- Sign up at https://vercel.com
- Install Vercel CLI: `npm install -g vercel`
- Login: `vercel login`

### 2. Azure Backend
- Azure VM deployed and running: **172.171.222.16**
- API accessible at: http://172.171.222.16/api
- Health check passing: http://172.171.222.16/health

### 3. Environment Setup
- Node.js 20+ installed
- Git repository connected to Vercel

## Environment Configuration

### Local Development (.env.local)
```bash
NEXT_PUBLIC_API_URL=http://localhost:8000/api
```

### Production (.env.production)
```bash
NEXT_PUBLIC_API_URL=http://172.171.222.16/api
```

## Deployment Steps

### Option 1: Deploy via Vercel CLI (Recommended)

#### Step 1: Navigate to web directory
```bash
cd web
```

#### Step 2: Build and test locally
```bash
npm run build
npm run start
# Test at http://localhost:3000
```

#### Step 3: Deploy to Vercel
```bash
# First deployment (creates project)
vercel

# Production deployment
vercel --prod
```

#### Step 4: Configure environment variables
```bash
# Set API URL
vercel env add NEXT_PUBLIC_API_URL
# Enter: http://172.171.222.16/api
# Select: Production
```

### Option 2: Deploy via Vercel Dashboard

#### Step 1: Connect Git Repository
1. Go to https://vercel.com/dashboard
2. Click "Import Project"
3. Select your GitHub repository
4. Set root directory to `web`

#### Step 2: Configure Build Settings
- **Framework Preset**: Next.js
- **Build Command**: `npm run build`
- **Output Directory**: `.next`
- **Install Command**: `npm install`

#### Step 3: Add Environment Variables
Go to Settings → Environment Variables:

| Key | Value | Environment |
|-----|-------|-------------|
| `NEXT_PUBLIC_API_URL` | `http://172.171.222.16/api` | Production |

#### Step 4: Deploy
Click "Deploy" - Vercel will automatically build and deploy.

## Post-Deployment Configuration

### 1. Update CORS on Azure Backend

SSH into Azure VM:
```bash
ssh kanadmin@172.171.222.16
```

Update backend CORS settings:
```bash
cd /opt/kanad
nano .env
```

Add your Vercel URL to CORS_ORIGINS:
```env
CORS_ORIGINS=["http://172.171.222.16","http://localhost:3000","https://your-app.vercel.app"]
```

Restart backend:
```bash
sudo systemctl restart kanad
```

### 2. Configure Custom Domain (Optional)

#### On Vercel:
1. Go to Project Settings → Domains
2. Add your domain: `kanad.yourdomain.com`
3. Follow DNS instructions

#### On Azure:
Update CORS to include your custom domain:
```env
CORS_ORIGINS=["https://kanad.yourdomain.com"]
```

### 3. Set up SSL/HTTPS on Azure (Recommended)

For production, enable HTTPS on the backend:

```bash
# SSH into Azure VM
ssh kanadmin@172.171.222.16

# Install Certbot
sudo apt-get install -y certbot python3-certbot-nginx

# Get SSL certificate (requires domain)
sudo certbot --nginx -d api.yourdomain.com
```

Update frontend to use HTTPS:
```bash
# On Vercel
vercel env add NEXT_PUBLIC_API_URL
# Enter: https://api.yourdomain.com/api
```

## Vercel Configuration (vercel.json)

The project includes a pre-configured `vercel.json`:

```json
{
  "buildCommand": "npm run build",
  "framework": "nextjs",
  "env": {
    "NEXT_PUBLIC_API_URL": "@kanad_api_url"
  },
  "headers": [
    {
      "source": "/:path*",
      "headers": [
        {"key": "X-Frame-Options", "value": "DENY"},
        {"key": "X-Content-Type-Options", "value": "nosniff"}
      ]
    }
  ]
}
```

## Testing Deployment

### 1. Health Check
```bash
# Test backend
curl http://172.171.222.16/health

# Test API
curl http://172.171.222.16/api/experiments/list?limit=5
```

### 2. Frontend Tests
Visit your Vercel URL and test:
- [ ] Landing page loads
- [ ] Dashboard accessible
- [ ] Can create experiment
- [ ] API calls work
- [ ] 3D visualization renders
- [ ] WebSocket connections work

### 3. Network Tests
Check browser console for:
- ✅ No CORS errors
- ✅ API calls succeed
- ✅ WebSocket connects
- ✅ Assets load properly

## Troubleshooting

### CORS Errors

**Symptom**: Console shows CORS policy errors

**Solution**:
```bash
# On Azure VM
cd /opt/kanad
nano .env

# Add Vercel URL to CORS_ORIGINS
CORS_ORIGINS=["https://your-app.vercel.app","http://172.171.222.16"]

# Restart
sudo systemctl restart kanad
```

### API Connection Failed

**Symptom**: "Failed to fetch" or timeout errors

**Check backend status**:
```bash
ssh kanadmin@172.171.222.16
sudo systemctl status kanad
sudo journalctl -u kanad -n 50
```

**Check connectivity**:
```bash
# From local machine
curl -v http://172.171.222.16/health
```

### Build Failures on Vercel

**Check build logs** on Vercel dashboard

**Common issues**:
1. **Missing dependencies**: Ensure `package.json` is correct
2. **TypeScript errors**: Run `npm run build` locally first
3. **Environment variables**: Check they're set in Vercel

### WebSocket Connection Issues

**Symptom**: Real-time updates don't work

**Solution**:
WebSockets may not work over HTTP. Options:
1. Enable HTTPS on Azure backend
2. Use polling fallback (already implemented)

## Performance Optimization

### 1. Enable Image Optimization
In `next.config.js`:
```javascript
module.exports = {
  images: {
    domains: ['172.171.222.16'],
  },
}
```

### 2. Enable Caching
Vercel automatically caches:
- Static pages
- API responses (with proper headers)
- Images

### 3. Monitor Performance
- Use Vercel Analytics
- Check Core Web Vitals
- Monitor API response times

## Cost Analysis

### Vercel Hosting
- **Free Tier**: 100 GB bandwidth, unlimited requests
- **Pro Tier**: $20/month (recommended for production)
  - 1 TB bandwidth
  - Advanced analytics
  - Custom domains
  - Team collaboration

### Total Monthly Cost
| Service | Tier | Cost |
|---------|------|------|
| Azure VM (Backend) | FX16-4mds_v2 | $1,805 |
| Vercel (Frontend) | Pro | $20 |
| **Total** | | **$1,825/month** |

**With $5,000 startup credits**: ~2.7 months of operation

## Monitoring & Maintenance

### Vercel Monitoring
- View deployments: `vercel ls`
- View logs: `vercel logs`
- Check analytics: Vercel Dashboard

### Backend Monitoring
```bash
# SSH into Azure VM
ssh kanadmin@172.171.222.16

# Check services
sudo systemctl status kanad nginx postgresql

# View logs
sudo journalctl -u kanad -f
```

## Rollback Procedure

### Rollback Frontend
```bash
# List deployments
vercel ls

# Promote previous deployment
vercel promote <deployment-url>
```

### Rollback Backend
```bash
# SSH into Azure VM
ssh kanadmin@172.171.222.16

# Revert code
cd /opt/kanad
git log --oneline
git reset --hard <previous-commit>

# Restart
sudo systemctl restart kanad
```

## Security Best Practices

### 1. Environment Variables
- ✅ Never commit `.env` files
- ✅ Use Vercel environment variables
- ✅ Rotate secrets regularly

### 2. CORS Configuration
- ✅ Specify exact origins (not `*`)
- ✅ Use HTTPS in production
- ✅ Limit allowed methods

### 3. Rate Limiting
Backend already implements:
- Session limits (5 per user)
- Request rate limiting
- Compute resource limits

## Scaling Strategy

### Frontend Scaling (Vercel)
- **Automatic**: Vercel scales edge functions automatically
- **CDN**: Global edge network included
- **No config needed**: Handles traffic spikes automatically

### Backend Scaling (Azure)
When needed, scale the VM:
```bash
# Scale up to FX32-8mds_v2 (672 GB RAM)
az vm deallocate --resource-group kanad-vm-rg --name kanad-fx16
az vm resize --resource-group kanad-vm-rg --name kanad-fx16 --size Standard_FX32-8mds_v2
az vm start --resource-group kanad-vm-rg --name kanad-fx16
```

## Next Steps After Deployment

1. **Set up monitoring**
   - Enable Vercel Analytics
   - Configure Azure Monitor
   - Set up error tracking (Sentry)

2. **Configure alerts**
   - Deployment notifications
   - Error rate alerts
   - Performance degradation alerts

3. **Document API URL**
   - Share with team
   - Update documentation
   - Configure in mobile apps (if any)

4. **User testing**
   - Test all features
   - Verify quantum calculations
   - Check 3D visualizations

## Support Resources

- **Vercel Docs**: https://vercel.com/docs
- **Next.js Docs**: https://nextjs.org/docs
- **Azure Docs**: https://docs.microsoft.com/azure
- **Project README**: See root directory

## Summary

You now have:
- ✅ Frontend deployed on Vercel
- ✅ Backend running on Azure VM
- ✅ CORS configured
- ✅ Environment variables set
- ✅ Production-ready infrastructure

**Frontend URL**: https://your-app.vercel.app
**Backend API**: http://172.171.222.16/api
**Health Check**: http://172.171.222.16/health

Total cost: **~$1,825/month** (~2.7 months on $5K credits)

🚀 **Ready for MVP deployment!**

---

*Generated: 2025-10-30*
*Frontend: Vercel*
*Backend: Azure FX16-4mds_v2*
*Status: ✅ Ready to Deploy*
