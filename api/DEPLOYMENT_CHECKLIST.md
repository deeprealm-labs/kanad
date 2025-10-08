# Kanad API - Deployment Checklist

## Pre-Deployment Steps

### 1. Environment Setup

- [ ] Activate Python virtual environment
  ```bash
  cd /home/mk/deeprealm/kanad
  source venv/bin/activate  # or your venv path
  ```

- [ ] Install API dependencies
  ```bash
  cd /home/mk/deeprealm/kanad/api
  pip install -r requirements.txt
  ```

- [ ] Install Kanad framework
  ```bash
  cd /home/mk/deeprealm/kanad
  pip install -e .
  ```

### 2. Database Setup

- [ ] Run cloud credentials migration
  ```bash
  cd /home/mk/deeprealm/kanad/api
  python migrations/add_cloud_credentials.py
  ```

- [ ] Verify database exists
  ```bash
  ls -la /home/mk/deeprealm/kanad/kanad.db
  ```

### 3. Configuration

- [ ] Create `.env` file (optional)
  ```bash
  cd /home/mk/deeprealm/kanad/api
  cat > .env << EOF
  DEBUG=True
  IBM_QUANTUM_TOKEN=your_token_here
  IBM_QUANTUM_CRN=your_crn_here
  BLUEQUBIT_API_KEY=your_key_here
  EOF
  ```

- [ ] Update CORS origins in `config.py` if needed
  ```python
  CORS_ORIGINS = [
      "http://localhost:3000",
      "https://your-frontend-domain.com"
  ]
  ```

### 4. Testing

- [ ] Start the server
  ```bash
  cd /home/mk/deeprealm/kanad/api
  uvicorn main:app --reload --port 8000
  ```

- [ ] Verify health endpoint
  ```bash
  curl http://localhost:8000/health
  ```

- [ ] Test info endpoint
  ```bash
  curl http://localhost:8000/api/v1/info | jq
  ```

- [ ] Run comprehensive tests
  ```bash
  python test_complete_api.py
  ```

### 5. Cloud Credentials Configuration

**Option A: Database Storage (Recommended)**

- [ ] Configure IBM Quantum credentials
  ```bash
  curl -X POST http://localhost:8000/api/v1/cloud-credentials/ibm \
    -H "Content-Type: application/json" \
    -d '{"crn": "your-crn", "api_key": "your-key"}'
  ```

- [ ] Configure BlueQubit credentials
  ```bash
  curl -X POST http://localhost:8000/api/v1/cloud-credentials/bluequbit \
    -H "Content-Type: application/json" \
    -d '{"api_token": "your-token"}'
  ```

- [ ] Verify credentials
  ```bash
  curl http://localhost:8000/api/v1/cloud-credentials/status
  ```

**Option B: Environment Variables (Fallback)**

- [ ] Set environment variables
  ```bash
  export IBM_QUANTUM_TOKEN="your-token"
  export IBM_QUANTUM_CRN="your-crn"
  export BLUEQUBIT_API_KEY="your-key"
  ```

## Production Deployment

### Security Checklist

- [ ] Enable HTTPS/TLS
  - Update Uvicorn command with SSL certificates
  - Or use reverse proxy (nginx/Apache)

- [ ] Implement credential encryption
  ```python
  # In cloud_credentials.py, add encryption
  from cryptography.fernet import Fernet
  ```

- [ ] Disable DEBUG mode
  ```python
  # In config.py
  DEBUG = False
  ```

- [ ] Set secure CORS origins
  ```python
  CORS_ORIGINS = ["https://your-frontend.com"]
  ```

- [ ] Implement rate limiting
  ```python
  from slowapi import Limiter
  # Add to main.py
  ```

- [ ] Add authentication/authorization
  ```python
  from fastapi.security import OAuth2PasswordBearer
  ```

### Performance Checklist

- [ ] Configure job queue workers
  ```python
  # In config.py
  MAX_CONCURRENT_JOBS = 4  # Based on server capacity
  ```

- [ ] Set appropriate timeouts
  ```python
  JOB_TIMEOUT_SECONDS = 3600  # 1 hour
  ```

- [ ] Enable database connection pooling
  ```python
  # In database.py
  pool_size=10, max_overflow=20
  ```

### Monitoring & Logging

- [ ] Configure structured logging
  ```python
  import logging.config
  # Setup JSON logging for production
  ```

- [ ] Set up monitoring
  - Application metrics (Prometheus/Grafana)
  - Error tracking (Sentry)
  - Uptime monitoring

- [ ] Configure log rotation
  ```bash
  # /etc/logrotate.d/kanad-api
  ```

### Docker Deployment (Optional)

- [ ] Create Dockerfile
  ```dockerfile
  FROM python:3.10-slim
  WORKDIR /app
  COPY requirements.txt .
  RUN pip install -r requirements.txt
  COPY . .
  CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
  ```

- [ ] Build image
  ```bash
  docker build -t kanad-api:latest .
  ```

- [ ] Run container
  ```bash
  docker run -p 8000:8000 \
    -e IBM_QUANTUM_TOKEN=$IBM_TOKEN \
    -v ./kanad.db:/app/kanad.db \
    kanad-api:latest
  ```

### Kubernetes Deployment (Optional)

- [ ] Create ConfigMap
  ```yaml
  apiVersion: v1
  kind: ConfigMap
  metadata:
    name: kanad-config
  ```

- [ ] Create Secret
  ```yaml
  apiVersion: v1
  kind: Secret
  metadata:
    name: kanad-secrets
  type: Opaque
  data:
    ibm-token: <base64-encoded>
  ```

- [ ] Create Deployment
  ```yaml
  apiVersion: apps/v1
  kind: Deployment
  metadata:
    name: kanad-api
  spec:
    replicas: 2
  ```

- [ ] Create Service
  ```yaml
  apiVersion: v1
  kind: Service
  metadata:
    name: kanad-api-service
  ```

## Post-Deployment Verification

- [ ] Test all endpoints
  ```bash
  python test_complete_api.py
  ```

- [ ] Create a test experiment
  ```bash
  curl -X POST $API_URL/api/v1/experiments/ \
    -H "Content-Type: application/json" \
    -d @test_experiment.json
  ```

- [ ] Verify cloud backends work
  - [ ] IBM Quantum job submission
  - [ ] BlueQubit job submission

- [ ] Check monitoring dashboards
  - [ ] API response times
  - [ ] Error rates
  - [ ] Job queue metrics

- [ ] Review logs
  ```bash
  tail -f /var/log/kanad-api.log
  ```

## Rollback Plan

If issues occur:

- [ ] Revert to previous version
  ```bash
  git checkout <previous-commit>
  ```

- [ ] Restore database backup
  ```bash
  cp kanad.db.backup kanad.db
  ```

- [ ] Restart services
  ```bash
  systemctl restart kanad-api
  ```

## Backup Strategy

- [ ] Database backups
  ```bash
  # Daily cron job
  0 2 * * * cp /path/to/kanad.db /backups/kanad.db.$(date +\%Y\%m\%d)
  ```

- [ ] Configuration backups
  ```bash
  tar -czf config-backup.tar.gz .env config.py
  ```

- [ ] Code repository
  ```bash
  git push origin main
  ```

## Support Information

**Documentation:**
- Implementation Report: `/home/mk/deeprealm/kanad/api/IMPLEMENTATION_REPORT.md`
- Completion Summary: `/home/mk/deeprealm/kanad/api/COMPLETION_SUMMARY.md`

**Endpoints:**
- API Info: `GET /api/v1/info`
- Health Check: `GET /health`
- API Docs: `GET /docs`
- ReDoc: `GET /redoc`

**Key Files:**
- Main App: `/home/mk/deeprealm/kanad/api/main.py`
- Config: `/home/mk/deeprealm/kanad/api/config.py`
- Database: `/home/mk/deeprealm/kanad/kanad.db`

**Contact:**
- GitHub Issues: [Create issue with logs]
- Documentation: [Link to docs]

---

**Last Updated:** 2025-10-08
**API Version:** 1.0.0
**Status:** ✅ Ready for Deployment
