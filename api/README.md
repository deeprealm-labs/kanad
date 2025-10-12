# Kanad API Backend

FastAPI backend server for the Kanad quantum chemistry framework.

## Features

- 🔬 **Quantum Chemistry Calculations**: HF, VQE, SQD solvers
- ☁️ **Cloud Backends**: IBM Quantum, BlueQubit support
- 📊 **Experiment Management**: Queue, track, and analyze experiments
- 🔄 **Real-time Updates**: Background job processing with progress tracking
- 💾 **Persistent Storage**: SQLite database for experiments and settings

## Quick Start

### 1. Install Dependencies

```bash
# Install API dependencies
cd api
pip install -r requirements.txt

# Kanad framework should already be installed
cd ..
pip install -e .
```

### 2. Set Environment Variables (Optional)

```bash
# IBM Quantum credentials
export IBM_API="your_ibm_api_token"
export IBM_CRN="your_ibm_crn"

# BlueQubit credentials
export BLUE_TOKEN="your_bluequbit_token"

# Server configuration
export API_HOST="0.0.0.0"
export API_PORT="8000"
export DEBUG="true"
```

### 3. Start the Server

```bash
# From the api directory
cd api
python main.py

# Or using uvicorn directly
uvicorn main:app --reload --host 0.0.0.0 --port 8000
```

The server will start at http://localhost:8000

### 4. View API Documentation

Open your browser and navigate to:
- Swagger UI: http://localhost:8000/docs
- ReDoc: http://localhost:8000/redoc

## API Endpoints

### Health Check
- `GET /health` - Server health status
- `GET /` - API information

### Molecules
- `POST /api/molecules/create` - Create molecule from SMILES or atoms
- `POST /api/molecules/validate-smiles` - Validate SMILES string
- `POST /api/molecules/bond-info` - Get bond type information

### Experiments
- `POST /api/experiments/submit` - Submit new experiment
- `GET /api/experiments/list` - List all experiments
- `GET /api/experiments/{id}` - Get experiment details
- `GET /api/experiments/{id}/status` - Get experiment status
- `GET /api/experiments/{id}/results` - Get experiment results
- `POST /api/experiments/{id}/cancel` - Cancel experiment
- `DELETE /api/experiments/{id}` - Delete experiment

### Jobs
- `GET /api/jobs/list` - List jobs in queue
- `GET /api/jobs/{id}` - Get job details
- `GET /api/jobs/{id}/status` - Get job status and progress
- `POST /api/jobs/{id}/cancel` - Cancel job
- `DELETE /api/jobs/{id}` - Delete job

### Settings
- `GET /api/settings/defaults` - Get default settings
- `PUT /api/settings/defaults` - Update default settings

### Library
- `GET /api/library/molecules` - Get molecule library
- `GET /api/library/molecules/{id}` - Get specific molecule

### Cloud
- `GET /api/cloud/backends` - List available backends
- `POST /api/cloud/credentials` - Store cloud credentials
- `GET /api/cloud/credentials/{provider}` - Check credential status

## Example Usage

### Submit a VQE Experiment

```bash
curl -X POST "http://localhost:8000/api/experiments/submit" \
  -H "Content-Type: application/json" \
  -d '{
    "molecule": {
      "smiles": "[H][H]",
      "basis": "sto-3g",
      "charge": 0,
      "multiplicity": 1
    },
    "configuration": {
      "method": "VQE",
      "ansatz": "hardware_efficient",
      "mapper": "jordan_wigner",
      "optimizer": "SLSQP",
      "max_iterations": 1000,
      "backend": "classical"
    },
    "execute_now": true
  }'
```

### Get Experiment Results

```bash
curl "http://localhost:8000/api/experiments/{experiment_id}/results"
```

## Database

The API uses SQLite for data persistence. The database file is created at:
```
kanad_experiments.db
```

Tables:
- `experiments` - Experiment records
- `jobs` - Job queue
- `user_settings` - User preferences
- `cloud_credentials` - Cloud provider credentials

## Architecture

```
api/
├── main.py                    # FastAPI app entry point
├── core/
│   ├── config.py             # Configuration settings
│   └── database.py           # Database operations
├── routes/
│   ├── health.py             # Health check
│   ├── molecules.py          # Molecule endpoints
│   ├── experiments.py        # Experiment endpoints
│   ├── jobs.py               # Job queue endpoints
│   ├── analysis.py           # Analysis endpoints
│   ├── settings.py           # Settings endpoints
│   ├── library.py            # Molecule library
│   └── cloud.py              # Cloud backend management
└── services/
    └── experiment_service.py # Experiment execution logic
```

## Frontend Integration

The API is designed to work with the Next.js frontend in `/web`.

Update the frontend's API URL:
```bash
# In /web/.env.local
NEXT_PUBLIC_API_URL=http://localhost:8000/api
```

## Development

### Running Tests

```bash
pytest api/tests/
```

### Code Style

```bash
black api/
flake8 api/
```

## Production Deployment

For production deployment:

1. Set `DEBUG=false`
2. Use a production-grade database (PostgreSQL)
3. Add authentication middleware
4. Set up HTTPS/TLS
5. Use a process manager (systemd, supervisor)
6. Add rate limiting
7. Set up monitoring and logging

## Troubleshooting

### Import Errors

Make sure the Kanad framework is installed:
```bash
pip install -e .
```

### Backend Connection Issues

Check that cloud credentials are properly set:
```bash
echo $IBM_API
echo $BLUE_TOKEN
```

### Database Locked

If you get "database is locked" errors, make sure only one server instance is running.

## License

Same as Kanad framework

## Support

For issues and questions, please open an issue on GitHub.
