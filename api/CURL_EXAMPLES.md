# Kanad API - cURL Examples

Complete collection of cURL commands for testing all API endpoints.

## Setup

```bash
# Set base URL variable
export API_URL="http://localhost:8000/api/v1"
```

## Health & Info

### Health Check
```bash
curl http://localhost:8000/health
```

### API Information
```bash
curl $API_URL/info
```

## Molecules

### Validate SMILES (Valid)
```bash
curl -X POST $API_URL/molecules/validate \
  -H "Content-Type: application/json" \
  -d '{"smiles": "CCO"}'
```

### Validate SMILES (Invalid)
```bash
curl -X POST $API_URL/molecules/validate \
  -H "Content-Type: application/json" \
  -d '{"smiles": "invalid_smiles"}'
```

### Get Molecule Library
```bash
curl $API_URL/molecules/library
```

### Get Molecule Categories
```bash
curl $API_URL/molecules/library/categories
```

### Get Specific Molecule
```bash
curl $API_URL/molecules/library/h2
```

### Filter Library by Category
```bash
curl "$API_URL/molecules/library?category=hydrocarbons"
```

## Settings

### Get Settings
```bash
curl $API_URL/settings/
```

### Update Settings
```bash
curl -X PUT $API_URL/settings/ \
  -H "Content-Type: application/json" \
  -d '{
    "method": "VQE",
    "ansatz": "hardware_efficient",
    "mapper": "bravyi_kitaev",
    "optimizer": "COBYLA",
    "backend": "classical",
    "circuit_optimization": true,
    "adaptive_vqe": false
  }'
```

### Reset Settings
```bash
curl -X DELETE $API_URL/settings/
```

## Experiments

### 1. Create H2 Experiment (VQE with UCC)
```bash
curl -X POST $API_URL/experiments/ \
  -H "Content-Type: application/json" \
  -d '{
    "name": "H2 VQE with UCC Ansatz",
    "molecule": {
      "smiles": "[H][H]",
      "basis": "sto-3g",
      "charge": 0,
      "multiplicity": 1
    },
    "configuration": {
      "method": "VQE",
      "ansatz": "ucc",
      "mapper": "jordan_wigner",
      "optimizer": "SLSQP",
      "backend": "classical",
      "max_iterations": 1000,
      "conv_threshold": 1e-6
    },
    "execute_immediately": true
  }'
```

### 2. Create Water Experiment (VQE with Hardware-Efficient)
```bash
curl -X POST $API_URL/experiments/ \
  -H "Content-Type: application/json" \
  -d '{
    "name": "H2O VQE Calculation",
    "molecule": {
      "smiles": "O",
      "basis": "sto-3g",
      "charge": 0,
      "multiplicity": 1
    },
    "configuration": {
      "method": "VQE",
      "ansatz": "hardware_efficient",
      "mapper": "jordan_wigner",
      "optimizer": "COBYLA",
      "backend": "classical",
      "max_iterations": 500
    },
    "execute_immediately": true
  }'
```

### 3. Create Ethanol Experiment (without immediate execution)
```bash
curl -X POST $API_URL/experiments/ \
  -H "Content-Type: application/json" \
  -d '{
    "name": "Ethanol VQE",
    "molecule": {
      "smiles": "CCO",
      "basis": "6-31g",
      "charge": 0,
      "multiplicity": 1
    },
    "configuration": {
      "method": "VQE",
      "ansatz": "ucc",
      "mapper": "bravyi_kitaev",
      "optimizer": "SLSQP",
      "backend": "classical"
    },
    "execute_immediately": false
  }'
```

### 4. Create HF Calculation
```bash
curl -X POST $API_URL/experiments/ \
  -H "Content-Type: application/json" \
  -d '{
    "name": "Methane HF",
    "molecule": {
      "smiles": "C",
      "basis": "sto-3g",
      "charge": 0,
      "multiplicity": 1
    },
    "configuration": {
      "method": "HF",
      "backend": "classical"
    },
    "execute_immediately": true
  }'
```

### List All Experiments
```bash
curl $API_URL/experiments/
```

### List Completed Experiments
```bash
curl "$API_URL/experiments/?status=completed"
```

### List with Pagination
```bash
curl "$API_URL/experiments/?limit=10&offset=0"
```

### Get Experiment Details
```bash
# Replace {id} with actual experiment ID
curl $API_URL/experiments/1
```

### Get Experiment Status
```bash
curl $API_URL/experiments/1/status
```

### Get Convergence Data
```bash
curl $API_URL/experiments/1/convergence
```

### Delete Experiment
```bash
curl -X DELETE $API_URL/experiments/1
```

## Queue Management

### Add Experiment to Queue
```bash
# First create experiment without immediate execution, then:
curl -X POST $API_URL/queue/ \
  -H "Content-Type: application/json" \
  -d '{
    "experiment_id": 3,
    "priority": 5
  }'
```

### Schedule Experiment for Later
```bash
curl -X POST $API_URL/queue/ \
  -H "Content-Type: application/json" \
  -d '{
    "experiment_id": 3,
    "priority": 10,
    "scheduled_time": "2025-10-08T20:00:00"
  }'
```

### List Queue
```bash
curl $API_URL/queue/
```

### List Queued Items Only
```bash
curl "$API_URL/queue/?status=queued"
```

### Get Queue Item
```bash
curl $API_URL/queue/1
```

### Update Queue Priority
```bash
curl -X PUT $API_URL/queue/1 \
  -H "Content-Type: application/json" \
  -d '{"priority": 100}'
```

### Pause Queue Item
```bash
curl -X PUT $API_URL/queue/1 \
  -H "Content-Type: application/json" \
  -d '{"status": "paused"}'
```

### Resume Queue Item
```bash
curl -X PUT $API_URL/queue/1 \
  -H "Content-Type: application/json" \
  -d '{"status": "queued"}'
```

### Execute Queue Item Manually
```bash
curl -X POST $API_URL/queue/1/execute
```

### Delete Queue Item
```bash
curl -X DELETE $API_URL/queue/1
```

## Advanced Examples

### Complex Molecule with Custom Configuration
```bash
curl -X POST $API_URL/experiments/ \
  -H "Content-Type: application/json" \
  -d '{
    "name": "Benzene VQE Study",
    "molecule": {
      "smiles": "c1ccccc1",
      "basis": "6-31g",
      "charge": 0,
      "multiplicity": 1
    },
    "configuration": {
      "method": "VQE",
      "ansatz": "hardware_efficient",
      "mapper": "bravyi_kitaev",
      "optimizer": "L-BFGS-B",
      "backend": "classical",
      "max_iterations": 2000,
      "conv_threshold": 1e-8
    },
    "execute_immediately": false
  }'
```

### Charged Molecule
```bash
curl -X POST $API_URL/experiments/ \
  -H "Content-Type: application/json" \
  -d '{
    "name": "Ammonium Ion",
    "molecule": {
      "smiles": "[NH4+]",
      "basis": "sto-3g",
      "charge": 1,
      "multiplicity": 1
    },
    "configuration": {
      "method": "VQE",
      "ansatz": "ucc",
      "mapper": "jordan_wigner",
      "optimizer": "SLSQP",
      "backend": "classical"
    },
    "execute_immediately": true
  }'
```

### Triplet State Molecule
```bash
curl -X POST $API_URL/experiments/ \
  -H "Content-Type: application/json" \
  -d '{
    "name": "O2 Triplet State",
    "molecule": {
      "smiles": "O=O",
      "basis": "sto-3g",
      "charge": 0,
      "multiplicity": 3
    },
    "configuration": {
      "method": "VQE",
      "ansatz": "ucc",
      "mapper": "jordan_wigner",
      "optimizer": "SLSQP",
      "backend": "classical"
    },
    "execute_immediately": true
  }'
```

## Polling for Results

### Bash Script to Wait for Completion
```bash
#!/bin/bash
EXPERIMENT_ID=1

while true; do
    STATUS=$(curl -s $API_URL/experiments/$EXPERIMENT_ID/status | jq -r '.status')
    PROGRESS=$(curl -s $API_URL/experiments/$EXPERIMENT_ID/status | jq -r '.progress')

    echo "Status: $STATUS, Progress: $PROGRESS%"

    if [ "$STATUS" = "completed" ] || [ "$STATUS" = "failed" ]; then
        echo "Experiment finished!"
        curl $API_URL/experiments/$EXPERIMENT_ID | jq
        break
    fi

    sleep 2
done
```

## Batch Operations

### Create Multiple Experiments
```bash
# H2
curl -X POST $API_URL/experiments/ -H "Content-Type: application/json" -d '{"name": "H2 Test", "molecule": {"smiles": "[H][H]", "basis": "sto-3g"}, "configuration": {"method": "VQE"}, "execute_immediately": true}'

# H2O
curl -X POST $API_URL/experiments/ -H "Content-Type: application/json" -d '{"name": "H2O Test", "molecule": {"smiles": "O", "basis": "sto-3g"}, "configuration": {"method": "VQE"}, "execute_immediately": true}'

# CH4
curl -X POST $API_URL/experiments/ -H "Content-Type: application/json" -d '{"name": "CH4 Test", "molecule": {"smiles": "C", "basis": "sto-3g"}, "configuration": {"method": "VQE"}, "execute_immediately": true}'
```

## Error Cases

### Invalid SMILES
```bash
curl -X POST $API_URL/experiments/ \
  -H "Content-Type: application/json" \
  -d '{
    "name": "Invalid Molecule",
    "molecule": {
      "smiles": "INVALID",
      "basis": "sto-3g"
    },
    "configuration": {
      "method": "VQE"
    }
  }'
```

### Invalid Configuration
```bash
curl -X POST $API_URL/experiments/ \
  -H "Content-Type: application/json" \
  -d '{
    "name": "Bad Config",
    "molecule": {
      "smiles": "O",
      "basis": "invalid_basis"
    },
    "configuration": {
      "method": "UNKNOWN_METHOD"
    }
  }'
```

## Testing Authentication (Future)

When authentication is implemented:

```bash
# Login
TOKEN=$(curl -X POST $API_URL/auth/login \
  -H "Content-Type: application/json" \
  -d '{"username": "user", "password": "pass"}' | jq -r '.access_token')

# Use token
curl $API_URL/experiments/ \
  -H "Authorization: Bearer $TOKEN"
```

## Tips

1. **Pretty Print JSON**: Pipe to `jq`
   ```bash
   curl $API_URL/experiments/ | jq
   ```

2. **Save Response**: Use `-o` flag
   ```bash
   curl $API_URL/experiments/1 -o experiment.json
   ```

3. **Show Headers**: Use `-v` flag
   ```bash
   curl -v $API_URL/experiments/
   ```

4. **Follow Redirects**: Use `-L` flag
   ```bash
   curl -L $API_URL/experiments/
   ```

5. **Timeout**: Use `--max-time` flag
   ```bash
   curl --max-time 30 $API_URL/experiments/1
   ```
