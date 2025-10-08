# Kanad Backend API - Usage Examples & Integration Guide

Complete guide for frontend developers integrating with the Kanad quantum chemistry backend.

---

## Quick Start

### 1. Start the API Server

```bash
cd /home/mk/deeprealm/kanad
python -m api.main
```

Server runs at: `http://localhost:8000`

### 2. Check API Health

```bash
curl http://localhost:8000/health
```

Expected response:
```json
{
  "status": "healthy",
  "job_queue_running": true,
  "queue_size": 0
}
```

### 3. View API Capabilities

```bash
curl http://localhost:8000/api/v1/info
```

---

## Example 1: Simple VQE Calculation (H2 Molecule)

### Frontend Code (JavaScript/React)

```javascript
const runVQECalculation = async () => {
  const experimentData = {
    name: "H2 VQE Ground State",
    molecule: {
      smiles: "[H][H]",
      basis: "sto-3g",
      charge: 0,
      multiplicity: 1
    },
    configuration: {
      method: "VQE",
      ansatz: "ucc",
      mapper: "jordan_wigner",
      optimizer: "SLSQP",
      backend: "classical",
      max_iterations: 1000,
      conv_threshold: 1e-6
    },
    execute_immediately: true
  };

  // Create experiment
  const response = await fetch('http://localhost:8000/api/v1/experiments', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(experimentData)
  });

  const { id } = await response.json();
  console.log(`Experiment created: ID ${id}`);

  // Poll for convergence (real-time energy updates)
  const pollInterval = setInterval(async () => {
    const convResponse = await fetch(
      `http://localhost:8000/api/v1/experiments/${id}/convergence`
    );
    const { status, convergence_data } = await convResponse.json();

    // Update energy plot
    updateEnergyGraph(convergence_data);

    // Check if complete
    if (status === 'completed' || status === 'failed') {
      clearInterval(pollInterval);

      // Get final results
      const resultsResponse = await fetch(
        `http://localhost:8000/api/v1/experiments/${id}`
      );
      const results = await resultsResponse.json();

      displayResults(results);
    }
  }, 2000); // Poll every 2 seconds
};
```

### cURL Example

```bash
# Create experiment
curl -X POST http://localhost:8000/api/v1/experiments \
  -H "Content-Type: application/json" \
  -d '{
    "name": "H2 VQE",
    "molecule": {
      "smiles": "[H][H]",
      "basis": "sto-3g"
    },
    "configuration": {
      "method": "VQE",
      "ansatz": "ucc",
      "mapper": "jordan_wigner",
      "optimizer": "SLSQP",
      "backend": "classical",
      "max_iterations": 1000
    },
    "execute_immediately": true
  }'

# Response:
# {"id": 1, "name": "H2 VQE", "status": "queued", ...}

# Check status
curl http://localhost:8000/api/v1/experiments/1/status

# Get convergence data (while running)
curl http://localhost:8000/api/v1/experiments/1/convergence

# Get full results (when complete)
curl http://localhost:8000/api/v1/experiments/1
```

---

## Example 2: Circuit Visualization

### Get Circuit in Multiple Formats

```javascript
const getCircuitVisualization = async (experimentId, format = 'json') => {
  const response = await fetch(
    `http://localhost:8000/api/v1/experiments/${experimentId}/circuit?format=${format}`
  );
  return await response.json();
};

// Usage
const jsonCircuit = await getCircuitVisualization(1, 'json');
const asciiCircuit = await getCircuitVisualization(1, 'ascii');
const qasmCircuit = await getCircuitVisualization(1, 'qasm');

console.log(jsonCircuit.gates);
console.log(asciiCircuit.ascii);
console.log(qasmCircuit.qasm);
```

### JSON Format Response

```json
{
  "n_qubits": 4,
  "n_electrons": 2,
  "ansatz_type": "ucc",
  "method": "VQE",
  "gates": [
    {"type": "x", "qubits": [0], "parameterized": false},
    {"type": "x", "qubits": [2], "parameterized": false},
    {"type": "ry", "qubits": [0], "parameter": "theta_0", "parameterized": true},
    {"type": "cx", "qubits": [0, 2], "parameterized": false}
  ],
  "metadata": {
    "gate_count": 24,
    "parameter_count": 8
  }
}
```

### Rendering Circuit (React Component Example)

```jsx
import React from 'react';

const CircuitVisualization = ({ gates, nQubits }) => {
  return (
    <div className="circuit-container">
      <h3>Quantum Circuit</h3>
      <svg width="800" height={nQubits * 60}>
        {gates.map((gate, idx) => (
          <Gate key={idx} gate={gate} position={idx} />
        ))}
      </svg>
    </div>
  );
};

const Gate = ({ gate, position }) => {
  const x = 100 + position * 50;
  const y = gate.qubits[0] * 60 + 30;

  if (gate.type === 'cx') {
    // Render CNOT gate
    return (
      <g>
        <circle cx={x} cy={y} r="10" fill="blue" />
        <circle cx={x} cy={gate.qubits[1] * 60 + 30} r="15"
                fill="none" stroke="blue" strokeWidth="2" />
        <line x1={x} y1={y} x2={x} y2={gate.qubits[1] * 60 + 30}
              stroke="blue" strokeWidth="2" />
      </g>
    );
  }

  return (
    <rect x={x - 15} y={y - 15} width="30" height="30"
          fill={gate.parameterized ? "orange" : "green"} />
  );
};
```

---

## Example 3: Experiment Report Generation

### Get Report in JSON Format

```javascript
const getExperimentReport = async (experimentId, format = 'json') => {
  const response = await fetch(
    `http://localhost:8000/api/v1/experiments/${experimentId}/report?format=${format}`
  );
  return await response.json();
};

// Usage
const report = await getExperimentReport(1, 'json');

console.log('Molecule:', report.molecule);
console.log('Method:', report.method);
console.log('Results:', report.results);
console.log('Energy:', report.results.energy, 'Ha');
console.log('Correlation:', report.results.correlation_energy, 'Ha');
```

### Report Response Structure

```json
{
  "experiment_id": 1,
  "name": "H2 VQE",
  "created_at": "2025-10-08T10:30:00",
  "completed_at": "2025-10-08T10:32:15",
  "execution_time": 135.5,

  "molecule": {
    "smiles": "[H][H]",
    "formula": "H2",
    "n_electrons": 2,
    "n_orbitals": 2,
    "basis": "sto-3g",
    "charge": 0,
    "multiplicity": 1
  },

  "method": {
    "type": "VQE",
    "ansatz": "ucc",
    "mapper": "jordan_wigner",
    "optimizer": "SLSQP",
    "backend": "classical",
    "max_iterations": 1000
  },

  "results": {
    "energy": -1.137283,
    "hf_energy": -1.116685,
    "correlation_energy": -0.020598,
    "converged": true,
    "iterations": 127
  },

  "convergence": {
    "data": [
      {"iteration": 1, "energy": -1.1167},
      {"iteration": 2, "energy": -1.1245},
      ...
    ],
    "final_iteration": 127
  },

  "analysis": {
    "bonding": {...},
    "properties": {...}
  }
}
```

### Get Markdown Report

```javascript
const report = await getExperimentReport(1, 'markdown');
console.log(report.content);

// Save to file or display in frontend
downloadFile('experiment_report.md', report.content);
```

---

## Example 4: SQD Multi-State Calculation

### Frontend Code

```javascript
const runSQDCalculation = async (smiles, nStates = 3) => {
  const experimentData = {
    name: `SQD ${nStates} States`,
    molecule: {
      smiles: smiles,
      basis: "sto-3g"
    },
    configuration: {
      method: "SQD",
      n_states: nStates,
      subspace_dim: 15,
      backend: "classical"
    },
    execute_immediately: true
  };

  const response = await fetch('http://localhost:8000/api/v1/experiments', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(experimentData)
  });

  const { id } = await response.json();

  // Wait for completion
  await waitForCompletion(id);

  // Get results
  const resultsResponse = await fetch(
    `http://localhost:8000/api/v1/experiments/${id}`
  );
  const results = await resultsResponse.json();

  return {
    groundStateEnergy: results.results.ground_state_energy,
    excitedEnergies: results.results.excited_state_energies,
    allEnergies: results.results.energies
  };
};

// Usage
const { groundStateEnergy, excitedEnergies } = await runSQDCalculation('O', 3);
console.log('Ground State:', groundStateEnergy, 'Ha');
console.log('1st Excited:', excitedEnergies[0], 'Ha');
console.log('2nd Excited:', excitedEnergies[1], 'Ha');
```

### cURL Example

```bash
curl -X POST http://localhost:8000/api/v1/experiments \
  -H "Content-Type: application/json" \
  -d '{
    "name": "H2O SQD Multi-State",
    "molecule": {
      "smiles": "O",
      "basis": "sto-3g"
    },
    "configuration": {
      "method": "SQD",
      "n_states": 3,
      "subspace_dim": 15
    },
    "execute_immediately": true
  }'
```

---

## Example 5: Excited States (UV-Vis Spectroscopy)

### Frontend Code

```javascript
const calculateExcitedStates = async (smiles, nStates = 5) => {
  const experimentData = {
    name: "Excited States (CIS)",
    molecule: {
      smiles: smiles,
      basis: "6-31g"
    },
    configuration: {
      method: "EXCITED_STATES",
      excited_method: "cis",
      n_states: nStates
    },
    execute_immediately: true
  };

  const response = await fetch('http://localhost:8000/api/v1/experiments', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(experimentData)
  });

  const { id } = await response.json();
  await waitForCompletion(id);

  const resultsResponse = await fetch(
    `http://localhost:8000/api/v1/experiments/${id}`
  );
  const results = await resultsResponse.json();

  return {
    groundState: results.results.ground_state_energy,
    excitationEnergies: results.results.excitation_energies_ev,
    oscillatorStrengths: results.results.oscillator_strengths,
    transitions: results.results.dominant_transitions
  };
};

// Usage
const spectrum = await calculateExcitedStates('[H][H]', 5);

// Display UV-Vis spectrum
spectrum.excitationEnergies.forEach((energy, i) => {
  const wavelength = 1240 / energy; // nm
  const intensity = spectrum.oscillatorStrengths[i];
  const transition = spectrum.transitions[i];

  console.log(`λ=${wavelength.toFixed(1)}nm (${energy.toFixed(2)}eV): f=${intensity.toFixed(4)} [${transition}]`);
});
```

### Response Format

```json
{
  "results": {
    "method": "ExcitedStates-CIS",
    "ground_state_energy": -1.116685,
    "excited_state_energies": [-1.05, -0.98, -0.92, -0.85],
    "excitation_energies_ev": [1.82, 3.71, 5.35, 7.24],
    "oscillator_strengths": [0.0, 0.421, 0.012, 0.325],
    "dominant_transitions": [
      "HOMO-0 → LUMO+0",
      "HOMO-1 → LUMO+0",
      "HOMO-0 → LUMO+1",
      "HOMO-1 → LUMO+1"
    ]
  }
}
```

---

## Example 6: Queue Management

### Get Queue Statistics

```javascript
const getQueueStats = async () => {
  const response = await fetch('http://localhost:8000/api/v1/queue/stats');
  return await response.json();
};

// Usage
const stats = await getQueueStats();

console.log('Queue:', stats.queue);
console.log('Active jobs:', stats.queue.running);
console.log('Pending jobs:', stats.queue.queued);
console.log('Success rate:', stats.experiments.success_rate, '%');
```

### Response

```json
{
  "queue": {
    "total": 15,
    "queued": 3,
    "running": 1,
    "completed": 10,
    "failed": 1,
    "scheduled": 0,
    "active_size": 3
  },
  "experiments": {
    "total": 50,
    "completed": 42,
    "failed": 3,
    "success_rate": 93.33
  }
}
```

### Dashboard Widget (React)

```jsx
const QueueStatsWidget = () => {
  const [stats, setStats] = useState(null);

  useEffect(() => {
    const fetchStats = async () => {
      const response = await fetch('http://localhost:8000/api/v1/queue/stats');
      const data = await response.json();
      setStats(data);
    };

    fetchStats();
    const interval = setInterval(fetchStats, 5000); // Update every 5s

    return () => clearInterval(interval);
  }, []);

  if (!stats) return <div>Loading...</div>;

  return (
    <div className="queue-stats">
      <h3>Job Queue Statistics</h3>
      <div className="stat-grid">
        <StatCard label="Running" value={stats.queue.running} color="blue" />
        <StatCard label="Queued" value={stats.queue.queued} color="yellow" />
        <StatCard label="Completed" value={stats.queue.completed} color="green" />
        <StatCard label="Failed" value={stats.queue.failed} color="red" />
      </div>
      <div className="success-rate">
        Success Rate: {stats.experiments.success_rate}%
      </div>
    </div>
  );
};
```

---

## Example 7: Real-time Convergence Monitoring

### Live Energy Plot (React)

```jsx
import React, { useEffect, useState } from 'react';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip } from 'recharts';

const ConvergencePlot = ({ experimentId }) => {
  const [convergenceData, setConvergenceData] = useState([]);
  const [status, setStatus] = useState('running');

  useEffect(() => {
    const pollConvergence = setInterval(async () => {
      const response = await fetch(
        `http://localhost:8000/api/v1/experiments/${experimentId}/convergence`
      );
      const data = await response.json();

      setConvergenceData(data.convergence_data);
      setStatus(data.status);

      if (data.status === 'completed' || data.status === 'failed') {
        clearInterval(pollConvergence);
      }
    }, 2000);

    return () => clearInterval(pollConvergence);
  }, [experimentId]);

  return (
    <div className="convergence-plot">
      <h3>Energy Convergence - {status}</h3>
      <LineChart width={800} height={400} data={convergenceData}>
        <CartesianGrid strokeDasharray="3 3" />
        <XAxis dataKey="iteration" label="Iteration" />
        <YAxis label="Energy (Ha)" />
        <Tooltip />
        <Line
          type="monotone"
          dataKey="energy"
          stroke="#8884d8"
          strokeWidth={2}
          dot={false}
        />
      </LineChart>
      {status === 'completed' && (
        <div className="final-energy">
          Final Energy: {convergenceData[convergenceData.length - 1]?.energy.toFixed(6)} Ha
        </div>
      )}
    </div>
  );
};
```

---

## Example 8: Batch Experiments

### Submit Multiple Experiments

```javascript
const runBatchExperiments = async (molecules) => {
  const experimentIds = [];

  for (const mol of molecules) {
    const experimentData = {
      name: `VQE - ${mol.name}`,
      molecule: {
        smiles: mol.smiles,
        basis: "sto-3g"
      },
      configuration: {
        method: "VQE",
        ansatz: "ucc",
        mapper: "jordan_wigner",
        optimizer: "SLSQP",
        backend: "classical"
      },
      execute_immediately: false // Don't execute yet
    };

    const response = await fetch('http://localhost:8000/api/v1/experiments', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(experimentData)
    });

    const { id } = await response.json();
    experimentIds.push(id);

    // Add to queue with priority
    await fetch('http://localhost:8000/api/v1/queue', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        experiment_id: id,
        priority: mol.priority || 0
      })
    });
  }

  return experimentIds;
};

// Usage
const molecules = [
  { name: 'H2', smiles: '[H][H]', priority: 1 },
  { name: 'H2O', smiles: 'O', priority: 2 },
  { name: 'CH4', smiles: 'C', priority: 0 }
];

const ids = await runBatchExperiments(molecules);
console.log('Submitted experiments:', ids);
```

---

## Example 9: Custom Configuration Options

### Hardware-Efficient Ansatz

```javascript
const experimentData = {
  name: "H2O Hardware-Efficient VQE",
  molecule: {
    smiles: "O",
    basis: "sto-3g"
  },
  configuration: {
    method: "VQE",
    ansatz: "hardware_efficient",  // Fewer parameters than UCC
    mapper: "jordan_wigner",
    optimizer: "COBYLA",
    backend: "classical",
    max_iterations: 500
  },
  execute_immediately: true
};
```

### Governance-Aware Ansatz

```javascript
const experimentData = {
  name: "NaCl Ionic Bond VQE",
  molecule: {
    smiles: "[Na+].[Cl-]",
    basis: "sto-3g"
  },
  configuration: {
    method: "VQE",
    ansatz: "governance",  // Automatically adapts to ionic bonding
    mapper: "bravyi_kitaev",
    optimizer: "SLSQP",
    backend: "classical"
  },
  execute_immediately: true
};
```

### IBM Quantum Hardware

```javascript
const experimentData = {
  name: "H2 on IBM Brisbane",
  molecule: {
    smiles: "[H][H]",
    basis: "sto-3g"
  },
  configuration: {
    method: "VQE",
    ansatz: "hardware_efficient",
    mapper: "jordan_wigner",
    optimizer: "COBYLA",
    backend: "ibm_quantum",
    backend_name: "ibm_brisbane",
    shots: 8192
  },
  execute_immediately: true
};
```

---

## Example 10: Error Handling

### Robust Frontend Implementation

```javascript
const runExperimentWithErrorHandling = async (experimentData) => {
  try {
    // Create experiment
    const response = await fetch('http://localhost:8000/api/v1/experiments', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(experimentData)
    });

    if (!response.ok) {
      const error = await response.json();
      throw new Error(error.detail || 'Failed to create experiment');
    }

    const { id } = await response.json();

    // Monitor progress
    return new Promise((resolve, reject) => {
      const pollInterval = setInterval(async () => {
        try {
          const statusResponse = await fetch(
            `http://localhost:8000/api/v1/experiments/${id}/status`
          );
          const { status, error_message } = await statusResponse.json();

          if (status === 'completed') {
            clearInterval(pollInterval);
            const resultsResponse = await fetch(
              `http://localhost:8000/api/v1/experiments/${id}`
            );
            const results = await resultsResponse.json();
            resolve(results);
          } else if (status === 'failed') {
            clearInterval(pollInterval);
            reject(new Error(error_message || 'Experiment failed'));
          }
        } catch (err) {
          clearInterval(pollInterval);
          reject(err);
        }
      }, 2000);

      // Timeout after 10 minutes
      setTimeout(() => {
        clearInterval(pollInterval);
        reject(new Error('Experiment timeout'));
      }, 600000);
    });

  } catch (error) {
    console.error('Experiment error:', error);
    throw error;
  }
};

// Usage
try {
  const results = await runExperimentWithErrorHandling(experimentData);
  console.log('Success:', results);
} catch (error) {
  console.error('Failed:', error.message);
  // Show error to user
}
```

---

## API Response Time Expectations

| Endpoint | Typical Response Time | Notes |
|----------|----------------------|-------|
| GET /experiments | < 50ms | Pagination helps |
| GET /experiments/{id} | < 50ms | Cached if accessed frequently |
| POST /experiments | < 100ms | Creates record + queues job |
| GET /queue/stats | < 30ms | Aggregation query |
| GET /circuit | < 100ms | Generates circuit on-demand |
| GET /report | < 200ms | Compiles comprehensive report |
| GET /convergence | < 50ms | Returns stored data |

**Computation Times (varies by molecule size):**
- H2 (HF): < 1 second
- H2 (VQE): 10-60 seconds
- H2O (VQE): 30-180 seconds
- Benzene (VQE): Minutes to hours

---

## Best Practices

### 1. Polling Strategy

```javascript
// Good: Adaptive polling
let pollInterval = 2000; // Start with 2 seconds
const maxPollInterval = 10000; // Max 10 seconds

const poll = async () => {
  // ... fetch convergence data ...

  // If data hasn't changed much, slow down polling
  if (dataUnchanged) {
    pollInterval = Math.min(pollInterval * 1.5, maxPollInterval);
  }
};
```

### 2. Circuit Caching

```javascript
// Cache circuit data (doesn't change after experiment completes)
const circuitCache = new Map();

const getCircuit = async (experimentId) => {
  if (circuitCache.has(experimentId)) {
    return circuitCache.get(experimentId);
  }

  const response = await fetch(
    `http://localhost:8000/api/v1/experiments/${experimentId}/circuit?format=json`
  );
  const circuit = await response.json();
  circuitCache.set(experimentId, circuit);
  return circuit;
};
```

### 3. Experiment Status Updates

```javascript
// Use status endpoint for lightweight polling
const checkStatus = async (experimentId) => {
  const response = await fetch(
    `http://localhost:8000/api/v1/experiments/${experimentId}/status`
  );
  return await response.json();
};

// Only fetch full results when complete
if (status === 'completed') {
  const fullResults = await fetch(
    `http://localhost:8000/api/v1/experiments/${experimentId}`
  ).then(r => r.json());
}
```

---

## Troubleshooting

### Experiment Fails Immediately

**Check:**
1. SMILES string validity
2. Basis set spelling
3. Method/ansatz/mapper combination compatibility

**Debug:**
```javascript
const response = await fetch(`/api/v1/experiments/${id}`);
const data = await response.json();
console.log(data.error_message);
```

### Circuit Endpoint Returns 400

**Reason:** Experiment not completed yet or method doesn't use circuits (HF)

**Solution:**
```javascript
const status = await checkStatus(experimentId);
if (status.status !== 'completed') {
  console.log('Wait for experiment to complete');
}

const config = experiment.configuration;
if (config.method === 'HF') {
  console.log('HF method does not use quantum circuits');
}
```

### Convergence Data Empty

**Reason:** Experiment hasn't started yet or method doesn't track convergence

**Solution:** Only VQE and SQD provide convergence data

---

## Additional Resources

- **OpenAPI Docs:** http://localhost:8000/docs
- **ReDoc:** http://localhost:8000/redoc
- **Health Check:** http://localhost:8000/health
- **API Info:** http://localhost:8000/api/v1/info

---

**Last Updated:** 2025-10-08
