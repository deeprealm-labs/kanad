# Qiskit Integration Examples

This directory contains examples demonstrating how to use the Kanad framework with Qiskit 2.x for quantum hardware execution.

## Prerequisites

```bash
# Install Qiskit 2.x
pip install qiskit>=2.2.0 qiskit-aer>=0.15.0

# For IBM Quantum hardware access
pip install qiskit-ibm-runtime>=0.30.0
```

## Examples

### 01_h2_vqe_comparison.py

Compares VQE performance across different backends:
- **Classical**: NumPy statevector (exact, fast)
- **Aer Statevector**: Qiskit Aer statevector simulator (exact)
- **Aer QASM**: Qiskit Aer QASM simulator (shot-based, with different shot counts)

**Run:**
```bash
python 01_h2_vqe_comparison.py
```

**Expected output:**
- Energy comparison table
- Accuracy analysis (error vs classical)
- Execution time comparison

### 02_hardware_execution_template.py

Template for running VQE on real IBM Quantum hardware.

**Setup:**
1. Create IBM Quantum account: https://quantum.ibm.com/
2. Get your API token
3. Save credentials:
```python
from qiskit_ibm_runtime import QiskitRuntimeService
QiskitRuntimeService.save_account(channel='ibm_quantum', token='YOUR_TOKEN')
```

**Run (simulation mode):**
```bash
python 02_hardware_execution_template.py
```

**Run on real hardware:**
```python
# Edit script to set use_real_hardware=True
result = run_on_ibm_hardware(backend_name='ibm_brisbane', use_real_hardware=True)
```

## Supported Backends

### Local Simulators
- `'classical'`: Custom NumPy simulator (exact, no Qiskit)
- `'aer_simulator_statevector'`: Aer statevector simulator (exact)
- `'aer_simulator'`: Aer QASM simulator (shot-based)

### IBM Quantum Hardware
- `'ibm_brisbane'`: 127-qubit Eagle processor
- `'ibm_kyoto'`: 127-qubit Eagle processor
- `'ibm_osaka'`: 127-qubit Eagle processor
- Check https://quantum.ibm.com/services/resources for available backends

## Usage Pattern

```python
from kanad.solvers.vqe_solver import VQESolver
from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper
from kanad.ansatze.hardware_efficient_ansatz import RealAmplitudesAnsatz

# Create VQE solver with Qiskit backend
solver = VQESolver(
    hamiltonian=your_hamiltonian,
    ansatz=RealAmplitudesAnsatz(n_qubits=4, n_electrons=2),
    mapper=JordanWignerMapper(),
    backend='aer_simulator',  # or 'ibm_brisbane' for hardware
    shots=1024,
    optimization_level=2
)

# Run VQE
result = solver.solve()
print(f"Energy: {result['energy']:.8f} Ha")
```

## Backend Options

### All Backends
- `backend`: Backend name (str)
- `shots`: Number of measurement shots (int, default=1024)
- `optimization_level`: Transpiler optimization level 0-3 (int, default=1)

### IBM Runtime (Hardware)
- `resilience_level`: Error mitigation level 0-2 (int, default=1)
  - 0: No mitigation
  - 1: Basic readout error mitigation
  - 2: Advanced mitigation (slower, more accurate)

Example:
```python
solver = VQESolver(
    ...,
    backend='ibm_brisbane',
    shots=8192,
    optimization_level=3,
    resilience_level=1  # Enable error mitigation
)
```

## Tips for Hardware Execution

1. **Start with simulation**: Always test with `'aer_simulator'` first
2. **Use shallow circuits**: Limit ansatz layers to 1-2 for NISQ hardware
3. **Increase shots**: Use 4096-8192 shots for better statistics
4. **Enable error mitigation**: Set `resilience_level=1` or `2`
5. **Monitor job queue**: Check https://quantum.ibm.com/jobs
6. **Compare with classical**: Validate results against NumPy simulation

## Troubleshooting

**Import Error: "No module named 'qiskit'"**
```bash
pip install qiskit>=2.2.0
```

**IBM Authentication Error**
```python
# Re-save credentials
from qiskit_ibm_runtime import QiskitRuntimeService
QiskitRuntimeService.save_account(channel='ibm_quantum', token='YOUR_TOKEN', overwrite=True)
```

**Large energy error on hardware**
- Increase `shots` to 8192 or more
- Enable error mitigation: `resilience_level=1`
- Reduce circuit depth: decrease `n_layers` in ansatz
- Try different backend with lower error rates

## Additional Resources

- [Qiskit Documentation](https://docs.quantum.ibm.com/)
- [IBM Quantum Experience](https://quantum.ibm.com/)
- [Kanad Documentation](../../README.md)
