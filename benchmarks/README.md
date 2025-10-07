# Kanad Benchmarking Suite

Phase 2: Framework benchmarking against industry standards

## Setup

### Main Kanad Environment
```bash
# Already set up with current packages
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

### Qiskit Nature Environment (Separate - uses deprecated packages)
```bash
# Create separate venv to avoid conflicts
python3 -m venv venv_qiskit_nature
source venv_qiskit_nature/bin/activate

# Install qiskit-nature with deprecated dependencies
pip install qiskit==0.45.3
pip install qiskit-nature==0.7.2
pip install pyscf==2.3.0
```

**Why separate?**
- qiskit-nature 0.7.2 requires older qiskit versions
- Conflicts with Kanad's qiskit-ibm-runtime requirements
- Keeps benchmarking isolated from production code

## Benchmark Categories

### 1. Energy Accuracy
Compare ground state energies against:
- Qiskit Nature VQE
- PySCF (classical reference)
- Published experimental values

**Molecules**: H2, LiH, H2O, N2

### 2. Performance & Speed
- Time to convergence
- Number of iterations
- Memory usage
- Scalability (qubits vs runtime)

### 3. Feature Completeness
- Supported basis sets
- Ansatz variety
- Mapper options
- Backend compatibility

### 4. Ease of Use
- Lines of code comparison
- API clarity
- Documentation quality
- Error messages

## Running Benchmarks

```bash
# Activate main Kanad environment
source venv/bin/activate

# Run Kanad benchmarks
python benchmarks/kanad_benchmarks.py

# Switch to qiskit-nature environment
deactivate
source venv_qiskit_nature/bin/activate

# Run qiskit-nature benchmarks
python benchmarks/qiskit_nature_benchmarks.py

# Compare results
python benchmarks/compare_results.py
```

## Expected Outputs

- `benchmarks/results/kanad_results.json`
- `benchmarks/results/qiskit_nature_results.json`
- `benchmarks/results/comparison_report.md`
- `benchmarks/results/performance_plots.png`

## Success Criteria

✅ Energy accuracy within 1 mHa (milliHartree) of qiskit-nature
✅ Competitive or faster convergence
✅ Cleaner API (fewer lines of code)
✅ Better error handling and validation
✅ Unique features: Governance protocols, bond-centric API
