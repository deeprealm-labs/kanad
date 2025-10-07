# Qiskit Nature Benchmarking Setup Guide

## Issue: Python 3.13 Incompatibility

qiskit-aer 0.13.3 and pyscf 2.3.0 fail to build on Python 3.13 due to C++ compilation issues.

## Solution Options

### Option 1: Use Python 3.11 (Recommended)

```bash
# Install Python 3.11 if not already installed
# macOS: brew install python@3.11
# Linux: apt install python3.11 python3.11-venv

# Create venv with Python 3.11
python3.11 -m venv venv_qiskit_nature
source venv_qiskit_nature/bin/activate

# Install packages
pip install qiskit==0.45.3 qiskit-aer==0.13.3 qiskit-nature==0.7.2 pyscf==2.3.0
```

### Option 2: Use Docker Container

```bash
cd benchmarks
docker build -t kanad-benchmark -f Dockerfile.qiskit_nature .
docker run -it -v $(pwd)/results:/results kanad-benchmark
```

```dockerfile
# benchmarks/Dockerfile.qiskit_nature
FROM python:3.11-slim

WORKDIR /benchmark

# Install system dependencies
RUN apt-get update && apt-get install -y \
    gcc g++ cmake ninja-build \
    libopenblas-dev libhdf5-dev \
    && rm -rf /var/lib/apt/lists/*

# Install qiskit-nature stack
RUN pip install --no-cache-dir \
    qiskit==0.45.3 \
    qiskit-aer==0.13.3 \
    qiskit-nature==0.7.2 \
    pyscf==2.3.0 \
    numpy scipy pandas matplotlib

COPY . /benchmark

CMD ["python", "qiskit_nature_benchmarks.py"]
```

### Option 3: Use Pre-built Wheels (if available)

```bash
# For macOS ARM64 (M1/M2/M3)
pip install qiskit-aer --only-binary :all:

# If not available, use Rosetta emulation:
arch -x86_64 pip install qiskit-aer==0.13.3
```

### Option 4: Skip Qiskit-Aer (Use Statevector Only)

```bash
# Install without qiskit-aer
pip install qiskit==0.45.3 qiskit-nature==0.7.2 pyscf==2.3.0

# Use BasicSimulator or StatevectorSimulator (slower, but works)
```

## Current Status

✅ Python 3.13 - Works for Kanad
❌ Python 3.13 - Fails for qiskit-aer/pyscf compilation
✅ Python 3.11 - Works for all packages

## Recommendation for User

**Use Python 3.11 in separate venv for benchmarking**

1. Install Python 3.11:
   ```bash
   brew install python@3.11  # macOS
   ```

2. Create benchmarking venv:
   ```bash
   python3.11 -m venv venv_qiskit_nature
   source venv_qiskit_nature/bin/activate
   pip install qiskit==0.45.3 qiskit-aer==0.13.3 qiskit-nature==0.7.2 pyscf==2.3.0
   ```

3. Run benchmarks:
   ```bash
   python benchmarks/qiskit_nature_benchmarks.py
   ```

4. Return to main Kanad development:
   ```bash
   deactivate
   source venv/bin/activate  # Main Kanad venv with Python 3.13
   ```

## Alternative: Use Cloud Benchmarking Service

If local compilation continues to fail, consider:
- Google Colab (free, Python 3.10)
- Azure Notebooks
- AWS SageMaker Notebooks
- Run benchmarks on separate machine with Python 3.11
