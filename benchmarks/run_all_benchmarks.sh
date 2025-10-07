#!/bin/bash
# Run All Benchmarks in Parallel
# This script runs all available benchmark suites simultaneously

echo "========================================================================"
echo "RUNNING ALL QUANTUM CHEMISTRY BENCHMARKS IN PARALLEL"
echo "========================================================================"
echo ""
echo "Benchmarks to run:"
echo "  1. Kanad (VQE)"
echo "  2. PySCF (Classical: HF, MP2, FCI)"
echo "  3. DFT (Classical: LDA, PBE, B3LYP)"
echo "  4. PennyLane QChem"
echo "  5. OpenFermion"
echo "  6. Qiskit Nature (requires Python 3.11 venv)"
echo ""

# Create results directory
mkdir -p benchmarks/results

# Timestamp
START_TIME=$(date +%s)
echo "Start time: $(date)"
echo ""

# Run benchmarks in parallel (background jobs)
echo "Starting benchmarks..."

# 1. Kanad (always available)
echo "  [1/6] Kanad VQE..."
python3 benchmarks/kanad_benchmarks.py > benchmarks/results/kanad_run.log 2>&1 &
PID_KANAD=$!

# 2. PySCF Classical
echo "  [2/6] PySCF Classical..."
python3 benchmarks/pyscf_benchmarks.py > benchmarks/results/pyscf_run.log 2>&1 &
PID_PYSCF=$!

# 3. DFT
echo "  [3/6] DFT (LDA, PBE, B3LYP)..."
python3 benchmarks/dft_benchmarks.py > benchmarks/results/dft_run.log 2>&1 &
PID_DFT=$!

# 4. PennyLane (if installed)
echo "  [4/6] PennyLane QChem..."
python3 benchmarks/pennylane_benchmarks.py > benchmarks/results/pennylane_run.log 2>&1 &
PID_PENNYLANE=$!

# 5. OpenFermion (if installed)
echo "  [5/6] OpenFermion..."
python3 benchmarks/openfermion_benchmarks.py > benchmarks/results/openfermion_run.log 2>&1 &
PID_OPENFERMION=$!

# 6. Qiskit Nature (requires Python 3.11 venv - skip if not available)
echo "  [6/6] Qiskit Nature (checking for venv)..."
if [ -d "venv_qiskit_nature" ]; then
    echo "    Found venv_qiskit_nature, activating..."
    source venv_qiskit_nature/bin/activate && \
    python benchmarks/qiskit_nature_benchmarks.py > benchmarks/results/qiskit_nature_run.log 2>&1 &
    PID_QISKIT=$!
    deactivate 2>/dev/null || true
else
    echo "    Skipping (venv_qiskit_nature not found)"
    PID_QISKIT=""
fi

echo ""
echo "All benchmarks launched in parallel!"
echo "Waiting for completion..."
echo ""

# Wait for all jobs with progress indication
wait_for_job() {
    local pid=$1
    local name=$2

    if [ ! -z "$pid" ]; then
        while kill -0 $pid 2>/dev/null; do
            echo -n "."
            sleep 5
        done
        wait $pid
        status=$?
        if [ $status -eq 0 ]; then
            echo " ✓ $name complete"
        else
            echo " ✗ $name failed (exit code: $status)"
        fi
    fi
}

echo -n "  Kanad VQE        "
wait_for_job $PID_KANAD "Kanad"

echo -n "  PySCF Classical  "
wait_for_job $PID_PYSCF "PySCF"

echo -n "  DFT              "
wait_for_job $PID_DFT "DFT"

echo -n "  PennyLane        "
wait_for_job $PID_PENNYLANE "PennyLane"

echo -n "  OpenFermion      "
wait_for_job $PID_OPENFERMION "OpenFermion"

if [ ! -z "$PID_QISKIT" ]; then
    echo -n "  Qiskit Nature    "
    wait_for_job $PID_QISKIT "Qiskit Nature"
fi

# Calculate elapsed time
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
MINUTES=$((ELAPSED / 60))
SECONDS=$((ELAPSED % 60))

echo ""
echo "========================================================================"
echo "✓ ALL BENCHMARKS COMPLETE"
echo "========================================================================"
echo ""
echo "Total time: ${MINUTES}m ${SECONDS}s"
echo ""

# Count results
RESULT_COUNT=$(ls benchmarks/results/*_results.json 2>/dev/null | wc -l | tr -d ' ')
echo "Results generated: $RESULT_COUNT JSON files"
echo ""

# List results
echo "Results files:"
ls -lh benchmarks/results/*_results.json 2>/dev/null | awk '{print "  " $9 " (" $5 ")"}'

echo ""
echo "Log files:"
ls -lh benchmarks/results/*_run.log 2>/dev/null | awk '{print "  " $9 " (" $5 ")"}'

echo ""
echo "Next step: Generate comprehensive comparison report"
echo "  python3 benchmarks/generate_comprehensive_report.py"
echo ""
