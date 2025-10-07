#!/bin/bash
#
# Run Failing Tests - Isolated Execution
# =====================================
#
# Runs only the 22 known failing tests to verify fixes
#

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

echo -e "${BLUE}╔════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║           Running Failing Tests (22 tests)                    ║${NC}"
echo -e "${BLUE}╚════════════════════════════════════════════════════════════════╝${NC}"
echo ""

# Track results
TOTAL=0
PASSED=0
FAILED=0

# Function to run a test
run_test() {
    local test_path=$1
    local test_name=$2

    ((TOTAL++))

    echo -e "${YELLOW}Testing: ${test_name}${NC}"

    if python3 -m pytest "$test_path" -v --tb=short 2>&1 | grep -q "PASSED"; then
        echo -e "${GREEN}✓ PASSED${NC}"
        ((PASSED++))
    else
        echo -e "${RED}✗ FAILED${NC}"
        ((FAILED++))
    fi
    echo ""
}

echo -e "${BLUE}═══════════════════════════════════════════════════════════════${NC}"
echo -e "${BLUE}VQE Solver Tests (12 tests)${NC}"
echo -e "${BLUE}═══════════════════════════════════════════════════════════════${NC}"
echo ""

run_test "tests/unit/test_vqe.py::TestVQESolver::test_vqe_solver_creation" "VQE Solver Creation"
run_test "tests/unit/test_vqe.py::TestVQESolver::test_vqe_energy_computation" "VQE Energy Computation"
run_test "tests/unit/test_vqe.py::TestVQESolver::test_vqe_optimization" "VQE Optimization"
run_test "tests/unit/test_vqe.py::TestVQESolver::test_vqe_with_ucc_ansatz" "VQE with UCC Ansatz"
run_test "tests/unit/test_vqe.py::TestVQESolver::test_vqe_callback" "VQE Callback"
run_test "tests/unit/test_vqe.py::TestVQESolver::test_vqe_energy_history" "VQE Energy History"
run_test "tests/unit/test_vqe.py::TestVQESolver::test_vqe_variance_computation" "VQE Variance"
run_test "tests/unit/test_vqe.py::TestVQEIntegration::test_vqe_with_different_optimizers" "VQE Different Optimizers"
run_test "tests/unit/test_vqe.py::TestVQEIntegration::test_vqe_parameter_initialization" "VQE Param Init"
run_test "tests/unit/test_vqe.py::TestVQEIntegration::test_vqe_reproducibility" "VQE Reproducibility"
run_test "tests/unit/test_vqe.py::TestVQEIntegration::test_vqe_converges_below_initial" "VQE Convergence"

echo -e "${BLUE}═══════════════════════════════════════════════════════════════${NC}"
echo -e "${BLUE}Qiskit Integration Tests (7 tests)${NC}"
echo -e "${BLUE}═══════════════════════════════════════════════════════════════${NC}"
echo ""

run_test "tests/unit/test_qiskit_integration.py::TestQiskitBackend::test_aer_simulator_init" "Aer Simulator Init"
run_test "tests/unit/test_qiskit_integration.py::TestQiskitBackend::test_aer_statevector_init" "Aer Statevector Init"
run_test "tests/unit/test_qiskit_integration.py::TestQiskitBackend::test_get_estimator" "Get Estimator"
run_test "tests/unit/test_qiskit_integration.py::TestQiskitBackend::test_backend_info" "Backend Info"
run_test "tests/unit/test_qiskit_integration.py::TestQiskitVQE::test_vqe_classical_backend" "VQE Classical Backend"
run_test "tests/unit/test_qiskit_integration.py::TestQiskitVQE::test_vqe_aer_backend" "VQE Aer Backend"
run_test "tests/unit/test_qiskit_integration.py::TestQiskitVQE::test_vqe_classical_vs_qiskit" "VQE Classical vs Qiskit"

echo -e "${BLUE}═══════════════════════════════════════════════════════════════${NC}"
echo -e "${BLUE}Import Tests (4 tests)${NC}"
echo -e "${BLUE}═══════════════════════════════════════════════════════════════${NC}"
echo ""

run_test "tests/unit/test_imports.py::test_main_import" "Main Import"
run_test "tests/unit/test_imports.py::test_solvers_import" "Solvers Import"
run_test "tests/unit/test_imports.py::test_backends_import" "Backends Import"
run_test "tests/unit/test_imports.py::test_ibm_solvers_import" "IBM Solvers Import"

# Summary
echo ""
echo -e "${BLUE}╔════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║                        Summary                                 ║${NC}"
echo -e "${BLUE}╚════════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo -e "${GREEN}Passed: $PASSED/$TOTAL${NC}"
echo -e "${RED}Failed: $FAILED/$TOTAL${NC}"
echo ""

if [ $FAILED -eq 0 ]; then
    echo -e "${GREEN}╔════════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${GREEN}║              All failing tests are now fixed! ✓                ║${NC}"
    echo -e "${GREEN}╚════════════════════════════════════════════════════════════════╝${NC}"
    exit 0
else
    echo -e "${RED}╔════════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${RED}║              $FAILED tests still failing                            ║${NC}"
    echo -e "${RED}╚════════════════════════════════════════════════════════════════╝${NC}"
    exit 1
fi
