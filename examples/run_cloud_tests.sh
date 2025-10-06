#!/bin/bash
#
# Kanad Cloud Tests Runner
# ========================
#
# Convenient script to run cloud quantum chemistry tests
#

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}╔════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║           Kanad Cloud Quantum Chemistry Tests                 ║${NC}"
echo -e "${BLUE}╚════════════════════════════════════════════════════════════════╝${NC}"
echo ""

# Check if we're in the right directory
if [ ! -f "examples/cloud_vqe_drug_molecules.py" ]; then
    echo -e "${RED}Error: Must run from kanad root directory${NC}"
    exit 1
fi

# Function to check environment variable
check_env() {
    local var_name=$1
    if [ -z "${!var_name}" ]; then
        echo -e "${RED}✗ $var_name not set${NC}"
        return 1
    else
        echo -e "${GREEN}✓ $var_name is set${NC}"
        return 0
    fi
}

# Function to run a test
run_test() {
    local test_name=$1
    local test_script=$2
    local backend=$3

    echo ""
    echo -e "${YELLOW}═══════════════════════════════════════════════════════════════${NC}"
    echo -e "${YELLOW}Running: $test_name${NC}"
    echo -e "${YELLOW}Backend: $backend${NC}"
    echo -e "${YELLOW}═══════════════════════════════════════════════════════════════${NC}"
    echo ""

    if python3 "$test_script" --backend "$backend"; then
        echo -e "${GREEN}✓ $test_name completed successfully${NC}"
        return 0
    else
        echo -e "${RED}✗ $test_name failed${NC}"
        return 1
    fi
}

# Parse command line arguments
BACKEND="bluequbit"
TEST="all"

while [[ $# -gt 0 ]]; do
    case $1 in
        --backend)
            BACKEND="$2"
            shift 2
            ;;
        --test)
            TEST="$2"
            shift 2
            ;;
        --help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --backend <ibm|bluequbit>    Choose backend (default: bluequbit)"
            echo "  --test <test_name>           Run specific test (default: all)"
            echo ""
            echo "Available tests:"
            echo "  drug_molecules    - Drug molecule VQE"
            echo "  materials         - Materials science"
            echo "  custom_workflow   - Custom workflow comparison"
            echo "  all              - Run all tests"
            echo ""
            echo "Examples:"
            echo "  $0 --backend bluequbit --test drug_molecules"
            echo "  $0 --backend ibm --test materials"
            echo "  $0 --test all"
            exit 0
            ;;
        *)
            echo -e "${RED}Unknown option: $1${NC}"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Check environment variables based on backend
echo ""
echo -e "${BLUE}Checking environment...${NC}"
echo ""

if [ "$BACKEND" == "bluequbit" ]; then
    if ! check_env "BLUE_TOKEN"; then
        echo ""
        echo -e "${YELLOW}To set BlueQubit token:${NC}"
        echo "  export BLUE_TOKEN=your_token_here"
        echo ""
        echo "Get your token from: https://app.bluequbit.io"
        exit 1
    fi
elif [ "$BACKEND" == "ibm" ]; then
    if ! check_env "IBM_API"; then
        echo ""
        echo -e "${YELLOW}To set IBM Quantum token:${NC}"
        echo "  export IBM_API=your_token_here"
        echo ""
        echo "Get your token from: https://quantum.ibm.com"
        exit 1
    fi
else
    echo -e "${RED}Invalid backend: $BACKEND${NC}"
    echo "Use 'ibm' or 'bluequbit'"
    exit 1
fi

# Check Python dependencies
echo ""
echo -e "${BLUE}Checking dependencies...${NC}"
echo ""

MISSING_DEPS=()

if ! python3 -c "import kanad" 2>/dev/null; then
    echo -e "${RED}✗ kanad not installed${NC}"
    MISSING_DEPS+=("kanad")
else
    echo -e "${GREEN}✓ kanad installed${NC}"
fi

if ! python3 -c "import rdkit" 2>/dev/null; then
    echo -e "${YELLOW}⚠ rdkit not installed (optional, needed for SMILES)${NC}"
    MISSING_DEPS+=("rdkit")
else
    echo -e "${GREEN}✓ rdkit installed${NC}"
fi

if [ ${#MISSING_DEPS[@]} -gt 0 ]; then
    echo ""
    echo -e "${YELLOW}To install missing dependencies:${NC}"
    echo "  pip install -e ."
    if [[ " ${MISSING_DEPS[@]} " =~ " rdkit " ]]; then
        echo "  pip install rdkit"
    fi
    exit 1
fi

# Create results directory
mkdir -p results
echo ""
echo -e "${GREEN}✓ Results will be saved to: results/${NC}"

# Run tests
SUCCESS_COUNT=0
FAIL_COUNT=0

echo ""
echo -e "${BLUE}╔════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║                    Running Tests                               ║${NC}"
echo -e "${BLUE}╚════════════════════════════════════════════════════════════════╝${NC}"

if [ "$TEST" == "all" ] || [ "$TEST" == "drug_molecules" ]; then
    if run_test "Drug Molecules VQE" "examples/cloud_vqe_drug_molecules.py" "$BACKEND"; then
        ((SUCCESS_COUNT++))
    else
        ((FAIL_COUNT++))
    fi
fi

if [ "$TEST" == "all" ] || [ "$TEST" == "materials" ]; then
    if run_test "Materials Science" "examples/cloud_materials_science.py" "$BACKEND"; then
        ((SUCCESS_COUNT++))
    else
        ((FAIL_COUNT++))
    fi
fi

if [ "$TEST" == "all" ] || [ "$TEST" == "custom_workflow" ]; then
    echo ""
    echo -e "${YELLOW}═══════════════════════════════════════════════════════════════${NC}"
    echo -e "${YELLOW}Running: Custom Workflow (uses all available backends)${NC}"
    echo -e "${YELLOW}═══════════════════════════════════════════════════════════════${NC}"
    echo ""

    if python3 "examples/cloud_custom_workflow.py"; then
        echo -e "${GREEN}✓ Custom Workflow completed successfully${NC}"
        ((SUCCESS_COUNT++))
    else
        echo -e "${RED}✗ Custom Workflow failed${NC}"
        ((FAIL_COUNT++))
    fi
fi

# Summary
echo ""
echo -e "${BLUE}╔════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║                        Summary                                 ║${NC}"
echo -e "${BLUE}╚═════════════════════════════════════════════════════════���══════╝${NC}"
echo ""
echo -e "${GREEN}Successful: $SUCCESS_COUNT${NC}"
echo -e "${RED}Failed: $FAIL_COUNT${NC}"
echo ""

if [ $FAIL_COUNT -eq 0 ]; then
    echo -e "${GREEN}╔════════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${GREEN}║              All tests completed successfully! ✓               ║${NC}"
    echo -e "${GREEN}╚════════════════════════════════════════════════════════════════╝${NC}"
    exit 0
else
    echo -e "${RED}╔════════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${RED}║              Some tests failed. Check logs above.              ║${NC}"
    echo -e "${RED}╚════════════════════════════════════════════════════════════════╝${NC}"
    exit 1
fi
