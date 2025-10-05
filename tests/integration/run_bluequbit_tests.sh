#!/bin/bash
# BlueQubit Integration Tests Runner

echo "================================================================================"
echo "                 BLUEQUBIT CLOUD BACKEND INTEGRATION TESTS"
echo "================================================================================"
echo ""

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Check environment
echo -e "${BLUE}[1/4] Checking Environment${NC}"
echo "────────────────────────────────────────────────────────────────────────────────"

# Activate virtual environment
if [ -f "env/bin/activate" ]; then
    source env/bin/activate
    echo -e "${GREEN}✓${NC} Virtual environment activated"
else
    echo -e "${RED}✗${NC} Virtual environment not found!"
    exit 1
fi

# Check Python version
PYTHON_VERSION=$(python --version 2>&1)
echo "  Python: $PYTHON_VERSION"

# Check if BlueQubit SDK is installed
if python -c "import bluequbit" 2>/dev/null; then
    BQ_VERSION=$(python -c "import bluequbit; print(bluequbit.__version__)" 2>/dev/null || echo "unknown")
    echo -e "${GREEN}✓${NC} BlueQubit SDK installed (version: $BQ_VERSION)"
else
    echo -e "${YELLOW}⚠${NC}  BlueQubit SDK not found. Installing..."
    pip install bluequbit -q
    echo -e "${GREEN}✓${NC} BlueQubit SDK installed"
fi

# Check for API token
echo ""
echo -e "${BLUE}[2/4] Checking API Token${NC}"
echo "────────────────────────────────────────────────────────────────────────────────"

if [ -f ".env" ]; then
    echo -e "${GREEN}✓${NC} .env file found"

    # Load token
    export $(grep -v '^#' .env | grep TOKEN | xargs)

    if [ -n "$TOKEN" ]; then
        echo -e "${GREEN}✓${NC} BlueQubit token loaded from .env"
        echo "  Token: ${TOKEN:0:10}... (${#TOKEN} chars)"
        export BLUEQUBIT_API_TOKEN=$TOKEN
    else
        echo -e "${YELLOW}⚠${NC}  No TOKEN found in .env file"
        echo ""
        echo "  To run cloud tests, add your BlueQubit token to .env:"
        echo "  TOKEN=<your-bluequbit-token>"
        echo ""
        echo "  Get your token at: https://app.bluequbit.io"
        echo ""
        echo "  Tests will be skipped if token is missing."
    fi
else
    echo -e "${YELLOW}⚠${NC}  .env file not found"
    echo "  Create .env with: TOKEN=<your-bluequbit-token>"
fi

# Run pytest
echo ""
echo -e "${BLUE}[3/4] Running Integration Tests${NC}"
echo "────────────────────────────────────────────────────────────────────────────────"
echo ""

# Set Python path
export PYTHONPATH="${PYTHONPATH}:$(pwd)"

# Run tests with different verbosity levels
if [ "$1" == "-v" ] || [ "$1" == "--verbose" ]; then
    # Verbose mode
    pytest tests/integration/test_bluequbit_backend.py -v -s \
        --tb=short \
        --color=yes \
        -m "not slow" \
        2>&1
elif [ "$1" == "-a" ] || [ "$1" == "--all" ]; then
    # Run all tests including slow ones
    echo -e "${YELLOW}Running ALL tests (including slow cloud execution tests)${NC}"
    echo ""
    pytest tests/integration/test_bluequbit_backend.py -v -s \
        --tb=short \
        --color=yes \
        2>&1
elif [ "$1" == "-f" ] || [ "$1" == "--fast" ]; then
    # Fast mode - skip slow tests
    echo -e "${YELLOW}Running FAST tests only (skip cloud execution)${NC}"
    echo ""
    pytest tests/integration/test_bluequbit_backend.py -v \
        --tb=line \
        --color=yes \
        -m "not slow" \
        2>&1
else
    # Default mode
    pytest tests/integration/test_bluequbit_backend.py \
        --tb=short \
        --color=yes \
        -m "not slow" \
        2>&1
fi

TEST_RESULT=$?

# Summary
echo ""
echo "================================================================================"
echo -e "${BLUE}[4/4] Test Summary${NC}"
echo "================================================================================"
echo ""

if [ $TEST_RESULT -eq 0 ]; then
    echo -e "${GREEN}✓ ALL TESTS PASSED${NC}"
    echo ""
    echo "BlueQubit integration is working correctly!"
    echo ""
elif [ $TEST_RESULT -eq 5 ]; then
    echo -e "${YELLOW}⚠ TESTS SKIPPED${NC}"
    echo ""
    echo "Most tests were skipped (likely due to missing API token)."
    echo ""
    echo "To run cloud tests:"
    echo "1. Get token from: https://app.bluequbit.io"
    echo "2. Add to .env: TOKEN=<your-token>"
    echo "3. Re-run tests: $0"
    echo ""
else
    echo -e "${RED}✗ SOME TESTS FAILED${NC}"
    echo ""
    echo "Check the output above for details."
    echo ""
fi

# Usage info
echo "Usage:"
echo "  $0           - Run standard tests (skip slow tests)"
echo "  $0 -v        - Verbose mode with detailed output"
echo "  $0 -a        - Run ALL tests including slow cloud execution"
echo "  $0 -f        - Fast mode (skip cloud execution tests)"
echo ""

echo "================================================================================"
exit $TEST_RESULT
