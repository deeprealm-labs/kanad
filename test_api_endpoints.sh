#!/bin/bash

# Kanad API Endpoint Testing Script
# Tests all new and existing endpoints

BASE_URL="http://localhost:8000"
API_V1="${BASE_URL}/api/v1"

echo "========================================="
echo "Kanad API Endpoint Testing"
echo "========================================="
echo ""

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Test function
test_endpoint() {
    local method=$1
    local endpoint=$2
    local data=$3
    local description=$4

    echo -e "${YELLOW}Testing:${NC} $description"
    echo -e "${YELLOW}Endpoint:${NC} $method $endpoint"

    if [ "$method" == "GET" ]; then
        response=$(curl -s -w "\n%{http_code}" "$endpoint")
    elif [ "$method" == "POST" ]; then
        response=$(curl -s -w "\n%{http_code}" -X POST "$endpoint" \
            -H "Content-Type: application/json" \
            -d "$data")
    fi

    http_code=$(echo "$response" | tail -n 1)
    body=$(echo "$response" | sed '$d')

    if [ "$http_code" -ge 200 ] && [ "$http_code" -lt 300 ]; then
        echo -e "${GREEN}✓ SUCCESS${NC} (HTTP $http_code)"
    else
        echo -e "${RED}✗ FAILED${NC} (HTTP $http_code)"
        echo "Response: $body"
    fi
    echo ""
}

# 1. Health Check
test_endpoint "GET" "${BASE_URL}/health" "" "Health Check"

# 2. API Info
test_endpoint "GET" "${API_V1}/info" "" "API Information and Capabilities"

# 3. Queue Statistics (NEW)
test_endpoint "GET" "${API_V1}/queue/stats" "" "Queue Statistics"

# 4. List Experiments
test_endpoint "GET" "${API_V1}/experiments?limit=10" "" "List Experiments"

# 5. List Queue
test_endpoint "GET" "${API_V1}/queue" "" "List Queue Items"

# 6. Create Simple HF Experiment
experiment_data='{
  "name": "Test H2 HF",
  "molecule": {
    "smiles": "[H][H]",
    "basis": "sto-3g"
  },
  "configuration": {
    "method": "HF"
  },
  "execute_immediately": true
}'

echo -e "${YELLOW}Creating test experiment...${NC}"
create_response=$(curl -s -X POST "${API_V1}/experiments" \
    -H "Content-Type: application/json" \
    -d "$experiment_data")

experiment_id=$(echo "$create_response" | grep -o '"id":[0-9]*' | grep -o '[0-9]*' | head -n 1)

if [ -n "$experiment_id" ]; then
    echo -e "${GREEN}✓ Experiment created with ID: $experiment_id${NC}"
    echo ""

    # Wait a bit for experiment to complete
    sleep 3

    # 7. Get Experiment Details
    test_endpoint "GET" "${API_V1}/experiments/${experiment_id}" "" "Get Experiment Details"

    # 8. Get Experiment Status
    test_endpoint "GET" "${API_V1}/experiments/${experiment_id}/status" "" "Get Experiment Status"

    # 9. Get Convergence Data
    test_endpoint "GET" "${API_V1}/experiments/${experiment_id}/convergence" "" "Get Convergence Data"

    # 10. Get Circuit Visualization - JSON (NEW)
    test_endpoint "GET" "${API_V1}/experiments/${experiment_id}/circuit?format=json" "" "Circuit Visualization (JSON)"

    # 11. Get Circuit Visualization - ASCII (NEW)
    test_endpoint "GET" "${API_V1}/experiments/${experiment_id}/circuit?format=ascii" "" "Circuit Visualization (ASCII)"

    # 12. Get Circuit Visualization - QASM (NEW)
    test_endpoint "GET" "${API_V1}/experiments/${experiment_id}/circuit?format=qasm" "" "Circuit Visualization (QASM)"

    # 13. Get Experiment Report - JSON (NEW)
    test_endpoint "GET" "${API_V1}/experiments/${experiment_id}/report?format=json" "" "Experiment Report (JSON)"

    # 14. Get Experiment Report - Markdown (NEW)
    test_endpoint "GET" "${API_V1}/experiments/${experiment_id}/report?format=markdown" "" "Experiment Report (Markdown)"

else
    echo -e "${RED}✗ Failed to create experiment${NC}"
    echo "Response: $create_response"
fi

echo ""
echo "========================================="
echo "Testing Complete"
echo "========================================="
echo ""
echo "For detailed API documentation, visit:"
echo "  - Swagger UI: ${BASE_URL}/docs"
echo "  - ReDoc: ${BASE_URL}/redoc"
echo ""
echo "New endpoints added:"
echo "  - GET /api/v1/experiments/{id}/circuit"
echo "  - GET /api/v1/experiments/{id}/report"
echo "  - GET /api/v1/queue/stats"
echo ""
