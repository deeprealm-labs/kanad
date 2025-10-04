#!/bin/bash
# Kanad API Server - Quick Start Script

echo "🚀 Starting Kanad Quantum Chemistry API Server..."
echo ""

# Activate virtual environment if exists
if [ -d "env" ]; then
    echo "Activating virtual environment..."
    source env/bin/activate
fi

echo "Server will be available at:"
echo "  - Web App: http://localhost:8000"
echo "  - API Docs: http://localhost:8000/api/docs"
echo ""
echo "Press Ctrl+C to stop the server"
echo ""

# Run the server
python -m uvicorn api.server:app --host 0.0.0.0 --port 8000 --reload
