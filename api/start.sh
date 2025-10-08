#!/bin/bash
# Kanad API Startup Script

echo "Starting Kanad Quantum Chemistry API..."

# Check if we're in the right directory
if [ ! -f "main.py" ]; then
    echo "Error: main.py not found. Please run this script from /home/mk/deeprealm/kanad/api/"
    exit 1
fi

# Check Python version
PYTHON_VERSION=$(python3 --version 2>&1 | grep -Po '(?<=Python )(.+)')
REQUIRED_VERSION="3.9"

if [ "$(printf '%s\n' "$REQUIRED_VERSION" "$PYTHON_VERSION" | sort -V | head -n1)" != "$REQUIRED_VERSION" ]; then
    echo "Error: Python 3.9+ required. Found: $PYTHON_VERSION"
    exit 1
fi

# Check if virtual environment exists
if [ ! -d "../venv" ]; then
    echo "Warning: No virtual environment found. Consider creating one:"
    echo "  python3 -m venv ../venv"
    echo "  source ../venv/bin/activate"
    echo "  pip install -r requirements.txt"
fi

# Create .env if it doesn't exist
if [ ! -f ".env" ]; then
    echo "Creating .env from .env.example..."
    cp .env.example .env
fi

# Start the server
echo "Starting server on http://localhost:8000"
echo "Documentation: http://localhost:8000/docs"
echo ""
echo "Press Ctrl+C to stop"
echo ""

python3 -m uvicorn api.main:app --reload --host 0.0.0.0 --port 8000
