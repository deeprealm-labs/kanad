#!/bin/bash

# Kanad API Server Startup Script

echo "🚀 Starting Kanad API Server..."
echo ""

# Check if in correct directory
if [ ! -d "api" ]; then
    echo "❌ Error: Please run this script from the kanad root directory"
    exit 1
fi

# Check if virtual environment exists
if [ ! -d "venv" ] && [ ! -d ".venv" ]; then
    echo "⚠️  Warning: No virtual environment found"
    echo "   Consider creating one: python -m venv venv"
    echo ""
fi

# Set environment variables (optional - you can also use .env file)
export API_HOST="${API_HOST:-0.0.0.0}"
export API_PORT="${API_PORT:-8000}"
export DEBUG="${DEBUG:-true}"

# Optional: Load credentials from environment
# export IBM_API="your_token"
# export IBM_CRN="your_crn"
# export BLUE_TOKEN="your_token"

echo "📋 Configuration:"
echo "   Host: $API_HOST"
echo "   Port: $API_PORT"
echo "   Debug: $DEBUG"
echo ""

# Start the server
cd api
echo "🔧 Installing/checking dependencies..."
pip install -q -r requirements.txt

echo ""
echo "✅ Server starting at http://localhost:$API_PORT"
echo "📚 API docs available at http://localhost:$API_PORT/docs"
echo ""
echo "Press Ctrl+C to stop the server"
echo ""

python main.py
