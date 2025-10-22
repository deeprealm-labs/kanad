#!/bin/bash

echo "🛑 Stopping API server..."
pkill -9 -f "uvicorn main:app"
sleep 2

echo "🧹 Clearing Python cache..."
find api -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null
find kanad -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null

echo "🚀 Starting API server..."
cd api
python3 -m uvicorn main:app --reload --port 8000 --log-level debug

