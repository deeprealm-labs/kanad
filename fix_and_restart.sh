#!/bin/bash

echo "🛑 Killing ALL API servers..."
pkill -9 -f "uvicorn.*main:app"
pkill -9 -f "restart_api.sh"
sleep 2

echo "🧹 Clearing ALL Python cache..."
find /home/mk/deeprealm/kanad -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null
find /home/mk/deeprealm/kanad -name "*.pyc" -delete 2>/dev/null
find /home/mk/deeprealm/kanad -name "*.pyo" -delete 2>/dev/null

echo "✅ Cache cleared"

echo "🚀 Starting CLEAN API server..."
cd /home/mk/deeprealm/kanad
source env/bin/activate
cd api
python3 -m uvicorn main:app --reload --port 8000 &

sleep 3

echo "✅ Server started. Check status:"
ss -tlnp | grep :8000
ps aux | grep "uvicorn.*main:app" | grep -v grep
