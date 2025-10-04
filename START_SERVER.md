# 🚀 How to Start Your Kanad Server

## Quick Start

```bash
./run_server.sh
```

Then open in your browser: **http://localhost:8000**

---

## What You'll See

### 1. Web Interface (http://localhost:8000)
- Beautiful purple gradient interface
- Quick preset buttons: H2, LiH, H2O, Na2
- Atom builder (add/remove atoms)
- Solver selection (HF, VQE, QPE, SQD)
- Optimization toggle
- Real-time results with energy and speedup stats

### 2. API Documentation (http://localhost:8000/api/docs)
- Interactive Swagger UI
- Try all API endpoints
- See request/response schemas
- Test with real data

### 3. API Endpoints
- `GET /api/health` - Check server status
- `GET /api/presets` - Get preset molecules
- `POST /api/compute` - Calculate energy
- `POST /api/optimize` - Analyze optimization
- `GET /api/solvers` - List available solvers

---

## Step-by-Step

### 1. Start the Server
```bash
# Option 1: Use the script
./run_server.sh

# Option 2: Manual
source env/bin/activate
python -m uvicorn api.server:app --host 0.0.0.0 --port 8000 --reload
```

You'll see:
```
🚀 Starting Kanad Quantum Chemistry API Server...

Activating virtual environment...
Server will be available at:
  - Web App: http://localhost:8000
  - API Docs: http://localhost:8000/api/docs

Press Ctrl+C to stop the server

INFO:     Serving static files from /home/mk/deeprealm/kanad/web
INFO:     Uvicorn running on http://0.0.0.0:8000
```

### 2. Open Your Browser
Go to: **http://localhost:8000**

### 3. Test with H2 Molecule
1. Click the **"H₂"** preset button
2. Select **"VQE (Accurate)"** solver
3. Check **"Enable Quantum Optimization"**
4. Click **"🚀 Compute Energy"**
5. See results in ~0.1 seconds!

### 4. Try the API
Open another terminal:
```bash
curl http://localhost:8000/api/health
```

You'll get:
```json
{
  "status": "healthy",
  "version": "1.0.0",
  "timestamp": 1234567890.123
}
```

---

## How It Works

### The Web App
- **Frontend**: Static HTML/CSS/JavaScript in `web/` folder
- **Served by**: FastAPI static file mounting
- **Calls**: Makes POST requests to `/api/compute`
- **Shows**: Energy, system info, optimization stats in real-time

### The API
- **Framework**: FastAPI (Python web framework)
- **Server**: Uvicorn (ASGI server)
- **Routes**: All prefixed with `/api/*`
- **Static files**: Mounted at `/` (root)

### File Structure
```
kanad/
├── api/
│   ├── __init__.py
│   └── server.py          ← FastAPI server
├── web/
│   ├── index.html         ← Web interface
│   ├── styles.css         ← Styling
│   └── app.js             ← JavaScript logic
└── run_server.sh          ← Startup script
```

---

## Troubleshooting

### Port Already in Use
If you see `Address already in use`, change the port:
```bash
python -m uvicorn api.server:app --port 8001
```

Then use: http://localhost:8001

### Module Not Found
Make sure you activated the environment:
```bash
source env/bin/activate
```

And installed dependencies:
```bash
pip install -r requirements.txt
```

### Web Page Shows API JSON
This means you're accessing `/api` instead of `/`.
Use: **http://localhost:8000** (not `/api`)

### Static Files Not Loading
Check that `web/` directory exists:
```bash
ls -la web/
# Should show: index.html, styles.css, app.js
```

---

## Example API Usage

### Compute H2 Energy
```bash
curl -X POST http://localhost:8000/api/compute \
  -H "Content-Type: application/json" \
  -d '{
    "atoms": [
      {"symbol": "H", "position": [0.0, 0.0, 0.0]},
      {"symbol": "H", "position": [0.74, 0.0, 0.0]}
    ],
    "bond_type": "covalent",
    "solver": "VQE",
    "optimize": true,
    "strategy": "balanced"
  }'
```

### Get Presets
```bash
curl http://localhost:8000/api/presets
```

### Analyze Optimization
```bash
curl -X POST http://localhost:8000/api/optimize \
  -H "Content-Type: application/json" \
  -d '{
    "atoms": [
      {"symbol": "Li", "position": [0.0, 0.0, 0.0]},
      {"symbol": "H", "position": [1.59, 0.0, 0.0]}
    ],
    "bond_type": "ionic",
    "strategy": "aggressive"
  }'
```

---

## What Happens When You Click "Compute"

1. **Web UI** (`app.js`) collects atom data
2. **JavaScript** makes POST to `/api/compute`
3. **FastAPI** receives request in `server.py`
4. **Backend** creates atoms, bond, Hamiltonian
5. **Solver** (VQE/QPE/SQD) computes energy
6. **Optimizer** (if enabled) reduces qubits
7. **Response** returns energy + system info
8. **Web UI** displays results with animations

All in 0.1-1 second! ⚡

---

## Next Steps

✅ **Local**: Server running on localhost:8000
🌐 **Deploy**: See [DEPLOYMENT_READY.md](DEPLOYMENT_READY.md)
📚 **Docs**: http://localhost:8000/api/docs
🔬 **Test**: Click presets and see quantum chemistry in action!

---

**Ready to go!** Start the server and open http://localhost:8000 🚀
