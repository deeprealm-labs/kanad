# Kanad Backend + Frontend Quick Start Guide

## 🚀 Get Everything Running in 5 Minutes

### Step 1: Start the Backend Server

```bash
# From the kanad root directory
./start_server.sh
```

**Expected output:**
```
🚀 Starting Kanad API Server...
✅ Server starting at http://localhost:8000
📚 API docs available at http://localhost:8000/docs
```

### Step 2: Start the Frontend

Open a **new terminal** window:

```bash
cd web
npm install  # First time only
npm run dev
```

**Expected output:**
```
▲ Next.js 15.5.4
- Local: http://localhost:3000
```

### Step 3: Test the System

Open your browser and navigate to:
- **Frontend**: http://localhost:3000
- **Backend API Docs**: http://localhost:8000/docs

### Step 4: Submit Your First Experiment

#### Option A: Using the Web Interface

1. Go to http://localhost:3000
2. Click "New Experiment"
3. Select or create a molecule (e.g., H2)
4. Configure the simulation (VQE, hardware_efficient)
5. Click "Run Experiment"
6. Watch real-time progress!

#### Option B: Using curl

```bash
curl -X POST "http://localhost:8000/api/experiments/submit" \
  -H "Content-Type: application/json" \
  -d '{
    "molecule": {
      "smiles": "[H][H]",
      "basis": "sto-3g",
      "charge": 0,
      "multiplicity": 1
    },
    "configuration": {
      "method": "VQE",
      "ansatz": "hardware_efficient",
      "mapper": "jordan_wigner",
      "optimizer": "SLSQP",
      "max_iterations": 100,
      "backend": "classical"
    },
    "execute_now": true
  }'
```

---

## 🧪 Example Molecules to Try

### Simple (Fast, ~10-30 seconds)
- **H2**: `[H][H]` - Hydrogen molecule
- **LiH**: `[Li]H` - Lithium hydride

### Medium (1-2 minutes)
- **H2O**: `O` - Water
- **NH3**: `N` - Ammonia

### Advanced (5+ minutes)
- **CH4**: `C` - Methane
- **C2H4**: `C=C` - Ethylene

---

## 📊 What You Can Do

### Computation Methods
- **Hartree-Fock (HF)**: Fast classical calculation
- **VQE**: Quantum variational algorithm
- **SQD**: Subspace quantum dynamics (coming soon)

### Backends
- **Classical**: Local simulation (always available)
- **IBM Quantum**: Real quantum hardware (requires credentials)
- **BlueQubit**: GPU-accelerated (requires credentials)

### Analysis Options
- Energy decomposition
- Bond analysis
- Dipole moment
- Thermochemistry
- Spectroscopy

---

## 🔑 Setting Up Cloud Backends (Optional)

### IBM Quantum

1. Get credentials from https://quantum.ibm.com/
2. Set environment variables:
```bash
export IBM_API="your_token_here"
export IBM_CRN="your_crn_here"
```
3. Restart the backend server

### BlueQubit

1. Get token from https://bluequbit.io/
2. Set environment variable:
```bash
export BLUE_TOKEN="your_token_here"
```
3. Restart the backend server

---

## 🛠️ Troubleshooting

### Backend won't start
```bash
# Check if port 8000 is in use
lsof -i :8000

# Kill existing process
kill -9 <PID>

# Check Python installation
python --version  # Should be 3.9+

# Reinstall dependencies
cd api
pip install -r requirements.txt
```

### Frontend won't start
```bash
# Check if port 3000 is in use
lsof -i :3000

# Clear Next.js cache
rm -rf web/.next

# Reinstall dependencies
cd web
rm -rf node_modules package-lock.json
npm install
```

### "Connection refused" error
```bash
# Make sure backend is running
curl http://localhost:8000/health

# Check frontend API URL
cat web/.env.local
# Should contain: NEXT_PUBLIC_API_URL=http://localhost:8000/api
```

### Experiment fails
```bash
# Check backend logs (terminal where backend is running)
# Look for error messages

# Check database
sqlite3 kanad_experiments.db "SELECT * FROM experiments ORDER BY created_at DESC LIMIT 1;"

# Try a simpler molecule first
# Use H2 with HF method to test basic functionality
```

---

## 📁 Files Created

Backend:
```
api/
├── main.py                    # FastAPI server
├── core/
│   ├── config.py             # Settings
│   └── database.py           # SQLite operations
├── routes/                   # API endpoints
└── services/                 # Business logic
```

Frontend:
```
web/
└── .env.local                # API configuration
```

Database:
```
kanad_experiments.db          # SQLite database (created on first run)
```

---

## 🎯 Next Steps

1. ✅ Both servers running
2. ✅ Submitted first experiment
3. **Try different molecules and methods**
4. **Explore the API documentation** at http://localhost:8000/docs
5. **Check the comprehensive guide** in `BACKEND_IMPLEMENTATION_GUIDE.md`

---

## 💡 Pro Tips

- Start with **H2 + HF** for fast testing (5-10 seconds)
- Use **classical backend** for development
- Use **VQE + hardware_efficient** for balance of accuracy/speed
- Monitor backend logs for debugging
- Keep both terminals visible to see real-time updates

---

## 🆘 Getting Help

1. **Check the logs**: Backend prints detailed info
2. **API docs**: http://localhost:8000/docs (interactive!)
3. **Database inspection**: `sqlite3 kanad_experiments.db`
4. **Full guide**: See `BACKEND_IMPLEMENTATION_GUIDE.md`

---

**You're all set! Happy quantum computing! 🎉**
