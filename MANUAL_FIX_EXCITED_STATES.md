# Manual Fix for Excited States - Step by Step

## The Problem

The code IS fixed, but Python bytecode cache is preventing the new code from running. Your terminal shows old messages that don't exist in the source code anymore.

## Manual Fix Steps

### Step 1: Kill ALL API Processes

```bash
cd /home/mk/deeprealm/kanad

# Kill everything
pkill -9 -f uvicorn
pkill -9 -f restart_api
pkill -9 -f "python.*multiprocessing"

# Wait for processes to die
sleep 3
```

### Step 2: Delete ALL Python Cache

```bash
# Delete all __pycache__ directories
find . -type d -name "__pycache__" -print0 | xargs -0 rm -rf

# Delete all .pyc files
find . -name "*.pyc" -delete

# Verify they're gone
find . -name "__pycache__" -o -name "*.pyc"
# Should show nothing
```

### Step 3: Verify the Code is Correct

```bash
# Check line 630-650 in experiment_service.py
sed -n '630,650p' /home/mk/deeprealm/kanad/api/services/experiment_service.py
```

You should see:
```python
# Get backend configuration
    backend_kwargs = {}

    # Get user-selected method FIRST - respect user's choice!
    # Support both snake_case (excited_method) and camelCase (excitedMethod)
    user_method = config.get('excited_method') or config.get('excitedMethod', 'cis')

    print(f"🔍 User selected excited states method: {user_method}")

    # Only configure backend_kwargs if method requires quantum execution
    if user_method == 'vqe' and backend_type in ['bluequbit', 'ibm_quantum']:
        print(f"⚠️  VQE excited states with quantum backend - requires many jobs")
```

**If you DON'T see this**, the file wasn't saved properly. Copy the code from the bottom of this document.

### Step 4: Restart API Fresh

```bash
cd /home/mk/deeprealm/kanad

# Activate environment
. env/bin/activate

# Start API manually (not with restart_api.sh to avoid multiple processes)
cd api
python3 -m uvicorn main:app --reload --port 8000
```

### Step 5: Test with CIS

1. Open web app: http://localhost:3000
2. Go to Settings
3. Select method: **Excited States**
4. Select ES Method: **CIS**
5. Set States: **3**
6. Backend: **BlueQubit** or **Classical**
7. Save settings
8. Run H2 experiment

###Step 6: Check Terminal Logs

You should now see:
```
🔬 Running Excited States calculation...
🔍 User selected excited states method: cis
📊 CIS is classical-only, switching backend to classical
🔍 DEBUG Excited States Config:
   config.get('excited_method'): None
   config.get('excitedMethod'): cis
   FINAL method: cis
   FINAL n_states: 3
   Backend type: classical
🔬 Starting CIS excited states calculation...
📊 Computing 3 excited states
⚙️ Running CIS solver...
✅ CIS calculation complete
```

**NO "using VQE method" message!**
**NO quantum jobs submitted!**

## If It Still Doesn't Work

The issue is that `config.get('excitedMethod')` is returning `None` because the field isn't being sent from the frontend.

Check if the BackendConfig model has the fields:

```bash
grep -A 5 "excit" /home/mk/deeprealm/kanad/api/routes/experiments.py
```

Should show:
```python
# Excited States specific fields
excited_method: Optional[str] = "cis"  # cis, tddft, vqe
excited_n_states: Optional[int] = 5  # Number of excited states to compute
```

## Complete Code for experiment_service.py (lines 630-690)

If the file doesn't have this, copy and paste:

```python
    # Get backend configuration
    backend_kwargs = {}

    # Get user-selected method FIRST - respect user's choice!
    # Support both snake_case (excited_method) and camelCase (excitedMethod)
    user_method = config.get('excited_method') or config.get('excitedMethod', 'cis')

    print(f"🔍 User selected excited states method: {user_method}")

    # Only configure backend_kwargs if method requires quantum execution
    if user_method == 'vqe' and backend_type in ['bluequbit', 'ibm_quantum']:
        print(f"⚠️  VQE excited states with quantum backend - requires many jobs")
        # Get backend credentials for VQE
        backend_kwargs = get_backend_kwargs(backend_type, config)
    elif backend_type in ['bluequbit', 'ibm_quantum'] and user_method != 'vqe':
        # User selected CIS or TDDFT with quantum backend - override to classical
        print(f"📊 {user_method.upper()} is classical-only, switching backend to classical")
        backend_type = 'classical'

    # Classical excited states (CIS/TDDFT) - create bond wrapper
    if molecule.n_atoms == 2:
        # Diatomic molecule - use bond API
        atom_1, atom_2 = molecule.atoms
        distance = atom_1.distance_to(atom_2)
        bond = BondFactory.create_bond(
            atom_1,
            atom_2,
            distance=distance,
            basis=molecule.basis
        )
        print(f"✅ Created bond for diatomic molecule: {atom_1.symbol}-{atom_2.symbol}")
    else:
        # Multi-atom molecule - create a minimal bond wrapper
        # CIS/TDDFT only need hamiltonian and molecule, not actual bond properties
        print(f"⚠️  Multi-atom molecule ({molecule.n_atoms} atoms) - creating wrapper for Excited States")

        # Create a simple wrapper class that mimics a bond
        class MoleculeBondWrapper:
            def __init__(self, molecule):
                self.molecule = molecule
                self.hamiltonian = molecule.hamiltonian
                self.atoms = molecule.atoms
                self.bond_type = "molecular"  # Not a traditional bond
                self.bond_order = None
                self.distance = None

        bond = MoleculeBondWrapper(molecule)
        print(f"✅ Created molecule wrapper for classical excited states")

    # Use the user-selected method (already read above)
    method = user_method

    # Support both snake_case (n_states) and camelCase (excitedNStates)
    n_states = config.get('n_states') or config.get('excitedNStates', 5)

    print(f"🔍 DEBUG Excited States Config:")
    print(f"   config.get('excited_method'): {config.get('excited_method')}")
    print(f"   config.get('excitedMethod'): {config.get('excitedMethod')}")
    print(f"   FINAL method: {method}")
    print(f"   FINAL n_states: {n_states}")
    print(f"   Backend type: {backend_type}")

    # Broadcast initial status to frontend
    if experiment_id:
        _broadcast_log_sync(experiment_id, f"🔬 Starting {method.upper()} excited states calculation...")
        _broadcast_log_sync(experiment_id, f"📊 Computing {n_states} excited states")
```

## Summary

The fix IS in the code, but Python is loading old bytecode. Follow the steps above carefully:

1. ✅ Kill everything
2. ✅ Delete all cache
3. ✅ Verify code is correct
4. ✅ Start API fresh (manually, one process only)
5. ✅ Test with CIS
6. ✅ Check logs show "User selected excited states method: cis"

The new logs with "🔍 User selected excited states method:" MUST appear, or the cache wasn't cleared properly.
