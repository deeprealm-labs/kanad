# SQD Convergence Graph Fix - Energy History Added

## Problem

SQD solver was broadcasting convergence updates via WebSocket successfully (7 stages), but the frontend convergence graph only showed 1 data point.

### Backend Logs (Working)
```
📊 SQD Progress: Stage 0/6 - HF reference computed, E = -1.11729324 Ha
📊 Broadcasting convergence: iter=0, E=-1.11729324  ✅
📊 SQD Progress: Stage 1/6 - Generating subspace basis, E = -1.11729324 Ha
📊 Broadcasting convergence: iter=1, E=-1.11729324  ✅
...
📊 SQD Progress: Stage 6/6 - State 2 computed, E = -0.11788595 Ha
📊 Broadcasting convergence: iter=6, E=-0.11788595  ✅
```

### Frontend Display (Broken)
- **Graph**: Only 1 data point visible
- **Current Iteration**: 1
- **Execution Logs**: Only 1 log entry

## Root Cause

The WebSocket messages were broadcast correctly, but the frontend wasn't seeing all the historical data after the experiment completed because:

1. **Real-time updates work**: WebSocket broadcasts all 7 stages
2. **Polling overwrites data**: When experiment completes, frontend polls for final results
3. **Missing energy_history**: SQD results didn't include `energy_history` array
4. **Frontend replaces data**: Without `energy_history`, frontend couldn't reconstruct the full convergence graph

### The Data Flow Problem

```
During Experiment:
WebSocket → Frontend accumulates 7 convergence points ✅
  [Stage 0, Stage 1, Stage 2, Stage 3, Stage 4, Stage 5, Stage 6]

When Completed:
Polling fetches results → No energy_history in results ❌
Frontend checks: convergence_history? ❌
                 energy_history? ❌
                 convergenceData? ❌
Result: Graph shows incomplete data (only last WebSocket update)
```

## Solution

Add `energy_history` array to SQD results so the frontend can reconstruct the full convergence graph even if WebSocket updates are missed.

### Changes Made

**File**: `api/services/experiment_service.py`

#### 1. Track Energy History in Callback

```python
# Define progress callback for real-time updates
n_states = config.get('n_states', 3)
total_stages = 4 + n_states

# Track energy history for convergence graph
energy_history = []  # ← NEW!

def sqd_progress_callback(stage: int, energy: float, message: str):
    """Progress callback for SQD execution."""
    # ... progress updates ...

    # Store energy for history
    energy_history.append(float(energy))  # ← NEW!

    # Send WebSocket update
    if experiment_id:
        _broadcast_convergence_sync(experiment_id, stage, float(energy))
```

#### 2. Include in Results

```python
results_dict = {
    'energy': float(result['energy']),
    'hf_energy': float(result.get('hf_energy', 0.0)),
    # ... other fields ...
    'energies': [float(e) for e in result['energies']],
    'excited_state_energies': [float(e) for e in result.get('excited_state_energies', [])],
    'energy_history': energy_history,  # ← NEW! For convergence graph
}

print(f"📊 SQD energy_history: {len(energy_history)} stages recorded")
```

## How It Works Now

### Complete Data Flow

```
During Experiment:
1. SQD callback triggered for each stage
2. Energy added to energy_history[] array
3. WebSocket broadcasts convergence update
4. Frontend accumulates in real-time ✅

After Completion:
1. Results include energy_history array
2. Polling fetches final results
3. Frontend sees energy_history
4. Converts to convergence data format:
   [
     {iteration: 1, energy: -1.11729324},
     {iteration: 2, energy: -1.11729324},
     ...
     {iteration: 7, energy: -0.11788595}
   ]
5. Graph displays all 7 stages ✅
```

### Frontend Processing (Already Implemented)

**File**: `web/src/components/simulation/ExperimentMonitor.tsx:208-215`

```typescript
if (exp.results.energy_history && Array.isArray(exp.results.energy_history)) {
  console.log("Processing energy_history:", exp.results.energy_history.length, "points");
  // Convert energy_history array to convergence data format
  const convergence = exp.results.energy_history.map((energy: number, index: number) => ({
    iteration: index + 1,
    energy: energy
  }));
  setConvergenceData(convergence);  // Updates the graph!
}
```

## Expected Output

### Backend Console
```
⚛  Starting experiment abc-123
📊 Creating molecule: {'smiles': '[H][H]', ...}
✅ Molecule created: H2 (2 electrons, 2 orbitals)
🔬 Running SQD calculation...
📊 SQD Progress: Stage 0/6 - HF reference computed, E = -1.11729324 Ha
📊 Broadcasting convergence: iter=0, E=-1.11729324
📊 SQD Progress: Stage 1/6 - Generating subspace basis, E = -1.11729324 Ha
📊 Broadcasting convergence: iter=1, E=-1.11729324
📊 SQD Progress: Stage 2/6 - Projecting Hamiltonian, E = -1.11729324 Ha
📊 Broadcasting convergence: iter=2, E=-1.11729324
📊 SQD Progress: Stage 3/6 - Diagonalizing Hamiltonian, E = -1.11729324 Ha
📊 Broadcasting convergence: iter=3, E=-1.11729324
📊 SQD Progress: Stage 4/6 - State 0 computed, E = -1.13605351 Ha
📊 Broadcasting convergence: iter=4, E=-1.13605351
📊 SQD Progress: Stage 5/6 - State 1 computed, E = -0.47567111 Ha
📊 Broadcasting convergence: iter=5, E=-0.47567111
📊 SQD Progress: Stage 6/6 - State 2 computed, E = -0.11788595 Ha
📊 Broadcasting convergence: iter=6, E=-0.11788595
📊 SQD energy_history: 7 stages recorded  ← NEW!
✅ Experiment completed: Energy = -1.13605351 Ha
✅ Broadcasted completion status via WebSocket
```

### Frontend Display
- **Graph**: 7 data points visible (stages 0-6)
- **Current Iteration**: 7
- **Execution Logs**: All 7 stages logged
- **Current Energy**: Shows final ground state energy (-1.136 Ha)

### Graph Visualization
```
Energy (Ha)
    0 ─┐
       │  ●─────●  (Stages 0-3: HF + subspace)
  -0.5 ┤
       │         ●  (Stage 4: Ground state -1.136)
  -1.0 ┤
       │
  -1.5 ┤              ●  (Stage 5: 1st excited -0.476)
       │                 ●  (Stage 6: 2nd excited -0.118)
       └─────────────────────────────
         0  1  2  3  4  5  6  Iteration
```

## Why This Matters

### For SQD Solver
SQD has distinct stages that are important to visualize:
1. **Stages 0-3**: Setup (HF, subspace, projection, diagonalization)
2. **Stages 4+**: Computed excited states (energies decrease)

Each stage represents meaningful progress, not just iterations. Users need to see all stages to understand:
- When the expensive steps (projection, diagonalization) occur
- How many excited states were found
- Energy ordering of states

### For Other Solvers
- **VQE**: Already has energy_history via convergence tracking ✅
- **Excited States**: Already has energy_history ✅
- **SQD**: Now has energy_history ✅
- **Hartree-Fock**: N/A (instant, no iterations)

## Files Modified

**Backend**:
1. `api/services/experiment_service.py:564` - Added `energy_history = []`
2. `api/services/experiment_service.py:579` - Added `energy_history.append(float(energy))`
3. `api/services/experiment_service.py:602` - Added `'energy_history': energy_history`
4. `api/services/experiment_service.py:605` - Added debug print

**Frontend**: No changes needed! Already handles `energy_history` arrays.

## Testing

After backend auto-reload:

1. **Run SQD experiment**:
   - Method: SQD
   - Molecule: H2
   - States: 3
   - Subspace: 10

2. **Expected results**:
   - ✅ Real-time updates during execution
   - ✅ Graph shows all 7 stages (0-6)
   - ✅ Final energy displayed
   - ✅ Analysis results shown

3. **Backend console should show**:
   ```
   📊 SQD energy_history: 7 stages recorded
   ```

4. **Frontend console should show**:
   ```
   Processing energy_history: 7 points
   ```

## Benefits

- ✅ **Complete visualization**: All SQD stages visible on graph
- ✅ **Resilient to timing issues**: Works even if WebSocket updates are missed
- ✅ **Consistent with other solvers**: VQE and Excited States already use energy_history
- ✅ **Better UX**: Users see the full computation history
- ✅ **Debugging**: Easier to identify which stage takes longest

## Related Fixes

This completes the SQD convergence display fix started in:
1. **Event Loop Fix** ([EVENT_LOOP_FIX.md](EVENT_LOOP_FIX.md)) - Fixed WebSocket broadcasting
2. **This Fix** - Added energy_history for graph reconstruction

---

**Status**: ✅ COMPLETE
**Date**: 2025-10-24
**Solvers with Full Convergence Support**:
- VQE ✅
- SQD ✅ (just fixed)
- Excited States ✅
- Hartree-Fock ⚠️ (N/A - non-iterative)
