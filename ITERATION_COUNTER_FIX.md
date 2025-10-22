# Iteration Counter Fix - Detailed Implementation

## Problem
The "Current Iteration" display was showing stale values from previous experiments. For example:
- Run SQD experiment #1 → completes with 6 iterations → shows "Current Iteration: 6" ✓
- Run SQD experiment #2 → shows "Current Iteration: 9" ✗ (accumulated from previous)

## Root Cause
The `ExperimentMonitor` component was not being fully remounted when switching between experiments. React was reusing the same component instance and only updating props, which meant:
1. State variables (like `currentIteration`) persisted across experiments
2. The useEffect dependency on `experimentId` wasn't reliably triggering state reset
3. Component kept accumulated state from previous experiments

## Solution Applied

### Fix 1: Force Component Remount with Key Prop
**File**: `web/src/app/dashboard/page.tsx` (line 268)

```typescript
<ExperimentMonitor
  key={currentExperimentId}  // ← Force new component instance for each experiment
  experimentId={currentExperimentId}
  experimentConfig={experimentConfig}
  onComplete={...}
  onBack={...}
  onQueueAnother={...}
/>
```

**Why This Works**:
- React's `key` prop identifies component instances
- When `key` changes, React **destroys** the old component and **creates** a new one
- New component instance = fresh state, no accumulated values
- This is the recommended React pattern for resetting component state

### Fix 2: Enhanced State Reset Logic
**File**: `web/src/components/simulation/ExperimentMonitor.tsx` (lines 66-94)

```typescript
// Reset state when experiment ID changes (new experiment)
useEffect(() => {
  console.log("🔄 COMPONENT MOUNT/UPDATE - Experiment ID:", experimentId);
  console.log("🔄 Current iteration before reset:", currentIteration);

  // Reset all state
  setCurrentIteration(0);        // ← Reset iteration counter
  setConvergenceData([]);        // ← Clear convergence data
  setProgress(0);                // ← Reset progress bar
  setLogs([]);                   // ← Clear logs
  setStatus("queued");           // ← Reset status
  setResults(null);              // ← Clear results
  setExperiment(null);           // ← Clear experiment data

  console.log("✅ State reset complete for experiment:", experimentId);

  // Close any existing WebSocket
  if (wsRef.current) {
    console.log("🔌 Closing previous WebSocket");
    wsRef.current.close();
    wsRef.current = null;
  }

  // Stop any existing polling
  if (pollingRef.current) {
    console.log("⏹️  Stopping previous polling");
    pollingRef.current();
    pollingRef.current = null;
  }
}, [experimentId]);
```

**Why This Helps**:
- Acts as safety net if component isn't remounted
- Ensures all state is reset even if component is reused
- Cleans up WebSocket and polling connections
- Comprehensive logging for debugging

### Fix 3: Set Correct Final Iteration Count
**File**: `web/src/components/simulation/ExperimentMonitor.tsx` (lines 241-263)

```typescript
} else if (exp.status === "completed") {
  setProgress(100);

  // Set final iteration count based on convergence data from results
  if (exp.results) {
    let finalIterations = 0;

    // Check for convergence_history (VQE/SQD)
    if (exp.results.convergence_history && Array.isArray(exp.results.convergence_history)) {
      finalIterations = exp.results.convergence_history.length;
      console.log("📊 Setting final iterations from convergence_history:", finalIterations);
    }
    // Check for energy_history (alternative format)
    else if (exp.results.energy_history && Array.isArray(exp.results.energy_history)) {
      finalIterations = exp.results.energy_history.length;
      console.log("📊 Setting final iterations from energy_history:", finalIterations);
    }

    if (finalIterations > 0) {
      console.log("✅ Setting currentIteration to:", finalIterations);
      setCurrentIteration(finalIterations);
    } else {
      console.log("⚠️  No convergence data found in results");
    }
  }

  if (previousStatus !== "completed") {
    addLog(`Experiment completed successfully`);
  }
}
```

**Why This Matters**:
- Ensures iteration count matches actual data when experiment completes
- Handles both `convergence_history` and `energy_history` formats
- Prevents showing incorrect iteration count due to polling race conditions

### Fix 4: Enhanced Logging for Debugging
**File**: `web/src/components/simulation/ExperimentMonitor.tsx` (lines 307-317)

```typescript
case "convergence":
  console.log("📊 WebSocket convergence update - iteration:", message.iteration, "energy:", message.energy);
  setConvergenceData((prev) => [
    ...prev,
    {
      iteration: message.iteration,
      energy: message.energy,
    },
  ]);
  console.log("✅ Setting currentIteration to:", message.iteration);
  setCurrentIteration(message.iteration);
```

**Purpose**:
- Track every iteration update via WebSocket
- Verify iteration numbers are sequential and correct
- Identify if WebSocket is sending wrong iteration numbers

## How to Verify the Fix

### Test Scenario 1: Sequential Experiments
1. Run SQD experiment on H2 molecule
2. Wait for completion → Should show "Current Iteration: 6"
3. Click "Queue Another" or go back to dashboard
4. Run another SQD experiment
5. **Expected**: Counter starts at 0, progresses to 6
6. **Bug if**: Counter starts at 6 or shows any number > 6

### Test Scenario 2: Different Methods
1. Run SQD experiment → completes with 6 iterations
2. Run VQE experiment with maxIterations=100
3. **Expected**: Counter starts at 0, progresses toward 100
4. **Bug if**: Counter starts at 6 or shows accumulated value

### Test Scenario 3: Browser Console Verification
Open browser console (F12 → Console tab) and look for:

**When starting new experiment:**
```
🔄 COMPONENT MOUNT/UPDATE - Experiment ID: <new-id>
🔄 Current iteration before reset: <previous-value>
✅ State reset complete for experiment: <new-id>
```

**During experiment (WebSocket updates):**
```
📊 WebSocket convergence update - iteration: 0 energy: -1.1167
✅ Setting currentIteration to: 0
📊 WebSocket convergence update - iteration: 1 energy: -1.1167
✅ Setting currentIteration to: 1
...
📊 WebSocket convergence update - iteration: 6 energy: -0.1683
✅ Setting currentIteration to: 6
```

**When experiment completes:**
```
📊 Setting final iterations from convergence_history: 6
✅ Setting currentIteration to: 6
```

### Test Scenario 4: Rapid Switching
1. Start experiment #1
2. Immediately go back and start experiment #2
3. **Expected**: Experiment #1 cancelled/cleaned up, #2 starts fresh
4. **Bug if**: Iteration counter or data from #1 appears in #2

## Technical Details

### React Component Lifecycle with Key Prop

**Without `key` prop** (old behavior):
```
Experiment A starts → Component mounts → State: iteration=0
...A runs... → State: iteration=6
Experiment B starts → Component updates (same instance) → State: iteration=6 (NOT RESET!)
```

**With `key` prop** (new behavior):
```
Experiment A starts → Component mounts (key="id-A") → State: iteration=0
...A runs... → State: iteration=6
Experiment B starts → Component unmounts (key changed)
                    → Component mounts (key="id-B") → State: iteration=0 (FRESH!)
```

### State Reset Timing

The useEffect runs at these times:
1. **Component mount**: Initial state setup
2. **experimentId changes**: State reset for new experiment
3. **Component unmount**: Cleanup (via return function)

### Why Both Fixes Are Needed

1. **Key prop** (primary fix):
   - Guarantees fresh component instance
   - Most reliable way to reset state in React
   - Follows React best practices

2. **useEffect reset** (secondary fix):
   - Safety net if component is reused somehow
   - Handles edge cases
   - Useful for debugging (logging)

## Known Limitations

### Doesn't Fix (Out of Scope)
1. Backend not using quantum circuits for SQD
2. Analysis data not displaying in frontend
3. EXCITED_STATES redirect UX

### Assumes
1. `experimentId` prop changes when new experiment starts
2. Parent component (`dashboard/page.tsx`) properly manages experiment IDs
3. WebSocket/polling provides correct iteration numbers

## Files Modified

1. **web/src/app/dashboard/page.tsx**
   - Added `key` prop to `<ExperimentMonitor>`

2. **web/src/components/simulation/ExperimentMonitor.tsx**
   - Enhanced state reset useEffect
   - Added logging throughout
   - Fixed final iteration count logic
   - Added WebSocket convergence logging

## Testing Checklist

- [ ] Run SQD experiment → completes with correct iteration count
- [ ] Run second SQD → starts from 0, not accumulated
- [ ] Run VQE after SQD → counter resets properly
- [ ] Run EXCITED_STATES → counter shows correct value
- [ ] Check browser console for reset logs
- [ ] Check browser console for convergence logs
- [ ] Verify graph shows correct number of points
- [ ] Verify "Iteration X / Y" display is correct

## Rollback Plan

If this fix causes issues, revert these changes:

```bash
# Revert page.tsx
git checkout HEAD -- web/src/app/dashboard/page.tsx

# Revert ExperimentMonitor.tsx
git checkout HEAD -- web/src/components/simulation/ExperimentMonitor.tsx

# Rebuild
cd web && npm run build
```

## Related Issues

This fix addresses the iteration counter bug but doesn't solve:
- **Analysis display issue** (separate problem - data exists in backend but not showing in UI)
- **SQD quantum backend issue** (SQD always runs classically, never uses quantum circuits)
- **EXCITED_STATES redirect** (confusing UX when quantum backend selected)

See `FIXES_APPLIED.md` for comprehensive status of all issues.
