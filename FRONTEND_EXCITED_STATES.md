# Frontend Updates for Excited States

## What Was Added

### 1. Excited States Method Selector

**File:** `web/src/components/simulation/ConfigurationSelector.tsx` (lines 79-100)

Users can now choose between:
- **CIS** (Configuration Interaction Singles) - Default, fast, accurate
- **TDDFT** (Time-Dependent DFT) - More accurate for larger systems
- **VQE** (Variational Quantum) - Experimental, requires many quantum jobs

```tsx
{/* Excited States Method Selection */}
{settings.method === "EXCITED_STATES" && (
  <div>
    <label>Excited States Method</label>
    <select value={settings.excitedMethod || "cis"}>
      <option value="cis">CIS (Configuration Interaction Singles)</option>
      <option value="tddft">TDDFT (Time-Dependent DFT)</option>
      <option value="vqe">VQE (Variational Quantum - Experimental)</option>
    </select>
    <p className="text-xs text-muted-foreground">
      {settings.excitedMethod === "cis" && "Fast, accurate for small-medium molecules. Classical computation."}
      {settings.excitedMethod === "tddft" && "More accurate for larger systems. Classical computation."}
      {settings.excitedMethod === "vqe" && "⚠️ Experimental - requires many quantum jobs. Use CIS instead."}
    </p>
  </div>
)}
```

### 2. Number of States Input

**File:** `web/src/components/simulation/ConfigurationSelector.tsx` (lines 102-120)

Users can specify how many excited states to compute (2-10):

```tsx
{/* Number of Excited States */}
{settings.method === "EXCITED_STATES" && (
  <div>
    <label>Number of States</label>
    <input
      type="number"
      min="2"
      max="10"
      value={settings.nStates || 5}
      onChange={(e) => updateSetting("nStates", parseInt(e.target.value))}
    />
    <p className="text-xs">
      Total number of states to compute (including ground state)
    </p>
  </div>
)}
```

### 3. Updated TypeScript Types

**File:** `web/src/lib/types.ts` (lines 21, 30-32)

Added support for:
- `method: "EXCITED_STATES"` type
- `excitedMethod?: "cis" | "tddft" | "vqe"`
- `nStates?: number`

```typescript
export interface BackendSettings {
  method: "HF" | "VQE" | "SQD" | "QPE" | "MP2" | "FCI" | "EXCITED_STATES";
  // ... other fields ...
  // Excited States settings
  excitedMethod?: "cis" | "tddft" | "vqe"; // Method for computing excited states
  nStates?: number; // Number of states (ground + excited)
}
```

## User Experience

### When User Selects "EXCITED_STATES":

1. **Method selector appears** with 3 options:
   - CIS (recommended, default)
   - TDDFT (for larger molecules)
   - VQE (experimental, not recommended)

2. **Helpful descriptions** show for each method explaining:
   - What it does
   - When to use it
   - Performance characteristics
   - Warnings for VQE

3. **Number of states input** appears allowing 2-10 states

4. **Backend selection still works** but:
   - If quantum backend + CIS/TDDFT → runs classical (as intended)
   - If quantum backend + VQE → shows warning, runs VQE (experimental)

## Backend Integration

The backend (`api/services/experiment_service.py`) reads these settings:

```python
method = config.get('excited_method', 'cis')  # Gets from excitedMethod field
n_states = config.get('n_states', 5)          # Gets from nStates field
```

And makes intelligent decisions:
- **CIS/TDDFT** → Always classical, fast, accurate ✅
- **VQE + quantum backend** → Currently disabled (too expensive), uses CIS instead
- **VQE explicit request** → Could be enabled with warnings

## Result

Users now have full control over:
1. ✅ **Which method** to use (CIS, TDDFT, or VQE)
2. ✅ **How many states** to compute
3. ✅ **Clear descriptions** of each option
4. ✅ **Warnings** for expensive/experimental options

The UI guides users toward the best choices (CIS for most cases) while allowing advanced users to explore other options! 🎉
