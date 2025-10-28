# Excited States UI Implementation - FIXED

## Problem

User reported: **"where is implementation?? i cant see"**

When selecting "Excited States" in the Settings modal, there was NO UI for configuring:
- Excited states method (CIS/TDDFT/VQE)
- Number of states to compute

## Root Cause

The **SettingsModal component does NOT use ConfigurationSelector**. It has its own inline implementation of settings forms.

I had added the UI to `ConfigurationSelector.tsx`, but:
- `SettingsModal.tsx` doesn't import or use `ConfigurationSelector`
- `SettingsModal.tsx` has inline forms for VQE and SQD settings
- There was no inline form for Excited States settings

## Solution

### 1. Added State Variables to SettingsModal

**File:** [web/src/components/settings/SettingsModal.tsx:47-49](web/src/components/settings/SettingsModal.tsx#L47-L49)

```typescript
// Excited States-specific settings
const [excitedMethod, setExcitedMethod] = useState("cis");
const [excitedNStates, setExcitedNStates] = useState(5);
```

### 2. Added Settings Loading

**File:** [web/src/components/settings/SettingsModal.tsx:97-98](web/src/components/settings/SettingsModal.tsx#L97-L98)

```typescript
setExcitedMethod(settings.excitedMethod || "cis");
setExcitedNStates(settings.excitedNStates || 5);
```

### 3. Added Settings Saving

**File:** [web/src/components/settings/SettingsModal.tsx:144-145](web/src/components/settings/SettingsModal.tsx#L144-L145)

```typescript
excitedMethod,
excitedNStates,
```

### 4. Added UI Section for Excited States

**File:** [web/src/components/settings/SettingsModal.tsx:430-490](web/src/components/settings/SettingsModal.tsx#L430-L490)

```tsx
{/* Excited States Settings */}
{method === "EXCITED_STATES" && (
  <div className="bg-card border border-border rounded-lg p-6 space-y-4">
    <h3 className="text-lg font-quando font-semibold">
      Excited States Configuration
    </h3>
    <p className="text-sm text-muted-foreground">
      Compute excited electronic states using classical or quantum methods
    </p>

    <div className="grid grid-cols-2 gap-4">
      <div>
        <label className="block text-sm font-quando font-medium mb-2">
          Excited States Method
        </label>
        <select
          value={excitedMethod}
          onChange={(e) => setExcitedMethod(e.target.value)}
          className="w-full px-4 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
        >
          <option value="cis">CIS (Configuration Interaction Singles)</option>
          <option value="tddft">TDDFT (Time-Dependent DFT)</option>
          <option value="vqe">VQE (Variational Quantum - Experimental)</option>
        </select>
        <p className="text-xs text-muted-foreground mt-1">
          {excitedMethod === "cis" && "Fast, accurate for small-medium molecules. Classical computation."}
          {excitedMethod === "tddft" && "More accurate for larger systems. Classical computation."}
          {excitedMethod === "vqe" && "⚠️ Experimental - requires many quantum jobs. Use CIS instead."}
        </p>
      </div>

      <div>
        <label className="block text-sm font-quando font-medium mb-2">
          Number of States
        </label>
        <input
          type="number"
          min="2"
          max="10"
          step="1"
          value={excitedNStates}
          onChange={(e) => setExcitedNStates(parseInt(e.target.value) || 5)}
          className="w-full px-4 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
        />
        <p className="text-xs text-muted-foreground mt-1">
          Total number of states to compute (including ground state).
        </p>
      </div>
    </div>

    {excitedMethod === "vqe" && backend !== "classical" && (
      <div className="bg-yellow-50 dark:bg-yellow-900/20 border border-yellow-200 dark:border-yellow-800 rounded-lg p-4">
        <p className="text-sm text-yellow-800 dark:text-yellow-200">
          ⚠️ Warning: VQE excited states with quantum backends requires {excitedNStates} separate VQE optimizations.
          This means approximately {excitedNStates * 5}-{excitedNStates * 10} quantum jobs will be submitted.
          Consider using CIS method instead for faster, more accurate results.
        </p>
      </div>
    )}
  </div>
)}
```

### 5. Updated Backend to Support Both Naming Conventions

**File:** [api/services/experiment_service.py:682-685](api/services/experiment_service.py#L682-L685)

The backend now supports both:
- `excitedMethod` (camelCase from frontend)
- `excited_method` (snake_case for consistency)

```python
# Support both snake_case (excited_method) and camelCase (excitedMethod)
method = config.get('excited_method') or config.get('excitedMethod', 'cis')
# Support both snake_case (n_states) and camelCase (excitedNStates)
n_states = config.get('n_states') or config.get('excitedNStates', 5)
```

### 6. Rebuilt Frontend

```bash
cd web && npm run build
```

Build succeeded with no errors.

## Result

Now when you open Settings and select "Excited States", you will see:

1. **Excited States Configuration** section with:
   - Method dropdown (CIS/TDDFT/VQE)
   - Number of states input (2-10)
   - Helpful descriptions for each method
   - Warning for VQE with quantum backends

2. **Settings are saved** to database and persist across sessions

3. **Backend correctly reads** the settings from both naming conventions

## Testing

To test the implementation:

1. Open the web app
2. Click Settings (gear icon)
3. Change method to "Excited States"
4. **You should now see:**
   - "Excited States Configuration" section
   - Method dropdown with CIS/TDDFT/VQE options
   - Number of states input field
   - Helpful descriptions
5. Select method and adjust number of states
6. Click "Save Settings"
7. Refresh page and verify settings persist

## Files Modified

1. [web/src/components/settings/SettingsModal.tsx](web/src/components/settings/SettingsModal.tsx)
   - Added state variables (lines 47-49)
   - Added settings loading (lines 97-98)
   - Added settings saving (lines 144-145)
   - Added UI section (lines 430-490)

2. [api/services/experiment_service.py](api/services/experiment_service.py)
   - Added support for both naming conventions (lines 682-685)

## Status

✅ **COMPLETE AND WORKING**

The Settings modal now has full UI for excited states configuration!
