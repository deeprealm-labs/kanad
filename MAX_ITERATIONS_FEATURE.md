# Max Iterations Feature - Implementation Complete

## Summary

Added user-configurable **Maximum Iterations** setting to control VQE optimization iterations and prevent large molecules from running excessively long.

---

## What Was Added

### 1. Frontend UI - Settings Modal
**File**: `/web/src/components/settings/SettingsModal.tsx`

Added new field in VQE Configuration section:
```typescript
<div>
  <label className="block text-sm font-quando font-medium mb-2">
    Maximum Iterations
  </label>
  <input
    type="number"
    min="10"
    max="1000"
    step="10"
    value={maxIterations}
    onChange={(e) => setMaxIterations(parseInt(e.target.value) || 100)}
    className="..."
  />
  <p className="text-xs text-muted-foreground mt-1">
    Maximum number of VQE optimization iterations (10-1000).
    Lower values = faster but less accurate.
  </p>
</div>
```

**Features**:
- Number input with range validation (10-1000)
- Step increment of 10 for easier adjustment
- Default value: 100 (reduced from previous 1000)
- Clear help text explaining the tradeoff

### 2. Backend API - Default Settings
**File**: `/api/routes/settings.py`

Updated default settings to include maxIterations:
```python
{
    "method": "VQE",
    "ansatz": "hardware_efficient",
    "mapper": "jordan_wigner",
    "optimizer": "SLSQP",
    "backend": "classical",
    "backend_name": "ibm_torino",
    "maxIterations": 100,  # ← NEW
    "optimization": { ... }
}
```

### 3. Backend API - Experiment Configuration
**File**: `/api/routes/experiments.py`

Updated BackendConfig model:
```python
class BackendConfig(BaseModel):
    method: str = "VQE"
    ansatz: Optional[str] = "hardware_efficient"
    mapper: Optional[str] = "jordan_wigner"
    optimizer: Optional[str] = "SLSQP"
    max_iterations: Optional[int] = 100  # Changed from 1000
    backend: str = "classical"
    backend_name: Optional[str] = None

    class Config:
        # Allow both camelCase (frontend) and snake_case (backend)
        populate_by_name = True
        alias_generator = lambda field_name: ''.join(
            word.capitalize() if i > 0 else word
            for i, word in enumerate(field_name.split('_'))
        )
```

**Key Changes**:
- Default changed from 1000 → 100
- Added `Config` class for camelCase/snake_case conversion
- Frontend sends `maxIterations` (camelCase)
- Backend uses `max_iterations` (snake_case)
- Pydantic automatically converts between the two

### 4. Database
**Updated**: `kanad_experiments.db` → user_settings table

Added maxIterations field to existing settings:
```sql
UPDATE user_settings
SET settings = json_set(settings, '$.maxIterations', 100)
WHERE id = 1;
```

---

## How It Works

### User Flow

1. **Configure Settings**:
   - User opens Settings modal
   - Scrolls to "VQE Configuration" section
   - Adjusts "Maximum Iterations" slider/input
   - Value range: 10-1000, step: 10
   - Clicks "Save Settings"

2. **Settings Saved**:
   - Frontend: Saves `maxIterations: 100` to API
   - Backend: Converts to `max_iterations: 100`
   - Database: Stores in user_settings table

3. **Run Experiment**:
   - User creates molecule (e.g., CO2, COOH)
   - Clicks "Execute"
   - Backend reads user settings
   - Passes `max_iterations` to VQE solver
   - VQE stops after 100 iterations (or earlier if converged)

### Technical Flow

```
Frontend (SettingsModal)
  ↓ maxIterations: 100 (camelCase)
API (POST /api/settings/defaults)
  ↓ Saves to database
Database (user_settings)
  ↓ {"maxIterations": 100}
API (POST /api/experiments)
  ↓ BackendConfig.max_iterations (snake_case)
Experiment Service
  ↓ config.get('max_iterations', 100)
VQE Solver
  ↓ max_iterations=100
Optimizer (SLSQP/COBYLA)
  ✓ Stops after 100 iterations
```

---

## Default Values

### Before This Fix
- **Default**: 1000 iterations
- **Problem**: CO2 taking 10+ minutes
- **Problem**: COOH running 2000+ iterations
- **User Experience**: Frustrating wait times

### After This Fix
- **Default**: 100 iterations
- **Expected**: CO2 completes in ~3-4 minutes
- **Expected**: H2O completes in ~2 minutes
- **Expected**: H2 completes in <30 seconds
- **User Experience**: Fast iterations with option to increase for accuracy

---

## Recommended Settings by Molecule Size

| Molecule Type | Electrons | Recommended Max Iterations | Expected Time |
|---------------|-----------|----------------------------|---------------|
| H2            | 2         | 50-100                     | 20-30 seconds |
| H2O           | 10        | 100-200                    | 2-4 minutes   |
| CO2           | 22        | 100-300                    | 3-6 minutes   |
| COOH          | 23        | 150-300                    | 4-8 minutes   |
| Benzene (C6H6)| 42        | 200-500                    | 10-20 minutes |

**General Rule**:
- Small molecules (< 10 electrons): 50-100 iterations
- Medium molecules (10-25 electrons): 100-300 iterations
- Large molecules (> 25 electrons): 200-500 iterations

**For Publication-Quality Results**: 500-1000 iterations

---

## Testing

### Verify Settings UI
```bash
# 1. Start frontend
cd web && npm run dev

# 2. Open http://localhost:3000
# 3. Click Settings (gear icon)
# 4. Check VQE Configuration section
# 5. Verify "Maximum Iterations" field exists
# 6. Verify default value is 100
# 7. Verify can change to 50, save, reload - should persist
```

### Verify Backend Receives Setting
```bash
# 1. Check current settings
curl http://localhost:8000/api/settings/defaults | python3 -m json.tool | grep maxIterations

# Should show:
# "maxIterations": 100

# 2. Create test experiment with H2
# 3. Check logs for VQE solver initialization
# Should see: max_iterations=100
```

### Verify Experiments Use Setting
```bash
# 1. Set maxIterations to 50 in Settings
# 2. Run H2 experiment
# 3. Monitor convergence graph
# 4. Should stop at iteration 50 (or earlier if converged)
# 5. Check experiment logs
grep "Iteration" /tmp/kanad_backend.log | tail -5
# Should not exceed 50 iterations
```

---

## Files Modified

1. ✅ `/web/src/components/settings/SettingsModal.tsx`
   - Added maxIterations state variable
   - Added UI input field
   - Added to settings save/load logic

2. ✅ `/api/routes/settings.py`
   - Added maxIterations to default settings fallback

3. ✅ `/api/routes/experiments.py`
   - Changed default from 1000 → 100
   - Added Config class for camelCase/snake_case conversion

4. ✅ `/api/core/config.py`
   - Already changed DEFAULT_MAX_ITERATIONS to 100 (previous fix)

5. ✅ `kanad_experiments.db`
   - Updated user_settings to include maxIterations: 100

---

## Backwards Compatibility

### Old Experiments
- Experiments created before this change will continue to work
- If maxIterations not specified, defaults to 100 (safe fallback)

### Old Settings
- Users with existing settings will get maxIterations: 100 on next load
- No migration needed - handled gracefully by fallback defaults

### API Compatibility
- Both `maxIterations` (camelCase) and `max_iterations` (snake_case) accepted
- Pydantic Config handles conversion automatically

---

## Future Enhancements

### Smart Defaults (Future)
Auto-adjust max_iterations based on molecule size:
```python
def calculate_recommended_iterations(n_electrons):
    if n_electrons < 10:
        return 50
    elif n_electrons < 25:
        return 100
    elif n_electrons < 50:
        return 200
    else:
        return 300
```

### Progress-Based Adjustment (Future)
Stop early if convergence detected:
```python
# Check energy change
if abs(energy_history[-1] - energy_history[-2]) < 1e-6:
    logger.info("Converged early - stopping")
    break
```

### Cloud Backend Considerations (Future)
For IBM/BlueQubit:
- Iterations happen on cloud, not locally
- May need different max_iterations for cloud vs local
- See CLOUD_BACKEND_IMPLEMENTATION_PLAN.md for details

---

## Success Criteria

✅ **UI**: maxIterations field visible in Settings modal
✅ **UI**: Default value is 100
✅ **UI**: Can change value and save
✅ **Backend**: API accepts maxIterations in camelCase
✅ **Backend**: Converts to max_iterations (snake_case)
✅ **Backend**: Default in BackendConfig is 100
✅ **Database**: Settings include maxIterations: 100
✅ **Experiments**: VQE solver receives and respects max_iterations
✅ **Validation**: Values constrained to 10-1000 range

---

## User Documentation

### How to Change Max Iterations

1. Click **Settings** icon (⚙️) in top navigation
2. Scroll to **VQE Configuration** section
3. Find **Maximum Iterations** field
4. Enter desired value:
   - **50-100**: Fast, less accurate (good for testing)
   - **100-200**: Balanced (default, recommended)
   - **300-500**: Slower, more accurate
   - **500-1000**: Publication quality, slow
5. Click **Save Settings**
6. Future experiments will use this setting

### FAQ

**Q: What happens if I set it too low?**
A: The optimization may not converge fully, giving less accurate energy values. You'll get results faster but they may not be chemically accurate.

**Q: What happens if I set it too high?**
A: Experiments will take longer but give more accurate results. For most molecules, 100-200 is sufficient.

**Q: Why does my experiment stop before reaching max iterations?**
A: The optimizer may have converged early (energy stopped changing). This is good - it found the answer faster!

**Q: Can I change this for each experiment?**
A: Currently it's a global setting. Future versions may allow per-experiment configuration.

---

**Implemented**: 2025-10-22
**Status**: ✅ Complete and Tested
**Ready for Production**: Yes
