# Basis Set Frontend Fix - Disable Unimplemented Options

## Problem

Frontend offered basis set options that weren't implemented in the backend, causing experiments to fail with:

```
❌ Experiment failed: Basis set 'cc-pvdz' not implemented yet
NotImplementedError: Basis set 'cc-pvdz' not implemented yet
```

## Root Cause

**Frontend** (`web/src/components/molecule/`)
```tsx
<select>
  <option value="sto-3g">STO-3G</option>
  <option value="6-31g">6-31G</option>
  <option value="6-31g*">6-31G*</option>       ← Offered but not implemented!
  <option value="cc-pvdz">cc-pVDZ</option>     ← Offered but not implemented!
  <option value="cc-pvtz">cc-pVTZ</option>     ← Offered but not implemented!
</select>
```

**Backend** (`kanad/core/integrals/basis_sets.py:647-652`)
```python
if self.basis_name == 'sto-3g':
    self._add_sto3g_functions(atom)
elif self.basis_name == '6-31g':
    self._add_6_31g_functions(atom)
else:
    raise NotImplementedError(f"Basis set '{self.basis_name}' not implemented yet")
```

Only **2 basis sets are implemented**:
- `sto-3g` ✅
- `6-31g` ✅

But **3 basis sets were offered in UI but not implemented**:
- `6-31g*` ❌
- `cc-pvdz` ❌
- `cc-pvtz` ❌

## Solution

Disabled unimplemented basis sets in frontend dropdowns so users can't select them.

### Changes Made

#### 1. MoleculeCreator.tsx

```tsx
<select value={basis} onChange={(e) => setBasis(e.target.value)}>
  <option value="sto-3g">STO-3G</option>
  <option value="6-31g">6-31G</option>
  <option value="6-31g*" disabled>6-31G* (not implemented)</option>
  <option value="cc-pvdz" disabled>cc-pVDZ (not implemented)</option>
  <option value="cc-pvtz" disabled>cc-pVTZ (not implemented)</option>
</select>
```

#### 2. MoleculeBuilder.tsx

Updated all 3 instances of basis set dropdowns (same change as above).

## Result

Users now see:
- **STO-3G** - clickable ✅
- **6-31G** - clickable ✅
- **6-31G* (not implemented)** - grayed out, not selectable ⚠️
- **cc-pVDZ (not implemented)** - grayed out, not selectable ⚠️
- **cc-pVTZ (not implemented)** - grayed out, not selectable ⚠️

## Testing Results

### Before Fix ❌
```
User selects: cc-pvdz
Backend error: NotImplementedError: Basis set 'cc-pvdz' not implemented yet
Experiment fails
```

### After Fix ✅
```
User can only select: sto-3g or 6-31g
Experiments complete successfully
```

### Verified Behavior

**STO-3G** (H2):
- 2 electrons, 2 orbitals
- Excited States: ground + 1 excited = 2 states total
- Works ✅

**6-31G** (H2):
- 2 electrons, 4 orbitals
- Excited States: ground + 3 excited = 4 states total
- Works ✅

## Files Modified

1. **web/src/components/molecule/MoleculeCreator.tsx:389-391**
   - Added `disabled` attribute to unimplemented basis sets

2. **web/src/components/molecule/MoleculeBuilder.tsx:197-199, 348-350, 508-510**
   - Added `disabled` attribute to unimplemented basis sets (3 locations)

## Future Work

To implement additional basis sets, add to `kanad/core/integrals/basis_sets.py`:

### Required Steps

1. **Add basis set data**:
   ```python
   def _add_6_31g_star_functions(self, atom: 'Atom') -> None:
       """Add 6-31G* basis functions (includes polarization)."""
       # Implement polarization functions (d orbitals on heavy atoms)
       pass
   ```

2. **Update build_basis() switch**:
   ```python
   elif self.basis_name == '6-31g*':
       self._add_6_31g_star_functions(atom)
   elif self.basis_name == 'cc-pvdz':
       self._add_cc_pvdz_functions(atom)
   ```

3. **Remove `disabled` from frontend dropdowns**

### Basis Set Implementation Priority

1. **6-31G*** - Adds d polarization functions (most requested)
2. **cc-pVDZ** - Correlation-consistent double-zeta
3. **cc-pVTZ** - Correlation-consistent triple-zeta (large, slow)

### Implementation Complexity

- **6-31G***: Moderate (add d functions to existing 6-31G)
- **cc-pVDZ**: High (completely different functional forms)
- **cc-pVTZ**: Very High (many more functions, very slow)

## Summary

**Issue**: Frontend offered 3 unimplemented basis sets causing experiment failures

**Fix**: Disabled unimplemented options in frontend UI

**Result**: Users can only select working basis sets (sto-3g, 6-31g)

**Status**: ✅ FIXED - No more basis set errors

---

**Date**: 2025-10-24
**Impact**: Prevents user confusion and failed experiments
**Future**: Implement 6-31G*, cc-pVDZ, cc-pVTZ basis sets in backend
