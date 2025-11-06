# Real Issues Audit - November 6, 2025

**Status:** INVESTIGATION IN PROGRESS

User is correct - previous claims of "all fixed" were premature. This document tracks ACTUAL unresolved issues.

---

## Issue #1: Governance Integration - PARTIALLY RESOLVED ⚠️

**Status:** SQD uses governance generators, but VQE ansatz consistency unclear

**Evidence:**
- SQD solver has validation filtering (lines 177-183, 240-246)
- Need to verify VQE ansatz uses same governance protocol

**Location:** kanad/solvers/vqe_solver.py

**Action Required:** Verify VQE ansatz construction uses governance protocol consistently

---

## Issue #2: Raman Quantum Path - NOT RESOLVED ❌

**Status:** Returns HF polarizability instead of using quantum 1-RDM

**Evidence:**
- File: kanad/analysis/raman_calculator.py
- Line 319: `alpha_classical = self._compute_polarizability(method='HF')`
- Line 339: `return alpha_classical` - Returns HF, not quantum
- Lines 323-329: TODO exists but not implemented

**Root Cause:** Finite-field quantum polarizability not implemented

**Action Required:**
1. Implement finite-field calculation: Apply ±F field
2. Recompute energy using quantum 1-RDM  
3. Extract α from energy curvature: α = -∂²E/∂F²

---

## Issue #3: Environment Effects - NOT RESOLVED ❌

**Status:** Multiple 0.0 fallbacks found

**Evidence:**
1. **pressure.py:325** - Returns 0.0 if no energy cached
2. **temperature.py:437** - Returns 0.0 if no excited states
3. **solvent.py:340, 367, 493, 534** - Multiple 0.0 returns
4. **ph_effects.py:110, 468, 490** - Multiple 0.0 returns

**Action Required:** Replace all 0.0 placeholders with proper calculations

---

## Issue #4: Spectroscopy Excited States - NOT RESOLVED ❌

**Status:** TD-DFT/QPE paths are placeholders

**Evidence:** Need to investigate spectroscopy modules

**Action Required:** 
1. Find placeholder energy offsets
2. Locate TD-DFT/QPE placeholder paths
3. Find excited-state Hessian TODO

---

## Issue #5: PDOS - NOT RESOLVED ❌

**Status:** Equal split among atoms, no orbital projection

**Evidence:** Need to search for equal split logic

**Action Required:** Implement orbital projection matrices for PDOS

---

## Issue #6: Fast VQE Expectation - NOT RESOLVED ❌

**Status:** Uses HF energy placeholder instead of SparsePauliOp expectation

**Evidence:** Need to check VQE solver expectation calculation

**Action Required:** Compute proper SparsePauliOp expectation value

---

## Next Steps

1. Complete investigation of all 6 issues
2. Create specific fix plans for each
3. Implement fixes systematically
4. Validate each fix with tests
5. NO premature celebration until ALL issues resolved

**No green checkmarks until actual validation shows issues are fixed.**

