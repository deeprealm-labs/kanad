# Application Modules Fix Plan

**Date:** November 6, 2025
**Status:** Analysis Complete - Ready to Fix

---

## Issues Identified

### 1. Materials Scout - CRITICAL ❌
**File:** [kanad/applications/materials_scout.py](kanad/applications/materials_scout.py)

**Problems:**
- **Line 769**: `total_gap += 2.0 * frac` - Default 2.0 eV for unknown materials
- **Lines 777-779**: Adds RANDOM NOISE to predictions (`np.random.normal`) - NON-DETERMINISTIC!

**Root Cause:**
The `_predict_bandgap()` method (line 745) is a fast screening filter, but:
1. Uses random noise (unacceptable for scientific computing)
2. Has limited material database (only 9 materials)
3. Defaults to arbitrary 2.0 eV

**Context:**
- Line 670 shows that `compute_dos_properties()` DOES compute real quantum HOMO-LUMO gaps
- So `_predict_bandgap()` is for initial screening, real calculation happens later
- BUT: Random noise and poor defaults make screening unreliable

**Solution:**
1. **REMOVE random noise entirely** (lines 777-779) - Use deterministic prediction
2. **Expand material database** with literature values
3. **Better default estimation** for unknown materials:
   - Use electronegativity difference to estimate band gap
   - Formula: `E_g ≈ C * |χ_A - χ_B|^n` where χ is electronegativity
   - Or use valence electron count correlation

**Impact:** HIGH - Affects all materials screening results

---

### 2. ADME Calculator - MEDIUM ⚠️
**File:** [kanad/analysis/adme_calculator.py](kanad/analysis/adme_calculator.py)

**Problems:**
- **Line 262**: `hb_contrib = -(desc.h_bond_donors * 0.5 + desc.h_bond_acceptors * 0.3)`
- **Line 265**: `aromatic_contrib = desc.aromatic_rings * 0.5`

**Root Cause:**
- Line 204 comment mentions "Wildman-Crippen method" - a validated literature approach
- BUT: The coefficients 0.5 and 0.3 are simplified/estimated
- Full Wildman-Crippen uses atom-type-specific contributions, not just counts

**Context:**
- Lines 269-276 show quantum corrections ARE being applied (HOMO-LUMO, polarizability)
- This is a hybrid QSPR + quantum approach
- References cite JCIM 2023, 2025 papers on quantum ML for ADME

**Solution Options:**

**Option A - Literature Validation:**
1. Find published Wildman-Crippen coefficients for H-bond groups
2. Add citations to code
3. Validate against known logP values (drugs, test set)

**Option B - Pure Quantum Approach:**
1. Rely more heavily on quantum descriptors (polarizability, dipole, HOMO-LUMO)
2. Use solvation free energy from quantum calculations
3. Fit quantum descriptors to experimental logP data

**Recommendation:** Option A (faster, validates current approach)
- Add proper citations
- Test against known compounds
- Document limitations

**Impact:** MEDIUM - Affects drug discovery predictions

---

### 3. Catalyst Optimizer - HIGH ⚠️
**File:** [kanad/applications/catalyst_optimizer.py](kanad/applications/catalyst_optimizer.py)

**Problems:**
- **Line 369**: `TOF = k_rate * 0.1` - Coverage factor of 0.1 (10% surface coverage)
- **Lines 382-385**: Elementary step barriers as fractions: 0.3, 1.0, 0.5, 0.2

**Root Cause:**
- Coverage factor should be computed from adsorption energies via Langmuir isotherm
- Elementary step barriers should come from transition state calculations or BEP relations

**Context:**
- Line 363-366 shows correct transition state theory (TST) is used: `k = (k_B*T/h) * exp(-E_act/k_B*T)`
- The math is correct, but input parameters are placeholders

**Solution:**

**For Coverage Factor (Line 369):**
```python
# BEFORE:
TOF = k_rate * 0.1  # Arbitrary 10% coverage

# AFTER: Use Langmuir isotherm
def _compute_coverage(self, E_ads, T=298):
    \"\"\"
    Compute surface coverage using Langmuir isotherm.

    θ = K*P / (1 + K*P)
    where K = exp(-ΔG_ads/RT)
    \"\"\"
    R = 1.987e-3  # kcal/(mol·K)
    P = 1.0  # atm (standard pressure)
    K = np.exp(-E_ads / (R * T))
    theta = K * P / (1 + K * P)
    return max(0.01, min(theta, 0.99))  # Clamp to 1-99%

TOF = k_rate * self._compute_coverage(E_ads, T)
```

**For Elementary Steps (Lines 382-385):**
```python
# BEFORE: Arbitrary fractions
elementary_steps=[
    {'name': 'Reactant adsorption', 'barrier': E_act * 0.3},
    ...
]

# AFTER: Use BEP (Brønsted-Evans-Polanyi) relations
def _estimate_elementary_barriers(self, E_act_rds, reaction_energy):
    \"\"\"
    Estimate barriers using BEP relations.

    E_barrier = E_0 + α * ΔE_rxn

    Typical α values:
    - Adsorption: α ≈ 0.2
    - Surface reactions: α ≈ 0.5-0.8
    - Desorption: α ≈ 0.2
    \"\"\"
    return [
        {'name': 'Reactant adsorption',
         'barrier': max(0.1, 0.2 * abs(reaction_energy))},  # BEP
        {'name': 'Bond activation',
         'barrier': E_act_rds},  # RDS (provided)
        {'name': 'Product formation',
         'barrier': max(0.1, 0.5 * abs(reaction_energy))},  # BEP
        {'name': 'Product desorption',
         'barrier': max(0.1, 0.2 * abs(reaction_energy))}   # BEP
    ]
```

**Impact:** HIGH - Affects catalysis predictions for chemical industry

---

## Priority Order

1. **Materials Scout** (CRITICAL) - Remove random noise, expand database
2. **Catalyst Optimizer** (HIGH) - Implement Langmuir coverage and BEP relations
3. **ADME Calculator** (MEDIUM) - Add literature citations and validation

---

## Implementation Steps

### Phase 1: Materials Scout Fix (1-2 hours)
1. Remove random noise (lines 777-779)
2. Expand material database from 9 to 30+ materials
3. Add electronegativity-based estimation for unknown materials
4. Test screening results for consistency

### Phase 2: Catalyst Optimizer Fix (2-3 hours)
1. Implement `_compute_coverage()` method with Langmuir isotherm
2. Implement `_estimate_elementary_barriers()` with BEP relations
3. Add adsorption energy calculation (or require as input)
4. Test against known catalytic systems

### Phase 3: ADME Calculator Validation (1-2 hours)
1. Research Wildman-Crippen original coefficients
2. Add citations to code
3. Create test set of known drugs (aspirin, ibuprofen, caffeine)
4. Validate logP predictions against experimental data
5. Document accuracy and limitations

---

## Expected Outcomes

### After Fixes:
- ✅ All predictions are deterministic (no random noise)
- ✅ Physics-based models (Langmuir, BEP, electronegativity trends)
- ✅ Literature-validated where possible
- ✅ Clearly documented limitations
- ✅ Test coverage for all fixes

### Remaining Limitations (ACCEPTABLE):
- ADME: Simplified QSPR model (not full atom-type Wildman-Crippen)
- Catalyst: BEP relations are approximate (real TS calculations would be better)
- Materials: Initial screening uses lookup table (real quantum calc in second pass)

**These limitations are acceptable because:**
1. They serve specific purposes (fast screening, initial estimates)
2. More accurate calculations ARE performed in second passes
3. Limitations are clearly documented
4. Results are deterministic and physics-based

---

**Ready to Implement:** Yes
**Estimated Time:** 4-7 hours total
**Risk:** Low - All changes are improvements over current placeholders
