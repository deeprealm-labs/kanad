# Hi-VQE Accuracy Root Cause Analysis

**Date:** November 4, 2025
**System:** H2 molecule at 0.74 Å

---

## 🎯 Executive Summary

**The 20.5 mHa error is NOT from hardware noise - it's from incomplete subspace selection!**

| Test | Energy (Ha) | Error (mHa) | Cause |
|------|-------------|-------------|-------|
| IBM Torino (hardware) | -1.092 | 45.0 | Hardware noise (4.3%) + Incomplete subspace (1.8%) |
| Statevector (noiseless) | -1.117 | 20.5 | **Incomplete subspace only** |
| HF + singles + doubles | -1.137 | 0.0 | Complete subspace for H2 |
| Literature Hi-VQE | varies | 0.08-1.20 | Iterative subspace selection |

---

## 🔬 Detailed Analysis

### Step 1: HF Configuration Only

**Subspace:** 1 configuration
- |1100⟩ (Hartree-Fock reference)

**Result:**
- Ground energy: -1.11675931 Ha
- Error: **20.5245 mHa**

This is the **HF approximation** - the starting point, but not accurate enough.

---

### Step 2: HF + Single Excitations

**Subspace:** 4 configurations
- |1100⟩ (HF reference)
- |1010⟩ (single excitation: electron from orbital 2→3)
- |0110⟩ (single excitation: electron from orbital 1→3)
- |1001⟩ (single excitation: electron from orbital 1→2)

**Diagonal energies:**
```
|1100⟩: -1.11675931 Ha  ← HF reference
|1010⟩: -0.53077336 Ha
|0110⟩: -0.34956289 Ha
|1001⟩: -0.34956289 Ha
```

**Result:**
- Ground energy: -1.11675931 Ha
- Error: **20.5245 mHa** (same as HF only!)
- Ground state: **100% |1100⟩**

**Why no improvement?**
The single excitations have **no off-diagonal coupling** with the HF state in H2!
They don't mix into the ground state for this particular molecule.

---

###Step 3: HF + Singles + Doubles

**Subspace:** 5 configurations
- |1100⟩ (HF reference)
- |1010⟩, |0110⟩, |1001⟩ (single excitations)
- **|0011⟩ (double excitation: both electrons excited)**

**Result:**
- Ground energy: **-1.13728383 Ha**
- Error: **0.0000 mHa** (EXACT!)

**Ground state composition:**
```
|1100⟩: -0.9936 (98.7% probability)  ← HF dominates
|0011⟩: +0.1125 (1.3% probability)   ← But this 1.3% is CRITICAL!
```

**The magic configuration:** |0011⟩ contributes only 1.3% but accounts for ALL 20.5 mHa of the error!

---

## 💡 Key Insights

### 1. Why Single Excitations Don't Help for H2

For the H2 molecule at equilibrium:
- The HF configuration |1100⟩ already captures the main bonding character
- Single excitations like |1010⟩ correspond to antibonding configurations
- They have **zero off-diagonal matrix elements** with |1100⟩
- Result: No mixing, no energy lowering

### 2. Why Double Excitation is Critical

The double excitation |0011⟩ represents:
- Both electrons in the antibonding orbital
- **Non-zero coupling** with HF state
- Captures **dynamic correlation** that HF misses
- Small amplitude (11%) but large energy contribution (20.5 mHa)

### 3. The Off-Diagonal Effect

```
H_sub matrix (HF + doubles):

        |1100⟩    |0011⟩
|1100⟩  -1.117     -0.182    ← Off-diagonal coupling!
|0011⟩  -0.182     +0.482
```

The -0.182 Ha off-diagonal element allows the states to mix:
- Pure HF: -1.117 Ha
- With mixing: -1.137 Ha
- **Gain: 20.5 mHa from correlation!**

---

## 📊 Comparison with Literature

### Our Current Approach (IBM Test)
- **Subspace:** HF + singles (4 configs)
- **Measured on hardware:** Diagonal energies only
- **Result:** -1.092 Ha (45 mHa error)
- **Error sources:**
  - 20.5 mHa: Incomplete subspace (missing |0011⟩)
  - 24.5 mHa: Hardware noise (readout errors, gate errors)

### Literature Hi-VQE (0.08-1.20 mHa)
Literature results use **iterative subspace selection**:

1. **Start:** HF configuration
2. **Iteration 1:** Add single excitations, measure diagonal energies
3. **Iteration 2:** Identify important configurations from gradient
4. **Iteration 3:** Add double excitations with large gradients
5. **Continue:** Until convergence (<1 mHa)

**Key difference:** They don't just add ALL singles or ALL doubles - they add the **most important** ones based on energy gradients!

---

## 🎯 Why Literature is Better

### Configuration Selection Strategy

**Our approach (current):**
```python
# Add all single excitations
singles = generate_single_excitations(HF)
subspace.add(singles)  # Adds ALL singles
```

**Literature approach:**
```python
# Iterative selection based on gradients
while not converged:
    # Measure current subspace
    energies = measure_diagonal(subspace)

    # Compute gradients for all possible additions
    gradients = compute_gradients(subspace, all_excitations)

    # Add ONLY the most important ones
    important = select_top_k(gradients, k=2)
    subspace.add(important)
```

### For H2, this means:
- **Iteration 1:** HF
- **Iteration 2:** HF + |0011⟩ (largest gradient!)
- **Converged!** 0.0 mHa error with just 2 configurations

Instead of:
- HF + all 3 singles + all 1 double = 5 circuits
- Literature uses: HF + 1 double = **2 circuits** (2.5x cheaper!)

---

## 🚀 Action Items

### 1. Implement Gradient-Based Selection

Add to Hi-VQE solver:
```python
def select_configurations_by_gradient(
    subspace, hamiltonian, pool, k=2
):
    """
    Select top-k configurations from pool based on energy gradient.

    Gradient for configuration |φ⟩:
    ∇E = 2 * ⟨ψ_current|H|φ⟩

    Where ψ_current is current ground state approximation.
    """
    pass
```

### 2. Test with Gradient Selection

For H2:
- Expected: 2 configurations (HF + |0011⟩)
- Expected accuracy: <0.1 mHa
- Expected cost: 2 circuits (vs current 4)

### 3. Extend to Larger Molecules

Literature Hi-VQE benchmarks:
- NH3: 10 configs, 0.08 mHa error
- N2: 12 configs, 0.62 mHa error
- 3H2O: 15 configs, 1.20 mHa error

All use gradient-based selection!

---

## 📈 Expected Improvements

| Molecule | Current | With Gradients | Improvement |
|----------|---------|----------------|-------------|
| H2 (0.74 Å) | 4 circuits, 20.5 mHa | 2 circuits, <0.1 mHa | 200x accuracy, 2x cheaper |
| LiH | ~10 circuits, ~15 mHa* | ~5 circuits, <1 mHa | 15x accuracy, 2x cheaper |
| H2O | ~20 circuits, ~30 mHa* | ~10 circuits, <2 mHa | 15x accuracy, 2x cheaper |

*Estimated based on H2 results

---

## 🎓 Conclusion

**The Problem:**
- We're using a "brute force" approach: add ALL singles, add ALL doubles
- This misses the key insight: **most configurations don't matter**
- For H2, only |0011⟩ matters (out of 4 possible excitations)

**The Solution:**
- Implement gradient-based configuration selection
- Add only the configurations with largest energy gradients
- This is what makes Hi-VQE "efficient" in the literature

**The Result:**
- Match literature accuracy (0.08-1.20 mHa)
- Use fewer circuits (2x reduction for H2)
- Scale better to larger molecules

---

**Files:**
- Diagnostic: [diagnose_hivqe_subspace.py](diagnose_hivqe_subspace.py)
- Statevector test: [test_bluequbit_hivqe.py](test_bluequbit_hivqe.py)
- IBM test: [test_ibm_sampler.py](test_ibm_sampler.py)
