# CRITICAL BUG REPORT: VQE Bond API Issue
**Date:** November 4, 2025
**Severity:** CRITICAL
**Status:** IDENTIFIED

---

## 🚨 PROBLEM SUMMARY

VQE gives DIFFERENT results depending on API used:

| API Path | Final Energy | HF Energy | Correlation | Recovery |
|----------|--------------|-----------|-------------|----------|
| **Direct (tests)** | -1.136 Ha | -1.117 Ha | -0.019 Ha | 94% ✅ |
| **Bond API (dashboard)** | -1.115 Ha | -1.117 Ha | +0.002 Ha | -7% ❌ |

**Bond API gives energy ABOVE HF!** This violates the variational principle for the Hamiltonian.

---

## 📊 TEST RESULTS

### Test 1: Direct Hamiltonian API (Working)
```python
# From test_all_optimizers.py
pauli_hamiltonian = openfermion_jordan_wigner(...)
ansatz = CovalentGovernanceAnsatz(n_qubits=4, n_electrons=2)
solver = VQESolver(hamiltonian=pauli_hamiltonian, ansatz=ansatz, ...)
```

**Result:**
- VQE Energy: -1.13605 Ha
- HF Energy: -1.11729 Ha
- Correlation: -0.01876 Ha (94% recovery) ✅

###Test 2: Bond API (Broken)
```python
# From test_bond_api_vqe.py
bond = BondFactory.create_bond(H, H, distance=0.74)
solver = VQESolver(bond=bond, ansatz_type='governance', ...)
```

**Result:**
- VQE Energy: -1.11518 Ha ❌ ABOVE HF!
- HF Energy: -1.11676 Ha
- Correlation: +0.00158 Ha (wrong sign!) ❌

---

## 🔍 ROOT CAUSE ANALYSIS

### Issue is NOT:
- ❌ Nuclear repulsion calculation (both use 0.715 Ha correctly for 0.74 Å)
- ❌ Parameter binding (verified working)
- ❌ Optimizer (COBYLA works in direct API)
- ❌ Ansatz selection (both use CovalentGovernanceAnsatz)

### Issue IS:
**Bond API creates DIFFERENT Hamiltonian or ansatz circuit!**

Evidence:
1. Same molecule (H2 @ 0.74 Å)
2. Same basis (STO-3G)
3. Same optimizer (COBYLA)
4. Same ansatz type (Covalent Governance)
5. **DIFFERENT energy landscapes!**

The Bond API Hamiltonian or ansatz must be constructed differently, causing:
- Energy to START near HF
- Optimization to make it WORSE (go above HF)
- Positive correlation energy (impossible for true ground state)

---

## 🎯 HYPOTHESIS

The Bond API likely has ONE of these issues:

### Hypothesis 1: Wrong Hamiltonian Construction
Bond.hamiltonian might be missing correlation terms or using wrong mapper

### Hypothesis 2: Wrong Ansatz Initialization
Bond API might initialize ansatz differently, causing it to diverge from ground state

### Hypothesis 3: Wrong Initial State
HF state preparation might be wrong in Bond API path

---

## 🔬 EVIDENCE

### Energy Variation During Optimization:

**Direct API (Good):**
```
Initial: -1.11082 Ha
Final:   -1.13605 Ha
Direction: DOWNWARD ✅
Improvement: -0.025 Ha ✅
```

**Bond API (Bad):**
```
Initial: -1.10894 Ha
Final:   -1.11518 Ha
Direction: UPWARD to -1.115 but still > HF ❌
Ended ABOVE HF energy ❌
```

---

## 🧪 NEXT STEPS

1. **Compare Hamiltonians:**
   - Print Pauli terms from direct API
   - Print Pauli terms from Bond API
   - Check if they match

2. **Compare Ansatz Circuits:**
   - Print circuit gates from direct API
   - Print circuit gates from Bond API
   - Check if they match

3. **Compare Initial States:**
   - Check HF state preparation
   - Verify qubit mapping

---

## 💡 TEMPORARY WORKAROUND

**For Users:**
Use direct Hamiltonian API instead of Bond API for VQE:

```python
# AVOID (broken):
bond = BondFactory.create_bond('H', 'H')
solver = VQESolver(bond=bond, ...)

# USE (working):
from kanad.core.hamiltonians.openfermion_jw import openfermion_jordan_wigner
hamiltonian = openfermion_jordan_wigner(...)
ansatz = CovalentGovernanceAnsatz(...)
solver = VQESolver(hamiltonian=hamiltonian, ansatz=ansatz, ...)
```

---

## 📋 ACTION ITEMS

- [ ] Debug Bond.hamiltonian.to_sparse_hamiltonian()
- [ ] Compare with openfermion_jordan_wigner() output
- [ ] Fix Bond API Hamiltonian construction
- [ ] Add regression test comparing both APIs
- [ ] Update dashboard to use fixed Bond API

---

**This explains why your dashboard shows HF energy!**
The Bond API is broken for VQE.
