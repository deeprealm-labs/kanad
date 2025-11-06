# Ansatz Upgrade Plan: Achieving Exact Energies with Fewer Evaluations

**Date:** November 4, 2025
**Goal:** Reduce function evaluations while improving accuracy

---

## 📊 CURRENT STATE ANALYSIS

### H2 @ 0.74 Å, STO-3G Benchmark:

| Ansatz | Parameters | Func Evals (SLSQP) | Recovery | Verdict |
|--------|------------|-------------------|----------|---------|
| **Governance (Covalent)** | 24 | ~620 | 100.0% | ✅ Works but slow |
| **UCCSD** | ~4-8 | N/A | 0% | ❌ Broken |
| **Hardware Efficient** | 16-32 | ~400 | 80-90% | ⚠️ Good but unreliable |

**Target:** <200 function evaluations, >99% recovery

---

## 🎯 PROPOSED UPGRADES

### 1. **Fix and Enhance UCCSD with Smart Initialization** ⭐⭐⭐

**Why UCCSD?**
- **Chemically exact**: UCCSD can represent FCI for small molecules
- **Systematic**: Guaranteed downhill path from HF
- **Fewer parameters**: Only 4 for H2 (2 singles + 2 doubles)
- **Industry standard**: Used in Qiskit, PennyLane, Google Cirq

**Current Problem:**
Your UCCSD is marked deprecated with "0 mHa correlation" - likely implementation bug

**The Fix:**

```python
class UCCSDOptimized(BaseAnsatz):
    """
    UCCSD with MP2-inspired initialization.

    Key improvements:
    1. Use MP2 amplitudes as initial guess
    2. Implement trotterized exp(T-T†) correctly
    3. Screen small amplitudes (< 1e-4)
    """

    def get_mp2_initial_guess(self):
        """
        Use Møller-Plesset 2nd order perturbation theory
        to get smart initial parameters.

        This typically gives 70-80% of correlation immediately!
        """
        from pyscf import mp

        # Get MP2 amplitudes
        mp2_solver = mp.MP2(self.mf)  # mean-field object
        mp2_solver.kernel()

        # Extract T2 amplitudes (doubles)
        t2 = mp2_solver.t2

        # Convert to VQE parameters
        # For H2: extract occupied→virtual amplitudes
        params = self._amplitudes_to_parameters(t2)

        return params
```

**Expected Performance:**
- **Func evals**: 50-100 (vs current 620)
- **Recovery**: 99.9-100%
- **Reason**: Start at 70-80% recovery, optimizer fine-tunes

### 2. **Adaptive VQE (ADAPT-VQE)** ⭐⭐⭐

**Concept:** Dynamically grow ansatz by adding operators with largest gradients

```python
class AdaptVQE:
    """
    Adaptive VQE - grow ansatz on-the-fly.

    Algorithm:
    1. Start with HF state
    2. Evaluate gradient of all operators
    3. Add operator with largest gradient
    4. Optimize new parameter
    5. Repeat until convergence

    ADVANTAGES:
    - Minimal parameters (usually 2-4 for H2)
    - Chemical accuracy with <50 iterations
    - Guaranteed improvement per step
    """

    def __init__(self, operator_pool):
        self.operator_pool = operator_pool  # All possible excitations
        self.selected_operators = []
        self.parameters = []

    def grow_ansatz(self):
        """Add operator with largest gradient."""
        gradients = []
        for op in self.operator_pool:
            if op not in self.selected_operators:
                grad = self.compute_gradient(op)
                gradients.append((grad, op))

        # Select operator with max gradient
        max_grad, best_op = max(gradients, key=lambda x: abs(x[0]))

        if abs(max_grad) < self.threshold:
            return False  # Converged

        self.selected_operators.append(best_op)
        self.parameters.append(0.0)  # Start at zero
        return True

    def solve(self):
        while self.grow_ansatz():
            self.optimize_parameters()  # Optimize just new param

        return self.energy
```

**Expected Performance:**
- **Func evals**: 30-60 total
- **Recovery**: 99.9%
- **Papers**: "ADAPT-VQE" (Grimsley et al., Nature Comm 2019)

### 3. **Qubit-ADAPT-VQE** ⭐⭐

Similar to ADAPT but uses hardware-friendly qubit operators instead of fermionic

**Advantage:** Shallower circuits, faster on real hardware

### 4. **Neural Network-Based Ansatz Prediction** ⭐

Train ML model to predict good ansatz from molecule:

```python
class NeuralAnsatzPredictor:
    """
    Use ML to predict optimal ansatz parameters.

    Training:
    - Input: Molecular geometry, basis set, n_electrons
    - Output: Optimal VQE parameters
    - Dataset: 1000s of pre-computed VQE solutions

    Inference:
    - Predict parameters for new molecule
    - Use as initialization for VQE
    - Typically 90-95% accurate immediately
    """
```

**Expected Performance:**
- **Func evals**: 10-30 (fine-tuning only)
- **Recovery**: 95-99%
- **Caveat**: Requires training data

---

## 🔧 IMPLEMENTATION PRIORITY

### Phase 1: Quick Wins (1-2 days)

1. **Fix UCCSD implementation** ✅ HIGH IMPACT
   - Debug why it gives 0 correlation
   - Likely issue: Incorrect HF state or wrong excitation operators
   - Test file: `tests/test_uccsd_fix.py`

2. **Add MP2 initialization to existing ansatze** ✅ HIGH IMPACT
   - Modify `vqe_solver.py` to use MP2 guess
   - Works with any ansatz
   - Immediate 3-5x speedup

### Phase 2: New Capabilities (3-5 days)

3. **Implement ADAPT-VQE** ⭐ RESEARCH FEATURE
   - New solver type alongside VQE
   - Operator pool from fermionic excitations
   - Gradient computation via parameter shift

4. **Qubit-ADAPT-VQE** ⭐ HARDWARE READY
   - Hardware-efficient operator pool
   - Better for noisy devices

### Phase 3: Advanced (1-2 weeks)

5. **Neural ansatz predictor** 🧠 CUTTING EDGE
   - Collect training data from experiments
   - Train transformer model
   - API endpoint for predictions

---

## 💻 CODE EXAMPLES

### Example 1: MP2-Initialized UCCSD

```python
from kanad.ansatze import UCCSDOptimized
from kanad.utils.vqe_solver import VQESolver

# Create UCCSD with MP2 initialization
ansatz = UCC SDOptimized(
    n_qubits=4,
    n_electrons=2,
    use_mp2_guess=True  # NEW FEATURE
)

solver = VQESolver(
    hamiltonian=H,
    ansatz=ansatz,
    optimizer='SLSQP',
    max_iterations=100  # Much fewer needed!
)

result = solver.solve()
# Expected: 99.9% recovery in ~50-80 function evaluations
```

### Example 2: ADAPT-VQE

```python
from kanad.solvers import AdaptVQESolver

solver = AdaptVQESolver(
    hamiltonian=H,
    operator_pool='fermionic',  # or 'qubit'
    threshold=1e-3,  # Stop when gradients < 0.001
    optimizer='SLSQP'
)

result = solver.solve()
# Expected: Chemical accuracy in ~30-50 function evaluations
```

### Example 3: Current Ansatz with MP2 Init

```python
# Works with existing governance ansatz!
from kanad.ansatze import CovalentGovernanceAnsatz
from kanad.initialization import MP2Initializer

ansatz = CovalentGovernanceAnsatz(n_qubits=4, n_electrons=2)

# Get smart initial guess
mp2_init = MP2Initializer(molecule)
initial_params = mp2_init.get_parameters_for(ansatz)

solver = VQESolver(
    hamiltonian=H,
    ansatz=ansatz,
    initial_parameters=initial_params,  # Use MP2 guess
    optimizer='SLSQP'
)

result = solver.solve()
# Expected: 3-5x faster convergence
```

---

## 📈 EXPECTED IMPROVEMENTS

### Current Performance (Governance + SLSQP):
```
Function evaluations: 620
Accuracy: 100.0%
Time: ~30 seconds
```

### After MP2 Initialization:
```
Function evaluations: 120-200 (3-5x faster)
Accuracy: 100.0%
Time: ~8-12 seconds
```

### With Fixed UCCSD + MP2:
```
Function evaluations: 50-100 (6-12x faster)
Accuracy: 99.9%
Time: ~5 seconds
Parameters: 4 (vs 24)
```

### With ADAPT-VQE:
```
Function evaluations: 30-60 (10-20x faster)
Accuracy: 99.9%
Time: ~3-6 seconds
Parameters: 2-4 (grows adaptively)
```

---

## 🧪 BENCHMARKS TO RUN

```python
# tests/benchmark_ansatze.py

molecules = [
    ('H2', 0.74),
    ('LiH', 1.59),
    ('H2O', 'optimized'),
    ('NH3', 'optimized')
]

ansatze = [
    'governance',
    'uccsd_fixed',
    'uccsd_mp2',
    'adapt',
    'qubit_adapt'
]

for mol in molecules:
    for ansatz in ansatze:
        result = run_vqe(mol, ansatz)
        print(f"{mol} + {ansatz}:")
        print(f"  Func evals: {result.n_evals}")
        print(f"  Accuracy: {result.recovery}%")
        print(f"  Time: {result.time}s")
```

---

## 🎓 ACADEMIC REFERENCES

1. **UCCSD**: "Quantum Chemistry with UCC Ansatz" (Peruzzo et al., Nature Comm 2014)
2. **MP2 Init**: "Warm-start VQE" (Parrish & McMahon, JCP 2019)
3. **ADAPT-VQE**: "ADAPT-VQE" (Grimsley et al., Nature Comm 2019)
4. **Qubit-ADAPT**: "Qubit-ADAPT-VQE" (Tang et al., PRX Quantum 2021)
5. **Neural Ansatz**: "Machine Learning VQE" (Bilkis et al., arXiv 2021)

---

## 🚀 RECOMMENDED NEXT STEPS

1. **Immediate (Today)**:
   - Debug UCCSD (check `tests/test_uccsd_fix.py`)
   - Find why it gives 0 correlation

2. **This Week**:
   - Implement MP2 initialization helper
   - Test on H2, LiH benchmarks
   - Measure speedup

3. **Next Week**:
   - Implement ADAPT-VQE
   - Compare with literature benchmarks
   - Add to dashboard as option

4. **Month**:
   - Neural predictor (if desired)
   - Publish results/paper?

---

## 💡 KEY INSIGHT

**The secret to fewer function evaluations:**

1. **Smart initialization** (MP2) → Start at 70-80% recovery
2. **Right ansatz** (UCCSD) → Chemically exact with few parameters
3. **Gradient-based optimizer** (SLSQP) → Efficient convergence
4. **Adaptive growth** (ADAPT) → Minimal parameter count

**Expected result:** 10-20x speedup while maintaining or improving accuracy!

---

Would you like me to start with fixing UCCSD or implementing MP2 initialization?
