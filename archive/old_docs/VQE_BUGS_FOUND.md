# VQE Optimization Bugs - ROOT CAUSE FOUND! 🔍

**Date:** November 4, 2025
**Status:** CRITICAL BUGS IDENTIFIED

---

## 🔴 BUG #1: VQE Solver Expects Wrong Hamiltonian Type

### The Problem:
VQE solver crashes when given a `SparsePauliOp` instead of `MolecularHamiltonian`

### Error:
```python
AttributeError: 'SparsePauliOp' object has no attribute 'n_orbitals'
```

### Location:
[kanad/utils/vqe_solver.py:558](kanad/utils/vqe_solver.py#L558)
```python
def _compute_energy_statevector(self, parameters):
    ...
    full_n_qubits = 2 * self.hamiltonian.n_orbitals  # ← CRASHES HERE!
```

### Root Cause:
- VQE solver's low-level API accepts `hamiltonian` parameter
- Documentation says it accepts "Hamiltonian object"
- Code assumes it's a `MolecularHamiltonian` with `.n_orbitals` attribute
- When tests pass `SparsePauliOp` directly, it crashes

### Why This Breaks Tests:
```python
# Test creates SparsePauliOp
pauli_hamiltonian = openfermion_jordan_wigner(...)  # Returns SparsePauliOp

# But VQE expects MolecularHamiltonian
solver = VQESolver(
    hamiltonian=pauli_hamiltonian,  # ← Wrong type!
    ansatz=ansatz,
    ...
)
```

### Fix Required:
1. **Option A:** Make VQE solver accept both types
   ```python
   if isinstance(self.hamiltonian, SparsePauliOp):
       # Use sparse method
       full_n_qubits = self.circuit.num_qubits
   else:
       # Use molecular hamiltonian method
       full_n_qubits = 2 * self.hamiltonian.n_orbitals
   ```

2. **Option B:** Force tests to use MolecularHamiltonian
   ```python
   # Create proper Molecular Hamiltonian object
   ham = MolecularHamiltonian(mol)
   solver = VQESolver(hamiltonian=ham, ...)
   ```

---

## 🔴 BUG #2: TwoLocalAnsatz Wrong Parameter Name

### The Problem:
TwoLocalAnsatz doesn't accept `reps` parameter

### Error:
```python
TypeError: TwoLocalAnsatz.__init__() got an unexpected keyword argument 'reps'
```

### Root Cause:
Test uses:
```python
ansatz = TwoLocalAnsatz(n_qubits=n_qubits, reps=2)  # ← reps not supported
```

But TwoLocalAnsatz expects different parameters

### Fix Required:
Check TwoLocalAnsatz signature and use correct parameter name

---

## 🟡 ISSUE #3: UCC Ansatz Still Used (Despite Deprecation)

### The Problem:
UCC ansatz is deprecated but still used in tests

### Warning:
```
DeprecationWarning: UCCAnsatz is deprecated and will be removed in v2.0.
Use CovalentGovernanceAnsatz (49x better on ionic molecules),
IonicGovernanceAnsatz, or HardwareEfficientAnsatz instead.
UCC family shows 0 mHa correlation in our implementation.
```

### Impact:
- Tests using UCC will always fail (0 correlation energy)
- Should use governance ansatze instead

---

## ✅ GOOD NEWS

### What Works:
1. ✅ Hamiltonian generation (`openfermion_jordan_wigner`) - WORKS
2. ✅ Ansatz creation (CovalentGovernance, TwoLocal) - WORKS
3. ✅ VQE solver initialization - WORKS
4. ✅ HF reference energy calculation - WORKS
5. ✅ Circuit binding and depth calculation - WORKS

### The Real Problem:
**VQE isn't broken** - the tests are using it incorrectly!

---

## 🔧 FIX PLAN

### Priority 1: Fix VQE Solver Type Handling
```python
# File: kanad/utils/vqe_solver.py

def _compute_energy_statevector(self, parameters):
    """
    Compute energy using statevector simulation.
    """
    # BEFORE (broken):
    full_n_qubits = 2 * self.hamiltonian.n_orbitals

    # AFTER (fixed):
    if isinstance(self.hamiltonian, SparsePauliOp):
        # Using sparse Pauli operator directly
        full_n_qubits = self.circuit.num_qubits
    elif hasattr(self.hamiltonian, 'n_orbitals'):
        # Using molecular Hamiltonian
        full_n_qubits = 2 * self.hamiltonian.n_orbitals
    else:
        raise TypeError(f"Unsupported Hamiltonian type: {type(self.hamiltonian)}")
```

### Priority 2: Fix TwoLocalAnsatz Usage
```python
# Check signature
from kanad.ansatze.two_local_ansatz import TwoLocalAnsatz
help(TwoLocalAnsatz.__init__)

# Use correct parameter name
ansatz = TwoLocalAnsatz(n_qubits=n_qubits, n_layers=2)  # or whatever it expects
```

### Priority 3: Update Tests to Use MolecularHamiltonian
```python
# RECOMMENDED APPROACH: Use high-level API
from api.services.experiment_service import create_molecule_from_config

molecule_config = {
    "smiles": "[H][H]",
    "basis": "sto-3g",
    "charge": 0,
    "multiplicity": 1
}

molecule = create_molecule_from_config(molecule_config)
ham = molecule.hamiltonian  # This is a MolecularHamiltonian

# Now VQE will work
solver = VQESolver(
    hamiltonian=ham,  # ← Correct type!
    ansatz=ansatz,
    optimizer='COBYLA',
    max_iterations=50
)
```

---

## 📊 TEST RESULTS

### Test Output Analysis:
```
✅ Hamiltonian generation: WORKS (15 Pauli terms)
✅ HF energy calculation: WORKS (-1.11675931 Ha)
✅ Nuclear repulsion: WORKS (0.71510434 Ha)
✅ Ansatz creation: WORKS (24 parameters for CovalentGovernance)
✅ Circuit binding: WORKS (4 qubits, depth 11)
❌ VQE optimization: CRASHES (wrong Hamiltonian type)
```

### Performance Stats:
- Hamiltonian build: ~0.1s
- HF calculation: ~0.05s
- Circuit construction: ~0.01s
- **VQE crash: immediate** (type error on first evaluation)

---

## 🎯 NEXT STEPS

1. **Implement Fix #1** (VQE type handling) - 30 min
2. **Fix TwoLocalAnsatz** parameter - 10 min
3. **Rerun tests** - 5 min
4. **If tests pass:** Write proper test suite - 1 hour
5. **If tests still fail:** Debug further

---

## 💡 KEY INSIGHT

**VQE optimization wasn't converging at iteration 0 because it was CRASHING, not converging!**

The original test output showed:
- Iterations: 0
- Correlation: 0.00 Ha

This looked like convergence failure, but it was actually:
1. VQE tries to evaluate energy
2. Crashes with AttributeError
3. Returns default/initial state
4. Reports 0 iterations

**Once we fix the type handling, VQE should work correctly!**

---

## 📝 DOCUMENTATION UPDATES NEEDED

After fixing:
1. Update VQE solver docstring to clarify accepted Hamiltonian types
2. Add type hints: `hamiltonian: Union[MolecularHamiltonian, SparsePauliOp]`
3. Add examples for both usage modes
4. Update test documentation

---

## ✅ CONFIDENCE LEVEL

**High (90%)**  - The bug is clearly identified and fixable.

**Why high confidence:**
- Error message is clear
- Root cause is obvious (type mismatch)
- Fix is straightforward
- Test setup is otherwise correct

**Expected outcome after fix:**
- VQE will run for 20-50 iterations
- Will recover correlation energy (-0.04 to -0.06 Ha for H2)
- All ansatze should work
- All optimizers should work

---

**Ready to implement the fix! 🚀**
