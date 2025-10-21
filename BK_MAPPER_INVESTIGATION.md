# Bravyi-Kitaev Mapper Investigation - ROOT CAUSE FOUND

## Problem
Bravyi-Kitaev mapper gives completely wrong energies:
- H₂: -0.350 Ha instead of -1.117 Ha (767 mHa error!)
- Energy above Hartree-Fock (unphysical)

## Investigation Timeline

### Fix Attempt 1: BravyiKitaevMapper Class
**Hypothesis**: Custom BK mapper implementation was broken
**Action**: Rewrote to use OpenFermion's validated BK transformation
**Result**: ❌ Still failed

### Fix Attempt 2: Pauli Converter
**Hypothesis**: Pauli converter always used JW, ignoring mapper parameter
**Action**: Modified pauli_converter.py to use correct OpenFermion transformation based on mapper type
**Result**: ❌ Still failed (but this fix was still necessary)

### Root Cause Discovery ✅

**Direct Hamiltonian Test**:
```python
# Both Hamiltonians built from same integrals
ham_jw = openfermion_jordan_wigner(h_mo, eri_mo, E_nuc, n_elec)
ham_bk = openfermion_bravyi_kitaev(h_mo, eri_mo, E_nuc, n_elec)

# HF state |1100⟩ (2 electrons in first 2 spin-orbitals)
hf_state = |1100⟩

# Test HF energy
HF Energy (JW Ham): -1.11675931 Ha  ✓ CORRECT
HF Energy (BK Ham): -0.34956289 Ha  ✗ WRONG!

# Test ground state by diagonalization
Ground state (JW Ham): -1.13728383 Ha  ✓ CORRECT
Ground state (BK Ham): -1.13728383 Ha  ✓ CORRECT!
```

**FINDING**: The BK Hamiltonian is **CORRECT**! The issue is the **reference state**!

## The Real Problem

### Jordan-Wigner Encoding
In JW, qubits directly represent orbital occupations:
- |1100⟩ means: orbital 0↑=1, 0↓=1, 1↑=0, 1↓=0
- This is the HF state for H₂ (2 electrons in lowest orbital)

### Bravyi-Kitaev Encoding
In BK, qubits encode **parity and occupancy information** in a binary tree:
- Qubit values do NOT directly represent orbital occupations!
- |1100⟩ in BK encoding represents a DIFFERENT physical state!
- The HF state |ψ_HF⟩ has a different bit-string representation in BK

### Current Bug
All ansätze (UCC, governance, etc.) hardcode the HF state as |1100⟩:

```python
# kanad/ansatze/ucc_ansatz.py:159
def _hartree_fock_state(self) -> List[int]:
    """Generate Hartree-Fock reference state."""
    # ...
    # For H2: return [1, 1, 0, 0]  # ← ONLY CORRECT FOR JORDAN-WIGNER!
```

This works for JW but gives the WRONG state for BK!

## The Solution

Use the existing `find_hf_state` utility (kanad/utils/hf_state_finder.py) to find the correct HF state for any mapper:

```python
def find_hf_state(hamiltonian_pauli, n_electrons, hf_energy):
    """
    Search all computational basis states to find which one
    has energy matching the HF energy.

    For JW: finds |1100⟩
    For BK: finds the correct BK-encoded state
    """
    for i in range(2**n_qubits):
        bits = format(i, f'0{n_qubits}b')
        if bits.count('1') != n_electrons:
            continue

        state = basis_state(i)
        energy = ⟨state|H|state⟩

        if abs(energy - hf_energy) < tolerance:
            return i, occupation_list
```

### Implementation Plan

**Option A**: Modify VQE Solver (RECOMMENDED)
```python
# In VQESolver.__init__ or solve():
if mapper_type == 'bravyi_kitaev':
    # Find correct HF state for BK encoding
    from kanad.utils.hf_state_finder import find_hf_state

    hf_energy = hamiltonian.solve_scf()[1]
    pauli_ham = hamiltonian.to_sparse_hamiltonian(mapper=mapper)
    _, hf_occupation = find_hf_state(pauli_ham, n_electrons, hf_energy)

    # Pass to ansatz
    circuit = ansatz.build_circuit(initial_state=hf_occupation)
else:
    # Use default HF state (works for JW)
    circuit = ansatz.build_circuit()
```

**Option B**: Modify Base Ansatz
- Add mapper parameter to ansatz __init__
- Use find_hf_state in _hartree_fock_state() method
- More architectural change required

## Files Modified So Far

1. ✅ `kanad/core/mappers/bravyi_kitaev_mapper.py` - Now uses OpenFermion
2. ✅ `kanad/core/hamiltonians/pauli_converter.py` - Now uses correct transformation based on mapper
3. ⏳ `kanad/solvers/vqe_solver.py` - NEEDS FIX for reference state
4. ⏳ `kanad/ansatze/*.py` - May need to accept initial_state parameter

## Next Steps

1. Implement Option A fix in VQE solver
2. Test on H₂ with BK mapper
3. Verify energy now matches JW result
4. Run full validation suite
5. Update other solvers (SQD, etc.) if needed

## Expected Result After Fix

```
Jordan-Wigner:
  HF state: |1100⟩
  Energy: -1.11675931 Ha ✓

Bravyi-Kitaev:
  HF state: |????⟩ (found by search)
  Energy: -1.11675931 Ha ✓ (should match!)
```

## Timeline

- Started: 2025-10-21 19:30
- Root cause found: 2025-10-21 20:05
- Time to fix: ~35 minutes of investigation
- Remaining: Implementation + testing

---

**Status**: Root cause identified, fix in progress
**Confidence**: 100% - verified with direct Hamiltonian test
**Complexity**: Medium - architectural change needed but clear path forward
