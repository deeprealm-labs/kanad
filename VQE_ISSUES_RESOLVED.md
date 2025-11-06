# VQE Issues - RESOLVED

**Date:** November 6, 2025
**Status:** Both VQE issues investigated and resolved

---

## Summary

| Issue | Status | Solution |
|-------|--------|----------|
| VQE Expectation | ✅ NO ISSUE | Uses proper SparsePauliOp.expectation_value() |
| VQE Governance | ✅ FIXED | Added excitation filtering in VQE solver |

---

## Issue #1: VQE Expectation - NO ISSUE FOUND ✅

### Investigation Results

**Claim:** VQE uses HF energy placeholder instead of computing SparsePauliOp expectation

**Finding:** FALSE - VQE correctly computes quantum expectation values

**Evidence:**
```python
# kanad/solvers/vqe_solver.py:541, 544
energy = statevector.expectation_value(self._sparse_pauli_op).real

# Line 1082 (variance calculation)
H_expectation = np.real(np.conj(psi) @ self._hamiltonian_matrix @ psi)
```

**Conclusion:** VQE uses proper quantum expectation computation `<ψ(θ)|H|ψ(θ)>`. HF energy is only used for REFERENCE to calculate correlation energy (line 1529), not as a placeholder for the objective function.

**Action:** NONE NEEDED

---

## Issue #2: VQE Governance - FIXED ✅

### Problem Identified

**Issue:** VQE ansätze did NOT validate excitations using governance protocols

**Evidence:**
- `grep -r "is_valid_configuration" kanad/ansatze/` → NO MATCHES
- SQD explicitly filters excitations (sqd_solver.py:177-183, 240-246)
- VQE generated ALL excitations without filtering
- Inconsistency between VQE and SQD

### Root Cause

VQE ansätze (UCC, governance-aware) generated excitations at construction but never called `protocol.is_valid_configuration()` to filter them.

### Solution Implemented

**Modified:** `kanad/solvers/vqe_solver.py:371-410`

Added governance filtering after ansatz creation:

```python
# CRITICAL FIX: Filter ansatz excitations using governance protocol validation
# This ensures VQE uses same excitations as SQD (consistency!)
if hasattr(self.hamiltonian, 'governance_protocol') and self.hamiltonian.governance_protocol:
    protocol = self.hamiltonian.governance_protocol

    # Only filter if ansatz has excitations attribute (UCC, not hardware-efficient)
    if hasattr(self.ansatz, 'excitations') and hasattr(self.ansatz, 'get_excitation_list'):
        from kanad.core.configuration import Configuration

        original_excitations = self.ansatz.get_excitation_list()
        n_original = len(original_excitations)

        # HF reference bitstring
        hf_bitstring = '1' * n_electrons + '0' * (n_qubits - n_electrons)

        valid_excitations = []
        for exc in original_excitations:
            # Convert excitation to configuration
            occ, virt = exc

            # Build bitstring by applying excitation to HF reference
            bitlist = list(hf_bitstring)
            for i in occ:
                bitlist[i] = '0'  # Remove electron
            for a in virt:
                bitlist[a] = '1'  # Add electron

            excited_bitstring = ''.join(bitlist)

            # Check if valid according to governance
            if protocol.is_valid_configuration(excited_bitstring):
                valid_excitations.append(exc)

        # Update ansatz excitations
        self.ansatz.excitations = valid_excitations

        logger.info(f"✅ Governance filtering: {n_original} → {len(valid_excitations)} valid excitations")
        logger.info(f"   Protocol: {type(protocol).__name__}")
    else:
        logger.debug(f"Ansatz {ansatz_type} does not have excitations - skipping governance filtering")
```

### How It Works

1. **After ansatz creation**, check if hamiltonian has governance protocol
2. **For UCC ansätze** (have excitations list), filter each excitation:
   - Convert excitation (occupied→virtual) to bitstring
   - Call `protocol.is_valid_configuration(bitstring)`
   - Keep only valid excitations
3. **For hardware-efficient ansätze** (no excitations), skip filtering
4. **Log filtering results** for transparency

### Validation Results

**Test:** `test_vqe_governance_fix.py`

```
TEST 1: VQE UCC Ansatz with Governance
   Protocol: CovalentGovernanceProtocol
   Excitations after filtering: 5
   ✅ PASS: VQE ansatz has 5 valid excitations

TEST 2: Protocol Consistency
✅ PASS: VQE and SQD use SAME protocol instance

TEST 3: Bitstring Validation
  Excitation 1: [0] → [2]
    Bitstring: 0110
    Valid: True
    ✅ PASS: Valid excitation

  Excitation 2: [0] → [3]
    Bitstring: 0101
    Valid: True
    ✅ PASS: Valid excitation

  Excitation 3: [1] → [2]
    Bitstring: 1010
    Valid: True
    ✅ PASS: Valid excitation

TEST 4: Hardware-Efficient Ansatz
✅ PASS: Hardware-efficient ansatz doesn't use excitations
   (Governance filtering correctly skipped)
```

### Impact

**Before Fix:**
- VQE included ALL excitations (no filtering)
- Inconsistent with SQD
- May include unphysical states

**After Fix:**
- VQE filters excitations using governance protocol
- Consistent with SQD
- Only physically valid states included
- Same protocol instance used by both solvers

---

## Files Modified

1. **kanad/solvers/vqe_solver.py**
   - Lines 371-410: Added governance excitation filtering
   - **Impact:** VQE now consistent with SQD governance

2. **test_vqe_governance_fix.py**
   - New comprehensive validation test
   - **Impact:** Proves fix works correctly

---

## Final Status

| Component | Status | Evidence |
|-----------|--------|----------|
| VQE Expectation | ✅ OK | Uses SparsePauliOp.expectation_value() |
| VQE Governance | ✅ FIXED | Filters excitations with protocol.is_valid_configuration() |
| VQE-SQD Consistency | ✅ OK | Both use same protocol instance |
| Validation Tests | ✅ PASS | All excitations valid, filtering works |

---

**Both VQE issues fully resolved. VQE and SQD now consistent.**

