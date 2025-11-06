# VQE Governance Issue - CONFIRMED

**Date:** November 6, 2025
**Status:** ISSUE FOUND - VQE ansätze do NOT validate excitations

---

## Investigation Results

### Issue #1: VQE Expectation - OK ✅

**Status:** No issue found - uses proper SparsePauliOp expectation

**Evidence:**
```python
# kanad/solvers/vqe_solver.py:541, 544
energy = statevector.expectation_value(self._sparse_pauli_op).real

# Line 1082
H_expectation = np.real(np.conj(psi) @ self._hamiltonian_matrix @ psi)
```

**Conclusion:** VQE correctly computes `<ψ(θ)|H|ψ(θ)>` using SparsePauliOp. HF energy is only used for reference (correlation energy calculation), not as placeholder.

**Action:** NONE NEEDED

---

## Issue #2: VQE Governance - NOT ENFORCED ❌

**Status:** CRITICAL - Ansätze do NOT call `is_valid_configuration()`

**Evidence:**

### Search Results:
```bash
grep -r "is_valid_configuration" kanad/ansatze/
# Result: NO MATCHES FOUND!
```

**Problem:** 
- SQD uses governance validation (lines 177-183, 240-246 in sqd_solver.py)
- VQE ansätze have governance protocols but DO NOT validate excitations
- UCC ansatz generates ALL excitations without checking validity

### UCC Ansatz (kanad/ansatze/ucc_ansatz.py:293-303)
```python
def get_excitation_list(self) -> List[Tuple[List[int], List[int]]]:
    """Get list of excitations."""
    return self.excitations  # NO VALIDATION!
```

### Governance-aware ansatz (kanad/ansatze/governance_aware_ansatz.py)
- Has protocols imported
- Line 335: Comment "Check governance (would use protocol.validate_orbital_hybridization)"
- BUT: NO ACTUAL VALIDATION IMPLEMENTED

---

## Root Cause

**VQE ansätze generate excitations at construction time, but never filter them using governance protocols.**

**Why this matters:**
1. SQD explicitly filters: `if governance_protocol.is_valid_configuration(bitstring)`
2. VQE never calls this validation
3. VQE may include unphysical excitations that SQD rejects
4. Inconsistent treatment between solvers

---

## Fix Plan

### Option 1: Add validation to UCC ansatz (Recommended)

Modify `UCCAnsatz.__init__()` to accept governance protocol and filter excitations:

```python
class UCCAnsatz(BaseAnsatz):
    def __init__(
        self,
        n_qubits: int,
        n_electrons: int,
        include_singles: bool = True,
        include_doubles: bool = True,
        governance_protocol: Optional[Any] = None  # ADD THIS
    ):
        super().__init__(n_qubits, n_electrons)
        self.include_singles = include_singles
        self.include_doubles = include_doubles
        self.governance_protocol = governance_protocol  # Store protocol
        
        # Generate all possible excitations
        all_excitations = []
        
        n_orb = n_qubits // 2
        n_occ = n_electrons // 2
        
        if include_singles:
            for i in range(n_occ):
                for a in range(n_occ, n_orb):
                    all_excitations.append(([i], [a]))
        
        if include_doubles:
            for i in range(n_occ):
                for j in range(i + 1, n_occ):
                    for a in range(n_occ, n_orb):
                        for b in range(a + 1, n_orb):
                            all_excitations.append(([i, j], [a, b]))
        
        # CRITICAL FIX: Filter excitations using governance protocol
        if governance_protocol is not None:
            valid_excitations = []
            for exc in all_excitations:
                # Convert excitation to bitstring
                bitstring = self._excitation_to_bitstring(exc, n_electrons, n_qubits)
                
                # Check if valid according to governance
                if governance_protocol.is_valid_configuration(bitstring):
                    valid_excitations.append(exc)
            
            logger.info(f"Governance filtering: {len(all_excitations)} → {len(valid_excitations)} excitations")
            self.excitations = valid_excitations
        else:
            self.excitations = all_excitations
    
    def _excitation_to_bitstring(self, excitation, n_electrons, n_qubits):
        """Convert excitation indices to bitstring configuration."""
        # Start with HF reference
        bitstring = '1' * n_electrons + '0' * (n_qubits - n_electrons)
        bitlist = list(bitstring)
        
        occ, virt = excitation
        # Apply excitation (spin-up and spin-down)
        for i in occ:
            bitlist[i] = '0'  # Remove electron from orbital i
            bitlist[i + n_qubits//2] = '0'  # Both spins
        for a in virt:
            bitlist[a] = '1'  # Add electron to orbital a
            bitlist[a + n_qubits//2] = '1'  # Both spins
        
        return ''.join(bitlist)
```

### Option 2: Add validation to VQE solver

Modify VQE solver to filter ansatz excitations after construction:

```python
# In VQESolver._init_ansatz() after creating ansatz
if hasattr(self.hamiltonian, 'governance_protocol') and self.hamiltonian.governance_protocol:
    # Filter ansatz excitations
    if hasattr(self.ansatz, 'excitations'):
        protocol = self.hamiltonian.governance_protocol
        valid_exc = []
        for exc in self.ansatz.excitations:
            bitstring = self._excitation_to_bitstring(exc)
            if protocol.is_valid_configuration(bitstring):
                valid_exc.append(exc)
        self.ansatz.excitations = valid_exc
        logger.info(f"Governance filtered: {len(valid_exc)} valid excitations")
```

---

## Priority

**HIGH** - This is an inconsistency between SQD and VQE that could lead to:
1. Different results from same Hamiltonian
2. VQE including unphysical states
3. Incorrect comparison between solvers

---

## Action Required

1. Implement excitation validation in UCC ansatz
2. Pass governance protocol from VQE solver to ansatz
3. Test VQE vs SQD consistency
4. Validate that both solvers use same excitations

