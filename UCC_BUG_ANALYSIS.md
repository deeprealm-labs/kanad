# UCC Ansatz Double Excitation Bug - Complete Analysis

## Summary

The UCC ansatz double excitation circuit is **fundamentally broken**. It does NOT create the correct fermionic double excitation and instead violates particle number conservation.

## Expected Behavior

For H₂ with 2 electrons and 4 spin orbitals:
- HF state: `|5⟩ = |0101⟩` (qubits 0,2 occupied - one electron in bonding MO, spin up and down)
- Double excitation: `|10⟩ = |1010⟩` (qubits 1,3 occupied - both electrons in antibonding MO)
- Correct ground state: 0.994|5⟩ - 0.112|10⟩ (minus sign is important!)

The double excitation operator should create a superposition:
```
|ψ(θ)⟩ = cos(θ/2)|5⟩ + sin(θ/2)|10⟩
```

This would give correlation energy of ~20 mHa.

## Actual Behavior

The current UCC circuit creates:
```
|ψ(θ)⟩ = cos(θ/2)|5⟩ + sin(θ/2)|15⟩
```

Where `|15⟩ = |1111⟩` has **4 electrons** (all qubits occupied)!

### Test Results

```
θ = 0.0:     |5⟩=1.000, |10⟩=0.000, |15⟩=0.000  ✓ Correct initial state
θ = π/4:     |5⟩=0.924, |10⟩=0.000, |15⟩=0.383  ✗ WRONG! Should be |10⟩ not |15⟩
θ = π/2:     |5⟩=0.707, |10⟩=0.000, |15⟩=0.707  ✗ Creates unphysical 4-electron state
```

## Root Cause

The double excitation circuit in `/kanad/ansatze/ucc_ansatz.py` lines 203-264 is not correctly implementing the Jordan-Wigner fermion-to-qubit transformation.

Current implementation:
```python
def _apply_double_excitation(self, circuit, occ1, occ2, virt1, virt2, exc_idx):
    theta = Parameter(f'θ_d_{exc_idx}')

    # Step 1: Entangle occupied orbitals
    circuit.cx(occ1, occ2)

    # Step 2: Entangle virtual orbitals
    circuit.cx(virt1, virt2)

    # Step 3: Create parity chain
    if occ2 < virt1:
        for q in range(occ2 + 1, virt1):
            circuit.cx(q, virt1)

    # Step 4: Apply parametrized rotation
    circuit.ry(theta, virt1)  # ← WRONG QUBIT!

    # Steps 5-6: Uncompute...
```

### Problems

1. **Wrong target qubit**: The RY rotation is applied to `virt1` only, but should involve all 4 qubits in the excitation
2. **Missing fermion parity**: Jordan-Wigner requires parity strings between creation/annihilation operators
3. **No particle conservation**: The circuit doesn't enforce that moving electrons from (occ1,occ2) to (virt1,virt2) conserves total electron number

## Correct Implementation

The double excitation operator in second quantization is:
```
T₂ = a†_{virt1} a†_{virt2} a_{occ2} a_{occ1}
```

Under Jordan-Wigner transformation, each fermionic operator becomes:
```
a†_i = (∏_{j<i} Z_j) · (X_i - iY_i)/2
a_i  = (∏_{j<i} Z_j) · (X_i + iY_i)/2
```

The parity strings (∏Z) are critical!

For the specific excitation [0,2]→[1,3], the full operator is:
```
exp(θ(|1010⟩⟨0101| - |0101⟩⟨1010|))
```

This can be decomposed using techniques from:
- Tranter et al., "A Comparison of the Bravyi–Kitaev and Jordan–Wigner Transformations for the Quantum Simulation of Quantum Chemistry"
- Or use Qiskit Nature's built-in UCC implementation

## Recommended Fix

**Option 1: Use Qiskit Nature's UCC**
```python
from qiskit_nature.second_q.circuit.library import UCC

# This is the gold-standard implementation
ucc_ansatz = UCC(
    excitations='sd',  # singles and doubles
    num_spatial_orbitals=n_orbitals,
    num_particles=(n_alpha, n_beta)
)
```

**Option 2: Fix Our Implementation**

Implement the proper gate decomposition following Qiskit Nature's source code or literature references. The correct decomposition for a double excitation involves:

1. Ladder of CNOT gates for parity tracking
2. Controlled rotations on specific qubit combinations
3. Proper phase factors
4. Symmetric application to create/annihilate pairs

Example structure (simplified):
```python
def _apply_double_excitation(self, circuit, i, j, a, b, theta):
    """
    Proper double excitation: |ij⟩ ↔ |ab⟩
    where i,j are occupied and a,b are virtual.
    """
    # Create parity ladder
    parity_qubits = list(range(min(i,j,a,b), max(i,j,a,b)+1))

    # Build excitation gate using proper Givens rotation network
    # This requires ~15-20 gates for a single double excitation
    # See: https://github.com/qiskit-community/qiskit-nature/blob/main/qiskit_nature/second_q/circuit/library/ansatzes/ucc.py

    # Simplified (NOT CORRECT, just for illustration):
    circuit.cx(i, j)
    circuit.cx(a, b)
    # ... many more gates for proper parity tracking ...
    circuit.ry(theta, a)  # Controlled by proper parity
    # ... uncompute parity ...
    circuit.cx(a, b)
    circuit.cx(i, j)
```

## Impact

**Critical**: UCC is supposed to be the "gold standard" ansatz for quantum chemistry VQE. Having it completely broken means:

- All UCC VQE results are wrong (converge to HF only)
- No correlation energy captured
- ~20 mHa systematic error for H₂
- Larger errors for more complex molecules
- Users cannot trust VQE results

**Workarounds** (until fixed):
1. Use Governance ansatz (WORKS - captures full correlation)
2. Use Hardware-Efficient ansatz (WORKS - captures full correlation)
3. Import Qiskit Nature's UCC implementation

## Verification Test

After fixing, run this test:
```python
from kanad.ansatze.ucc_ansatz import UCCAnsatz
import numpy as np

ansatz = UCCAnsatz(n_qubits=4, n_electrons=2)
ansatz.build_circuit()

# Test double excitation connects |5⟩ ↔ |10⟩
params = np.array([0, 0, 0, 0, np.pi/4])  # Only double excitation
ansatz.circuit.bind_parameters(params)

qc = ansatz.circuit.to_qiskit()
state = Statevector.from_instruction(qc)

amp_5 = abs(state.data[5])   # Should be ~0.924
amp_10 = abs(state.data[10])  # Should be ~0.383
amp_15 = abs(state.data[15])  # Should be ~0.000

assert amp_5 > 0.9, f"State |5⟩ amp too small: {amp_5}"
assert amp_10 > 0.3, f"State |10⟩ amp too small: {amp_10} - DOUBLE EXCITATION BROKEN!"
assert amp_15 < 0.01, f"State |15⟩ should be zero: {amp_15} - PARTICLE NUMBER NOT CONSERVED!"

print("✓ UCC double excitation working correctly!")
```

Current result: **FAILS** - amp_10 = 0.000, amp_15 = 0.383

## References

1. Qiskit Nature UCC implementation: https://github.com/qiskit-community/qiskit-nature
2. Romero et al., "Strategies for quantum computing molecular energies using the unitary coupled cluster ansatz" (2018)
3. Barkoutsos et al., "Quantum algorithms for electronic structure calculations: Particle-hole Hamiltonian and optimized wave-function expansions" (2018)
