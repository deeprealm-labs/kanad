# Phase 1: Fix Strategy - Strict Scientific Approach

## Problem Analysis

The VQESolver was refactored from a **component-based API** to a **bond-based API**:

### Old API (Tests Use This)
```python
solver = VQESolver(
    hamiltonian=hamiltonian,  # Direct Hamiltonian
    ansatz=ansatz,            # Direct ansatz object
    mapper=mapper,            # Direct mapper object
    optimizer='COBYLA'
)
```

### New API (Current Implementation)
```python
bond = BondFactory.create_bond('H', 'H')
solver = VQESolver(
    bond=bond,                # Bond object
    ansatz_type='ucc',        # String identifier
    mapper_type='jordan_wigner',
    optimizer_method='COBYLA'
)
```

## Decision: Support Both APIs ✅

**Rationale:**
1. **Unit Testing** requires direct component injection
2. **User API** benefits from high-level bond interface
3. **Backward Compatibility** maintains existing code
4. **Scientific Rigor** allows testing specific configurations

## Implementation Strategy

### Option 1: Dual-Mode VQESolver (RECOMMENDED)

Make VQESolver detect which API is being used and adapt:

```python
class VQESolver:
    def __init__(
        self,
        bond: Optional['BaseBond'] = None,
        # New API parameters
        ansatz_type: str = 'ucc',
        mapper_type: str = 'jordan_wigner',
        # Old API parameters (for testing)
        hamiltonian: Optional['Hamiltonian'] = None,
        ansatz: Optional['BaseAnsatz'] = None,
        mapper: Optional['BaseMapper'] = None,
        # Common parameters
        optimizer: str = 'SLSQP',
        max_iterations: int = 1000,
        **kwargs
    ):
        # Detect which API is being used
        if bond is not None:
            # New high-level API
            self._init_from_bond(bond, ansatz_type, mapper_type, **kwargs)
        elif hamiltonian is not None and ansatz is not None:
            # Old low-level API (for testing)
            self._init_from_components(hamiltonian, ansatz, mapper, **kwargs)
        else:
            raise ValueError("Must provide either 'bond' or 'hamiltonian+ansatz'")
```

**Advantages:**
- ✅ Backward compatible
- ✅ Tests don't need modification
- ✅ Supports both use cases
- ✅ Clean separation of concerns

**Disadvantages:**
- More complex __init__
- Two code paths to maintain

### Option 2: Separate Classes

Create two separate classes:

```python
class VQESolver:  # High-level, user-facing
    def __init__(self, bond, ansatz_type='ucc', ...):
        ...

class VQESolverCore:  # Low-level, for testing
    def __init__(self, hamiltonian, ansatz, mapper, ...):
        ...
```

**Advantages:**
- Clean separation
- Each class has single responsibility

**Disadvantages:**
- Code duplication
- Need to maintain both
- Tests need updating anyway

### Option 3: Update All Tests

Update all tests to use new API:

```python
# Convert this:
solver = VQESolver(hamiltonian=ham, ansatz=ans, mapper=map)

# To this:
from kanad.bonds import BondFactory
bond = BondFactory.create_bond('H', 'H')
solver = VQESolver(bond)
```

**Advantages:**
- Single API
- Cleaner codebase

**Disadvantages:**
- Loses fine-grained testing control
- Can't test specific component combinations
- Not suitable for unit testing

## DECISION: Option 1 - Dual-Mode VQESolver ✅

**Why:**
1. **Scientific Testing Requires Precision** - Need to test exact components
2. **Backward Compatibility** - Don't break existing code
3. **Best of Both Worlds** - High-level for users, low-level for testing

## Implementation Plan

### Step 1: Modify VQESolver __init__

```python
class VQESolver(BaseSolver):
    def __init__(
        self,
        bond: Optional['BaseBond'] = None,
        # High-level API (bond-based)
        ansatz_type: str = 'ucc',
        mapper_type: str = 'jordan_wigner',
        # Low-level API (component-based, for testing)
        hamiltonian: Optional['Hamiltonian'] = None,
        ansatz: Optional['BaseAnsatz'] = None,
        mapper: Optional['BaseMapper'] = None,
        # Common parameters
        optimizer: str = 'SLSQP',
        max_iterations: int = 1000,
        backend: str = 'statevector',
        **kwargs
    ):
        # Auto-detect API mode
        if bond is not None:
            # High-level API
            if hamiltonian is not None or ansatz is not None:
                raise ValueError("Cannot mix bond-based and component-based API")
            self._init_from_bond(bond, ansatz_type, mapper_type, **kwargs)
        elif hamiltonian is not None:
            # Low-level API (for testing)
            if bond is not None:
                raise ValueError("Cannot mix bond-based and component-based API")
            if ansatz is None:
                raise ValueError("ansatz required when using component-based API")
            self._init_from_components(hamiltonian, ansatz, mapper, **kwargs)
        else:
            raise ValueError("Must provide either 'bond' or 'hamiltonian'")

        # Common initialization
        self.optimizer = optimizer
        self.max_iterations = max_iterations
        self.backend = backend
        self.energy_history = []
```

### Step 2: Add _init_from_components method

```python
def _init_from_components(
    self,
    hamiltonian,
    ansatz,
    mapper,
    **kwargs
):
    """Initialize from individual components (for testing)."""
    self.hamiltonian = hamiltonian
    self.ansatz = ansatz
    self.mapper = mapper if mapper is not None else JordanWignerMapper()

    # Set molecule and bond to None (not available in this mode)
    self.molecule = None
    self.bond = None

    # Disable analysis features (need bond for that)
    self.enable_analysis = False
    self.enable_optimization = False

    logger.info("VQE initialized in component mode (testing)")
```

### Step 3: Verify All Tests Pass

Run tests systematically:

```bash
# Test each VQE test individually
pytest tests/unit/test_vqe.py::TestVQESolver::test_vqe_solver_creation -v
pytest tests/unit/test_vqe.py::TestVQESolver::test_vqe_energy_computation -v
# ... etc
```

### Step 4: Add Strict Scientific Validation

After fixing, add validation tests:

```python
def test_vqe_energy_accuracy_h2():
    """
    STRICT TEST: VQE must match FCI energy within 1 mHa.

    Reference: PySCF FCI
    Molecule: H2 at 0.74 Å
    Basis: STO-3G
    Expected FCI: -1.137 Ha
    Tolerance: 0.001 Ha (1 mHa)
    """
    # Create H2
    h1 = Atom('H', [0.0, 0.0, 0.0])
    h2 = Atom('H', [0.0, 0.0, 0.74])
    molecule = Molecule([h1, h2])

    # Build Hamiltonian
    representation = LCAORepresentation(molecule)
    hamiltonian = CovalentHamiltonian(molecule, representation, 'sto-3g')

    # Get FCI reference
    from pyscf import fci
    fci_energy = fci.FCI(hamiltonian.mf).kernel()[0]

    # VQE with UCC
    ansatz = UCCAnsatz(n_qubits=4, n_electrons=2, include_singles=True, include_doubles=True)
    mapper = JordanWignerMapper()

    solver = VQESolver(
        hamiltonian=hamiltonian,
        ansatz=ansatz,
        mapper=mapper,
        optimizer='SLSQP',
        max_iterations=200
    )

    result = solver.solve()
    vqe_energy = result['energy']

    error = abs(vqe_energy - fci_energy)

    # STRICT REQUIREMENT: <1 mHa error
    assert error < 0.001, f"VQE error {error*1000:.3f} mHa exceeds 1 mHa tolerance"

    # Log for validation report
    logger.info(f"H2 Energy Validation:")
    logger.info(f"  FCI: {fci_energy:.6f} Ha")
    logger.info(f"  VQE: {vqe_energy:.6f} Ha")
    logger.info(f"  Error: {error*1000:.3f} mHa")
```

## Scientific Validation Criteria

### Energy Accuracy (STRICT)

| Molecule | Basis | Method | Max Error | Rationale |
|----------|-------|--------|-----------|-----------|
| H2 | STO-3G | VQE-UCC | 1 mHa | Chemical accuracy |
| H2 | 6-31G | VQE-UCC | 1 mHa | Chemical accuracy |
| LiH | STO-3G | VQE-UCC | 5 mHa | Larger system |
| H2O | STO-3G | VQE-UCC | 10 mHa | Polyatomic |

### Properties (STRICT)

| Property | Molecule | Max Error | Rationale |
|----------|----------|-----------|-----------|
| Dipole | H2O | 5% | Standard for QC |
| Bond length | H2 | 0.01 Å | Spectroscopic |
| Frequency | H2O | 10% | Harmonic approx |

### Numerical Stability (STRICT)

- Energy must be real: `assert np.isreal(energy)`
- Energy must be finite: `assert np.isfinite(energy)`
- Gradient descent: Final ≤ Initial (within tolerance)
- Reproducibility: Same params → same result (±1e-6)

## Acceptance Criteria

### Must Pass (BLOCKING)
- [ ] All 344 unit tests pass
- [ ] Energy accuracy within specified tolerances
- [ ] All properties validated against references
- [ ] Numerical stability verified
- [ ] No regressions in existing tests

### Should Pass (IMPORTANT)
- [ ] Code coverage >95%
- [ ] All edge cases tested
- [ ] Documentation updated
- [ ] Performance benchmarks run

### Nice to Have
- [ ] Additional validation molecules
- [ ] Performance improvements
- [ ] Better error messages

## Next Steps

1. ✅ Implement dual-mode VQESolver
2. ✅ Run all VQE tests
3. ✅ Add strict validation tests
4. ✅ Move to Qiskit integration tests
5. ✅ Fix import tests
6. ✅ Achieve 100% pass rate

**Let's implement the dual-mode VQESolver now!**
