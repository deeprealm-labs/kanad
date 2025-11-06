# Remaining Work Plan - Hi-VQE Full Integration

**Date:** November 4, 2025
**Goal:** Complete Hi-VQE integration across all components and prepare for cloud deployment

---

## Phase 1: Core Integration (CURRENT)

### Task 1.1: Active Space Integration with Hamiltonian Construction ⚡ HIGH PRIORITY

**Problem:** Currently, active space selection works independently, but Hamiltonians are built with full orbital space. This causes qubit count mismatches.

**Solution:** Integrate active space directly into Hamiltonian construction.

**Files to Modify:**
1. `kanad/core/hamiltonians/covalent_hamiltonian.py`
2. `kanad/core/hamiltonians/ionic_hamiltonian.py`
3. `kanad/core/hamiltonians/metallic_hamiltonian.py`
4. `kanad/core/hamiltonians/openfermion_jw.py`

**Changes Needed:**
- Add `frozen_orbitals` and `active_orbitals` parameters to Hamiltonian constructors
- Build fermionic operators only for active space
- Add frozen core energy contribution
- Update `to_sparse_hamiltonian()` to use active space

**Expected Impact:**
- LiH: 12 → 10 qubits (works properly)
- BeH: Works with qubit reduction
- All molecules: Consistent qubit count between active space and Hamiltonian

### Task 1.2: Update All Hamiltonian Types

**Hamiltonian Types to Update:**
1. ✅ CovalentHamiltonian (governance-aware)
2. ✅ IonicHamiltonian (governance-aware)
3. ✅ MetallicHamiltonian (governance-aware)
4. ⚠️ Standard molecular Hamiltonian (MolecularHamiltonian in molecule.py)

**Ensure all support:**
- Active space construction
- Frozen core energy
- Correct qubit count
- Governance protocol integration

### Task 1.3: VQE Solver Hi-VQE Mode Integration

**File:** `kanad/utils/vqe_solver.py`

**Add Hi-VQE mode:**
```python
def solve(self, mode='standard'):
    """
    Solve VQE problem.

    Args:
        mode: 'standard' (Pauli measurements) or 'hivqe' (configuration sampling)
    """
    if mode == 'hivqe':
        return self._solve_hivqe()
    else:
        return self._solve_standard()
```

**Implementation:**
1. Configuration sampling from quantum state
2. Subspace expansion with governance guidance
3. Classical diagonalization
4. Parameter optimization to sample important configs better
5. Iteration until convergence

### Task 1.4: Circuit Preparation from Configurations

**Purpose:** Prepare quantum circuits from sampled configurations for next iteration

**Implementation:**
- Convert configuration amplitudes to circuit parameters
- Use sampled configs to initialize next VQE iteration
- Integrate with existing ansatz construction

---

## Phase 2: Governance-Guided Excitations

### Task 2.1: Extend Governance Protocols with Excitation Methods

**Files to Modify:**
1. `kanad/governance/protocols/covalent_protocol.py`
2. `kanad/governance/protocols/ionic_protocol.py`
3. `kanad/governance/protocols/metallic_protocol.py`

**Methods to Add:**
```python
class CovalentGovernanceProtocol:
    def generate_single_excitations(self, config: str) -> List[str]:
        """Generate only bonding-relevant single excitations."""
        # Only excitations within/between bonding pairs
        # Preserve hybridization character
        # Maintain spin pairing

    def generate_double_excitations(self, config: str) -> List[str]:
        """Generate only MO-preserving double excitations."""
        # Only excitations that preserve bonding/antibonding character

    def is_physical_excitation(self, config_i: str, config_j: str) -> bool:
        """Validate excitation preserves physical constraints."""
        # Check spin symmetry
        # Check hybridization preservation
        # Check no long-range entanglement (covalent)
```

**Expected Impact:**
- Reduce excitation count by 5-10x
- Faster subspace growth
- More efficient convergence

### Task 2.2: Integrate with ConfigurationSubspace

**File:** `kanad/core/configuration.py`

**Update methods to use governance:**
- `generate_single_excitations()` - call protocol if available
- `generate_double_excitations()` - call protocol if available
- `filter_configs()` - use protocol validation

---

## Phase 3: Cloud Backend Integration

### Task 3.1: IBM Quantum Cloud Support

**File:** `kanad/utils/vqe_solver.py`

**Add cloud backend support:**
```python
def _solve_hivqe_cloud(self, backend='ibm_brisbane'):
    """
    Run Hi-VQE on IBM Quantum cloud backend.

    Key optimizations:
    - Shot budget allocation
    - Error mitigation (readout, gate errors)
    - Job batching
    - Result caching
    """
```

**Optimizations Needed:**
1. **Shot Budget Allocation:**
   - Allocate more shots to important configurations
   - Adaptive sampling based on amplitude estimates

2. **Error Mitigation:**
   - Readout error mitigation
   - Gate error calibration
   - Noise-aware circuit compilation

3. **Job Management:**
   - Batch multiple measurements
   - Queue management
   - Result retrieval and caching

### Task 3.2: Configuration Sampling from Real Backend

**File:** `kanad/core/configuration.py`

**Update for real quantum hardware:**
```python
def sample_configurations_from_backend(
    circuit: QuantumCircuit,
    backend: IBMBackend,
    shots: int = 1000,
    n_electrons: int = None,
    error_mitigation: bool = True
) -> List[Tuple[str, int]]:
    """
    Sample configurations from real quantum backend.

    Handles:
    - Measurement error mitigation
    - Shot noise
    - Invalid configuration filtering
    """
```

### Task 3.3: Noise-Aware Circuit Preparation

**Considerations:**
- Circuit depth minimization
- Native gate sets for IBM hardware
- Qubit connectivity constraints
- Transpilation optimization

---

## Phase 4: Testing & Validation

### Task 4.1: Comprehensive Molecule Test Suite

**Test molecules:**
1. ✅ H2 (2 electrons, 4 qubits)
2. ✅ H2O (10 electrons, 14 qubits → 12 with active space)
3. ⚠️ LiH (4 electrons, 12 qubits → 10 with active space)
4. ⚠️ BeH (5 electrons, needs testing)
5. ⚠️ NH3 (10 electrons, 16 qubits → 14 with active space)
6. ⚠️ CH4 (larger test case)

**Metrics to Track:**
- Energy accuracy vs FCI
- Convergence iterations
- Measurement reduction
- Subspace size
- Wall-clock time

### Task 4.2: Benchmarking vs Literature

**Compare against:**
- Qunova's Hi-VQE results
- Standard VQE results
- Classical methods (CCSD, FCI)

**Generate publishable stats:**
- Error bars
- Scaling analysis
- Noise resilience

---

## Implementation Priority

### Week 1: Core Integration (Current Week)
1. ✅ Active space implementation
2. ✅ Configuration sampling
3. ✅ Classical diagonalization
4. ⚠️ **Active space + Hamiltonian integration** (DO NOW)
5. ⚠️ **VQE solver Hi-VQE mode** (DO NOW)

### Week 2: Governance & Optimization
1. Governance-guided excitations
2. Smart configuration filtering
3. Subspace pruning optimization

### Week 3: Cloud Integration
1. IBM Quantum backend support
2. Error mitigation
3. Shot budget optimization
4. Job batching

### Week 4: Testing & Publishing
1. Comprehensive benchmarks
2. Comparison with literature
3. Documentation
4. Publication preparation

---

## Immediate Next Steps (This Session)

### 1. Fix Active Space + Hamiltonian Integration ⚡
- Modify all Hamiltonian types to accept frozen/active orbitals
- Test with LiH, BeH, H2O
- Verify qubit count matches

### 2. Add Hi-VQE Mode to VQE Solver
- Implement `_solve_hivqe()` method
- Integrate configuration sampling
- Add classical diagonalization loop
- Test with H2

### 3. Test Full Pipeline
- H2 with Hi-VQE mode
- H2O with active space + Hi-VQE
- LiH with active space + Hi-VQE

### 4. Document & Prepare for Governance Excitations
- Clean up code
- Add docstrings
- Prepare governance protocol extensions

---

## Success Criteria

### Must Have (Production Ready):
- ✅ 100-1000x measurement reduction
- ✅ 2-10 iteration convergence
- ⚠️ Active space working for all molecules
- ⚠️ Hi-VQE mode in VQE solver
- ⚠️ Cloud backend support

### Should Have (Enhanced):
- Governance-guided excitations (5-10x subspace reduction)
- Error mitigation for noisy backends
- Adaptive shot allocation

### Nice to Have (Future):
- Multi-reference configurations
- Excited state calculations
- Dynamic correlation recovery

---

**Let's start with Task 1.1: Active Space + Hamiltonian Integration!**
