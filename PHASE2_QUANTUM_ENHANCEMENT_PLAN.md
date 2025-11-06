# Phase 2: Quantum Hardware Enhancement Plan

**Date:** November 6, 2025
**Status:** In Progress
**Goal:** Maximize quantum execution across all framework components

---

## 📊 CURRENT STATE

Based on comprehensive codebase analysis:

- **Quantum Infrastructure:** ✅ **Excellent** (VQE, SQD, backends ready)
- **Current Quantum Readiness:** **60%**
- **Target Quantum Readiness:** **95%**
- **Main Work:** Wire existing quantum solvers to applications (not building new algorithms)

---

## 🎯 PHASE 2 PRIORITIES

### Priority 1: Auto-Select SPSA for Cloud Backends ⚡ **1 HOUR**

**Impact:** CRITICAL - 20x reduction in quantum jobs for cloud execution

**Problem:**
- Gradient-based optimizers (SLSQP, L-BFGS-B) use ~40 function evaluations per iteration
- Each function eval = 1 quantum job on IBM/BlueQubit
- Users must manually know to select SPSA for cloud backends

**Solution:**
- Auto-detect when backend is 'ibm' or 'bluequbit'
- Automatically switch to SPSA (2 evaluations/iteration)
- Log warning to inform user

**File:** [kanad/solvers/vqe_solver.py](kanad/solvers/vqe_solver.py) Line 1286-1351

**Implementation:**
```python
def _solve_standard_vqe(self, initial_parameters):
    """Standard VQE with smart optimizer selection."""

    # Auto-select SPSA for cloud backends
    if self.backend in ['ibm', 'bluequbit'] and self.optimizer_method not in ['SPSA', 'COBYLA']:
        logger.warning(f"☁️ Cloud backend '{self.backend}' detected")
        logger.warning(f"📉 Switching from {self.optimizer_method} to SPSA")
        logger.warning(f"   Reason: SPSA uses 2 evals/iter vs ~40 for gradient-based optimizers")
        logger.warning(f"   Expected speedup: 20x fewer quantum jobs")
        self.optimizer_method = 'SPSA'
        self.max_iter = min(self.max_iter, 50)  # SPSA converges faster
```

**Testing:**
```python
# Test with IBM backend
solver = VQESolver(bond=h2_bond, backend='ibm', optimizer='SLSQP')
# Should auto-switch to SPSA and log warning
```

**Estimated Time:** 1 hour
**Expected Result:** Users can use cloud backends without knowing optimizer details

---

### Priority 2: Enable SQD on Quantum Hardware ⚡ **2-3 DAYS**

**Impact:** HIGH - SQD is more noise-resistant than VQE, perfect for NISQ hardware

**Current State:**
- SQD works perfectly on statevector
- Architecture supports quantum backends
- Missing: Quantum Hamiltonian projection implementation

**Problem:**
```python
# kanad/solvers/sqd_solver.py Line 288-295
H_matrix = self.hamiltonian.to_matrix(n_qubits=n_qubits, use_mo_basis=True)
H_sub = np.zeros((n_basis, n_basis), dtype=complex)
for i in range(n_basis):
    for j in range(i, n_basis):
        H_sub[i, j] = np.vdot(basis[i], H_matrix @ basis[j])  # Classical only
```

**Solution:** Implement quantum Hadamard test for `<ψ_i|H|ψ_j>`

**Implementation:**
1. Add quantum Hamiltonian projection method
2. Use EstimatorV2 (IBM) or sampler (BlueQubit)
3. Fallback to statevector if backend not available

**Benefits:**
- No optimization loop (single diagonalization)
- Returns ground + excited states in one shot
- Better for noisy hardware than VQE

**Estimated Time:** 2-3 days
**Expected Result:** `SQDSolver(backend='ibm')` works on real hardware

---

### Priority 3: Add Quantum UV-Vis Spectroscopy ⚡ **2-3 DAYS**

**Impact:** MEDIUM-HIGH - First production quantum UV-Vis calculator

**Current State:**
- UV-Vis uses classical TD-DFT/CIS
- ExcitedStatesSolver with SQD exists but not integrated
- Just needs wiring

**File:** [kanad/analysis/spectroscopy.py](kanad/analysis/spectroscopy.py) Line 59-87

**Implementation:**
```python
def compute_excitations(
    self,
    n_states: int = 5,
    method: str = 'TDA',  # Add 'quantum_sqd' option
    backend: str = 'statevector',  # NEW parameter
    verbose: bool = True
) -> Dict[str, Any]:
    """Compute excited states."""

    if method == 'quantum_sqd':
        from kanad.solvers import ExcitedStatesSolver

        excited_solver = ExcitedStatesSolver(
            bond=self.molecule.bond,
            method='sqd',
            n_states=n_states,
            backend=backend,
            subspace_dim=15
        )

        return excited_solver.solve()
```

**Testing:**
```python
# Test quantum UV-Vis
uv_vis = UVVisCalculator(molecule=benzene)
spectrum = uv_vis.compute_spectrum(method='quantum_sqd', backend='ibm')
```

**Estimated Time:** 2-3 days
**Expected Result:** Quantum UV-Vis spectroscopy available as an option

---

### Priority 4: Integrate Quantum into Drug Discovery ⚡ **1 WEEK**

**Impact:** CRITICAL - Market differentiator, delivers promised accuracy

**Current State:**
```python
# kanad/applications/drug_discovery.py Line 393
E_complex = E_ligand - 0.010  # Placeholder
```

**Problem:** Drug discovery claims quantum advantage but uses placeholders

**Solution:**
```python
def _quantum_binding(self, ligand, target, pH, temperature, solvent):
    """Real quantum binding affinity calculation."""
    from kanad.solvers import SQDSolver
    from kanad.bonds import BondFactory

    # 1. Create ligand bond
    ligand_bond = BondFactory.create_from_molecule(ligand)

    # 2. Run quantum SQD
    solver = SQDSolver(
        bond=ligand_bond,
        subspace_dim=12,
        backend=self.backend,  # Use IBM/BlueQubit
        **self.backend_kwargs
    )
    result = solver.solve(n_states=1)

    # 3. Compute binding energy
    E_ligand = result['ground_state_energy']
    Delta_G = (E_ligand + protein_correction) * 627.509  # Ha → kcal/mol

    return BindingResult(affinity=Delta_G, method='quantum_sqd', ...)
```

**Benefits:**
- Delivers promised <1 kcal/mol accuracy
- Differentiates from SwissADME (force fields)
- Validates "10x better accuracy" claims

**Estimated Time:** 1 week
**Expected Result:** Real quantum binding affinity calculations

---

### Priority 5: Test on Real Hardware ⚡ **ONGOING**

**Testing Plan:**

1. **VQE with SPSA on IBM Quantum:**
   - H2 molecule with IBM backend
   - Verify auto-switch to SPSA
   - Compare cost: SLSQP vs SPSA

2. **SQD on IBM Quantum:**
   - H2 excited states
   - Verify accuracy vs statevector
   - Test error mitigation

3. **Quantum UV-Vis:**
   - Benzene molecule
   - Compare quantum vs classical excitation energies
   - Validate against experiment

4. **Drug Discovery Binding:**
   - Small ligand (caffeine, aspirin)
   - Quantum binding affinity
   - Compare to experimental data

---

## 📈 SUCCESS METRICS

| Metric | Current | Target | Status |
|--------|---------|--------|--------|
| **Quantum Readiness** | 60% | 95% | In Progress |
| **Cloud Backend Efficiency** | 40 jobs/iter | 2 jobs/iter | Pending |
| **SQD Hardware Support** | Statevector only | IBM + BlueQubit | Pending |
| **Quantum UV-Vis** | Not available | Production-ready | Pending |
| **Drug Discovery Quantum** | Placeholder | Real calculations | Pending |

---

## 🗓️ TIMELINE

### Week 1 (Nov 6-12)
- ✅ Phase 1 Cleanup (DONE)
- ⏳ Priority 1: Auto-select SPSA (1 hour)
- ⏳ Priority 2: SQD quantum support (2-3 days)
- ⏳ Priority 3: Quantum UV-Vis (2-3 days)

### Week 2 (Nov 13-19)
- Priority 4: Drug Discovery integration (1 week)
- Testing on IBM Quantum hardware

### Week 3-4 (Nov 20-Dec 3)
- Advanced features (pH-dependent, metabolism)
- Performance optimization
- Documentation

---

## 🎯 QUICK WINS FIRST

**Today's Goals:**
1. Auto-select SPSA for cloud backends (1 hour) ← **START HERE**
2. Begin SQD quantum implementation (start design)
3. Plan quantum UV-Vis integration

**Why this order?**
- SPSA auto-selection has IMMEDIATE impact on all cloud users
- SQD needs careful design but high value
- UV-Vis is mostly wiring, good for momentum

---

## 📝 IMPLEMENTATION NOTES

### Auto-Select SPSA
- **File:** kanad/solvers/vqe_solver.py
- **Lines to modify:** 1286-1293
- **Testing:** H2 with IBM backend, verify SPSA auto-selected

### SQD Quantum Support
- **File:** kanad/solvers/sqd_solver.py
- **New methods:** `_project_hamiltonian_quantum()`, `_hadamard_test()`
- **Backend integration:** Use EstimatorV2 for IBM, sampler for BlueQubit

### Quantum UV-Vis
- **File:** kanad/analysis/spectroscopy.py
- **Change:** Add `method='quantum_sqd'` option
- **Integration:** Call ExcitedStatesSolver with SQD backend

### Drug Discovery
- **File:** kanad/applications/drug_discovery.py
- **Function:** `_quantum_binding()` Line 360-420
- **Change:** Replace placeholder with SQDSolver

---

## 🚀 NEXT STEPS

**Immediate Actions:**
1. Implement SPSA auto-selection
2. Test with IBM backend
3. Commit and move to SQD enhancement

**Let's start with Priority 1 - it's quick and has huge impact!**

---

*Generated: November 6, 2025*
*Phase 2: Quantum Enhancement - IN PROGRESS*
