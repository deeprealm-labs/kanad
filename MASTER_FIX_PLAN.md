# 🔧 MASTER FIX PLAN - Systematic Framework Completion

**Date:** November 6, 2025
**Status:** 🚧 **IN PROGRESS - SYSTEMATIC FIXES**
**Goal:** Fix all critical issues, validate comprehensively, ensure production readiness

---

## Executive Summary

Based on comprehensive investigation, we have **5 critical issues** blocking proper quantum functionality. This plan addresses them systematically with clear priorities, testing, and validation.

**Approach:** Bottom-up fixes - start with foundational issues, test incrementally, validate thoroughly.

---

## Issue Hierarchy & Dependencies

```
Critical Issues (Must Fix):
1. ✅ Density Matrix Extraction (DONE - Phase 1)
   ↓ (blocks)
2. ⏳ Quantum Property Calculations (NMR, Raman, Dipole)
   ↓ (needs)
3. ⏳ Governance Integration in Solvers
   ↓ (enables)
4. ⏳ Error Mitigation Automation
   ↓ (improves)
5. ⏳ Environment Effects Completion

High Priority (Should Fix):
6. ⏳ Gradient-based Optimization
7. ⏳ Active Space Validation
8. ⏳ Open-Shell Support
... (8 total)

Medium Priority (Nice to Have):
9. ⏳ Excited States
10. ⏳ PDOS Improvements
... (6 total)
```

---

## PHASE 1: Foundation (COMPLETED ✅)

### What We Fixed
- ✅ Added `get_density_matrix()` to Hamiltonians
- ✅ Property calculator extracts from Hamiltonians
- ✅ Removed hardcoded uniform density fallbacks
- ✅ Validated with tests (H2, LiH)

**Time:** 1 hour
**Files Modified:** 4 files, +80 lines

---

## PHASE 2: Quantum Properties (CURRENT - 4-6 hours)

### Critical Issue #2: NMR/Raman Using HF Fallback

**Goal:** Make quantum properties actually use quantum corrections

### 2.1 Improve NMR Quantum Corrections (2 hours)

**Current Problem:**
```python
# nmr_calculator.py:429
correlation_factor = correlation_energy * 10.0  # ❌ Linear scaling, arbitrary
```

**Fix:**
```python
# Better correlation correction for NMR shielding
def _compute_quantum_nmr_correction(self, correlation_energy, hf_energy, atom_type):
    """
    Compute bonding-aware NMR shielding correction.

    Correlation effects on NMR:
    - Typical: 5-20 ppm for light atoms
    - Varies by atom type and bonding
    """
    # Percentage correlation
    corr_fraction = abs(correlation_energy / hf_energy)

    # Atom-specific scaling (from literature)
    scaling = {
        'H': 15.0,  # ppm per 1% correlation
        'C': 25.0,
        'N': 30.0,
        'O': 35.0,
        'F': 40.0
    }

    atom_scale = scaling.get(atom_type, 20.0)

    # Bonding correction (from bond type)
    bond_factor = 1.0
    if hasattr(self, 'bond_type'):
        if self.bond_type == 'covalent':
            bond_factor = 1.2  # More delocalized
        elif self.bond_type == 'ionic':
            bond_factor = 0.8  # More localized

    return corr_fraction * atom_scale * bond_factor * 100  # Convert to ppm
```

**Test:**
```python
# test_quantum_nmr_corrections.py
def test_h2o_nmr_quantum_vs_classical():
    """Validate quantum NMR shows correlation effects."""
    molecule = create_h2o()
    nmr = NMRCalculator(molecule)

    classical = nmr.compute_chemical_shifts(method='HF')
    quantum = nmr.compute_chemical_shifts(method='sqd')

    # Should be different but not wildly (within 20 ppm)
    assert not np.allclose(classical['shifts'], quantum['shifts'], atol=1.0)
    assert all(abs(q - c) < 20.0 for q, c in zip(quantum['shifts'], classical['shifts']))

    # H atoms should have different shifts
    assert len(set(quantum['shifts'])) > 1
```

### 2.2 Fix Raman Polarizability (2 hours)

**Current Problem:**
```python
# raman_calculator.py:166
alpha_iso = n_electrons * 0.8  # ❌ Hardcoded formula
```

**Fix:**
```python
def _compute_polarizability_from_density(self, density_matrix, mo_energies):
    """
    Compute polarizability from quantum density matrix.

    α = 2 Σ_{occ,virt} |⟨occ|μ|virt⟩|² / (E_virt - E_occ)

    Where μ is dipole operator
    """
    n_occ = self.molecule.n_electrons // 2

    # Get dipole integrals in MO basis
    dipole_x = self._compute_dipole_integrals('x')
    dipole_y = self._compute_dipole_integrals('y')
    dipole_z = self._compute_dipole_integrals('z')

    alpha = 0.0
    for i in range(n_occ):  # Occupied
        for a in range(n_occ, len(mo_energies)):  # Virtual
            delta_E = mo_energies[a] - mo_energies[i]
            if delta_E < 0.1:  # Skip near-degenerate
                continue

            # Transition dipole moment
            mu_x = dipole_x[i, a]
            mu_y = dipole_y[i, a]
            mu_z = dipole_z[i, a]
            mu_sq = mu_x**2 + mu_y**2 + mu_z**2

            alpha += 2.0 * mu_sq / delta_E

    return alpha / 3.0  # Isotropic average
```

**Test:**
```python
def test_h2_raman_quantum_vs_classical():
    """Validate quantum Raman polarizability."""
    bond = BondFactory.create_bond('H', 'H', distance=0.74)
    raman = RamanCalculator(bond)

    classical = raman.compute_raman_spectrum(method='HF')
    quantum = raman.compute_raman_spectrum(method='sqd')

    # Quantum should be within 2x of classical (not 1500x!)
    ratio = quantum['intensities'][0] / classical['intensities'][0]
    assert 0.5 < ratio < 2.0, f"Ratio {ratio} should be near 1.0"
```

---

## PHASE 3: Governance Integration (4-6 hours)

### Critical Issue #3: Governance Not Used in Subspace Construction

**Current Problem:**
```python
# sqd_solver.py:293-312
def _get_governance_protocol(self):
    # ✅ Gets protocol
    # ❌ But never uses it for subspace!
```

**Fix:**
```python
def _generate_governance_guided_subspace(self, n_basis_states):
    """
    Generate subspace using governance protocol.

    Governance selects important configurations based on bonding type:
    - Covalent: Shared electron pairs, low-spin states
    - Ionic: Charge-separated, local excitations
    - Metallic: Delocalized, band-like states
    """
    protocol = self._get_governance_protocol()

    if protocol is None:
        # Fallback to random
        return self._generate_random_subspace(n_basis_states)

    # Get governance-filtered configurations
    important_configs = protocol.select_important_configurations(
        n_qubits=self.n_qubits,
        n_electrons=self.hamiltonian.n_electrons,
        max_configs=n_basis_states
    )

    # Build circuits for these configurations
    basis_circuits = []
    for config in important_configs:
        circuit = self._build_circuit_for_configuration(config)
        basis_circuits.append(circuit)

    return basis_circuits
```

**Integration Points:**
1. SQDSolver: Use in `_generate_basis_states()` (line 334)
2. VQESolver: Use in ansatz construction (line 245)
3. Validate speedup: Should see 5-10x reduction in subspace size

**Test:**
```python
def test_governance_subspace_reduction():
    """Validate governance reduces subspace size."""
    bond = BondFactory.create_bond('H', 'H', distance=0.74)

    # Without governance
    solver_no_gov = SQDSolver(bond, use_governance=False, subspace_dim=100)
    result_no_gov = solver_no_gov.solve()

    # With governance
    solver_gov = SQDSolver(bond, use_governance=True, subspace_dim=20)  # 5x smaller
    result_gov = solver_gov.solve()

    # Should get similar energy with much smaller subspace
    assert abs(result_gov['energies'][0] - result_no_gov['energies'][0]) < 1e-3
    assert solver_gov.subspace_dim < solver_no_gov.subspace_dim / 3
```

---

## PHASE 4: Error Mitigation Automation (2-3 hours)

### Critical Issue #4: Manual Error Mitigation Configuration

**Current Problem:**
```python
# Error mitigation exists but not automatic
# Must manually configure ZNE, readout error, etc.
```

**Fix:**
```python
# backends/ibm/backend.py
def _auto_configure_error_mitigation(self, backend_name):
    """
    Automatically configure error mitigation based on backend.

    - Simulators: No mitigation needed
    - Real hardware: Full mitigation stack
    """
    if 'simulator' in backend_name.lower() or 'statevector' in backend_name.lower():
        return {
            'zne_enabled': False,
            'readout_mitigation': False,
            'dynamical_decoupling': False
        }

    # Real hardware - enable all mitigations
    return {
        'zne_enabled': True,
        'zne_extrapolation': 'exponential',
        'zne_scale_factors': [1, 3, 5],
        'readout_mitigation': True,
        'readout_method': 'tensored',
        'dynamical_decoupling': True,
        'dd_sequence': 'XY4'
    }
```

**Test:**
```python
def test_auto_error_mitigation():
    """Validate automatic error mitigation."""
    from kanad.backends.ibm import IBMQuantumBackend

    # Simulator - no mitigation
    backend_sim = IBMQuantumBackend('aer_simulator')
    assert not backend_sim.mitigation_config['zne_enabled']

    # Real hardware - full mitigation
    backend_hw = IBMQuantumBackend('ibm_kyoto')
    assert backend_hw.mitigation_config['zne_enabled']
    assert backend_hw.mitigation_config['readout_mitigation']
```

---

## PHASE 5: Environment Effects (2-3 hours)

### Critical Issue #5: Placeholder Environment Implementations

**Current Problems:**
1. Temperature effects: Boltzmann factors not computed
2. Pressure effects: Volume changes not propagated
3. pH effects: Protonation states hardcoded

**Fixes:**

### 5.1 Temperature Modulator
```python
# environment/temperature_modulator.py
def compute_thermal_population(self, energies, temperature=298.15):
    """
    Compute Boltzmann populations at temperature T.

    P_i = exp(-E_i/kT) / Σ_j exp(-E_j/kT)
    """
    k_B = 8.617333e-5  # eV/K
    beta = 1.0 / (k_B * temperature)

    # Boltzmann factors
    exp_factors = np.exp(-beta * np.array(energies))
    Z = np.sum(exp_factors)  # Partition function

    populations = exp_factors / Z
    return populations
```

### 5.2 Pressure Modulator
```python
def compute_volume_change(self, pressure, bulk_modulus):
    """
    Compute volume change under pressure.

    ΔV/V₀ = -P/B
    where B is bulk modulus
    """
    volume_fraction = 1.0 - pressure / bulk_modulus
    return volume_fraction
```

### 5.3 pH Modulator
```python
def determine_protonation_state(self, molecule, pH):
    """
    Determine protonation state based on pKa values.

    Uses Henderson-Hasselbalch equation:
    pH = pKa + log([A-]/[HA])
    """
    protonation_sites = self._identify_ionizable_groups(molecule)

    protonated = {}
    for site, pKa in protonation_sites.items():
        # Fraction protonated
        delta_pH = pH - pKa
        fraction = 1.0 / (1.0 + 10**delta_pH)
        protonated[site] = fraction > 0.5

    return protonated
```

---

## PHASE 6: High Priority Fixes (6-8 hours)

### 6.1 Gradient-Based Optimization
### 6.2 Active Space Validation
### 6.3 Open-Shell Support
### 6.4 Metallic Hamiltonian Improvements
### 6.5 Correlation Methods
### 6.6 Basis Set Expansion
### 6.7 Configuration Explorer Completion
### 6.8 Spectroscopy Excited States

*(Details in separate sections below)*

---

## PHASE 7: Comprehensive Testing (4-6 hours)

### 7.1 Unit Tests for Each Fix
- Test each fixed function individually
- Validate inputs/outputs
- Check edge cases

### 7.2 Integration Tests
- Test solver + property calculator
- Test governance + error mitigation
- Test environment + solvers

### 7.3 Validation Tests
- Compare to known exact results
- Validate against literature values
- Check physical constraints

### 7.4 Performance Tests
- Measure governance speedup
- Validate error mitigation improvement
- Check memory usage

---

## PHASE 8: Production Validation (2-3 hours)

### 8.1 Benchmark Suite
```python
# tests/benchmarks/benchmark_suite.py
class BenchmarkSuite:
    """Comprehensive benchmark suite for production validation."""

    def benchmark_h2_dissociation(self):
        """H2 dissociation curve - exact results known."""
        distances = np.linspace(0.5, 3.0, 10)
        energies = []

        for d in distances:
            bond = BondFactory.create_bond('H', 'H', distance=d)
            solver = VQESolver(bond, backend='statevector', use_governance=True)
            result = solver.solve()
            energies.append(result['energy'])

        # Compare to exact FCI
        exact_energies = self._get_exact_h2_energies(distances)
        errors = np.abs(np.array(energies) - exact_energies)

        assert np.max(errors) < 1e-3, "H2 energies should match FCI"

    def benchmark_h2o_properties(self):
        """H2O properties - compare to experiment."""
        molecule = create_h2o()

        # Dipole moment
        prop_calc = PropertyCalculator(molecule)
        dipole = prop_calc.compute_quantum_molecular_properties(method='sqd')
        assert 1.8 < dipole['dipole_magnitude'] < 1.9, "H2O dipole ~1.85 D"

        # NMR shifts
        nmr = NMRCalculator(molecule)
        shifts = nmr.compute_chemical_shifts(method='sqd')
        # O should be downfield from H
        assert shifts['O'] > shifts['H']
```

### 8.2 Accuracy Validation
- NMR shifts within 10 ppm of experiment
- Raman intensities within 2x of experiment
- Dipole moments within 0.1 D
- Bond energies within 1 kcal/mol

### 8.3 Performance Validation
- Governance speedup 5-10x (measured)
- Error mitigation improves fidelity 2-3x
- Memory usage < 8 GB for 20-qubit systems

---

## Timeline & Effort Estimate

| Phase | Tasks | Estimated Time | Priority |
|-------|-------|---------------|----------|
| ✅ Phase 1 | Foundation fixes | 1 hour | CRITICAL |
| 🔄 Phase 2 | Quantum properties | 4-6 hours | CRITICAL |
| ⏳ Phase 3 | Governance integration | 4-6 hours | CRITICAL |
| ⏳ Phase 4 | Error mitigation | 2-3 hours | CRITICAL |
| ⏳ Phase 5 | Environment effects | 2-3 hours | CRITICAL |
| ⏳ Phase 6 | High priority fixes | 6-8 hours | HIGH |
| ⏳ Phase 7 | Comprehensive testing | 4-6 hours | CRITICAL |
| ⏳ Phase 8 | Production validation | 2-3 hours | CRITICAL |
| **TOTAL** | | **25-36 hours** | **3-5 days** |

---

## Success Criteria

### Must Have (Blocking Release)
- [ ] All 5 critical issues fixed
- [ ] Comprehensive test suite passing (>95%)
- [ ] Validation benchmarks meet accuracy targets
- [ ] No hardcoded fallbacks or placeholders
- [ ] Performance metrics validated

### Should Have (High Priority)
- [ ] 8 high-priority issues addressed
- [ ] Gradient-based optimization working
- [ ] Open-shell systems supported
- [ ] Active space validated

### Nice to Have (Future)
- [ ] 6 medium-priority issues improved
- [ ] Excited states working
- [ ] Full correlation methods
- [ ] Extended basis sets

---

## Risk Mitigation

### Risk 1: Fixes Break Existing Features
**Mitigation:** Test incrementally, keep backups, use git branches

### Risk 2: Performance Degrades
**Mitigation:** Benchmark before/after, profile bottlenecks

### Risk 3: Accuracy Doesn't Improve
**Mitigation:** Validate against exact results, compare to literature

### Risk 4: Timeline Overruns
**Mitigation:** Prioritize ruthlessly, defer nice-to-haves

---

## Execution Strategy

### Day 1 (8 hours)
- Morning: Phase 2 (NMR/Raman fixes)
- Afternoon: Phase 3 (Governance integration)
- Evening: Test Phase 2 & 3

### Day 2 (8 hours)
- Morning: Phase 4 (Error mitigation)
- Afternoon: Phase 5 (Environment effects)
- Evening: Phase 6 (Start high priority)

### Day 3 (8 hours)
- Morning: Phase 6 (Finish high priority)
- Afternoon: Phase 7 (Comprehensive testing)
- Evening: Phase 8 (Production validation)

### Day 4-5 (Flex)
- Buffer for issues
- Additional validation
- Documentation updates

---

## Status Tracking

| Issue | Status | Owner | Completion |
|-------|--------|-------|------------|
| Density Matrix | ✅ DONE | Claude | 100% |
| NMR Corrections | 🔄 IN PROGRESS | | 0% |
| Raman Polarizability | ⏳ PENDING | | 0% |
| Governance Integration | ⏳ PENDING | | 0% |
| Error Mitigation | ⏳ PENDING | | 0% |
| Environment Effects | ⏳ PENDING | | 0% |
| High Priority Fixes | ⏳ PENDING | | 0% |
| Testing Campaign | ⏳ PENDING | | 0% |
| Production Validation | ⏳ PENDING | | 0% |

---

**Date:** November 6, 2025
**Status:** 🚧 **PHASE 2 STARTING**
**Next:** Fix NMR quantum corrections
