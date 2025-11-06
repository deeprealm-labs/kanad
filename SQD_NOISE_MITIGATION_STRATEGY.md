# Noise Mitigation Strategy for SQD and Krylov-SQD

## Executive Summary

Our IBM Torino test showed **124 mHa error** for H2 molecule. Research shows this can scale poorly to larger molecules without mitigation. However, multiple proven strategies can reduce this to **<1 mHa (chemical accuracy)** even for larger systems.

## Problem: Noise Scaling with Molecule Size

| System | Qubits | Without Mitigation | Target |
|--------|--------|-------------------|---------|
| H2 | 4 | 124 mHa | <1 mHa |
| LiH | 8 | ~500 mHa (est.) | <5 mHa |
| H2O | 14 | ~2000 mHa (est.) | <10 mHa |
| N2 | 20 | ~10000 mHa (est.) | <50 mHa |

**Challenge**: Noise grows exponentially with system size without mitigation.

---

## Mitigation Strategies (Ranked by Impact)

### 1. Symmetry-Based Error Correction ⭐⭐⭐⭐⭐
**Impact**: 10-20x error reduction
**Cost**: ~0% overhead (classical post-processing)
**Status**: ✅ Already partially implemented in Kanad

**How it works**:
- SQD samples can violate physical symmetries due to noise
- Particle number conservation: reject/correct samples with wrong electron count
- Spin symmetry: reject/correct samples with wrong total spin
- Spatial symmetry: use molecular point group symmetries

**Implementation**:
```python
# Already exists in kanad.governance.protocols.covalent_protocol
from kanad.governance.protocols import CovalentGovernanceProtocol

protocol = CovalentGovernanceProtocol()
# Automatically enforces:
# - Particle number conservation
# - Spin symmetry
# - Spatial symmetries
```

**Next Steps**:
- Extend to molecular point group symmetries (C2v, D2h, etc.)
- Add iterative refinement loop for noisy samples
- Implement configuration recovery from invalid samples

---

### 2. Reference-State Error Mitigation (REM) ⭐⭐⭐⭐⭐
**Impact**: 10-100x error reduction
**Cost**: 2x shots (calibration + measurement)
**Status**: 🔨 Needs implementation

**How it works**:
- Measure known reference states (HF, simple excitations) on quantum hardware
- Compare to exact theoretical values
- Use discrepancy to build error mitigation map
- Apply correction to all measurements

**References**:
- Paper: "Reference-State Error Mitigation" (J. Chem. Theory Comput. 2023)
- Demonstrated on IBM hardware: up to 100x improvement
- Works with circuit depth >1000 gates

**Implementation Plan**:
```python
class ReferenceStateErrorMitigator:
    def __init__(self, hamiltonian, reference_states):
        """
        Args:
            reference_states: Dict of known states and their exact energies
                e.g., {'HF': -1.1373, 'LUMO': -0.4234, ...}
        """
        self.reference_map = {}

    def calibrate(self, backend):
        """Measure reference states on hardware, build error map"""
        for state_name, exact_energy in self.reference_states.items():
            measured = self._measure_on_hardware(state_name, backend)
            self.reference_map[state_name] = exact_energy - measured

    def mitigate(self, measured_energy, state_type='general'):
        """Apply error correction based on reference map"""
        correction = self._interpolate_correction(state_type)
        return measured_energy + correction
```

**Next Steps**:
1. Implement `ReferenceStateErrorMitigator` class
2. Integrate with SQDSolver and KrylovSQDSolver
3. Test on H2, LiH on IBM hardware

---

### 3. Multi-basis Measurements ⭐⭐⭐⭐
**Impact**: 5-10x error reduction, 3x fewer shots
**Cost**: More complex circuit compilation
**Status**: 🔨 Needs implementation

**How it works**:
- Instead of measuring in computational basis only (Z-basis)
- Rotate to X, Y bases for some measurements
- Pauli grouping: measure commuting operators together
- Reduces total measurement shots by ~3x

**References**:
- Paper: "Efficient and noise resilient measurements" (npj Quantum 2020)
- Low-rank factorization: 1000x measurement reduction for large molecules

**Implementation Plan**:
```python
class MultiBasesMeasurementOptimizer:
    def optimize_measurement_bases(self, hamiltonian):
        """
        Group Hamiltonian terms by commuting Pauli operators
        Returns: List of measurement bases and associated terms
        """
        # 1. Parse Hamiltonian into Pauli terms
        # 2. Build commutation graph
        # 3. Color graph (group commuting terms)
        # 4. For each group, find optimal measurement basis
        pass

    def measure_multi_basis(self, circuit, bases, backend):
        """Execute measurements in multiple bases"""
        pass
```

---

### 4. Sample-based Krylov (SKQD) ⭐⭐⭐⭐
**Impact**: Better convergence + noise resilience
**Cost**: Moderate (hybrid quantum-classical)
**Status**: 🔨 Major implementation needed

**How it works**:
- Combines Krylov subspace methods with sample-based approach
- Builds Krylov space K_m(H, |ψ₀⟩) using noisy samples
- Self-consistent iteration: sample → project → refine → repeat
- Guaranteed convergence even with significant noise

**References**:
- Paper: "Sample-based Krylov Quantum Diagonalization" (arXiv:2501.09702, 2025)
- IBM/U.Tokyo: Demonstrated on 56 qubits

**Implementation Plan**:
This is a major upgrade to our current Krylov-SQD:
1. Replace state preparation with sampling from ansatz circuits
2. Add iterative refinement loop
3. Implement Hadamard test for H|ψ⟩ on real hardware
4. Add convergence checks

---

### 5. M3 Readout Error Mitigation ⭐⭐⭐
**Impact**: 2-5x error reduction
**Cost**: 15% overhead (calibration)
**Status**: ✅ Available in `kanad.error_mitigation.lite_mitigation`

**How it works**:
- Matrix-free approach to readout error correction
- Calibrate measurement errors on 4-qubit neighborhoods
- Apply Bayesian correction to measurement outcomes
- Amortized cost ~0% with caching

**Already implemented!**
```python
from kanad.error_mitigation import create_lite_mitigator

mitigator = create_lite_mitigator(n_qubits=4)
mitigator.calibrate(backend)  # One-time cost
corrected_counts = mitigator.apply(noisy_counts)
```

---

### 6. Zero-Noise Extrapolation (ZNE) ⭐⭐⭐
**Impact**: 3-5x error reduction
**Cost**: 3x shots (multiple noise levels)
**Status**: 🔨 Needs implementation

**How it works**:
- Run circuit at multiple noise levels (stretch gates)
- Extrapolate to zero-noise limit
- Works well for coherent errors

**References**:
- Standard technique in Qiskit Runtime
- "Error mitigation extends computational reach" (Nature 2019)

---

### 7. Dynamical Decoupling (DD) ⭐⭐
**Impact**: 1.5-2x error reduction
**Cost**: Longer circuit execution time
**Status**: ⚠️ Attempted but may not be enabled on IBM

**How it works**:
- Insert Pauli pulse sequences during idle times
- Decouples qubits from environmental noise
- Reduces decoherence errors

**Current Status**:
```python
# Already tried in test_sqd_ibm_torino.py
try:
    sampler.options.dynamical_decoupling.enable = True
except:
    pass  # May not be available
```

---

## Recommended Implementation Priority

### Phase 1: Quick Wins (1-2 weeks) 🎯
1. **Enhance symmetry filtering** in existing SQD
   - Add point group symmetries
   - Implement configuration recovery
   - Target: 10-20x error reduction

2. **Integrate M3 mitigation** with SQD
   - Already have LiteMitigator
   - Add to diagonal measurement step
   - Target: additional 2-3x improvement

**Expected Result**: 124 mHa → ~5-10 mHa for H2

---

### Phase 2: Advanced Mitigation (2-4 weeks) 🚀
3. **Reference-State Error Mitigation**
   - Implement REM for SQD diagonal measurements
   - Calibrate on HF + low-lying excitations
   - Target: 5-10x additional improvement

4. **Multi-basis Measurements**
   - Pauli grouping for Hamiltonian
   - Reduce measurement shots by 3x
   - Lower noise via fewer measurements

**Expected Result**: 5-10 mHa → <1 mHa for H2, <10 mHa for H2O

---

### Phase 3: Research Frontier (1-2 months) 🔬
5. **Sample-based Krylov (SKQD)**
   - Major architectural change
   - Replace current Krylov-SQD
   - Enable Hadamard test on real hardware
   - Target: best possible accuracy + convergence

6. **Zero-Noise Extrapolation**
   - Add as option to all solvers
   - Combine with other techniques

**Expected Result**: Chemical accuracy (<1 mHa) on systems up to 20 qubits

---

## Projected Performance

### H2 Molecule (4 qubits)
| Method | Error (mHa) | Shots | Runtime |
|--------|-------------|-------|---------|
| Current (no mitigation) | 124.0 | 32K | 18s |
| + Symmetry filtering | 12.0 | 32K | 18s |
| + M3 readout mitigation | 4.0 | 40K | 20s |
| + REM | 0.8 | 80K | 25s |
| + Multi-basis | 0.4 | 30K | 22s |

### LiH Molecule (8 qubits)
| Method | Error (mHa) | Shots | Runtime |
|--------|-------------|-------|---------|
| Baseline (estimated) | ~500 | 100K | 60s |
| Phase 1 complete | ~50 | 100K | 65s |
| Phase 2 complete | ~5 | 120K | 75s |
| Phase 3 complete | <1 | 80K | 70s |

### H2O Molecule (14 qubits)
| Method | Error (mHa) | Shots | Runtime |
|--------|-------------|-------|---------|
| Baseline (estimated) | ~2000 | 500K | 300s |
| Phase 1 complete | ~200 | 500K | 320s |
| Phase 2 complete | ~20 | 400K | 290s |
| Phase 3 complete | ~5 | 300K | 280s |

---

## Cost-Benefit Analysis

### Phase 1 (Symmetry + M3)
- **Implementation time**: 1-2 weeks
- **Expected improvement**: 10-30x error reduction
- **Cost**: <10% overhead
- **ROI**: ⭐⭐⭐⭐⭐ Excellent

### Phase 2 (REM + Multi-basis)
- **Implementation time**: 2-4 weeks
- **Expected improvement**: Additional 5-10x
- **Cost**: 2x shots for calibration, amortized across runs
- **ROI**: ⭐⭐⭐⭐⭐ Excellent

### Phase 3 (SKQD + ZNE)
- **Implementation time**: 1-2 months
- **Expected improvement**: Reaches research frontier
- **Cost**: 3x shots (ZNE), moderate complexity
- **ROI**: ⭐⭐⭐⭐ Good (research impact)

---

## Conclusions

1. **124 mHa error is fixable**: Multiple proven strategies exist
2. **Quick wins available**: Symmetry + M3 can achieve 10-30x improvement in 2 weeks
3. **Chemical accuracy achievable**: <1 mHa for small molecules, <10 mHa for medium
4. **Scaling to larger molecules**: Phase 2-3 essential for molecules beyond 10 qubits
5. **Kanad is well-positioned**: Already has governance system for symmetry filtering

## Next Steps

**Immediate (this week)**:
1. Enhance symmetry filtering in CovalentGovernanceProtocol
2. Integrate M3 mitigation with SQD diagonal measurements
3. Test on H2, measure improvement

**Short-term (2-4 weeks)**:
4. Implement Reference-State Error Mitigation
5. Add multi-basis measurement optimization
6. Test on LiH, H2O

**Medium-term (1-2 months)**:
7. Research and implement SKQD
8. Add ZNE option
9. Publish results, compare to IBM's official SQD addon

---

## References

1. IBM SQD Documentation: https://docs.quantum.ibm.com/guides/qiskit-addons-sqd
2. Sample-based Krylov QD (2025): https://arxiv.org/abs/2501.09702
3. Reference-State Error Mitigation (2023): JCTC doi:10.1021/acs.jctc.2c00807
4. Multireference Error Mitigation (2025): arXiv:2503.05967
5. Efficient Measurements (2020): npj Quantum Information doi:10.1038/s41534-020-00341-7
6. KQD on 56 Qubits: Nature Communications (IBM/U.Tokyo)
