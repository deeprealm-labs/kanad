# 🔍 Hardcoding and Mock Investigation Report

**Date:** Investigation Report  
**Scope:** `kanad/analysis` and `kanad/applications` modules  
**Status:** ⚠️ **MULTIPLE HARDCODED VALUES AND PLACEHOLDERS FOUND**

---

## Executive Summary

The investigation revealed **significant hardcoding and placeholder implementations** across both analysis and application modules. While the infrastructure is functional, many calculations use:
- Hardcoded constants instead of computed values
- Placeholder implementations marked with TODO/FIXME
- Random values for approximations
- Simplified formulas instead of full quantum calculations

**Critical Issues Found:**
- **NMR Calculator**: Uses hardcoded fallback values (-50 ppm constant)
- **Raman Calculator**: Uses approximate formulas instead of quantum density extraction
- **Property Calculator**: Missing density matrix extraction from quantum states
- **Applications**: Multiple placeholder values for formation energies, properties, etc.

---

## 🔴 CRITICAL ISSUES - Analysis Module

### 1. **NMR Calculator** (`kanad/analysis/nmr_calculator.py`)

**Issue:** Hardcoded fallback values when density matrix unavailable

**Location:** Lines 229, 398-401

```python
# Line 229: Placeholder comment
# This is a placeholder for demonstration

# Lines 398-401: Hardcoded fallback
logger.warning("No density matrix available, using approximate method")
n_orbitals = sum(1 if atom.symbol == 'H' else 5 for atom in self.atoms)
rdm1 = np.eye(2 * n_orbitals) * 0.5  # Uniform distribution (hardcoded!)
```

**Problem:**
- Uses uniform density matrix (0.5) when quantum density unavailable
- Falls back to classical approximation instead of extracting from quantum state
- Results in constant -50 ppm chemical shifts (from validation reports)

**Impact:** CRITICAL - Quantum NMR produces incorrect values

---

### 2. **Raman Calculator** (`kanad/analysis/raman_calculator.py`)

**Issue:** Simplified polarizability formulas instead of quantum extraction

**Location:** Lines 166, 171, 263

```python
# Line 166: Hardcoded empirical factor
alpha_iso = n_electrons * 0.8  # Empirical factor (hardcoded!)

# Line 171: Hardcoded anisotropy
alpha_parallel = alpha_iso * (1 + 0.5 * R)  # Along bond (hardcoded 0.5!)

# Line 263: Hardcoded correlation factor
correlation_factor = 1.0 + (abs(correlation_energy) / abs(hf_energy)) * 0.5  # Hardcoded 0.5!
```

**Problem:**
- Uses classical formula `alpha_iso = n_electrons * 0.8` instead of quantum density
- Hardcoded anisotropy factor (0.5) and correlation scaling (0.5)
- Not extracting polarizability from quantum density matrix

**Impact:** CRITICAL - Quantum Raman shows 1500x error vs classical

---

### 3. **Property Calculator** (`kanad/analysis/property_calculator.py`)

**Issue:** Missing density matrix extraction from quantum solvers

**Location:** Lines 748, 763, 843

```python
# Line 748: TODO comment
# TODO: Implement proper density matrix extraction from SQD eigenvector
density_matrix = None  # Will use HF as fallback

# Line 763: TODO comment  
# TODO: Implement proper density matrix extraction from VQE state
density_matrix = None  # Will use HF as fallback

# Line 843: Placeholder note
# This is a placeholder implementation. Full quantum polarizability
# requires computing the response to electric fields at the quantum level
```

**Problem:**
- Density matrix extraction from quantum solvers not implemented
- Falls back to HF density instead of using quantum state
- Polarizability uses classical method even when quantum backend specified

**Impact:** HIGH - Quantum properties not truly quantum

---

### 4. **DOS Calculator** (`kanad/analysis/dos_calculator.py`)

**Issue:** Placeholder PDOS implementation

**Location:** Lines 271, 281

```python
# Line 271: Placeholder comment
# For now, return placeholder

# Line 281: Hardcoded equal distribution
# Placeholder: split DOS equally among atoms
pdos_dict[idx] = total_result['dos'] / n_atoms  # Equal split (hardcoded!)
```

**Problem:**
- Projected DOS (PDOS) splits DOS equally among atoms
- Not using actual orbital projections from PySCF

**Impact:** MEDIUM - PDOS values are incorrect

---

### 5. **Spectroscopy Calculator** (`kanad/analysis/spectroscopy.py`)

**Issue:** Placeholder excited state energies

**Location:** Lines 634, 640, 1022

```python
# Line 634: Hardcoded placeholder energy
energies.append(E_ground + 0.1)  # Placeholder

# Line 640: Hardcoded placeholder energy
energies.append(E_ground + 0.1 * i)  # Placeholder

# Line 1022: TODO comment
# TODO: Implement excited state Hessian calculation
```

**Problem:**
- Uses hardcoded energy offsets (0.1 Ha) when excited states fail
- Excited state Hessian not implemented for vibronic spectra

**Impact:** MEDIUM - Excited state calculations incomplete

---

### 6. **Configuration Explorer** (`kanad/analysis/configuration_explorer.py`)

**Issue:** Placeholder energy calculations

**Location:** Lines 466, 710, 730, 736-737

```python
# Line 466: Placeholder comment
# Placeholder - actual implementation uses VQE/SQD/Hi-VQE

# Line 710: Placeholder comment
# Placeholder - full implementation would:

# Line 730: Placeholder comment
# Placeholder

# Line 736-737: Placeholder geometry optimization
# Placeholder - full implementation uses quantum gradients
logger.debug("Geometry optimization (placeholder)")
```

**Problem:**
- Energy calculation uses fallback methods
- Geometry optimization not fully implemented
- Missing quantum gradient calculations

**Impact:** MEDIUM - Configuration exploration incomplete

---

## 🟡 MODERATE ISSUES - Applications Module

### 7. **Alloy Designer** (`kanad/applications/alloy_designer.py`)

**Issue:** Hardcoded formation energies and properties

**Location:** Lines 320, 501, 522, 532, 600, 603

```python
# Line 320: Placeholder comment
# Compute properties (placeholder - would use quantum solver)

# Line 501: Placeholder comment
# Placeholder - full implementation would use gradient-based optimization

# Line 522: Placeholder comment
# Placeholder - would use actual quantum calculation

# Line 532: Hardcoded formation energy
Delta_H = -5.0  # kJ/mol (placeholder for ternary+)

# Line 600: Hardcoded entropy
S_vib = 0.1  # kJ/(mol·K) (placeholder)

# Line 603: Placeholder comment
# Volume (placeholder)
```

**Problem:**
- Formation energies use hardcoded values (-5.0 kJ/mol)
- Vibrational entropy hardcoded (0.1 kJ/(mol·K))
- Properties computed with placeholders instead of quantum calculations

**Impact:** MEDIUM - Alloy screening uses simplified models

---

### 8. **Catalyst Optimizer** (`kanad/applications/catalyst_optimizer.py`)

**Issue:** Placeholder transition state and properties

**Location:** Lines 356, 472, 474, 514, 593

```python
# Line 356: Placeholder comment
# Placeholder activation energy

# Line 472: Placeholder comment
# Placeholder TS result

# Line 474: Hardcoded geometry
'geometry': np.zeros((10, 3)),  # Placeholder

# Line 514: Placeholder comment
# Placeholder animation

# Line 593: Placeholder comment
# Placeholder - would use actual quantum calculation
```

**Problem:**
- Transition state geometry is zeros array (hardcoded)
- Activation energy calculations use placeholders
- Reaction animation not implemented

**Impact:** MEDIUM - Catalyst optimization incomplete

---

### 9. **Drug Discovery** (`kanad/applications/drug_discovery.py`)

**Issue:** Placeholder toxicity and metabolite predictions

**Location:** Lines 98, 510, 525, 542, 582, 620, 629, 634, 641, 713

```python
# Line 98: Placeholder comment
# Toxicity (placeholder - need ML models)

# Line 510: Hardcoded pose
pose=np.zeros((10, 3)),  # Placeholder for atomic coordinates

# Line 525: Placeholder comment
# Placeholder - would use trained ML model

# Line 542: Placeholder comment
# Placeholder

# Line 582: Placeholder comment
# Placeholder - full implementation would use:

# Line 620: Placeholder comment
# For now, return placeholder

# Line 629: Hardcoded return
return True  # Placeholder

# Line 634: Hardcoded score
score = 0.8  # Placeholder

# Line 641: Hardcoded SMILES
return "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin placeholder

# Line 713: Placeholder comment
# Placeholder - real implementation would use geometry
```

**Problem:**
- Toxicity predictions need ML models (not implemented)
- Metabolite predictions use hardcoded SMILES (Aspirin)
- Binding poses use zeros array
- Many functions return placeholder values

**Impact:** MEDIUM - Drug discovery features incomplete

---

### 10. **Materials Scout** (`kanad/applications/materials_scout.py`)

**Issue:** Random corrections instead of proper predictions

**Location:** Lines 777, 779

```python
# Line 777: Random correction
correction = np.random.normal(0, 0.05)  # ±0.05 eV std

# Line 779: Random correction
correction = np.random.normal(0, 0.2)  # ±0.2 eV std
```

**Problem:**
- Uses random normal distributions for corrections
- Not using actual quantum/ML predictions
- Results vary randomly between runs

**Impact:** MEDIUM - Material screening non-deterministic

---

## 🔵 MINOR ISSUES - Hardcoded Constants

### Physical Constants (Acceptable)
These are standard physical constants and are acceptable:
- NMR gyromagnetic ratios (lines 75-80 in nmr_calculator.py)
- Reference shielding constants (lines 94-101 in nmr_calculator.py)
- Conversion factors (Bohr_to_A, Ha_to_J, etc.)

### Empirical Parameters (Needs Review)
These should be configurable or validated:
- `alpha_iso = n_electrons * 0.8` (raman_calculator.py:166)
- `correlation_factor *= 0.5` (raman_calculator.py:263)
- `Delta_H = -5.0` kJ/mol (alloy_designer.py:532)
- `S_vib = 0.1` kJ/(mol·K) (alloy_designer.py:600)

---

## 📊 Summary Statistics

| Module | Critical Issues | Moderate Issues | Minor Issues | Total |
|--------|----------------|-----------------|--------------|-------|
| **NMR Calculator** | 1 | 0 | 0 | 1 |
| **Raman Calculator** | 1 | 0 | 0 | 1 |
| **Property Calculator** | 1 | 0 | 0 | 1 |
| **DOS Calculator** | 0 | 1 | 0 | 1 |
| **Spectroscopy** | 0 | 1 | 0 | 1 |
| **Config Explorer** | 0 | 1 | 0 | 1 |
| **Alloy Designer** | 0 | 1 | 0 | 1 |
| **Catalyst Optimizer** | 0 | 1 | 0 | 1 |
| **Drug Discovery** | 0 | 1 | 0 | 1 |
| **Materials Scout** | 0 | 1 | 0 | 1 |
| **TOTAL** | **3** | **7** | **0** | **10** |

---

## 🎯 Recommendations

### Priority 1 (CRITICAL - Fix Immediately)
1. **NMR Calculator**: Implement proper density matrix extraction from quantum solvers
2. **Raman Calculator**: Replace hardcoded formulas with quantum density-based calculations
3. **Property Calculator**: Implement density matrix extraction from SQD/VQE states

### Priority 2 (HIGH - Fix Soon)
4. **DOS Calculator**: Implement proper PDOS using orbital projections
5. **Spectroscopy**: Implement proper excited state calculations
6. **Configuration Explorer**: Complete geometry optimization with quantum gradients

### Priority 3 (MEDIUM - Fix When Convenient)
7. **Alloy Designer**: Replace placeholder energies with quantum calculations
8. **Catalyst Optimizer**: Implement proper transition state finding
9. **Drug Discovery**: Implement ML models for toxicity and metabolites
10. **Materials Scout**: Replace random corrections with proper predictions

---

## 📝 Notes

- **Infrastructure is solid**: The framework for quantum calculations exists and works
- **Main issue**: Property extraction from quantum states is incomplete
- **Validation reports confirm**: Quantum NMR and Raman produce incorrect values (from QUANTUM_VALUES_VALIDATION_REPORT.md)
- **Test coverage**: Unit tests pass, but actual quantum values are wrong

---

## 🔗 Related Files

- `QUANTUM_VALUES_VALIDATION_REPORT.md` - Confirms quantum values are incorrect
- `CRITICAL_ISSUES_REPORT.md` - Details quantum implementation issues
- `kanad/analysis/nmr_calculator.py` - NMR implementation
- `kanad/analysis/raman_calculator.py` - Raman implementation
- `kanad/analysis/property_calculator.py` - Property calculations

---

**Report Generated:** Investigation complete
**Next Steps:** Implement proper density matrix extraction from quantum solvers

