# Basis Sets in Kanad

## Overview

Kanad provides comprehensive basis set support through the `BasisSetRegistry`. Users can specify basis sets when creating molecular systems, allowing full control over computational accuracy vs. speed tradeoffs.

## Available Basis Sets

Kanad supports **50+ basis sets** including:

### Built-in Basis Sets
- **sto-3g**: Minimal basis (3 Gaussians per STO) - fastest, least accurate
- **6-31g**: Split-valence basis - good balance of speed and accuracy

### PySCF Basis Sets (when PySCF is available)

#### Minimal Basis Sets
- `sto-3g`, `sto-6g`, `minao`

#### Split-Valence Basis Sets
- `3-21g`, `6-31g`, `6-311g`

#### Polarization Basis Sets
- `6-31g*`, `6-31g**` (d on heavy atoms, p on H)
- `6-31g(d)`, `6-31g(d,p)` (alternative notation)
- `6-311g*`, `6-311g**`

#### Diffuse Function Basis Sets
- `6-31+g` (diffuse on heavy atoms)
- `6-31++g` (diffuse on all atoms)
- `6-31+g*`, `6-31++g**` (with polarization)

#### Correlation-Consistent Basis Sets
- `cc-pvdz` (double-zeta)
- `cc-pvtz` (triple-zeta)
- `cc-pvqz` (quadruple-zeta)
- `cc-pv5z` (quintuple-zeta)
- `aug-cc-pvdz`, `aug-cc-pvtz`, `aug-cc-pvqz` (augmented with diffuse)

#### Def2 Basis Sets (Ahlrichs)
- `def2-svp` (split-valence + polarization)
- `def2-svpd` (with diffuse)
- `def2-tzvp` (triple-zeta + polarization)
- `def2-tzvpd`, `def2-tzvpp`
- `def2-qzvp` (quadruple-zeta)

#### Periodic System Basis Sets
- `gth-dzvp`, `gth-tzvp` (for periodic boundary conditions)

## Usage

### Basic Usage

```python
from kanad.core.atom import Atom
from kanad.bonds.covalent_bond import CovalentBond

# Create atoms
H1 = Atom('H', position=[0.0, 0.0, 0.0])
H2 = Atom('H', position=[0.74, 0.0, 0.0])

# Create bond with specified basis set
bond = CovalentBond(H1, H2, basis='6-31g')

# Compute energy
result = bond.compute_energy(method='HF')
print(f"Energy with 6-31g: {result['energy']:.6f} Ha")
```

### Comparing Different Basis Sets

```python
from kanad.core.atom import Atom
from kanad.bonds.covalent_bond import CovalentBond

# Create atoms
H1 = Atom('H', position=[0.0, 0.0, 0.0])
H2 = Atom('H', position=[0.74, 0.0, 0.0])

# Compare energies with different basis sets
basis_sets = ['sto-3g', '6-31g', 'cc-pvdz', 'cc-pvtz']

for basis in basis_sets:
    bond = CovalentBond(H1, H2, basis=basis)
    result = bond.compute_energy(method='HF')
    print(f"{basis:15} : {result['energy']:12.8f} Ha")
```

Expected output:
```
sto-3g          :  -1.11675931 Ha
6-31g           :  -1.12675532 Ha
cc-pvdz         :  -1.13152341 Ha
cc-pvtz         :  -1.13362214 Ha
```

### Using Basis Set Registry

```python
from kanad.core.integrals.basis_registry import BasisSetRegistry

# List all available basis sets
available = BasisSetRegistry.list_available_basis_sets()
print(f"Available: {len(available)} basis sets")

# Get information about a specific basis
info = BasisSetRegistry.get_info('cc-pvdz')
print(info)
# Output: {'name': 'cc-pvdz', 'available': True,
#          'description': 'Correlation-consistent double-zeta', 'source': 'PySCF'}

# Get recommendations for different purposes
print(BasisSetRegistry.recommend_basis('minimal'))     # → 'sto-3g'
print(BasisSetRegistry.recommend_basis('general'))     # → '6-31g*'
print(BasisSetRegistry.recommend_basis('accurate'))    # → 'cc-pvtz'
print(BasisSetRegistry.recommend_basis('correlation')) # → 'cc-pvdz'
print(BasisSetRegistry.recommend_basis('diffuse'))     # → 'aug-cc-pvdz'
```

## Validation

The basis set is validated when creating bonds:

```python
from kanad.bonds.covalent_bond import CovalentBond

# This will work
bond = CovalentBond(H1, H2, basis='sto-3g')

# This will raise ValueError
try:
    bond = CovalentBond(H1, H2, basis='invalid-basis')
except ValueError as e:
    print(f"Error: {e}")
    # Output: Basis set 'invalid-basis' not available.
    #         Available basis sets: 3-21g, 4-31g, 6-21g, ...
```

## Choosing a Basis Set

### For Quick Calculations
- **sto-3g**: Fastest, reasonable for qualitative predictions
- Use for: Initial structure optimization, quick screening

### For General Purpose
- **6-31g** or **6-31g***: Good balance of speed and accuracy
- Use for: Most routine calculations, geometry optimizations

### For High Accuracy
- **cc-pvdz** or **cc-pvtz**: Systematically improvable
- Use for: Benchmark calculations, publication-quality results

### For Correlated Methods (MP2, CCSD)
- **cc-pvdz** or higher: Designed for correlation
- Use for: Post-Hartree-Fock methods

### For Anions and Excited States
- **aug-cc-pvdz**: Includes diffuse functions
- Use for: Negative ions, Rydberg states, polarizabilities

### For Periodic Systems
- **gth-dzvp**: Optimized for plane waves
- Use for: Solids, surfaces, crystals

## Energy Convergence

Larger basis sets generally give lower (more negative) energies:

```
Basis Size:    sto-3g < 6-31g < cc-pvdz < cc-pvtz < cc-pvqz
Energy:        -1.117  < -1.127 < -1.132  < -1.134  < -1.135 Ha
Cost:          Fast    → Medium → Slow    → Slower  → Slowest
```

The energy converges toward the **complete basis set (CBS) limit**.

## Basis Set Superposition Error (BSSE)

For weakly bound systems (van der Waals, hydrogen bonds), use **counterpoise correction**:

```python
# TODO: Implement counterpoise correction in future version
# This accounts for BSSE in interaction energies
```

## Element Support

### Built-in Basis Sets
- **sto-3g**: All elements H-Ar (Z=1-18)
- **6-31g**: H, C, N, O, F only (common organic elements)

### PySCF Basis Sets
- Most basis sets support all elements up to Z=86 (Rn)
- Some specialized basis sets may have limited element coverage

**Important**: If using 6-31g with elements like Li, Na, use PySCF integration or switch to sto-3g.

## Performance Notes

### Computational Cost Scaling

For a molecule with N basis functions:
- **Integrals**: O(N⁴) - most expensive step
- **SCF iteration**: O(N³)
- **Memory**: O(N²) for matrices, O(N⁴) for full ERI tensor

Typical sizes:
- **sto-3g**: 1 function per H, 5 functions per C,N,O
- **6-31g**: 2 functions per H, 9 functions per C,N,O
- **cc-pvdz**: 5 functions per H, 14 functions per C,N,O
- **cc-pvtz**: 14 functions per H, 30 functions per C,N,O

### Memory Usage

For H₂O with different basis sets:
- **sto-3g**: 7 functions → ~0.01 MB
- **6-31g**: 13 functions → ~0.1 MB
- **cc-pvdz**: 24 functions → ~1 MB
- **cc-pvtz**: 58 functions → ~30 MB

## Testing Basis Sets

To verify basis set functionality:

```bash
# Run basis set validation suite
python tests/validation/06_basis_set_validation.py
```

This tests:
1. Basis set registry functionality
2. H₂ with multiple basis sets
3. LiH with sto-3g
4. Basis propagation through components
5. Invalid basis rejection
6. Energy comparison

## References

1. Hehre, W. J.; Ditchfield, R.; Pople, J. A. *J. Chem. Phys.* **1972**, 56, 2257. (6-31G basis)
2. Dunning, T. H. *J. Chem. Phys.* **1989**, 90, 1007. (cc-pVDZ basis)
3. Weigend, F.; Ahlrichs, R. *Phys. Chem. Chem. Phys.* **2005**, 7, 3297. (def2 basis sets)

## Future Enhancements

- [ ] Automatic basis set extrapolation (CBS limit)
- [ ] Counterpoise correction for BSSE
- [ ] Effective core potentials (ECPs) for heavy elements
- [ ] Custom basis set definition from files
- [ ] Basis set optimization for specific systems
