"""
Test to diagnose normalization issue
"""
import numpy as np
from kanad.core.integrals.basis_sets import BasisSet, GaussianPrimitive
from kanad.core.atom import Atom
from kanad.core.integrals.overlap import OverlapIntegrals

# Create two H atoms at H2 bond length
atom1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
atom2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))  # Experimental bond length

# Build basis set
basis = BasisSet('sto-3g')
basis.build_basis([atom1, atom2])

print("="*70)
print("NORMALIZATION DIAGNOSTIC TEST")
print("="*70)
print()

print("H2 molecule at 0.74 Å bond length")
print()

# Check basis functions
print(f"Number of basis functions: {len(basis.basis_functions)}")
print()

# Print primitive data for first H atom
bf1 = basis.basis_functions[0]  # 1s on atom 1
print("Basis function 1 (H1 1s):")
print(f"  Number of primitives: {len(bf1.primitives)}")
for i, prim in enumerate(bf1.primitives):
    print(f"  Primitive {i}:")
    print(f"    Exponent: {prim.exponent:.8f}")
    print(f"    Coefficient: {prim.coefficient:.8f}")
    print(f"    Angular momentum: {prim.angular_momentum}")

    # Test self-overlap (should be 1 if normalized)
    self_overlap = OverlapIntegrals.overlap_primitive(prim, prim)
    print(f"    Self-overlap <φ|φ>: {self_overlap:.8f} (should be ~1)")
print()

# Compute overlap matrix
print("Computing overlap matrix S:")
n_basis = len(basis.basis_functions)
S = np.zeros((n_basis, n_basis))

for i in range(n_basis):
    for j in range(n_basis):
        S[i, j] = OverlapIntegrals.overlap_contracted(
            basis.basis_functions[i],
            basis.basis_functions[j]
        )

print(f"S = ")
print(S)
print()

print("Diagonal elements (should be 1.0):")
for i in range(n_basis):
    print(f"  S[{i},{i}] = {S[i,i]:.8f}")
print()

print("Off-diagonal elements (orbital overlap):")
print(f"  S[0,1] = {S[0,1]:.8f}  (H1-H2 overlap at 0.74 Å)")
print(f"  Expected: ~0.66 for H2 at 0.74 Å")
print()

# Check if issue is in primitive normalization
print("="*70)
print("CHECKING PRIMITIVE NORMALIZATION")
print("="*70)
print()

# Get first primitive of first basis function
prim = bf1.primitives[0]
α = prim.exponent
coeff = prim.coefficient

# Manual normalization check for s-orbital (l=0)
# For a 3D Gaussian: N = (2α/π)^(3/4)
manual_norm = (2 * α / np.pi) ** 0.75

print(f"Exponent α = {α:.8f}")
print(f"Coefficient from STO-3G data = {coeff:.8f}")
print(f"Manual normalization (2α/π)^(3/4) = {manual_norm:.8f}")
print()

# Compute what normalization we're actually using
actual_norm = prim._normalization_constant()
print(f"Actual _normalization_constant() = {actual_norm:.8f}")
print()

# The issue: let's check the evaluate function
center = prim.center
point = center  # Evaluate at center
value_at_center = prim.evaluate(point)
print(f"Value at center (should include normalization): {value_at_center:.8f}")
print()

# Analytical overlap for two identical s-type Gaussians at same point
# <φ|φ> = ∫ φ² dr = N² ∫ exp(-2αr²) dr = N² * (π/(2α))^(3/2)
# For normalized: N² * (π/(2α))^(3/2) = 1
# So: N = (2α/π)^(3/4)

analytical_self_overlap = actual_norm**2 * (np.pi / (2*α))**(3/2)
print(f"Analytical self-overlap with our norm: {analytical_self_overlap:.8f}")
print()

# What coefficient should give us proper normalization?
# In STO-3G, coefficients are contraction coefficients
# The contraction c_i * φ_i needs to be normalized
# For a contracted function: <Σ c_i φ_i | Σ c_j φ_j> = 1

print("="*70)
print("CHECKING CONTRACTED FUNCTION")
print("="*70)
print()

# Overlap of contracted function with itself
bf1_self = OverlapIntegrals.overlap_contracted(bf1, bf1)
print(f"Contracted basis function self-overlap: {bf1_self:.8f}")
print(f"Should be 1.0 for normalized basis function")
print()

# Let's see what the STO-3G contraction should give
print("STO-3G contraction details:")
total = 0.0
for i, prim_i in enumerate(bf1.primitives):
    for j, prim_j in enumerate(bf1.primitives):
        overlap_ij = OverlapIntegrals.overlap_primitive(prim_i, prim_j)
        contrib = prim_i.coefficient * prim_j.coefficient * overlap_ij
        total += contrib
        print(f"  c[{i}] * c[{j}] * S[{i},{j}] = {prim_i.coefficient:.6f} * {prim_j.coefficient:.6f} * {overlap_ij:.6f} = {contrib:.6f}")

print(f"  Total = {total:.8f}")
print()

print("DIAGNOSIS:")
if abs(S[0,0] - 1.0) > 0.1:
    print("❌ PROBLEM: Diagonal elements >> 1.0")
    print("   Root cause: Primitives are not properly normalized")
    print("   The normalization constant calculation is incorrect")
else:
    print("✓ Diagonal elements close to 1.0")
