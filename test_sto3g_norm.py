"""
Test STO-3G coefficient normalization

The issue: STO-3G coefficients from literature are pre-normalized
for a specific convention. We need to renormalize them for our
Gaussian primitives.
"""
import numpy as np

# STO-3G data for H (from our code)
sto3g_H = {
    's': [
        (3.42525091, 0.15432897),
        (0.62391373, 0.53532814),
        (0.16885540, 0.44463454)
    ]
}

print("="*70)
print("STO-3G COEFFICIENT ANALYSIS")
print("="*70)
print()

# The STO-3G coefficients are designed such that when you compute:
# χ(r) = Σ d_i × (2α_i/π)^(3/4) × exp(-α_i r²)
#
# This is normalized.
#
# But we're computing:
# χ(r) = Σ c_i × N_i × exp(-α_i r²)
#
# Where N_i = (2α_i/π)^(3/4) is our normalization
# And c_i are the coefficients
#
# So we need: c_i × N_i = d_i
# Therefore: c_i = d_i / N_i

print("STO-3G H 1s orbital:")
print()

total_norm = 0.0
for i, (exp, d_coeff) in enumerate(sto3g_H['s']):
    # Our normalization
    N_i = (2 * exp / np.pi) ** 0.75

    # The coefficient we should use
    c_i = d_coeff / N_i

    # But wait - let's check if the d_coefficients are already meant
    # to give a normalized contraction

    print(f"Primitive {i}:")
    print(f"  α = {exp:.8f}")
    print(f"  d (from STO-3G) = {d_coeff:.8f}")
    print(f"  N (our norm) = {N_i:.8f}")
    print(f"  d/N = {d_coeff/N_i:.8f}")
    print()

# The question: are the d coefficients normalized such that:
# Σᵢⱼ dᵢ dⱼ Sᵢⱼ = 1  where Sᵢⱼ = ∫ φᵢ φⱼ dr
#
# Let's compute the overlap matrix for normalized primitives

print("="*70)
print("OVERLAP MATRIX FOR NORMALIZED PRIMITIVES")
print("="*70)
print()

def gaussian_overlap_ss(alpha_a, alpha_b, R_AB):
    """
    Overlap between two s-type Gaussians.

    <φ_a|φ_b> = (π/(α_a + α_b))^(3/2) × exp(-α_a α_b R² / (α_a + α_b))

    For normalized Gaussians: multiply by N_a × N_b
    """
    gamma = alpha_a + alpha_b
    K = np.exp(-alpha_a * alpha_b * R_AB**2 / gamma)
    return (np.pi / gamma)**(3/2) * K

# For H-H at same position (self-overlap)
R = 0.0  # same atom

S_unnorm = np.zeros((3, 3))
S_norm = np.zeros((3, 3))

for i in range(3):
    for j in range(3):
        alpha_i = sto3g_H['s'][i][0]
        alpha_j = sto3g_H['s'][j][0]

        S_unnorm[i,j] = gaussian_overlap_ss(alpha_i, alpha_j, R)

        # With our normalization
        N_i = (2 * alpha_i / np.pi) ** 0.75
        N_j = (2 * alpha_j / np.pi) ** 0.75
        S_norm[i,j] = N_i * N_j * S_unnorm[i,j]

print("Overlap matrix (unnormalized primitives):")
print(S_unnorm)
print()

print("Overlap matrix (normalized primitives):")
print(S_norm)
print()

# Now compute contracted overlap with d coefficients
d_coeffs = np.array([d for (_, d) in sto3g_H['s']])

overlap_with_d = d_coeffs.T @ S_unnorm @ d_coeffs
overlap_with_d_normalized = d_coeffs.T @ S_norm @ d_coeffs

print(f"Contracted overlap (d coefficients, unnormalized primitives): {overlap_with_d:.8f}")
print(f"Contracted overlap (d coefficients, normalized primitives): {overlap_with_d_normalized:.8f}")
print()

print("CONCLUSION:")
print("The STO-3G d coefficients are designed to work with UNNORMALIZED primitives!")
print("When we normalize our primitives, we must renormalize the contraction coefficients.")
print()

# The solution: renormalize the contraction
# We need: Σᵢⱼ c_i c_j <φ_i|φ_j> = 1
# Where <φ_i|φ_j> uses our normalized primitives

print("="*70)
print("RENORMALIZING CONTRACTION COEFFICIENTS")
print("="*70)
print()

# Start with d coefficients
c_unnormalized = d_coeffs.copy()

# Compute normalization factor for the contraction
N_contraction = np.sqrt(c_unnormalized.T @ S_norm @ c_unnormalized)
print(f"Contraction normalization factor: {N_contraction:.8f}")

# Renormalized coefficients
c_normalized = c_unnormalized / N_contraction

print(f"\nRenormalized contraction coefficients:")
for i, c in enumerate(c_normalized):
    print(f"  c[{i}] = {c:.8f} (was {c_unnormalized[i]:.8f})")

# Verify
final_overlap = c_normalized.T @ S_norm @ c_normalized
print(f"\nFinal contracted overlap: {final_overlap:.8f}")
print("✓ Should be 1.0!" if abs(final_overlap - 1.0) < 1e-6 else "✗ Still wrong")
