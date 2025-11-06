"""
Diagnose Hi-VQE Subspace Selection

Understand why single-iteration Hi-VQE gives -1.116 Ha
instead of exact -1.137 Ha, and what configurations are needed.
"""
import numpy as np
from kanad.bonds import BondFactory
from kanad.core.configuration import ConfigurationSubspace
from kanad.governance.protocols.covalent_protocol import CovalentGovernanceProtocol
from kanad.core.classical_solver import SubspaceHamiltonianBuilder

print("="*80)
print("HI-VQE SUBSPACE ANALYSIS - H2 MOLECULE")
print("="*80)

# Setup
bond = BondFactory.create_bond('H', 'H', distance=0.74)
hamiltonian = bond.hamiltonian.to_sparse_hamiltonian(mapper='jordan_wigner')

protocol = CovalentGovernanceProtocol()
subspace = ConfigurationSubspace(4, 2, protocol=protocol)

# Add HF configuration
hf_config = subspace.get_hf_configuration()
subspace.add_config(hf_config)

print(f"\n📌 Step 1: HF Configuration Only")
print(f"   Subspace size: {len(subspace)}")
print(f"   Config: {hf_config}")

# Build Hamiltonian for HF only
builder = SubspaceHamiltonianBuilder(hamiltonian)
H_sub = builder.project_fast(subspace)
eigenvalues = np.linalg.eigvalsh(H_sub)
print(f"   Ground energy: {eigenvalues[0]:.8f} Ha")
print(f"   Expected (FCI): -1.13728383 Ha")
print(f"   Error: {abs(eigenvalues[0] - (-1.13728383))*1000:.4f} mHa")

# Add single excitations
single_excs = subspace.generate_single_excitations(hf_config)
subspace.add_configs(single_excs)

print(f"\n📌 Step 2: HF + Single Excitations")
print(f"   Subspace size: {len(subspace)}")
for i, config in enumerate(subspace):
    print(f"   Config {i}: {config}")

# Build Hamiltonian
H_sub = builder.project_fast(subspace)
eigenvalues = np.linalg.eigvalsh(H_sub)
eigenvectors = np.linalg.eigh(H_sub)[1]
ground_state = eigenvectors[:, 0]

print(f"\n   Ground energy: {eigenvalues[0]:.8f} Ha")
print(f"   Expected (FCI): -1.13728383 Ha")
print(f"   Error: {abs(eigenvalues[0] - (-1.13728383))*1000:.4f} mHa")

print(f"\n   Ground state composition:")
for i, (config, amp) in enumerate(zip(subspace, ground_state)):
    prob = amp**2
    if abs(prob) > 0.01:
        print(f"     {config}: {amp:+.6f} (prob={prob:.6f})")

# Add double excitations
print(f"\n📌 Step 3: Adding Double Excitations")
double_excs = subspace.generate_double_excitations(hf_config)
print(f"   Generated {len(double_excs)} double excitations")

subspace.add_configs(double_excs)
print(f"   Subspace size: {len(subspace)}")

# Build Hamiltonian
H_sub = builder.project_fast(subspace)
eigenvalues = np.linalg.eigvalsh(H_sub)
eigenvectors = np.linalg.eigh(H_sub)[1]
ground_state = eigenvectors[:, 0]

print(f"\n   Ground energy: {eigenvalues[0]:.8f} Ha")
print(f"   Expected (FCI): -1.13728383 Ha")
print(f"   Error: {abs(eigenvalues[0] - (-1.13728383))*1000:.4f} mHa")

print(f"\n   Ground state composition:")
for i, (config, amp) in enumerate(zip(subspace, ground_state)):
    prob = amp**2
    if abs(prob) > 0.01:
        print(f"     {config}: {amp:+.6f} (prob={prob:.6f})")

# Now test the diagonal energies for HF + single excitations
print(f"\n📌 Step 4: Diagonal Energy Analysis (HF + Singles)")
print(f"   This is what we measure on quantum hardware:")

# Reset to HF + singles
subspace2 = ConfigurationSubspace(4, 2, protocol=protocol)
subspace2.add_config(hf_config)
subspace2.add_configs(single_excs)

# Compute diagonal energies
H_sub2 = builder.project_fast(subspace2)
diagonal_energies = [H_sub2[i, i] for i in range(len(subspace2))]

print(f"\n   Diagonal energies ⟨config|H|config⟩:")
for config, diag_e in zip(subspace2, diagonal_energies):
    print(f"     {config}: {diag_e:.8f} Ha")

# Show off-diagonal elements
print(f"\n   Off-diagonal elements ⟨i|H|j⟩:")
for i in range(len(subspace2)):
    for j in range(i+1, len(subspace2)):
        off_diag = H_sub2[i, j]
        if abs(off_diag) > 0.001:
            print(f"     ⟨{subspace2.configurations[i]}|H|{subspace2.configurations[j]}⟩ = {off_diag:.6f} Ha")

# Key insight: The off-diagonal elements are what enable mixing
# Without them, we just get the HF diagonal energy
print(f"\n🔍 KEY INSIGHT:")
print(f"   The diagonal energy of HF is: {diagonal_energies[0]:.8f} Ha")
print(f"   But with off-diagonal mixing, we get: {eigenvalues[0]:.8f} Ha")
print(f"   The off-diagonal coupling lowers the energy!")

print(f"\n💡 CONCLUSION:")
print(f"   - With HF + singles (4 configs): {eigenvalues[0]:.8f} Ha (20.5 mHa error)")
print(f"   - With HF + singles + doubles: {eigenvalues[0]:.8f} Ha")
print(f"   - Literature Hi-VQE achieves 0.08-1.20 mHa error")
print(f"   - Need to include more important configurations (iterative selection)")

print("="*80)
