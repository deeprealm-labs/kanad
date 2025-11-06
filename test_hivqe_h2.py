"""Test full Hi-VQE workflow with H2"""
import numpy as np
from kanad.bonds import BondFactory
from kanad.core.active_space import get_governance_active_space
from kanad.core.configuration import ConfigurationSubspace
from kanad.core.classical_solver import compute_subspace_energy, get_important_configurations
from kanad.governance.protocols.covalent_protocol import CovalentGovernanceProtocol

print("="*80)
print("HI-VQE FULL WORKFLOW TEST - H2")
print("="*80)

# Step 1: Create H2 molecule
print(f"\n{'─'*80}")
print("Step 1: Molecule Setup")
print(f"{'─'*80}")

bond = BondFactory.create_bond('H', 'H', distance=0.74)
molecule = bond.molecule

print(f"\nH2 Molecule:")
print(f"  Atoms: {molecule.n_atoms}")
print(f"  Electrons: {molecule.n_electrons}")

# Step 2: Get governance-aware active space
print(f"\n{'─'*80}")
print("Step 2: Governance-Aware Active Space")
print(f"{'─'*80}")

protocol = CovalentGovernanceProtocol()
frozen, active, n_active_electrons = get_governance_active_space(molecule, protocol)

print(f"\nActive Space:")
print(f"  Frozen orbitals: {frozen}")
print(f"  Active orbitals: {active}")
print(f"  Active electrons: {n_active_electrons}")
print(f"  Qubits: {len(active) * 2}")

n_qubits = len(active) * 2

# Step 3: Build Hamiltonian
print(f"\n{'─'*80}")
print("Step 3: Build Hamiltonian")
print(f"{'─'*80}")

# Use bond's Hamiltonian and convert to sparse Pauli
hamiltonian = bond.hamiltonian.to_sparse_hamiltonian(mapper='jordan_wigner')

print(f"\nHamiltonian:")
print(f"  Pauli terms: {len(hamiltonian)}")
print(f"  Qubits: {hamiltonian.num_qubits}")

# Step 4: Initialize configuration subspace
print(f"\n{'─'*80}")
print("Step 4: Configuration Subspace Initialization")
print(f"{'─'*80}")

subspace = ConfigurationSubspace(
    n_qubits=n_qubits,
    n_electrons=n_active_electrons,
    protocol=protocol
)

# Start with HF configuration
hf_config = subspace.get_hf_configuration()
subspace.add_config(hf_config)

print(f"\nInitial subspace:")
print(f"  HF configuration: {hf_config}")
print(f"  Subspace size: {len(subspace)}")

# Step 5: Classical diagonalization (HF only)
print(f"\n{'─'*80}")
print("Step 5: Classical Diagonalization (HF only)")
print(f"{'─'*80}")

energy_hf, amplitudes_hf = compute_subspace_energy(hamiltonian, subspace, use_fast=True)

print(f"\nHF Subspace Energy: {energy_hf:.8f} Ha")
print(f"  (This is the HF energy from our configuration-based approach)")

# Step 6: Expand subspace with single excitations
print(f"\n{'─'*80}")
print("Step 6: Subspace Expansion (Single Excitations)")
print(f"{'─'*80}")

single_excs = subspace.generate_single_excitations(hf_config)
added = subspace.add_configs(single_excs)

print(f"\nGenerated {len(single_excs)} single excitations")
print(f"Added {added} to subspace")
print(f"Subspace size: {len(subspace)}")

# Step 7: Classical diagonalization (HF + singles)
print(f"\n{'─'*80}")
print("Step 7: Classical Diagonalization (HF + Singles)")
print(f"{'─'*80}")

energy_cis, amplitudes_cis = compute_subspace_energy(hamiltonian, subspace, use_fast=True)

print(f"\nCIS Energy (HF + Singles): {energy_cis:.8f} Ha")
print(f"  vs HF: {energy_hf:.8f} Ha")
print(f"  Improvement: {energy_hf - energy_cis:.8f} Ha ({(energy_hf - energy_cis)*627.5:.2f} kcal/mol)")

# Step 8: Expand with double excitations
print(f"\n{'─'*80}")
print("Step 8: Subspace Expansion (Double Excitations)")
print(f"{'─'*80}")

double_excs = subspace.generate_double_excitations(hf_config)
added = subspace.add_configs(double_excs)

print(f"\nGenerated {len(double_excs)} double excitations from HF")
print(f"Added {added} to subspace")
print(f"Subspace size: {len(subspace)}")

# Step 9: Final classical diagonalization (HF + singles + doubles)
print(f"\n{'─'*80}")
print("Step 9: Classical Diagonalization (HF + Singles + Doubles)")
print(f"{'─'*80}")

energy_cisd, amplitudes_cisd = compute_subspace_energy(hamiltonian, subspace, use_fast=True)

print(f"\nCISD Energy (HF + S + D): {energy_cisd:.8f} Ha")
print(f"  vs CIS: {energy_cis:.8f} Ha")
print(f"  Improvement: {energy_cis - energy_cisd:.8f} Ha ({(energy_cis - energy_cisd)*627.5:.2f} kcal/mol)")

# Step 10: Analyze amplitudes
print(f"\n{'─'*80}")
print("Step 10: Configuration Analysis")
print(f"{'─'*80}")

important = get_important_configurations(subspace, amplitudes_cisd, threshold=0.05)

print(f"\nImportant configurations (|amplitude| > 0.05):")
for config, amp in important[:10]:  # Show top 10
    print(f"  {config}: {amp:+.6f} (|amp|²={abs(amp)**2:.6f})")

# Step 11: Summary
print(f"\n{'='*80}")
print("HI-VQE SUMMARY")
print(f"{'='*80}")

print(f"\nEnergy Progression:")
print(f"  HF only:         {energy_hf:.8f} Ha")
print(f"  + Singles:       {energy_cis:.8f} Ha (Δ = {energy_hf - energy_cis:.6f} Ha)")
print(f"  + Doubles:       {energy_cisd:.8f} Ha (Δ = {energy_cis - energy_cisd:.6f} Ha)")

print(f"\nSubspace:")
print(f"  Final size: {len(subspace)} configurations")
print(f"  vs Full CI: {2**n_qubits} configurations")
print(f"  Reduction: {2**n_qubits / len(subspace):.1f}x smaller")

print(f"\nKey Hi-VQE Advantage:")
print(f"  ✅ Only 1 Z measurement per iteration (not {len(hamiltonian)} Pauli measurements!)")
print(f"  ✅ Exact energy in subspace (no measurement noise)")
print(f"  ✅ Subspace grows intelligently (physics-guided excitations)")

print(f"\n🎯 Hi-VQE workflow complete!")
