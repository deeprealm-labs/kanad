"""
Example: Ionic Bonding with Governance-Aware Ansätze

Shows how to model ionic systems (NaCl, LiF, etc.) using IonicGovernanceAnsatz
WITHOUT relying on buggy charged molecule Hamiltonians.

The key insight: Ionic character is captured through ANSATZ STRUCTURE,
not through molecular charge!
"""

import numpy as np
from kanad.core.molecule import Molecule
from kanad.core.atom import Atom
from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper
from kanad.ansatze.governance_aware_ansatz import IonicGovernanceAnsatz, CovalentGovernanceAnsatz
from kanad.ansatze.two_local_ansatz import TwoLocalAnsatz
from kanad.solvers.vqe_solver import VQESolver
from pyscf import fci

print("="*70)
print("IONIC BONDING EXAMPLE: LiH")
print("="*70)

print("""
LiH has significant ionic character (Li⁺H⁻-like).
Electronegativity: Li = 0.98, H = 2.20
ΔEN = 1.22 (< 1.7, so partially ionic)

We'll compare three approaches:
1. TwoLocal (hardware-efficient, no physics)
2. CovalentGovernanceAnsatz (wrong physics for LiH!)
3. IonicGovernanceAnsatz (correct physics for LiH)
""")

# Build LiH as NEUTRAL molecule
atoms = [
    Atom('Li', [0.0, 0.0, 0.0]),
    Atom('H', [0.0, 0.0, 1.5])
]

mol = Molecule(atoms, charge=0, spin=0, basis='sto-3g')  # Neutral!
_ = mol.n_electrons

# Get reference energies
hf_energy = mol._hamiltonian.hf_energy
cisolver = fci.FCI(mol._hamiltonian.mf)
fci_energy = cisolver.kernel()[0]
correlation_energy = fci_energy - hf_energy

print(f"\nLiH Reference Energies:")
print(f"  HF:  {hf_energy:.6f} Ha")
print(f"  FCI: {fci_energy:.6f} Ha")
print(f"  Correlation: {correlation_energy*1000:.3f} mHa")

mapper = JordanWignerMapper()
n_qubits = 2 * mol.n_orbitals

# ============================================================================
# Test 1: TwoLocal (baseline, no physics)
# ============================================================================
print(f"\n" + "="*70)
print("TEST 1: TwoLocal Ansatz (no physics)")
print("="*70)

twolocal = TwoLocalAnsatz(
    n_qubits=n_qubits,
    n_electrons=mol.n_electrons,
    n_layers=2,
    rotation_gates='ry',
    entanglement='linear'
)

vqe_twolocal = VQESolver(
    mol._hamiltonian,
    twolocal,
    mapper,
    optimizer='COBYLA',
    max_iterations=150
)

result_twolocal = vqe_twolocal.solve()

corr_captured_twolocal = (result_twolocal['energy'] - hf_energy) / correlation_energy * 100

print(f"  Energy: {result_twolocal['energy']:.6f} Ha")
print(f"  Error from FCI: {(result_twolocal['energy'] - fci_energy)*1000:.3f} mHa")
print(f"  Correlation captured: {corr_captured_twolocal:.1f}%")
print(f"  Parameters: {len(twolocal.parameters)}")

# ============================================================================
# Test 2: CovalentGovernanceAnsatz (wrong physics!)
# ============================================================================
print(f"\n" + "="*70)
print("TEST 2: CovalentGovernanceAnsatz (WRONG physics for ionic LiH)")
print("="*70)

covalent = CovalentGovernanceAnsatz(
    n_qubits=n_qubits,
    n_electrons=mol.n_electrons,
    n_layers=2,
    hybridization='sp'
)

vqe_covalent = VQESolver(
    mol._hamiltonian,
    covalent,
    mapper,
    optimizer='COBYLA',
    max_iterations=150
)

result_covalent = vqe_covalent.solve()

corr_captured_covalent = (result_covalent['energy'] - hf_energy) / correlation_energy * 100

print(f"  Energy: {result_covalent['energy']:.6f} Ha")
print(f"  Error from FCI: {(result_covalent['energy'] - fci_energy)*1000:.3f} mHa")
print(f"  Correlation captured: {corr_captured_covalent:.1f}%")
print(f"  Parameters: {len(covalent.parameters)}")

# ============================================================================
# Test 3: IonicGovernanceAnsatz (CORRECT physics!)
# ============================================================================
print(f"\n" + "="*70)
print("TEST 3: IonicGovernanceAnsatz (CORRECT physics for ionic LiH)")
print("="*70)

ionic = IonicGovernanceAnsatz(
    n_qubits=n_qubits,
    n_electrons=mol.n_electrons,
    n_layers=2
)

# Build circuit with ionization threshold
# ΔEN = 1.22 suggests significant ionic character
circuit_ionic = ionic.build_circuit(ionization_threshold=1.2)

vqe_ionic = VQESolver(
    mol._hamiltonian,
    ionic,
    mapper,
    optimizer='COBYLA',
    max_iterations=150
)

result_ionic = vqe_ionic.solve()

corr_captured_ionic = (result_ionic['energy'] - hf_energy) / correlation_energy * 100

print(f"  Energy: {result_ionic['energy']:.6f} Ha")
print(f"  Error from FCI: {(result_ionic['energy'] - fci_energy)*1000:.3f} mHa")
print(f"  Correlation captured: {corr_captured_ionic:.1f}%")
print(f"  Parameters: {len(ionic.parameters)}")

# ============================================================================
# Summary
# ============================================================================
print(f"\n" + "="*70)
print("SUMMARY")
print("="*70)

results = [
    ('TwoLocal (no physics)', result_twolocal['energy'], len(twolocal.parameters), corr_captured_twolocal),
    ('Covalent (wrong physics)', result_covalent['energy'], len(covalent.parameters), corr_captured_covalent),
    ('Ionic (correct physics)', result_ionic['energy'], len(ionic.parameters), corr_captured_ionic),
]

print(f"\n{'Ansatz':<30} {'Energy (Ha)':<15} {'Error (mHa)':<15} {'Params':<10} {'Corr %'}")
print("-"*85)
print(f"{'FCI (exact)':<30} {fci_energy:<15.6f} {0.0:<15.3f} {'-':<10} {'100.0'}")
for name, energy, n_params, corr_pct in results:
    error = (energy - fci_energy) * 1000
    print(f"{name:<30} {energy:<15.6f} {error:<15.3f} {n_params:<10} {corr_pct:.1f}%")

print(f"\n" + "="*70)
print("KEY INSIGHT")
print("="*70)
print("""
For ionic systems like LiH, NaCl, LiF:
- Use NEUTRAL charge (charge=0)
- Use IonicGovernanceAnsatz to encode ionic physics
- Ansatz structure (localized orbitals, charge-transfer gates)
  captures ionic bonding WITHOUT needing charged molecules

This avoids the Qiskit Nature charged molecule bug AND is more
physically motivated (ionic character emerges from wavefunction
structure, not formal charges).
""")

print("="*70)
