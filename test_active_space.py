"""Test governance-aware active space selection"""
from kanad.bonds import BondFactory
from kanad.core.active_space import get_governance_active_space
from kanad.governance.protocols.covalent_protocol import CovalentGovernanceProtocol

print("="*80)
print("GOVERNANCE-AWARE ACTIVE SPACE TEST")
print("="*80)

# Test 1: H2O (should reduce from 14 to 12 qubits)
print(f"\n{'─'*80}")
print("Test 1: H2O Active Space Reduction")
print(f"{'─'*80}")

# Create H2O molecule
h2o = BondFactory.create_molecule(['O', 'H', 'H'], geometry='water')
print(f"\nH2O Molecule:")
print(f"  Atoms: {h2o.n_atoms}")
print(f"  Total electrons: {h2o.n_electrons}")

# Get governance-guided active space
protocol = CovalentGovernanceProtocol()
frozen, active, n_active_electrons = get_governance_active_space(h2o, protocol)

print(f"\nGovernance-Guided Active Space:")
print(f"  Frozen orbitals: {frozen} (count: {len(frozen)})")
print(f"  Active orbitals: {active} (count: {len(active)})")
print(f"  Active electrons: {n_active_electrons}")
print(f"  Qubits with governance: {len(active) * 2}")

# Verify expectations
assert len(frozen) == 1, f"Expected 1 frozen orbital (O 1s), got {len(frozen)}"
assert len(active) == 6, f"Expected 6 active orbitals, got {len(active)}"
assert len(active) * 2 == 12, f"Expected 12 qubits, got {len(active) * 2}"
assert n_active_electrons == 8, f"Expected 8 active electrons, got {n_active_electrons}"

print(f"\n✅ H2O Test PASSED:")
print(f"   Qubit reduction: 14 → 12 qubits (2 qubit savings)")
print(f"   Core orbital frozen: O 1s")
print(f"   Active orbitals: O 2s, 2p (3) + H 1s (2)")

# Test 2: H2 (no core orbitals to freeze)
print(f"\n{'─'*80}")
print("Test 2: H2 Active Space (No Core Orbitals)")
print(f"{'─'*80}")

bond_h2 = BondFactory.create_bond('H', 'H', distance=0.74)
h2 = bond_h2.molecule

print(f"\nH2 Molecule:")
print(f"  Atoms: {h2.n_atoms}")
print(f"  Total electrons: {h2.n_electrons}")

frozen_h2, active_h2, n_active_electrons_h2 = get_governance_active_space(h2, protocol)

print(f"\nActive Space:")
print(f"  Frozen orbitals: {frozen_h2}")
print(f"  Active orbitals: {active_h2}")
print(f"  Active electrons: {n_active_electrons_h2}")

assert len(frozen_h2) == 0, f"H2 should have no frozen orbitals, got {len(frozen_h2)}"
assert len(active_h2) == 2, f"H2 should have 2 orbitals active, got {len(active_h2)}"
assert n_active_electrons_h2 == 2, f"H2 should have 2 electrons active, got {n_active_electrons_h2}"

print(f"\n✅ H2 Test PASSED: No core orbitals frozen (as expected)")

# Test 3: NH3 (should also freeze N 1s core)
print(f"\n{'─'*80}")
print("Test 3: NH3 Active Space Reduction")
print(f"{'─'*80}")

nh3 = BondFactory.create_molecule(['N', 'H', 'H', 'H'], geometry='auto')
print(f"\nNH3 Molecule:")
print(f"  Atoms: {nh3.n_atoms}")
print(f"  Total electrons: {nh3.n_electrons}")

frozen_nh3, active_nh3, n_active_electrons_nh3 = get_governance_active_space(nh3, protocol)

print(f"\nGovernance-Guided Active Space:")
print(f"  Frozen orbitals: {frozen_nh3}")
print(f"  Active orbitals: {active_nh3}")
print(f"  Active electrons: {n_active_electrons_nh3}")
print(f"  Qubits with governance: {len(active_nh3) * 2}")

assert len(frozen_nh3) == 1, f"Expected 1 frozen orbital (N 1s), got {len(frozen_nh3)}"
assert len(active_nh3) == 7, f"Expected 7 active orbitals, got {len(active_nh3)}"  # N has 4 valence (2s,2p), H has 1 each (3)
assert n_active_electrons_nh3 == 8, f"Expected 8 active electrons, got {n_active_electrons_nh3}"  # 10 total - 2 frozen

total_orbitals = 8  # N: 5 orbitals (1s + 2s + 2p), H: 1 each (3)
qubit_reduction = total_orbitals * 2 - len(active_nh3) * 2
print(f"\n✅ NH3 Test PASSED:")
print(f"   Qubit reduction: {total_orbitals * 2} → {len(active_nh3) * 2} qubits ({qubit_reduction} qubit savings)")

print(f"\n{'='*80}")
print("ALL TESTS PASSED ✅")
print(f"{'='*80}")
print("\nKey Results:")
print(f"  ✅ H2O: 14 → 12 qubits (freeze O 1s)")
print(f"  ✅ H2: All orbitals active (no core)")
print(f"  ✅ NH3: Freeze N 1s (similar to O)")
print(f"\n🎯 Governance protocols successfully guide active space selection!")
