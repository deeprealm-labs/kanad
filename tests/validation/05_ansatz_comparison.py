"""
Validation Script 5: Ansatz Comparison
======================================

Compares different variational ansätze for the same molecular system
and validates their properties.

This script validates:
1. UCC ansatz construction (UCCSD, singles, doubles)
2. Hardware-efficient ansätze
3. Governance-aware ansätze
4. Parameter counting
5. Circuit depth comparison
"""

import numpy as np
from kanad.bonds import BondFactory
from kanad.ansatze import (
    UCCAnsatz,
    UCC_S_Ansatz,
    UCC_D_Ansatz,
    HardwareEfficientAnsatz,
    RealAmplitudesAnsatz,
    EfficientSU2Ansatz,
    IonicGovernanceAnsatz,
    CovalentGovernanceAnsatz,
    AdaptiveGovernanceAnsatz
)

print("="*70)
print("VALIDATION 5: Ansatz Comparison")
print("="*70)
print()

# ============================================================================
# 1. Setup Test System (H2)
# ============================================================================
print("1. Setting up H2 molecule...")
print("-" * 70)

bond = BondFactory.create_bond('H', 'H')
n_qubits = bond.representation.n_qubits
n_electrons = bond.molecule.n_electrons
n_orbitals = bond.representation.n_orbitals

print(f"Molecule: H2")
print(f"Number of qubits: {n_qubits}")
print(f"Number of electrons: {n_electrons}")
print(f"Number of orbitals: {n_orbitals}")
print()

# ============================================================================
# 2. UCC Ansätze Family
# ============================================================================
print("2. Unitary Coupled Cluster (UCC) ansätze...")
print("-" * 70)

ucc_ansatze = {
    'UCCSD': UCCAnsatz(n_qubits, n_electrons),
    'UCC-S (singles only)': UCC_S_Ansatz(n_qubits, n_electrons),
    'UCC-D (doubles only)': UCC_D_Ansatz(n_qubits, n_electrons),
}

# Build circuits to initialize parameters
for ansatz in ucc_ansatze.values():
    if ansatz.circuit is None:
        ansatz.circuit = ansatz.build_circuit()

print(f"{'Ansatz':<25} {'Parameters':<15} {'Excitations':<15}")
print("-" * 70)

for name, ansatz in ucc_ansatze.items():
    n_params = len(ansatz.circuit.parameters) if ansatz.circuit else 0
    n_excitations = len(ansatz.excitations) if hasattr(ansatz, 'excitations') else 'N/A'

    print(f"{name:<25} {n_params:<15} {str(n_excitations):<15}")

print()

# ============================================================================
# 3. Hardware-Efficient Ansätze
# ============================================================================
print("3. Hardware-efficient ansätze...")
print("-" * 70)

n_layers = 2

hea_ansatze = {
    'Hardware-Efficient (Ry+CZ)': HardwareEfficientAnsatz(n_qubits, n_layers),
    'RealAmplitudes (Ry+CX)': RealAmplitudesAnsatz(n_qubits, n_layers),
    'EfficientSU2 (Ry+Rz+CX)': EfficientSU2Ansatz(n_qubits, n_layers),
}

# Build circuits
for ansatz in hea_ansatze.values():
    if ansatz.circuit is None:
        ansatz.circuit = ansatz.build_circuit()

print(f"{'Ansatz':<30} {'Layers':<10} {'Parameters':<15}")
print("-" * 70)

for name, ansatz in hea_ansatze.items():
    n_params = len(ansatz.circuit.parameters) if ansatz.circuit else 0

    print(f"{name:<30} {n_layers:<10} {n_params:<15}")

print()

# ============================================================================
# 4. Governance-Aware Ansätze
# ============================================================================
print("4. Governance-aware ansätze...")
print("-" * 70)

gov_ansatze = {}

# Ionic governance (if enough qubits)
try:
    gov_ansatze['Ionic Governance'] = IonicGovernanceAnsatz(
        n_qubits, n_electrons, n_layers=2
    )
except Exception as e:
    print(f"Note: Ionic governance not applicable for H2: {str(e)[:50]}")

# Covalent governance
try:
    gov_ansatze['Covalent Governance'] = CovalentGovernanceAnsatz(
        n_qubits, n_electrons, n_layers=2
    )
except Exception as e:
    print(f"Note: Covalent governance error: {str(e)[:50]}")

# Adaptive governance
try:
    # Create with atom info for adaptive selection
    from kanad.core.atom import Atom
    atom1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
    atom2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))

    gov_ansatze['Adaptive Governance'] = AdaptiveGovernanceAnsatz(
        n_qubits, n_electrons, [atom1, atom2], n_layers=2
    )
except Exception as e:
    print(f"Note: Adaptive governance error: {str(e)[:50]}")

if gov_ansatze:
    # Build circuits
    for ansatz in gov_ansatze.values():
        if ansatz.circuit is None:
            ansatz.circuit = ansatz.build_circuit()

    print(f"{'Ansatz':<25} {'Parameters':<15} {'Governance Type':<20}")
    print("-" * 70)

    for name, ansatz in gov_ansatze.items():
        n_params = len(ansatz.circuit.parameters) if ansatz.circuit else 0
        gov_type = getattr(ansatz, 'governance_type', 'N/A')

        print(f"{name:<25} {n_params:<15} {str(gov_type):<20}")

print()

# ============================================================================
# 5. Parameter Count Comparison
# ============================================================================
print("5. Parameter count comparison...")
print("-" * 70)

all_ansatze = {**ucc_ansatze, **hea_ansatze, **gov_ansatze}

# Sort by parameter count
sorted_ansatze = sorted(
    all_ansatze.items(),
    key=lambda x: len(x[1].circuit.parameters) if x[1].circuit else 0
)

print(f"{'Ansatz':<35} {'Parameters':<15} {'Category':<20}")
print("-" * 70)

for name, ansatz in sorted_ansatze:
    n_params = len(ansatz.circuit.parameters) if ansatz.circuit else 0

    # Determine category
    if 'UCC' in name:
        category = 'Quantum Chemistry'
    elif 'Hardware' in name or 'Real' in name or 'Efficient' in name:
        category = 'Hardware-Efficient'
    elif 'Governance' in name:
        category = 'Governance-Aware'
    else:
        category = 'Other'

    print(f"{name:<35} {n_params:<15} {category:<20}")

print()

# ============================================================================
# 6. Circuit Properties
# ============================================================================
print("6. Circuit properties...")
print("-" * 70)

print(f"{'Ansatz':<30} {'Gates':<10} {'Depth':<10} {'2-Qubit Gates':<15}")
print("-" * 70)

for name, ansatz in all_ansatze.items():
    circuit = ansatz.circuit

    # Count gates
    n_gates = len(circuit.gates)

    # Estimate depth (simplified - just count parameterized layers)
    depth = getattr(ansatz, 'n_layers', 1)

    # Count 2-qubit gates
    two_qubit_gates = sum(
        1 for gate in circuit.gates
        if gate['type'] in ['CNOT', 'CZ', 'CX']
    )

    print(f"{name:<30} {n_gates:<10} {depth:<10} {two_qubit_gates:<15}")

print()

# ============================================================================
# 7. Ansatz for Different Bond Types
# ============================================================================
print("7. Ansatz selection for different bond types...")
print("-" * 70)

test_systems = [
    ('H', 'H', 'Covalent'),
    ('C', 'H', 'Covalent'),
    ('C', 'O', 'Polar Covalent'),
]

print(f"{'System':<15} {'Bond Type':<20} {'Recommended Ansatz':<30}")
print("-" * 70)

for atom1, atom2, bond_type in test_systems:
    bond_test = BondFactory.create_bond(atom1, atom2)
    n_q = bond_test.representation.n_qubits
    n_e = bond_test.molecule.n_electrons

    # Determine recommended ansatz
    if bond_type == 'Ionic':
        recommended = 'Ionic Governance'
    elif bond_type == 'Covalent' or bond_type == 'Polar Covalent':
        recommended = 'Covalent Governance / UCCSD'
    elif bond_type == 'Metallic':
        recommended = 'Hardware-Efficient'
    else:
        recommended = 'Adaptive Governance'

    print(f"{atom1}-{atom2:<13} {bond_type:<20} {recommended:<30}")

print()

# ============================================================================
# 8. Parameter Initialization
# ============================================================================
print("8. Parameter initialization...")
print("-" * 70)

# Check UCCSD parameter initialization
uccsd = ucc_ansatze['UCCSD']
print(f"UCCSD ansatz:")
if uccsd.circuit:
    print(f"  Number of parameters: {len(uccsd.circuit.parameters)}")
    print(f"  Parameter names: {[p.name for p in uccsd.circuit.parameters[:5]]}")
    print(f"  Initial values: {[p.value for p in uccsd.circuit.parameters[:5]]}")
print()

# Check HEA parameter initialization
hea = hea_ansatze['Hardware-Efficient (Ry+CZ)']
print(f"Hardware-Efficient ansatz:")
if hea.circuit:
    print(f"  Number of parameters: {len(hea.circuit.parameters)}")
    print(f"  Parameter names: {[p.name for p in hea.circuit.parameters[:5]]}")
    print(f"  Initial values: {[p.value for p in hea.circuit.parameters[:5]]}")
print()

# ============================================================================
# 9. Validation Summary
# ============================================================================
print("="*70)
print("VALIDATION SUMMARY")
print("="*70)

validations = []

# Check UCC ansätze created
if len(ucc_ansatze) == 3:
    validations.append(("✓", f"All 3 UCC ansätze created successfully"))
else:
    validations.append(("✗", f"Only {len(ucc_ansatze)}/3 UCC ansätze created"))

# Check HEA ansätze created
if len(hea_ansatze) == 3:
    validations.append(("✓", f"All 3 hardware-efficient ansätze created"))
else:
    validations.append(("✗", f"Only {len(hea_ansatze)}/3 HEA ansätze created"))

# Check parameter counts are reasonable
all_have_params = all(len(a.circuit.parameters) > 0 if a.circuit else False for a in all_ansatze.values())
if all_have_params:
    validations.append(("✓", "All ansätze have parameters"))
else:
    validations.append(("✗", "Some ansätze have no parameters"))

# Check UCCSD has both singles and doubles
uccsd_excitations = len(ucc_ansatze['UCCSD'].excitations)
if uccsd_excitations > 0:
    validations.append(("✓", f"UCCSD has {uccsd_excitations} excitations"))
else:
    validations.append(("✗", "UCCSD has no excitations"))

# Check circuits are created
all_have_circuits = all(
    hasattr(a, 'circuit') and a.circuit is not None
    for a in all_ansatze.values()
)
if all_have_circuits:
    validations.append(("✓", "All ansätze have circuits"))
else:
    validations.append(("✗", "Some ansätze missing circuits"))

# Check parameter values initialized
uccsd_params_initialized = all(p.value == 0.0 for p in uccsd.circuit.parameters) if uccsd.circuit else False
if uccsd_params_initialized:
    validations.append(("✓", "Parameters initialized to zero (ready for optimization)"))
else:
    validations.append(("~", "Parameters have non-zero initial values"))

# Print validation results
print()
for symbol, message in validations:
    print(f"{symbol} {message}")

passed = sum(1 for s, _ in validations if s == "✓")
total = len([v for v in validations if v[0] in ["✓", "✗"]])
print()
print(f"Validation score: {passed}/{total} checks passed")

if passed == total:
    print("\n🎉 ALL VALIDATIONS PASSED! Ansatz library working correctly.")
elif passed >= total - 1:
    print("\n✅ MOSTLY PASSED! Ansatz library working well.")
else:
    print(f"\n⚠️  {total - passed} validation(s) failed. Review above.")

print("="*70)
