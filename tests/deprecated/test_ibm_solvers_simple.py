"""
Simple test to verify IBM solvers can be imported and initialized.
No actual quantum execution - just structure validation.
"""

import numpy as np
from kanad.core.atom import Atom
from kanad.bonds.bond_factory import BondFactory
from kanad.ansatze.ucc_ansatz import UCCAnsatz

print("=" * 80)
print("IBM SOLVERS - STRUCTURE VALIDATION")
print("=" * 80)

# Create H2 molecule
print("\n1. Creating H2 molecule...")
H1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
H2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))
h2_bond = BondFactory.create_bond(H1, H2)
hamiltonian = h2_bond.hamiltonian
print(f"   ✓ H2 molecule: {hamiltonian.n_orbitals} orbitals, {h2_bond.molecule.n_electrons} electrons")

# Create ansatz
print("\n2. Creating UCCSD ansatz...")
ansatz = UCCAnsatz(n_qubits=4, n_electrons=2, include_singles=True, include_doubles=True)
print(f"   ✓ UCCSD ansatz created")

# Test imports
print("\n3. Testing solver imports...")
try:
    from kanad.backends.ibm import IBMVQESolver, IBMQPESolver, IBMSQDSolver, IBMRuntimeBackend
    print("   ✓ All solvers import successfully")
    print("      - IBMRuntimeBackend")
    print("      - IBMVQESolver")
    print("      - IBMQPESolver")
    print("      - IBMSQDSolver")
except ImportError as e:
    print(f"   ❌ Import failed: {e}")
    exit(1)

# Test Hamiltonian conversion
print("\n4. Testing sparse Hamiltonian conversion...")
try:
    sparse_h = hamiltonian.to_sparse_hamiltonian()
    print(f"   ✓ Sparse Hamiltonian: {len(sparse_h)} Pauli terms")
    print(f"   ✓ Compression: {(2**(2*hamiltonian.n_orbitals))**2} matrix elements → {len(sparse_h)} terms")
except Exception as e:
    print(f"   ❌ Conversion failed: {e}")
    exit(1)

# Test circuit creation
print("\n5. Testing ansatz circuit...")
try:
    circuit = ansatz.build_circuit()
    qiskit_circuit = circuit.to_qiskit()
    print(f"   ✓ Circuit: {qiskit_circuit.num_qubits} qubits, {qiskit_circuit.depth()} depth")
    print(f"   ✓ Parameters: {circuit.get_num_parameters()}")
except Exception as e:
    print(f"   ❌ Circuit creation failed: {e}")
    exit(1)

print("\n" + "=" * 80)
print("VALIDATION SUMMARY")
print("=" * 80)
print("✅ All structural components working correctly:")
print("   ✓ Molecule creation")
print("   ✓ Hamiltonian construction")
print("   ✓ Sparse Pauli conversion (fast!)")
print("   ✓ Ansatz circuit generation")
print("   ✓ IBM solver imports")
print()
print("📝 Next step: Run solvers on IBM Quantum")
print("   Use test_ibm_vqe.py with your credentials")
print("=" * 80)
