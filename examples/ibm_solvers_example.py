"""
Example: Using IBM Quantum solvers with modular backend architecture.

This demonstrates the clean new structure:
  kanad/backends/ibm/
    ├── backend.py (IBMRuntimeBackend)
    ├── vqe_solver.py (IBMVQESolver)
    ├── qpe_solver.py (IBMQPESolver)
    └── sqd_solver.py (IBMSQDSolver)
"""

import numpy as np
from kanad.core.atom import Atom
from kanad.bonds.bond_factory import BondFactory
from kanad.ansatze.ucc_ansatz import UCCAnsatz

# New modular imports
from kanad.backends.ibm import (
    IBMRuntimeBackend,
    IBMVQESolver,
    IBMQPESolver,
    IBMSQDSolver
)

print("=" * 80)
print("IBM QUANTUM SOLVERS EXAMPLE")
print("=" * 80)

# Create H2 molecule
print("\n1. Creating H2 molecule...")
H1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
H2 = Atom('H', position=np.array([0.74, 0.0, 0.0]))
h2_bond = BondFactory.create_bond(H1, H2)
hamiltonian = h2_bond.hamiltonian

print(f"   Molecule: H2 (0.74 Å)")
print(f"   Electrons: {h2_bond.molecule.n_electrons}")
print(f"   Orbitals: {hamiltonian.n_orbitals}")

# Create ansatz
print("\n2. Creating UCCSD ansatz...")
ansatz = UCCAnsatz(
    n_qubits=4,
    n_electrons=2,
    include_singles=True,
    include_doubles=True
)
print(f"   Ansatz: UCCSD")
print(f"   Qubits: 4")

# Option 1: VQE Solver
print("\n3. VQE Solver Example")
print("-" * 80)
print("""
# Usage:
from kanad.backends.ibm import IBMVQESolver

vqe = IBMVQESolver(
    hamiltonian=hamiltonian,
    ansatz=ansatz,
    backend_name='ibm_torino',  # or 'ibm_simulator_statevector'
    token='your_ibm_token',
    instance='your_crn',
    optimizer='COBYLA',
    max_iterations=100,
    shots=4096,
    optimization_level=3,
    resilience_level=1
)

result = vqe.solve()
print(f"Ground state energy: {result['energy']:.6f} Ha")
""")

# Option 2: QPE Solver
print("\n4. QPE Solver Example")
print("-" * 80)
print("""
# Usage:
from kanad.backends.ibm import IBMQPESolver

qpe = IBMQPESolver(
    hamiltonian=hamiltonian,
    n_ancilla=6,  # Phase precision: 2^(-6)
    backend_name='ibm_torino',
    token='your_ibm_token',
    instance='your_crn',
    shots=8192,
    optimization_level=3,
    resilience_level=2  # Higher for QPE
)

result = qpe.solve()
print(f"Ground state energy: {result['energy']:.6f} Ha")
print(f"Precision: ±{result['precision']:.2e}")
""")

# Option 3: SQD Solver
print("\n5. SQD Solver Example")
print("-" * 80)
print("""
# Usage:
from kanad.backends.ibm import IBMSQDSolver

sqd = IBMSQDSolver(
    hamiltonian=hamiltonian,
    ansatz=ansatz,
    n_samples=10,  # Number of quantum states to sample
    backend_name='ibm_torino',
    token='your_ibm_token',
    instance='your_crn',
    shots=8192,
    optimization_level=3,
    resilience_level=1
)

result = sqd.solve()
print(f"Ground state energy: {result['energy']:.6f} Ha")
print(f"All eigenvalues: {result['eigenvalues']}")
""")

# Backend reuse
print("\n6. Backend Reuse Example")
print("-" * 80)
print("""
# Create backend once, reuse for multiple solvers:
from kanad.backends.ibm import IBMRuntimeBackend, IBMVQESolver, IBMQPESolver

backend = IBMRuntimeBackend(
    backend_name='ibm_torino',
    token='your_token',
    instance='your_crn',
    shots=4096,
    optimization_level=3
)

# Use same backend for different solvers
vqe = IBMVQESolver(hamiltonian, ansatz, backend=backend)
qpe = IBMQPESolver(hamiltonian, n_ancilla=6, backend=backend)

vqe_result = vqe.solve()
qpe_result = qpe.solve()
""")

print("\n" + "=" * 80)
print("MODULAR ARCHITECTURE BENEFITS")
print("=" * 80)
print("✓ Clean organization: kanad/backends/ibm/, /bluequbit/, etc.")
print("✓ Backend-specific optimizations in each solver")
print("✓ Easy to add new providers (IonQ, Braket, etc.)")
print("✓ Backward compatible with existing code")
print("✓ Simple imports: from kanad.backends.ibm import IBMVQESolver")
print("=" * 80)
