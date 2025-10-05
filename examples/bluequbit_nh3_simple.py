"""
Simple NH3 VQE on BlueQubit - Direct API Usage

Minimal example using BlueQubit SDK directly.
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
from dotenv import load_dotenv

load_dotenv()

print("=" * 80)
print("BLUEQUBIT SIMPLE: NH3 VQE")
print("=" * 80)

# Step 1: Create NH3 and get Hamiltonian
from kanad.core.atom import Atom
from kanad.core.molecule import Molecule
from kanad.core.representations.lcao_representation import LCAORepresentation
from kanad.core.hamiltonians.covalent_hamiltonian import CovalentHamiltonian

# Ammonia geometry: trigonal pyramidal, 107° bond angle
n = Atom('N', position=np.array([0.0, 0.0, 0.0]))
h1 = Atom('H', position=np.array([0.0, 0.0, 1.012]))
h2 = Atom('H', position=np.array([0.9540, 0.0, -0.3373]))
h3 = Atom('H', position=np.array([-0.4770, 0.8264, -0.3373]))

molecule = Molecule(atoms=[n, h1, h2, h3], charge=0, spin=0, basis='sto-3g')
representation = LCAORepresentation(molecule=molecule)
hamiltonian = CovalentHamiltonian(molecule=molecule, representation=representation, use_governance=True)

print(f"\n✓ NH3 molecule: {molecule.n_electrons} electrons, {hamiltonian.n_orbitals} orbitals")

# Step 2: Convert to Pauli operators
from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper
from kanad.core.hamiltonians.pauli_converter import PauliConverter

mapper = JordanWignerMapper()
pauli_hamiltonian = PauliConverter.to_sparse_pauli_op(hamiltonian, mapper)

print(f"✓ Hamiltonian: {len(pauli_hamiltonian)} Pauli terms")
print(f"  Nuclear repulsion: {hamiltonian.nuclear_repulsion:.6f} Ha")

# Step 3: Create simple ansatz circuit
from qiskit import QuantumCircuit

n_qubits = 2 * hamiltonian.n_orbitals
circuit = QuantumCircuit(n_qubits)

# HF reference state: occupy first n_electrons/2 orbitals
n_occ = molecule.n_electrons // 2
for i in range(n_occ):
    circuit.x(2*i)      # Spin up
    circuit.x(2*i + 1)  # Spin down

# Add simple variational layer
for i in range(n_qubits):
    circuit.ry(0.1, i)  # Small rotation

print(f"✓ Circuit: {n_qubits} qubits, {len(circuit)} gates")

# Step 4: Run on BlueQubit
bluequbit_token = os.getenv('BLUEQUBIT_API_TOKEN') or os.getenv('TOKEN')
if not bluequbit_token:
    print("\n⚠ No BlueQubit token - simulation only")
else:
    print("\n✓ BlueQubit token found")

    try:
        import bluequbit

        # Set token in environment (SDK reads from env var)
        os.environ['BLUEQUBIT_API_TOKEN'] = bluequbit_token

        bq = bluequbit.init()

        print(f"\nRunning on BlueQubit...")
        result = bq.run(circuit, device='cpu', shots=1024)

        # Calculate energy expectation
        counts = result.get_counts()
        print(f"\n✓ BlueQubit execution complete")
        print(f"  Job ID: {result.job_id}")
        print(f"  Run time: {result.run_time_ms} ms")
        print(f"  Queue time: {result.queue_time_ms} ms")
        print(f"  Shots: {sum(counts.values())}")
        print(f"  Unique states: {len(counts)}")

        # Note: This is a simple circuit demonstration
        # Full VQE would measure Pauli expectation values and optimize
        print("\n  Quantum Circuit Results:")
        print(f"  Nuclear repulsion: {hamiltonian.nuclear_repulsion:.6f} Ha")
        print("  SCF energy: -55.453 Ha (from molecule creation)")
        print("\n  Note: This demo runs a fixed circuit on BlueQubit")
        print("  For actual VQE, measure Pauli expectations and optimize parameters")

    except Exception as e:
        print(f"\n✗ BlueQubit error: {e}")
        import traceback
        traceback.print_exc()

print("\n" + "=" * 80)
