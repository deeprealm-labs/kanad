#!/usr/bin/env python3
"""
Phase 1: Core Foundations - Simple Working Tests
================================================

Tests what actually exists in the codebase, discovering the API as we go.
"""

import numpy as np
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

print("="*80)
print("PHASE 1: CORE FOUNDATIONS - DISCOVERING ACTUAL API")
print("="*80)

results = {'passed': 0, 'failed': 0, 'skipped': 0}

# ==============================================================================
# Test 1: Basic Atom Properties
# ==============================================================================

print("\n[1] Atom Properties")
print("-"*80)

from kanad.core.atom import Atom

try:
    h = Atom('H', position=np.array([0.0, 0.0, 0.0]))
    print(f"  Symbol: {h.symbol}")
    print(f"  Atomic number: {h.atomic_number}")
    print(f"  Atomic mass: {h.atomic_mass} amu")
    print(f"  Position: {h.position}")

    assert h.atomic_number == 1, "H should have Z=1"
    assert abs(h.atomic_mass - 1.008) < 0.1, "H mass ~1.008 amu"

    print(f"  ✅ Atom properties working")
    results['passed'] += 1

except Exception as e:
    print(f"  ❌ Failed: {e}")
    results['failed'] += 1

# ==============================================================================
# Test 2: H2 Hamiltonian Construction
# ==============================================================================

print("\n[2] H2 Covalent Hamiltonian")
print("-"*80)

from kanad.core.representations.base_representation import Molecule
from kanad.core.representations.lcao_representation import LCAORepresentation
from kanad.core.hamiltonians.covalent_hamiltonian import CovalentHamiltonian

try:
    h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
    h2 = Atom('H', position=np.array([0.0, 0.0, 0.74]))
    molecule = Molecule([h1, h2])

    representation = LCAORepresentation(molecule)
    hamiltonian = CovalentHamiltonian(molecule, representation, basis_name='sto-3g')

    print(f"  Orbitals: {hamiltonian.n_orbitals}")
    print(f"  Electrons: {hamiltonian.n_electrons}")
    print(f"  Nuclear repulsion: {hamiltonian.nuclear_repulsion:.6f} Ha")

    # Check h_core exists
    print(f"  h_core shape: {hamiltonian.h_core.shape}")
    print(f"  h_core[0,0]: {hamiltonian.h_core[0,0]:.6f}")

    # Check ERI exists
    print(f"  ERI shape: {hamiltonian.eri.shape}")
    print(f"  ERI[0,0,0,0]: {hamiltonian.eri[0,0,0,0]:.6f}")

    # Basic sanity checks
    assert hamiltonian.n_orbitals == 2, "H2 should have 2 orbitals"
    assert hamiltonian.n_electrons == 2, "H2 should have 2 electrons"
    assert 0.3 < hamiltonian.nuclear_repulsion < 0.4, f"Nuc repulsion should be ~0.378 Ha"

    print(f"  ✅ Hamiltonian construction working")
    results['passed'] += 1

except Exception as e:
    print(f"  ❌ Failed: {e}")
    import traceback
    traceback.print_exc()
    results['failed'] += 1

# ==============================================================================
# Test 3: Jordan-Wigner Mapping
# ==============================================================================

print("\n[3] Jordan-Wigner Mapper")
print("-"*80)

from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper
from kanad.core.hamiltonians.pauli_converter import PauliConverter

try:
    mapper = JordanWignerMapper()

    # Map hamiltonian to qubits
    pauli_op = PauliConverter.to_sparse_pauli_op(
        hamiltonian,
        mapper,
        use_qiskit_nature=True  # Use working implementation
    )

    print(f"  Qubits: {pauli_op.num_qubits}")
    print(f"  Pauli terms: {len(pauli_op)}")

    # Check energy eigenvalue
    H_matrix = pauli_op.to_matrix()
    eigenvalues = np.linalg.eigvalsh(H_matrix)
    ground_state_energy = eigenvalues[0]

    print(f"  Ground state energy: {ground_state_energy:.6f} Ha")
    print(f"  Expected (FCI): ~-1.137 Ha")

    error = abs(ground_state_energy - (-1.137))
    print(f"  Error: {error:.6f} Ha")

    if error < 0.05:
        print(f"  ✅ Mapper producing correct energies!")
        results['passed'] += 1
    else:
        print(f"  ⚠️  Energy off by {error:.3f} Ha (might be basis set or method)")
        results['failed'] += 1

except Exception as e:
    print(f"  ❌ Failed: {e}")
    import traceback
    traceback.print_exc()
    results['failed'] += 1

# ==============================================================================
# Test 4: VQE Solver
# ==============================================================================

print("\n[4] VQE Solver (CPU)")
print("-"*80)

from kanad.ansatze.hardware_efficient_ansatz import RealAmplitudesAnsatz
from kanad.solvers.vqe_solver import VQESolver

try:
    n_qubits = hamiltonian.n_orbitals * 2
    ansatz = RealAmplitudesAnsatz(
        n_qubits=n_qubits,
        n_electrons=hamiltonian.n_electrons,
        n_layers=1
    )

    solver = VQESolver(
        hamiltonian=hamiltonian,
        ansatz=ansatz,
        mapper=mapper,
        backend='aer_simulator',
        max_iterations=20  # Quick test
    )

    result = solver.solve()

    print(f"  VQE Energy: {result['energy']:.6f} Ha")
    print(f"  Iterations: {result['iterations']}")
    print(f"  Converged: {result['converged']}")

    # Check if energy is negative (should be for stable molecule)
    if result['energy'] < 0:
        print(f"  ✅ VQE running and producing negative energy")
        results['passed'] += 1
    else:
        print(f"  ⚠️  VQE energy positive (unexpected)")
        results['failed'] += 1

except Exception as e:
    print(f"  ❌ Failed: {e}")
    import traceback
    traceback.print_exc()
    results['failed'] += 1

# ==============================================================================
# Test 5: Metallic Bond (Band Theory)
# ==============================================================================

print("\n[5] Metallic Bond - Band Theory")
print("-"*80)

from kanad.bonds.metallic_bond import MetallicBond

try:
    # Simple Na chain
    na_atoms = [Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(5)]
    bond = MetallicBond(na_atoms, lattice_type='1d_chain')

    result = bond.compute_energy()

    print(f"  Cohesive energy: {result['energy']:.6f} eV")
    print(f"  Number of bands: {len(result.get('band_energies', []))}")

    # Energy should be negative (cohesive)
    if result['energy'] < 0:
        print(f"  ✅ Band theory working, negative cohesive energy")
        results['passed'] += 1
    else:
        print(f"  ⚠️  Positive energy (unexpected for cohesive)")
        results['failed'] += 1

except Exception as e:
    print(f"  ❌ Failed: {e}")
    import traceback
    traceback.print_exc()
    results['failed'] += 1

# ==============================================================================
# Test 6: GPU Backend (cuQuantum)
# ==============================================================================

print("\n[6] GPU Backend - cuQuantum TensorNet")
print("-"*80)

try:
    from cuquantum import tensornet
    import cupy as cp

    # Get circuit
    circuit = ansatz.build_circuit()
    params = ansatz.initialize_parameters('small_random')

    # Convert to Qiskit and bind
    qiskit_circuit = circuit.to_qiskit()
    param_dict = {p: params[i] for i, p in enumerate(qiskit_circuit.parameters)}
    bound_circuit = qiskit_circuit.assign_parameters(param_dict)

    # Run on GPU
    converter = tensornet.CircuitToEinsum(bound_circuit, backend='cupy')
    expression, operands = converter.state_vector()
    sv = tensornet.contract(expression, *operands)

    print(f"  Statevector shape: {sv.shape}")
    print(f"  Statevector norm: {cp.linalg.norm(sv.reshape(-1)):.6f}")

    assert abs(cp.linalg.norm(sv.reshape(-1)) - 1.0) < 1e-6, "Statevector not normalized"

    print(f"  ✅ GPU backend (cuQuantum) working!")
    results['passed'] += 1

except ImportError:
    print(f"  ⚠️  cuQuantum or CuPy not available - skipping")
    results['skipped'] += 1

except Exception as e:
    print(f"  ❌ Failed: {e}")
    import traceback
    traceback.print_exc()
    results['failed'] += 1

# ==============================================================================
# Test 7: IonicHamiltonian (Known Issues)
# ==============================================================================

print("\n[7] Ionic Hamiltonian (Known Bugs)")
print("-"*80)

from kanad.core.hamiltonians.ionic_hamiltonian import IonicHamiltonian

try:
    li = Atom('Li', position=np.array([0.0, 0.0, 0.0]))
    f = Atom('F', position=np.array([0.0, 0.0, 1.56]))
    molecule_lif = Molecule([li, f])

    representation_lif = LCAORepresentation(molecule_lif)

    # Try to create IonicHamiltonian
    ionic_ham = IonicHamiltonian(molecule_lif, representation_lif)

    print(f"  Orbitals: {ionic_ham.n_orbitals}")
    print(f"  Electrons: {ionic_ham.n_electrons}")
    print(f"  Nuclear repulsion: {ionic_ham.nuclear_repulsion:.6f} Ha")

    # Known issue: energy might be positive
    print(f"  h_core[0,0]: {ionic_ham.h_core[0,0]:.6f}")

    if ionic_ham.h_core[0,0] > 0:
        print(f"  ⚠️  h_core positive (KNOWN BUG - missing electron-nucleus terms)")
        results['failed'] += 1
    else:
        print(f"  ✅ h_core negative (good!)")
        results['passed'] += 1

except ValueError as e:
    if "basis not available" in str(e):
        print(f"  ⚠️  Li/F basis missing (KNOWN ISSUE)")
        results['failed'] += 1
    else:
        raise

except Exception as e:
    print(f"  ❌ Failed: {e}")
    import traceback
    traceback.print_exc()
    results['failed'] += 1

# ==============================================================================
# Summary
# ==============================================================================

print("\n" + "="*80)
print("PHASE 1 SUMMARY")
print("="*80)

total = results['passed'] + results['failed'] + results['skipped']
print(f"\nResults: {results['passed']}/{total} passed")
print(f"  ✅ Passed: {results['passed']}")
print(f"  ❌ Failed: {results['failed']}")
print(f"  ⚠️  Skipped: {results['skipped']}")

if results['failed'] == 0:
    print("\n🎉 All tests passed! Core foundations are solid.")
    sys.exit(0)
else:
    print(f"\n⚠️  {results['failed']} tests failed - need fixes before proceeding")
    sys.exit(1)
