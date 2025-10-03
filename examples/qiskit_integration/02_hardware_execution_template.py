"""
IBM Quantum Hardware Execution Template

This script shows how to run VQE on real IBM Quantum hardware.

PREREQUISITES:
1. IBM Quantum account: https://quantum.ibm.com/
2. Save credentials:
   from qiskit_ibm_runtime import QiskitRuntimeService
   QiskitRuntimeService.save_account(channel='ibm_quantum', token='YOUR_TOKEN')

3. Check available backends:
   service = QiskitRuntimeService()
   print(service.backends())
"""

import numpy as np
from kanad.core.atom import Atom
from kanad.core.representations.base_representation import Molecule
from kanad.core.representations.lcao_representation import LCAORepresentation
from kanad.core.hamiltonians.covalent_hamiltonian import CovalentHamiltonian
from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper
from kanad.ansatze.hardware_efficient_ansatz import RealAmplitudesAnsatz
from kanad.solvers.vqe_solver import VQESolver


def run_on_ibm_hardware(backend_name='ibm_brisbane', use_real_hardware=False):
    """
    Run VQE on IBM Quantum hardware.

    Args:
        backend_name: IBM backend name (e.g., 'ibm_brisbane', 'ibm_kyoto')
        use_real_hardware: If False, runs on simulator for testing
    """

    print("=" * 70)
    print("IBM Quantum Hardware VQE Example")
    print("=" * 70)

    if not use_real_hardware:
        print("\n⚠️  Running in SIMULATION mode")
        print("Set use_real_hardware=True to run on actual quantum hardware")
        backend_name = 'aer_simulator'

    # Create H2 molecule (small enough for NISQ hardware)
    print("\n1. Creating H2 molecule...")
    h1 = Atom('H', position=np.array([0.0, 0.0, 0.0]))
    h2 = Atom('H', position=np.array([0.0, 0.0, 1.4]))
    molecule = Molecule([h1, h2])

    # Build Hamiltonian
    print("2. Building Hamiltonian...")
    representation = LCAORepresentation(molecule)
    hamiltonian = CovalentHamiltonian(
        molecule,
        representation,
        basis_name='sto-3g'
    )

    # Hardware-efficient ansatz (shallow circuit for NISQ)
    print("3. Creating hardware-efficient ansatz...")
    ansatz = RealAmplitudesAnsatz(
        n_qubits=hamiltonian.n_orbitals,
        n_electrons=hamiltonian.n_electrons,
        n_layers=1,  # Shallow for hardware
        entanglement='linear'
    )

    mapper = JordanWignerMapper()

    print(f"   Qubits required: {hamiltonian.n_orbitals}")
    print(f"   Circuit depth: {ansatz.get_circuit_depth()}")
    print(f"   Parameters: {ansatz.get_num_parameters()}")

    # Create VQE solver with hardware backend
    print(f"\n4. Initializing VQE with backend: {backend_name}")

    backend_options = {}
    if use_real_hardware:
        # Hardware-specific options
        backend_options = {
            'optimization_level': 3,  # Maximum transpilation optimization
        }

    solver = VQESolver(
        hamiltonian=hamiltonian,
        ansatz=ansatz,
        mapper=mapper,
        backend=backend_name,
        optimizer='COBYLA',  # Gradient-free (better for noisy hardware)
        max_iterations=50,  # Limit iterations to save queue time
        shots=8192,  # More shots for better statistics
        **backend_options
    )

    print(f"\n{solver}")

    # Run VQE
    print("\n5. Running VQE optimization...")
    if use_real_hardware:
        print("   ⏳ This may take a while due to hardware queue...")
        print("   Check job status at: https://quantum.ibm.com/jobs")

    result = solver.solve()

    # Print results
    print("\n" + "=" * 70)
    print("RESULTS")
    print("=" * 70)
    print(f"Ground state energy: {result['energy']:.8f} Ha")
    print(f"Converged: {result['converged']}")
    print(f"Iterations: {result['iterations']}")
    print(f"Optimizer message: {result['optimizer_message']}")

    # Compare with classical result (for validation)
    print("\n6. Validating with classical simulation...")
    classical_solver = VQESolver(
        hamiltonian=hamiltonian,
        ansatz=ansatz,
        mapper=mapper,
        backend='classical',
        optimizer='COBYLA',
        max_iterations=50
    )

    classical_result = classical_solver.solve(
        initial_parameters=result['parameters']  # Use same final parameters
    )

    print(f"Classical energy: {classical_result['energy']:.8f} Ha")
    error_mha = (result['energy'] - classical_result['energy']) * 1000
    print(f"Error: {error_mha:.4f} mHa")

    if use_real_hardware:
        print("\n✅ Hardware execution complete!")
        print(f"   Energy difference (hardware vs classical): {error_mha:.4f} mHa")
        if abs(error_mha) < 10:
            print("   ✓ Results are in good agreement!")
        else:
            print("   ⚠️  Large error - may need error mitigation or more shots")
    else:
        print("\n✓ Simulation complete!")
        print("   To run on real hardware:")
        print("   1. Set use_real_hardware=True")
        print("   2. Ensure IBM Quantum credentials are configured")

    return result


if __name__ == "__main__":
    # Run in simulation mode by default (safe for testing)
    result = run_on_ibm_hardware(use_real_hardware=False)

    # Uncomment to run on real hardware:
    # result = run_on_ibm_hardware(backend_name='ibm_brisbane', use_real_hardware=True)
