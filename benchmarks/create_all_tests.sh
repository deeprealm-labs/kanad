#!/bin/bash
# Generate all 12 individual test scripts

MOLS=("H2:4:2" "LiH:12:4" "HeH+:4:2")
ANSATZE=("HardwareEfficient" "TwoLocal" "CovalentGov" "IonicGov")

for mol_spec in "${MOLS[@]}"; do
    IFS=':' read -r mol nq ne <<< "$mol_spec"
    
    # Set geometry
    if [ "$mol" = "H2" ]; then
        atoms="[Atom('H', [0,0,0]), Atom('H', [0,0,0.735])]"
        charge=0
    elif [ "$mol" = "LiH" ]; then
        atoms="[Atom('Li', [0,0,0]), Atom('H', [0,0,1.595])]"
        charge=0
    elif [ "$mol" = "HeH+" ]; then
        atoms="[Atom('He', [0,0,0]), Atom('H', [0,0,0.775])]"
        charge=1
    fi
    
    for ansatz in "${ANSATZE[@]}"; do
        filename="test_${mol,,}_${ansatz,,}.py"
        
        # Set ansatz params
        if [ "$ansatz" = "HardwareEfficient" ]; then
            ansatz_class="HardwareEfficientAnsatz"
            ansatz_import="from kanad.ansatze.hardware_efficient_ansatz import HardwareEfficientAnsatz"
            ansatz_init="HardwareEfficientAnsatz(n_qubits=$nq, n_electrons=$ne, n_layers=2, entanglement='linear')"
        elif [ "$ansatz" = "TwoLocal" ]; then
            ansatz_class="TwoLocalAnsatz"
            ansatz_import="from kanad.ansatze.two_local_ansatz import TwoLocalAnsatz"
            ansatz_init="TwoLocalAnsatz(n_qubits=$nq, n_electrons=$ne, n_layers=2)"
        elif [ "$ansatz" = "CovalentGov" ]; then
            ansatz_class="CovalentGovernanceAnsatz"
            ansatz_import="from kanad.ansatze.governance_aware_ansatz import CovalentGovernanceAnsatz"
            ansatz_init="CovalentGovernanceAnsatz(n_qubits=$nq, n_electrons=$ne, n_layers=2)"
        elif [ "$ansatz" = "IonicGov" ]; then
            ansatz_class="IonicGovernanceAnsatz"
            ansatz_import="from kanad.ansatze.governance_aware_ansatz import IonicGovernanceAnsatz"
            ansatz_init="IonicGovernanceAnsatz(n_qubits=$nq, n_electrons=$ne, n_layers=2)"
        fi
        
        cat > "benchmarks/$filename" << EOF
"""Test $mol with $ansatz ansatz"""
import time
import numpy as np
from kanad.core.atom import Atom
from kanad.core.molecule import Molecule
from kanad.core.hamiltonians.openfermion_jw import openfermion_jordan_wigner
from kanad.solvers.vqe_solver import VQESolver
$ansatz_import

mol = Molecule(atoms=$atoms, charge=$charge)

# FCI
mo_e, C = mol.hamiltonian.compute_molecular_orbitals()
h_mo = C.T @ mol.hamiltonian.h_core @ C
eri_mo = np.einsum('pi,qj,pqrs,rk,sl->ijkl', C, C, mol.hamiltonian.eri, C, C)
pauli_h = openfermion_jordan_wigner(h_mo, eri_mo, mol.hamiltonian.nuclear_repulsion, $ne // 2)
E_fci = np.linalg.eigh(pauli_h.to_matrix())[0][0]

# VQE
ansatz = $ansatz_init
solver = VQESolver(hamiltonian=mol.hamiltonian, molecule=mol, ansatz=ansatz, optimizer='SLSQP', max_iterations=150)
t0 = time.time()
result = solver.solve()
dt = time.time() - t0

E_vqe = result['energy']
err = (E_vqe - E_fci) * 1000

print(f"$mol/$ansatz/SLSQP: E={E_vqe:.6f} Ha, error={err:+.3f} mHa, time={dt:.1f}s")
assert abs(err) < 50, f"Error too large: {err} mHa"
print("✓ PASS")
EOF
    done
done

echo "Created 12 test scripts"
