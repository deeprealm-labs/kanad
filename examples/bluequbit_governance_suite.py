"""
BlueQubit Governance Suite: diversified bonding types, ansatzes, and solvers

Runs a curated set of medium-to-larger molecules on BlueQubit CPU device,
validates governance, and tries multiple ansatzes/solvers.

Note: VQE integration via framework backend adapter is pending; for now,
we focus on governance validation and circuit metrics, and run small fixed
circuits on BlueQubit CPU as smoke tests.
"""

import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import numpy as np
import time
from dotenv import load_dotenv

from kanad.core.atom import Atom
from kanad.core.molecule import Molecule
from kanad.core.representations.lcao_representation import LCAORepresentation
from kanad.core.representations.second_quantization import SecondQuantizationRepresentation
from kanad.core.hamiltonians.covalent_hamiltonian import CovalentHamiltonian
from kanad.core.hamiltonians.ionic_hamiltonian import IonicHamiltonian
from kanad.core.hamiltonians.metallic_hamiltonian import MetallicHamiltonian
from kanad.ansatze.governance_aware_ansatz import (
    CovalentGovernanceAnsatz,
    IonicGovernanceAnsatz,
    AdaptiveGovernanceAnsatz,
)
from kanad.ansatze.ucc_ansatz import UCCAnsatz


def header(t: str):
    print("\n" + "=" * 80)
    print(t)
    print("=" * 80)


def governance_summary(title: str, validation: dict):
    header(title)
    for check in validation.get('checks', []):
        status = "✓" if check['passed'] else "✗"
        print(f"[{status}] {check['name']}: {check.get('message', '')}")
    if validation.get('all_checks_passed'):
        print("\nGovernance validation PASSED")
    else:
        print("\nGovernance validation has warnings/failures")


def smoke_run_on_bluequbit(circuit_qubits: int, label: str, shots: int = 128):
    token = os.getenv('BLUEQUBIT_API_TOKEN') or os.getenv('TOKEN')
    if not token:
        print(f"\n⚠ BLUEQUBIT_API_TOKEN not set - skipping {label} smoke run")
        return
    try:
        import bluequbit
        os.environ['BLUEQUBIT_API_TOKEN'] = token
        bq = bluequbit.init()

        # Minimal circuit: prepare HF over given qubits
        from qiskit import QuantumCircuit
        qc = QuantumCircuit(circuit_qubits)
        for i in range(min(2, circuit_qubits)):
            qc.x(i)
        for i in range(circuit_qubits):
            qc.ry(0.05, i)

        print(f"\nBlueQubit smoke: {label} ({circuit_qubits} qubits)")
        # Simple retry for transient API issues
        last_err = None
        for attempt in range(3):
            try:
                result = bq.run(qc, device='cpu', shots=shots)
                counts = result.get_counts()
                print(f"  Job: {result.job_id}, states: {len(counts)}, shots: {sum(counts.values())}")
                break
            except Exception as e:
                last_err = e
                print(f"  Attempt {attempt+1} failed: {e}")
                time.sleep(2.0)
        else:
            print(f"  ✗ Smoke run failed after retries: {last_err}")
        time.sleep(1.0)
    except Exception as e:
        print(f"\n✗ BlueQubit smoke error for {label}: {e}")


def run_covalent_systems():
    # CH4 (tetrahedral)
    c = Atom('C', position=np.array([0.0, 0.0, 0.0]))
    Hs = [
        Atom('H', position=np.array([0.629, 0.629, 0.629])),
        Atom('H', position=np.array([0.629, -0.629, -0.629])),
        Atom('H', position=np.array([-0.629, 0.629, -0.629])),
        Atom('H', position=np.array([-0.629, -0.629, 0.629])),
    ]
    mol_ch4 = Molecule([c] + Hs, charge=0, spin=0, basis='sto-3g')
    rep_ch4 = LCAORepresentation(mol_ch4)
    H_ch4 = CovalentHamiltonian(molecule=mol_ch4, representation=rep_ch4, use_governance=True)
    validation = H_ch4.validate_with_governance()
    governance_summary("CH4 Governance", validation)

    # Ansatz choices
    ansatzes = [
        ("CovalentGovernanceAnsatz", CovalentGovernanceAnsatz(2 * H_ch4.n_orbitals, mol_ch4.n_electrons, n_layers=1)),
        ("AdaptiveGovernanceAnsatz", AdaptiveGovernanceAnsatz(2 * H_ch4.n_orbitals, mol_ch4.n_electrons, n_layers=1, bonding_type='auto')),
        ("UCCSD", UCCAnsatz(2 * H_ch4.n_orbitals, mol_ch4.n_electrons, include_singles=True, include_doubles=True)),
    ]
    print(f"\nCH4 Circuit scales:")
    for name, a in ansatzes:
        circuit = a.build_circuit()
        print(f"  {name:26s} qubits={circuit.n_qubits if hasattr(circuit,'n_qubits') else 2*H_ch4.n_orbitals}, gates={len(circuit.gates) if hasattr(circuit,'gates') else 'N/A'}")

    smoke_run_on_bluequbit(2 * H_ch4.n_orbitals, "CH4")

    # NH3 (ammonia)
    n = Atom('N', position=np.array([0.0, 0.0, 0.0]))
    h1 = Atom('H', position=np.array([0.0, 0.0, 1.012]))
    h2 = Atom('H', position=np.array([0.9540, 0.0, -0.3373]))
    h3 = Atom('H', position=np.array([-0.4770, 0.8264, -0.3373]))
    mol_nh3 = Molecule([n, h1, h2, h3], charge=0, spin=0, basis='sto-3g')
    rep_nh3 = LCAORepresentation(mol_nh3)
    H_nh3 = CovalentHamiltonian(molecule=mol_nh3, representation=rep_nh3, use_governance=True)
    validation = H_nh3.validate_with_governance()
    governance_summary("NH3 Governance", validation)

    ansatzes = [
        ("CovalentGovernanceAnsatz", CovalentGovernanceAnsatz(2 * H_nh3.n_orbitals, mol_nh3.n_electrons, n_layers=1)),
        ("AdaptiveGovernanceAnsatz", AdaptiveGovernanceAnsatz(2 * H_nh3.n_orbitals, mol_nh3.n_electrons, n_layers=1, bonding_type='auto')),
        ("UCCSD", UCCAnsatz(2 * H_nh3.n_orbitals, mol_nh3.n_electrons, include_singles=True, include_doubles=True)),
    ]
    print(f"\nNH3 Circuit scales:")
    for name, a in ansatzes:
        circuit = a.build_circuit()
        print(f"  {name:26s} qubits={circuit.n_qubits if hasattr(circuit,'n_qubits') else 2*H_nh3.n_orbitals}, gates={len(circuit.gates) if hasattr(circuit,'gates') else 'N/A'}")

    smoke_run_on_bluequbit(2 * H_nh3.n_orbitals, "NH3")


def run_ionic_systems():
    # LiH (ionic)
    Li = Atom('Li', position=np.array([0.0, 0.0, 0.0]))
    H = Atom('H', position=np.array([1.60, 0.0, 0.0]))
    mol_lih = Molecule([Li, H], charge=0, spin=0)
    rep_lih = SecondQuantizationRepresentation(molecule=mol_lih)
    H_lih = IonicHamiltonian(molecule=mol_lih, representation=rep_lih, use_governance=True)
    validation = H_lih.validate_with_governance()
    governance_summary("LiH Governance", validation)

    ansatzes = [
        ("IonicGovernanceAnsatz", IonicGovernanceAnsatz(2 * H_lih.n_orbitals, H_lih.n_electrons, n_layers=2)),
        ("AdaptiveGovernanceAnsatz", AdaptiveGovernanceAnsatz(2 * H_lih.n_orbitals, H_lih.n_electrons, n_layers=1, bonding_type='ionic')),
        ("UCCSD", UCCAnsatz(2 * H_lih.n_orbitals, H_lih.n_electrons, include_singles=True, include_doubles=True)),
    ]
    print(f"\nLiH Circuit scales:")
    for name, a in ansatzes:
        circuit = a.build_circuit()
        print(f"  {name:26s} qubits={circuit.n_qubits if hasattr(circuit,'n_qubits') else 2*H_lih.n_orbitals}, gates={len(circuit.gates) if hasattr(circuit,'gates') else 'N/A'}")

    smoke_run_on_bluequbit(2 * H_lih.n_orbitals, "LiH")


def run_metallic_systems():
    # Simple 1D metallic chain (4 atoms)
    positions = [
        np.array([0.0, 0.0, 0.0]),
        np.array([1.2, 0.0, 0.0]),
        np.array([2.4, 0.0, 0.0]),
        np.array([3.6, 0.0, 0.0]),
    ]
    atoms = [Atom('Na', position=p) for p in positions]
    mol_metal = Molecule(atoms)
    H_met = MetallicHamiltonian(molecule=mol_metal, lattice_type='1d_chain', hopping_parameter=-1.0, onsite_energy=0.0, hubbard_u=0.0)

    # Report tight-binding size
    print("\nMetallic chain (tight-binding):")
    print(f"  Sites: {H_met.n_sites}, orbitals: {H_met.n_orbitals}")

    # Hardware-efficient circuits are suitable here
    try:
        from kanad.ansatze.hardware_efficient_ansatz import HardwareEfficientAnsatz
        he = HardwareEfficientAnsatz(n_qubits=H_met.n_orbitals, n_electrons=H_met.n_electrons, n_layers=2, entanglement='linear')
        circuit = he.build_circuit()
        print(f"  HardwareEfficientAnsatz: qubits={circuit.n_qubits}, gates={len(circuit.gates)}")
        smoke_run_on_bluequbit(H_met.n_orbitals, "Metal chain")
    except Exception as e:
        print(f"  Note: hardware-efficient ansatz unavailable: {e}")

    # CO2 (linear) — governance via covalent Hamiltonian; capacity check only
    c = Atom('C', position=np.array([0.0, 0.0, 0.0]))
    o1 = Atom('O', position=np.array([1.16, 0.0, 0.0]))
    o2 = Atom('O', position=np.array([-1.16, 0.0, 0.0]))
    mol_co2 = Molecule([c, o1, o2], charge=0, spin=0, basis='sto-3g')
    rep_co2 = LCAORepresentation(mol_co2)
    H_co2 = CovalentHamiltonian(molecule=mol_co2, representation=rep_co2, use_governance=True)
    validation = H_co2.validate_with_governance()
    governance_summary("CO2 Governance", validation)

    print("\nCO2 capacity estimate:")
    q_co2 = 2 * H_co2.n_orbitals
    print(f"  Estimated qubits: {q_co2}")
    if q_co2 > 30:
        print("  ⚠ Likely exceeds BlueQubit ~30-qubit CPU capacity")
    else:
        smoke_run_on_bluequbit(q_co2, "CO2")


def main():
    load_dotenv()
    header("BLUEQUBIT GOVERNANCE SUITE")
    run_covalent_systems()
    run_ionic_systems()
    run_metallic_systems()


if __name__ == "__main__":
    main()


