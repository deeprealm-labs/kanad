"""
Kanad: Governance-Driven Multi-Representation Quantum Chemistry Framework

A quantum chemistry framework for molecular simulation with:
- Multiple bonding types (ionic, covalent, metallic) with governance protocols
- Comprehensive Hamiltonians for all bond types
- Multiple ansatze (UCC, hardware-efficient, governance-aware)
- Multiple solvers (VQE, QPE, SQD, FCI)
- Cloud backends (IBM Quantum, BlueQubit, cuQuantum)
"""

__version__ = "0.1.0"
__author__ = "Kanad Framework Team"

# Core Modules
from kanad.core.molecule import Molecule

# Bond Classes
from kanad.bonds.bond_factory import BondFactory
from kanad.bonds.ionic_bond import IonicBond
from kanad.bonds.covalent_bond import CovalentBond
from kanad.bonds.metallic_bond import MetallicBond

# Hamiltonians
from kanad.core.hamiltonians.molecular_hamiltonian import MolecularHamiltonian
from kanad.core.hamiltonians.ionic_hamiltonian import IonicHamiltonian
from kanad.core.hamiltonians.covalent_hamiltonian import CovalentHamiltonian
from kanad.core.hamiltonians.metallic_hamiltonian import MetallicHamiltonian

# Ansatze
from kanad.ansatze.base_ansatz import BaseAnsatz
from kanad.ansatze.ucc_ansatz import UCCAnsatz
from kanad.ansatze.hardware_efficient_ansatz import HardwareEfficientAnsatz
from kanad.ansatze.governance_aware_ansatz import (
    IonicGovernanceAnsatz,
    CovalentGovernanceAnsatz,
    AdaptiveGovernanceAnsatz
)

# Solvers
from kanad.solvers.vqe_solver import VQESolver
from kanad.solvers.qpe_solver import QPESolver
from kanad.solvers.sqd_solver import SQDSolver
from kanad.solvers.fci_solver import FCISolver

# Backends
from kanad.backends.qiskit_backend import QiskitBackend
from kanad.backends.ibm import IBMRuntimeBackend, IBMVQESolver, IBMQPESolver, IBMSQDSolver
from kanad.backends.bluequbit import BlueQubitBackend

# Governance
from kanad.governance.protocols.base_protocol import BaseGovernanceProtocol, BondingType, GovernanceRule
from kanad.governance.protocols.ionic_protocol import IonicGovernanceProtocol
from kanad.governance.protocols.covalent_protocol import CovalentGovernanceProtocol
from kanad.governance.protocols.metallic_protocol import MetallicGovernanceProtocol

__all__ = [
    # Core
    'Molecule',

    # Bonds
    'BondFactory',
    'IonicBond',
    'CovalentBond',
    'MetallicBond',

    # Hamiltonians
    'MolecularHamiltonian',
    'IonicHamiltonian',
    'CovalentHamiltonian',
    'MetallicHamiltonian',

    # Ansatze
    'BaseAnsatz',
    'UCCAnsatz',
    'HardwareEfficientAnsatz',
    'IonicGovernanceAnsatz',
    'CovalentGovernanceAnsatz',
    'AdaptiveGovernanceAnsatz',

    # Solvers
    'VQESolver',
    'QPESolver',
    'SQDSolver',
    'FCISolver',

    # Backends
    'QiskitBackend',
    'IBMRuntimeBackend',
    'IBMVQESolver',
    'IBMQPESolver',
    'IBMSQDSolver',
    'BlueQubitBackend',

    # Governance
    'BaseGovernanceProtocol',
    'BondingType',
    'GovernanceRule',
    'IonicGovernanceProtocol',
    'CovalentGovernanceProtocol',
    'MetallicGovernanceProtocol',
]
