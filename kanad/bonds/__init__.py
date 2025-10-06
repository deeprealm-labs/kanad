"""
Bond creation and management.

Provides user-friendly interface for creating and analyzing chemical bonds.

This module serves as the PRIMARY INTERFACE for the Kanad framework.
It exposes all core components needed for quantum chemistry calculations:
- Bond creation and management
- Analysis tools (energy, bonding, properties)
- Optimization tools (geometry, circuits, orbitals)
- All necessary core framework components
"""

# Bond Classes (Primary Interface)
from kanad.bonds.bond_factory import BondFactory, BondType
from kanad.bonds.base_bond import BaseBond
from kanad.bonds.ionic_bond import IonicBond
from kanad.bonds.covalent_bond import CovalentBond
from kanad.bonds.metallic_bond import MetallicBond

# Core Framework Components (for solver access)
from kanad.core.atom import Atom
from kanad.core.representations.base_representation import Molecule

# Hamiltonians
from kanad.core.hamiltonians.molecular_hamiltonian import MolecularHamiltonian
from kanad.core.hamiltonians.ionic_hamiltonian import IonicHamiltonian
from kanad.core.hamiltonians.covalent_hamiltonian import CovalentHamiltonian
from kanad.core.hamiltonians.metallic_hamiltonian import MetallicHamiltonian

# Ansatze (all types)
from kanad.ansatze.base_ansatz import BaseAnsatz
from kanad.ansatze.ucc_ansatz import UCCAnsatz
from kanad.ansatze.hardware_efficient_ansatz import HardwareEfficientAnsatz
from kanad.ansatze.governance_aware_ansatz import (
    CovalentGovernanceAnsatz,
    IonicGovernanceAnsatz,
    AdaptiveGovernanceAnsatz
)

# Mappers
from kanad.core.mappers.jordan_wigner_mapper import JordanWignerMapper
from kanad.core.mappers.bravyi_kitaev_mapper import BravyiKitaevMapper
# Note: ParityMapper not implemented, using available mappers

# Analysis Module (full exposure)
from kanad.analysis import (
    EnergyAnalyzer,
    BondingAnalyzer,
    CorrelationAnalyzer,
    PropertyCalculator,
    BondLengthScanner,
    ThermochemistryCalculator,
    FrequencyCalculator,
    UVVisCalculator,
    ExcitedStateSolver,
    DOSCalculator,
    UncertaintyAnalyzer
)

# Optimization Module (full exposure)
from kanad.optimization import (
    QuantumOptimizer,
    OrbitalOptimizer,
    CircuitOptimizer,
    AdaptiveOptimizer,
    GeometryOptimizer
)

__all__ = [
    # Primary Bond Interface
    'BondFactory',
    'BondType',
    'BaseBond',
    'IonicBond',
    'CovalentBond',
    'MetallicBond',

    # Core Components
    'Atom',
    'Molecule',

    # Hamiltonians
    'MolecularHamiltonian',
    'IonicHamiltonian',
    'CovalentHamiltonian',
    'MetallicHamiltonian',

    # Ansatze
    'BaseAnsatz',
    'UCCAnsatz',
    'HardwareEfficientAnsatz',
    'CovalentGovernanceAnsatz',
    'IonicGovernanceAnsatz',
    'AdaptiveGovernanceAnsatz',

    # Mappers
    'JordanWignerMapper',
    'BravyiKitaevMapper',

    # Analysis Tools
    'EnergyAnalyzer',
    'BondingAnalyzer',
    'CorrelationAnalyzer',
    'PropertyCalculator',
    'BondLengthScanner',
    'ThermochemistryCalculator',
    'FrequencyCalculator',
    'UVVisCalculator',
    'ExcitedStateSolver',
    'DOSCalculator',
    'UncertaintyAnalyzer',

    # Optimization Tools
    'QuantumOptimizer',
    'OrbitalOptimizer',
    'CircuitOptimizer',
    'AdaptiveOptimizer',
    'GeometryOptimizer',
]
