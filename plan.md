# **Kanad: Governance-Driven Multi-Representation Quantum Chemistry Framework**
## **Complete Architectural Blueprint**

A quantum chemistry framework that allows you to create molecules with different different kind of bondings, bonding is governed by protocol with bonding representative 

---

## **PHASE 0: Scientific Foundations & Theory**

### **Core Physics Principles**

#### **1. Bonding Representations (Mathematical Foundation)**

```python
"""
IONIC BONDING
- Physics: Electron transfer between atoms
- Wavefunction: |ψ⟩ = a†_i a_j |0⟩  (localized transfer)
- Operator: T_ij = a†_i a_j
- Entanglement: Minimal (2-qubit)
"""

"""
COVALENT BONDING  
- Physics: Orbital hybridization (LCAO)
- Wavefunction: |ψ_bond⟩ = α|φ_A⟩ + β|φ_B⟩
- Operator: H = Σ_AB (α*α c†_A c_A + β*β c†_B c_B + αβ* c†_A c_B + α*β c†_B c_A)
- Entanglement: 2-qubit Bell-like states
"""

"""
METALLIC BONDING
- Physics: Delocalized band structure
- Wavefunction: |ψ_k⟩ = (1/√N) Σ_j e^(ikj) |j⟩  (Bloch states)
- Operator: H_k = Σ_k ε_k c†_k c_k
- Entanglement: Multi-qubit GHZ-like states
"""
```

---

## **PHASE 1: Core Architecture**

### **Directory Structure**

```
kanad/
├── core/
│   ├── constants/
│   │   ├── __init__.py
│   │   ├── physical_constants.py      # ℏ, e, m_e, etc.
│   │   ├── atomic_data.py             # Periodic table, radii, EN
│   │   └── conversion_factors.py      # Hartree, eV, Bohr, etc.
│   │
│   ├── integrals/
│   │   ├── __init__.py
│   │   ├── one_electron.py            # Kinetic + nuclear attraction
│   │   ├── two_electron.py            # Electron-electron repulsion
│   │   ├── overlap.py                 # ⟨φ_i|φ_j⟩ integrals
│   │   └── basis_sets.py              # STO-3G, 6-31G, etc.
│   │
│   ├── representations/
│   │   ├── __init__.py
│   │   ├── base_representation.py     # Abstract base class
│   │   ├── second_quantization.py     # Fock space (ionic)
│   │   ├── lcao_representation.py     # Hybrid orbitals (covalent)
│   │   └── bloch_representation.py    # k-space (metallic)
│   │
│   ├── hamiltonians/
│   │   ├── __init__.py
│   │   ├── molecular_hamiltonian.py   # H = T + V_ne + V_ee
│   │   ├── ionic_hamiltonian.py       # Transfer Hamiltonian
│   │   ├── covalent_hamiltonian.py    # Hybridization Hamiltonian
│   │   └── metallic_hamiltonian.py    # Band Hamiltonian
│   │
│   └── mappers/
│       ├── __init__.py
│       ├── base_mapper.py             # Abstract mapper interface
│       ├── jordan_wigner.py           # For ionic (localized)
│       ├── hybrid_orbital_mapper.py   # For covalent (paired)
│       └── momentum_space_mapper.py   # For metallic (delocalized)
│
├── governance/
│   ├── __init__.py
│   ├── protocols/
│   │   ├── __init__.py
│   │   ├── base_protocol.py           # Abstract governance
│   │   ├── ionic_protocol.py          # Transfer rules
│   │   ├── covalent_protocol.py       # Hybridization rules
│   │   └── metallic_protocol.py       # Delocalization rules
│   │
│   ├── operators/
│   │   ├── __init__.py
│   │   ├── ionic_operators.py         # a†_i a_j gates
│   │   ├── covalent_operators.py      # Hybridization gates
│   │   └── metallic_operators.py      # Band formation gates
│   │
│   ├── constraints/
│   │   ├── __init__.py
│   │   ├── symmetry_constraints.py    # Point group, spin, etc.
│   │   ├── physical_constraints.py    # Pauli exclusion, fermion parity
│   │   └── bonding_constraints.py     # Bond-specific rules
│   │
│   └── engine.py                      # Governance execution engine
│
├── ansatze/
│   ├── __init__.py
│   ├── base_ansatz.py                 # Abstract ansatz class
│   ├── ionic_ansatz.py                # Minimal entanglement
│   ├── covalent_ansatz.py             # Bell-pair based
│   ├── metallic_ansatz.py             # Highly entangled
│   └── adaptive_ansatz.py             # Governance-driven selection
│
├── analysis/
│   ├── __init__.py
│   ├── energy_estimation.py           # VQE, QPE protocols
│   ├── bonding_analysis.py            # Bond orders, overlap, etc.
│   ├── spectroscopy.py                # Excitation energies
│   └── visualization.py               # Orbital plots, band diagrams
│
├── bonds/
│   ├── __init__.py
│   ├── bond_factory.py                # Bond creation interface
│   ├── ionic_bond.py                  # Ionic bond class
│   ├── covalent_bond.py               # Covalent bond class
│   └── metallic_bond.py               # Metallic bond class
│
├── research/
│   ├── __init__.py
│   ├── workflows/                     # Pre-built research pipelines
│   ├── experiments/                   # Custom experiment templates
│   └── notebooks/                     # Tutorial notebooks
│
└── tests/
    ├── unit/
    ├── integration/
    └── validation/                    # Against known results
```

---

## **PHASE 2: Implementation Roadmap**

### **Step 1: Physical Constants & Atomic Data** 

```python
# kanad/core/constants/physical_constants.py

from typing import Final
import numpy as np
from dataclasses import dataclass

@dataclass(frozen=True)
class PhysicalConstants:
    """Fundamental physical constants in SI units."""
    
    # Universal constants
    SPEED_OF_LIGHT: Final[float] = 299792458.0  # m/s
    PLANCK_CONSTANT: Final[float] = 6.62607015e-34  # J⋅s
    HBAR: Final[float] = 1.054571817e-34  # J⋅s
    
    # Atomic constants
    ELEMENTARY_CHARGE: Final[float] = 1.602176634e-19  # C
    ELECTRON_MASS: Final[float] = 9.1093837015e-31  # kg
    PROTON_MASS: Final[float] = 1.67262192369e-27  # kg
    BOHR_RADIUS: Final[float] = 5.29177210903e-11  # m
    
    # Atomic units (Hartree atomic units)
    HARTREE_TO_EV: Final[float] = 27.211386245988  # eV
    BOHR_TO_ANGSTROM: Final[float] = 0.529177210903  # Å
    
    # Energy conversion
    EV_TO_JOULE: Final[float] = 1.602176634e-19  # J
    KCAL_MOL_TO_HARTREE: Final[float] = 0.0015936  # E_h

# Singleton instance
CONSTANTS = PhysicalConstants()
```

```python
# kanad/core/constants/atomic_data.py

from enum import Enum
from dataclasses import dataclass
from typing import Dict

@dataclass
class AtomicProperties:
    """Properties of an atomic element."""
    symbol: str
    atomic_number: int
    atomic_mass: float  # amu
    covalent_radius: float  # Angstrom
    electronegativity: float  # Pauling scale
    valence_electrons: int
    ionization_energy: float  # eV

class PeriodicTable:
    """Complete periodic table with bonding-relevant properties."""
    
    ELEMENTS: Dict[str, AtomicProperties] = {
        'H': AtomicProperties('H', 1, 1.008, 0.31, 2.20, 1, 13.598),
        'C': AtomicProperties('C', 6, 12.011, 0.76, 2.55, 4, 11.260),
        'N': AtomicProperties('N', 7, 14.007, 0.71, 3.04, 5, 14.534),
        'O': AtomicProperties('O', 8, 15.999, 0.66, 3.44, 6, 13.618),
        'Na': AtomicProperties('Na', 11, 22.990, 1.66, 0.93, 1, 5.139),
        'Cl': AtomicProperties('Cl', 17, 35.45, 1.02, 3.16, 7, 12.968),
        'Fe': AtomicProperties('Fe', 26, 55.845, 1.32, 1.83, 8, 7.902),
        # Add all elements...
    }
    
    @classmethod
    def get_element(cls, symbol: str) -> AtomicProperties:
        return cls.ELEMENTS[symbol]
    
    @classmethod
    def electronegativity_difference(cls, sym1: str, sym2: str) -> float:
        """Calculate EN difference for bonding classification."""
        en1 = cls.ELEMENTS[sym1].electronegativity
        en2 = cls.ELEMENTS[sym2].electronegativity
        return abs(en1 - en2)
```

---

### **Step 2: Integral Computation Engine** 

```python
# kanad/core/integrals/one_electron.py

import numpy as np
from typing import Tuple, List
from ..constants import CONSTANTS

class OneElectronIntegrals:
    """Compute kinetic and nuclear attraction integrals."""
    
    def __init__(self, molecule: 'Molecule', basis_set: str = 'sto-3g'):
        self.molecule = molecule
        self.basis = BasisSet(basis_set)
        self.n_orbitals = len(molecule.atoms) * self.basis.n_functions
        
    def compute_overlap(self) -> np.ndarray:
        """
        Compute overlap matrix S_ij = ⟨φ_i|φ_j⟩
        
        Returns:
            S: (n_orbitals, n_orbitals) overlap matrix
        """
        S = np.zeros((self.n_orbitals, self.n_orbitals))
        
        for i in range(self.n_orbitals):
            for j in range(self.n_orbitals):
                S[i, j] = self._overlap_integral(i, j)
        
        return S
    
    def compute_kinetic(self) -> np.ndarray:
        """
        Compute kinetic energy matrix T_ij = ⟨φ_i|-½∇²|φ_j⟩
        
        Returns:
            T: (n_orbitals, n_orbitals) kinetic matrix
        """
        T = np.zeros((self.n_orbitals, self.n_orbitals))
        
        for i in range(self.n_orbitals):
            for j in range(self.n_orbitals):
                T[i, j] = self._kinetic_integral(i, j)
        
        return T
    
    def compute_nuclear_attraction(self) -> np.ndarray:
        """
        Compute nuclear-electron attraction V_ij = ⟨φ_i|Σ_A -Z_A/r_iA|φ_j⟩
        
        Returns:
            V: (n_orbitals, n_orbitals) nuclear attraction matrix
        """
        V = np.zeros((self.n_orbitals, self.n_orbitals))
        
        for i in range(self.n_orbitals):
            for j in range(self.n_orbitals):
                for atom in self.molecule.atoms:
                    V[i, j] += self._nuclear_integral(i, j, atom)
        
        return V
    
    def compute_core_hamiltonian(self) -> np.ndarray:
        """
        Compute core Hamiltonian H_core = T + V_ne
        
        Returns:
            H_core: (n_orbitals, n_orbitals) core Hamiltonian
        """
        T = self.compute_kinetic()
        V = self.compute_nuclear_attraction()
        return T + V
```

```python
# kanad/core/integrals/two_electron.py

class TwoElectronIntegrals:
    """Compute electron-electron repulsion integrals."""
    
    def __init__(self, molecule: 'Molecule', basis_set: str = 'sto-3g'):
        self.molecule = molecule
        self.basis = BasisSet(basis_set)
        self.n_orbitals = len(molecule.atoms) * self.basis.n_functions
        
    def compute_eri(self) -> np.ndarray:
        """
        Compute electron repulsion integrals (ij|kl) = ⟨ij|r₁₂⁻¹|kl⟩
        
        Returns:
            ERI: (n_orbitals, n_orbitals, n_orbitals, n_orbitals) tensor
        """
        ERI = np.zeros((self.n_orbitals,) * 4)
        
        # Use 8-fold symmetry: (ij|kl) = (ji|kl) = (ij|lk) = (ji|lk)
        #                                = (kl|ij) = (lk|ij) = (kl|ji) = (lk|ji)
        
        for i in range(self.n_orbitals):
            for j in range(i + 1):
                for k in range(self.n_orbitals):
                    for l in range(k + 1):
                        if self._quartet_index(i, j, k, l):
                            value = self._eri_integral(i, j, k, l)
                            self._fill_symmetric_eri(ERI, i, j, k, l, value)
        
        return ERI
    
    def compute_eri_compressed(self) -> Dict[Tuple[int, ...], float]:
        """
        Compute ERIs in compressed format (only unique values).
        Uses chemist's notation (ij|kl).
        """
        eri_dict = {}
        
        for i in range(self.n_orbitals):
            for j in range(i + 1):
                ij = self._compound_index(i, j)
                for k in range(self.n_orbitals):
                    for l in range(k + 1):
                        kl = self._compound_index(k, l)
                        if ij >= kl:
                            value = self._eri_integral(i, j, k, l)
                            if abs(value) > 1e-12:  # Sparsity threshold
                                eri_dict[(i, j, k, l)] = value
        
        return eri_dict
```

---

### **Phase 3: Representation Layer** 

```python
# kanad/core/representations/base_representation.py

from abc import ABC, abstractmethod
from typing import List, Dict, Any
import numpy as np

class BaseRepresentation(ABC):
    """Abstract base class for quantum representations of bonding."""
    
    def __init__(self, molecule: 'Molecule'):
        self.molecule = molecule
        self.n_qubits = None
        self.n_electrons = molecule.n_electrons
        
    @abstractmethod
    def build_hamiltonian(self) -> 'Hamiltonian':
        """Build the Hamiltonian in this representation."""
        pass
    
    @abstractmethod
    def get_ground_state(self) -> np.ndarray:
        """Get the expected ground state wavefunction."""
        pass
    
    @abstractmethod
    def compute_observables(self, state: np.ndarray) -> Dict[str, float]:
        """Compute physical observables from quantum state."""
        pass
    
    @abstractmethod
    def to_qubit_operator(self) -> 'QubitOperator':
        """Map representation to qubit operators."""
        pass
```

```python
# kanad/core/representations/lcao_representation.py

class LCAORepresentation(BaseRepresentation):
    """
    Linear Combination of Atomic Orbitals representation for covalent bonding.
    
    Theory:
        |ψ_bond⟩ = Σ_i c_i |φ_i⟩
        
        where |φ_i⟩ are atomic orbitals and c_i are mixing coefficients.
        
    For a simple 2-center bond:
        |ψ_bond⟩ = c_A |φ_A⟩ + c_B |φ_B⟩
        
    Hybridization types:
        - sp: linear (180°)
        - sp²: trigonal planar (120°)
        - sp³: tetrahedral (109.5°)
    """
    
    def __init__(self, molecule: 'Molecule', hybridization: str = 'sp3'):
        super().__init__(molecule)
        self.hybridization = hybridization
        self.hybrid_orbitals = self._construct_hybrid_orbitals()
        
    def _construct_hybrid_orbitals(self) -> List['HybridOrbital']:
        """Construct hybridized atomic orbitals based on geometry."""
        hybrids = []
        
        for atom in self.molecule.atoms:
            if atom.element == 'C':
                if self.hybridization == 'sp3':
                    hybrids.extend(self._sp3_hybrids(atom))
                elif self.hybridization == 'sp2':
                    hybrids.extend(self._sp2_hybrids(atom))
                elif self.hybridization == 'sp':
                    hybrids.extend(self._sp_hybrids(atom))
            # Handle other elements...
            
        return hybrids
    
    def _sp3_hybrids(self, atom: 'Atom') -> List['HybridOrbital']:
        """
        Construct sp³ hybrid orbitals: |sp³⟩ = ½(|s⟩ + |px⟩ + |py⟩ + |pz⟩)
        
        Returns 4 tetrahedral hybrid orbitals.
        """
        s_orbital = atom.get_orbital('2s')
        px_orbital = atom.get_orbital('2px')
        py_orbital = atom.get_orbital('2py')
        pz_orbital = atom.get_orbital('2pz')
        
        # Tetrahedral geometry coefficients
        hybrids = [
            HybridOrbital([
                (0.5, s_orbital),
                (0.5, px_orbital),
                (0.5, py_orbital),
                (0.5, pz_orbital)
            ], direction=[1, 1, 1]),
            # ... 3 more hybrids
        ]
        
        return hybrids
    
    def build_hamiltonian(self) -> 'CovalentHamiltonian':
        """
        Build Hamiltonian in LCAO basis:
        
        H = Σ_i ε_i c†_i c_i + Σ_ij t_ij (c†_i c_j + h.c.) + Σ_ijkl V_ijkl c†_i c†_j c_k c_l
        
        where:
            ε_i: orbital energies
            t_ij: hopping/hybridization integrals
            V_ijkl: two-electron integrals in hybrid basis
        """
        # Compute integrals in hybrid orbital basis
        S = self._overlap_matrix()
        H_core = self._core_hamiltonian()
        ERI = self._two_electron_integrals()
        
        return CovalentHamiltonian(
            orbital_energies=np.diag(H_core),
            hopping_matrix=H_core,
            eri_tensor=ERI,
            hybrid_orbitals=self.hybrid_orbitals
        )
    
    def get_bonding_orbitals(self, atom_i: int, atom_j: int) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute bonding and antibonding molecular orbitals.
        
        |ψ_bonding⟩ = N_+ (|φ_i⟩ + |φ_j⟩)
        |ψ_antibonding⟩ = N_- (|φ_i⟩ - |φ_j⟩)
        
        Returns:
            (bonding_MO, antibonding_MO)
        """
        hybrid_i = self.hybrid_orbitals[atom_i]
        hybrid_j = self.hybrid_orbitals[atom_j]
        
        S_ij = self._overlap(hybrid_i, hybrid_j)
        
        # Normalization factors
        N_plus = 1.0 / np.sqrt(2 * (1 + S_ij))
        N_minus = 1.0 / np.sqrt(2 * (1 - S_ij))
        
        bonding = N_plus * (hybrid_i + hybrid_j)
        antibonding = N_minus * (hybrid_i - hybrid_j)
        
        return bonding, antibonding
```

```python
# kanad/core/representations/bloch_representation.py

class BlochRepresentation(BaseRepresentation):
    """
    Bloch wave representation for metallic bonding.
    
    Theory:
        |ψ_k⟩ = (1/√N) Σ_j e^(ikR_j) |j⟩
        
        where k is the crystal momentum and R_j are lattice positions.
        
    Band structure:
        ε(k) = ε_0 - 2t cos(ka)
        
        where t is the hopping parameter and a is the lattice constant.
    """
    
    def __init__(self, molecule: 'Molecule', lattice_type: str = '1d'):
        super().__init__(molecule)
        self.lattice_type = lattice_type
        self.k_points = self._generate_k_mesh()
        self.bands = None
        
    def _generate_k_mesh(self) -> np.ndarray:
        """Generate k-point mesh in Brillouin zone."""
        if self.lattice_type == '1d':
            # 1D chain: k ∈ [-π/a, π/a]
            n_k = 50
            return np.linspace(-np.pi, np.pi, n_k)
        elif self.lattice_type == '2d':
            # 2D lattice
            n_k = 20
            kx = np.linspace(-np.pi, np.pi, n_k)
            ky = np.linspace(-np.pi, np.pi, n_k)
            return np.meshgrid(kx, ky)
        # Add 3D...
    
    def build_hamiltonian(self) -> 'MetallicHamiltonian':
        """
        Build tight-binding Hamiltonian in k-space:
        
        H(k) = Σ_k ε(k) c†_k c_k
        
        where ε(k) is the band dispersion.
        """
        # Tight-binding parameters
        hopping_t = self._compute_hopping_parameter()
        on_site_ε = self._compute_on_site_energy()
        
        # Build H(k) for each k-point
        hamiltonians_k = []
        for k in self.k_points:
            H_k = self._hamiltonian_at_k(k, hopping_t, on_site_ε)
            hamiltonians_k.append(H_k)
        
        return MetallicHamiltonian(
            k_points=self.k_points,
            hamiltonians_k=hamiltonians_k,
            hopping_parameter=hopping_t,
            lattice=self.lattice_type
        )
    
    def _hamiltonian_at_k(self, k: float, t: float, ε_0: float) -> np.ndarray:
        """
        Compute H(k) for a specific k-point.
        
        For 1D chain:
            H(k) = ε_0 - 2t cos(ka)
        """
        if self.lattice_type == '1d':
            a = self._lattice_constant()
            energy = ε_0 - 2 * t * np.cos(k * a)
            return np.array([[energy]])
        # Add 2D, 3D...
    
    def compute_band_structure(self) -> Dict[str, np.ndarray]:
        """
        Compute full band structure ε_n(k).
        
        Returns:
            bands: Dictionary with k-points and energy eigenvalues
        """
        bands = {'k_points': self.k_points, 'energies': []}
        
        for k in self.k_points:
            H_k = self._hamiltonian_at_k(k, self.hopping_parameter, self.on_site_energy)
            eigenvalues = np.linalg.eigvalsh(H_k)
            bands['energies'].append(eigenvalues)
        
        bands['energies'] = np.array(bands['energies'])
        self.bands = bands
        return bands
    
    def compute_fermi_surface(self) -> Dict[str, Any]:
        """
        Compute Fermi surface: locus of k where ε(k) = E_F.
        """
        if self.bands is None:
            self.compute_band_structure()
        
        E_F = self._fermi_energy()
        
        fermi_surface = {
            'k_fermi': [],
            'fermi_energy': E_F,
            'dos_at_fermi': self._dos_at_energy(E_F)
        }
        
        # Find k-points where bands cross E_F
        for i, energy in enumerate(self.bands['energies']):
            if np.any(np.abs(energy - E_F) < 0.01):  # Tolerance
                fermi_surface['k_fermi'].append(self.k_points[i])
        
        return fermi_surface
    
    def to_real_space(self) -> 'SecondQuantization':
        """
        Transform from k-space back to real space via Fourier transform.
        
        |j⟩ = (1/√N) Σ_k e^(-ikR_j) |k⟩
        """
        # Inverse Fourier transform
        real_space_ham = self._inverse_fourier_transform()
        return SecondQuantization(real_space_ham)
```

---

### **Phase 4: Governance Protocol Layer** 

This is THE INNOVATION - let me detail this carefully:

```python
# kanad/governance/protocols/base_protocol.py

from abc import ABC, abstractmethod
from typing import List, Dict, Any, Tuple
from enum import Enum

class BondingType(Enum):
    IONIC = "ionic"
    COVALENT = "covalent"
    METALLIC = "metallic"

class GovernanceRule:
    """A single governance rule for circuit construction."""
    
    def __init__(
        self,
        name: str,
        condition: callable,
        action: callable,
        priority: int = 0
    ):
        self.name = name
        self.condition = condition  # (state, context) -> bool
        self.action = action  # (state, context) -> modified_state
        self.priority = priority
        
    def applies(self, state: 'QuantumState', context: Dict) -> bool:
        return self.condition(state, context)
    
    def execute(self, state: 'QuantumState', context: Dict) -> 'QuantumState':
        return self.action(state, context)

class BaseGovernanceProtocol(ABC):
    """
    Abstract base class for governance protocols.
    
    Governance protocols DICTATE:
        1. Which operators are physically valid
        2. What circuit topology is appropriate
        3. How entanglement should be structured
        4. Which symmetries must be preserved
    """
    
    def __init__(self, bond_type: BondingType):
        self.bond_type = bond_type
        self.rules: List[GovernanceRule] = []
        self._initialize_rules()
        
    @abstractmethod
    def _initialize_rules(self):
        """Initialize bonding-specific governance rules."""
        pass
    
    @abstractmethod
    def validate_operator(self, operator: 'QuantumOperator') -> bool:
        """Check if operator is physically valid for this bonding type."""
        pass
    
    @abstractmethod
    def construct_ansatz(self, representation: 'BaseRepresentation') -> 'QuantumCircuit':
        """Construct physically appropriate ansatz circuit."""
        pass
    
    @abstractmethod
    def enforce_constraints(self, circuit: 'QuantumCircuit') -> 'QuantumCircuit':
        """Apply physical constraints to circuit."""
        pass
    
    def apply_governance(
        self,
        initial_circuit: 'QuantumCircuit',
        context: Dict
    ) -> 'QuantumCircuit':
        """
        Apply all governance rules in priority order.
        """
        circuit = initial_circuit
        state = circuit.get_state()
        
        # Sort rules by priority
        sorted_rules = sorted(self.rules, key=lambda r: r.priority, reverse=True)
        
        for rule in sorted_rules:
            if rule.applies(state, context):
                state = rule.execute(state, context)
                circuit = self._state_to_circuit(state)
        
        return circuit
```

```python
# kanad/governance/protocols/covalent_protocol.py

class CovalentGovernanceProtocol(BaseGovernanceProtocol):
    """
    Governance protocol for covalent bonding.
    
    PHYSICAL PRINCIPLES:
        1. Orbital hybridization must occur before bonding
        2. Bonding/antibonding pairs must be created
        3. Electron pairs must be entangled (Bell states)
        4. Symmetry: σ/π bond character must be preserved
        
    CIRCUIT REQUIREMENTS:
        1. Paired qubit entanglement (not chained)
        2. Variational orbital rotations
        3. Symmetry-adapted basis
    """
    
    def __init__(self):
        super().__init__(BondingType.COVALENT)
        
    def _initialize_rules(self):
        """Initialize covalent bonding rules."""
        
        # Rule 1: Hybridization before bonding
        self.rules.append(GovernanceRule(
            name="hybridization_first",
            condition=lambda state, ctx: not state.is_hybridized,
            action=self._apply_hybridization,
            priority=100
        ))
        
        # Rule 2: Create bonding/antibonding pairs
        self.rules.append(GovernanceRule(
            name="molecular_orbital_formation",
            condition=lambda state, ctx: state.is_hybridized and not state.has_mo_pairs,
            action=self._create_mo_pairs,
            priority=90
        ))
        
        # Rule 3: Entangle electron pairs
        self.rules.append(GovernanceRule(
            name="electron_pair_entanglement",
            condition=lambda state, ctx: state.has_mo_pairs and not state.is_paired,
            action=self._entangle_pairs,
            priority=80
        ))
        
        # Rule 4: Preserve spin symmetry
        self.rules.append(GovernanceRule(
            name="spin_symmetry",
            condition=lambda state, ctx: True,  # Always check
            action=self._enforce_spin_symmetry,
            priority=70
        ))
        
    def validate_operator(self, operator: 'QuantumOperator') -> bool:
        """
        Validate if operator is appropriate for covalent bonding.
        
        ALLOWED:
            - Orbital rotation operators (single-qubit rotations)
            - Pairing operators (2-qubit gates on bonding pairs)
            - Symmetric excitations
            
        FORBIDDEN:
            - Long-range entanglement (>2 qubits for single bond)
            - Asymmetric excitations
            - Bare transfer operators (must go through MO formation)
        """
        if isinstance(operator, TransferOperator):
            # Raw transfer not allowed - must use MO formation
            return False
        
        if isinstance(operator, EntanglingGate):
            # Check if it's a valid pairing gate
            if len(operator.qubits) != 2:
                return False
            if not self._are_bonding_pair(operator.qubits[0], operator.qubits[1]):
                return False
        
        return True
    
    def construct_ansatz(self, representation: LCAORepresentation) -> 'QuantumCircuit':
        """
        Construct covalent bonding ansatz.
        
        CIRCUIT STRUCTURE:
            1. Hybrid orbital preparation (single-qubit rotations)
            2. Bonding/antibonding MO formation (Givens rotations)
            3. Electron pairing (Bell state preparation)
            4. Variational correlation (small entanglement)
        """
        circuit = QuantumCircuit(representation.n_qubits)
        
        # Step 1: Prepare hybrid orbitals
        for i, hybrid in enumerate(representation.hybrid_orbitals):
            # Rotate atomic orbital basis to hybrid basis
            theta, phi = self._hybridization_angles(hybrid)
            circuit.ry(theta, i)
            circuit.rz(phi, i)
        
        # Step 2: Form molecular orbitals
        for bond in representation.bonds:
            i, j = bond.atom_indices
            # Givens rotation creates bonding/antibonding combination
            theta = self._bonding_angle(bond)
            circuit.givens(theta, i, j)
        
        # Step 3: Create electron pairs (Bell states)
        for pair in representation.electron_pairs:
            i, j = pair
            # |00⟩ + |11⟩ (bonding configuration)
            circuit.h(i)
            circuit.cx(i, j)
        
        # Step 4: Variational correlation
        for bond in representation.bonds:
            i, j = bond.atom_indices
            # Small amplitude excitations
            circuit.rxx(Parameter(f'θ_rxx_{i}{j}'), i, j)
            circuit.ryy(Parameter(f'θ_ryy_{i}{j}'), i, j)
            circuit.rzz(Parameter(f'θ_rzz_{i}{j}'), i, j)
        
        return circuit
    
    def _apply_hybridization(self, state: 'QuantumState', context: Dict) -> 'QuantumState':
        """Apply orbital hybridization transformation."""
        hybridization_type = context.get('hybridization', 'sp3')
        
        if hybridization_type == 'sp3':
            # |hybrid⟩ = ½(|s⟩ + |px⟩ + |py⟩ + |pz⟩)
            mixing_matrix = np.array([
                [0.5, 0.5, 0.5, 0.5],
                [0.5, -0.5, -0.5, 0.5],
                [0.5, -0.5, 0.5, -0.5],
                [0.5, 0.5, -0.5, -0.5]
            ])
        elif hybridization_type == 'sp2':
            # Trigonal planar
            mixing_matrix = self._sp2_mixing_matrix()
        # ... etc
        
        state.apply_transformation(mixing_matrix)
        state.is_hybridized = True
        return state
    
    def _create_mo_pairs(self, state: 'QuantumState', context: Dict) -> 'QuantumState':
        """Create bonding/antibonding molecular orbital pairs."""
        for bond in context['bonds']:
            i, j = bond.atom_indices
            overlap = bond.overlap
            
            # Normalization
            N_plus = 1.0 / np.sqrt(2 * (1 + overlap))
            N_minus = 1.0 / np.sqrt(2 * (1 - overlap))
            
            # Bonding: |ψ_+⟩ = N_+ (|i⟩ + |j⟩)
            # Antibonding: |ψ_-⟩ = N_- (|i⟩ - |j⟩)
            
            state.create_mo_pair(i, j, N_plus, N_minus)
        
        state.has_mo_pairs = True
        return state
    
    def _entangle_pairs(self, state: 'QuantumState', context: Dict) -> 'QuantumState':
        """Entangle electron pairs in bonding orbitals."""
        for pair in context['electron_pairs']:
            i, j = pair
            # Create Bell state: (|01⟩ + |10⟩)/√2 for opposite spins
            state.apply_bell_pair(i, j)
        
        state.is_paired = True
        return state
```

```python
# kanad/governance/protocols/metallic_protocol.py

class MetallicGovernanceProtocol(BaseGovernanceProtocol):
    """
    Governance protocol for metallic bonding.
    
    PHYSICAL PRINCIPLES:
        1. Electrons must be delocalized over ALL atoms
        2. Wavefunction must have Bloch wave character
        3. Multi-qubit entanglement (GHZ-like states)
        4. Translation symmetry must be preserved
        
    CIRCUIT REQUIREMENTS:
        1. Quantum Fourier Transform (real→k space)
        2. Highly entangled ansatz
        3. Periodic boundary conditions
    """
    
    def __init__(self):
        super().__init__(BondingType.METALLIC)
        
    def _initialize_rules(self):
        """Initialize metallic bonding rules."""
        
        # Rule 1: Transform to momentum space
        self.rules.append(GovernanceRule(
            name="momentum_space_transform",
            condition=lambda state, ctx: not state.in_k_space,
            action=self._apply_qft,
            priority=100
        ))
        
        # Rule 2: Create delocalized state
        self.rules.append(GovernanceRule(
            name="delocalization",
            condition=lambda state, ctx: state.in_k_space and state.is_localized,
            action=self._delocalize,
            priority=90
        ))
        
        # Rule 3: Build GHZ-like entanglement
        self.rules.append(GovernanceRule(
            name="collective_entanglement",
            condition=lambda state, ctx: not state.is_collectively_entangled,
            action=self._create_ghz_state,
            priority=80
        ))
        
    def construct_ansatz(self, representation: BlochRepresentation) -> 'QuantumCircuit':
        """
        Construct metallic bonding ansatz.
        
        CIRCUIT STRUCTURE:
            1. QFT (transform to k-space)
            2. k-dependent rotations (band energies)
            3. Multi-qubit entanglement (collective state)
            4. Inverse QFT (back to real space for measurements)
        """
        circuit = QuantumCircuit(representation.n_qubits)
        n = representation.n_qubits
        
        # Step 1: Quantum Fourier Transform
        circuit = self._apply_qft_circuit(circuit, n)
        
        # Step 2: Apply k-dependent phase rotations
        for k_idx in range(n):
            k = representation.k_points[k_idx]
            energy_k = representation.band_energy(k)
            circuit.rz(2 * energy_k, k_idx)
        
        # Step 3: Create collective entanglement (GHZ state)
        # |GHZ⟩ = (|000...0⟩ + |111...1⟩)/√2
        circuit.h(0)
        for i in range(1, n):
            circuit.cx(0, i)
        
        # Step 4: Variational layers with full connectivity
        for layer in range(context.get('n_layers', 3)):
            # All-to-all entanglement
            for i in range(n):
                for j in range(i + 1, n):
                    circuit.rzz(Parameter(f'θ_{layer}_{i}_{j}'), i, j)
            
            # Single-qubit rotations
            for i in range(n):
                circuit.ry(Parameter(f'φ_{layer}_{i}'), i)
        
        # Step 5: Inverse QFT for measurement
        circuit = self._apply_inverse_qft_circuit(circuit, n)
        
        return circuit
    
    def _apply_qft(self, state: 'QuantumState', context: Dict) -> 'QuantumState':
        """Apply Quantum Fourier Transform: |j⟩ → |k⟩"""
        n = state.n_qubits
        qft_matrix = self._qft_matrix(n)
        state.apply_transformation(qft_matrix)
        state.in_k_space = True
        return state
    
    def _delocalize(self, state: 'QuantumState', context: Dict) -> 'QuantumState':
        """Create delocalized Bloch state: |ψ_k⟩ = (1/√N) Σ_j e^(ikj) |j⟩"""
        n = state.n_qubits
        
        # Equal superposition over all sites
        delocalization_op = np.ones((n, n)) / np.sqrt(n)
        
        state.apply_transformation(delocalization_op)
        state.is_localized = False
        return state
    
    def _create_ghz_state(self, state: 'QuantumState', context: Dict) -> 'QuantumState':
        """Create GHZ-like maximally entangled state."""
        n = state.n_qubits
        
        # |GHZ⟩ = (|0⟩^⊗n + |1⟩^⊗n) / √2
        ghz_state = np.zeros(2**n)
        ghz_state[0] = 1/np.sqrt(2)  # |000...0⟩
        ghz_state[-1] = 1/np.sqrt(2)  # |111...1⟩
        
        state.set_statevector(ghz_state)
        state.is_collectively_entangled = True
        return state
```

# **Kanad Framework - Complete Implementation Plan (Phases 5-10)**

---

## **PHASE 5: Hamiltonian Builders** 

### **Bond-Specific Hamiltonian Construction**

```python
# kanad/core/hamiltonians/molecular_hamiltonian.py

from abc import ABC, abstractmethod
from typing import List, Dict, Tuple
import numpy as np
from scipy.linalg import eigh

class MolecularHamiltonian(ABC):
    """
    Base class for molecular Hamiltonians.
    
    General form:
        H = H_core + H_ee + H_nn
        
    where:
        H_core = T + V_ne (kinetic + nuclear-electron)
        H_ee = electron-electron repulsion
        H_nn = nuclear-nuclear repulsion (constant)
    """
    
    def __init__(
        self,
        n_orbitals: int,
        n_electrons: int,
        nuclear_repulsion: float = 0.0
    ):
        self.n_orbitals = n_orbitals
        self.n_electrons = n_electrons
        self.nuclear_repulsion = nuclear_repulsion
        
        # Integral matrices (to be filled by subclasses)
        self.h_core = None  # (n_orbitals, n_orbitals)
        self.eri = None     # (n_orbitals, n_orbitals, n_orbitals, n_orbitals)
        
    @abstractmethod
    def compute_energy(self, density_matrix: np.ndarray) -> float:
        """Compute total energy given density matrix."""
        pass
    
    @abstractmethod
    def to_second_quantized(self) -> 'FermionicOperator':
        """Convert to second-quantized form."""
        pass
    
    def get_fci_energy(self) -> Tuple[float, np.ndarray]:
        """
        Compute exact Full Configuration Interaction energy.
        Used for validation and benchmarking.
        """
        from pyscf import fci
        # FCI solver implementation
        pass
    
    def get_hf_energy(self) -> Tuple[float, np.ndarray]:
        """
        Compute Hartree-Fock energy and orbitals.
        Used as initial guess for VQE.
        """
        # Self-consistent field iteration
        pass
```

```python
# kanad/core/hamiltonians/ionic_hamiltonian.py

class IonicHamiltonian(MolecularHamiltonian):
    """
    Hamiltonian for ionic bonding.
    
    KEY FEATURE: Dominated by electron transfer terms
    
    H_ionic = Σ_i ε_i n_i + Σ_ij t_ij a†_i a_j + U Σ_i n_i↑ n_i↓
    
    where:
        ε_i: site energies (electronegativity difference drives transfer)
        t_ij: hopping amplitude (small, nearest-neighbor only)
        U: on-site Coulomb repulsion
    """
    
    def __init__(
        self,
        donor_atom: 'Atom',
        acceptor_atom: 'Atom',
        integrals: 'OneElectronIntegrals'
    ):
        super().__init__(
            n_orbitals=2,  # Simplified 2-orbital model
            n_electrons=2,
            nuclear_repulsion=self._compute_nuclear_repulsion(donor_atom, acceptor_atom)
        )
        
        self.donor = donor_atom
        self.acceptor = acceptor_atom
        self.integrals = integrals
        
        # Build Hamiltonian matrices
        self._build_hamiltonian()
        
    def _build_hamiltonian(self):
        """
        Build ionic Hamiltonian emphasizing charge transfer.
        """
        # Site energies (electronegativity difference)
        en_donor = self.donor.properties.electronegativity
        en_acceptor = self.acceptor.properties.electronegativity
        
        # Donor has HIGHER energy (easier to lose electron)
        # Acceptor has LOWER energy (wants to gain electron)
        ε_donor = -en_donor  # More negative = lower energy
        ε_acceptor = -en_acceptor
        
        # Core Hamiltonian
        self.h_core = np.array([
            [ε_donor, self._compute_transfer_integral()],
            [self._compute_transfer_integral(), ε_acceptor]
        ])
        
        # Two-electron integrals
        self.eri = self._compute_eri_ionic()
        
    def _compute_transfer_integral(self) -> float:
        """
        Compute electron transfer integral t_DA.
        
        t_DA = ⟨φ_D|H|φ_A⟩
        
        This is SMALL for ionic bonds (weak overlap).
        """
        overlap = self.integrals.compute_overlap_between(
            self.donor.orbital,
            self.acceptor.orbital
        )
        
        # Transfer integral scales with overlap
        # For ionic: small overlap → small t
        distance = np.linalg.norm(self.donor.position - self.acceptor.position)
        t = overlap * np.exp(-distance / 2.0)  # Exponential decay
        
        return t
    
    def _compute_eri_ionic(self) -> np.ndarray:
        """
        Compute two-electron repulsion integrals.
        
        For ionic bonding, the dominant term is:
            U = ⟨ii|ii⟩ (on-site repulsion)
        """
        eri = np.zeros((2, 2, 2, 2))
        
        # On-site repulsion (Hubbard U)
        U_donor = self._hubbard_u(self.donor)
        U_acceptor = self._hubbard_u(self.acceptor)
        
        eri[0, 0, 0, 0] = U_donor
        eri[1, 1, 1, 1] = U_acceptor
        
        # Inter-site repulsion (much smaller)
        V_inter = self._inter_site_coulomb()
        eri[0, 0, 1, 1] = V_inter
        eri[1, 1, 0, 0] = V_inter
        
        return eri
    
    def to_second_quantized(self) -> 'FermionicOperator':
        """
        Convert to fermionic operators.
        
        H = ε_D a†_D a_D + ε_A a†_A a_A + t(a†_D a_A + a†_A a_D)
          + U_D n_D↑ n_D↓ + U_A n_A↑ n_A↓
        """
        from qiskit_nature.second_q.operators import FermionicOp
        
        ops = []
        
        # Site energies
        ops.append((f'+_0 -_0', self.h_core[0, 0]))  # Donor
        ops.append((f'+_1 -_1', self.h_core[1, 1]))  # Acceptor
        
        # Transfer terms (hopping)
        t = self.h_core[0, 1]
        ops.append((f'+_0 -_1', t))  # A → D transfer
        ops.append((f'+_1 -_0', t))  # D → A transfer
        
        # On-site repulsion (spin sectors)
        U_D = self.eri[0, 0, 0, 0]
        U_A = self.eri[1, 1, 1, 1]
        ops.append((f'+_0 -_0 +_2 -_2', U_D))  # Donor up-down
        ops.append((f'+_1 -_1 +_3 -_3', U_A))  # Acceptor up-down
        
        return FermionicOp(ops, num_spin_orbitals=4)
    
    def compute_charge_transfer(self, density_matrix: np.ndarray) -> float:
        """
        Compute amount of charge transferred from donor to acceptor.
        
        Δq = ρ_AA - ρ_DD
        
        For full transfer (Na+ Cl-): Δq ≈ 1.0
        """
        charge_donor = density_matrix[0, 0]
        charge_acceptor = density_matrix[1, 1]
        
        return charge_acceptor - charge_donor
```

```python
# kanad/core/hamiltonians/covalent_hamiltonian.py

class CovalentHamiltonian(MolecularHamiltonian):
    """
    Hamiltonian for covalent bonding.
    
    KEY FEATURE: Orbital hybridization and MO formation
    
    H_covalent = Σ_μν h_μν c†_μ c_ν + ½ Σ_μνλσ (μν|λσ) c†_μ c†_ν c_σ c_λ
    
    where μ,ν run over HYBRID ORBITALS, not atomic orbitals.
    
    This gives bonding/antibonding splitting:
        E_bonding = ⟨φ_+ | H | φ_+⟩
        E_antibonding = ⟨φ_- | H | φ_-⟩
    """
    
    def __init__(
        self,
        atoms: List['Atom'],
        bonds: List['CovalentBond'],
        hybridization: str = 'sp3'
    ): 
        self.atoms = atoms
        self.bonds = bonds
        self.hybridization = hybridization
        
        # Generate hybrid orbitals
        self.hybrid_orbitals = self._generate_hybrids()
        
        super().__init__(
            n_orbitals=len(self.hybrid_orbitals),
            n_electrons=sum(atom.n_valence for atom in atoms),
            nuclear_repulsion=self._compute_nuclear_repulsion(atoms)
        )
        
        self._build_hamiltonian()
        
    def _generate_hybrids(self) -> List['HybridOrbital']:
        """Generate hybridized atomic orbitals."""
        hybrids = []
        
        for atom in self.atoms:
            if self.hybridization == 'sp3':
                # Carbon: 1 s + 3 p → 4 sp³ hybrids
                s_coeff = 0.5
                p_coeff = 0.5
                
                hybrids.extend([
                    HybridOrbital(atom, [
                        (s_coeff, '2s'),
                        (p_coeff, '2px'),
                        (p_coeff, '2py'),
                        (p_coeff, '2pz')
                    ], direction=[1, 1, 1]),
                    # ... 3 more tetrahedrally oriented
                ])
                
            elif self.hybridization == 'sp2':
                # Trigonal planar: 1 s + 2 p → 3 sp² + 1 pure p
                s_coeff = 1/np.sqrt(3)
                p_coeff = np.sqrt(2/3)
                
                hybrids.extend([
                    HybridOrbital(atom, [
                        (s_coeff, '2s'),
                        (p_coeff, '2px'),
                    ], direction=[1, 0, 0]),
                    # ... 2 more at 120°
                    HybridOrbital(atom, [(1.0, '2pz')], direction=[0, 0, 1])  # π orbital
                ])
                
        return hybrids
    
    def _build_hamiltonian(self):
        """
        Build Hamiltonian in hybrid orbital basis.
        """
        n = len(self.hybrid_orbitals)
        
        # Transform atomic integrals to hybrid basis
        U = self._hybridization_matrix()  # AO → HO transformation
        
        # h_hybrid = U† h_atomic U
        h_atomic = self._atomic_core_hamiltonian()
        self.h_core = U.T @ h_atomic @ U
        
        # ERI transformation (4-index)
        eri_atomic = self._atomic_eri()
        self.eri = self._transform_eri(eri_atomic, U)
        
    def _hybridization_matrix(self) -> np.ndarray:
        """
        Construct transformation matrix from atomic to hybrid orbitals.
        
        For sp³:
            |sp³_1⟩ = ½(|s⟩ + |px⟩ + |py⟩ + |pz⟩)
            |sp³_2⟩ = ½(|s⟩ + |px⟩ - |py⟩ - |pz⟩)
            |sp³_3⟩ = ½(|s⟩ - |px⟩ + |py⟩ - |pz⟩)
            |sp³_4⟩ = ½(|s⟩ - |px⟩ - |py⟩ + |pz⟩)
        """
        if self.hybridization == 'sp3':
            return 0.5 * np.array([
                [ 1,  1,  1,  1],
                [ 1,  1, -1, -1],
                [ 1, -1,  1, -1],
                [ 1, -1, -1,  1]
            ])
        # Add sp², sp, etc.
    
    def compute_molecular_orbitals(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Diagonalize to get bonding/antibonding molecular orbitals.
        
        Returns:
            energies: MO energies
            coefficients: MO coefficients in hybrid basis
        """
        # Solve generalized eigenvalue problem: H C = S C E
        S = self._overlap_matrix()
        energies, coefficients = eigh(self.h_core, S)
        
        return energies, coefficients
    
    def get_bonding_antibonding_split(self, bond_idx: int) -> Dict[str, float]:
        """
        Compute bonding/antibonding energy splitting for a specific bond.
        
        Δε = E_antibonding - E_bonding
        
        Larger splitting → stronger bond
        """
        bond = self.bonds[bond_idx]
        i, j = bond.orbital_indices
        
        energies, coeffs = self.compute_molecular_orbitals()
        
        # Identify bonding MO (lower energy, same sign)
        # and antibonding MO (higher energy, opposite sign)
        bonding_idx = np.argmin(energies)
        antibonding_idx = np.argmax(energies)
        
        return {
            'bonding_energy': energies[bonding_idx],
            'antibonding_energy': energies[antibonding_idx],
            'splitting': energies[antibonding_idx] - energies[bonding_idx],
            'bonding_character': coeffs[:, bonding_idx],
            'antibonding_character': coeffs[:, antibonding_idx]
        }
    
    def to_second_quantized(self) -> 'FermionicOperator':
        """
        Express in second-quantized form with HYBRID orbitals.
        
        H = Σ_μν h_μν c†_μ c_ν + ½ Σ_μνλσ (μν|λσ) c†_μ c†_ν c_σ c_λ
        """
        from qiskit_nature.second_q.operators import FermionicOp
        
        ops = []
        n = self.n_orbitals
        
        # One-body terms (in hybrid basis)
        for i in range(n):
            for j in range(n):
                if abs(self.h_core[i, j]) > 1e-10:
                    ops.append((f'+_{i} -_{j}', self.h_core[i, j]))
        
        # Two-body terms
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    for l in range(n):
                        if abs(self.eri[i, j, k, l]) > 1e-10:
                            ops.append((
                                f'+_{i} +_{j} -_{l} -_{k}',
                                0.5 * self.eri[i, j, k, l]
                            ))
        
        return FermionicOp(ops, num_spin_orbitals=2*n)
```

```python
# kanad/core/hamiltonians/metallic_hamiltonian.py

class MetallicHamiltonian(MolecularHamiltonian):
    """
    Hamiltonian for metallic bonding.
    
    KEY FEATURE: k-space representation with band structure
    
    H_metallic = Σ_k ε(k) c†_k c_k
    
    where:
        ε(k) is the band dispersion
        k is crystal momentum
        
    Real space equivalent (tight-binding):
        H = -t Σ_⟨ij⟩ (c†_i c_j + h.c.) + ε_0 Σ_i c†_i c_i
    """
    
    def __init__(
        self,
        lattice_type: str,
        n_sites: int,
        hopping_parameter: float,
        on_site_energy: float = 0.0
    ):
        self.lattice_type = lattice_type
        self.n_sites = n_sites
        self.hopping_t = hopping_parameter
        self.on_site_ε = on_site_energy
        
        super().__init__(
            n_orbitals=n_sites,
            n_electrons=n_sites,  # Half-filling for simplicity
            nuclear_repulsion=0.0  # Ignore for now
        )
        
        # Generate k-space mesh
        self.k_points = self._generate_k_mesh()
        self.bands = None
        
        self._build_hamiltonian()
        
    def _build_hamiltonian(self):
        """
        Build tight-binding Hamiltonian.
        
        Real space:
            H_ij = -t if ⟨ij⟩ are nearest neighbors
            H_ii = ε_0
        """
        n = self.n_sites
        self.h_core = np.zeros((n, n))
        
        # On-site energies
        for i in range(n):
            self.h_core[i, i] = self.on_site_ε
        
        # Hopping terms (nearest-neighbor)
        if self.lattice_type == '1d_chain':
            for i in range(n - 1):
                self.h_core[i, i+1] = -self.hopping_t
                self.h_core[i+1, i] = -self.hopping_t
            # Periodic boundary conditions
            self.h_core[0, n-1] = -self.hopping_t
            self.h_core[n-1, 0] = -self.hopping_t
            
        elif self.lattice_type == '2d_square':
            # Square lattice connectivity
            for i in range(n):
                neighbors = self._get_2d_neighbors(i)
                for j in neighbors:
                    self.h_core[i, j] = -self.hopping_t
        
        # Two-electron terms (ignore for now, or use Hubbard U)
        self.eri = np.zeros((n, n, n, n))
    
    def compute_band_structure(self) -> Dict[str, np.ndarray]:
        """
        Compute ε(k) for all k-points.
        
        For 1D chain:
            ε(k) = ε_0 - 2t cos(ka)
        """
        energies = []
        
        for k in self.k_points:
            if self.lattice_type == '1d_chain':
                # Analytical dispersion
                a = 1.0  # Lattice constant
                ε_k = self.on_site_ε - 2 * self.hopping_t * np.cos(k * a)
                energies.append(ε_k)
            else:
                # Numerical diagonalization
                H_k = self._hamiltonian_at_k(k)
                ε_k = np.linalg.eigvalsh(H_k)
                energies.append(ε_k)
        
        self.bands = {
            'k_points': self.k_points,
            'energies': np.array(energies)
        }
        
        return self.bands
    
    def _hamiltonian_at_k(self, k: Union[float, np.ndarray]) -> np.ndarray:
        """
        Construct H(k) by Fourier transforming H(r).
        
        H(k) = Σ_R H(R) e^(ik·R)
        """
        n = self.n_sites
        H_k = np.zeros((n, n), dtype=complex)
        
        for i in range(n):
            for j in range(n):
                R_ij = self._lattice_vector(i, j)
                phase = np.exp(1j * np.dot(k, R_ij))
                H_k[i, j] = self.h_core[i, j] * phase
        
        return H_k
    
    def compute_fermi_energy(self, filling: float = 0.5) -> float:
        """
        Compute Fermi energy for given electron filling.
        
        E_F is the energy where integrated DOS equals n_electrons.
        """
        if self.bands is None:
            self.compute_band_structure()
        
        # Sort all eigenvalues
        all_energies = np.sort(self.bands['energies'].flatten())
        
        # Fermi energy at filling fraction
        idx_fermi = int(filling * len(all_energies))
        return all_energies[idx_fermi]
    
    def compute_dos(self, energy_range: Tuple[float, float], n_points: int = 1000) -> Dict:
        """
        Compute density of states ρ(E).
        
        ρ(E) = (1/N) Σ_k δ(E - ε_k)
        """
        if self.bands is None:
            self.compute_band_structure()
        
        energies = np.linspace(energy_range[0], energy_range[1], n_points)
        dos = np.zeros(n_points)
        
        # Gaussian broadening
        σ = 0.1  # Broadening parameter
        
        for ε_k in self.bands['energies'].flatten():
            dos += np.exp(-(energies - ε_k)**2 / (2*σ**2)) / (σ * np.sqrt(2*np.pi))
        
        dos /= len(self.bands['energies'])
        
        return {'energies': energies, 'dos': dos}
    
    def to_second_quantized(self) -> 'FermionicOperator':
        """
        Express in k-space second quantization.
        
        H = Σ_k ε(k) c†_k c_k
        """
        from qiskit_nature.second_q.operators import FermionicOp
        
        if self.bands is None:
            self.compute_band_structure()
        
        ops = []
        for k_idx, ε_k in enumerate(self.bands['energies']):
            ops.append((f'+_{k_idx} -_{k_idx}', ε_k))
        
        return FermionicOp(ops, num_spin_orbitals=len(self.k_points))
```

---

## **PHASE 6: Mapper Implementation** 

### **Custom Qubit Mappings for Each Bonding Type**

```python
# kanad/core/mappers/base_mapper.py

from abc import ABC, abstractmethod
from qiskit_nature.second_q.operators import FermionicOp
from qiskit.quantum_info import SparsePauliOp

class BaseMapper(ABC):
    """
    Abstract base class for fermionic-to-qubit mappers.
    
    Each bonding type needs a DIFFERENT mapping strategy:
        - Ionic: Jordan-Wigner (local, sequential)
        - Covalent: Paired mapping (orbital-centric)
        - Metallic: Momentum-space mapping (collective)
    """
    
    @abstractmethod
    def map(self, fermionic_op: FermionicOp) -> SparsePauliOp:
        """Map fermionic operator to qubit operator."""
        pass
    
    @abstractmethod
    def inverse_map(self, qubit_op: SparsePauliOp) -> FermionicOp:
        """Map qubit operator back to fermionic."""
        pass
    
    @abstractmethod
    def n_qubits(self, n_spin_orbitals: int) -> int:
        """Compute number of qubits needed."""
        pass
```

```python
# kanad/core/mappers/hybrid_orbital_mapper.py

class HybridOrbitalMapper(BaseMapper):
    """
    Custom mapper for covalent bonding with hybrid orbitals.
    
    KEY INNOVATION: Maps bonding/antibonding PAIRS to qubits
    
    Instead of:
        |orbital_i⟩ → |qubit_i⟩ (Jordan-Wigner)
        
    We use:
        |bonding⟩, |antibonding⟩ → |qubit_pair⟩
        
    This naturally encodes the molecular orbital structure.
    
    Encoding:
        |00⟩ → both empty
        |01⟩ → antibonding occupied
        |10⟩ → bonding occupied  (ground state)
        |11⟩ → both occupied (excited state)
    """
    
    def __init__(self, molecular_orbitals: List[Tuple[int, int]]):
        """
        Args:
            molecular_orbitals: List of (bonding, antibonding) orbital pairs
        """
        self.mo_pairs = molecular_orbitals
        self.n_pairs = len(molecular_orbitals)
        
    def n_qubits(self, n_spin_orbitals: int) -> int:
        """
        Each MO pair → 2 qubits
        
        For N hybrid orbitals → N/2 bonds → N qubits
        """
        return n_spin_orbitals  # Same as JW, but different interpretation
    
    def map(self, fermionic_op: FermionicOp) -> SparsePauliOp:
        """
        Map fermionic operators in MO basis to Pauli operators.
        
        Key transformations:
            c†_bonding c_bonding → ½(I - Z_i)  (bonding occupancy)
            c†_antibonding c_antibonding → ½(I - Z_j)  (antibonding occupancy)
            c†_bonding c_antibonding → X_i X_j + Y_i Y_j  (excitation)
        """
        pauli_sum = SparsePauliOp([])
        
        for term, coeff in fermionic_op.items():
            pauli_term = self._map_term(term)
            pauli_sum += coeff * pauli_term
        
        return pauli_sum
    
    def _map_term(self, fermionic_term: str) -> SparsePauliOp:
        """
        Map individual fermionic term to Pauli string.
        
        Examples:
            '+_0 -_0' (bonding occupancy) → (I - Z_0) / 2
            '+_1 -_1' (antibonding) → (I - Z_1) / 2
            '+_0 -_1' (bonding→antibonding excitation) → (X_0 X_1 + Y_0 Y_1) / 2
        """
        # Parse fermionic term
        creators, annihilators = self._parse_term(fermionic_term)
        
        if creators == annihilators:
            # Number operator: n_i = c†_i c_i
            i = creators[0]
            return 0.5 * (SparsePauliOp(['I']) - SparsePauliOp(['Z'], coeffs=[1.0], paulis_at=[i]))
        
        elif len(creators) == 1 and len(annihilators) == 1:
            # Excitation operator: c†_i c_j
            i, j = creators[0], annihilators[0]
            
            if self._are_mo_pair(i, j):
                # Bonding-antibonding excitation
                return 0.5 * (
                    self._pauli_string('X', [i, j]) +
                    self._pauli_string('Y', [i, j])
                )
            else:
                # Inter-orbital transfer (use JW-like)
                return self._jordan_wigner_excitation(i, j)
        
        # Handle more complex terms...
    
    def _are_mo_pair(self, i: int, j: int) -> bool:
        """Check if orbitals i and j form a bonding/antibonding pair."""
        return (i, j) in self.mo_pairs or (j, i) in self.mo_pairs
    
    def create_excitation_operator(
        self,
        bonding_idx: int,
        antibonding_idx: int
    ) -> SparsePauliOp:
        """
        Create excitation from bonding to antibonding MO.
        
        This is the KEY covalent excitation operator:
            |bonding occupied⟩ → |antibonding occupied⟩
        """
        i, j = bonding_idx, antibonding_idx
        
        # RXX + RYY gate in Pauli form
        excitation = (
            SparsePauliOp.from_list([('XX', 0.25)], num_qubits=self.n_qubits) +
            SparsePauliOp.from_list([('YY', 0.25)], num_qubits=self.n_qubits)
        )
        
        return excitation
```

```python
# kanad/core/mappers/momentum_space_mapper.py

class MomentumSpaceMapper(BaseMapper):
    """
    Custom mapper for metallic bonding in k-space.
    
    KEY INNOVATION: Encodes momentum eigenstates
    
    Instead of position basis |r⟩, use momentum basis |k⟩:
        |k⟩ = (1/√N) Σ_r e^(ikr) |r⟩
        
    This is achieved via Quantum Fourier Transform on qubits.
    
    Encoding:
        Computational basis |j⟩ represents k_j in Brillouin zone
        
    Advantages:
        - Band structure directly accessible
        - Translation symmetry manifest
        - Collective excitations natural
    """
    
    def __init__(self, n_k_points: int, lattice_constant: float = 1.0):
        self.n_k = n_k_points
        self.a = lattice_constant
        self.k_points = self._generate_k_mesh()
        
    def _generate_k_mesh(self) -> np.ndarray:
        """Generate k-points in Brillouin zone [-π/a, π/a]."""
        return np.linspace(-np.pi/self.a, np.pi/self.a, self.n_k, endpoint=False)
    
    def n_qubits(self, n_spin_orbitals: int) -> int:
        """Number of qubits = log2(n_k_points)."""
        return int(np.ceil(np.log2(self.n_k)))
    
    def map(self, fermionic_op: FermionicOp) -> SparsePauliOp:
        """
        Map fermionic operators in k-space to Pauli operators.
        
        c†_k c_k → ½(I - Z_k)  (occupation at k)
        
        But k is encoded in binary, so we need to map k-index to qubits.
        """
        pauli_sum = SparsePauliOp([])
        
        for term, coeff in fermionic_op.items():
            k_idx = self._extract_k_index(term)
            pauli_term = self._map_k_occupation(k_idx)
            pauli_sum += coeff * pauli_term
        
        return pauli_sum
    
    def _map_k_occupation(self, k_idx: int) -> SparsePauliOp:
        """
        Map occupation at k-point to Pauli operator.
        
        Since k is encoded in binary across multiple qubits,
        we need to create a projector onto |k⟩ state.
        """
        n_q = self.n_qubits(self.n_k)
        binary_k = format(k_idx, f'0{n_q}b')
        
        # Create projector |k⟩⟨k| = (1/2)^n ∏_i (I + (-1)^{b_i} Z_i)
        projector = SparsePauliOp(['I' * n_q], coeffs=[1.0])
        
        for i, bit in enumerate(binary_k):
            if bit == '0':
                projector = 0.5 * (projector + SparsePauliOp(['Z'], coeffs=[1.0], paulis_at=[i]))
            else:
                projector = 0.5 * (projector - SparsePauliOp(['Z'], coeffs=[1.0], paulis_at=[i]))
        
        return projector
    
    def create_collective_excitation(self) -> SparsePauliOp:
        """
        Create collective excitation across all k-states.
        
        This is analogous to creating a plasmon excitation.
        """
        n_q = self.n_qubits(self.n_k)
        
        # Superposition of all k-states
        collective = SparsePauliOp(['X' * n_q], coeffs=[1.0 / np.sqrt(self.n_k)])
        
        return collective
    
    def apply_qft_encoding(self, circuit: 'QuantumCircuit') -> 'QuantumCircuit':
        """
        Apply Quantum Fourier Transform to convert real→k space.
        
        This should be applied BEFORE measurement to get k-space properties.
        """
        n_q = circuit.num_qubits
        
        # Standard QFT circuit
        for j in range(n_q):
            circuit.h(j)
            for k in range(j + 1, n_q):
                circuit.cp(2 * np.pi / (2 ** (k - j + 1)), k, j)
        
        # Reverse qubit order (bit-reversal)
        for i in range(n_q // 2):
            circuit.swap(i, n_q - i - 1)
        
        return circuit
```

---

## **PHASE 7: Energy Estimation & Analysis** 

```python
# kanad/analysis/energy_estimation.py

from qiskit.algorithms import VQE, NumPyMinimumEigensolver
from qiskit.algorithms.optimizers import SLSQP, COBYLA, SPSA
from typing import Dict, Any, Optional

class GovernanceVQE:
    """
    VQE with governance-aware circuit construction.
    
    INNOVATION: Ansatz is built by governance protocol, not user
    """
    
    def __init__(
        self,
        hamiltonian: 'MolecularHamiltonian',
        governance_protocol: 'BaseGovernanceProtocol',
        backend: Optional['Backend'] = None
    ):
        self.hamiltonian = hamiltonian
        self.governance = governance_protocol
        self.backend = backend or Aer.get_backend('statevector_simulator')
        
        # Build governed ansatz
        self.ansatz = self.governance.construct_ansatz(
            self.hamiltonian.representation
        )
        
    def compute_energy(
        self,
        optimizer: str = 'SLSQP',
        initial_point: Optional[np.ndarray] = None
    ) -> Dict[str, Any]:
        """
        Run VQE with governance constraints.
        """
        # Map Hamiltonian to qubits
        qubit_ham = self._map_hamiltonian()
        
        # Select optimizer
        if optimizer == 'SLSQP':
            opt = SLSQP(maxiter=1000)
        elif optimizer == 'COBYLA':
            opt = COBYLA(maxiter=1000)
        elif optimizer == 'SPSA':
            opt = SPSA(maxiter=1000)
        
        # Create VQE instance
        vqe = VQE(
            ansatz=self.ansatz,
            optimizer=opt,
            quantum_instance=self.backend,
            initial_point=initial_point
        )
        
        # Run optimization
        result = vqe.compute_minimum_eigenvalue(qubit_ham)
        
        # Extract and analyze results
        return {
            'energy': result.eigenvalue.real,
            'optimal_parameters': result.optimal_parameters,
            'optimizer_evals': result.optimizer_evals,
            'optimal_circuit': self.ansatz.bind_parameters(result.optimal_parameters),
            'governance_metrics': self._analyze_governance(result)
        }
    
    def _map_hamiltonian(self) -> SparsePauliOp:
        """Map Hamiltonian using governance-appropriate mapper."""
        bond_type = self.governance.bond_type
        
        if bond_type == BondingType.IONIC:
            mapper = JordanWignerMapper()
        elif bond_type == BondingType.COVALENT:
            mapper = HybridOrbitalMapper(self.hamiltonian.mo_pairs)
        elif bond_type == BondingType.METALLIC:
            mapper = MomentumSpaceMapper(self.hamiltonian.n_k_points)
        
        fermionic_op = self.hamiltonian.to_second_quantized()
        return mapper.map(fermionic_op)
    
    def _analyze_governance(self, vqe_result) -> Dict[str, float]:
        """
        Analyze how well the result obeys governance rules.
        """
        optimal_circuit = vqe_result.optimal_circuit
        
        metrics = {
            'constraint_violations': self.governance.count_violations(optimal_circuit),
            'entanglement_structure': self._analyze_entanglement(optimal_circuit),
            'symmetry_preservation': self._check_symmetries(optimal_circuit)
        }
        
        return metrics
```

```python
# kanad/analysis/bonding_analysis.py

class BondingAnalyzer:
    """
    Analyze bonding characteristics from quantum states.
    """
    
    def __init__(self, bond: 'Bond', hamiltonian: 'MolecularHamiltonian'):
        self.bond = bond
        self.hamiltonian = hamiltonian
        
    def compute_bond_order(self, density_matrix: np.ndarray) -> float:
        """
        Compute bond order from density matrix.
        
        Bond order = (n_bonding - n_antibonding) / 2
        
        Where n is the occupation number.
        """
        if isinstance(self.hamiltonian, CovalentHamiltonian):
            bonding_idx, antibonding_idx = self.bond.mo_indices
            
            n_bonding = density_matrix[bonding_idx, bonding_idx]
            n_antibonding = density_matrix[antibonding_idx, antibonding_idx]
            
            return (n_bonding - n_antibonding) / 2.0
        
        # For other bond types...
    
    def compute_bond_length(self) -> float:
        """Compute equilibrium bond length."""
        return np.linalg.norm(
            self.bond.atom_i.position - self.bond.atom_j.position
        )
    
    def compute_bond_energy(self, total_energy: float) -> float:
        """
        Extract bond dissociation energy.
        
        D_e = E(separated) - E(bonded)
        """
        # Run calculation with atoms separated
        separated_energy = self._compute_separated_energy()
        
        return separated_energy - total_energy
    
    def analyze_bonding_character(self, wavefunction: np.ndarray) -> Dict[str, float]:
        """
        Analyze ionic vs covalent character.
        
        Returns:
            {'ionic': 0.3, 'covalent': 0.7}
        """
        if isinstance(self.bond, CovalentBond):
            # Compute wavefunction overlap
            overlap = self._compute_overlap(wavefunction)
            
            # Higher overlap → more covalent
            # Lower overlap → more ionic
            
            ionic_character = 1.0 - overlap
            covalent_character = overlap
            
            return {
                'ionic': ionic_character,
                'covalent': covalent_character,
                'polarity': self._compute_polarity(wavefunction)
            }
    
    def compute_dipole_moment(self, wavefunction: np.ndarray) -> np.ndarray:
        """
        Compute electric dipole moment μ = Σ q_i r_i
        """
        charges = self._compute_partial_charges(wavefunction)
        positions = [atom.position for atom in self.bond.atoms]
        
        dipole = np.sum([q * r for q, r in zip(charges, positions)], axis=0)
        return dipole
    
    def compute_mulliken_charges(self, density_matrix: np.ndarray) -> Dict[int, float]:
        """
        Compute Mulliken population analysis.
        
        q_A = Z_A - Σ_μ∈A Σ_ν (P S)_μν
        """
        S = self.hamiltonian._overlap_matrix()
        PS = density_matrix @ S
        
        charges = {}
        for atom_idx, atom in enumerate(self.bond.atoms):
            # Sum over orbitals on this atom
            orbital_indices = self._get_atom_orbitals(atom_idx)
            population = sum(PS[i, i] for i in orbital_indices)
            charges[atom_idx] = atom.atomic_number - population
        
        return charges
```

---

## **PHASE 8: Bond Factory & Research Interface** 

```python
# kanad/bonds/bond_factory.py

from typing import Union, List, Dict, Any
from enum import Enum

class BondType(Enum):
    IONIC = "ionic"
    COVALENT = "covalent"  
    METALLIC = "metallic"
    AUTO = "auto"  # Determine from electronegativity

class BondFactory:
    """
    Factory for creating bonds with automatic governance.
    
    USER-FACING API - This is what researchers interact with!
    """
    
    @staticmethod
    def create_bond(
        atom_1: Union[str, 'Atom'],
        atom_2: Union[str, 'Atom'],
        bond_type: BondType = BondType.AUTO,
        **kwargs
    ) -> 'Bond':
        """
        Create a bond between two atoms.
        
        Args:
            atom_1: First atom (symbol or Atom object)
            atom_2: Second atom
            bond_type: Type of bond (or AUTO to determine)
            **kwargs: Additional parameters
            
        Returns:
            Bond object with appropriate governance
            
        Example:
            >>> bond = BondFactory.create_bond('Na', 'Cl')
            >>> print(bond.type)  # 'ionic'
            >>> energy = bond.compute_energy()
        """
        # Convert strings to Atom objects
        if isinstance(atom_1, str):
            atom_1 = Atom(atom_1)
        if isinstance(atom_2, str):
            atom_2 = Atom(atom_2)
        
        # Auto-determine bond type
        if bond_type == BondType.AUTO:
            bond_type = BondFactory._determine_bond_type(atom_1, atom_2)
        
        # Create appropriate bond with governance
        if bond_type == BondType.IONIC:
            return IonicBond(atom_1, atom_2, **kwargs)
        elif bond_type == BondType.COVALENT:
            return CovalentBond(atom_1, atom_2, **kwargs)
        elif bond_type == BondType.METALLIC:
            return MetallicBond(atom_1, atom_2, **kwargs)
    
    @staticmethod
    def _determine_bond_type(atom_1: 'Atom', atom_2: 'Atom') -> BondType:
        """
        Determine bond type from electronegativity difference.
        
        Rules:
            ΔEN > 1.7 → Ionic
            0.4 < ΔEN ≤ 1.7 → Polar covalent
            ΔEN ≤ 0.4 → Nonpolar covalent
            Both metals → Metallic
        """
        en1 = atom_1.properties.electronegativity
        en2 = atom_2.properties.electronegativity
        delta_en = abs(en1 - en2)
        
        if atom_1.is_metal() and atom_2.is_metal():
            return BondType.METALLIC
        elif delta_en > 1.7:
            return BondType.IONIC
        else:
            return BondType.COVALENT
    
    @staticmethod
    def create_molecule(
        atoms: List[Union[str, 'Atom']],
        geometry: Union[str, np.ndarray] = 'auto',
        bond_types: Optional[List[BondType]] = None
    ) -> 'Molecule':
        """
        Create a molecule with multiple bonds.
        
        Example:
            >>> mol = BondFactory.create_molecule(
            ...     atoms=['H', 'O', 'H'],
            ...     geometry='water'
            ... )
            >>> mol.analyze_all_bonds()
        """
        # Parse atoms
        atom_objects = [Atom(a) if isinstance(a, str) else a for a in atoms]
        
        # Determine geometry
        if geometry == 'auto':
            positions = self._optimize_geometry(atom_objects)
        elif isinstance(geometry, str):
            positions = self._load_geometry_template(geometry)
        else:
            positions = geometry
        
        # Assign positions
        for atom, pos in zip(atom_objects, positions):
            atom.position = pos
        
        # Create bonds
        bonds = []
        connectivity = self._determine_connectivity(atom_objects)
        
        for (i, j), bond_type in zip(connectivity, bond_types or [BondType.AUTO]*len(connectivity)):
            bond = BondFactory.create_bond(
                atom_objects[i],
                atom_objects[j],
                bond_type=bond_type
            )
            bonds.append(bond)
        
        return Molecule(atoms=atom_objects, bonds=bonds)
```

```python
# kanad/bonds/ionic_bond.py

class IonicBond:
    """
    Ionic bond with automatic governance.
    """
    
    def __init__(
        self,
        donor: 'Atom',
        acceptor: 'Atom',
        distance: Optional[float] = None
    ):
        self.donor = donor
        self.acceptor = acceptor
        self.distance = distance or self._equilibrium_distance()
        
        # Create governance protocol
        self.governance = IonicGovernanceProtocol()
        
        # Build Hamiltonian
        self.hamiltonian = IonicHamiltonian(donor, acceptor)
        
        # Create mapper
        self.mapper = JordanWignerMapper()
        
    def compute_energy(
        self,
        method: str = 'VQE',
        backend: Optional['Backend'] = None
    ) -> Dict[str, Any]:
        """
        Compute bond energy using specified method.
        
        Args:
            method: 'VQE', 'QPE', or 'exact'
            backend: Quantum backend (default: statevector simulator)
            
        Returns:
            Dictionary with energy and analysis
        """
        if method == 'VQE':
            vqe = GovernanceVQE(
                self.hamiltonian,
                self.governance,
                backend=backend
            )
            result = vqe.compute_energy()
            
        elif method == 'exact':
            result = {'energy': self.hamiltonian.get_fci_energy()[0]}
        
        # Add bonding analysis
        result['bond_analysis'] = self.analyze()
        
        return result
    
    def analyze(self) -> Dict[str, Any]:
        """
        Analyze ionic bond properties.
        """
        analyzer = BondingAnalyzer(self, self.hamiltonian)
        
        return {
            'bond_type': 'ionic',
            'charge_transfer': self._compute_charge_transfer(),
            'ionic_character': 1.0,
            'covalent_character': 0.0,
            'dipole_moment': analyzer.compute_dipole_moment(),
            'bond_length': self.distance,
            'lattice_energy': self._compute_lattice_energy()
        }
    
    def _compute_charge_transfer(self) -> float:
        """Amount of electron transferred (0 to 1)."""
        # From Mulliken analysis or density matrix
        pass
    
    def visualize_circuit(self) -> None:
        """Display the governed quantum circuit."""
        ansatz = self.governance.construct_ansatz(self.hamiltonian.representation)
        print(ansatz.draw('mpl'))
```

```python
# kanad/bonds/covalent_bond.py

class CovalentBond:
    """
    Covalent bond with orbital hybridization.
    """
    
    def __init__(
        self,
        atom_1: 'Atom',
        atom_2: 'Atom',
        hybridization: str = 'sp3',
        bond_order: int = 1
    ):
        self.atom_1 = atom_1
        self.atom_2 = atom_2
        self.hybridization = hybridization
        self.bond_order = bond_order
        
        # Governance
        self.governance = CovalentGovernanceProtocol()
        
        # Build Hamiltonian in hybrid orbital basis
        self.hamiltonian = CovalentHamiltonian(
            atoms=[atom_1, atom_2],
            bonds=[self],
            hybridization=hybridization
        )
        
        # Custom mapper for MO pairs
        self.mapper = HybridOrbitalMapper(self.hamiltonian.mo_pairs)
        
    def compute_energy(
        self,
        method: str = 'VQE',
        backend: Optional['Backend'] = None
    ) -> Dict[str, Any]:
        """Compute covalent bond energy."""
        vqe = GovernanceVQE(self.hamiltonian, self.governance, backend)
        result = vqe.compute_energy()
        
        result['bond_analysis'] = self.analyze()
        return result
    
    def analyze(self) -> Dict[str, Any]:
        """
        Analyze covalent bond properties.
        """
        analyzer = BondingAnalyzer(self, self.hamiltonian)
        
        # Get MO energies
        mo_data = self.hamiltonian.get_bonding_antibonding_split(0)
        
        return {
            'bond_type': 'covalent',
            'hybridization': self.hybridization,
            'bond_order': analyzer.compute_bond_order(),
            'bonding_energy': mo_data['bonding_energy'],
            'antibonding_energy': mo_data['antibonding_energy'],
            'mo_splitting': mo_data['splitting'],
            'overlap': self._compute_overlap(),
            'ionic_character': analyzer.analyze_bonding_character()['ionic'],
            'covalent_character': analyzer.analyze_bonding_character()['covalent']
        }
    
    def visualize_orbitals(self) -> None:
        """
        Visualize bonding and antibonding molecular orbitals.
        """
        mo_data = self.hamiltonian.get_bonding_antibonding_split(0)
        
        # 3D orbital plot
        fig = plt.figure(figsize=(12, 5))
        
        # Bonding MO
        ax1 = fig.add_subplot(121, projection='3d')
        self._plot_orbital(mo_data['bonding_character'], ax1, title='Bonding')
        
        # Antibonding MO
        ax2 = fig.add_subplot(122, projection='3d')
        self._plot_orbital(mo_data['antibonding_character'], ax2, title='Antibonding')
        
        plt.show()
```

```python
# kanad/bonds/metallic_bond.py

class MetallicBond:
    """
    Metallic bond with delocalized electrons.
    """
    
    def __init__(
        self,
        atoms: List['Atom'],
        lattice_type: str = '1d_chain'
    ):
        self.atoms = atoms
        self.lattice_type = lattice_type
        self.n_atoms = len(atoms)
        
        # Governance
        self.governance = MetallicGovernanceProtocol()
        
        # Build band Hamiltonian
        self.hamiltonian = MetallicHamiltonian(
            lattice_type=lattice_type,
            n_sites=self.n_atoms,
            hopping_parameter=self._compute_hopping(),
            on_site_energy=self._compute_on_site()
        )
        
        # k-space mapper
        self.mapper = MomentumSpaceMapper(n_k_points=self.n_atoms)
        
    def compute_energy(
        self,
        method: str = 'VQE',
        backend: Optional['Backend'] = None
    ) -> Dict[str, Any]:
        """Compute metallic system energy."""
        vqe = GovernanceVQE(self.hamiltonian, self.governance, backend)
        result = vqe.compute_energy()
        
        result['bond_analysis'] = self.analyze()
        return result
    
    def analyze(self) -> Dict[str, Any]:
        """
        Analyze metallic bonding properties.
        """
        # Compute band structure
        bands = self.hamiltonian.compute_band_structure()
        
        # Fermi surface
        fermi_data = self.hamiltonian.compute_fermi_surface()
        
        # DOS
        dos_data = self.hamiltonian.compute_dos(
            energy_range=(bands['energies'].min(), bands['energies'].max())
        )
        
        return {
            'bond_type': 'metallic',
            'lattice_type': self.lattice_type,
            'n_atoms': self.n_atoms,
            'band_structure': bands,
            'fermi_energy': fermi_data['fermi_energy'],
            'dos_at_fermi': fermi_data['dos_at_fermi'],
            'bandwidth': bands['energies'].max() - bands['energies'].min(),
            'conductivity': self._estimate_conductivity(dos_data)
        }
    
    def visualize_bands(self) -> None:
        """
        Visualize electronic band structure.
        """
        bands = self.hamiltonian.compute_band_structure()
        
        plt.figure(figsize=(10, 6))
        plt.plot(bands['k_points'], bands['energies'], 'b-', linewidth=2)
        plt.axhline(self.hamiltonian.compute_fermi_energy(), 
                   color='r', linestyle='--', label='Fermi level')
        plt.xlabel('k (momentum)')
        plt.ylabel('Energy')
        plt.title(f'{self.lattice_type} Band Structure')
        plt.legend()
        plt.grid(True)
        plt.show()
    
    def visualize_fermi_surface(self) -> None:
        """Visualize Fermi surface in k-space."""
        fermi_data = self.hamiltonian.compute_fermi_surface()
        # Plot Fermi surface...
```

---

## **PHASE 9: Testing & Validation Framework** 

```python
# tests/validation/test_known_molecules.py

import pytest
import numpy as np
from kanad.bonds import BondFactory

class TestKnownMolecules:
    """
    Validate against experimentally known values.
    """
    
    def test_h2_molecule(self):
        """
        Test H2 covalent bond.
        
        Known values:
            - Bond length: 0.74 Å
            - Dissociation energy: 4.52 eV
            - FCI energy: -1.174 Ha (STO-3G)
        """
        bond = BondFactory.create_bond('H', 'H', bond_type='covalent')
        result = bond.compute_energy(method='exact')
        
        # Check energy (within chemical accuracy ~1.6 mHa)
        assert abs(result['energy'] + 1.174) < 0.002
        
        # Check bond length
        assert abs(result['bond_analysis']['bond_length'] - 0.74) < 0.05
    
    def test_nacl_ionic(self):
        """
        Test NaCl ionic bond.
        
        Known values:
            - Lattice energy: -787 kJ/mol
            - Bond length: 2.36 Å
            - Charge transfer: ~0.85 e
        """
        bond = BondFactory.create_bond('Na', 'Cl', bond_type='ionic')
        result = bond.compute_energy()
        
        # Check charge transfer
        assert result['bond_analysis']['charge_transfer'] > 0.7
        
        # Check ionic character
        assert result['bond_analysis']['ionic_character'] > 0.9
    
    def test_metallic_sodium(self):
        """
        Test metallic sodium lattice.
        
        Known values:
            - BCC structure
            - Cohesive energy: -1.11 eV/atom
            - Fermi energy: 3.1 eV
        """
        atoms = [Atom('Na') for _ in range(8)]
        bond = MetallicBond(atoms, lattice_type='bcc')
        result = bond.compute_energy()
        
        # Check Fermi energy
        E_F = result['bond_analysis']['fermi_energy']
        assert abs(E_F - 3.1) < 0.5
```

```python
# tests/integration/test_governance_protocols.py

class TestGovernanceProtocols:
    """
    Test that governance protocols enforce physical constraints.
    """
    
    def test_covalent_creates_bell_pairs(self):
        """Covalent governance must create Bell-pair entanglement."""
        bond = CovalentBond(Atom('C'), Atom('C'))
        ansatz = bond.governance.construct_ansatz(bond.hamiltonian.representation)
        
        # Check circuit creates Bell states
        assert self._has_bell_pair_structure(ansatz)
    
    def test_metallic_creates_ghz_states(self):
        """Metallic governance must create collective entanglement."""
        bond = MetallicBond([Atom('Na')] * 4)
        ansatz = bond.governance.construct_ansatz(bond.hamiltonian.representation)
        
        # Check for GHZ-like entanglement
        assert self._has_collective_entanglement(ansatz)
    
    def test_ionic_minimal_entanglement(self):
        """Ionic governance should minimize entanglement."""
        bond = IonicBond(Atom('Na'), Atom('Cl'))
        ansatz = bond.governance.construct_ansatz(bond.hamiltonian.representation)
        
        # Compute entanglement entropy
        entropy = self._compute_entanglement_entropy(ansatz)
        assert entropy < 1.0  # Should be close to 0
```

---

## **PHASE 10: Research Workflows & Documentation** 

```python
# kanad/research/workflows/bond_comparison.py

class BondComparisonWorkflow:
    """
    Compare different bonding types for research.
    """
    
    @staticmethod
    def compare_bond_types(atom_1: str, atom_2: str) -> Dict[str, Any]:
        """
        Compare ionic vs covalent vs metallic for same atom pair.
        
        Example:
            >>> results = BondComparisonWorkflow.compare_bond_types('Na', 'Cl')
            >>> print(results['best_bond_type'])  # 'ionic'
        """
        results = {}
        
        for bond_type in ['ionic', 'covalent', 'metallic']:
            try:
                bond = BondFactory.create_bond(atom_1, atom_2, bond_type=bond_type)
                energy_data = bond.compute_energy()
                
                results[bond_type] = {
                    'energy': energy_data['energy'],
                    'circuit_depth': energy_data['optimal_circuit'].depth(),
                    'n_parameters': len(energy_data['optimal_parameters']),
                    'analysis': energy_data['bond_analysis']
                }
            except Exception as e:
                results[bond_type] = {'error': str(e)}
        
        # Determine most stable
        energies = {k: v['energy'] for k, v in results.items() if 'energy' in v}
        best = min(energies, key=energies.get)
        
        results['best_bond_type'] = best
        results['energy_differences'] = {
            k: v - energies[best] for k, v in energies.items()
        }
        
        return results
```

```python
# examples/01_simple_ionic_bond.py

"""
Example 1: Simple Ionic Bond (NaCl)
====================================

This example demonstrates how to create and analyze an ionic bond.
"""

from kanad.bonds import BondFactory
import matplotlib.pyplot as plt

# Create NaCl bond
print("Creating Na-Cl ionic bond...")
bond = BondFactory.create_bond('Na', 'Cl')

print(f"Bond type automatically determined: {bond.__class__.__name__}")

# Compute energy
print("\nComputing bond energy with VQE...")
result = bond.compute_energy(method='VQE')

print(f"Ground state energy: {result['energy']:.6f} Hartree")
print(f"Charge transfer: {result['bond_analysis']['charge_transfer']:.3f} e")
print(f"Ionic character: {result['bond_analysis']['ionic_character']:.1%}")

# Visualize circuit
print("\nQuantum circuit structure:")
bond.visualize_circuit()

# Compare with exact
print("\nComparing with exact FCI calculation...")
exact_result = bond.compute_energy(method='exact')
error = abs(result['energy'] - exact_result['energy'])
print(f"VQE error: {error * 1000:.2f} mHartree")
```

```python
# examples/02_covalent_bond_analysis.py

"""
Example 2: Covalent Bond Analysis (H2)
=======================================

Demonstrates orbital hybridization and MO formation.
"""

from kanad.bonds import BondFactory
import numpy as np

# Create H2 molecule
print("Creating H2 covalent bond...")
bond = BondFactory.create_bond('H', 'H', bond_type='covalent')

# Compute energy
result = bond.compute_energy()

print(f"\nBond Energy: {result['energy']:.6f} Hartree")
print(f"Bond Order: {result['bond_analysis']['bond_order']:.2f}")
print(f"MO Splitting: {result['bond_analysis']['mo_splitting']:.6f} Hartree")

# Visualize molecular orbitals
print("\nVisualizing bonding and antibonding orbitals...")
bond.visualize_orbitals()

# Analyze bonding character
character = result['bond_analysis']
print(f"\nBonding Character:")
print(f"  Covalent: {character['covalent_character']:.1%}")
print(f"  Ionic: {character['ionic_character']:.1%}")
```

```python
# examples/03_metallic_band_structure.py

"""
Example 3: Metallic Band Structure
===================================

Demonstrates delocalized bonding and band formation.
"""

from kanad.bonds import MetallicBond, Atom

# Create sodium chain
print("Creating metallic sodium chain (8 atoms)...")
atoms = [Atom('Na') for _ in range(8)]
bond = MetallicBond(atoms, lattice_type='1d_chain')

# Compute energy
result = bond.compute_energy()

print(f"\nTotal Energy: {result['energy']:.6f} Hartree")
print(f"Fermi Energy: {result['bond_analysis']['fermi_energy']:.6f} Hartree")
print(f"Bandwidth: {result['bond_analysis']['bandwidth']:.6f} Hartree")
print(f"DOS at Fermi level: {result['bond_analysis']['dos_at_fermi']:.4f}")

# Visualize band structure
print("\nPlotting band structure...")
bond.visualize_bands()

# Visualize Fermi surface
bond.visualize_fermi_surface()
```

---

## **FINAL ARCHITECTURE SUMMARY**

### **What Makes This Framework Innovative**

1. **Multi-Representation**: Not just second quantization
   - Ionic: Standard fermion operators
   - Covalent: LCAO/hybrid orbital basis
   - Metallic: Momentum-space (Bloch waves)

2. **Governance Layer**: Physical principles enforce circuit design
   - Ionic: Minimal entanglement (local transfer)
   - Covalent: Bell pairs (orbital hybridization)
   - Metallic: GHZ states (collective delocalization)

3. **Custom Mappers**: Optimized for each bonding type
   - Ionic: Jordan-Wigner
   - Covalent: MO-pair encoding
   - Metallic: k-space with QFT

4. **Automatic Bond Classification**: Users just specify atoms

5. **Research-Ready**: Pre-built workflows for common experiments

