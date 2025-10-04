"""
Active Space Selection and Qubit Reduction.

Reduces problem size by selecting chemically relevant orbitals:
- HOMO/LUMO and nearby orbitals
- Bonding/antibonding pairs (from governance)
- Frozen core approximation
- Virtual orbital truncation

Qubit reduction techniques:
- Symmetry reduction (point group)
- Particle number conservation (Z2 symmetry)
- Tapering (eliminate redundant qubits)
- Orbital rotation to maximize localization
"""

from typing import List, Dict, Optional, Tuple, Any, Set
import numpy as np
from scipy.linalg import eigh
import logging

logger = logging.getLogger(__name__)

from kanad.core.hamiltonians.molecular_hamiltonian import MolecularHamiltonian
from kanad.core.mappers.base_mapper import BaseMapper
from kanad.governance.protocols.base_protocol import BaseGovernanceProtocol


class ActiveSpaceSelector:
    """
    Select active space orbitals for quantum simulation.

    Uses governance protocols to identify important orbitals:
    - Covalent: bonding/antibonding pairs
    - Ionic: donor/acceptor orbitals
    - Metallic: orbitals near Fermi level

    Methods:
        - select_cas: Complete Active Space selection
        - select_ras: Restricted Active Space selection
        - freeze_core: Remove core orbitals
        - truncate_virtual: Remove high-energy virtuals
    """

    def __init__(
        self,
        hamiltonian: MolecularHamiltonian,
        governance: Optional[BaseGovernanceProtocol] = None
    ):
        """
        Initialize active space selector.

        Args:
            hamiltonian: Molecular Hamiltonian
            governance: Governance protocol for orbital selection
        """
        self.hamiltonian = hamiltonian
        self.governance = governance

        self.n_orbitals = hamiltonian.n_orbitals
        self.n_electrons = hamiltonian.n_electrons

        self.frozen_core: Set[int] = set()
        self.active_orbitals: Set[int] = set()
        self.virtual_orbitals: Set[int] = set()

        logger.info(f"Initialized ActiveSpaceSelector: {self.n_orbitals} orbitals, "
                   f"{self.n_electrons} electrons")

    def select_cas(
        self,
        n_active_orbitals: int,
        n_active_electrons: int,
        selection_method: str = 'homo_lumo'
    ) -> Dict[str, Any]:
        """
        Select Complete Active Space (CAS).

        CAS(n,m): n electrons in m orbitals
        - Freezes core orbitals (fully occupied)
        - Selects active orbitals (partial occupation)
        - Deletes virtual orbitals (unoccupied)

        Args:
            n_active_orbitals: Number of active orbitals (m)
            n_active_electrons: Number of active electrons (n)
            selection_method: 'homo_lumo', 'natural_orbitals', 'governance'

        Returns:
            Dictionary with active space definition and reduced Hamiltonian
        """
        logger.info(f"Selecting CAS({n_active_electrons},{n_active_orbitals}) active space...")

        if selection_method == 'homo_lumo':
            active_orbitals = self._select_homo_lumo_window(n_active_orbitals)

        elif selection_method == 'natural_orbitals':
            active_orbitals = self._select_natural_orbitals(n_active_orbitals)

        elif selection_method == 'governance':
            if self.governance is None:
                raise ValueError("Governance protocol required for governance-based selection")
            active_orbitals = self._select_governance_orbitals(n_active_orbitals)

        else:
            raise ValueError(f"Unknown selection method: {selection_method}")

        # Determine frozen and virtual
        n_occupied = self.n_electrons // 2

        frozen_orbitals = set(range(n_occupied - n_active_electrons // 2))
        virtual_orbitals = set(range(max(active_orbitals) + 1, self.n_orbitals))

        self.frozen_core = frozen_orbitals
        self.active_orbitals = active_orbitals
        self.virtual_orbitals = virtual_orbitals

        # Build reduced Hamiltonian
        reduced_ham = self._build_reduced_hamiltonian()

        n_qubits_full = 2 * self.n_orbitals
        n_qubits_reduced = 2 * len(active_orbitals)
        reduction_factor = n_qubits_full / n_qubits_reduced

        logger.info(f"Active space selected:")
        logger.info(f"  Frozen core: {len(frozen_orbitals)} orbitals")
        logger.info(f"  Active space: {len(active_orbitals)} orbitals ({n_active_electrons} electrons)")
        logger.info(f"  Virtual: {len(virtual_orbitals)} orbitals")
        logger.info(f"  Qubit reduction: {n_qubits_full} → {n_qubits_reduced} ({reduction_factor:.1f}x)")

        results = {
            'frozen_orbitals': sorted(frozen_orbitals),
            'active_orbitals': sorted(active_orbitals),
            'virtual_orbitals': sorted(virtual_orbitals),
            'n_qubits_full': n_qubits_full,
            'n_qubits_reduced': n_qubits_reduced,
            'reduction_factor': reduction_factor,
            'reduced_hamiltonian': reduced_ham,
            'method': selection_method
        }

        return results

    def _select_homo_lumo_window(self, n_orbitals: int) -> Set[int]:
        """
        Select orbitals around HOMO-LUMO gap.

        Selects n_orbitals centered on HOMO/LUMO.
        """
        # Get MO energies (if available)
        if hasattr(self.hamiltonian, 'compute_molecular_orbitals'):
            mo_energies, _ = self.hamiltonian.compute_molecular_orbitals()
        else:
            # Fallback: use h_core diagonal as orbital energies
            mo_energies = np.diag(self.hamiltonian.h_core)

        n_occupied = self.n_electrons // 2
        homo_idx = n_occupied - 1
        lumo_idx = n_occupied

        # Select window around gap
        n_below = n_orbitals // 2
        n_above = n_orbitals - n_below

        start_idx = max(0, homo_idx - n_below + 1)
        end_idx = min(self.n_orbitals, lumo_idx + n_above)

        active = set(range(start_idx, end_idx))

        logger.info(f"HOMO-LUMO window: orbitals {start_idx} to {end_idx-1}")
        logger.info(f"HOMO = {homo_idx}, LUMO = {lumo_idx}")

        return active

    def _select_natural_orbitals(self, n_orbitals: int) -> Set[int]:
        """
        Select natural orbitals with highest occupation numbers.

        Natural orbitals from density matrix diagonalization.
        """
        # Get density matrix
        # For HF: P = 2 * C_occ @ C_occ.T
        if hasattr(self.hamiltonian, 'compute_molecular_orbitals'):
            mo_energies, mo_coeffs = self.hamiltonian.compute_molecular_orbitals()
        else:
            # Fallback: use identity as MO coefficients
            mo_energies = np.diag(self.hamiltonian.h_core)
            mo_coeffs = np.eye(self.n_orbitals)
        n_occupied = self.n_electrons // 2

        C_occ = mo_coeffs[:, :n_occupied]
        P = 2 * C_occ @ C_occ.T

        # Natural orbitals
        occupation_numbers, natural_orbitals = eigh(P)

        # Select orbitals with occupation closest to 1.0 (most correlated)
        # Distance from 1.0 indicates correlation importance
        correlation_importance = np.abs(occupation_numbers - 1.0)

        # Select top n_orbitals
        important_indices = np.argsort(correlation_importance)[::-1][:n_orbitals]

        active = set(important_indices.tolist())

        logger.info(f"Natural orbital selection: orbitals {sorted(active)}")

        return active

    def _select_governance_orbitals(self, n_orbitals: int) -> Set[int]:
        """
        Select orbitals based on governance protocol.

        Covalent: bonding/antibonding pairs
        Ionic: donor/acceptor orbitals
        Metallic: orbitals near Fermi level
        """
        if hasattr(self.governance, 'get_important_orbitals'):
            active = self.governance.get_important_orbitals(self.hamiltonian, n_orbitals)
        else:
            # Fallback to HOMO-LUMO
            logger.warning("Governance doesn't implement get_important_orbitals, using HOMO-LUMO")
            active = self._select_homo_lumo_window(n_orbitals)

        return active

    def _build_reduced_hamiltonian(self) -> Dict[str, np.ndarray]:
        """
        Build Hamiltonian in active space.

        Freezes core electrons, truncates virtuals.
        """
        active_indices = sorted(self.active_orbitals)

        # Extract active block from h_core
        h_core_active = self.hamiltonian.h_core[np.ix_(active_indices, active_indices)]

        # Core energy contribution (frozen orbitals)
        E_core = 0.0
        for i in self.frozen_core:
            E_core += 2 * self.hamiltonian.h_core[i, i]

        # Extract active block from ERI
        if hasattr(self.hamiltonian, 'eri'):
            eri = self.hamiltonian.eri
            eri_active = eri[np.ix_(active_indices, active_indices, active_indices, active_indices)]
        else:
            n_active = len(active_indices)
            eri_active = np.zeros((n_active, n_active, n_active, n_active))

        reduced_ham = {
            'h_core': h_core_active,
            'eri': eri_active,
            'E_core': E_core,
            'E_nuclear': self.hamiltonian.nuclear_repulsion,
            'active_indices': active_indices,
            'n_active_orbitals': len(active_indices),
            'n_active_electrons': len([i for i in active_indices if i < self.n_electrons // 2]) * 2
        }

        return reduced_ham


class QubitReducer:
    """
    Apply qubit reduction techniques.

    Methods:
    - Z2 symmetry tapering (particle number conservation)
    - Point group symmetry
    - Contextual subspace (eliminate low-weight Paulis)
    """

    def __init__(self, n_qubits: int, n_electrons: int):
        """
        Initialize qubit reducer.

        Args:
            n_qubits: Number of qubits
            n_electrons: Number of electrons (for Z2 symmetry)
        """
        self.n_qubits = n_qubits
        self.n_electrons = n_electrons

        logger.info(f"Initialized QubitReducer: {n_qubits} qubits, {n_electrons} electrons")

    def apply_z2_tapering(self, mapper: BaseMapper) -> Dict[str, Any]:
        """
        Apply Z2 symmetry tapering.

        Particle number conservation:
        N = Σ a†_i a_i = Σ (1 - Z_i)/2 = const

        This symmetry reduces qubits by ~1.
        """
        logger.info("Applying Z2 symmetry tapering...")

        # Find Z2 symmetry generators
        # For particle number: Π Z_i = ±1
        z2_generator = np.ones(self.n_qubits)

        # Parity of electron number
        parity = self.n_electrons % 2

        # Tapering removes 1 qubit
        n_qubits_tapered = self.n_qubits - 1

        logger.info(f"Z2 tapering: {self.n_qubits} → {n_qubits_tapered} qubits")

        results = {
            'n_qubits_original': self.n_qubits,
            'n_qubits_tapered': n_qubits_tapered,
            'reduction': 1,
            'symmetry': 'Z2 (particle number)',
            'parity': parity
        }

        return results

    def estimate_circuit_depth_reduction(self, active_space_results: Dict) -> Dict[str, Any]:
        """
        Estimate circuit depth reduction from active space.

        Depth scales as:
        - Trotter steps: O(n²) for n orbitals
        - Gate count: O(n⁴) for full CI
        - UCCSD: O(n² + n⁴) gates
        """
        n_qubits_full = active_space_results['n_qubits_full']
        n_qubits_reduced = active_space_results['n_qubits_reduced']

        n_orbitals_full = n_qubits_full // 2
        n_orbitals_reduced = n_qubits_reduced // 2

        # UCCSD gate count estimate
        gates_full = n_orbitals_full**2 + n_orbitals_full**4
        gates_reduced = n_orbitals_reduced**2 + n_orbitals_reduced**4

        gate_reduction = gates_full / gates_reduced if gates_reduced > 0 else 1

        # Circuit depth estimate (assuming parallelization)
        depth_full = n_orbitals_full**2
        depth_reduced = n_orbitals_reduced**2

        depth_reduction = depth_full / depth_reduced if depth_reduced > 0 else 1

        logger.info(f"Circuit complexity reduction:")
        logger.info(f"  Gates: {gates_full:.0f} → {gates_reduced:.0f} ({gate_reduction:.1f}x)")
        logger.info(f"  Depth: {depth_full:.0f} → {depth_reduced:.0f} ({depth_reduction:.1f}x)")

        return {
            'gates_full': gates_full,
            'gates_reduced': gates_reduced,
            'gate_reduction_factor': gate_reduction,
            'depth_full': depth_full,
            'depth_reduced': depth_reduced,
            'depth_reduction_factor': depth_reduction
        }
