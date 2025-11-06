"""
Configuration-based quantum state representation for Hi-VQE.

Hi-VQE uses configuration sampling instead of measuring all Pauli terms:
- Sample bitstrings from quantum state (Z basis measurement only)
- Filter valid configurations (correct electron count, spin, etc.)
- Build configuration subspace
- Diagonalize Hamiltonian in subspace (classical!)
- Generate excitations from important configurations

This reduces measurements from ~15,000 Pauli terms to 1 Z measurement per iteration.
"""

import numpy as np
from typing import List, Tuple, Set, Optional, Dict
import logging
from collections import Counter

logger = logging.getLogger(__name__)


class Configuration:
    """
    Represents a single Slater determinant configuration.

    A configuration is a bitstring like |110010⟩ representing occupation of spin orbitals.
    """

    def __init__(self, bitstring: str, n_qubits: int = None):
        """
        Initialize configuration from bitstring.

        Args:
            bitstring: Binary string like '110010' (1 = occupied, 0 = empty)
            n_qubits: Total number of qubits (pads bitstring if needed)
        """
        self.bitstring = bitstring

        # Pad if necessary
        if n_qubits is not None and len(bitstring) < n_qubits:
            self.bitstring = bitstring + '0' * (n_qubits - len(bitstring))

        self.n_qubits = len(self.bitstring)
        self.n_electrons = self.bitstring.count('1')

    def __str__(self) -> str:
        return f"|{self.bitstring}⟩"

    def __repr__(self) -> str:
        return f"Configuration('{self.bitstring}', n_electrons={self.n_electrons})"

    def __eq__(self, other) -> bool:
        return self.bitstring == other.bitstring

    def __hash__(self) -> int:
        return hash(self.bitstring)

    def to_int(self) -> int:
        """Convert bitstring to integer."""
        return int(self.bitstring, 2)

    @classmethod
    def from_int(cls, value: int, n_qubits: int):
        """Create configuration from integer."""
        bitstring = format(value, f'0{n_qubits}b')
        return cls(bitstring, n_qubits)


class ConfigurationSubspace:
    """
    Manages a subspace of configurations for Hi-VQE.

    The subspace grows as we sample more configurations and generate excitations.
    """

    def __init__(self, n_qubits: int, n_electrons: int, protocol=None):
        """
        Initialize configuration subspace.

        Args:
            n_qubits: Number of qubits (spin orbitals)
            n_electrons: Number of electrons (for filtering)
            protocol: Optional governance protocol for physics-based filtering
        """
        self.n_qubits = n_qubits
        self.n_electrons = n_electrons
        self.protocol = protocol

        # Track configurations as set (for fast membership testing)
        self._configs: Set[Configuration] = set()

        # Track order for indexing
        self._config_list: List[Configuration] = []

        logger.info(f"ConfigurationSubspace initialized: {n_qubits} qubits, {n_electrons} electrons")

    def add_config(self, config: Configuration) -> bool:
        """
        Add configuration to subspace (if valid).

        Returns:
            True if added (new), False if already present
        """
        # Check if already in subspace
        if config in self._configs:
            return False

        # Validate configuration
        if not self.is_valid(config):
            return False

        # Add to subspace
        self._configs.add(config)
        self._config_list.append(config)

        return True

    def add_configs(self, configs: List[Configuration]) -> int:
        """
        Add multiple configurations.

        Returns:
            Number of new configurations added
        """
        count = 0
        for config in configs:
            if self.add_config(config):
                count += 1

        if count > 0:
            logger.debug(f"Added {count} new configurations (total: {len(self)})")

        return count

    def is_valid(self, config: Configuration) -> bool:
        """
        Check if configuration is valid.

        Basic checks:
        - Correct number of electrons

        Governance checks (if protocol provided):
        - Spin symmetry
        - Charge conservation
        - Physics-specific rules
        """
        # Basic check: correct electron count
        if config.n_electrons != self.n_electrons:
            return False

        # Governance check (if protocol provided)
        if self.protocol is not None:
            if hasattr(self.protocol, 'is_valid_configuration'):
                if not self.protocol.is_valid_configuration(config.bitstring):
                    return False

        return True

    def get_hf_configuration(self) -> Configuration:
        """
        Get Hartree-Fock configuration (lowest orbitals occupied).

        For n_electrons=4, n_qubits=8: |11110000⟩
        """
        bitstring = '1' * self.n_electrons + '0' * (self.n_qubits - self.n_electrons)
        return Configuration(bitstring, self.n_qubits)

    def generate_single_excitations(self, config: Configuration) -> List[Configuration]:
        """
        Generate single excitations from a configuration.

        Single excitation: Move one electron from occupied to virtual orbital.

        If governance protocol is provided, generate only physically meaningful excitations.
        """
        excitations = []

        # Use governance protocol if available
        if self.protocol is not None and hasattr(self.protocol, 'generate_single_excitations'):
            bitstrings = self.protocol.generate_single_excitations(config.bitstring)
            excitations = [Configuration(bs, self.n_qubits) for bs in bitstrings]
            return excitations

        # Otherwise, generate all single excitations
        bits = list(config.bitstring)

        for i in range(self.n_qubits):
            if bits[i] == '1':  # Occupied orbital
                for j in range(self.n_qubits):
                    if bits[j] == '0':  # Virtual orbital
                        # Create excitation
                        new_bits = bits.copy()
                        new_bits[i] = '0'
                        new_bits[j] = '1'
                        excitations.append(Configuration(''.join(new_bits), self.n_qubits))

        return excitations

    def generate_double_excitations(self, config: Configuration) -> List[Configuration]:
        """
        Generate double excitations from a configuration.

        Double excitation: Move two electrons from occupied to virtual orbitals.
        """
        excitations = []

        # Use governance protocol if available
        if self.protocol is not None and hasattr(self.protocol, 'generate_double_excitations'):
            bitstrings = self.protocol.generate_double_excitations(config.bitstring)
            excitations = [Configuration(bs, self.n_qubits) for bs in bitstrings]
            return excitations

        # Otherwise, generate all double excitations
        bits = list(config.bitstring)
        occupied = [i for i, b in enumerate(bits) if b == '1']
        virtual = [i for i, b in enumerate(bits) if b == '0']

        for i in range(len(occupied)):
            for j in range(i + 1, len(occupied)):
                for a in range(len(virtual)):
                    for b in range(a + 1, len(virtual)):
                        # Create double excitation
                        new_bits = bits.copy()
                        new_bits[occupied[i]] = '0'
                        new_bits[occupied[j]] = '0'
                        new_bits[virtual[a]] = '1'
                        new_bits[virtual[b]] = '1'
                        excitations.append(Configuration(''.join(new_bits), self.n_qubits))

        return excitations

    def prune(self, amplitudes: np.ndarray, threshold: float = 1e-4) -> int:
        """
        Prune configurations with low amplitudes.

        Args:
            amplitudes: Amplitudes of configurations in subspace
            threshold: Remove configurations with |amplitude| < threshold

        Returns:
            Number of configurations removed
        """
        assert len(amplitudes) == len(self._config_list), "Amplitude count mismatch"

        # Find configurations to keep
        keep_indices = [i for i, amp in enumerate(amplitudes) if abs(amp) >= threshold]

        # Remove low-amplitude configs
        removed = len(self) - len(keep_indices)

        if removed > 0:
            self._config_list = [self._config_list[i] for i in keep_indices]
            self._configs = set(self._config_list)
            logger.info(f"Pruned {removed} low-amplitude configurations (threshold={threshold})")

        return removed

    def __len__(self) -> int:
        return len(self._config_list)

    def __getitem__(self, idx: int) -> Configuration:
        return self._config_list[idx]

    def __iter__(self):
        return iter(self._config_list)

    @property
    def configs(self) -> List[Configuration]:
        """Get list of configurations."""
        return self._config_list


def sample_configurations_from_statevector(
    statevector: np.ndarray,
    n_shots: int = 1000,
    n_electrons: int = None
) -> List[Tuple[str, int]]:
    """
    Sample configurations from statevector (simulates Z basis measurement).

    Args:
        statevector: Quantum state vector
        n_shots: Number of samples
        n_electrons: Filter by electron count (optional)

    Returns:
        List of (bitstring, count) tuples
    """
    n_qubits = int(np.log2(len(statevector)))

    # Compute probabilities
    probabilities = np.abs(statevector) ** 2

    # Sample indices
    indices = np.random.choice(len(statevector), size=n_shots, p=probabilities)

    # Convert to bitstrings
    bitstrings = [format(idx, f'0{n_qubits}b') for idx in indices]

    # Filter by electron count if provided
    if n_electrons is not None:
        bitstrings = [bs for bs in bitstrings if bs.count('1') == n_electrons]

    # Count occurrences
    counts = Counter(bitstrings)

    return list(counts.items())


def sample_configurations_from_counts(
    counts: Dict[str, int],
    n_electrons: int = None
) -> List[Tuple[str, int]]:
    """
    Extract configurations from measurement counts (real quantum backend).

    Args:
        counts: Dictionary of {bitstring: count} from quantum measurement
        n_electrons: Filter by electron count (optional)

    Returns:
        List of (bitstring, count) tuples
    """
    filtered = []

    for bitstring, count in counts.items():
        # Filter by electron count
        if n_electrons is None or bitstring.count('1') == n_electrons:
            filtered.append((bitstring, count))

    return filtered
