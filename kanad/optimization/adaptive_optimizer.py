"""
Adaptive Optimizer - Dynamic active space optimization.

Adapts active space during simulation:
- ADAPT-VQE style orbital selection
- Dynamic expansion based on gradient
- Error-based truncation
- Iterative improvement
"""

from typing import Dict, Any, List, Optional, Callable
import numpy as np
import logging

logger = logging.getLogger(__name__)

from kanad.core.hamiltonians.molecular_hamiltonian import MolecularHamiltonian
from kanad.solvers.active_space import ActiveSpaceSelector


class AdaptiveOptimizer:
    """
    Adaptive active space optimization.

    Features:
    - **ADAPT-VQE**: Grow ansatz iteratively based on gradient
    - **Dynamic Expansion**: Add orbitals when needed
    - **Error-Based Truncation**: Remove orbitals below threshold
    - **Importance Sampling**: Prioritize important orbital pairs

    Example:
        >>> optimizer = AdaptiveOptimizer(hamiltonian)
        >>> result = optimizer.adaptive_vqe(
        ...     energy_solver=vqe_solver,
        ...     max_orbitals=10,
        ...     gradient_threshold=1e-3
        ... )
    """

    def __init__(
        self,
        hamiltonian: MolecularHamiltonian,
        initial_active_orbitals: Optional[List[int]] = None
    ):
        """
        Initialize adaptive optimizer.

        Args:
            hamiltonian: Molecular Hamiltonian
            initial_active_orbitals: Starting active space (default: HOMO-LUMO)
        """
        self.hamiltonian = hamiltonian
        self.n_orbitals = hamiltonian.n_orbitals
        self.n_electrons = hamiltonian.n_electrons

        # Initialize with HOMO-LUMO if not specified
        if initial_active_orbitals is None:
            n_occ = self.n_electrons // 2
            initial_active_orbitals = [n_occ - 1, n_occ]  # HOMO, LUMO

        self.active_orbitals = set(initial_active_orbitals)

        logger.info(f"AdaptiveOptimizer initialized: starting with {len(self.active_orbitals)} orbitals")

    def adaptive_vqe(
        self,
        energy_solver: Callable,
        max_orbitals: int = 10,
        gradient_threshold: float = 1e-3,
        max_iterations: int = 20,
        energy_threshold: float = 1e-5
    ) -> Dict[str, Any]:
        """
        ADAPT-VQE style adaptive optimization.

        Iteratively:
        1. Compute energy gradient for all orbital excitations
        2. Add orbital with largest gradient
        3. Re-optimize energy
        4. Repeat until convergence

        Args:
            energy_solver: Function that computes energy for given active space
            max_orbitals: Maximum active space size
            gradient_threshold: Stop when max gradient < threshold
            max_iterations: Maximum iterations
            energy_threshold: Energy convergence threshold

        Returns:
            Optimization results with final active space
        """
        logger.info("="*70)
        logger.info("ADAPTIVE VQE OPTIMIZATION")
        logger.info("="*70)

        energies = []
        active_space_sizes = []

        energy_old = float('inf')

        for iteration in range(max_iterations):
            logger.info(f"\n--- Iteration {iteration + 1} ---")
            logger.info(f"Active orbitals: {sorted(self.active_orbitals)}")

            # Compute energy with current active space
            energy = energy_solver(list(self.active_orbitals))
            energies.append(energy)
            active_space_sizes.append(len(self.active_orbitals))

            logger.info(f"Energy: {energy:.6f} Ha")

            # Check energy convergence
            delta_energy = abs(energy - energy_old)
            if delta_energy < energy_threshold:
                logger.info(f"Energy converged (ΔE = {delta_energy:.2e} Ha)")
                break

            # Stop if max orbitals reached
            if len(self.active_orbitals) >= max_orbitals:
                logger.info(f"Maximum active space size reached ({max_orbitals})")
                break

            # Compute gradients for all possible orbital additions
            gradients = self._compute_orbital_gradients(energy_solver)

            # Find orbital with maximum gradient
            max_gradient_orbital, max_gradient = self._select_next_orbital(gradients)

            logger.info(f"Max gradient: {max_gradient:.6f} (orbital {max_gradient_orbital})")

            # Check gradient threshold
            if max_gradient < gradient_threshold:
                logger.info(f"Gradient below threshold ({gradient_threshold})")
                break

            # Add orbital
            self.active_orbitals.add(max_gradient_orbital)
            logger.info(f"Added orbital {max_gradient_orbital}")

            energy_old = energy

        # Final summary
        logger.info("\n" + "="*70)
        logger.info("ADAPTIVE VQE COMPLETE")
        logger.info("="*70)
        logger.info(f"Final active space: {sorted(self.active_orbitals)}")
        logger.info(f"Final energy: {energies[-1]:.6f} Ha")
        logger.info(f"Energy improvement: {energies[0] - energies[-1]:.6f} Ha")
        logger.info(f"Iterations: {len(energies)}")
        logger.info("="*70)

        results = {
            'final_active_orbitals': sorted(self.active_orbitals),
            'final_energy': energies[-1],
            'initial_energy': energies[0],
            'energy_improvement': energies[0] - energies[-1],
            'energies': energies,
            'active_space_sizes': active_space_sizes,
            'iterations': len(energies),
            'converged': delta_energy < energy_threshold
        }

        return results

    def importance_sampling(
        self,
        n_samples: int = 1000,
        temperature: float = 1.0
    ) -> Dict[str, Any]:
        """
        Sample important orbital excitations based on energy contributions.

        Uses Boltzmann sampling to select excitations:
        P(excitation) ∝ exp(-ΔE / kT)

        Args:
            n_samples: Number of samples
            temperature: Sampling temperature (higher = more exploration)

        Returns:
            Sampled orbital pairs and their importance
        """
        logger.info("Importance sampling for orbital selection...")

        # Compute energy contributions for all orbital pairs
        orbital_pairs = []
        energy_contributions = []

        for i in range(self.n_orbitals):
            for j in range(i + 1, self.n_orbitals):
                # Estimate energy contribution (simplified)
                # In practice: compute <ij||ij> integrals
                contribution = self._estimate_pair_energy(i, j)
                orbital_pairs.append((i, j))
                energy_contributions.append(abs(contribution))

        energy_contributions = np.array(energy_contributions)

        # Boltzmann probabilities
        probabilities = np.exp(-energy_contributions / temperature)
        probabilities /= probabilities.sum()

        # Sample orbital pairs
        sample_indices = np.random.choice(
            len(orbital_pairs),
            size=n_samples,
            p=probabilities,
            replace=True
        )

        # Count occurrences
        sampled_orbitals = set()
        orbital_importance = {}

        for idx in sample_indices:
            i, j = orbital_pairs[idx]
            sampled_orbitals.add(i)
            sampled_orbitals.add(j)

            orbital_importance[i] = orbital_importance.get(i, 0) + 1
            orbital_importance[j] = orbital_importance.get(j, 0) + 1

        # Normalize importance
        total_counts = sum(orbital_importance.values())
        for orbital in orbital_importance:
            orbital_importance[orbital] /= total_counts

        # Select top orbitals
        top_orbitals = sorted(orbital_importance.items(), key=lambda x: x[1], reverse=True)

        logger.info(f"Sampled {len(sampled_orbitals)} important orbitals")
        logger.info(f"Top 5 orbitals: {top_orbitals[:5]}")

        results = {
            'sampled_orbitals': sorted(sampled_orbitals),
            'orbital_importance': orbital_importance,
            'top_orbitals': top_orbitals[:10],
            'n_samples': n_samples,
            'temperature': temperature
        }

        return results

    def _compute_orbital_gradients(
        self,
        energy_solver: Callable
    ) -> Dict[int, float]:
        """
        Compute energy gradient for each possible orbital addition.

        Gradient ≈ ΔE when adding orbital

        Args:
            energy_solver: Energy computation function

        Returns:
            Dictionary mapping orbital → gradient
        """
        gradients = {}

        # Current energy
        E_current = energy_solver(list(self.active_orbitals))

        # Candidate orbitals (not in active space)
        candidates = set(range(self.n_orbitals)) - self.active_orbitals

        for orbital in candidates:
            # Add orbital temporarily
            test_active = list(self.active_orbitals) + [orbital]

            # Compute energy
            E_test = energy_solver(test_active)

            # Gradient = energy lowering from adding orbital
            gradient = abs(E_current - E_test)
            gradients[orbital] = gradient

        return gradients

    def _select_next_orbital(
        self,
        gradients: Dict[int, float]
    ) -> tuple[int, float]:
        """
        Select orbital with maximum gradient.

        Args:
            gradients: Orbital gradients

        Returns:
            Tuple of (orbital_index, gradient_value)
        """
        max_orbital = max(gradients.items(), key=lambda x: x[1])
        return max_orbital

    def _estimate_pair_energy(self, i: int, j: int) -> float:
        """
        Estimate energy contribution of orbital pair (i,j).

        Uses two-electron integral <ij||ij>

        Args:
            i: First orbital index
            j: Second orbital index

        Returns:
            Estimated energy contribution
        """
        # Access ERI if available
        if hasattr(self.hamiltonian, 'eri'):
            eri = self.hamiltonian.eri
            # Antisymmetrized integral <ij||ij> = <ij|ij> - <ij|ji>
            if i < eri.shape[0] and j < eri.shape[1]:
                energy = eri[i, j, i, j] - eri[i, j, j, i]
                return energy

        # Fallback: random estimate
        return np.random.rand()

    def prune_active_space(
        self,
        energy_solver: Callable,
        error_threshold: float = 1e-4
    ) -> Dict[str, Any]:
        """
        Remove unimportant orbitals from active space.

        Tests each orbital removal and keeps if error < threshold.

        Args:
            energy_solver: Energy computation function
            error_threshold: Maximum allowed error

        Returns:
            Pruning results
        """
        logger.info("Pruning active space...")

        # Reference energy
        E_ref = energy_solver(list(self.active_orbitals))

        removed_orbitals = []

        for orbital in sorted(self.active_orbitals):
            # Try removing orbital
            test_active = [orb for orb in self.active_orbitals if orb != orbital]

            if len(test_active) < 2:
                # Keep at least 2 orbitals
                continue

            # Compute energy without orbital
            E_test = energy_solver(test_active)

            error = abs(E_test - E_ref)

            if error < error_threshold:
                # Safe to remove
                self.active_orbitals.remove(orbital)
                removed_orbitals.append(orbital)
                logger.info(f"  Removed orbital {orbital} (error = {error:.2e})")

        logger.info(f"Pruned {len(removed_orbitals)} orbitals")

        results = {
            'removed_orbitals': removed_orbitals,
            'final_active_orbitals': sorted(self.active_orbitals),
            'n_removed': len(removed_orbitals)
        }

        return results
