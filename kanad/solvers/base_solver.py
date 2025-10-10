"""
Base Solver Class for Kanad Framework.

All solvers inherit from this class and work with the bonds module interface.
This provides consistent API and automatic integration with analysis and optimization.
"""

from typing import Dict, Any, Optional
from abc import ABC, abstractmethod
import logging
import numpy as np

logger = logging.getLogger(__name__)


class BaseSolver(ABC):
    """
    Base class for all quantum chemistry solvers.

    Design Philosophy:
    1. Solvers take a BOND as input (not raw Hamiltonians)
    2. Solvers automatically integrate analysis tools
    3. Solvers automatically integrate optimization tools
    4. Solvers provide rich, comprehensive results

    This makes solvers the "one-stop-shop" for users.
    """

    def __init__(
        self,
        bond: 'BaseBond',
        enable_analysis: bool = True,
        enable_optimization: bool = True,
        **kwargs
    ):
        """
        Initialize base solver.

        Args:
            bond: Bond object (from BondFactory)
            enable_analysis: Enable automatic analysis (default: True)
            enable_optimization: Enable automatic optimization (default: True)
            **kwargs: Solver-specific parameters
        """
        self.bond = bond
        self.enable_analysis = enable_analysis
        self.enable_optimization = enable_optimization

        # Extract components from bond
        self.hamiltonian = bond.hamiltonian
        self.molecule = bond.molecule
        self.atoms = bond.atoms

        # Initialize analysis tools if enabled
        if enable_analysis:
            self._init_analysis_tools()

        # Initialize optimization tools if enabled
        if enable_optimization:
            self._init_optimization_tools()

        # Storage for results
        self.results = {}

        logger.info(f"Initialized {self.__class__.__name__} for {bond.bond_type} bond")

    def _init_analysis_tools(self):
        """Initialize analysis tools from bonds module."""
        from kanad.bonds import (
            EnergyAnalyzer,
            BondingAnalyzer,
            PropertyCalculator
        )

        self.energy_analyzer = EnergyAnalyzer(self.hamiltonian)
        self.bonding_analyzer = BondingAnalyzer(self.hamiltonian)

        # PropertyCalculator may require PySCF molecule - initialize if available
        try:
            self.property_calculator = PropertyCalculator(self.hamiltonian)
        except (AttributeError, TypeError):
            logger.debug("PropertyCalculator initialization skipped (requires PySCF molecule)")
            self.property_calculator = None

        logger.debug("Analysis tools initialized")

    def _init_optimization_tools(self):
        """Initialize optimization tools from bonds module."""
        from kanad.bonds import (
            CircuitOptimizer,
            OrbitalOptimizer
        )

        # CircuitOptimizer for quantum circuits (if available)
        self.circuit_optimizer = CircuitOptimizer() if hasattr(self, 'circuit') else None

        # OrbitalOptimizer requires MO coefficients - defer initialization
        self.orbital_optimizer = None

        logger.debug("Optimization tools initialized")

    @abstractmethod
    def solve(self, **kwargs) -> Dict[str, Any]:
        """
        Solve for ground state energy and properties.

        Must be implemented by subclass.

        Returns:
            Dictionary with comprehensive results including:
                - energy: Ground state energy (Hartree)
                - converged: Convergence status
                - iterations: Number of iterations
                - analysis: Detailed analysis (if enabled)
                - optimization_stats: Optimization statistics (if enabled)
        """
        pass

    def _add_analysis_to_results(self, energy: float, density_matrix: Optional[np.ndarray] = None):
        """
        Add automatic analysis to results.

        Args:
            energy: Computed energy
            density_matrix: Density matrix (if available)
        """
        if not self.enable_analysis:
            return

        analysis = {}

        # Energy decomposition
        try:
            # EnergyAnalyzer.decompose_energy() takes only density_matrix
            analysis['energy_components'] = self.energy_analyzer.decompose_energy(
                density_matrix if density_matrix is not None else np.eye(self.hamiltonian.n_orbitals)
            )
        except Exception as e:
            logger.warning(f"Energy decomposition failed: {e}")
            analysis['energy_components'] = None

        # Bonding analysis
        try:
            # BondingAnalyzer: some methods need density_matrix, others don't
            bonding_info = {}

            if hasattr(self.bonding_analyzer, 'analyze_bonding_type'):
                bonding_info['bond_type'] = self.bonding_analyzer.analyze_bonding_type()

            if hasattr(self.bonding_analyzer, 'analyze_bond_orders') and density_matrix is not None:
                bonding_info['bond_orders'] = self.bonding_analyzer.analyze_bond_orders(density_matrix)

            analysis['bonding'] = bonding_info if bonding_info else None
        except Exception as e:
            logger.warning(f"Bonding analysis failed: {e}")
            analysis['bonding'] = None

        # Molecular properties
        try:
            if self.property_calculator is not None:
                analysis['properties'] = self.property_calculator.calculate_properties(
                    self.molecule,
                    self.hamiltonian,
                    density_matrix if density_matrix is not None else None
                )
            else:
                analysis['properties'] = None
        except Exception as e:
            logger.warning(f"Property calculation failed: {e}")
            analysis['properties'] = None

        self.results['analysis'] = analysis

    def _add_optimization_stats(self):
        """Add optimization statistics to results."""
        if not self.enable_optimization:
            return

        opt_stats = {}

        # Circuit optimization stats (if applicable)
        if hasattr(self, 'circuit') and hasattr(self, 'circuit_optimizer') and self.circuit_optimizer:
            opt_stats['circuit'] = {
                'gates_before': getattr(self, '_gates_before_opt', None),
                'gates_after': getattr(self, '_gates_after_opt', None),
                'depth_before': getattr(self, '_depth_before_opt', None),
                'depth_after': getattr(self, '_depth_after_opt', None),
            }

        # Orbital optimization stats
        opt_stats['orbitals'] = {
            'localization_applied': getattr(self, '_orbital_localization', False),
            'rotation_applied': getattr(self, '_orbital_rotation', False),
        }

        self.results['optimization_stats'] = opt_stats

    def get_reference_energy(self) -> float:
        """
        Get Hartree-Fock reference energy for comparison.

        Returns:
            HF energy (Hartree)
        """
        try:
            density_matrix, hf_energy = self.hamiltonian.solve_scf(
                max_iterations=100,
                conv_tol=1e-8,
                use_diis=True
            )
            return hf_energy
        except Exception as e:
            logger.warning(f"Could not compute HF reference: {e}")
            return None

    def validate_results(self) -> Dict[str, Any]:
        """
        Validate solver results against known checks.

        Returns:
            Validation report
        """
        validation = {
            'passed': True,
            'checks': []
        }

        # Check 1: Energy is real
        if 'energy' in self.results:
            energy = self.results['energy']
            is_real = np.isreal(energy) and not np.isnan(energy) and not np.isinf(energy)
            validation['checks'].append({
                'name': 'energy_is_real',
                'passed': is_real,
                'message': f"Energy is {'valid' if is_real else 'invalid'}: {energy}"
            })
            validation['passed'] = validation['passed'] and is_real

        # Check 2: Energy below HF (for correlated methods)
        if 'energy' in self.results and hasattr(self, '_is_correlated') and self._is_correlated:
            hf_energy = self.get_reference_energy()
            if hf_energy is not None:
                # Use 1e-5 Ha (10 μHa) tolerance for VQE numerical precision
                # Previous 1e-6 was too strict and caused false positives
                below_hf = self.results['energy'] <= hf_energy + 1e-5
                validation['checks'].append({
                    'name': 'energy_below_hf',
                    'passed': below_hf,
                    'message': f"Energy {'≤' if below_hf else '>'} HF ({self.results['energy']:.6f} vs {hf_energy:.6f})"
                })
                if not below_hf:
                    validation['passed'] = False
                    logger.warning(f"Correlated method energy ({self.results['energy']:.6f}) above HF ({hf_energy:.6f})!")

        # Check 3: Convergence
        if 'converged' in self.results:
            converged = self.results['converged']
            validation['checks'].append({
                'name': 'convergence',
                'passed': converged,
                'message': f"Solver {'converged' if converged else 'did not converge'}"
            })
            if not converged:
                logger.info("Solver did not fully converge (may need more iterations or different optimizer)")

        return validation

    def print_summary(self):
        """Print human-readable summary of results."""
        print("=" * 80)
        print(f"{self.__class__.__name__} RESULTS")
        print("=" * 80)

        # System info
        print(f"\nSystem: {'-'.join([a.symbol for a in self.atoms])}")
        print(f"Bond Type: {self.bond.bond_type}")
        print(f"Electrons: {self.molecule.n_electrons}")
        print(f"Orbitals: {self.hamiltonian.n_orbitals}")

        # Energy
        if 'energy' in self.results:
            print(f"\nGround State Energy: {self.results['energy']:.8f} Hartree")

        # Convergence
        if 'converged' in self.results:
            status = "✓ Converged" if self.results['converged'] else "✗ Not Converged"
            print(f"Status: {status}")

        if 'iterations' in self.results:
            print(f"Iterations: {self.results['iterations']}")

        # Correlation energy
        if 'correlation_energy' in self.results:
            print(f"Correlation Energy: {self.results['correlation_energy']:.8f} Hartree")

        # Analysis
        if 'analysis' in self.results and self.results['analysis']:
            print("\n" + "-" * 80)
            print("ANALYSIS")
            print("-" * 80)

            if self.results['analysis'].get('bonding'):
                print("\nBonding Analysis:")
                bonding = self.results['analysis']['bonding']
                for key, value in bonding.items():
                    if isinstance(value, (int, float)):
                        print(f"  {key}: {value:.4f}")

        # Validation
        validation = self.validate_results()
        if not validation['passed']:
            print("\n" + "-" * 80)
            print("⚠ VALIDATION WARNINGS")
            print("-" * 80)
            for check in validation['checks']:
                if not check['passed']:
                    print(f"✗ {check['name']}: {check['message']}")

        print("=" * 80)
