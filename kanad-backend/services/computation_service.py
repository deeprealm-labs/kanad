"""
Computation Service - Core integration with Kanad framework.

This service provides the bridge between the FastAPI backend and the
Kanad quantum chemistry framework.
"""

from typing import Dict, Any, List, Optional, Callable
import numpy as np
import logging
from datetime import datetime

# Import Kanad framework components
from kanad.core.atom import Atom
from kanad.core.molecule import Molecule
from kanad.bonds import BondFactory
from kanad.solvers import VQESolver
from kanad.analysis import (
    EnergyAnalyzer,
    BondingAnalyzer,
    PropertyCalculator,
    ThermochemistryCalculator,
    FrequencyCalculator,
    UVVisCalculator,
    ExcitedStateSolver,
    VibronicCalculator,
    UncertaintyAnalyzer,
    BondLengthScanner,
    DOSCalculator
)
from kanad.io import from_smiles, from_xyz, to_xyz

logger = logging.getLogger(__name__)


class ComputationService:
    """
    Service for quantum chemistry computations using Kanad.

    Provides high-level API for molecule creation, computation,
    and analysis.
    """

    def __init__(self):
        """Initialize computation service."""
        logger.info("ComputationService initialized")

    async def create_molecule_from_atoms(
        self,
        atoms: List[Dict[str, Any]],
        basis: str = 'sto-3g',
        charge: int = 0,
        multiplicity: int = 1
    ) -> Dict[str, Any]:
        """
        Create molecule from list of atoms.

        Args:
            atoms: List of dicts with 'element' and 'position' keys
            basis: Basis set name
            charge: Molecular charge
            multiplicity: Spin multiplicity (2S+1)

        Returns:
            Dict with molecule object and metadata
        """
        logger.info(f"Creating molecule from {len(atoms)} atoms, basis={basis}")

        try:
            # Convert to Kanad Atom objects
            kanad_atoms = []
            for atom_data in atoms:
                element = atom_data['element']
                position = np.array(atom_data['position'])
                kanad_atoms.append(Atom(element, position))

            # Create Molecule
            spin = multiplicity - 1  # Convert multiplicity to spin
            molecule = Molecule(
                atoms=kanad_atoms,
                charge=charge,
                spin=spin,
                basis=basis
            )

            # Generate Lewis structure (simplified - would use RDKit in production)
            lewis_structure = self._generate_lewis_structure(kanad_atoms)

            result = {
                "molecule": molecule,
                "formula": molecule.formula,
                "n_electrons": molecule.n_electrons,
                "n_orbitals": molecule.n_orbitals,
                "n_qubits": 2 * molecule.n_orbitals,
                "lewis_structure": lewis_structure,
                "geometry": self._serialize_geometry(kanad_atoms)
            }

            logger.info(f"Molecule created: {molecule.formula}, {molecule.n_electrons} electrons")
            return result

        except Exception as e:
            logger.error(f"Failed to create molecule: {e}")
            raise

    async def create_molecule_from_smiles(
        self,
        smiles: str,
        basis: str = 'sto-3g',
        optimize_geometry: bool = False
    ) -> Dict[str, Any]:
        """
        Create molecule from SMILES string.

        Args:
            smiles: SMILES notation
            basis: Basis set
            optimize_geometry: Whether to optimize geometry

        Returns:
            Dict with molecule and metadata
        """
        logger.info(f"Creating molecule from SMILES: {smiles}")

        try:
            # Parse SMILES using Kanad's parser
            molecule = from_smiles(smiles, basis=basis)

            # Optionally optimize geometry
            if optimize_geometry:
                logger.info("Optimizing geometry...")
                from kanad.optimization import GeometryOptimizer
                optimizer = GeometryOptimizer(molecule)
                optimized_result = optimizer.optimize()
                molecule = optimized_result['molecule']

            result = {
                "molecule": molecule,
                "formula": molecule.formula,
                "smiles": smiles,
                "n_electrons": molecule.n_electrons,
                "n_orbitals": molecule.n_orbitals,
                "n_qubits": 2 * molecule.n_orbitals,
                "geometry": self._serialize_geometry(molecule.atoms)
            }

            return result

        except Exception as e:
            logger.error(f"Failed to create molecule from SMILES: {e}")
            raise

    async def run_computation(
        self,
        molecule: Molecule,
        method: str = 'VQE',
        config: Dict[str, Any] = None,
        progress_callback: Optional[Callable] = None
    ) -> Dict[str, Any]:
        """
        Run quantum chemistry computation.

        Args:
            molecule: Kanad Molecule object
            method: Computation method ('HF', 'VQE', 'MP2', etc.)
            config: Method-specific configuration
            progress_callback: Async callback for progress updates

        Returns:
            Computation results
        """
        config = config or {}
        logger.info(f"Running {method} computation on {molecule.formula}")

        try:
            if method.upper() == 'HF':
                return await self._run_hf(molecule, progress_callback)
            elif method.upper() == 'VQE':
                return await self._run_vqe(molecule, config, progress_callback)
            elif method.upper() == 'MP2':
                return await self._run_mp2(molecule, progress_callback)
            elif method.upper() == 'SQD':
                return await self._run_sqd(molecule, config, progress_callback)
            elif method.upper() == 'EXCITED_STATES':
                return await self._run_excited_states(molecule, config, progress_callback)
            else:
                raise ValueError(f"Unknown method: {method}")

        except Exception as e:
            logger.error(f"Computation failed: {e}")
            raise

    async def _run_hf(
        self,
        molecule: Molecule,
        progress_callback: Optional[Callable] = None
    ) -> Dict[str, Any]:
        """Run Hartree-Fock calculation."""
        if progress_callback:
            await progress_callback(10, "Building Hamiltonian...")

        hamiltonian = molecule.hamiltonian
        hf_energy = hamiltonian.hf_energy

        if progress_callback:
            await progress_callback(100, "HF complete")

        return {
            'method': 'HF',
            'energy': hf_energy,
            'hf_energy': hf_energy,
            'correlation_energy': 0.0,
            'converged': hamiltonian._scf_converged,
            'n_iterations': getattr(hamiltonian, '_scf_iterations', 0)
        }

    async def _run_vqe(
        self,
        molecule: Molecule,
        config: Dict[str, Any],
        progress_callback: Optional[Callable] = None
    ) -> Dict[str, Any]:
        """Run VQE calculation."""
        logger.info("Initializing VQE solver...")

        if progress_callback:
            await progress_callback(5, "Initializing VQE...")

        # Extract config
        ansatz_type = config.get('ansatz', 'hardware_efficient')
        mapper_type = config.get('mapper', 'jordan_wigner')
        optimizer = config.get('optimizer', 'SLSQP')
        max_iterations = config.get('max_iterations', 1000)
        backend_type = config.get('backend', {}).get('type', 'classical')

        # Create bond from molecule (for governance protocols)
        # For multi-atom molecules, use the molecule directly
        bond = molecule  # VQESolver can accept Molecule too

        # Create VQE solver
        from kanad.solvers.vqe_solver import VQESolver

        solver = VQESolver(
            bond=bond,
            ansatz_type=ansatz_type,
            mapper_type=mapper_type,
            optimizer=optimizer,
            max_iterations=max_iterations,
            backend=backend_type
        )

        if progress_callback:
            await progress_callback(10, "VQE initialized, starting optimization...")

        # Define progress callback wrapper
        def vqe_callback(iteration, energy, params=None):
            if progress_callback:
                progress = 10 + int(80 * iteration / max_iterations)
                import asyncio
                asyncio.create_task(
                    progress_callback(progress, f"Iteration {iteration}/{max_iterations}: E={energy:.6f} Ha")
                )

        # Run VQE
        result = solver.solve(callback=vqe_callback)

        if progress_callback:
            await progress_callback(100, "VQE complete!")

        return {
            'method': 'VQE',
            'energy': result['energy'],
            'hf_energy': result.get('hf_energy'),
            'correlation_energy': result.get('correlation_energy'),
            'converged': result.get('converged'),
            'n_iterations': result.get('iterations'),
            'convergence_history': self._format_convergence_history(result.get('energy_history', []))
        }

    async def _run_mp2(
        self,
        molecule: Molecule,
        progress_callback: Optional[Callable] = None
    ) -> Dict[str, Any]:
        """Run MP2 calculation."""
        if progress_callback:
            await progress_callback(10, "Running HF reference...")

        # Get HF energy first
        hamiltonian = molecule.hamiltonian
        hf_energy = hamiltonian.hf_energy

        if progress_callback:
            await progress_callback(50, "Computing MP2 correlation...")

        # MP2 correlation
        from kanad.core.correlation import MP2Solver
        mp2 = MP2Solver(hamiltonian)
        mp2_energy = mp2.compute_energy()

        if progress_callback:
            await progress_callback(100, "MP2 complete")

        return {
            'method': 'MP2',
            'energy': mp2_energy,
            'hf_energy': hf_energy,
            'correlation_energy': mp2_energy - hf_energy,
            'converged': True,
            'n_iterations': 1
        }

    async def run_analysis(
        self,
        molecule: Molecule,
        results: Dict[str, Any],
        analysis_requests: Dict[str, bool]
    ) -> Dict[str, Any]:
        """
        Run requested analyses on computation results.

        Args:
            molecule: Molecule object
            results: Computation results
            analysis_requests: Dict of analysis flags

        Returns:
            Analysis results
        """
        logger.info("Running analyses...")
        analysis = {}

        try:
            hamiltonian = molecule.hamiltonian

            # Energy decomposition
            if analysis_requests.get('energy_decomposition'):
                logger.info("Running energy decomposition...")
                analyzer = EnergyAnalyzer(hamiltonian)
                decomp = analyzer.decompose_energy()
                analysis['energy_decomposition'] = decomp

            # Bond analysis
            if analysis_requests.get('bond_analysis'):
                logger.info("Running bond analysis...")
                analyzer = BondingAnalyzer(hamiltonian)
                bonds = {
                    'bond_orders': analyzer.compute_bond_order(),
                    'homo_lumo_gap': analyzer.compute_homo_lumo_gap()
                }
                analysis['bond_analysis'] = bonds

            # Molecular properties
            if analysis_requests.get('dipole_moment'):
                logger.info("Computing dipole moment...")
                calculator = PropertyCalculator(hamiltonian)
                dipole = calculator.compute_dipole_moment()
                analysis['dipole_moment'] = {
                    'magnitude': np.linalg.norm(dipole),
                    'direction': dipole.tolist()
                }

            # Thermochemistry
            if analysis_requests.get('thermochemistry'):
                logger.info("Computing thermochemistry...")
                thermo = ThermochemistryCalculator(molecule, results['energy'])
                thermo_data = thermo.compute_all()
                analysis['thermochemistry'] = thermo_data

            # Spectroscopy
            if analysis_requests.get('spectroscopy'):
                logger.info("Computing UV-Vis spectrum...")
                spectro = UVVisCalculator(molecule)
                spectrum = spectro.compute_spectrum()
                analysis['spectroscopy'] = spectrum

            return analysis

        except Exception as e:
            logger.error(f"Analysis failed: {e}")
            return analysis  # Return partial results

    def _generate_lewis_structure(self, atoms: List[Atom]) -> Dict[str, Any]:
        """
        Generate Lewis structure (simplified version).

        In production, this would use RDKit or similar.
        """
        # Simple bond detection based on distances
        bonds = []
        for i in range(len(atoms)):
            for j in range(i + 1, len(atoms)):
                dist = np.linalg.norm(atoms[i].position - atoms[j].position)
                # Simple heuristic: bond if distance < sum of covalent radii * 1.2
                expected_dist = (atoms[i].covalent_radius + atoms[j].covalent_radius) * 1.2
                if dist < expected_dist:
                    bonds.append([i, j])

        return {
            'bonds': bonds,
            'lone_pairs': {},  # Would calculate from valence electrons
            'formal_charges': {}
        }

    def _serialize_geometry(self, atoms: List[Atom]) -> Dict[str, Any]:
        """Serialize molecular geometry to dict."""
        return {
            'atoms': [
                {
                    'element': atom.symbol,
                    'position': atom.position.tolist()
                }
                for atom in atoms
            ]
        }

    def _format_convergence_history(self, energy_history: np.ndarray) -> List[Dict]:
        """Format convergence history for API response."""
        if len(energy_history) == 0:
            return []

        return [
            {'iteration': i, 'energy': float(e)}
            for i, e in enumerate(energy_history)
        ]

    # ===== MISSING SOLVERS =====

    async def _run_sqd(
        self,
        molecule: Molecule,
        config: Dict[str, Any],
        progress_callback: Optional[Callable] = None
    ) -> Dict[str, Any]:
        """
        Run Subspace Quantum Diagonalization (SQD).

        SQD solver for efficient ground state calculations.
        """
        logger.info("Initializing SQD solver...")

        if progress_callback:
            await progress_callback(5, "Initializing SQD...")

        try:
            from kanad.solvers.sqd_solver import SQDSolver

            # Extract config
            active_space = config.get('advanced', {}).get('active_space')
            max_iterations = config.get('max_iterations', 1000)

            # Create SQD solver
            solver = SQDSolver(molecule)

            # Apply active space if specified
            if active_space:
                from kanad.solvers.active_space import ActiveSpaceSelector
                selector = ActiveSpaceSelector(molecule)
                active_orbitals = selector.select_cas(
                    n_electrons=active_space['n_electrons'],
                    n_orbitals=active_space['n_orbitals']
                )
                solver.set_active_space(active_orbitals)

            if progress_callback:
                await progress_callback(10, "Running SQD calculation...")

            # Run SQD
            result = solver.solve()

            if progress_callback:
                await progress_callback(100, "SQD complete!")

            return {
                'method': 'SQD',
                'energy': result['energy'],
                'hf_energy': result.get('hf_energy'),
                'correlation_energy': result.get('correlation_energy', result['energy'] - result.get('hf_energy', 0)),
                'converged': True,
                'n_iterations': result.get('iterations', 1),
                'subspace_dimension': result.get('subspace_dim', 0)
            }

        except Exception as e:
            logger.error(f"SQD solver failed: {e}")
            raise

    async def _run_excited_states(
        self,
        molecule: Molecule,
        config: Dict[str, Any],
        progress_callback: Optional[Callable] = None
    ) -> Dict[str, Any]:
        """
        Run excited states calculation.

        Computes multiple excited electronic states.
        """
        logger.info("Initializing excited states solver...")

        if progress_callback:
            await progress_callback(5, "Initializing excited states calculation...")

        try:
            from kanad.solvers.excited_states_solver import ExcitedStatesSolver

            # Get number of excited states to compute
            n_states = config.get('n_excited_states', 5)

            # Create solver
            solver = ExcitedStatesSolver(molecule)

            if progress_callback:
                await progress_callback(10, f"Computing {n_states} excited states...")

            # Compute excited states
            result = solver.solve(n_states=n_states)

            if progress_callback:
                await progress_callback(100, "Excited states calculation complete!")

            # Format results
            excited_energies = result.get('excited_energies', [])
            oscillator_strengths = result.get('oscillator_strengths', [])

            return {
                'method': 'EXCITED_STATES',
                'energy': result['ground_state_energy'],
                'hf_energy': result.get('hf_energy'),
                'correlation_energy': 0.0,
                'converged': True,
                'n_iterations': 1,
                'excited_states': [
                    {
                        'state_number': i + 1,
                        'energy': float(e),
                        'excitation_energy': float(e - result['ground_state_energy']),
                        'oscillator_strength': float(oscillator_strengths[i]) if i < len(oscillator_strengths) else 0.0
                    }
                    for i, e in enumerate(excited_energies)
                ]
            }

        except Exception as e:
            logger.error(f"Excited states solver failed: {e}")
            raise

    # ===== ENHANCED ANALYSIS =====

    async def run_enhanced_analysis(
        self,
        molecule: Molecule,
        results: Dict[str, Any],
        analysis_requests: Dict[str, bool]
    ) -> Dict[str, Any]:
        """
        Run enhanced analysis with all available tools.

        Includes spectroscopy, vibrational analysis, uncertainty, bond scanning, and DOS.
        """
        logger.info("Running enhanced analyses...")
        analysis = await self.run_analysis(molecule, results, analysis_requests)

        try:
            hamiltonian = molecule.hamiltonian

            # Spectroscopy analysis (UV-Vis, excited states)
            if analysis_requests.get('spectroscopy'):
                logger.info("Running spectroscopy analysis...")
                try:
                    uv_vis_calc = UVVisCalculator(molecule)
                    uv_spectrum = uv_vis_calc.compute_spectrum()

                    excited_solver = ExcitedStateSolver(molecule)
                    excited_states = excited_solver.solve(n_states=5)

                    analysis['spectroscopy'] = {
                        'uv_vis_spectrum': uv_spectrum,
                        'excited_states': excited_states
                    }
                except Exception as e:
                    logger.warning(f"Spectroscopy analysis failed: {e}")

            # Vibrational analysis
            if analysis_requests.get('vibrational'):
                logger.info("Running vibrational analysis...")
                try:
                    freq_calc = FrequencyCalculator(molecule)
                    frequencies = freq_calc.compute_frequencies()
                    analysis['vibrational'] = {
                        'frequencies': frequencies.get('frequencies', []),
                        'normal_modes': frequencies.get('normal_modes', []),
                        'zero_point_energy': frequencies.get('zpe', 0.0)
                    }
                except Exception as e:
                    logger.warning(f"Vibrational analysis failed: {e}")

            # Uncertainty analysis
            if analysis_requests.get('uncertainty'):
                logger.info("Computing uncertainty estimates...")
                try:
                    uncertainty_analyzer = UncertaintyAnalyzer(molecule, results)
                    uncertainty = uncertainty_analyzer.estimate_errors()
                    analysis['uncertainty'] = {
                        'energy_error': uncertainty.get('energy_error', 0.0),
                        'statistical_error': uncertainty.get('statistical_error', 0.0),
                        'systematic_error': uncertainty.get('systematic_error', 0.0)
                    }
                except Exception as e:
                    logger.warning(f"Uncertainty analysis failed: {e}")

            # Bond scanning (potential energy surface scan)
            if analysis_requests.get('bond_scan'):
                logger.info("Running bond scan...")
                try:
                    bond_length_scanner = BondLengthScanner(molecule)
                    scan_data = bond_length_scanner.scan()
                    analysis['bond_scan'] = {
                        'scanned_bonds': scan_data.get('bonds', []),
                        'energies': scan_data.get('energies', []),
                        'distances': scan_data.get('distances', [])
                    }
                except Exception as e:
                    logger.warning(f"Bond scan failed: {e}")

            # Density of States (for metallurgy/periodic systems)
            if analysis_requests.get('dos'):
                logger.info("Computing density of states...")
                try:
                    dos_calculator = DOSCalculator(hamiltonian)
                    dos_data = dos_calculator.compute_dos()
                    analysis['dos'] = {
                        'energies': dos_data.get('energies', []),
                        'dos_values': dos_data.get('dos', []),
                        'fermi_energy': dos_data.get('fermi_energy', 0.0)
                    }
                except Exception as e:
                    logger.warning(f"DOS calculation failed: {e}")

            return analysis

        except Exception as e:
            logger.error(f"Enhanced analysis failed: {e}")
            return analysis  # Return partial results

    # ===== OPTIMIZATION TOOLS =====

    async def run_geometry_optimization(
        self,
        molecule: Molecule,
        progress_callback: Optional[Callable] = None
    ) -> Dict[str, Any]:
        """Run geometry optimization."""
        logger.info("Starting geometry optimization...")

        if progress_callback:
            await progress_callback(5, "Initializing geometry optimizer...")

        try:
            from kanad.optimization import GeometryOptimizer

            optimizer = GeometryOptimizer(molecule)

            if progress_callback:
                await progress_callback(10, "Optimizing geometry...")

            result = optimizer.optimize()

            if progress_callback:
                await progress_callback(100, "Geometry optimization complete!")

            return {
                'optimized_molecule': result['molecule'],
                'final_energy': result['energy'],
                'optimized_geometry': self._serialize_geometry(result['molecule'].atoms),
                'n_steps': result.get('n_steps', 0),
                'converged': result.get('converged', True)
            }

        except Exception as e:
            logger.error(f"Geometry optimization failed: {e}")
            raise

    async def run_orbital_optimization(
        self,
        molecule: Molecule,
        progress_callback: Optional[Callable] = None
    ) -> Dict[str, Any]:
        """Run orbital optimization."""
        logger.info("Starting orbital optimization...")

        try:
            from kanad.optimization.orbital_optimizer import OrbitalOptimizer

            optimizer = OrbitalOptimizer(molecule.hamiltonian)
            result = optimizer.optimize()

            return {
                'optimized_orbitals': result['orbitals'],
                'energy_improvement': result.get('energy_improvement', 0.0)
            }

        except Exception as e:
            logger.error(f"Orbital optimization failed: {e}")
            raise

    async def run_circuit_optimization(
        self,
        circuit: Any,
        progress_callback: Optional[Callable] = None
    ) -> Dict[str, Any]:
        """Run quantum circuit optimization."""
        logger.info("Starting circuit optimization...")

        try:
            from kanad.optimization.circuit_optimizer import CircuitOptimizer

            optimizer = CircuitOptimizer(circuit)
            result = optimizer.optimize()

            return {
                'optimized_circuit': result['circuit'],
                'depth_reduction': result.get('depth_reduction', 0),
                'gate_count_reduction': result.get('gate_count_reduction', 0)
            }

        except Exception as e:
            logger.error(f"Circuit optimization failed: {e}")
            raise

    async def run_adaptive_vqe(
        self,
        molecule: Molecule,
        config: Dict[str, Any],
        progress_callback: Optional[Callable] = None
    ) -> Dict[str, Any]:
        """Run Adaptive VQE."""
        logger.info("Starting Adaptive VQE...")

        try:
            from kanad.optimization.adaptive_optimizer import AdaptiveOptimizer

            optimizer = AdaptiveOptimizer(molecule)
            result = optimizer.run_adaptive_vqe()

            return {
                'energy': result['energy'],
                'circuit': result['circuit'],
                'n_parameters': result.get('n_parameters', 0)
            }

        except Exception as e:
            logger.error(f"Adaptive VQE failed: {e}")
            raise

    # ===== PERIODIC/CRYSTAL SYSTEMS =====

    async def create_periodic_system(
        self,
        lattice_vectors: List[List[float]],
        basis_atoms: List[Dict[str, Any]],
        basis: str = 'sto-3g'
    ) -> Dict[str, Any]:
        """
        Create periodic system (crystal/material).

        Args:
            lattice_vectors: 3x3 lattice vectors
            basis_atoms: Atoms in unit cell
            basis: Basis set

        Returns:
            Periodic system data
        """
        logger.info("Creating periodic system...")

        try:
            from kanad.core.lattice import Lattice
            from kanad.core.hamiltonians.periodic_hamiltonian import PeriodicHamiltonian
            from kanad.io.crystal_builder import CrystalBuilder

            # Create lattice
            lattice = Lattice(np.array(lattice_vectors))

            # Create crystal structure
            crystal = CrystalBuilder.from_lattice_and_atoms(
                lattice,
                basis_atoms,
                basis=basis
            )

            # Create periodic Hamiltonian
            hamiltonian = PeriodicHamiltonian(crystal)

            result = {
                'crystal': crystal,
                'hamiltonian': hamiltonian,
                'lattice_parameters': {
                    'a': lattice.a,
                    'b': lattice.b,
                    'c': lattice.c,
                    'alpha': lattice.alpha,
                    'beta': lattice.beta,
                    'gamma': lattice.gamma
                },
                'space_group': crystal.get_space_group(),
                'n_atoms': crystal.n_atoms
            }

            logger.info(f"Periodic system created: {crystal.n_atoms} atoms")
            return result

        except Exception as e:
            logger.error(f"Failed to create periodic system: {e}")
            raise

    async def analyze_band_structure(
        self,
        periodic_hamiltonian: Any,
        k_path: List[List[float]] = None
    ) -> Dict[str, Any]:
        """
        Compute band structure for periodic system.

        Args:
            periodic_hamiltonian: PeriodicHamiltonian object
            k_path: k-space path for band structure

        Returns:
            Band structure data
        """
        logger.info("Computing band structure...")

        try:
            # Compute band structure along high-symmetry path
            bands = periodic_hamiltonian.compute_band_structure(k_path=k_path)

            # Compute Fermi energy
            fermi_energy = periodic_hamiltonian.compute_fermi_energy()

            # Compute band gap
            band_gap = periodic_hamiltonian.compute_band_gap()

            result = {
                'bands': bands,
                'fermi_energy': fermi_energy,
                'band_gap': band_gap,
                'is_metal': band_gap < 0.01,  # Small gap ~ metallic
                'is_semiconductor': 0.01 <= band_gap <= 3.0,
                'is_insulator': band_gap > 3.0
            }

            logger.info(f"Band structure computed: Eg = {band_gap:.3f} eV")
            return result

        except Exception as e:
            logger.error(f"Band structure calculation failed: {e}")
            raise
