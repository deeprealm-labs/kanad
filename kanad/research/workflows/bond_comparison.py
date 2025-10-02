"""
Bond comparison workflow for research.

Compares different bonding types for the same atom pair.
"""

from typing import Dict, Any, List
import numpy as np

from kanad.bonds import BondFactory, BondType


class BondComparisonWorkflow:
    """
    Compare different bonding types for research purposes.

    Helps researchers understand which bonding model
    is most appropriate for a given atom pair.
    """

    @staticmethod
    def compare_bond_types(
        atom_1: str,
        atom_2: str,
        methods: List[str] = ['HF'],
        bond_types: List[str] = None
    ) -> Dict[str, Any]:
        """
        Compare ionic vs covalent vs metallic for same atom pair.

        Args:
            atom_1: First atom symbol
            atom_2: Second atom symbol
            methods: Computational methods to use (default: ['HF'])
            bond_types: Bond types to try (default: all three)

        Returns:
            Dictionary with comparison results

        Example:
            >>> results = BondComparisonWorkflow.compare_bond_types('Na', 'Cl')
            >>> print(results['best_bond_type'])  # 'ionic'
            >>> print(results['energy_differences'])
        """
        if bond_types is None:
            bond_types = ['ionic', 'covalent']  # Metallic requires multiple atoms

        results = {}

        for bond_type in bond_types:
            try:
                # Create bond
                bond = BondFactory.create_bond(
                    atom_1,
                    atom_2,
                    bond_type=bond_type
                )

                # Compute energy with first method
                method = methods[0] if methods else 'HF'
                energy_data = bond.compute_energy(method=method)

                results[bond_type] = {
                    'success': True,
                    'energy': energy_data['energy'],
                    'method': energy_data['method'],
                    'converged': energy_data['converged'],
                    'bond_analysis': energy_data['bond_analysis'],
                    'bond_object': bond
                }

            except Exception as e:
                results[bond_type] = {
                    'success': False,
                    'error': str(e)
                }

        # Determine most stable (lowest energy)
        successful_results = {
            k: v for k, v in results.items()
            if v.get('success', False)
        }

        if successful_results:
            energies = {k: v['energy'] for k, v in successful_results.items()}
            best = min(energies, key=energies.get)

            results['best_bond_type'] = best
            results['energy_differences'] = {
                k: (v - energies[best]) * 27.211  # Convert to eV
                for k, v in energies.items()
            }
        else:
            results['best_bond_type'] = None
            results['energy_differences'] = {}

        # Add quick info
        results['quick_info'] = BondFactory.quick_bond_info(atom_1, atom_2)

        return results

    @staticmethod
    def analyze_bond_character(
        atom_1: str,
        atom_2: str
    ) -> Dict[str, Any]:
        """
        Analyze bonding character without full calculation.

        Fast analysis based on atomic properties only.

        Args:
            atom_1: First atom symbol
            atom_2: Second atom symbol

        Returns:
            Dictionary with bonding character analysis

        Example:
            >>> analysis = BondComparisonWorkflow.analyze_bond_character('H', 'Cl')
            >>> print(analysis['predicted_type'])  # 'covalent'
            >>> print(analysis['ionic_percentage'])  # 18%
        """
        info = BondFactory.quick_bond_info(atom_1, atom_2)

        # Calculate bonding percentages
        delta_en = info['electronegativity_difference']

        # Empirical formulas
        ionic_percentage = (1.0 - np.exp(-0.25 * delta_en**2)) * 100
        covalent_percentage = 100 - ionic_percentage

        return {
            'atom_1': atom_1,
            'atom_2': atom_2,
            'predicted_type': info['predicted_type'],
            'electronegativity_difference': delta_en,
            'ionic_percentage': ionic_percentage,
            'covalent_percentage': covalent_percentage,
            'bond_polarity': 'polar' if delta_en > 0.4 else 'nonpolar',
            'estimated_bond_length': info['estimated_bond_length'],
            'rationale': info['rationale']
        }

    @staticmethod
    def scan_bond_lengths(
        atom_1: str,
        atom_2: str,
        bond_type: str = 'auto',
        distances: List[float] = None,
        method: str = 'HF'
    ) -> Dict[str, Any]:
        """
        Perform potential energy surface scan.

        Computes energy at different bond lengths to find
        equilibrium geometry and dissociation curve.

        Args:
            atom_1: First atom symbol
            atom_2: Second atom symbol
            bond_type: Bond type (default: auto-detect)
            distances: Bond lengths to scan (Angstroms)
            method: Computational method

        Returns:
            Dictionary with scan results

        Example:
            >>> scan = BondComparisonWorkflow.scan_bond_lengths('H', 'H')
            >>> print(scan['equilibrium_distance'])  # ~0.74 Å
            >>> print(scan['dissociation_energy'])  # eV
        """
        # Default distance range
        if distances is None:
            # Estimate equilibrium distance
            info = BondFactory.quick_bond_info(atom_1, atom_2)
            r_eq = info['estimated_bond_length']

            # Scan from 0.5 * r_eq to 3.0 * r_eq
            distances = np.linspace(0.5 * r_eq, 3.0 * r_eq, 15)

        energies = []
        converged_flags = []

        for distance in distances:
            try:
                bond = BondFactory.create_bond(
                    atom_1,
                    atom_2,
                    bond_type=bond_type,
                    distance=distance
                )

                result = bond.compute_energy(method=method, max_iterations=50)

                energies.append(result['energy'])
                converged_flags.append(result['converged'])

            except Exception as e:
                energies.append(np.nan)
                converged_flags.append(False)

        energies = np.array(energies)

        # Find equilibrium (minimum energy)
        valid_indices = ~np.isnan(energies)
        if np.any(valid_indices):
            min_idx = np.nanargmin(energies)
            equilibrium_distance = distances[min_idx]
            equilibrium_energy = energies[min_idx]

            # Estimate dissociation energy (E_∞ - E_min)
            # Approximate E_∞ as energy at largest distance
            e_infinity = energies[valid_indices][-1]
            dissociation_energy = (e_infinity - equilibrium_energy) * 27.211  # eV

        else:
            equilibrium_distance = None
            equilibrium_energy = None
            dissociation_energy = None

        return {
            'atom_1': atom_1,
            'atom_2': atom_2,
            'bond_type': bond_type,
            'distances': distances,
            'energies': energies,
            'converged': converged_flags,
            'equilibrium_distance': equilibrium_distance,
            'equilibrium_energy': equilibrium_energy,
            'dissociation_energy': dissociation_energy,
            'method': method
        }
