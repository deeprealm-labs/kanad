#!/usr/bin/env python3
"""
Scientific Test 1: 1D Chain Physics and Edge Cases

This test rigorously validates the 1D tight-binding chain implementation
against known theoretical results and explores edge cases.

Scientific Focus:
- Tight-binding theory validation
- Periodic vs non-periodic boundary conditions
- Finite-size effects
- Band structure accuracy
- Metallic character detection
- Edge cases and numerical stability
"""

import numpy as np
import sys
import traceback
from typing import Dict, Any, List, Tuple
import matplotlib.pyplot as plt

# Add the project root to the path
sys.path.insert(0, '/home/mk/deeprealm/kanad')

from kanad.core.atom import Atom
from kanad.bonds.metallic_bond import MetallicBond
from kanad.core.hamiltonians.metallic_hamiltonian import MetallicHamiltonian

class Test1DChainPhysics:
    """Comprehensive test of 1D chain physics."""
    
    def __init__(self):
        self.test_results = {}
        self.critical_issues = []
        self.warnings = []
        self.passed_tests = 0
        self.failed_tests = 0
        
    def test_tight_binding_matrix_construction(self):
        """Test 1.1: Tight-binding matrix construction accuracy."""
        print("="*70)
        print("TEST 1.1: TIGHT-BINDING MATRIX CONSTRUCTION")
        print("="*70)
        
        issues = []
        
        # Test different system sizes
        for n in [2, 3, 4, 5, 6, 7, 8, 10, 16]:
            print(f"\nTesting N={n} atoms:")
            
            atoms = [Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(n)]
            bond = MetallicBond(atoms, hopping_parameter=-1.0, periodic=True)
            H = bond.hamiltonian.h_tight_binding
            
            # Check matrix properties
            print(f"  Matrix shape: {H.shape}")
            print(f"  Matrix is symmetric: {np.allclose(H, H.T)}")
            print(f"  Matrix is real: {np.allclose(H, H.real)}")
            print(f"  Diagonal elements: {np.diag(H)}")
            print(f"  Off-diagonal elements: {H[np.triu_indices_from(H, k=1)]}")
            
            # Check specific properties
            if not np.allclose(H, H.T):
                issues.append(f"N={n}: Matrix not symmetric")
            
            if not np.allclose(H, H.real):
                issues.append(f"N={n}: Matrix has imaginary parts")
            
            # Check periodic boundary conditions
            if n > 2:
                if abs(H[0, n-1] - (-1.0)) > 1e-10:
                    issues.append(f"N={n}: PBC not properly implemented (H[0,{n-1}]={H[0, n-1]})")
                if abs(H[n-1, 0] - (-1.0)) > 1e-10:
                    issues.append(f"N={n}: PBC not properly implemented (H[{n-1},0]={H[n-1, 0]})")
            
            # Check nearest-neighbor hopping
            for i in range(n-1):
                if abs(H[i, i+1] - (-1.0)) > 1e-10:
                    issues.append(f"N={n}: Nearest-neighbor hopping incorrect (H[{i},{i+1}]={H[i, i+1]})")
                if abs(H[i+1, i] - (-1.0)) > 1e-10:
                    issues.append(f"N={n}: Nearest-neighbor hopping incorrect (H[{i+1},{i}]={H[i+1, i]})")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ All tight-binding matrices are correct")
            self.passed_tests += 1
        
        return issues
    
    def test_eigenvalue_accuracy(self):
        """Test 1.2: Eigenvalue accuracy against analytical results."""
        print("="*70)
        print("TEST 1.2: EIGENVALUE ACCURACY")
        print("="*70)
        
        issues = []
        
        # Test known analytical results
        test_cases = [
            (2, [-1.0, 1.0]),  # N=2: eigenvalues are ±t
            (3, [-2.0, 1.0, 1.0]),  # N=3: -2t, t, t
            (4, [-2.0, 0.0, 0.0, 2.0]),  # N=4: -2t, 0, 0, 2t
        ]
        
        for n, expected_eigenvalues in test_cases:
            print(f"\nTesting N={n} (expected: {expected_eigenvalues}):")
            
            atoms = [Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(n)]
            bond = MetallicBond(atoms, hopping_parameter=-1.0, periodic=True)
            eigenvalues = np.linalg.eigvalsh(bond.hamiltonian.h_tight_binding)
            
            print(f"  Computed: {eigenvalues}")
            print(f"  Expected: {expected_eigenvalues}")
            
            # Check if eigenvalues match (allowing for numerical precision)
            if len(eigenvalues) != len(expected_eigenvalues):
                issues.append(f"N={n}: Wrong number of eigenvalues ({len(eigenvalues)} vs {len(expected_eigenvalues)})")
            else:
                for i, (computed, expected) in enumerate(zip(eigenvalues, expected_eigenvalues)):
                    if abs(computed - expected) > 1e-10:
                        issues.append(f"N={n}: Eigenvalue {i} incorrect ({computed} vs {expected})")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ All eigenvalues are correct")
            self.passed_tests += 1
        
        return issues
    
    def test_band_structure_accuracy(self):
        """Test 1.3: Band structure accuracy against analytical dispersion."""
        print("="*70)
        print("TEST 1.3: BAND STRUCTURE ACCURACY")
        print("="*70)
        
        issues = []
        
        # Test 1D chain band structure
        atoms = [Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(4)]
        bond = MetallicBond(atoms, hopping_parameter=-1.0, periodic=True)
        
        band_structure = bond.get_band_structure()
        k_points = band_structure['k_points']
        energies = band_structure['energies']
        
        print(f"k-points: {len(k_points)} points from {k_points[0]:.3f} to {k_points[-1]:.3f}")
        print(f"Energy shape: {energies.shape}")
        
        # Test dispersion relation E(k) = -2t cos(k)
        t = -1.0
        max_error = 0.0
        
        for i, k in enumerate(k_points):
            computed_E = energies[i, 0]
            expected_E = 2 * t * np.cos(k)
            error = abs(computed_E - expected_E)
            max_error = max(max_error, error)
            
            if error > 1e-10:
                issues.append(f"k={k:.6f}: E={computed_E:.6f}, expected={expected_E:.6f}, error={error:.2e}")
        
        print(f"Maximum dispersion error: {max_error:.2e}")
        
        # Test specific k-points
        k_zero_idx = np.argmin(np.abs(k_points))
        k_pi_idx = np.argmin(np.abs(k_points - np.pi))
        
        E_k0 = energies[k_zero_idx, 0]
        E_kpi = energies[k_pi_idx, 0]
        
        expected_E_k0 = 2 * t * np.cos(k_points[k_zero_idx])
        expected_E_kpi = 2 * t * np.cos(k_points[k_pi_idx])
        
        print(f"E(k≈0): {E_k0:.6f} eV, expected: {expected_E_k0:.6f} eV")
        print(f"E(k≈π): {E_kpi:.6f} eV, expected: {expected_E_kpi:.6f} eV")
        
        if abs(E_k0 - expected_E_k0) > 1e-10:
            issues.append(f"E(k≈0) incorrect: {E_k0} vs {expected_E_k0}")
        if abs(E_kpi - expected_E_kpi) > 1e-10:
            issues.append(f"E(k≈π) incorrect: {E_kpi} vs {expected_E_kpi}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues[:5]:  # Show first 5 issues
                print(f"  - {issue}")
            if len(issues) > 5:
                print(f"  ... and {len(issues) - 5} more issues")
            self.failed_tests += 1
        else:
            print(f"\n✅ Band structure is accurate")
            self.passed_tests += 1
        
        return issues
    
    def test_metallic_character_physics(self):
        """Test 1.4: Metallic character detection physics."""
        print("="*70)
        print("TEST 1.4: METALLIC CHARACTER DETECTION PHYSICS")
        print("="*70)
        
        issues = []
        
        # Test different system sizes
        test_sizes = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
        
        for n in test_sizes:
            atoms = [Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(n)]
            bond = MetallicBond(atoms, hopping_parameter=-1.0, periodic=True)
            
            # Get eigenvalues and check metallic character
            eigenvalues = np.linalg.eigvalsh(bond.hamiltonian.h_tight_binding)
            is_metallic = bond.hamiltonian.is_metallic()
            fermi_energy = bond.hamiltonian.get_fermi_energy(eigenvalues)
            dos_fermi = bond.hamiltonian.compute_dos_at_fermi(eigenvalues)
            
            # Expected metallic character based on physics
            # For 1D chain with 1 electron per site:
            # - Odd N: metallic (states at Fermi level)
            # - Even N: depends on degeneracies
            
            # Check if there are states at Fermi level
            states_at_fermi = np.abs(eigenvalues - fermi_energy) < 1e-12
            expected_metallic = np.sum(states_at_fermi) > 0
            
            if is_metallic != expected_metallic:
                issues.append(f"N={n}: Metallic={is_metallic}, expected={expected_metallic}, "
                            f"DOS={dos_fermi}, states_at_fermi={np.sum(states_at_fermi)}")
            
            # Check Fermi energy calculation
            n_electrons = n
            n_bands_occupied = int(np.ceil(n_electrons / 2.0))
            
            if n_electrons % 2 == 0 and n_bands_occupied < len(eigenvalues):
                expected_fermi = (eigenvalues[n_bands_occupied - 1] + eigenvalues[n_bands_occupied]) / 2.0
            else:
                expected_fermi = eigenvalues[n_bands_occupied - 1]
            
            if abs(fermi_energy - expected_fermi) > 1e-10:
                issues.append(f"N={n}: Fermi energy incorrect: {fermi_energy} vs {expected_fermi}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues[:10]:  # Show first 10 issues
                print(f"  - {issue}")
            if len(issues) > 10:
                print(f"  ... and {len(issues) - 10} more issues")
            self.failed_tests += 1
        else:
            print(f"\n✅ Metallic character detection is correct")
            self.passed_tests += 1
        
        return issues
    
    def test_periodic_vs_non_periodic(self):
        """Test 1.5: Periodic vs non-periodic boundary conditions."""
        print("="*70)
        print("TEST 1.5: PERIODIC VS NON-PERIODIC BOUNDARY CONDITIONS")
        print("="*70)
        
        issues = []
        
        # Test different system sizes
        for n in [3, 4, 5, 6, 8, 10]:
            print(f"\nTesting N={n}:")
            
            atoms = [Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(n)]
            
            # Periodic boundary conditions
            bond_pbc = MetallicBond(atoms, hopping_parameter=-1.0, periodic=True)
            H_pbc = bond_pbc.hamiltonian.h_tight_binding
            eigenvalues_pbc = np.linalg.eigvalsh(H_pbc)
            energy_pbc = bond_pbc.compute_energy(method='tight_binding')['energy']
            
            # Non-periodic boundary conditions
            bond_npbc = MetallicBond(atoms, hopping_parameter=-1.0, periodic=False)
            H_npbc = bond_npbc.hamiltonian.h_tight_binding
            eigenvalues_npbc = np.linalg.eigvalsh(H_npbc)
            energy_npbc = bond_npbc.compute_energy(method='tight_binding')['energy']
            
            print(f"  PBC energy: {energy_pbc:.6f} eV")
            print(f"  NPBC energy: {energy_npbc:.6f} eV")
            print(f"  Energy difference: {abs(energy_pbc - energy_npbc):.6f} eV")
            
            # Check that PBC and NPBC give different results
            if abs(energy_pbc - energy_npbc) < 1e-10:
                issues.append(f"N={n}: PBC and NPBC give identical energies")
            
            # Check that PBC has more hopping terms
            pbc_hopping_terms = np.sum(np.abs(H_pbc) > 1e-10) - n  # Subtract diagonal
            npbc_hopping_terms = np.sum(np.abs(H_npbc) > 1e-10) - n  # Subtract diagonal
            
            if pbc_hopping_terms <= npbc_hopping_terms:
                issues.append(f"N={n}: PBC should have more hopping terms than NPBC")
            
            print(f"  PBC hopping terms: {pbc_hopping_terms}")
            print(f"  NPBC hopping terms: {npbc_hopping_terms}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ PBC vs NPBC behavior is correct")
            self.passed_tests += 1
        
        return issues
    
    def test_finite_size_effects(self):
        """Test 1.6: Finite-size effects and scaling."""
        print("="*70)
        print("TEST 1.6: FINITE-SIZE EFFECTS AND SCALING")
        print("="*70)
        
        issues = []
        
        # Test energy per atom scaling
        sizes = [2, 4, 6, 8, 10, 12, 16, 20, 24, 32]
        energies_per_atom = []
        
        for n in sizes:
            atoms = [Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(n)]
            bond = MetallicBond(atoms, hopping_parameter=-1.0, periodic=True)
            result = bond.compute_energy(method='tight_binding')
            energy_per_atom = result['energy'] / n
            energies_per_atom.append(energy_per_atom)
            
            print(f"N={n:2d}: E/N={energy_per_atom:8.4f} eV/atom")
        
        # Check if energy per atom converges to expected value
        # For large N, E/N should approach -2|t| = -2 eV
        expected_large_N = -2.0  # eV/atom
        
        # Use the last few points to estimate large-N limit
        last_energies = energies_per_atom[-3:]
        avg_large_N = np.mean(last_energies)
        
        print(f"\nLarge-N limit estimate: {avg_large_N:.4f} eV/atom")
        print(f"Expected large-N limit: {expected_large_N:.4f} eV/atom")
        print(f"Error: {abs(avg_large_N - expected_large_N):.4f} eV/atom")
        
        if abs(avg_large_N - expected_large_N) > 0.1:  # 0.1 eV tolerance
            issues.append(f"Large-N limit incorrect: {avg_large_N:.4f} vs {expected_large_N:.4f}")
        
        # Check if energy per atom is monotonically decreasing (more negative)
        for i in range(1, len(energies_per_atom)):
            if energies_per_atom[i] > energies_per_atom[i-1]:
                issues.append(f"Energy per atom not monotonically decreasing: "
                            f"N={sizes[i-1]}->{sizes[i]}: {energies_per_atom[i-1]:.4f}->{energies_per_atom[i]:.4f}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Finite-size effects are correct")
            self.passed_tests += 1
        
        return issues
    
    def test_numerical_stability(self):
        """Test 1.7: Numerical stability and edge cases."""
        print("="*70)
        print("TEST 1.7: NUMERICAL STABILITY AND EDGE CASES")
        print("="*70)
        
        issues = []
        
        # Test extreme hopping parameters
        extreme_hopping = [1e-10, 1e-5, 1e-2, 1e2, 1e5, 1e10]
        
        for t in extreme_hopping:
            print(f"\nTesting hopping parameter t={t}:")
            
            try:
                atoms = [Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(4)]
                bond = MetallicBond(atoms, hopping_parameter=t, periodic=True)
                
                # Check if matrix is well-conditioned
                H = bond.hamiltonian.h_tight_binding
                cond_num = np.linalg.cond(H)
                print(f"  Condition number: {cond_num:.2e}")
                
                if cond_num > 1e12:
                    issues.append(f"t={t}: Matrix ill-conditioned (cond={cond_num:.2e})")
                
                # Check eigenvalues
                eigenvalues = np.linalg.eigvalsh(H)
                if np.any(np.isnan(eigenvalues)) or np.any(np.isinf(eigenvalues)):
                    issues.append(f"t={t}: Invalid eigenvalues (NaN or Inf)")
                
                # Check energy computation
                result = bond.compute_energy(method='tight_binding')
                if np.isnan(result['energy']) or np.isinf(result['energy']):
                    issues.append(f"t={t}: Invalid energy (NaN or Inf)")
                
            except Exception as e:
                issues.append(f"t={t}: Exception during computation: {e}")
        
        # Test very large systems
        large_sizes = [50, 100, 200]
        
        for n in large_sizes:
            print(f"\nTesting large system N={n}:")
            
            try:
                atoms = [Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(n)]
                bond = MetallicBond(atoms, hopping_parameter=-1.0, periodic=True)
                
                # Check memory usage and computation time
                import time
                start_time = time.time()
                result = bond.compute_energy(method='tight_binding')
                end_time = time.time()
                
                print(f"  Computation time: {end_time - start_time:.4f} seconds")
                print(f"  Energy: {result['energy']:.4f} eV")
                
                if end_time - start_time > 10.0:  # More than 10 seconds
                    issues.append(f"N={n}: Computation too slow ({end_time - start_time:.2f}s)")
                
            except Exception as e:
                issues.append(f"N={n}: Exception during computation: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Numerical stability is good")
            self.passed_tests += 1
        
        return issues
    
    def run_all_tests(self):
        """Run all 1D chain physics tests."""
        print("="*70)
        print("SCIENTIFIC TEST 1: 1D CHAIN PHYSICS AND EDGE CASES")
        print("="*70)
        
        # Run all tests
        self.test_tight_binding_matrix_construction()
        self.test_eigenvalue_accuracy()
        self.test_band_structure_accuracy()
        self.test_metallic_character_physics()
        self.test_periodic_vs_non_periodic()
        self.test_finite_size_effects()
        self.test_numerical_stability()
        
        # Print summary
        print("\n" + "="*70)
        print("TEST 1 SUMMARY")
        print("="*70)
        
        print(f"Passed tests: {self.passed_tests}")
        print(f"Failed tests: {self.failed_tests}")
        print(f"Total issues found: {len(self.critical_issues)}")
        
        if self.critical_issues:
            print(f"\n🚨 CRITICAL ISSUES FOUND:")
            for i, issue in enumerate(self.critical_issues, 1):
                print(f"  {i}. {issue}")
        else:
            print(f"\n✅ NO CRITICAL ISSUES FOUND")
        
        return {
            'passed': self.passed_tests,
            'failed': self.failed_tests,
            'issues': self.critical_issues
        }

def main():
    """Run the 1D chain physics test suite."""
    tester = Test1DChainPhysics()
    results = tester.run_all_tests()
    return results

if __name__ == "__main__":
    main()
