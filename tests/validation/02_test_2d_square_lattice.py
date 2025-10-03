#!/usr/bin/env python3
"""
Scientific Test 2: 2D Square Lattice Physics and Band Structure

This test rigorously validates the 2D square lattice implementation
against known theoretical results and explores 2D-specific physics.

Scientific Focus:
- 2D tight-binding theory validation
- Square lattice band structure
- 2D Brillouin zone physics
- Van Hove singularities
- 2D metallic character
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

class Test2DSquareLattice:
    """Comprehensive test of 2D square lattice physics."""
    
    def __init__(self):
        self.test_results = {}
        self.critical_issues = []
        self.warnings = []
        self.passed_tests = 0
        self.failed_tests = 0
        
    def test_2d_lattice_construction(self):
        """Test 2.1: 2D square lattice construction."""
        print("="*70)
        print("TEST 2.1: 2D SQUARE LATTICE CONSTRUCTION")
        print("="*70)
        
        issues = []
        
        # Test different 2D lattice sizes
        test_sizes = [(2, 2), (3, 3), (4, 4), (2, 3), (3, 4), (4, 2)]
        
        for nx, ny in test_sizes:
            print(f"\nTesting {nx}×{ny} lattice:")
            
            # Create 2D square lattice
            atoms = []
            for i in range(nx):
                for j in range(ny):
                    atoms.append(Atom('Na', position=np.array([i*3.0, j*3.0, 0.0])))
            
            try:
                bond = MetallicBond(atoms, lattice_type='2d_square', hopping_parameter=-1.0, periodic=True)
                H = bond.hamiltonian.h_tight_binding
                
                print(f"  Total atoms: {len(atoms)}")
                print(f"  Matrix shape: {H.shape}")
                print(f"  Matrix is symmetric: {np.allclose(H, H.T)}")
                print(f"  Matrix is real: {np.allclose(H, H.real)}")
                
                # Check matrix properties
                if not np.allclose(H, H.T):
                    issues.append(f"{nx}×{ny}: Matrix not symmetric")
                
                if not np.allclose(H, H.real):
                    issues.append(f"{nx}×{ny}: Matrix has imaginary parts")
                
                # Check that matrix has correct size
                if H.shape != (len(atoms), len(atoms)):
                    issues.append(f"{nx}×{ny}: Matrix size incorrect: {H.shape} vs ({len(atoms)}, {len(atoms)})")
                
                # Check diagonal elements (should be zero for tight-binding)
                if not np.allclose(np.diag(H), 0.0):
                    issues.append(f"{nx}×{ny}: Diagonal elements not zero: {np.diag(H)}")
                
            except Exception as e:
                issues.append(f"{nx}×{ny}: Exception during construction: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ All 2D lattice constructions are correct")
            self.passed_tests += 1
        
        return issues
    
    def test_2d_band_structure_physics(self):
        """Test 2.2: 2D band structure physics."""
        print("="*70)
        print("TEST 2.2: 2D BAND STRUCTURE PHYSICS")
        print("="*70)
        
        issues = []
        
        # Test 2D square lattice band structure
        atoms = [Atom('Na', position=np.array([i*3.0, 0.0, 0.0])) for i in range(4)]  # 1D for now
        bond = MetallicBond(atoms, lattice_type='2d_square', hopping_parameter=-1.0, periodic=True)
        
        try:
            band_structure = bond.get_band_structure()
            print(f"Band structure keys: {list(band_structure.keys())}")
            
            # Check if 2D band structure is returned
            if 'kx' in band_structure and 'ky' in band_structure:
                kx = band_structure['kx']
                ky = band_structure['ky']
                energies = band_structure['energies']
                
                print(f"kx points: {len(kx)} from {kx[0]:.3f} to {kx[-1]:.3f}")
                print(f"ky points: {len(ky)} from {ky[0]:.3f} to {ky[-1]:.3f}")
                print(f"Energy shape: {energies.shape}")
                
                # Test 2D dispersion relation: E(k) = -2t[cos(kx) + cos(ky)]
                t = -1.0
                max_error = 0.0
                
                for i, kx_val in enumerate(kx):
                    for j, ky_val in enumerate(ky):
                        computed_E = energies[i, j]
                        expected_E = 2 * t * (np.cos(kx_val) + np.cos(ky_val))
                        error = abs(computed_E - expected_E)
                        max_error = max(max_error, error)
                        
                        if error > 1e-10:
                            issues.append(f"k=({kx_val:.3f}, {ky_val:.3f}): E={computed_E:.6f}, expected={expected_E:.6f}, error={error:.2e}")
                
                print(f"Maximum 2D dispersion error: {max_error:.2e}")
                
                # Test specific k-points
                # Γ point (0, 0)
                gamma_idx = (len(kx)//2, len(ky)//2)
                E_gamma = energies[gamma_idx[0], gamma_idx[1]]
                expected_E_gamma = 2 * t * (np.cos(0) + np.cos(0))  # = 4t = -4 eV
                
                print(f"E(Γ): {E_gamma:.6f} eV, expected: {expected_E_gamma:.6f} eV")
                
                if abs(E_gamma - expected_E_gamma) > 1e-10:
                    issues.append(f"E(Γ) incorrect: {E_gamma} vs {expected_E_gamma}")
                
                # M point (π, π)
                M_idx = (-1, -1)  # Last point
                E_M = energies[M_idx[0], M_idx[1]]
                expected_E_M = 2 * t * (np.cos(np.pi) + np.cos(np.pi))  # = -4t = +4 eV
                
                print(f"E(M): {E_M:.6f} eV, expected: {expected_E_M:.6f} eV")
                
                if abs(E_M - expected_E_M) > 1e-10:
                    issues.append(f"E(M) incorrect: {E_M} vs {expected_E_M}")
                
                # X point (π, 0)
                X_idx = (-1, len(ky)//2)
                E_X = energies[X_idx[0], X_idx[1]]
                expected_E_X = 2 * t * (np.cos(np.pi) + np.cos(0))  # = 0 eV
                
                print(f"E(X): {E_X:.6f} eV, expected: {expected_E_X:.6f} eV")
                
                if abs(E_X - expected_E_X) > 1e-10:
                    issues.append(f"E(X) incorrect: {E_X} vs {expected_E_X}")
                
            else:
                issues.append("2D band structure not implemented - missing kx, ky keys")
                
        except Exception as e:
            issues.append(f"Exception during 2D band structure test: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues[:5]:  # Show first 5 issues
                print(f"  - {issue}")
            if len(issues) > 5:
                print(f"  ... and {len(issues) - 5} more issues")
            self.failed_tests += 1
        else:
            print(f"\n✅ 2D band structure is correct")
            self.passed_tests += 1
        
        return issues
    
    def test_2d_metallic_character(self):
        """Test 2.3: 2D metallic character detection."""
        print("="*70)
        print("TEST 2.3: 2D METALLIC CHARACTER DETECTION")
        print("="*70)
        
        issues = []
        
        # Test different 2D system sizes
        test_sizes = [(2, 2), (3, 3), (4, 4), (2, 3), (3, 4)]
        
        for nx, ny in test_sizes:
            print(f"\nTesting {nx}×{ny} lattice:")
            
            # Create 2D square lattice
            atoms = []
            for i in range(nx):
                for j in range(ny):
                    atoms.append(Atom('Na', position=np.array([i*3.0, j*3.0, 0.0])))
            
            try:
                bond = MetallicBond(atoms, lattice_type='2d_square', hopping_parameter=-1.0, periodic=True)
                
                # Get eigenvalues and check metallic character
                eigenvalues = np.linalg.eigvalsh(bond.hamiltonian.h_tight_binding)
                is_metallic = bond.hamiltonian.is_metallic()
                fermi_energy = bond.hamiltonian.get_fermi_energy(eigenvalues)
                dos_fermi = bond.hamiltonian.compute_dos_at_fermi(eigenvalues)
                
                print(f"  Total atoms: {len(atoms)}")
                print(f"  Eigenvalues: {eigenvalues}")
                print(f"  Is metallic: {is_metallic}")
                print(f"  Fermi energy: {fermi_energy:.4f} eV")
                print(f"  DOS at Fermi: {dos_fermi}")
                
                # Check if metallic character makes sense
                n_electrons = len(atoms)  # 1 electron per atom
                n_bands_occupied = int(np.ceil(n_electrons / 2.0))
                
                # For 2D systems, metallic character depends on band filling
                # Check if there are states at Fermi level
                states_at_fermi = np.abs(eigenvalues - fermi_energy) < 1e-12
                expected_metallic = np.sum(states_at_fermi) > 0
                
                if is_metallic != expected_metallic:
                    issues.append(f"{nx}×{ny}: Metallic={is_metallic}, expected={expected_metallic}, "
                                f"DOS={dos_fermi}, states_at_fermi={np.sum(states_at_fermi)}")
                
            except Exception as e:
                issues.append(f"{nx}×{ny}: Exception during metallic character test: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ 2D metallic character detection is correct")
            self.passed_tests += 1
        
        return issues
    
    def test_2d_energy_scaling(self):
        """Test 2.4: 2D energy scaling with system size."""
        print("="*70)
        print("TEST 2.4: 2D ENERGY SCALING")
        print("="*70)
        
        issues = []
        
        # Test different 2D system sizes
        test_sizes = [(2, 2), (3, 3), (4, 4), (5, 5), (6, 6)]
        energies_per_atom = []
        
        for nx, ny in test_sizes:
            print(f"\nTesting {nx}×{ny} lattice:")
            
            # Create 2D square lattice
            atoms = []
            for i in range(nx):
                for j in range(ny):
                    atoms.append(Atom('Na', position=np.array([i*3.0, j*3.0, 0.0])))
            
            try:
                bond = MetallicBond(atoms, lattice_type='2d_square', hopping_parameter=-1.0, periodic=True)
                result = bond.compute_energy(method='tight_binding')
                
                energy_per_atom = result['energy'] / len(atoms)
                energies_per_atom.append(energy_per_atom)
                
                print(f"  Total atoms: {len(atoms)}")
                print(f"  Total energy: {result['energy']:.4f} eV")
                print(f"  Energy per atom: {energy_per_atom:.4f} eV/atom")
                
            except Exception as e:
                issues.append(f"{nx}×{ny}: Exception during energy computation: {e}")
        
        # Check energy per atom scaling
        if len(energies_per_atom) > 1:
            print(f"\nEnergy per atom scaling:")
            for i, (nx, ny) in enumerate(test_sizes):
                print(f"  {nx}×{ny}: {energies_per_atom[i]:.4f} eV/atom")
            
            # Check if energy per atom is reasonable
            # For 2D square lattice, energy per atom should be around -4|t| = -4 eV
            expected_2d = -4.0  # eV/atom
            
            # Use the last few points to estimate large-N limit
            last_energies = energies_per_atom[-2:]
            avg_large_N = np.mean(last_energies)
            
            print(f"\nLarge-N limit estimate: {avg_large_N:.4f} eV/atom")
            print(f"Expected 2D limit: {expected_2d:.4f} eV/atom")
            print(f"Error: {abs(avg_large_N - expected_2d):.4f} eV/atom")
            
            if abs(avg_large_N - expected_2d) > 0.5:  # 0.5 eV tolerance
                issues.append(f"2D large-N limit incorrect: {avg_large_N:.4f} vs {expected_2d:.4f}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ 2D energy scaling is correct")
            self.passed_tests += 1
        
        return issues
    
    def test_2d_brillouin_zone_symmetry(self):
        """Test 2.5: 2D Brillouin zone symmetry."""
        print("="*70)
        print("TEST 2.5: 2D BRILLOUIN ZONE SYMMETRY")
        print("="*70)
        
        issues = []
        
        # Test 2D square lattice band structure
        atoms = [Atom('Na', position=np.array([i*3.0, 0.0, 0.0])) for i in range(4)]  # 1D for now
        bond = MetallicBond(atoms, lattice_type='2d_square', hopping_parameter=-1.0, periodic=True)
        
        try:
            band_structure = bond.get_band_structure()
            
            if 'kx' in band_structure and 'ky' in band_structure:
                kx = band_structure['kx']
                ky = band_structure['ky']
                energies = band_structure['energies']
                
                # Test symmetry properties
                # E(kx, ky) = E(-kx, ky) = E(kx, -ky) = E(-kx, -ky)
                
                n_kx, n_ky = len(kx), len(ky)
                
                # Test E(kx, ky) = E(-kx, ky)
                for i in range(n_kx):
                    for j in range(n_ky):
                        kx_val = kx[i]
                        ky_val = ky[j]
                        
                        # Find corresponding -kx point
                        kx_neg_idx = np.argmin(np.abs(kx + kx_val))
                        E_pos = energies[i, j]
                        E_neg = energies[kx_neg_idx, j]
                        
                        if abs(E_pos - E_neg) > 1e-10:
                            issues.append(f"E({kx_val:.3f}, {ky_val:.3f}) ≠ E({-kx_val:.3f}, {ky_val:.3f}): {E_pos:.6f} vs {E_neg:.6f}")
                
                # Test E(kx, ky) = E(kx, -ky)
                for i in range(n_kx):
                    for j in range(n_ky):
                        kx_val = kx[i]
                        ky_val = ky[j]
                        
                        # Find corresponding -ky point
                        ky_neg_idx = np.argmin(np.abs(ky + ky_val))
                        E_pos = energies[i, j]
                        E_neg = energies[i, ky_neg_idx]
                        
                        if abs(E_pos - E_neg) > 1e-10:
                            issues.append(f"E({kx_val:.3f}, {ky_val:.3f}) ≠ E({kx_val:.3f}, {-ky_val:.3f}): {E_pos:.6f} vs {E_neg:.6f}")
                
                print(f"Tested {n_kx}×{n_ky} k-points for symmetry")
                
            else:
                issues.append("2D band structure not implemented - cannot test symmetry")
                
        except Exception as e:
            issues.append(f"Exception during symmetry test: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues[:5]:  # Show first 5 issues
                print(f"  - {issue}")
            if len(issues) > 5:
                print(f"  ... and {len(issues) - 5} more issues")
            self.failed_tests += 1
        else:
            print(f"\n✅ 2D Brillouin zone symmetry is correct")
            self.passed_tests += 1
        
        return issues
    
    def test_2d_edge_cases(self):
        """Test 2.6: 2D edge cases and numerical stability."""
        print("="*70)
        print("TEST 2.6: 2D EDGE CASES AND NUMERICAL STABILITY")
        print("="*70)
        
        issues = []
        
        # Test extreme hopping parameters
        extreme_hopping = [1e-10, 1e-5, 1e-2, 1e2, 1e5, 1e10]
        
        for t in extreme_hopping:
            print(f"\nTesting hopping parameter t={t}:")
            
            try:
                atoms = [Atom('Na', position=np.array([i*3.0, 0.0, 0.0])) for i in range(4)]
                bond = MetallicBond(atoms, lattice_type='2d_square', hopping_parameter=t, periodic=True)
                
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
        
        # Test very large 2D systems
        large_sizes = [(10, 10), (15, 15), (20, 20)]
        
        for nx, ny in large_sizes:
            print(f"\nTesting large 2D system {nx}×{ny}:")
            
            try:
                atoms = []
                for i in range(nx):
                    for j in range(ny):
                        atoms.append(Atom('Na', position=np.array([i*3.0, j*3.0, 0.0])))
                
                bond = MetallicBond(atoms, lattice_type='2d_square', hopping_parameter=-1.0, periodic=True)
                
                # Check memory usage and computation time
                import time
                start_time = time.time()
                result = bond.compute_energy(method='tight_binding')
                end_time = time.time()
                
                print(f"  Total atoms: {len(atoms)}")
                print(f"  Computation time: {end_time - start_time:.4f} seconds")
                print(f"  Energy: {result['energy']:.4f} eV")
                
                if end_time - start_time > 30.0:  # More than 30 seconds
                    issues.append(f"{nx}×{ny}: Computation too slow ({end_time - start_time:.2f}s)")
                
            except Exception as e:
                issues.append(f"{nx}×{ny}: Exception during computation: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ 2D numerical stability is good")
            self.passed_tests += 1
        
        return issues
    
    def run_all_tests(self):
        """Run all 2D square lattice tests."""
        print("="*70)
        print("SCIENTIFIC TEST 2: 2D SQUARE LATTICE PHYSICS")
        print("="*70)
        
        # Run all tests
        self.test_2d_lattice_construction()
        self.test_2d_band_structure_physics()
        self.test_2d_metallic_character()
        self.test_2d_energy_scaling()
        self.test_2d_brillouin_zone_symmetry()
        self.test_2d_edge_cases()
        
        # Print summary
        print("\n" + "="*70)
        print("TEST 2 SUMMARY")
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
    """Run the 2D square lattice test suite."""
    tester = Test2DSquareLattice()
    results = tester.run_all_tests()
    return results

if __name__ == "__main__":
    main()
