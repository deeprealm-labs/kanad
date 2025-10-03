#!/usr/bin/env python3
"""
Scientific Test 3: 3D Cubic Lattice and High-Dimensional Physics

This test rigorously validates the 3D cubic lattice implementation
and explores high-dimensional physics effects.

Scientific Focus:
- 3D tight-binding theory validation
- Cubic lattice band structure
- 3D Brillouin zone physics
- High-dimensional scaling
- 3D metallic character
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

class Test3DCubicLattice:
    """Comprehensive test of 3D cubic lattice physics."""
    
    def __init__(self):
        self.test_results = {}
        self.critical_issues = []
        self.warnings = []
        self.passed_tests = 0
        self.failed_tests = 0
        
    def test_3d_lattice_construction(self):
        """Test 3.1: 3D cubic lattice construction."""
        print("="*70)
        print("TEST 3.1: 3D CUBIC LATTICE CONSTRUCTION")
        print("="*70)
        
        issues = []
        
        # Test different 3D lattice sizes
        test_sizes = [(2, 2, 2), (3, 3, 3), (2, 2, 3), (2, 3, 4)]
        
        for nx, ny, nz in test_sizes:
            print(f"\nTesting {nx}×{ny}×{nz} lattice:")
            
            # Create 3D cubic lattice
            atoms = []
            for i in range(nx):
                for j in range(ny):
                    for k in range(nz):
                        atoms.append(Atom('Na', position=np.array([i*3.0, j*3.0, k*3.0])))
            
            try:
                bond = MetallicBond(atoms, lattice_type='3d_cubic', hopping_parameter=-1.0, periodic=True)
                H = bond.hamiltonian.h_tight_binding
                
                print(f"  Total atoms: {len(atoms)}")
                print(f"  Matrix shape: {H.shape}")
                print(f"  Matrix is symmetric: {np.allclose(H, H.T)}")
                print(f"  Matrix is real: {np.allclose(H, H.real)}")
                
                # Check matrix properties
                if not np.allclose(H, H.T):
                    issues.append(f"{nx}×{ny}×{nz}: Matrix not symmetric")
                
                if not np.allclose(H, H.real):
                    issues.append(f"{nx}×{ny}×{nz}: Matrix has imaginary parts")
                
                # Check that matrix has correct size
                if H.shape != (len(atoms), len(atoms)):
                    issues.append(f"{nx}×{ny}×{nz}: Matrix size incorrect: {H.shape} vs ({len(atoms)}, {len(atoms)})")
                
                # Check diagonal elements (should be zero for tight-binding)
                if not np.allclose(np.diag(H), 0.0):
                    issues.append(f"{nx}×{ny}×{nz}: Diagonal elements not zero: {np.diag(H)}")
                
            except Exception as e:
                issues.append(f"{nx}×{ny}×{nz}: Exception during construction: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ All 3D lattice constructions are correct")
            self.passed_tests += 1
        
        return issues
    
    def test_3d_band_structure_physics(self):
        """Test 3.2: 3D band structure physics."""
        print("="*70)
        print("TEST 3.2: 3D BAND STRUCTURE PHYSICS")
        print("="*70)
        
        issues = []
        
        # Test 3D cubic lattice band structure
        atoms = [Atom('Na', position=np.array([i*3.0, 0.0, 0.0])) for i in range(4)]  # 1D for now
        bond = MetallicBond(atoms, lattice_type='3d_cubic', hopping_parameter=-1.0, periodic=True)
        
        try:
            band_structure = bond.get_band_structure()
            print(f"Band structure keys: {list(band_structure.keys())}")
            
            # Check if 3D band structure is returned
            if 'kx' in band_structure and 'ky' in band_structure and 'kz' in band_structure:
                kx = band_structure['kx']
                ky = band_structure['ky']
                kz = band_structure['kz']
                energies = band_structure['energies']
                
                print(f"kx points: {len(kx)} from {kx[0]:.3f} to {kx[-1]:.3f}")
                print(f"ky points: {len(ky)} from {ky[0]:.3f} to {ky[-1]:.3f}")
                print(f"kz points: {len(kz)} from {kz[0]:.3f} to {kz[-1]:.3f}")
                print(f"Energy shape: {energies.shape}")
                
                # Test 3D dispersion relation: E(k) = -2t[cos(kx) + cos(ky) + cos(kz)]
                t = -1.0
                max_error = 0.0
                
                for i, kx_val in enumerate(kx):
                    for j, ky_val in enumerate(ky):
                        for k, kz_val in enumerate(kz):
                            computed_E = energies[i, j, k]
                            expected_E = 2 * t * (np.cos(kx_val) + np.cos(ky_val) + np.cos(kz_val))
                            error = abs(computed_E - expected_E)
                            max_error = max(max_error, error)
                            
                            if error > 1e-10:
                                issues.append(f"k=({kx_val:.3f}, {ky_val:.3f}, {kz_val:.3f}): E={computed_E:.6f}, expected={expected_E:.6f}, error={error:.2e}")
                
                print(f"Maximum 3D dispersion error: {max_error:.2e}")
                
                # Test specific k-points
                # Γ point (0, 0, 0)
                gamma_idx = (len(kx)//2, len(ky)//2, len(kz)//2)
                E_gamma = energies[gamma_idx[0], gamma_idx[1], gamma_idx[2]]
                expected_E_gamma = 2 * t * (np.cos(0) + np.cos(0) + np.cos(0))  # = 6t = -6 eV
                
                print(f"E(Γ): {E_gamma:.6f} eV, expected: {expected_E_gamma:.6f} eV")
                
                if abs(E_gamma - expected_E_gamma) > 1e-10:
                    issues.append(f"E(Γ) incorrect: {E_gamma} vs {expected_E_gamma}")
                
                # R point (π, π, π)
                R_idx = (-1, -1, -1)  # Last point
                E_R = energies[R_idx[0], R_idx[1], R_idx[2]]
                expected_E_R = 2 * t * (np.cos(np.pi) + np.cos(np.pi) + np.cos(np.pi))  # = -6t = +6 eV
                
                print(f"E(R): {E_R:.6f} eV, expected: {expected_E_R:.6f} eV")
                
                if abs(E_R - expected_E_R) > 1e-10:
                    issues.append(f"E(R) incorrect: {E_R} vs {expected_E_R}")
                
                # X point (π, 0, 0)
                X_idx = (-1, len(ky)//2, len(kz)//2)
                E_X = energies[X_idx[0], X_idx[1], X_idx[2]]
                expected_E_X = 2 * t * (np.cos(np.pi) + np.cos(0) + np.cos(0))  # = -2t = +2 eV
                
                print(f"E(X): {E_X:.6f} eV, expected: {expected_E_X:.6f} eV")
                
                if abs(E_X - expected_E_X) > 1e-10:
                    issues.append(f"E(X) incorrect: {E_X} vs {expected_E_X}")
                
            else:
                issues.append("3D band structure not implemented - missing kx, ky, kz keys")
                
        except Exception as e:
            issues.append(f"Exception during 3D band structure test: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues[:5]:  # Show first 5 issues
                print(f"  - {issue}")
            if len(issues) > 5:
                print(f"  ... and {len(issues) - 5} more issues")
            self.failed_tests += 1
        else:
            print(f"\n✅ 3D band structure is correct")
            self.passed_tests += 1
        
        return issues
    
    def test_3d_energy_scaling(self):
        """Test 3.3: 3D energy scaling with system size."""
        print("="*70)
        print("TEST 3.3: 3D ENERGY SCALING")
        print("="*70)
        
        issues = []
        
        # Test different 3D system sizes
        test_sizes = [(2, 2, 2), (3, 3, 3), (4, 4, 4), (5, 5, 5)]
        energies_per_atom = []
        
        for nx, ny, nz in test_sizes:
            print(f"\nTesting {nx}×{ny}×{nz} lattice:")
            
            # Create 3D cubic lattice
            atoms = []
            for i in range(nx):
                for j in range(ny):
                    for k in range(nz):
                        atoms.append(Atom('Na', position=np.array([i*3.0, j*3.0, k*3.0])))
            
            try:
                bond = MetallicBond(atoms, lattice_type='3d_cubic', hopping_parameter=-1.0, periodic=True)
                result = bond.compute_energy(method='tight_binding')
                
                energy_per_atom = result['energy'] / len(atoms)
                energies_per_atom.append(energy_per_atom)
                
                print(f"  Total atoms: {len(atoms)}")
                print(f"  Total energy: {result['energy']:.4f} eV")
                print(f"  Energy per atom: {energy_per_atom:.4f} eV/atom")
                
            except Exception as e:
                issues.append(f"{nx}×{ny}×{nz}: Exception during energy computation: {e}")
        
        # Check energy per atom scaling
        if len(energies_per_atom) > 1:
            print(f"\nEnergy per atom scaling:")
            for i, (nx, ny, nz) in enumerate(test_sizes):
                print(f"  {nx}×{ny}×{nz}: {energies_per_atom[i]:.4f} eV/atom")
            
            # Check if energy per atom is reasonable
            # For 3D cubic lattice, energy per atom should be around -6|t| = -6 eV
            expected_3d = -6.0  # eV/atom
            
            # Use the last few points to estimate large-N limit
            last_energies = energies_per_atom[-2:]
            avg_large_N = np.mean(last_energies)
            
            print(f"\nLarge-N limit estimate: {avg_large_N:.4f} eV/atom")
            print(f"Expected 3D limit: {expected_3d:.4f} eV/atom")
            print(f"Error: {abs(avg_large_N - expected_3d):.4f} eV/atom")
            
            if abs(avg_large_N - expected_3d) > 1.0:  # 1 eV tolerance
                issues.append(f"3D large-N limit incorrect: {avg_large_N:.4f} vs {expected_3d:.4f}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ 3D energy scaling is correct")
            self.passed_tests += 1
        
        return issues
    
    def test_3d_metallic_character(self):
        """Test 3.4: 3D metallic character detection."""
        print("="*70)
        print("TEST 3.4: 3D METALLIC CHARACTER DETECTION")
        print("="*70)
        
        issues = []
        
        # Test different 3D system sizes
        test_sizes = [(2, 2, 2), (3, 3, 3), (4, 4, 4), (2, 2, 3), (2, 3, 4)]
        
        for nx, ny, nz in test_sizes:
            print(f"\nTesting {nx}×{ny}×{nz} lattice:")
            
            # Create 3D cubic lattice
            atoms = []
            for i in range(nx):
                for j in range(ny):
                    for k in range(nz):
                        atoms.append(Atom('Na', position=np.array([i*3.0, j*3.0, k*3.0])))
            
            try:
                bond = MetallicBond(atoms, lattice_type='3d_cubic', hopping_parameter=-1.0, periodic=True)
                
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
                
                # For 3D systems, metallic character depends on band filling
                # Check if there are states at Fermi level
                states_at_fermi = np.abs(eigenvalues - fermi_energy) < 1e-12
                expected_metallic = np.sum(states_at_fermi) > 0
                
                if is_metallic != expected_metallic:
                    issues.append(f"{nx}×{ny}×{nz}: Metallic={is_metallic}, expected={expected_metallic}, "
                                f"DOS={dos_fermi}, states_at_fermi={np.sum(states_at_fermi)}")
                
            except Exception as e:
                issues.append(f"{nx}×{ny}×{nz}: Exception during metallic character test: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ 3D metallic character detection is correct")
            self.passed_tests += 1
        
        return issues
    
    def test_3d_brillouin_zone_symmetry(self):
        """Test 3.5: 3D Brillouin zone symmetry."""
        print("="*70)
        print("TEST 3.5: 3D BRILLOUIN ZONE SYMMETRY")
        print("="*70)
        
        issues = []
        
        # Test 3D cubic lattice band structure
        atoms = [Atom('Na', position=np.array([i*3.0, 0.0, 0.0])) for i in range(4)]  # 1D for now
        bond = MetallicBond(atoms, lattice_type='3d_cubic', hopping_parameter=-1.0, periodic=True)
        
        try:
            band_structure = bond.get_band_structure()
            
            if 'kx' in band_structure and 'ky' in band_structure and 'kz' in band_structure:
                kx = band_structure['kx']
                ky = band_structure['ky']
                kz = band_structure['kz']
                energies = band_structure['energies']
                
                # Test symmetry properties
                # E(kx, ky, kz) = E(-kx, ky, kz) = E(kx, -ky, kz) = E(kx, ky, -kz)
                
                n_kx, n_ky, n_kz = len(kx), len(ky), len(kz)
                
                # Test E(kx, ky, kz) = E(-kx, ky, kz)
                for i in range(n_kx):
                    for j in range(n_ky):
                        for k in range(n_kz):
                            kx_val = kx[i]
                            ky_val = ky[j]
                            kz_val = kz[k]
                            
                            # Find corresponding -kx point
                            kx_neg_idx = np.argmin(np.abs(kx + kx_val))
                            E_pos = energies[i, j, k]
                            E_neg = energies[kx_neg_idx, j, k]
                            
                            if abs(E_pos - E_neg) > 1e-10:
                                issues.append(f"E({kx_val:.3f}, {ky_val:.3f}, {kz_val:.3f}) ≠ E({-kx_val:.3f}, {ky_val:.3f}, {kz_val:.3f}): {E_pos:.6f} vs {E_neg:.6f}")
                
                print(f"Tested {n_kx}×{n_ky}×{n_kz} k-points for symmetry")
                
            else:
                issues.append("3D band structure not implemented - cannot test symmetry")
                
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
            print(f"\n✅ 3D Brillouin zone symmetry is correct")
            self.passed_tests += 1
        
        return issues
    
    def test_3d_edge_cases(self):
        """Test 3.6: 3D edge cases and numerical stability."""
        print("="*70)
        print("TEST 3.6: 3D EDGE CASES AND NUMERICAL STABILITY")
        print("="*70)
        
        issues = []
        
        # Test extreme hopping parameters
        extreme_hopping = [1e-10, 1e-5, 1e-2, 1e2, 1e5, 1e10]
        
        for t in extreme_hopping:
            print(f"\nTesting hopping parameter t={t}:")
            
            try:
                atoms = [Atom('Na', position=np.array([i*3.0, 0.0, 0.0])) for i in range(4)]
                bond = MetallicBond(atoms, lattice_type='3d_cubic', hopping_parameter=t, periodic=True)
                
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
        
        # Test very large 3D systems
        large_sizes = [(5, 5, 5), (6, 6, 6), (7, 7, 7)]
        
        for nx, ny, nz in large_sizes:
            print(f"\nTesting large 3D system {nx}×{ny}×{nz}:")
            
            try:
                atoms = []
                for i in range(nx):
                    for j in range(ny):
                        for k in range(nz):
                            atoms.append(Atom('Na', position=np.array([i*3.0, j*3.0, k*3.0])))
                
                bond = MetallicBond(atoms, lattice_type='3d_cubic', hopping_parameter=-1.0, periodic=True)
                
                # Check memory usage and computation time
                import time
                start_time = time.time()
                result = bond.compute_energy(method='tight_binding')
                end_time = time.time()
                
                print(f"  Total atoms: {len(atoms)}")
                print(f"  Computation time: {end_time - start_time:.4f} seconds")
                print(f"  Energy: {result['energy']:.4f} eV")
                
                if end_time - start_time > 60.0:  # More than 60 seconds
                    issues.append(f"{nx}×{ny}×{nz}: Computation too slow ({end_time - start_time:.2f}s)")
                
            except Exception as e:
                issues.append(f"{nx}×{ny}×{nz}: Exception during computation: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ 3D numerical stability is good")
            self.passed_tests += 1
        
        return issues
    
    def run_all_tests(self):
        """Run all 3D cubic lattice tests."""
        print("="*70)
        print("SCIENTIFIC TEST 3: 3D CUBIC LATTICE PHYSICS")
        print("="*70)
        
        # Run all tests
        self.test_3d_lattice_construction()
        self.test_3d_band_structure_physics()
        self.test_3d_energy_scaling()
        self.test_3d_metallic_character()
        self.test_3d_brillouin_zone_symmetry()
        self.test_3d_edge_cases()
        
        # Print summary
        print("\n" + "="*70)
        print("TEST 3 SUMMARY")
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
    """Run the 3D cubic lattice test suite."""
    tester = Test3DCubicLattice()
    results = tester.run_all_tests()
    return results

if __name__ == "__main__":
    main()
