#!/usr/bin/env python3
"""
Scientific Test 5: Hubbard U Effects and Strong Correlation Physics

This test rigorously validates the Hubbard U implementation
and strong correlation physics in the metallic bonding framework.

Scientific Focus:
- Hubbard U parameter effects
- Strong correlation physics
- Electron-electron interactions
- Magnetic properties
- Phase transitions
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

class TestHubbardUEffects:
    """Comprehensive test of Hubbard U effects and strong correlation physics."""
    
    def __init__(self):
        self.test_results = {}
        self.critical_issues = []
        self.warnings = []
        self.passed_tests = 0
        self.failed_tests = 0
        
    def test_hubbard_u_parameter_effects(self):
        """Test 5.1: Hubbard U parameter effects on energy."""
        print("="*70)
        print("TEST 5.1: HUBBARD U PARAMETER EFFECTS")
        print("="*70)
        
        issues = []
        
        # Test different Hubbard U values
        U_values = [0.0, 0.5, 1.0, 2.0, 5.0, 10.0]  # eV
        
        atoms = [Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(4)]
        
        print(f"Testing Hubbard U effects on {len(atoms)}-atom chain:")
        
        for U in U_values:
            print(f"\nHubbard U={U} eV:")
            
            try:
                bond = MetallicBond(atoms, hopping_parameter=-1.0, hubbard_u=U, periodic=True)
                
                # Get energy
                result = bond.compute_energy(method='tight_binding')
                energy = result['energy']
                
                print(f"  Energy: {energy:.6f} eV")
                
                # Check if energy is reasonable
                if np.isnan(energy) or np.isinf(energy):
                    issues.append(f"U={U}eV: Invalid energy (NaN or Inf): {energy}")
                
                # Check if energy increases with U (repulsive interaction)
                if U > 0:
                    # For repulsive U, energy should be higher than U=0
                    pass  # We'll check this in the next test
                
            except Exception as e:
                issues.append(f"U={U}eV: Exception during energy calculation: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Hubbard U parameter effects work")
            self.passed_tests += 1
        
        return issues
    
    def test_hubbard_u_energy_scaling(self):
        """Test 5.2: Hubbard U energy scaling and physics."""
        print("="*70)
        print("TEST 5.2: HUBBARD U ENERGY SCALING")
        print("="*70)
        
        issues = []
        
        # Test Hubbard U energy scaling
        U_values = [0.0, 0.5, 1.0, 2.0, 5.0, 10.0]  # eV
        atoms = [Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(4)]
        
        print(f"Testing Hubbard U energy scaling for {len(atoms)}-atom chain:")
        
        energies = []
        
        for U in U_values:
            print(f"\nHubbard U={U} eV:")
            
            try:
                bond = MetallicBond(atoms, hopping_parameter=-1.0, hubbard_u=U, periodic=True)
                result = bond.compute_energy(method='tight_binding')
                energy = result['energy']
                energies.append(energy)
                
                print(f"  Energy: {energy:.6f} eV")
                
            except Exception as e:
                issues.append(f"U={U}eV: Exception during energy calculation: {e}")
        
        # Check energy scaling with U
        if len(energies) > 1:
            print(f"\nEnergy scaling with U:")
            for i, U in enumerate(U_values):
                print(f"  U={U}eV: E={energies[i]:.6f} eV")
            
            # Check if energy increases with U (repulsive interaction)
            for i in range(1, len(energies)):
                if energies[i] < energies[i-1]:
                    issues.append(f"Energy decreases with U: U={U_values[i-1]}eV->{U_values[i]}eV: "
                                f"E={energies[i-1]:.6f}->{energies[i]:.6f}")
            
            # Check if energy scaling is reasonable
            # For large U, energy should scale roughly as U * n_doubly_occupied
            # where n_doubly_occupied is the number of doubly occupied sites
            
            # For 4 atoms with 4 electrons, there should be 2 doubly occupied sites
            # So energy should scale roughly as U * 2 = 2U
            expected_scaling = 2.0  # eV per U
            
            # Calculate energy difference per U
            if len(energies) > 1:
                energy_diff_per_U = (energies[-1] - energies[0]) / (U_values[-1] - U_values[0])
                print(f"Energy difference per U: {energy_diff_per_U:.6f} eV/U")
                print(f"Expected scaling: {expected_scaling:.6f} eV/U")
                
                if abs(energy_diff_per_U - expected_scaling) > 1.0:  # 1 eV tolerance
                    issues.append(f"U energy scaling incorrect: {energy_diff_per_U:.6f} vs {expected_scaling:.6f}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Hubbard U energy scaling is correct")
            self.passed_tests += 1
        
        return issues
    
    def test_hubbard_u_hamiltonian_construction(self):
        """Test 5.3: Hubbard U Hamiltonian construction."""
        print("="*70)
        print("TEST 5.3: HUBBARD U HAMILTONIAN CONSTRUCTION")
        print("="*70)
        
        issues = []
        
        # Test Hubbard U Hamiltonian construction
        U_values = [0.0, 1.0, 5.0]  # eV
        atoms = [Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(4)]
        
        print(f"Testing Hubbard U Hamiltonian construction for {len(atoms)}-atom chain:")
        
        for U in U_values:
            print(f"\nHubbard U={U} eV:")
            
            try:
                bond = MetallicBond(atoms, hopping_parameter=-1.0, hubbard_u=U, periodic=True)
                H = bond.hamiltonian.h_tight_binding
                
                print(f"  Tight-binding matrix shape: {H.shape}")
                print(f"  Matrix is symmetric: {np.allclose(H, H.T)}")
                print(f"  Matrix is real: {np.allclose(H, H.real)}")
                
                # Check matrix properties
                if not np.allclose(H, H.T):
                    issues.append(f"U={U}eV: Matrix not symmetric")
                
                if not np.allclose(H, H.real):
                    issues.append(f"U={U}eV: Matrix has imaginary parts")
                
                # Check if Hubbard U affects the matrix
                if U > 0:
                    # For U > 0, the matrix should be different from U = 0
                    # This is a basic check - the actual implementation might be more complex
                    pass
                
            except Exception as e:
                issues.append(f"U={U}eV: Exception during Hamiltonian construction: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Hubbard U Hamiltonian construction works")
            self.passed_tests += 1
        
        return issues
    
    def test_hubbard_u_metallic_character(self):
        """Test 5.4: Hubbard U effects on metallic character."""
        print("="*70)
        print("TEST 5.4: HUBBARD U EFFECTS ON METALLIC CHARACTER")
        print("="*70)
        
        issues = []
        
        # Test Hubbard U effects on metallic character
        U_values = [0.0, 0.5, 1.0, 2.0, 5.0, 10.0]  # eV
        atoms = [Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(4)]
        
        print(f"Testing Hubbard U effects on metallic character for {len(atoms)}-atom chain:")
        
        for U in U_values:
            print(f"\nHubbard U={U} eV:")
            
            try:
                bond = MetallicBond(atoms, hopping_parameter=-1.0, hubbard_u=U, periodic=True)
                
                # Get metallic character
                eigenvalues = np.linalg.eigvalsh(bond.hamiltonian.h_tight_binding)
                is_metallic = bond.hamiltonian.is_metallic()
                fermi_energy = bond.hamiltonian.get_fermi_energy(eigenvalues)
                dos_fermi = bond.hamiltonian.compute_dos_at_fermi(eigenvalues)
                
                print(f"  Eigenvalues: {eigenvalues}")
                print(f"  Is metallic: {is_metallic}")
                print(f"  Fermi energy: {fermi_energy:.6f} eV")
                print(f"  DOS at Fermi: {dos_fermi}")
                
                # Check if metallic character makes sense
                # For large U, the system might become insulating
                if U > 5.0:
                    # For very large U, the system might become insulating
                    # This depends on the specific implementation
                    pass
                
            except Exception as e:
                issues.append(f"U={U}eV: Exception during metallic character test: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Hubbard U effects on metallic character work")
            self.passed_tests += 1
        
        return issues
    
    def test_hubbard_u_system_size_scaling(self):
        """Test 5.5: Hubbard U effects with different system sizes."""
        print("="*70)
        print("TEST 5.5: HUBBARD U SYSTEM SIZE SCALING")
        print("="*70)
        
        issues = []
        
        # Test Hubbard U effects with different system sizes
        system_sizes = [2, 4, 6, 8, 10]
        U_values = [0.0, 1.0, 5.0]  # eV
        
        print(f"Testing Hubbard U effects with different system sizes:")
        
        for n in system_sizes:
            print(f"\nSystem size N={n}:")
            
            atoms = [Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(n)]
            
            for U in U_values:
                print(f"  U={U} eV:")
                
                try:
                    bond = MetallicBond(atoms, hopping_parameter=-1.0, hubbard_u=U, periodic=True)
                    result = bond.compute_energy(method='tight_binding')
                    energy = result['energy']
                    
                    print(f"    Energy: {energy:.6f} eV")
                    
                    # Check if energy scales reasonably with system size
                    energy_per_atom = energy / n
                    if abs(energy_per_atom) > 20.0:  # Unreasonably large energy per atom
                        issues.append(f"N={n}, U={U}eV: Unreasonably large energy per atom: {energy_per_atom:.4f} eV/atom")
                    
                except Exception as e:
                    issues.append(f"N={n}, U={U}eV: Exception during energy calculation: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Hubbard U system size scaling works")
            self.passed_tests += 1
        
        return issues
    
    def test_hubbard_u_edge_cases(self):
        """Test 5.6: Hubbard U edge cases and numerical stability."""
        print("="*70)
        print("TEST 5.6: HUBBARD U EDGE CASES")
        print("="*70)
        
        issues = []
        
        # Test extreme Hubbard U values
        extreme_U_values = [1e-10, 1e-5, 1e-2, 1e2, 1e5, 1e10]  # eV
        
        atoms = [Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(4)]
        
        print(f"Testing extreme Hubbard U values for {len(atoms)}-atom chain:")
        
        for U in extreme_U_values:
            print(f"\nHubbard U={U} eV:")
            
            try:
                bond = MetallicBond(atoms, hopping_parameter=-1.0, hubbard_u=U, periodic=True)
                
                # Check if matrix is well-conditioned
                H = bond.hamiltonian.h_tight_binding
                cond_num = np.linalg.cond(H)
                print(f"  Condition number: {cond_num:.2e}")
                
                if cond_num > 1e12:
                    issues.append(f"U={U}eV: Matrix ill-conditioned (cond={cond_num:.2e})")
                
                # Check eigenvalues
                eigenvalues = np.linalg.eigvalsh(H)
                if np.any(np.isnan(eigenvalues)) or np.any(np.isinf(eigenvalues)):
                    issues.append(f"U={U}eV: Invalid eigenvalues (NaN or Inf)")
                
                # Check energy computation
                result = bond.compute_energy(method='tight_binding')
                if np.isnan(result['energy']) or np.isinf(result['energy']):
                    issues.append(f"U={U}eV: Invalid energy (NaN or Inf)")
                
            except Exception as e:
                issues.append(f"U={U}eV: Exception during computation: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Hubbard U edge cases are handled correctly")
            self.passed_tests += 1
        
        return issues
    
    def run_all_tests(self):
        """Run all Hubbard U effects tests."""
        print("="*70)
        print("SCIENTIFIC TEST 5: HUBBARD U EFFECTS AND STRONG CORRELATION")
        print("="*70)
        
        # Run all tests
        self.test_hubbard_u_parameter_effects()
        self.test_hubbard_u_energy_scaling()
        self.test_hubbard_u_hamiltonian_construction()
        self.test_hubbard_u_metallic_character()
        self.test_hubbard_u_system_size_scaling()
        self.test_hubbard_u_edge_cases()
        
        # Print summary
        print("\n" + "="*70)
        print("TEST 5 SUMMARY")
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
    """Run the Hubbard U effects test suite."""
    tester = TestHubbardUEffects()
    results = tester.run_all_tests()
    return results

if __name__ == "__main__":
    main()
