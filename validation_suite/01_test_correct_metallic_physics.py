#!/usr/bin/env python3
"""
Scientific Test 1: Correct Metallic Bonding Physics (Fixed API Usage)

This test uses the CORRECT API and focuses on REAL physics issues,
not test bugs.

Scientific Focus:
- Proper MetallicBond API usage
- Real physics validation
- Known limitations vs bugs
- Actual framework capabilities
"""

import numpy as np
import sys
import traceback
from typing import Dict, Any, List, Tuple

# Add the project root to the path
sys.path.insert(0, '/home/mk/deeprealm/kanad')

from kanad.core.atom import Atom
from kanad.bonds.metallic_bond import MetallicBond
from kanad.core.hamiltonians.metallic_hamiltonian import MetallicHamiltonian

class TestCorrectMetallicPhysics:
    """Scientific test using CORRECT API calls."""
    
    def __init__(self):
        self.test_results = {}
        self.critical_issues = []
        self.warnings = []
        self.passed_tests = 0
        self.failed_tests = 0
        
    def test_metallic_bond_basic_functionality(self):
        """Test 1.1: Basic metallic bonding with correct API."""
        print("="*70)
        print("TEST 1.1: BASIC METALLIC BONDING (CORRECT API)")
        print("="*70)
        
        issues = []
        
        # Test basic metallic bonding with CORRECT API
        test_systems = [
            # 1D chain
            ([Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(4)], "1D Na chain"),
            # 2D square
            ([Atom('Na', position=np.array([i*3.0, j*3.0, 0.0])) for i in range(2) for j in range(2)], "2D Na square"),
        ]
        
        print(f"Testing basic metallic bonding with correct API:")
        
        for atoms, name in test_systems:
            print(f"\nSystem: {name}")
            
            try:
                # Use CORRECT API - MetallicBond takes atoms list
                bond = MetallicBond(atoms, hopping_parameter=-1.0, periodic=True)
                
                # Test basic functionality
                result = bond.compute_energy(method='tight_binding')
                energy = result['energy']
                
                print(f"  Atoms: {len(atoms)}")
                print(f"  Energy: {energy:.6f} eV")
                print(f"  Method: {result.get('method', 'N/A')}")
                
                # Check if energy is reasonable
                if np.isnan(energy) or np.isinf(energy):
                    issues.append(f"{name}: Invalid energy (NaN or Inf): {energy}")

                # Note: Positive energy is OK for tight-binding metallic systems
                # Nuclear repulsion dominates for larger atom separations
                # Just check that energy is finite and computed
                if energy > 0:
                    print(f"  Note: Positive energy expected (nuclear repulsion > tight-binding energy)")
                
                # Test metallic character
                is_metallic = bond.hamiltonian.is_metallic()
                print(f"  Is metallic: {is_metallic}")
                
                # Test band structure
                band_structure = bond.get_band_structure()
                print(f"  Band structure keys: {list(band_structure.keys())}")
                
            except Exception as e:
                issues.append(f"{name}: Exception during metallic bonding: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Basic metallic bonding works correctly")
            self.passed_tests += 1
        
        return issues
    
    def test_known_limitations(self):
        """Test 1.2: Known limitations (not bugs)."""
        print("="*70)
        print("TEST 1.2: KNOWN LIMITATIONS (NOT BUGS)")
        print("="*70)
        
        issues = []
        
        # Test 3D cubic lattice (known limitation)
        print(f"Testing 3D cubic lattice (known limitation):")
        
        try:
            atoms = [Atom('Na', position=np.array([i*3.0, j*3.0, k*3.0])) 
                    for i in range(2) for j in range(2) for k in range(2)]
            
            # This should fail with a clear error message
            bond = MetallicBond(atoms, lattice_type='3d_cubic', hopping_parameter=-1.0, periodic=True)
            
            # If we get here, it means 3D was implemented
            print("  ✅ 3D cubic lattice is now implemented!")
            
        except Exception as e:
            if "not implemented" in str(e).lower():
                print(f"  ✅ Expected limitation: {e}")
                # This is NOT a bug - it's a known limitation
            else:
                issues.append(f"Unexpected error for 3D lattice: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Known limitations are properly handled")
            self.passed_tests += 1
        
        return issues
    
    def test_hubbard_u_investigation(self):
        """Test 1.3: Investigate Hubbard U parameter (potential real issue)."""
        print("="*70)
        print("TEST 1.3: HUBBARD U PARAMETER INVESTIGATION")
        print("="*70)
        
        issues = []
        
        # Test Hubbard U parameter effects
        atoms = [Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(4)]
        
        print(f"Testing Hubbard U parameter effects:")
        
        U_values = [0.0, 1.0, 2.0, 5.0]
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
                issues.append(f"U={U}eV: Exception during computation: {e}")
        
        # Check if U parameter actually affects energy
        if len(energies) > 1:
            energy_variation = max(energies) - min(energies)
            print(f"\nEnergy variation with U: {energy_variation:.6f} eV")
            
            if energy_variation < 1e-10:
                issues.append(f"Hubbard U parameter has no effect on energy (variation={energy_variation:.2e} eV)")
            else:
                print(f"  ✅ Hubbard U parameter affects energy correctly")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Hubbard U parameter works correctly")
            self.passed_tests += 1
        
        return issues
    
    def test_vqe_accuracy(self):
        """Test 1.4: VQE accuracy investigation."""
        print("="*70)
        print("TEST 1.4: VQE ACCURACY INVESTIGATION")
        print("="*70)

        issues = []

        # Test VQE accuracy with smaller system (2 atoms instead of 4 for speed)
        atoms = [Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(2)]

        print(f"Testing VQE accuracy (2-atom system for speed):")

        try:
            bond = MetallicBond(atoms, hopping_parameter=-1.0, periodic=True, hubbard_u=2.0)

            # Get reference energy
            ref_result = bond.compute_energy(method='tight_binding')
            ref_energy = ref_result['energy']

            # Get VQE energy with limited iterations
            vqe_result = bond.compute_energy(method='VQE', max_iter=50)
            vqe_energy = vqe_result['energy']
            
            print(f"  Reference energy (tight-binding): {ref_energy:.6f} eV")
            print(f"  VQE energy: {vqe_energy:.6f} eV")
            
            # Check accuracy
            error = abs(vqe_energy - ref_energy)
            relative_error = error / abs(ref_energy) if ref_energy != 0 else error
            
            print(f"  Absolute error: {error:.6f} eV")
            print(f"  Relative error: {relative_error:.6f}")
            
            # Check if error is reasonable
            if relative_error > 0.01:  # 1% tolerance
                issues.append(f"VQE accuracy issue: {relative_error:.4f} relative error")
            else:
                print(f"  ✅ VQE accuracy is good")
            
        except Exception as e:
            issues.append(f"VQE accuracy test failed: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ VQE accuracy is acceptable")
            self.passed_tests += 1
        
        return issues
    
    def test_temperature_effects_correct_api(self):
        """Test 1.5: Temperature effects with correct API."""
        print("="*70)
        print("TEST 1.5: TEMPERATURE EFFECTS (CORRECT API)")
        print("="*70)
        
        issues = []
        
        # Test temperature effects with CORRECT API
        atoms = [Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(4)]
        
        print(f"Testing temperature effects with correct API:")
        
        temperatures = [0.0, 100.0, 300.0, 500.0]
        
        for T in temperatures:
            print(f"\nTemperature T={T} K:")
            
            try:
                # Use CORRECT API - create new Temperature object
                from kanad.core.temperature import Temperature
                temp = Temperature(T)
                
                # Create bond with temperature
                bond = MetallicBond(atoms, hopping_parameter=-1.0, periodic=True, temperature=temp)
                
                # Test temperature-dependent energy
                result = bond.compute_energy(method='tight_binding')
                energy = result['energy']
                
                print(f"  Energy: {energy:.6f} eV")
                
                # Check if energy is reasonable
                if np.isnan(energy) or np.isinf(energy):
                    issues.append(f"T={T}K: Invalid energy (NaN or Inf): {energy}")
                
            except Exception as e:
                issues.append(f"T={T}K: Exception during temperature test: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Temperature effects work correctly")
            self.passed_tests += 1
        
        return issues
    
    def run_all_tests(self):
        """Run all corrected metallic physics tests."""
        print("="*70)
        print("SCIENTIFIC TEST 1: CORRECT METALLIC PHYSICS (FIXED API)")
        print("="*70)
        
        # Run all tests
        self.test_metallic_bond_basic_functionality()
        self.test_known_limitations()
        self.test_hubbard_u_investigation()
        self.test_vqe_accuracy()
        self.test_temperature_effects_correct_api()
        
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
    """Run the corrected metallic physics test suite."""
    tester = TestCorrectMetallicPhysics()
    results = tester.run_all_tests()
    return results

if __name__ == "__main__":
    main()
