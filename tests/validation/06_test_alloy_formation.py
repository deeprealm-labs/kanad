#!/usr/bin/env python3
"""
Scientific Test 6: Alloy Formation and Mixing Thermodynamics

This test rigorously validates the alloy formation and mixing
thermodynamics in the metallic bonding framework.

Scientific Focus:
- Alloy formation thermodynamics
- Mixing entropy calculations
- Mixing enthalpy calculations
- Free energy of mixing
- Alloy stability predictions
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

class TestAlloyFormation:
    """Comprehensive test of alloy formation and mixing thermodynamics."""
    
    def __init__(self):
        self.test_results = {}
        self.critical_issues = []
        self.warnings = []
        self.passed_tests = 0
        self.failed_tests = 0
        
    def test_alloy_formation_basic(self):
        """Test 6.1: Basic alloy formation functionality."""
        print("="*70)
        print("TEST 6.1: BASIC ALLOY FORMATION FUNCTIONALITY")
        print("="*70)
        
        issues = []
        
        # Test different alloy compositions
        alloy_compositions = [
            ([Atom('Na', position=np.array([0.0, 0, 0])), Atom('K', position=np.array([3.0, 0, 0]))], "Na-K"),
            ([Atom('Na', position=np.array([0.0, 0, 0])), Atom('Na', position=np.array([3.0, 0, 0])), 
             Atom('K', position=np.array([6.0, 0, 0])), Atom('K', position=np.array([9.0, 0, 0]))], "Na2K2"),
            ([Atom('Cu', position=np.array([0.0, 0, 0])), Atom('Zn', position=np.array([3.0, 0, 0]))], "Cu-Zn"),
        ]
        
        print(f"Testing basic alloy formation:")
        
        for atoms, name in alloy_compositions:
            print(f"\nAlloy: {name}")
            
            try:
                bond = MetallicBond(atoms, hopping_parameter=-1.0, periodic=True)
                
                # Get basic properties
                result = bond.compute_energy(method='tight_binding')
                energy = result['energy']
                
                print(f"  Atoms: {[atom.symbol for atom in atoms]}")
                print(f"  Energy: {energy:.6f} eV")
                
                # Check if energy is reasonable
                if np.isnan(energy) or np.isinf(energy):
                    issues.append(f"{name}: Invalid energy (NaN or Inf): {energy}")
                
                # Check if energy is negative (bound system)
                if energy > 0:
                    issues.append(f"{name}: Positive energy (unbound system): {energy}")
                
            except Exception as e:
                issues.append(f"{name}: Exception during alloy formation: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Basic alloy formation works")
            self.passed_tests += 1
        
        return issues
    
    def test_mixing_entropy_calculation(self):
        """Test 6.2: Mixing entropy calculation."""
        print("="*70)
        print("TEST 6.2: MIXING ENTROPY CALCULATION")
        print("="*70)
        
        issues = []
        
        # Test mixing entropy for different compositions
        test_compositions = [
            ([Atom('Na', position=np.array([0.0, 0, 0])), Atom('K', position=np.array([3.0, 0, 0]))], "Na-K", 0.5),
            ([Atom('Na', position=np.array([0.0, 0, 0])), Atom('Na', position=np.array([3.0, 0, 0])), 
             Atom('K', position=np.array([6.0, 0, 0])), Atom('K', position=np.array([9.0, 0, 0]))], "Na2K2", 0.5),
            ([Atom('Na', position=np.array([0.0, 0, 0])), Atom('Na', position=np.array([3.0, 0, 0])), 
             Atom('Na', position=np.array([6.0, 0, 0])), Atom('K', position=np.array([9.0, 0, 0]))], "Na3K", 0.25),
        ]
        
        print(f"Testing mixing entropy calculation:")
        
        for atoms, name, expected_x in test_compositions:
            print(f"\nAlloy: {name} (x={expected_x})")
            
            try:
                bond = MetallicBond(atoms, hopping_parameter=-1.0, periodic=True)
                
                # Get mixing entropy
                mixing_entropy = bond.get_mixing_entropy()
                
                print(f"  Mixing entropy: {mixing_entropy:.6f} eV/K")
                
                # Check mixing entropy properties
                # 1. Mixing entropy should be non-negative
                if mixing_entropy < 0:
                    issues.append(f"{name}: Negative mixing entropy: {mixing_entropy}")
                
                # 2. Mixing entropy should be finite
                if np.isnan(mixing_entropy) or np.isinf(mixing_entropy):
                    issues.append(f"{name}: Invalid mixing entropy (NaN or Inf): {mixing_entropy}")
                
                # 3. For ideal mixing, S_mix = -k_B * [x*ln(x) + (1-x)*ln(1-x)]
                # where x is the mole fraction of one component
                k_B = 8.617e-5  # eV/K
                expected_entropy = -k_B * (expected_x * np.log(expected_x) + (1 - expected_x) * np.log(1 - expected_x))
                
                print(f"  Expected entropy: {expected_entropy:.6f} eV/K")
                print(f"  Error: {abs(mixing_entropy - expected_entropy):.6f} eV/K")
                
                if abs(mixing_entropy - expected_entropy) > 1e-6:  # 1 μeV/K tolerance
                    issues.append(f"{name}: Mixing entropy incorrect: {mixing_entropy:.6f} vs {expected_entropy:.6f}")
                
            except Exception as e:
                issues.append(f"{name}: Exception during mixing entropy calculation: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Mixing entropy calculation is correct")
            self.passed_tests += 1
        
        return issues
    
    def test_mixing_enthalpy_calculation(self):
        """Test 6.3: Mixing enthalpy calculation."""
        print("="*70)
        print("TEST 6.3: MIXING ENTHALPY CALCULATION")
        print("="*70)
        
        issues = []
        
        # Test mixing enthalpy for different compositions
        test_compositions = [
            ([Atom('Na', position=np.array([0.0, 0, 0])), Atom('K', position=np.array([3.0, 0, 0]))], "Na-K"),
            ([Atom('Na', position=np.array([0.0, 0, 0])), Atom('Na', position=np.array([3.0, 0, 0])), 
             Atom('K', position=np.array([6.0, 0, 0])), Atom('K', position=np.array([9.0, 0, 0]))], "Na2K2"),
            ([Atom('Cu', position=np.array([0.0, 0, 0])), Atom('Zn', position=np.array([3.0, 0, 0]))], "Cu-Zn"),
        ]
        
        print(f"Testing mixing enthalpy calculation:")
        
        for atoms, name in test_compositions:
            print(f"\nAlloy: {name}")
            
            try:
                bond = MetallicBond(atoms, hopping_parameter=-1.0, periodic=True)
                
                # Get mixing enthalpy
                mixing_enthalpy = bond.get_mixing_enthalpy()
                
                print(f"  Mixing enthalpy: {mixing_enthalpy:.6f} eV")
                
                # Check mixing enthalpy properties
                # 1. Mixing enthalpy should be finite
                if np.isnan(mixing_enthalpy) or np.isinf(mixing_enthalpy):
                    issues.append(f"{name}: Invalid mixing enthalpy (NaN or Inf): {mixing_enthalpy}")
                
                # 2. Mixing enthalpy can be positive or negative
                # Positive: endothermic mixing (unfavorable)
                # Negative: exothermic mixing (favorable)
                
            except Exception as e:
                issues.append(f"{name}: Exception during mixing enthalpy calculation: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Mixing enthalpy calculation works")
            self.passed_tests += 1
        
        return issues
    
    def test_free_energy_of_mixing(self):
        """Test 6.4: Free energy of mixing calculation."""
        print("="*70)
        print("TEST 6.4: FREE ENERGY OF MIXING CALCULATION")
        print("="*70)
        
        issues = []
        
        # Test free energy of mixing for different compositions
        test_compositions = [
            ([Atom('Na', position=np.array([0.0, 0, 0])), Atom('K', position=np.array([3.0, 0, 0]))], "Na-K"),
            ([Atom('Na', position=np.array([0.0, 0, 0])), Atom('Na', position=np.array([3.0, 0, 0])), 
             Atom('K', position=np.array([6.0, 0, 0])), Atom('K', position=np.array([9.0, 0, 0]))], "Na2K2"),
            ([Atom('Cu', position=np.array([0.0, 0, 0])), Atom('Zn', position=np.array([3.0, 0, 0]))], "Cu-Zn"),
        ]
        
        print(f"Testing free energy of mixing calculation:")
        
        for atoms, name in test_compositions:
            print(f"\nAlloy: {name}")
            
            try:
                bond = MetallicBond(atoms, hopping_parameter=-1.0, periodic=True)
                
                # Get free energy of mixing
                free_energy_mixing = bond.get_free_energy_of_mixing()
                
                print(f"  Free energy of mixing: {free_energy_mixing:.6f} eV")
                
                # Check free energy of mixing properties
                # 1. Free energy of mixing should be finite
                if np.isnan(free_energy_mixing) or np.isinf(free_energy_mixing):
                    issues.append(f"{name}: Invalid free energy of mixing (NaN or Inf): {free_energy_mixing}")
                
                # 2. Free energy of mixing can be positive or negative
                # Negative: mixing is thermodynamically favorable
                # Positive: mixing is thermodynamically unfavorable
                
            except Exception as e:
                issues.append(f"{name}: Exception during free energy of mixing calculation: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Free energy of mixing calculation works")
            self.passed_tests += 1
        
        return issues
    
    def test_alloy_stability_predictions(self):
        """Test 6.5: Alloy stability predictions."""
        print("="*70)
        print("TEST 6.5: ALLOY STABILITY PREDICTIONS")
        print("="*70)
        
        issues = []
        
        # Test alloy stability for different compositions
        test_compositions = [
            ([Atom('Na', position=np.array([0.0, 0, 0])), Atom('K', position=np.array([3.0, 0, 0]))], "Na-K"),
            ([Atom('Na', position=np.array([0.0, 0, 0])), Atom('Na', position=np.array([3.0, 0, 0])), 
             Atom('K', position=np.array([6.0, 0, 0])), Atom('K', position=np.array([9.0, 0, 0]))], "Na2K2"),
            ([Atom('Cu', position=np.array([0.0, 0, 0])), Atom('Zn', position=np.array([3.0, 0, 0]))], "Cu-Zn"),
        ]
        
        print(f"Testing alloy stability predictions:")
        
        for atoms, name in test_compositions:
            print(f"\nAlloy: {name}")
            
            try:
                bond = MetallicBond(atoms, hopping_parameter=-1.0, periodic=True)
                
                # Get alloy stability
                is_stable = bond.is_alloy_stable()
                
                print(f"  Is stable: {is_stable}")
                
                # Check if stability prediction makes sense
                # This depends on the specific implementation
                # For now, just check that it returns a boolean
                if not isinstance(is_stable, bool):
                    issues.append(f"{name}: Stability prediction not boolean: {is_stable}")
                
            except Exception as e:
                issues.append(f"{name}: Exception during stability prediction: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Alloy stability predictions work")
            self.passed_tests += 1
        
        return issues
    
    def test_alloy_edge_cases(self):
        """Test 6.6: Alloy edge cases and numerical stability."""
        print("="*70)
        print("TEST 6.6: ALLOY EDGE CASES")
        print("="*70)
        
        issues = []
        
        # Test edge cases
        edge_cases = [
            ([Atom('Na', position=np.array([0.0, 0, 0]))], "Single atom"),
            ([Atom('Na', position=np.array([0.0, 0, 0])), Atom('Na', position=np.array([3.0, 0, 0]))], "Pure Na"),
            ([Atom('Na', position=np.array([0.0, 0, 0])), Atom('K', position=np.array([3.0, 0, 0])), 
             Atom('Cu', position=np.array([6.0, 0, 0])), Atom('Zn', position=np.array([9.0, 0, 0]))], "Multi-component"),
        ]
        
        print(f"Testing alloy edge cases:")
        
        for atoms, name in edge_cases:
            print(f"\nCase: {name}")
            
            try:
                bond = MetallicBond(atoms, hopping_parameter=-1.0, periodic=True)
                
                # Get basic properties
                result = bond.compute_energy(method='tight_binding')
                energy = result['energy']
                
                print(f"  Atoms: {[atom.symbol for atom in atoms]}")
                print(f"  Energy: {energy:.6f} eV")
                
                # Check if energy is reasonable
                if np.isnan(energy) or np.isinf(energy):
                    issues.append(f"{name}: Invalid energy (NaN or Inf): {energy}")
                
            except Exception as e:
                issues.append(f"{name}: Exception during edge case test: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Alloy edge cases are handled correctly")
            self.passed_tests += 1
        
        return issues
    
    def run_all_tests(self):
        """Run all alloy formation tests."""
        print("="*70)
        print("SCIENTIFIC TEST 6: ALLOY FORMATION AND MIXING THERMODYNAMICS")
        print("="*70)
        
        # Run all tests
        self.test_alloy_formation_basic()
        self.test_mixing_entropy_calculation()
        self.test_mixing_enthalpy_calculation()
        self.test_free_energy_of_mixing()
        self.test_alloy_stability_predictions()
        self.test_alloy_edge_cases()
        
        # Print summary
        print("\n" + "="*70)
        print("TEST 6 SUMMARY")
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
    """Run the alloy formation test suite."""
    tester = TestAlloyFormation()
    results = tester.run_all_tests()
    return results

if __name__ == "__main__":
    main()
