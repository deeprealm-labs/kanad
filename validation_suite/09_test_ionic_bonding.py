#!/usr/bin/env python3
"""
Scientific Test 9: Ionic Bonding and Crystal Structures

This test rigorously validates ionic bonding implementations
with crystal structures and ionic compounds.

Scientific Focus:
- Ionic bond formation
- Crystal structure geometry
- Madelung constants
- Lattice energy calculations
- Ionic radii and coordination
- Crystal field effects
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
from kanad.bonds.ionic_bond import IonicBond
from kanad.core.hamiltonians.ionic_hamiltonian import IonicHamiltonian

class TestIonicBonding:
    """Comprehensive test of ionic bonding and crystal structures."""
    
    def __init__(self):
        self.test_results = {}
        self.critical_issues = []
        self.warnings = []
        self.passed_tests = 0
        self.failed_tests = 0
        
    def test_simple_ionic_compounds(self):
        """Test 9.1: Simple ionic compounds."""
        print("="*70)
        print("TEST 9.1: SIMPLE IONIC COMPOUNDS")
        print("="*70)
        
        issues = []
        
        # Test simple ionic compounds
        test_compounds = [
            # NaCl (rock salt structure)
            ([Atom('Na', position=np.array([0.0, 0.0, 0.0])),
              Atom('Cl', position=np.array([2.82, 0.0, 0.0]))], "NaCl"),
            
            # CsCl (cesium chloride structure)
            ([Atom('Cs', position=np.array([0.0, 0.0, 0.0])),
              Atom('Cl', position=np.array([3.57, 0.0, 0.0]))], "CsCl"),
            
            # CaF2 (fluorite structure)
            ([Atom('Ca', position=np.array([0.0, 0.0, 0.0])),
              Atom('F', position=np.array([2.36, 0.0, 0.0])),
              Atom('F', position=np.array([-2.36, 0.0, 0.0]))], "CaF2"),
        ]
        
        print(f"Testing simple ionic compounds:")
        
        for atoms, name in test_compounds:
            print(f"\nCompound: {name}")
            
            try:
                bond = IonicBond(atoms)
                
                # Get basic properties
                result = bond.compute_energy(method='ionic')
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
                issues.append(f"{name}: Exception during ionic bonding: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Simple ionic compounds work")
            self.passed_tests += 1
        
        return issues
    
    def test_madelung_constants(self):
        """Test 9.2: Madelung constants calculation."""
        print("="*70)
        print("TEST 9.2: MADELUNG CONSTANTS CALCULATION")
        print("="*70)
        
        issues = []
        
        # Test Madelung constants for different crystal structures
        test_structures = [
            # NaCl (rock salt) - Madelung constant ≈ 1.748
            ([Atom('Na', position=np.array([0.0, 0.0, 0.0])),
              Atom('Cl', position=np.array([2.82, 0.0, 0.0]))], "NaCl", 1.748),
            
            # CsCl (cesium chloride) - Madelung constant ≈ 1.763
            ([Atom('Cs', position=np.array([0.0, 0.0, 0.0])),
              Atom('Cl', position=np.array([3.57, 0.0, 0.0]))], "CsCl", 1.763),
            
            # ZnS (zinc blende) - Madelung constant ≈ 1.638
            ([Atom('Zn', position=np.array([0.0, 0.0, 0.0])),
              Atom('S', position=np.array([2.34, 0.0, 0.0]))], "ZnS", 1.638),
        ]
        
        print(f"Testing Madelung constants calculation:")
        
        for atoms, name, expected_madelung in test_structures:
            print(f"\nStructure: {name}")
            
            try:
                bond = IonicBond(atoms)
                
                # Get Madelung constant
                madelung = bond.get_madelung_constant()
                
                print(f"  Computed Madelung constant: {madelung:.6f}")
                print(f"  Expected Madelung constant: {expected_madelung:.6f}")
                print(f"  Error: {abs(madelung - expected_madelung):.6f}")
                
                # Check if Madelung constant is reasonable
                if np.isnan(madelung) or np.isinf(madelung):
                    issues.append(f"{name}: Invalid Madelung constant (NaN or Inf): {madelung}")
                
                # Check if Madelung constant is close to expected value
                if abs(madelung - expected_madelung) > 0.1:  # 0.1 tolerance
                    issues.append(f"{name}: Madelung constant incorrect: {madelung:.6f} vs {expected_madelung:.6f}")
                
            except Exception as e:
                issues.append(f"{name}: Exception during Madelung constant calculation: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Madelung constants calculation is correct")
            self.passed_tests += 1
        
        return issues
    
    def test_lattice_energy_calculations(self):
        """Test 9.3: Lattice energy calculations."""
        print("="*70)
        print("TEST 9.3: LATTICE ENERGY CALCULATIONS")
        print("="*70)
        
        issues = []
        
        # Test lattice energy calculations
        test_compounds = [
            # NaCl - lattice energy ≈ -787 kJ/mol
            ([Atom('Na', position=np.array([0.0, 0.0, 0.0])),
              Atom('Cl', position=np.array([2.82, 0.0, 0.0]))], "NaCl", -787),
            
            # CsCl - lattice energy ≈ -657 kJ/mol
            ([Atom('Cs', position=np.array([0.0, 0.0, 0.0])),
              Atom('Cl', position=np.array([3.57, 0.0, 0.0]))], "CsCl", -657),
            
            # CaF2 - lattice energy ≈ -2633 kJ/mol
            ([Atom('Ca', position=np.array([0.0, 0.0, 0.0])),
              Atom('F', position=np.array([2.36, 0.0, 0.0])),
              Atom('F', position=np.array([-2.36, 0.0, 0.0]))], "CaF2", -2633),
        ]
        
        print(f"Testing lattice energy calculations:")
        
        for atoms, name, expected_lattice_energy in test_compounds:
            print(f"\nCompound: {name}")
            
            try:
                bond = IonicBond(atoms)
                
                # Get lattice energy
                lattice_energy = bond.get_lattice_energy()
                
                print(f"  Computed lattice energy: {lattice_energy:.6f} eV")
                print(f"  Expected lattice energy: {expected_lattice_energy:.6f} eV")
                print(f"  Error: {abs(lattice_energy - expected_lattice_energy):.6f} eV")
                
                # Check if lattice energy is reasonable
                if np.isnan(lattice_energy) or np.isinf(lattice_energy):
                    issues.append(f"{name}: Invalid lattice energy (NaN or Inf): {lattice_energy}")
                
                # Check if lattice energy is negative (exothermic)
                if lattice_energy > 0:
                    issues.append(f"{name}: Positive lattice energy (endothermic): {lattice_energy}")
                
            except Exception as e:
                issues.append(f"{name}: Exception during lattice energy calculation: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Lattice energy calculations work")
            self.passed_tests += 1
        
        return issues
    
    def test_ionic_radii_and_coordination(self):
        """Test 9.4: Ionic radii and coordination."""
        print("="*70)
        print("TEST 9.4: IONIC RADII AND COORDINATION")
        print("="*70)
        
        issues = []
        
        # Test ionic radii and coordination
        test_compounds = [
            # NaCl - Na+ (6-coordinate), Cl- (6-coordinate)
            ([Atom('Na', position=np.array([0.0, 0.0, 0.0])),
              Atom('Cl', position=np.array([2.82, 0.0, 0.0]))], "NaCl"),
            
            # CsCl - Cs+ (8-coordinate), Cl- (8-coordinate)
            ([Atom('Cs', position=np.array([0.0, 0.0, 0.0])),
              Atom('Cl', position=np.array([3.57, 0.0, 0.0]))], "CsCl"),
            
            # CaF2 - Ca2+ (8-coordinate), F- (4-coordinate)
            ([Atom('Ca', position=np.array([0.0, 0.0, 0.0])),
              Atom('F', position=np.array([2.36, 0.0, 0.0])),
              Atom('F', position=np.array([-2.36, 0.0, 0.0]))], "CaF2"),
        ]
        
        print(f"Testing ionic radii and coordination:")
        
        for atoms, name in test_compounds:
            print(f"\nCompound: {name}")
            
            try:
                bond = IonicBond(atoms)
                
                # Get ionic radii
                ionic_radii = bond.get_ionic_radii()
                
                print(f"  Ionic radii: {ionic_radii}")
                
                # Check if ionic radii are reasonable
                if not isinstance(ionic_radii, dict):
                    issues.append(f"{name}: Ionic radii not a dictionary: {type(ionic_radii)}")
                
                for ion, radius in ionic_radii.items():
                    if radius <= 0 or radius > 5.0:  # Unreasonable ionic radius
                        issues.append(f"{name}: Unreasonable ionic radius {ion}: {radius}")
                
                # Get coordination numbers
                coordination = bond.get_coordination_numbers()
                
                print(f"  Coordination numbers: {coordination}")
                
                # Check if coordination numbers are reasonable
                if not isinstance(coordination, dict):
                    issues.append(f"{name}: Coordination numbers not a dictionary: {type(coordination)}")
                
                for ion, coord in coordination.items():
                    if coord <= 0 or coord > 12:  # Unreasonable coordination number
                        issues.append(f"{name}: Unreasonable coordination number {ion}: {coord}")
                
            except Exception as e:
                issues.append(f"{name}: Exception during ionic radii test: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Ionic radii and coordination work")
            self.passed_tests += 1
        
        return issues
    
    def test_crystal_field_effects(self):
        """Test 9.5: Crystal field effects."""
        print("="*70)
        print("TEST 9.5: CRYSTAL FIELD EFFECTS")
        print("="*70)
        
        issues = []
        
        # Test crystal field effects
        test_compounds = [
            # FeCl3 - Fe3+ with crystal field splitting
            ([Atom('Fe', position=np.array([0.0, 0.0, 0.0])),
              Atom('Cl', position=np.array([2.2, 0.0, 0.0])),
              Atom('Cl', position=np.array([-2.2, 0.0, 0.0])),
              Atom('Cl', position=np.array([0.0, 2.2, 0.0]))], "FeCl3"),
            
            # CuSO4 - Cu2+ with crystal field splitting
            ([Atom('Cu', position=np.array([0.0, 0.0, 0.0])),
              Atom('S', position=np.array([2.3, 0.0, 0.0])),
              Atom('O', position=np.array([-1.6, 1.6, 0.0])),
              Atom('O', position=np.array([-1.6, -1.6, 0.0])),
              Atom('O', position=np.array([1.6, 1.6, 0.0])),
              Atom('O', position=np.array([1.6, -1.6, 0.0]))], "CuSO4"),
        ]
        
        print(f"Testing crystal field effects:")
        
        for atoms, name in test_compounds:
            print(f"\nCompound: {name}")
            
            try:
                bond = IonicBond(atoms)
                
                # Get crystal field splitting
                crystal_field_splitting = bond.get_crystal_field_splitting()
                
                print(f"  Crystal field splitting: {crystal_field_splitting:.6f} eV")
                
                # Check if crystal field splitting is reasonable
                if np.isnan(crystal_field_splitting) or np.isinf(crystal_field_splitting):
                    issues.append(f"{name}: Invalid crystal field splitting (NaN or Inf): {crystal_field_splitting}")
                
                # Check if crystal field splitting is positive
                if crystal_field_splitting < 0:
                    issues.append(f"{name}: Negative crystal field splitting: {crystal_field_splitting}")
                
                # Get d-orbital energies
                d_orbital_energies = bond.get_d_orbital_energies()
                
                print(f"  d-orbital energies: {d_orbital_energies}")
                
                # Check if d-orbital energies are reasonable
                if not isinstance(d_orbital_energies, dict):
                    issues.append(f"{name}: d-orbital energies not a dictionary: {type(d_orbital_energies)}")
                
            except Exception as e:
                issues.append(f"{name}: Exception during crystal field test: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Crystal field effects work")
            self.passed_tests += 1
        
        return issues
    
    def test_ionic_edge_cases(self):
        """Test 9.6: Ionic bonding edge cases."""
        print("="*70)
        print("TEST 9.6: IONIC BONDING EDGE CASES")
        print("="*70)
        
        issues = []
        
        # Test edge cases
        edge_cases = [
            # Single ion
            ([Atom('Na', position=np.array([0.0, 0.0, 0.0]))], "Single Na+"),
            
            # Two ions far apart
            ([Atom('Na', position=np.array([0.0, 0.0, 0.0])),
              Atom('Cl', position=np.array([10.0, 0.0, 0.0]))], "Na+ and Cl- far apart"),
            
            # Very large crystal
            ([Atom('Na', position=np.array([i*2.82, 0.0, 0.0])) for i in range(10)] +
             [Atom('Cl', position=np.array([i*2.82 + 1.41, 0.0, 0.0])) for i in range(10)], "Na10Cl10"),
        ]
        
        print(f"Testing ionic bonding edge cases:")
        
        for atoms, name in edge_cases:
            print(f"\nCase: {name}")
            
            try:
                bond = IonicBond(atoms)
                
                # Get basic properties
                result = bond.compute_energy(method='ionic')
                energy = result['energy']
                
                print(f"  Atoms: {[atom.symbol for atom in atoms]}")
                print(f"  Energy: {energy:.6f} eV")
                
                # Check if energy is reasonable
                if np.isnan(energy) or np.isinf(energy):
                    issues.append(f"{name}: Invalid energy (NaN or Inf): {energy}")
                
            except Exception as e:
                issues.append(f"{name}: Exception during ionic edge case test: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Ionic bonding edge cases are handled correctly")
            self.passed_tests += 1
        
        return issues
    
    def run_all_tests(self):
        """Run all ionic bonding tests."""
        print("="*70)
        print("SCIENTIFIC TEST 9: IONIC BONDING AND CRYSTAL STRUCTURES")
        print("="*70)
        
        # Run all tests
        self.test_simple_ionic_compounds()
        self.test_madelung_constants()
        self.test_lattice_energy_calculations()
        self.test_ionic_radii_and_coordination()
        self.test_crystal_field_effects()
        self.test_ionic_edge_cases()
        
        # Print summary
        print("\n" + "="*70)
        print("TEST 9 SUMMARY")
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
    """Run the ionic bonding test suite."""
    tester = TestIonicBonding()
    results = tester.run_all_tests()
    return results

if __name__ == "__main__":
    main()
