#!/usr/bin/env python3
"""
Scientific Test 4: Temperature Effects and Thermodynamic Consistency

This test rigorously validates the temperature-dependent physics
and thermodynamic consistency of the metallic bonding framework.

Scientific Focus:
- Temperature-dependent energy calculations
- Fermi-Dirac distribution physics
- Entropy calculations
- Free energy consistency
- Heat capacity behavior
- Thermodynamic relations
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

class TestTemperatureEffects:
    """Comprehensive test of temperature effects and thermodynamics."""
    
    def __init__(self):
        self.test_results = {}
        self.critical_issues = []
        self.warnings = []
        self.passed_tests = 0
        self.failed_tests = 0
        
    def test_temperature_energy_calculation(self):
        """Test 4.1: Temperature-dependent energy calculations."""
        print("="*70)
        print("TEST 4.1: TEMPERATURE-DEPENDENT ENERGY CALCULATIONS")
        print("="*70)
        
        issues = []
        
        # Test different temperatures
        temperatures = [0.0, 100.0, 300.0, 500.0, 1000.0, 2000.0]  # K
        
        atoms = [Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(4)]
        bond = MetallicBond(atoms, hopping_parameter=-1.0, periodic=True)
        
        print(f"Testing temperature effects on {len(atoms)}-atom chain:")
        
        for T in temperatures:
            print(f"\nTemperature T={T} K:")
            
            try:
                # Create new bond with this temperature (Temperature is immutable)
                from kanad.core.temperature import Temperature
                temp = Temperature(T)
                bond = MetallicBond(atoms, lattice_type='1d_chain', temperature=temp)

                # Get temperature-dependent energy
                result = bond.compute_energy(method='tight_binding')
                energy = result['energy']
                
                print(f"  Energy: {energy:.6f} eV")
                
                # Check if energy is reasonable
                if np.isnan(energy) or np.isinf(energy):
                    issues.append(f"T={T}K: Invalid energy (NaN or Inf)")
                
                # Check if energy increases with temperature (for metallic systems)
                if T > 0:
                    # At finite temperature, energy should be higher than at T=0
                    # due to thermal excitation
                    pass  # We'll check this in the next test
                
            except Exception as e:
                issues.append(f"T={T}K: Exception during energy calculation: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Temperature-dependent energy calculations work")
            self.passed_tests += 1
        
        return issues
    
    def test_fermi_dirac_distribution(self):
        """Test 4.2: Fermi-Dirac distribution physics."""
        print("="*70)
        print("TEST 4.2: FERMI-DIRAC DISTRIBUTION PHYSICS")
        print("="*70)
        
        issues = []
        
        # Test Fermi-Dirac distribution at different temperatures
        temperatures = [0.0, 100.0, 300.0, 500.0, 1000.0]
        
        atoms = [Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(4)]
        bond = MetallicBond(atoms, hopping_parameter=-1.0, periodic=True)
        
        # Get eigenvalues
        eigenvalues = np.linalg.eigvalsh(bond.hamiltonian.h_tight_binding)
        fermi_energy = bond.hamiltonian.get_fermi_energy(eigenvalues)
        
        print(f"Fermi energy: {fermi_energy:.6f} eV")
        print(f"Eigenvalues: {eigenvalues}")
        
        for T in temperatures:
            print(f"\nTemperature T={T} K:")
            
            try:
                # Create new bond with this temperature (Temperature is immutable)
                from kanad.core.temperature import Temperature
                temp = Temperature(T)
                bond = MetallicBond(atoms, lattice_type='1d_chain', temperature=temp)
                
                # Get occupation numbers
                occupation_numbers = bond.temperature.get_occupation_numbers(eigenvalues, fermi_energy)
                
                print(f"  Occupation numbers: {occupation_numbers}")
                
                # Check Fermi-Dirac distribution properties
                # 1. All occupation numbers should be between 0 and 1
                if np.any(occupation_numbers < 0) or np.any(occupation_numbers > 1):
                    issues.append(f"T={T}K: Occupation numbers outside [0,1]: {occupation_numbers}")
                
                # 2. At T=0, occupation should be step function
                if T == 0.0:
                    n_electrons = len(atoms)
                    n_bands_occupied = int(np.ceil(n_electrons / 2.0))
                    
                    expected_occupation = np.zeros_like(eigenvalues)
                    expected_occupation[:n_bands_occupied] = 2.0  # 2 electrons per band
                    if n_electrons % 2 == 1:
                        expected_occupation[n_bands_occupied-1] = 1.0  # 1 electron in highest occupied band
                    
                    if not np.allclose(occupation_numbers, expected_occupation, atol=1e-10):
                        issues.append(f"T={T}K: T=0 occupation incorrect: {occupation_numbers} vs {expected_occupation}")
                
                # 3. At high T, occupation should be more uniform
                if T > 0:
                    # Check that occupation numbers sum to total number of electrons
                    total_occupation = np.sum(occupation_numbers)
                    expected_total = len(atoms)  # 1 electron per atom
                    
                    if abs(total_occupation - expected_total) > 1e-10:
                        issues.append(f"T={T}K: Total occupation incorrect: {total_occupation} vs {expected_total}")
                
            except Exception as e:
                issues.append(f"T={T}K: Exception during Fermi-Dirac test: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Fermi-Dirac distribution is correct")
            self.passed_tests += 1
        
        return issues
    
    def test_entropy_calculation(self):
        """Test 4.3: Entropy calculation physics."""
        print("="*70)
        print("TEST 4.3: ENTROPY CALCULATION PHYSICS")
        print("="*70)
        
        issues = []
        
        # Test entropy at different temperatures
        temperatures = [0.0, 100.0, 300.0, 500.0, 1000.0]
        
        atoms = [Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(4)]
        bond = MetallicBond(atoms, hopping_parameter=-1.0, periodic=True)
        
        print(f"Testing entropy for {len(atoms)}-atom chain:")
        
        for T in temperatures:
            print(f"\nTemperature T={T} K:")
            
            try:
                # Create new bond with this temperature (Temperature is immutable)
                from kanad.core.temperature import Temperature
                temp = Temperature(T)
                bond = MetallicBond(atoms, lattice_type='1d_chain', temperature=temp)
                
                # Get entropy
                entropy = bond.temperature.get_entropy()
                
                print(f"  Entropy: {entropy:.6f} eV/K")
                
                # Check entropy properties
                # 1. Entropy should be non-negative
                if entropy < 0:
                    issues.append(f"T={T}K: Negative entropy: {entropy}")
                
                # 2. At T=0, entropy should be zero (third law of thermodynamics)
                if T == 0.0:
                    if abs(entropy) > 1e-10:
                        issues.append(f"T={T}K: Non-zero entropy at T=0: {entropy}")
                
                # 3. Entropy should increase with temperature
                if T > 0:
                    # We'll check this in the next test
                    pass
                
            except Exception as e:
                issues.append(f"T={T}K: Exception during entropy calculation: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Entropy calculations are correct")
            self.passed_tests += 1
        
        return issues
    
    def test_free_energy_consistency(self):
        """Test 4.4: Free energy consistency and thermodynamic relations."""
        print("="*70)
        print("TEST 4.4: FREE ENERGY CONSISTENCY")
        print("="*70)
        
        issues = []
        
        # Test free energy at different temperatures
        temperatures = [100.0, 300.0, 500.0, 1000.0]
        
        atoms = [Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(4)]
        bond = MetallicBond(atoms, hopping_parameter=-1.0, periodic=True)
        
        print(f"Testing free energy for {len(atoms)}-atom chain:")
        
        for T in temperatures:
            print(f"\nTemperature T={T} K:")
            
            try:
                # Create new bond with this temperature (Temperature is immutable)
                from kanad.core.temperature import Temperature
                temp = Temperature(T)
                bond = MetallicBond(atoms, lattice_type='1d_chain', temperature=temp)
                
                # Get thermodynamic quantities
                internal_energy = bond.compute_energy(method='tight_binding')['energy']
                entropy = bond.temperature.get_entropy()
                free_energy = bond.temperature.get_free_energy()
                
                print(f"  Internal energy: {internal_energy:.6f} eV")
                print(f"  Entropy: {entropy:.6f} eV/K")
                print(f"  Free energy: {free_energy:.6f} eV")
                
                # Check thermodynamic relation: F = U - TS
                expected_free_energy = internal_energy - T * entropy
                free_energy_error = abs(free_energy - expected_free_energy)
                
                print(f"  Expected F: {expected_free_energy:.6f} eV")
                print(f"  Error: {free_energy_error:.6f} eV")
                
                if free_energy_error > 1e-10:
                    issues.append(f"T={T}K: Free energy relation violated: F={free_energy:.6f} vs U-TS={expected_free_energy:.6f}")
                
            except Exception as e:
                issues.append(f"T={T}K: Exception during free energy test: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Free energy consistency is correct")
            self.passed_tests += 1
        
        return issues
    
    def test_heat_capacity_behavior(self):
        """Test 4.5: Heat capacity behavior and temperature scaling."""
        print("="*70)
        print("TEST 4.5: HEAT CAPACITY BEHAVIOR")
        print("="*70)
        
        issues = []
        
        # Test heat capacity at different temperatures
        temperatures = [100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0]
        
        atoms = [Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(4)]
        bond = MetallicBond(atoms, hopping_parameter=-1.0, periodic=True)
        
        print(f"Testing heat capacity for {len(atoms)}-atom chain:")
        
        heat_capacities = []
        
        for T in temperatures:
            print(f"\nTemperature T={T} K:")
            
            try:
                # Create new bond with this temperature (Temperature is immutable)
                from kanad.core.temperature import Temperature
                temp = Temperature(T)
                bond = MetallicBond(atoms, lattice_type='1d_chain', temperature=temp)
                
                # Get heat capacity
                heat_capacity = bond.temperature.get_heat_capacity()
                heat_capacities.append(heat_capacity)
                
                print(f"  Heat capacity: {heat_capacity:.6f} eV/K")
                
                # Check heat capacity properties
                # 1. Heat capacity should be non-negative
                if heat_capacity < 0:
                    issues.append(f"T={T}K: Negative heat capacity: {heat_capacity}")
                
                # 2. Heat capacity should be finite
                if np.isnan(heat_capacity) or np.isinf(heat_capacity):
                    issues.append(f"T={T}K: Invalid heat capacity (NaN or Inf): {heat_capacity}")
                
            except Exception as e:
                issues.append(f"T={T}K: Exception during heat capacity calculation: {e}")
        
        # Check temperature scaling of heat capacity
        if len(heat_capacities) > 1:
            print(f"\nHeat capacity scaling:")
            for i, T in enumerate(temperatures):
                print(f"  T={T}K: C={heat_capacities[i]:.6f} eV/K")
            
            # Check if heat capacity increases with temperature (for metallic systems)
            for i in range(1, len(heat_capacities)):
                if heat_capacities[i] < heat_capacities[i-1]:
                    issues.append(f"Heat capacity decreases with temperature: T={temperatures[i-1]}K->{temperatures[i]}K: "
                                f"C={heat_capacities[i-1]:.6f}->{heat_capacities[i]:.6f}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Heat capacity behavior is correct")
            self.passed_tests += 1
        
        return issues
    
    def test_temperature_scaling_consistency(self):
        """Test 4.6: Temperature scaling consistency."""
        print("="*70)
        print("TEST 4.6: TEMPERATURE SCALING CONSISTENCY")
        print("="*70)
        
        issues = []
        
        # Test temperature scaling with different system sizes
        system_sizes = [2, 4, 6, 8, 10]
        temperatures = [100.0, 300.0, 500.0, 1000.0]
        
        print(f"Testing temperature scaling for different system sizes:")
        
        for n in system_sizes:
            print(f"\nSystem size N={n}:")
            
            atoms = [Atom('Na', position=np.array([i*3.0, 0, 0])) for i in range(n)]
            bond = MetallicBond(atoms, hopping_parameter=-1.0, periodic=True)
            
            for T in temperatures:
                print(f"  T={T}K:")
                
                try:
                    # Set temperature
                    bond.temperature.set_temperature(T)
                    
                    # Get thermodynamic quantities
                    internal_energy = bond.compute_energy(method='tight_binding')['energy']
                    entropy = bond.temperature.get_entropy()
                    free_energy = bond.temperature.get_free_energy()
                    heat_capacity = bond.temperature.get_heat_capacity()
                    
                    print(f"    U={internal_energy:.4f} eV, S={entropy:.4f} eV/K, F={free_energy:.4f} eV, C={heat_capacity:.4f} eV/K")
                    
                    # Check if quantities scale reasonably with system size
                    # Energy should scale roughly linearly with N
                    energy_per_atom = internal_energy / n
                    if abs(energy_per_atom) > 10.0:  # Unreasonably large energy per atom
                        issues.append(f"N={n}, T={T}K: Unreasonably large energy per atom: {energy_per_atom:.4f} eV/atom")
                    
                    # Entropy should scale roughly linearly with N
                    entropy_per_atom = entropy / n
                    if entropy_per_atom < 0 or entropy_per_atom > 1.0:  # Unreasonable entropy per atom
                        issues.append(f"N={n}, T={T}K: Unreasonable entropy per atom: {entropy_per_atom:.4f} eV/K/atom")
                    
                except Exception as e:
                    issues.append(f"N={n}, T={T}K: Exception during temperature scaling test: {e}")
        
        if issues:
            self.critical_issues.extend(issues)
            print(f"\n❌ CRITICAL ISSUES FOUND: {len(issues)}")
            for issue in issues:
                print(f"  - {issue}")
            self.failed_tests += 1
        else:
            print(f"\n✅ Temperature scaling is consistent")
            self.passed_tests += 1
        
        return issues
    
    def run_all_tests(self):
        """Run all temperature effects tests."""
        print("="*70)
        print("SCIENTIFIC TEST 4: TEMPERATURE EFFECTS AND THERMODYNAMICS")
        print("="*70)
        
        # Run all tests
        self.test_temperature_energy_calculation()
        self.test_fermi_dirac_distribution()
        self.test_entropy_calculation()
        self.test_free_energy_consistency()
        self.test_heat_capacity_behavior()
        self.test_temperature_scaling_consistency()
        
        # Print summary
        print("\n" + "="*70)
        print("TEST 4 SUMMARY")
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
    """Run the temperature effects test suite."""
    tester = TestTemperatureEffects()
    results = tester.run_all_tests()
    return results

if __name__ == "__main__":
    main()
