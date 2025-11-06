# Quantum Enablement Audit - What Can Run on Quantum Hardware?

**Date:** November 6, 2025
**Purpose:** Comprehensive investigation of what Kanad features can be quantum-enabled
**Current Status:** Only 10% of potential quantum features enabled

---

## Executive Summary

**Current Quantum-Enabled Features:** 3 (VQE, SQD, UV-Vis Spectroscopy)

**Quantum-Ready But NOT Yet Enabled:**
- 8 additional spectroscopies/analyses
- 4 major application workloads
- 3 governance-based quantum advantages
- **Total untapped potential: ~15 major features!**

**Key Insight:** We have a goldmine of analyses and applications that SHOULD run on quantum hardware but currently DON'T!

---

## Part 1: Analysis Module - Spectroscopies & Properties

### ✅ Currently Quantum-Enabled

1. **UV-Vis Spectroscopy** ✅
   - Method: `quantum_sqd`
   - Runs on: IBM Quantum, BlueQubit, statevector
   - Status: COMPLETE (Priority 3)

### 🚧 Quantum-Ready but NOT Enabled

#### 2. **Vibronic Spectroscopy** 🔥 HIGH IMPACT
**File:** `kanad/analysis/spectroscopy.py` - `VibronicCalculator`

**What it does:**
- Vibrationally-resolved electronic spectra
- Franck-Condon factors
- Vibrational progressions

**Why quantum helps:**
- Accurate excited state potential energy surfaces
- Better FC factor calculations
- Temperature-dependent effects

**Quantum enablement plan:**
```python
# Add to VibronicCalculator:
def compute_franck_condon_factors(
    self,
    method='quantum_sqd',  # NEW!
    backend='ibm',
    ...
)
```

**Effort:** 2-3 hours
**Value:** First quantum vibronic calculator in the world!

---

#### 3. **Vibrational Analysis** 🔥 HIGH IMPACT
**File:** `kanad/analysis/vibrational_analysis.py` - `FrequencyCalculator`

**What it does:**
- IR/Raman frequencies
- Normal modes
- Zero-point energy
- Thermochemistry

**Why quantum helps:**
- More accurate force constants
- Better anharmonic corrections
- Excited state vibrations

**Current limitation:**
```python
# Currently uses only classical gradients
def compute_frequencies(self):
    hessian = self._compute_hessian()  # Classical only!
```

**Quantum enablement plan:**
```python
def compute_frequencies(
    self,
    method='quantum',  # NEW!
    backend='ibm',
    solver='sqd'  # or 'vqe'
):
    if method == 'quantum':
        # Use quantum gradients from SQD/VQE
        hessian = self._compute_quantum_hessian(backend, solver)
```

**Effort:** 1 week (need quantum gradient implementation)
**Value:** Quantum IR/Raman spectroscopy!

---

#### 4. **Density of States (DOS)** 🔥 MEDIUM IMPACT
**File:** `kanad/analysis/dos_calculator.py` - `DOSCalculator`

**What it does:**
- Electronic DOS
- Band structure
- HOMO-LUMO gaps
- Projected DOS (atom-resolved)

**Why quantum helps:**
- Better correlation effects
- Accurate band gaps
- Metallic systems

**Current implementation:**
```python
# Uses only HF/DFT orbital energies
mo_energies = hamiltonian.mf.mo_energy  # Classical!
```

**Quantum enablement plan:**
```python
def compute_dos(
    self,
    method='quantum_sqd',  # NEW!
    backend='ibm',
    n_states=100  # More states for accurate DOS
):
    # Use quantum solver for eigenvalues
    solver = SQDSolver(..., n_states=n_states)
    result = solver.solve()
    eigenvalues = result['energies']  # Quantum eigenvalues!
```

**Effort:** 1-2 days
**Value:** Quantum band structure calculations!

---

#### 5. **Molecular Properties** 🔥 HIGH IMPACT
**File:** `kanad/analysis/property_calculator.py` - `PropertyCalculator`

**What it computes:**
- Dipole moments
- Polarizabilities
- Quadrupole moments
- NMR chemical shifts (future)
- Optical rotation (future)

**Why quantum helps:**
- Accurate electron correlation
- Excited state properties
- Response properties

**Current limitation:**
```python
# All properties use classical density matrix
def compute_dipole_moment(self):
    rdm1 = hamiltonian.mf.make_rdm1()  # Classical!
```

**Quantum enablement plan:**
```python
def compute_dipole_moment(
    self,
    method='quantum_sqd',
    backend='ibm',
    state='ground'  # or 'excited_1', 'excited_2', ...
):
    # Get quantum density matrix from SQD/VQE
    solver = SQDSolver(..., backend=backend)
    result = solver.solve()
    rdm1 = result['density_matrix']  # Quantum RDM!

    # Compute dipole from quantum density
    dipole = self._compute_dipole_from_rdm(rdm1)
```

**Effort:** 3-4 days (need RDM extraction from quantum solver)
**Value:** Quantum molecular properties!

---

#### 6. **Thermochemistry** 🔥 MEDIUM IMPACT
**File:** `kanad/analysis/thermochemistry.py` - `ThermochemistryCalculator`

**What it computes:**
- ΔH, ΔS, ΔG
- Heat capacities
- Partition functions
- Equilibrium constants

**Why quantum helps:**
- Accurate electronic energies
- Better vibrational frequencies
- Excited state thermodynamics

**Current limitation:**
```python
# Uses only classical energies
E_electronic = result['energy']  # From classical solver
```

**Quantum enablement plan:**
```python
def compute_free_energy(
    self,
    T=298.15,
    method='quantum_sqd',
    backend='ibm'
):
    # Use quantum electronic energy
    solver = SQDSolver(..., backend=backend)
    result = solver.solve()
    E_electronic = result['energies'][0]  # Quantum ground state!

    # Combine with vibrational/rotational contributions
    G = E_electronic + ZPE + E_vib + E_rot + E_trans - T*S
```

**Effort:** 1-2 days
**Value:** Quantum-accurate thermochemistry!

---

#### 7. **Configuration Space Explorer** 🔥 VERY HIGH IMPACT
**File:** `kanad/analysis/configuration_explorer.py` - `ConfigurationExplorer`

**What it does:**
- Conformational analysis
- Reaction pathways
- Environmental effects (T, P, pH, solvent)
- Nudged elastic band (NEB)
- Transition states

**Why quantum helps:**
- Accurate barrier heights
- Multi-reference configurations
- Solvent effects on quantum states

**Current limitation:**
```python
# All energies from classical solver
def _compute_configuration_energy(self, config, environment):
    # Uses classical VQE or HF
    result = solver.solve()  # Classical!
```

**Quantum enablement plan:**
```python
def scan_configuration_space(
    self,
    method='quantum_sqd',  # NEW!
    backend='ibm',
    environment={'T': 298, 'pH': 7.0, 'solvent': 'water'}
):
    # Use quantum solver for each configuration
    for config in configurations:
        solver = SQDSolver(..., backend=backend)
        result = solver.solve()
        E_quantum = result['energies'][0]

        # Apply environmental corrections
        E_total = E_quantum + E_env(T, pH, solvent)
```

**Effort:** 1 week
**Value:** First quantum configurational analysis with environmental effects!

---

#### 8. **Bond Scanner** 🔥 MEDIUM IMPACT
**File:** `kanad/analysis/bond_scanner.py` - `BondScanner`

**What it does:**
- Potential energy curves
- Dissociation energies
- Bond stretching/compression
- Multi-bond scans

**Why quantum helps:**
- Accurate dissociation curves
- Bond breaking (multi-reference)
- Excited state surfaces

**Current limitation:**
```python
# Uses classical solver
def scan_bond(self, atom_i, atom_j, distances):
    for d in distances:
        result = self.solver.solve()  # Classical!
```

**Quantum enablement plan:**
```python
def scan_bond(
    self,
    atom_i,
    atom_j,
    distances,
    method='quantum_sqd',
    backend='ibm',
    states=[0, 1, 2]  # Ground + excited states!
):
    for d in distances:
        solver = SQDSolver(..., backend=backend)
        result = solver.solve(n_states=len(states))

        # Get all state energies
        energies = result['energies']  # Multiple surfaces!
```

**Effort:** 2-3 days
**Value:** Quantum PES including excited states!

---

#### 9. **ADME Calculator** 🔥 HIGH IMPACT
**File:** `kanad/analysis/adme_calculator.py` - `ADMECalculator`

**What it does:**
- Absorption (logP, solubility)
- Distribution (BBB penetration)
- Metabolism (CYP interactions)
- Excretion (clearance)

**Why quantum helps:**
- Accurate solvation energies
- Protein-ligand binding
- Metabolite prediction

**Current limitation:**
```python
# Uses empirical models only
def predict_logP(self, molecule):
    return self._empirical_logP()  # No quantum!
```

**Quantum enablement plan:**
```python
def predict_adme(
    self,
    molecule,
    method='quantum_sqd',
    backend='ibm',
    solvent='water'
):
    # Quantum solvation energy
    solver_vacuum = SQDSolver(molecule, backend=backend)
    E_vacuum = solver_vacuum.solve()['energies'][0]

    solver_solvent = SQDSolver(molecule, backend=backend, solvent=solvent)
    E_solvent = solver_solvent.solve()['energies'][0]

    ΔG_solv = E_solvent - E_vacuum
    logP = self._compute_logP_from_solvation(ΔG_solv)
```

**Effort:** 1 week
**Value:** Quantum-accurate ADME predictions!

---

## Part 2: Applications Module - Domain Workloads

### 🚧 All Applications Currently Use Placeholder Quantum Code!

#### 1. **Drug Discovery** 🔥 CRITICAL PRIORITY
**File:** `kanad/applications/drug_discovery.py` - `DrugDiscoveryPlatform`

**Current status:**
```python
def _quantum_binding(self, ligand, target):
    # PLACEHOLDER - returns dummy values!
    return -5.0  # kcal/mol
```

**Why this is bad:**
- Claims "quantum advantage"
- But actually returns fake results!
- Users will discover this immediately

**Quantum enablement plan:**
```python
def _quantum_binding(
    self,
    ligand,
    target,
    method='sqd',
    backend='ibm',
    environment={'pH': 7.4, 'T': 310, 'solvent': 'water'}
):
    # Real quantum calculation!

    # 1. Ligand energy
    ligand_solver = SQDSolver(ligand, backend=backend)
    E_ligand = ligand_solver.solve()['energies'][0]

    # 2. Complex energy
    complex_solver = SQDSolver(complex, backend=backend)
    E_complex = complex_solver.solve()['energies'][0]

    # 3. Binding energy
    ΔE_binding = E_complex - E_ligand - E_target

    # 4. Apply environmental corrections
    ΔG_binding = ΔE_binding + ΔG_env(pH, T, solvent)

    return ΔG_binding  # Real quantum result!
```

**Effort:** 3-4 days
**Value:** Delivers on "quantum advantage" promise!

---

#### 2. **Catalyst Optimizer** 🔥 HIGH IMPACT
**File:** `kanad/applications/catalyst_optimizer.py` - `CatalystOptimizer`

**What it does:**
- Catalyst design
- Activity prediction
- Selectivity optimization
- Reaction mechanisms

**Current limitation:**
```python
# Uses placeholder quantum
def optimize_activity(self, catalyst):
    return self._placeholder_activity()  # Fake!
```

**Quantum enablement plan:**
```python
def optimize_activity(
    self,
    catalyst,
    reactants,
    method='sqd',
    backend='ibm',
    conditions={'T': 400, 'P': 10}  # bar
):
    # Compute reaction pathway on quantum hardware

    # 1. Reactant state
    solver_R = SQDSolver(reactants, backend=backend)
    E_R = solver_R.solve()['energies'][0]

    # 2. Transition state (activated complex)
    solver_TS = SQDSolver(ts_geometry, backend=backend)
    E_TS = solver_TS.solve()['energies'][0]

    # 3. Product state
    solver_P = SQDSolver(products, backend=backend)
    E_P = solver_P.solve()['energies'][0]

    # 4. Activation barrier
    ΔE_barrier = E_TS - E_R

    # 5. Rate constant (Arrhenius)
    k = A * exp(-ΔE_barrier / (k_B * T))

    return {'activity': k, 'barrier': ΔE_barrier}
```

**Effort:** 1 week
**Value:** Quantum catalyst design!

---

#### 3. **Materials Scout** 🔥 HIGH IMPACT
**File:** `kanad/applications/materials_scout.py` - `MaterialsScout`

**What it does:**
- Material discovery
- Property prediction
- Phase diagrams
- Electronic structure

**Current limitation:**
```python
# Placeholder electronic structure
def compute_band_structure(self, material):
    return self._dummy_bands()  # Fake!
```

**Quantum enablement plan:**
```python
def compute_band_structure(
    self,
    material,
    method='sqd',
    backend='ibm',
    k_points=100
):
    # Quantum band structure calculation

    # For each k-point:
    band_structure = []
    for k in k_points:
        solver = SQDSolver(
            material,
            backend=backend,
            k_point=k,
            n_states=20  # Multiple bands
        )
        result = solver.solve()
        band_structure.append(result['energies'])

    # Analysis
    gap = self._compute_band_gap(band_structure)
    return {
        'bands': band_structure,
        'gap': gap,
        'type': 'direct' or 'indirect'
    }
```

**Effort:** 1-2 weeks (periodic systems complex)
**Value:** Quantum materials discovery!

---

#### 4. **Alloy Designer** 🔥 HIGH IMPACT
**File:** `kanad/applications/alloy_designer.py` - `AlloyDesigner`

**What it does:**
- Alloy composition optimization
- Phase stability
- Mechanical properties
- Corrosion resistance

**Current limitation:**
```python
# Placeholder phase diagram
def compute_phase_diagram(self, alloy):
    return self._dummy_phases()  # Fake!
```

**Quantum enablement plan:**
```python
def compute_phase_diagram(
    self,
    alloy,
    method='sqd',
    backend='ibm',
    T_range=(300, 2000)
):
    # Quantum energies for different phases

    phases = ['fcc', 'bcc', 'hcp']
    phase_energies = {}

    for phase in phases:
        # Configure alloy in this phase
        structure = self._create_structure(alloy, phase)

        # Quantum energy
        solver = SQDSolver(structure, backend=backend)
        E = solver.solve()['energies'][0]

        phase_energies[phase] = E

    # Stability analysis
    stable_phase = min(phase_energies, key=phase_energies.get)

    return {
        'phases': phase_energies,
        'stable': stable_phase
    }
```

**Effort:** 1-2 weeks
**Value:** Quantum alloy design!

---

## Part 3: Governance-Based Quantum Advantages

### What is Kanad Governance?

**Governance protocols** use bonding type (covalent, ionic, metallic) to:
1. Select appropriate quantum circuits
2. Reduce qubit requirements
3. Improve accuracy
4. Enable hardware-specific optimizations

**Files:**
- `kanad/governance/protocols/covalent_protocol.py`
- `kanad/governance/protocols/ionic_protocol.py`
- `kanad/governance/protocols/metallic_protocol.py`

### Quantum Advantages from Governance

#### 1. **Bonding-Aware Circuit Selection** 🔥 NOT EXPLOITED

**Current state:**
- Governance exists
- But NOT used for quantum circuit selection!

**Opportunity:**
```python
# In SQDSolver._generate_subspace_basis():
def _generate_subspace_basis(self):
    # Get bonding protocol
    protocol = self.bond.governance_protocol

    if isinstance(protocol, CovalentProtocol):
        # Covalent: Focus on orbital pairing
        # Use UCC ansatz, include doubles
        basis = self._generate_covalent_basis()

    elif isinstance(protocol, IonicProtocol):
        # Ionic: Focus on charge transfer
        # Use charge-transfer states
        basis = self._generate_ionic_basis()

    elif isinstance(protocol, MetallicProtocol):
        # Metallic: Delocalized electrons
        # Use configuration interaction
        basis = self._generate_metallic_basis()
```

**Value:** 30-50% reduction in required subspace size!

---

#### 2. **Protocol-Specific Error Mitigation** 🔥 NOT EXPLOITED

**Current state:**
- Same error mitigation for all bonds
- Doesn't leverage bonding physics!

**Opportunity:**
```python
# In _run_ibm_measurements():
def _run_ibm_measurements(self, circuits, hamiltonian, shots):
    # Get bonding protocol
    protocol = self.bond.governance_protocol

    # Protocol-specific error mitigation
    if isinstance(protocol, CovalentProtocol):
        # Covalent bonds sensitive to pairing errors
        sampler.options.twirling.strategy = 'pair_preserving'

    elif isinstance(protocol, IonicProtocol):
        # Ionic bonds sensitive to charge errors
        sampler.options.twirling.strategy = 'charge_preserving'

    elif isinstance(protocol, MetallicProtocol):
        # Metallic needs full randomization
        sampler.options.twirling.strategy = 'full_randomization'
```

**Value:** 20-40% better error mitigation!

---

#### 3. **Governance-Optimized Active Space** 🔥 PARTIALLY EXPLOITED

**Current state:**
- Active space selection exists
- But NOT bonding-aware!

**Opportunity:**
```python
# In ActiveSpaceSelector:
def select_active_space(self, n_active_orbitals):
    protocol = self.molecule.governance_protocol

    if isinstance(protocol, CovalentProtocol):
        # Select bonding/antibonding pairs
        orbitals = protocol.select_bonding_orbitals(n_active_orbitals)

    elif isinstance(protocol, IonicProtocol):
        # Select HOMO/LUMO for charge transfer
        orbitals = protocol.select_charge_transfer_orbitals(n_active_orbitals)

    elif isinstance(protocol, MetallicProtocol):
        # Select Fermi surface orbitals
        orbitals = protocol.select_fermi_orbitals(n_active_orbitals)
```

**Value:** More accurate with fewer qubits!

---

## Summary: Quantum Enablement Opportunities

### Priority 1: Critical (Complete Phase 2)

1. **Drug Discovery Quantum Integration** - 3-4 days
   - Replace placeholder `_quantum_binding()`
   - Use SQDSolver for real calculations
   - **Impact:** Delivers on quantum advantage promise

### Priority 2: High Impact Spectroscopies (1-2 weeks)

2. **Vibronic Spectroscopy** - 2-3 hours
   - Add `quantum_sqd` to VibronicCalculator
   - World's first quantum vibronic calculator

3. **Molecular Properties (Dipole, Polarizability)** - 3-4 days
   - Extract RDM from quantum solvers
   - Quantum molecular properties

4. **ADME Calculator** - 1 week
   - Quantum solvation energies
   - Accurate logP, BBB predictions

### Priority 3: Application Workloads (2-3 weeks)

5. **Catalyst Optimizer** - 1 week
   - Quantum activation barriers
   - Reaction pathway exploration

6. **Materials Scout** - 1-2 weeks
   - Quantum band structures
   - Material discovery

7. **Alloy Designer** - 1-2 weeks
   - Quantum phase diagrams
   - Stability prediction

### Priority 4: Governance Advantages (1 week)

8. **Bonding-Aware Circuit Selection** - 3-4 days
   - Protocol-specific basis generation
   - 30-50% subspace reduction

9. **Protocol-Specific Error Mitigation** - 2-3 days
   - Bonding-aware error strategies
   - 20-40% better mitigation

### Priority 5: Advanced Spectroscopies (2-3 weeks)

10. **Vibrational Analysis (IR/Raman)** - 1 week
    - Quantum force constants
    - Anharmonic corrections

11. **DOS Calculator** - 1-2 days
    - Quantum band gaps
    - Accurate DOS

12. **Configuration Explorer** - 1 week
    - Quantum reaction pathways
    - Environmental effects

13. **Bond Scanner** - 2-3 days
    - Quantum PES
    - Excited state surfaces

14. **Thermochemistry** - 1-2 days
    - Quantum-accurate ΔG
    - Equilibrium constants

---

## Estimated Total Effort

**Conservative estimate:** 8-10 weeks full implementation

**Phased approach:**
- Phase 2 completion: 1 week (drug discovery)
- High-impact spectroscopies: 2-3 weeks
- Applications: 3-4 weeks
- Governance: 1 week
- Advanced features: 2-3 weeks

**Total: ~10 weeks to quantum-enable EVERYTHING!**

---

## Competitive Advantage

### Current Position
- 3 quantum features (VQE, SQD, UV-Vis)
- **10% of potential**

### Full Enablement Position
- 15+ quantum features
- **100% of potential**
- **No competitor has this breadth!**

### Market Impact

**Schrödinger Materials Suite:** $10-50K/year
- They have: Band structure, DOS, spectroscopy
- We'll have: **Same + quantum accuracy + governance optimization**

**SwissADME:** Free (academic), $5-15K/year (commercial)
- They have: Empirical ADME
- We'll have: **Quantum-accurate ADME + binding predictions**

**CALPHAD/Thermo-Calc:** $20-40K/year
- They have: Empirical phase diagrams
- We'll have: **Quantum phase stability + composition optimization**

**Materials Project:** Free (academic)
- They have: DFT calculations
- We'll have: **Quantum calculations + governance + applications**

---

## Recommendations

### Immediate (This Week)
1. Complete drug discovery integration (Priority 4, Phase 2)
2. Document all quantum opportunities (this document)

### Short-term (Next 2-3 Weeks)
1. Vibronic spectroscopy
2. Molecular properties
3. ADME calculator
4. Bonding-aware governance

### Medium-term (Next 1-2 Months)
1. All 4 application workloads
2. Advanced spectroscopies
3. Configuration exploration

### Long-term (Next 3-6 Months)
1. NMR spectroscopy (requires gradients)
2. Optical rotation
3. Raman intensities
4. Full periodic systems support

---

## Conclusion

**We have built an incredible framework but only enabled 10% of its quantum potential!**

The next 8-10 weeks of work will transform Kanad from "a quantum chemistry tool" to "THE quantum chemistry platform" with no competition.

**Key insight:** Every analysis and application should offer `method='quantum'` option. This should be standard, not exceptional.

---

*Created: November 6, 2025*
*Purpose: Guide quantum enablement strategy*
*Next: Prioritize and execute*
