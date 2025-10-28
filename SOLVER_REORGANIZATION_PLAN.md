# Solver Reorganization Plan - Kanad Platform

## Current Problem

**User Confusion**: Current structure mixes ground state and excited state methods
- "VQE" sounds like it's only for ground states
- "Excited States" is a separate category but uses VQE internally
- Users can't easily configure VQE settings for excited states
- Hardcoded values (max_iterations=100) override user settings

**Current Structure** (Confusing):
```
Methods:
  ├── HF (Hartree-Fock)
  ├── VQE (Variational Quantum Eigensolver) ← Ground state only?
  ├── SQD (Subspace Quantum Diagonalization)
  └── EXCITED_STATES
       ├── CIS (classical)
       ├── TDDFT (classical)
       └── VQE (quantum) ← Wait, VQE again? How do I configure it?
```

---

## New Proposed Structure

### **Two Main Categories**:

#### 1. **Ground State Calculator**
*Calculate the lowest energy state of a molecule*

**Available Solvers**:
- **VQE** (Variational Quantum Eigensolver) - Quantum
- **HF** (Hartree-Fock) - Classical
- **QPE** (Quantum Phase Estimation) - Quantum [Future]
- **SQD** (Subspace Quantum Diagonalization) - Hybrid

**Configuration Panel** (changes based on selected solver):
```typescript
If VQE selected:
  - Ansatz: UCC, Hardware Efficient, etc.
  - Mapper: Jordan-Wigner, Parity, etc.
  - Optimizer: COBYLA, SLSQP, etc.
  - Max Iterations: [slider 10-1000]
  - Backend: Classical, BlueQubit, IBM

If HF selected:
  - SCF Convergence: [input]
  - Max SCF Iterations: [input]

If SQD selected:
  - N States: [input]
  - Solver: Davidson, Lanczos
  - Max Iterations: [input]
```

---

#### 2. **Excited State Calculator**
*Calculate higher energy electronic states*

**Available Solvers**:
- **VQE** (Orthogonally-Constrained VQE) - Quantum
- **CIS** (Configuration Interaction Singles) - Classical
- **TDDFT** (Time-Dependent DFT) - Classical
- **SQD** (Subspace Quantum Diagonalization) - Hybrid

**Common Settings** (all solvers):
- **N States**: Number of excited states to compute [1-10]

**Solver-Specific Configuration**:
```typescript
If VQE selected:
  └── VQE Configuration (same as Ground State)
       - Ansatz: UCC, Hardware Efficient
       - Mapper: Jordan-Wigner, Parity
       - Optimizer: COBYLA, SLSQP
       - Max Iterations: [slider 10-1000]
       - Penalty Weight: [slider 0.1-10.0]
       - Backend: Classical, BlueQubit, IBM

If CIS selected:
  - Max Singles: [input]
  - Davidson Convergence: [input]

If TDDFT selected:
  - Functional: B3LYP, PBE, etc.
  - TDA: [checkbox] Tamm-Dancoff Approximation

If SQD selected:
  - Solver: Davidson, Lanczos
  - Max Iterations: [input]
```

---

## Implementation Plan

### Phase 1: Backend Restructuring

#### 1.1 Create New API Schema

**File**: `api/routes/experiments.py`

```python
class MethodCategory(str, Enum):
    GROUND_STATE = "ground_state"
    EXCITED_STATE = "excited_state"

class GroundStateSolver(str, Enum):
    VQE = "vqe"
    HF = "hf"
    QPE = "qpe"
    SQD = "sqd"

class ExcitedStateSolver(str, Enum):
    VQE = "vqe"
    CIS = "cis"
    TDDFT = "tddft"
    SQD = "sqd"

class ExperimentConfig(BaseModel):
    # Main category
    category: MethodCategory

    # Ground State fields
    ground_solver: Optional[GroundStateSolver] = None
    vqe_config: Optional[VQEConfig] = None
    hf_config: Optional[HFConfig] = None
    sqd_config: Optional[SQDConfig] = None

    # Excited State fields
    excited_solver: Optional[ExcitedStateSolver] = None
    n_excited_states: Optional[int] = None
    excited_vqe_config: Optional[VQEConfig] = None  # Same as ground VQE but separate
    cis_config: Optional[CISConfig] = None
    tddft_config: Optional[TDDFTConfig] = None

class VQEConfig(BaseModel):
    ansatz: str = "hardware_efficient"
    mapper: str = "jordan_wigner"
    optimizer: str = "COBYLA"
    max_iterations: int = 100  # User-configurable, NOT hardcoded
    backend: str = "classical"
    backend_name: Optional[str] = None
    shots: Optional[int] = None
```

#### 1.2 Update Experiment Service

**File**: `api/services/experiment_service.py`

```python
async def execute_experiment(experiment_id: str, config: Dict):
    category = config.get('category')

    if category == 'ground_state':
        solver = config.get('ground_solver')

        if solver == 'vqe':
            vqe_config = config.get('vqe_config', {})
            result = execute_ground_vqe(
                molecule,
                ansatz=vqe_config['ansatz'],
                mapper=vqe_config['mapper'],
                optimizer=vqe_config['optimizer'],
                max_iterations=vqe_config['max_iterations'],  # FROM CONFIG!
                backend=vqe_config['backend'],
                ...
            )

        elif solver == 'hf':
            result = execute_hf(molecule, ...)

        elif solver == 'sqd':
            result = execute_sqd(molecule, ...)

    elif category == 'excited_state':
        solver = config.get('excited_solver')
        n_states = config.get('n_excited_states', 3)

        if solver == 'vqe':
            vqe_config = config.get('excited_vqe_config', {})
            result = execute_excited_vqe(
                molecule,
                n_states=n_states,
                ansatz=vqe_config['ansatz'],
                optimizer=vqe_config['optimizer'],
                max_iterations=vqe_config['max_iterations'],  # FROM CONFIG!
                penalty_weight=vqe_config.get('penalty_weight', 5.0),
                ...
            )

        elif solver == 'cis':
            result = execute_cis(molecule, n_states=n_states, ...)

        elif solver == 'tddft':
            result = execute_tddft(molecule, n_states=n_states, ...)
```

#### 1.3 Remove Hardcoded Values

**File**: `kanad/solvers/excited_states_solver.py`

**Changes**:
- ✅ Remove `max_iterations = getattr(self, '_max_iterations', 100)` default
- ✅ Require explicit `max_iterations` parameter
- ✅ Remove all other hardcoded solver parameters

---

### Phase 2: Frontend Restructuring

#### 2.1 Settings Modal Redesign

**File**: `web/src/components/settings/SettingsModal.tsx`

**New Structure**:
```typescript
<SettingsModal>
  {/* Step 1: Choose Category */}
  <CategorySelector>
    <Option value="ground_state">
      Ground State Calculator
      <Description>Find the lowest energy configuration</Description>
    </Option>
    <Option value="excited_state">
      Excited State Calculator
      <Description>Calculate higher energy states</Description>
    </Option>
  </CategorySelector>

  {/* Step 2: Choose Solver (based on category) */}
  {category === 'ground_state' && (
    <SolverSelector>
      <Option value="vqe">VQE - Variational Quantum</Option>
      <Option value="hf">Hartree-Fock</Option>
      <Option value="qpe">QPE - Quantum Phase</Option>
      <Option value="sqd">SQD - Subspace Diag</Option>
    </SolverSelector>
  )}

  {category === 'excited_state' && (
    <>
      <NStatesInput label="Number of Excited States" />
      <SolverSelector>
        <Option value="vqe">VQE - Quantum</Option>
        <Option value="cis">CIS - Classical</Option>
        <Option value="tddft">TDDFT - Classical</Option>
        <Option value="sqd">SQD - Hybrid</Option>
      </SolverSelector>
    </>
  )}

  {/* Step 3: Configure Selected Solver */}
  {(solver === 'vqe') && <VQEConfigPanel />}
  {(solver === 'hf') && <HFConfigPanel />}
  {(solver === 'cis') && <CISConfigPanel />}
  {(solver === 'tddft') && <TDDFTConfigPanel />}
  {(solver === 'sqd') && <SQDConfigPanel />}
</SettingsModal>
```

#### 2.2 VQE Configuration Component (Reusable)

**File**: `web/src/components/settings/VQEConfigPanel.tsx`

```typescript
export const VQEConfigPanel = ({ onChange, defaults }) => {
  return (
    <div className="space-y-4">
      <Select label="Ansatz" value={config.ansatz} onChange={...}>
        <Option value="hardware_efficient">Hardware Efficient</Option>
        <Option value="uccsd">UCCSD</Option>
        <Option value="ucc">UCC</Option>
      </Select>

      <Select label="Mapper" value={config.mapper} onChange={...}>
        <Option value="jordan_wigner">Jordan-Wigner</Option>
        <Option value="parity">Parity</Option>
        <Option value="bravyi_kitaev">Bravyi-Kitaev</Option>
      </Select>

      <Select label="Optimizer" value={config.optimizer} onChange={...}>
        <Option value="COBYLA">COBYLA</Option>
        <Option value="SLSQP">SLSQP</Option>
        <Option value="L-BFGS-B">L-BFGS-B</Option>
      </Select>

      <Slider
        label="Max Iterations"
        min={10}
        max={1000}
        value={config.max_iterations}
        onChange={...}
      />

      <Select label="Backend" value={config.backend} onChange={...}>
        <Option value="classical">Classical Simulator</Option>
        <Option value="bluequbit">BlueQubit Cloud</Option>
        <Option value="ibm">IBM Quantum</Option>
      </Select>

      {config.backend === 'bluequbit' && (
        <Select label="BlueQubit Device" ...>
          <Option value="cpu">CPU</Option>
          <Option value="gpu">GPU</Option>
        </Select>
      )}
    </div>
  );
};
```

---

### Phase 3: UI/UX Improvements

#### 3.1 Professional Method Cards

```typescript
<MethodCard category="ground_state" solver="vqe">
  <Icon>🔬</Icon>
  <Title>VQE Ground State</Title>
  <Description>
    Variational Quantum Eigensolver finds the ground state energy
    using a quantum-classical hybrid approach.
  </Description>
  <Badges>
    <Badge color="purple">Quantum</Badge>
    <Badge color="green">NISQ-Ready</Badge>
  </Badges>
  <EstimatedTime>~2-5 minutes</EstimatedTime>
</MethodCard>
```

#### 3.2 Configuration Summary Panel

Show user's final configuration before submitting:

```
┌─────────────────────────────────────┐
│ Experiment Configuration Summary    │
├─────────────────────────────────────┤
│ Category: Excited State Calculator  │
│ Solver: VQE (Quantum)               │
│ N States: 3                         │
│                                     │
│ VQE Settings:                       │
│  ├─ Ansatz: Hardware Efficient      │
│  ├─ Optimizer: COBYLA               │
│  ├─ Max Iterations: 50              │ ← USER SET, NOT 100!
│  ├─ Penalty Weight: 5.0             │
│  └─ Backend: BlueQubit (CPU)        │
│                                     │
│ Estimated Cost: $0.00               │
│ Estimated Time: ~3 minutes          │
└─────────────────────────────────────┘
```

---

## Migration Strategy

### Step 1: Backward Compatibility
Keep old API working while adding new structure:

```python
# Support both old and new formats
if 'method' in config:
    # Old format: method="EXCITED_STATES"
    category = map_old_method_to_category(config['method'])
    solver = map_old_method_to_solver(config['method'])
else:
    # New format: category="excited_state", solver="vqe"
    category = config['category']
    solver = config[f'{category}_solver']
```

### Step 2: Frontend Feature Flag
```typescript
const USE_NEW_SOLVER_UI = true;  // Toggle

{USE_NEW_SOLVER_UI ? (
  <NewSettingsModal />
) : (
  <OldSettingsModal />
)}
```

### Step 3: Gradual Migration
1. Week 1: Backend supports both old & new
2. Week 2: Frontend shows new UI (with fallback)
3. Week 3: Test thoroughly
4. Week 4: Remove old code

---

## Testing Plan

### Backend Tests
```bash
# Test ground state VQE with custom iterations
POST /api/experiments/submit
{
  "category": "ground_state",
  "ground_solver": "vqe",
  "vqe_config": {
    "max_iterations": 25,  # Should use 25, NOT 100!
    "ansatz": "ucc"
  }
}

# Test excited state VQE with custom iterations
POST /api/experiments/submit
{
  "category": "excited_state",
  "excited_solver": "vqe",
  "n_excited_states": 3,
  "excited_vqe_config": {
    "max_iterations": 30,  # Should use 30, NOT 100!
    "penalty_weight": 10.0
  }
}
```

### Frontend Tests
- [ ] Ground State: Select VQE, set iterations to 25, verify in logs
- [ ] Excited State: Select VQE, set iterations to 30, verify in logs
- [ ] Switch between solvers, verify correct config panels show
- [ ] Submit experiment, verify summary shows correct values

---

## Benefits

1. **Clear Mental Model**: Users understand Ground vs Excited
2. **No Hardcoded Values**: Everything configurable
3. **Solver-Specific UI**: Each solver shows only relevant options
4. **Extensible**: Easy to add new solvers (QPE, CCSD, etc.)
5. **Professional**: Matches scientific software UX standards

---

## Timeline

### Week 1: Backend (Current Sprint)
- [x] Remove hardcoded max_iterations from excited_states_solver.py
- [ ] Create new API schema (MethodCategory, etc.)
- [ ] Update experiment_service.py to handle both categories
- [ ] Add backward compatibility layer

### Week 2: Frontend
- [ ] Design new SettingsModal component tree
- [ ] Implement CategorySelector
- [ ] Create reusable VQEConfigPanel
- [ ] Build SolverSelector components
- [ ] Add configuration summary panel

### Week 3: Testing & Polish
- [ ] Integration testing
- [ ] Fix bugs
- [ ] UI/UX polish
- [ ] Performance optimization

### Week 4: Migration & Docs
- [ ] Remove old code
- [ ] Update documentation
- [ ] User guide with examples
- [ ] Video tutorials

