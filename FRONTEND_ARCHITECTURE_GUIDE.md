# Kanad Frontend Architecture Guide
**For: kanad-frontend-architect agent**

---

## 1. Project Overview

**Project Name**: Kanad Quantum Chemistry Framework - Web GUI
**Framework**: Next.js 15 (App Router) + TypeScript + Tailwind CSS 4 + shadcn/ui
**Backend API**: FastAPI at http://localhost:8000/api
**Location**: `/home/mk/deeprealm/kanad/web/`

---

## 2. Design Language & Theme

### Visual Identity
- **Primary Brand Color**: Dark Saffron/Orange (#FF8C00 or similar)
- **Color Scheme**:
  - Black backgrounds for branding/navigation
  - White/light backgrounds for content areas
  - High contrast, minimalist aesthetic
- **Artistic Theme**: The artistic image in `/public/image.webp` (figure with orange clouds) represents innovation and exploration

### Typography
1. **Logo Font**: "Bietro DEMO" (located at `/public/fonts/Bietro DEMO.otf`)
   - Use for "deeprealm" and "kanad" logos everywhere
2. **Body Font**: "Quando" (Google Font or similar elegant serif)
3. **UI Font**: Clean modern sans-serif (Inter, Geist, or system fonts)

### Design References
- **Home Page**: `/public/home.jpg`
- **Dashboard**: `/public/dashboard.jpg`

---

## 3. Backend API Documentation

### Base URL
```
http://localhost:8000/api
```

### Complete API Endpoints (36 total)

#### Authentication (`/api/auth`)
```typescript
POST /api/auth/register
  Body: { email, password, name?, institution?, field? }
  Returns: { user, access_token }

POST /api/auth/login
  Body: { email, password }
  Returns: { access_token, token_type: "bearer" }

POST /api/auth/refresh
  Headers: { Authorization: "Bearer <token>" }
  Returns: { access_token }
```

#### Molecules (`/api/molecules`)
```typescript
POST /api/molecules/create
  Body: {
    method: "atoms" | "smiles" | "library" | "xyz",
    data: { atoms, basis, charge, multiplicity } | { smiles, basis }
  }
  Returns: {
    molecule_id, formula, geometry, lewis_structure,
    n_electrons, n_orbitals, n_qubits, preview
  }

GET /api/molecules
  Returns: { molecules: [...] }

GET /api/molecules/{id}
  Returns: { molecule details }
```

#### Simulations (`/api/simulations`)
```typescript
POST /api/simulations/configure
  Body: {
    molecule_id: string,
    method: "HF" | "VQE" | "MP2" | "SQD" | "EXCITED_STATES",
    ansatz?: "ucc" | "hardware_efficient" | "governance" | "two_local" | "ucc_correct_double",
    mapper?: "jordan_wigner" | "bravyi_kitaev" | "parity" | "hybrid_orbital",
    optimizer: "SLSQP" | "COBYLA" | "L-BFGS-B",
    max_iterations: number,
    backend: {
      type: "classical" | "ibm_quantum" | "bluequbit",
      backend_name?: string,
      use_user_credentials?: boolean
    },
    analysis: {
      energy_decomposition: boolean,
      bond_analysis: boolean,
      dipole_moment: boolean,
      thermochemistry: boolean,
      spectroscopy: boolean,
      vibrational: boolean,
      uncertainty: boolean,
      bond_scan: boolean,
      dos: boolean
    },
    optimization: {
      geometry: boolean,
      orbitals: boolean,
      circuit: boolean,
      adaptive: boolean
    },
    advanced: {
      active_space?: { n_electrons, n_orbitals },
      frozen_core: boolean,
      symmetry: string
    }
  }
  Returns: {
    simulation_id,
    preview: { n_qubits, n_parameters, circuit_depth, estimated_time, estimated_cost },
    terms_and_conditions
  }

POST /api/simulations/{id}/submit
  Body: { accepted_terms: true, user_credentials?: {...} }
  Returns: { job_id, status: "queued" }
```

#### Jobs (`/api/jobs`)
```typescript
GET /api/jobs
  Returns: { jobs: [{ job_id, molecule_name, method, status, progress, created_at, backend }] }

GET /api/jobs/{id}/status
  Returns: {
    job_id, status, progress, current_iteration, max_iterations,
    current_energy, best_energy, message
  }

GET /api/jobs/{id}/results
  Returns: {
    job_id, status, molecule,
    results: { method, energy, hf_energy, correlation_energy, convergence_history },
    analysis: { energy_decomposition, bond_analysis, dipole_moment, thermochemistry, ... },
    llm_report: { summary, key_findings, interpretation, recommendations }
  }

DELETE /api/jobs/{id}
  Returns: { message: "Job cancelled" }

GET /api/jobs/{id}/export?format=json|csv|xyz
  Returns: File download

GET /api/jobs/{id}/report?format=pdf|html|markdown
  Returns: File download (LLM-generated report)
```

#### Analysis (`/api/analysis`)
```typescript
POST /api/analysis/energy-decomposition
  Body: { molecule_id }
  Returns: { kinetic, nuclear_attraction, electron_repulsion }

POST /api/analysis/bond-analysis
  Body: { molecule_id }
  Returns: { bonds: [{ atoms, order, length }], homo_lumo_gap }

POST /api/analysis/vibrational
  Body: { molecule_id }
  Returns: { frequencies, normal_modes, zero_point_energy }
```

#### Cloud Providers (`/api/cloud`)
```typescript
POST /api/cloud/credentials
  Body: {
    provider: "ibm_quantum" | "bluequbit",
    credentials: { ibm_api?, ibm_crn?, blue_token? }
  }
  Returns: { message: "Credentials stored securely" }

GET /api/cloud/backends
  Returns: {
    providers: [{
      name, status,
      backends: [{ name, qubits, queue_depth, online, estimated_wait }]
    }]
  }
```

#### Library (`/api/library`)
```typescript
GET /api/library
  Returns: {
    categories: [{
      name: "Organic" | "Inorganic" | "Proteins" | "Crystals",
      molecules: [{ id, name, formula, smiles }]
    }]
  }
```

#### Settings (`/api/settings`)
```typescript
GET /api/settings/defaults
  Returns: { computation, optimization, analysis, cloud }

PUT /api/settings/defaults
  Body: { defaults object }
  Returns: { message: "Settings updated" }
```

#### Domain-Specific: Metallurgy (`/api/metallurgy`)
```typescript
POST /api/metallurgy/crystal-structure
  Body: { lattice_vectors, basis_atoms, basis }
  Returns: { band_structure, dos, fermi_energy, band_gap }

POST /api/metallurgy/alloy-properties
  Body: { composition, crystal_structure }
  Returns: { properties }
```

#### Domain-Specific: Bioscience (`/api/bioscience`)
```typescript
POST /api/bioscience/protein-ligand
  Body: { protein_pdb, ligand_smiles }
  Returns: { binding_energy, binding_sites }

POST /api/bioscience/drug-properties
  Body: { molecule_smiles }
  Returns: { adme_properties }
```

#### Domain-Specific: Chemical Engineering (`/api/chemical-engineering`)
```typescript
POST /api/chemical-engineering/reaction-pathway
  Body: { reactants, products, catalyst? }
  Returns: { pathway, activation_energy, transition_states }

POST /api/chemical-engineering/catalyst-screening
  Body: { reaction, candidates }
  Returns: { ranked_catalysts }
```

#### Batch Scheduling (`/api/schedules`)
```typescript
POST /api/schedules/create
  Body: {
    name,
    experiments: [{ molecule_id, method, basis, ansatz, mapper }],
    execution_mode: "sequential" | "parallel",
    priority: "low" | "normal" | "high"
  }
  Returns: { schedule_id, total_experiments, estimated_total_time, job_ids }

GET /api/schedules/{id}/progress
  Returns: { completed, total, current_job, progress }
```

#### User Profile (`/api/user`)
```typescript
GET /api/user/profile
  Returns: { user_id, email, name, institution, field, stats }

GET /api/user/history
  Returns: { jobs: [...] }
```

---

## 4. Kanad Framework Capabilities

### Quantum Chemistry Methods
1. **Hartree-Fock (HF)**: Classical mean-field method
2. **VQE (Variational Quantum Eigensolver)**: Hybrid quantum-classical
3. **MP2 (Møller-Plesset)**: Post-HF perturbation theory
4. **SQD (Subspace Quantum Diagonalization)**: Advanced diagonalization
5. **Excited States**: Multiple electronic states

### Ansätze (VQE Circuit Structures)
1. **UCC**: Unitary Coupled Cluster
2. **Hardware-Efficient**: Optimized for quantum hardware
3. **Governance**: Ionic/Covalent/Metallic bonding-aware
4. **Two-Local**: Two-qubit gate structure
5. **UCC Correct Double**: Improved UCC variant

### Mappers (Fermion-to-Qubit)
1. **Jordan-Wigner**: Standard transformation
2. **Bravyi-Kitaev**: Reduced qubit requirements
3. **Parity**: Alternative mapping
4. **Hybrid Orbital**: For metallurgy applications

### Analysis Tools
1. **Energy Decomposition**: Kinetic, nuclear, electron repulsion
2. **Bond Analysis**: Bond orders, lengths, HOMO-LUMO gap
3. **Dipole Moment**: Molecular polarity
4. **Thermochemistry**: Enthalpy, entropy, Gibbs free energy
5. **Spectroscopy**: UV-Vis, excited states
6. **Vibrational Analysis**: Frequencies, normal modes
7. **Uncertainty Analysis**: Error estimation
8. **Bond Scanning**: Potential energy surfaces
9. **DOS**: Density of states (for materials)

### Optimization Features
1. **Geometry Optimization**: Molecular structure
2. **Orbital Optimization**: Molecular orbitals
3. **Circuit Optimization**: Quantum circuit depth/gate count
4. **Adaptive VQE**: Dynamic ansatz construction

### Supported Backends
1. **Classical Simulator**: Local CPU/GPU
2. **IBM Quantum**: Real quantum hardware (requires API key)
3. **BlueQubit**: Cloud quantum computing (requires token)

---

## 5. Page Structure & Flow

### Pages to Build

#### 1. **Home Page** (`/`)
**Reference**: `/public/home.jpg`

**Layout**:
```
┌─────────────────────────────────────────┐
│  Left (50%)         │  Right (50%)      │
│  White BG           │  Black BG         │
│                     │                   │
│  welcome            │  [Artistic Image] │
│  to the place for   │                   │
│  creation           │                   │
│  exploration        │                   │
│  &                  │  deeprealm (vert) │
│  invention          │                   │
│                     │  kanad (logo)     │
│  [email input]      │                   │
│  [Go] [Google]      │                   │
└─────────────────────────────────────────┘
```

**Features**:
- Email input + "Go" button for sign-up
- Google OAuth integration
- Smooth transition to dashboard after login
- Responsive: Stack vertically on mobile

**Components to Build**:
- `HeroSection.tsx`
- `AuthForm.tsx`
- `GoogleAuthButton.tsx`

---

#### 2. **Dashboard Page** (`/dashboard`)
**Reference**: `/public/dashboard.jpg`

**Layout**:
```
┌──────────────────────────────────────────────────┐
│ Sidebar │  Main Content Area                     │
│ (Black) │  (White)                    [Settings] │
│         │                                        │
│ [User]  │  click anywhere to start               │
│         │                                        │
│ Docs    │  [Dynamic Content Based on Mode]       │
│ Tuts    │                                        │
│         │                                        │
│ kanad   │                                        │
└──────────────────────────────────────────────────┘
```

**Sidebar Components**:
- User profile (avatar, name, role)
- Navigation: Docs, Tutorials, Library
- "kanad" logo at bottom (Bietro DEMO font)

**Main Area States**:
1. **Welcome State**: "click anywhere to start"
2. **Molecule Builder**: Create/import molecules
3. **Simulation Config**: Configure computations
4. **Job Monitor**: Track running jobs
5. **Results Viewer**: View and analyze results

**Settings Modal** (gear icon):
- Backend configuration (API URL)
- Cloud credentials (IBM Quantum, BlueQubit)
- Default computation settings
- User preferences

---

#### 3. **Experimentation Workflow** (Multi-step Windows)

As per API_BUILD_PLAN.md, the experimentation should follow these steps:

**Step 1: Molecule Creation**
```
┌─────────────────────────────────────────┐
│  Create Molecule                        │
│  ────────────────────────────────────── │
│                                         │
│  [ ] From Atoms                         │
│  [ ] From SMILES                        │
│  [ ] From Library                       │
│  [ ] Upload XYZ/PDB                     │
│                                         │
│  [Input fields based on selection]      │
│                                         │
│  [Preview: Lewis Structure]             │
│                                         │
│  [ Cancel ] [ Next: Configure ]         │
└─────────────────────────────────────────┘
```

**Step 2: Simulation Configuration**
```
┌─────────────────────────────────────────┐
│  Configure Simulation                   │
│  ────────────────────────────────────── │
│                                         │
│  Method:  [VQE ▼]                       │
│  Ansatz:  [Hardware-Efficient ▼]        │
│  Mapper:  [Jordan-Wigner ▼]             │
│  Optimizer: [SLSQP ▼]                   │
│  Max Iterations: [1000]                 │
│                                         │
│  Backend: [ Classical ] [ IBM ] [ BQ ]  │
│                                         │
│  Analysis Options:                      │
│  ☑ Energy Decomposition                 │
│  ☑ Bond Analysis                        │
│  ☐ Thermochemistry                      │
│  ☐ Vibrational                          │
│                                         │
│  Circuit Preview:                       │
│  • Qubits: 12                           │
│  • Parameters: 48                       │
│  • Depth: 156                           │
│  • Est. Time: ~2 hours                  │
│                                         │
│  [ Back ] [ Preview ] [ Submit ]        │
└─────────────────────────────────────────┘
```

**Step 3: Job Submission**
```
┌─────────────────────────────────────────┐
│  Submit Computation                     │
│  ────────────────────────────────────── │
│                                         │
│  Summary:                               │
│  • Molecule: H2O                        │
│  • Method: VQE                          │
│  • Backend: IBM Quantum (ibm_torino)    │
│  • Estimated Cost: 1,200 credits        │
│                                         │
│  ☑ I accept the terms and conditions    │
│                                         │
│  Cloud Credentials (if not saved):      │
│  API Key: [____________]                │
│                                         │
│  [ Cancel ] [ Submit Job ]              │
└─────────────────────────────────────────┘
```

**Step 4: Real-Time Monitoring**
```
┌─────────────────────────────────────────┐
│  Job #12345 - Running                   │
│  ────────────────────────────────────── │
│                                         │
│  Progress: [████████░░] 80%             │
│                                         │
│  Current Iteration: 800 / 1000          │
│  Current Energy: -1.13647 Ha            │
│  Best Energy: -1.13682 Ha               │
│                                         │
│  Real-time Logs:                        │
│  ┌───────────────────────────────────┐ │
│  │ [09:15:23] Iteration 798: -1.1365 │ │
│  │ [09:15:24] Iteration 799: -1.1364 │ │
│  │ [09:15:25] Iteration 800: -1.1365 │ │
│  │ ...                               │ │
│  └───────────────────────────────────┘ │
│                                         │
│  Convergence Chart:                     │
│  [Live updating line chart]             │
│                                         │
│  [ Cancel Job ] [ View Details ]        │
└─────────────────────────────────────────┘
```

**Step 5: Results Display**
```
┌─────────────────────────────────────────┐
│  Results - H2O VQE Calculation          │
│  ────────────────────────────────────── │
│                                         │
│  ┌ Energy Results ──────────────────┐  │
│  │ Ground State: -76.0267 Ha        │  │
│  │ HF Energy: -76.0109 Ha           │  │
│  │ Correlation: -0.0158 Ha          │  │
│  │ Converged: ✓ (985 iterations)    │  │
│  └──────────────────────────────────┘  │
│                                         │
│  ┌ Energy Decomposition ─────────────┐ │
│  │ [Pie chart: Kinetic, Nuclear,    │ │
│  │  Electron Repulsion]             │ │
│  └──────────────────────────────────┘  │
│                                         │
│  ┌ Bond Analysis ────────────────────┐ │
│  │ O-H1: 0.96 Å (order: 0.95)       │ │
│  │ O-H2: 0.96 Å (order: 0.95)       │ │
│  │ HOMO-LUMO Gap: 12.8 eV           │ │
│  └──────────────────────────────────┘  │
│                                         │
│  ┌ Convergence History ──────────────┐ │
│  │ [Line chart showing energy vs    │ │
│  │  iteration]                      │ │
│  └──────────────────────────────────┘  │
│                                         │
│  ┌ AI Report ────────────────────────┐ │
│  │ [LLM-generated summary]          │ │
│  │ Key Findings:                    │ │
│  │ • Strong correlation energy...   │ │
│  │ • Bond lengths consistent with...│ │
│  └──────────────────────────────────┘  │
│                                         │
│  [ Export JSON ] [ Export PDF ]         │
│  [ Download XYZ ] [ New Calculation ]   │
└─────────────────────────────────────────┘
```

---

## 6. UI Components to Build (shadcn/ui)

### Core Components
- `Button`, `Input`, `Select`, `Checkbox`, `Radio`
- `Card`, `Dialog`, `Sheet`, `Popover`
- `Table`, `Tabs`, `Accordion`
- `Progress`, `Badge`, `Avatar`
- `Tooltip`, `DropdownMenu`, `ContextMenu`

### Charts & Visualizations (shadcn charts)
- `LineChart` - Convergence history
- `PieChart` - Energy decomposition
- `BarChart` - Bond orders, property comparisons
- `ScatterPlot` - PES scans
- `HeatMap` - DOS visualization

### Custom Components
1. **MoleculeViewer3D**
   - Use Three.js or React Three Fiber
   - Display molecular structure
   - Interactive rotation, zoom
   - Highlight bonds, orbitals

2. **LewisStructureDiagram**
   - SVG-based 2D molecular structure
   - Show bonds, lone pairs, formal charges

3. **ConvergenceMonitor**
   - Real-time updating chart
   - WebSocket connection for live logs
   - Pause/resume, export data

4. **QuantumCircuitPreview**
   - Display quantum circuit diagram
   - Show qubits, gates, depth
   - Highlight ansatz structure

5. **CloudBackendSelector**
   - Grid of available backends
   - Show queue depth, estimated wait
   - Status indicators (online/offline)

6. **AnalysisResultsPanel**
   - Tabbed interface for different analyses
   - Charts, tables, visualizations
   - Export functionality

7. **JobQueueManager**
   - List of running/completed jobs
   - Filter, sort, search
   - Batch operations

---

## 7. State Management

### Recommended: Zustand or React Context

**Stores to Create**:

```typescript
// authStore.ts
interface AuthStore {
  user: User | null;
  token: string | null;
  login: (email, password) => Promise<void>;
  logout: () => void;
  isAuthenticated: boolean;
}

// moleculeStore.ts
interface MoleculeStore {
  molecules: Molecule[];
  currentMolecule: Molecule | null;
  createMolecule: (data) => Promise<Molecule>;
  deleteMolecule: (id) => Promise<void>;
}

// simulationStore.ts
interface SimulationStore {
  config: SimulationConfig;
  setMethod: (method) => void;
  setAnsatz: (ansatz) => void;
  setBackend: (backend) => void;
  submit: () => Promise<string>; // Returns job_id
}

// jobStore.ts
interface JobStore {
  jobs: Job[];
  currentJob: Job | null;
  fetchJobs: () => Promise<void>;
  cancelJob: (id) => Promise<void>;
  subscribeToJob: (id) => void; // WebSocket
}

// settingsStore.ts
interface SettingsStore {
  apiUrl: string;
  cloudCredentials: CloudCredentials;
  defaults: DefaultSettings;
  updateSettings: (settings) => Promise<void>;
}
```

---

## 8. WebSocket Integration

For **real-time job monitoring**:

```typescript
// hooks/useJobMonitor.ts
import { useEffect, useState } from 'react';

export function useJobMonitor(jobId: string) {
  const [status, setStatus] = useState<JobStatus>();
  const [logs, setLogs] = useState<string[]>([]);

  useEffect(() => {
    const ws = new WebSocket(`ws://localhost:8000/api/jobs/${jobId}/logs`);

    ws.onmessage = (event) => {
      const data = JSON.parse(event.data);
      setStatus(data.status);
      setLogs(prev => [...prev, data.log]);
    };

    return () => ws.close();
  }, [jobId]);

  return { status, logs };
}
```

---

## 9. Responsive Design

### Breakpoints (Tailwind)
- `sm`: 640px - Tablet portrait
- `md`: 768px - Tablet landscape
- `lg`: 1024px - Desktop
- `xl`: 1280px - Large desktop
- `2xl`: 1536px - Extra large

### Mobile Adaptations
- **Home Page**: Stack left/right vertically
- **Dashboard**: Hamburger menu for sidebar, full-width main area
- **Workflow**: Full-screen modals instead of side-by-side
- **Charts**: Responsive sizing, touch-friendly interactions

---

## 10. Deployment & Environment

### Environment Variables
```env
NEXT_PUBLIC_API_URL=http://localhost:8000/api
NEXT_PUBLIC_WS_URL=ws://localhost:8000/api
NEXTAUTH_SECRET=...
GOOGLE_CLIENT_ID=...
GOOGLE_CLIENT_SECRET=...
```

### Build Commands
```bash
npm run dev      # Development server (Turbopack)
npm run build    # Production build
npm run start    # Production server
```

---

## 11. Key Implementation Notes

### 1. **Font Integration**
```css
/* globals.css */
@font-face {
  font-family: 'Bietro';
  src: url('/fonts/Bietro DEMO.otf') format('opentype');
}

@import url('https://fonts.googleapis.com/css2?family=Quando&display=swap');

.logo {
  font-family: 'Bietro', sans-serif;
  color: #FF8C00; /* Dark saffron */
}

body {
  font-family: 'Quando', serif;
}
```

### 2. **Authentication Flow**
- Email/password login → JWT token
- Google OAuth → Backend verification → JWT token
- Store token in localStorage/cookies
- Include in Authorization header: `Bearer <token>`

### 3. **Error Handling**
```typescript
try {
  const response = await fetch('/api/endpoint');
  if (!response.ok) {
    throw new Error(`HTTP ${response.status}`);
  }
  const data = await response.json();
} catch (error) {
  toast.error('Operation failed: ' + error.message);
}
```

### 4. **Loading States**
- Skeleton loaders for data fetching
- Progress bars for long operations
- Spinner overlays for form submissions

### 5. **Accessibility**
- ARIA labels for interactive elements
- Keyboard navigation support
- Screen reader compatibility
- High contrast mode

---

## 12. Priority Features (MVP)

### Phase 1: Core Functionality ⭐
1. ✅ Home page with authentication
2. ✅ Dashboard with sidebar navigation
3. ✅ Molecule creation (atoms, SMILES)
4. ✅ VQE configuration (basic options)
5. ✅ Job submission and monitoring
6. ✅ Results display (energy, convergence chart)

### Phase 2: Enhanced Features ⭐⭐
7. ⏳ Settings modal (cloud credentials, defaults)
8. ⏳ Complete analysis visualizations
9. ⏳ Molecule library integration
10. ⏳ Export functionality (JSON, PDF, XYZ)
11. ⏳ Real-time WebSocket monitoring

### Phase 3: Advanced Features ⭐⭐⭐
12. ⏳ 3D molecular visualization
13. ⏳ Domain-specific tools (metallurgy, bioscience, chem eng)
14. ⏳ Batch scheduling
15. ⏳ User profile & history
16. ⏳ Docs & tutorials

---

## 13. Testing Strategy

### Unit Tests (Jest + React Testing Library)
- Component rendering
- User interactions
- State management
- API integration

### E2E Tests (Playwright)
- Complete workflow: Login → Create Molecule → Run VQE → View Results
- Settings configuration
- Error handling

---

## 14. Performance Optimization

1. **Code Splitting**: Dynamic imports for large components
2. **Image Optimization**: Next.js Image component
3. **Data Fetching**: SWR or React Query for caching
4. **Lazy Loading**: Charts, 3D viewer loaded on demand
5. **Memoization**: React.memo for heavy computations

---

## 15. File Structure

```
/home/mk/deeprealm/kanad/web/
├── public/
│   ├── fonts/
│   │   └── Bietro DEMO.otf
│   ├── home.jpg
│   ├── dashboard.jpg
│   └── image.webp
├── src/
│   ├── app/
│   │   ├── layout.tsx
│   │   ├── page.tsx (Home)
│   │   ├── dashboard/
│   │   │   └── page.tsx
│   │   ├── login/
│   │   └── api/ (if using Next.js API routes)
│   ├── components/
│   │   ├── ui/ (shadcn components)
│   │   ├── molecule/
│   │   │   ├── MoleculeBuilder.tsx
│   │   │   ├── MoleculeViewer3D.tsx
│   │   │   └── LewisStructure.tsx
│   │   ├── simulation/
│   │   │   ├── ConfigWizard.tsx
│   │   │   ├── MethodSelector.tsx
│   │   │   └── BackendSelector.tsx
│   │   ├── job/
│   │   │   ├── JobMonitor.tsx
│   │   │   ├── JobList.tsx
│   │   │   └── ResultsViewer.tsx
│   │   ├── analysis/
│   │   │   ├── EnergyChart.tsx
│   │   │   ├── BondTable.tsx
│   │   │   └── ConvergencePlot.tsx
│   │   ├── settings/
│   │   │   └── SettingsModal.tsx
│   │   └── layout/
│   │       ├── Sidebar.tsx
│   │       └── Header.tsx
│   ├── lib/
│   │   ├── api/ (API client functions)
│   │   ├── utils.ts
│   │   └── constants.ts
│   ├── hooks/
│   │   ├── useAuth.ts
│   │   ├── useJobMonitor.ts
│   │   └── useApi.ts
│   ├── store/ (Zustand stores)
│   │   ├── authStore.ts
│   │   ├── moleculeStore.ts
│   │   └── jobStore.ts
│   └── types/
│       └── api.ts (TypeScript interfaces)
├── components.json (shadcn config)
├── .mcp.json (shadcn MCP)
└── package.json
```

---

## 16. Success Criteria

The frontend is complete when:
- ✅ Home page matches `/public/home.jpg` design
- ✅ Dashboard matches `/public/dashboard.jpg` design
- ✅ User can create molecule from SMILES
- ✅ User can configure and submit VQE calculation
- ✅ User can monitor job progress in real-time
- ✅ User can view results with charts and analysis
- ✅ Settings modal works for cloud credentials
- ✅ All API endpoints integrated correctly
- ✅ Responsive design works on mobile/tablet
- ✅ "kanad" logo uses Bietro DEMO font everywhere
- ✅ Clean, professional, high-performance UI

---

**Ready to build! 🚀**
