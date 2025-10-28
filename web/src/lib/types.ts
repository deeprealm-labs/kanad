// API Types for Kanad Frontend

// ===== Molecule Types =====
export interface Atom {
  symbol: string;
  x: number;
  y: number;
  z: number;
}

export interface MoleculeConfig {
  smiles?: string;
  atoms?: Atom[];
  basis: string;
  charge: number;
  multiplicity: number;
}

// ===== Backend Configuration Types =====
export interface BackendSettings {
  method: "HF" | "VQE" | "SQD" | "QPE" | "MP2" | "FCI" | "EXCITED_STATES";
  ansatz?: "ucc" | "uccsd" | "hardware_efficient" | "governance" | "hea";
  mapper?: "jordan_wigner" | "bravyi_kitaev" | "hybrid_orbital" | "parity";
  hamiltonian?: "covalent" | "ionic" | "metallic" | "custom";
  optimizer?: "SLSQP" | "COBYLA" | "L-BFGS-B" | "Powell" | "Nelder-Mead" | "CG" | "BFGS" | "TNC" | "ADAM";
  backend: "classical" | "ibm_quantum" | "bluequbit";
  backendName?: string; // e.g., "ibm_torino", "ibm_brisbane"
  bluequbitDevice?: "cpu" | "gpu" | "mps.cpu" | "mps.gpu" | "pauli-path"; // BlueQubit device selection
  maxIterations?: number; // Maximum VQE iterations (default: 100)
  // Excited States settings
  excitedMethod?: "cis" | "sqd" | "tddft" | "vqe"; // Method for computing excited states
  nStates?: number; // Number of states (ground + excited)
  // SQD-specific settings
  subspaceDim?: number;
  circuitDepth?: number;
  optimization?: {
    geometry: boolean;
    orbitals: boolean;
    circuit: boolean;
    adaptive: boolean;
  };
}

// ===== Analysis Configuration =====
export interface AnalysisConfig {
  energyDecomposition?: boolean;
  bondAnalysis?: boolean;
  dipoleMoment?: boolean;
  polarizability?: boolean;
  thermochemistry?: boolean;
  spectroscopy?: boolean;
}

// ===== Experiment Types =====
export type ExperimentStatus = "queued" | "running" | "completed" | "failed" | "cancelled";

export interface ExperimentResults {
  energy: number;
  dipoleMoment?: number;
  bondLengths?: number[];
  converged: boolean;
  iterations: number;
  convergenceData?: ConvergencePoint[];
  polarizability?: number;
  thermochemistry?: any;
  spectroscopy?: any;
}

export interface ConvergencePoint {
  iteration: number;
  energy: number;
  timestamp?: string;
}

export interface Experiment {
  id: string;
  status: ExperimentStatus;
  molecule: MoleculeConfig;
  backendSettings: BackendSettings;
  analysis?: AnalysisConfig;
  results?: ExperimentResults;
  method: string;
  backend: string;
  timestamp: string;
  startedAt?: string;
  completedAt?: string;
  errorMessage?: string;
}

// ===== Queue Types =====
export interface QueuedJob {
  id: string;
  name: string;
  molecule: MoleculeConfig;
  method: string;
  backend: string;
  scheduledTime?: string;
  priority: number;
  status: "queued" | "scheduled" | "running" | "paused";
  createdAt: string;
  backendSettings: BackendSettings;
  analysis?: AnalysisConfig;
}

// ===== API Request/Response Types =====

// Experiment Submission
export interface CreateExperimentRequest {
  molecule: MoleculeConfig;
  backendSettings: BackendSettings;
  analysis?: AnalysisConfig;
  executeNow?: boolean; // true for immediate execution, false to add to queue
}

export interface CreateExperimentResponse {
  experimentId: string;
  status: ExperimentStatus;
  message: string;
}

// Get Experiments
export interface GetExperimentsRequest {
  status?: ExperimentStatus;
  limit?: number;
  offset?: number;
  search?: string;
}

export interface GetExperimentsResponse {
  experiments: Experiment[];
  total: number;
  limit: number;
  offset: number;
}

// Get Single Experiment
export interface GetExperimentResponse {
  experiment: Experiment;
}

// Update Experiment
export interface UpdateExperimentRequest {
  status?: ExperimentStatus;
}

// Queue Management
export interface GetQueueResponse {
  queue: QueuedJob[];
  total: number;
}

export interface UpdateQueueItemRequest {
  priority?: number;
  status?: "queued" | "scheduled" | "running" | "paused";
  scheduledTime?: string;
}

// Settings
export interface UserSettings {
  method: string;
  ansatz: string;
  mapper: string;
  optimizer: string;
  backend: string;
  backendName?: string;
  bluequbitDevice?: string;
  maxIterations?: number;
  optimization: {
    geometry: boolean;
    orbitals: boolean;
    circuit: boolean;
    adaptive: boolean;
  };
}

export interface GetSettingsResponse {
  settings: UserSettings;
}

export interface UpdateSettingsRequest {
  settings: UserSettings;
}

// SMILES Validation
export interface ValidateSmilesRequest {
  smiles: string;
}

export interface ValidateSmilesResponse {
  valid: boolean;
  message?: string;
  molecularFormula?: string;
  molecularWeight?: number;
}

// Export
export type ExportFormat = "json" | "csv";

export interface ExportExperimentRequest {
  experimentId: string;
  format: ExportFormat;
}

// Error Response
export interface ApiError {
  error: string;
  message: string;
  statusCode: number;
  details?: any;
}

// WebSocket Message Types
export interface WSMessage {
  type: "status" | "progress" | "convergence" | "log" | "result" | "error";
  experimentId: string;
  data: any;
}

export interface WSStatusUpdate {
  status: ExperimentStatus;
  progress?: number;
}

export interface WSConvergenceUpdate {
  iteration: number;
  energy: number;
  timestamp: string;
}

export interface WSLogUpdate {
  message: string;
  timestamp: string;
  level: "info" | "warning" | "error";
}

export interface WSResultUpdate {
  results: ExperimentResults;
}

// Circuit Data
export interface CircuitData {
  diagram: string; // ASCII or visual representation
  depth: number;
  gates: number;
  parameters: number;
  qubits: number;
  operations?: Array<{
    gate: string;
    qubits: number[];
    parameters?: number[];
  }>;
}

// Experiment Report
export interface ExperimentReport {
  experiment_id: string;
  molecule: {
    smiles?: string;
    formula?: string;
    atoms: Atom[];
    charge: number;
    multiplicity: number;
  };
  configuration: BackendSettings;
  results: {
    energy: number;
    hf_energy?: number;
    correlation_energy?: number;
    dipole_moment?: number;
    convergence: {
      converged: boolean;
      iterations: number;
      final_gradient?: number;
    };
    properties?: {
      [key: string]: any;
    };
    // SQD-specific fields
    subspace_dim?: number;
    circuit_depth?: number;
    energies?: number[];
  };
  convergence_data?: ConvergencePoint[];
  circuit?: CircuitData;
  timestamp: string;
  duration?: number;
}

// Configuration Options
export interface ConfigurationOptions {
  methods: string[];
  ansatze: string[];
  mappers: string[];
  hamiltonians: string[];
  optimizers: string[];
  basis_sets: string[];
  backends: string[];
}
