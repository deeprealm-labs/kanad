// API Types for Kanad Backend

export interface User {
  user_id: string;
  email: string;
  name?: string;
  institution?: string;
  field?: string;
}

export interface Atom {
  element: string;
  position: [number, number, number];
}

export interface Molecule {
  molecule_id: string;
  name?: string;
  formula: string;
  atoms?: Atom[];
  bonds?: number[][];
  charge: number;
  multiplicity: number;
  basis: string;
  n_electrons?: number;
  n_orbitals?: number;
  n_qubits?: number;
  geometry_optimized?: boolean;
}

export type ComputationMethod = "HF" | "VQE" | "MP2" | "SQD" | "EXCITED_STATES";
export type AnsatzType = "ucc" | "hardware_efficient" | "governance" | "two_local" | "ucc_correct_double";
export type MapperType = "jordan_wigner" | "bravyi_kitaev" | "parity" | "hybrid_orbital";
export type BackendType = "classical" | "ibm_quantum" | "bluequbit";
export type JobStatus = "queued" | "running" | "completed" | "failed" | "cancelled";

export interface SimulationConfig {
  simulation_id?: string;
  molecule_id: string;
  method: ComputationMethod;
  ansatz?: AnsatzType;
  mapper?: MapperType;
  optimizer: string;
  max_iterations: number;
  convergence_threshold?: number;
  backend: {
    type: BackendType;
    backend_name?: string;
    use_user_credentials?: boolean;
  };
  analysis: {
    energy_decomposition: boolean;
    bond_analysis: boolean;
    dipole_moment: boolean;
    thermochemistry: boolean;
    spectroscopy: boolean;
    vibrational: boolean;
    uncertainty: boolean;
    bond_scan: boolean;
    dos: boolean;
  };
  optimization?: {
    geometry: boolean;
    orbitals: boolean;
    circuit: boolean;
    adaptive: boolean;
  };
  advanced?: {
    active_space?: { n_electrons: number; n_orbitals: number };
    frozen_core: boolean;
    symmetry: string;
  };
}

export interface Job {
  job_id: string;
  user_id?: string;
  simulation_id?: string;
  molecule_id?: string;
  molecule_name?: string;
  method?: string;
  status: JobStatus;
  progress: number;
  created_at: string;
  started_at?: string;
  completed_at?: string;
  backend?: string;
  cloud_job_id?: string;
  current_iteration?: number;
  max_iterations?: number;
  current_energy?: number;
  best_energy?: number;
  message?: string;
}

export interface JobResults {
  job_id: string;
  status: string;
  molecule?: Molecule;
  results: {
    method: string;
    energy: number;
    hf_energy: number;
    correlation_energy: number;
    n_iterations: number;
    converged: boolean;
    convergence_history: Array<{ iteration: number; energy: number }>;
  };
  analysis?: {
    energy_decomposition?: {
      kinetic: number;
      nuclear_attraction: number;
      electron_repulsion: number;
    };
    bond_analysis?: {
      bonds: Array<{ atoms: number[]; order: number; length: number }>;
      homo_lumo_gap: number;
    };
    dipole_moment?: {
      magnitude: number;
      direction: [number, number, number];
    };
    thermochemistry?: {
      enthalpy: number;
      entropy: number;
      gibbs_free_energy: number;
    };
  };
  llm_report?: {
    summary: string;
    key_findings: string[];
    interpretation: string;
    recommendations: string[];
  };
}

export interface CloudCredentials {
  ibm_api?: string;
  ibm_crn?: string;
  blue_token?: string;
}

export interface Settings {
  computation: {
    default_basis: string;
    default_method: ComputationMethod;
    default_ansatz: AnsatzType;
    default_mapper: MapperType;
    default_optimizer: string;
  };
  optimization: {
    geometry_optimization: boolean;
    orbital_optimization: boolean;
    circuit_optimization: boolean;
    adaptive_vqe: boolean;
  };
  analysis: {
    auto_analyze: boolean;
    default_analyses: string[];
  };
  cloud: {
    default_backend: BackendType;
    auto_select_backend: boolean;
  };
}
