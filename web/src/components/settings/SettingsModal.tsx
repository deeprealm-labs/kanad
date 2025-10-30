"use client";

import { useState, useEffect } from "react";
import { X, Save, Check } from "lucide-react";
import * as api from "@/lib/api";
import { useToast } from "@/components/ui/toast";

interface SettingsModalProps {
  isOpen: boolean;
  onClose: () => void;
  onSave?: () => void;
}

const DEFAULT_SETTINGS = {
  method: "VQE",
  ansatz: "hardware_efficient",
  mapper: "jordan_wigner",
  optimizer: "COBYLA",  // Default to COBYLA - better for cloud backends
  backend: "classical",
  backendName: "ibm_torino",
  maxIterations: 100,
  optimization: {
    geometry: false,
    orbitals: false,
    circuit: true,
    adaptive: false,
  },
};

export default function SettingsModal({ isOpen, onClose, onSave }: SettingsModalProps) {
  // Method: VQE, HF, SQD (ground state only - excited states in Advanced Analysis)
  const [method, setMethod] = useState("VQE");

  // Advanced Analysis (optional - runs during experiment if enabled)
  const [advancedAnalysisEnabled, setAdvancedAnalysisEnabled] = useState(false);
  const [advancedAnalysisProfile, setAdvancedAnalysisProfile] = useState<string | null>(null);

  // Common settings
  const [ansatz, setAnsatz] = useState(DEFAULT_SETTINGS.ansatz);
  const [mapper, setMapper] = useState(DEFAULT_SETTINGS.mapper);
  const [optimizer, setOptimizer] = useState(DEFAULT_SETTINGS.optimizer);
  const [backend, setBackend] = useState(DEFAULT_SETTINGS.backend);
  const [backendName, setBackendName] = useState(DEFAULT_SETTINGS.backendName);
  const [bluequbitDevice, setBluequbitDevice] = useState("gpu");
  const [hamiltonian, setHamiltonian] = useState("molecular");
  const [basisSet, setBasisSet] = useState("sto-3g");
  const [maxIterations, setMaxIterations] = useState(DEFAULT_SETTINGS.maxIterations);

  // SQD-specific settings
  const [subspaceDim, setSubspaceDim] = useState(10);
  const [circuitDepth, setCircuitDepth] = useState(3);
  const [nStates, setNStates] = useState(3);

  const [optimization, setOptimization] = useState(DEFAULT_SETTINGS.optimization);
  const [analysis, setAnalysis] = useState({
    energy_decomposition: true,
    bond_analysis: true,
    dipole_moment: true,
    polarizability: false,
    thermochemistry: false,
    spectroscopy: false,
    vibrational: false,
  });
  const [configOptions, setConfigOptions] = useState<any>(null);
  const toast = useToast();

  // Load configuration options from API
  useEffect(() => {
    const loadConfig = async () => {
      try {
        const options = await api.getConfigurationOptions();
        setConfigOptions(options);
      } catch (error) {
        console.error("Failed to load configuration options:", error);
      }
    };
    if (isOpen) {
      loadConfig();
    }
  }, [isOpen]);

  // Load settings from API and localStorage on mount
  useEffect(() => {
    const loadSettings = async () => {
      try {
        const settings = await api.getSettings();

        setMethod(settings.method || DEFAULT_SETTINGS.method);
        setAnsatz(settings.ansatz || DEFAULT_SETTINGS.ansatz);
        setMapper(settings.mapper || DEFAULT_SETTINGS.mapper);
        setOptimizer(settings.optimizer || DEFAULT_SETTINGS.optimizer);
        setBackend(settings.backend || DEFAULT_SETTINGS.backend);
        setBackendName(settings.backendName || DEFAULT_SETTINGS.backendName);
        setBluequbitDevice(settings.bluequbitDevice || "gpu");
        setHamiltonian(settings.hamiltonian || "molecular");
        setBasisSet(settings.basisSet || "sto-3g");
        setMaxIterations(settings.maxIterations || DEFAULT_SETTINGS.maxIterations);
        setSubspaceDim(settings.subspaceDim || 10);
        setCircuitDepth(settings.circuitDepth || 3);
        setNStates(settings.nStates || 3);
        setOptimization(settings.optimization || DEFAULT_SETTINGS.optimization);
        setAnalysis(settings.analysis || analysis);

        // Load advanced analysis settings
        setAdvancedAnalysisEnabled(settings.advancedAnalysisEnabled || false);
        setAdvancedAnalysisProfile(settings.advancedAnalysisProfile || null);
      } catch (error) {
        console.error("Failed to load settings from API:", error);
        // Fallback to localStorage
        const saved = localStorage.getItem("kanad_settings");
        if (saved) {
          try {
            const settings = JSON.parse(saved);

            setMethod(settings.method || DEFAULT_SETTINGS.method);
            setAnsatz(settings.ansatz || DEFAULT_SETTINGS.ansatz);
            setMapper(settings.mapper || DEFAULT_SETTINGS.mapper);
            setOptimizer(settings.optimizer || DEFAULT_SETTINGS.optimizer);
            setBackend(settings.backend || DEFAULT_SETTINGS.backend);
            setBackendName(settings.backendName || DEFAULT_SETTINGS.backendName);
            setBluequbitDevice(settings.bluequbitDevice || "gpu");
            setHamiltonian(settings.hamiltonian || "molecular");
            setBasisSet(settings.basisSet || "sto-3g");
            setMaxIterations(settings.maxIterations || DEFAULT_SETTINGS.maxIterations);
            setSubspaceDim(settings.subspaceDim || 10);
            setCircuitDepth(settings.circuitDepth || 3);
            setNStates(settings.nStates || 3);
            setOptimization(settings.optimization || DEFAULT_SETTINGS.optimization);
            setAnalysis(settings.analysis || analysis);
            setAdvancedAnalysisEnabled(settings.advancedAnalysisEnabled || false);
            setAdvancedAnalysisProfile(settings.advancedAnalysisProfile || null);
          } catch (e) {
            console.error("Failed to parse localStorage settings:", e);
          }
        }
      }
    };
    if (isOpen) {
      loadSettings();
    }
  }, [isOpen]);

  const handleSave = async () => {
    const settings = {
      method,
      ansatz,
      mapper,
      optimizer,
      backend,
      backendName,
      bluequbitDevice,
      hamiltonian,
      basisSet,
      maxIterations,
      subspaceDim,
      circuitDepth,
      nStates,
      optimization,
      analysis,
      advancedAnalysisEnabled,
      advancedAnalysisProfile,
    };

    try {
      await api.updateSettings(settings);
      console.log("Saved settings to API:", settings);
      toast.success("Settings saved successfully");
    } catch (error) {
      console.error("Failed to save settings to API:", error);
      toast.error("Failed to save settings to database, saved locally instead");
      // Fallback to localStorage
      localStorage.setItem("kanad_settings", JSON.stringify(settings));
      console.log("Saved settings to localStorage:", settings);
    }

    // Notify parent that settings were saved
    if (onSave) {
      onSave();
    }

    onClose();
  };

  if (!isOpen) return null;

  return (
    <div className="fixed inset-0 z-50 flex items-center justify-center backdrop-blur-sm bg-black/20">
      <div className="bg-background border border-border rounded-lg w-full max-w-4xl max-h-[90vh] overflow-hidden flex flex-col">
        {/* Header */}
        <div className="flex items-center justify-between p-6 border-b border-border">
          <h2 className="text-2xl font-quando font-bold">Settings</h2>
          <button
            onClick={onClose}
            className="p-2 hover:bg-accent rounded-md transition"
          >
            <X className="w-5 h-5" />
          </button>
        </div>

        {/* Content */}
        <div className="flex-1 overflow-auto p-6 space-y-6">
          {/* Method Selection */}
          <div>
            <h3 className="text-lg font-quando font-semibold mb-3">
              Method
            </h3>
            <select
              value={method}
              onChange={(e) => setMethod(e.target.value)}
              className="w-full px-4 py-3 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
            >
              <option value="VQE">VQE - Variational Quantum Eigensolver</option>
              <option value="HF">HF - Hartree-Fock (Classical)</option>
              <option value="SQD">SQD - Subspace Quantum Diagonalization</option>
            </select>
            <p className="text-xs text-muted-foreground mt-2">
              All methods compute ground state. For excited states, enable Advanced Analysis below.
            </p>
          </div>

          {/* VQE Settings */}
          {method === "VQE" && (
            <div className="bg-card border border-border rounded-lg p-6 space-y-4">
              <h3 className="text-lg font-quando font-semibold">
                VQE Configuration
              </h3>

              <div className="grid grid-cols-2 gap-4">
                <div>
                  <label className="block text-sm font-quando font-medium mb-2">
                    Ansatz
                  </label>
                  <select
                    value={ansatz}
                    onChange={(e) => setAnsatz(e.target.value)}
                    className="w-full px-4 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                  >
                    {configOptions?.ansatze?.map((a: any) => (
                      <option key={a.value} value={a.value}>
                        {a.label} - {a.description}
                      </option>
                    )) || (
                      <>
                        <option value="ucc">UCC (Higher accuracy)</option>
                        <option value="hardware_efficient">Hardware-Efficient (Faster)</option>
                      </>
                    )}
                  </select>
                </div>

                <div>
                  <label className="block text-sm font-quando font-medium mb-2">
                    Qubit Mapper
                  </label>
                  <select
                    value={mapper}
                    onChange={(e) => setMapper(e.target.value)}
                    className="w-full px-4 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                  >
                    {configOptions?.mappers?.map((m: any) => (
                      <option key={m.value} value={m.value}>
                        {m.label} - {m.description}
                      </option>
                    )) || (
                      <>
                        <option value="jordan_wigner">Jordan-Wigner (Standard)</option>
                        <option value="bravyi_kitaev">Bravyi-Kitaev (Reduced gates)</option>
                      </>
                    )}
                  </select>
                </div>

                <div>
                  <label className="block text-sm font-quando font-medium mb-2">
                    Optimizer
                  </label>
                  <select
                    value={optimizer}
                    onChange={(e) => setOptimizer(e.target.value)}
                    className="w-full px-4 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                  >
                    {configOptions?.optimizers?.map((opt: any) => (
                      <option key={opt.value} value={opt.value}>
                        {opt.label} - {opt.description}
                      </option>
                    )) || (
                      <>
                        <option value="COBYLA">COBYLA (Recommended for cloud - ~1-3x iterations)</option>
                        <option value="POWELL">POWELL (Fast, ~1-2x iterations)</option>
                        <option value="SLSQP">SLSQP (Gradient-based - ~40x iterations!)</option>
                        <option value="L-BFGS-B">L-BFGS-B (Gradient-based - ~40x iterations!)</option>
                      </>
                    )}
                  </select>
                  {backend !== "classical" && (optimizer === "SLSQP" || optimizer === "L-BFGS-B") && (
                    <p className="text-xs text-yellow-600 dark:text-yellow-400 mt-1">
                      ⚠️ Warning: {optimizer} uses ~40 function evaluations per iteration.
                      With {maxIterations || 100} max iterations, this means ~{(maxIterations || 100) * 40} quantum jobs!
                      Consider COBYLA or POWELL for fewer jobs.
                    </p>
                  )}
                </div>

                <div>
                  <label className="block text-sm font-quando font-medium mb-2">
                    Hamiltonian
                  </label>
                  <select
                    value={hamiltonian}
                    onChange={(e) => setHamiltonian(e.target.value)}
                    className="w-full px-4 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                  >
                    {configOptions?.hamiltonians?.map((h: any) => (
                      <option key={h.value} value={h.value}>
                        {h.label} - {h.description}
                      </option>
                    )) || (
                      <>
                        <option value="molecular">Molecular</option>
                        <option value="ionic">Ionic</option>
                        <option value="covalent">Covalent</option>
                      </>
                    )}
                  </select>
                </div>

                <div>
                  <label className="block text-sm font-quando font-medium mb-2">
                    Basis Set
                  </label>
                  <select
                    value={basisSet}
                    onChange={(e) => setBasisSet(e.target.value)}
                    className="w-full px-4 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                  >
                    {configOptions?.basis_sets?.map((b: any) => (
                      <option key={b.value} value={b.value}>
                        {b.label}
                      </option>
                    )) || (
                      <>
                        <option value="sto-3g">STO-3G</option>
                        <option value="6-31g">6-31G</option>
                      </>
                    )}
                  </select>
                </div>

                <div>
                  <label className="block text-sm font-quando font-medium mb-2">
                    Maximum Iterations
                  </label>
                  <input
                    type="number"
                    min="10"
                    max="1000"
                    step="10"
                    value={maxIterations}
                    onChange={(e) => setMaxIterations(parseInt(e.target.value) || 100)}
                    className="w-full px-4 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                  />
                  <p className="text-xs text-muted-foreground mt-1">
                    Maximum number of VQE optimization iterations (10-1000). Lower values = faster but less accurate.
                  </p>
                </div>
              </div>
            </div>
          )}

          {/* SQD Settings */}
          {method === "SQD" && (
            <div className="bg-card border border-border rounded-lg p-6 space-y-4">
              <h3 className="text-lg font-quando font-semibold">
                SQD Configuration
              </h3>
              <p className="text-sm text-muted-foreground">
                Subspace Quantum Diagonalization - Finds multiple energy eigenstates with lower circuit depth
              </p>

              <div className="grid grid-cols-2 gap-4">
                <div>
                  <label className="block text-sm font-quando font-medium mb-2">
                    Subspace Dimension
                  </label>
                  <input
                    type="number"
                    min="2"
                    max="20"
                    step="1"
                    value={subspaceDim}
                    onChange={(e) => setSubspaceDim(parseInt(e.target.value) || 10)}
                    className="w-full px-4 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                  />
                  <p className="text-xs text-muted-foreground mt-1">
                    Size of quantum subspace (2-20). Higher = more accurate but slower.
                  </p>
                </div>

                <div>
                  <label className="block text-sm font-quanto font-medium mb-2">
                    Circuit Depth
                  </label>
                  <input
                    type="number"
                    min="1"
                    max="10"
                    step="1"
                    value={circuitDepth}
                    onChange={(e) => setCircuitDepth(parseInt(e.target.value) || 3)}
                    className="w-full px-4 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                  />
                  <p className="text-xs text-muted-foreground mt-1">
                    Quantum circuit layers (1-10).
                  </p>
                </div>

                <div>
                  <label className="block text-sm font-quando font-medium mb-2">
                    Number of States
                  </label>
                  <input
                    type="number"
                    min="1"
                    max="10"
                    step="1"
                    value={nStates}
                    onChange={(e) => setNStates(parseInt(e.target.value) || 3)}
                    className="w-full px-4 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                  />
                  <p className="text-xs text-muted-foreground mt-1">
                    Number of energy states to compute.
                  </p>
                </div>
              </div>
            </div>
          )}

          {/* Backend Selection */}
          <div>
            <h3 className="text-lg font-quando font-semibold mb-3">
              Backend Selection
            </h3>
            <div className="space-y-3">
              <label className="flex items-center p-4 border border-border rounded-lg hover:bg-accent cursor-pointer transition">
                <input
                  type="radio"
                  name="backend"
                  value="classical"
                  checked={backend === "classical"}
                  onChange={(e) => setBackend(e.target.value)}
                  className="mr-3"
                />
                <div className="flex-1">
                  <div className="font-quando font-medium">
                    Classical Simulation (Local)
                  </div>
                  <div className="text-xs text-muted-foreground mt-1">
                    Run on your local machine
                  </div>
                </div>
              </label>

              {/* IBM Quantum - Now with real-time visualization */}
              <label className="flex items-center p-4 border border-border rounded-lg hover:bg-accent cursor-pointer transition">
                <input
                  type="radio"
                  name="backend"
                  value="ibm_quantum"
                  checked={backend === "ibm_quantum"}
                  onChange={(e) => setBackend(e.target.value)}
                  className="mr-3"
                />
                <div className="flex-1">
                  <div className="font-quando font-medium">
                    IBM Quantum (Batch Mode + Live Updates)
                  </div>
                  <div className="text-xs text-muted-foreground mt-1">
                    Real quantum hardware. Jobs run in batch mode with live convergence visualization.
                  </div>
                </div>
              </label>

              <label className="flex items-center p-4 border border-border rounded-lg hover:bg-accent cursor-pointer transition">
                <input
                  type="radio"
                  name="backend"
                  value="bluequbit"
                  checked={backend === "bluequbit"}
                  onChange={(e) => setBackend(e.target.value)}
                  className="mr-3"
                />
                <div className="flex-1">
                  <div className="font-quando font-medium">
                    BlueQubit (GPU-accelerated)
                  </div>
                  <div className="text-xs text-muted-foreground mt-1">
                    Fast cloud simulation (requires account credits)
                  </div>
                  {backend === "bluequbit" && (
                    <select
                      value={bluequbitDevice}
                      onChange={(e) => setBluequbitDevice(e.target.value)}
                      className="mt-2 w-full px-3 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando text-sm"
                      onClick={(e) => e.stopPropagation()}
                    >
                      <option value="gpu">GPU (36 qubits, fastest)</option>
                      <option value="cpu">CPU (34 qubits)</option>
                      <option value="mps.gpu">MPS GPU (40+ qubits)</option>
                      <option value="mps.cpu">MPS CPU (40+ qubits)</option>
                      <option value="pauli-path">Pauli-Path (50+ qubits, requires config)</option>
                    </select>
                  )}
                </div>
              </label>
            </div>
          </div>

          {/* Advanced Analysis */}
          <div className="bg-card border border-border rounded-lg p-6">
            <h3 className="text-lg font-quando font-semibold mb-3">
              Advanced Analysis
            </h3>
            <p className="text-sm text-muted-foreground mb-4">
              Run domain-specific advanced analysis during experiment (excited states, ADME, NMR, etc.)
            </p>

            <label className="flex items-center cursor-pointer mb-4">
              <input
                type="checkbox"
                checked={advancedAnalysisEnabled}
                onChange={(e) => setAdvancedAnalysisEnabled(e.target.checked)}
                className="mr-3"
              />
              <span className="font-quando text-sm font-medium">Enable Advanced Analysis</span>
            </label>

            {advancedAnalysisEnabled && (
              <div>
                <label className="block text-sm font-quando font-medium mb-2">
                  Analysis Profile
                </label>
                <select
                  value={advancedAnalysisProfile || ""}
                  onChange={(e) => setAdvancedAnalysisProfile(e.target.value || null)}
                  className="w-full px-4 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                >
                  <option value="">Select a profile...</option>
                  <option value="drug_discovery">Drug Discovery - ADME Properties</option>
                  <option value="spectroscopy">Spectroscopy - UV-Vis, NMR, Excited States</option>
                  <option value="materials">Materials Science - Density of States</option>
                  <option value="catalysis">Catalysis - Transition State Analysis</option>
                </select>
                <p className="text-xs text-muted-foreground mt-2">
                  {advancedAnalysisProfile === "drug_discovery" && "📊 Molecular weight, LogP, BBB permeability, drug-likeness"}
                  {advancedAnalysisProfile === "spectroscopy" && "🌈 UV-Vis spectra, excited states, vibrational frequencies"}
                  {advancedAnalysisProfile === "materials" && "⚡ Electronic density of states, band structure"}
                  {advancedAnalysisProfile === "catalysis" && "🔬 Vibrational frequencies, transition state verification"}
                </p>
              </div>
            )}
          </div>

          {/* Optimization Settings */}
          <div>
            <h3 className="text-lg font-quando font-semibold mb-3">
              Optimization Settings
            </h3>
            <div className="grid grid-cols-2 gap-4">
              {Object.entries(optimization).map(([key, value]) => (
                <label key={key} className="flex items-center cursor-pointer">
                  <input
                    type="checkbox"
                    checked={value}
                    onChange={(e) =>
                      setOptimization({
                        ...optimization,
                        [key]: e.target.checked,
                      })
                    }
                    className="mr-3"
                  />
                  <span className="font-quando text-sm">
                    {key
                      .replace(/([A-Z])/g, " $1")
                      .replace(/^./, (str) => str.toUpperCase())}
                    {" optimization"}
                  </span>
                </label>
              ))}
            </div>
          </div>

          {/* Analysis Configuration */}
          {/* <div className="bg-card border border-border rounded-lg p-6">
            <h3 className="text-lg font-quando font-semibold mb-3">
              Analysis Options
            </h3>
            <p className="text-sm text-muted-foreground mb-4">
              Select which analyses to perform after computation
            </p>
            <div className="grid grid-cols-2 gap-3">
              <label className="flex items-center cursor-pointer">
                <input
                  type="checkbox"
                  checked={analysis.energy_decomposition}
                  onChange={(e) =>
                    setAnalysis({ ...analysis, energy_decomposition: e.target.checked })
                  }
                  className="mr-3"
                />
                <span className="font-quando text-sm">Energy Decomposition</span>
              </label>
              <label className="flex items-center cursor-pointer">
                <input
                  type="checkbox"
                  checked={analysis.bond_analysis}
                  onChange={(e) =>
                    setAnalysis({ ...analysis, bond_analysis: e.target.checked })
                  }
                  className="mr-3"
                />
                <span className="font-quando text-sm">Bond Analysis</span>
              </label>
              <label className="flex items-center cursor-pointer">
                <input
                  type="checkbox"
                  checked={analysis.dipole_moment}
                  onChange={(e) =>
                    setAnalysis({ ...analysis, dipole_moment: e.target.checked })
                  }
                  className="mr-3"
                />
                <span className="font-quando text-sm">Dipole Moment</span>
              </label>
              <label className="flex items-center cursor-pointer">
                <input
                  type="checkbox"
                  checked={analysis.polarizability}
                  onChange={(e) =>
                    setAnalysis({ ...analysis, polarizability: e.target.checked })
                  }
                  className="mr-3"
                />
                <span className="font-quando text-sm">Polarizability</span>
              </label>
              <label className="flex items-center cursor-pointer">
                <input
                  type="checkbox"
                  checked={analysis.thermochemistry}
                  onChange={(e) =>
                    setAnalysis({ ...analysis, thermochemistry: e.target.checked })
                  }
                  className="mr-3"
                />
                <span className="font-quando text-sm">Thermochemistry</span>
              </label>
              <label className="flex items-center cursor-pointer">
                <input
                  type="checkbox"
                  checked={analysis.spectroscopy}
                  onChange={(e) =>
                    setAnalysis({ ...analysis, spectroscopy: e.target.checked })
                  }
                  className="mr-3"
                />
                <span className="font-quando text-sm">Spectroscopy</span>
              </label>
              <label className="flex items-center cursor-pointer">
                <input
                  type="checkbox"
                  checked={analysis.vibrational}
                  onChange={(e) =>
                    setAnalysis({ ...analysis, vibrational: e.target.checked })
                  }
                  className="mr-3"
                />
                <span className="font-quando text-sm">Vibrational Analysis</span>
              </label>
            </div>
          </div> */}
        </div>

        {/* Footer */}
        <div className="border-t border-border p-6 flex gap-4">
          <button
            onClick={onClose}
            className="flex-1 px-6 py-3 border border-border rounded-lg hover:bg-accent transition font-quando"
          >
            Cancel
          </button>
          <button
            onClick={handleSave}
            className="flex-1 px-6 py-3 bg-brand-orange text-white rounded-lg hover:bg-brand-orange-dark transition font-quando flex items-center justify-center gap-2"
          >
            <Save className="w-5 h-5" />
            Save Settings
          </button>
        </div>
      </div>
    </div>
  );
}
