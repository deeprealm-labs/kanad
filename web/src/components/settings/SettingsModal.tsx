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
  // Method: VQE, SQD, KRYLOV_SQD, EXCITED_STATES, HF
  const [method, setMethod] = useState("VQE");

  // VQE Mode: standard or hivqe
  const [vqeMode, setVqeMode] = useState("standard");

  // Application-based configuration
  const [application, setApplication] = useState<string>("quantum-chemistry");
  const [applicationConfig, setApplicationConfig] = useState<any>({});

  // Common settings
  const [ansatz, setAnsatz] = useState(DEFAULT_SETTINGS.ansatz);
  const [mapper, setMapper] = useState(DEFAULT_SETTINGS.mapper);
  const [optimizer, setOptimizer] = useState(DEFAULT_SETTINGS.optimizer);
  const [backend, setBackend] = useState(DEFAULT_SETTINGS.backend);
  const [backendName, setBackendName] = useState(DEFAULT_SETTINGS.backendName);
  const [bluequbitDevice, setBluequbitDevice] = useState("gpu");
  const [hamiltonian, setHamiltonian] = useState("molecular");
  const [maxIterations, setMaxIterations] = useState(DEFAULT_SETTINGS.maxIterations);

  // VQE Advanced Options
  const [useActiveSpace, setUseActiveSpace] = useState(false);
  const [hivqeMaxIterations, setHivqeMaxIterations] = useState(10);
  const [hivqeSubspaceThreshold, setHivqeSubspaceThreshold] = useState(0.05);

  // SQD-specific settings
  const [subspaceDim, setSubspaceDim] = useState(10);
  const [circuitDepth, setCircuitDepth] = useState(3);
  const [nStates, setNStates] = useState(3);

  // Krylov-SQD specific settings
  const [krylovDim, setKrylovDim] = useState(15);

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
        setVqeMode(settings.vqe_mode || "standard");
        setAnsatz(settings.ansatz || DEFAULT_SETTINGS.ansatz);
        setMapper(settings.mapper || DEFAULT_SETTINGS.mapper);
        setOptimizer(settings.optimizer || DEFAULT_SETTINGS.optimizer);
        setBackend(settings.backend || DEFAULT_SETTINGS.backend);
        setBackendName(settings.backendName || DEFAULT_SETTINGS.backendName);
        setBluequbitDevice(settings.bluequbitDevice || "gpu");
        setHamiltonian(settings.hamiltonian || "molecular");
        setMaxIterations(settings.maxIterations || DEFAULT_SETTINGS.maxIterations);
        setUseActiveSpace(settings.use_active_space || false);
        setHivqeMaxIterations(settings.hivqe_max_iterations || 10);
        setHivqeSubspaceThreshold(settings.hivqe_subspace_threshold || 0.05);
        setSubspaceDim(settings.subspaceDim || 10);
        setCircuitDepth(settings.circuitDepth || 3);
        setNStates(settings.nStates || 3);
        setKrylovDim(settings.krylov_dim || 15);
        setOptimization(settings.optimization || DEFAULT_SETTINGS.optimization);
        setAnalysis(settings.analysis || analysis);

        // Load application configuration
        setApplication(settings.application || "quantum-chemistry");
        setApplicationConfig(settings.applicationConfig || {});
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
            setMaxIterations(settings.maxIterations || DEFAULT_SETTINGS.maxIterations);
            setSubspaceDim(settings.subspaceDim || 10);
            setCircuitDepth(settings.circuitDepth || 3);
            setNStates(settings.nStates || 3);
            setOptimization(settings.optimization || DEFAULT_SETTINGS.optimization);
            setAnalysis(settings.analysis || analysis);
            setApplication(settings.application || "quantum-chemistry");
            setApplicationConfig(settings.applicationConfig || {});
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
      vqe_mode: vqeMode,
      ansatz,
      mapper,
      optimizer,
      backend,
      backendName,
      bluequbitDevice,
      hamiltonian,
      maxIterations,
      use_active_space: useActiveSpace,
      hivqe_max_iterations: hivqeMaxIterations,
      hivqe_subspace_threshold: hivqeSubspaceThreshold,
      subspaceDim,
      circuitDepth,
      nStates,
      krylov_dim: krylovDim,
      optimization,
      analysis,
      application,
      applicationConfig,
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
      <div className="bg-background border border-border rounded-lg w-full max-w-6xl max-h-[90vh] overflow-hidden flex flex-col">
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
          {/* Application Focus - NOW AT THE TOP */}
          <div className="bg-gradient-to-br from-brand-orange/5 to-brand-orange/10 border-2 border-brand-orange/30 rounded-xl p-6">
            <div className="flex items-center gap-3 mb-4">
              <div className="w-10 h-10 rounded-lg bg-brand-orange/20 flex items-center justify-center">
                <svg className="w-6 h-6 text-brand-orange" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 12l2 2 4-4m5.618-4.016A11.955 11.955 0 0112 2.944a11.955 11.955 0 01-8.618 3.04A12.02 12.02 0 003 9c0 5.591 3.824 10.29 9 11.622 5.176-1.332 9-6.03 9-11.622 0-1.042-.133-2.052-.382-3.016z" />
                </svg>
              </div>
              <div>
                <h3 className="text-xl font-quando font-bold text-brand-orange">
                  Application Domain
                </h3>
                <p className="text-xs text-muted-foreground mt-0.5">
                  Configure your research focus and analysis priorities
                </p>
              </div>
            </div>

            <div className="space-y-4">
              <div>
                <label className="block text-sm font-quando font-semibold mb-2 text-foreground">
                  Select Your Research Domain
                </label>
                <select
                  value={application}
                  onChange={(e) => setApplication(e.target.value)}
                  className="w-full px-4 py-3 border-2 border-brand-orange/30 bg-background rounded-lg focus:outline-none focus:ring-2 focus:ring-brand-orange focus:border-brand-orange font-quando text-base shadow-sm hover:border-brand-orange/50 transition-colors"
                >
                  <option value="quantum-chemistry">Quantum Chemistry - General molecular calculations</option>
                  <option value="drug-discovery">Drug Discovery - ADME, pharmacokinetics, binding affinity</option>
                  <option value="materials-science">Materials Science - Band structure, density of states</option>
                  <option value="catalysis">Catalysis - Reaction mechanisms, transition states</option>
                  <option value="spectroscopy">Spectroscopy - UV-Vis, NMR, excited states</option>
                  <option value="molecular-dynamics">Molecular Dynamics - Time-dependent simulations</option>
                </select>
              </div>
            </div>
          </div>

          {/* Application Workflow Configuration */}
          <div>
            <h3 className="text-lg font-quando font-semibold mb-3">
              Application Workflow
            </h3>
            <p className="text-sm text-muted-foreground mb-4">
              Configure your experiment workflow based on real-world use cases
            </p>

            <div className="space-y-4">
              {/* Drug Discovery Workflow */}
              {application === "drug-discovery" && (
                <div className="space-y-4">
                  {/* Use Case Selection */}
                  <div className="bg-white dark:bg-card border border-brand-orange/30 rounded-xl p-5 shadow-sm">
                    <div className="flex items-center gap-2 mb-3">
                      <svg className="w-5 h-5 text-brand-orange" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 5H7a2 2 0 00-2 2v12a2 2 0 002 2h10a2 2 0 002-2V7a2 2 0 00-2-2h-2M9 5a2 2 0 002 2h2a2 2 0 002-2M9 5a2 2 0 012-2h2a2 2 0 012 2" />
                      </svg>
                      <h4 className="text-base font-quando font-bold text-brand-orange">Drug Discovery Use Case</h4>
                    </div>
                    <div>
                      <label className="block text-sm font-quando font-medium mb-2">Select Research Goal</label>
                      <select
                        value={applicationConfig?.useCase || "lead-optimization"}
                        onChange={(e) => setApplicationConfig({...applicationConfig, useCase: e.target.value})}
                        className="w-full px-4 py-3 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                      >
                        <option value="lead-optimization">Lead Optimization - Improve drug candidates</option>
                        <option value="target-screening">Target Screening - Find binding molecules</option>
                        <option value="adme-prediction">ADME Prediction - Pharmacokinetic properties</option>
                        <option value="toxicity-assessment">Toxicity Assessment - Safety profiling</option>
                      </select>
                    </div>
                  </div>

                  {/* Workflow Parameters */}
                  <div className="bg-white dark:bg-card border border-brand-orange/30 rounded-xl p-5 shadow-sm">
                    <h5 className="text-sm font-quando font-bold text-brand-orange mb-3">Workflow Parameters</h5>
                    <div className="grid grid-cols-2 gap-4">
                      <div>
                        <label className="block text-sm font-quando font-medium mb-2">Target Protein (Optional)</label>
                        <input
                          type="text"
                          placeholder="e.g., EGFR, ACE2"
                          value={applicationConfig?.targetProtein || ""}
                          onChange={(e) => setApplicationConfig({...applicationConfig, targetProtein: e.target.value})}
                          className="w-full px-3 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando text-sm"
                        />
                      </div>
                      <div>
                        <label className="block text-sm font-quando font-medium mb-2">Binding Affinity Threshold</label>
                        <input
                          type="number"
                          placeholder="-8.0"
                          step="0.1"
                          value={applicationConfig?.bindingThreshold || ""}
                          onChange={(e) => setApplicationConfig({...applicationConfig, bindingThreshold: parseFloat(e.target.value)})}
                          className="w-full px-3 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quanto text-sm"
                        />
                        <p className="text-xs text-muted-foreground mt-1">kcal/mol (more negative = stronger binding)</p>
                      </div>
                    </div>
                  </div>

                  {/* Analysis Pipeline */}
                  <div className="bg-white dark:bg-card border border-brand-orange/30 rounded-xl p-5 shadow-sm">
                    <h5 className="text-sm font-quando font-bold text-brand-orange mb-3">Analysis Pipeline</h5>
                    <div className="space-y-2">
                      {[
                        { key: 'molecularDocking', label: 'Molecular Docking', desc: 'Protein-ligand binding poses' },
                        { key: 'admeAnalysis', label: 'ADME Analysis', desc: 'Drug-like properties assessment' },
                        { key: 'toxicityPrediction', label: 'Toxicity Prediction', desc: 'Safety and side effects' }
                      ].map(({ key, label, desc }) => (
                        <label key={key} className="flex items-center p-3 bg-muted/30 rounded-lg hover:bg-muted/50 transition-colors cursor-pointer">
                          <input
                            type="checkbox"
                            checked={applicationConfig?.[key] !== false}
                            onChange={(e) => setApplicationConfig({...applicationConfig, [key]: e.target.checked})}
                            className="mr-3 w-4 h-4 text-brand-orange focus:ring-brand-orange rounded"
                          />
                          <div className="flex-1">
                            <div className="font-quando font-medium text-sm">{label}</div>
                            <div className="text-xs text-muted-foreground">{desc}</div>
                          </div>
                        </label>
                      ))}
                    </div>
                  </div>
                </div>
              )}

              {/* Materials Science Workflow */}
              {application === "materials-science" && (
                <div className="space-y-4">
                  <div className="bg-white dark:bg-card border border-brand-orange/30 rounded-xl p-5 shadow-sm">
                    <div className="flex items-center gap-2 mb-3">
                      <svg className="w-5 h-5 text-brand-orange" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 3v2m6-2v2M9 19v2m6-2v2M5 9H3m2 6H3m18-6h-2m2 6h-2M7 19h10a2 2 0 002-2V7a2 2 0 00-2-2H7a2 2 0 00-2 2v10a2 2 0 002 2zM9 9h6v6H9V9z" />
                      </svg>
                      <h4 className="text-base font-quando font-bold text-brand-orange">Materials Science Use Case</h4>
                    </div>
                    <div>
                      <label className="block text-sm font-quando font-medium mb-2">Select Research Goal</label>
                      <select
                        value={applicationConfig?.useCase || "electronic-structure"}
                        onChange={(e) => setApplicationConfig({...applicationConfig, useCase: e.target.value})}
                        className="w-full px-4 py-3 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                      >
                        <option value="electronic-structure">Electronic Structure - Band gaps, DOS</option>
                        <option value="optical-properties">Optical Properties - Absorption spectra</option>
                        <option value="mechanical-properties">Mechanical Properties - Elastic constants</option>
                        <option value="thermal-properties">Thermal Properties - Phonon dispersion</option>
                      </select>
                    </div>
                  </div>

                  <div className="bg-white dark:bg-card border border-brand-orange/30 rounded-xl p-5 shadow-sm">
                    <h5 className="text-sm font-quando font-bold text-brand-orange mb-3">Material Properties</h5>
                    <div className="grid grid-cols-2 gap-4">
                      <div>
                        <label className="block text-sm font-quando font-medium mb-2">Target Band Gap (eV)</label>
                        <input
                          type="number"
                          placeholder="1.5"
                          step="0.1"
                          value={applicationConfig?.targetBandGap || ""}
                          onChange={(e) => setApplicationConfig({...applicationConfig, targetBandGap: parseFloat(e.target.value)})}
                          className="w-full px-3 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando text-sm"
                        />
                      </div>
                      <div>
                        <label className="block text-sm font-quando font-medium mb-2">Crystal System</label>
                        <select
                          value={applicationConfig?.crystalSystem || "cubic"}
                          onChange={(e) => setApplicationConfig({...applicationConfig, crystalSystem: e.target.value})}
                          className="w-full px-3 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando text-sm"
                        >
                          <option value="cubic">Cubic</option>
                          <option value="hexagonal">Hexagonal</option>
                          <option value="tetragonal">Tetragonal</option>
                          <option value="orthorhombic">Orthorhombic</option>
                        </select>
                      </div>
                    </div>
                  </div>
                </div>
              )}

              {/* Catalysis Workflow */}
              {application === "catalysis" && (
                <div className="space-y-4">
                  <div className="bg-white dark:bg-card border border-brand-orange/30 rounded-xl p-5 shadow-sm">
                    <div className="flex items-center gap-2 mb-3">
                      <svg className="w-5 h-5 text-brand-orange" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M13 10V3L4 14h7v7l9-11h-7z" />
                      </svg>
                      <h4 className="text-base font-quando font-bold text-brand-orange">Catalysis Use Case</h4>
                    </div>
                    <div>
                      <label className="block text-sm font-quando font-medium mb-2">Select Research Goal</label>
                      <select
                        value={applicationConfig?.useCase || "reaction-mechanism"}
                        onChange={(e) => setApplicationConfig({...applicationConfig, useCase: e.target.value})}
                        className="w-full px-4 py-3 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                      >
                        <option value="reaction-mechanism">Reaction Mechanism - Pathways and intermediates</option>
                        <option value="catalyst-screening">Catalyst Screening - Find optimal catalysts</option>
                        <option value="activation-energy">Activation Energy - Barrier heights</option>
                        <option value="selectivity-analysis">Selectivity Analysis - Product distribution</option>
                      </select>
                    </div>
                  </div>

                  <div className="bg-white dark:bg-card border border-brand-orange/30 rounded-xl p-5 shadow-sm">
                    <h5 className="text-sm font-quando font-bold text-brand-orange mb-3">Reaction Conditions</h5>
                    <div className="grid grid-cols-2 gap-4">
                      <div>
                        <label className="block text-sm font-quando font-medium mb-2">Temperature (K)</label>
                        <input
                          type="number"
                          placeholder="298"
                          value={applicationConfig?.temperature || ""}
                          onChange={(e) => setApplicationConfig({...applicationConfig, temperature: parseInt(e.target.value)})}
                          className="w-full px-3 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando text-sm"
                        />
                      </div>
                      <div>
                        <label className="block text-sm font-quando font-medium mb-2">Pressure (atm)</label>
                        <input
                          type="number"
                          placeholder="1.0"
                          step="0.1"
                          value={applicationConfig?.pressure || ""}
                          onChange={(e) => setApplicationConfig({...applicationConfig, pressure: parseFloat(e.target.value)})}
                          className="w-full px-3 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quanto text-sm"
                        />
                      </div>
                    </div>
                  </div>
                </div>
              )}

              {/* Spectroscopy Workflow */}
              {application === "spectroscopy" && (
                <div className="space-y-4">
                  <div className="bg-white dark:bg-card border border-brand-orange/30 rounded-xl p-5 shadow-sm">
                    <div className="flex items-center gap-2 mb-3">
                      <svg className="w-5 h-5 text-brand-orange" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M7 12l3-3 3 3 4-4M8 21l4-4 4 4M3 4h18M4 4h16v12a1 1 0 01-1 1H5a1 1 0 01-1-1V4z" />
                      </svg>
                      <h4 className="text-base font-quando font-bold text-brand-orange">Spectroscopy Use Case</h4>
                    </div>
                    <div>
                      <label className="block text-sm font-quando font-medium mb-2">Select Spectroscopy Type</label>
                      <select
                        value={applicationConfig?.useCase || "uv-vis"}
                        onChange={(e) => setApplicationConfig({...applicationConfig, useCase: e.target.value})}
                        className="w-full px-4 py-3 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                      >
                        <option value="uv-vis">UV-Vis Spectroscopy - Electronic transitions</option>
                        <option value="nmr">NMR Spectroscopy - Chemical shifts</option>
                        <option value="ir-raman">IR/Raman - Vibrational modes</option>
                        <option value="excited-states">Excited States - TD-DFT analysis</option>
                      </select>
                    </div>
                  </div>

                  <div className="bg-white dark:bg-card border border-brand-orange/30 rounded-xl p-5 shadow-sm">
                    <h5 className="text-sm font-quando font-bold text-brand-orange mb-3">Spectral Parameters</h5>
                    <div className="grid grid-cols-2 gap-4">
                      <div>
                        <label className="block text-sm font-quando font-medium mb-2">Wavelength Range (nm)</label>
                        <div className="flex gap-2">
                          <input
                            type="number"
                            placeholder="200"
                            value={applicationConfig?.wavelengthMin || ""}
                            onChange={(e) => setApplicationConfig({...applicationConfig, wavelengthMin: parseInt(e.target.value)})}
                            className="w-1/2 px-3 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando text-sm"
                          />
                          <input
                            type="number"
                            placeholder="800"
                            value={applicationConfig?.wavelengthMax || ""}
                            onChange={(e) => setApplicationConfig({...applicationConfig, wavelengthMax: parseInt(e.target.value)})}
                            className="w-1/2 px-3 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando text-sm"
                          />
                        </div>
                      </div>
                      <div>
                        <label className="block text-sm font-quando font-medium mb-2">Number of Excited States</label>
                        <input
                          type="number"
                          placeholder="10"
                          value={applicationConfig?.nExcitedStates || ""}
                          onChange={(e) => setApplicationConfig({...applicationConfig, nExcitedStates: parseInt(e.target.value)})}
                          className="w-full px-3 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando text-sm"
                        />
                      </div>
                    </div>
                  </div>
                </div>
              )}

              {/* Quantum Chemistry - General */}
              {application === "quantum-chemistry" && (
                <div className="bg-white dark:bg-card border border-brand-orange/30 rounded-xl p-5 shadow-sm">
                  <div className="flex items-center gap-2 mb-3">
                    <svg className="w-5 h-5 text-brand-orange" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                      <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9.663 17h4.673M12 3v1m6.364 1.636l-.707.707M21 12h-1M4 12H3m3.343-5.657l-.707-.707m2.828 9.9a5 5 0 117.072 0l-.548.547A3.374 3.374 0 0014 18.469V19a2 2 0 11-4 0v-.531c0-.895-.356-1.754-.988-2.386l-.548-.547z" />
                    </svg>
                    <h4 className="text-base font-quando font-bold text-brand-orange">General Quantum Chemistry</h4>
                  </div>
                  <p className="text-sm text-muted-foreground">
                    Standard molecular calculations for ground state energy, geometry optimization, and molecular properties.
                    Configure specific application domains above for specialized workflows.
                  </p>
                </div>
              )}

              {/* Molecular Dynamics */}
              {application === "molecular-dynamics" && (
                <div className="space-y-4">
                  <div className="bg-white dark:bg-card border border-brand-orange/30 rounded-xl p-5 shadow-sm">
                    <div className="flex items-center gap-2 mb-3">
                      <svg className="w-5 h-5 text-brand-orange" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                        <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M13 10V3L4 14h7v7l9-11h-7z" />
                      </svg>
                      <h4 className="text-base font-quando font-bold text-brand-orange">Molecular Dynamics Simulation</h4>
                    </div>
                    <div className="grid grid-cols-2 gap-4">
                      <div>
                        <label className="block text-sm font-quando font-medium mb-2">Simulation Time (ps)</label>
                        <input
                          type="number"
                          placeholder="1000"
                          value={applicationConfig?.simulationTime || ""}
                          onChange={(e) => setApplicationConfig({...applicationConfig, simulationTime: parseInt(e.target.value)})}
                          className="w-full px-3 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando text-sm"
                        />
                      </div>
                      <div>
                        <label className="block text-sm font-quando font-medium mb-2">Temperature (K)</label>
                        <input
                          type="number"
                          placeholder="300"
                          value={applicationConfig?.temperature || ""}
                          onChange={(e) => setApplicationConfig({...applicationConfig, temperature: parseInt(e.target.value)})}
                          className="w-full px-3 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando text-sm"
                        />
                      </div>
                    </div>
                  </div>
                </div>
              )}
            </div>
          </div>

          
          {/* Solver Method Selection */}
          <div>
            <h3 className="text-lg font-quando font-semibold mb-3">
              Solver Method
            </h3>
            <div className="grid grid-cols-2 gap-3">
              {configOptions?.methods?.map((methodOption: any) => (
                <label key={methodOption.value} className="flex items-center p-4 border-2 border-border rounded-lg hover:bg-accent cursor-pointer transition">
                  <input
                    type="radio"
                    name="method"
                    value={methodOption.value}
                    checked={method === methodOption.value}
                    onChange={(e) => setMethod(e.target.value)}
                    className="mr-3"
                  />
                  <div>
                    <div className="font-quando font-semibold">{methodOption.label}</div>
                    <div className="text-xs text-muted-foreground">{methodOption.description}</div>
                  </div>
                </label>
              )) || (
                /* Fallback if config not loaded */
                <>
                  <label className="flex items-center p-4 border-2 border-border rounded-lg hover:bg-accent cursor-pointer transition">
                    <input
                      type="radio"
                      name="method"
                      value="VQE"
                      checked={method === "VQE"}
                      onChange={(e) => setMethod(e.target.value)}
                      className="mr-3"
                    />
                    <div>
                      <div className="font-quando font-semibold">VQE</div>
                      <div className="text-xs text-muted-foreground">Variational Quantum Eigensolver</div>
                    </div>
                  </label>
                  <label className="flex items-center p-4 border-2 border-border rounded-lg hover:bg-accent cursor-pointer transition">
                    <input
                      type="radio"
                      name="method"
                      value="SQD"
                      checked={method === "SQD"}
                      onChange={(e) => setMethod(e.target.value)}
                      className="mr-3"
                    />
                    <div>
                      <div className="font-quando font-semibold">SQD</div>
                      <div className="text-xs text-muted-foreground">Subspace Quantum Diagonalization</div>
                    </div>
                  </label>
                </>
              )}
            </div>
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
                    VQE Mode
                  </label>
                  <select
                    value={vqeMode}
                    onChange={(e) => setVqeMode(e.target.value)}
                    className="w-full px-4 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                  >
                    {configOptions?.vqe_modes?.map((mode: any) => (
                      <option key={mode.value} value={mode.value}>
                        {mode.label}
                      </option>
                    )) || (
                      <>
                        <option value="standard">Standard VQE</option>
                        <option value="hivqe">Hi-VQE (1000x measurement reduction)</option>
                      </>
                    )}
                  </select>
                  {vqeMode === "hivqe" && backend !== "classical" && (
                    <p className="text-xs text-green-600 dark:text-green-400 mt-1">
                      🚀 Hi-VQE: 99.98% cost savings on quantum hardware (1000x fewer measurements)
                    </p>
                  )}
                  {vqeMode === "hivqe" && backend === "classical" && (
                    <p className="text-xs text-yellow-600 dark:text-yellow-400 mt-1">
                      ℹ️ Hi-VQE provides no benefit on classical simulator (use for real hardware)
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

              {/* Advanced VQE Options */}
              <div className="border-t border-border pt-4">
                <h4 className="text-md font-quando font-semibold mb-3">Advanced Options</h4>

                <div className="space-y-4">
                  {/* Active Space Reduction */}
                  <label className="flex items-center p-3 border border-border rounded-lg hover:bg-accent cursor-pointer transition">
                    <input
                      type="checkbox"
                      checked={useActiveSpace}
                      onChange={(e) => setUseActiveSpace(e.target.checked)}
                      className="mr-3"
                    />
                    <div className="flex-1">
                      <div className="font-quando font-medium">Active Space Reduction</div>
                      <div className="text-xs text-muted-foreground">
                        Freeze core orbitals to reduce qubit count by ~17% (better for large molecules)
                      </div>
                    </div>
                  </label>

                  {/* Hi-VQE Specific Options */}
                  {vqeMode === "hivqe" && (
                    <div className="bg-green-50 dark:bg-green-900/20 border border-green-200 dark:border-green-800 rounded-lg p-4 space-y-3">
                      <div className="font-quando font-medium text-green-800 dark:text-green-200">
                        Hi-VQE Advanced Settings
                      </div>

                      <div>
                        <label className="block text-sm font-quando font-medium mb-2">
                          Subspace Iterations
                        </label>
                        <input
                          type="number"
                          min="5"
                          max="20"
                          value={hivqeMaxIterations}
                          onChange={(e) => setHivqeMaxIterations(parseInt(e.target.value) || 10)}
                          className="w-full px-4 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                        />
                        <p className="text-xs text-muted-foreground mt-1">
                          Number of subspace expansion iterations (5-20, default: 10)
                        </p>
                      </div>

                      <div>
                        <label className="block text-sm font-quando font-medium mb-2">
                          Amplitude Threshold
                        </label>
                        <input
                          type="number"
                          min="0.01"
                          max="0.1"
                          step="0.01"
                          value={hivqeSubspaceThreshold}
                          onChange={(e) => setHivqeSubspaceThreshold(parseFloat(e.target.value) || 0.05)}
                          className="w-full px-4 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                        />
                        <p className="text-xs text-muted-foreground mt-1">
                          Configuration importance threshold (0.01-0.1, default: 0.05)
                        </p>
                      </div>
                    </div>
                  )}
                </div>
              </div>
            </div>
          )}

          {/* Krylov-SQD Settings */}
          {method === "KRYLOV_SQD" && (
            <div className="bg-card border border-border rounded-lg p-6 space-y-4">
              <h3 className="text-lg font-quando font-semibold">
                Krylov-SQD Configuration
              </h3>
              <p className="text-sm text-muted-foreground">
                Krylov Subspace Quantum Diagonalization - 10-20x more efficient than standard SQD. Computes ground + excited states.
              </p>

              <div className="grid grid-cols-2 gap-4">
                <div>
                  <label className="block text-sm font-quando font-medium mb-2">
                    Krylov Dimension
                  </label>
                  <input
                    type="number"
                    min="10"
                    max="30"
                    value={krylovDim}
                    onChange={(e) => setKrylovDim(parseInt(e.target.value) || 15)}
                    className="w-full px-4 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                  />
                  <p className="text-xs text-muted-foreground mt-1">
                    Subspace dimension (10-30, default: 15). Much smaller than standard SQD (50-100).
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
                    value={nStates}
                    onChange={(e) => setNStates(parseInt(e.target.value) || 3)}
                    className="w-full px-4 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                  />
                  <p className="text-xs text-muted-foreground mt-1">
                    Number of eigenvalues to compute (ground + excited states)
                  </p>
                </div>
              </div>

              <div className="bg-blue-50 dark:bg-blue-900/20 border border-blue-200 dark:border-blue-800 rounded-lg p-3">
                <p className="text-sm text-blue-800 dark:text-blue-200">
                  ⚠️ Note: Krylov-SQD currently only supports diatomic molecules (2 atoms). Use VQE for larger systems.
                </p>
              </div>
            </div>
          )}

          {/* EXCITED_STATES Settings */}
          {method === "EXCITED_STATES" && (
            <div className="bg-card border border-border rounded-lg p-6 space-y-4">
              <h3 className="text-lg font-quando font-semibold">
                Excited States Configuration
              </h3>
              <p className="text-sm text-muted-foreground">
                Calculate excited state energies and properties for electronic excitations.
              </p>

              <div className="grid grid-cols-2 gap-4">
                <div>
                  <label className="block text-sm font-quando font-medium mb-2">
                    Number of Excited States
                  </label>
                  <input
                    type="number"
                    min="1"
                    max="10"
                    value={nStates}
                    onChange={(e) => setNStates(parseInt(e.target.value) || 3)}
                    className="w-full px-4 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                  />
                </div>
              </div>
            </div>
          )}

          {/* Old HI-VQE section removed - now integrated into VQE */}
          {method === "HI-VQE_OLD" && (
            <div className="bg-card border border-border rounded-lg p-6 space-y-4">
              <h3 className="text-lg font-quando font-semibold">
                HI-VQE Configuration
              </h3>
              <p className="text-sm text-muted-foreground">
                Hardware-Independent VQE with error mitigation - Combines VQE optimization with noise suppression techniques
              </p>

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
                    Optimizer
                  </label>
                  <select
                    value={optimizer}
                    onChange={(e) => setOptimizer(e.target.value)}
                    className="w-full px-4 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                  >
                    <option value="COBYLA">COBYLA (Recommended)</option>
                    <option value="POWELL">POWELL</option>
                    <option value="SLSQP">SLSQP</option>
                  </select>
                </div>

                <div>
                  <label className="block text-sm font-quando font-medium mb-2">
                    Error Mitigation Level
                  </label>
                  <select
                    value={applicationConfig?.errorMitigation || "standard"}
                    onChange={(e) => setApplicationConfig({...applicationConfig, errorMitigation: e.target.value})}
                    className="w-full px-4 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                  >
                    <option value="none">None</option>
                    <option value="standard">Standard (Readout Error)</option>
                    <option value="advanced">Advanced (Full ZNE)</option>
                  </select>
                  <p className="text-xs text-muted-foreground mt-1">
                    Higher levels provide better accuracy but increase computation time
                  </p>
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
                  <label className="block text-sm font-quando font-medium mb-2">
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

          {/* Error Mitigation - Only for Quantum Backends */}
          {backend !== "classical" && (
            <div className="bg-card border border-border rounded-lg p-6 space-y-4">
              <h3 className="text-lg font-quando font-semibold">
                Error Mitigation
              </h3>
              <p className="text-sm text-muted-foreground">
                Configure error mitigation strategies for quantum hardware execution
              </p>

              <div className="grid grid-cols-2 gap-4">
                <div>
                  <label className="block text-sm font-quando font-medium mb-2">
                    Mitigation Level
                  </label>
                  <select
                    value={applicationConfig?.errorMitigationLevel || "standard"}
                    onChange={(e) => setApplicationConfig({...applicationConfig, errorMitigationLevel: e.target.value})}
                    className="w-full px-4 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                  >
                    <option value="none">None - No mitigation</option>
                    <option value="standard">Standard - Basic error suppression</option>
                    <option value="advanced">Advanced - Full error mitigation suite</option>
                  </select>
                  <p className="text-xs text-muted-foreground mt-1">
                    Higher levels provide better accuracy but require more circuit executions
                  </p>
                </div>

                <div>
                  <label className="flex items-center cursor-pointer p-4 border border-border rounded-lg hover:bg-accent transition">
                    <input
                      type="checkbox"
                      checked={applicationConfig?.enableZNE || false}
                      onChange={(e) => setApplicationConfig({...applicationConfig, enableZNE: e.target.checked})}
                      className="mr-3"
                    />
                    <div>
                      <span className="font-quando text-sm font-medium block">Zero-Noise Extrapolation (ZNE)</span>
                      <span className="text-xs text-muted-foreground">Extrapolate to zero-noise limit</span>
                    </div>
                  </label>
                </div>

                <div>
                  <label className="flex items-center cursor-pointer p-4 border border-border rounded-lg hover:bg-accent transition">
                    <input
                      type="checkbox"
                      checked={applicationConfig?.enableReadoutMitigation || false}
                      onChange={(e) => setApplicationConfig({...applicationConfig, enableReadoutMitigation: e.target.checked})}
                      className="mr-3"
                    />
                    <div>
                      <span className="font-quando text-sm font-medium block">Readout Error Mitigation</span>
                      <span className="text-xs text-muted-foreground">Correct measurement errors</span>
                    </div>
                  </label>
                </div>

                <div>
                  <label className="flex items-center cursor-pointer p-4 border border-border rounded-lg hover:bg-accent transition">
                    <input
                      type="checkbox"
                      checked={applicationConfig?.enableDynamicalDecoupling || false}
                      onChange={(e) => setApplicationConfig({...applicationConfig, enableDynamicalDecoupling: e.target.checked})}
                      className="mr-3"
                    />
                    <div>
                      <span className="font-quando text-sm font-medium block">Dynamical Decoupling</span>
                      <span className="text-xs text-muted-foreground">Suppress decoherence errors</span>
                    </div>
                  </label>
                </div>
              </div>
            </div>
          )}

          {/* Active Space Reduction */}
          <div className="bg-card border border-border rounded-lg p-6 space-y-4">
            <h3 className="text-lg font-quando font-semibold">
              Active Space Reduction
            </h3>
            <p className="text-sm text-muted-foreground">
              Reduce computational complexity by limiting the number of orbitals in the calculation
            </p>

            <div className="grid grid-cols-2 gap-4">
              <div>
                <label className="flex items-center cursor-pointer">
                  <input
                    type="checkbox"
                    checked={applicationConfig?.enableActiveSpace || false}
                    onChange={(e) => setApplicationConfig({...applicationConfig, enableActiveSpace: e.target.checked})}
                    className="mr-3"
                  />
                  <span className="font-quando text-sm font-medium">Enable Active Space Reduction</span>
                </label>
              </div>
            </div>

            {applicationConfig?.enableActiveSpace && (
              <div className="grid grid-cols-2 gap-4 mt-4">
                <div>
                  <label className="block text-sm font-quando font-medium mb-2">
                    Active Electrons
                  </label>
                  <input
                    type="number"
                    min="2"
                    max="20"
                    step="2"
                    value={applicationConfig?.activeElectrons || 4}
                    onChange={(e) => setApplicationConfig({...applicationConfig, activeElectrons: parseInt(e.target.value) || 4})}
                    className="w-full px-4 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                  />
                  <p className="text-xs text-muted-foreground mt-1">
                    Number of electrons in active space (must be even)
                  </p>
                </div>

                <div>
                  <label className="block text-sm font-quando font-medium mb-2">
                    Active Orbitals
                  </label>
                  <input
                    type="number"
                    min="2"
                    max="20"
                    step="1"
                    value={applicationConfig?.activeOrbitals || 4}
                    onChange={(e) => setApplicationConfig({...applicationConfig, activeOrbitals: parseInt(e.target.value) || 4})}
                    className="w-full px-4 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                  />
                  <p className="text-xs text-muted-foreground mt-1">
                    Number of orbitals in active space
                  </p>
                </div>

                <div>
                  <label className="flex items-center cursor-pointer p-4 border border-border rounded-lg hover:bg-accent transition">
                    <input
                      type="checkbox"
                      checked={applicationConfig?.frozenCore || false}
                      onChange={(e) => setApplicationConfig({...applicationConfig, frozenCore: e.target.checked})}
                      className="mr-3"
                    />
                    <div>
                      <span className="font-quando text-sm font-medium block">Frozen Core Approximation</span>
                      <span className="text-xs text-muted-foreground">Freeze core orbitals</span>
                    </div>
                  </label>
                </div>

                <div>
                  <label className="block text-sm font-quando font-medium mb-2">
                    Virtual Orbital Truncation
                  </label>
                  <input
                    type="number"
                    min="0"
                    max="50"
                    step="1"
                    value={applicationConfig?.virtualTruncation || 0}
                    onChange={(e) => setApplicationConfig({...applicationConfig, virtualTruncation: parseInt(e.target.value) || 0})}
                    className="w-full px-4 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                  />
                  <p className="text-xs text-muted-foreground mt-1">
                    Number of high-energy virtual orbitals to remove (0 = none)
                  </p>
                </div>
              </div>
            )}
          </div>

          {/* Optimization Settings */}
          <div className="bg-card border border-border rounded-lg p-6 space-y-4">
            <h3 className="text-lg font-quando font-semibold">
              Optimization Settings
            </h3>
            <p className="text-sm text-muted-foreground">
              Configure convergence criteria and optimization parameters
            </p>

            <div className="grid grid-cols-2 gap-4">
              <div>
                <label className="block text-sm font-quando font-medium mb-2">
                  Maximum Iterations
                </label>
                <input
                  type="number"
                  min="10"
                  max="500"
                  step="10"
                  value={maxIterations}
                  onChange={(e) => setMaxIterations(parseInt(e.target.value) || 100)}
                  className="w-full px-4 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                />
                <p className="text-xs text-muted-foreground mt-1">
                  Maximum optimization iterations
                </p>
              </div>

              <div>
                <label className="block text-sm font-quando font-medium mb-2">
                  Convergence Threshold
                </label>
                <input
                  type="number"
                  min="0.00001"
                  max="0.01"
                  step="0.00001"
                  value={applicationConfig?.convergenceThreshold || 0.001}
                  onChange={(e) => setApplicationConfig({...applicationConfig, convergenceThreshold: parseFloat(e.target.value) || 0.001})}
                  className="w-full px-4 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                />
                <p className="text-xs text-muted-foreground mt-1">
                  Energy convergence tolerance (Ha)
                </p>
              </div>

              <div>
                <label className="block text-sm font-quando font-medium mb-2">
                  Gradient Method
                </label>
                <select
                  value={applicationConfig?.gradientMethod || "finite-diff"}
                  onChange={(e) => setApplicationConfig({...applicationConfig, gradientMethod: e.target.value})}
                  className="w-full px-4 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                >
                  <option value="finite-diff">Finite Difference</option>
                  <option value="parameter-shift">Parameter Shift Rule</option>
                  <option value="spsa">SPSA (Stochastic)</option>
                </select>
                <p className="text-xs text-muted-foreground mt-1">
                  Method for computing gradients
                </p>
              </div>

              <div>
                <label className="flex items-center cursor-pointer p-4 border border-border rounded-lg hover:bg-accent transition">
                  <input
                    type="checkbox"
                    checked={applicationConfig?.adaptiveConvergence || false}
                    onChange={(e) => setApplicationConfig({...applicationConfig, adaptiveConvergence: e.target.checked})}
                    className="mr-3"
                  />
                  <div>
                    <span className="font-quando text-sm font-medium block">Adaptive Convergence</span>
                    <span className="text-xs text-muted-foreground">Adjust thresholds dynamically</span>
                  </div>
                </label>
              </div>
            </div>

            {/* Legacy optimization checkboxes */}
            <div className="pt-4 border-t border-border">
              <h4 className="text-sm font-quando font-medium mb-3">Additional Optimizations</h4>
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
