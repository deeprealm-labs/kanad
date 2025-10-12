"use client";

import { useState, useEffect } from "react";
import { X, Save } from "lucide-react";
import * as api from "@/lib/api";

interface SettingsModalProps {
  isOpen: boolean;
  onClose: () => void;
  onSave?: () => void;
}

const DEFAULT_SETTINGS = {
  method: "VQE",
  ansatz: "hardware_efficient",
  mapper: "jordan_wigner",
  optimizer: "SLSQP",
  backend: "classical",
  backendName: "ibm_torino",
  optimization: {
    geometry: false,
    orbitals: false,
    circuit: true,
    adaptive: false,
  },
};

export default function SettingsModal({ isOpen, onClose, onSave }: SettingsModalProps) {
  const [method, setMethod] = useState(DEFAULT_SETTINGS.method);
  const [ansatz, setAnsatz] = useState(DEFAULT_SETTINGS.ansatz);
  const [mapper, setMapper] = useState(DEFAULT_SETTINGS.mapper);
  const [optimizer, setOptimizer] = useState(DEFAULT_SETTINGS.optimizer);
  const [backend, setBackend] = useState(DEFAULT_SETTINGS.backend);
  const [backendName, setBackendName] = useState(DEFAULT_SETTINGS.backendName);
  const [hamiltonian, setHamiltonian] = useState("molecular");
  const [basisSet, setBasisSet] = useState("sto-3g");

  const [optimization, setOptimization] = useState(DEFAULT_SETTINGS.optimization);
  const [configOptions, setConfigOptions] = useState<any>(null);

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
        setHamiltonian(settings.hamiltonian || "molecular");
        setBasisSet(settings.basisSet || "sto-3g");
        setOptimization(settings.optimization || DEFAULT_SETTINGS.optimization);
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
            setHamiltonian(settings.hamiltonian || "molecular");
            setBasisSet(settings.basisSet || "sto-3g");
            setOptimization(settings.optimization || DEFAULT_SETTINGS.optimization);
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
      hamiltonian,
      basisSet,
      optimization,
    };

    try {
      await api.updateSettings(settings);
      console.log("Saved settings to API:", settings);
    } catch (error) {
      console.error("Failed to save settings to API:", error);
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
          {/* Computation Method */}
          <div>
            <h3 className="text-lg font-quando font-semibold mb-3">
              Computation Method
            </h3>
            <select
              value={method}
              onChange={(e) => setMethod(e.target.value)}
              className="w-full px-4 py-3 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
            >
              {configOptions?.methods?.map((m: any) => (
                <option key={m.value} value={m.value}>
                  {m.label} - {m.description}
                </option>
              )) || (
                <>
                  <option value="HF">Hartree-Fock (HF)</option>
                  <option value="VQE">Variational Quantum Eigensolver (VQE)</option>
                  <option value="SQD">Subspace Quantum Diagonalization</option>
                  <option value="EXCITED_STATES">Excited States</option>
                </>
              )}
            </select>
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
                        <option value="SLSQP">SLSQP</option>
                        <option value="COBYLA">COBYLA</option>
                        <option value="L-BFGS-B">L-BFGS-B</option>
                      </>
                    )}
                  </select>
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
                    IBM Quantum (Cloud)
                  </div>
                  <div className="text-xs text-muted-foreground mt-1">
                    Run on real quantum hardware
                  </div>
                  {backend === "ibm_quantum" && (
                    <select
                      value={backendName}
                      onChange={(e) => setBackendName(e.target.value)}
                      className="mt-2 w-full px-3 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando text-sm"
                      onClick={(e) => e.stopPropagation()}
                    >
                      {configOptions?.ibm_backends?.map((b: any) => (
                        <option key={b.value} value={b.value}>
                          {b.label} ({b.qubits} qubits)
                        </option>
                      )) || (
                        <>
                          <option value="ibm_torino">ibm_torino (133 qubits)</option>
                          <option value="ibm_brisbane">ibm_brisbane (127 qubits)</option>
                          <option value="ibm_kyoto">ibm_kyoto (127 qubits)</option>
                        </>
                      )}
                    </select>
                  )}
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
                    10x faster simulation on GPUs
                  </div>
                </div>
              </label>
            </div>
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
