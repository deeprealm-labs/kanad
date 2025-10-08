"use client";

import { useState } from "react";
import { X, Save } from "lucide-react";

interface SettingsModalProps {
  isOpen: boolean;
  onClose: () => void;
}

export default function SettingsModal({ isOpen, onClose }: SettingsModalProps) {
  const [method, setMethod] = useState("VQE");
  const [ansatz, setAnsatz] = useState("hardware_efficient");
  const [mapper, setMapper] = useState("jordan_wigner");
  const [optimizer, setOptimizer] = useState("SLSQP");
  const [backend, setBackend] = useState("classical");
  const [backendName, setBackendName] = useState("ibm_torino");

  const [optimization, setOptimization] = useState({
    geometry: false,
    orbitals: false,
    circuit: true,
    adaptive: false,
  });

  const handleSave = () => {
    // TODO: Save to localStorage or backend
    const settings = {
      method,
      ansatz,
      mapper,
      optimizer,
      backend,
      backendName,
      optimization,
    };
    console.log("Saved settings:", settings);
    onClose();
  };

  if (!isOpen) return null;

  return (
    <div className="fixed inset-0 z-50 flex items-center justify-center bg-black bg-opacity-50">
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
              <option value="HF">Hartree-Fock (HF)</option>
              <option value="VQE">Variational Quantum Eigensolver (VQE)</option>
              <option value="MP2">Møller-Plesset 2nd Order (MP2)</option>
              <option value="FCI">Full Configuration Interaction (FCI)</option>
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
                    <option value="ucc">UCC (Higher accuracy)</option>
                    <option value="hardware_efficient">
                      Hardware-Efficient (Faster)
                    </option>
                    <option value="governance">
                      Governance (Bonding-aware)
                    </option>
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
                    <option value="jordan_wigner">
                      Jordan-Wigner (Standard)
                    </option>
                    <option value="bravyi_kitaev">
                      Bravyi-Kitaev (Reduced gates)
                    </option>
                    <option value="hybrid_orbital">
                      Hybrid Orbital (Advanced)
                    </option>
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
                    <option value="SLSQP">SLSQP</option>
                    <option value="COBYLA">COBYLA</option>
                    <option value="L-BFGS-B">L-BFGS-B</option>
                    <option value="ADAM">ADAM</option>
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
                      <option value="ibm_torino">
                        ibm_torino (133 qubits)
                      </option>
                      <option value="ibm_brisbane">
                        ibm_brisbane (127 qubits)
                      </option>
                      <option value="ibm_kyoto">ibm_kyoto (127 qubits)</option>
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
