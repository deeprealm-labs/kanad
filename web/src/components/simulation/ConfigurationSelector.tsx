"use client";

import { useState, useEffect } from "react";
import * as api from "@/lib/api";
import type { BackendSettings } from "@/lib/types";

interface ConfigurationSelectorProps {
  settings: BackendSettings;
  onChange: (settings: BackendSettings) => void;
}

export default function ConfigurationSelector({
  settings,
  onChange,
}: ConfigurationSelectorProps) {
  const [options, setOptions] = useState<any>(null);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    loadOptions();
  }, []);

  const loadOptions = async () => {
    try {
      const opts = await api.getConfigurationOptions();
      setOptions(opts);
    } catch (error) {
      console.error("Failed to load configuration options:", error);
    } finally {
      setLoading(false);
    }
  };

  const updateSetting = (key: keyof BackendSettings, value: any) => {
    onChange({ ...settings, [key]: value });
  };

  if (loading) {
    return (
      <div className="text-sm text-muted-foreground">
        Loading configuration options...
      </div>
    );
  }

  // Helper function to format option names
  const formatOptionName = (name: string) => {
    if (!name || typeof name !== 'string') return '';
    // Replace underscores with spaces and capitalize words
    return name
      .split("_")
      .map((word) => word.charAt(0).toUpperCase() + word.slice(1))
      .join(" ");
  };

  return (
    <div className="space-y-4">
      {/* Method Selection */}
      <div>
        <label className="block text-sm font-quando font-medium mb-2">
          Method
        </label>
        <select
          value={settings.method}
          onChange={(e) => updateSetting("method", e.target.value)}
          className="w-full px-4 py-2 border border-border rounded-lg bg-background focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando text-sm"
        >
          {options?.methods?.map((method: any) => (
            <option key={method.value} value={method.value}>
              {method.label}
            </option>
          ))}
        </select>
        <p className="text-xs text-muted-foreground mt-1">
          Choose the computational method for ground state or excited states
        </p>
      </div>

      {/* Excited States Method Selection */}
      {settings.method === "EXCITED_STATES" && (
        <div>
          <label className="block text-sm font-quando font-medium mb-2">
            Excited States Method
          </label>
          <select
            value={settings.excitedMethod || "cis"}
            onChange={(e) => updateSetting("excitedMethod", e.target.value)}
            className="w-full px-4 py-2 border border-border rounded-lg bg-background focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando text-sm"
          >
            <option value="cis">CIS (Configuration Interaction Singles - Recommended)</option>
            <option value="sqd">SQD (Subspace Quantum Diagonalization)</option>
            <option value="tddft">TDDFT (Time-Dependent DFT)</option>
          </select>
          <p className="text-xs text-muted-foreground mt-1">
            {settings.excitedMethod === "cis" && "Fast, accurate for small-medium molecules. Classical computation."}
            {settings.excitedMethod === "sqd" && "Quantum subspace method - excellent for excited states. Lower circuit depth."}
            {settings.excitedMethod === "tddft" && "More accurate for larger systems. Classical computation."}
          </p>
        </div>
      )}

      {/* Number of Excited States */}
      {settings.method === "EXCITED_STATES" && (
        <div>
          <label className="block text-sm font-quando font-medium mb-2">
            Number of States
          </label>
          <input
            type="number"
            min="2"
            max="10"
            value={settings.nStates || 5}
            onChange={(e) => updateSetting("nStates", parseInt(e.target.value))}
            className="w-full px-4 py-2 border border-border rounded-lg bg-background focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando text-sm"
          />
          <p className="text-xs text-muted-foreground mt-1">
            Total number of states to compute (including ground state)
          </p>
        </div>
      )}

      {/* SQD Settings for Excited States */}
      {settings.method === "EXCITED_STATES" && settings.excitedMethod === "sqd" && (
        <div className="grid grid-cols-2 gap-4">
          <div>
            <label className="block text-sm font-quando font-medium mb-2">
              Subspace Dimension
            </label>
            <input
              type="number"
              min="4"
              max="20"
              value={settings.subspaceDim || 10}
              onChange={(e) => updateSetting("subspaceDim", parseInt(e.target.value))}
              className="w-full px-4 py-2 border border-border rounded-lg bg-background focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando text-sm"
            />
            <p className="text-xs text-muted-foreground mt-1">
              Larger = more accurate but slower
            </p>
          </div>
          <div>
            <label className="block text-sm font-quando font-medium mb-2">
              Circuit Depth
            </label>
            <input
              type="number"
              min="1"
              max="5"
              value={settings.circuitDepth || 3}
              onChange={(e) => updateSetting("circuitDepth", parseInt(e.target.value))}
              className="w-full px-4 py-2 border border-border rounded-lg bg-background focus:outline-none focus:ring-2 focus:ring-brand-orange font-quanto text-sm"
            />
            <p className="text-xs text-muted-foreground mt-1">
              Quantum circuit complexity
            </p>
          </div>
        </div>
      )}


      {/* Ansatz Selection (for VQE/SQD) */}
      {(settings.method === "VQE" || settings.method === "SQD") && (
        <div>
          <label className="block text-sm font-quando font-medium mb-2">
            Ansatz
          </label>
          <select
            value={settings.ansatz}
            onChange={(e) => updateSetting("ansatz", e.target.value)}
            className="w-full px-4 py-2 border border-border rounded-lg bg-background focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando text-sm"
          >
            {options?.ansatze?.map((ansatz: any) => (
              <option key={ansatz.value} value={ansatz.value}>
                {ansatz.label}
              </option>
            ))}
          </select>
          <p className="text-xs text-muted-foreground mt-1">
            Quantum circuit ansatz for variational algorithms
          </p>
        </div>
      )}

      {/* Mapper Selection */}
      {(settings.method === "VQE" || settings.method === "SQD") && (
        <div>
          <label className="block text-sm font-quando font-medium mb-2">
            Qubit Mapper
          </label>
          <select
            value={settings.mapper}
            onChange={(e) => updateSetting("mapper", e.target.value)}
            className="w-full px-4 py-2 border border-border rounded-lg bg-background focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando text-sm"
          >
            {options?.mappers?.map((mapper: any) => (
              <option key={mapper.value} value={mapper.value}>
                {mapper.label}
              </option>
            ))}
          </select>
          <p className="text-xs text-muted-foreground mt-1">
            Fermion-to-qubit mapping transformation
          </p>
        </div>
      )}

      {/* Hamiltonian Type */}
      <div>
        <label className="block text-sm font-quando font-medium mb-2">
          Hamiltonian Type
        </label>
        <select
          value={settings.hamiltonian || "covalent"}
          onChange={(e) => updateSetting("hamiltonian", e.target.value)}
          className="w-full px-4 py-2 border border-border rounded-lg bg-background focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando text-sm"
        >
          {options?.hamiltonians?.map((ham: any) => (
            <option key={ham.value} value={ham.value}>
              {ham.label}
            </option>
          ))}
        </select>
        <p className="text-xs text-muted-foreground mt-1">
          Bond type governance for molecular Hamiltonian
        </p>
      </div>

      {/* Optimizer Selection */}
      {(settings.method === "VQE" || settings.method === "SQD") && (
        <div>
          <label className="block text-sm font-quando font-medium mb-2">
            Optimizer
          </label>
          <select
            value={settings.optimizer}
            onChange={(e) => updateSetting("optimizer", e.target.value)}
            className="w-full px-4 py-2 border border-border rounded-lg bg-background focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando text-sm"
          >
            {options?.optimizers?.map((opt: any) => (
              <option key={opt.value} value={opt.value}>
                {opt.label}
              </option>
            ))}
          </select>
          <p className="text-xs text-muted-foreground mt-1">
            {settings.optimizer === "COBYLA" && "✅ COBYLA: ~3 function evals/iter (BEST for quantum backends)"}
            {settings.optimizer === "Powell" && "✅ Powell: ~5 function evals/iter (good for quantum backends)"}
            {settings.optimizer === "Nelder-Mead" && "⚠️ Nelder-Mead: ~10 function evals/iter (moderate cost)"}
            {settings.optimizer === "SLSQP" && "⚠️ SLSQP: gradient-based, many evals/iter (HIGH cost on quantum!)"}
            {settings.optimizer === "L-BFGS-B" && "⚠️ L-BFGS-B: gradient-based, many evals/iter (HIGH cost on quantum!)"}
            {settings.optimizer === "BFGS" && "⚠️ BFGS: gradient-based, many evals/iter (HIGH cost on quantum!)"}
            {settings.optimizer === "CG" && "⚠️ CG: gradient-based, many evals/iter (HIGH cost on quantum!)"}
            {settings.optimizer === "TNC" && "⚠️ TNC: gradient-based, many evals/iter (HIGH cost on quantum!)"}
          </p>
        </div>
      )}

      {/* Max Iterations */}
      {(settings.method === "VQE" || settings.method === "SQD") && (
        <div>
          <label className="block text-sm font-quando font-medium mb-2">
            Max Iterations
          </label>
          <input
            type="number"
            min="10"
            max="500"
            value={settings.maxIterations || 100}
            onChange={(e) => updateSetting("maxIterations", parseInt(e.target.value))}
            className="w-full px-4 py-2 border border-border rounded-lg bg-background focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando text-sm"
          />
          <p className="text-xs text-muted-foreground mt-1">
            {settings.optimizer === "COBYLA" && `✅ ~${(settings.maxIterations || 100) * 3} total function evaluations`}
            {settings.optimizer === "Powell" && `✅ ~${(settings.maxIterations || 100) * 5} total function evaluations`}
            {settings.optimizer === "Nelder-Mead" && `⚠️ ~${(settings.maxIterations || 100) * 10} total function evaluations`}
            {settings.optimizer === "SLSQP" && `⚠️ Gradient-based: HIGH cost on quantum backends!`}
            {settings.optimizer === "L-BFGS-B" && `⚠️ Gradient-based: HIGH cost on quantum backends!`}
            {settings.optimizer === "BFGS" && `⚠️ Gradient-based: HIGH cost on quantum backends!`}
            {settings.optimizer === "CG" && `⚠️ Gradient-based: HIGH cost on quantum backends!`}
            {settings.optimizer === "TNC" && `⚠️ Gradient-based: HIGH cost on quantum backends!`}
          </p>
        </div>
      )}

      {/* Backend Selection */}
      <div>
        <label className="block text-sm font-quando font-medium mb-2">
          Backend
        </label>
        <select
          value={settings.backend}
          onChange={(e) => updateSetting("backend", e.target.value)}
          className="w-full px-4 py-2 border border-border rounded-lg bg-background focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando text-sm"
        >
          {options?.backends?.map((backend: any) => (
            <option key={backend.value} value={backend.value}>
              {backend.label}
            </option>
          ))}
        </select>
        <p className="text-xs text-muted-foreground mt-1">
          Execution backend (classical simulator or quantum hardware)
        </p>
      </div>

      {/* Backend Name (for IBM/BlueQubit) */}
      {(settings.backend === "ibm_quantum" ||
        settings.backend === "bluequbit") && (
        <div>
          <label className="block text-sm font-quando font-medium mb-2">
            {settings.backend === "ibm_quantum"
              ? "IBM Backend Name"
              : "BlueQubit Instance"}
          </label>
          <input
            type="text"
            value={settings.backendName || ""}
            onChange={(e) => updateSetting("backendName", e.target.value)}
            placeholder={
              settings.backend === "ibm_quantum"
                ? "e.g., ibm_torino, ibm_brisbane"
                : "e.g., bluequbit-gpu"
            }
            className="w-full px-4 py-2 border border-border rounded-lg bg-background focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando text-sm"
          />
          <p className="text-xs text-muted-foreground mt-1">
            Specific backend instance name
          </p>
        </div>
      )}
    </div>
  );
}
