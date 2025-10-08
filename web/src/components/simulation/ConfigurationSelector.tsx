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
          {options?.methods?.map((method: string) => (
            <option key={method} value={method}>
              {method}
            </option>
          ))}
        </select>
        <p className="text-xs text-muted-foreground mt-1">
          Choose the computational method for ground state or excited states
        </p>
      </div>

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
            {options?.ansatze?.map((ansatz: string) => (
              <option key={ansatz} value={ansatz}>
                {ansatz.toUpperCase().replace("_", "-")}
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
            {options?.mappers?.map((mapper: string) => (
              <option key={mapper} value={mapper}>
                {mapper
                  .split("_")
                  .map((w) => w.charAt(0).toUpperCase() + w.slice(1))
                  .join("-")}
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
          {options?.hamiltonians?.map((ham: string) => (
            <option key={ham} value={ham}>
              {ham.charAt(0).toUpperCase() + ham.slice(1)}
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
            {options?.optimizers?.map((opt: string) => (
              <option key={opt} value={opt}>
                {opt}
              </option>
            ))}
          </select>
          <p className="text-xs text-muted-foreground mt-1">
            Classical optimizer for parameter optimization
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
          {options?.backends?.map((backend: string) => (
            <option key={backend} value={backend}>
              {backend
                .split("_")
                .map((w) => w.charAt(0).toUpperCase() + w.slice(1))
                .join(" ")}
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
