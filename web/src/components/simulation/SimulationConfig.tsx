"use client";

import { useState, useEffect } from "react";
import { ArrowLeft, Zap, Cpu, Sparkles, Database } from "lucide-react";
import * as api from "@/lib/api";

type ComputationMethod = "HF" | "VQE" | "MP2" | "FCI";
type BackendType = "classical" | "ibm_quantum" | "bluequbit";

interface SimulationConfigProps {
  onBack: () => void;
  onComplete: (config: SimulationConfiguration) => void;
}

interface SimulationConfiguration {
  method: ComputationMethod;
  ansatz?: string;
  mapper?: string;
  optimizer?: string;
  maxIterations?: number;
  backend: BackendType;
  backendName?: string;
  analysis: {
    energyDecomposition: boolean;
    bondAnalysis: boolean;
    dipoleMoment: boolean;
    polarizability: boolean;
    thermochemistry: boolean;
    spectroscopy: boolean;
    vibrational: boolean;
  };
  optimization: {
    geometry: boolean;
    orbitals: boolean;
    circuit: boolean;
    adaptive: boolean;
  };
}

export default function SimulationConfig({
  onBack,
  onComplete,
}: SimulationConfigProps) {
  const [method, setMethod] = useState<ComputationMethod>("VQE");
  const [ansatz, setAnsatz] = useState("hardware_efficient");
  const [mapper, setMapper] = useState("jordan_wigner");
  const [optimizer, setOptimizer] = useState("SLSQP");
  const [maxIterations, setMaxIterations] = useState(1000);
  const [backend, setBackend] = useState<BackendType>("classical");
  const [backendName, setBackendName] = useState("ibm_torino");

  // Configuration options from API
  const [configOptions, setConfigOptions] = useState<any>(null);

  const [analysis, setAnalysis] = useState({
    energyDecomposition: true,
    bondAnalysis: true,
    dipoleMoment: true,
    polarizability: false,
    thermochemistry: true,
    spectroscopy: false,
    vibrational: false,
  });

  const [optimization, setOptimization] = useState({
    geometry: false,
    orbitals: false,
    circuit: true,
    adaptive: false,
  });

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
    loadConfig();
  }, []);

  const handleSubmit = () => {
    const config: SimulationConfiguration = {
      method,
      ansatz: method === "VQE" ? ansatz : undefined,
      mapper: method === "VQE" ? mapper : undefined,
      optimizer: method === "VQE" ? optimizer : undefined,
      maxIterations: method === "VQE" ? maxIterations : undefined,
      backend,
      backendName: backend !== "classical" ? backendName : undefined,
      analysis,
      optimization,
    };
    onComplete(config);
  };

  return (
    <div className="h-full p-8 overflow-auto bg-background">
      <div className="max-w-4xl mx-auto">
        {/* Header */}
        <div className="flex items-center mb-8">
          <button
            onClick={onBack}
            className="mr-4 p-2 hover:bg-accent rounded-md transition"
          >
            <ArrowLeft className="w-5 h-5" />
          </button>
          <h1 className="text-3xl font-quando font-bold">
            Configure Computation
          </h1>
        </div>

        {/* Method Selection */}
        <div className="mb-8">
          <label className="block text-sm font-quando font-medium mb-3">
            Computation Method
          </label>
          <div className="grid grid-cols-4 gap-4">
            <button
              onClick={() => setMethod("HF")}
              className={`p-4 border-2 rounded-lg transition font-quando ${
                method === "HF"
                  ? "border-brand-orange bg-orange-50 dark:bg-orange-950"
                  : "border-border hover:border-muted-foreground"
              }`}
            >
              <div className="mb-2 flex justify-center">
                <Zap className="w-8 h-8" />
              </div>
              <div className="text-sm font-medium">Hartree-Fock</div>
              <div className="text-xs text-muted-foreground mt-1">Fast</div>
            </button>

            <button
              onClick={() => setMethod("VQE")}
              className={`p-4 border-2 rounded-lg transition font-quando ${
                method === "VQE"
                  ? "border-brand-orange bg-orange-50 dark:bg-orange-950"
                  : "border-border hover:border-muted-foreground"
              }`}
            >
              <div className="mb-2 flex justify-center">
                <Cpu className="w-8 h-8" />
              </div>
              <div className="text-sm font-medium">VQE</div>
              <div className="text-xs text-muted-foreground mt-1">Quantum</div>
            </button>

            <button
              onClick={() => setMethod("MP2")}
              className={`p-4 border-2 rounded-lg transition font-quando ${
                method === "MP2"
                  ? "border-brand-orange bg-orange-50 dark:bg-orange-950"
                  : "border-border hover:border-muted-foreground"
              }`}
            >
              <div className="mb-2 flex justify-center">
                <Sparkles className="w-8 h-8" />
              </div>
              <div className="text-sm font-medium">MP2</div>
              <div className="text-xs text-muted-foreground mt-1">
                Accurate
              </div>
            </button>

            <button
              onClick={() => setMethod("FCI")}
              className={`p-4 border-2 rounded-lg transition font-quando ${
                method === "FCI"
                  ? "border-brand-orange bg-orange-50 dark:bg-orange-950"
                  : "border-border hover:border-muted-foreground"
              }`}
            >
              <div className="mb-2 flex justify-center">
                <Database className="w-8 h-8" />
              </div>
              <div className="text-sm font-medium">FCI</div>
              <div className="text-xs text-muted-foreground mt-1">Exact</div>
            </button>
          </div>
        </div>

        {/* VQE Configuration */}
        {method === "VQE" && (
          <div className="space-y-6 bg-card border border-border rounded-lg p-6 mb-8">
            <h3 className="text-lg font-quando font-semibold">
              VQE Configuration
            </h3>

            <div className="grid grid-cols-2 gap-6">
              <div>
                <label className="block text-sm font-quando font-medium mb-2">
                  Ansatz
                </label>
                <select
                  value={ansatz}
                  onChange={(e) => setAnsatz(e.target.value)}
                  className="w-full px-4 py-3 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                >
                  {configOptions?.ansatze?.map((a: any) => (
                    <option key={a.value} value={a.value}>
                      {a.label} - {a.description}
                    </option>
                  )) || (
                    <>
                      <option value="ucc">UCC (Higher accuracy)</option>
                      <option value="hardware_efficient">Hardware-Efficient (Faster)</option>
                      <option value="governance">Governance (Bonding-aware)</option>
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
                  className="w-full px-4 py-3 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                >
                  {configOptions?.mappers?.map((m: any) => (
                    <option key={m.value} value={m.value}>
                      {m.label} - {m.description}
                    </option>
                  )) || (
                    <>
                      <option value="jordan_wigner">Jordan-Wigner (Standard)</option>
                      <option value="bravyi_kitaev">Bravyi-Kitaev (Reduced gates)</option>
                      <option value="hybrid_orbital">Hybrid Orbital (Advanced)</option>
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
                  className="w-full px-4 py-3 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
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
                  Max Iterations
                </label>
                <input
                  type="number"
                  value={maxIterations}
                  onChange={(e) => setMaxIterations(parseInt(e.target.value))}
                  className="w-full px-4 py-3 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                />
              </div>
            </div>
          </div>
        )}

        {/* Backend Selection */}
        <div className="mb-8">
          <label className="block text-sm font-quando font-medium mb-3">
            Backend Selection
          </label>
          <div className="space-y-3">
            <label className="flex items-center p-4 border border-border rounded-lg hover:bg-accent cursor-pointer transition">
              <input
                type="radio"
                name="backend"
                value="classical"
                checked={backend === "classical"}
                onChange={(e) => setBackend(e.target.value as BackendType)}
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
                onChange={(e) => setBackend(e.target.value as BackendType)}
                className="mr-3"
              />
              <div className="flex-1">
                <div className="font-quando font-medium">IBM Quantum (Cloud)</div>
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
                onChange={(e) => setBackend(e.target.value as BackendType)}
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

        {/* Analysis Options */}
        <div className="mb-8 bg-card border border-border rounded-lg p-6">
          <h3 className="text-lg font-quando font-semibold mb-4">
            Analysis & Properties
          </h3>
          <div className="grid grid-cols-2 gap-4">
            {Object.entries(analysis).map(([key, value]) => (
              <label key={key} className="flex items-center cursor-pointer">
                <input
                  type="checkbox"
                  checked={value}
                  onChange={(e) =>
                    setAnalysis({ ...analysis, [key]: e.target.checked })
                  }
                  className="mr-3"
                />
                <span className="font-quando text-sm">
                  {key
                    .replace(/([A-Z])/g, " $1")
                    .replace(/^./, (str) => str.toUpperCase())}
                </span>
              </label>
            ))}
          </div>
        </div>

        {/* Optimization Settings */}
        <div className="mb-8 bg-card border border-border rounded-lg p-6">
          <h3 className="text-lg font-quando font-semibold mb-4">
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

        {/* Action Buttons */}
        <div className="flex gap-4">
          <button
            onClick={onBack}
            className="flex-1 px-6 py-3 border border-border rounded-lg hover:bg-accent transition font-quando"
          >
            Back
          </button>
          <button
            onClick={handleSubmit}
            className="flex-1 px-6 py-3 bg-brand-orange text-white rounded-lg hover:bg-brand-orange-dark transition font-quando"
          >
            Review Configuration →
          </button>
        </div>
      </div>
    </div>
  );
}
