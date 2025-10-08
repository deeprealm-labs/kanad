"use client";

import { useState } from "react";
import { ArrowLeft, CheckCircle2, Settings } from "lucide-react";
import ConfigurationSelector from "./ConfigurationSelector";
import type { BackendSettings } from "@/lib/types";

interface PreviewWindowProps {
  molecule: any;
  backendSettings: BackendSettings;
  onBack: () => void;
  onExecute: (config: any) => void;
  onQueue?: (config: any) => void;
}

export default function PreviewWindow({
  molecule,
  backendSettings: initialSettings,
  onBack,
  onExecute,
  onQueue,
}: PreviewWindowProps) {
  const [acceptedTerms, setAcceptedTerms] = useState(false);
  const [showAdvancedSettings, setShowAdvancedSettings] = useState(false);
  const [backendSettings, setBackendSettings] = useState<BackendSettings>(initialSettings);
  const [analysis, setAnalysis] = useState({
    energyDecomposition: true,
    bondAnalysis: true,
    dipoleMoment: true,
    polarizability: false,
    thermochemistry: true,
    spectroscopy: false,
  });

  const handleExecute = () => {
    if (acceptedTerms) {
      const config = { molecule, backendSettings, analysis };
      console.log("PreviewWindow.handleExecute - config:", config);
      console.log("PreviewWindow.handleExecute - backendSettings:", backendSettings);
      onExecute(config);
    }
  };

  const handleQueue = () => {
    if (acceptedTerms && onQueue) {
      onQueue({ molecule, backendSettings, analysis });
    }
  };

  return (
    <div className="h-full flex flex-col bg-background">
      {/* Header */}
      <div className="flex items-center justify-between p-6 border-b border-border">
        <div className="flex items-center gap-4">
          <button
            onClick={onBack}
            className="p-2 hover:bg-accent rounded-md transition"
          >
            <ArrowLeft className="w-5 h-5" />
          </button>
          <h1 className="text-2xl font-quando font-bold">
            Review Configuration
          </h1>
        </div>
      </div>

      {/* Main Content - 2 Columns */}
      <div className="flex-1 grid grid-cols-2 gap-6 p-6 overflow-hidden">
        {/* Left Column: Summary */}
        <div className="flex flex-col h-full overflow-auto space-y-6">
          {/* Molecular Configuration */}
          <div className="bg-card border border-border rounded-lg p-6">
            <h3 className="text-lg font-quando font-semibold mb-4">
              Molecular Configuration
            </h3>
            <div className="space-y-3 text-sm">
              <div className="flex justify-between">
                <span className="text-muted-foreground">SMILES:</span>
                <span className="font-mono">
                  {molecule.smiles || "From atoms"}
                </span>
              </div>
              <div className="flex justify-between">
                <span className="text-muted-foreground">Basis Set:</span>
                <span className="font-quando">{molecule.basis}</span>
              </div>
              <div className="flex justify-between">
                <span className="text-muted-foreground">Charge:</span>
                <span className="font-quando">{molecule.charge}</span>
              </div>
              <div className="flex justify-between">
                <span className="text-muted-foreground">Multiplicity:</span>
                <span className="font-quando">{molecule.multiplicity}</span>
              </div>
              {molecule.atoms && molecule.atoms.length > 0 && (
                <div>
                  <div className="text-muted-foreground mb-2">Atoms:</div>
                  <div className="bg-muted rounded p-2 font-mono text-xs max-h-32 overflow-auto">
                    {molecule.atoms.map((atom: any, i: number) => (
                      <div key={i}>
                        {atom.symbol} {atom.x.toFixed(4)} {atom.y.toFixed(4)}{" "}
                        {atom.z.toFixed(4)}
                      </div>
                    ))}
                  </div>
                </div>
              )}
            </div>
          </div>

          {/* Backend Configuration */}
          <div className="bg-card border border-border rounded-lg p-6">
            <div className="flex items-center justify-between mb-4">
              <h3 className="text-lg font-quando font-semibold">
                Backend Configuration
              </h3>
              <button
                onClick={() => setShowAdvancedSettings(!showAdvancedSettings)}
                className="flex items-center gap-2 px-3 py-1 text-sm border border-border rounded-lg hover:bg-accent transition font-quando"
              >
                <Settings className="w-4 h-4" />
                {showAdvancedSettings ? "Hide" : "Show"} Settings
              </button>
            </div>

            {showAdvancedSettings ? (
              <ConfigurationSelector
                settings={backendSettings}
                onChange={setBackendSettings}
              />
            ) : (
              <div className="space-y-3 text-sm">
                <div className="flex justify-between">
                  <span className="text-muted-foreground">Method:</span>
                  <span className="font-quando">
                    {backendSettings?.method || "VQE"}
                  </span>
                </div>
                {backendSettings?.ansatz && (
                  <div className="flex justify-between">
                    <span className="text-muted-foreground">Ansatz:</span>
                    <span className="font-quando">
                      {backendSettings.ansatz.toUpperCase().replace("_", "-")}
                    </span>
                  </div>
                )}
                {backendSettings?.mapper && (
                  <div className="flex justify-between">
                    <span className="text-muted-foreground">Mapper:</span>
                    <span className="font-quando">
                      {backendSettings.mapper
                        .split("_")
                        .map((w) => w.charAt(0).toUpperCase() + w.slice(1))
                        .join("-")}
                    </span>
                  </div>
                )}
                {backendSettings?.hamiltonian && (
                  <div className="flex justify-between">
                    <span className="text-muted-foreground">Hamiltonian:</span>
                    <span className="font-quando">
                      {backendSettings.hamiltonian.charAt(0).toUpperCase() +
                        backendSettings.hamiltonian.slice(1)}
                    </span>
                  </div>
                )}
                <div className="flex justify-between">
                  <span className="text-muted-foreground">Backend:</span>
                  <span className="font-quando">
                    {backendSettings?.backend
                      ?.split("_")
                      .map((w) => w.charAt(0).toUpperCase() + w.slice(1))
                      .join(" ") || "Classical"}
                  </span>
                </div>
                {backendSettings?.optimizer && (
                  <div className="flex justify-between">
                    <span className="text-muted-foreground">Optimizer:</span>
                    <span className="font-quando">
                      {backendSettings.optimizer}
                    </span>
                  </div>
                )}
              </div>
            )}
          </div>

          {/* Analysis Properties */}
          <div className="bg-card border border-border rounded-lg p-6">
            <h3 className="text-lg font-quando font-semibold mb-4">
              Analysis Properties
            </h3>
            <div className="space-y-2">
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
        </div>

        {/* Right Column: Circuit Preview */}
        <div className="flex flex-col h-full">
          <div className="bg-card border border-border rounded-lg p-6 h-full flex flex-col">
            <h3 className="text-lg font-quando font-semibold mb-4">
              Quantum Circuit Preview
            </h3>

            {/* Circuit Diagram */}
            <div className="flex-1 bg-muted rounded-lg p-6 overflow-auto">
              <div className="font-mono text-sm space-y-2">
                {backendSettings.method === "VQE" || backendSettings.method === "SQD" ? (
                  <>
                    <div className="text-muted-foreground mb-4">
                      Circuit will be generated when experiment starts
                    </div>
                    <div className="space-y-2">
                      <div className="text-xs text-muted-foreground">Preview based on configuration:</div>
                      <div className="mt-2">
                        <span className="text-muted-foreground">Ansatz: </span>
                        <span className="font-quando">{backendSettings.ansatz?.toUpperCase() || "UCC"}</span>
                      </div>
                      <div>
                        <span className="text-muted-foreground">Mapper: </span>
                        <span className="font-quando">
                          {backendSettings.mapper?.split("_").map(w => w.charAt(0).toUpperCase() + w.slice(1)).join("-") || "Jordan-Wigner"}
                        </span>
                      </div>
                      {molecule.atoms && (
                        <div>
                          <span className="text-muted-foreground">Atoms: </span>
                          <span className="font-quando">{molecule.atoms.length}</span>
                        </div>
                      )}
                      <div className="mt-4 p-3 bg-background rounded border border-border">
                        <div className="text-xs text-muted-foreground mb-2">
                          Estimated Circuit Properties:
                        </div>
                        <div className="text-xs space-y-1">
                          <div>Qubits: ~{molecule.atoms ? molecule.atoms.length * 2 : 4}-{molecule.atoms ? molecule.atoms.length * 4 : 8}</div>
                          <div>Parameters: Depends on ansatz</div>
                          <div>Gates: Generated at runtime</div>
                        </div>
                      </div>
                    </div>
                  </>
                ) : (
                  <div className="text-muted-foreground">
                    {backendSettings.method === "HF"
                      ? "Hartree-Fock method does not use quantum circuits"
                      : "Circuit preview available after experiment creation"}
                  </div>
                )}
              </div>
            </div>

            {/* Estimates */}
            <div className="mt-6 space-y-2 text-sm">
              <div className="flex justify-between">
                <span className="text-muted-foreground">Method:</span>
                <span className="font-quando">{backendSettings.method || "VQE"}</span>
              </div>
              <div className="flex justify-between">
                <span className="text-muted-foreground">Backend:</span>
                <span className="font-quando">
                  {backendSettings.backend?.split("_").map(w => w.charAt(0).toUpperCase() + w.slice(1)).join(" ") || "Classical"}
                </span>
              </div>
              <div className="flex justify-between">
                <span className="text-muted-foreground">
                  Estimated runtime:
                </span>
                <span className="font-quando">
                  {backendSettings.method === "HF" ? "< 1 second" : "~30-60 seconds"}
                </span>
              </div>
            </div>
          </div>
        </div>
      </div>

      {/* Bottom Bar: Terms & Execute */}
      <div className="border-t border-border p-6">
        <div className="max-w-4xl mx-auto space-y-4">
          {/* Terms and Conditions */}
          <div className="bg-card border border-border rounded-lg p-4">
            <label className="flex items-start cursor-pointer">
              <input
                type="checkbox"
                checked={acceptedTerms}
                onChange={(e) => setAcceptedTerms(e.target.checked)}
                className="mr-3 mt-1"
              />
              <div className="text-sm font-quando">
                <p className="font-medium mb-1">
                  I accept the terms and conditions
                </p>
                <p className="text-muted-foreground text-xs">
                  By using IBM Quantum or cloud backends, you agree to their
                  respective terms of service. Results will be stored securely
                  and used for computation only. This experiment will consume
                  quantum credits as estimated above.
                </p>
              </div>
            </label>
          </div>

          {/* Action Buttons */}
          <div className="flex gap-4">
            <button
              onClick={onBack}
              className="px-6 py-3 border border-border rounded-lg hover:bg-accent transition font-quando"
            >
              Back
            </button>
            {onQueue && (
              <button
                onClick={handleQueue}
                disabled={!acceptedTerms}
                className="flex-1 px-6 py-3 border-2 border-brand-orange text-brand-orange rounded-lg hover:bg-brand-orange hover:text-white transition font-quando disabled:opacity-50 disabled:cursor-not-allowed flex items-center justify-center gap-2"
              >
                <CheckCircle2 className="w-5 h-5" />
                Add to Queue
              </button>
            )}
            <button
              onClick={handleExecute}
              disabled={!acceptedTerms}
              className="flex-1 px-6 py-3 bg-brand-orange text-white rounded-lg hover:bg-brand-orange-dark transition font-quando disabled:opacity-50 disabled:cursor-not-allowed flex items-center justify-center gap-2"
            >
              <CheckCircle2 className="w-5 h-5" />
              Execute Now
            </button>
          </div>
        </div>
      </div>
    </div>
  );
}
