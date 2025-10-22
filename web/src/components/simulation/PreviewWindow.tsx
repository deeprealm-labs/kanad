"use client";

import { useState, useEffect, useMemo } from "react";
import { ArrowLeft, CheckCircle2, Settings } from "lucide-react";
import type { BackendSettings } from "@/lib/types";
import * as api from "@/lib/api";

interface PreviewWindowProps {
  molecule: any;
  backendSettings: BackendSettings;
  onBack: () => void;
  onExecute: (config: any) => void;
  onQueue?: (config: any) => void;
  onRefreshSettings?: () => Promise<void>;
  onOpenSettings?: () => void;
}

export default function PreviewWindow({
  molecule,
  backendSettings,
  onBack,
  onExecute,
  onQueue,
  onRefreshSettings,
  onOpenSettings,
}: PreviewWindowProps) {
  const [acceptedTerms, setAcceptedTerms] = useState(false);
  const [isLoadingSettings, setIsLoadingSettings] = useState(true);
  const [circuitPreview, setCircuitPreview] = useState<any>(null);
  const [isLoadingCircuit, setIsLoadingCircuit] = useState(false);

  // Memoize stringified versions to prevent infinite loops
  const moleculeKey = useMemo(() => JSON.stringify(molecule), [molecule]);
  const settingsKey = useMemo(() => JSON.stringify(backendSettings), [backendSettings]);

  // Reload settings when component mounts and wait for them
  useEffect(() => {
    const loadSettings = async () => {
      if (onRefreshSettings) {
        await onRefreshSettings();
      }
      setIsLoadingSettings(false);
    };
    loadSettings();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []); // Only run on mount - onRefreshSettings reference changes don't matter

  // Load circuit preview when molecule or settings change
  useEffect(() => {
    const loadCircuit = async () => {
      if (!molecule || !backendSettings) return;

      setIsLoadingCircuit(true);
      try {
        const response = await api.getCircuitPreview(molecule, backendSettings);
        if (response.success) {
          setCircuitPreview(response.preview);
        }
      } catch (error) {
        console.error("Failed to load circuit preview:", error);
        setCircuitPreview(null);
      } finally {
        setIsLoadingCircuit(false);
      }
    };

    loadCircuit();
  }, [moleculeKey, settingsKey]);
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
  const [showSettings, setShowSettings] = useState(false);

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

      {/* Main Content - 2 Columns with proper overflow */}
      <div className="flex-1 flex gap-6 p-6 overflow-hidden">
        {/* Left Column: Summary - Scrollable */}
        <div className="flex-1 flex flex-col overflow-y-auto space-y-6 pr-2">
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
              {onOpenSettings && (
                <button
                  onClick={onOpenSettings}
                  className="flex items-center gap-2 px-3 py-1 text-sm border border-border rounded-lg hover:bg-accent transition font-quando"
                >
                  <Settings className="w-4 h-4" />
                  Edit Settings
                </button>
              )}
            </div>

            <div className="space-y-3 text-sm">
                <div className="flex justify-between">
                  <span className="text-muted-foreground">Method:</span>
                  <span className="font-quando">
                    {backendSettings?.method || "VQE"}
                  </span>
                </div>
                {/* Only show ansatz for VQE */}
                {backendSettings?.method === "VQE" && backendSettings?.ansatz && (
                  <div className="flex justify-between">
                    <span className="text-muted-foreground">Ansatz:</span>
                    <span className="font-quando">
                      {backendSettings.ansatz.toUpperCase().replace("_", "-")}
                    </span>
                  </div>
                )}
                {/* Only show mapper for VQE */}
                {backendSettings?.method === "VQE" && backendSettings?.mapper && (
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
                {/* Only show optimizer for VQE */}
                {backendSettings?.method === "VQE" && backendSettings?.optimizer && (
                  <div className="flex justify-between">
                    <span className="text-muted-foreground">Optimizer:</span>
                    <span className="font-quando">
                      {backendSettings.optimizer}
                    </span>
                  </div>
                )}
              </div>
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

        {/* Right Column: Circuit Preview - 2x width, Scrollable */}
        <div className="flex-[2] flex flex-col overflow-hidden">
          <div className="bg-card border border-border rounded-lg p-6 h-full flex flex-col overflow-hidden">
            <h3 className="text-lg font-quando font-semibold mb-4 flex-shrink-0">
              Quantum Circuit Preview
            </h3>

            {/* Circuit Content - Scrollable */}
            <div className="flex-1 overflow-y-auto">
              <div className="space-y-4">
                {isLoadingCircuit ? (
                  <div className="text-muted-foreground">Loading circuit preview...</div>
                ) : circuitPreview && circuitPreview.circuit_diagram ? (
                  <>
                    {/* Note for non-VQE methods */}
                    {circuitPreview.note && (
                      <div className="p-3 bg-blue-50 dark:bg-blue-900/20 border border-blue-200 dark:border-blue-800 rounded text-sm text-blue-700 dark:text-blue-300">
                        {circuitPreview.note}
                      </div>
                    )}

                    {/* ASCII Circuit Diagram - Horizontally scrollable */}
                    <div className="bg-background rounded border border-border overflow-x-auto max-w-full">
                      <div className="p-4 min-w-max">
                        <pre className="text-xs leading-relaxed whitespace-pre font-mono">
                          {circuitPreview.circuit_diagram}
                        </pre>
                      </div>
                    </div>

                    {/* Circuit Statistics */}
                    {circuitPreview.statistics && (
                      <div className="p-4 bg-background rounded border border-border">
                        <div className="text-sm font-semibold text-muted-foreground mb-3">
                          Circuit Properties
                        </div>
                        <div className="grid grid-cols-2 gap-3 text-sm">
                          <div>
                            <span className="text-muted-foreground">Qubits: </span>
                            <span className="font-quando font-semibold">{circuitPreview.statistics.n_qubits}</span>
                          </div>
                          <div>
                            <span className="text-muted-foreground">Depth: </span>
                            <span className="font-quando font-semibold">{circuitPreview.statistics.depth}</span>
                          </div>
                          <div>
                            <span className="text-muted-foreground">Parameters: </span>
                            <span className="font-quando font-semibold">{circuitPreview.n_parameters || 0}</span>
                          </div>
                          <div>
                            <span className="text-muted-foreground">Total Gates: </span>
                            <span className="font-quando font-semibold">{circuitPreview.statistics.total_gates || 0}</span>
                          </div>
                        </div>
                        {circuitPreview.statistics.gate_counts && Object.keys(circuitPreview.statistics.gate_counts).length > 0 && (
                          <div className="mt-3 pt-3 border-t border-border">
                            <div className="text-sm text-muted-foreground mb-2">Gate Distribution:</div>
                            <div className="flex flex-wrap gap-2">
                              {Object.entries(circuitPreview.statistics.gate_counts).map(([gate, count]) => (
                                <span key={gate} className="px-2 py-1 bg-muted rounded text-sm font-mono">
                                  {gate.toUpperCase()}: <span className="font-semibold">{String(count)}</span>
                                </span>
                              ))}
                            </div>
                          </div>
                        )}
                      </div>
                    )}

                    {/* Configuration Info - only for VQE */}
                    {backendSettings.method === "VQE" && (
                      <div className="p-4 bg-muted rounded border border-border">
                        <div className="text-sm font-semibold mb-3">VQE Configuration</div>
                        <div className="space-y-2 text-sm">
                          <div>
                            <span className="text-muted-foreground">Ansatz: </span>
                            <span className="font-quando font-semibold">{backendSettings.ansatz?.toUpperCase().replace("_", "-") || "UCC"}</span>
                          </div>
                          <div>
                            <span className="text-muted-foreground">Mapper: </span>
                            <span className="font-quando font-semibold">
                              {backendSettings.mapper?.split("_").map(w => w.charAt(0).toUpperCase() + w.slice(1)).join("-") || "Jordan-Wigner"}
                            </span>
                          </div>
                        </div>
                      </div>
                    )}
                  </>
                ) : (
                  <div className="text-muted-foreground p-8 text-center">
                    Unable to generate circuit preview
                  </div>
                )}
              </div>
            </div>

            {/* Estimates - Fixed at bottom */}
            <div className="mt-4 p-4 bg-muted rounded border border-border flex-shrink-0">
              <div className="grid grid-cols-3 gap-4 text-sm">
                <div>
                  <div className="text-muted-foreground text-xs mb-1">Method</div>
                  <div className="font-quando font-semibold">{backendSettings.method || "VQE"}</div>
                </div>
                <div>
                  <div className="text-muted-foreground text-xs mb-1">Backend</div>
                  <div className="font-quando font-semibold">
                    {backendSettings.backend?.split("_").map(w => w.charAt(0).toUpperCase() + w.slice(1)).join(" ") || "Classical"}
                  </div>
                </div>
                <div>
                  <div className="text-muted-foreground text-xs mb-1">Est. Runtime</div>
                  <div className="font-quando font-semibold">
                    {backendSettings.method === "HF" ? "< 1s" : "~30-60s"}
                  </div>
                </div>
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
