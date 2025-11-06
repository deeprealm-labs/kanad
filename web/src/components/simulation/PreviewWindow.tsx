"use client";

import { useState, useEffect, useMemo } from "react";
import { ArrowLeft, CheckCircle2, Settings, Eye, X } from "lucide-react";
import type { BackendSettings } from "@/lib/types";
import * as api from "@/lib/api";
import QuantumCircuitViewer from "@/components/circuit/QuantumCircuitViewer";

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
  const [showCircuitModal, setShowCircuitModal] = useState(false);

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

  return (
    <div className="h-full flex flex-col bg-background">
      {/* Compact Header */}
      <div className="flex items-center gap-3 px-6 py-4 border-b border-border bg-card">
        <button
          onClick={onBack}
          className="p-1.5 hover:bg-accent rounded-md transition"
        >
          <ArrowLeft className="w-4 h-4" />
        </button>
        <div>
          <h1 className="text-xl font-quando font-bold">Experiment Preview</h1>
          <p className="text-xs text-muted-foreground">Review your configuration before execution</p>
        </div>
      </div>

      {/* Main Content */}
      <div className="flex-1 overflow-y-auto p-8">
        <div className="max-w-5xl mx-auto space-y-6">
          {/* Configuration Summary - Single Panel */}
          <div className="bg-card border border-border rounded-lg p-6">
            <h2 className="text-lg font-quando font-semibold mb-6">Configuration Summary</h2>

            <div className="grid grid-cols-2 gap-x-12 gap-y-6">
              {/* Left Column - Molecule */}
              <div className="space-y-4">
                <h3 className="text-sm font-quando font-semibold text-brand-orange uppercase tracking-wide">Molecule</h3>
                <div className="space-y-3 text-sm">
                  <div className="flex justify-between">
                    <span className="text-muted-foreground">Formula:</span>
                    <span className="font-mono font-semibold">
                      {molecule.smiles || (molecule.atoms ? `${molecule.atoms.length} atoms` : "N/A")}
                    </span>
                  </div>
                  <div className="flex justify-between">
                    <span className="text-muted-foreground">Basis Set:</span>
                    <span className="font-quando">{molecule.basis}</span>
                  </div>
                  <div className="flex justify-between">
                    <span className="text-muted-foreground">Charge:</span>
                    <span className="font-mono">{molecule.charge}</span>
                  </div>
                  <div className="flex justify-between">
                    <span className="text-muted-foreground">Multiplicity:</span>
                    <span className="font-mono">{molecule.multiplicity}</span>
                  </div>
                </div>
              </div>

              {/* Right Column - Quantum Method */}
              <div className="space-y-4">
                <h3 className="text-sm font-quando font-semibold text-brand-orange uppercase tracking-wide">Quantum Method</h3>
                <div className="space-y-3 text-sm">
                  <div className="flex justify-between">
                    <span className="text-muted-foreground">Solver:</span>
                    <span className="font-quando font-semibold">
                      {backendSettings?.method || "VQE"}
                    </span>
                  </div>
                  {backendSettings?.method === "VQE" && backendSettings?.ansatz && (
                    <div className="flex justify-between">
                      <span className="text-muted-foreground">Ansatz:</span>
                      <span className="font-quando">
                        {backendSettings.ansatz.toUpperCase().replace("_", "-")}
                      </span>
                    </div>
                  )}
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
                  <div className="flex justify-between">
                    <span className="text-muted-foreground">Backend:</span>
                    <span className="font-quando">
                      {backendSettings?.backend
                        ?.split("_")
                        .map((w) => w.charAt(0).toUpperCase() + w.slice(1))
                        .join(" ") || "Classical"}
                    </span>
                  </div>
                  {backendSettings?.method === "VQE" && backendSettings?.optimizer && (
                    <div className="flex justify-between">
                      <span className="text-muted-foreground">Optimizer:</span>
                      <span className="font-quando">{backendSettings.optimizer}</span>
                    </div>
                  )}
                </div>
              </div>
            </div>

            {/* Circuit Preview Button */}
            <div className="mt-6 pt-6 border-t border-border">
              <button
                onClick={() => setShowCircuitModal(true)}
                disabled={isLoadingCircuit}
                className="w-full flex items-center justify-center gap-3 px-6 py-4 bg-brand-orange text-white rounded-lg hover:bg-brand-orange-dark transition font-quando font-semibold disabled:opacity-50 disabled:cursor-not-allowed"
              >
                <Eye className="w-5 h-5" />
                {isLoadingCircuit ? "Loading Circuit..." : "View Quantum Circuit Preview"}
              </button>
              {circuitPreview && circuitPreview.statistics && (
                <div className="mt-3 flex justify-center gap-6 text-xs text-muted-foreground">
                  <span>Qubits: <strong className="text-foreground">{circuitPreview.n_qubits || circuitPreview.statistics.n_qubits}</strong></span>
                  <span>Depth: <strong className="text-foreground">{circuitPreview.statistics.depth}</strong></span>
                  <span>Gates: <strong className="text-foreground">{circuitPreview.statistics.total_gates}</strong></span>
                  {circuitPreview.n_parameters > 0 && (
                    <span>Parameters: <strong className="text-foreground">{circuitPreview.n_parameters}</strong></span>
                  )}
                </div>
              )}
            </div>
          </div>
        </div>
      </div>

      {/* Bottom Bar: Terms & Execute */}
      <div className="border-t border-border p-6 bg-card">
        <div className="max-w-5xl mx-auto space-y-4">
          {/* Terms and Conditions */}
          <div className="bg-background border border-border rounded-lg p-4">
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
                  and used for computation only.
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

      {/* Circuit Preview Modal */}
      {showCircuitModal && (
        <div className="fixed inset-0 z-50 flex items-center justify-center bg-black/60 backdrop-blur-sm p-4">
          <div className="bg-card border border-border rounded-lg w-full max-w-6xl max-h-[90vh] flex flex-col shadow-2xl">
            {/* Modal Header */}
            <div className="flex items-center justify-between p-6 border-b border-border">
              <div>
                <h2 className="text-xl font-quando font-bold">Quantum Circuit</h2>
                <p className="text-sm text-muted-foreground mt-1">
                  Generated circuit for {backendSettings?.method || "VQE"} solver
                </p>
              </div>
              <button
                onClick={() => setShowCircuitModal(false)}
                className="p-2 hover:bg-accent rounded-lg transition"
              >
                <X className="w-5 h-5" />
              </button>
            </div>

            {/* Modal Content */}
            <div className="flex-1 overflow-y-auto p-6">
              {isLoadingCircuit ? (
                <div className="flex items-center justify-center h-64">
                  <div className="text-muted-foreground">Loading circuit preview...</div>
                </div>
              ) : circuitPreview && circuitPreview.circuit_diagram ? (
                <div className="space-y-4">
                  {circuitPreview.note && (
                    <div className="p-3 bg-blue-50 dark:bg-blue-900/20 border border-blue-200 dark:border-blue-800 rounded text-sm text-blue-700 dark:text-blue-300">
                      {circuitPreview.note}
                    </div>
                  )}

                  {console.log('PreviewWindow: circuitPreview.circuit_image =', circuitPreview.circuit_image ? `${circuitPreview.circuit_image.substring(0, 50)}... (${circuitPreview.circuit_image.length} chars)` : 'undefined/null')}

                  <QuantumCircuitViewer
                    circuitDiagram={circuitPreview.circuit_diagram}
                    stats={circuitPreview.statistics}
                    gates={circuitPreview.gates}
                    circuitImage={circuitPreview.circuit_image}
                    className="w-full"
                  />
                </div>
              ) : (
                <div className="flex items-center justify-center h-64">
                  <div className="text-muted-foreground">Unable to generate circuit preview</div>
                </div>
              )}
            </div>

            {/* Modal Footer */}
            <div className="flex items-center justify-between p-6 border-t border-border bg-background">
              {circuitPreview && circuitPreview.statistics && (
                <div className="flex gap-6 text-sm">
                  <div>
                    <span className="text-muted-foreground">Qubits:</span>
                    <span className="ml-2 font-mono font-semibold">{circuitPreview.n_qubits || circuitPreview.statistics.n_qubits}</span>
                  </div>
                  <div>
                    <span className="text-muted-foreground">Depth:</span>
                    <span className="ml-2 font-mono font-semibold">{circuitPreview.statistics.depth}</span>
                  </div>
                  <div>
                    <span className="text-muted-foreground">Gates:</span>
                    <span className="ml-2 font-mono font-semibold">{circuitPreview.statistics.total_gates}</span>
                  </div>
                  {circuitPreview.n_parameters > 0 && (
                    <div>
                      <span className="text-muted-foreground">Parameters:</span>
                      <span className="ml-2 font-mono font-semibold">{circuitPreview.n_parameters}</span>
                    </div>
                  )}
                </div>
              )}
              <button
                onClick={() => setShowCircuitModal(false)}
                className="px-6 py-2 bg-brand-orange text-white rounded-lg hover:bg-brand-orange-dark transition font-quando"
              >
                Close
              </button>
            </div>
          </div>
        </div>
      )}
    </div>
  );
}
