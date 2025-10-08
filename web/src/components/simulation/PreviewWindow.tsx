"use client";

import { useState } from "react";
import { ArrowLeft, CheckCircle2 } from "lucide-react";

interface PreviewWindowProps {
  molecule: any;
  backendSettings: any;
  onBack: () => void;
  onExecute: (config: any) => void;
}

export default function PreviewWindow({
  molecule,
  backendSettings,
  onBack,
  onExecute,
}: PreviewWindowProps) {
  const [acceptedTerms, setAcceptedTerms] = useState(false);
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
      onExecute({ molecule, backendSettings, analysis });
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
            <h3 className="text-lg font-quando font-semibold mb-4">
              Backend Configuration
            </h3>
            <div className="space-y-3 text-sm">
              <div className="flex justify-between">
                <span className="text-muted-foreground">Method:</span>
                <span className="font-quando">
                  {backendSettings?.method || "VQE"}
                </span>
              </div>
              <div className="flex justify-between">
                <span className="text-muted-foreground">Ansatz:</span>
                <span className="font-quando">
                  {backendSettings?.ansatz || "Hardware-Efficient"}
                </span>
              </div>
              <div className="flex justify-between">
                <span className="text-muted-foreground">Mapper:</span>
                <span className="font-quando">
                  {backendSettings?.mapper || "Jordan-Wigner"}
                </span>
              </div>
              <div className="flex justify-between">
                <span className="text-muted-foreground">Backend:</span>
                <span className="font-quando">
                  {backendSettings?.backend || "Classical"}
                </span>
              </div>
              <div className="flex justify-between">
                <span className="text-muted-foreground">Optimizer:</span>
                <span className="font-quando">
                  {backendSettings?.optimizer || "SLSQP"}
                </span>
              </div>
            </div>
            <div className="mt-4 text-xs text-muted-foreground font-quando">
              Configure in Settings to change these options
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

        {/* Right Column: Circuit Preview */}
        <div className="flex flex-col h-full">
          <div className="bg-card border border-border rounded-lg p-6 h-full flex flex-col">
            <h3 className="text-lg font-quando font-semibold mb-4">
              Quantum Circuit Preview
            </h3>

            {/* Circuit Diagram */}
            <div className="flex-1 bg-muted rounded-lg p-6 overflow-auto">
              <div className="font-mono text-sm space-y-2">
                <div>q0: ──H──●──────●──Ry(θ₁)──●──────</div>
                <div>q1: ──H──X──●───┼──────────X──Ry(θ₂)</div>
                <div>q2: ──H─────X───●──Ry(θ₃)──────────</div>
                <div>q3: ──H──────────┼──────────●──Ry(θ₄)</div>
                <div className="mt-4 text-muted-foreground">
                  Circuit depth: 8 | Parameters: 12 | Gates: 24
                </div>
              </div>
            </div>

            {/* Estimates */}
            <div className="mt-6 space-y-2 text-sm">
              <div className="flex justify-between">
                <span className="text-muted-foreground">Estimated qubits:</span>
                <span className="font-quando">12</span>
              </div>
              <div className="flex justify-between">
                <span className="text-muted-foreground">
                  Estimated runtime:
                </span>
                <span className="font-quando">~30 seconds</span>
              </div>
              <div className="flex justify-between">
                <span className="text-muted-foreground">Queue position:</span>
                <span className="font-quando">-</span>
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
              className="flex-1 px-6 py-3 border border-border rounded-lg hover:bg-accent transition font-quando"
            >
              Back
            </button>
            <button
              onClick={handleExecute}
              disabled={!acceptedTerms}
              className="flex-1 px-6 py-3 bg-brand-orange text-white rounded-lg hover:bg-brand-orange-dark transition font-quando disabled:opacity-50 disabled:cursor-not-allowed flex items-center justify-center gap-2"
            >
              <CheckCircle2 className="w-5 h-5" />
              Execute Experiment
            </button>
          </div>
        </div>
      </div>
    </div>
  );
}
