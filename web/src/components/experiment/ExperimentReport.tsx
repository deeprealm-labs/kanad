"use client";

import { useState, useEffect } from "react";
import { X, Download, FileText, Loader2 } from "lucide-react";
import {
  LineChart,
  Line,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  ResponsiveContainer,
} from "recharts";
import * as api from "@/lib/api";
import type { ExperimentReport as ReportType } from "@/lib/types";
import { useToast } from "@/components/ui/toast";

interface ExperimentReportProps {
  experimentId: string;
  onClose: () => void;
}

export default function ExperimentReport({
  experimentId,
  onClose,
}: ExperimentReportProps) {
  const [report, setReport] = useState<ReportType | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const toast = useToast();

  useEffect(() => {
    loadReport();
  }, [experimentId]);

  const loadReport = async () => {
    try {
      setLoading(true);
      const data = await api.getExperimentReport(experimentId);
      setReport(data);
      setError(null);
    } catch (err: any) {
      console.error("Failed to load report:", err);
      setError(err.message || "Failed to load experiment report");
    } finally {
      setLoading(false);
    }
  };

  const handleDownloadJSON = async () => {
    try {
      await api.downloadExperiment(experimentId, "json");
      toast.success("Report downloaded successfully!");
    } catch (err: any) {
      toast.error(err.message || "Failed to download report");
    }
  };

  const handleDownloadPDF = () => {
    toast.info("PDF export coming soon!");
  };

  if (loading) {
    return (
      <div className="fixed inset-0 z-50 flex items-center justify-center backdrop-blur-sm bg-black/50">
        <div className="bg-background border border-border rounded-lg p-12">
          <Loader2 className="w-12 h-12 animate-spin text-brand-orange mx-auto mb-4" />
          <p className="text-muted-foreground font-quando">Loading report...</p>
        </div>
      </div>
    );
  }

  if (error || !report) {
    return (
      <div className="fixed inset-0 z-50 flex items-center justify-center backdrop-blur-sm bg-black/50">
        <div className="bg-background border border-border rounded-lg p-12 max-w-md">
          <div className="text-center">
            <FileText className="w-12 h-12 text-red-500 mx-auto mb-4" />
            <h3 className="text-lg font-quando font-bold mb-2">
              Failed to Load Report
            </h3>
            <p className="text-sm text-muted-foreground mb-6">{error}</p>
            <button
              onClick={onClose}
              className="px-6 py-2 bg-brand-orange text-white rounded-lg hover:bg-brand-orange-dark transition font-quando"
            >
              Close
            </button>
          </div>
        </div>
      </div>
    );
  }

  return (
    <div className="fixed inset-0 z-50 flex items-center justify-center backdrop-blur-sm bg-black/50 p-6">
      <div className="bg-background border border-border rounded-lg w-full max-w-6xl max-h-[90vh] flex flex-col">
        {/* Header */}
        <div className="flex items-center justify-between p-6 border-b border-border">
          <div>
            <h2 className="text-2xl font-quando font-bold">Experiment Report</h2>
            <p className="text-sm text-muted-foreground">
              ID: {report.experiment_id}
            </p>
          </div>
          <div className="flex items-center gap-2">
            <button
              onClick={handleDownloadJSON}
              className="flex items-center gap-2 px-4 py-2 border border-border rounded-lg hover:bg-accent transition font-quando text-sm"
            >
              <Download className="w-4 h-4" />
              JSON
            </button>
            <button
              onClick={handleDownloadPDF}
              className="flex items-center gap-2 px-4 py-2 border border-border rounded-lg hover:bg-accent transition font-quando text-sm"
            >
              <Download className="w-4 h-4" />
              PDF
            </button>
            <button
              onClick={onClose}
              className="p-2 hover:bg-accent rounded-lg transition"
            >
              <X className="w-5 h-5" />
            </button>
          </div>
        </div>

        {/* Content - Scrollable */}
        <div className="flex-1 overflow-auto p-6 space-y-6">
          {/* Molecule Information */}
          <div className="bg-card border border-border rounded-lg p-6">
            <h3 className="text-lg font-quando font-semibold mb-4">
              Molecular Structure
            </h3>
            <div className="grid grid-cols-2 gap-6">
              <div className="space-y-3 text-sm">
                <div className="flex justify-between">
                  <span className="text-muted-foreground">SMILES:</span>
                  <span className="font-mono">
                    {report.molecule.smiles || "N/A"}
                  </span>
                </div>
                <div className="flex justify-between">
                  <span className="text-muted-foreground">Formula:</span>
                  <span className="font-mono">
                    {report.molecule.formula || "N/A"}
                  </span>
                </div>
                <div className="flex justify-between">
                  <span className="text-muted-foreground">Charge:</span>
                  <span className="font-quando">{report.molecule.charge}</span>
                </div>
                <div className="flex justify-between">
                  <span className="text-muted-foreground">Multiplicity:</span>
                  <span className="font-quando">
                    {report.molecule.multiplicity}
                  </span>
                </div>
              </div>

              {report.molecule.atoms && report.molecule.atoms.length > 0 && (
                <div>
                  <div className="text-sm text-muted-foreground mb-2">
                    Atomic Coordinates:
                  </div>
                  <div className="bg-muted rounded p-3 font-mono text-xs max-h-48 overflow-auto">
                    {report.molecule.atoms.map((atom: any, i: number) => (
                      <div key={i}>
                        {atom.symbol} {atom.x.toFixed(6)} {atom.y.toFixed(6)}{" "}
                        {atom.z.toFixed(6)}
                      </div>
                    ))}
                  </div>
                </div>
              )}
            </div>
          </div>

          {/* Calculation Parameters */}
          <div className="bg-card border border-border rounded-lg p-6">
            <h3 className="text-lg font-quando font-semibold mb-4">
              Calculation Parameters
            </h3>
            <div className="grid grid-cols-3 gap-4 text-sm">
              <div className="flex justify-between">
                <span className="text-muted-foreground">Method:</span>
                <span className="font-quando font-medium">
                  {report.configuration.method}
                </span>
              </div>
              <div className="flex justify-between">
                <span className="text-muted-foreground">Ansatz:</span>
                <span className="font-quando font-medium">
                  {report.configuration.ansatz || "N/A"}
                </span>
              </div>
              <div className="flex justify-between">
                <span className="text-muted-foreground">Mapper:</span>
                <span className="font-quando font-medium">
                  {report.configuration.mapper || "N/A"}
                </span>
              </div>
              <div className="flex justify-between">
                <span className="text-muted-foreground">Backend:</span>
                <span className="font-quando font-medium">
                  {report.configuration.backend}
                </span>
              </div>
              <div className="flex justify-between">
                <span className="text-muted-foreground">Optimizer:</span>
                <span className="font-quando font-medium">
                  {report.configuration.optimizer || "N/A"}
                </span>
              </div>
              {report.configuration.hamiltonian && (
                <div className="flex justify-between">
                  <span className="text-muted-foreground">Hamiltonian:</span>
                  <span className="font-quando font-medium">
                    {report.configuration.hamiltonian}
                  </span>
                </div>
              )}
            </div>
          </div>

          {/* Results */}
          <div className="bg-card border border-border rounded-lg p-6">
            <h3 className="text-lg font-quando font-semibold mb-4">
              Energy Results
            </h3>
            <div className="grid grid-cols-2 gap-4">
              <div className="bg-gradient-to-br from-brand-orange/10 to-brand-orange/5 rounded-lg p-6 border border-brand-orange/20">
                <div className="text-xs text-muted-foreground mb-2">
                  Ground State Energy
                </div>
                <div className="text-3xl font-quando font-bold text-brand-orange">
                  {report.results.energy.toFixed(8)}
                  <span className="text-sm font-normal ml-2">Ha</span>
                </div>
              </div>

              {report.results.hf_energy !== undefined && (
                <div className="bg-muted rounded-lg p-6">
                  <div className="text-xs text-muted-foreground mb-2">
                    Hartree-Fock Energy
                  </div>
                  <div className="text-2xl font-quando font-bold">
                    {report.results.hf_energy.toFixed(8)}
                    <span className="text-sm font-normal ml-2">Ha</span>
                  </div>
                </div>
              )}

              {report.results.correlation_energy !== undefined && (
                <div className="bg-muted rounded-lg p-6">
                  <div className="text-xs text-muted-foreground mb-2">
                    Correlation Energy
                  </div>
                  <div className="text-2xl font-quando font-bold">
                    {report.results.correlation_energy.toFixed(8)}
                    <span className="text-sm font-normal ml-2">Ha</span>
                  </div>
                </div>
              )}

              {report.results.dipole_moment !== undefined && (
                <div className="bg-muted rounded-lg p-6">
                  <div className="text-xs text-muted-foreground mb-2">
                    Dipole Moment
                  </div>
                  <div className="text-2xl font-quando font-bold">
                    {report.results.dipole_moment.toFixed(6)}
                    <span className="text-sm font-normal ml-2">D</span>
                  </div>
                </div>
              )}
            </div>
          </div>

          {/* Convergence Information */}
          <div className="bg-card border border-border rounded-lg p-6">
            <h3 className="text-lg font-quando font-semibold mb-4">
              Convergence Analysis
            </h3>
            <div className="grid grid-cols-3 gap-4 mb-6">
              <div className="bg-muted rounded-lg p-4">
                <div className="text-xs text-muted-foreground mb-1">Status</div>
                <div
                  className={`text-xl font-quando font-bold ${
                    report.results.convergence.converged
                      ? "text-green-600"
                      : "text-red-600"
                  }`}
                >
                  {report.results.convergence.converged
                    ? "Converged"
                    : "Not Converged"}
                </div>
              </div>
              <div className="bg-muted rounded-lg p-4">
                <div className="text-xs text-muted-foreground mb-1">
                  Iterations
                </div>
                <div className="text-xl font-quando font-bold">
                  {report.results.convergence.iterations}
                </div>
              </div>
              {report.duration && (
                <div className="bg-muted rounded-lg p-4">
                  <div className="text-xs text-muted-foreground mb-1">
                    Duration
                  </div>
                  <div className="text-xl font-quando font-bold">
                    {(report.duration / 1000).toFixed(2)}s
                  </div>
                </div>
              )}
            </div>

            {/* Convergence Graph */}
            {report.convergence_data && report.convergence_data.length > 0 && (
              <div className="h-64">
                <ResponsiveContainer width="100%" height="100%">
                  <LineChart data={report.convergence_data}>
                    <CartesianGrid
                      strokeDasharray="3 3"
                      stroke="hsl(var(--border))"
                    />
                    <XAxis
                      dataKey="iteration"
                      label={{
                        value: "Iteration",
                        position: "insideBottom",
                        offset: -5,
                      }}
                      tick={{ fontSize: 11 }}
                      stroke="hsl(var(--muted-foreground))"
                    />
                    <YAxis
                      label={{
                        value: "Energy (Ha)",
                        angle: -90,
                        position: "insideLeft",
                      }}
                      tick={{ fontSize: 11 }}
                      stroke="hsl(var(--muted-foreground))"
                      domain={["auto", "auto"]}
                    />
                    <Tooltip
                      contentStyle={{
                        backgroundColor: "hsl(var(--card))",
                        border: "1px solid hsl(var(--border))",
                        borderRadius: "8px",
                        fontSize: "12px",
                      }}
                      formatter={(value: any) => [
                        value.toFixed(8) + " Ha",
                        "Energy",
                      ]}
                    />
                    <Line
                      type="monotone"
                      dataKey="energy"
                      stroke="#ea580c"
                      strokeWidth={2}
                      dot={{ r: 2 }}
                      animationDuration={300}
                    />
                  </LineChart>
                </ResponsiveContainer>
              </div>
            )}
          </div>

          {/* Circuit Information */}
          {report.circuit && (
            <div className="bg-card border border-border rounded-lg p-6">
              <h3 className="text-lg font-quando font-semibold mb-4">
                Quantum Circuit
              </h3>
              <div className="grid grid-cols-4 gap-4 mb-4">
                <div className="bg-muted rounded-lg p-3 text-center">
                  <div className="text-xs text-muted-foreground mb-1">Qubits</div>
                  <div className="text-xl font-quando font-bold">
                    {report.circuit.qubits}
                  </div>
                </div>
                <div className="bg-muted rounded-lg p-3 text-center">
                  <div className="text-xs text-muted-foreground mb-1">Depth</div>
                  <div className="text-xl font-quando font-bold">
                    {report.circuit.depth}
                  </div>
                </div>
                <div className="bg-muted rounded-lg p-3 text-center">
                  <div className="text-xs text-muted-foreground mb-1">Gates</div>
                  <div className="text-xl font-quando font-bold">
                    {report.circuit.gates}
                  </div>
                </div>
                <div className="bg-muted rounded-lg p-3 text-center">
                  <div className="text-xs text-muted-foreground mb-1">
                    Parameters
                  </div>
                  <div className="text-xl font-quando font-bold">
                    {report.circuit.parameters}
                  </div>
                </div>
              </div>
              <div className="bg-muted rounded-lg p-4 font-mono text-xs overflow-auto">
                <pre className="whitespace-pre">{report.circuit.diagram}</pre>
              </div>
            </div>
          )}

          {/* Additional Properties */}
          {report.results.properties &&
            Object.keys(report.results.properties).length > 0 && (
              <div className="bg-card border border-border rounded-lg p-6">
                <h3 className="text-lg font-quando font-semibold mb-4">
                  Additional Properties
                </h3>
                <div className="grid grid-cols-2 gap-4 text-sm">
                  {Object.entries(report.results.properties).map(
                    ([key, value]) => (
                      <div key={key} className="flex justify-between">
                        <span className="text-muted-foreground">
                          {key.replace(/_/g, " ").replace(/\b\w/g, (l) =>
                            l.toUpperCase()
                          )}
                          :
                        </span>
                        <span className="font-quando font-medium">
                          {typeof value === "number"
                            ? value.toFixed(6)
                            : String(value)}
                        </span>
                      </div>
                    )
                  )}
                </div>
              </div>
            )}
        </div>

        {/* Footer */}
        <div className="border-t border-border p-6">
          <div className="flex justify-between items-center text-sm text-muted-foreground">
            <div>
              Generated: {new Date(report.timestamp).toLocaleString()}
            </div>
            <button
              onClick={onClose}
              className="px-6 py-2 bg-brand-orange text-white rounded-lg hover:bg-brand-orange-dark transition font-quando"
            >
              Close
            </button>
          </div>
        </div>
      </div>
    </div>
  );
}
