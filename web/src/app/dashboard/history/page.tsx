"use client";

import { useState, useEffect } from "react";
import {
  ArrowLeft,
  Search,
  Download,
  Trash2,
  Filter,
  CheckCircle2,
  Clock,
  Play,
  FileText,
  Loader2,
  Ban,
  XCircle,
} from "lucide-react";
import Link from "next/link";
import * as api from "@/lib/api";
import ExperimentReport from "@/components/experiment/ExperimentReport";
import { useToast } from "@/components/ui/toast";

export default function HistoryPage() {
  const [experiments, setExperiments] = useState<any[]>([]);
  const [searchTerm, setSearchTerm] = useState("");
  const [filterStatus, setFilterStatus] = useState<string>("all");
  const [selectedExperiment, setSelectedExperiment] = useState<any>(null);
  const [showReport, setShowReport] = useState(false);
  const [loading, setLoading] = useState(true);
  const toast = useToast();

  useEffect(() => {
    loadExperiments();
  }, []);

  const loadExperiments = async () => {
    try {
      setLoading(true);
      const response = await api.getExperiments({ limit: 100 });
      setExperiments(response.experiments);
    } catch (error: any) {
      console.error("Failed to load experiments from API:", error);

      // Fallback to localStorage
      const saved = localStorage.getItem("kanad_experiments");
      if (saved) {
        try {
          setExperiments(JSON.parse(saved));
        } catch (e) {
          console.error("Failed to parse localStorage experiments:", e);
        }
      }

      if (error.statusCode !== 503) {
        toast.error("Failed to load experiments. Using local data.");
      }
    } finally {
      setLoading(false);
    }
  };

  const filteredExperiments = experiments.filter((exp) => {
    const matchesSearch =
      exp.molecule?.smiles?.toLowerCase().includes(searchTerm.toLowerCase()) ||
      exp.method?.toLowerCase().includes(searchTerm.toLowerCase()) ||
      exp.backend?.toLowerCase().includes(searchTerm.toLowerCase());

    const matchesFilter =
      filterStatus === "all" || exp.status === filterStatus;

    return matchesSearch && matchesFilter;
  });

  const handleDelete = async (id: string) => {
    try {
      await api.deleteExperiment(id);
      const updated = experiments.filter((e) => e.id !== id);
      setExperiments(updated);
      if (selectedExperiment?.id === id) {
        setSelectedExperiment(null);
      }
      toast.success("Experiment deleted successfully");
    } catch (error: any) {
      console.error("Failed to delete experiment:", error);

      // Fallback to localStorage
      const updated = experiments.filter((e) => e.id !== id);
      setExperiments(updated);
      localStorage.setItem("kanad_experiments", JSON.stringify(updated));
      if (selectedExperiment?.id === id) {
        setSelectedExperiment(null);
      }

      toast.warning("Experiment deleted locally (offline mode)");
    }
  };

  const handleDownload = async (id: string) => {
    try {
      await api.downloadExperiment(id, "json");
      toast.success("Experiment exported successfully!");
    } catch (error: any) {
      toast.error(error.message || "Failed to export experiment");
    }
  };

  const getStatusIcon = (status: string) => {
    switch (status) {
      case "completed":
        return <CheckCircle2 className="w-5 h-5 text-green-600" />;
      case "running":
        return <Play className="w-5 h-5 text-blue-600" />;
      case "queued":
        return <Clock className="w-5 h-5 text-yellow-600" />;
      case "cancelled":
        return <Ban className="w-5 h-5 text-orange-600" />;
      case "failed":
        return <XCircle className="w-5 h-5 text-red-600" />;
      default:
        return null;
    }
  };

  if (loading) {
    return (
      <div className="h-full flex items-center justify-center bg-background">
        <div className="text-center">
          <Loader2 className="w-12 h-12 animate-spin text-brand-orange mx-auto mb-4" />
          <p className="text-muted-foreground font-quando">Loading experiments...</p>
        </div>
      </div>
    );
  }

  return (
    <div className="h-full flex flex-col bg-background">
      {showReport && selectedExperiment && (
        <ExperimentReport
          experimentId={selectedExperiment.id}
          onClose={() => setShowReport(false)}
        />
      )}

      {/* Header */}
      <div className="flex items-center justify-between px-6 py-4 border-b border-border">
        <div className="flex items-center gap-4">
          <Link
            href="/dashboard"
            className="p-2 hover:bg-accent rounded-md transition"
          >
            <ArrowLeft className="w-5 h-5" />
          </Link>
          <h1 className="text-2xl font-quando font-bold">
            History
          </h1>
        </div>
        <div className="flex items-center gap-3">
          <div className="relative">
            <Search className="absolute left-3 top-1/2 -translate-y-1/2 w-4 h-4 text-muted-foreground" />
            <input
              type="text"
              placeholder="Search experiments..."
              value={searchTerm}
              onChange={(e) => setSearchTerm(e.target.value)}
              className="pl-10 pr-4 py-2 border border-border rounded-lg bg-background focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando text-sm"
            />
          </div>
          <select
            value={filterStatus}
            onChange={(e) => setFilterStatus(e.target.value)}
            className="px-4 py-2 border border-border rounded-lg bg-background focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando text-sm"
          >
            <option value="all">All Status</option>
            <option value="completed">Completed</option>
            <option value="running">Running</option>
            <option value="queued">Queued</option>
            <option value="cancelled">Cancelled</option>
            <option value="failed">Failed</option>
          </select>
        </div>
      </div>

      {/* Main Content */}
      <div className="flex-1 grid grid-cols-2 gap-6 p-6 min-h-0">
        {/* Left: Experiments List */}
        <div className="flex flex-col min-h-0">
          <div className="mb-4">
            <p className="text-sm text-muted-foreground font-quando">
              {filteredExperiments.length} experiment
              {filteredExperiments.length !== 1 ? "s" : ""} found
            </p>
          </div>

          <div className="flex-1 overflow-auto space-y-3 min-h-0">
            {filteredExperiments.length === 0 ? (
              <div className="bg-card border border-border rounded-lg p-12 flex flex-col items-center justify-center text-center">
                <Search className="w-16 h-16 text-muted-foreground mb-4" />
                <p className="text-lg font-quando font-semibold mb-2">
                  No experiments found
                </p>
                <p className="text-sm text-muted-foreground">
                  Try adjusting your search or filters
                </p>
              </div>
            ) : (
              filteredExperiments.map((exp) => (
                <div
                  key={exp.id}
                  onClick={() => setSelectedExperiment(exp)}
                  className={`bg-card border rounded-lg p-4 cursor-pointer transition ${
                    selectedExperiment?.id === exp.id
                      ? "border-brand-orange"
                      : "border-border hover:border-brand-orange/50"
                  }`}
                >
                  <div className="flex items-start justify-between mb-3">
                    <div className="flex items-center gap-3">
                      {getStatusIcon(exp.status)}
                      <div>
                        <h3 className="font-quando font-semibold">
                          {exp.molecule?.smiles || "Custom Molecule"}
                        </h3>
                        <p className="text-xs text-muted-foreground">
                          {new Date(exp.timestamp).toLocaleString()}
                        </p>
                      </div>
                    </div>
                  </div>

                  <div className="grid grid-cols-2 gap-2 text-xs">
                    <div>
                      <span className="text-muted-foreground">Method:</span>
                      <span className="ml-2 font-quando font-medium">
                        {exp.method}
                      </span>
                    </div>
                    <div>
                      <span className="text-muted-foreground">Backend:</span>
                      <span className="ml-2 font-quando font-medium">
                        {exp.backend}
                      </span>
                    </div>
                  </div>

                  {exp.results && (
                    <div className="mt-3 pt-3 border-t border-border">
                      <div className="flex items-center justify-between">
                        <span className="text-xs text-muted-foreground">
                          Energy:
                        </span>
                        <span className="font-mono text-sm font-bold">
                          {exp.results.energy.toFixed(6)} Ha
                        </span>
                      </div>
                    </div>
                  )}
                </div>
              ))
            )}
          </div>
        </div>

        {/* Right: Experiment Details */}
        <div className="bg-card border border-border rounded-lg p-6 flex flex-col min-h-0">
          {selectedExperiment ? (
            <>
              <div className="flex items-center justify-between mb-6">
                <h2 className="text-xl font-quando font-bold">
                  Experiment Details
                </h2>
                <div className="flex gap-2">
                  {selectedExperiment.status === "completed" && (
                    <button
                      onClick={() => setShowReport(true)}
                      className="flex items-center gap-2 px-3 py-2 bg-brand-orange text-white rounded-lg hover:bg-brand-orange-dark transition font-quando text-sm"
                    >
                      <FileText className="w-4 h-4" />
                      View Report
                    </button>
                  )}
                  <button
                    onClick={() => handleDownload(selectedExperiment.id)}
                    className="p-2 border border-border rounded-lg hover:bg-accent transition"
                    title="Download JSON"
                  >
                    <Download className="w-4 h-4" />
                  </button>
                  <button
                    onClick={() => handleDelete(selectedExperiment.id)}
                    className="p-2 border border-border rounded-lg hover:bg-red-100 dark:hover:bg-red-900 transition text-red-600"
                    title="Delete"
                  >
                    <Trash2 className="w-4 h-4" />
                  </button>
                </div>
              </div>

              <div className="flex-1 overflow-auto space-y-6 min-h-0">
                {/* Status */}
                <div>
                  <h3 className="text-sm font-quando font-semibold mb-2 text-muted-foreground">
                    Status
                  </h3>
                  <div className="flex items-center gap-2">
                    {getStatusIcon(selectedExperiment.status)}
                    <span className="font-quando capitalize">
                      {selectedExperiment.status}
                    </span>
                  </div>
                </div>

                {/* Molecule */}
                <div>
                  <h3 className="text-sm font-quando font-semibold mb-2 text-muted-foreground">
                    Molecule
                  </h3>
                  <div className="bg-muted rounded-lg p-4 space-y-2 text-sm">
                    <div className="flex justify-between">
                      <span>SMILES:</span>
                      <span className="font-mono">
                        {selectedExperiment.molecule?.smiles || "N/A"}
                      </span>
                    </div>
                    <div className="flex justify-between">
                      <span>Basis Set:</span>
                      <span className="font-quando">
                        {selectedExperiment.molecule?.basis || "sto-3g"}
                      </span>
                    </div>
                    <div className="flex justify-between">
                      <span>Charge:</span>
                      <span className="font-quando">
                        {selectedExperiment.molecule?.charge || 0}
                      </span>
                    </div>
                    <div className="flex justify-between">
                      <span>Multiplicity:</span>
                      <span className="font-quando">
                        {selectedExperiment.molecule?.multiplicity || 1}
                      </span>
                    </div>
                  </div>
                </div>

                {/* Configuration */}
                <div>
                  <h3 className="text-sm font-quando font-semibold mb-2 text-muted-foreground">
                    Configuration
                  </h3>
                  <div className="bg-muted rounded-lg p-4 space-y-2 text-sm">
                    <div className="flex justify-between">
                      <span>Method:</span>
                      <span className="font-quando">
                        {selectedExperiment.method}
                      </span>
                    </div>
                    <div className="flex justify-between">
                      <span>Backend:</span>
                      <span className="font-quando">
                        {selectedExperiment.backend}
                      </span>
                    </div>
                    {selectedExperiment.backendSettings?.ansatz && (
                      <div className="flex justify-between">
                        <span>Ansatz:</span>
                        <span className="font-quando">
                          {selectedExperiment.backendSettings.ansatz}
                        </span>
                      </div>
                    )}
                  </div>
                </div>

                {/* Results */}
                {selectedExperiment.results && (
                  <div>
                    <h3 className="text-sm font-quando font-semibold mb-2 text-muted-foreground">
                      Results
                    </h3>
                    <div className="space-y-3">
                      <div className="bg-muted rounded-lg p-4">
                        <div className="text-xs text-muted-foreground mb-1">
                          Ground State Energy
                        </div>
                        <div className="text-2xl font-quando font-bold">
                          {selectedExperiment.results.energy.toFixed(6)} Ha
                        </div>
                        <div className="text-xs text-muted-foreground mt-1">
                          Converged in {selectedExperiment.results.iterations}{" "}
                          iterations
                        </div>
                      </div>

                      {selectedExperiment.results.dipoleMoment !== undefined && (
                        <div className="bg-muted rounded-lg p-4">
                          <div className="text-xs text-muted-foreground mb-1">
                            Dipole Moment
                          </div>
                          <div className="text-xl font-quando font-bold">
                            {selectedExperiment.results.dipoleMoment.toFixed(4)}{" "}
                            D
                          </div>
                        </div>
                      )}

                      <div className="bg-muted rounded-lg p-4">
                        <div className="text-xs text-muted-foreground mb-1">
                          Convergence Status
                        </div>
                        <div
                          className={`text-xl font-quando font-bold ${
                            selectedExperiment.results.converged
                              ? "text-green-600"
                              : "text-red-600"
                          }`}
                        >
                          {selectedExperiment.results.converged
                            ? "Converged"
                            : "Not Converged"}
                        </div>
                      </div>
                    </div>
                  </div>
                )}

                {/* Timestamp */}
                <div>
                  <h3 className="text-sm font-quando font-semibold mb-2 text-muted-foreground">
                    Timestamp
                  </h3>
                  <p className="text-sm font-quando">
                    {new Date(selectedExperiment.timestamp).toLocaleString()}
                  </p>
                </div>
              </div>
            </>
          ) : (
            <div className="flex-1 flex items-center justify-center text-center">
              <div>
                <Filter className="w-16 h-16 text-muted-foreground mx-auto mb-4" />
                <p className="text-lg font-quando font-semibold mb-2">
                  Select an experiment
                </p>
                <p className="text-sm text-muted-foreground">
                  Click on an experiment to view details
                </p>
              </div>
            </div>
          )}
        </div>
      </div>
    </div>
  );
}
