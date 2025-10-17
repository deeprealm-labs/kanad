"use client";

import { useState, useEffect, useRef } from "react";
import {
  CheckCircle2,
  XCircle,
  Loader2,
  Clock,
  Download,
  ArrowLeft,
  Plus,
  Ban,
} from "lucide-react";
import {
  LineChart,
  Line,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  ResponsiveContainer,
} from "recharts";
import { useToast } from "@/components/ui/toast";
import * as api from "@/lib/api";
import type { Experiment, ConvergencePoint, ExperimentStatus } from "@/lib/types";
import CancelConfirmationDialog from "@/components/experiment/CancelConfirmationDialog";
import AnalysisResults from "@/components/simulation/AnalysisResults";

interface ExperimentMonitorProps {
  experimentId?: string | null;
  experimentConfig: any;
  onComplete: () => void;
  onBack: () => void;
  onQueueAnother?: () => void;
}

export default function ExperimentMonitor({
  experimentId,
  experimentConfig,
  onComplete,
  onBack,
  onQueueAnother,
}: ExperimentMonitorProps) {
  const [status, setStatus] = useState<ExperimentStatus>("queued");
  const [progress, setProgress] = useState(0);
  const [logs, setLogs] = useState<string[]>([
    "Initializing experiment...",
    "Validating molecular configuration...",
  ]);
  const [startTime] = useState(Date.now());
  const [elapsedTime, setElapsedTime] = useState(0);
  const [results, setResults] = useState<any>(null);
  const [convergenceData, setConvergenceData] = useState<ConvergencePoint[]>([]);
  const [currentIteration, setCurrentIteration] = useState(0);
  const [experiment, setExperiment] = useState<Experiment | null>(null);
  const [circuitData, setCircuitData] = useState<any>(null);
  const [showCancelDialog, setShowCancelDialog] = useState(false);
  const [isCancelling, setIsCancelling] = useState(false);

  const wsRef = useRef<WebSocket | null>(null);
  const pollingRef = useRef<(() => void) | null>(null);
  const previousStatusRef = useRef<ExperimentStatus>("queued");
  const toast = useToast();

  // Monitor experiment progress
  useEffect(() => {
    let timer: NodeJS.Timeout | null = null;

    // Only run timer if experiment is not completed/failed/cancelled
    if (status !== "completed" && status !== "failed" && status !== "cancelled") {
      timer = setInterval(() => {
        setElapsedTime(Math.floor((Date.now() - startTime) / 1000));
      }, 1000);
    }

    if (experimentId) {
      // Try WebSocket connection first
      try {
        const ws = api.createWebSocket(experimentId);
        wsRef.current = ws;

        ws.onopen = () => {
          console.log("WebSocket connected");
          addLog("Connected to real-time updates");
        };

        ws.onmessage = (event) => {
          try {
            const message = JSON.parse(event.data);
            handleWebSocketMessage(message);
          } catch (error) {
            console.error("Failed to parse WebSocket message:", error);
          }
        };

        ws.onerror = () => {
          // WebSocket not available, fall back to polling
          startPolling();
        };

        ws.onclose = () => {
          // Connection closed, polling will handle updates
        };
      } catch (error) {
        // WebSocket not supported, fall back to polling
        startPolling();
      }
    } else {
      // No experiment ID - this shouldn't happen
      console.error("ExperimentMonitor: No experiment ID provided");
      setStatus("failed");
      addLog("Error: No experiment ID provided");
    }

    return () => {
      if (timer) {
        clearInterval(timer);
      }
      if (wsRef.current) {
        wsRef.current.close();
      }
      if (pollingRef.current) {
        pollingRef.current();
      }
    };
  }, [experimentId, startTime, status]);

  const startPolling = () => {
    if (!experimentId) return;

    const stopPolling = api.pollExperimentStatus(
      experimentId,
      async (exp) => {
        const previousStatus = previousStatusRef.current;

        // Update refs first
        previousStatusRef.current = exp.status;

        setExperiment(exp);
        setStatus(exp.status);

        // Fetch circuit data when experiment starts or is running
        if ((exp.status === "running" || exp.status === "completed") && !circuitData) {
          try {
            const circuit = await api.getExperimentCircuit(experimentId);
            setCircuitData(circuit);
          } catch (error) {
            console.error("Failed to fetch circuit data:", error);
          }
        }

        if (exp.results) {
          console.log("Received results:", exp.results);
          setResults(exp.results);

          // Handle convergence data from different sources
          if (exp.results.convergence_history && Array.isArray(exp.results.convergence_history)) {
            console.log("Processing convergence_history:", exp.results.convergence_history.length, "points");
            setConvergenceData(exp.results.convergence_history);
          } else if (exp.results.energy_history && Array.isArray(exp.results.energy_history)) {
            console.log("Processing energy_history:", exp.results.energy_history.length, "points");
            // Convert energy_history array to convergence data format
            const convergence = exp.results.energy_history.map((energy: number, index: number) => ({
              iteration: index + 1,
              energy: energy
            }));
            setConvergenceData(convergence);
          } else if (exp.convergenceData) {
            console.log("Using exp.convergenceData");
            setConvergenceData(exp.convergenceData);
          } else if (exp.results.convergenceData) {
            console.log("Using exp.results.convergenceData");
            setConvergenceData(exp.results.convergenceData);
          } else {
            console.log("No convergence data found in:", Object.keys(exp.results));
          }
        }

        // Update progress based on convergence data or status
        if (exp.status === "running") {
          if (exp.convergenceData && exp.convergenceData.length > 0) {
            const maxIterations = 100; // default max
            const currentProgress = Math.min(
              90,
              10 + (exp.convergenceData.length / maxIterations) * 80
            );
            setProgress(currentProgress);
            setCurrentIteration(exp.convergenceData.length);
          } else {
            setProgress(Math.min(90, progress + 5));
          }
        } else if (exp.status === "completed") {
          setProgress(100);
          if (previousStatus !== "completed") {
            addLog(`Experiment completed successfully`);
          }
        } else if (exp.status === "failed" && previousStatus !== "failed") {
          addLog(`Experiment failed`);
        }

        // Only log status changes, not every poll
        if (previousStatus !== exp.status) {
          addLog(`Status: ${exp.status}`);
        }

        // Stop polling AFTER processing all data if experiment is done
        if (exp.status === "completed" || exp.status === "failed" || exp.status === "cancelled") {
          console.log("Stopping polling - experiment done, results processed");

          // Reset cancelling flag if experiment was cancelled
          if (exp.status === "cancelled") {
            setIsCancelling(false);
          }

          if (pollingRef.current) {
            pollingRef.current();
            pollingRef.current = null;
          }
        }
      },
      2000
    );

    pollingRef.current = stopPolling;
  };

  const handleWebSocketMessage = (message: any) => {
    switch (message.type) {
      case "status":
        setStatus(message.data.status);
        if (message.data.progress !== undefined) {
          setProgress(message.data.progress);
        }
        addLog(`Status updated: ${message.data.status}`);
        break;

      case "progress":
        setProgress(message.data.progress);
        break;

      case "convergence":
        setConvergenceData((prev) => [
          ...prev,
          {
            iteration: message.data.iteration,
            energy: message.data.energy,
          },
        ]);
        setCurrentIteration(message.data.iteration);
        break;

      case "log":
        addLog(message.data.message);
        break;

      case "result":
        setResults(message.data.results);
        setStatus("completed");
        setProgress(100);
        addLog("Experiment completed successfully!");
        toast.success("Experiment completed!");
        break;

      case "error":
        setStatus("failed");
        addLog(`Error: ${message.data.message}`);
        toast.error(message.data.message);
        break;
    }
  };

  const addLog = (message: string) => {
    setLogs((prev) => [...prev, message]);
  };

  const formatTime = (seconds: number) => {
    const mins = Math.floor(seconds / 60);
    const secs = seconds % 60;
    return `${mins}:${secs.toString().padStart(2, "0")}`;
  };

  const getStatusIcon = () => {
    switch (status) {
      case "queued":
        return <Clock className="w-6 h-6 text-yellow-500" />;
      case "running":
        return <Loader2 className="w-6 h-6 text-blue-500 animate-spin" />;
      case "completed":
        return <CheckCircle2 className="w-6 h-6 text-green-500" />;
      case "failed":
        return <XCircle className="w-6 h-6 text-red-500" />;
      case "cancelled":
        return <Ban className="w-6 h-6 text-orange-500" />;
    }
  };

  const getStatusText = () => {
    switch (status) {
      case "queued":
        return "Queued";
      case "running":
        return "Running";
      case "completed":
        return "Completed";
      case "failed":
        return "Failed";
      case "cancelled":
        return "Cancelled";
    }
  };

  const handleCancelClick = () => {
    setShowCancelDialog(true);
  };

  const handleCancelConfirm = async () => {
    if (!experimentId) return;

    try {
      setIsCancelling(true);
      const response = await api.cancelExperiment(experimentId);
      setShowCancelDialog(false);
      toast.success(response.message || "Experiment is being cancelled...");
      addLog("Cancellation requested - waiting for experiment to stop");

      // DO NOT stop polling here - let it continue until backend confirms cancellation
      // The polling callback will detect the "cancelled" status and stop automatically
    } catch (error: any) {
      console.error("Failed to cancel experiment:", error);
      toast.error(error.message || "Failed to cancel experiment");
      setIsCancelling(false);
    }
    // Note: setIsCancelling(false) will be called when polling detects "cancelled" status
  };

  const handleExport = async () => {
    if (!experimentId) {
      toast.warning("Cannot export simulated results");
      return;
    }

    try {
      await api.downloadExperiment(experimentId, "json");
      toast.success("Results exported successfully!");
    } catch (error: any) {
      toast.error(error.message || "Failed to export results");
    }
  };

  return (
    <div className="h-full flex flex-col bg-background">
      {/* Cancel Confirmation Dialog */}
      <CancelConfirmationDialog
        isOpen={showCancelDialog}
        onClose={() => setShowCancelDialog(false)}
        onConfirm={handleCancelConfirm}
        experimentName={experiment?.molecule?.smiles || experimentConfig?.molecule?.smiles || "Custom Experiment"}
        isRunning={status === "running"}
        isLoading={isCancelling}
      />

      {/* Header - Fixed height */}
      <div className="flex items-center justify-between px-6 py-4 border-b border-border">
        <div className="flex items-center gap-4">
          <button
            onClick={onBack}
            className="p-2 hover:bg-accent rounded-md transition"
          >
            <ArrowLeft className="w-5 h-5" />
          </button>
          <h1 className="text-xl font-quando font-bold">Experiment Monitor</h1>
        </div>
        <div className="flex items-center gap-6">
          <div className="flex items-center gap-2">
            {getStatusIcon()}
            <div>
              <div className="text-sm font-quando font-semibold">
                {getStatusText()}
              </div>
              <div className="text-xs text-muted-foreground font-quando">
                {formatTime(elapsedTime)}
              </div>
            </div>
          </div>
          <div className="flex items-center gap-2">
            {(status === "running" || status === "queued") && experimentId && (
              <button
                onClick={handleCancelClick}
                disabled={isCancelling}
                className="flex items-center gap-2 px-4 py-2 text-sm bg-red-600 text-white rounded-lg hover:bg-red-700 transition font-quando disabled:opacity-50 disabled:cursor-not-allowed"
              >
                {isCancelling ? (
                  <>
                    <Loader2 className="w-4 h-4 animate-spin" />
                    Cancelling...
                  </>
                ) : (
                  <>
                    <Ban className="w-4 h-4" />
                    Cancel
                  </>
                )}
              </button>
            )}
            {onQueueAnother && status === "running" && (
              <button
                onClick={onQueueAnother}
                className="flex items-center gap-2 px-4 py-2 text-sm border-2 border-brand-orange text-brand-orange rounded-lg hover:bg-brand-orange hover:text-white transition font-quando"
              >
                <Plus className="w-4 h-4" />
                Queue Another
              </button>
            )}
            {status === "completed" && (
              <button
                onClick={handleExport}
                className="flex items-center gap-2 px-4 py-2 text-sm border border-border rounded-lg hover:bg-accent transition font-quando"
              >
                <Download className="w-4 h-4" />
                Export
              </button>
            )}
          </div>
        </div>
      </div>

      {/* Dashboard Grid - No scrolling, everything fits */}
      <div className="flex-1 p-6 grid grid-cols-3 gap-4 min-h-0">
        {/* Left Column - Configuration & Progress */}
        <div className="flex flex-col gap-4 min-h-0">
          {/* Configuration Card */}
          <div className="bg-card border border-border rounded-lg p-4">
            <h3 className="text-sm font-quando font-semibold mb-3">
              Configuration
            </h3>
            <div className="space-y-2 text-xs">
              <div className="flex justify-between">
                <span className="text-muted-foreground">Method:</span>
                <span className="font-quando font-medium">
                  {experimentConfig?.backendSettings?.method || "VQE"}
                </span>
              </div>
              <div className="flex justify-between">
                <span className="text-muted-foreground">Backend:</span>
                <span className="font-quando font-medium">
                  {experimentConfig?.backendSettings?.backend || "Classical"}
                </span>
              </div>
              <div className="flex justify-between">
                <span className="text-muted-foreground">Ansatz:</span>
                <span className="font-quando font-medium">
                  {experimentConfig?.backendSettings?.ansatz || "HEA"}
                </span>
              </div>
              <div className="flex justify-between">
                <span className="text-muted-foreground">Basis:</span>
                <span className="font-quando font-medium">
                  {experimentConfig?.molecule?.basis || "STO-3G"}
                </span>
              </div>
            </div>
          </div>

          {/* Progress Card */}
          <div className="bg-card border border-border rounded-lg p-4">
            <div className="flex justify-between items-center mb-2">
              <h3 className="text-sm font-quando font-semibold">Progress</h3>
              <span className="text-xs font-quando text-muted-foreground">
                {progress}%
              </span>
            </div>
            <div className="w-full bg-muted rounded-full h-2 overflow-hidden">
              <div
                className="h-full bg-brand-orange transition-all duration-300 ease-out"
                style={{ width: `${progress}%` }}
              />
            </div>
            {status === "running" && currentIteration > 0 && (
              <div className="mt-3 text-xs text-muted-foreground font-quando">
                Iteration {currentIteration} / 42
              </div>
            )}
          </div>

          {/* Results & Analysis */}
          <div className="bg-card border border-border rounded-lg p-4 flex-1 min-h-0 flex flex-col overflow-hidden">
            <h3 className="text-sm font-quando font-semibold mb-3 flex-shrink-0">
              {status === "completed" ? "Analysis Results" : "Live Metrics"}
            </h3>
            <div className="flex-1 overflow-y-auto min-h-0">
              {status === "completed" && results ? (
                <AnalysisResults results={results} />
              ) : (
                <div className="space-y-3">
                  <div className="bg-muted rounded-lg p-3">
                    <div className="text-xs text-muted-foreground mb-1">
                      Current Energy
                    </div>
                    <div className="text-lg font-quando font-bold">
                      {convergenceData.length > 0
                        ? convergenceData[convergenceData.length - 1].energy.toFixed(6)
                        : "---"}{" "}
                      <span className="text-sm font-normal">Ha</span>
                    </div>
                  </div>
                  {convergenceData.length > 0 && (
                    <div className="bg-muted rounded-lg p-3">
                      <div className="text-xs text-muted-foreground mb-1">
                        Current Iteration
                      </div>
                      <div className="text-lg font-quando font-bold">
                        {convergenceData.length}
                      </div>
                    </div>
                  )}
                </div>
              )}
            </div>
          </div>
        </div>

        {/* Center Column - Convergence Chart */}
        <div className="bg-card border border-border rounded-lg p-4 min-h-0 flex flex-col">
          <h3 className="text-sm font-quando font-semibold mb-3">
            Energy Convergence
          </h3>
          <div className="flex-1 min-h-0">
            <ResponsiveContainer width="100%" height="100%">
              <LineChart data={convergenceData}>
                <CartesianGrid strokeDasharray="3 3" stroke="hsl(var(--border))" />
                <XAxis
                  dataKey="iteration"
                  label={{ value: "Iteration", position: "insideBottom", offset: -5 }}
                  tick={{ fontSize: 11 }}
                  stroke="hsl(var(--muted-foreground))"
                />
                <YAxis
                  label={{ value: "Energy (Ha)", angle: -90, position: "insideLeft" }}
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
                  formatter={(value: any) => [value.toFixed(6) + " Ha", "Energy"]}
                />
                <Line
                  type="monotone"
                  dataKey="energy"
                  stroke="#ea580c"
                  strokeWidth={2}
                  dot={false}
                  animationDuration={300}
                />
              </LineChart>
            </ResponsiveContainer>
          </div>
        </div>

        {/* Right Column - Logs */}
        <div className="bg-card border border-border rounded-lg p-4 min-h-0 flex flex-col">
          <h3 className="text-sm font-quando font-semibold mb-3">
            Execution Logs
          </h3>
          <div className="flex-1 bg-muted rounded-lg p-3 font-mono text-xs overflow-auto min-h-0">
            <div className="space-y-1">
              {logs.map((log, i) => (
                <div key={i} className="text-foreground">
                  <span className="text-muted-foreground mr-2">
                    [{formatTime(Math.floor((i + 1) * 1.5))}]
                  </span>
                  {log}
                </div>
              ))}
            </div>
          </div>
        </div>
      </div>

      {/* Footer Actions - Only show when completed */}
      {status === "completed" && (
        <div className="border-t border-border px-6 py-4">
          <div className="flex gap-4 justify-end">
            <button
              onClick={onBack}
              className="px-6 py-2.5 text-sm border border-border rounded-lg hover:bg-accent transition font-quando"
            >
              New Experiment
            </button>
            <button
              onClick={onComplete}
              className="px-6 py-2.5 text-sm bg-brand-orange text-white rounded-lg hover:bg-brand-orange-dark transition font-quando"
            >
              View in Dashboard
            </button>
          </div>
        </div>
      )}
    </div>
  );
}
