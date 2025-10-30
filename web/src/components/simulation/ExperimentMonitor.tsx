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
  BarChart,
  Bar,
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
  onCancel?: () => void;
}

export default function ExperimentMonitor({
  experimentId,
  experimentConfig,
  onComplete,
  onBack,
  onQueueAnother,
  onCancel,
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
  const wsConnectedRef = useRef<boolean>(false);
  const pollingRef = useRef<(() => void) | null>(null);
  const previousStatusRef = useRef<ExperimentStatus>("queued");
  const toast = useToast();

  // Initialize or update when experiment ID changes
  useEffect(() => {
    console.log("🔄 COMPONENT MOUNT/UPDATE - Experiment ID:", experimentId);
    console.log("🔄 Current state - iteration:", currentIteration, "convergence points:", convergenceData.length);

    // DON'T reset state if we're remounting with the same experiment
    // Only reset if it's truly a NEW experiment (different ID)
    if (!experimentId) {
      console.log("❌ No experiment ID - cannot initialize");
      return;
    }

    // Close any existing WebSocket
    if (wsRef.current) {
      console.log("🔌 Closing previous WebSocket");
      wsRef.current.close();
      wsRef.current = null;
    }

    // Stop any existing polling
    if (pollingRef.current) {
      console.log("⏹️  Stopping previous polling");
      pollingRef.current();
      pollingRef.current = null;
    }

    // Fetch experiment data to restore state (if reopening an experiment)
    console.log("📥 Fetching experiment data to restore state...");
    api.getExperiment(experimentId).then((response) => {
      const exp = response.experiment || response;
      console.log("📥 Loaded experiment:", exp.id, "status:", exp.status);

      // Update experiment state
      setExperiment(exp);
      setStatus(exp.status);

      // Restore progress from job data if available
      if (exp.job && typeof exp.job.progress === 'number') {
        console.log("📊 Restoring progress from job:", exp.job.progress);
        setProgress(exp.job.progress);
        if (exp.job.current_iteration) {
          setCurrentIteration(exp.job.current_iteration);
        }
      }

      // Restore convergence data from results if available
      if (exp.results) {
        setResults(exp.results);

        // Only restore convergence data if we don't have it yet
        if (convergenceData.length === 0) {
          if (exp.results.convergence_history && Array.isArray(exp.results.convergence_history)) {
            console.log("📊 Restoring convergence_history:", exp.results.convergence_history.length, "points");
            setConvergenceData(exp.results.convergence_history);
            if (exp.results.convergence_history.length > 0) {
              setCurrentIteration(exp.results.convergence_history.length);
            }
          } else if (exp.results.energy_history && Array.isArray(exp.results.energy_history)) {
            console.log("📊 Restoring energy_history:", exp.results.energy_history.length, "points");
            const convergence = exp.results.energy_history.map((energy: number, index: number) => ({
              iteration: index + 1,
              energy: energy
            }));
            setConvergenceData(convergence);
            setCurrentIteration(convergence.length);
          }
        } else {
          console.log("📊 Keeping existing convergence data:", convergenceData.length, "points");
        }
      }

      console.log("✅ State restored for experiment:", experimentId);
    }).catch((error) => {
      console.error("❌ Failed to fetch experiment data:", error);
      // Reset to clean state if fetch fails
      setCurrentIteration(0);
      setConvergenceData([]);
      setProgress(0);
      setLogs(["Initializing experiment..."]);
      setStatus("queued");
      setResults(null);
    });
  }, [experimentId]);

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
      // Reset WebSocket connection tracking
      wsConnectedRef.current = false;

      // Try WebSocket connection first
      try {
        console.log("🔌 Attempting WebSocket connection for experiment:", experimentId);
        const ws = api.createWebSocket(experimentId);
        wsRef.current = ws;

        ws.onopen = () => {
          console.log("✅ WebSocket connected successfully");
          wsConnectedRef.current = true;
          addLog("✅ Connected to real-time updates");
        };

        ws.onmessage = (event) => {
          try {
            const message = JSON.parse(event.data);
            console.log("📨 WebSocket message received:", message);
            handleWebSocketMessage(message);
          } catch (error) {
            console.error("❌ Failed to parse WebSocket message:", error);
          }
        };

        ws.onerror = (error) => {
          console.error("❌ WebSocket error:", error);
          // Don't immediately fall back to polling - wait for onclose
          // Sometimes errors are transient
        };

        ws.onclose = (event) => {
          console.log("🔌 WebSocket closed:", event.code, event.reason, "wasConnected:", wsConnectedRef.current);

          // Only show warning and start polling if WebSocket never connected
          // Give a small delay to ensure onopen event had time to fire
          setTimeout(() => {
            if (!wsConnectedRef.current && !pollingRef.current) {
              console.log("⚠️  WebSocket never connected, falling back to polling");
              addLog("⚠️  WebSocket failed, using polling");
              startPolling();
            } else {
              console.log("✅ WebSocket closed normally after successful connection");
            }
          }, 100);
        };
      } catch (error) {
        console.error("❌ WebSocket creation failed:", error);
        addLog("⚠️  WebSocket not available, using polling");
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
          console.log("📊 Received results:", exp.results);
          console.log("📊 Has analysis in results?", !!exp.results.analysis);
          if (exp.results.analysis) {
            console.log("📊 Analysis keys:", Object.keys(exp.results.analysis));
          }
          setResults(exp.results);

          // Handle convergence data from different sources
          // Only replace if we don't already have WebSocket data
          setConvergenceData((prevData) => {
            if (prevData.length > 0) {
              console.log("Keeping existing WebSocket convergence data:", prevData.length, "points");
              return prevData;
            }

            if (exp.results.convergence_history && Array.isArray(exp.results.convergence_history)) {
              console.log("Processing convergence_history:", exp.results.convergence_history.length, "points");
              return exp.results.convergence_history;
            } else if (exp.results.energy_history && Array.isArray(exp.results.energy_history)) {
              console.log("Processing energy_history:", exp.results.energy_history.length, "points");
              // Convert energy_history array to convergence data format
              const convergence = exp.results.energy_history.map((energy: number, index: number) => ({
                iteration: index + 1,
                energy: energy
              }));
              return convergence;
            } else if (exp.convergenceData) {
              console.log("Using exp.convergenceData");
              return exp.convergenceData;
            } else if (exp.results.convergenceData) {
              console.log("Using exp.results.convergenceData");
              return exp.results.convergenceData;
            } else {
              console.log("No convergence data found in:", Object.keys(exp.results));
              return prevData;
            }
          });
        }

        // Update progress from job data (if available) or fallback to convergence-based calculation
        if (exp.status === "running") {
          // PRIORITY 1: Use job progress if available (most accurate)
          if (exp.job && typeof exp.job.progress === 'number') {
            console.log("📊 Using job.progress:", exp.job.progress);
            setProgress(Math.min(90, exp.job.progress));

            // Update iteration from job data
            if (exp.job.current_iteration) {
              setCurrentIteration((prev) => Math.max(prev, exp.job.current_iteration));
            }
          }
          // FALLBACK: Calculate from convergence data if job data not available
          else if (exp.convergenceData && exp.convergenceData.length > 0) {
            console.log("⚠️  No job.progress, calculating from convergence data");
            // Get actual max iterations from experiment config
            const method = exp.config?.backendSettings?.method;
            let maxIterations = 100;

            if (method === "VQE") {
              maxIterations = exp.config?.backendSettings?.maxIterations || 100;
            } else if (method === "SQD") {
              maxIterations = 6; // SQD has 7 stages (0-6)
            } else if (method === "EXCITED_STATES") {
              maxIterations = exp.config?.backendSettings?.nStates || 5;
            }

            const currentProgress = Math.min(
              90,
              10 + (exp.convergenceData.length / maxIterations) * 80
            );
            setProgress(currentProgress);
            // Only update iteration if it's higher (don't go backwards with stale data)
            setCurrentIteration((prev) => Math.max(prev, exp.convergenceData.length));
          } else {
            console.log("⚠️  No progress data available, incrementing slowly");
            setProgress(Math.min(90, progress + 5));
          }
        } else if (exp.status === "completed") {
          setProgress(100);
          // Set final iteration count based on convergence data from results
          if (exp.results) {
            let finalIterations = 0;
            if (exp.results.convergence_history && Array.isArray(exp.results.convergence_history)) {
              finalIterations = exp.results.convergence_history.length;
              console.log("📊 Setting final iterations from convergence_history:", finalIterations);
            } else if (exp.results.energy_history && Array.isArray(exp.results.energy_history)) {
              finalIterations = exp.results.energy_history.length;
              console.log("📊 Setting final iterations from energy_history:", finalIterations);
            }
            if (finalIterations > 0) {
              console.log("✅ Setting currentIteration to:", finalIterations);
              setCurrentIteration(finalIterations);
            } else {
              console.log("⚠️  No convergence data found in results");
            }
          }
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

          // Reset cancelling flag and notify parent if experiment was cancelled
          if (exp.status === "cancelled") {
            setIsCancelling(false);
            // Notify parent component so it can refresh experiment list
            if (onCancel) {
              onCancel();
            }
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
        setStatus(message.status);
        if (message.progress !== undefined) {
          setProgress(message.progress);
        }
        addLog(`Status updated: ${message.status}`);

        // If experiment completed, fetch final results
        if (message.status === 'completed' && experimentId) {
          console.log("✅ Experiment completed - fetching final results");
          api.getExperiment(experimentId).then((response) => {
            const exp = response.experiment || response;
            console.log("📊 Received results:", exp.results);
            console.log("📊 Has analysis in results?", !!exp.results?.analysis);
            console.log("📊 Has advanced_analysis in results?", !!exp.results?.advanced_analysis);
            if (exp.results?.analysis) {
              console.log("📊 Analysis keys:", Object.keys(exp.results.analysis));
            }
            if (exp.results?.advanced_analysis) {
              console.log("📊 Advanced analysis profile:", exp.results.advanced_analysis.profile);
              console.log("📊 Advanced analysis status:", exp.results.advanced_analysis.status);
              console.log("📊 Advanced analysis results:", exp.results.advanced_analysis.results);
            }

            if (exp.results) {
              setResults(exp.results);

              // Handle convergence data - only use if we don't have WebSocket data
              setConvergenceData((prevData) => {
                if (prevData.length > 0) {
                  console.log("Keeping existing WebSocket data on completion:", prevData.length, "points");
                  return prevData;
                }

                if (exp.results.convergence_history && Array.isArray(exp.results.convergence_history)) {
                  setCurrentIteration(exp.results.convergence_history.length);
                  return exp.results.convergence_history;
                } else if (exp.results.energy_history && Array.isArray(exp.results.energy_history)) {
                  const convergence = exp.results.energy_history.map((energy: number, index: number) => ({
                    iteration: index + 1,
                    energy: energy
                  }));
                  setCurrentIteration(convergence.length);
                  return convergence;
                }
                return prevData;
              });
            }
          }).catch((error) => {
            console.error("❌ Failed to fetch final results:", error);
          });
        }
        break;

      case "progress":
        setProgress(message.progress);
        break;

      case "convergence":
        console.log("📊 WebSocket convergence update - iteration:", message.iteration, "energy:", message.energy);
        setConvergenceData((prev) => [
          ...prev,
          {
            iteration: message.iteration,
            energy: message.energy,
          },
        ]);
        console.log("✅ Setting currentIteration to:", message.iteration);
        setCurrentIteration(message.iteration);

        // Add log for convergence point
        // For SQD (low iteration count), log every stage. For VQE, log every 10 iterations
        const shouldLog = message.iteration < 10 || message.iteration % 10 === 0;
        if (shouldLog) {
          const estIter = message.is_optimizer_iteration
            ? message.iteration
            : Math.floor(message.iteration / 40);
          addLog(`Function eval ${message.iteration} (~iter ${estIter}): E = ${message.energy.toFixed(8)} Ha`);
        }
        break;

      case "log":
        addLog(message.message);
        break;

      case "result":
        setResults(message.results);
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
          <div>
            <h1 className="text-xl font-quando font-bold">Experiment Monitor</h1>
            {(experiment?.molecule || experimentConfig?.molecule) && (
              <div className="flex items-center gap-4 text-xs text-muted-foreground mt-1">
                <span className="font-mono">
                  {experiment?.molecule?.smiles || experimentConfig?.molecule?.smiles || "Custom Molecule"}
                </span>
                <span>•</span>
                <span>{experiment?.method || experimentConfig?.method || "VQE"}</span>
                <span>•</span>
                <span>{experiment?.backend || experimentConfig?.backend || "Classical"}</span>
              </div>
            )}
          </div>
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

              {/* VQE-specific fields */}
              {experimentConfig?.backendSettings?.method === "VQE" && (
                <div className="flex justify-between">
                  <span className="text-muted-foreground">Ansatz:</span>
                  <span className="font-quando font-medium">
                    {experimentConfig?.backendSettings?.ansatz || "HEA"}
                  </span>
                </div>
              )}

              {/* SQD-specific fields */}
              {experimentConfig?.backendSettings?.method === "SQD" && (
                <>
                  <div className="flex justify-between">
                    <span className="text-muted-foreground">Subspace:</span>
                    <span className="font-quando font-medium">
                      {experimentConfig?.backendSettings?.subspaceDim || 10}
                    </span>
                  </div>
                  <div className="flex justify-between">
                    <span className="text-muted-foreground">States:</span>
                    <span className="font-quando font-medium">
                      {experimentConfig?.backendSettings?.nStates || 3}
                    </span>
                  </div>
                </>
              )}

              {/* Excited States-specific fields */}
              {experimentConfig?.backendSettings?.method === "EXCITED_STATES" && (
                <>
                  <div className="flex justify-between">
                    <span className="text-muted-foreground">ES Method:</span>
                    <span className="font-quando font-medium">
                      {(experimentConfig?.backendSettings?.excited_method || experimentConfig?.backendSettings?.excitedMethod)?.toUpperCase() || "CIS"}
                    </span>
                  </div>
                  <div className="flex justify-between">
                    <span className="text-muted-foreground">States:</span>
                    <span className="font-quando font-medium">
                      {experimentConfig?.backendSettings?.nStates || experimentConfig?.backendSettings?.excited_n_states || experimentConfig?.backendSettings?.excitedNStates || 5}
                    </span>
                  </div>
                </>
              )}

              <div className="flex justify-between">
                <span className="text-muted-foreground">Basis:</span>
                <span className="font-quanto font-medium">
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
                {(() => {
                  const method = experimentConfig?.backendSettings?.method;
                  let maxIterations;

                  if (method === "VQE") {
                    maxIterations = experimentConfig?.backendSettings?.maxIterations || 100;
                  } else if (method === "SQD") {
                    maxIterations = 6; // SQD has 7 stages (0-6)
                  } else if (method === "EXCITED_STATES") {
                    // If using quantum backend, EXCITED_STATES redirects to SQD (7 stages)
                    const backend = experimentConfig?.backendSettings?.backend;
                    if (backend === "bluequbit" || backend === "ibm_quantum") {
                      maxIterations = 6; // Redirected to SQD
                    } else {
                      maxIterations = experimentConfig?.backendSettings?.nStates || 5;
                    }
                  } else {
                    maxIterations = 100;
                  }

                  return `Iteration ${currentIteration} / ${maxIterations}`;
                })()}
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
            {(() => {
              const method = experimentConfig?.backendSettings?.method;
              if (method === "SQD") return "Energy Spectrum";
              if (method === "EXCITED_STATES") return "Excited States Spectrum";
              return "Energy Convergence";
            })()}
          </h3>
          <div className="flex-1 min-h-0">
            <ResponsiveContainer width="100%" height="100%">
              {(() => {
                const method = experimentConfig?.backendSettings?.method;

                // For SQD and Excited States, use BarChart to show energy levels
                if (method === "SQD" || method === "EXCITED_STATES") {
                  // For SQD: filter to show only eigenvalues (stages 4+)
                  // For Excited States: show all data (already contains only states)
                  let displayData = convergenceData;
                  let stateOffset = 0;

                  if (method === "SQD" && convergenceData.length > 4) {
                    // SQD has 7 stages (0-6), we only want eigenvalues (4-6)
                    displayData = convergenceData.filter(d => d.iteration >= 4);
                    stateOffset = 4;
                  } else if (method === "EXCITED_STATES") {
                    // Excited States already sends only state data (iter 1, 2, 3...)
                    stateOffset = 1; // States are labeled starting from 1
                  }

                  return (
                    <BarChart data={displayData}>
                      <CartesianGrid strokeDasharray="3 3" stroke="hsl(var(--border))" />
                      <XAxis
                        dataKey="iteration"
                        label={{
                          value: "Energy State",
                          position: "insideBottom",
                          offset: -5
                        }}
                        tick={{ fontSize: 11 }}
                        stroke="hsl(var(--muted-foreground))"
                        tickFormatter={(value) => `S${value - stateOffset}`}
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
                        labelFormatter={(label) => `State ${label - stateOffset}`}
                      />
                      <Bar
                        dataKey="energy"
                        fill="#ea580c"
                        radius={[4, 4, 0, 0]}
                        animationDuration={300}
                      />
                    </BarChart>
                  );
                }

                // For VQE, use traditional LineChart
                return (
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
                );
              })()}
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
