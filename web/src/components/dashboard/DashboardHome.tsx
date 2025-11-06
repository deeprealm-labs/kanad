"use client";

import { useState, useEffect } from "react";
import { Play, History, Clock, CheckCircle2, TrendingUp, Loader2, XCircle } from "lucide-react";
import * as api from "@/lib/api";
import { useToast } from "@/components/ui/toast";
import CancelConfirmationDialog from "@/components/experiment/CancelConfirmationDialog";

interface DashboardHomeProps {
  onNewExperiment: () => void;
  onViewExperiment?: (experimentId: string, config: any) => void;
  experiments: any[];
  onRefresh?: () => void;
}

export default function DashboardHome({
  onNewExperiment,
  onViewExperiment,
  experiments = [],
  onRefresh,
}: DashboardHomeProps) {
  const [experimentStats, setExperimentStats] = useState<any>(null);
  const [queueStats, setQueueStats] = useState<any>(null);
  const [loadingStats, setLoadingStats] = useState(false);
  const [showCancelDialog, setShowCancelDialog] = useState(false);
  const [experimentToCancel, setExperimentToCancel] = useState<any>(null);
  const [isCancelling, setIsCancelling] = useState(false);
  const toast = useToast();

  useEffect(() => {
    loadAllStats();
  }, []);

  const loadAllStats = async () => {
    try {
      setLoadingStats(true);
      const [expStats, qStats] = await Promise.all([
        api.getExperimentStatistics(),
        api.getQueueStatistics()
      ]);
      setExperimentStats(expStats);
      setQueueStats(qStats);
    } catch (error) {
      console.error("Failed to load stats:", error);
    } finally {
      setLoadingStats(false);
    }
  };

  // Use stats from API if available, otherwise fallback to local experiments array
  const totalCount = experimentStats?.total ?? experiments.length;
  const completedCount = experimentStats?.completed ?? experiments.filter((e) => e.status === "completed").length;
  const runningCount = experimentStats?.running ?? experiments.filter((e) => e.status === "running").length;
  const queuedCount = experimentStats?.queued ?? queueStats?.queued ?? experiments.filter((e) => e.status === "queued").length;

  const recentExperiments = experiments.slice(0, 10);
  const runningExperiments = experiments.filter((e) => e.status === "running");

  const handleCancelClick = (exp: any, event: React.MouseEvent) => {
    event.stopPropagation();
    setExperimentToCancel(exp);
    setShowCancelDialog(true);
  };

  const handleCancelConfirm = async () => {
    if (!experimentToCancel) return;

    try {
      setIsCancelling(true);
      await api.cancelExperiment(experimentToCancel.id);
      setShowCancelDialog(false);
      setExperimentToCancel(null);
      toast.success("Experiment cancelled successfully");

      // Refresh data via callback (better than page reload)
      if (onRefresh) {
        onRefresh();
      }
      await loadAllStats();
    } catch (error: any) {
      console.error("Failed to cancel experiment:", error);
      toast.error(error.message || "Failed to cancel experiment");
    } finally {
      setIsCancelling(false);
    }
  };

  return (
    <div className="h-full flex flex-col bg-background overflow-hidden">
      {/* Cancel Confirmation Dialog */}
      {experimentToCancel && (
        <CancelConfirmationDialog
          isOpen={showCancelDialog}
          onClose={() => {
            setShowCancelDialog(false);
            setExperimentToCancel(null);
          }}
          onConfirm={handleCancelConfirm}
          experimentName={experimentToCancel.molecule?.smiles || experimentToCancel.name || "Custom Experiment"}
          isRunning={experimentToCancel.status === "running"}
          isLoading={isCancelling}
        />
      )}

      {/* Header */}
      <div className="px-8 pt-4 pb-3 border-b border-border flex-shrink-0">
        <div>
          <h1 className="text-2xl font-quando font-bold tracking-tight">Dashboard</h1>
          <p className="text-muted-foreground font-quando text-xs mt-1">
            Welcome back! Here&apos;s an overview of your quantum experiments.
          </p>
        </div>
      </div>

      {/* Main Content Area - Two Column Layout */}
      <div className="flex-1 flex gap-6 p-6 overflow-hidden">

        {/* LEFT COLUMN: Stats and Quick Actions */}
        <div className="w-80 flex-shrink-0 flex flex-col gap-4 overflow-y-auto">
          {/* Stats Cards */}
          <div className="space-y-3">
            <h2 className="text-lg font-quando font-bold mb-3">Statistics</h2>

            <div className="bg-card border border-border rounded-lg p-4">
              <div className="flex items-center justify-between">
                <div>
                  <p className="text-xs text-muted-foreground font-quando uppercase mb-1">Total</p>
                  <p className="text-2xl font-quando font-bold">{totalCount}</p>
                </div>
                <TrendingUp className="w-6 h-6 text-muted-foreground opacity-50" />
              </div>
            </div>

            <div className="bg-card border border-border rounded-lg p-4">
              <div className="flex items-center justify-between">
                <div>
                  <p className="text-xs text-muted-foreground font-quando uppercase mb-1">Completed</p>
                  <p className="text-2xl font-quando font-bold">{completedCount}</p>
                </div>
                <CheckCircle2 className="w-6 h-6 text-muted-foreground opacity-50" />
              </div>
            </div>

            <div className="bg-card border border-border rounded-lg p-4">
              <div className="flex items-center justify-between">
                <div>
                  <p className="text-xs text-muted-foreground font-quando uppercase mb-1">Running</p>
                  <p className="text-2xl font-quando font-bold">{runningCount}</p>
                </div>
                <Play className="w-6 h-6 text-muted-foreground opacity-50" />
              </div>
            </div>

            <div className="bg-card border border-border rounded-lg p-4">
              <div className="flex items-center justify-between">
                <div>
                  <p className="text-xs text-muted-foreground font-quando uppercase mb-1">In Queue</p>
                  <p className="text-2xl font-quando font-bold">{queuedCount}</p>
                </div>
                <Clock className="w-6 h-6 text-muted-foreground opacity-50" />
              </div>
            </div>
          </div>

          {/* Quick Action Buttons */}
          <div className="space-y-3">
            <h2 className="text-lg font-quando font-bold mb-3">Quick Actions</h2>

            <button
              onClick={onNewExperiment}
              className="w-full bg-brand-orange text-white rounded-lg p-4 hover:bg-brand-orange-dark transition-all flex flex-col items-center justify-center gap-2 group"
            >
              <Play className="w-8 h-8 group-hover:scale-110 transition-transform" />
              <div className="text-center">
                <p className="text-base font-quando font-bold">New Experiment</p>
                <p className="text-xs opacity-90">Create calculation</p>
              </div>
            </button>

            <a
              href="/dashboard/queue"
              className="w-full bg-card border border-border rounded-lg p-4 hover:border-brand-orange transition-all flex flex-col items-center justify-center gap-2 group block"
            >
              <Clock className="w-8 h-8 group-hover:scale-110 transition-transform" />
              <div className="text-center">
                <p className="text-base font-quando font-bold">Job Queue</p>
                <p className="text-xs text-muted-foreground">Manage queue</p>
              </div>
            </a>
          </div>
        </div>

        {/* RIGHT COLUMN: Recent Experiments */}
        <div className="flex-1 flex flex-col overflow-hidden">
          <div className="flex items-center justify-between mb-4">
            <h2 className="text-2xl font-quando font-bold">Recent Experiments</h2>
            <a
              href="/dashboard/history"
              className="text-sm text-brand-orange hover:underline font-quando font-medium"
            >
              View All →
            </a>
          </div>

          {experiments.length === 0 ? (
            <div className="flex-1 flex items-center justify-center">
              <div className="text-center">
                <History className="w-16 h-16 text-muted-foreground mb-4 mx-auto" />
                <p className="text-xl font-quando font-semibold mb-2">No experiments yet</p>
                <p className="text-sm text-muted-foreground mb-6">
                  Start your first quantum chemistry experiment
                </p>
                <button
                  onClick={onNewExperiment}
                  className="px-6 py-3 bg-brand-orange text-white rounded-lg hover:bg-brand-orange-dark transition-all font-quando font-semibold"
                >
                  Create Experiment
                </button>
              </div>
            </div>
          ) : (
            <div className="flex-1 space-y-2 overflow-y-auto pr-2">
              {recentExperiments.map((exp) => (
                <div
                  key={exp.id}
                  onClick={() => exp.status === "completed" && onViewExperiment && onViewExperiment(exp.id, exp)}
                  className={`bg-card border border-border rounded-lg p-3 hover:border-brand-orange transition-all hover:shadow-sm ${
                    exp.status === "completed" ? "cursor-pointer" : ""
                  }`}
                >
                  <div className="flex items-center justify-between gap-4">
                    <div className="flex-1 min-w-0">
                      <div className="flex items-center gap-2 mb-1">
                        <h3 className="font-quando font-semibold text-sm truncate">
                          {exp.molecule?.smiles || exp.name || "Custom Molecule"}
                        </h3>
                        <span
                          className={`px-2 py-0.5 text-xs rounded font-quando flex-shrink-0 ${
                            exp.status === "completed"
                              ? "bg-green-100 text-green-800 dark:bg-green-900 dark:text-green-200"
                              : exp.status === "running"
                              ? "bg-blue-100 text-blue-800 dark:bg-blue-900 dark:text-blue-200"
                              : exp.status === "cancelled"
                              ? "bg-orange-100 text-orange-800 dark:bg-orange-900 dark:text-orange-200"
                              : exp.status === "failed"
                              ? "bg-red-100 text-red-800 dark:bg-red-900 dark:text-red-200"
                              : "bg-yellow-100 text-yellow-800 dark:bg-yellow-900 dark:text-yellow-200"
                          }`}
                        >
                          {exp.status}
                        </span>
                      </div>
                      <div className="flex gap-4 text-xs text-muted-foreground font-quando">
                        <span>{exp.method || "VQE"}</span>
                        <span>{exp.backend || "Classical"}</span>
                        <span>{new Date(exp.timestamp).toLocaleDateString()}</span>
                      </div>
                    </div>
                    {(exp.results?.energy || exp.energy) && (
                      <div className="text-right flex-shrink-0">
                        <p className="text-xs text-muted-foreground">Energy</p>
                        <p className="font-mono font-bold text-sm">
                          {(exp.results?.energy || exp.energy).toFixed(6)} Ha
                        </p>
                      </div>
                    )}
                  </div>
                </div>
              ))}
            </div>
          )}
        </div>
      </div>
    </div>
  );
}
