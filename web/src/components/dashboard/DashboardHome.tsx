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
}

export default function DashboardHome({
  onNewExperiment,
  onViewExperiment,
  experiments = [],
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

  const recentExperiments = experiments.slice(0, 3);
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

      // Refresh data
      window.location.reload();
    } catch (error: any) {
      console.error("Failed to cancel experiment:", error);
      toast.error(error.message || "Failed to cancel experiment");
    } finally {
      setIsCancelling(false);
    }
  };

  return (
    <div className="h-full flex flex-col bg-background p-6 space-y-6">
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
      <div>
        <h1 className="text-3xl font-quando font-bold mb-2">Dashboard</h1>
        <p className="text-muted-foreground font-quando">
          Welcome back! Here&apos;s an overview of your quantum experiments.
        </p>
      </div>

     

      {/* Quick Stats */}
      <div className="grid grid-cols-4 gap-4">
        <div className="bg-card border border-border rounded-lg p-6">
          <div className="flex items-center justify-between">
            <div>
              <p className="text-sm text-muted-foreground font-quando mb-1">
                Total Experiments
              </p>
              <p className="text-3xl font-quando font-bold">{totalCount}</p>
            </div>
            <TrendingUp className="w-8 h-8 " />
          </div>
        </div>

        <div className="bg-card border border-border rounded-lg p-6">
          <div className="flex items-center justify-between">
            <div>
              <p className="text-sm text-muted-foreground font-quando mb-1">
                Completed
              </p>
              <p className="text-3xl font-quando font-bold ">
                {completedCount}
              </p>
            </div>
            <CheckCircle2 className="w-8 h-8 " />
          </div>
        </div>

        <div className="bg-card border border-border rounded-lg p-6">
          <div className="flex items-center justify-between">
            <div>
              <p className="text-sm text-muted-foreground font-quando mb-1">
                Running
              </p>
              <p className="text-3xl font-quando font-bold ">
                {runningCount}
              </p>
            </div>
            <Play className="w-8 h-8 " />
          </div>
        </div>

        <div className="bg-card border border-border rounded-lg p-6">
          <div className="flex items-center justify-between">
            <div>
              <p className="text-sm text-muted-foreground font-quando mb-1">
                In Queue
              </p>
              <p className="text-3xl font-quando font-bold ">
                {queuedCount}
              </p>
            </div>
            <Clock className="w-8 h-8 " />
          </div>
        </div>
      </div>

      {/* Quick Actions */}
      <div className="grid grid-cols-2 gap-4">
        <button
          onClick={onNewExperiment}
          className="bg-brand-orange text-white rounded-lg p-8 hover:bg-brand-orange-dark transition flex flex-col items-center justify-center gap-3"
        >
          <Play className="w-12 h-12" />
          <div>
            <p className="text-xl font-quando font-bold">New Experiment</p>
            <p className="text-sm opacity-90">Create a new quantum calculation</p>
          </div>
        </button>

        <a
          href="/dashboard/queue"
          className="bg-card border-2 border-border rounded-lg p-8 hover:border-brand-orange transition flex flex-col items-center justify-center gap-3"
        >
          <Clock className="w-12 h-12" />
          <div>
            <p className="text-xl font-quando font-bold">Job Queue</p>
            <p className="text-sm text-muted-foreground">
              Manage and schedule experiments
            </p>
          </div>
        </a>
      </div>

      {/* Running Experiments */}
      {runningExperiments.length > 0 && (
        <div>
          <h2 className="text-xl font-quando font-bold mb-4">Running Experiments</h2>
          <div className="space-y-3">
            {runningExperiments.map((exp) => (
              <div
                key={exp.id}
                onClick={() => onViewExperiment && onViewExperiment(exp.id, exp)}
                className="bg-card border-2 border-blue-500 rounded-lg p-4 hover:border-brand-orange transition cursor-pointer"
              >
                <div className="flex items-center justify-between">
                  <div className="flex-1">
                    <div className="flex items-center gap-3 mb-2">
                      <Loader2 className="w-5 h-5 text-blue-500 animate-spin" />
                      <h3 className="font-quando font-semibold">
                        {exp.molecule?.smiles || exp.name || "Custom Molecule"}
                      </h3>
                      <span className="px-2 py-1 text-xs rounded font-quando bg-blue-100 text-blue-800 dark:bg-blue-900 dark:text-blue-200">
                        Running
                      </span>
                    </div>
                    <div className="flex gap-6 text-xs text-muted-foreground font-quando">
                      <span>Method: {exp.method || "VQE"}</span>
                      <span>Backend: {exp.backend || "Classical"}</span>
                      {exp.startedAt && (
                        <span>
                          Started: {new Date(exp.startedAt).toLocaleTimeString()}
                        </span>
                      )}
                    </div>
                  </div>
                  <div className="flex items-center gap-3">
                    <button
                      onClick={(e) => handleCancelClick(exp, e)}
                      className="p-2 border-2 border-red-600 text-red-600 rounded-lg hover:bg-red-600 hover:text-white transition"
                      title="Cancel experiment"
                    >
                      <XCircle className="w-5 h-5" />
                    </button>
                    <div className="text-right">
                      <p className="text-sm text-blue-600 font-quando font-medium">
                        Click to monitor
                      </p>
                    </div>
                  </div>
                </div>
              </div>
            ))}
          </div>
        </div>
      )}

      {/* Recent Experiments */}
      <div className="flex-1 min-h-0">
        <div className="flex items-center justify-between mb-4">
          <h2 className="text-xl font-quando font-bold">Recent Experiments</h2>
          <a
            href="/dashboard/history"
            className="text-sm text-brand-orange hover:underline font-quando"
          >
            View All
          </a>
        </div>

        {experiments.length === 0 ? (
          <div className="bg-card border border-border rounded-lg p-12 flex flex-col items-center justify-center text-center">
            <History className="w-16 h-16 text-muted-foreground mb-4" />
            <p className="text-lg font-quando font-semibold mb-2">
              No experiments yet
            </p>
            <p className="text-sm text-muted-foreground mb-6">
              Start your first quantum chemistry experiment
            </p>
            <button
              onClick={onNewExperiment}
              className="px-6 py-3 bg-brand-orange text-white rounded-lg hover:bg-brand-orange-dark transition font-quando"
            >
              Create Experiment
            </button>
          </div>
        ) : (
          <div className="space-y-3 overflow-auto max-h-64">
            {recentExperiments.map((exp) => (
              <div
                key={exp.id}
                onClick={() => exp.status === "completed" && onViewExperiment && onViewExperiment(exp.id, exp)}
                className={`bg-card border border-border rounded-lg p-4 hover:border-brand-orange transition ${
                  exp.status === "completed" ? "cursor-pointer" : ""
                }`}
              >
                <div className="flex items-center justify-between">
                  <div className="flex-1">
                    <div className="flex items-center gap-3 mb-2">
                      <h3 className="font-quando font-semibold">
                        {exp.molecule?.smiles || exp.name || "Custom Molecule"}
                      </h3>
                      <span
                        className={`px-2 py-1 text-xs rounded font-quando ${
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
                    <div className="flex gap-6 text-xs text-muted-foreground font-quando">
                      <span>Method: {exp.method || "VQE"}</span>
                      <span>Backend: {exp.backend || "Classical"}</span>
                      <span>
                        {new Date(exp.timestamp).toLocaleDateString()}
                      </span>
                    </div>
                  </div>
                  {exp.results && exp.results.energy && (
                    <div className="text-right">
                      <p className="text-sm text-muted-foreground">Energy</p>
                      <p className="font-mono font-bold text-lg">
                        {exp.results.energy.toFixed(6)} Ha
                      </p>
                    </div>
                  )}
                  {exp.energy && !exp.results && (
                    <div className="text-right">
                      <p className="text-sm text-muted-foreground">Energy</p>
                      <p className="font-mono font-bold text-lg">
                        {exp.energy.toFixed(6)} Ha
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
  );
}
