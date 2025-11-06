"use client";

import { useState, useEffect } from "react";
import { ArrowLeft, CheckCircle2, XCircle, Loader2, Clock, Download } from "lucide-react";
import Link from "next/link";
import * as api from "@/lib/api";
import { useParams } from "next/navigation";
import CampaignAnimation from "@/components/campaign/CampaignAnimation";
import CampaignResults from "@/components/campaign/CampaignResults";

export default function CampaignMonitorPage() {
  const params = useParams();
  const campaignId = params.id as string;

  const [campaign, setCampaign] = useState<any>(null);
  const [experiments, setExperiments] = useState<any[]>([]);
  const [progress, setProgress] = useState<any>({
    total: 0,
    completed: 0,
    running: 0,
    pending: 0,
    failed: 0,
    percentage: 0
  });
  const [loading, setLoading] = useState(true);
  const [logs, setLogs] = useState<string[]>(["Initializing campaign..."]);
  const [elapsedTime, setElapsedTime] = useState(0);
  const [startTime] = useState(Date.now());

  useEffect(() => {
    if (!campaignId) return;

    // Timer for elapsed time
    const timer = setInterval(() => {
      if (campaign?.status !== "completed" && campaign?.status !== "failed") {
        setElapsedTime(Math.floor((Date.now() - startTime) / 1000));
      }
    }, 1000);

    // Start polling for campaign progress
    const cleanup = api.pollCampaignProgress(campaignId, (data) => {
      setCampaign(data.campaign);
      setExperiments(data.experiments || []);
      setProgress({
        total: data.progress?.total || 0,
        completed: data.progress?.completed || 0,
        running: data.progress?.running || 0,
        pending: data.progress?.pending || 0,
        failed: data.progress?.failed || 0,
        percentage: data.progress?.percentage || 0
      });
      setLoading(false);

      // Add log for status changes
      if (data.campaign?.status === "completed" && campaign?.status !== "completed") {
        addLog("Campaign completed successfully!");
      }
    });

    return () => {
      clearInterval(timer);
      cleanup();
    };
  }, [campaignId, campaign?.status, startTime]);

  const addLog = (message: string) => {
    setLogs((prev) => [...prev, message]);
  };

  const formatTime = (seconds: number) => {
    const mins = Math.floor(seconds / 60);
    const secs = seconds % 60;
    return `${mins}:${secs.toString().padStart(2, "0")}`;
  };

  const getStatusIcon = () => {
    switch (campaign?.status) {
      case "completed":
        return <CheckCircle2 className="w-6 h-6 text-green-500" />;
      case "running":
        return <Loader2 className="w-6 h-6 text-blue-500 animate-spin" />;
      case "failed":
        return <XCircle className="w-6 h-6 text-red-500" />;
      default:
        return <Clock className="w-6 h-6 text-yellow-500" />;
    }
  };

  const getStatusText = () => {
    switch (campaign?.status) {
      case "completed":
        return "Completed";
      case "running":
        return "Running";
      case "failed":
        return "Failed";
      default:
        return "Queued";
    }
  };

  // Determine campaign state
  const isInProgress = campaign?.status === "queued" || campaign?.status === "running";
  const isCompleted = campaign?.status === "completed";

  if (loading) {
    return (
      <div className="min-h-screen bg-background flex items-center justify-center">
        <Loader2 className="w-8 h-8 animate-spin" />
      </div>
    );
  }

  return (
    <div className="h-screen flex flex-col bg-background">
      {/* Header */}
      <div className="flex items-center justify-between px-6 py-4 border-b border-border">
        <div className="flex items-center gap-4">
          <Link href="/dashboard/queue">
            <button className="p-2 hover:bg-accent rounded-md transition">
              <ArrowLeft className="w-5 h-5" />
            </button>
          </Link>
          <div>
            <h1 className="text-xl font-quando font-bold">
              {campaign?.name || "Campaign Monitor"}
            </h1>
            <p className="text-xs text-muted-foreground mt-1">
              {campaign?.description || "Sequential experiment execution"}
            </p>
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
            {isCompleted && (
              <>
                <button className="flex items-center gap-2 px-4 py-2 text-sm border border-border rounded-lg hover:bg-accent transition font-quando">
                  <Download className="w-4 h-4" />
                  Export All
                </button>
                <Link href="/dashboard/queue">
                  <button className="px-4 py-2 text-sm border border-border rounded-lg hover:bg-accent transition font-quando">
                    Back to Queue
                  </button>
                </Link>
                <Link href="/dashboard/history">
                  <button className="px-4 py-2 text-sm bg-brand-orange text-white rounded-lg hover:bg-brand-orange-dark transition font-quando">
                    View All Results
                  </button>
                </Link>
              </>
            )}
          </div>
        </div>
      </div>

      {/* Main Content - Conditional Rendering */}
      <div className="flex-1 overflow-hidden">
        {isInProgress ? (
          // DURING EXECUTION: Show animation, progress bar, and logs
          <div className="h-full grid grid-cols-3 gap-6 p-6">
            {/* Left & Center: Animation + Progress */}
            <div className="col-span-2 flex flex-col gap-6">
              {/* Animation */}
              <div className="bg-card border border-border rounded-lg p-8 flex-1 flex items-center justify-center">
                <CampaignAnimation status={campaign?.status || "queued"} progress={progress} />
              </div>

              {/* Progress Stats */}
              <div className="grid grid-cols-4 gap-4">
                <div className="bg-card border border-border rounded-lg p-4">
                  <div className="text-xs text-muted-foreground font-quando mb-1">
                    Total
                  </div>
                  <div className="text-2xl font-quando font-bold">
                    {progress.total}
                  </div>
                </div>
                <div className="bg-card border border-border rounded-lg p-4">
                  <div className="text-xs text-muted-foreground font-quando mb-1">
                    Completed
                  </div>
                  <div className="text-2xl font-quando font-bold text-green-600">
                    {progress.completed}
                  </div>
                </div>
                <div className="bg-card border border-border rounded-lg p-4">
                  <div className="text-xs text-muted-foreground font-quando mb-1">
                    Running
                  </div>
                  <div className="text-2xl font-quando font-bold text-blue-600">
                    {progress.running}
                  </div>
                </div>
                <div className="bg-card border border-border rounded-lg p-4">
                  <div className="text-xs text-muted-foreground font-quando mb-1">
                    Pending
                  </div>
                  <div className="text-2xl font-quando font-bold text-gray-600">
                    {progress.pending}
                  </div>
                </div>
              </div>

              {/* Progress Bar */}
              <div className="bg-card border border-border rounded-lg p-6">
                <div className="flex justify-between items-center mb-3">
                  <h3 className="text-lg font-quando font-bold">Campaign Progress</h3>
                  <span className="text-sm font-quando font-semibold text-muted-foreground">
                    {progress.percentage.toFixed(0)}%
                  </span>
                </div>
                <div className="w-full bg-muted rounded-full h-4 overflow-hidden">
                  <div
                    className="h-full bg-gradient-to-r from-brand-orange to-brand-yellow transition-all duration-300 ease-out"
                    style={{ width: `${progress.percentage}%` }}
                  />
                </div>
                <div className="mt-4 text-sm text-muted-foreground font-quando">
                  {progress.completed} of {progress.total} experiments completed
                </div>
              </div>
            </div>

            {/* Right: Logs */}
            <div className="bg-card border border-border rounded-lg p-6 flex flex-col min-h-0">
              <h3 className="text-lg font-quando font-bold mb-4">
                Execution Logs
              </h3>
              <div className="flex-1 bg-muted rounded-lg p-4 font-mono text-xs overflow-auto min-h-0">
                <div className="space-y-1">
                  {logs.map((log, i) => (
                    <div key={i} className="text-foreground">
                      <span className="text-muted-foreground mr-2">
                        [{formatTime(Math.floor((i + 1) * 2))}]
                      </span>
                      {log}
                    </div>
                  ))}
                  {experiments.map((exp, i) => (
                    <div key={exp.id} className="text-foreground">
                      <span className="text-muted-foreground mr-2">
                        [{formatTime(Math.floor((i + logs.length + 1) * 2))}]
                      </span>
                      Experiment {i + 1}: {exp.status}
                      {exp.status === "completed" && exp.results?.energy &&
                        ` - E = ${exp.results.energy.toFixed(6)} Ha`
                      }
                    </div>
                  ))}
                </div>
              </div>
            </div>
          </div>
        ) : isCompleted ? (
          // AFTER COMPLETION: Show beautiful results with comparison
          <CampaignResults campaign={campaign} experiments={experiments} />
        ) : (
          // FAILED: Show error message
          <div className="h-full flex items-center justify-center">
            <div className="text-center">
              <div className="mb-4">
                <XCircle className="w-16 h-16 text-red-500 mx-auto" />
              </div>
              <h2 className="text-2xl font-quando font-bold mb-2">
                Campaign Failed
              </h2>
              <p className="text-muted-foreground mb-6">
                The campaign encountered an error during execution.
              </p>
              <Link href="/dashboard/queue">
                <button className="px-6 py-2.5 bg-brand-orange text-white rounded-lg hover:bg-brand-orange-dark transition font-quando">
                  Back to Queue
                </button>
              </Link>
            </div>
          </div>
        )}
      </div>

      {/* Footer Actions - Only show when completed */}
      {isCompleted && (
        <div className="border-t border-border px-6 py-4">
          <div className="flex gap-4 justify-end">
            <Link href="/dashboard/queue">
              <button className="px-6 py-2.5 text-sm border border-border rounded-lg hover:bg-accent transition font-quando">
                Back to Queue
              </button>
            </Link>
            <Link href="/dashboard/history">
              <button className="px-6 py-2.5 text-sm bg-brand-orange text-white rounded-lg hover:bg-brand-orange-dark transition font-quando">
                View All Results
              </button>
            </Link>
          </div>
        </div>
      )}
    </div>
  );
}
