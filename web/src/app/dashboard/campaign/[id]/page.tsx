"use client";

import { useState, useEffect } from "react";
import { ArrowLeft, CheckCircle2, XCircle, Loader2, Clock } from "lucide-react";
import Link from "next/link";
import * as api from "@/lib/api";
import { useParams } from "next/navigation";

export default function CampaignMonitorPage() {
  const params = useParams();
  const campaignId = params.id as string;

  const [campaign, setCampaign] = useState<any>(null);
  const [experiments, setExperiments] = useState<any[]>([]);
  const [progress, setProgress] = useState<any>(null);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    if (!campaignId) return;

    // Start polling for campaign progress
    const cleanup = api.pollCampaignProgress(campaignId, (data) => {
      setCampaign(data.campaign);
      setExperiments(data.experiments || []);
      setProgress(data.progress);
      setLoading(false);
    });

    return cleanup;
  }, [campaignId]);

  const getStatusColor = (status: string) => {
    switch (status) {
      case "completed":
        return "text-green-600 bg-green-100 dark:bg-green-900/20";
      case "running":
        return "text-blue-600 bg-blue-100 dark:bg-blue-900/20";
      case "failed":
        return "text-red-600 bg-red-100 dark:bg-red-900/20";
      case "pending":
        return "text-gray-600 bg-gray-100 dark:bg-gray-900/20";
      default:
        return "text-gray-600 bg-gray-100 dark:bg-gray-900/20";
    }
  };

  const getStatusIcon = (status: string) => {
    switch (status) {
      case "completed":
        return <CheckCircle2 className="w-5 h-5" />;
      case "running":
        return <Loader2 className="w-5 h-5 animate-spin" />;
      case "failed":
        return <XCircle className="w-5 h-5" />;
      default:
        return <Clock className="w-5 h-5" />;
    }
  };

  if (loading) {
    return (
      <div className="min-h-screen bg-background flex items-center justify-center">
        <Loader2 className="w-8 h-8 animate-spin" />
      </div>
    );
  }

  return (
    <div className="min-h-screen bg-background p-8">
      <div className="max-w-6xl mx-auto">
        {/* Header */}
        <div className="flex items-center gap-4 mb-8">
          <Link href="/dashboard/queue">
            <button className="p-2 hover:bg-accent rounded-lg transition">
              <ArrowLeft className="w-6 h-6" />
            </button>
          </Link>
          <div className="flex-1">
            <h1 className="text-3xl font-quando font-bold">
              {campaign?.name || "Campaign"}
            </h1>
            <p className="text-muted-foreground mt-1">
              {campaign?.description || "Sequential experiment execution"}
            </p>
          </div>
          <div className={`px-4 py-2 rounded-lg ${getStatusColor(campaign?.status || "")}`}>
            <span className="font-quando font-semibold capitalize">
              {campaign?.status || "Unknown"}
            </span>
          </div>
        </div>

        {/* Progress Overview */}
        <div className="grid grid-cols-4 gap-4 mb-8">
          <div className="bg-card border border-border rounded-lg p-6">
            <p className="text-sm text-muted-foreground font-quando mb-1">
              Total Experiments
            </p>
            <p className="text-3xl font-quando font-bold">
              {progress?.total || 0}
            </p>
          </div>
          <div className="bg-card border border-border rounded-lg p-6">
            <p className="text-sm text-muted-foreground font-quando mb-1">
              Completed
            </p>
            <p className="text-3xl font-quando font-bold text-green-600">
              {progress?.completed || 0}
            </p>
          </div>
          <div className="bg-card border border-border rounded-lg p-6">
            <p className="text-sm text-muted-foreground font-quando mb-1">
              Failed
            </p>
            <p className="text-3xl font-quando font-bold text-red-600">
              {progress?.failed || 0}
            </p>
          </div>
          <div className="bg-card border border-border rounded-lg p-6">
            <p className="text-sm text-muted-foreground font-quando mb-1">
              Progress
            </p>
            <p className="text-3xl font-quando font-bold">
              {Math.round(progress?.percentage || 0)}%
            </p>
          </div>
        </div>

        {/* Progress Bar */}
        <div className="bg-card border border-border rounded-lg p-6 mb-8">
          <div className="flex items-center justify-between mb-3">
            <h3 className="text-lg font-quando font-semibold">
              Overall Progress
            </h3>
            <span className="text-sm text-muted-foreground">
              {progress?.completed + progress?.failed || 0} of {progress?.total || 0}
            </span>
          </div>
          <div className="w-full bg-gray-200 dark:bg-gray-800 rounded-full h-3">
            <div
              className="bg-brand-orange h-3 rounded-full transition-all duration-500"
              style={{ width: `${progress?.percentage || 0}%` }}
            />
          </div>
        </div>

        {/* Experiment List */}
        <div className="space-y-3">
          <h2 className="text-xl font-quando font-semibold mb-4">
            Experiment Sequence
          </h2>

          {experiments.map((exp, index) => (
            <div
              key={exp.id}
              className="bg-card border border-border rounded-lg p-4 flex items-center gap-4"
            >
              <div className="w-8 h-8 rounded-full bg-brand-orange/10 text-brand-orange flex items-center justify-center font-semibold flex-shrink-0">
                {index + 1}
              </div>

              <div className="flex-1">
                <h3 className="font-quando font-semibold">
                  Experiment {exp.id.slice(0, 8)}
                </h3>
                <div className="flex gap-4 text-sm text-muted-foreground mt-1">
                  <span>Method: {exp.method}</span>
                  <span>Backend: {exp.backend}</span>
                </div>
              </div>

              <div
                className={`flex items-center gap-2 px-3 py-1.5 rounded-lg ${getStatusColor(
                  exp.status
                )}`}
              >
                {getStatusIcon(exp.status)}
                <span className="font-quando text-sm capitalize">
                  {exp.status}
                </span>
              </div>
            </div>
          ))}
        </div>

        {/* Complete Action */}
        {campaign?.status === "completed" && (
          <div className="bg-card border border-green-600 rounded-lg p-6 mt-8">
            <div className="flex items-center gap-4">
              <CheckCircle2 className="w-8 h-8 text-green-600" />
              <div className="flex-1">
                <h3 className="text-lg font-quando font-semibold text-green-600">
                  Campaign Completed!
                </h3>
                <p className="text-sm text-muted-foreground mt-1">
                  All experiments have finished executing
                </p>
              </div>
              <Link href="/dashboard/history">
                <button className="px-6 py-3 bg-brand-orange text-white rounded-lg hover:bg-brand-orange-dark transition font-quando">
                  View Results
                </button>
              </Link>
            </div>
          </div>
        )}

        {/* Failed Status */}
        {campaign?.status === "failed" && (
          <div className="bg-card border border-red-600 rounded-lg p-6 mt-8">
            <div className="flex items-center gap-4">
              <XCircle className="w-8 h-8 text-red-600" />
              <div className="flex-1">
                <h3 className="text-lg font-quando font-semibold text-red-600">
                  Campaign Failed
                </h3>
                <p className="text-sm text-muted-foreground mt-1">
                  Some experiments encountered errors
                </p>
              </div>
            </div>
          </div>
        )}
      </div>
    </div>
  );
}
