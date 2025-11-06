"use client";

import { useState } from "react";
import { LineChart, Line, BarChart, Bar, CartesianGrid, XAxis, YAxis, Tooltip, ResponsiveContainer, Legend } from "recharts";
import { CheckCircle2, XCircle, TrendingDown, Activity } from "lucide-react";
import Link from "next/link";

interface CampaignResultsProps {
  campaign: any;
  experiments: any[];
}

export default function CampaignResults({ campaign, experiments }: CampaignResultsProps) {
  const [activeTab, setActiveTab] = useState<"overview" | "comparison">("overview");

  // Calculate statistics
  const completedExperiments = experiments.filter(e => e.status === "completed");
  const failedExperiments = experiments.filter(e => e.status === "failed");
  const successRate = (completedExperiments.length / experiments.length) * 100;

  // Extract energies for comparison
  const energyComparison = completedExperiments
    .map((exp, i) => ({
      index: i + 1,
      id: exp.id.slice(0, 8),
      energy: exp.results?.energy || 0,
      method: exp.method,
      converged: exp.results?.converged
    }))
    .filter(e => e.energy !== 0);

  return (
    <div className="h-full flex flex-col bg-background">
      {/* Tabs Navigation */}
      <div className="flex gap-4 px-6 py-4 border-b border-border">
        <button
          onClick={() => setActiveTab("overview")}
          className={`px-4 py-2 text-sm font-quando font-semibold rounded-lg transition ${
            activeTab === "overview"
              ? "bg-brand-orange text-white"
              : "text-muted-foreground hover:text-foreground hover:bg-accent"
          }`}
        >
          <div className="flex items-center gap-2">
            <Activity className="w-4 h-4" />
            Campaign Overview
          </div>
        </button>
        <button
          onClick={() => setActiveTab("comparison")}
          className={`px-4 py-2 text-sm font-quando font-semibold rounded-lg transition ${
            activeTab === "comparison"
              ? "bg-brand-orange text-white"
              : "text-muted-foreground hover:text-foreground hover:bg-accent"
          }`}
        >
          <div className="flex items-center gap-2">
            <TrendingDown className="w-4 h-4" />
            Energy Comparison
          </div>
        </button>
      </div>

      {/* Tab Content */}
      <div className="flex-1 overflow-y-auto p-6">
        {activeTab === "overview" ? (
          <OverviewTab
            campaign={campaign}
            experiments={experiments}
            completedExperiments={completedExperiments}
            failedExperiments={failedExperiments}
            successRate={successRate}
          />
        ) : (
          <ComparisonTab energyComparison={energyComparison} />
        )}
      </div>
    </div>
  );
}

// Overview Tab
function OverviewTab({
  campaign,
  experiments,
  completedExperiments,
  failedExperiments,
  successRate
}: {
  campaign: any;
  experiments: any[];
  completedExperiments: any[];
  failedExperiments: any[];
  successRate: number;
}) {
  return (
    <div className="grid grid-cols-2 gap-6">
      {/* Campaign Statistics */}
      <div className="bg-card border border-border rounded-lg p-6">
        <h3 className="text-lg font-quando font-bold mb-6">Campaign Statistics</h3>
        <div className="space-y-4">
          <div className="bg-gradient-to-r from-brand-orange/10 to-brand-orange/5 rounded-lg p-4 border border-brand-orange/20">
            <div className="text-sm text-muted-foreground mb-1">Total Experiments</div>
            <div className="text-3xl font-quando font-bold text-brand-orange">
              {experiments.length}
            </div>
          </div>

          <div className="grid grid-cols-2 gap-3">
            <div className="bg-gradient-to-br from-green-50 to-green-100 dark:from-green-900/30 dark:to-green-900/20 border border-green-200 dark:border-green-800 rounded-lg p-4">
              <div className="flex items-center gap-2 mb-2">
                <CheckCircle2 className="w-4 h-4 text-green-600" />
                <div className="text-sm text-green-700 dark:text-green-300">Completed</div>
              </div>
              <div className="text-2xl font-quanto font-bold text-green-900 dark:text-green-100">
                {completedExperiments.length}
              </div>
            </div>

            <div className="bg-gradient-to-br from-red-50 to-red-100 dark:from-red-900/30 dark:to-red-900/20 border border-red-200 dark:border-red-800 rounded-lg p-4">
              <div className="flex items-center gap-2 mb-2">
                <XCircle className="w-4 h-4 text-red-600" />
                <div className="text-sm text-red-700 dark:text-red-300">Failed</div>
              </div>
              <div className="text-2xl font-quanto font-bold text-red-900 dark:text-red-100">
                {failedExperiments.length}
              </div>
            </div>
          </div>

          <div className="bg-muted rounded-lg p-4">
            <div className="text-sm text-muted-foreground mb-1">Success Rate</div>
            <div className="text-2xl font-quanto font-bold">
              {successRate.toFixed(1)}%
            </div>
            <div className="w-full bg-background rounded-full h-2 mt-2">
              <div
                className="h-2 bg-green-500 rounded-full transition-all duration-300"
                style={{ width: `${successRate}%` }}
              />
            </div>
          </div>

          {campaign?.created_at && (
            <div className="bg-muted rounded-lg p-4">
              <div className="text-sm text-muted-foreground mb-1">Execution Time</div>
              <div className="text-xl font-quanto font-bold">
                {campaign.completed_at
                  ? `${Math.round((new Date(campaign.completed_at).getTime() - new Date(campaign.created_at).getTime()) / 1000 / 60)} min`
                  : "In Progress"}
              </div>
            </div>
          )}
        </div>
      </div>

      {/* Experiment List */}
      <div className="bg-card border border-border rounded-lg p-6">
        <h3 className="text-lg font-quando font-bold mb-6">Experiment Results</h3>
        <div className="space-y-3 max-h-[600px] overflow-y-auto">
          {experiments.map((exp, index) => (
            <div
              key={exp.id}
              className={`border rounded-lg p-4 ${
                exp.status === "completed"
                  ? "border-green-200 dark:border-green-800 bg-green-50 dark:bg-green-900/10"
                  : exp.status === "failed"
                  ? "border-red-200 dark:border-red-800 bg-red-50 dark:bg-red-900/10"
                  : "border-border bg-background"
              }`}
            >
              <div className="flex items-start justify-between gap-3">
                <div className="flex items-center gap-3">
                  <div className="w-8 h-8 rounded-full bg-brand-orange/20 text-brand-orange flex items-center justify-center font-semibold text-sm flex-shrink-0">
                    {index + 1}
                  </div>
                  <div>
                    <div className="font-quanto font-semibold text-sm">
                      Experiment {exp.id.slice(0, 8)}
                    </div>
                    <div className="text-xs text-muted-foreground mt-0.5">
                      {exp.method} • {exp.backend}
                    </div>
                  </div>
                </div>

                <div className="flex flex-col items-end gap-1">
                  {exp.status === "completed" && exp.results?.energy && (
                    <div className="text-xs">
                      <span className="text-muted-foreground">E:</span>{" "}
                      <span className="font-mono font-semibold">{exp.results.energy.toFixed(6)} Ha</span>
                    </div>
                  )}
                  <div
                    className={`flex items-center gap-1 px-2 py-0.5 rounded text-xs font-quando font-semibold ${
                      exp.status === "completed"
                        ? "text-green-700 bg-green-100 dark:text-green-300 dark:bg-green-900/30"
                        : exp.status === "failed"
                        ? "text-red-700 bg-red-100 dark:text-red-300 dark:bg-red-900/30"
                        : "text-gray-700 bg-gray-100 dark:text-gray-300 dark:bg-gray-900/30"
                    }`}
                  >
                    {exp.status === "completed" && <CheckCircle2 className="w-3 h-3" />}
                    {exp.status === "failed" && <XCircle className="w-3 h-3" />}
                    <span className="capitalize">{exp.status}</span>
                  </div>
                </div>
              </div>

              {exp.status === "completed" && (
                <Link href={`/dashboard/history/${exp.id}`}>
                  <button className="mt-3 text-xs text-brand-orange hover:text-brand-yellow font-quando font-semibold transition">
                    View Details →
                  </button>
                </Link>
              )}
            </div>
          ))}
        </div>
      </div>
    </div>
  );
}

// Comparison Tab
function ComparisonTab({ energyComparison }: { energyComparison: any[] }) {
  // Find lowest energy
  const lowestEnergy = energyComparison.length > 0
    ? Math.min(...energyComparison.map(e => e.energy))
    : 0;

  return (
    <div className="space-y-6">
      {/* Energy Bar Chart */}
      <div className="bg-card border border-border rounded-lg p-6">
        <h3 className="text-lg font-quando font-bold mb-6">Energy Comparison Across Experiments</h3>
        <ResponsiveContainer width="100%" height={400}>
          <BarChart data={energyComparison}>
            <CartesianGrid strokeDasharray="3 3" stroke="hsl(var(--border))" />
            <XAxis
              dataKey="id"
              label={{ value: "Experiment ID", position: "insideBottom", offset: -5 }}
              tick={{ fontSize: 11 }}
              stroke="hsl(var(--muted-foreground))"
            />
            <YAxis
              label={{ value: "Ground State Energy (Ha)", angle: -90, position: "insideLeft" }}
              tick={{ fontSize: 11 }}
              stroke="hsl(var(--muted-foreground))"
            />
            <Tooltip
              contentStyle={{
                backgroundColor: "hsl(var(--card))",
                border: "1px solid hsl(var(--border))",
                borderRadius: "8px",
                fontSize: "12px",
              }}
              formatter={(value: any, name: string, props: any) => [
                `${value.toFixed(8)} Ha`,
                `Energy (${props.payload.method})`
              ]}
              labelFormatter={(label) => `Experiment ${label}`}
            />
            <Bar dataKey="energy" fill="#ea580c" radius={[4, 4, 0, 0]} />
          </BarChart>
        </ResponsiveContainer>
      </div>

      {/* Best Result Card */}
      {energyComparison.length > 0 && (
        <div className="bg-gradient-to-r from-green-50 to-green-100 dark:from-green-900/30 dark:to-green-900/20 border-2 border-green-500 rounded-lg p-6">
          <h3 className="text-lg font-quando font-bold text-green-700 dark:text-green-300 mb-4 flex items-center gap-2">
            <TrendingDown className="w-5 h-5" />
            Lowest Energy Found
          </h3>
          <div className="grid grid-cols-3 gap-4">
            <div>
              <div className="text-sm text-muted-foreground mb-1">Experiment</div>
              <div className="font-quanto font-semibold">
                {energyComparison.find(e => e.energy === lowestEnergy)?.id}
              </div>
            </div>
            <div>
              <div className="text-sm text-muted-foreground mb-1">Energy</div>
              <div className="text-2xl font-quanto font-bold text-green-700 dark:text-green-300">
                {lowestEnergy.toFixed(8)} <span className="text-base">Ha</span>
              </div>
            </div>
            <div>
              <div className="text-sm text-muted-foreground mb-1">Method</div>
              <div className="font-quanto font-semibold">
                {energyComparison.find(e => e.energy === lowestEnergy)?.method}
              </div>
            </div>
          </div>
        </div>
      )}

      {/* Method Performance Table */}
      <div className="bg-card border border-border rounded-lg p-6">
        <h3 className="text-lg font-quando font-bold mb-6">Detailed Results</h3>
        <div className="overflow-x-auto">
          <table className="w-full">
            <thead>
              <tr className="border-b border-border">
                <th className="text-left py-3 px-4 text-sm font-quando font-semibold">Experiment</th>
                <th className="text-left py-3 px-4 text-sm font-quando font-semibold">Method</th>
                <th className="text-right py-3 px-4 text-sm font-quando font-semibold">Energy (Ha)</th>
                <th className="text-center py-3 px-4 text-sm font-quando font-semibold">Converged</th>
              </tr>
            </thead>
            <tbody>
              {energyComparison.map((exp) => (
                <tr key={exp.id} className="border-b border-border hover:bg-accent transition">
                  <td className="py-3 px-4 font-mono text-sm">{exp.id}</td>
                  <td className="py-3 px-4 font-quanto text-sm">{exp.method}</td>
                  <td className="py-3 px-4 font-mono text-sm text-right">
                    {exp.energy.toFixed(8)}
                  </td>
                  <td className="py-3 px-4 text-center">
                    {exp.converged ? (
                      <CheckCircle2 className="w-4 h-4 text-green-600 mx-auto" />
                    ) : (
                      <XCircle className="w-4 h-4 text-red-600 mx-auto" />
                    )}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      </div>
    </div>
  );
}
