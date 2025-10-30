"use client";

import { useState, useEffect } from "react";
import { X, TrendingDown, Award, Zap } from "lucide-react";
import {
  LineChart,
  Line,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  Legend,
  ResponsiveContainer,
} from "recharts";
import * as api from "@/lib/api";
import { useToast } from "@/components/ui/toast";

interface ExperimentComparisonProps {
  experimentIds: string[];
  onClose: () => void;
}

const COLORS = [
  "#ea580c", // brand-orange
  "#3b82f6", // blue
  "#10b981", // green
  "#f59e0b", // amber
  "#8b5cf6", // purple
  "#ef4444", // red
  "#06b6d4", // cyan
  "#ec4899", // pink
  "#84cc16", // lime
  "#6366f1", // indigo
];

export default function ExperimentComparison({
  experimentIds,
  onClose,
}: ExperimentComparisonProps) {
  const [comparisonData, setComparisonData] = useState<any>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const toast = useToast();

  useEffect(() => {
    loadComparisonData();
  }, [experimentIds]);

  const loadComparisonData = async () => {
    try {
      setLoading(true);
      const data = await api.compareExperiments(experimentIds);
      setComparisonData(data);
      setError(null);
    } catch (err: any) {
      console.error("Failed to load comparison:", err);
      setError(err.message || "Failed to load comparison data");
    } finally {
      setLoading(false);
    }
  };

  // Merge convergence data from all experiments for overlay chart
  const getMergedConvergenceData = () => {
    if (!comparisonData?.experiments) return [];

    const maxLength = Math.max(
      ...comparisonData.experiments.map(
        (exp: any) => exp.convergence_history?.length || exp.energy_history?.length || 0
      )
    );

    const mergedData = [];
    for (let i = 0; i < maxLength; i++) {
      const point: any = { iteration: i + 1 };
      comparisonData.experiments.forEach((exp: any, idx: number) => {
        const convergence = exp.convergence_history || exp.energy_history?.map((e: number, i: number) => ({ iteration: i + 1, energy: e })) || [];
        if (convergence[i]) {
          point[`exp_${idx}`] = convergence[i].energy;
        }
      });
      mergedData.push(point);
    }
    return mergedData;
  };

  if (loading) {
    return (
      <div className="fixed inset-0 z-50 flex items-center justify-center backdrop-blur-sm bg-black/50">
        <div className="bg-background border border-border rounded-lg p-12">
          <div className="w-12 h-12 border-4 border-brand-orange border-t-transparent rounded-full animate-spin mx-auto mb-4"></div>
          <p className="text-muted-foreground font-quando">Loading comparison...</p>
        </div>
      </div>
    );
  }

  if (error || !comparisonData) {
    return (
      <div className="fixed inset-0 z-50 flex items-center justify-center backdrop-blur-sm bg-black/50">
        <div className="bg-background border border-border rounded-lg p-12 max-w-md">
          <div className="text-center">
            <X className="w-12 h-12 text-red-500 mx-auto mb-4" />
            <h3 className="text-lg font-quando font-bold mb-2">
              Failed to Load Comparison
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

  const { experiments, comparison_stats } = comparisonData;
  const convergenceData = getMergedConvergenceData();

  return (
    <div className="fixed inset-0 z-50 flex items-center justify-center backdrop-blur-sm bg-black/50 p-6">
      <div className="bg-background border border-border rounded-lg w-full max-w-7xl max-h-[90vh] flex flex-col">
        {/* Header */}
        <div className="flex items-center justify-between p-6 border-b border-border">
          <div>
            <h2 className="text-2xl font-quando font-bold">Experiment Comparison</h2>
            <p className="text-sm text-muted-foreground">
              Comparing {experiments.length} experiments
            </p>
          </div>
          <button
            onClick={onClose}
            className="p-2 hover:bg-accent rounded-lg transition"
          >
            <X className="w-5 h-5" />
          </button>
        </div>

        {/* Content - Scrollable */}
        <div className="flex-1 overflow-auto p-6 space-y-6">
          {/* Comparison Stats */}
          <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
            <div className="bg-gradient-to-br from-green-500/10 to-green-500/5 border border-green-500/20 rounded-lg p-4">
              <div className="flex items-center gap-2 text-green-600 mb-2">
                <Award className="w-5 h-5" />
                <span className="text-xs font-quando font-medium">Best Energy</span>
              </div>
              <div className="text-2xl font-quando font-bold">
                {comparison_stats.lowest_energy?.toFixed(6)}
                <span className="text-sm font-normal ml-1">Ha</span>
              </div>
              {comparison_stats.best_experiment && (
                <div className="text-xs text-muted-foreground mt-1">
                  {comparison_stats.best_experiment.molecule}
                </div>
              )}
            </div>

            <div className="bg-card border border-border rounded-lg p-4">
              <div className="flex items-center gap-2 text-muted-foreground mb-2">
                <TrendingDown className="w-5 h-5" />
                <span className="text-xs font-quando font-medium">Energy Range</span>
              </div>
              <div className="text-2xl font-quando font-bold">
                {comparison_stats.energy_range?.toFixed(6)}
                <span className="text-sm font-normal ml-1">Ha</span>
              </div>
            </div>

            <div className="bg-card border border-border rounded-lg p-4">
              <div className="flex items-center gap-2 text-muted-foreground mb-2">
                <Zap className="w-5 h-5" />
                <span className="text-xs font-quando font-medium">Average Energy</span>
              </div>
              <div className="text-2xl font-quando font-bold">
                {comparison_stats.average_energy?.toFixed(6)}
                <span className="text-sm font-normal ml-1">Ha</span>
              </div>
            </div>

            <div className="bg-card border border-border rounded-lg p-4">
              <div className="text-xs font-quando font-medium text-muted-foreground mb-2">
                Status
              </div>
              <div className="text-2xl font-quando font-bold">
                {comparison_stats.completed}/{comparison_stats.total_compared}
              </div>
              <div className="text-xs text-muted-foreground">completed</div>
            </div>
          </div>

          {/* Convergence Overlay Chart */}
          {convergenceData.length > 0 && (
            <div className="bg-card border border-border rounded-lg p-6">
              <h3 className="text-lg font-quando font-semibold mb-4">
                Convergence Comparison
              </h3>
              <div className="h-64">
                <ResponsiveContainer width="100%" height="100%">
                  <LineChart data={convergenceData}>
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
                    />
                    <Legend />
                    {experiments.map((exp: any, idx: number) => (
                      <Line
                        key={exp.id}
                        type="monotone"
                        dataKey={`exp_${idx}`}
                        name={exp.molecule?.smiles || `Exp ${idx + 1}`}
                        stroke={COLORS[idx % COLORS.length]}
                        strokeWidth={2}
                        dot={false}
                      />
                    ))}
                  </LineChart>
                </ResponsiveContainer>
              </div>
            </div>
          )}

          {/* Side-by-Side Comparison Table */}
          <div className="bg-card border border-border rounded-lg p-6 overflow-x-auto">
            <h3 className="text-lg font-quando font-semibold mb-4">
              Detailed Comparison
            </h3>
            <table className="w-full text-sm">
              <thead>
                <tr className="border-b border-border">
                  <th className="text-left py-2 px-3 font-quando font-medium text-muted-foreground">
                    Property
                  </th>
                  {experiments.map((exp: any, idx: number) => (
                    <th
                      key={exp.id}
                      className="text-right py-2 px-3 font-quando font-medium"
                      style={{ color: COLORS[idx % COLORS.length] }}
                    >
                      Exp {idx + 1}
                    </th>
                  ))}
                </tr>
              </thead>
              <tbody>
                <tr className="border-b border-border/50">
                  <td className="py-2 px-3 text-muted-foreground">Molecule</td>
                  {experiments.map((exp: any) => (
                    <td key={exp.id} className="text-right py-2 px-3 font-mono text-xs">
                      {exp.molecule?.smiles || "N/A"}
                    </td>
                  ))}
                </tr>
                <tr className="border-b border-border/50">
                  <td className="py-2 px-3 text-muted-foreground">Method</td>
                  {experiments.map((exp: any) => (
                    <td key={exp.id} className="text-right py-2 px-3 font-quando">
                      {exp.method}
                    </td>
                  ))}
                </tr>
                <tr className="border-b border-border/50">
                  <td className="py-2 px-3 text-muted-foreground">Energy (Ha)</td>
                  {experiments.map((exp: any) => (
                    <td key={exp.id} className="text-right py-2 px-3 font-mono font-bold">
                      {exp.energy?.toFixed(6) || "N/A"}
                    </td>
                  ))}
                </tr>
                <tr className="border-b border-border/50">
                  <td className="py-2 px-3 text-muted-foreground">HF Energy (Ha)</td>
                  {experiments.map((exp: any) => (
                    <td key={exp.id} className="text-right py-2 px-3 font-mono">
                      {exp.hf_energy?.toFixed(6) || "N/A"}
                    </td>
                  ))}
                </tr>
                <tr className="border-b border-border/50">
                  <td className="py-2 px-3 text-muted-foreground">Correlation (Ha)</td>
                  {experiments.map((exp: any) => (
                    <td key={exp.id} className="text-right py-2 px-3 font-mono">
                      {exp.correlation_energy?.toFixed(6) || "N/A"}
                    </td>
                  ))}
                </tr>
                <tr className="border-b border-border/50">
                  <td className="py-2 px-3 text-muted-foreground">Iterations</td>
                  {experiments.map((exp: any) => (
                    <td key={exp.id} className="text-right py-2 px-3 font-quando">
                      {exp.iterations || "N/A"}
                    </td>
                  ))}
                </tr>
                <tr className="border-b border-border/50">
                  <td className="py-2 px-3 text-muted-foreground">Converged</td>
                  {experiments.map((exp: any) => (
                    <td key={exp.id} className="text-right py-2 px-3">
                      <span className={exp.converged ? "text-green-600" : "text-red-600"}>
                        {exp.converged ? "Yes" : "No"}
                      </span>
                    </td>
                  ))}
                </tr>
                <tr>
                  <td className="py-2 px-3 text-muted-foreground">Backend</td>
                  {experiments.map((exp: any) => (
                    <td key={exp.id} className="text-right py-2 px-3 font-quando">
                      {exp.backend}
                    </td>
                  ))}
                </tr>
              </tbody>
            </table>
          </div>
        </div>

        {/* Footer */}
        <div className="border-t border-border p-6">
          <button
            onClick={onClose}
            className="w-full px-6 py-3 bg-brand-orange text-white rounded-lg hover:bg-brand-orange-dark transition font-quando"
          >
            Close Comparison
          </button>
        </div>
      </div>
    </div>
  );
}
