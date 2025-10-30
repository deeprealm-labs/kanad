"use client";

import { useState, useEffect } from "react";
import { X, Beaker, Sparkles, Cpu, Pill, FlaskConical, Loader2, CheckCircle2, XCircle } from "lucide-react";
import * as api from "@/lib/api";
import { useToast } from "@/components/ui/toast";

interface AdvancedAnalysisModalProps {
  experimentId: string;
  experimentName?: string;
  onClose: () => void;
}

interface AnalysisProfile {
  name: string;
  display_name: string;
  description: string;
  modules: string[];
}

const PROFILE_ICONS = {
  drug_discovery: Pill,
  spectroscopy: Sparkles,
  materials: Cpu,
  catalysis: FlaskConical,
};

const PROFILE_COLORS = {
  drug_discovery: "pink",
  spectroscopy: "purple",
  materials: "green",
  catalysis: "orange",
};

export default function AdvancedAnalysisModal({
  experimentId,
  experimentName,
  onClose,
}: AdvancedAnalysisModalProps) {
  const [profiles, setProfiles] = useState<AnalysisProfile[]>([]);
  const [selectedProfile, setSelectedProfile] = useState<string | null>(null);
  const [isLoading, setIsLoading] = useState(false);
  const [isRunning, setIsRunning] = useState(false);
  const [results, setResults] = useState<any>(null);
  const [error, setError] = useState<string | null>(null);
  const toast = useToast();

  useEffect(() => {
    loadProfiles();
  }, []);

  const loadProfiles = async () => {
    setIsLoading(true);
    try {
      const response = await api.getAnalysisProfiles();
      console.log("Profiles response:", response);

      // Handle different response formats
      let profilesData = [];
      if (Array.isArray(response)) {
        // Already an array
        profilesData = response;
      } else if (response && response.profiles) {
        // Response has a "profiles" field
        if (Array.isArray(response.profiles)) {
          // profiles is an array
          profilesData = response.profiles;
        } else if (typeof response.profiles === 'object') {
          // profiles is an object - convert to array with keys
          profilesData = Object.entries(response.profiles).map(([key, profile]: [string, any]) => ({
            name: key,
            display_name: profile.name || key,
            description: profile.description || '',
            modules: profile.analyses || [],
          }));
        }
      } else if (response && typeof response === 'object') {
        // Response is an object - convert to array
        profilesData = Object.entries(response).map(([key, profile]: [string, any]) => ({
          name: key,
          display_name: profile.name || key,
          description: profile.description || '',
          modules: profile.analyses || [],
        }));
      }

      console.log("Profiles data:", profilesData);
      setProfiles(Array.isArray(profilesData) ? profilesData : []);
    } catch (err: any) {
      console.error("Failed to load profiles:", err);
      toast.error("Failed to load analysis profiles");
      setProfiles([]);
    } finally {
      setIsLoading(false);
    }
  };

  const runAnalysis = async () => {
    if (!selectedProfile) return;

    setIsRunning(true);
    setError(null);
    setResults(null);

    try {
      const response = await api.runAdvancedAnalysis(experimentId, selectedProfile);
      console.log("Analysis results:", response);

      if (response.status === "success" || response.status === "completed") {
        setResults(response.results);
        toast.success("Analysis completed successfully!");
      } else if (response.status === "partial") {
        setResults(response.results);
        toast.warning("Analysis completed with some failures");
      } else {
        setError(response.error || "Analysis failed");
        toast.error("Analysis failed");
      }
    } catch (err: any) {
      console.error("Analysis failed:", err);
      setError(err.message || "Failed to run analysis");
      toast.error("Failed to run analysis");
    } finally {
      setIsRunning(false);
    }
  };

  const getProfileColor = (profileName: string) => {
    return PROFILE_COLORS[profileName as keyof typeof PROFILE_COLORS] || "gray";
  };

  const getProfileIcon = (profileName: string) => {
    const IconComponent = PROFILE_ICONS[profileName as keyof typeof PROFILE_ICONS] || Beaker;
    return IconComponent;
  };

  const renderResults = () => {
    if (!results) return null;

    return (
      <div className="space-y-4">
        <h3 className="text-lg font-quando font-semibold">Analysis Results</h3>

        {Object.entries(results).map(([moduleName, moduleResults]: [string, any]) => {
          const hasError = moduleResults?.error;
          const hasData = moduleResults && !hasError;

          return (
            <div
              key={moduleName}
              className={`border rounded-lg p-4 ${
                hasError
                  ? "border-red-300 bg-red-50 dark:border-red-800 dark:bg-red-900/20"
                  : "border-border bg-card"
              }`}
            >
              <div className="flex items-start justify-between mb-2">
                <h4 className="font-quando font-semibold capitalize flex items-center gap-2">
                  {hasError ? (
                    <XCircle className="w-4 h-4 text-red-600" />
                  ) : (
                    <CheckCircle2 className="w-4 h-4 text-green-600" />
                  )}
                  {moduleName.replace(/_/g, " ")}
                </h4>
              </div>

              {hasError ? (
                <p className="text-sm text-red-600 dark:text-red-400">{moduleResults.error}</p>
              ) : (
                <div className="space-y-2 text-sm">
                  {/* ADME Results */}
                  {moduleName === "adme" && moduleResults.adme_properties && (
                    <div className="grid grid-cols-2 gap-2">
                      <div className="bg-muted rounded p-2">
                        <div className="text-muted-foreground mb-1 text-xs">Molecular Weight</div>
                        <div className="font-semibold">
                          {moduleResults.descriptors?.molecular_weight?.toFixed(2) || "N/A"} g/mol
                        </div>
                      </div>
                      <div className="bg-muted rounded p-2">
                        <div className="text-muted-foreground mb-1 text-xs">LogP</div>
                        <div className="font-semibold">
                          {moduleResults.adme_properties?.logP?.toFixed(2) || "N/A"}
                        </div>
                      </div>
                      <div className="bg-muted rounded p-2">
                        <div className="text-muted-foreground mb-1 text-xs">BBB Permeability</div>
                        <div className="font-semibold">
                          {moduleResults.adme_properties?.bbb_permeability || "N/A"}
                        </div>
                      </div>
                      <div className="bg-muted rounded p-2">
                        <div className="text-muted-foreground mb-1 text-xs">Drug Likeness</div>
                        <div className="font-semibold">
                          {moduleResults.adme_properties?.drug_likeness || "N/A"}
                        </div>
                      </div>
                    </div>
                  )}

                  {/* Energy Analysis Results */}
                  {moduleName === "energy" && moduleResults.total_energy !== undefined && (
                    <div className="grid grid-cols-2 gap-2">
                      <div className="bg-blue-50 dark:bg-blue-900/20 border border-blue-200 dark:border-blue-800 rounded p-2">
                        <div className="text-blue-700 dark:text-blue-300 mb-1 text-xs">Total Energy</div>
                        <div className="font-mono font-semibold">
                          {moduleResults.total_energy.toFixed(6)} Ha
                        </div>
                      </div>
                      {moduleResults.electronic_energy !== undefined && (
                        <div className="bg-muted rounded p-2">
                          <div className="text-muted-foreground mb-1 text-xs">Electronic Energy</div>
                          <div className="font-mono font-semibold">
                            {moduleResults.electronic_energy.toFixed(6)} Ha
                          </div>
                        </div>
                      )}
                      {moduleResults.nuclear_repulsion !== undefined && (
                        <div className="bg-muted rounded p-2">
                          <div className="text-muted-foreground mb-1 text-xs">Nuclear Repulsion</div>
                          <div className="font-mono font-semibold">
                            {moduleResults.nuclear_repulsion.toFixed(6)} Ha
                          </div>
                        </div>
                      )}
                      {moduleResults.energy_per_electron !== undefined && (
                        <div className="bg-muted rounded p-2">
                          <div className="text-muted-foreground mb-1 text-xs">Energy per Electron</div>
                          <div className="font-mono font-semibold">
                            {moduleResults.energy_per_electron.toFixed(4)} Ha/e⁻
                          </div>
                        </div>
                      )}
                    </div>
                  )}

                  {/* Correlation Analysis Results */}
                  {moduleName === "correlation" && moduleResults.correlation_energy !== undefined && (
                    <div className="grid grid-cols-2 gap-2">
                      <div className="bg-muted rounded p-2">
                        <div className="text-muted-foreground mb-1 text-xs">HF Energy</div>
                        <div className="font-mono font-semibold">
                          {moduleResults.hf_energy.toFixed(6)} Ha
                        </div>
                      </div>
                      <div className="bg-muted rounded p-2">
                        <div className="text-muted-foreground mb-1 text-xs">VQE Energy</div>
                        <div className="font-mono font-semibold">
                          {moduleResults.vqe_energy.toFixed(6)} Ha
                        </div>
                      </div>
                      <div className="bg-purple-50 dark:bg-purple-900/20 border border-purple-200 dark:border-purple-800 rounded p-2">
                        <div className="text-purple-700 dark:text-purple-300 mb-1 text-xs">Correlation Energy</div>
                        <div className="font-mono font-semibold">
                          {moduleResults.correlation_energy.toFixed(6)} Ha
                        </div>
                      </div>
                      {moduleResults.correlation_percentage !== undefined && (
                        <div className="bg-muted rounded p-2">
                          <div className="text-muted-foreground mb-1 text-xs">Correlation %</div>
                          <div className="font-mono font-semibold">
                            {moduleResults.correlation_percentage.toFixed(2)}%
                          </div>
                        </div>
                      )}
                    </div>
                  )}

                  {/* Generic JSON display for other modules */}
                  {moduleName !== "adme" && moduleName !== "energy" && moduleName !== "correlation" && (
                    <pre className="bg-muted rounded p-2 text-xs overflow-x-auto max-h-48 overflow-y-auto">
                      {JSON.stringify(moduleResults, null, 2)}
                    </pre>
                  )}
                </div>
              )}
            </div>
          );
        })}
      </div>
    );
  };

  return (
    <div className="fixed inset-0 z-50 flex items-center justify-center bg-black/50 backdrop-blur-sm">
      <div className="bg-background border border-border rounded-lg shadow-xl w-full max-w-3xl max-h-[90vh] overflow-hidden flex flex-col">
        {/* Header */}
        <div className="flex items-center justify-between p-6 border-b border-border">
          <div>
            <h2 className="text-xl font-quando font-bold">Advanced Analysis</h2>
            {experimentName && (
              <p className="text-sm text-muted-foreground mt-1">{experimentName}</p>
            )}
          </div>
          <button
            onClick={onClose}
            className="p-2 hover:bg-accent rounded-lg transition"
            disabled={isRunning}
          >
            <X className="w-5 h-5" />
          </button>
        </div>

        {/* Content */}
        <div className="flex-1 overflow-y-auto p-6">
          {isLoading ? (
            <div className="flex items-center justify-center py-12">
              <Loader2 className="w-8 h-8 animate-spin text-brand-orange" />
            </div>
          ) : results ? (
            renderResults()
          ) : (
            <>
              <p className="text-sm text-muted-foreground mb-6">
                Select an analysis profile to run advanced domain-specific analysis on this experiment.
              </p>

              <div className="space-y-3">
                {Array.isArray(profiles) && profiles.length > 0 ? (
                  profiles.map((profile) => {
                    const Icon = getProfileIcon(profile.name);
                    const color = getProfileColor(profile.name);
                    const isSelected = selectedProfile === profile.name;

                  return (
                    <button
                      key={profile.name}
                      onClick={() => setSelectedProfile(profile.name)}
                      className={`w-full text-left p-4 rounded-lg border-2 transition ${
                        isSelected
                          ? `border-${color}-500 bg-${color}-50 dark:bg-${color}-900/20`
                          : "border-border hover:border-muted-foreground bg-card"
                      }`}
                      disabled={isRunning}
                    >
                      <div className="flex items-start gap-3">
                        <div
                          className={`p-2 rounded-lg ${
                            isSelected
                              ? `bg-${color}-100 dark:bg-${color}-800/30`
                              : "bg-muted"
                          }`}
                        >
                          <Icon
                            className={`w-5 h-5 ${
                              isSelected
                                ? `text-${color}-600 dark:text-${color}-400`
                                : "text-muted-foreground"
                            }`}
                          />
                        </div>
                        <div className="flex-1">
                          <h3 className="font-quando font-semibold mb-1">
                            {profile.display_name}
                          </h3>
                          <p className="text-sm text-muted-foreground mb-2">
                            {profile.description}
                          </p>
                          <div className="flex flex-wrap gap-1">
                            {Array.isArray(profile.modules) && profile.modules.map((module) => (
                              <span
                                key={module}
                                className="text-xs px-2 py-1 bg-muted rounded"
                              >
                                {module}
                              </span>
                            ))}
                          </div>
                        </div>
                      </div>
                    </button>
                  );
                })
                ) : (
                  <div className="text-center py-8 text-muted-foreground">
                    {isLoading ? "Loading profiles..." : "No analysis profiles available"}
                  </div>
                )}
              </div>

              {error && (
                <div className="mt-4 p-4 bg-red-50 dark:bg-red-900/20 border border-red-300 dark:border-red-800 rounded-lg">
                  <p className="text-sm text-red-600 dark:text-red-400">{error}</p>
                </div>
              )}
            </>
          )}
        </div>

        {/* Footer */}
        <div className="flex items-center justify-end gap-3 p-6 border-t border-border">
          {results ? (
            <>
              <button
                onClick={() => {
                  setResults(null);
                  setSelectedProfile(null);
                }}
                className="px-4 py-2 border border-border rounded-lg hover:bg-accent transition font-quando"
              >
                Run Another
              </button>
              <button
                onClick={onClose}
                className="px-4 py-2 bg-brand-orange text-white rounded-lg hover:bg-brand-orange/90 transition font-quando"
              >
                Done
              </button>
            </>
          ) : (
            <>
              <button
                onClick={onClose}
                className="px-4 py-2 border border-border rounded-lg hover:bg-accent transition font-quando"
                disabled={isRunning}
              >
                Cancel
              </button>
              <button
                onClick={runAnalysis}
                disabled={!selectedProfile || isRunning}
                className="px-4 py-2 bg-brand-orange text-white rounded-lg hover:bg-brand-orange/90 transition font-quando disabled:opacity-50 disabled:cursor-not-allowed flex items-center gap-2"
              >
                {isRunning ? (
                  <>
                    <Loader2 className="w-4 h-4 animate-spin" />
                    Running...
                  </>
                ) : (
                  "Run Analysis"
                )}
              </button>
            </>
          )}
        </div>
      </div>
    </div>
  );
}
