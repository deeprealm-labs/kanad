"use client";

import { useState, useEffect } from "react";
import MoleculeCreator from "@/components/molecule/MoleculeCreator";
import LewisStructureView from "@/components/molecule/LewisStructureView";
import PreviewWindow from "@/components/simulation/PreviewWindow";
import ExperimentMonitor from "@/components/simulation/ExperimentMonitor";
import DashboardHome from "@/components/dashboard/DashboardHome";
import { useToast } from "@/components/ui/toast";
import { FullPageLoader } from "@/components/ui/loading";
import * as api from "@/lib/api";
import type { Experiment, BackendSettings } from "@/lib/types";

type WorkflowStep = "home" | "create" | "review" | "preview" | "running";

export default function DashboardPage() {
  const [currentStep, setCurrentStep] = useState<WorkflowStep>("home");
  const [molecule, setMolecule] = useState<any>(null);
  const [experimentConfig, setExperimentConfig] = useState<any>(null);
  const [experiments, setExperiments] = useState<Experiment[]>([]);
  const [currentExperimentId, setCurrentExperimentId] = useState<string | null>(null);
  const [loading, setLoading] = useState(false);
  const [backendSettings, setBackendSettings] = useState<BackendSettings>({
    method: "VQE",
    ansatz: "hardware_efficient",
    mapper: "jordan_wigner",
    backend: "classical",
    optimizer: "SLSQP",
  });

  const toast = useToast();

  // Load experiments from API
  useEffect(() => {
    loadExperiments();
    loadSettings();
  }, []);

  const loadExperiments = async () => {
    try {
      // Try to load from API
      const response = await api.getExperiments({ limit: 10 });
      setExperiments(response.experiments);
    } catch (error: any) {
      console.error("Failed to load experiments from API:", error);

      // Fallback to localStorage if API is not available
      const saved = localStorage.getItem("kanad_experiments");
      if (saved) {
        try {
          setExperiments(JSON.parse(saved));
        } catch (e) {
          console.error("Failed to parse localStorage experiments:", e);
        }
      }

      // Show toast only if it's not a connection error
      if (error.statusCode !== 503) {
        toast.error("Failed to load experiments. Using local data.");
      }
    }
  };

  const loadSettings = async () => {
    try {
      const settings = await api.getSettings();
      setBackendSettings(settings as any);
    } catch (error) {
      console.error("Failed to load settings from API:", error);

      // Fallback to localStorage
      const saved = localStorage.getItem("kanad_settings");
      if (saved) {
        try {
          setBackendSettings(JSON.parse(saved));
        } catch (e) {
          console.error("Failed to parse localStorage settings:", e);
        }
      }
    }
  };

  // Save experiment
  const saveExperiment = async (exp: Partial<Experiment>) => {
    try {
      setLoading(true);

      // Create experiment via API
      const response = await api.createExperiment({
        molecule: exp.molecule!,
        configuration: exp.backendSettings!,
        analysis: exp.analysis,
        execute_now: false, // Just save, don't execute
      });

      toast.success("Experiment saved successfully!");

      // Reload experiments
      await loadExperiments();
    } catch (error: any) {
      console.error("Failed to save experiment:", error);
      toast.error(error.message || "Failed to save experiment");

      // Fallback to localStorage
      const newExperiment: Experiment = {
        id: Date.now().toString(),
        timestamp: new Date().toISOString(),
        status: "completed",
        method: exp.backendSettings?.method || "VQE",
        backend: exp.backendSettings?.backend || "classical",
        ...exp,
      } as Experiment;

      const updated = [newExperiment, ...experiments];
      setExperiments(updated);
      localStorage.setItem("kanad_experiments", JSON.stringify(updated));
    } finally {
      setLoading(false);
    }
  };

  // Add job to queue
  const addToQueue = async (config: any) => {
    try {
      setLoading(true);

      // Submit to API (execute_now: false means add to queue)
      const response = await api.createExperiment({
        molecule: config.molecule,
        configuration: config.backendSettings || backendSettings,
        analysis: config.analysis,
        execute_now: false,
      });

      toast.success("Experiment added to queue!");
      setCurrentStep("home");
    } catch (error: any) {
      console.error("Failed to add to queue:", error);
      toast.error(error.message || "Failed to add experiment to queue");

      // Fallback to localStorage
      const saved = localStorage.getItem("kanad_queue");
      const queue = saved ? JSON.parse(saved) : [];

      const newJob = {
        id: Date.now().toString(),
        name: config.molecule?.smiles || "Custom Molecule",
        molecule: config.molecule,
        method: config.backendSettings?.method || backendSettings.method,
        backend: config.backendSettings?.backend || backendSettings.backend,
        priority: queue.length + 1,
        status: "queued" as const,
        createdAt: new Date().toISOString(),
        backendSettings: config.backendSettings || backendSettings,
        analysis: config.analysis,
      };

      const updated = [...queue, newJob];
      localStorage.setItem("kanad_queue", JSON.stringify(updated));
      toast.success("Experiment added to queue (offline mode)");
      setCurrentStep("home");
    } finally {
      setLoading(false);
    }
  };

  // Execute experiment
  const executeExperiment = async (config: any) => {
    try {
      setLoading(true);

      // Debug logging
      console.log("executeExperiment - config:", config);
      console.log("executeExperiment - config.backendSettings:", config.backendSettings);
      console.log("executeExperiment - dashboard backendSettings:", backendSettings);

      // Ensure we have valid backend settings with defaults
      const settings = config.backendSettings || backendSettings || {
        method: "VQE",
        ansatz: "ucc",
        mapper: "jordan_wigner",
        backend: "classical",
        optimizer: "SLSQP",
      };

      console.log("executeExperiment - using settings:", settings);

      // Submit to API with execute_now: true
      const response = await api.createExperiment({
        molecule: config.molecule,
        configuration: settings,
        analysis: config.analysis,
        execute_now: true,
      });

      console.log("Experiment submitted successfully:", response);

      // Get the experiment ID from response
      const experimentId = response.experiment_id || response.experimentId || response.id;

      if (!experimentId) {
        throw new Error("No experiment ID returned from API");
      }

      setCurrentExperimentId(experimentId);
      setExperimentConfig(config);
      setCurrentStep("running");
      toast.success("Experiment started!");
    } catch (error: any) {
      console.error("Failed to execute experiment:", error);
      toast.error(error.message || "Failed to start experiment");

      // Don't fall back to offline mode - show error
      setLoading(false);
      return;
    } finally {
      setLoading(false);
    }
  };

  if (loading) {
    return <FullPageLoader message="Processing..." />;
  }

  if (currentStep === "create") {
    return (
      <MoleculeCreator
        onComplete={(mol) => {
          setMolecule(mol);
          setCurrentStep("review");
        }}
        onCancel={() => setCurrentStep("home")}
      />
    );
  }

  if (currentStep === "review" && molecule) {
    return (
      <LewisStructureView
        molecule={molecule}
        onExecute={() => setCurrentStep("preview")}
      />
    );
  }

  if (currentStep === "preview" && molecule) {
    return (
      <PreviewWindow
        molecule={molecule}
        backendSettings={backendSettings}
        onBack={() => setCurrentStep("review")}
        onExecute={(config) => {
          executeExperiment(config);
        }}
        onQueue={(config) => {
          addToQueue(config);
        }}
        onRefreshSettings={loadSettings}
        onOpenSettings={() => setCurrentStep("review")}
      />
    );
  }

  if (currentStep === "running") {
    return (
      <ExperimentMonitor
        experimentId={currentExperimentId}
        experimentConfig={experimentConfig}
        onComplete={async () => {
          await loadExperiments();
          setCurrentStep("home");
        }}
        onBack={() => setCurrentStep("home")}
        onQueueAnother={() => setCurrentStep("create")}
      />
    );
  }

  // Home view
  return (
  
    <DashboardHome
      onNewExperiment={() => setCurrentStep("create")}
      onViewExperiment={(experimentId, config) => {
        setCurrentExperimentId(experimentId);
        setExperimentConfig(config);
        setCurrentStep("running");
      }}
      experiments={experiments}
    />
       

  );
}
