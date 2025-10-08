"use client";

import { useState } from "react";
import MoleculeCreator from "@/components/molecule/MoleculeCreator";
import LewisStructureView from "@/components/molecule/LewisStructureView";
import PreviewWindow from "@/components/simulation/PreviewWindow";

type WorkflowStep = "idle" | "create" | "review" | "preview" | "running";

export default function DashboardPage() {
  const [currentStep, setCurrentStep] = useState<WorkflowStep>("idle");
  const [molecule, setMolecule] = useState<any>(null);

  // Mock backend settings - in production, this comes from settings
  const backendSettings = {
    method: "VQE",
    ansatz: "Hardware-Efficient",
    mapper: "Jordan-Wigner",
    backend: "Classical",
    optimizer: "SLSQP",
  };

  if (currentStep === "create") {
    return (
      <MoleculeCreator
        onComplete={(mol) => {
          setMolecule(mol);
          setCurrentStep("review");
        }}
        onCancel={() => setCurrentStep("idle")}
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
          console.log("Executing experiment:", config);
          // TODO: Submit to backend and move to running step
          setCurrentStep("running");
        }}
      />
    );
  }

  if (currentStep === "running") {
    return (
      <div className="h-full flex items-center justify-center">
        <div className="text-center">
          <div className="text-2xl font-quando font-bold mb-4">
            Experiment Running...
          </div>
          <div className="text-muted-foreground font-quando">
            Job monitoring view coming soon
          </div>
          <button
            onClick={() => setCurrentStep("idle")}
            className="mt-6 px-6 py-3 border border-border rounded-lg hover:bg-accent transition font-quando"
          >
            Back to Dashboard
          </button>
        </div>
      </div>
    );
  }

  return (
    <div
      onClick={() => setCurrentStep("create")}
      className="flex items-center justify-center h-full cursor-pointer hover:bg-accent transition-colors"
    >
      <h1 className="text-4xl font-quando text-muted-foreground">
        click anywhere to start
      </h1>
    </div>
  );
}
