"use client";

import { useState } from "react";
import MoleculeBuilder from "@/components/molecule/MoleculeBuilder";
import SimulationConfig from "@/components/simulation/SimulationConfig";

type WorkflowStep = "idle" | "molecule" | "config" | "monitor" | "results";

export default function DashboardPage() {
  const [currentStep, setCurrentStep] = useState<WorkflowStep>("idle");

  if (currentStep === "molecule") {
    return (
      <MoleculeBuilder
        onClose={() => setCurrentStep("idle")}
        onComplete={() => setCurrentStep("config")}
      />
    );
  }

  if (currentStep === "config") {
    return (
      <SimulationConfig
        onBack={() => setCurrentStep("molecule")}
        onComplete={(config) => {
          console.log("Configuration:", config);
          // TODO: Submit job and move to monitor step
          setCurrentStep("idle");
        }}
      />
    );
  }

  return (
    <div
      onClick={() => setCurrentStep("molecule")}
      className="flex items-center justify-center h-full cursor-pointer hover:bg-accent transition-colors"
    >
      <h1 className="text-4xl font-quando text-muted-foreground">
        click anywhere to start
      </h1>
    </div>
  );
}
