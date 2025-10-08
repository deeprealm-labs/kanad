"use client";

import { useState } from "react";
import { ZoomIn, ZoomOut, RotateCcw } from "lucide-react";

interface LewisStructureViewProps {
  molecule: any;
  onExecute: () => void;
}

export default function LewisStructureView({
  molecule,
  onExecute,
}: LewisStructureViewProps) {
  const [zoom, setZoom] = useState(1);

  const handleZoomIn = () => setZoom(Math.min(zoom + 0.2, 3));
  const handleZoomOut = () => setZoom(Math.max(zoom - 0.2, 0.5));
  const handleReset = () => setZoom(1);

  // Simple example - in production, this would generate actual Lewis structure
  const renderLewisStructure = () => {
    if (molecule.smiles === "Cl" || molecule.smiles === "HCl" || molecule.smiles === "[H]Cl") {
      return (
        <div className="relative" style={{ transform: `scale(${zoom})`, transition: "transform 0.2s" }}>
          {/* H atom */}
          <div className="absolute" style={{ left: "200px", top: "200px" }}>
            <div className="relative">
              <div className="absolute -top-8 left-1/2 -translate-x-1/2 text-sm text-muted-foreground">
                6
              </div>
              <div className="absolute -top-6 right-0 w-6 h-6 bg-background border border-border rounded-full flex items-center justify-center text-xs">
                δ+
              </div>
              <div className="text-8xl font-bold">H</div>
              <div className="absolute -bottom-8 left-1/2 -translate-x-1/2 text-sm text-muted-foreground">
                4
              </div>
            </div>
          </div>

          {/* Bond */}
          <div className="absolute" style={{ left: "310px", top: "240px", width: "80px", height: "8px", background: "currentColor" }}></div>

          {/* Cl atom */}
          <div className="absolute" style={{ left: "400px", top: "200px" }}>
            <div className="relative">
              <div className="absolute -top-8 left-1/2 -translate-x-1/2 text-sm text-muted-foreground">
                7
              </div>
              <div className="absolute -top-6 left-0 w-6 h-6 bg-background border border-border rounded-full flex items-center justify-center text-xs">
                δ−
              </div>
              <div className="flex items-center gap-2">
                <div className="text-8xl font-bold">Cl</div>
                <div className="flex flex-col gap-1">
                  <div className="flex gap-1">
                    <div className="w-2 h-2 bg-foreground rounded-full"></div>
                    <div className="w-2 h-2 bg-foreground rounded-full"></div>
                  </div>
                  <div className="flex gap-1">
                    <div className="w-2 h-2 bg-foreground rounded-full"></div>
                    <div className="w-2 h-2 bg-foreground rounded-full"></div>
                  </div>
                  <div className="flex gap-1">
                    <div className="w-2 h-2 bg-foreground rounded-full"></div>
                    <div className="w-2 h-2 bg-foreground rounded-full"></div>
                  </div>
                </div>
              </div>
              <div className="absolute -bottom-8 left-8 text-sm text-muted-foreground flex gap-2">
                <span>17</span>
                <span>35</span>
                <span>11</span>
              </div>
            </div>
          </div>
        </div>
      );
    }

    // Generic molecule display
    return (
      <div className="text-center" style={{ transform: `scale(${zoom})`, transition: "transform 0.2s" }}>
        <div className="text-6xl font-bold font-mono mb-4">
          {molecule.smiles || molecule.atoms?.map((a: any) => a.symbol).join("-")}
        </div>
        <div className="text-muted-foreground font-quando">
          Basis: {molecule.basis} | Charge: {molecule.charge} | Multiplicity: {molecule.multiplicity}
        </div>
      </div>
    );
  };

  return (
    <div className="h-full flex items-center justify-center relative">
      {/* Zoom Controls */}
      <div className="absolute left-8 top-1/2 -translate-y-1/2 flex flex-col gap-2">
        <button
          onClick={handleReset}
          className="w-10 h-10 bg-card border border-border rounded-full flex items-center justify-center hover:bg-accent transition"
          title="Reset zoom"
        >
          <RotateCcw className="w-5 h-5" />
        </button>
        <button
          onClick={handleZoomIn}
          className="w-10 h-10 bg-card border border-border rounded-full flex items-center justify-center hover:bg-accent transition"
          title="Zoom in"
        >
          <ZoomIn className="w-5 h-5" />
        </button>
        <button
          onClick={handleZoomOut}
          className="w-10 h-10 bg-card border border-border rounded-full flex items-center justify-center hover:bg-accent transition"
          title="Zoom out"
        >
          <ZoomOut className="w-5 h-5" />
        </button>
      </div>

      {/* Lewis Structure */}
      <div className="flex items-center justify-center">
        {renderLewisStructure()}
      </div>

      {/* Execute Button */}
      <button
        onClick={onExecute}
        className="absolute bottom-8 right-8 px-6 py-3 bg-background border-2 border-border rounded-lg hover:border-brand-orange transition font-quando flex items-center gap-2"
      >
        <span>execute now</span>
      </button>
    </div>
  );
}
