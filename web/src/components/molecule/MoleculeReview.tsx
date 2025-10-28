"use client";

import { useState } from "react";
import { ArrowLeft, Play, Settings, Download, Share2, Info } from "lucide-react";
import Molecule3DViewer from "./Molecule3DViewer";
import { createMoleculeData } from "@/utils/moleculeUtils";

interface MoleculeReviewProps {
  molecule: any;
  onExecute: () => void;
  onBack: () => void;
}

export default function MoleculeReview({
  molecule,
  onExecute,
  onBack,
}: MoleculeReviewProps) {
  const [selectedAtom, setSelectedAtom] = useState<any>(null);
  const [showDetails, setShowDetails] = useState(false);

  const moleculeData = createMoleculeData(
    molecule.smiles,
    molecule.atoms,
    molecule.xyzData,
    molecule.name
  );

  const handleAtomClick = (atom: any) => {
    setSelectedAtom(atom);
    setShowDetails(true);
  };

  const getMoleculeInfo = () => {
    const info = [];
    
    if (molecule.smiles) {
      info.push({ label: "SMILES", value: molecule.smiles });
    }
    
    if (molecule.atoms && molecule.atoms.length > 0) {
      info.push({ label: "Atoms", value: molecule.atoms.length.toString() });
    }
    
    if (molecule.basis) {
      info.push({ label: "Basis Set", value: molecule.basis.toUpperCase() });
    }
    
    if (molecule.charge !== undefined) {
      info.push({ label: "Charge", value: molecule.charge.toString() });
    }
    
    if (molecule.multiplicity !== undefined) {
      info.push({ label: "Multiplicity", value: molecule.multiplicity.toString() });
    }

    return info;
  };

  const moleculeInfo = getMoleculeInfo();

  return (
    <div className="h-full flex flex-col bg-background">
      {/* Header */}
      <div className="flex items-center justify-between px-6 py-4 border-b border-border">
        <div className="flex items-center gap-4">
          <button
            onClick={onBack}
            className="p-2 hover:bg-accent rounded-md transition"
          >
            <ArrowLeft className="w-5 h-5" />
          </button>
          <div>
            <h1 className="text-xl font-quando font-bold">Molecule Review</h1>
            <p className="text-sm text-muted-foreground font-quando">
              {molecule.name || "Custom Molecule"}
            </p>
          </div>
        </div>
        
        <div className="flex items-center gap-2">
          <button
            onClick={() => setShowDetails(!showDetails)}
            className="p-2 hover:bg-accent rounded-md transition"
            title="Show Details"
          >
            <Info className="w-5 h-5" />
          </button>
          <button
            onClick={onExecute}
            className="px-4 py-2 bg-brand-orange text-white rounded-lg hover:bg-brand-orange-dark transition font-quando flex items-center gap-2"
          >
            <Play className="w-4 h-4" />
            Execute Experiment
          </button>
        </div>
      </div>

      {/* Main Content */}
      <div className="flex-1 grid grid-cols-1 lg:grid-cols-3 gap-6 p-6 overflow-hidden min-h-0">
        {/* Left: 3D Visualization */}
        <div className="lg:col-span-2 flex flex-col min-h-0">
          <div className="flex items-center justify-between mb-4">
            <h2 className="text-lg font-quando font-semibold">3D Structure</h2>
            <div className="flex items-center gap-2 text-sm text-muted-foreground font-quando">
              <span>Interactive 3D Model</span>
            </div>
          </div>
          
          <div className="flex-1 min-h-0 h-full">
            <Molecule3DViewer
              molecule={moleculeData}
              height="500px"
              showControls={true}
              onAtomClick={handleAtomClick}
              className="border border-border rounded-lg"
            />
          </div>
        </div>

        {/* Right: Details Panel */}
        <div className="flex flex-col min-h-0">
          <h2 className="text-lg font-quando font-semibold mb-4">Molecule Details</h2>
          
          <div className="space-y-4 flex-1 overflow-auto min-h-0">
            {/* Basic Info */}
            <div className="bg-card border border-border rounded-lg p-4">
              <h3 className="text-sm font-quando font-semibold mb-3">Configuration</h3>
              <div className="space-y-2">
                {moleculeInfo.map((info, index) => (
                  <div key={index} className="flex justify-between text-sm">
                    <span className="text-muted-foreground font-quando">{info.label}:</span>
                    <span className="font-mono font-medium">{info.value}</span>
                  </div>
                ))}
              </div>
            </div>

            {/* Atom Details */}
            {selectedAtom && (
              <div className="bg-card border border-border rounded-lg p-4">
                <h3 className="text-sm font-quando font-semibold mb-3">Selected Atom</h3>
                <div className="space-y-2 text-sm">
                  <div className="flex justify-between">
                    <span className="text-muted-foreground font-quando">Element:</span>
                    <span className="font-mono font-medium">{selectedAtom.elem}</span>
                  </div>
                  <div className="flex justify-between">
                    <span className="text-muted-foreground font-quando">Index:</span>
                    <span className="font-mono font-medium">{selectedAtom.index}</span>
                  </div>
                  <div className="flex justify-between">
                    <span className="text-muted-foreground font-quando">X:</span>
                    <span className="font-mono font-medium">{selectedAtom.x?.toFixed(3)}</span>
                  </div>
                  <div className="flex justify-between">
                    <span className="text-muted-foreground font-quando">Y:</span>
                    <span className="font-mono font-medium">{selectedAtom.y?.toFixed(3)}</span>
                  </div>
                  <div className="flex justify-between">
                    <span className="text-muted-foreground font-quando">Z:</span>
                    <span className="font-mono font-medium">{selectedAtom.z?.toFixed(3)}</span>
                  </div>
                </div>
              </div>
            )}

            {/* Extended Details */}
            {showDetails && (
              <div className="bg-card border border-border rounded-lg p-4">
                <h3 className="text-sm font-quando font-semibold mb-3">Extended Information</h3>
                <div className="space-y-3 text-sm">
                  <div>
                    <span className="text-muted-foreground font-quando block mb-1">Description:</span>
                    <p className="text-muted-foreground">
                      This molecule has been prepared for quantum chemistry calculations. 
                      The 3D structure shows the optimized geometry that will be used 
                      for computational analysis.
                    </p>
                  </div>
                  
                  {molecule.atoms && molecule.atoms.length > 0 && (
                    <div>
                      <span className="text-muted-foreground font-quando block mb-1">Atom List:</span>
                      <div className="grid grid-cols-2 gap-1 text-xs font-mono">
                        {molecule.atoms.map((atom: any, index: number) => (
                          <div key={index} className="flex justify-between">
                            <span>{atom.symbol}</span>
                            <span className="text-muted-foreground">
                              ({atom.x.toFixed(2)}, {atom.y.toFixed(2)}, {atom.z.toFixed(2)})
                            </span>
                          </div>
                        ))}
                      </div>
                    </div>
                  )}
                </div>
              </div>
            )}

            {/* Actions */}
            <div className="bg-card border border-border rounded-lg p-4">
              <h3 className="text-sm font-quando font-semibold mb-3">Actions</h3>
              <div className="space-y-2">
                <button
                  onClick={() => {
                    // TODO: Implement download functionality
                    console.log("Download molecule data");
                  }}
                  className="w-full px-3 py-2 text-sm border border-border rounded hover:bg-accent transition font-quando flex items-center gap-2"
                >
                  <Download className="w-4 h-4" />
                  Download Structure
                </button>
                <button
                  onClick={() => {
                    // TODO: Implement share functionality
                    console.log("Share molecule");
                  }}
                  className="w-full px-3 py-2 text-sm border border-border rounded hover:bg-accent transition font-quando flex items-center gap-2"
                >
                  <Share2 className="w-4 h-4" />
                  Share Molecule
                </button>
              </div>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}
