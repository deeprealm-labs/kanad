"use client";

import { useState, useEffect } from "react";
import { ZoomIn, ZoomOut, RotateCcw } from "lucide-react";

// Declare global SmilesDrawer
declare global {
  interface Window {
    SmilesDrawer: any;
  }
}

interface LewisStructureViewProps {
  molecule: any;
  onExecute: () => void;
}

export default function LewisStructureView({
  molecule,
  onExecute,
}: LewisStructureViewProps) {
  const [zoom, setZoom] = useState(1);
  const [isReady, setIsReady] = useState(false);

  // Check if SmilesDrawer is loaded
  useEffect(() => {
    const checkLibrary = () => {
      if (typeof window !== "undefined" && window.SmilesDrawer) {
        setIsReady(true);
      } else {
        setTimeout(checkLibrary, 100);
      }
    };
    checkLibrary();
  }, []);

  // Draw molecule when ready - using SVG instead
  useEffect(() => {
    // Handle case where molecule is created from atoms without SMILES
    if (!molecule.smiles && molecule.atoms && molecule.atoms.length > 0) {
      const svg = document.getElementById("molecule-svg");
      if (!svg) return;

      // Display atoms in a grid layout
      const atomList = molecule.atoms.map((atom: any) =>
        `<div style="display: inline-block; margin: 10px; padding: 15px 20px; background: ${atom.color}; border: 2px solid hsl(var(--border)); border-radius: 50%; font-size: 24px; font-weight: bold; color: #1a1a1a; min-width: 60px; text-align: center;">${atom.symbol}</div>`
      ).join('');

      svg.innerHTML = `
        <div style="display: flex; align-items: center; justify-content: center; height: 100%; flex-direction: column; font-family: var(--font-quando); padding: 40px;">
          <div style="font-size: 24px; font-weight: bold; margin-bottom: 30px;">
            Molecular Structure
          </div>
          <div style="display: flex; flex-wrap: wrap; justify-content: center; gap: 10px; max-width: 600px; margin-bottom: 20px;">
            ${atomList}
          </div>
          <div style="font-size: 14px; color: var(--muted-foreground); margin-top: 20px;">
            ${molecule.atoms.length} atom${molecule.atoms.length !== 1 ? 's' : ''} | Charge: ${molecule.charge} | Multiplicity: ${molecule.multiplicity}
          </div>
          <div style="font-size: 12px; color: var(--muted-foreground); margin-top: 10px;">
            Basis Set: ${molecule.basis || 'sto-3g'}
          </div>
        </div>
      `;
      return;
    }

    if (isReady && molecule.smiles) {
      console.log("Drawing molecule with SMILES:", molecule.smiles);

      // Clean the SMILES string
      const cleanSmiles = window.SmilesDrawer.clean(molecule.smiles);
      console.log("Cleaned SMILES:", cleanSmiles);

      // Wait for SVG container to be in DOM
      setTimeout(() => {
        const svg = document.getElementById("molecule-svg");
        if (!svg) {
          console.error("SVG not found");
          return;
        }

        // Special handling for very simple molecules that library can't handle
        const simpleMolecules: Record<string, string> = {
          'C': 'CH₄ (Methane)',
          'c': 'CH₄ (Methane)',
          'N': 'NH₃ (Ammonia)',
          'O': 'H₂O (Water)',
          'Cl': 'HCl (Hydrogen Chloride)',
          '[H]Cl': 'HCl (Hydrogen Chloride)',
          'HCl': 'HCl (Hydrogen Chloride)',
        };

        // Check if it's a very simple molecule or too short (likely to fail)
        if (simpleMolecules[cleanSmiles] || cleanSmiles.length <= 2) {
          const displayName = simpleMolecules[cleanSmiles] || `Single atom: ${cleanSmiles}`;
          svg.innerHTML = `
            <div style="display: flex; align-items: center; justify-content: center; height: 100%; flex-direction: column; font-family: var(--font-quando);">
              <div style="font-size: 48px; font-weight: bold; margin-bottom: 20px;">
                ${displayName}
              </div>
              <div style="font-size: 18px; color: var(--muted-foreground);">
                SMILES: ${cleanSmiles}
              </div>
              <div style="font-size: 14px; color: var(--muted-foreground); margin-top: 10px;">
                (Structural formula for simple molecules)
              </div>
            </div>
          `;
          return;
        }

        console.log("SVG found");

        // Custom options with website aesthetics
        const options = {
          width: 800,
          height: 600,
          bondThickness: 1.5,
          bondLength: 20,
          shortBondLength: 0.85,
          bondSpacing: 0.25,
          atomVisualization: 'default',
          isomeric: true,
          debug: false,
          terminalCarbons: false,
          explicitHydrogens: false,
          overlapSensitivity: 0.42,
          overlapResolutionIterations: 3,
          compactDrawing: true,
          fontFamily: 'var(--font-quando), Arial, Helvetica, sans-serif',
          fontSizeLarge: 11,
          fontSizeSmall: 9,
          padding: 40,
          themes: {
            light: {
              C: '#1a1a1a',
              O: '#ef4444',
              N: '#3b82f6',
              F: '#10b981',
              CL: '#10b981',
              BR: '#a16207',
              I: '#6b21a8',
              P: '#ea580c',
              S: '#eab308',
              B: '#f97316',
              SI: '#78716c',
              H: '#525252',
              BACKGROUND: '#ffffff'
            }
          }
        };

        const drawer = new window.SmilesDrawer.SvgDrawer(options);

        window.SmilesDrawer.parse(
          cleanSmiles,
          (tree: any) => {
            console.log("Parse successful, tree:", tree);
            try {
              // Create a new SVG element for the library to use
              const svgElement = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
              svgElement.setAttribute('xmlns', 'http://www.w3.org/2000/svg');

              // Draw using SVG - pass the actual element, not string ID
              drawer.draw(tree, svgElement, "light");

              // Clear container and append the drawn SVG
              svg.innerHTML = "";
              svg.appendChild(svgElement);

              console.log("Draw successful");
            } catch (err) {
              console.error("Draw error:", err);
              console.error("Tree that caused error:", tree);

              // Fallback display
              svg.innerHTML = `
                <div style="display: flex; align-items: center; justify-content: center; height: 100%; flex-direction: column; font-family: var(--font-quando);">
                  <div style="font-size: 36px; font-weight: bold; margin-bottom: 20px;">
                    ${cleanSmiles}
                  </div>
                  <div style="font-size: 14px; color: var(--muted-foreground);">
                    (Structure visualization not available)
                  </div>
                </div>
              `;
            }
          },
          (err: any) => {
            console.error("Parse error:", err);
            // Show error in UI
            svg.innerHTML = `<div style="padding: 20px; color: red; font-family: var(--font-quando);">Error parsing SMILES: ${err}</div>`;
          }
        );
      }, 300);
    }
  }, [isReady, molecule]);

  const handleZoomIn = () => setZoom(Math.min(zoom + 0.2, 3));
  const handleZoomOut = () => setZoom(Math.max(zoom - 0.2, 0.5));
  const handleReset = () => setZoom(1);

  return (
    <div className="h-full flex items-center justify-center relative bg-background overflow-hidden">
      {/* Zoom Controls - with high z-index to stay on top */}
      <div className="absolute left-8 top-1/2 -translate-y-1/2 flex flex-col gap-2 z-50">
        <button
          onClick={handleReset}
          className="w-10 h-10 bg-card border border-border rounded-full flex items-center justify-center hover:bg-accent transition shadow-lg"
          title="Reset zoom"
        >
          <RotateCcw className="w-5 h-5" />
        </button>
        <button
          onClick={handleZoomIn}
          className="w-10 h-10 bg-card border border-border rounded-full flex items-center justify-center hover:bg-accent transition shadow-lg"
          title="Zoom in"
        >
          <ZoomIn className="w-5 h-5" />
        </button>
        <button
          onClick={handleZoomOut}
          className="w-10 h-10 bg-card border border-border rounded-full flex items-center justify-center hover:bg-accent transition shadow-lg"
          title="Zoom out"
        >
          <ZoomOut className="w-5 h-5" />
        </button>
      </div>

      {/* Molecular Structure SVG - with lower z-index */}
      <div className="flex items-center justify-center z-10" style={{ transform: `scale(${zoom})`, transition: "transform 0.2s", transformOrigin: "center" }}>
        <div
          id="molecule-svg"
          style={{
            width: "800px",
            height: "600px",
            background: "transparent"
          }}
        />
      </div>

      {/* Execute Button - with high z-index */}
      <button
        onClick={onExecute}
        className="absolute bottom-8 right-8 px-6 py-3 bg-background border-2 border-border rounded-lg hover:border-brand-orange transition font-quando flex items-center gap-2 z-50 shadow-lg"
      >
        <span>execute now</span>
      </button>
    </div>
  );
}
