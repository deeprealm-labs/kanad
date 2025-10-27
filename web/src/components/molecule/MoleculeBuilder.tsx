"use client";

import { useState } from "react";
import { Microscope, Atom, Library, FileUp, X } from "lucide-react";
import { moleculeLibrary } from "@/data/molecule-library";

type InputMethod = "smiles" | "atoms" | "library" | "xyz";

interface MoleculeBuilderProps {
  onClose: () => void;
  onComplete?: () => void;
}

interface AtomCoordinate {
  element: string;
  x: number;
  y: number;
  z: number;
}

export default function MoleculeBuilder({
  onClose,
  onComplete,
}: MoleculeBuilderProps) {
  const [method, setMethod] = useState<InputMethod>("smiles");
  const [smiles, setSmiles] = useState("");
  const [basis, setBasis] = useState("sto-3g");
  const [charge, setCharge] = useState(0);
  const [multiplicity, setMultiplicity] = useState(1);
  const [atoms, setAtoms] = useState<AtomCoordinate[]>([]);
  const [xyzFile, setXyzFile] = useState<File | null>(null);
  const [isDragging, setIsDragging] = useState(false);

  const handleCreate = async () => {
    // API call will be implemented
    console.log("Creating molecule:", { method, smiles, basis, charge, multiplicity, atoms });
  };

  const handleFileUpload = (file: File) => {
    if (file.name.endsWith(".xyz")) {
      setXyzFile(file);
      // Parse XYZ file
      const reader = new FileReader();
      reader.onload = (e) => {
        const content = e.target?.result as string;
        parseXYZFile(content);
      };
      reader.readAsText(file);
    }
  };

  const parseXYZFile = (content: string) => {
    // Simple XYZ parser
    const lines = content.trim().split("\n");
    if (lines.length < 3) return;

    const numAtoms = parseInt(lines[0]);
    // line 1 is comment
    const atomLines = lines.slice(2, 2 + numAtoms);

    const parsedAtoms: AtomCoordinate[] = atomLines.map((line) => {
      const parts = line.trim().split(/\s+/);
      return {
        element: parts[0],
        x: parseFloat(parts[1]),
        y: parseFloat(parts[2]),
        z: parseFloat(parts[3]),
      };
    });

    setAtoms(parsedAtoms);
  };

  const handleDragOver = (e: React.DragEvent) => {
    e.preventDefault();
    setIsDragging(true);
  };

  const handleDragLeave = (e: React.DragEvent) => {
    e.preventDefault();
    setIsDragging(false);
  };

  const handleDrop = (e: React.DragEvent) => {
    e.preventDefault();
    setIsDragging(false);
    const file = e.dataTransfer.files[0];
    if (file) handleFileUpload(file);
  };

  return (
    <div className="h-full p-8 overflow-auto bg-background">
      <div className="max-w-4xl mx-auto">
        {/* Header */}
        <div className="flex items-center justify-between mb-8">
          <h1 className="text-3xl font-quando font-bold">Create Molecule</h1>
          <button
            onClick={onClose}
            className="text-muted-foreground hover:text-foreground transition"
          >
            <X className="w-6 h-6" />
          </button>
        </div>

        {/* Method Selection */}
        <div className="mb-8">
          <label className="block text-sm font-quando font-medium mb-3">
            Input Method
          </label>
          <div className="grid grid-cols-4 gap-4">
            <button
              onClick={() => setMethod("smiles")}
              className={`p-4 border-2 rounded-lg transition font-quando ${
                method === "smiles"
                  ? "border-brand-orange bg-orange-50 dark:bg-orange-950"
                  : "border-border hover:border-muted-foreground"
              }`}
            >
              <div className="mb-2 flex justify-center">
                <Microscope className="w-8 h-8" />
              </div>
              <div className="text-sm">From SMILES</div>
            </button>
            <button
              onClick={() => setMethod("atoms")}
              className={`p-4 border-2 rounded-lg transition font-quando ${
                method === "atoms"
                  ? "border-brand-orange bg-orange-50 dark:bg-orange-950"
                  : "border-border hover:border-muted-foreground"
              }`}
            >
              <div className="mb-2 flex justify-center">
                <Atom className="w-8 h-8" />
              </div>
              <div className="text-sm">From Atoms</div>
            </button>
            <button
              onClick={() => setMethod("library")}
              className={`p-4 border-2 rounded-lg transition font-quando ${
                method === "library"
                  ? "border-brand-orange bg-orange-50 dark:bg-orange-950"
                  : "border-border hover:border-muted-foreground"
              }`}
            >
              <div className="mb-2 flex justify-center">
                <Library className="w-8 h-8" />
              </div>
              <div className="text-sm">From Library</div>
            </button>
            <button
              onClick={() => setMethod("xyz")}
              className={`p-4 border-2 rounded-lg transition font-quando ${
                method === "xyz"
                  ? "border-brand-orange bg-orange-50 dark:bg-orange-950"
                  : "border-border hover:border-muted-foreground"
              }`}
            >
              <div className="mb-2 flex justify-center">
                <FileUp className="w-8 h-8" />
              </div>
              <div className="text-sm">Upload XYZ</div>
            </button>
          </div>
        </div>

        {/* SMILES Input */}
        {method === "smiles" && (
          <div className="space-y-6 bg-card border border-border rounded-lg p-6">
            <div>
              <label className="block text-sm font-quando font-medium mb-2">
                SMILES String
              </label>
              <input
                type="text"
                value={smiles}
                onChange={(e) => setSmiles(e.target.value)}
                placeholder="e.g., CCO (Ethanol), c1ccccc1 (Benzene)"
                className="w-full px-4 py-3 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-mono"
              />
              <p className="text-xs text-muted-foreground mt-1 font-quando">
                Enter a valid SMILES notation for your molecule
              </p>
            </div>

            <div className="grid grid-cols-3 gap-4">
              <div>
                <label className="block text-sm font-quando font-medium mb-2">
                  Basis Set
                </label>
                <select
                  value={basis}
                  onChange={(e) => setBasis(e.target.value)}
                  className="w-full px-4 py-3 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                >
                  <option value="sto-3g">STO-3G</option>
                  <option value="6-31g">6-31G</option>
                  <option value="6-31g*" disabled>6-31G* (not implemented)</option>
                  <option value="cc-pvdz" disabled>cc-pVDZ (not implemented)</option>
                  <option value="cc-pvtz" disabled>cc-pVTZ (not implemented)</option>
                </select>
              </div>

              <div>
                <label className="block text-sm font-quando font-medium mb-2">
                  Charge
                </label>
                <input
                  type="number"
                  value={charge}
                  onChange={(e) => setCharge(parseInt(e.target.value))}
                  className="w-full px-4 py-3 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                />
              </div>

              <div>
                <label className="block text-sm font-quando font-medium mb-2">
                  Multiplicity
                </label>
                <input
                  type="number"
                  value={multiplicity}
                  onChange={(e) => setMultiplicity(parseInt(e.target.value))}
                  min={1}
                  className="w-full px-4 py-3 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                />
              </div>
            </div>

            {/* Preview Area */}
            {smiles && (
              <div className="border-t border-border pt-6">
                <h3 className="text-sm font-quando font-medium mb-3">
                  Molecule Preview
                </h3>
                <div className="bg-muted border border-border rounded-lg p-8 text-center">
                  <p className="text-muted-foreground font-quando">
                    Molecular structure visualization will appear here
                  </p>
                  <p className="text-xs text-muted-foreground mt-2 font-mono">
                    SMILES: {smiles}
                  </p>
                </div>
              </div>
            )}
          </div>
        )}

        {/* Atoms Input */}
        {method === "atoms" && (
          <div className="space-y-6 bg-card border border-border rounded-lg p-6">
            <div>
              <label className="block text-sm font-quando font-medium mb-2">
                Atomic Coordinates
              </label>
              <div className="space-y-3">
                {atoms.map((atom, index) => (
                  <div key={index} className="flex gap-2 items-center">
                    <select
                      value={atom.element}
                      onChange={(e) => {
                        const newAtoms = [...atoms];
                        newAtoms[index].element = e.target.value;
                        setAtoms(newAtoms);
                      }}
                      className="w-20 px-2 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando text-sm"
                    >
                      <option value="H">H</option>
                      <option value="C">C</option>
                      <option value="N">N</option>
                      <option value="O">O</option>
                      <option value="F">F</option>
                      <option value="S">S</option>
                      <option value="Cl">Cl</option>
                      <option value="Li">Li</option>
                      <option value="Fe">Fe</option>
                      <option value="Ti">Ti</option>
                    </select>
                    <input
                      type="number"
                      step="0.01"
                      value={atom.x}
                      onChange={(e) => {
                        const newAtoms = [...atoms];
                        newAtoms[index].x = parseFloat(e.target.value) || 0;
                        setAtoms(newAtoms);
                      }}
                      placeholder="X"
                      className="flex-1 px-3 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-mono text-sm"
                    />
                    <input
                      type="number"
                      step="0.01"
                      value={atom.y}
                      onChange={(e) => {
                        const newAtoms = [...atoms];
                        newAtoms[index].y = parseFloat(e.target.value) || 0;
                        setAtoms(newAtoms);
                      }}
                      placeholder="Y"
                      className="flex-1 px-3 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-mono text-sm"
                    />
                    <input
                      type="number"
                      step="0.01"
                      value={atom.z}
                      onChange={(e) => {
                        const newAtoms = [...atoms];
                        newAtoms[index].z = parseFloat(e.target.value) || 0;
                        setAtoms(newAtoms);
                      }}
                      placeholder="Z"
                      className="flex-1 px-3 py-2 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-mono text-sm"
                    />
                    <button
                      onClick={() => {
                        const newAtoms = atoms.filter((_, i) => i !== index);
                        setAtoms(newAtoms);
                      }}
                      className="p-2 text-muted-foreground hover:text-destructive transition"
                    >
                      <X className="w-4 h-4" />
                    </button>
                  </div>
                ))}
              </div>
              <button
                onClick={() => {
                  setAtoms([...atoms, { element: "H", x: 0, y: 0, z: 0 }]);
                }}
                className="mt-3 w-full px-4 py-2 border border-border rounded-md hover:bg-accent transition font-quando text-sm"
              >
                + Add Atom
              </button>
            </div>

            <div className="grid grid-cols-3 gap-4">
              <div>
                <label className="block text-sm font-quando font-medium mb-2">
                  Basis Set
                </label>
                <select
                  value={basis}
                  onChange={(e) => setBasis(e.target.value)}
                  className="w-full px-4 py-3 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                >
                  <option value="sto-3g">STO-3G</option>
                  <option value="6-31g">6-31G</option>
                  <option value="6-31g*" disabled>6-31G* (not implemented)</option>
                  <option value="cc-pvdz" disabled>cc-pVDZ (not implemented)</option>
                  <option value="cc-pvtz" disabled>cc-pVTZ (not implemented)</option>
                </select>
              </div>

              <div>
                <label className="block text-sm font-quando font-medium mb-2">
                  Charge
                </label>
                <input
                  type="number"
                  value={charge}
                  onChange={(e) => setCharge(parseInt(e.target.value))}
                  className="w-full px-4 py-3 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                />
              </div>

              <div>
                <label className="block text-sm font-quando font-medium mb-2">
                  Multiplicity
                </label>
                <input
                  type="number"
                  value={multiplicity}
                  onChange={(e) => setMultiplicity(parseInt(e.target.value))}
                  min={1}
                  className="w-full px-4 py-3 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                />
              </div>
            </div>

            {atoms.length > 0 && (
              <div className="border-t border-border pt-6">
                <h3 className="text-sm font-quando font-medium mb-3">
                  Molecule Summary
                </h3>
                <div className="bg-muted border border-border rounded-lg p-4">
                  <p className="text-sm text-muted-foreground font-quando">
                    Total atoms: {atoms.length}
                  </p>
                  <p className="text-xs text-muted-foreground mt-2 font-mono">
                    {atoms.map((a) => a.element).join(" ")}
                  </p>
                </div>
              </div>
            )}
          </div>
        )}

        {/* Library */}
        {method === "library" && (
          <div className="bg-card border border-border rounded-lg p-6 space-y-6">
            {moleculeLibrary.map((category) => (
              <div key={category.name}>
                <h3 className="text-lg font-quando font-semibold mb-2">
                  {category.name}
                </h3>
                <p className="text-sm text-muted-foreground mb-4 font-quando">
                  {category.description}
                </p>
                <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-3">
                  {category.molecules.map((molecule) => (
                    <button
                      key={molecule.id}
                      onClick={() => {
                        setSmiles(molecule.smiles);
                        setBasis(molecule.basis);
                        setCharge(molecule.charge);
                        setMultiplicity(molecule.multiplicity);
                        setMethod("smiles"); // Switch to SMILES view to show selection
                      }}
                      className="p-4 border border-border rounded-lg hover:border-brand-orange hover:bg-accent transition text-left"
                    >
                      <div className="font-quando font-medium mb-1">
                        {molecule.name}
                      </div>
                      <div className="text-xs text-muted-foreground mb-2 font-mono">
                        {molecule.formula}
                      </div>
                      <div className="text-xs text-muted-foreground font-quando">
                        {molecule.description}
                      </div>
                    </button>
                  ))}
                </div>
              </div>
            ))}
          </div>
        )}

        {/* XYZ Upload */}
        {method === "xyz" && (
          <div className="bg-card border border-border rounded-lg p-6 space-y-6">
            <div
              onDragOver={handleDragOver}
              onDragLeave={handleDragLeave}
              onDrop={handleDrop}
              className={`border-2 border-dashed rounded-lg p-12 text-center transition ${
                isDragging
                  ? "border-brand-orange bg-orange-50 dark:bg-orange-950"
                  : "border-border"
              }`}
            >
              <input
                type="file"
                accept=".xyz"
                onChange={(e) => {
                  const file = e.target.files?.[0];
                  if (file) handleFileUpload(file);
                }}
                className="hidden"
                id="xyz-upload"
              />
              <label htmlFor="xyz-upload" className="cursor-pointer">
                <div className="mb-4 flex justify-center">
                  <FileUp className="w-16 h-16 text-muted-foreground" />
                </div>
                <p className="text-muted-foreground font-quando mb-2">
                  {xyzFile
                    ? `Selected: ${xyzFile.name}`
                    : "Drag and drop your XYZ file here"}
                </p>
                <p className="text-xs text-muted-foreground font-quando">
                  or click to browse
                </p>
              </label>
            </div>

            {atoms.length > 0 && (
              <div className="space-y-4">
                <div>
                  <label className="block text-sm font-quando font-medium mb-2">
                    Parsed Atoms ({atoms.length})
                  </label>
                  <div className="bg-muted border border-border rounded-lg p-4 max-h-64 overflow-auto">
                    <div className="space-y-1 font-mono text-xs">
                      {atoms.map((atom, index) => (
                        <div key={index} className="text-muted-foreground">
                          {atom.element.padEnd(3)} {atom.x.toFixed(4).padStart(10)}{" "}
                          {atom.y.toFixed(4).padStart(10)}{" "}
                          {atom.z.toFixed(4).padStart(10)}
                        </div>
                      ))}
                    </div>
                  </div>
                </div>

                <div className="grid grid-cols-3 gap-4">
                  <div>
                    <label className="block text-sm font-quando font-medium mb-2">
                      Basis Set
                    </label>
                    <select
                      value={basis}
                      onChange={(e) => setBasis(e.target.value)}
                      className="w-full px-4 py-3 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                    >
                      <option value="sto-3g">STO-3G</option>
                      <option value="6-31g">6-31G</option>
                      <option value="6-31g*">6-31G*</option>
                      <option value="cc-pvdz">cc-pVDZ</option>
                      <option value="cc-pvtz">cc-pVTZ</option>
                    </select>
                  </div>

                  <div>
                    <label className="block text-sm font-quando font-medium mb-2">
                      Charge
                    </label>
                    <input
                      type="number"
                      value={charge}
                      onChange={(e) => setCharge(parseInt(e.target.value))}
                      className="w-full px-4 py-3 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                    />
                  </div>

                  <div>
                    <label className="block text-sm font-quando font-medium mb-2">
                      Multiplicity
                    </label>
                    <input
                      type="number"
                      value={multiplicity}
                      onChange={(e) => setMultiplicity(parseInt(e.target.value))}
                      min={1}
                      className="w-full px-4 py-3 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                    />
                  </div>
                </div>
              </div>
            )}
          </div>
        )}

        {/* Action Buttons */}
        <div className="flex gap-4 mt-8">
          <button
            onClick={onClose}
            className="flex-1 px-6 py-3 border border-border rounded-lg hover:bg-accent transition font-quando"
          >
            Cancel
          </button>
          <button
            onClick={() => {
              handleCreate();
              if (onComplete) onComplete();
            }}
            disabled={
              (method === "smiles" && !smiles) ||
              (method === "atoms" && atoms.length === 0)
            }
            className="flex-1 px-6 py-3 bg-brand-orange text-white rounded-lg hover:bg-brand-orange-dark transition font-quando disabled:opacity-50 disabled:cursor-not-allowed"
          >
            Continue to Configuration →
          </button>
        </div>
      </div>
    </div>
  );
}
