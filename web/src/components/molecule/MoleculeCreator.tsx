"use client";

import { useState } from "react";
import { Search, X, Plus } from "lucide-react";
import { moleculeLibrary } from "@/data/molecule-library";

interface AtomData {
  symbol: string;
  name: string;
  atomicNumber: number;
  massNumber: number;
  color: string;
}

interface DroppedAtom {
  id: string;
  symbol: string;
  atomicNumber: number;
  massNumber: number;
  color: string;
  x: number;
  y: number;
  z: number;
}

interface MoleculeCreatorProps {
  onComplete: (molecule: any) => void;
  onCancel: () => void;
}

// Periodic table with colors based on element groups
const periodicTable: AtomData[] = [
  { symbol: "H", name: "Hydrogen", atomicNumber: 1, massNumber: 1, color: "#e8f4f8" },
  { symbol: "He", name: "Helium", atomicNumber: 2, massNumber: 4, color: "#d0e8f2" },
  { symbol: "Li", name: "Lithium", atomicNumber: 3, massNumber: 7, color: "#f4d0e8" },
  { symbol: "Be", name: "Beryllium", atomicNumber: 4, massNumber: 9, color: "#f4d0e8" },
  { symbol: "B", name: "Boron", atomicNumber: 5, massNumber: 11, color: "#e8f4d0" },
  { symbol: "C", name: "Carbon", atomicNumber: 6, massNumber: 12, color: "#e8f4f8" },
  { symbol: "N", name: "Nitrogen", atomicNumber: 7, massNumber: 14, color: "#e8f4f8" },
  { symbol: "O", name: "Oxygen", atomicNumber: 8, massNumber: 16, color: "#e8f4f8" },
  { symbol: "F", name: "Fluorine", atomicNumber: 9, massNumber: 19, color: "#f8f4d0" },
  { symbol: "Ne", name: "Neon", atomicNumber: 10, massNumber: 20, color: "#d0e8f2" },
  { symbol: "Na", name: "Sodium", atomicNumber: 11, massNumber: 23, color: "#f4d0e8" },
  { symbol: "Mg", name: "Magnesium", atomicNumber: 12, massNumber: 24, color: "#f4d0e8" },
  { symbol: "Al", name: "Aluminum", atomicNumber: 13, massNumber: 27, color: "#d0e8e8" },
  { symbol: "Si", name: "Silicon", atomicNumber: 14, massNumber: 28, color: "#e8f4d0" },
  { symbol: "P", name: "Phosphorus", atomicNumber: 15, massNumber: 31, color: "#e8f4f8" },
  { symbol: "S", name: "Sulfur", atomicNumber: 16, massNumber: 32, color: "#e8f4f8" },
  { symbol: "Cl", name: "Chlorine", atomicNumber: 17, massNumber: 35, color: "#f8f4d0" },
  { symbol: "Ar", name: "Argon", atomicNumber: 18, massNumber: 40, color: "#d0e8f2" },
  { symbol: "Fe", name: "Iron", atomicNumber: 26, massNumber: 56, color: "#f0d8e8" },
  { symbol: "Cu", name: "Copper", atomicNumber: 29, massNumber: 64, color: "#f0d8e8" },
  { symbol: "Zn", name: "Zinc", atomicNumber: 30, massNumber: 65, color: "#f0d8e8" },
  { symbol: "Ti", name: "Titanium", atomicNumber: 22, massNumber: 48, color: "#f0d8e8" },
];

export default function MoleculeCreator({
  onComplete,
  onCancel,
}: MoleculeCreatorProps) {
  const [searchTerm, setSearchTerm] = useState("");
  const [smilesInput, setSmilesInput] = useState("");
  const [smiles, setSmiles] = useState("");
  const [droppedAtoms, setDroppedAtoms] = useState<DroppedAtom[]>([]);
  const [basis, setBasis] = useState("sto-3g");
  const [charge, setCharge] = useState(0);
  const [multiplicity, setMultiplicity] = useState(1);
  const [viewMode, setViewMode] = useState<"atoms" | "molecules">("atoms");

  const filteredAtoms = periodicTable.filter(
    (atom) =>
      atom.symbol.toLowerCase().includes(searchTerm.toLowerCase()) ||
      atom.name.toLowerCase().includes(searchTerm.toLowerCase())
  );

  const filteredMolecules = moleculeLibrary.flatMap(cat => cat.molecules).filter(
    (mol) =>
      mol.name.toLowerCase().includes(searchTerm.toLowerCase()) ||
      mol.formula.toLowerCase().includes(searchTerm.toLowerCase())
  );

  const handleDragStart = (e: React.DragEvent, atom: AtomData) => {
    e.dataTransfer.setData("atom", JSON.stringify(atom));
  };

  const handleDrop = (e: React.DragEvent) => {
    e.preventDefault();
    const atomData = JSON.parse(e.dataTransfer.getData("atom")) as AtomData;
    const newAtom: DroppedAtom = {
      id: Math.random().toString(36).substr(2, 9),
      symbol: atomData.symbol,
      atomicNumber: atomData.atomicNumber,
      massNumber: atomData.massNumber,
      color: atomData.color,
      x: 0,
      y: 0,
      z: 0,
    };
    setDroppedAtoms([...droppedAtoms, newAtom]);
  };

  const handleDragOver = (e: React.DragEvent) => {
    e.preventDefault();
  };

  const removeAtom = (id: string) => {
    setDroppedAtoms(droppedAtoms.filter((a) => a.id !== id));
  };

  const handleAddSmiles = () => {
    if (smilesInput.trim()) {
      setSmiles(smilesInput);
      setDroppedAtoms([]); // Clear atoms when SMILES is added
    }
  };

  const handleLoadMolecule = (molecule: any) => {
    setSmilesInput(molecule.smiles);
    setSmiles(molecule.smiles);
    setBasis(molecule.basis);
    setCharge(molecule.charge);
    setMultiplicity(molecule.multiplicity);
    setDroppedAtoms([]); // Clear atoms when molecule is loaded
  };

  const handleSubmit = () => {
    if (smiles || droppedAtoms.length > 0) {
      onComplete({
        smiles,
        atoms: droppedAtoms,
        basis,
        charge,
        multiplicity,
      });
    }
  };

  return (
    <div className="h-full flex flex-col bg-background overflow-hidden">
      {/* Top Bar */}
      <div className="flex items-center justify-between px-6 py-4 border-b border-border">
        <h1 className="text-xl font-quando font-bold">Create Molecule</h1>
        <button
          onClick={onCancel}
          className="p-2 hover:bg-accent rounded-md transition"
        >
          <X className="w-5 h-5" />
        </button>
      </div>

      {/* Main Content - Fixed Height, No Scrolling */}
      <div className="flex-1 grid grid-cols-3 gap-4 p-4 overflow-hidden min-h-0">
        {/* Left: Library & Search */}
        <div className="flex flex-col h-full min-h-0">
          <div className="flex items-center justify-between mb-2">
            <h2 className="text-base font-quando font-semibold">Library</h2>
            <div className="flex gap-1">
              <button
                onClick={() => setViewMode("atoms")}
                className={`px-3 py-1 text-xs rounded ${
                  viewMode === "atoms"
                    ? "bg-brand-orange text-white"
                    : "bg-card border border-border"
                }`}
              >
                Atoms
              </button>
              <button
                onClick={() => setViewMode("molecules")}
                className={`px-3 py-1 text-xs rounded ${
                  viewMode === "molecules"
                    ? "bg-brand-orange text-white"
                    : "bg-card border border-border"
                }`}
              >
                Molecules
              </button>
            </div>
          </div>

          {/* Search */}
          <div className="relative mb-2">
            <Search className="absolute left-2 top-2 w-4 h-4 text-muted-foreground" />
            <input
              type="text"
              value={searchTerm}
              onChange={(e) => setSearchTerm(e.target.value)}
              placeholder={`Search ${viewMode}...`}
              className="w-full pl-8 pr-3 py-1.5 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando text-xs"
            />
          </div>

          {/* Library Content */}
          <div className="flex-1 overflow-auto border border-border rounded-lg p-2 min-h-0">
            {viewMode === "atoms" ? (
              <div className="grid grid-cols-3 gap-1.5">
                {filteredAtoms.map((atom) => (
                  <div
                    key={atom.symbol}
                    draggable
                    onDragStart={(e) => handleDragStart(e, atom)}
                    style={{ backgroundColor: atom.color }}
                    className="border border-border rounded p-1.5 cursor-move hover:border-brand-orange transition"
                  >
                    <div className="text-center">
                      <div className="text-[9px] text-muted-foreground leading-none">
                        {atom.atomicNumber}
                      </div>
                      <div className="text-lg font-bold font-quando leading-none my-0.5">
                        {atom.symbol}
                      </div>
                      <div className="text-[9px] text-muted-foreground leading-none">
                        {atom.massNumber}
                      </div>
                    </div>
                  </div>
                ))}
              </div>
            ) : (
              <div className="space-y-1.5">
                {filteredMolecules.map((mol) => (
                  <button
                    key={mol.id}
                    onClick={() => handleLoadMolecule(mol)}
                    className="w-full p-2 border border-border rounded hover:border-brand-orange hover:bg-accent transition text-left"
                  >
                    <div className="font-quando font-medium text-sm">
                      {mol.name}
                    </div>
                    <div className="text-[10px] text-muted-foreground">
                      {mol.formula}
                    </div>
                  </button>
                ))}
              </div>
            )}
          </div>
        </div>

        {/* Center: Drop Zone */}
        <div className="flex flex-col h-full min-h-0">
          <h2 className="text-base font-quando font-semibold mb-2">
            Build Molecule
          </h2>

          {/* Drop Zone */}
          <div
            onDrop={handleDrop}
            onDragOver={handleDragOver}
            className="flex-1 border-2 border-dashed border-border rounded-lg p-3 overflow-auto min-h-0 hover:border-brand-orange transition"
          >
            {!smiles && droppedAtoms.length === 0 ? (
              <div className="h-full flex items-center justify-center text-center text-muted-foreground font-quando">
                <div>
                  <p className="mb-1 text-sm">Drop atoms here</p>
                  <p className="text-xs">or use SMILES below</p>
                </div>
              </div>
            ) : smiles ? (
              <div className="h-full flex items-center justify-center">
                <div className="text-center">
                  <div className="text-xs text-muted-foreground mb-1">SMILES</div>
                  <div className="text-2xl font-mono font-bold">{smiles}</div>
                </div>
              </div>
            ) : (
              <div className="space-y-1.5">
                {droppedAtoms.map((atom) => (
                  <div
                    key={atom.id}
                    className="flex items-center gap-2 p-2 bg-card border border-border rounded"
                  >
                    <div
                      style={{ backgroundColor: atom.color }}
                      className="w-10 h-10 border border-border rounded flex-shrink-0 flex flex-col items-center justify-center"
                    >
                      <div className="text-[8px] text-muted-foreground leading-none">
                        {atom.atomicNumber}
                      </div>
                      <div className="text-sm font-bold leading-none my-0.5">
                        {atom.symbol}
                      </div>
                      <div className="text-[8px] text-muted-foreground leading-none">
                        {atom.massNumber}
                      </div>
                    </div>
                    <div className="flex gap-1 flex-1">
                      <input
                        type="number"
                        step="0.01"
                        value={atom.x}
                        onChange={(e) => {
                          const updated = droppedAtoms.map((a) =>
                            a.id === atom.id
                              ? { ...a, x: parseFloat(e.target.value) || 0 }
                              : a
                          );
                          setDroppedAtoms(updated);
                        }}
                        placeholder="X"
                        className="w-full px-2 py-1 border border-input bg-background rounded text-xs font-mono"
                      />
                      <input
                        type="number"
                        step="0.01"
                        value={atom.y}
                        onChange={(e) => {
                          const updated = droppedAtoms.map((a) =>
                            a.id === atom.id
                              ? { ...a, y: parseFloat(e.target.value) || 0 }
                              : a
                          );
                          setDroppedAtoms(updated);
                        }}
                        placeholder="Y"
                        className="w-full px-2 py-1 border border-input bg-background rounded text-xs font-mono"
                      />
                      <input
                        type="number"
                        step="0.01"
                        value={atom.z}
                        onChange={(e) => {
                          const updated = droppedAtoms.map((a) =>
                            a.id === atom.id
                              ? { ...a, z: parseFloat(e.target.value) || 0 }
                              : a
                          );
                          setDroppedAtoms(updated);
                        }}
                        placeholder="Z"
                        className="w-full px-2 py-1 border border-input bg-background rounded text-xs font-mono"
                      />
                    </div>
                    <button
                      onClick={() => removeAtom(atom.id)}
                      className="p-1 hover:text-destructive transition flex-shrink-0"
                    >
                      <X className="w-4 h-4" />
                    </button>
                  </div>
                ))}
              </div>
            )}
          </div>

          {/* SMILES Input */}
          <div className="mt-2">
            <label className="block text-xs font-quando font-medium mb-1">
              Or enter SMILES
            </label>
            <div className="flex gap-1">
              <input
                type="text"
                value={smilesInput}
                onChange={(e) => setSmilesInput(e.target.value)}
                placeholder="e.g., CCO, c1ccccc1"
                className="flex-1 px-3 py-1.5 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-mono text-xs"
              />
              <button
                onClick={handleAddSmiles}
                disabled={!smilesInput.trim()}
                className="px-3 py-1.5 bg-brand-orange text-white rounded-md hover:bg-brand-orange-dark transition text-xs font-quando disabled:opacity-50 disabled:cursor-not-allowed flex items-center gap-1"
              >
                <Plus className="w-3 h-3" />
                Add
              </button>
            </div>
          </div>
        </div>

        {/* Right: Configuration */}
        <div className="flex flex-col h-full min-h-0">
          <h2 className="text-base font-quando font-semibold mb-2">
            Configuration
          </h2>

          <div className="space-y-3 flex-1 overflow-auto min-h-0">
            <div>
              <label className="block text-xs font-quando font-medium mb-1">
                Basis Set
              </label>
              <select
                value={basis}
                onChange={(e) => setBasis(e.target.value)}
                className="w-full px-3 py-1.5 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando text-xs"
              >
                <option value="sto-3g">STO-3G</option>
                <option value="6-31g">6-31G</option>
                <option value="6-31g*">6-31G*</option>
                <option value="cc-pvdz">cc-pVDZ</option>
                <option value="cc-pvtz">cc-pVTZ</option>
              </select>
            </div>

            <div>
              <label className="block text-xs font-quando font-medium mb-1">
                Charge
              </label>
              <input
                type="number"
                value={charge}
                onChange={(e) => setCharge(parseInt(e.target.value))}
                className="w-full px-3 py-1.5 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando text-xs"
              />
            </div>

            <div>
              <label className="block text-xs font-quando font-medium mb-1">
                Multiplicity
              </label>
              <input
                type="number"
                value={multiplicity}
                onChange={(e) => setMultiplicity(parseInt(e.target.value))}
                min={1}
                className="w-full px-3 py-1.5 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando text-xs"
              />
            </div>

            <div>
              <label className="block text-xs font-quando font-medium mb-1">
                Upload XYZ File
              </label>
              <input
                type="file"
                accept=".xyz"
                className="w-full px-3 py-1.5 border border-input bg-background rounded-md font-quando text-xs file:mr-2 file:px-3 file:py-1 file:rounded file:border-0 file:bg-brand-orange file:text-white file:text-xs"
              />
            </div>
          </div>

          <div className="mt-3 pt-3 border-t border-border">
            <button
              onClick={handleSubmit}
              disabled={!smiles && droppedAtoms.length === 0}
              className="w-full px-4 py-2 bg-brand-orange text-white rounded-lg hover:bg-brand-orange-dark transition font-quando disabled:opacity-50 disabled:cursor-not-allowed text-sm"
            >
              Create Molecule →
            </button>
          </div>
        </div>
      </div>
    </div>
  );
}
