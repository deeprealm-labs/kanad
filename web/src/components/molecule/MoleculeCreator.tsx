"use client";

import { useState, useRef, useEffect } from "react";
import { Search, X, Plus } from "lucide-react";
import { moleculeLibrary } from "@/data/molecule-library";
import Molecule3DViewer from "./Molecule3DViewer";
import { createMoleculeData, processFileUpload, AtomData as UtilAtomData } from "@/utils/moleculeUtils";
import * as api from "@/lib/api";

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
  const [xyzFile, setXyzFile] = useState<File | null>(null);
  const [xyzData, setXyzData] = useState<string>("");
  const [isProcessingFile, setIsProcessingFile] = useState(false);
  const [basisSets, setBasisSets] = useState<any[]>([]);
  const fileInputRef = useRef<HTMLInputElement>(null);

  // Load basis sets from API
  useEffect(() => {
    const loadBasisSets = async () => {
      try {
        const config = await api.getConfigurationOptions();
        setBasisSets(config.basis_sets || []);
      } catch (error) {
        console.error("Failed to load basis sets:", error);
        // Fallback to minimal basis sets
        setBasisSets([
          { value: "sto-3g", label: "STO-3G", category: "minimal" },
          { value: "6-31g", label: "6-31G", category: "split_valence" },
        ]);
      }
    };
    loadBasisSets();
  }, []);

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
    setSmiles(""); // Clear SMILES when atoms are added
    setSmilesInput(""); // Clear SMILES input
    setXyzFile(null); // Clear XYZ file
    setXyzData(""); // Clear XYZ data
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
      setXyzFile(null); // Clear XYZ file when SMILES is added
      setXyzData(""); // Clear XYZ data when SMILES is added
    }
  };

  const handleLoadMolecule = (molecule: any) => {
    setSmilesInput(molecule.smiles);
    setSmiles(molecule.smiles);
    setBasis(molecule.basis);
    setCharge(molecule.charge);
    setMultiplicity(molecule.multiplicity);
    setDroppedAtoms([]); // Clear atoms when molecule is loaded
    setXyzFile(null); // Clear XYZ file
    setXyzData(""); // Clear XYZ data
  };

  const handleFileUpload = async (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (!file) return;

    if (!file.name.endsWith('.xyz')) {
      alert('Please upload an XYZ file');
      return;
    }

    setIsProcessingFile(true);
    try {
      const content = await file.text();
      console.log("XYZ file content:", content);
      
      // Parse the XYZ file to get atoms
      const lines = content.trim().split('\n');
      const numAtoms = parseInt(lines[0]);
      const atomLines = lines.slice(2, 2 + numAtoms);
      
      const parsedAtoms = atomLines.map((line) => {
        const parts = line.trim().split(/\s+/);
        return {
          symbol: parts[0],
          x: parseFloat(parts[1]),
          y: parseFloat(parts[2]),
          z: parseFloat(parts[3]),
          atomicNumber: getAtomicNumber(parts[0]),
        };
      });

      console.log("Parsed atoms:", parsedAtoms);
      
      setXyzFile(file);
      setSmiles(''); // Clear SMILES when XYZ is loaded
      setDroppedAtoms([]); // Clear atoms when XYZ is loaded
      setSmilesInput(''); // Clear SMILES input
      
      // Store XYZ data for the 3D viewer
      setXyzData(content);
    } catch (error) {
      console.error('Error processing XYZ file:', error);
      alert('Error processing XYZ file: ' + (error as Error).message);
    } finally {
      setIsProcessingFile(false);
    }
  };

  const getAtomicNumber = (symbol: string): number => {
    const atomicNumbers: Record<string, number> = {
      'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8,
      'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15,
      'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Ti': 22, 'Fe': 26,
      'Cu': 29, 'Zn': 30, 'Br': 35, 'I': 53,
    };
    return atomicNumbers[symbol] || 0;
  };

  const handleSubmit = () => {
    if (smiles || droppedAtoms.length > 0 || xyzFile) {
      const moleculeData = createMoleculeData(
        smiles,
        droppedAtoms.map(atom => ({
          symbol: atom.symbol,
          x: atom.x,
          y: atom.y,
          z: atom.z,
          atomicNumber: atom.atomicNumber,
        })),
        xyzData, // Pass XYZ data to viewer
        xyzFile?.name.replace('.xyz', '')
      );

      onComplete({
        ...moleculeData,
        xyzFile,
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
                      <div className="text-[9px] text-black leading-none">
                        {atom.atomicNumber}
                      </div>
                      <div className="text-lg font-bold font-quando leading-none my-0.5 text-black">
                        {atom.symbol}
                      </div>
                      <div className="text-[9px] text-black leading-none">
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
                      <div className="text-[8px] text-black leading-none">
                        {atom.atomicNumber}
                      </div>
                      <div className="text-sm font-bold leading-none my-0.5 text-black">
                        {atom.symbol}
                      </div>
                      <div className="text-[8px] text-black leading-none">
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

        {/* Right: Configuration & 3D Preview */}
        <div className="flex flex-col h-full min-h-0">
          <h2 className="text-base font-quando font-semibold mb-2">
            Configuration & Preview
          </h2>

          {/* 3D Visualization - Now at the top and larger */}
          <div className="mb-3">
            {(() => {
              const moleculeData = {
                atoms: droppedAtoms.map(atom => ({
                  symbol: atom.symbol,
                  x: atom.x,
                  y: atom.y,
                  z: atom.z,
                  atomicNumber: atom.atomicNumber,
                })),
                smiles: smiles,
                xyzData: xyzData,
              };

              console.log("=== MOLECULE CREATOR DEBUG ===");
              console.log("Dropped atoms:", droppedAtoms);
              console.log("SMILES:", smiles);
              console.log("XYZ data:", xyzData);
              console.log("Molecule data being passed to 3D viewer:", moleculeData);

              return (
                <Molecule3DViewer
                  molecule={moleculeData}
                  height="300px"
                  showControls={true}
                />
              );
            })()}
          </div>

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
                {basisSets.length > 0 ? (
                  <>
                    {/* Group basis sets by category */}
                    {['minimal', 'split_valence', 'polarization', 'diffuse', 'correlation_consistent', 'def2', 'periodic'].map((category) => {
                      const categoryBases = basisSets.filter((b) => b.category === category);
                      if (categoryBases.length === 0) return null;

                      const categoryLabels: Record<string, string> = {
                        minimal: 'Minimal Basis',
                        split_valence: 'Split-Valence',
                        polarization: 'Polarization',
                        diffuse: 'Diffuse Functions',
                        correlation_consistent: 'Correlation-Consistent',
                        def2: 'Def2 Series',
                        periodic: 'Periodic Systems'
                      };

                      return (
                        <optgroup key={category} label={categoryLabels[category] || category}>
                          {categoryBases.map((b: any) => (
                            <option key={b.value} value={b.value}>
                              {b.label}
                            </option>
                          ))}
                        </optgroup>
                      );
                    })}
                  </>
                ) : (
                  <>
                    <option value="sto-3g">STO-3G</option>
                    <option value="6-31g">6-31G</option>
                  </>
                )}
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
                ref={fileInputRef}
                type="file"
                accept=".xyz"
                onChange={handleFileUpload}
                disabled={isProcessingFile}
                className="w-full px-3 py-1.5 border border-input bg-background rounded-md font-quando text-xs file:mr-2 file:px-3 file:py-1 file:rounded file:border-0 file:bg-brand-orange file:text-white file:text-xs disabled:opacity-50"
              />
              {xyzFile && (
                <div className="mt-1 text-xs text-muted-foreground font-quando">
                  ✓ {xyzFile.name}
                </div>
              )}
              {isProcessingFile && (
                <div className="mt-1 text-xs text-brand-orange font-quando">
                  Processing file...
                </div>
              )}
            </div>
          </div>

          <div className="mt-3 pt-3 border-t border-border">
            <button
              onClick={handleSubmit}
              disabled={!smiles && droppedAtoms.length === 0 && !xyzFile}
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
