"use client";

import { useState, useRef } from "react";
import { ZoomIn, ZoomOut, RotateCcw, Download } from "lucide-react";

interface CircuitStats {
  n_qubits: number;
  depth: number;
  total_gates: number;
  gate_counts: Record<string, number>;
}

interface StructuredGate {
  type: string;
  qubits: number[];
  params?: any[];
}

interface QuantumCircuitViewerProps {
  circuitDiagram: string;
  stats: CircuitStats;
  gates?: StructuredGate[];  // Optional structured gate data from backend
  circuitImage?: string;  // Optional base64 PNG from Matplotlib
  className?: string;
}

interface Gate {
  qubit: number;
  column: number;
  type: string;
  label: string;
}

interface CNOTGate {
  control: number;
  target: number;
  column: number;
}

// Professional color scheme matching research papers
const GATE_COLORS: Record<string, { bg: string; text: string; border: string }> = {
  // Pauli Gates
  'H': { bg: '#FFD700', text: '#000000', border: '#000000' },      // Gold - Hadamard
  'X': { bg: '#4169E1', text: '#FFFFFF', border: '#000080' },      // Royal Blue - Pauli-X
  'Y': { bg: '#32CD32', text: '#FFFFFF', border: '#006400' },      // Lime Green - Pauli-Y
  'Z': { bg: '#FF6347', text: '#FFFFFF', border: '#8B0000' },      // Tomato - Pauli-Z
  'I': { bg: '#F5F5F5', text: '#000000', border: '#808080' },      // White Smoke - Identity

  // Rotation Gates
  'RX': { bg: '#9370DB', text: '#FFFFFF', border: '#4B0082' },     // Medium Purple - Rotation-X
  'RY': { bg: '#FF69B4', text: '#FFFFFF', border: '#C71585' },     // Hot Pink - Rotation-Y
  'RZ': { bg: '#20B2AA', text: '#FFFFFF', border: '#008B8B' },     // Light Sea Green - Rotation-Z
  'Rx': { bg: '#9370DB', text: '#FFFFFF', border: '#4B0082' },
  'Ry': { bg: '#FF69B4', text: '#FFFFFF', border: '#C71585' },
  'Rz': { bg: '#20B2AA', text: '#FFFFFF', border: '#008B8B' },

  // Phase Gates
  'S': { bg: '#87CEEB', text: '#000000', border: '#4682B4' },      // Sky Blue - S gate
  'T': { bg: '#DDA0DD', text: '#000000', border: '#8B008B' },      // Plum - T gate
  'S†': { bg: '#B0C4DE', text: '#000000', border: '#6495ED' },     // Light Steel Blue - S dagger
  'T†': { bg: '#BA55D3', text: '#FFFFFF', border: '#8B008B' },     // Medium Orchid - T dagger
  'SDAGGER': { bg: '#B0C4DE', text: '#000000', border: '#6495ED' },
  'TDAGGER': { bg: '#BA55D3', text: '#FFFFFF', border: '#8B008B' },

  // Multi-qubit Gates
  'CNOT': { bg: '#FF4500', text: '#FFFFFF', border: '#8B0000' },   // Orange Red - CNOT
  'CX': { bg: '#FF4500', text: '#FFFFFF', border: '#8B0000' },     // Same as CNOT
  'CZ': { bg: '#DC143C', text: '#FFFFFF', border: '#8B0000' },     // Crimson - CZ
  'CY': { bg: '#FF1493', text: '#FFFFFF', border: '#8B008B' },     // Deep Pink - CY
  'SWAP': { bg: '#FF8C00', text: '#FFFFFF', border: '#FF4500' },   // Dark Orange - SWAP
  'TOFFOLI': { bg: '#8B4513', text: '#FFFFFF', border: '#654321' }, // Saddle Brown - Toffoli
  'CCX': { bg: '#8B4513', text: '#FFFFFF', border: '#654321' },    // Same as Toffoli

  // Special Gates
  'U': { bg: '#9932CC', text: '#FFFFFF', border: '#4B0082' },      // Dark Orchid - Universal
  'U1': { bg: '#9932CC', text: '#FFFFFF', border: '#4B0082' },
  'U2': { bg: '#8A2BE2', text: '#FFFFFF', border: '#4B0082' },
  'U3': { bg: '#9400D3', text: '#FFFFFF', border: '#4B0082' },
  'P': { bg: '#6A5ACD', text: '#FFFFFF', border: '#483D8B' },      // Slate Blue - Phase
  'PHASE': { bg: '#6A5ACD', text: '#FFFFFF', border: '#483D8B' },
};

const DEFAULT_COLORS = { bg: '#D3D3D3', text: '#000000', border: '#696969' };

export default function QuantumCircuitViewer({
  circuitDiagram,
  stats,
  gates: structuredGates,
  circuitImage,
  className = "",
}: QuantumCircuitViewerProps) {
  const [scale, setScale] = useState(1.0);
  const [hoveredGate, setHoveredGate] = useState<number | null>(null);
  const svgRef = useRef<SVGSVGElement>(null);

  // Export circuit as PNG
  const exportAsPNG = () => {
    if (circuitImage) {
      // If we have a base64 image, download that
      const link = document.createElement('a');
      link.href = `data:image/png;base64,${circuitImage}`;
      link.download = 'quantum_circuit.png';
      link.click();
    } else if (svgRef.current) {
      // Convert SVG to PNG
      const svgElement = svgRef.current;
      const svgData = new XMLSerializer().serializeToString(svgElement);
      const canvas = document.createElement('canvas');
      const ctx = canvas.getContext('2d');
      const img = new Image();

      canvas.width = svgElement.viewBox.baseVal.width * 2; // 2x for better quality
      canvas.height = svgElement.viewBox.baseVal.height * 2;

      img.onload = () => {
        if (ctx) {
          ctx.scale(2, 2);
          ctx.drawImage(img, 0, 0);
          canvas.toBlob((blob) => {
            if (blob) {
              const url = URL.createObjectURL(blob);
              const link = document.createElement('a');
              link.href = url;
              link.download = 'quantum_circuit.png';
              link.click();
              URL.revokeObjectURL(url);
            }
          });
        }
      };

      img.src = 'data:image/svg+xml;base64,' + btoa(svgData);
    }
  };

  // Export circuit as SVG
  const exportAsSVG = () => {
    if (svgRef.current) {
      const svgData = new XMLSerializer().serializeToString(svgRef.current);
      const blob = new Blob([svgData], { type: 'image/svg+xml;charset=utf-8' });
      const url = URL.createObjectURL(blob);
      const link = document.createElement('a');
      link.href = url;
      link.download = 'quantum_circuit.svg';
      link.click();
      URL.revokeObjectURL(url);
    }
  };

  // If we have a Matplotlib-generated image, display that instead
  if (circuitImage) {
    return (
      <div className={`relative ${className}`}>
        {/* Compact Toolbar */}
        <div className="absolute top-1 right-1 z-10 flex gap-1 bg-white/95 backdrop-blur-sm rounded-md p-1 border border-gray-300 shadow-sm">
          <button
            onClick={() => setScale(prev => Math.min(2.5, prev + 0.2))}
            className="p-1.5 hover:bg-gray-100 rounded text-gray-700 transition-colors"
            title="Zoom In"
          >
            <ZoomIn className="w-3.5 h-3.5" />
          </button>
          <button
            onClick={() => setScale(prev => Math.max(0.6, prev - 0.2))}
            className="p-1.5 hover:bg-gray-100 rounded text-gray-700 transition-colors"
            title="Zoom Out"
          >
            <ZoomOut className="w-3.5 h-3.5" />
          </button>
          <button
            onClick={() => setScale(1.0)}
            className="p-1.5 hover:bg-gray-100 rounded text-gray-700 transition-colors"
            title="Reset Zoom"
          >
            <RotateCcw className="w-3.5 h-3.5" />
          </button>
          <div className="h-4 w-px bg-gray-300 my-1"></div>
          <button
            onClick={exportAsPNG}
            className="p-1.5 hover:bg-brand-orange/10 hover:text-brand-orange rounded text-gray-700 transition-colors"
            title="Download as PNG"
          >
            <Download className="w-3.5 h-3.5" />
          </button>
          <div className="px-1.5 py-1.5 text-[10px] font-mono text-gray-600 flex items-center">
            {Math.round(scale * 100)}%
          </div>
        </div>

        {/* Matplotlib Circuit Image */}
        <div className="overflow-auto border border-gray-300 rounded-md bg-white" style={{ maxHeight: '450px' }}>
          <img
            src={`data:image/png;base64,${circuitImage}`}
            alt="Quantum Circuit Diagram"
            style={{
              width: `${scale * 100}%`,
              height: 'auto',
              display: 'block'
            }}
          />
        </div>

        {/* Compact Legend */}
        <div className="mt-2 p-2 bg-gray-50 rounded-md border border-gray-200">
          <div className="flex flex-wrap gap-2 items-center">
            <span className="text-[10px] font-semibold text-gray-600 mr-1">Gates:</span>
            {Object.entries(stats.gate_counts || {}).map(([gate, count]) => {
              const colors = GATE_COLORS[gate] || DEFAULT_COLORS;
              return (
                <div key={gate} className="flex items-center gap-1">
                  <div
                    className="w-5 h-5 rounded flex items-center justify-center text-[9px] font-bold border"
                    style={{
                      backgroundColor: colors.bg,
                      color: colors.text,
                      borderColor: colors.border,
                      borderWidth: '1.5px'
                    }}
                  >
                    {gate}
                  </div>
                  <span className="text-[10px] text-gray-600">×{count}</span>
                </div>
              );
            })}
          </div>
        </div>

        {/* Compact Info Bar */}
        <div className="mt-1 flex gap-3 text-[10px] text-gray-600 font-mono">
          <span>{stats.n_qubits} qubits</span>
          <span>•</span>
          <span>depth {stats.depth}</span>
          <span>•</span>
          <span>{stats.total_gates} gates</span>
        </div>
      </div>
    );
  }

  // Parse structured gates (preferred method)
  const parseStructuredGates = () => {
    if (!structuredGates || structuredGates.length === 0) {
      return null;
    }

    const gates: Gate[] = [];
    const cnotGates: CNOTGate[] = [];
    let columnIdx = 0;

    // Process gates in order - each gate gets its own column
    structuredGates.forEach((gate) => {
      const gateType = gate.type.toUpperCase();

      if (gate.qubits.length === 1) {
        // Single-qubit gate
        const qubit = gate.qubits[0];
        let label = gateType;

        // Add parameter notation if present - display actual values like Qiskit
        if (gate.params && gate.params.length > 0) {
          const params = gate.params.map((p: number) => {
            // Format parameter value to 3 decimal places
            return typeof p === 'number' ? p.toFixed(3) : p;
          }).join(',');

          // For rotation gates, use Qiskit-style notation: Rz(θ)
          if (['RX', 'RY', 'RZ', 'Rx', 'Ry', 'Rz'].includes(gateType)) {
            label = `${gateType}(${params})`;
          } else if (['U1', 'U2', 'U3', 'P', 'PHASE'].includes(gateType)) {
            label = `${gateType}(${params})`;
          } else {
            label = `${gateType}(${params})`;
          }
        }

        gates.push({
          qubit,
          column: columnIdx,
          type: gateType,
          label,
        });
      } else if (gate.qubits.length === 2) {
        // Two-qubit gate (CNOT, CZ, etc.)
        const [control, target] = gate.qubits;
        cnotGates.push({
          control,
          target,
          column: columnIdx,
        });
      }

      columnIdx++;
    });

    return {
      gates,
      cnotGates,
      nQubits: stats.n_qubits,
      maxColumn: columnIdx - 1,
    };
  };

  // Parse circuit from ASCII (fallback method)
  const parseCircuit = () => {
    const lines = circuitDiagram.trim().split('\n');
    const gates: Gate[] = [];
    const cnotGates: CNOTGate[] = [];

    // First pass: determine maximum number of columns by finding all gate positions
    let maxPatterns = 0;
    const allPatterns: string[][] = [];

    lines.forEach((line, qubitIdx) => {
      const parts = line.split(':');
      if (parts.length < 2) {
        allPatterns.push([]);
        return;
      }

      const gatesSection = parts[1];
      const patterns = gatesSection.match(/\[(.*?)\]|──●──|──⊕──|──│──|─────/g) || [];
      allPatterns.push(patterns);
      maxPatterns = Math.max(maxPatterns, patterns.length);
    });

    // Second pass: identify columns with actual gates (not just spacers)
    const columnMap: number[] = []; // Maps pattern index to actual column number
    let currentCol = 0;

    for (let patIdx = 0; patIdx < maxPatterns; patIdx++) {
      let hasGate = false;
      for (let qubitIdx = 0; qubitIdx < allPatterns.length; qubitIdx++) {
        const pattern = allPatterns[qubitIdx][patIdx];
        if (pattern && pattern !== '─────' && pattern !== '──│──') {
          hasGate = true;
          break;
        }
      }
      if (hasGate) {
        columnMap[patIdx] = currentCol;
        currentCol++;
      } else {
        columnMap[patIdx] = -1; // Spacer, no column
      }
    }

    // Third pass: create gates with proper column indices
    allPatterns.forEach((patterns, qubitIdx) => {
      patterns.forEach((pattern, patIdx) => {
        const colIdx = columnMap[patIdx];
        if (colIdx === -1) return; // Skip spacers

        if (pattern.startsWith('[')) {
          // Single-qubit gate
          const content = pattern.slice(1, -1);
          const gateType = content.split('(')[0];
          gates.push({
            qubit: qubitIdx,
            column: colIdx,
            type: gateType,
            label: content,
          });
        } else if (pattern === '──●──') {
          // Control
          cnotGates.push({
            control: qubitIdx,
            target: -1,
            column: colIdx
          });
        } else if (pattern === '──⊕──') {
          // Target
          const cnot = cnotGates.find(g => g.column === colIdx && g.target === -1);
          if (cnot) cnot.target = qubitIdx;
        }
      });
    });

    return {
      gates,
      cnotGates: cnotGates.filter(g => g.target !== -1),
      nQubits: lines.length,
      maxColumn: Math.max(
        ...gates.map(g => g.column),
        ...cnotGates.map(g => g.column),
        0
      )
    };
  };

  // Use structured gates if available, otherwise parse ASCII
  const parsedData = parseStructuredGates() || parseCircuit();
  const { gates, cnotGates, nQubits, maxColumn } = parsedData;

  // Compact Qiskit-style layout - much tighter spacing
  const GATE_SIZE = 32;           // Reduced from 40
  const WIRE_SPACING = 40;        // Reduced from 50 - tighter vertical spacing
  const COLUMN_WIDTH = 50;        // Reduced from 70 - gates closer horizontally
  const MARGIN_LEFT = 40;         // Reduced from 50
  const MARGIN_TOP = 15;          // Reduced from 20
  const MARGIN_RIGHT = 15;        // Reduced from 20
  const MARGIN_BOTTOM = 15;       // Reduced from 20

  const svgWidth = MARGIN_LEFT + (maxColumn + 1) * COLUMN_WIDTH + MARGIN_RIGHT;
  const svgHeight = MARGIN_TOP + nQubits * WIRE_SPACING + MARGIN_BOTTOM;

  const handleZoomIn = () => setScale(prev => Math.min(2.5, prev + 0.2));
  const handleZoomOut = () => setScale(prev => Math.max(0.6, prev - 0.2));
  const handleReset = () => setScale(1.0);

  return (
    <div className={`relative ${className}`}>
      {/* Compact Toolbar */}
      <div className="absolute top-1 right-1 z-10 flex gap-1 bg-white/95 backdrop-blur-sm rounded-md p-1 border border-gray-300 shadow-sm">
        <button
          onClick={handleZoomIn}
          className="p-1.5 hover:bg-gray-100 rounded text-gray-700 transition-colors"
          title="Zoom In"
        >
          <ZoomIn className="w-3.5 h-3.5" />
        </button>
        <button
          onClick={handleZoomOut}
          className="p-1.5 hover:bg-gray-100 rounded text-gray-700 transition-colors"
          title="Zoom Out"
        >
          <ZoomOut className="w-3.5 h-3.5" />
        </button>
        <button
          onClick={handleReset}
          className="p-1.5 hover:bg-gray-100 rounded text-gray-700 transition-colors"
          title="Reset Zoom"
        >
          <RotateCcw className="w-3.5 h-3.5" />
        </button>
        <div className="h-4 w-px bg-gray-300 my-1"></div>
        <button
          onClick={exportAsPNG}
          className="p-1.5 hover:bg-brand-orange/10 hover:text-brand-orange rounded text-gray-700 transition-colors"
          title="Download as PNG"
        >
          <Download className="w-3.5 h-3.5" />
        </button>
        <div className="px-1.5 py-1.5 text-[10px] font-mono text-gray-600 flex items-center">
          {Math.round(scale * 100)}%
        </div>
      </div>

      {/* SVG Circuit - Grid Layout */}
      <div className="overflow-auto border border-gray-300 rounded-md bg-white" style={{ maxHeight: '450px' }}>
        <svg
          ref={svgRef}
          width={svgWidth * scale}
          height={svgHeight * scale}
          viewBox={`0 0 ${svgWidth} ${svgHeight}`}
        >
          {/* White background */}
          <rect width={svgWidth} height={svgHeight} fill="#FFFFFF" />

          {/* Qubit wires - horizontal grid lines */}
          {Array.from({ length: nQubits }).map((_, i) => {
            const y = MARGIN_TOP + i * WIRE_SPACING;
            return (
              <g key={`wire-${i}`}>
                {/* Wire line - thinner for compact look */}
                <line
                  x1={MARGIN_LEFT}
                  y1={y}
                  x2={svgWidth - MARGIN_RIGHT}
                  y2={y}
                  stroke="#000000"
                  strokeWidth="1"
                />
                {/* Qubit label - smaller font */}
                <text
                  x={MARGIN_LEFT - 10}
                  y={y}
                  textAnchor="end"
                  dominantBaseline="middle"
                  fill="#000000"
                  fontSize="10"
                  fontFamily="sans-serif"
                  fontWeight="500"
                >
                  q{i}
                </text>
              </g>
            );
          })}

          {/* CNOT gates - compact style */}
          {cnotGates.map((cnot, idx) => {
            const x = MARGIN_LEFT + cnot.column * COLUMN_WIDTH;
            const y1 = MARGIN_TOP + cnot.control * WIRE_SPACING;
            const y2 = MARGIN_TOP + cnot.target * WIRE_SPACING;

            return (
              <g key={`cnot-${idx}`}>
                {/* Vertical connection - thinner */}
                <line
                  x1={x}
                  y1={y1}
                  x2={x}
                  y2={y2}
                  stroke="#000000"
                  strokeWidth="1.5"
                />
                {/* Control dot - smaller */}
                <circle
                  cx={x}
                  cy={y1}
                  r="4"
                  fill="#000000"
                />
                {/* Target circle with cross - smaller */}
                <circle
                  cx={x}
                  cy={y2}
                  r="10"
                  fill="none"
                  stroke="#000000"
                  strokeWidth="1.5"
                />
                <line
                  x1={x}
                  y1={y2 - 10}
                  x2={x}
                  y2={y2 + 10}
                  stroke="#000000"
                  strokeWidth="1.5"
                />
                <line
                  x1={x - 10}
                  y1={y2}
                  x2={x + 10}
                  y2={y2}
                  stroke="#000000"
                  strokeWidth="1.5"
                />
              </g>
            );
          })}

          {/* Single-qubit gates */}
          {gates.map((gate, idx) => {
            const x = MARGIN_LEFT + gate.column * COLUMN_WIDTH;
            const y = MARGIN_TOP + gate.qubit * WIRE_SPACING;
            const colors = GATE_COLORS[gate.type] || DEFAULT_COLORS;
            const isHovered = hoveredGate === idx;

            return (
              <g
                key={`gate-${idx}`}
                onMouseEnter={() => setHoveredGate(idx)}
                onMouseLeave={() => setHoveredGate(null)}
                style={{ cursor: 'pointer' }}
              >
                {/* Gate box - centered on wire */}
                <rect
                  x={x - GATE_SIZE / 2}
                  y={y - GATE_SIZE / 2}
                  width={GATE_SIZE}
                  height={GATE_SIZE}
                  fill={colors.bg}
                  stroke={colors.border}
                  strokeWidth={isHovered ? "3" : "2"}
                  rx="2"
                  style={{
                    transition: 'all 0.15s ease',
                    filter: isHovered ? 'brightness(1.1) drop-shadow(0px 2px 4px rgba(0,0,0,0.2))' : 'none'
                  }}
                />
                {/* Gate label - centered, show full label with params if available */}
                <text
                  x={x}
                  y={y}
                  textAnchor="middle"
                  dominantBaseline="central"
                  fill={colors.text}
                  fontSize={gate.label && gate.label.length > 6 ? "8" : (gate.label && gate.label.length > 3 ? "9" : (isHovered ? "11" : "10"))}
                  fontFamily="sans-serif"
                  fontWeight="600"
                  style={{
                    transition: 'font-size 0.15s ease',
                    pointerEvents: 'none'
                  }}
                >
                  {gate.label || gate.type}
                </text>
                {/* Tooltip on hover */}
                {isHovered && (
                  <g>
                    <rect
                      x={x - 50}
                      y={y - GATE_SIZE / 2 - 30}
                      width="100"
                      height="22"
                      fill="#1F2937"
                      rx="4"
                      opacity="0.95"
                    />
                    <text
                      x={x}
                      y={y - GATE_SIZE / 2 - 19}
                      textAnchor="middle"
                      fill="#FFFFFF"
                      fontSize="11"
                      fontFamily="sans-serif"
                      fontWeight="500"
                      style={{ pointerEvents: 'none' }}
                    >
                      {gate.label || gate.type}
                    </text>
                  </g>
                )}
              </g>
            );
          })}
        </svg>
      </div>

      {/* Compact Legend */}
      <div className="mt-2 p-2 bg-gray-50 rounded-md border border-gray-200">
        <div className="flex flex-wrap gap-2 items-center">
          <span className="text-[10px] font-semibold text-gray-600 mr-1">Gates:</span>
          {Object.entries(stats.gate_counts || {}).map(([gate, count]) => {
            const colors = GATE_COLORS[gate] || DEFAULT_COLORS;
            return (
              <div key={gate} className="flex items-center gap-1">
                <div
                  className="w-5 h-5 rounded flex items-center justify-center text-[9px] font-bold border"
                  style={{
                    backgroundColor: colors.bg,
                    color: colors.text,
                    borderColor: colors.border,
                    borderWidth: '1.5px'
                  }}
                >
                  {gate}
                </div>
                <span className="text-[10px] text-gray-600">×{count}</span>
              </div>
            );
          })}
        </div>
      </div>

      {/* Compact Info Bar */}
      <div className="mt-1 flex gap-3 text-[10px] text-gray-600 font-mono">
        <span>{nQubits} qubits</span>
        <span>•</span>
        <span>depth {maxColumn + 1}</span>
        <span>•</span>
        <span>{stats.total_gates} gates</span>
      </div>
    </div>
  );
}
