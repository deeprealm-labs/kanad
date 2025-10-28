// Utility functions for molecule data processing

export interface AtomData {
  symbol: string;
  x: number;
  y: number;
  z: number;
  atomicNumber?: number;
}

export interface MoleculeData {
  smiles?: string;
  atoms?: AtomData[];
  xyzData?: string;
  name?: string;
  formula?: string;
}

/**
 * Parse XYZ file content and extract atom coordinates
 */
export function parseXYZFile(content: string): AtomData[] {
  const lines = content.trim().split('\n');
  
  if (lines.length < 3) {
    throw new Error('Invalid XYZ file format');
  }

  const numAtoms = parseInt(lines[0]);
  if (isNaN(numAtoms) || numAtoms <= 0) {
    throw new Error('Invalid number of atoms in XYZ file');
  }

  const atoms: AtomData[] = [];
  const atomLines = lines.slice(2, 2 + numAtoms);

  for (const line of atomLines) {
    const parts = line.trim().split(/\s+/);
    if (parts.length < 4) {
      throw new Error(`Invalid atom line: ${line}`);
    }

    const symbol = parts[0];
    const x = parseFloat(parts[1]);
    const y = parseFloat(parts[2]);
    const z = parseFloat(parts[3]);

    if (isNaN(x) || isNaN(y) || isNaN(z)) {
      throw new Error(`Invalid coordinates in line: ${line}`);
    }

    atoms.push({
      symbol,
      x,
      y,
      z,
      atomicNumber: getAtomicNumber(symbol),
    });
  }

  return atoms;
}

/**
 * Convert atoms array to XYZ format string
 */
export function atomsToXYZ(atoms: AtomData[], comment?: string): string {
  const header = `${atoms.length}`;
  const commentLine = comment || 'Generated molecule';
  const atomLines = atoms.map(atom => 
    `${atom.symbol} ${atom.x.toFixed(6)} ${atom.y.toFixed(6)} ${atom.z.toFixed(6)}`
  );

  return [header, commentLine, ...atomLines].join('\n');
}

/**
 * Get atomic number from element symbol
 */
export function getAtomicNumber(symbol: string): number {
  const atomicNumbers: Record<string, number> = {
    'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8,
    'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15,
    'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Ti': 22, 'Fe': 26,
    'Cu': 29, 'Zn': 30, 'Br': 35, 'I': 53,
  };
  
  return atomicNumbers[symbol] || 0;
}

/**
 * Get element color based on atomic number (CPK coloring)
 */
export function getElementColor(atomicNumber: number): string {
  const colors: Record<number, string> = {
    1: '#FFFFFF',   // H - White
    6: '#909090',   // C - Dark grey
    7: '#3050F8',   // N - Blue
    8: '#FF0D0D',   // O - Red
    9: '#90E050',   // F - Light green
    15: '#FF8000',  // P - Orange
    16: '#FFFF30',  // S - Yellow
    17: '#1FF01F',  // Cl - Green
    35: '#A62929',  // Br - Dark red
    53: '#940094',  // I - Dark violet
  };
  
  return colors[atomicNumber] || '#C0C0C0'; // Default silver
}

/**
 * Validate SMILES string format
 */
export function validateSMILES(smiles: string): boolean {
  if (!smiles || typeof smiles !== 'string') {
    return false;
  }

  // Basic SMILES validation - check for common patterns
  const validChars = /^[A-Za-z0-9@+\-\[\]()=#\\\/\s]*$/;
  if (!validChars.test(smiles)) {
    return false;
  }

  // Check for balanced brackets
  let bracketCount = 0;
  for (const char of smiles) {
    if (char === '[') bracketCount++;
    if (char === ']') bracketCount--;
    if (bracketCount < 0) return false;
  }

  return bracketCount === 0;
}

/**
 * Generate simple 3D coordinates for common molecules
 * This is a fallback when 3Dmol.js can't parse SMILES
 */
export function generateSimpleCoordinates(smiles: string): AtomData[] {
  const normalized = smiles.trim().toUpperCase();
  
  // Common molecules with known structures
  const molecules: Record<string, AtomData[]> = {
    'H2': [
      { symbol: 'H', x: 0, y: 0, z: 0, atomicNumber: 1 },
      { symbol: 'H', x: 0.74, y: 0, z: 0, atomicNumber: 1 },
    ],
    'H2O': [
      { symbol: 'O', x: 0, y: 0, z: 0, atomicNumber: 8 },
      { symbol: 'H', x: 0.96, y: 0.24, z: 0, atomicNumber: 1 },
      { symbol: 'H', x: -0.24, y: 0.96, z: 0, atomicNumber: 1 },
    ],
    'CH4': [
      { symbol: 'C', x: 0, y: 0, z: 0, atomicNumber: 6 },
      { symbol: 'H', x: 1.09, y: 1.09, z: 1.09, atomicNumber: 1 },
      { symbol: 'H', x: -1.09, y: -1.09, z: 1.09, atomicNumber: 1 },
      { symbol: 'H', x: -1.09, y: 1.09, z: -1.09, atomicNumber: 1 },
      { symbol: 'H', x: 1.09, y: -1.09, z: -1.09, atomicNumber: 1 },
    ],
    'NH3': [
      { symbol: 'N', x: 0, y: 0, z: 0, atomicNumber: 7 },
      { symbol: 'H', x: 1.01, y: 0, z: 0, atomicNumber: 1 },
      { symbol: 'H', x: -0.34, y: 0.94, z: 0, atomicNumber: 1 },
      { symbol: 'H', x: -0.34, y: -0.47, z: 0.82, atomicNumber: 1 },
    ],
    'CO2': [
      { symbol: 'C', x: 0, y: 0, z: 0, atomicNumber: 6 },
      { symbol: 'O', x: 1.16, y: 0, z: 0, atomicNumber: 8 },
      { symbol: 'O', x: -1.16, y: 0, z: 0, atomicNumber: 8 },
    ],
  };

  return molecules[normalized] || [];
}

/**
 * Process file upload and extract molecule data
 */
export async function processFileUpload(file: File): Promise<MoleculeData> {
  return new Promise((resolve, reject) => {
    const reader = new FileReader();
    
    reader.onload = (e) => {
      try {
        const content = e.target?.result as string;
        
        if (file.name.endsWith('.xyz')) {
          const atoms = parseXYZFile(content);
          resolve({
            atoms,
            xyzData: content,
            name: file.name.replace('.xyz', ''),
          });
        } else {
          reject(new Error('Unsupported file format. Please upload an XYZ file.'));
        }
      } catch (error) {
        reject(error);
      }
    };
    
    reader.onerror = () => {
      reject(new Error('Failed to read file'));
    };
    
    reader.readAsText(file);
  });
}

/**
 * Create molecule data from various inputs
 */
export function createMoleculeData(
  smiles?: string,
  atoms?: AtomData[],
  xyzData?: string,
  name?: string
): MoleculeData {
  const molecule: MoleculeData = {};

  if (smiles && validateSMILES(smiles)) {
    molecule.smiles = smiles;
  }

  if (atoms && atoms.length > 0) {
    molecule.atoms = atoms;
    if (!molecule.xyzData) {
      molecule.xyzData = atomsToXYZ(atoms, name);
    }
  }

  if (xyzData) {
    molecule.xyzData = xyzData;
    try {
      molecule.atoms = parseXYZFile(xyzData);
    } catch (error) {
      console.warn('Failed to parse XYZ data:', error);
    }
  }

  if (name) {
    molecule.name = name;
  }

  return molecule;
}
