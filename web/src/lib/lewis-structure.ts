// Lewis Structure Generator
// Follows VSEPR theory and chemistry laws

interface Atom {
  element: string;
  x: number;
  y: number;
  valenceElectrons: number;
  bonds: Bond[];
  lonePairs: number;
  formalCharge: number;
  partialCharge?: string; // δ+, δ-, etc.
}

interface Bond {
  type: "single" | "double" | "triple";
  to: number; // index of connected atom
  electrons: number;
}

interface LewisStructure {
  atoms: Atom[];
  totalElectrons: number;
  geometry: string;
  bondAngles: number[];
}

// Valence electrons for each element
const VALENCE_ELECTRONS: Record<string, number> = {
  H: 1,
  C: 4,
  N: 5,
  O: 6,
  F: 7,
  Cl: 7,
  Br: 7,
  I: 7,
  S: 6,
  P: 5,
  B: 3,
  Si: 4,
  Li: 1,
  Na: 1,
  K: 1,
  Mg: 2,
  Ca: 2,
  Fe: 2, // simplified
  Cu: 1, // simplified
  Zn: 2,
  Ti: 4, // simplified
  Al: 3,
  He: 2,
  Ne: 8,
  Ar: 8,
};

// Electronegativity values (Pauling scale)
const ELECTRONEGATIVITY: Record<string, number> = {
  H: 2.20,
  C: 2.55,
  N: 3.04,
  O: 3.44,
  F: 3.98,
  Cl: 3.16,
  Br: 2.96,
  I: 2.66,
  S: 2.58,
  P: 2.19,
  B: 2.04,
  Si: 1.90,
  Li: 0.98,
  Na: 0.93,
  K: 0.82,
  Mg: 1.31,
  Ca: 1.00,
  Fe: 1.83,
  Cu: 1.90,
  Zn: 1.65,
  Ti: 1.54,
  Al: 1.61,
};

export class LewisStructureGenerator {
  /**
   * Generate Lewis structure from SMILES
   */
  static fromSMILES(smiles: string): LewisStructure | null {
    // Parse common SMILES patterns
    const normalized = smiles.trim().toUpperCase();

    // Diatomic molecules
    if (normalized === "O=C=O" || normalized === "CO2") {
      return this.generateCO2();
    }
    if (normalized === "O" || normalized === "H2O") {
      return this.generateH2O();
    }
    if (normalized === "[H][H]" || normalized === "H2") {
      return this.generateH2();
    }
    if (normalized === "N#N" || normalized === "N2") {
      return this.generateN2();
    }
    if (normalized === "CCO" || normalized === "C2H5OH") {
      return this.generateEthanol();
    }
    if (normalized === "C1CCCCC1" || normalized === "C6H6") {
      return this.generateBenzene();
    }
    if (normalized === "[H]CL" || normalized === "HCL" || normalized === "CL") {
      return this.generateHCl();
    }
    if (normalized === "N") {
      return this.generateNH3();
    }
    if (normalized === "C") {
      return this.generateCH4();
    }
    if (normalized === "CC(=O)O") {
      return this.generateAceticAcid();
    }

    // Default: try to parse simple molecule
    return this.generateFromFormula(smiles);
  }

  /**
   * Generate H-Cl molecule (from image)
   */
  static generateHCl(): LewisStructure {
    const atoms: Atom[] = [
      {
        element: "H",
        x: 200,
        y: 200,
        valenceElectrons: 1,
        bonds: [{ type: "single", to: 1, electrons: 2 }],
        lonePairs: 0,
        formalCharge: 0,
        partialCharge: "δ+",
      },
      {
        element: "Cl",
        x: 400,
        y: 200,
        valenceElectrons: 7,
        bonds: [{ type: "single", to: 0, electrons: 2 }],
        lonePairs: 3, // 3 lone pairs (6 electrons)
        formalCharge: 0,
        partialCharge: "δ−",
      },
    ];

    return {
      atoms,
      totalElectrons: 8,
      geometry: "linear",
      bondAngles: [180],
    };
  }

  /**
   * Generate H2O molecule
   */
  static generateH2O(): LewisStructure {
    const atoms: Atom[] = [
      {
        element: "H",
        x: 150,
        y: 250,
        valenceElectrons: 1,
        bonds: [{ type: "single", to: 1, electrons: 2 }],
        lonePairs: 0,
        formalCharge: 0,
      },
      {
        element: "O",
        x: 250,
        y: 200,
        valenceElectrons: 6,
        bonds: [
          { type: "single", to: 0, electrons: 2 },
          { type: "single", to: 2, electrons: 2 },
        ],
        lonePairs: 2, // 2 lone pairs (4 electrons)
        formalCharge: 0,
      },
      {
        element: "H",
        x: 350,
        y: 250,
        valenceElectrons: 1,
        bonds: [{ type: "single", to: 1, electrons: 2 }],
        lonePairs: 0,
        formalCharge: 0,
      },
    ];

    return {
      atoms,
      totalElectrons: 8,
      geometry: "bent",
      bondAngles: [104.5],
    };
  }

  /**
   * Generate H2 molecule
   */
  static generateH2(): LewisStructure {
    const atoms: Atom[] = [
      {
        element: "H",
        x: 200,
        y: 200,
        valenceElectrons: 1,
        bonds: [{ type: "single", to: 1, electrons: 2 }],
        lonePairs: 0,
        formalCharge: 0,
      },
      {
        element: "H",
        x: 300,
        y: 200,
        valenceElectrons: 1,
        bonds: [{ type: "single", to: 0, electrons: 2 }],
        lonePairs: 0,
        formalCharge: 0,
      },
    ];

    return {
      atoms,
      totalElectrons: 2,
      geometry: "linear",
      bondAngles: [180],
    };
  }

  /**
   * Generate NH3 molecule
   */
  static generateNH3(): LewisStructure {
    const atoms: Atom[] = [
      {
        element: "N",
        x: 250,
        y: 200,
        valenceElectrons: 5,
        bonds: [
          { type: "single", to: 1, electrons: 2 },
          { type: "single", to: 2, electrons: 2 },
          { type: "single", to: 3, electrons: 2 },
        ],
        lonePairs: 1, // 1 lone pair
        formalCharge: 0,
      },
      {
        element: "H",
        x: 150,
        y: 250,
        valenceElectrons: 1,
        bonds: [{ type: "single", to: 0, electrons: 2 }],
        lonePairs: 0,
        formalCharge: 0,
      },
      {
        element: "H",
        x: 250,
        y: 300,
        valenceElectrons: 1,
        bonds: [{ type: "single", to: 0, electrons: 2 }],
        lonePairs: 0,
        formalCharge: 0,
      },
      {
        element: "H",
        x: 350,
        y: 250,
        valenceElectrons: 1,
        bonds: [{ type: "single", to: 0, electrons: 2 }],
        lonePairs: 0,
        formalCharge: 0,
      },
    ];

    return {
      atoms,
      totalElectrons: 8,
      geometry: "trigonal pyramidal",
      bondAngles: [107],
    };
  }

  /**
   * Generate CH4 molecule
   */
  static generateCH4(): LewisStructure {
    const atoms: Atom[] = [
      {
        element: "C",
        x: 250,
        y: 200,
        valenceElectrons: 4,
        bonds: [
          { type: "single", to: 1, electrons: 2 },
          { type: "single", to: 2, electrons: 2 },
          { type: "single", to: 3, electrons: 2 },
          { type: "single", to: 4, electrons: 2 },
        ],
        lonePairs: 0,
        formalCharge: 0,
      },
      {
        element: "H",
        x: 150,
        y: 150,
        valenceElectrons: 1,
        bonds: [{ type: "single", to: 0, electrons: 2 }],
        lonePairs: 0,
        formalCharge: 0,
      },
      {
        element: "H",
        x: 350,
        y: 150,
        valenceElectrons: 1,
        bonds: [{ type: "single", to: 0, electrons: 2 }],
        lonePairs: 0,
        formalCharge: 0,
      },
      {
        element: "H",
        x: 150,
        y: 250,
        valenceElectrons: 1,
        bonds: [{ type: "single", to: 0, electrons: 2 }],
        lonePairs: 0,
        formalCharge: 0,
      },
      {
        element: "H",
        x: 350,
        y: 250,
        valenceElectrons: 1,
        bonds: [{ type: "single", to: 0, electrons: 2 }],
        lonePairs: 0,
        formalCharge: 0,
      },
    ];

    return {
      atoms,
      totalElectrons: 8,
      geometry: "tetrahedral",
      bondAngles: [109.5],
    };
  }

  /**
   * Generate CO2 molecule
   */
  static generateCO2(): LewisStructure {
    const atoms: Atom[] = [
      {
        element: "O",
        x: 150,
        y: 200,
        valenceElectrons: 6,
        bonds: [{ type: "double", to: 1, electrons: 4 }],
        lonePairs: 2,
        formalCharge: 0,
      },
      {
        element: "C",
        x: 250,
        y: 200,
        valenceElectrons: 4,
        bonds: [
          { type: "double", to: 0, electrons: 4 },
          { type: "double", to: 2, electrons: 4 },
        ],
        lonePairs: 0,
        formalCharge: 0,
      },
      {
        element: "O",
        x: 350,
        y: 200,
        valenceElectrons: 6,
        bonds: [{ type: "double", to: 1, electrons: 4 }],
        lonePairs: 2,
        formalCharge: 0,
      },
    ];

    return {
      atoms,
      totalElectrons: 16,
      geometry: "linear",
      bondAngles: [180],
    };
  }

  /**
   * Generate N2 molecule
   */
  static generateN2(): LewisStructure {
    const atoms: Atom[] = [
      {
        element: "N",
        x: 200,
        y: 200,
        valenceElectrons: 5,
        bonds: [{ type: "triple", to: 1, electrons: 6 }],
        lonePairs: 1,
        formalCharge: 0,
      },
      {
        element: "N",
        x: 300,
        y: 200,
        valenceElectrons: 5,
        bonds: [{ type: "triple", to: 0, electrons: 6 }],
        lonePairs: 1,
        formalCharge: 0,
      },
    ];

    return {
      atoms,
      totalElectrons: 10,
      geometry: "linear",
      bondAngles: [180],
    };
  }

  /**
   * Generate Ethanol (C2H5OH)
   */
  static generateEthanol(): LewisStructure {
    // Simplified for now
    const atoms: Atom[] = [
      {
        element: "C",
        x: 150,
        y: 200,
        valenceElectrons: 4,
        bonds: [
          { type: "single", to: 1, electrons: 2 },
          { type: "single", to: 2, electrons: 2 },
        ],
        lonePairs: 0,
        formalCharge: 0,
      },
      {
        element: "C",
        x: 250,
        y: 200,
        valenceElectrons: 4,
        bonds: [
          { type: "single", to: 0, electrons: 2 },
          { type: "single", to: 3, electrons: 2 },
        ],
        lonePairs: 0,
        formalCharge: 0,
      },
      {
        element: "H",
        x: 100,
        y: 150,
        valenceElectrons: 1,
        bonds: [{ type: "single", to: 0, electrons: 2 }],
        lonePairs: 0,
        formalCharge: 0,
      },
      {
        element: "O",
        x: 350,
        y: 200,
        valenceElectrons: 6,
        bonds: [
          { type: "single", to: 1, electrons: 2 },
          { type: "single", to: 4, electrons: 2 },
        ],
        lonePairs: 2,
        formalCharge: 0,
      },
      {
        element: "H",
        x: 400,
        y: 180,
        valenceElectrons: 1,
        bonds: [{ type: "single", to: 3, electrons: 2 }],
        lonePairs: 0,
        formalCharge: 0,
      },
    ];

    return {
      atoms,
      totalElectrons: 20,
      geometry: "various",
      bondAngles: [],
    };
  }

  /**
   * Generate Benzene (simplified)
   */
  static generateBenzene(): LewisStructure {
    const centerX = 250;
    const centerY = 200;
    const radius = 80;
    const atoms: Atom[] = [];

    // Create hexagon of carbons
    for (let i = 0; i < 6; i++) {
      const angle = (i * 60 - 90) * (Math.PI / 180);
      atoms.push({
        element: "C",
        x: centerX + radius * Math.cos(angle),
        y: centerY + radius * Math.sin(angle),
        valenceElectrons: 4,
        bonds: [],
        lonePairs: 0,
        formalCharge: 0,
      });
    }

    return {
      atoms,
      totalElectrons: 30,
      geometry: "planar hexagonal",
      bondAngles: [120],
    };
  }

  /**
   * Generate Acetic Acid
   */
  static generateAceticAcid(): LewisStructure {
    const atoms: Atom[] = [
      {
        element: "C",
        x: 150,
        y: 200,
        valenceElectrons: 4,
        bonds: [{ type: "single", to: 1, electrons: 2 }],
        lonePairs: 0,
        formalCharge: 0,
      },
      {
        element: "C",
        x: 250,
        y: 200,
        valenceElectrons: 4,
        bonds: [
          { type: "single", to: 0, electrons: 2 },
          { type: "double", to: 2, electrons: 4 },
          { type: "single", to: 3, electrons: 2 },
        ],
        lonePairs: 0,
        formalCharge: 0,
      },
      {
        element: "O",
        x: 250,
        y: 100,
        valenceElectrons: 6,
        bonds: [{ type: "double", to: 1, electrons: 4 }],
        lonePairs: 2,
        formalCharge: 0,
      },
      {
        element: "O",
        x: 350,
        y: 200,
        valenceElectrons: 6,
        bonds: [
          { type: "single", to: 1, electrons: 2 },
          { type: "single", to: 4, electrons: 2 },
        ],
        lonePairs: 2,
        formalCharge: 0,
      },
      {
        element: "H",
        x: 400,
        y: 180,
        valenceElectrons: 1,
        bonds: [{ type: "single", to: 3, electrons: 2 }],
        lonePairs: 0,
        formalCharge: 0,
      },
    ];

    return {
      atoms,
      totalElectrons: 24,
      geometry: "various",
      bondAngles: [],
    };
  }

  /**
   * Fallback: Generate from simple formula
   */
  static generateFromFormula(formula: string): LewisStructure | null {
    // Very basic parser for simple molecules
    // Returns null if can't parse
    return null;
  }

  /**
   * Calculate partial charges based on electronegativity difference
   */
  static calculatePartialCharge(
    element1: string,
    element2: string
  ): [string, string] {
    const en1 = ELECTRONEGATIVITY[element1] || 2.5;
    const en2 = ELECTRONEGATIVITY[element2] || 2.5;
    const diff = Math.abs(en1 - en2);

    if (diff < 0.4) {
      return ["", ""]; // Non-polar
    }

    if (en1 > en2) {
      return ["δ−", "δ+"];
    } else {
      return ["δ+", "δ−"];
    }
  }
}
