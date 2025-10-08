export interface LibraryMolecule {
  id: string;
  name: string;
  formula: string;
  smiles: string;
  description: string;
  basis: string;
  charge: number;
  multiplicity: number;
}

export interface MoleculeCategory {
  name: string;
  description: string;
  molecules: LibraryMolecule[];
}

export const moleculeLibrary: MoleculeCategory[] = [
  {
    name: "Small Molecules",
    description: "Common small molecules for testing and benchmarking",
    molecules: [
      {
        id: "h2",
        name: "Hydrogen",
        formula: "H₂",
        smiles: "[H][H]",
        description: "Simplest diatomic molecule",
        basis: "sto-3g",
        charge: 0,
        multiplicity: 1,
      },
      {
        id: "h2o",
        name: "Water",
        formula: "H₂O",
        smiles: "O",
        description: "Essential solvent molecule",
        basis: "sto-3g",
        charge: 0,
        multiplicity: 1,
      },
      {
        id: "co2",
        name: "Carbon Dioxide",
        formula: "CO₂",
        smiles: "O=C=O",
        description: "Linear triatomic molecule",
        basis: "sto-3g",
        charge: 0,
        multiplicity: 1,
      },
      {
        id: "nh3",
        name: "Ammonia",
        formula: "NH₃",
        smiles: "N",
        description: "Pyramidal molecule with lone pair",
        basis: "sto-3g",
        charge: 0,
        multiplicity: 1,
      },
      {
        id: "ch4",
        name: "Methane",
        formula: "CH₄",
        smiles: "C",
        description: "Tetrahedral hydrocarbon",
        basis: "sto-3g",
        charge: 0,
        multiplicity: 1,
      },
      {
        id: "lih",
        name: "Lithium Hydride",
        formula: "LiH",
        smiles: "[LiH]",
        description: "Simple ionic molecule",
        basis: "sto-3g",
        charge: 0,
        multiplicity: 1,
      },
    ],
  },
  {
    name: "Metallurgy",
    description: "Metallic bonds and compounds for materials science",
    molecules: [
      {
        id: "fe_co",
        name: "Iron-Carbon Bond",
        formula: "FeC",
        smiles: "[Fe]=C",
        description: "Steel composition modeling",
        basis: "sto-3g",
        charge: 0,
        multiplicity: 1,
      },
      {
        id: "tio2",
        name: "Titanium Dioxide",
        formula: "TiO₂",
        smiles: "O=[Ti]=O",
        description: "Photocatalyst and pigment",
        basis: "sto-3g",
        charge: 0,
        multiplicity: 1,
      },
      {
        id: "feo",
        name: "Iron Oxide",
        formula: "FeO",
        smiles: "[Fe]=O",
        description: "Common corrosion product",
        basis: "sto-3g",
        charge: 0,
        multiplicity: 1,
      },
    ],
  },
  {
    name: "Bioscience",
    description: "Biological molecules and amino acids",
    molecules: [
      {
        id: "caffeine",
        name: "Caffeine",
        formula: "C₈H₁₀N₄O₂",
        smiles: "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        description: "Common stimulant compound",
        basis: "sto-3g",
        charge: 0,
        multiplicity: 1,
      },
      {
        id: "alanine",
        name: "Alanine",
        formula: "C₃H₇NO₂",
        smiles: "CC(N)C(=O)O",
        description: "Simple amino acid",
        basis: "sto-3g",
        charge: 0,
        multiplicity: 1,
      },
      {
        id: "glycine",
        name: "Glycine",
        formula: "C₂H₅NO₂",
        smiles: "NCC(=O)O",
        description: "Simplest amino acid",
        basis: "sto-3g",
        charge: 0,
        multiplicity: 1,
      },
    ],
  },
  {
    name: "Pharmaceutical",
    description: "Drug molecules and medicinal compounds",
    molecules: [
      {
        id: "aspirin",
        name: "Aspirin",
        formula: "C₉H₈O₄",
        smiles: "CC(=O)OC1=CC=CC=C1C(=O)O",
        description: "Common pain reliever",
        basis: "sto-3g",
        charge: 0,
        multiplicity: 1,
      },
      {
        id: "ethanol",
        name: "Ethanol",
        formula: "C₂H₆O",
        smiles: "CCO",
        description: "Common alcohol",
        basis: "sto-3g",
        charge: 0,
        multiplicity: 1,
      },
      {
        id: "acetic_acid",
        name: "Acetic Acid",
        formula: "C₂H₄O₂",
        smiles: "CC(=O)O",
        description: "Weak organic acid",
        basis: "sto-3g",
        charge: 0,
        multiplicity: 1,
      },
    ],
  },
  {
    name: "Organic Chemistry",
    description: "Common organic molecules and hydrocarbons",
    molecules: [
      {
        id: "benzene",
        name: "Benzene",
        formula: "C₆H₆",
        smiles: "c1ccccc1",
        description: "Aromatic hydrocarbon",
        basis: "sto-3g",
        charge: 0,
        multiplicity: 1,
      },
      {
        id: "formaldehyde",
        name: "Formaldehyde",
        formula: "CH₂O",
        smiles: "C=O",
        description: "Simplest aldehyde",
        basis: "sto-3g",
        charge: 0,
        multiplicity: 1,
      },
      {
        id: "acetone",
        name: "Acetone",
        formula: "C₃H₆O",
        smiles: "CC(=O)C",
        description: "Common ketone solvent",
        basis: "sto-3g",
        charge: 0,
        multiplicity: 1,
      },
    ],
  },
];
