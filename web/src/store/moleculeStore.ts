import { create } from "zustand";
import { api } from "@/lib/api/client";
import type { Molecule } from "@/types/api";

interface MoleculeStore {
  molecules: Molecule[];
  currentMolecule: Molecule | null;
  isLoading: boolean;
  error: string | null;

  fetchMolecules: () => Promise<void>;
  createMolecule: (data: any) => Promise<Molecule>;
  setCurrentMolecule: (molecule: Molecule | null) => void;
  deleteMolecule: (id: string) => Promise<void>;
}

export const useMoleculeStore = create<MoleculeStore>((set, get) => ({
  molecules: [],
  currentMolecule: null,
  isLoading: false,
  error: null,

  fetchMolecules: async () => {
    set({ isLoading: true, error: null });
    try {
      const response = await api.molecules.list();
      set({ molecules: response.molecules, isLoading: false });
    } catch (error: any) {
      set({ error: error.message, isLoading: false });
    }
  },

  createMolecule: async (data: any) => {
    set({ isLoading: true, error: null });
    try {
      const molecule = await api.molecules.create(data);
      set((state) => ({
        molecules: [...state.molecules, molecule],
        currentMolecule: molecule,
        isLoading: false,
      }));
      return molecule;
    } catch (error: any) {
      set({ error: error.message, isLoading: false });
      throw error;
    }
  },

  setCurrentMolecule: (molecule: Molecule | null) => {
    set({ currentMolecule: molecule });
  },

  deleteMolecule: async (id: string) => {
    try {
      await api.molecules.delete(id);
      set((state) => ({
        molecules: state.molecules.filter((m) => m.molecule_id !== id),
        currentMolecule:
          state.currentMolecule?.molecule_id === id ? null : state.currentMolecule,
      }));
    } catch (error: any) {
      set({ error: error.message });
      throw error;
    }
  },
}));
