import { create } from "zustand";
import { api } from "@/lib/api/client";
import type { User } from "@/types/api";

interface AuthStore {
  user: User | null;
  token: string | null;
  isAuthenticated: boolean;
  isLoading: boolean;
  error: string | null;

  login: (email: string, password: string) => Promise<void>;
  register: (email: string, password: string, name?: string) => Promise<void>;
  logout: () => void;
  setUser: (user: User | null) => void;
}

export const useAuthStore = create<AuthStore>((set) => ({
  user: null,
  token: typeof window !== "undefined" ? localStorage.getItem("auth_token") : null,
  isAuthenticated: false,
  isLoading: false,
  error: null,

  login: async (email: string, password: string) => {
    set({ isLoading: true, error: null });
    try {
      const response = await api.auth.login(email, password);
      localStorage.setItem("auth_token", response.access_token);
      set({
        token: response.access_token,
        user: response.user,
        isAuthenticated: true,
        isLoading: false,
      });
    } catch (error: any) {
      set({ error: error.message, isLoading: false });
      throw error;
    }
  },

  register: async (email: string, password: string, name?: string) => {
    set({ isLoading: true, error: null });
    try {
      const response = await api.auth.register({ email, password, name });
      localStorage.setItem("auth_token", response.access_token);
      set({
        token: response.access_token,
        user: response.user,
        isAuthenticated: true,
        isLoading: false,
      });
    } catch (error: any) {
      set({ error: error.message, isLoading: false });
      throw error;
    }
  },

  logout: () => {
    localStorage.removeItem("auth_token");
    set({
      user: null,
      token: null,
      isAuthenticated: false,
    });
  },

  setUser: (user: User | null) => {
    set({ user, isAuthenticated: !!user });
  },
}));
