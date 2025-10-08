import { create } from "zustand";
import { api } from "@/lib/api/client";
import type { Job, JobResults } from "@/types/api";

interface JobStore {
  jobs: Job[];
  currentJob: Job | null;
  currentResults: JobResults | null;
  isLoading: boolean;
  error: string | null;

  fetchJobs: () => Promise<void>;
  getJobStatus: (id: string) => Promise<Job>;
  getJobResults: (id: string) => Promise<JobResults>;
  cancelJob: (id: string) => Promise<void>;
  setCurrentJob: (job: Job | null) => void;
  updateJobProgress: (jobId: string, progress: number) => void;
}

export const useJobStore = create<JobStore>((set, get) => ({
  jobs: [],
  currentJob: null,
  currentResults: null,
  isLoading: false,
  error: null,

  fetchJobs: async () => {
    set({ isLoading: true, error: null });
    try {
      const response = await api.jobs.list();
      set({ jobs: response.jobs, isLoading: false });
    } catch (error: any) {
      set({ error: error.message, isLoading: false });
    }
  },

  getJobStatus: async (id: string) => {
    try {
      const job = await api.jobs.getStatus(id);
      set((state) => ({
        jobs: state.jobs.map((j) => (j.job_id === id ? job : j)),
        currentJob: state.currentJob?.job_id === id ? job : state.currentJob,
      }));
      return job;
    } catch (error: any) {
      set({ error: error.message });
      throw error;
    }
  },

  getJobResults: async (id: string) => {
    set({ isLoading: true, error: null });
    try {
      const results = await api.jobs.getResults(id);
      set({ currentResults: results, isLoading: false });
      return results;
    } catch (error: any) {
      set({ error: error.message, isLoading: false });
      throw error;
    }
  },

  cancelJob: async (id: string) => {
    try {
      await api.jobs.cancel(id);
      set((state) => ({
        jobs: state.jobs.map((j) =>
          j.job_id === id ? { ...j, status: "cancelled" as const } : j
        ),
      }));
    } catch (error: any) {
      set({ error: error.message });
      throw error;
    }
  },

  setCurrentJob: (job: Job | null) => {
    set({ currentJob: job });
  },

  updateJobProgress: (jobId: string, progress: number) => {
    set((state) => ({
      jobs: state.jobs.map((j) => (j.job_id === jobId ? { ...j, progress } : j)),
      currentJob:
        state.currentJob?.job_id === jobId
          ? { ...state.currentJob, progress }
          : state.currentJob,
    }));
  },
}));
