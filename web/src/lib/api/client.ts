// API Client for Kanad Backend

const API_URL = process.env.NEXT_PUBLIC_API_URL || "http://localhost:8000/api";

export class APIError extends Error {
  constructor(public status: number, message: string) {
    super(message);
    this.name = "APIError";
  }
}

export async function apiCall<T = any>(
  endpoint: string,
  options?: RequestInit
): Promise<T> {
  const token = typeof window !== "undefined" ? localStorage.getItem("auth_token") : null;

  const headers: Record<string, string> = {
    "Content-Type": "application/json",
  };

  if (token) {
    headers.Authorization = `Bearer ${token}`;
  }

  // Merge with provided headers
  if (options?.headers) {
    const optHeaders = options.headers as Record<string, string>;
    Object.assign(headers, optHeaders);
  }

  const response = await fetch(`${API_URL}${endpoint}`, {
    ...options,
    headers,
  });

  if (!response.ok) {
    const error = await response.json().catch(() => ({ detail: "Unknown error" }));
    throw new APIError(response.status, error.detail || error.message || "API Error");
  }

  return response.json();
}

export const api = {
  // Auth
  auth: {
    register: (data: { email: string; password: string; name?: string }) =>
      apiCall("/auth/register", {
        method: "POST",
        body: JSON.stringify(data),
      }),
    login: (email: string, password: string) =>
      apiCall("/auth/login", {
        method: "POST",
        body: JSON.stringify({ email, password }),
      }),
    refresh: () => apiCall("/auth/refresh", { method: "POST" }),
  },

  // Molecules
  molecules: {
    create: (data: any) =>
      apiCall("/molecules/create", {
        method: "POST",
        body: JSON.stringify(data),
      }),
    list: () => apiCall("/molecules"),
    get: (id: string) => apiCall(`/molecules/${id}`),
    delete: (id: string) =>
      apiCall(`/molecules/${id}`, { method: "DELETE" }),
  },

  // Simulations
  simulations: {
    configure: (config: any) =>
      apiCall("/simulations/configure", {
        method: "POST",
        body: JSON.stringify(config),
      }),
    submit: (id: string, data: any) =>
      apiCall(`/simulations/${id}/submit`, {
        method: "POST",
        body: JSON.stringify(data),
      }),
  },

  // Jobs
  jobs: {
    list: () => apiCall("/jobs"),
    getStatus: (id: string) => apiCall(`/jobs/${id}/status`),
    getResults: (id: string) => apiCall(`/jobs/${id}/results`),
    cancel: (id: string) => apiCall(`/jobs/${id}`, { method: "DELETE" }),
    export: (id: string, format: string) =>
      fetch(`${API_URL}/jobs/${id}/export?format=${format}`).then((r) => r.blob()),
  },

  // Library
  library: {
    get: () => apiCall("/library"),
  },

  // Cloud
  cloud: {
    storeCredentials: (provider: string, credentials: any) =>
      apiCall("/cloud/credentials", {
        method: "POST",
        body: JSON.stringify({ provider, credentials }),
      }),
    getBackends: () => apiCall("/cloud/backends"),
  },

  // Settings
  settings: {
    getDefaults: () => apiCall("/settings/defaults"),
    updateDefaults: (settings: any) =>
      apiCall("/settings/defaults", {
        method: "PUT",
        body: JSON.stringify(settings),
      }),
  },

  // User
  user: {
    getProfile: () => apiCall("/user/profile"),
    getHistory: () => apiCall("/user/history"),
  },
};
