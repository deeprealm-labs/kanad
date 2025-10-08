// Kanad API Client Library
// This module provides all API functions for communicating with the FastAPI backend

import type {
  CreateExperimentRequest,
  CreateExperimentResponse,
  GetExperimentsRequest,
  GetExperimentsResponse,
  GetExperimentResponse,
  UpdateExperimentRequest,
  GetQueueResponse,
  UpdateQueueItemRequest,
  GetSettingsResponse,
  UpdateSettingsRequest,
  ValidateSmilesRequest,
  ValidateSmilesResponse,
  ExportFormat,
  ApiError,
  Experiment,
  QueuedJob,
  UserSettings,
} from "./types";

// ===== Configuration =====
const API_BASE_URL = process.env.NEXT_PUBLIC_API_URL || "http://localhost:8000";
const API_TIMEOUT = 30000; // 30 seconds
const MAX_RETRIES = 3;
const RETRY_DELAY = 1000; // 1 second

// ===== Utility Functions =====

/**
 * Custom fetch with timeout support
 */
async function fetchWithTimeout(
  url: string,
  options: RequestInit = {},
  timeout = API_TIMEOUT
): Promise<Response> {
  const controller = new AbortController();
  const timeoutId = setTimeout(() => controller.abort(), timeout);

  try {
    const response = await fetch(url, {
      ...options,
      signal: controller.signal,
    });
    clearTimeout(timeoutId);
    return response;
  } catch (error) {
    clearTimeout(timeoutId);
    throw error;
  }
}

/**
 * Retry logic for failed requests
 */
async function fetchWithRetry(
  url: string,
  options: RequestInit = {},
  retries = MAX_RETRIES
): Promise<Response> {
  try {
    return await fetchWithTimeout(url, options);
  } catch (error: any) {
    if (retries > 0 && (error.name === "AbortError" || error.name === "TypeError")) {
      await new Promise((resolve) => setTimeout(resolve, RETRY_DELAY));
      return fetchWithRetry(url, options, retries - 1);
    }
    throw error;
  }
}

/**
 * Handle API response and errors
 */
async function handleResponse<T>(response: Response): Promise<T> {
  if (!response.ok) {
    let errorData: ApiError;
    try {
      errorData = await response.json();
    } catch {
      errorData = {
        error: "Unknown Error",
        message: response.statusText || "An unexpected error occurred",
        statusCode: response.status,
      };
    }
    throw new ApiError(
      errorData.message || "API request failed",
      errorData.statusCode,
      errorData.details
    );
  }

  // Handle empty responses (204 No Content)
  if (response.status === 204) {
    return {} as T;
  }

  try {
    return await response.json();
  } catch {
    throw new ApiError("Invalid JSON response from server", response.status);
  }
}

/**
 * Get default headers for API requests
 */
function getHeaders(): HeadersInit {
  const headers: HeadersInit = {
    "Content-Type": "application/json",
    Accept: "application/json",
  };

  // Add authentication token if available
  if (typeof window !== "undefined") {
    const token = localStorage.getItem("kanad_auth_token");
    if (token) {
      headers["Authorization"] = `Bearer ${token}`;
    }
  }

  return headers;
}

// ===== Custom Error Class =====
export class ApiError extends Error {
  constructor(
    message: string,
    public statusCode: number = 500,
    public details?: any
  ) {
    super(message);
    this.name = "ApiError";
  }
}

// ===== API Client Functions =====

// ----- Experiments -----

/**
 * Create a new experiment
 * @param request - Experiment configuration
 * @returns Created experiment details
 */
export async function createExperiment(
  request: CreateExperimentRequest
): Promise<CreateExperimentResponse> {
  // Debug logging
  console.log("createExperiment called with:", JSON.stringify(request, null, 2));

  // Ensure backendSettings exists
  if (!request.backendSettings) {
    console.error("Missing backendSettings in request:", request);
    throw new Error("backendSettings is required. Received: " + JSON.stringify(request));
  }

  // Transform frontend request to backend schema
  const backendRequest = {
    name: request.molecule?.smiles || "Custom Experiment",
    molecule: {
      smiles: request.molecule?.smiles,
      atoms: request.molecule?.atoms,
      basis: request.molecule?.basis || "sto-3g",
      charge: request.molecule?.charge || 0,
      multiplicity: request.molecule?.multiplicity || 1,
    },
    configuration: {
      method: request.backendSettings.method || "VQE",
      ansatz: request.backendSettings.ansatz || "ucc",
      mapper: request.backendSettings.mapper || "jordan_wigner",
      optimizer: request.backendSettings.optimizer || "SLSQP",
      backend: request.backendSettings.backend || "classical",
      backend_name: request.backendSettings.backendName,
    },
    execute_immediately: request.executeNow || false,
  };

  const response = await fetchWithRetry(`${API_BASE_URL}/api/v1/experiments`, {
    method: "POST",
    headers: getHeaders(),
    body: JSON.stringify(backendRequest),
  });

  const result: any = await handleResponse(response);

  // Transform backend response to frontend schema
  return {
    experimentId: result.id?.toString() || "",
    status: result.status || "queued",
    message: result.message || "Experiment created",
  };
}

/**
 * Get list of experiments with optional filters
 * @param params - Query parameters for filtering experiments
 * @returns List of experiments
 */
export async function getExperiments(
  params: GetExperimentsRequest = {}
): Promise<GetExperimentsResponse> {
  const queryParams = new URLSearchParams();

  if (params.status) queryParams.append("status", params.status);
  if (params.limit) queryParams.append("limit", params.limit.toString());
  if (params.offset) queryParams.append("offset", params.offset.toString());
  if (params.search) queryParams.append("search", params.search);

  const url = `${API_BASE_URL}/api/v1/experiments?${queryParams.toString()}`;
  const response = await fetchWithRetry(url, {
    method: "GET",
    headers: getHeaders(),
  });

  const result: any = await handleResponse(response);

  // Transform backend response to frontend schema
  return {
    experiments: result.experiments.map((exp: any) => ({
      id: exp.id?.toString() || "",
      name: exp.name || exp.smiles || "Unknown",
      smiles: exp.smiles,
      molecule: exp.molecule_data,
      configuration: exp.configuration,
      status: exp.status,
      energy: exp.energy,
      hf_energy: exp.hf_energy,
      correlation_energy: exp.correlation_energy,
      results: exp.results,
      convergenceData: exp.convergence_data,
      error: exp.error_message,
      timestamp: exp.created_at,
      startedAt: exp.started_at,
      completedAt: exp.completed_at,
      method: exp.configuration?.method || "VQE",
      backend: exp.configuration?.backend || "classical",
    })),
    total: result.total || 0,
    offset: result.offset || 0,
    limit: result.limit || 100,
  };
}

/**
 * Get a single experiment by ID
 * @param experimentId - Experiment ID
 * @returns Experiment details
 */
export async function getExperiment(
  experimentId: string
): Promise<GetExperimentResponse> {
  const response = await fetchWithRetry(
    `${API_BASE_URL}/api/v1/experiments/${experimentId}`,
    {
      method: "GET",
      headers: getHeaders(),
    }
  );

  const exp: any = await handleResponse(response);

  // Transform backend response to frontend schema
  return {
    experiment: {
      id: exp.id?.toString() || "",
      name: exp.name || exp.smiles || "Unknown",
      smiles: exp.smiles,
      molecule: exp.molecule_data,
      configuration: exp.configuration,
      status: exp.status,
      energy: exp.energy,
      hf_energy: exp.hf_energy,
      correlation_energy: exp.correlation_energy,
      results: exp.results,
      convergenceData: exp.convergence_data,
      error: exp.error_message,
      timestamp: exp.created_at,
      startedAt: exp.started_at,
      completedAt: exp.completed_at,
      method: exp.configuration?.method || "VQE",
      backend: exp.configuration?.backend || "classical",
    },
  };
}

/**
 * Update an experiment (e.g., cancel a running experiment)
 * @param experimentId - Experiment ID
 * @param request - Update data
 */
export async function updateExperiment(
  experimentId: string,
  request: UpdateExperimentRequest
): Promise<void> {
  const response = await fetchWithRetry(
    `${API_BASE_URL}/api/v1/experiments/${experimentId}`,
    {
      method: "PATCH",
      headers: getHeaders(),
      body: JSON.stringify(request),
    }
  );

  return handleResponse<void>(response);
}

/**
 * Cancel a running or queued experiment
 * @param experimentId - Experiment ID to cancel
 * @returns Cancellation confirmation
 */
export async function cancelExperiment(
  experimentId: string
): Promise<{ status: string; message: string }> {
  const response = await fetchWithRetry(
    `${API_BASE_URL}/api/v1/experiments/${experimentId}/cancel`,
    {
      method: "PATCH",
      headers: getHeaders(),
    }
  );

  return handleResponse<{ status: string; message: string }>(response);
}

/**
 * Delete an experiment
 * @param experimentId - Experiment ID
 */
export async function deleteExperiment(experimentId: string): Promise<void> {
  const response = await fetchWithRetry(
    `${API_BASE_URL}/api/v1/experiments/${experimentId}`,
    {
      method: "DELETE",
      headers: getHeaders(),
    }
  );

  return handleResponse<void>(response);
}

/**
 * Export experiment results
 * @param experimentId - Experiment ID
 * @param format - Export format (json or csv)
 * @returns Blob containing the exported data
 */
export async function exportExperiment(
  experimentId: string,
  format: ExportFormat = "json"
): Promise<Blob> {
  const response = await fetchWithRetry(
    `${API_BASE_URL}/api/v1/experiments/${experimentId}/export?format=${format}`,
    {
      method: "GET",
      headers: {
        ...getHeaders(),
        Accept: format === "json" ? "application/json" : "text/csv",
      },
    }
  );

  if (!response.ok) {
    const errorData = await response.json();
    throw new ApiError(errorData.message, response.status);
  }

  return await response.blob();
}

/**
 * Download experiment export as file
 * @param experimentId - Experiment ID
 * @param format - Export format
 * @param filename - Optional filename
 */
export async function downloadExperiment(
  experimentId: string,
  format: ExportFormat = "json",
  filename?: string
): Promise<void> {
  const blob = await exportExperiment(experimentId, format);
  const url = window.URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename || `experiment-${experimentId}.${format}`;
  document.body.appendChild(a);
  a.click();
  window.URL.revokeObjectURL(url);
  document.body.removeChild(a);
}

/**
 * Get circuit data for an experiment
 * @param experimentId - Experiment ID
 * @returns Circuit visualization data
 */
export async function getExperimentCircuit(
  experimentId: string
): Promise<any> {
  const response = await fetchWithRetry(
    `${API_BASE_URL}/api/v1/experiments/${experimentId}/circuit`,
    {
      method: "GET",
      headers: getHeaders(),
    }
  );

  return handleResponse(response);
}

/**
 * Get detailed report for an experiment
 * @param experimentId - Experiment ID
 * @returns Detailed experiment report
 */
export async function getExperimentReport(
  experimentId: string
): Promise<any> {
  const response = await fetchWithRetry(
    `${API_BASE_URL}/api/v1/experiments/${experimentId}/report`,
    {
      method: "GET",
      headers: getHeaders(),
    }
  );

  return handleResponse(response);
}

/**
 * Get available solver/ansatz/mapper options from backend
 * @returns Available configuration options
 */
export async function getConfigurationOptions(): Promise<{
  methods: string[];
  ansatze: string[];
  mappers: string[];
  hamiltonians: string[];
  optimizers: string[];
  basis_sets: string[];
  backends: string[];
}> {
  try {
    const response = await fetchWithRetry(`${API_BASE_URL}/api/v1/info`, {
      method: "GET",
      headers: getHeaders(),
    });

    const data: any = await handleResponse(response);
    return {
      methods: data.capabilities?.methods || ["VQE", "HF", "SQD", "QPE"],
      ansatze: data.capabilities?.ansatze || ["ucc", "uccsd", "hardware_efficient", "governance", "hea"],
      mappers: data.capabilities?.mappers || ["jordan_wigner", "bravyi_kitaev", "hybrid_orbital", "parity"],
      hamiltonians: data.capabilities?.hamiltonians || ["covalent", "ionic", "metallic", "custom"],
      optimizers: data.capabilities?.optimizers || ["SLSQP", "COBYLA", "L-BFGS-B", "ADAM", "POWELL"],
      basis_sets: data.capabilities?.basis_sets || ["sto-3g", "6-31g", "6-31g*", "6-31g**", "cc-pvdz"],
      backends: data.capabilities?.backends || ["classical", "ibm_quantum", "bluequbit"],
    };
  } catch (error) {
    console.error("Failed to fetch configuration options:", error);
    // Return fallback options
    return {
      methods: ["VQE", "HF", "SQD", "QPE"],
      ansatze: ["ucc", "uccsd", "hardware_efficient", "governance", "hea"],
      mappers: ["jordan_wigner", "bravyi_kitaev", "hybrid_orbital", "parity"],
      hamiltonians: ["covalent", "ionic", "metallic", "custom"],
      optimizers: ["SLSQP", "COBYLA", "L-BFGS-B", "ADAM", "POWELL"],
      basis_sets: ["sto-3g", "6-31g", "6-31g*", "6-31g**", "cc-pvdz", "cc-pvtz"],
      backends: ["classical", "ibm_quantum", "bluequbit"],
    };
  }
}

/**
 * Get queue statistics
 * @returns Queue statistics and counts
 */
export async function getQueueStatistics(): Promise<{
  total: number;
  running: number;
  queued: number;
  scheduled: number;
  paused: number;
}> {
  try {
    const response = await fetchWithRetry(`${API_BASE_URL}/api/v1/queue/stats`, {
      method: "GET",
      headers: getHeaders(),
    });

    return handleResponse(response);
  } catch (error) {
    console.error("Failed to fetch queue statistics:", error);
    // Return fallback data
    return {
      total: 0,
      running: 0,
      queued: 0,
      scheduled: 0,
      paused: 0,
    };
  }
}

// ----- Queue Management -----

/**
 * Get the job queue
 * @returns List of queued jobs
 */
export async function getQueue(): Promise<GetQueueResponse> {
  const response = await fetchWithRetry(`${API_BASE_URL}/api/v1/queue`, {
    method: "GET",
    headers: getHeaders(),
  });

  return handleResponse<GetQueueResponse>(response);
}

/**
 * Update a queue item (priority, status, schedule)
 * @param jobId - Job ID
 * @param request - Update data
 */
export async function updateQueueItem(
  jobId: string,
  request: UpdateQueueItemRequest
): Promise<void> {
  const response = await fetchWithRetry(`${API_BASE_URL}/api/v1/queue/${jobId}`, {
    method: "PATCH",
    headers: getHeaders(),
    body: JSON.stringify(request),
  });

  return handleResponse<void>(response);
}

/**
 * Delete a job from the queue
 * @param jobId - Job ID
 */
export async function deleteQueueItem(jobId: string): Promise<void> {
  const response = await fetchWithRetry(`${API_BASE_URL}/api/v1/queue/${jobId}`, {
    method: "DELETE",
    headers: getHeaders(),
  });

  return handleResponse<void>(response);
}

/**
 * Reorder queue items
 * @param jobIds - Array of job IDs in new order
 */
export async function reorderQueue(jobIds: string[]): Promise<void> {
  const response = await fetchWithRetry(`${API_BASE_URL}/api/v1/queue/reorder`, {
    method: "POST",
    headers: getHeaders(),
    body: JSON.stringify({ jobIds }),
  });

  return handleResponse<void>(response);
}

// ----- Settings -----

/**
 * Get user settings
 * @returns User settings
 */
export async function getSettings(): Promise<GetSettingsResponse> {
  const response = await fetchWithRetry(`${API_BASE_URL}/api/v1/settings`, {
    method: "GET",
    headers: getHeaders(),
  });

  return handleResponse<GetSettingsResponse>(response);
}

/**
 * Update user settings
 * @param request - Settings data
 */
export async function updateSettings(
  request: UpdateSettingsRequest
): Promise<void> {
  const response = await fetchWithRetry(`${API_BASE_URL}/api/v1/settings`, {
    method: "PUT",
    headers: getHeaders(),
    body: JSON.stringify(request),
  });

  return handleResponse<void>(response);
}

// ----- SMILES Validation -----

/**
 * Validate a SMILES string
 * @param smiles - SMILES string to validate
 * @returns Validation result
 */
export async function validateSmiles(
  smiles: string
): Promise<ValidateSmilesResponse> {
  const response = await fetchWithRetry(
    `${API_BASE_URL}/api/v1/validate/smiles`,
    {
      method: "POST",
      headers: getHeaders(),
      body: JSON.stringify({ smiles } as ValidateSmilesRequest),
    }
  );

  return handleResponse<ValidateSmilesResponse>(response);
}

// ----- WebSocket Connection -----

/**
 * Create a WebSocket connection for real-time updates
 * @param experimentId - Experiment ID to monitor
 * @returns WebSocket instance
 */
export function createWebSocket(experimentId: string): WebSocket {
  const wsUrl = API_BASE_URL.replace("http", "ws");
  const token = localStorage.getItem("kanad_auth_token");
  const url = `${wsUrl}/ws/experiments/${experimentId}${token ? `?token=${token}` : ""}`;

  return new WebSocket(url);
}

/**
 * Poll for experiment status (alternative to WebSocket)
 * @param experimentId - Experiment ID
 * @param onUpdate - Callback for updates
 * @param interval - Polling interval in milliseconds
 * @returns Stop function
 */
export function pollExperimentStatus(
  experimentId: string,
  onUpdate: (experiment: Experiment) => void,
  interval = 2000
): () => void {
  let isActive = true;

  const poll = async () => {
    while (isActive) {
      try {
        const { experiment } = await getExperiment(experimentId);
        onUpdate(experiment);

        // Stop polling if experiment is complete or failed
        if (experiment.status === "completed" || experiment.status === "failed") {
          isActive = false;
          break;
        }
      } catch (error) {
        console.error("Error polling experiment status:", error);
      }

      await new Promise((resolve) => setTimeout(resolve, interval));
    }
  };

  poll();

  return () => {
    isActive = false;
  };
}

// ----- Cloud Credentials -----

/**
 * Save IBM Quantum credentials
 * @param crn - IBM Cloud Resource Name
 * @param apiKey - IBM Quantum API key
 */
export async function saveIBMCredentials(
  crn: string,
  apiKey: string
): Promise<{ provider: string; configured: boolean; message: string }> {
  const response = await fetchWithRetry(
    `${API_BASE_URL}/api/v1/cloud-credentials/ibm`,
    {
      method: "POST",
      headers: getHeaders(),
      body: JSON.stringify({ crn, api_key: apiKey }),
    }
  );

  return handleResponse(response);
}

/**
 * Save BlueQubit credentials
 * @param apiToken - BlueQubit API token
 */
export async function saveBlueQubitCredentials(
  apiToken: string
): Promise<{ provider: string; configured: boolean; message: string }> {
  const response = await fetchWithRetry(
    `${API_BASE_URL}/api/v1/cloud-credentials/bluequbit`,
    {
      method: "POST",
      headers: getHeaders(),
      body: JSON.stringify({ api_token: apiToken }),
    }
  );

  return handleResponse(response);
}

/**
 * Get cloud credentials status
 * @returns Status of IBM and BlueQubit credentials
 */
export async function getCloudCredentialsStatus(): Promise<{
  ibm: { configured: boolean };
  bluequbit: { configured: boolean };
}> {
  const response = await fetchWithRetry(
    `${API_BASE_URL}/api/v1/cloud-credentials/status`,
    {
      method: "GET",
      headers: getHeaders(),
    }
  );

  return handleResponse(response);
}

/**
 * Delete IBM Quantum credentials
 */
export async function deleteIBMCredentials(): Promise<{ message: string }> {
  const response = await fetchWithRetry(
    `${API_BASE_URL}/api/v1/cloud-credentials/ibm`,
    {
      method: "DELETE",
      headers: getHeaders(),
    }
  );

  return handleResponse(response);
}

/**
 * Delete BlueQubit credentials
 */
export async function deleteBlueQubitCredentials(): Promise<{ message: string }> {
  const response = await fetchWithRetry(
    `${API_BASE_URL}/api/v1/cloud-credentials/bluequbit`,
    {
      method: "DELETE",
      headers: getHeaders(),
    }
  );

  return handleResponse(response);
}

// ----- Health Check -----

/**
 * Check API health
 * @returns Health status
 */
export async function healthCheck(): Promise<{ status: string; message: string }> {
  try {
    const response = await fetchWithTimeout(`${API_BASE_URL}/health`, {
      method: "GET",
    }, 5000);

    return handleResponse(response);
  } catch (error) {
    throw new ApiError("API is unreachable", 503);
  }
}

// ----- Migration Helper -----

/**
 * Migrate localStorage data to backend
 * This helper function migrates existing localStorage data to the API
 */
export async function migrateLocalStorageToAPI(): Promise<{
  experiments: number;
  queue: number;
  settings: boolean;
}> {
  const results = {
    experiments: 0,
    queue: 0,
    settings: false,
  };

  try {
    // Migrate experiments
    const experimentsData = localStorage.getItem("kanad_experiments");
    if (experimentsData) {
      const experiments: Experiment[] = JSON.parse(experimentsData);
      for (const exp of experiments) {
        try {
          await createExperiment({
            molecule: exp.molecule,
            backendSettings: exp.backendSettings,
            analysis: exp.analysis,
            executeNow: false,
          });
          results.experiments++;
        } catch (error) {
          console.error("Failed to migrate experiment:", exp.id, error);
        }
      }
    }

    // Migrate queue
    const queueData = localStorage.getItem("kanad_queue");
    if (queueData) {
      const queue: QueuedJob[] = JSON.parse(queueData);
      for (const job of queue) {
        try {
          await createExperiment({
            molecule: job.molecule,
            backendSettings: job.backendSettings,
            analysis: job.analysis,
            executeNow: false,
          });
          results.queue++;
        } catch (error) {
          console.error("Failed to migrate queue item:", job.id, error);
        }
      }
    }

    // Migrate settings
    const settingsData = localStorage.getItem("kanad_settings");
    if (settingsData) {
      try {
        const settings: UserSettings = JSON.parse(settingsData);
        await updateSettings({ settings });
        results.settings = true;
      } catch (error) {
        console.error("Failed to migrate settings:", error);
      }
    }

    return results;
  } catch (error) {
    console.error("Migration failed:", error);
    throw error;
  }
}

// Export all as default for convenience
export default {
  // Experiments
  createExperiment,
  getExperiments,
  getExperiment,
  updateExperiment,
  cancelExperiment,
  deleteExperiment,
  exportExperiment,
  downloadExperiment,
  getExperimentCircuit,
  getExperimentReport,

  // Queue
  getQueue,
  updateQueueItem,
  deleteQueueItem,
  reorderQueue,
  getQueueStatistics,

  // Settings
  getSettings,
  updateSettings,

  // Configuration
  getConfigurationOptions,

  // Validation
  validateSmiles,

  // Cloud Credentials
  saveIBMCredentials,
  saveBlueQubitCredentials,
  getCloudCredentialsStatus,
  deleteIBMCredentials,
  deleteBlueQubitCredentials,

  // WebSocket
  createWebSocket,
  pollExperimentStatus,

  // Health
  healthCheck,

  // Migration
  migrateLocalStorageToAPI,
};
