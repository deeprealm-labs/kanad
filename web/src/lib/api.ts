// Kanad API Client Library - Matches Backend Exactly
// This module provides all API functions for communicating with the FastAPI backend

const API_BASE_URL = process.env.NEXT_PUBLIC_API_URL || "http://localhost:8000/api";
const API_TIMEOUT = 30000;

class APIError extends Error {
  constructor(public status: number, message: string) {
    super(message);
    this.name = "APIError";
  }
}

async function apiCall<T = any>(endpoint: string, options?: RequestInit): Promise<T> {
  const controller = new AbortController();
  const timeoutId = setTimeout(() => controller.abort(), API_TIMEOUT);

  try {
    const response = await fetch(`${API_BASE_URL}${endpoint}`, {
      ...options,
      signal: controller.signal,
      headers: {
        "Content-Type": "application/json",
        ...options?.headers,
      },
    });

    clearTimeout(timeoutId);

    if (!response.ok) {
      const error = await response.json().catch(() => ({ detail: "Unknown error" }));

      // Handle FastAPI validation errors (422)
      if (response.status === 422 && Array.isArray(error.detail)) {
        const validationErrors = error.detail.map((e: any) =>
          `${e.loc.join('.')}: ${e.msg}`
        ).join(', ');
        throw new APIError(response.status, validationErrors);
      }

      throw new APIError(response.status, error.detail || error.message || "API Error");
    }

    return response.json();
  } catch (error: any) {
    clearTimeout(timeoutId);
    if (error.name === "AbortError") {
      throw new APIError(408, "Request timeout");
    }
    throw error;
  }
}

// ===== EXPERIMENTS =====

export async function submitExperiment(data: any) {
  return apiCall("/experiments/submit", {
    method: "POST",
    body: JSON.stringify(data),
  });
}

export async function listExperiments(status?: string, limit = 50, offset = 0) {
  const params = new URLSearchParams();
  if (status) params.append("status", status);
  params.append("limit", limit.toString());
  params.append("offset", offset.toString());

  return apiCall(`/experiments/list?${params}`);
}

// Backwards compatible wrapper that accepts both object and positional params
export function getExperiments(statusOrOptions?: string | { status?: string; limit?: number; offset?: number }, limit = 50, offset = 0) {
  if (typeof statusOrOptions === 'object') {
    return listExperiments(statusOrOptions?.status, statusOrOptions?.limit || 50, statusOrOptions?.offset || 0);
  }
  return listExperiments(statusOrOptions, limit, offset);
}

export const createExperiment = submitExperiment;

export async function getExperiment(experimentId: string) {
  return apiCall(`/experiments/${experimentId}`);
}

export async function getExperimentStatus(experimentId: string) {
  return apiCall(`/experiments/${experimentId}/status`);
}

export async function getExperimentResults(experimentId: string) {
  return apiCall(`/experiments/${experimentId}/results`);
}

export async function getExperimentReport(experimentId: string) {
  // Report endpoint not yet implemented on backend
  // Return experiment data as report for now
  const experiment = await getExperiment(experimentId);
  return {
    experiment: experiment.experiment || experiment,
  };
}

export async function cancelExperiment(experimentId: string) {
  return apiCall(`/experiments/${experimentId}/cancel`, { method: "POST" });
}

export async function deleteExperiment(experimentId: string) {
  return apiCall(`/experiments/${experimentId}`, { method: "DELETE" });
}

export async function downloadExperiment(experimentId: string, format: string) {
  // Note: Export endpoint not yet implemented on backend
  // For now, fetch the experiment data and create a download
  const experiment = await getExperiment(experimentId);
  const dataStr = JSON.stringify(experiment, null, 2);
  const blob = new Blob([dataStr], { type: 'application/json' });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = `experiment-${experimentId}.${format}`;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(url);
}

export function pollExperimentStatus(
  experimentId: string,
  callback: (experiment: any) => void,
  interval: number = 2000
): () => void {
  const poll = async () => {
    try {
      const response = await getExperiment(experimentId);
      // Extract the experiment object from the response
      const experiment = response.experiment || response;
      callback(experiment);
    } catch (error) {
      console.error("Polling error:", error);
    }
  };

  poll(); // Initial call
  const intervalId = setInterval(poll, interval);

  // Return cleanup function
  return () => clearInterval(intervalId);
}

export async function getExperimentCircuit(experimentId: string) {
  // Note: Circuit endpoint not yet implemented on backend
  // Return mock circuit data for now
  return {
    gates: [],
    depth: 0,
    qubits: 0,
  };
}

export async function getExperimentStatistics() {
  // Get ALL experiments (no limit) to calculate accurate stats
  const response = await getExperiments({ limit: 1000 }); // Large limit to get all
  const experiments = response.experiments || [];

  return {
    total: experiments.length,
    completed: experiments.filter((e: any) => e.status === "completed").length,
    running: experiments.filter((e: any) => e.status === "running").length,
    failed: experiments.filter((e: any) => e.status === "failed").length,
    pending: experiments.filter((e: any) => e.status === "pending").length,
    queued: experiments.filter((e: any) => e.status === "queued").length,
  };
}

// ===== JOBS =====

export async function listJobs(status?: string) {
  const params = status ? `?status=${status}` : "";
  return apiCall(`/jobs/list${params}`);
}

export async function getQueue() {
  // Get all pending/queued EXPERIMENTS (not jobs)
  const response = await getExperiments({ limit: 1000 });
  const experiments = response.experiments || [];

  // Filter for queued/pending experiments (not yet executed, not part of a campaign)
  const queue = experiments
    .filter((e: any) =>
      (e.status === "pending" || e.status === "queued") &&
      !e.campaign_id
    )
    .map((e: any) => ({
      id: e.id,
      name: `${e.molecule?.atoms?.[0]?.symbol || "Custom"} Molecule`,
      molecule: e.molecule,
      method: e.method || "VQE",
      backend: e.backend || "classical",
      status: e.status,
      createdAt: e.created_at,
    }));

  return { queue };
}

export async function getQueueStatistics() {
  // Return statistics based on job list
  const jobs = await listJobs();
  const jobList = jobs.jobs || [];

  return {
    total: jobList.length,
    running: jobList.filter((j: any) => j.status === "running").length,
    pending: jobList.filter((j: any) => j.status === "pending").length,
    completed: jobList.filter((j: any) => j.status === "completed").length,
    failed: jobList.filter((j: any) => j.status === "failed").length,
  };
}

export async function getJob(jobId: string) {
  return apiCall(`/jobs/${jobId}`);
}

export async function getJobStatus(jobId: string) {
  return apiCall(`/jobs/${jobId}/status`);
}

export async function cancelJob(jobId: string) {
  return apiCall(`/jobs/${jobId}/cancel`, { method: "POST" });
}

export async function deleteJob(jobId: string) {
  return apiCall(`/jobs/${jobId}`, { method: "DELETE" });
}

// ===== CAMPAIGNS =====

export async function createCampaign(data: {
  name: string;
  description?: string;
  experiment_ids: string[];
}) {
  return apiCall("/campaigns/create", {
    method: "POST",
    body: JSON.stringify(data),
  });
}

export async function executeCampaign(campaignId: string) {
  return apiCall(`/campaigns/${campaignId}/execute`, {
    method: "POST",
  });
}

export async function getCampaign(campaignId: string) {
  return apiCall(`/campaigns/${campaignId}`);
}

export async function listCampaigns(limit: number = 100, offset: number = 0) {
  return apiCall(`/campaigns/list?limit=${limit}&offset=${offset}`);
}

export async function cancelCampaign(campaignId: string) {
  return apiCall(`/campaigns/${campaignId}/cancel`, {
    method: "POST",
  });
}

export function pollCampaignProgress(
  campaignId: string,
  callback: (campaign: any) => void,
  interval: number = 2000
): () => void {
  const poll = async () => {
    try {
      const response = await getCampaign(campaignId);
      callback(response);
    } catch (error) {
      console.error("Campaign polling error:", error);
    }
  };

  poll(); // Initial call
  const intervalId = setInterval(poll, interval);

  return () => clearInterval(intervalId);
}

// ===== MOLECULES =====

export async function createMolecule(data: any) {
  return apiCall("/molecules/create", {
    method: "POST",
    body: JSON.stringify(data),
  });
}

export async function validateSmiles(smiles: string) {
  return apiCall("/molecules/validate-smiles", {
    method: "POST",
    body: JSON.stringify({ smiles }),
  });
}

export async function getBondInfo(atom1: string, atom2: string) {
  return apiCall(`/molecules/bond-info?atom_1=${atom1}&atom_2=${atom2}`, {
    method: "POST",
  });
}

// ===== SETTINGS =====

export async function getSettings() {
  const response = await apiCall("/settings/defaults");
  return response.settings; // Extract settings from response
}

export async function updateSettings(settings: any) {
  return apiCall("/settings/defaults", {
    method: "PUT",
    body: JSON.stringify({ settings }),
  });
}

// ===== LIBRARY =====

export async function getMoleculeLibrary() {
  return apiCall("/library/molecules");
}

export async function getLibraryMolecule(moleculeId: string) {
  return apiCall(`/library/molecules/${moleculeId}`);
}

// ===== CLOUD =====

export async function getBackends() {
  return apiCall("/cloud/backends");
}

export async function storeCredentials(provider: string, credentials: any) {
  return apiCall("/cloud/credentials", {
    method: "POST",
    body: JSON.stringify({ provider, credentials }),
  });
}

export async function getCredentials(provider: string) {
  return apiCall(`/cloud/credentials/${provider}`);
}

export async function getCloudCredentialsStatus() {
  return apiCall("/cloud/credentials/status");
}

export async function saveIBMCredentials(crn: string, apiKey: string) {
  return apiCall("/cloud/credentials", {
    method: "POST",
    body: JSON.stringify({
      provider: "ibm",
      credentials: {
        crn: crn,
        api_token: apiKey
      }
    })
  });
}

export async function saveBlueQubitCredentials(token: string) {
  return apiCall("/cloud/credentials", {
    method: "POST",
    body: JSON.stringify({
      provider: "bluequbit",
      credentials: {
        api_token: token
      }
    })
  });
}

export async function deleteIBMCredentials() {
  return apiCall("/cloud/credentials/ibm", {
    method: "DELETE"
  });
}

export async function deleteBlueQubitCredentials() {
  return apiCall("/cloud/credentials/bluequbit", {
    method: "DELETE"
  });
}

export async function getCloudBackends() {
  return apiCall("/cloud/backends");
}

// ===== CONFIGURATION =====

export async function getConfigurationOptions() {
  return apiCall("/configuration/options");
}

// ===== CIRCUITS =====

export async function getCircuitPreview(molecule: any, configuration: any) {
  return apiCall("/circuits/preview", {
    method: "POST",
    body: JSON.stringify({ molecule, configuration }),
  });
}

// ===== WEBSOCKET =====

export function createWebSocket(experimentId: string): WebSocket {
  // Convert HTTP URL to WebSocket URL
  const wsUrl = API_BASE_URL.replace('http://', 'ws://').replace('https://', 'wss://');
  const fullWsUrl = `${wsUrl}/ws/experiments/${experimentId}`;
  console.log("🔌 Creating WebSocket connection to:", fullWsUrl);
  return new WebSocket(fullWsUrl);
}

// ===== HEALTH =====

export async function healthCheck() {
  return fetch(`${API_BASE_URL.replace('/api', '')}/health`).then(r => r.json());
}

export { APIError };
