// Kanad API Client Library - Matches Backend Exactly
// This module provides all API functions for communicating with the FastAPI backend

const API_BASE_URL = process.env.NEXT_PUBLIC_API_URL || "http://localhost:8000/api";
const WS_BASE_URL = process.env.NEXT_PUBLIC_WS_URL || "ws://localhost:8000/api";
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

  // Get access token from localStorage for authenticated requests
  const accessToken = typeof window !== 'undefined' ? localStorage.getItem('kanad_access_token') : null;

  try {
    const headers: Record<string, string> = {
      "Content-Type": "application/json",
      ...(options?.headers as Record<string, string>),
    };

    // Add Authorization header if token exists
    if (accessToken) {
      headers["Authorization"] = `Bearer ${accessToken}`;
    }

    const response = await fetch(`${API_BASE_URL}${endpoint}`, {
      ...options,
      signal: controller.signal,
      headers,
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

      throw new APIError(response.status, error.detail || error.message || error.error || "API Error");
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

// ===== ANALYSIS =====

export async function getBondAnalysis(experimentId: string) {
  return apiCall(`/analysis/${experimentId}/bond-analysis`);
}

export async function getEnergyDecomposition(experimentId: string) {
  return apiCall(`/analysis/${experimentId}/energy-decomposition`);
}

export async function getFullAnalysis(experimentId: string) {
  return apiCall(`/analysis/${experimentId}/full-analysis`);
}

export async function compareExperiments(experimentIds: string[]) {
  return apiCall("/analysis/compare", {
    method: "POST",
    body: JSON.stringify({ experiment_ids: experimentIds }),
  });
}

// ===== ADVANCED ANALYSIS (AnalysisService) =====

export async function getAnalysisProfiles() {
  return apiCall("/analysis/profiles");
}

export async function getAnalysisModules() {
  return apiCall("/analysis/modules");
}

export async function runAdvancedAnalysis(experimentId: string, profileName: string) {
  return apiCall(`/analysis/${experimentId}/run-analysis`, {
    method: "POST",
    body: JSON.stringify({ profile: profileName }),
  });
}

export async function getAnalysisResults(experimentId: string) {
  return apiCall(`/analysis/${experimentId}/results`);
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

export async function getCampaignResults(campaignId: string) {
  return apiCall(`/campaigns/${campaignId}/results`);
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
  // Use dedicated WebSocket URL (bypasses Vercel proxy for production)
  const fullWsUrl = `${WS_BASE_URL}/ws/experiments/${experimentId}`;
  console.log("🔌 Creating WebSocket connection to:", fullWsUrl);
  return new WebSocket(fullWsUrl);
}

// ===== AUTHENTICATION =====

export async function register(data: {
  email: string;
  password: string;
  full_name: string;
  access_key: string;
}) {
  return apiCall("/auth/register", {
    method: "POST",
    body: JSON.stringify(data),
  });
}

export async function verifyOTP(data: {
  email: string;
  otp: string;
}) {
  return apiCall("/auth/verify-email", {
    method: "POST",
    body: JSON.stringify(data),
  });
}

export async function login(data: {
  email: string;
  password: string;
}) {
  return apiCall("/auth/login", {
    method: "POST",
    body: JSON.stringify(data),
  });
}

export async function refreshToken(refreshToken: string) {
  return apiCall("/auth/refresh", {
    method: "POST",
    body: JSON.stringify({ refresh_token: refreshToken }),
  });
}

export async function logout(refreshToken: string) {
  return apiCall("/auth/logout", {
    method: "POST",
    body: JSON.stringify({ refresh_token: refreshToken }),
  });
}

export async function requestPasswordReset(email: string) {
  return apiCall("/auth/request-password-reset", {
    method: "POST",
    body: JSON.stringify({ email }),
  });
}

export async function resetPassword(data: {
  email: string;
  otp: string;
  new_password: string;
}) {
  return apiCall("/auth/reset-password", {
    method: "POST",
    body: JSON.stringify(data),
  });
}

export async function googleAuth(credential: string) {
  return apiCall("/auth/google", {
    method: "POST",
    body: JSON.stringify({ id_token: credential }),
  });
}

export async function getAuthStatus() {
  return apiCall("/auth/status");
}

// ===== USER PROFILE & ACCOUNT MANAGEMENT =====

export async function getUserProfile() {
  return apiCall("/users/me");
}

export async function updateUserProfile(data: {
  full_name?: string;
  avatar_url?: string;
}) {
  return apiCall("/users/me", {
    method: "PUT",
    body: JSON.stringify(data),
  });
}

export async function changePassword(data: {
  current_password: string;
  new_password: string;
}) {
  return apiCall("/users/change-password", {
    method: "POST",
    body: JSON.stringify(data),
  });
}

export async function getUserStatistics() {
  return apiCall("/users/me/statistics");
}

export async function getUserSessions() {
  return apiCall("/users/me/sessions");
}

export async function revokeSession(sessionId: number) {
  return apiCall(`/users/me/sessions/${sessionId}`, {
    method: "DELETE",
  });
}

export async function deleteAccount() {
  return apiCall("/users/me", {
    method: "DELETE",
  });
}

// ===== ADMIN: ACCESS KEY MANAGEMENT =====

export async function createAccessKey(data: {
  description?: string;
  max_uses: number;
  expires_in_days?: number;
}) {
  return apiCall("/admin/keys", {
    method: "POST",
    body: JSON.stringify(data),
  });
}

export async function listAccessKeys(activeOnly: boolean = false) {
  return apiCall(`/admin/keys?active_only=${activeOnly}`);
}

export async function getAccessKey(keyId: number) {
  return apiCall(`/admin/keys/${keyId}`);
}

export async function deactivateAccessKey(keyId: number) {
  return apiCall(`/admin/keys/${keyId}/deactivate`, {
    method: "PUT",
  });
}

export async function deleteAccessKey(keyId: number) {
  return apiCall(`/admin/keys/${keyId}`, {
    method: "DELETE",
  });
}

// ===== ADMIN: USER MANAGEMENT =====

export async function listUsers(params?: {
  role?: string;
  verified_only?: boolean;
  active_only?: boolean;
  skip?: number;
  limit?: number;
}) {
  const query = new URLSearchParams();
  if (params?.role) query.append("role", params.role);
  if (params?.verified_only) query.append("verified_only", "true");
  if (params?.active_only !== undefined) query.append("active_only", params.active_only.toString());
  if (params?.skip) query.append("skip", params.skip.toString());
  if (params?.limit) query.append("limit", params.limit.toString());

  return apiCall(`/admin/users?${query.toString()}`);
}

export async function getUser(userId: number) {
  return apiCall(`/admin/users/${userId}`);
}

export async function updateUser(userId: number, data: {
  role?: string;
  is_active?: boolean;
  is_verified?: boolean;
}) {
  return apiCall(`/admin/users/${userId}`, {
    method: "PUT",
    body: JSON.stringify(data),
  });
}

export async function createUser(data: {
  email: string;
  full_name: string;
  role?: string;
  send_credentials?: boolean;
}) {
  return apiCall("/admin/users", {
    method: "POST",
    body: JSON.stringify(data),
  });
}

export async function deleteUser(userId: number) {
  return apiCall(`/admin/users/${userId}`, {
    method: "DELETE",
  });
}

// ===== ADMIN: LIVE MONITORING =====

export async function getLiveExperiments() {
  return apiCall("/admin/experiments/live");
}

export async function getRecentExperiments(limit: number = 50) {
  return apiCall(`/admin/experiments/recent?limit=${limit}`);
}

// ===== ADMIN: USAGE TRACKING =====

export async function getUsageByBackend(days: number = 30) {
  return apiCall(`/admin/usage/by-backend?days=${days}`);
}

export async function getUsageByUser(days: number = 30, limit: number = 10) {
  return apiCall(`/admin/usage/by-user?days=${days}&limit=${limit}`);
}

export async function getUsageTimeline(days: number = 30) {
  return apiCall(`/admin/usage/timeline?days=${days}`);
}

// ===== ADMIN: SYSTEM STATISTICS =====

export async function getSystemStats() {
  return apiCall("/admin/stats/overview");
}

export async function getGrowthStats(days: number = 30) {
  return apiCall(`/admin/stats/growth?days=${days}`);
}

// ===== HEALTH =====

export async function healthCheck() {
  return fetch(`${API_BASE_URL.replace('/api', '')}/health`).then(r => r.json());
}

export { APIError };
