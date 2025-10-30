"use client";

import { useState, useEffect } from "react";
import { useAuth } from "@/contexts/AuthContext";
import { useToast } from "@/components/ui/toast";
import * as api from "@/lib/api";
import {
  BarChart3,
  Users,
  Key,
  Activity,
  TrendingUp,
  Shield,
  Clock,
  CheckCircle,
  XCircle,
  Plus,
  Trash2,
  Eye,
  RefreshCw,
  Download,
  Filter,
  Search,
} from "lucide-react";

type Tab = "overview" | "keys" | "users" | "monitoring" | "analytics";

interface AccessKey {
  id: number;
  key: string;
  description?: string;
  max_uses: number;
  used_count: number;
  is_active: boolean;
  expires_at?: string;
  created_by?: number;
  created_at: string;
}

interface User {
  id: number;
  email: string;
  full_name: string;
  role: string;
  is_verified: boolean;
  is_active: boolean;
  has_google_auth: boolean;
  access_key_id?: number;
  created_at: string;
  last_login?: string;
  experiments_count: number;
  campaigns_count: number;
}

interface LiveExperiment {
  id: number;
  name: string;
  user_email: string;
  status: string;
  progress: number;
  backend?: string;
  method?: string;
  started_at?: string;
  running_time_seconds?: number;
}

interface SystemStats {
  users: {
    total: number;
    active: number;
    verified: number;
    new_today: number;
  };
  experiments: {
    total: number;
    running: number;
    completed: number;
    started_today: number;
    success_rate: number;
  };
  campaigns: {
    total: number;
  };
  access_keys: {
    total: number;
    active: number;
  };
}

interface UsageByBackend {
  backend: string;
  experiment_count: number;
  average_progress: number;
}

interface UsageByUser {
  user_email: string;
  user_name: string;
  total_experiments: number;
  completed_experiments: number;
  success_rate: number;
}

// Simple password check - In production, use proper authentication
const ADMIN_PASSWORD = "kanad_admin_2024";

export default function AdminDashboard() {
  const { user } = useAuth();
  const toast = useToast();
  const [activeTab, setActiveTab] = useState<Tab>("overview");
  const [loading, setLoading] = useState(false);

  // Password protection
  const [isAuthenticated, setIsAuthenticated] = useState(false);
  const [password, setPassword] = useState("");
  const [passwordError, setPasswordError] = useState("");

  // Overview data
  const [systemStats, setSystemStats] = useState<SystemStats | null>(null);

  // Access Keys data
  const [accessKeys, setAccessKeys] = useState<AccessKey[]>([]);
  const [showCreateKeyModal, setShowCreateKeyModal] = useState(false);
  const [showActiveKeysOnly, setShowActiveKeysOnly] = useState(false);

  // Users data
  const [users, setUsers] = useState<User[]>([]);
  const [userFilters, setUserFilters] = useState({
    role: "",
    verified_only: false,
    active_only: true,
  });
  const [selectedUser, setSelectedUser] = useState<User | null>(null);

  // Monitoring data
  const [liveExperiments, setLiveExperiments] = useState<LiveExperiment[]>([]);
  const [autoRefresh, setAutoRefresh] = useState(true);

  // Analytics data
  const [usageByBackend, setUsageByBackend] = useState<UsageByBackend[]>([]);
  const [usageByUser, setUsageByUser] = useState<UsageByUser[]>([]);
  const [analyticsDays, setAnalyticsDays] = useState(30);

  // Check session storage for admin auth
  useEffect(() => {
    const adminAuth = sessionStorage.getItem("admin_authenticated");
    if (adminAuth === "true") {
      setIsAuthenticated(true);
    }
  }, []);

  // Check if user is admin role and show helpful message
  useEffect(() => {
    if (user && user.role !== "admin" && isAuthenticated) {
      toast.error("Your account doesn't have admin role. Please logout and login again if you were recently promoted.");
    }
  }, [user, toast, isAuthenticated]);

  // Load data based on active tab - MOVED BEFORE EARLY RETURN
  useEffect(() => {
    if (!isAuthenticated) return;

    const loadData = async () => {
      setLoading(true);
      try {
        switch (activeTab) {
          case "overview":
            await loadOverview();
            break;
          case "keys":
            await loadAccessKeys();
            break;
          case "users":
            await loadUsers();
            break;
          case "monitoring":
            await loadLiveExperiments();
            break;
          case "analytics":
            await loadAnalytics();
            break;
        }
      } catch (error: any) {
        toast.error(error.message || "Failed to load data");
      } finally {
        setLoading(false);
      }
    };

    loadData();
  }, [activeTab, isAuthenticated]);

  // Auto-refresh for monitoring tab - MOVED BEFORE EARLY RETURN
  useEffect(() => {
    if (activeTab === "monitoring" && autoRefresh && isAuthenticated) {
      const interval = setInterval(() => {
        loadLiveExperiments();
      }, 5000); // Refresh every 5 seconds

      return () => clearInterval(interval);
    }
  }, [activeTab, autoRefresh, isAuthenticated]);

  const handlePasswordSubmit = (e: React.FormEvent) => {
    e.preventDefault();
    if (password === ADMIN_PASSWORD) {
      sessionStorage.setItem("admin_authenticated", "true");
      setIsAuthenticated(true);
      setPasswordError("");
      toast.success("Admin access granted");
    } else {
      setPasswordError("Incorrect password");
      setPassword("");
    }
  };

  const loadOverview = async () => {
    const stats = await api.getSystemStats();
    setSystemStats(stats);
  };

  const loadAccessKeys = async () => {
    const keys = await api.listAccessKeys(showActiveKeysOnly);
    setAccessKeys(keys);
  };

  const loadUsers = async () => {
    const userList = await api.listUsers(userFilters);
    setUsers(userList);
  };

  const loadLiveExperiments = async () => {
    const experiments = await api.getLiveExperiments();
    setLiveExperiments(experiments);
  };

  const loadAnalytics = async () => {
    const [backendData, userData] = await Promise.all([
      api.getUsageByBackend(analyticsDays),
      api.getUsageByUser(analyticsDays, 10),
    ]);
    setUsageByBackend(backendData);
    setUsageByUser(userData);
  };

  const formatDuration = (seconds?: number) => {
    if (!seconds) return "N/A";
    const hours = Math.floor(seconds / 3600);
    const minutes = Math.floor((seconds % 3600) / 60);
    const secs = seconds % 60;
    if (hours > 0) return `${hours}h ${minutes}m`;
    if (minutes > 0) return `${minutes}m ${secs}s`;
    return `${secs}s`;
  };

  const formatDate = (dateString?: string) => {
    if (!dateString) return "Never";
    return new Date(dateString).toLocaleString();
  };

  const formatRelativeDate = (dateString: string) => {
    const date = new Date(dateString);
    const now = new Date();
    const diffMs = now.getTime() - date.getTime();
    const diffDays = Math.floor(diffMs / (1000 * 60 * 60 * 24));

    if (diffDays === 0) return "Today";
    if (diffDays === 1) return "Yesterday";
    if (diffDays < 7) return `${diffDays} days ago`;
    if (diffDays < 30) return `${Math.floor(diffDays / 7)} weeks ago`;
    return `${Math.floor(diffDays / 30)} months ago`;
  };

  // Show password prompt if not authenticated
  if (!isAuthenticated) {
    return (
      <div className="min-h-screen bg-background flex items-center justify-center p-4">
        <div className="w-full max-w-md">
          <div className="bg-card border-2 border-border rounded-lg p-8">
            <div className="flex items-center justify-center mb-6">
              <Shield className="w-16 h-16 text-brand-orange" />
            </div>
            <h1 className="text-3xl font-quando font-bold text-center mb-2">
              Admin Panel
            </h1>
            <p className="text-muted-foreground text-center mb-8">
              Enter admin password to continue
            </p>

            <form onSubmit={handlePasswordSubmit} className="space-y-4">
              <div>
                <label className="block text-sm font-quando font-semibold mb-2">
                  Password
                </label>
                <input
                  type="password"
                  value={password}
                  onChange={(e) => {
                    setPassword(e.target.value);
                    setPasswordError("");
                  }}
                  className="w-full px-4 py-3 border-2 border-border rounded-md bg-background focus:border-brand-orange focus:outline-none font-quando"
                  placeholder="Enter admin password"
                  autoFocus
                />
                {passwordError && (
                  <p className="text-red-500 text-sm mt-2">{passwordError}</p>
                )}
              </div>

              <button
                type="submit"
                className="w-full px-4 py-3 bg-brand-orange text-white rounded-md hover:bg-brand-orange/90 transition font-quando font-semibold"
              >
                Access Admin Panel
              </button>
            </form>

            <div className="mt-6 pt-6 border-t border-border">
              <p className="text-xs text-muted-foreground text-center">
                This area is restricted to authorized administrators only.
                <br />
                All access attempts are logged.
              </p>
            </div>
          </div>

          {/* Kanad Branding */}
          <div className="mt-8 text-center">
            <h2 className="font-bietro text-5xl text-brand-orange">kanad</h2>
            <p className="text-xs text-muted-foreground mt-2">
              Quantum Chemistry Platform
            </p>
          </div>
        </div>
      </div>
    );
  }

  // ========== RENDER FUNCTIONS ==========

  const renderOverview = () => {
    if (!systemStats) return <div>Loading statistics...</div>;

    return (
      <div className="space-y-6">
        {/* Stats Grid */}
        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-6">
          {/* Users Card */}
          <div className="bg-card border border-border rounded-lg p-6">
            <div className="flex items-center justify-between mb-4">
              <Users className="w-8 h-8 text-brand-orange" />
              <span className="text-sm text-muted-foreground">
                +{systemStats.users.new_today} today
              </span>
            </div>
            <h3 className="text-2xl font-quando font-bold mb-1">
              {systemStats.users.total}
            </h3>
            <p className="text-muted-foreground">Total Users</p>
            <div className="mt-4 flex gap-4 text-sm">
              <div>
                <span className="text-green-500">{systemStats.users.active}</span>
                <span className="text-muted-foreground"> active</span>
              </div>
              <div>
                <span className="text-blue-500">{systemStats.users.verified}</span>
                <span className="text-muted-foreground"> verified</span>
              </div>
            </div>
          </div>

          {/* Experiments Card */}
          <div className="bg-card border border-border rounded-lg p-6">
            <div className="flex items-center justify-between mb-4">
              <Activity className="w-8 h-8 text-brand-orange" />
              <span className="text-sm text-muted-foreground">
                +{systemStats.experiments.started_today} today
              </span>
            </div>
            <h3 className="text-2xl font-quando font-bold mb-1">
              {systemStats.experiments.total}
            </h3>
            <p className="text-muted-foreground">Total Experiments</p>
            <div className="mt-4 flex gap-4 text-sm">
              <div>
                <span className="text-yellow-500">{systemStats.experiments.running}</span>
                <span className="text-muted-foreground"> running</span>
              </div>
              <div>
                <span className="text-green-500">
                  {systemStats.experiments.success_rate}%
                </span>
                <span className="text-muted-foreground"> success</span>
              </div>
            </div>
          </div>

          {/* Access Keys Card */}
          <div className="bg-card border border-border rounded-lg p-6">
            <div className="flex items-center justify-between mb-4">
              <Key className="w-8 h-8 text-brand-orange" />
              <Shield className="w-5 h-5 text-muted-foreground" />
            </div>
            <h3 className="text-2xl font-quando font-bold mb-1">
              {systemStats.access_keys.total}
            </h3>
            <p className="text-muted-foreground">Access Keys</p>
            <div className="mt-4 text-sm">
              <span className="text-green-500">{systemStats.access_keys.active}</span>
              <span className="text-muted-foreground"> active keys</span>
            </div>
          </div>

          {/* Campaigns Card */}
          <div className="bg-card border border-border rounded-lg p-6">
            <div className="flex items-center justify-between mb-4">
              <TrendingUp className="w-8 h-8 text-brand-orange" />
              <BarChart3 className="w-5 h-5 text-muted-foreground" />
            </div>
            <h3 className="text-2xl font-quando font-bold mb-1">
              {systemStats.campaigns.total}
            </h3>
            <p className="text-muted-foreground">Total Campaigns</p>
          </div>
        </div>

        {/* Quick Actions */}
        <div className="bg-card border border-border rounded-lg p-6">
          <h3 className="text-lg font-quando font-bold mb-4">Quick Actions</h3>
          <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
            <button
              onClick={() => setActiveTab("keys")}
              className="px-4 py-3 border border-border rounded-md hover:bg-accent transition font-quando text-left"
            >
              <Key className="w-5 h-5 mb-2 text-brand-orange" />
              <div className="font-semibold">Manage Access Keys</div>
              <div className="text-sm text-muted-foreground">
                Create and manage early access keys
              </div>
            </button>
            <button
              onClick={() => setActiveTab("users")}
              className="px-4 py-3 border border-border rounded-md hover:bg-accent transition font-quando text-left"
            >
              <Users className="w-5 h-5 mb-2 text-brand-orange" />
              <div className="font-semibold">Manage Users</div>
              <div className="text-sm text-muted-foreground">
                View and manage user accounts
              </div>
            </button>
            <button
              onClick={() => setActiveTab("monitoring")}
              className="px-4 py-3 border border-border rounded-md hover:bg-accent transition font-quanto text-left"
            >
              <Activity className="w-5 h-5 mb-2 text-brand-orange" />
              <div className="font-semibold">Live Monitoring</div>
              <div className="text-sm text-muted-foreground">
                Monitor running experiments in real-time
              </div>
            </button>
          </div>
        </div>
      </div>
    );
  };

  const renderAccessKeys = () => {
    return (
      <div className="space-y-6">
        {/* Header */}
        <div className="flex items-center justify-between">
          <div>
            <h2 className="text-2xl font-quando font-bold">Access Keys Management</h2>
            <p className="text-muted-foreground mt-1">
              Create and manage early access registration keys
            </p>
          </div>
          <button
            onClick={() => setShowCreateKeyModal(true)}
            className="px-4 py-2 bg-brand-orange text-white rounded-md hover:bg-brand-orange/90 transition font-quando flex items-center gap-2"
          >
            <Plus className="w-4 h-4" />
            Create Key
          </button>
        </div>

        {/* Filters */}
        <div className="flex items-center gap-4">
          <label className="flex items-center gap-2 cursor-pointer">
            <input
              type="checkbox"
              checked={showActiveKeysOnly}
              onChange={(e) => {
                setShowActiveKeysOnly(e.target.checked);
                setLoading(true);
                api.listAccessKeys(e.target.checked).then((keys) => {
                  setAccessKeys(keys);
                  setLoading(false);
                });
              }}
              className="rounded border-border"
            />
            <span className="text-sm font-quando">Show active keys only</span>
          </label>
          <button
            onClick={loadAccessKeys}
            className="ml-auto px-3 py-2 border border-border rounded-md hover:bg-accent transition"
            title="Refresh"
          >
            <RefreshCw className="w-4 h-4" />
          </button>
        </div>

        {/* Keys Table */}
        <div className="bg-card border border-border rounded-lg overflow-hidden">
          <div className="overflow-x-auto">
            <table className="w-full">
              <thead className="bg-accent border-b border-border">
                <tr>
                  <th className="px-6 py-3 text-left text-xs font-quando font-semibold uppercase tracking-wider">
                    Key
                  </th>
                  <th className="px-6 py-3 text-left text-xs font-quando font-semibold uppercase tracking-wider">
                    Description
                  </th>
                  <th className="px-6 py-3 text-left text-xs font-quando font-semibold uppercase tracking-wider">
                    Usage
                  </th>
                  <th className="px-6 py-3 text-left text-xs font-quando font-semibold uppercase tracking-wider">
                    Status
                  </th>
                  <th className="px-6 py-3 text-left text-xs font-quando font-semibold uppercase tracking-wider">
                    Expires
                  </th>
                  <th className="px-6 py-3 text-left text-xs font-quando font-semibold uppercase tracking-wider">
                    Created
                  </th>
                  <th className="px-6 py-3 text-right text-xs font-quando font-semibold uppercase tracking-wider">
                    Actions
                  </th>
                </tr>
              </thead>
              <tbody className="divide-y divide-border">
                {accessKeys.length === 0 ? (
                  <tr>
                    <td colSpan={7} className="px-6 py-8 text-center text-muted-foreground">
                      No access keys found
                    </td>
                  </tr>
                ) : (
                  accessKeys.map((key) => (
                    <tr key={key.id} className="hover:bg-accent/50 transition">
                      <td className="px-6 py-4">
                        <code className="text-xs bg-accent px-2 py-1 rounded">
                          {key.key}
                        </code>
                      </td>
                      <td className="px-6 py-4 text-sm">
                        {key.description || <span className="text-muted-foreground italic">No description</span>}
                      </td>
                      <td className="px-6 py-4 text-sm">
                        <div className="flex items-center gap-2">
                          <span className={key.used_count >= key.max_uses ? "text-red-500" : ""}>
                            {key.used_count} / {key.max_uses}
                          </span>
                          <div className="w-20 h-2 bg-border rounded-full overflow-hidden">
                            <div
                              className="h-full bg-brand-orange transition-all"
                              style={{
                                width: `${(key.used_count / key.max_uses) * 100}%`,
                              }}
                            />
                          </div>
                        </div>
                      </td>
                      <td className="px-6 py-4">
                        {key.is_active ? (
                          <span className="inline-flex items-center gap-1 px-2 py-1 bg-green-500/10 text-green-500 rounded text-xs font-semibold">
                            <CheckCircle className="w-3 h-3" />
                            Active
                          </span>
                        ) : (
                          <span className="inline-flex items-center gap-1 px-2 py-1 bg-red-500/10 text-red-500 rounded text-xs font-semibold">
                            <XCircle className="w-3 h-3" />
                            Inactive
                          </span>
                        )}
                      </td>
                      <td className="px-6 py-4 text-sm text-muted-foreground">
                        {key.expires_at ? formatRelativeDate(key.expires_at) : "Never"}
                      </td>
                      <td className="px-6 py-4 text-sm text-muted-foreground">
                        {formatRelativeDate(key.created_at)}
                      </td>
                      <td className="px-6 py-4 text-right">
                        <div className="flex items-center justify-end gap-2">
                          {key.is_active && (
                            <button
                              onClick={async () => {
                                if (confirm("Deactivate this access key?")) {
                                  try {
                                    await api.deactivateAccessKey(key.id);
                                    toast.success("Access key deactivated");
                                    loadAccessKeys();
                                  } catch (error: any) {
                                    toast.error(error.message);
                                  }
                                }
                              }}
                              className="p-2 hover:bg-accent rounded transition"
                              title="Deactivate"
                            >
                              <XCircle className="w-4 h-4 text-orange-500" />
                            </button>
                          )}
                          {key.used_count === 0 && (
                            <button
                              onClick={async () => {
                                if (confirm("Permanently delete this access key?")) {
                                  try {
                                    await api.deleteAccessKey(key.id);
                                    toast.success("Access key deleted");
                                    loadAccessKeys();
                                  } catch (error: any) {
                                    toast.error(error.message);
                                  }
                                }
                              }}
                              className="p-2 hover:bg-accent rounded transition"
                              title="Delete"
                            >
                              <Trash2 className="w-4 h-4 text-red-500" />
                            </button>
                          )}
                        </div>
                      </td>
                    </tr>
                  ))
                )}
              </tbody>
            </table>
          </div>
        </div>

        {/* Create Key Modal */}
        {showCreateKeyModal && (
          <CreateAccessKeyModal
            onClose={() => setShowCreateKeyModal(false)}
            onSuccess={() => {
              setShowCreateKeyModal(false);
              loadAccessKeys();
            }}
          />
        )}
      </div>
    );
  };

  const renderUsers = () => {
    return (
      <div className="space-y-6">
        {/* Header */}
        <div className="flex items-center justify-between">
          <div>
            <h2 className="text-2xl font-quando font-bold">User Management</h2>
            <p className="text-muted-foreground mt-1">
              View and manage all user accounts
            </p>
          </div>
        </div>

        {/* Filters */}
        <div className="flex flex-wrap items-center gap-4">
          <select
            value={userFilters.role}
            onChange={(e) => {
              setUserFilters({ ...userFilters, role: e.target.value });
              setLoading(true);
              api
                .listUsers({ ...userFilters, role: e.target.value || undefined })
                .then((userList) => {
                  setUsers(userList);
                  setLoading(false);
                });
            }}
            className="px-3 py-2 border border-border rounded-md bg-background font-quando"
          >
            <option value="">All Roles</option>
            <option value="admin">Admin</option>
            <option value="user">User</option>
            <option value="viewer">Viewer</option>
          </select>

          <label className="flex items-center gap-2 cursor-pointer">
            <input
              type="checkbox"
              checked={userFilters.verified_only}
              onChange={(e) => {
                const newFilters = { ...userFilters, verified_only: e.target.checked };
                setUserFilters(newFilters);
                setLoading(true);
                api.listUsers(newFilters).then((userList) => {
                  setUsers(userList);
                  setLoading(false);
                });
              }}
              className="rounded border-border"
            />
            <span className="text-sm font-quando">Verified only</span>
          </label>

          <label className="flex items-center gap-2 cursor-pointer">
            <input
              type="checkbox"
              checked={userFilters.active_only}
              onChange={(e) => {
                const newFilters = { ...userFilters, active_only: e.target.checked };
                setUserFilters(newFilters);
                setLoading(true);
                api.listUsers(newFilters).then((userList) => {
                  setUsers(userList);
                  setLoading(false);
                });
              }}
              className="rounded border-border"
            />
            <span className="text-sm font-quando">Active only</span>
          </label>

          <button
            onClick={loadUsers}
            className="ml-auto px-3 py-2 border border-border rounded-md hover:bg-accent transition"
            title="Refresh"
          >
            <RefreshCw className="w-4 h-4" />
          </button>
        </div>

        {/* Users Table */}
        <div className="bg-card border border-border rounded-lg overflow-hidden">
          <div className="overflow-x-auto">
            <table className="w-full">
              <thead className="bg-accent border-b border-border">
                <tr>
                  <th className="px-6 py-3 text-left text-xs font-quando font-semibold uppercase tracking-wider">
                    User
                  </th>
                  <th className="px-6 py-3 text-left text-xs font-quando font-semibold uppercase tracking-wider">
                    Role
                  </th>
                  <th className="px-6 py-3 text-left text-xs font-quando font-semibold uppercase tracking-wider">
                    Status
                  </th>
                  <th className="px-6 py-3 text-left text-xs font-quando font-semibold uppercase tracking-wider">
                    Auth Type
                  </th>
                  <th className="px-6 py-3 text-left text-xs font-quando font-semibold uppercase tracking-wider">
                    Activity
                  </th>
                  <th className="px-6 py-3 text-left text-xs font-quando font-semibold uppercase tracking-wider">
                    Usage
                  </th>
                  <th className="px-6 py-3 text-right text-xs font-quando font-semibold uppercase tracking-wider">
                    Actions
                  </th>
                </tr>
              </thead>
              <tbody className="divide-y divide-border">
                {users.length === 0 ? (
                  <tr>
                    <td colSpan={7} className="px-6 py-8 text-center text-muted-foreground">
                      No users found
                    </td>
                  </tr>
                ) : (
                  users.map((u) => (
                    <tr key={u.id} className="hover:bg-accent/50 transition">
                      <td className="px-6 py-4">
                        <div>
                          <div className="font-semibold">{u.full_name}</div>
                          <div className="text-sm text-muted-foreground">{u.email}</div>
                        </div>
                      </td>
                      <td className="px-6 py-4">
                        <span
                          className={`inline-flex px-2 py-1 rounded text-xs font-semibold uppercase ${
                            u.role === "admin"
                              ? "bg-purple-500/10 text-purple-500"
                              : u.role === "user"
                                ? "bg-blue-500/10 text-blue-500"
                                : "bg-gray-500/10 text-gray-500"
                          }`}
                        >
                          {u.role}
                        </span>
                      </td>
                      <td className="px-6 py-4">
                        <div className="flex flex-col gap-1">
                          {u.is_active ? (
                            <span className="inline-flex items-center gap-1 text-xs text-green-500">
                              <CheckCircle className="w-3 h-3" />
                              Active
                            </span>
                          ) : (
                            <span className="inline-flex items-center gap-1 text-xs text-red-500">
                              <XCircle className="w-3 h-3" />
                              Inactive
                            </span>
                          )}
                          {u.is_verified && (
                            <span className="inline-flex items-center gap-1 text-xs text-blue-500">
                              <CheckCircle className="w-3 h-3" />
                              Verified
                            </span>
                          )}
                        </div>
                      </td>
                      <td className="px-6 py-4 text-sm">
                        {u.has_google_auth ? (
                          <span className="text-blue-500">Google OAuth</span>
                        ) : (
                          <span className="text-muted-foreground">Email/Password</span>
                        )}
                      </td>
                      <td className="px-6 py-4 text-sm text-muted-foreground">
                        <div>Joined: {formatRelativeDate(u.created_at)}</div>
                        {u.last_login && (
                          <div className="text-xs">
                            Last login: {formatRelativeDate(u.last_login)}
                          </div>
                        )}
                      </td>
                      <td className="px-6 py-4 text-sm">
                        <div>{u.experiments_count} experiments</div>
                        <div className="text-xs text-muted-foreground">
                          {u.campaigns_count} campaigns
                        </div>
                      </td>
                      <td className="px-6 py-4 text-right">
                        <button
                          onClick={() => setSelectedUser(u)}
                          className="p-2 hover:bg-accent rounded transition"
                          title="View Details"
                        >
                          <Eye className="w-4 h-4" />
                        </button>
                      </td>
                    </tr>
                  ))
                )}
              </tbody>
            </table>
          </div>
        </div>

        {/* User Details Modal */}
        {selectedUser && (
          <UserDetailsModal
            user={selectedUser}
            onClose={() => setSelectedUser(null)}
            onUpdate={() => {
              setSelectedUser(null);
              loadUsers();
            }}
          />
        )}
      </div>
    );
  };

  const renderMonitoring = () => {
    return (
      <div className="space-y-6">
        {/* Header */}
        <div className="flex items-center justify-between">
          <div>
            <h2 className="text-2xl font-quando font-bold">Live Experiment Monitoring</h2>
            <p className="text-muted-foreground mt-1">
              Real-time monitoring of running quantum experiments
            </p>
          </div>
          <div className="flex items-center gap-4">
            <label className="flex items-center gap-2 cursor-pointer">
              <input
                type="checkbox"
                checked={autoRefresh}
                onChange={(e) => setAutoRefresh(e.target.checked)}
                className="rounded border-border"
              />
              <span className="text-sm font-quando">Auto-refresh (5s)</span>
            </label>
            <button
              onClick={loadLiveExperiments}
              className="px-3 py-2 border border-border rounded-md hover:bg-accent transition"
              title="Refresh Now"
            >
              <RefreshCw className="w-4 h-4" />
            </button>
          </div>
        </div>

        {/* Live Experiments */}
        {liveExperiments.length === 0 ? (
          <div className="bg-card border border-border rounded-lg p-12 text-center">
            <Activity className="w-16 h-16 text-muted-foreground mx-auto mb-4" />
            <h3 className="text-lg font-quando font-bold mb-2">No Running Experiments</h3>
            <p className="text-muted-foreground">
              There are currently no experiments running on the platform
            </p>
          </div>
        ) : (
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
            {liveExperiments.map((exp) => (
              <div key={exp.id} className="bg-card border border-border rounded-lg p-6">
                <div className="flex items-start justify-between mb-4">
                  <div className="flex-1">
                    <h3 className="font-quando font-bold text-lg mb-1">{exp.name}</h3>
                    <p className="text-sm text-muted-foreground">{exp.user_email}</p>
                  </div>
                  <span className="px-3 py-1 bg-yellow-500/10 text-yellow-500 rounded text-xs font-semibold">
                    {exp.status}
                  </span>
                </div>

                {/* Progress Bar */}
                <div className="mb-4">
                  <div className="flex items-center justify-between text-sm mb-2">
                    <span className="text-muted-foreground">Progress</span>
                    <span className="font-semibold">{exp.progress.toFixed(1)}%</span>
                  </div>
                  <div className="w-full h-3 bg-border rounded-full overflow-hidden">
                    <div
                      className="h-full bg-brand-orange transition-all duration-300"
                      style={{ width: `${exp.progress}%` }}
                    />
                  </div>
                </div>

                {/* Details Grid */}
                <div className="grid grid-cols-2 gap-4 text-sm">
                  <div>
                    <span className="text-muted-foreground">Backend:</span>
                    <span className="ml-2 font-semibold">{exp.backend || "N/A"}</span>
                  </div>
                  <div>
                    <span className="text-muted-foreground">Method:</span>
                    <span className="ml-2 font-semibold">{exp.method || "N/A"}</span>
                  </div>
                  <div>
                    <span className="text-muted-foreground">Started:</span>
                    <span className="ml-2">{exp.started_at ? formatRelativeDate(exp.started_at) : "N/A"}</span>
                  </div>
                  <div>
                    <span className="text-muted-foreground">Runtime:</span>
                    <span className="ml-2">{formatDuration(exp.running_time_seconds)}</span>
                  </div>
                </div>

                {/* Experiment ID */}
                <div className="mt-4 pt-4 border-t border-border">
                  <span className="text-xs text-muted-foreground">Experiment ID: </span>
                  <code className="text-xs bg-accent px-2 py-1 rounded">{exp.id}</code>
                </div>
              </div>
            ))}
          </div>
        )}
      </div>
    );
  };

  const renderAnalytics = () => {
    return (
      <div className="space-y-6">
        {/* Header */}
        <div className="flex items-center justify-between">
          <div>
            <h2 className="text-2xl font-quando font-bold">Usage Analytics</h2>
            <p className="text-muted-foreground mt-1">
              Platform usage statistics and trends
            </p>
          </div>
          <select
            value={analyticsDays}
            onChange={(e) => {
              setAnalyticsDays(Number(e.target.value));
              setLoading(true);
              Promise.all([
                api.getUsageByBackend(Number(e.target.value)),
                api.getUsageByUser(Number(e.target.value), 10),
              ]).then(([backendData, userData]) => {
                setUsageByBackend(backendData);
                setUsageByUser(userData);
                setLoading(false);
              });
            }}
            className="px-3 py-2 border border-border rounded-md bg-background font-quando"
          >
            <option value={7}>Last 7 days</option>
            <option value={30}>Last 30 days</option>
            <option value={90}>Last 90 days</option>
            <option value={365}>Last year</option>
          </select>
        </div>

        {/* Backend Usage */}
        <div className="bg-card border border-border rounded-lg p-6">
          <h3 className="text-lg font-quando font-bold mb-4">Usage by Backend</h3>
          {usageByBackend.length === 0 ? (
            <p className="text-muted-foreground text-center py-8">No data available</p>
          ) : (
            <div className="space-y-4">
              {usageByBackend.map((item, index) => (
                <div key={index}>
                  <div className="flex items-center justify-between text-sm mb-2">
                    <span className="font-semibold">{item.backend}</span>
                    <div className="flex items-center gap-4 text-muted-foreground">
                      <span>{item.experiment_count} experiments</span>
                      <span>{item.average_progress.toFixed(1)}% avg progress</span>
                    </div>
                  </div>
                  <div className="w-full h-3 bg-border rounded-full overflow-hidden">
                    <div
                      className="h-full bg-brand-orange"
                      style={{
                        width: `${
                          (item.experiment_count /
                            Math.max(...usageByBackend.map((b) => b.experiment_count))) *
                          100
                        }%`,
                      }}
                    />
                  </div>
                </div>
              ))}
            </div>
          )}
        </div>

        {/* Top Users */}
        <div className="bg-card border border-border rounded-lg p-6">
          <h3 className="text-lg font-quando font-bold mb-4">Top Users by Activity</h3>
          {usageByUser.length === 0 ? (
            <p className="text-muted-foreground text-center py-8">No data available</p>
          ) : (
            <div className="overflow-x-auto">
              <table className="w-full">
                <thead className="border-b border-border">
                  <tr>
                    <th className="px-4 py-3 text-left text-xs font-quando font-semibold uppercase">
                      User
                    </th>
                    <th className="px-4 py-3 text-right text-xs font-quando font-semibold uppercase">
                      Total
                    </th>
                    <th className="px-4 py-3 text-right text-xs font-quando font-semibold uppercase">
                      Completed
                    </th>
                    <th className="px-4 py-3 text-right text-xs font-quando font-semibold uppercase">
                      Success Rate
                    </th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-border">
                  {usageByUser.map((item, index) => (
                    <tr key={index}>
                      <td className="px-4 py-3">
                        <div>
                          <div className="font-semibold">{item.user_name}</div>
                          <div className="text-xs text-muted-foreground">
                            {item.user_email}
                          </div>
                        </div>
                      </td>
                      <td className="px-4 py-3 text-right font-semibold">
                        {item.total_experiments}
                      </td>
                      <td className="px-4 py-3 text-right text-green-500">
                        {item.completed_experiments}
                      </td>
                      <td className="px-4 py-3 text-right">
                        <span
                          className={`font-semibold ${
                            item.success_rate >= 80
                              ? "text-green-500"
                              : item.success_rate >= 50
                                ? "text-yellow-500"
                                : "text-red-500"
                          }`}
                        >
                          {item.success_rate.toFixed(1)}%
                        </span>
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          )}
        </div>
      </div>
    );
  };

  // All authentication checks passed, show admin dashboard
  return (
    <div className="min-h-screen bg-background">
      {/* Page Header */}
      <div className="border-b border-border bg-card">
        <div className="max-w-7xl mx-auto px-6 py-6">
          <div className="flex items-center gap-3 mb-2">
            <Shield className="w-8 h-8 text-brand-orange" />
            <h1 className="text-3xl font-quando font-bold">Admin Dashboard</h1>
          </div>
          <p className="text-muted-foreground">
            Platform administration and monitoring tools
          </p>
        </div>
      </div>

      {/* Tabs */}
      <div className="border-b border-border bg-card sticky top-0 z-10">
        <div className="max-w-7xl mx-auto px-6">
          <div className="flex gap-1">
            <button
              onClick={() => setActiveTab("overview")}
              className={`px-6 py-3 font-quando font-semibold transition border-b-2 ${
                activeTab === "overview"
                  ? "border-brand-orange text-brand-orange"
                  : "border-transparent text-muted-foreground hover:text-foreground"
              }`}
            >
              <BarChart3 className="w-4 h-4 inline mr-2" />
              Overview
            </button>
            <button
              onClick={() => setActiveTab("keys")}
              className={`px-6 py-3 font-quando font-semibold transition border-b-2 ${
                activeTab === "keys"
                  ? "border-brand-orange text-brand-orange"
                  : "border-transparent text-muted-foreground hover:text-foreground"
              }`}
            >
              <Key className="w-4 h-4 inline mr-2" />
              Access Keys
            </button>
            <button
              onClick={() => setActiveTab("users")}
              className={`px-6 py-3 font-quando font-semibold transition border-b-2 ${
                activeTab === "users"
                  ? "border-brand-orange text-brand-orange"
                  : "border-transparent text-muted-foreground hover:text-foreground"
              }`}
            >
              <Users className="w-4 h-4 inline mr-2" />
              Users
            </button>
            <button
              onClick={() => setActiveTab("monitoring")}
              className={`px-6 py-3 font-quando font-semibold transition border-b-2 ${
                activeTab === "monitoring"
                  ? "border-brand-orange text-brand-orange"
                  : "border-transparent text-muted-foreground hover:text-foreground"
              }`}
            >
              <Activity className="w-4 h-4 inline mr-2" />
              Live Monitoring
            </button>
            <button
              onClick={() => setActiveTab("analytics")}
              className={`px-6 py-3 font-quando font-semibold transition border-b-2 ${
                activeTab === "analytics"
                  ? "border-brand-orange text-brand-orange"
                  : "border-transparent text-muted-foreground hover:text-foreground"
              }`}
            >
              <TrendingUp className="w-4 h-4 inline mr-2" />
              Analytics
            </button>
          </div>
        </div>
      </div>

      {/* Content */}
      <div className="max-w-7xl mx-auto px-6 py-8">
        {loading ? (
          <div className="flex items-center justify-center py-12">
            <RefreshCw className="w-8 h-8 text-brand-orange animate-spin" />
          </div>
        ) : (
          <>
            {activeTab === "overview" && renderOverview()}
            {activeTab === "keys" && renderAccessKeys()}
            {activeTab === "users" && renderUsers()}
            {activeTab === "monitoring" && renderMonitoring()}
            {activeTab === "analytics" && renderAnalytics()}
          </>
        )}
      </div>
    </div>
  );
}

// ========== MODALS ==========

function CreateAccessKeyModal({
  onClose,
  onSuccess,
}: {
  onClose: () => void;
  onSuccess: () => void;
}) {
  const toast = useToast();
  const [formData, setFormData] = useState({
    description: "",
    max_uses: 1,
    expires_in_days: 30,
  });
  const [createdKey, setCreatedKey] = useState<string | null>(null);
  const [loading, setLoading] = useState(false);

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    setLoading(true);

    try {
      const result = await api.createAccessKey({
        description: formData.description || undefined,
        max_uses: formData.max_uses,
        expires_in_days: formData.expires_in_days || undefined,
      });
      setCreatedKey(result.key);
      toast.success("Access key created successfully");
    } catch (error: any) {
      toast.error(error.message || "Failed to create access key");
    } finally {
      setLoading(false);
    }
  };

  if (createdKey) {
    return (
      <div className="fixed inset-0 bg-black/50 flex items-center justify-center z-50 p-4">
        <div className="bg-card border border-border rounded-lg max-w-md w-full p-6">
          <h3 className="text-xl font-quando font-bold mb-4">Access Key Created</h3>
          <p className="text-sm text-muted-foreground mb-4">
            Save this key now. You won't be able to see it again!
          </p>
          <div className="bg-accent p-4 rounded-lg mb-6">
            <code className="text-sm break-all">{createdKey}</code>
          </div>
          <div className="flex gap-3">
            <button
              onClick={() => {
                navigator.clipboard.writeText(createdKey);
                toast.success("Key copied to clipboard");
              }}
              className="flex-1 px-4 py-2 border border-border rounded-md hover:bg-accent transition font-quando"
            >
              Copy Key
            </button>
            <button
              onClick={() => {
                onSuccess();
                onClose();
              }}
              className="flex-1 px-4 py-2 bg-brand-orange text-white rounded-md hover:bg-brand-orange/90 transition font-quando"
            >
              Done
            </button>
          </div>
        </div>
      </div>
    );
  }

  return (
    <div className="fixed inset-0 bg-black/50 flex items-center justify-center z-50 p-4">
      <div className="bg-card border border-border rounded-lg max-w-md w-full p-6">
        <h3 className="text-xl font-quando font-bold mb-4">Create Access Key</h3>
        <form onSubmit={handleSubmit} className="space-y-4">
          <div>
            <label className="block text-sm font-quando font-semibold mb-2">
              Description (optional)
            </label>
            <input
              type="text"
              value={formData.description}
              onChange={(e) =>
                setFormData({ ...formData, description: e.target.value })
              }
              className="w-full px-3 py-2 border border-border rounded-md bg-background"
              placeholder="e.g., Beta testers batch 1"
            />
          </div>

          <div>
            <label className="block text-sm font-quando font-semibold mb-2">
              Max Uses *
            </label>
            <input
              type="number"
              min="1"
              max="1000"
              value={formData.max_uses}
              onChange={(e) =>
                setFormData({ ...formData, max_uses: Number(e.target.value) })
              }
              className="w-full px-3 py-2 border border-border rounded-md bg-background"
              required
            />
            <p className="text-xs text-muted-foreground mt-1">
              Number of users who can register with this key
            </p>
          </div>

          <div>
            <label className="block text-sm font-quando font-semibold mb-2">
              Expires in (days)
            </label>
            <input
              type="number"
              min="1"
              max="365"
              value={formData.expires_in_days}
              onChange={(e) =>
                setFormData({ ...formData, expires_in_days: Number(e.target.value) })
              }
              className="w-full px-3 py-2 border border-border rounded-md bg-background"
            />
            <p className="text-xs text-muted-foreground mt-1">
              Leave empty for no expiration
            </p>
          </div>

          <div className="flex gap-3 pt-4">
            <button
              type="button"
              onClick={onClose}
              className="flex-1 px-4 py-2 border border-border rounded-md hover:bg-accent transition font-quando"
              disabled={loading}
            >
              Cancel
            </button>
            <button
              type="submit"
              className="flex-1 px-4 py-2 bg-brand-orange text-white rounded-md hover:bg-brand-orange/90 transition font-quando disabled:opacity-50"
              disabled={loading}
            >
              {loading ? "Creating..." : "Create Key"}
            </button>
          </div>
        </form>
      </div>
    </div>
  );
}

function UserDetailsModal({
  user,
  onClose,
  onUpdate,
}: {
  user: User;
  onClose: () => void;
  onUpdate: () => void;
}) {
  const toast = useToast();
  const [editing, setEditing] = useState(false);
  const [formData, setFormData] = useState({
    role: user.role,
    is_active: user.is_active,
    is_verified: user.is_verified,
  });

  const handleUpdate = async () => {
    try {
      await api.updateUser(user.id, formData);
      toast.success("User updated successfully");
      onUpdate();
    } catch (error: any) {
      toast.error(error.message || "Failed to update user");
    }
  };

  const handleDelete = async () => {
    if (
      !confirm(
        `Permanently delete user ${user.email}? This will also delete all their experiments and campaigns.`
      )
    ) {
      return;
    }

    try {
      await api.deleteUser(user.id);
      toast.success("User deleted successfully");
      onUpdate();
    } catch (error: any) {
      toast.error(error.message || "Failed to delete user");
    }
  };

  return (
    <div className="fixed inset-0 bg-black/50 flex items-center justify-center z-50 p-4">
      <div className="bg-card border border-border rounded-lg max-w-2xl w-full p-6 max-h-[90vh] overflow-y-auto">
        <div className="flex items-center justify-between mb-6">
          <h3 className="text-xl font-quando font-bold">User Details</h3>
          <button
            onClick={onClose}
            className="p-2 hover:bg-accent rounded transition"
          >
            <XCircle className="w-5 h-5" />
          </button>
        </div>

        <div className="space-y-6">
          {/* Basic Info */}
          <div className="grid grid-cols-2 gap-4">
            <div>
              <label className="text-sm text-muted-foreground">Email</label>
              <div className="font-semibold">{user.email}</div>
            </div>
            <div>
              <label className="text-sm text-muted-foreground">Full Name</label>
              <div className="font-semibold">{user.full_name}</div>
            </div>
            <div>
              <label className="text-sm text-muted-foreground">User ID</label>
              <div className="font-mono text-sm">{user.id}</div>
            </div>
            <div>
              <label className="text-sm text-muted-foreground">Auth Type</label>
              <div>
                {user.has_google_auth ? (
                  <span className="text-blue-500">Google OAuth</span>
                ) : (
                  <span>Email/Password</span>
                )}
              </div>
            </div>
          </div>

          {/* Editable Fields */}
          <div className="border-t border-border pt-6">
            <div className="flex items-center justify-between mb-4">
              <h4 className="font-quando font-bold">Account Settings</h4>
              <button
                onClick={() => setEditing(!editing)}
                className="text-sm text-brand-orange hover:underline"
              >
                {editing ? "Cancel" : "Edit"}
              </button>
            </div>

            <div className="grid grid-cols-2 gap-4">
              <div>
                <label className="text-sm text-muted-foreground block mb-2">Role</label>
                {editing ? (
                  <select
                    value={formData.role}
                    onChange={(e) =>
                      setFormData({ ...formData, role: e.target.value })
                    }
                    className="w-full px-3 py-2 border border-border rounded-md bg-background"
                  >
                    <option value="admin">Admin</option>
                    <option value="user">User</option>
                    <option value="viewer">Viewer</option>
                  </select>
                ) : (
                  <div className="font-semibold capitalize">{user.role}</div>
                )}
              </div>

              <div className="space-y-2">
                <label className="flex items-center gap-2 cursor-pointer">
                  <input
                    type="checkbox"
                    checked={editing ? formData.is_active : user.is_active}
                    onChange={(e) =>
                      setFormData({ ...formData, is_active: e.target.checked })
                    }
                    disabled={!editing}
                    className="rounded border-border"
                  />
                  <span className="text-sm">Active</span>
                </label>
                <label className="flex items-center gap-2 cursor-pointer">
                  <input
                    type="checkbox"
                    checked={editing ? formData.is_verified : user.is_verified}
                    onChange={(e) =>
                      setFormData({ ...formData, is_verified: e.target.checked })
                    }
                    disabled={!editing}
                    className="rounded border-border"
                  />
                  <span className="text-sm">Verified</span>
                </label>
              </div>
            </div>

            {editing && (
              <button
                onClick={handleUpdate}
                className="mt-4 px-4 py-2 bg-brand-orange text-white rounded-md hover:bg-brand-orange/90 transition font-quando"
              >
                Save Changes
              </button>
            )}
          </div>

          {/* Activity Stats */}
          <div className="border-t border-border pt-6">
            <h4 className="font-quando font-bold mb-4">Activity Statistics</h4>
            <div className="grid grid-cols-2 gap-4">
              <div>
                <label className="text-sm text-muted-foreground">Experiments</label>
                <div className="text-2xl font-bold text-brand-orange">
                  {user.experiments_count}
                </div>
              </div>
              <div>
                <label className="text-sm text-muted-foreground">Campaigns</label>
                <div className="text-2xl font-bold text-brand-orange">
                  {user.campaigns_count}
                </div>
              </div>
              <div>
                <label className="text-sm text-muted-foreground">Joined</label>
                <div className="text-sm">
                  {new Date(user.created_at).toLocaleDateString()}
                </div>
              </div>
              <div>
                <label className="text-sm text-muted-foreground">Last Login</label>
                <div className="text-sm">
                  {user.last_login
                    ? new Date(user.last_login).toLocaleDateString()
                    : "Never"}
                </div>
              </div>
            </div>
          </div>

          {/* Danger Zone */}
          <div className="border-t border-red-500/20 pt-6">
            <h4 className="font-quando font-bold text-red-500 mb-4">Danger Zone</h4>
            <button
              onClick={handleDelete}
              className="px-4 py-2 bg-red-500 text-white rounded-md hover:bg-red-600 transition font-quando flex items-center gap-2"
            >
              <Trash2 className="w-4 h-4" />
              Delete User
            </button>
            <p className="text-xs text-muted-foreground mt-2">
              This action cannot be undone. All user data will be permanently deleted.
            </p>
          </div>
        </div>
      </div>
    </div>
  );
}
