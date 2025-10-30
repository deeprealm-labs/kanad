"use client";

import { useState, useEffect } from "react";
import { X, User, Lock, BarChart3, LogOut, Trash2, Save } from "lucide-react";
import * as api from "@/lib/api";
import { useToast } from "@/components/ui/toast";
import { useAuth } from "@/contexts/AuthContext";

interface UserSettingsModalProps {
  isOpen: boolean;
  onClose: () => void;
}

type Tab = "profile" | "account" | "usage";

export default function UserSettingsModal({ isOpen, onClose }: UserSettingsModalProps) {
  const [activeTab, setActiveTab] = useState<Tab>("profile");
  const toast = useToast();
  const { user, logout } = useAuth();

  // Profile state
  const [fullName, setFullName] = useState("");
  const [email, setEmail] = useState("");
  const [isGoogleUser, setIsGoogleUser] = useState(false);
  const [profileLoading, setProfileLoading] = useState(false);

  // Password state
  const [currentPassword, setCurrentPassword] = useState("");
  const [newPassword, setNewPassword] = useState("");
  const [confirmPassword, setConfirmPassword] = useState("");
  const [passwordLoading, setPasswordLoading] = useState(false);

  // Statistics state
  const [statistics, setStatistics] = useState<any>(null);
  const [statsLoading, setStatsLoading] = useState(false);

  // Sessions state
  const [sessions, setSessions] = useState<any[]>([]);
  const [sessionsLoading, setSessionsLoading] = useState(false);

  // Load profile data
  useEffect(() => {
    if (isOpen && activeTab === "profile") {
      loadProfile();
    }
  }, [isOpen, activeTab]);

  // Load statistics
  useEffect(() => {
    if (isOpen && activeTab === "usage") {
      loadStatistics();
    }
  }, [isOpen, activeTab]);

  // Load sessions
  useEffect(() => {
    if (isOpen && activeTab === "account") {
      loadSessions();
    }
  }, [isOpen, activeTab]);

  const loadProfile = async () => {
    try {
      const profile = await api.getUserProfile();
      setFullName(profile.full_name || "");
      setEmail(profile.email || "");
      setIsGoogleUser(!!profile.google_id);
    } catch (error: any) {
      toast.error(error.message || "Failed to load profile");
    }
  };

  const loadStatistics = async () => {
    setStatsLoading(true);
    try {
      const stats = await api.getUserStatistics();
      setStatistics(stats);
    } catch (error: any) {
      toast.error(error.message || "Failed to load statistics");
    } finally {
      setStatsLoading(false);
    }
  };

  const loadSessions = async () => {
    setSessionsLoading(true);
    try {
      const sessionsList = await api.getUserSessions();
      setSessions(sessionsList);
    } catch (error: any) {
      toast.error(error.message || "Failed to load sessions");
    } finally {
      setSessionsLoading(false);
    }
  };

  const handleUpdateProfile = async () => {
    setProfileLoading(true);
    try {
      await api.updateUserProfile({
        full_name: fullName,
      });
      toast.success("Profile updated successfully");
    } catch (error: any) {
      toast.error(error.message || "Failed to update profile");
    } finally {
      setProfileLoading(false);
    }
  };

  const handleChangePassword = async () => {
    if (newPassword !== confirmPassword) {
      toast.error("Passwords do not match");
      return;
    }

    if (newPassword.length < 8) {
      toast.error("Password must be at least 8 characters");
      return;
    }

    setPasswordLoading(true);
    try {
      await api.changePassword({
        current_password: currentPassword,
        new_password: newPassword,
      });
      toast.success("Password changed successfully");
      setCurrentPassword("");
      setNewPassword("");
      setConfirmPassword("");
    } catch (error: any) {
      toast.error(error.message || "Failed to change password");
    } finally {
      setPasswordLoading(false);
    }
  };

  const handleRevokeSession = async (sessionId: number) => {
    try {
      await api.revokeSession(sessionId);
      toast.success("Session revoked successfully");
      loadSessions();
    } catch (error: any) {
      toast.error(error.message || "Failed to revoke session");
    }
  };

  const handleDeleteAccount = async () => {
    if (!confirm("Are you sure you want to deactivate your account? This action cannot be undone.")) {
      return;
    }

    try {
      await api.deleteAccount();
      toast.success("Account deactivated successfully");

      // Account is now deactivated and all sessions are revoked
      // Clear local storage and redirect without calling logout API
      localStorage.removeItem("kanad_user");
      localStorage.removeItem("kanad_access_token");
      localStorage.removeItem("kanad_refresh_token");

      onClose();
      window.location.href = "/";
    } catch (error: any) {
      toast.error(error.message || "Failed to deactivate account");
    }
  };

  if (!isOpen) return null;

  const tabs = [
    { id: "profile" as Tab, label: "Profile", icon: User },
    { id: "account" as Tab, label: "Account", icon: Lock },
    { id: "usage" as Tab, label: "Usage", icon: BarChart3 },
  ];

  return (
    <div className="fixed inset-0 bg-black/50 backdrop-blur-sm flex items-center justify-center z-50 p-4">
      <div className="bg-background border border-border rounded-lg w-full max-w-3xl max-h-[90vh] overflow-hidden flex flex-col">
        {/* Header */}
        <div className="flex items-center justify-between p-6 border-b border-border">
          <h2 className="text-2xl font-quando font-bold">User Settings</h2>
          <button
            onClick={onClose}
            className="p-2 hover:bg-accent rounded-md transition"
          >
            <X className="w-5 h-5" />
          </button>
        </div>

        {/* Tabs */}
        <div className="flex border-b border-border px-6">
          {tabs.map((tab) => {
            const Icon = tab.icon;
            return (
              <button
                key={tab.id}
                onClick={() => setActiveTab(tab.id)}
                className={`flex items-center gap-2 px-4 py-3 font-quando border-b-2 transition ${
                  activeTab === tab.id
                    ? "border-brand-orange text-brand-orange"
                    : "border-transparent text-muted-foreground hover:text-foreground"
                }`}
              >
                <Icon className="w-4 h-4" />
                {tab.label}
              </button>
            );
          })}
        </div>

        {/* Content */}
        <div className="flex-1 overflow-y-auto p-6">
          {activeTab === "profile" && (
            <div className="space-y-6">
              <div>
                <label className="block text-sm font-quando font-semibold mb-2">
                  Full Name
                </label>
                <input
                  type="text"
                  value={fullName}
                  onChange={(e) => setFullName(e.target.value)}
                  className="w-full px-4 py-2 bg-accent border border-border rounded-md font-quando focus:outline-none focus:ring-2 focus:ring-brand-orange"
                  placeholder="Enter your full name"
                />
              </div>

              <div>
                <label className="block text-sm font-quando font-semibold mb-2">
                  Email
                </label>
                <input
                  type="email"
                  value={email}
                  disabled
                  className="w-full px-4 py-2 bg-muted border border-border rounded-md font-quando opacity-60 cursor-not-allowed"
                />
                <p className="text-xs text-muted-foreground mt-1 font-quando">
                  Email cannot be changed
                </p>
              </div>

              {isGoogleUser && (
                <div className="bg-blue-500/10 border border-blue-500/20 rounded-md p-4">
                  <p className="text-sm font-quando text-blue-600 dark:text-blue-400">
                    You're signed in with Google OAuth
                  </p>
                </div>
              )}

              <button
                onClick={handleUpdateProfile}
                disabled={profileLoading}
                className="flex items-center gap-2 px-4 py-2 bg-brand-orange text-white rounded-md hover:bg-brand-orange/90 transition font-quando disabled:opacity-50"
              >
                <Save className="w-4 h-4" />
                {profileLoading ? "Saving..." : "Save Changes"}
              </button>
            </div>
          )}

          {activeTab === "account" && (
            <div className="space-y-8">
              {/* Change Password */}
              {!isGoogleUser && (
                <div className="space-y-4">
                  <h3 className="text-lg font-quando font-bold">Change Password</h3>

                  <div>
                    <label className="block text-sm font-quando font-semibold mb-2">
                      Current Password
                    </label>
                    <input
                      type="password"
                      value={currentPassword}
                      onChange={(e) => setCurrentPassword(e.target.value)}
                      className="w-full px-4 py-2 bg-accent border border-border rounded-md font-quando focus:outline-none focus:ring-2 focus:ring-brand-orange"
                    />
                  </div>

                  <div>
                    <label className="block text-sm font-quando font-semibold mb-2">
                      New Password
                    </label>
                    <input
                      type="password"
                      value={newPassword}
                      onChange={(e) => setNewPassword(e.target.value)}
                      className="w-full px-4 py-2 bg-accent border border-border rounded-md font-quando focus:outline-none focus:ring-2 focus:ring-brand-orange"
                    />
                  </div>

                  <div>
                    <label className="block text-sm font-quando font-semibold mb-2">
                      Confirm New Password
                    </label>
                    <input
                      type="password"
                      value={confirmPassword}
                      onChange={(e) => setConfirmPassword(e.target.value)}
                      className="w-full px-4 py-2 bg-accent border border-border rounded-md font-quando focus:outline-none focus:ring-2 focus:ring-brand-orange"
                    />
                  </div>

                  <button
                    onClick={handleChangePassword}
                    disabled={passwordLoading || !currentPassword || !newPassword}
                    className="flex items-center gap-2 px-4 py-2 bg-brand-orange text-white rounded-md hover:bg-brand-orange/90 transition font-quando disabled:opacity-50"
                  >
                    <Lock className="w-4 h-4" />
                    {passwordLoading ? "Changing..." : "Change Password"}
                  </button>
                </div>
              )}

              {/* Active Sessions */}
              <div className="space-y-4">
                <h3 className="text-lg font-quando font-bold">Active Sessions</h3>

                {sessionsLoading ? (
                  <p className="text-muted-foreground font-quando">Loading sessions...</p>
                ) : sessions.length === 0 ? (
                  <p className="text-muted-foreground font-quando">No active sessions</p>
                ) : (
                  <div className="space-y-2">
                    {sessions.map((session) => (
                      <div
                        key={session.id}
                        className="flex items-center justify-between p-4 bg-accent border border-border rounded-md"
                      >
                        <div>
                          <p className="font-quando font-semibold">
                            {session.is_current ? "Current Session" : "Session"}
                          </p>
                          <p className="text-xs text-muted-foreground font-quando">
                            Created: {new Date(session.created_at).toLocaleString()}
                          </p>
                        </div>
                        {!session.is_current && (
                          <button
                            onClick={() => handleRevokeSession(session.id)}
                            className="px-3 py-1 text-sm bg-red-500/10 text-red-600 dark:text-red-400 border border-red-500/20 rounded-md hover:bg-red-500/20 transition font-quando"
                          >
                            Revoke
                          </button>
                        )}
                      </div>
                    ))}
                  </div>
                )}
              </div>

              {/* Delete Account */}
              <div className="space-y-4 pt-4 border-t border-border">
                <h3 className="text-lg font-quando font-bold text-red-600 dark:text-red-400">
                  Danger Zone
                </h3>
                <p className="text-sm text-muted-foreground font-quando">
                  Deactivating your account will disable access to all experiments and data.
                </p>
                <button
                  onClick={handleDeleteAccount}
                  className="flex items-center gap-2 px-4 py-2 bg-red-500/10 text-red-600 dark:text-red-400 border border-red-500/20 rounded-md hover:bg-red-500/20 transition font-quando"
                >
                  <Trash2 className="w-4 h-4" />
                  Deactivate Account
                </button>
              </div>
            </div>
          )}

          {activeTab === "usage" && (
            <div className="space-y-6">
              {statsLoading ? (
                <p className="text-muted-foreground font-quando">Loading statistics...</p>
              ) : statistics ? (
                <div className="grid grid-cols-2 gap-4">
                  <div className="p-4 bg-accent border border-border rounded-md">
                    <p className="text-sm text-muted-foreground font-quando">
                      Total Experiments
                    </p>
                    <p className="text-3xl font-quando font-bold mt-2">
                      {statistics.total_experiments}
                    </p>
                  </div>

                  <div className="p-4 bg-accent border border-border rounded-md">
                    <p className="text-sm text-muted-foreground font-quando">
                      Completed
                    </p>
                    <p className="text-3xl font-quando font-bold mt-2 text-green-600 dark:text-green-400">
                      {statistics.completed_experiments}
                    </p>
                  </div>

                  <div className="p-4 bg-accent border border-border rounded-md">
                    <p className="text-sm text-muted-foreground font-quando">
                      Failed
                    </p>
                    <p className="text-3xl font-quando font-bold mt-2 text-red-600 dark:text-red-400">
                      {statistics.failed_experiments}
                    </p>
                  </div>

                  <div className="p-4 bg-accent border border-border rounded-md">
                    <p className="text-sm text-muted-foreground font-quando">
                      Total Campaigns
                    </p>
                    <p className="text-3xl font-quando font-bold mt-2">
                      {statistics.total_campaigns}
                    </p>
                  </div>

                  <div className="p-4 bg-accent border border-border rounded-md col-span-2">
                    <p className="text-sm text-muted-foreground font-quando">
                      Account Age
                    </p>
                    <p className="text-3xl font-quando font-bold mt-2">
                      {statistics.account_age_days} days
                    </p>
                  </div>

                  {statistics.last_experiment_date && (
                    <div className="p-4 bg-accent border border-border rounded-md col-span-2">
                      <p className="text-sm text-muted-foreground font-quando">
                        Last Experiment
                      </p>
                      <p className="text-lg font-quando font-semibold mt-2">
                        {new Date(statistics.last_experiment_date).toLocaleString()}
                      </p>
                    </div>
                  )}
                </div>
              ) : (
                <p className="text-muted-foreground font-quando">No statistics available</p>
              )}
            </div>
          )}
        </div>
      </div>
    </div>
  );
}
