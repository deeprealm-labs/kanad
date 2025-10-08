"use client";

import { useState, useEffect } from "react";
import { Save, Check, X, Key, Cloud, Trash2 } from "lucide-react";
import * as api from "@/lib/api";

export default function SettingsPage() {
  // IBM Quantum credentials
  const [ibmCrn, setIbmCrn] = useState("");
  const [ibmApiKey, setIbmApiKey] = useState("");
  const [ibmConfigured, setIbmConfigured] = useState(false);
  const [ibmSaving, setIbmSaving] = useState(false);
  const [ibmMessage, setIbmMessage] = useState<{ type: "success" | "error"; text: string } | null>(null);

  // BlueQubit credentials
  const [blueQubitToken, setBlueQubitToken] = useState("");
  const [blueQubitConfigured, setBlueQubitConfigured] = useState(false);
  const [blueQubitSaving, setBlueQubitSaving] = useState(false);
  const [blueQubitMessage, setBlueQubitMessage] = useState<{ type: "success" | "error"; text: string } | null>(null);

  useEffect(() => {
    loadCredentialsStatus();
  }, []);

  const loadCredentialsStatus = async () => {
    try {
      const status = await api.getCloudCredentialsStatus();
      setIbmConfigured(status.ibm.configured);
      setBlueQubitConfigured(status.bluequbit.configured);
    } catch (error) {
      console.error("Failed to load credentials status:", error);
    }
  };

  const handleSaveIBM = async () => {
    if (!ibmCrn || !ibmApiKey) {
      setIbmMessage({ type: "error", text: "Please enter both CRN and API Key" });
      return;
    }

    setIbmSaving(true);
    setIbmMessage(null);

    try {
      const response = await api.saveIBMCredentials(ibmCrn, ibmApiKey);
      setIbmMessage({ type: "success", text: response.message });
      setIbmConfigured(true);
      // Clear form for security
      setIbmCrn("");
      setIbmApiKey("");

      // Clear success message after 3 seconds
      setTimeout(() => setIbmMessage(null), 3000);
    } catch (error: any) {
      setIbmMessage({ type: "error", text: error.message || "Failed to save credentials" });
    } finally {
      setIbmSaving(false);
    }
  };

  const handleDeleteIBM = async () => {
    if (!confirm("Are you sure you want to delete IBM Quantum credentials?")) {
      return;
    }

    try {
      await api.deleteIBMCredentials();
      setIbmConfigured(false);
      setIbmMessage({ type: "success", text: "IBM Quantum credentials deleted" });
      setTimeout(() => setIbmMessage(null), 3000);
    } catch (error: any) {
      setIbmMessage({ type: "error", text: error.message || "Failed to delete credentials" });
    }
  };

  const handleSaveBlueQubit = async () => {
    if (!blueQubitToken) {
      setBlueQubitMessage({ type: "error", text: "Please enter API Token" });
      return;
    }

    setBlueQubitSaving(true);
    setBlueQubitMessage(null);

    try {
      const response = await api.saveBlueQubitCredentials(blueQubitToken);
      setBlueQubitMessage({ type: "success", text: response.message });
      setBlueQubitConfigured(true);
      // Clear form for security
      setBlueQubitToken("");

      // Clear success message after 3 seconds
      setTimeout(() => setBlueQubitMessage(null), 3000);
    } catch (error: any) {
      setBlueQubitMessage({ type: "error", text: error.message || "Failed to save credentials" });
    } finally {
      setBlueQubitSaving(false);
    }
  };

  const handleDeleteBlueQubit = async () => {
    if (!confirm("Are you sure you want to delete BlueQubit credentials?")) {
      return;
    }

    try {
      await api.deleteBlueQubitCredentials();
      setBlueQubitConfigured(false);
      setBlueQubitMessage({ type: "success", text: "BlueQubit credentials deleted" });
      setTimeout(() => setBlueQubitMessage(null), 3000);
    } catch (error: any) {
      setBlueQubitMessage({ type: "error", text: error.message || "Failed to delete credentials" });
    }
  };

  return (
    <div className="h-full overflow-auto p-6 bg-background">
      <div className="max-w-4xl mx-auto">
        {/* Header */}
        <div className="mb-8">
          <h1 className="text-3xl font-quando font-bold mb-2">Settings</h1>
          <p className="text-muted-foreground font-quando">
            Configure cloud credentials for quantum backends
          </p>
        </div>

        {/* IBM Quantum Section */}
        <div className="bg-card border border-border rounded-lg p-6 mb-6">
          <div className="flex items-center justify-between mb-4">
            <div className="flex items-center gap-3">
              <div className="w-10 h-10 bg-blue-500/10 rounded-lg flex items-center justify-center">
                <Cloud className="w-6 h-6 text-blue-500" />
              </div>
              <div>
                <h2 className="text-xl font-quando font-semibold">IBM Quantum</h2>
                <p className="text-sm text-muted-foreground font-quando">
                  Configure access to IBM Quantum hardware and simulators
                </p>
              </div>
            </div>
            {ibmConfigured && (
              <div className="flex items-center gap-2">
                <div className="flex items-center gap-2 px-3 py-1 bg-green-500/10 rounded-full">
                  <Check className="w-4 h-4 text-green-500" />
                  <span className="text-sm font-quando text-green-500">Configured</span>
                </div>
                <button
                  onClick={handleDeleteIBM}
                  className="p-2 hover:bg-red-500/10 rounded-lg transition"
                  title="Delete credentials"
                >
                  <Trash2 className="w-4 h-4 text-red-500" />
                </button>
              </div>
            )}
          </div>

          <div className="space-y-4">
            <div>
              <label className="block text-sm font-quando font-medium mb-2">
                Cloud Resource Name (CRN)
              </label>
              <input
                type="text"
                value={ibmCrn}
                onChange={(e) => setIbmCrn(e.target.value)}
                placeholder="crn:v1:bluemix:public:quantum-computing:..."
                className="w-full px-4 py-2 border border-border rounded-lg bg-background focus:outline-none focus:ring-2 focus:ring-brand-orange font-mono text-sm"
              />
              <p className="text-xs text-muted-foreground mt-1 font-quando">
                Found in IBM Quantum dashboard under your service instance
              </p>
            </div>

            <div>
              <label className="block text-sm font-quando font-medium mb-2 flex items-center gap-2">
                <Key className="w-4 h-4" />
                API Key
              </label>
              <input
                type="password"
                value={ibmApiKey}
                onChange={(e) => setIbmApiKey(e.target.value)}
                placeholder="Enter your IBM Quantum API key"
                className="w-full px-4 py-2 border border-border rounded-lg bg-background focus:outline-none focus:ring-2 focus:ring-brand-orange font-mono text-sm"
              />
              <p className="text-xs text-muted-foreground mt-1 font-quando">
                Your API key is stored securely and never displayed
              </p>
            </div>

            {ibmMessage && (
              <div
                className={`p-3 rounded-lg flex items-center gap-2 ${
                  ibmMessage.type === "success"
                    ? "bg-green-500/10 text-green-500"
                    : "bg-red-500/10 text-red-500"
                }`}
              >
                {ibmMessage.type === "success" ? (
                  <Check className="w-4 h-4" />
                ) : (
                  <X className="w-4 h-4" />
                )}
                <span className="text-sm font-quando">{ibmMessage.text}</span>
              </div>
            )}

            <button
              onClick={handleSaveIBM}
              disabled={ibmSaving}
              className="flex items-center gap-2 px-6 py-2 bg-brand-orange text-white rounded-lg hover:bg-brand-orange-dark transition font-quando disabled:opacity-50 disabled:cursor-not-allowed"
            >
              <Save className="w-4 h-4" />
              {ibmSaving ? "Saving..." : "Save IBM Credentials"}
            </button>
          </div>
        </div>

        {/* BlueQubit Section */}
        <div className="bg-card border border-border rounded-lg p-6">
          <div className="flex items-center justify-between mb-4">
            <div className="flex items-center gap-3">
              <div className="w-10 h-10 bg-purple-500/10 rounded-lg flex items-center justify-center">
                <Cloud className="w-6 h-6 text-purple-500" />
              </div>
              <div>
                <h2 className="text-xl font-quando font-semibold">BlueQubit</h2>
                <p className="text-sm text-muted-foreground font-quando">
                  Configure access to BlueQubit GPU-accelerated simulators
                </p>
              </div>
            </div>
            {blueQubitConfigured && (
              <div className="flex items-center gap-2">
                <div className="flex items-center gap-2 px-3 py-1 bg-green-500/10 rounded-full">
                  <Check className="w-4 h-4 text-green-500" />
                  <span className="text-sm font-quando text-green-500">Configured</span>
                </div>
                <button
                  onClick={handleDeleteBlueQubit}
                  className="p-2 hover:bg-red-500/10 rounded-lg transition"
                  title="Delete credentials"
                >
                  <Trash2 className="w-4 h-4 text-red-500" />
                </button>
              </div>
            )}
          </div>

          <div className="space-y-4">
            <div>
              <label className="block text-sm font-quando font-medium mb-2 flex items-center gap-2">
                <Key className="w-4 h-4" />
                API Token
              </label>
              <input
                type="password"
                value={blueQubitToken}
                onChange={(e) => setBlueQubitToken(e.target.value)}
                placeholder="Enter your BlueQubit API token"
                className="w-full px-4 py-2 border border-border rounded-lg bg-background focus:outline-none focus:ring-2 focus:ring-brand-orange font-mono text-sm"
              />
              <p className="text-xs text-muted-foreground mt-1 font-quando">
                Get your API token from the BlueQubit dashboard
              </p>
            </div>

            {blueQubitMessage && (
              <div
                className={`p-3 rounded-lg flex items-center gap-2 ${
                  blueQubitMessage.type === "success"
                    ? "bg-green-500/10 text-green-500"
                    : "bg-red-500/10 text-red-500"
                }`}
              >
                {blueQubitMessage.type === "success" ? (
                  <Check className="w-4 h-4" />
                ) : (
                  <X className="w-4 h-4" />
                )}
                <span className="text-sm font-quando">{blueQubitMessage.text}</span>
              </div>
            )}

            <button
              onClick={handleSaveBlueQubit}
              disabled={blueQubitSaving}
              className="flex items-center gap-2 px-6 py-2 bg-brand-orange text-white rounded-lg hover:bg-brand-orange-dark transition font-quando disabled:opacity-50 disabled:cursor-not-allowed"
            >
              <Save className="w-4 h-4" />
              {blueQubitSaving ? "Saving..." : "Save BlueQubit Credentials"}
            </button>
          </div>
        </div>

        {/* Information Footer */}
        <div className="mt-6 p-4 bg-muted rounded-lg">
          <p className="text-sm text-muted-foreground font-quando">
            <strong>Note:</strong> Credentials are stored securely in the database. They are required to submit jobs to cloud quantum backends. You can update or delete them at any time.
          </p>
        </div>
      </div>
    </div>
  );
}
