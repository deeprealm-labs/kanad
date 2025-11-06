"use client";

import { useState, useEffect } from "react";
import {
  ArrowLeft,
  Plus,
  Play,
  Trash2,
  GripVertical,
  Loader2,
  CheckCircle,
} from "lucide-react";
import Link from "next/link";
import * as api from "@/lib/api";
import { useRouter } from "next/navigation";

interface QueuedExperiment {
  id: string;
  name: string;
  molecule: any;
  method: string;
  backend: string;
  status: string;
  createdAt: string;
}

export default function QueuePage() {
  const [queue, setQueue] = useState<QueuedExperiment[]>([]);
  const [loading, setLoading] = useState(true);
  const [campaignName, setCampaignName] = useState("");
  const [showCampaignDialog, setShowCampaignDialog] = useState(false);
  const [isCreatingCampaign, setIsCreatingCampaign] = useState(false);
  const [isExecuting, setIsExecuting] = useState(false);
  const [draggedIndex, setDraggedIndex] = useState<number | null>(null);
  const router = useRouter();

  useEffect(() => {
    loadQueue();
  }, []);

  const loadQueue = async () => {
    try {
      setLoading(true);
      const response = await api.getQueue();
      setQueue(response.queue || []);
    } catch (error) {
      console.error("Failed to load queue:", error);
    } finally {
      setLoading(false);
    }
  };

  const handleCreateAndExecuteCampaign = async () => {
    if (queue.length === 0) {
      alert("No experiments in queue!");
      return;
    }

    if (!campaignName.trim()) {
      alert("Please enter a campaign name");
      return;
    }

    try {
      setIsCreatingCampaign(true);

      // Create campaign
      const campaignResponse = await api.createCampaign({
        name: campaignName,
        description: `Sequential execution of ${queue.length} experiments`,
        experiment_ids: queue.map((exp) => exp.id),
      });

      const campaignId = campaignResponse.campaign_id;

      // Execute campaign
      await api.executeCampaign(campaignId);

      // Navigate to campaign monitor
      router.push(`/dashboard/campaign/${campaignId}`);
    } catch (error: any) {
      console.error("Failed to execute campaign:", error);
      alert(`Failed to execute campaign: ${error.message || "Unknown error"}`);
    } finally {
      setIsCreatingCampaign(false);
      setShowCampaignDialog(false);
    }
  };

  const handleQuickExecute = async () => {
    if (queue.length === 0) {
      alert("No experiments in queue!");
      return;
    }

    const name = `Campaign ${new Date().toLocaleString()}`;

    try {
      setIsExecuting(true);

      // Create and execute campaign
      const campaignResponse = await api.createCampaign({
        name,
        description: `Sequential execution of ${queue.length} experiments`,
        experiment_ids: queue.map((exp) => exp.id),
      });

      const campaignId = campaignResponse.campaign_id;
      await api.executeCampaign(campaignId);

      // Navigate to campaign monitor
      router.push(`/dashboard/campaign/${campaignId}`);
    } catch (error: any) {
      console.error("Failed to execute campaign:", error);
      alert(`Failed to execute campaign: ${error.message || "Unknown error"}`);
    } finally {
      setIsExecuting(false);
    }
  };

  const handleRemoveFromQueue = async (experimentId: string) => {
    try {
      await api.deleteExperiment(experimentId);
      setQueue(queue.filter((exp) => exp.id !== experimentId));
    } catch (error) {
      console.error("Failed to remove experiment:", error);
    }
  };

  const handleDragStart = (index: number) => {
    setDraggedIndex(index);
  };

  const handleDragOver = (e: React.DragEvent, index: number) => {
    e.preventDefault();
    if (draggedIndex === null || draggedIndex === index) return;

    const newQueue = [...queue];
    const draggedItem = newQueue[draggedIndex];
    newQueue.splice(draggedIndex, 1);
    newQueue.splice(index, 0, draggedItem);

    setQueue(newQueue);
    setDraggedIndex(index);
  };

  const handleDragEnd = () => {
    setDraggedIndex(null);
  };

  if (loading) {
    return (
      <div className="min-h-screen bg-background flex items-center justify-center">
        <Loader2 className="w-8 h-8 animate-spin" />
      </div>
    );
  }

  return (
    <div className="h-full flex flex-col bg-background overflow-hidden">
      {/* Header */}
      <div className="px-8 pt-4 pb-3 border-b border-border flex-shrink-0">
        <div className="flex items-center justify-between">
          <div>
            <h1 className="text-2xl font-quando font-bold tracking-tight">Job Queue</h1>
            <p className="text-muted-foreground font-quando text-xs mt-1">
              Manage and execute experiments sequentially
            </p>
          </div>

          <Link href="/dashboard">
            <button className="px-6 py-3 bg-brand-orange text-white rounded-lg hover:bg-brand-orange-dark transition font-quando flex items-center gap-2">
              <Plus className="w-5 h-5" />
              New Experiment
            </button>
          </Link>
        </div>
      </div>

      {/* Main Content */}
      <div className="flex-1 flex flex-col p-6 overflow-hidden">
        <div className="flex-1 overflow-y-auto space-y-6">
          {/* Stats */}
          <div className="grid grid-cols-4 gap-4 mb-6">
            <div className="bg-card border border-border rounded-lg p-6">
              <p className="text-sm text-muted-foreground font-quando mb-1">
                Total in Queue
              </p>
              <p className="text-3xl font-quando font-bold">{queue.length}</p>
            </div>
            <div className="bg-card border border-border rounded-lg p-6">
              <p className="text-sm text-muted-foreground font-quando mb-1">
                Status
              </p>
              <p className="text-lg font-quando font-bold text-blue-600">
                {queue.length > 0 ? "Ready" : "Empty"}
              </p>
            </div>
            <div className="bg-card border border-border rounded-lg p-6">
              <p className="text-sm text-muted-foreground font-quando mb-1">
                Estimated Time
              </p>
              <p className="text-lg font-quando font-bold">
                ~{queue.length * 2} min
              </p>
            </div>
            <div className="bg-card border border-border rounded-lg p-6">
              <p className="text-sm text-muted-foreground font-quando mb-1">
                Mode
              </p>
              <p className="text-lg font-quando font-bold">Sequential</p>
            </div>
          </div>

          {/* Execute Buttons */}
          {queue.length > 0 && (
            <div className="bg-card border border-border rounded-lg p-6 mb-8">
              <div className="flex items-center justify-between">
                <div>
                  <h3 className="text-lg font-quando font-semibold mb-2">
                    Execute Campaign
                  </h3>
                  <p className="text-sm text-muted-foreground">
                    Run all {queue.length} experiments sequentially as one batch
                  </p>
                </div>
                <div className="flex gap-3">
                  <button
                    onClick={() => setShowCampaignDialog(true)}
                    className="px-6 py-3 border border-border rounded-lg hover:bg-accent transition font-quando"
                    disabled={isExecuting || isCreatingCampaign}
                  >
                    Name & Execute
                  </button>
                  <button
                    onClick={handleQuickExecute}
                    className="px-6 py-3 bg-brand-orange text-white rounded-lg hover:bg-brand-orange-dark transition font-quando flex items-center gap-2"
                    disabled={isExecuting || isCreatingCampaign}
                  >
                    {isExecuting ? (
                      <>
                        <Loader2 className="w-5 h-5 animate-spin" />
                        Starting...
                      </>
                    ) : (
                      <>
                        <Play className="w-5 h-5" />
                        Quick Execute
                      </>
                    )}
                  </button>
                </div>
              </div>
            </div>
          )}

          {/* Queue List */}
          {queue.length === 0 ? (
            <div className="bg-card border border-border rounded-lg p-12 flex flex-col items-center justify-center text-center">
              <CheckCircle className="w-16 h-16 text-muted-foreground mb-4" />
              <h3 className="text-xl font-quando font-semibold mb-2">
                No experiments in queue
              </h3>
              <p className="text-muted-foreground mb-6">
                Start creating experiments to build your queue
              </p>
              <Link href="/dashboard">
                <button className="px-6 py-3 bg-brand-orange text-white rounded-lg hover:bg-brand-orange-dark transition font-quando">
                  Create Experiment
                </button>
              </Link>
            </div>
          ) : (
            <div className="space-y-3">
              <div className="flex items-center justify-between mb-4">
                <h2 className="text-xl font-quando font-semibold">
                  Queued Experiments ({queue.length})
                </h2>
                <p className="text-sm text-muted-foreground">
                  Drag to reorder execution sequence
                </p>
              </div>

              {queue.map((exp, index) => (
                <div
                  key={exp.id}
                  draggable
                  onDragStart={() => handleDragStart(index)}
                  onDragOver={(e) => handleDragOver(e, index)}
                  onDragEnd={handleDragEnd}
                  className={`bg-card border border-border rounded-lg p-4 flex items-center gap-4 hover:border-brand-orange transition cursor-move ${
                    draggedIndex === index ? "opacity-50" : ""
                  }`}
                >
                  <div className="flex items-center gap-3">
                    <GripVertical className="w-5 h-5 text-muted-foreground" />
                    <div className="w-8 h-8 rounded-full bg-brand-orange/10 text-brand-orange flex items-center justify-center font-semibold">
                      {index + 1}
                    </div>
                  </div>

                  <div className="flex-1">
                    <h3 className="font-quando font-semibold">{exp.name || "Unnamed Experiment"}</h3>
                    <div className="flex gap-4 text-sm text-muted-foreground mt-1">
                      <span>Method: {exp.method}</span>
                      <span>Backend: {exp.backend}</span>
                    </div>
                  </div>

                  <button
                    onClick={() => handleRemoveFromQueue(exp.id)}
                    className="p-2 hover:bg-red-100 dark:hover:bg-red-900/20 rounded-lg transition text-red-600"
                  >
                    <Trash2 className="w-5 h-5" />
                  </button>
                </div>
              ))}
            </div>
          )}
        </div>
      </div>

      {/* Campaign Name Dialog */}
      {showCampaignDialog && (
        <div className="fixed inset-0 z-50 flex items-center justify-center backdrop-blur-sm bg-black/20">
          <div className="bg-background border border-border rounded-lg w-full max-w-md p-6">
            <h3 className="text-xl font-quando font-semibold mb-4">
              Name Your Campaign
            </h3>
            <p className="text-sm text-muted-foreground mb-4">
              Give this batch execution a memorable name
            </p>
            <input
              type="text"
              value={campaignName}
              onChange={(e) => setCampaignName(e.target.value)}
              placeholder="e.g., H2 Bond Length Scan"
              className="w-full px-4 py-3 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando mb-6"
              onKeyPress={(e) => {
                if (e.key === "Enter") handleCreateAndExecuteCampaign();
              }}
            />
            <div className="flex gap-3">
              <button
                onClick={() => setShowCampaignDialog(false)}
                className="flex-1 px-4 py-2 border border-border rounded-lg hover:bg-accent transition font-quando"
                disabled={isCreatingCampaign}
              >
                Cancel
              </button>
              <button
                onClick={handleCreateAndExecuteCampaign}
                className="flex-1 px-4 py-2 bg-brand-orange text-white rounded-lg hover:bg-brand-orange-dark transition font-quando flex items-center justify-center gap-2"
                disabled={isCreatingCampaign}
              >
                {isCreatingCampaign ? (
                  <>
                    <Loader2 className="w-4 h-4 animate-spin" />
                    Creating...
                  </>
                ) : (
                  <>
                    <Play className="w-4 h-4" />
                    Execute
                  </>
                )}
              </button>
            </div>
          </div>
        </div>
      )}
    </div>
  );
}
