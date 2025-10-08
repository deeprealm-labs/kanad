"use client";

import { useState, useEffect } from "react";
import {
  ArrowLeft,
  Plus,
  Calendar,
  Clock,
  Play,
  Pause,
  Trash2,
  MoveUp,
  MoveDown,
} from "lucide-react";
import Link from "next/link";
import * as api from "@/lib/api";

interface QueuedJob {
  id: string;
  name: string;
  molecule: any;
  method: string;
  backend: string;
  scheduledTime?: string;
  priority: number;
  status: "queued" | "scheduled" | "running" | "paused";
  createdAt: string;
}

export default function QueuePage() {
  const [queue, setQueue] = useState<QueuedJob[]>([]);
  const [showScheduleModal, setShowScheduleModal] = useState(false);
  const [selectedJob, setSelectedJob] = useState<QueuedJob | null>(null);
  const [scheduleDate, setScheduleDate] = useState("");
  const [scheduleTime, setScheduleTime] = useState("");
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    loadQueue();
  }, []);

  const loadQueue = async () => {
    try {
      setLoading(true);
      const response = await api.getQueue();
      setQueue(response.queue || []);
    } catch (error) {
      console.error("Failed to load queue from API:", error);

      // Fallback to localStorage
      const saved = localStorage.getItem("kanad_queue");
      if (saved) {
        try {
          setQueue(JSON.parse(saved));
        } catch (e) {
          console.error("Failed to parse localStorage queue:", e);
        }
      }
    } finally {
      setLoading(false);
    }
  };

  const saveQueue = (updatedQueue: QueuedJob[]) => {
    setQueue(updatedQueue);
    localStorage.setItem("kanad_queue", JSON.stringify(updatedQueue));
  };

  const handleSchedule = () => {
    if (selectedJob && scheduleDate && scheduleTime) {
      const scheduled = `${scheduleDate}T${scheduleTime}`;
      const updated = queue.map((job) =>
        job.id === selectedJob.id
          ? { ...job, scheduledTime: scheduled, status: "scheduled" as const }
          : job
      );
      saveQueue(updated);
      setShowScheduleModal(false);
      setSelectedJob(null);
      setScheduleDate("");
      setScheduleTime("");
    }
  };

  const handleDelete = (id: string) => {
    const updated = queue.filter((job) => job.id !== id);
    saveQueue(updated);
  };

  const handlePriorityChange = (id: string, direction: "up" | "down") => {
    const index = queue.findIndex((job) => job.id === id);
    if (index === -1) return;

    const updated = [...queue];
    if (direction === "up" && index > 0) {
      [updated[index], updated[index - 1]] = [updated[index - 1], updated[index]];
    } else if (direction === "down" && index < queue.length - 1) {
      [updated[index], updated[index + 1]] = [updated[index + 1], updated[index]];
    }

    // Update priority numbers
    updated.forEach((job, i) => {
      job.priority = i + 1;
    });

    saveQueue(updated);
  };

  const handleTogglePause = (id: string) => {
    const updated = queue.map((job) =>
      job.id === id
        ? {
            ...job,
            status:
              job.status === "paused"
                ? ("queued" as const)
                : ("paused" as const),
          }
        : job
    );
    saveQueue(updated);
  };

  const getStatusColor = (status: string) => {
    switch (status) {
      case "running":
        return "bg-blue-100 text-blue-800 dark:bg-blue-900 dark:text-blue-200";
      case "scheduled":
        return "bg-purple-100 text-purple-800 dark:bg-purple-900 dark:text-purple-200";
      case "paused":
        return "bg-gray-100 text-gray-800 dark:bg-gray-900 dark:text-gray-200";
      default:
        return "bg-yellow-100 text-yellow-800 dark:bg-yellow-900 dark:text-yellow-200";
    }
  };

  return (
    <div className="h-full flex flex-col bg-background">
      {/* Header */}
      <div className="flex items-center justify-between px-6 py-4 border-b border-border">
        <div className="flex items-center gap-4">
          <Link
            href="/dashboard"
            className="p-2 hover:bg-accent rounded-md transition"
          >
            <ArrowLeft className="w-5 h-5" />
          </Link>
          <h1 className="text-2xl font-quando font-bold">Job Queue</h1>
        </div>
        <Link
          href="/dashboard"
          className="flex items-center gap-2 px-4 py-2 bg-brand-orange text-white rounded-lg hover:bg-brand-orange-dark transition font-quando"
        >
          <Plus className="w-5 h-5" />
          New Experiment
        </Link>
      </div>

      {/* Stats Bar */}
      <div className="grid grid-cols-4 gap-4 px-6 py-4 border-b border-border">
        <div className="bg-card border border-border rounded-lg p-4">
          <p className="text-xs text-muted-foreground font-quando mb-1">
            Total in Queue
          </p>
          <p className="text-2xl font-quando font-bold">{queue.length}</p>
        </div>
        <div className="bg-card border border-border rounded-lg p-4">
          <p className="text-xs text-muted-foreground font-quando mb-1">
            Running
          </p>
          <p className="text-2xl font-quando font-bold text-blue-600">
            {queue.filter((j) => j.status === "running").length}
          </p>
        </div>
        <div className="bg-card border border-border rounded-lg p-4">
          <p className="text-xs text-muted-foreground font-quando mb-1">
            Scheduled
          </p>
          <p className="text-2xl font-quando font-bold text-purple-600">
            {queue.filter((j) => j.status === "scheduled").length}
          </p>
        </div>
        <div className="bg-card border border-border rounded-lg p-4">
          <p className="text-xs text-muted-foreground font-quando mb-1">
            Paused
          </p>
          <p className="text-2xl font-quando font-bold text-gray-600">
            {queue.filter((j) => j.status === "paused").length}
          </p>
        </div>
      </div>

      {/* Queue List */}
      <div className="flex-1 overflow-auto p-6">
        {queue.length === 0 ? (
          <div className="bg-card border border-border rounded-lg p-12 flex flex-col items-center justify-center text-center">
            <Clock className="w-16 h-16 text-muted-foreground mb-4" />
            <p className="text-lg font-quando font-semibold mb-2">
              No jobs in queue
            </p>
            <p className="text-sm text-muted-foreground mb-6">
              Start creating experiments to build your queue
            </p>
            <Link
              href="/dashboard"
              className="px-6 py-3 bg-brand-orange text-white rounded-lg hover:bg-brand-orange-dark transition font-quando"
            >
              Create Experiment
            </Link>
          </div>
        ) : (
          <div className="space-y-3">
            {queue.map((job, index) => (
              <div
                key={job.id}
                className="bg-card border border-border rounded-lg p-4 hover:border-brand-orange/50 transition"
              >
                <div className="flex items-center gap-4">
                  {/* Priority Number */}
                  <div className="flex flex-col items-center gap-1">
                    <button
                      onClick={() => handlePriorityChange(job.id, "up")}
                      disabled={index === 0}
                      className="p-1 hover:bg-accent rounded disabled:opacity-30"
                    >
                      <MoveUp className="w-4 h-4" />
                    </button>
                    <div className="text-lg font-quando font-bold text-muted-foreground">
                      #{job.priority}
                    </div>
                    <button
                      onClick={() => handlePriorityChange(job.id, "down")}
                      disabled={index === queue.length - 1}
                      className="p-1 hover:bg-accent rounded disabled:opacity-30"
                    >
                      <MoveDown className="w-4 h-4" />
                    </button>
                  </div>

                  {/* Job Details */}
                  <div className="flex-1">
                    <div className="flex items-center gap-3 mb-2">
                      <h3 className="font-quando font-semibold">
                        {job.molecule?.smiles || "Custom Molecule"}
                      </h3>
                      <span
                        className={`px-2 py-1 text-xs rounded font-quando ${getStatusColor(
                          job.status
                        )}`}
                      >
                        {job.status}
                      </span>
                    </div>
                    <div className="flex gap-6 text-xs text-muted-foreground">
                      <span>Method: {job.method}</span>
                      <span>Backend: {job.backend}</span>
                      {job.scheduledTime && (
                        <span className="flex items-center gap-1">
                          <Calendar className="w-3 h-3" />
                          {new Date(job.scheduledTime).toLocaleString()}
                        </span>
                      )}
                    </div>
                  </div>

                  {/* Actions */}
                  <div className="flex gap-2">
                    <button
                      onClick={() => {
                        setSelectedJob(job);
                        setShowScheduleModal(true);
                      }}
                      className="p-2 border border-border rounded-lg hover:bg-accent transition"
                      title="Schedule"
                    >
                      <Calendar className="w-4 h-4" />
                    </button>
                    <button
                      onClick={() => handleTogglePause(job.id)}
                      className="p-2 border border-border rounded-lg hover:bg-accent transition"
                      title={job.status === "paused" ? "Resume" : "Pause"}
                    >
                      {job.status === "paused" ? (
                        <Play className="w-4 h-4" />
                      ) : (
                        <Pause className="w-4 h-4" />
                      )}
                    </button>
                    <button
                      onClick={() => handleDelete(job.id)}
                      className="p-2 border border-border rounded-lg hover:bg-red-100 dark:hover:bg-red-900 transition text-red-600"
                      title="Delete"
                    >
                      <Trash2 className="w-4 h-4" />
                    </button>
                  </div>
                </div>
              </div>
            ))}
          </div>
        )}
      </div>

      {/* Schedule Modal */}
      {showScheduleModal && (
        <div className="fixed inset-0 z-50 flex items-center justify-center backdrop-blur-sm bg-black/20">
          <div className="bg-background border border-border rounded-lg w-full max-w-md p-6">
            <h2 className="text-xl font-quando font-bold mb-4">
              Schedule Experiment
            </h2>

            <div className="space-y-4 mb-6">
              <div>
                <label className="block text-sm font-quando font-medium mb-2">
                  Date
                </label>
                <input
                  type="date"
                  value={scheduleDate}
                  onChange={(e) => setScheduleDate(e.target.value)}
                  className="w-full px-4 py-2 border border-border rounded-lg bg-background focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                />
              </div>

              <div>
                <label className="block text-sm font-quando font-medium mb-2">
                  Time
                </label>
                <input
                  type="time"
                  value={scheduleTime}
                  onChange={(e) => setScheduleTime(e.target.value)}
                  className="w-full px-4 py-2 border border-border rounded-lg bg-background focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                />
              </div>

              {selectedJob && (
                <div className="bg-muted rounded-lg p-3 text-sm">
                  <p className="font-quando">
                    <strong>Molecule:</strong>{" "}
                    {selectedJob.molecule?.smiles || "Custom"}
                  </p>
                  <p className="font-quando">
                    <strong>Method:</strong> {selectedJob.method}
                  </p>
                </div>
              )}
            </div>

            <div className="flex gap-3">
              <button
                onClick={() => {
                  setShowScheduleModal(false);
                  setSelectedJob(null);
                }}
                className="flex-1 px-4 py-2 border border-border rounded-lg hover:bg-accent transition font-quando"
              >
                Cancel
              </button>
              <button
                onClick={handleSchedule}
                disabled={!scheduleDate || !scheduleTime}
                className="flex-1 px-4 py-2 bg-brand-orange text-white rounded-lg hover:bg-brand-orange-dark transition font-quando disabled:opacity-50"
              >
                Schedule
              </button>
            </div>
          </div>
        </div>
      )}
    </div>
  );
}
