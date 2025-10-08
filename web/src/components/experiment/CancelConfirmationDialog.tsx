"use client";

import { AlertTriangle, Loader2 } from "lucide-react";

interface CancelConfirmationDialogProps {
  isOpen: boolean;
  onClose: () => void;
  onConfirm: () => void;
  experimentName: string;
  isRunning: boolean;
  isLoading?: boolean;
}

export default function CancelConfirmationDialog({
  isOpen,
  onClose,
  onConfirm,
  experimentName,
  isRunning,
  isLoading = false,
}: CancelConfirmationDialogProps) {
  if (!isOpen) return null;

  return (
    <div className="fixed inset-0 z-50 flex items-center justify-center backdrop-blur-sm bg-black/40">
      <div className="bg-background border-2 border-red-500 rounded-lg w-full max-w-md p-6 shadow-2xl">
        <div className="flex items-start gap-4 mb-4">
          <div className="p-2 bg-red-100 dark:bg-red-900/30 rounded-full">
            <AlertTriangle className="w-6 h-6 text-red-600" />
          </div>
          <div className="flex-1">
            <h2 className="text-xl font-quando font-bold mb-2">
              Cancel Experiment
            </h2>
            <p className="text-sm text-muted-foreground">
              Are you sure you want to cancel this experiment?
            </p>
          </div>
        </div>

        <div className="bg-muted rounded-lg p-4 mb-4">
          <p className="font-quando font-semibold text-sm mb-1">
            Experiment:
          </p>
          <p className="font-mono text-sm break-all">{experimentName}</p>
        </div>

        {isRunning && (
          <div className="bg-yellow-50 dark:bg-yellow-900/20 border border-yellow-200 dark:border-yellow-800 rounded-lg p-3 mb-4">
            <div className="flex items-start gap-2">
              <AlertTriangle className="w-4 h-4 text-yellow-600 mt-0.5 flex-shrink-0" />
              <p className="text-xs text-yellow-800 dark:text-yellow-200">
                This experiment is currently running. Cancellation will stop it
                immediately and may save partial results if available.
              </p>
            </div>
          </div>
        )}

        <div className="bg-red-50 dark:bg-red-900/20 border border-red-200 dark:border-red-800 rounded-lg p-3 mb-6">
          <p className="text-xs text-red-800 dark:text-red-200 font-quando">
            This action cannot be undone.
          </p>
        </div>

        <div className="flex gap-3">
          <button
            onClick={onClose}
            disabled={isLoading}
            className="flex-1 px-4 py-2.5 border-2 border-border rounded-lg hover:bg-accent transition font-quando font-medium disabled:opacity-50"
          >
            Keep Running
          </button>
          <button
            onClick={onConfirm}
            disabled={isLoading}
            className="flex-1 px-4 py-2.5 bg-red-600 text-white rounded-lg hover:bg-red-700 transition font-quando font-medium disabled:opacity-50 flex items-center justify-center gap-2"
          >
            {isLoading ? (
              <>
                <Loader2 className="w-4 h-4 animate-spin" />
                Cancelling...
              </>
            ) : (
              "Cancel Experiment"
            )}
          </button>
        </div>
      </div>
    </div>
  );
}
