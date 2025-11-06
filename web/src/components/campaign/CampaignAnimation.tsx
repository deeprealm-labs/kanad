"use client";

/**
 * Campaign Animation Component
 *
 * Shows a visual representation of multiple experiments running in sequence
 */

export default function CampaignAnimation({
  status,
  progress
}: {
  status: string;
  progress: { total: number; completed: number; running: number; pending: number };
}) {
  return (
    <div className="relative w-full h-full flex flex-col items-center justify-center">
      {/* Sequential Experiment Flow */}
      <div className="relative w-full max-w-2xl">
        {/* Timeline Line */}
        <div className="absolute top-1/2 left-0 right-0 h-1 bg-muted -translate-y-1/2" />
        <div
          className="absolute top-1/2 left-0 h-1 bg-gradient-to-r from-brand-orange to-brand-yellow -translate-y-1/2 transition-all duration-500"
          style={{ width: `${(progress.completed / progress.total) * 100}%` }}
        />

        {/* Experiment Nodes */}
        <div className="relative flex justify-between items-center py-12">
          {Array.from({ length: Math.min(progress.total, 5) }).map((_, i) => {
            const isCompleted = i < progress.completed;
            const isRunning = i === progress.completed && progress.running > 0;
            const isPending = i > progress.completed;

            return (
              <div key={i} className="flex flex-col items-center gap-3">
                {/* Node */}
                <div
                  className={`w-16 h-16 rounded-full border-4 flex items-center justify-center font-quanto font-bold text-lg transition-all duration-300 ${
                    isCompleted
                      ? "bg-green-500 border-green-600 text-white shadow-lg shadow-green-500/50"
                      : isRunning
                      ? "bg-brand-orange border-brand-yellow text-white animate-pulse shadow-lg shadow-brand-orange/50"
                      : "bg-muted border-border text-muted-foreground"
                  }`}
                >
                  {i + 1}
                </div>

                {/* Label */}
                <div className="text-xs font-quanto text-center">
                  <div className={`font-semibold ${isCompleted ? "text-green-600" : isRunning ? "text-brand-orange" : "text-muted-foreground"}`}>
                    Exp {i + 1}
                  </div>
                  <div className="text-muted-foreground text-[10px] mt-0.5">
                    {isCompleted ? "Complete" : isRunning ? "Running" : "Pending"}
                  </div>
                </div>

                {/* Running indicator */}
                {isRunning && (
                  <div className="absolute -top-8 animate-bounce">
                    <div className="w-3 h-3 bg-brand-orange rounded-full" />
                  </div>
                )}
              </div>
            );
          })}
        </div>

        {/* Show "+N more" if there are more than 5 experiments */}
        {progress.total > 5 && (
          <div className="text-center text-sm text-muted-foreground mt-4">
            + {progress.total - 5} more experiments
          </div>
        )}
      </div>

      {/* Status text */}
      <div className="absolute bottom-0 text-center">
        <p className="text-lg font-quando font-semibold text-foreground">
          {status === "queued" && "Preparing Campaign..."}
          {status === "running" && `Executing Experiment ${progress.completed + 1} of ${progress.total}`}
          {status === "completed" && "Campaign Complete!"}
          {status === "failed" && "Campaign Failed"}
        </p>
        <p className="text-sm text-muted-foreground mt-2">
          {status === "running" && `${progress.completed} completed, ${progress.pending} remaining`}
        </p>
      </div>
    </div>
  );
}
