"use client";

import { Loader2 } from "lucide-react";

// ===== Full Page Loader =====
export function FullPageLoader({ message }: { message?: string }) {
  return (
    <div className="fixed inset-0 z-50 flex items-center justify-center bg-background/80 backdrop-blur-sm">
      <div className="flex flex-col items-center gap-4">
        <Loader2 className="w-12 h-12 text-brand-orange animate-spin" />
        {message && <p className="text-sm font-quando text-muted-foreground">{message}</p>}
      </div>
    </div>
  );
}

// ===== Spinner =====
export function Spinner({ size = "default", className = "" }: { size?: "sm" | "default" | "lg"; className?: string }) {
  const sizeClasses = {
    sm: "w-4 h-4",
    default: "w-6 h-6",
    lg: "w-8 h-8",
  };

  return (
    <Loader2 className={`${sizeClasses[size]} text-brand-orange animate-spin ${className}`} />
  );
}

// ===== Loading Button =====
export function LoadingButton({
  loading,
  children,
  disabled,
  ...props
}: {
  loading: boolean;
  children: React.ReactNode;
  disabled?: boolean;
  [key: string]: any;
}) {
  return (
    <button {...props} disabled={disabled || loading}>
      {loading ? (
        <div className="flex items-center gap-2">
          <Spinner size="sm" />
          <span>Loading...</span>
        </div>
      ) : (
        children
      )}
    </button>
  );
}

// ===== Skeleton Loader =====
export function Skeleton({ className = "" }: { className?: string }) {
  return (
    <div className={`animate-pulse bg-muted rounded ${className}`} />
  );
}

// ===== Card Skeleton =====
export function CardSkeleton() {
  return (
    <div className="bg-card border border-border rounded-lg p-4 space-y-3">
      <Skeleton className="h-5 w-3/4" />
      <Skeleton className="h-4 w-1/2" />
      <div className="grid grid-cols-2 gap-2 mt-4">
        <Skeleton className="h-4 w-full" />
        <Skeleton className="h-4 w-full" />
      </div>
    </div>
  );
}

// ===== Table Skeleton =====
export function TableSkeleton({ rows = 5 }: { rows?: number }) {
  return (
    <div className="space-y-3">
      {Array.from({ length: rows }).map((_, i) => (
        <div key={i} className="bg-card border border-border rounded-lg p-4">
          <div className="flex items-center gap-4">
            <Skeleton className="h-10 w-10 rounded-full" />
            <div className="flex-1 space-y-2">
              <Skeleton className="h-4 w-3/4" />
              <Skeleton className="h-3 w-1/2" />
            </div>
            <Skeleton className="h-8 w-20" />
          </div>
        </div>
      ))}
    </div>
  );
}

// ===== Loading Overlay =====
export function LoadingOverlay({
  visible,
  message,
}: {
  visible: boolean;
  message?: string;
}) {
  if (!visible) return null;

  return (
    <div className="absolute inset-0 z-10 flex items-center justify-center bg-background/80 backdrop-blur-sm rounded-lg">
      <div className="flex flex-col items-center gap-3">
        <Spinner size="lg" />
        {message && <p className="text-sm font-quando text-muted-foreground">{message}</p>}
      </div>
    </div>
  );
}

// ===== Progress Bar =====
export function ProgressBar({
  progress,
  message,
  showPercentage = true,
}: {
  progress: number;
  message?: string;
  showPercentage?: boolean;
}) {
  return (
    <div className="w-full space-y-2">
      {message && <p className="text-sm font-quando text-muted-foreground">{message}</p>}
      <div className="relative w-full bg-muted rounded-full h-2 overflow-hidden">
        <div
          className="h-full bg-brand-orange transition-all duration-300 ease-out"
          style={{ width: `${Math.min(100, Math.max(0, progress))}%` }}
        />
      </div>
      {showPercentage && (
        <p className="text-xs text-right font-quando text-muted-foreground">
          {Math.round(progress)}%
        </p>
      )}
    </div>
  );
}
