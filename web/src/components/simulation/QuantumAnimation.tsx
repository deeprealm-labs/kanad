"use client";

import Lottie from "lottie-react";

/**
 * Quantum Animation Component
 *
 * Uses a quantum/science-themed Lottie animation to show experiment progress.
 * Fallback to CSS animation if Lottie animation is not available.
 */

export default function QuantumAnimation({ status }: { status: string }) {
  // Use a simple quantum-inspired CSS animation
  // In production, you would replace this with a proper Lottie JSON animation

  return (
    <div className="relative w-full h-full flex items-center justify-center">
      {/* Animated Quantum Orbitals */}
      <div className="relative w-64 h-64">
        {/* Central nucleus */}
        <div className="absolute top-1/2 left-1/2 transform -translate-x-1/2 -translate-y-1/2">
          <div className="w-12 h-12 bg-brand-orange rounded-full animate-pulse shadow-lg shadow-brand-orange/50" />
        </div>

        {/* Orbital rings */}
        <div className="absolute inset-0 animate-spin-slow">
          <div className="absolute top-1/2 left-1/2 transform -translate-x-1/2 -translate-y-1/2 w-40 h-40 border-2 border-brand-orange/30 rounded-full" />
          <div className="absolute top-2 left-1/2 transform -translate-x-1/2 w-3 h-3 bg-blue-500 rounded-full shadow-lg shadow-blue-500/50" />
        </div>

        <div className="absolute inset-0 animate-spin-slow-reverse" style={{ animationDelay: '0.5s' }}>
          <div className="absolute top-1/2 left-1/2 transform -translate-x-1/2 -translate-y-1/2 w-52 h-52 border-2 border-purple-500/30 rounded-full" />
          <div className="absolute bottom-2 left-1/2 transform -translate-x-1/2 w-3 h-3 bg-purple-500 rounded-full shadow-lg shadow-purple-500/50" />
        </div>

        <div className="absolute inset-0 animate-spin-medium">
          <div className="absolute top-1/2 left-1/2 transform -translate-x-1/2 -translate-y-1/2 w-48 h-48 border-2 border-green-500/30 rounded-full" />
          <div className="absolute top-1/2 right-2 transform -translate-y-1/2 w-3 h-3 bg-green-500 rounded-full shadow-lg shadow-green-500/50" />
        </div>

        {/* Particle effects */}
        <div className="absolute inset-0 animate-pulse-slow">
          <div className="absolute top-1/4 left-1/4 w-1 h-1 bg-white rounded-full opacity-60" />
          <div className="absolute top-3/4 left-1/3 w-1 h-1 bg-white rounded-full opacity-60" />
          <div className="absolute top-1/3 left-3/4 w-1 h-1 bg-white rounded-full opacity-60" />
          <div className="absolute top-2/3 left-2/3 w-1 h-1 bg-white rounded-full opacity-60" />
        </div>
      </div>

      {/* Status text */}
      <div className="absolute bottom-0 text-center">
        <p className="text-lg font-quando font-semibold text-foreground">
          {status === "queued" && "Preparing Quantum Circuit..."}
          {status === "running" && "Optimizing Quantum Parameters..."}
          {status === "completed" && "Experiment Complete!"}
          {status === "failed" && "Experiment Failed"}
        </p>
        <p className="text-sm text-muted-foreground mt-2">
          {status === "running" && "This may take a few minutes"}
        </p>
      </div>

      {/* Add custom animations to global CSS */}
      <style jsx>{`
        @keyframes spin-slow {
          from { transform: rotate(0deg); }
          to { transform: rotate(360deg); }
        }
        @keyframes spin-slow-reverse {
          from { transform: rotate(360deg); }
          to { transform: rotate(0deg); }
        }
        @keyframes spin-medium {
          from { transform: rotate(0deg); }
          to { transform: rotate(360deg); }
        }
        @keyframes pulse-slow {
          0%, 100% { opacity: 0.3; }
          50% { opacity: 0.8; }
        }

        .animate-spin-slow {
          animation: spin-slow 8s linear infinite;
        }
        .animate-spin-slow-reverse {
          animation: spin-slow-reverse 10s linear infinite;
        }
        .animate-spin-medium {
          animation: spin-medium 6s linear infinite;
        }
        .animate-pulse-slow {
          animation: pulse-slow 3s ease-in-out infinite;
        }
      `}</style>
    </div>
  );
}
