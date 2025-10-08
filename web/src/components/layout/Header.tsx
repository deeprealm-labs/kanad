"use client";

import { useState } from "react";
import { Settings } from "lucide-react";
import { ThemeToggle } from "@/components/theme-toggle";

export default function Header() {
  const [showSettings, setShowSettings] = useState(false);

  return (
    <header className="bg-background  px-8 py-6">
      <div className="flex justify-end gap-3">
        <ThemeToggle />
        <button
          onClick={() => setShowSettings(!showSettings)}
          className="flex items-center gap-2 px-4 py-2 border border-border rounded-md hover:bg-accent transition font-quando"
          title="Settings"
        >
          <Settings className="w-5 h-5" />
          <span className="text-sm">settings</span>
        </button>
      </div>
    </header>
  );
}
