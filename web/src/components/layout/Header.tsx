"use client";

import { useState } from "react";
import { Settings, LogOut, User as UserIcon } from "lucide-react";
import { ThemeToggle } from "@/components/theme-toggle";
import SettingsModal from "@/components/settings/SettingsModal";
import UserSettingsModal from "@/components/settings/UserSettingsModal";
import { useAuth } from "@/contexts/AuthContext";
import { useRouter } from "next/navigation";

export default function Header() {
  const [showSettings, setShowSettings] = useState(false);
  const [showUserSettings, setShowUserSettings] = useState(false);
  const { logout, user } = useAuth();
  const router = useRouter();

  const handleLogout = async () => {
    await logout();
    router.push("/");
  };

  return (
    <>
      <header className="bg-background px-8 py-6">
        <div className="flex justify-between items-center">
          {/* User Info */}
          <div className="flex items-center gap-3">
            <div className="w-8 h-8 bg-brand-orange rounded-full flex items-center justify-center text-white font-quando font-bold">
              {user?.full_name?.charAt(0).toUpperCase() || "U"}
            </div>
            <div>
              <p className="text-sm font-quando font-semibold">{user?.full_name || "User"}</p>
              <p className="text-xs text-muted-foreground font-quando">{user?.email}</p>
            </div>
          </div>

          {/* Actions */}
          <div className="flex gap-3">
            <ThemeToggle />
            <button
              onClick={() => setShowUserSettings(true)}
              className="flex items-center gap-2 px-4 py-2 border border-border rounded-md hover:bg-accent transition font-quando"
              title="Account"
            >
              <UserIcon className="w-5 h-5" />
              <span className="text-sm">account</span>
            </button>
            <button
              onClick={() => setShowSettings(true)}
              className="flex items-center gap-2 px-4 py-2 border border-border rounded-md hover:bg-accent transition font-quando"
              title="Settings"
            >
              <Settings className="w-5 h-5" />
              <span className="text-sm">settings</span>
            </button>
            <button
              onClick={handleLogout}
              className="flex items-center gap-2 px-4 py-2 bg-red-500/10 text-red-600 dark:text-red-400 border border-red-500/20 rounded-md hover:bg-red-500/20 transition font-quando"
              title="Logout"
            >
              <LogOut className="w-5 h-5" />
              <span className="text-sm">logout</span>
            </button>
          </div>
        </div>
      </header>

      <SettingsModal
        isOpen={showSettings}
        onClose={() => setShowSettings(false)}
      />

      <UserSettingsModal
        isOpen={showUserSettings}
        onClose={() => setShowUserSettings(false)}
      />
    </>
  );
}
