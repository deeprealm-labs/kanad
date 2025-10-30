"use client";

import React, { createContext, useContext, useState, useEffect } from "react";
import { login as apiLogin, logout as apiLogout, refreshToken, googleAuth } from "@/lib/api";

interface User {
  id: string;
  email: string;
  full_name: string;
  is_verified: boolean;
  role: string;
}

interface AuthContextType {
  user: User | null;
  accessToken: string | null;
  isAuthenticated: boolean;
  isLoading: boolean;
  login: (email: string, password: string) => Promise<void>;
  logout: () => Promise<void>;
  loginWithGoogle: (credential: string) => Promise<void>;
  setUser: (user: User | null) => void;
  setAccessToken: (token: string | null) => void;
}

const AuthContext = createContext<AuthContextType | undefined>(undefined);

export function AuthProvider({ children }: { children: React.ReactNode }) {
  const [user, setUser] = useState<User | null>(null);
  const [accessToken, setAccessToken] = useState<string | null>(null);
  const [isLoading, setIsLoading] = useState(true);

  // Load auth state from localStorage on mount
  useEffect(() => {
    const storedUser = localStorage.getItem("kanad_user");
    const storedAccessToken = localStorage.getItem("kanad_access_token");
    const storedRefreshToken = localStorage.getItem("kanad_refresh_token");

    if (storedUser && storedAccessToken) {
      setUser(JSON.parse(storedUser));
      setAccessToken(storedAccessToken);
    }

    setIsLoading(false);

    // Set up token refresh interval (refresh every 14 minutes, token expires in 15)
    if (storedRefreshToken) {
      const refreshInterval = setInterval(async () => {
        try {
          const response = await refreshToken(storedRefreshToken);
          setAccessToken(response.access_token);
          localStorage.setItem("kanad_access_token", response.access_token);
        } catch (error) {
          console.error("Token refresh failed:", error);
          // If refresh fails, log out the user
          await handleLogout();
        }
      }, 14 * 60 * 1000); // 14 minutes

      return () => clearInterval(refreshInterval);
    }
  }, []);

  const handleLogin = async (email: string, password: string) => {
    const response = await apiLogin({ email, password });

    const userData: User = {
      id: response.user.id,
      email: response.user.email,
      full_name: response.user.full_name,
      is_verified: response.user.is_verified,
      role: response.user.role,
    };

    setUser(userData);
    setAccessToken(response.access_token);

    // Store in localStorage
    localStorage.setItem("kanad_user", JSON.stringify(userData));
    localStorage.setItem("kanad_access_token", response.access_token);
    localStorage.setItem("kanad_refresh_token", response.refresh_token);
  };

  const handleGoogleLogin = async (credential: string) => {
    const response = await googleAuth(credential);

    const userData: User = {
      id: response.user.id,
      email: response.user.email,
      full_name: response.user.full_name,
      is_verified: response.user.is_verified,
      role: response.user.role,
    };

    setUser(userData);
    setAccessToken(response.access_token);

    // Store in localStorage
    localStorage.setItem("kanad_user", JSON.stringify(userData));
    localStorage.setItem("kanad_access_token", response.access_token);
    localStorage.setItem("kanad_refresh_token", response.refresh_token);
  };

  const handleLogout = async () => {
    const storedRefreshToken = localStorage.getItem("kanad_refresh_token");

    if (storedRefreshToken) {
      try {
        await apiLogout(storedRefreshToken);
      } catch (error) {
        console.error("Logout API call failed:", error);
      }
    }

    setUser(null);
    setAccessToken(null);

    // Clear localStorage
    localStorage.removeItem("kanad_user");
    localStorage.removeItem("kanad_access_token");
    localStorage.removeItem("kanad_refresh_token");
  };

  return (
    <AuthContext.Provider
      value={{
        user,
        accessToken,
        isAuthenticated: !!user,
        isLoading,
        login: handleLogin,
        logout: handleLogout,
        loginWithGoogle: handleGoogleLogin,
        setUser,
        setAccessToken,
      }}
    >
      {children}
    </AuthContext.Provider>
  );
}

export function useAuth() {
  const context = useContext(AuthContext);
  if (context === undefined) {
    throw new Error("useAuth must be used within an AuthProvider");
  }
  return context;
}
