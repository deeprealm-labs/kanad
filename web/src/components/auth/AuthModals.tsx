"use client";

import { useState, useEffect } from "react";
import { register, verifyOTP } from "@/lib/api";
import { useAuth } from "@/contexts/AuthContext";
import { useRouter } from "next/navigation";

// Google Sign-In types
declare global {
  interface Window {
    google?: {
      accounts: {
        id: {
          initialize: (config: any) => void;
          renderButton: (element: HTMLElement, config: any) => void;
          prompt: () => void;
        };
      };
    };
  }
}

interface AuthModalsProps {
  showLogin: boolean;
  showRegister: boolean;
  onClose: () => void;
}

export function AuthModals({ showLogin, showRegister, onClose }: AuthModalsProps) {
  const { login, loginWithGoogle } = useAuth();
  const router = useRouter();

  // Login state
  const [loginEmail, setLoginEmail] = useState("");
  const [loginPassword, setLoginPassword] = useState("");
  const [loginError, setLoginError] = useState("");
  const [loginLoading, setLoginLoading] = useState(false);

  // Register state
  const [registerEmail, setRegisterEmail] = useState("");
  const [registerPassword, setRegisterPassword] = useState("");
  const [registerFullName, setRegisterFullName] = useState("");
  const [registerAccessKey, setRegisterAccessKey] = useState("");
  const [registerError, setRegisterError] = useState("");
  const [registerLoading, setRegisterLoading] = useState(false);

  // OTP verification state
  const [showOtpVerification, setShowOtpVerification] = useState(false);
  const [otpEmail, setOtpEmail] = useState("");
  const [otp, setOtp] = useState("");
  const [otpError, setOtpError] = useState("");
  const [otpLoading, setOtpLoading] = useState(false);

  const handleLogin = async (e: React.FormEvent) => {
    e.preventDefault();
    setLoginError("");
    setLoginLoading(true);

    try {
      await login(loginEmail, loginPassword);
      onClose();
      router.push("/dashboard");
    } catch (error: any) {
      setLoginError(error.message || "Login failed");
    } finally {
      setLoginLoading(false);
    }
  };

  const handleRegister = async (e: React.FormEvent) => {
    e.preventDefault();
    setRegisterError("");
    setRegisterLoading(true);

    try {
      const response = await register({
        email: registerEmail,
        password: registerPassword,
        full_name: registerFullName,
        access_key: registerAccessKey,
      });

      // Show OTP verification modal
      setOtpEmail(registerEmail);
      setShowOtpVerification(true);
    } catch (error: any) {
      // Better error message formatting
      let errorMessage = "Registration failed";
      if (error.message) {
        errorMessage = error.message;
      } else if (error.detail) {
        errorMessage = error.detail;
      }
      setRegisterError(errorMessage);
    } finally {
      setRegisterLoading(false);
    }
  };

  const handleVerifyOtp = async (e: React.FormEvent) => {
    e.preventDefault();
    setOtpError("");
    setOtpLoading(true);

    try {
      await verifyOTP({ email: otpEmail, otp });

      // Auto-login after verification
      await login(otpEmail, registerPassword);

      setShowOtpVerification(false);
      onClose();
      router.push("/dashboard");
    } catch (error: any) {
      setOtpError(error.message || "OTP verification failed");
    } finally {
      setOtpLoading(false);
    }
  };

  const handleGoogleSignIn = async (credential: string) => {
    try {
      await loginWithGoogle(credential);
      onClose();
      router.push("/dashboard");
    } catch (error: any) {
      setLoginError(error.message || "Google sign-in failed");
    }
  };

  // Initialize Google Sign-In when login modal is shown
  useEffect(() => {
    if (!showLogin || typeof window === "undefined") return;

    const initializeGoogleSignIn = () => {
      if (window.google) {
        const clientId = process.env.NEXT_PUBLIC_GOOGLE_CLIENT_ID;
        if (!clientId) {
          console.error("Google Client ID not configured");
          return;
        }

        window.google.accounts.id.initialize({
          client_id: clientId,
          callback: (response: any) => {
            handleGoogleSignIn(response.credential);
          },
        });

        const googleButtonDiv = document.getElementById("google-signin-button");
        if (googleButtonDiv) {
          window.google.accounts.id.renderButton(googleButtonDiv, {
            theme: "outline",
            size: "large",
            width: 400,
            text: "signin_with",
          });
        }
      }
    };

    // Wait for Google script to load
    const checkGoogle = setInterval(() => {
      if (window.google) {
        clearInterval(checkGoogle);
        initializeGoogleSignIn();
      }
    }, 100);

    return () => clearInterval(checkGoogle);
  }, [showLogin]);

  if (!showLogin && !showRegister && !showOtpVerification) {
    return null;
  }

  return (
    <>
      {/* Backdrop */}
      <div
        className="fixed inset-0 bg-black/60 backdrop-blur-sm z-50"
        onClick={onClose}
      />

      {/* Login Modal */}
      {showLogin && (
        <div className="fixed inset-0 z-50 flex items-center justify-center p-4">
          <div
            className="bg-background border border-border rounded-lg shadow-2xl max-w-md w-full p-8"
            onClick={(e) => e.stopPropagation()}
          >
            <h2 className="text-3xl font-quando font-bold text-brand-orange mb-6">
              Sign In
            </h2>

            <form onSubmit={handleLogin} className="space-y-4">
              <div>
                <label className="block text-sm font-quando mb-2">Email</label>
                <input
                  type="email"
                  value={loginEmail}
                  onChange={(e) => setLoginEmail(e.target.value)}
                  className="w-full px-4 py-3 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                  placeholder="email@example.com"
                  required
                />
              </div>

              <div>
                <label className="block text-sm font-quando mb-2">Password</label>
                <input
                  type="password"
                  value={loginPassword}
                  onChange={(e) => setLoginPassword(e.target.value)}
                  className="w-full px-4 py-3 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                  placeholder="••••••••"
                  required
                />
              </div>

              {loginError && (
                <div className="bg-red-50 dark:bg-red-900/20 text-red-600 dark:text-red-400 px-4 py-3 rounded-md text-sm">
                  {loginError}
                </div>
              )}

              <button
                type="submit"
                disabled={loginLoading}
                className="w-full bg-brand-orange text-white py-3 rounded-md hover:bg-brand-orange/90 transition font-quando disabled:opacity-50"
              >
                {loginLoading ? "Signing in..." : "Sign In"}
              </button>
            </form>

            <div className="mt-6">
              <div className="relative">
                <div className="absolute inset-0 flex items-center">
                  <div className="w-full border-t border-border"></div>
                </div>
                <div className="relative flex justify-center text-sm">
                  <span className="px-2 bg-background text-muted-foreground">
                    Or continue with
                  </span>
                </div>
              </div>

              {/* Google Sign-In Button */}
              <div id="google-signin-button" className="mt-4 w-full flex justify-center"></div>
            </div>

            <button
              onClick={onClose}
              className="absolute top-4 right-4 text-muted-foreground hover:text-foreground"
            >
              ✕
            </button>
          </div>
        </div>
      )}

      {/* Register Modal */}
      {showRegister && !showOtpVerification && (
        <div className="fixed inset-0 z-50 flex items-center justify-center p-4">
          <div
            className="bg-background border border-border rounded-lg shadow-2xl max-w-md w-full p-8 max-h-[90vh] overflow-y-auto"
            onClick={(e) => e.stopPropagation()}
          >
            <h2 className="text-3xl font-quando font-bold text-brand-orange mb-6">
              Create Account
            </h2>

            <form onSubmit={handleRegister} className="space-y-4">
              <div>
                <label className="block text-sm font-quando mb-2">Full Name</label>
                <input
                  type="text"
                  value={registerFullName}
                  onChange={(e) => setRegisterFullName(e.target.value)}
                  className="w-full px-4 py-3 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                  placeholder="John Doe"
                  required
                />
              </div>

              <div>
                <label className="block text-sm font-quando mb-2">Email</label>
                <input
                  type="email"
                  value={registerEmail}
                  onChange={(e) => setRegisterEmail(e.target.value)}
                  className="w-full px-4 py-3 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                  placeholder="email@example.com"
                  required
                />
              </div>

              <div>
                <label className="block text-sm font-quando mb-2">Password</label>
                <input
                  type="password"
                  value={registerPassword}
                  onChange={(e) => setRegisterPassword(e.target.value)}
                  className="w-full px-4 py-3 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                  placeholder="••••••••"
                  required
                  minLength={8}
                />
                <p className="text-xs text-muted-foreground mt-1">
                  At least 8 characters, include uppercase, lowercase, and number
                </p>
              </div>

              <div>
                <label className="block text-sm font-quando mb-2">Access Key</label>
                <input
                  type="text"
                  value={registerAccessKey}
                  onChange={(e) => setRegisterAccessKey(e.target.value)}
                  className="w-full px-4 py-3 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-quando"
                  placeholder="ACCESS-KEY-2024"
                  required
                />
                <p className="text-xs text-muted-foreground mt-1">
                  Request an access key to join the platform
                </p>
              </div>

              {registerError && (
                <div className="bg-red-50 dark:bg-red-900/20 text-red-600 dark:text-red-400 px-4 py-3 rounded-md text-sm">
                  {registerError}
                </div>
              )}

              <button
                type="submit"
                disabled={registerLoading}
                className="w-full bg-brand-orange text-white py-3 rounded-md hover:bg-brand-orange/90 transition font-quando disabled:opacity-50"
              >
                {registerLoading ? "Creating account..." : "Create Account"}
              </button>
            </form>

            <button
              onClick={onClose}
              className="absolute top-4 right-4 text-muted-foreground hover:text-foreground"
            >
              ✕
            </button>
          </div>
        </div>
      )}

      {/* OTP Verification Modal */}
      {showOtpVerification && (
        <div className="fixed inset-0 z-50 flex items-center justify-center p-4">
          <div
            className="bg-background border border-border rounded-lg shadow-2xl max-w-md w-full p-8"
            onClick={(e) => e.stopPropagation()}
          >
            <h2 className="text-3xl font-quando font-bold text-brand-orange mb-6">
              Verify Email
            </h2>

            <p className="text-muted-foreground mb-6 font-quando">
              We&apos;ve sent a verification code to <strong>{otpEmail}</strong>. Please enter it below.
            </p>

            <form onSubmit={handleVerifyOtp} className="space-y-4">
              <div>
                <label className="block text-sm font-quando mb-2">Verification Code</label>
                <input
                  type="text"
                  value={otp}
                  onChange={(e) => setOtp(e.target.value)}
                  className="w-full px-4 py-3 border border-input bg-background rounded-md focus:outline-none focus:ring-2 focus:ring-brand-orange font-mono text-2xl text-center tracking-widest"
                  placeholder="000000"
                  required
                  maxLength={6}
                />
              </div>

              {otpError && (
                <div className="bg-red-50 dark:bg-red-900/20 text-red-600 dark:text-red-400 px-4 py-3 rounded-md text-sm">
                  {otpError}
                </div>
              )}

              <button
                type="submit"
                disabled={otpLoading}
                className="w-full bg-brand-orange text-white py-3 rounded-md hover:bg-brand-orange/90 transition font-quando disabled:opacity-50"
              >
                {otpLoading ? "Verifying..." : "Verify Email"}
              </button>
            </form>

            <button
              onClick={() => {
                setShowOtpVerification(false);
                onClose();
              }}
              className="absolute top-4 right-4 text-muted-foreground hover:text-foreground"
            >
              ✕
            </button>
          </div>
        </div>
      )}
    </>
  );
}
