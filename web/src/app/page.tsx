"use client";

import { useState, useEffect } from "react";
import Image from "next/image";
import { AuthModals } from "@/components/auth/AuthModals";
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

export default function Home() {
  const [showLoginModal, setShowLoginModal] = useState(false);
  const [showRegisterModal, setShowRegisterModal] = useState(false);
  const { isAuthenticated, loginWithGoogle } = useAuth();
  const router = useRouter();

  // Redirect to dashboard if authenticated (using useEffect to avoid render issues)
  useEffect(() => {
    if (isAuthenticated) {
      router.push("/dashboard");
    }
  }, [isAuthenticated, router]);

  // Initialize Google Sign-In for home page button
  useEffect(() => {
    // Don't initialize if already authenticated
    if (typeof window === "undefined" || isAuthenticated) return;

    const initializeGoogleSignIn = () => {
      if (window.google) {
        window.google.accounts.id.initialize({
          client_id: "336249296791-njvoeikug1vurqk1u6i1uqej4jufpcnf.apps.googleusercontent.com",
          callback: async (response: any) => {
            try {
              await loginWithGoogle(response.credential);
              router.push("/dashboard");
            } catch (error) {
              console.error("Google Sign-In failed:", error);
              setShowLoginModal(true);
            }
          },
        });
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
  }, [loginWithGoogle, router, isAuthenticated]);

  const handleGoogleSignIn = () => {
    if (window.google) {
      window.google.accounts.id.prompt();
    }
  };

  return (
    <div className="min-h-screen flex flex-col md:flex-row">
      {/* Left Side - Theme-Aware Background */}
      <div className="w-full md:w-1/2 bg-background flex items-center justify-center p-8 md:p-16">
        <div className="max-w-lg w-full space-y-8">
          {/* Welcome Text */}
          <div className="space-y-2">
            <h1 className="text-4xl md:text-5xl lg:text-6xl font-quando leading-tight">
              <span className="block">welcome</span>
              <span className="block text-2xl md:text-3xl lg:text-4xl text-muted-foreground">
                to the place for
              </span>
            </h1>
            <div className="text-4xl md:text-5xl lg:text-6xl font-quando leading-tight">
              <span className="block">creation</span>
              <span className="block">exploration</span>
              <span className="block flex items-center gap-2">
                <span>&</span>
              </span>
              <span className="block">invention</span>
            </div>
          </div>

          {/* Authentication Buttons */}
          <div className="space-y-4">
            {/* Sign In / Register Buttons */}
            <div className="flex gap-3">
              <button
                onClick={() => setShowLoginModal(true)}
                className="flex-1 px-6 py-4 border-2 border-border rounded-md hover:bg-accent transition font-quando text-lg"
              >
                Sign In
              </button>
              <button
                onClick={() => setShowRegisterModal(true)}
                className="flex-1 px-6 py-4 bg-black text-white rounded-md hover:bg-brand-orange/90 transition font-quando text-lg"
              >
                Create Account
              </button>
            </div>

            {/* Google Sign In Button - Temporarily disabled until OAuth origins propagate */}
            {/*
            <button
              onClick={handleGoogleSignIn}
              className="w-full px-6 py-4 border-2 border-border rounded-md hover:bg-accent transition font-quando text-lg flex items-center justify-center gap-3"
            >
              <svg className="w-5 h-5" viewBox="0 0 24 24">
                <path fill="#4285F4" d="M22.56 12.25c0-.78-.07-1.53-.2-2.25H12v4.26h5.92c-.26 1.37-1.04 2.53-2.21 3.31v2.77h3.57c2.08-1.92 3.28-4.74 3.28-8.09z"/>
                <path fill="#34A853" d="M12 23c2.97 0 5.46-.98 7.28-2.66l-3.57-2.77c-.98.66-2.23 1.06-3.71 1.06-2.86 0-5.29-1.93-6.16-4.53H2.18v2.84C3.99 20.53 7.7 23 12 23z"/>
                <path fill="#FBBC05" d="M5.84 14.09c-.22-.66-.35-1.36-.35-2.09s.13-1.43.35-2.09V7.07H2.18C1.43 8.55 1 10.22 1 12s.43 3.45 1.18 4.93l2.85-2.22.81-.62z"/>
                <path fill="#EA4335" d="M12 5.38c1.62 0 3.06.56 4.21 1.64l3.15-3.15C17.45 2.09 14.97 1 12 1 7.7 1 3.99 3.47 2.18 7.07l3.66 2.84c.87-2.6 3.3-4.53 6.16-4.53z"/>
              </svg>
              Sign in with Google
            </button>
            */}

          </div>
        </div>
      </div>

      {/* Right Side - Black Background with Artistic Image */}
      <div className="w-full md:w-1/2 bg-black flex items-center justify-center relative p-4 overflow-hidden min-h-screen bg-white">
        {/* Artistic Image */}
        {/* Artistic Image */}
        <div className="relative w-full h-full flex items-center justify-center ">
          <Image
            src="/image.webp"
            alt="Kanad - Creation and Exploration"
            fill
            className="object-cover rounded-xl"
            priority
          />

          {/* Vertical "deeprealm" text */}
          {/* <div className="absolute left-0 top-1/3 -translate-y-1/2 text-white text-sm tracking-widest rotate-90">
            deeprealm
          </div> */}

          {/* Kanad Logo - Much Bigger, Covers Image Width, Centered */}
          <div className="absolute bottom-12 left-0 right-0 px-8 text-center">
            <h1 className="font-bietro text-[clamp(3rem,12vw,16rem)] text-brand-orange font-bold leading-none relative w-full whitespace-nowrap">
              kanad
            </h1>
          </div>
        </div>
      </div>

      {/* Authentication Modals */}
      <AuthModals
        showLogin={showLoginModal}
        showRegister={showRegisterModal}
        onClose={() => {
          setShowLoginModal(false);
          setShowRegisterModal(false);
        }}
      />
    </div>
  );
}
