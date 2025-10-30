import type { Metadata } from "next";
import { Quando } from "next/font/google";
import "./globals.css";
import { ThemeProvider } from "@/components/theme-provider";
import { ToastProvider } from "@/components/ui/toast";
import { AuthProvider } from "@/contexts/AuthContext";
import Script from "next/script";

const quando = Quando({
  weight: "400",
  subsets: ["latin"],
  variable: "--font-quando",
  display: "swap",
});

export const metadata: Metadata = {
  title: "Kanad - Quantum Chemistry Platform",
  description: "The place for creation, exploration, and invention in quantum chemistry",
};

export default function RootLayout({
  children,
}: Readonly<{
  children: React.ReactNode;
}>) {
  return (
    <html lang="en" suppressHydrationWarning>
      <head>
        <Script
          src="/smiles-drawer.min.js"
          strategy="beforeInteractive"
        />
        <Script
          src="https://accounts.google.com/gsi/client"
          strategy="beforeInteractive"
        />
      </head>
      <body className={`${quando.variable} antialiased`} style={{ fontFamily: 'var(--font-quando), serif' }}>
        <ThemeProvider
          attribute="class"
          defaultTheme="light"
          enableSystem
          disableTransitionOnChange
        >
          <AuthProvider>
            <ToastProvider>
              {children}
            </ToastProvider>
          </AuthProvider>
        </ThemeProvider>
      </body>
    </html>
  );
}
