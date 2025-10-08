import type { Metadata } from "next";
import { Quando } from "next/font/google";
import "./globals.css";
import { ThemeProvider } from "@/components/theme-provider";

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
      <body className={`${quando.variable} antialiased`} style={{ fontFamily: 'var(--font-quando), serif' }}>
        <ThemeProvider
          attribute="class"
          defaultTheme="light"
          enableSystem
          disableTransitionOnChange
        >
          {children}
        </ThemeProvider>
      </body>
    </html>
  );
}
