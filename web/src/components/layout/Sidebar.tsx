"use client";

import Link from "next/link";
import { useState } from "react";

export default function Sidebar() {
  const [isOpen, setIsOpen] = useState(false);

  return (
    <>
      {/* Mobile Hamburger Button */}
      <button
        onClick={() => setIsOpen(!isOpen)}
        className="md:hidden fixed top-4 left-4 z-50 p-2 bg-black text-white rounded-md"
      >
        <svg className="w-6 h-6" fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 6h16M4 12h16M4 18h16" />
        </svg>
      </button>

      {/* Sidebar */}
      <aside
        className={`
          ${isOpen ? "translate-x-0" : "-translate-x-full"}
          md:translate-x-0
          fixed md:relative
          z-40
          w-72
          h-screen
          bg-black dark:bg-black
       31 -            text-white dark:text-white
          flex
          flex-col
          transition-transform
          duration-300
          border-r border-border
          rounded-lg
        `}
      >
        {/* Navigation */}
        <nav className="flex-1 p-4 pt-6">
          <div className="space-y-1">
            <a
              href="/dashboard"
              className="block px-4 py-2 text-gray-400 hover:text-white hover:bg-gray-800 rounded-md transition font-quando"
            >
              Dashboard
            </a>
            <a
              href="/dashboard/history"
              className="block px-4 py-2 text-gray-400 hover:text-white hover:bg-gray-800 rounded-md transition font-quando"
            >
             History
            </a>
            <a
              href="/dashboard/queue"
              className="block px-4 py-2 text-gray-400 hover:text-white hover:bg-gray-800 rounded-md transition font-quando"
            >
              Job Queue
            </a>
            <a
              href="/dashboard/backend"
              className="block px-4 py-2 text-gray-400 hover:text-white hover:bg-gray-800 rounded-md transition font-quando"
            >
              Backend
            </a>
            <div className="border-t border-gray-800 my-2"></div>
            <a
              href="#"
              className="block px-4 py-2 text-gray-400 hover:text-white hover:bg-gray-800 rounded-md transition font-quando"
            >
              Docs
            </a>
            <a
              href="#"
              className="block px-4 py-2 text-gray-400 hover:text-white hover:bg-gray-800 rounded-md transition font-quando"
            >
              Tutorials
            </a>
          </div>
        </nav>

        {/* Kanad Logo at Bottom - Much Bigger, Covers Sidebar Width */}
        <div className="p-4 border-t border-gray-800">
          <h1 className="font-bietro text-7xl text-brand-orange leading-none ">
            kanad
          </h1>
        </div>
      </aside>

      {/* Overlay for mobile */}
      {isOpen && (
        <div
          onClick={() => setIsOpen(false)}
          className="md:hidden fixed inset-0 bg-black bg-opacity-50 z-30"
        />
      )}
    </>
  );
}
