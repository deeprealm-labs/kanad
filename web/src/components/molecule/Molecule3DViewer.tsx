"use client";

import { useEffect, useRef, useState, useCallback } from "react";
import { RotateCcw, ZoomIn, ZoomOut, Maximize2, Minimize2, Settings, Atom } from "lucide-react";

interface AtomData {
  symbol: string;
  x: number;
  y: number;
  z: number;
  atomicNumber?: number;
}

interface Molecule3DViewerProps {
  molecule?: {
    smiles?: string;
    atoms?: AtomData[];
    xyzData?: string;
  };
  width?: number | string;
  height?: number | string;
  className?: string;
  showControls?: boolean;
  onAtomClick?: (atom: any) => void;
}

// Declare global 3Dmol
declare global {
  interface Window {
    $3Dmol: any;
  }
}

export default function Molecule3DViewer({
  molecule,
  width = "100%",
  height = "400px",
  className = "",
  showControls = true,
  onAtomClick,
}: Molecule3DViewerProps) {
  const viewerRef = useRef<HTMLDivElement>(null);
  const viewerInstanceRef = useRef<any>(null);
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [isFullscreen, setIsFullscreen] = useState(false);
  const [viewStyle, setViewStyle] = useState<"stick" | "sphere" | "cartoon" | "surface">("stick");

  // Load 3Dmol.js dynamically
  useEffect(() => {
    const load3Dmol = async () => {
      if (typeof window !== "undefined" && !window.$3Dmol) {
        console.log("Loading 3Dmol.js library...");

        // Check if script already exists
        const existingScript = document.querySelector('script[src*="3Dmol"]');
        if (existingScript) {
          existingScript.remove();
        }

        const script = document.createElement("script");
        script.src = "https://3dmol.org/build/3Dmol-min.js";
        script.async = true;
        script.onload = () => {
          console.log("✓ 3Dmol.js loaded successfully");
          // Re-initialize viewer after library loads
          setTimeout(() => {
            initializeViewer();
          }, 200);
        };
        script.onerror = () => {
          console.error("❌ Failed to load 3Dmol.js");
          setError("Failed to load 3Dmol.js library");
        };
        document.head.appendChild(script);
      } else if (window.$3Dmol) {
        console.log("✓ 3Dmol.js already loaded");
        // Library already loaded, initialize immediately
        setTimeout(() => {
          initializeViewer();
        }, 100);
      }
    };

    load3Dmol();
  }, []);

  // Initialize viewer
  const initializeViewer = useCallback(() => {
    if (!viewerRef.current) {
      console.log("❌ Viewer ref not available");
      return;
    }

    if (!window.$3Dmol) {
      console.log("❌ 3Dmol.js not loaded");
      return;
    }

    try {
      console.log("Initializing 3Dmol viewer...");

      // Clear existing viewer
      if (viewerInstanceRef.current) {
        try {
          viewerInstanceRef.current.clear();
        } catch (e) {
          console.log("Could not clear existing viewer:", e);
        }
      }

      // Create new viewer with proper configuration
      viewerInstanceRef.current = window.$3Dmol.createViewer(viewerRef.current, {
        backgroundColor: "white",
        defaultcolors: window.$3Dmol.rasmolElementColors,
        antialias: true,
        quality: 'medium'
      });

      console.log("✓ Viewer instance created");

      // Force immediate render to ensure canvas is visible
      viewerInstanceRef.current.render();

      // Add click handler for atoms (only if the method exists)
      if (onAtomClick && typeof viewerInstanceRef.current.addClickHandler === 'function') {
        try {
          viewerInstanceRef.current.addClickHandler({
            onAtom: (atom: any) => {
              onAtomClick(atom);
            },
          });
        } catch (err) {
          console.log("Click handler not supported, continuing without it");
        }
      }

      console.log("✓ 3Dmol viewer initialized successfully");

      // Force initial render and resize
      setTimeout(() => {
        if (viewerInstanceRef.current && viewerRef.current) {
          // Ensure canvas is properly sized
          const rect = viewerRef.current.getBoundingClientRect();
          if (rect.width > 0 && rect.height > 0) {
            viewerInstanceRef.current.resize();
            viewerInstanceRef.current.render();
            console.log("✓ Canvas resized and rendered:", rect.width, "x", rect.height);
          } else {
            console.log("❌ Canvas has invalid dimensions:", rect.width, "x", rect.height);
          }
        }
      }, 100);

    } catch (err) {
      console.error("❌ Error initializing 3Dmol viewer:", err);
      setError("Failed to initialize 3D viewer");
    }
  }, [onAtomClick]);

  // Convert SMILES to 3D coordinates using a simple conversion API or fallback
  const convertSMILESto3D = async (smiles: string): Promise<string | null> => {
    console.log(`🔄 Converting SMILES "${smiles}" to 3D structure...`);

    try {
      // Try using PubChem's API to convert SMILES to SDF format
      console.log("Trying PubChem API...");
      const pubchemUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/${encodeURIComponent(smiles)}/SDF`;
      console.log("PubChem URL:", pubchemUrl);

      const response = await fetch(pubchemUrl, {
        method: 'GET',
        mode: 'cors'
      });

      console.log("PubChem response status:", response.status);

      if (response.ok) {
        const sdf = await response.text();
        console.log("✓ Converted SMILES to SDF using PubChem, length:", sdf.length);
        console.log("SDF preview:", sdf.substring(0, 200));
        return sdf;
      } else {
        console.log("PubChem returned non-OK status:", response.status, response.statusText);
      }
    } catch (err) {
      console.log("❌ PubChem conversion failed:", err);
    }

    // Fallback: Use NIH/NCI CACTUS service
    try {
      console.log("Trying CACTUS API...");
      const cactusUrl = `https://cactus.nci.nih.gov/chemical/structure/${encodeURIComponent(smiles)}/sdf`;
      console.log("CACTUS URL:", cactusUrl);

      const response = await fetch(cactusUrl, {
        method: 'GET',
        mode: 'cors'
      });

      console.log("CACTUS response status:", response.status);

      if (response.ok) {
        const sdf = await response.text();
        console.log("✓ Converted SMILES to SDF using CACTUS, length:", sdf.length);
        console.log("SDF preview:", sdf.substring(0, 200));
        return sdf;
      } else {
        console.log("CACTUS returned non-OK status:", response.status, response.statusText);
      }
    } catch (err) {
      console.log("❌ CACTUS conversion failed:", err);
    }

    console.log("❌ All SMILES conversion methods failed");
    return null;
  };

  // Load molecule data
  const loadMolecule = useCallback(async () => {
    if (!viewerInstanceRef.current) {
      console.log("No viewer instance available");
      return;
    }

    console.log("=== LOADING MOLECULE ===");
    console.log("Molecule data:", molecule);
    setIsLoading(true);
    setError(null);

    try {
      // Clear previous molecule
      viewerInstanceRef.current.clear();

      let modelAdded = false;

      // Try atoms array first (for live updates)
      if (molecule?.atoms && molecule.atoms.length > 0) {
        console.log("Processing atoms array:", molecule.atoms);
        try {
          // Convert atoms to XYZ format
          const xyzContent = `${molecule.atoms.length}\nGenerated from atoms\n${molecule.atoms
            .map((atom) => `${atom.symbol} ${atom.x.toFixed(6)} ${atom.y.toFixed(6)} ${atom.z.toFixed(6)}`)
            .join("\n")}`;

          console.log("Adding atoms model (XYZ format):", xyzContent);
          viewerInstanceRef.current.addModel(xyzContent, "xyz");
          modelAdded = true;
          console.log("✓ Successfully loaded molecule from atoms array:", molecule.atoms.length, "atoms");
        } catch (err) {
          console.error("Failed to load atoms array:", err);
        }
      }

      // Try XYZ data next
      if (!modelAdded && molecule?.xyzData && molecule.xyzData.trim()) {
        console.log("Processing XYZ data:", molecule.xyzData.substring(0, 200) + "...");
        try {
          console.log("Adding XYZ model");
          viewerInstanceRef.current.addModel(molecule.xyzData, "xyz");
          modelAdded = true;
          console.log("✓ Successfully loaded molecule from XYZ data");
        } catch (err) {
          console.error("Failed to load XYZ data:", err);
        }
      }

      // Try SMILES if no atoms or XYZ (convert to SDF first)
      if (!modelAdded && molecule?.smiles && molecule.smiles.trim()) {
        console.log("Processing SMILES:", molecule.smiles);
        try {
          // Convert SMILES to SDF format for 3Dmol.js
          const sdf = await convertSMILESto3D(molecule.smiles);

          if (sdf) {
            viewerInstanceRef.current.addModel(sdf, "sdf");
            modelAdded = true;
            console.log("✓ Successfully loaded molecule from SMILES (converted to SDF)");
          } else {
            console.log("Failed to convert SMILES to 3D structure");
            setError("Unable to generate 3D structure from SMILES. Try using an XYZ file or manually placing atoms.");
          }
        } catch (err) {
          console.error("Failed to load SMILES:", err);
          setError("Failed to convert SMILES to 3D structure");
        }
      }

      if (modelAdded) {
        console.log("Applying view style:", viewStyle);
        // Apply style
        applyViewStyle();

        // Zoom to fit
        viewerInstanceRef.current.zoomTo();

        // Force render with proper canvas update
        viewerInstanceRef.current.render();

        // Ensure canvas is visible and properly sized
        setTimeout(() => {
          if (viewerInstanceRef.current && viewerRef.current) {
            const canvas = viewerRef.current.querySelector('canvas');
            if (canvas) {
              console.log("Canvas found, dimensions:", canvas.width, "x", canvas.height);
              canvas.style.display = 'block';
              canvas.style.visibility = 'visible';
            }
            viewerInstanceRef.current.resize();
            viewerInstanceRef.current.render();
          }
        }, 50);

        setTimeout(() => {
          if (viewerInstanceRef.current) {
            viewerInstanceRef.current.render();
            console.log("Final render completed");
          }
        }, 200);

        console.log("✓ Molecule rendered successfully");
      } else {
        console.log("❌ No molecule data to display");
        setError("No valid molecule data provided");
      }
    } catch (err) {
      console.error("Error loading molecule:", err);
      setError("Failed to load molecule: " + (err as Error).message);
    } finally {
      setIsLoading(false);
    }
  }, [molecule, viewStyle]);

  // Apply visualization style
  const applyViewStyle = useCallback(() => {
    if (!viewerInstanceRef.current) {
      console.log("No viewer instance for applying style");
      return;
    }

    console.log("Applying view style:", viewStyle);
    
    try {
      const style = {};
      switch (viewStyle) {
        case "stick":
          Object.assign(style, { stick: { radius: 0.1 } });
          break;
        case "sphere":
          Object.assign(style, { sphere: { scale: 0.3 } });
          break;
        case "cartoon":
          Object.assign(style, { cartoon: {} });
          break;
        case "surface":
          Object.assign(style, { surface: { opacity: 0.7 } });
          break;
      }

      console.log("Style object:", style);
      viewerInstanceRef.current.setStyle({}, style);
      viewerInstanceRef.current.render();
      console.log("✓ Style applied successfully");
    } catch (err) {
      console.error("Error applying style:", err);
    }
  }, [viewStyle]);

  // Initialize viewer when 3Dmol is ready
  useEffect(() => {
    if (window.$3Dmol) {
      initializeViewer();
    }
  }, [initializeViewer]);

  // Load molecule when data changes
  useEffect(() => {
    if (viewerInstanceRef.current && molecule) {
      console.log("Molecule data changed, loading molecule...");
      loadMolecule();
    } else if (molecule && !viewerInstanceRef.current) {
      console.log("Molecule data available but viewer not ready, will retry...");
      // Retry after a delay if viewer isn't ready yet
      const retryTimer = setTimeout(() => {
        if (viewerInstanceRef.current && molecule) {
          console.log("Retrying molecule load...");
          loadMolecule();
        }
      }, 500);
      return () => clearTimeout(retryTimer);
    }
  }, [molecule, loadMolecule]);

  // Apply style when viewStyle changes
  useEffect(() => {
    if (viewerInstanceRef.current) {
      applyViewStyle();
    }
  }, [viewStyle, applyViewStyle]);

  // Handle container resize
  useEffect(() => {
    if (!viewerRef.current || !viewerInstanceRef.current) return;

    const resizeObserver = new ResizeObserver(() => {
      if (viewerInstanceRef.current) {
        setTimeout(() => {
          viewerInstanceRef.current?.resize();
          viewerInstanceRef.current?.render();
        }, 100);
      }
    });

    resizeObserver.observe(viewerRef.current);

    return () => {
      resizeObserver.disconnect();
    };
  }, []);

  const handleReset = () => {
    if (viewerInstanceRef.current) {
      viewerInstanceRef.current.zoomTo();
      viewerInstanceRef.current.render();
    }
  };

  const handleZoomIn = () => {
    if (viewerInstanceRef.current) {
      viewerInstanceRef.current.zoom(1.2);
      viewerInstanceRef.current.render();
    }
  };

  const handleZoomOut = () => {
    if (viewerInstanceRef.current) {
      viewerInstanceRef.current.zoom(0.8);
      viewerInstanceRef.current.render();
    }
  };

  const toggleFullscreen = () => {
    if (!viewerRef.current) return;

    if (!isFullscreen) {
      // Enter fullscreen
      if (viewerRef.current.requestFullscreen) {
        viewerRef.current.requestFullscreen();
      } else if ((viewerRef.current as any).webkitRequestFullscreen) {
        (viewerRef.current as any).webkitRequestFullscreen();
      } else if ((viewerRef.current as any).msRequestFullscreen) {
        (viewerRef.current as any).msRequestFullscreen();
      }
    } else {
      // Exit fullscreen
      if (document.exitFullscreen) {
        document.exitFullscreen();
      } else if ((document as any).webkitExitFullscreen) {
        (document as any).webkitExitFullscreen();
      } else if ((document as any).msExitFullscreen) {
        (document as any).msExitFullscreen();
      }
    }
  };

  // Handle fullscreen change events
  useEffect(() => {
    const handleFullscreenChange = () => {
      setIsFullscreen(!!document.fullscreenElement);

      // Resize viewer when entering/exiting fullscreen
      setTimeout(() => {
        if (viewerInstanceRef.current) {
          viewerInstanceRef.current.resize();
          viewerInstanceRef.current.render();
        }
      }, 100);
    };

    document.addEventListener('fullscreenchange', handleFullscreenChange);
    document.addEventListener('webkitfullscreenchange', handleFullscreenChange);
    document.addEventListener('msfullscreenchange', handleFullscreenChange);

    return () => {
      document.removeEventListener('fullscreenchange', handleFullscreenChange);
      document.removeEventListener('webkitfullscreenchange', handleFullscreenChange);
      document.removeEventListener('msfullscreenchange', handleFullscreenChange);
    };
  }, []);

  const containerStyle = {
    width: typeof width === "number" ? `${width}px` : width,
    height: typeof height === "number" ? `${height}px` : height,
  };

  return (
    <div className={`relative bg-card border border-border rounded-lg overflow-hidden ${className}`}>
      {/* Header */}
      <div className="flex items-center justify-between px-3 py-2 border-b border-border bg-muted/30">
        <div className="flex items-center gap-2">
          <Atom className="w-4 h-4 text-brand-orange" />
          <span className="text-sm font-quando font-medium">
            3D Molecule Viewer
          </span>
        </div>

        {showControls && (
          <div className="flex items-center gap-1">
            {/* Style Controls */}
            <select
              value={viewStyle}
              onChange={(e) => setViewStyle(e.target.value as any)}
              className="px-2 py-1 text-xs border border-input bg-background rounded-md font-quando hover:border-brand-orange/50 focus:border-brand-orange focus:outline-none transition-all"
            >
              <option value="stick">Stick</option>
              <option value="sphere">Sphere</option>
              <option value="cartoon">Cartoon</option>
              <option value="surface">Surface</option>
            </select>

            {/* View Controls */}
            <button
              onClick={handleReset}
              className="p-1 hover:bg-brand-orange/10 hover:text-brand-orange rounded transition-all"
              title="Reset View"
            >
              <RotateCcw className="w-3.5 h-3.5" />
            </button>
            <button
              onClick={handleZoomOut}
              className="p-1 hover:bg-brand-orange/10 hover:text-brand-orange rounded transition-all"
              title="Zoom Out"
            >
              <ZoomOut className="w-3.5 h-3.5" />
            </button>
            <button
              onClick={handleZoomIn}
              className="p-1 hover:bg-brand-orange/10 hover:text-brand-orange rounded transition-all"
              title="Zoom In"
            >
              <ZoomIn className="w-3.5 h-3.5" />
            </button>
            <button
              onClick={toggleFullscreen}
              className="p-1 hover:bg-brand-orange/10 hover:text-brand-orange rounded transition-all"
              title={isFullscreen ? "Exit Fullscreen" : "Fullscreen"}
            >
              {isFullscreen ? <Minimize2 className="w-3.5 h-3.5" /> : <Maximize2 className="w-3.5 h-3.5" />}
            </button>
          </div>
        )}
      </div>

      {/* Viewer Container */}
      <div className="relative bg-background">
        <div
          ref={viewerRef}
          style={containerStyle}
          className="w-full"
        />
        
        
        {/* Loading Overlay */}
        {isLoading && (
          <div className="absolute inset-0 bg-background/80 flex items-center justify-center">
            <div className="flex items-center gap-2 text-sm font-quando">
              <div className="w-4 h-4 border-2 border-brand-orange border-t-transparent rounded-full animate-spin" />
              Loading molecule...
            </div>
          </div>
        )}

        {/* Error Overlay */}
        {error && (
          <div className="absolute inset-0 bg-background/80 flex items-center justify-center">
            <div className="text-center p-4">
              <div className="text-sm font-quando text-destructive mb-2">
                {error}
              </div>
              <button
                onClick={() => {
                  setError(null);
                  if (molecule) loadMolecule();
                }}
                className="px-3 py-1 text-xs bg-brand-orange text-white rounded hover:bg-brand-orange-dark transition font-quando"
              >
                Retry
              </button>
            </div>
          </div>
        )}

        {/* Empty State */}
        {!molecule && !isLoading && !error && (
          <div className="absolute inset-0 flex items-center justify-center">
            <div className="text-center text-muted-foreground font-quando">
              <Atom className="w-8 h-8 mx-auto mb-2 opacity-50" />
              <div className="text-sm">No molecule data</div>
              <div className="text-xs">Add SMILES, upload XYZ, or build with atoms</div>
            </div>
          </div>
        )}

        {/* Debug Info */}
        {process.env.NODE_ENV === 'development' && (
          <div className="absolute top-2 left-2 bg-black/80 text-white text-xs p-2 rounded font-mono">
            <div>3Dmol: {window.$3Dmol ? '✓' : '✗'}</div>
            <div>Viewer: {viewerInstanceRef.current ? '✓' : '✗'}</div>
            <div>Atoms: {molecule?.atoms?.length || (molecule?.xyzData ? molecule.xyzData.trim().split('\n')[0] : '0')}</div>
            <div>SMILES: {molecule?.smiles || 'none'}</div>
            <div>XYZ: {molecule?.xyzData ? '✓' : '✗'}</div>
            <div>Style: {viewStyle}</div>
          </div>
        )}
      </div>

    </div>
  );
}
