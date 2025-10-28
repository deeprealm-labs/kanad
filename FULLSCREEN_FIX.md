# Fullscreen Mode Fix

## Problem
Fullscreen mode showed blank white screen because it was trying to reuse the same viewer reference that was already attached to the main container.

## Solution
Changed from custom fullscreen modal to **native browser fullscreen API**.

### Changes Made:

**File:** `web/src/components/molecule/Molecule3DViewer.tsx`

### Before (Broken):
```tsx
// Custom fullscreen modal with duplicate viewer ref
{isFullscreen && (
  <div className="fixed inset-0 z-50 bg-background">
    <div ref={viewerRef} className="w-full h-full" /> // ❌ Ref already used!
  </div>
)}
```

### After (Fixed):
```tsx
const toggleFullscreen = () => {
  if (!viewerRef.current) return;

  if (!isFullscreen) {
    // Use native Fullscreen API
    viewerRef.current.requestFullscreen();
  } else {
    document.exitFullscreen();
  }
};

// Listen for fullscreen changes and resize viewer
useEffect(() => {
  const handleFullscreenChange = () => {
    setIsFullscreen(!!document.fullscreenElement);
    
    // Resize viewer after transition
    setTimeout(() => {
      if (viewerInstanceRef.current) {
        viewerInstanceRef.current.resize();
        viewerInstanceRef.current.render();
      }
    }, 100);
  };

  document.addEventListener('fullscreenchange', handleFullscreenChange);
  return () => document.removeEventListener('fullscreenchange', handleFullscreenChange);
}, []);
```

## Benefits:

1. ✅ Uses native browser fullscreen (better UX)
2. ✅ No duplicate viewer instances
3. ✅ Auto-resizes when entering/exiting fullscreen
4. ✅ Works across all browsers (Chrome, Firefox, Safari)
5. ✅ User can exit with Esc key (browser standard)

## How It Works:

1. Click fullscreen button → calls `viewerRef.current.requestFullscreen()`
2. Browser puts the viewer div in fullscreen mode
3. Fullscreen change event fires → updates state
4. Viewer auto-resizes to fill screen
5. Press Esc or click minimize → exits fullscreen

## Test:

1. Refresh page
2. Load a molecule (SMILES, XYZ, or atoms)
3. Click fullscreen button (⛶)
4. ✅ Should see molecule in fullscreen!
5. Can still rotate, zoom, change view styles
6. Press Esc or click minimize to exit

All fixed! 🎉
