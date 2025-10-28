# 3D Visualization Debugging Guide

## 🔍 How to Debug the 3D Molecule Viewer

### Step 1: Open Browser Console
1. Press **F12** or right-click → Inspect
2. Go to **Console** tab
3. Clear console (🗑️ icon)

### Step 2: Load a Molecule

#### For SMILES (e.g., "CCO"):
You should see these logs in order:

```
1. User inputs SMILES "CCO"

2. === MOLECULE CREATOR DEBUG ===
   Dropped atoms: []
   SMILES: CCO
   XYZ data:
   Molecule data being passed to 3D viewer: {atoms: [], smiles: "CCO", xyzData: ""}

3. Molecule data changed, loading molecule...

4. === LOADING MOLECULE ===
   Molecule data: {atoms: [], smiles: "CCO", xyzData: ""}

5. Processing SMILES: CCO
   🔄 Converting SMILES "CCO" to 3D structure...
   Trying PubChem API...
   PubChem URL: https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/CCO/SDF
   PubChem response status: 200
   ✓ Converted SMILES to SDF using PubChem, length: 1234
   SDF preview: C2H6O...

6. ✓ Successfully loaded molecule from SMILES (converted to SDF)

7. Applying view style: stick
   Style object: {stick: {radius: 0.1}}
   ✓ Style applied successfully

8. ✓ Molecule rendered successfully
```

**If you see all ✓ symbols → Everything is working!**

---

#### For XYZ File:
You should see:

```
1. XYZ file content: 3
   Water molecule
   O 0.000000 0.000000 0.000000
   ...

2. Parsed atoms: [{symbol: "O", x: 0, y: 0, z: 0, atomicNumber: 8}, ...]

3. === MOLECULE CREATOR DEBUG ===
   XYZ data: 3
   Water molecule
   ...

4. Molecule data changed, loading molecule...

5. === LOADING MOLECULE ===

6. Processing XYZ data: 3...
   Adding XYZ model

7. ✓ Successfully loaded molecule from XYZ data

8. Applying view style: sphere

9. ✓ Molecule rendered successfully
```

---

### Common Error Patterns

#### ❌ Error 1: "Maximum update depth exceeded"
**Cause:** Infinite loop (FIXED in latest version)
**Solution:** Update code to remove `error` from dependency array

#### ❌ Error 2: "No viewer instance available"
**Logs show:**
```
❌ Viewer ref not available
```
or
```
❌ 3Dmol.js not loaded
```

**Solution:**
- Wait for library to load (takes 1-2 seconds)
- Check internet connection
- Refresh page

#### ❌ Error 3: "Failed to convert SMILES"
**Logs show:**
```
❌ PubChem conversion failed: TypeError: Failed to fetch
❌ CACTUS conversion failed: TypeError: Failed to fetch
❌ All SMILES conversion methods failed
```

**Causes:**
- No internet connection
- CORS error (browser blocking request)
- External API down

**Solutions:**
- Check internet connection
- Use XYZ file instead
- Use manual atom placement
- Wait and retry

#### ❌ Error 4: Blank canvas (no errors in console)
**Logs show:**
```
✓ 3Dmol.js loaded successfully
✓ Viewer instance created
✓ Successfully loaded molecule
✓ Molecule rendered successfully
```

But canvas is blank.

**Causes:**
- Canvas has invalid dimensions
- 3Dmol.js rendering issue
- Browser incompatibility

**Solutions:**
- Check console for: "Canvas resized and rendered: X x Y"
- If dimensions are 0x0 → CSS issue
- Try different browser
- Click "Reset View" button

---

### Debug Panel (Development Mode)

In the top-left corner of 3D viewer, you'll see:

```
3Dmol: ✓     ← Library loaded
Viewer: ✓    ← Viewer initialized
Atoms: 3     ← Number of atoms
SMILES: CCO  ← SMILES string (or "none")
XYZ: ✓       ← XYZ data present
Style: stick ← Current view style
```

**All should show ✓ when working!**

---

### Step-by-Step SMILES Debugging

1. **Enter SMILES:** Type "CCO" and click Add
2. **Check Console:**
   ```
   🔄 Converting SMILES "CCO" to 3D structure...
   ```

3. **Check API Call:**
   ```
   Trying PubChem API...
   PubChem URL: https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/CCO/SDF
   ```

4. **Check Response:**
   ```
   PubChem response status: 200
   ```
   - **200** = Success ✅
   - **404** = Not found ❌
   - **503** = Service unavailable ❌

5. **Check Conversion:**
   ```
   ✓ Converted SMILES to SDF using PubChem, length: 1234
   SDF preview: C2H6O...
   ```

6. **Check Loading:**
   ```
   ✓ Successfully loaded molecule from SMILES (converted to SDF)
   ```

7. **Check Rendering:**
   ```
   ✓ Molecule rendered successfully
   ```

**If ANY step fails, check the error message and follow solutions above.**

---

### Network Tab Debugging

1. Open DevTools → **Network** tab
2. Upload XYZ or enter SMILES
3. Look for requests to:
   - `pubchem.ncbi.nlm.nih.gov`
   - `cactus.nci.nih.gov`

4. Click on request to see:
   - **Status:** Should be 200
   - **Preview:** Should show SDF data
   - **Response:** Full SDF content

**If status is not 200:**
- Check CORS headers
- Try different browser
- Check if API is blocked by firewall/proxy

---

### Testing Checklist

- [ ] 3Dmol.js loads (check console for "✓ 3Dmol.js loaded")
- [ ] Viewer initializes (check console for "✓ Viewer instance created")
- [ ] XYZ file uploads successfully
- [ ] XYZ data parses correctly
- [ ] Atoms display in 3D viewer
- [ ] SMILES converts to SDF
- [ ] SMILES molecule displays
- [ ] Manual atoms work
- [ ] Can rotate molecule
- [ ] Can zoom in/out
- [ ] View styles work
- [ ] No console errors
- [ ] No infinite loops

---

### Quick Fixes

#### Fix 1: Force Reload Molecule
In development mode, click the 🔄 button in viewer header.

#### Fix 2: Refresh Page
Sometimes viewer gets stuck, full page refresh helps.

#### Fix 3: Clear Browser Cache
Shift + F5 or Ctrl + Shift + R

#### Fix 4: Try Different Input Method
- SMILES not working? → Try XYZ file
- XYZ not working? → Try manual atoms
- Manual atoms not working? → Try SMILES

---

### Success Indicators

✅ **All working when you see:**
1. No ❌ symbols in console
2. Debug panel shows all ✓
3. Molecule visible in canvas
4. Can rotate with mouse
5. Can change view styles
6. No error overlays

---

### Getting More Help

1. **Check all console logs** (F12)
2. **Screenshot error messages**
3. **Note what you tried:**
   - Input method (SMILES/XYZ/atoms)
   - Input value (what SMILES? which XYZ file?)
   - Browser and version
4. **Review documentation:**
   - `XYZ_VISUALIZATION_FIXES.md`
   - `INFINITE_LOOP_FIX.md`
   - `TESTING_3D_VISUALIZATION.md`

---

**Happy Debugging! 🔍**
