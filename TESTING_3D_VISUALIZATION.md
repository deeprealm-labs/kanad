# Testing 3D Molecule Visualization

## Quick Start

### 1. Start Development Server
```bash
cd web
npm run dev
```

Navigate to: `http://localhost:3000/dashboard`

## Test Scenarios

### Test 1: Manual Atom Placement ✓
**Steps:**
1. Click "New Experiment" button
2. Drag Hydrogen (H) from periodic table to build zone
3. Drag Oxygen (O) to build zone
4. Drag another Hydrogen (H)
5. Set coordinates:
   - O: x=0, y=0, z=0
   - H: x=0, y=0, z=1.159
   - H: x=0, y=0, z=-1.159
6. **Expected**: See water molecule in 3D preview (right panel)
7. Rotate, zoom, pan the 3D view
8. Try different view styles (Stick, Sphere, Cartoon, Surface)

**Success Criteria:**
- ✓ Atoms appear in 3D viewer
- ✓ Can rotate molecule with mouse drag
- ✓ Coordinates update molecule position
- ✓ View style changes work

---

### Test 2: SMILES Input ✓
**Steps:**
1. Click "New Experiment"
2. In center panel, find "Or enter SMILES" input
3. Test these SMILES strings:

   **Water (H2O):**
   ```
   O
   ```

   **Ethanol (C2H5OH):**
   ```
   CCO
   ```

   **Benzene (C6H6):**
   ```
   c1ccccc1
   ```

   **Methane (CH4):**
   ```
   C
   ```

4. Click "Add" button
5. **Expected**:
   - Loading spinner appears
   - After 1-2 seconds, 3D structure appears in right panel
   - Can interact with 3D molecule

**Success Criteria:**
- ✓ SMILES converts to 3D structure
- ✓ Loading state shows
- ✓ Error handling for invalid SMILES
- ✓ Molecule is properly centered

---

### Test 3: XYZ File Upload ✓
**Steps:**
1. Click "New Experiment"
2. Find "Upload XYZ File" button in configuration panel
3. Click and select `water.xyz` from project root
4. **Expected**:
   - "Processing file..." message appears
   - 3D water molecule renders in preview
   - File name shows "✓ water.xyz"

**Test with this XYZ content:**
```xyz
3
Water molecule
O 0.000000 0.000000 0.000000
H 0.000000 0.000000 1.159076
H 0.000000 0.000000 -1.159076
```

**Success Criteria:**
- ✓ File uploads successfully
- ✓ Coordinates parsed correctly
- ✓ 3D structure renders
- ✓ Can still adjust configuration

---

### Test 4: Molecule Library ✓
**Steps:**
1. Click "New Experiment"
2. In left panel, click "Molecules" tab
3. Click on any pre-configured molecule (e.g., "Water", "Ethanol")
4. **Expected**:
   - SMILES auto-fills
   - 3D structure appears
   - Configuration pre-fills

**Success Criteria:**
- ✓ Library molecules load
- ✓ 3D preview works
- ✓ Configuration updates

---

### Test 5: Real-time Updates ✓
**Steps:**
1. Click "New Experiment"
2. Drag 2 Hydrogen atoms to build zone
3. Set coordinates:
   - H1: x=0, y=0, z=0
   - H2: x=0.74, y=0, z=0
4. **Watch 3D preview update in real-time**
5. Adjust H2 x-coordinate to 1.5
6. **Expected**: Molecule stretches in 3D view immediately

**Success Criteria:**
- ✓ Updates happen instantly
- ✓ No lag or delay
- ✓ Smooth transitions

---

### Test 6: Interactive Controls ✓
**Steps:**
1. Create any molecule (atoms, SMILES, or XYZ)
2. Test all controls in 3D viewer header:
   - **View Style dropdown**: Switch between Stick/Sphere/Cartoon/Surface
   - **Reset button** (↻): Resets camera position
   - **Zoom Out** (−): Zooms out
   - **Zoom In** (+): Zooms in
   - **Fullscreen** (⛶): Opens fullscreen view
3. **Mouse interactions**:
   - Click & drag to rotate
   - Scroll to zoom
   - Right-click & drag to pan

**Success Criteria:**
- ✓ All buttons work
- ✓ Mouse controls respond
- ✓ Fullscreen mode works
- ✓ View styles change appearance

---

### Test 7: Complete Workflow ✓
**Steps:**
1. Start at Dashboard Home
2. Click "New Experiment" → **MoleculeCreator** opens
3. Input molecule (any method)
4. See 3D preview update
5. Configure settings (basis set, charge, multiplicity)
6. Click "Create Molecule →" → **MoleculeReview** opens
7. **Expected**: Large 3D viewer with full details
8. Verify:
   - 3D structure is prominent
   - Can rotate and inspect
   - Molecule details show on right
   - "Execute Experiment" button ready

**Success Criteria:**
- ✓ Smooth workflow transition
- ✓ 3D viewer on both screens
- ✓ No Lewis structure anywhere
- ✓ Beautiful, cohesive design

---

## Visual Quality Checks

### Design System Integration ✓
- [ ] Orange brand color used for accents
- [ ] Gradient backgrounds present
- [ ] Rounded corners (rounded-xl)
- [ ] Shadow effects visible
- [ ] Smooth transitions on hover
- [ ] Font: Quando used correctly
- [ ] Spacing consistent with rest of app

### 3D Visualization Quality ✓
- [ ] Atoms colored correctly (CPK colors)
- [ ] Bonds rendered between atoms
- [ ] Smooth rotation/zoom
- [ ] Anti-aliasing enabled (smooth edges)
- [ ] Good lighting
- [ ] Canvas renders properly (no blank boxes)

---

## Error Handling Tests

### Test Invalid SMILES ✓
**Input:** `INVALID123`
**Expected:** Error message with retry option

### Test Empty State ✓
**Steps:** Open MoleculeCreator without adding anything
**Expected:** "No molecule data" message with icon

### Test Network Failure ✓
**Steps:**
1. Disconnect internet
2. Try SMILES input
**Expected:** Error message about conversion failure

### Test Corrupt XYZ File ✓
**Create file with:**
```
not a valid xyz file
```
**Expected:** Parse error message

---

## Performance Checks

### Load Time ✓
- [ ] 3Dmol.js loads within 2 seconds
- [ ] Initial viewer ready quickly
- [ ] No visible lag on interaction

### Responsiveness ✓
- [ ] Coordinate changes update <100ms
- [ ] SMILES conversion <3s
- [ ] File upload processes <500ms
- [ ] Smooth 60fps rotation

---

## Browser Compatibility

Test in:
- [ ] Chrome/Edge (primary)
- [ ] Firefox
- [ ] Safari (may have CORS issues with external APIs)

---

## Known Issues & Workarounds

### Issue 1: SMILES Conversion Fails
**Cause:** PubChem/CACTUS API down or CORS issue
**Workaround:** Use XYZ file upload or manual atom placement

### Issue 2: 3Dmol.js Loading Slowly
**Cause:** CDN latency
**Workaround:** Wait for library to load, fallback 2D view shows

### Issue 3: Canvas Not Visible
**Cause:** Browser extension blocking canvas
**Workaround:** Disable extensions or try different browser

---

## Success Criteria Summary

### Must Have ✓
- ✓ 3D visualization works for atoms, SMILES, XYZ
- ✓ Real-time updates functional
- ✓ Interactive controls work
- ✓ Beautiful design matches frontend
- ✓ No Lewis structure references
- ✓ Build succeeds with no errors

### Nice to Have ✓
- ✓ Smooth animations
- ✓ Error handling graceful
- ✓ Loading states clear
- ✓ Fullscreen mode
- ✓ Multiple view styles

---

## Reporting Issues

If you find bugs, document:
1. Steps to reproduce
2. Expected behavior
3. Actual behavior
4. Browser/OS
5. Console errors (F12 → Console)

---

## Next Steps

After testing:
1. ✓ Verify all test cases pass
2. ✓ Check design quality
3. ✓ Test on multiple browsers
4. ✓ Document any issues found
5. ✓ Deploy to production when ready

---

## Quick Test Commands

```bash
# Build and check for errors
npm run build

# Start dev server
npm run dev

# Lint check
npm run lint
```

---

## File Locations for Testing

- **Test XYZ file**: `/home/mk/deeprealm/kanad/water.xyz`
- **Test XYZ file 2**: `/home/mk/deeprealm/kanad/web/test_h2.xyz`
- **Component**: `web/src/components/molecule/Molecule3DViewer.tsx`
- **Creator**: `web/src/components/molecule/MoleculeCreator.tsx`
- **Review**: `web/src/components/molecule/MoleculeReview.tsx`
