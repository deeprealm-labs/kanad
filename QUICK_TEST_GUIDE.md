# Quick Test Guide - 3D Molecule Visualization

## 🚀 Quick Start

```bash
cd /home/mk/deeprealm/kanad/web
npm run dev
```

Open browser: `http://localhost:3000/dashboard`

---

## ✅ Quick Tests

### Test 1: XYZ File (30 seconds)
1. Click "New Experiment"
2. Upload `water.xyz` file
3. ✅ Should see 3D water molecule instantly
4. Try rotating with mouse
5. Try different view styles (Stick, Sphere, Cartoon)

### Test 2: SMILES (30 seconds)
1. Click "New Experiment"  
2. Type `CCO` in SMILES field
3. Click "Add"
4. ✅ Should see ethanol molecule in 3D (after 1-2s loading)

### Test 3: Manual Atoms (1 minute)
1. Click "New Experiment"
2. Drag H, O, H to build zone
3. Set coords: O(0,0,0), H(0,0,1.16), H(0,0.76,-0.50)
4. ✅ Should see real-time 3D preview update

---

## 🔍 Console Checks

Open Browser Console (F12):

### Should See These Logs:
```
✓ 3Dmol.js loaded successfully
✓ Viewer instance created  
✓ 3Dmol viewer initialized successfully
✓ Canvas resized and rendered: 400 x 300
✓ Successfully loaded molecule from XYZ data
✓ Molecule rendered successfully
```

### Check 3Dmol is Loaded:
Type in console:
```javascript
window.$3Dmol
```
Should return: `{createViewer: ƒ, ...}`

---

## ❌ If Something's Wrong

### Blank Canvas?
- Check console for "❌" errors
- Try refreshing page
- Click "Reset View" button

### "No molecule data"?
- File didn't upload
- Check file is valid .xyz format
- Try manual atoms instead

### SMILES not working?
- Check internet connection
- Try XYZ file instead
- May be external API issue

---

## 📊 What Should Work

- [x] XYZ file upload → 3D display
- [x] SMILES input → 3D display  
- [x] Manual atoms → 3D display
- [x] Mouse rotation/zoom
- [x] View style changes
- [x] Fullscreen mode
- [x] Real-time coordinate updates

---

## 🎯 Expected UX

1. **Upload XYZ**: Instant 3D molecule
2. **Type SMILES**: 1-2s loading → 3D molecule
3. **Drag atoms**: Real-time 3D preview
4. **Smooth rotation**: 60fps interaction
5. **Beautiful design**: Orange accents, clean layout

---

## 📞 Getting Help

If issues persist:
1. Check console (F12) for errors
2. Review `XYZ_VISUALIZATION_FIXES.md`
3. Try different browser (Chrome recommended)
4. Check `water.xyz` file is not corrupted

---

**Expected Result: Beautiful, interactive 3D molecules! 🎉**
