# Kanad Frontend Quick Start Guide

**Time to first screen**: ~2 minutes

---

## 1. Start the Frontend (Development Mode)

```bash
cd /home/mk/deeprealm/kanad/web
npm run dev
```

**Expected output**:
```
  ▲ Next.js 15.5.4 (Turbopack)
  - Local:        http://localhost:3000
  - Environments: .env.local

 ✓ Starting...
 ✓ Ready in 1.2s
```

**Open**: http://localhost:3000

---

## 2. What You'll See

### Home Page (/)
- Split-screen design
- Left: Welcome text + email input + "Go" button
- Right: Artistic image with "kanad" logo (orange, Bietro font)
- Click "Go" → Goes to dashboard

### Dashboard (/dashboard)
- Black sidebar on left:
  - User profile
  - Navigation (Docs, Tutorials)
  - Kanad logo at bottom
- White main area:
  - "click anywhere to start" message
  - Quick action buttons
- Settings gear icon top-right

---

## 3. Current Status

### ✅ Working
- Home page (pixel-perfect match to design)
- Dashboard layout (sidebar + header)
- Fonts (Bietro DEMO, Quando)
- Brand colors (#FF8C00 orange)
- Responsive design (mobile hamburger menu)
- TypeScript compilation
- Build system (Next.js 15 + Turbopack)

### ⏳ Not Yet Implemented
- Molecule Builder
- Simulation Config Wizard
- Job Monitor
- Results Viewer
- Settings Modal
- Authentication (login/register)

---

## 4. File Locations

### Pages
- Home: `/home/mk/deeprealm/kanad/web/src/app/page.tsx`
- Dashboard: `/home/mk/deeprealm/kanad/web/src/app/dashboard/page.tsx`

### Components
- Sidebar: `/home/mk/deeprealm/kanad/web/src/components/layout/Sidebar.tsx`
- Header: `/home/mk/deeprealm/kanad/web/src/components/layout/Header.tsx`

### API Client
- `/home/mk/deeprealm/kanad/web/src/lib/api/client.ts`

### State Stores (Zustand)
- `/home/mk/deeprealm/kanad/web/src/store/authStore.ts`
- `/home/mk/deeprealm/kanad/web/src/store/moleculeStore.ts`
- `/home/mk/deeprealm/kanad/web/src/store/jobStore.ts`

### Types
- `/home/mk/deeprealm/kanad/web/src/types/api.ts`

---

## 5. Connect to Backend API

### Start Backend (Separate Terminal)
```bash
cd /home/mk/deeprealm/kanad/kanad-backend
source ../env/bin/activate
export DATABASE_URL="sqlite:///./kanad_dev.db"
uvicorn api.main:app --reload --port 8000
```

### Configure Frontend
Create `/home/mk/deeprealm/kanad/web/.env.local`:
```env
NEXT_PUBLIC_API_URL=http://localhost:8000/api
NEXT_PUBLIC_WS_URL=ws://localhost:8000/api
```

---

## 6. Build for Production

```bash
cd /home/mk/deeprealm/kanad/web
npm run build
npm run start
```

Production server runs on: http://localhost:3000

---

## 7. Next Development Steps

### Phase 2: Build Core Features

**1. Molecule Builder** (2-3 hours)
- Create `/src/components/molecule/MoleculeBuilder.tsx`
- Tabs for: Atoms, SMILES, Library, Upload
- Connect to `POST /api/molecules/create`

**2. Simulation Config Wizard** (3-4 hours)
- Create `/src/components/simulation/ConfigWizard.tsx`
- Multi-step form (Method → VQE Config → Backend → Analysis → Submit)
- Connect to `POST /api/simulations/configure` → `POST /api/simulations/{id}/submit`

**3. Job Monitor** (3-4 hours)
- Create `/src/components/job/JobMonitor.tsx`
- Real-time progress, logs, convergence chart
- Connect to `GET /api/jobs/{id}/status` + WebSocket for logs

**4. Results Viewer** (4-5 hours)
- Create `/src/components/job/ResultsViewer.tsx`
- Energy cards, charts, AI report, export buttons
- Connect to `GET /api/jobs/{id}/results`

**5. Settings Modal** (2-3 hours)
- Create `/src/components/settings/SettingsModal.tsx`
- Cloud credentials, defaults, user profile
- Connect to `POST /api/cloud/credentials`, `GET/PUT /api/settings/defaults`

---

## 8. Design Reference Images

**View design mockups**:
- Home: `/home/mk/deeprealm/kanad/web/public/home.jpg`
- Dashboard: `/home/mk/deeprealm/kanad/web/public/dashboard.jpg`
- Hero image: `/home/mk/deeprealm/kanad/web/public/image.webp`

---

## 9. Key Architecture Decisions

### State Management: Zustand
```typescript
import { useAuthStore } from "@/store/authStore";

const { user, login, logout } = useAuthStore();
```

### API Calls
```typescript
import { api } from "@/lib/api/client";

const molecule = await api.molecules.create({ ... });
```

### Styling
- Tailwind CSS 4 for utility classes
- shadcn/ui components (install as needed)
- Custom fonts: Bietro DEMO (logo), Quando (body)

---

## 10. Troubleshooting

### Port 3000 already in use
```bash
# Kill existing process
lsof -ti:3000 | xargs kill -9

# Or use different port
npm run dev -- -p 3001
```

### Fonts not loading
- Check `/public/fonts/Bietro DEMO.otf` exists
- Check browser console for font errors
- Clear browser cache

### API connection errors
- Ensure backend is running on port 8000
- Check `.env.local` has correct `NEXT_PUBLIC_API_URL`
- Check CORS settings in backend

---

## 11. Documentation

**Full Documentation**:
- Architecture: `/home/mk/deeprealm/kanad/FRONTEND_ARCHITECTURE_GUIDE.md`
- Progress: `/home/mk/deeprealm/kanad/FRONTEND_IMPLEMENTATION_PROGRESS.md`
- Backend API: `/home/mk/deeprealm/kanad/API_BUILD_PLAN.md`
- Backend Status: `/home/mk/deeprealm/kanad/BACKEND_COMPLETE_SUMMARY.md`

---

## Success! 🎉

You now have a working Next.js frontend with:
- Beautiful home page matching design
- Functional dashboard layout
- Complete API client
- State management ready
- Type-safe TypeScript
- Production-ready build system

**Next**: Build the Molecule Builder component and start creating molecules!

---

**Questions?** Check the documentation or inspect the code - everything is well-commented and typed!
