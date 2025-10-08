# Kanad Frontend - API Integration Complete

**Date**: October 8, 2025
**Status**: Core API Integration Layer Complete
**Backend API URL**: http://localhost:8000

---

## Overview

The Next.js frontend has been successfully updated to connect with the FastAPI backend. All core infrastructure is in place, with graceful fallback to localStorage when the API is unavailable.

---

## What Has Been Completed

### 1. Core API Infrastructure

#### `/web/src/lib/types.ts`
- **Complete TypeScript type definitions** for all API requests and responses
- Types for: Experiments, Queue, Settings, Molecules, Backend Configuration, WebSocket messages
- Fully typed interfaces matching backend API schema

#### `/web/src/lib/api.ts`
- **Comprehensive API client library** with all functions
- Features:
  - Automatic retry logic (3 retries with exponential backoff)
  - Request timeout handling (30s default)
  - Error handling with custom ApiError class
  - JWT authentication support (token from localStorage)
  - WebSocket support for real-time updates
  - Polling fallback when WebSocket unavailable
  - Export/download functionality
  - Migration helper for localStorage to API

**API Functions Implemented:**
```typescript
// Experiments
- createExperiment(request)
- getExperiments(params)
- getExperiment(experimentId)
- updateExperiment(experimentId, request)
- deleteExperiment(experimentId)
- exportExperiment(experimentId, format)
- downloadExperiment(experimentId, format, filename)

// Queue
- getQueue()
- updateQueueItem(jobId, request)
- deleteQueueItem(jobId)
- reorderQueue(jobIds)

// Settings
- getSettings()
- updateSettings(request)

// Validation
- validateSmiles(smiles)

// Real-time
- createWebSocket(experimentId)
- pollExperimentStatus(experimentId, onUpdate, interval)

// Health
- healthCheck()

// Migration
- migrateLocalStorageToAPI()
```

### 2. UI Components

#### `/web/src/components/ui/toast.tsx`
- **Toast notification system** for user feedback
- Support for success, error, warning, and info messages
- Auto-dismiss with configurable duration
- Beautiful animations and theming
- Context-based API via `useToast()` hook

#### `/web/src/components/ui/loading.tsx`
- **Loading components** for better UX
- Components:
  - `FullPageLoader` - Full-screen loading overlay
  - `Spinner` - Configurable spinner (sm/default/lg)
  - `LoadingButton` - Button with loading state
  - `Skeleton` - Skeleton loader for content
  - `CardSkeleton` - Card-shaped skeleton
  - `TableSkeleton` - Table row skeletons
  - `LoadingOverlay` - Overlay for specific sections
  - `ProgressBar` - Progress bar with percentage

### 3. Updated Pages & Components

#### `/web/src/app/layout.tsx`
- Added `ToastProvider` to root layout
- All pages now have access to toast notifications

#### `/web/src/app/dashboard/page.tsx`
- **Fully integrated with API**
- Features:
  - Load experiments from API with localStorage fallback
  - Load settings from API with localStorage fallback
  - Create experiments via API (executeNow: true/false)
  - Add experiments to queue via API
  - Real-time experiment monitoring
  - Error handling with toast notifications
  - Loading states with FullPageLoader

#### `/web/src/components/simulation/ExperimentMonitor.tsx`
- **Real-time experiment monitoring**
- Features:
  - WebSocket support for live updates
  - Polling fallback (2-second interval)
  - Real-time convergence chart updates
  - Live logs streaming
  - Progress tracking
  - Export functionality via API
  - Simulation mode when offline
  - Toast notifications for status changes

### 4. Environment Configuration

#### `/web/.env.example`
```bash
NEXT_PUBLIC_API_URL=http://localhost:8000
```

#### `/web/.env.local`
```bash
NEXT_PUBLIC_API_URL=http://localhost:8000
```

---

## Backend Endpoints Required

The frontend expects the following backend API endpoints:

### Experiments
- `POST /api/experiments` - Create new experiment
- `GET /api/experiments?status=&limit=&offset=&search=` - List experiments
- `GET /api/experiments/{id}` - Get single experiment
- `PATCH /api/experiments/{id}` - Update experiment
- `DELETE /api/experiments/{id}` - Delete experiment
- `GET /api/experiments/{id}/export?format=json|csv` - Export results

### Queue
- `GET /api/queue` - Get queue
- `PATCH /api/queue/{id}` - Update queue item
- `DELETE /api/queue/{id}` - Delete queue item
- `POST /api/queue/reorder` - Reorder queue

### Settings
- `GET /api/settings` - Get user settings
- `PUT /api/settings` - Update user settings

### Validation
- `POST /api/validate/smiles` - Validate SMILES string

### WebSocket
- `WS /ws/experiments/{id}?token=` - Real-time experiment updates

### Health
- `GET /health` - Health check

---

## Request/Response Schemas

### Create Experiment
**Request:**
```json
{
  "molecule": {
    "smiles": "H2O",
    "basis": "sto-3g",
    "charge": 0,
    "multiplicity": 1
  },
  "backendSettings": {
    "method": "VQE",
    "ansatz": "hardware_efficient",
    "mapper": "jordan_wigner",
    "backend": "classical",
    "optimizer": "SLSQP"
  },
  "analysis": {
    "energyDecomposition": true,
    "bondAnalysis": true,
    "dipoleMoment": true
  },
  "executeNow": true
}
```

**Response:**
```json
{
  "experimentId": "exp_123456",
  "status": "queued",
  "message": "Experiment created successfully"
}
```

### Get Experiments
**Response:**
```json
{
  "experiments": [
    {
      "id": "exp_123456",
      "status": "completed",
      "molecule": { ... },
      "backendSettings": { ... },
      "results": {
        "energy": -1.137283,
        "dipoleMoment": 0.0,
        "converged": true,
        "iterations": 42,
        "convergenceData": [
          { "iteration": 1, "energy": -1.05 },
          { "iteration": 2, "energy": -1.10 }
        ]
      },
      "timestamp": "2025-10-08T10:30:00Z",
      "completedAt": "2025-10-08T10:32:15Z"
    }
  ],
  "total": 1,
  "limit": 10,
  "offset": 0
}
```

### WebSocket Messages
The frontend handles the following WebSocket message types:

```typescript
// Status Update
{
  "type": "status",
  "experimentId": "exp_123",
  "data": {
    "status": "running",
    "progress": 45
  }
}

// Convergence Update
{
  "type": "convergence",
  "experimentId": "exp_123",
  "data": {
    "iteration": 10,
    "energy": -1.125,
    "timestamp": "2025-10-08T10:30:15Z"
  }
}

// Log Update
{
  "type": "log",
  "experimentId": "exp_123",
  "data": {
    "message": "Iteration 10/42 completed",
    "timestamp": "2025-10-08T10:30:15Z",
    "level": "info"
  }
}

// Result Update
{
  "type": "result",
  "experimentId": "exp_123",
  "data": {
    "results": {
      "energy": -1.137283,
      "dipoleMoment": 0.0,
      "converged": true,
      "iterations": 42
    }
  }
}

// Error
{
  "type": "error",
  "experimentId": "exp_123",
  "data": {
    "message": "Circuit execution failed"
  }
}
```

---

## How to Run

### 1. Start Backend (FastAPI)
```bash
cd /home/mk/deeprealm/kanad/api
python -m uvicorn main:app --reload --host 0.0.0.0 --port 8000
```

### 2. Start Frontend (Next.js)
```bash
cd /home/mk/deeprealm/kanad/web
npm run dev
```

The frontend will be available at `http://localhost:3000`

---

## Features Implemented

### ✅ Core Features
- [x] API client library with TypeScript types
- [x] Error handling with retry logic
- [x] Toast notification system
- [x] Loading states and skeletons
- [x] Real-time WebSocket support
- [x] Polling fallback mechanism
- [x] Export/download functionality
- [x] Graceful localStorage fallback
- [x] JWT authentication support (prepared)

### ✅ Dashboard Integration
- [x] Load experiments from API
- [x] Load settings from API
- [x] Create experiments via API
- [x] Add to queue via API
- [x] Real-time experiment monitoring
- [x] Error handling with toasts
- [x] Loading states

### ✅ Experiment Monitor
- [x] WebSocket real-time updates
- [x] Polling fallback
- [x] Live convergence chart
- [x] Log streaming
- [x] Progress tracking
- [x] Export results button
- [x] Toast notifications

---

## Remaining Work (Optional Enhancements)

The following pages still use localStorage primarily but can be updated to use the API:

### 1. History Page (`/web/src/app/dashboard/history/page.tsx`)
- Currently uses localStorage
- **How to update**: Replace `localStorage.getItem("kanad_experiments")` with `api.getExperiments()`
- Add search/filter parameters to API call
- Add loading states with `CardSkeleton`
- Implement delete via `api.deleteExperiment(id)`
- Implement export via `api.downloadExperiment(id, format)`

### 2. Queue Page (`/web/src/app/dashboard/queue/page.tsx`)
- Currently uses localStorage
- **How to update**: Replace `localStorage.getItem("kanad_queue")` with `api.getQueue()`
- Use `api.updateQueueItem()` for priority/status changes
- Use `api.deleteQueueItem()` for deletion
- Use `api.reorderQueue()` for reordering

### 3. Settings Modal (`/web/src/components/settings/SettingsModal.tsx`)
- Currently uses localStorage
- **How to update**: Load settings with `api.getSettings()` in useEffect
- Save settings with `api.updateSettings()` instead of localStorage
- Add loading state while saving
- Show toast on success/error

### 4. Preview Window (`/web/src/components/simulation/PreviewWindow.tsx`)
- Already receives callbacks from dashboard
- **Optional**: Add SMILES validation using `api.validateSmiles()`
- Show validation feedback before submission

---

## Migration Strategy

If you have existing data in localStorage, you can migrate it to the API:

```typescript
import { migrateLocalStorageToAPI } from "@/lib/api";

// Run once when backend is ready
const results = await migrateLocalStorageToAPI();
console.log(`Migrated ${results.experiments} experiments, ${results.queue} queue items, settings: ${results.settings}`);
```

---

## Error Handling

The frontend handles errors gracefully:

1. **API Unavailable (503)**: Silently falls back to localStorage, no toast shown
2. **Network Error**: Shows toast "Failed to connect to API", uses localStorage
3. **Server Error (500)**: Shows toast with error message
4. **Validation Error (400)**: Shows toast with validation details
5. **Authentication Error (401)**: Would redirect to login (when auth is implemented)

---

## Backward Compatibility

The frontend maintains full backward compatibility:

- **With API**: All features work with real-time updates
- **Without API**: All features work with localStorage (offline mode)
- **Mixed Mode**: Uses API when available, localStorage as fallback

---

## Code Examples

### Using Toast Notifications
```typescript
import { useToast } from "@/components/ui/toast";

const toast = useToast();

// Success
toast.success("Experiment completed!");

// Error
toast.error("Failed to load data");

// Warning
toast.warning("Connection lost, using local data");

// Info
toast.info("Connecting to WebSocket...");
```

### Using API Functions
```typescript
import * as api from "@/lib/api";

// Create experiment
try {
  const response = await api.createExperiment({
    molecule: { smiles: "H2O", basis: "sto-3g", charge: 0, multiplicity: 1 },
    backendSettings: { method: "VQE", backend: "classical" },
    executeNow: true,
  });
  console.log("Experiment ID:", response.experimentId);
} catch (error) {
  console.error("Failed:", error.message);
}

// Get experiments with filters
const { experiments } = await api.getExperiments({
  status: "completed",
  limit: 20,
  offset: 0,
  search: "H2O"
});

// Real-time monitoring (WebSocket)
const ws = api.createWebSocket("exp_123");
ws.onmessage = (event) => {
  const message = JSON.parse(event.data);
  // Handle message
};

// Polling fallback
const stopPolling = api.pollExperimentStatus(
  "exp_123",
  (experiment) => {
    console.log("Status:", experiment.status);
  },
  2000 // poll every 2 seconds
);
```

### Using Loading Components
```typescript
import { FullPageLoader, Spinner, LoadingButton, CardSkeleton } from "@/components/ui/loading";

// Full page loading
{loading && <FullPageLoader message="Processing..." />}

// Inline spinner
<Spinner size="lg" />

// Loading button
<LoadingButton loading={isSubmitting} className="...">
  Submit
</LoadingButton>

// Skeleton while loading
{loading ? <CardSkeleton /> : <ActualCard data={data} />}
```

---

## Testing the Integration

### 1. Test Without Backend
```bash
# Frontend only
cd /home/mk/deeprealm/kanad/web
npm run dev
```
- All features should work with localStorage
- No error toasts for connection failures
- Graceful fallback to offline mode

### 2. Test With Backend
```bash
# Terminal 1: Backend
cd /home/mk/deeprealm/kanad/api
python -m uvicorn main:app --reload --port 8000

# Terminal 2: Frontend
cd /home/mk/deeprealm/kanad/web
npm run dev
```
- All features should use API
- Real-time updates via WebSocket
- Toast notifications for success/errors

### 3. Test WebSocket Fallback
- Disable WebSocket in backend → Should fall back to polling
- Frontend should handle gracefully and still get updates

---

## Backend Development Checklist

The backend agent should implement:

- [ ] All API endpoints listed above
- [ ] WebSocket support for real-time updates
- [ ] SMILES validation endpoint
- [ ] Export functionality (JSON/CSV)
- [ ] Pagination for experiments list
- [ ] Search/filter support
- [ ] Queue management
- [ ] Settings persistence (per-user)
- [ ] Error responses matching frontend expectations
- [ ] CORS configuration for `http://localhost:3000`
- [ ] Authentication (JWT) - optional for MVP

---

## Summary

The Kanad frontend is now **fully prepared** for backend integration with:

1. ✅ Complete API client library (`/web/src/lib/api.ts`)
2. ✅ Full TypeScript type definitions (`/web/src/lib/types.ts`)
3. ✅ Toast notification system
4. ✅ Loading states and skeletons
5. ✅ Real-time WebSocket + polling support
6. ✅ Export/download functionality
7. ✅ Graceful fallback to localStorage
8. ✅ Dashboard fully integrated with API
9. ✅ Experiment monitor with real-time updates
10. ✅ Error handling and retry logic

**The frontend will work seamlessly whether the backend is running or not**, providing a robust user experience in all scenarios.

---

## Questions?

Review this document and the implementation files:
- `/web/src/lib/api.ts` - API client library
- `/web/src/lib/types.ts` - Type definitions
- `/web/src/app/dashboard/page.tsx` - Dashboard with API integration
- `/web/src/components/simulation/ExperimentMonitor.tsx` - Real-time monitoring
- `/web/src/components/ui/toast.tsx` - Toast notifications
- `/web/src/components/ui/loading.tsx` - Loading components

The backend agent can use this document as a specification for implementing the FastAPI endpoints.
