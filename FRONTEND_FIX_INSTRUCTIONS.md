# Frontend 404 Errors - FIXED

## Problem
The frontend was calling routes without the `/api` prefix:
- ❌ `http://localhost:8000/experiments/submit`
- ❌ `http://localhost:8000/settings/defaults`
- ❌ `http://localhost:8000/circuits/preview`

But the backend routes are:
- ✅ `http://localhost:8000/api/experiments/submit`
- ✅ `http://localhost:8000/api/settings/defaults`
- ✅ `http://localhost:8000/api/circuits/preview`

## Root Cause
The `web/.env.local` file had the wrong API URL:
```
NEXT_PUBLIC_API_URL=http://localhost:8000  ← WRONG
```

## Fix Applied
Updated `web/.env.local` to include `/api` prefix:
```
NEXT_PUBLIC_API_URL=http://localhost:8000/api  ← CORRECT
```

## Action Required: Restart Next.js Dev Server

**You MUST restart the Next.js dev server for the change to take effect:**

1. Go to the terminal where Next.js is running
2. Press `Ctrl+C` to stop the server
3. Restart it:
   ```bash
   cd web
   npm run dev
   ```

## Verification

After restarting, check the browser console. You should see requests to:
- ✅ `http://localhost:8000/api/experiments/submit`
- ✅ `http://localhost:8000/api/settings/defaults`
- ✅ `http://localhost:8000/api/circuits/preview`

All 404 errors should be gone!

## Backend Status

The backend is running correctly:
- ✅ Server: `http://localhost:8000`
- ✅ Health: `http://localhost:8000/health`
- ✅ All API routes working at `/api/*`

---

**File modified:** `web/.env.local`
**Change:** `NEXT_PUBLIC_API_URL` updated from `http://localhost:8000` to `http://localhost:8000/api`
**Next step:** Restart Next.js dev server (Ctrl+C then `npm run dev`)
