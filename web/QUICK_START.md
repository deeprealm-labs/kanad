# Quick Start Guide - Frontend API Integration

## For Frontend Developers

### Installation
```bash
cd /home/mk/deeprealm/kanad/web
npm install
```

### Run Development Server
```bash
npm run dev
```

Frontend will be available at http://localhost:3000

### Environment Setup
The `.env.local` file is already configured:
```bash
NEXT_PUBLIC_API_URL=http://localhost:8000
```

To change the API URL, edit `/web/.env.local`

### Using the API Client

```typescript
import * as api from "@/lib/api";
import { useToast } from "@/components/ui/toast";

function MyComponent() {
  const toast = useToast();

  const handleSubmit = async () => {
    try {
      const response = await api.createExperiment({
        molecule: { smiles: "H2O", basis: "sto-3g", charge: 0, multiplicity: 1 },
        backendSettings: { method: "VQE", backend: "classical" },
        executeNow: true
      });

      toast.success("Experiment started!");
      console.log("Experiment ID:", response.experimentId);
    } catch (error) {
      toast.error(error.message);
    }
  };

  return <button onClick={handleSubmit}>Start Experiment</button>;
}
```

### Key Files to Know
- `/web/src/lib/api.ts` - All API functions
- `/web/src/lib/types.ts` - TypeScript types
- `/web/src/components/ui/toast.tsx` - Toast notifications
- `/web/src/components/ui/loading.tsx` - Loading components

## For Backend Developers

### Start Backend Server
```bash
cd /home/mk/deeprealm/kanad/api
python -m uvicorn main:app --reload --host 0.0.0.0 --port 8000
```

Backend should be available at http://localhost:8000

### Required Endpoints
See `/web/API_INTEGRATION_COMPLETE.md` for complete specifications.

**Minimum endpoints needed:**
- `POST /api/experiments` - Create experiment
- `GET /api/experiments` - List experiments
- `GET /api/experiments/{id}` - Get experiment
- `WS /ws/experiments/{id}` - Real-time updates
- `GET /health` - Health check

### CORS Configuration
Add to your FastAPI app:
```python
from fastapi.middleware.cors import CORSMiddleware

app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:3000"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
```

### WebSocket Message Format
```json
{
  "type": "convergence",
  "experimentId": "exp_123",
  "data": {
    "iteration": 10,
    "energy": -1.125
  }
}
```

Types: `status`, `progress`, `convergence`, `log`, `result`, `error`

## Testing the Integration

### Test 1: Frontend Only (Offline Mode)
```bash
# Only start frontend
cd /home/mk/deeprealm/kanad/web
npm run dev
```
✅ Should work with localStorage fallback
✅ No error toasts about API

### Test 2: Full Stack
```bash
# Terminal 1: Backend
cd /home/mk/deeprealm/kanad/api
python -m uvicorn main:app --reload --port 8000

# Terminal 2: Frontend
cd /home/mk/deeprealm/kanad/web
npm run dev
```
✅ Should use API for all operations
✅ Real-time updates via WebSocket
✅ Toast notifications on success/error

### Test 3: Health Check
```bash
# Check if backend is running
curl http://localhost:8000/health

# Should return:
{"status": "ok", "message": "API is healthy"}
```

## Troubleshooting

### "Failed to connect to API"
- Check if backend is running on port 8000
- Check CORS configuration in backend
- Verify `NEXT_PUBLIC_API_URL` in `.env.local`

### "WebSocket connection failed"
- Backend should have WebSocket support at `/ws/experiments/{id}`
- Frontend will automatically fall back to polling
- No action needed, polling works fine

### "Cannot find module '@/lib/api'"
- Run `npm install` to ensure all dependencies are installed
- Check that TypeScript paths are configured in `tsconfig.json`

## Available Scripts

```bash
npm run dev       # Start development server
npm run build     # Build for production
npm run start     # Start production server
npm run lint      # Run ESLint
```

## Documentation

- **Complete API Guide**: `/web/API_INTEGRATION_COMPLETE.md`
- **File Changes**: `/web/INTEGRATION_SUMMARY.md`
- **This Guide**: `/web/QUICK_START.md`

## Support

For questions about:
- **Frontend**: Check `/web/src/lib/api.ts` for API functions
- **Backend**: Check `/web/API_INTEGRATION_COMPLETE.md` for endpoint specs
- **Types**: Check `/web/src/lib/types.ts` for all interfaces

## Status

🟢 **Frontend**: Ready for backend integration
🟡 **Backend**: Needs to be implemented
🟢 **Fallback**: Works offline with localStorage

---

**Happy coding!**
