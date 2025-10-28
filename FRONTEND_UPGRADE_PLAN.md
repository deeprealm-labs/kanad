# Frontend Upgrade Plan - Kanad Quantum Chemistry Platform

## Current Critical Issues Identified

### 1. **Progress Bar Stuck at 0%**
**Root Cause**: Frontend fetches `experiment` object which doesn't include `job.progress`
- Backend updates `job.progress` correctly (verified in logs: `progress=80.0%`)
- Frontend polls `/api/experiments/{id}` which returns experiment data
- **Missing**: Job progress is in separate `jobs` table, not included in experiment response

**Fix Required**:
- Modify `/api/experiments/{id}` endpoint to include related job data
- OR: Frontend should fetch job data separately `/api/jobs/{id}`
- OR: WebSocket should broadcast job progress updates

---

### 2. **Cancellation Not Working**
**Symptoms**: Click cancel → shows "Cancelling..." → experiment continues running
**Investigation Needed**:
- Check `/api/experiments/{id}/cancel` endpoint implementation
- Verify cancellation signal reaches running experiment process
- Check if multiprocessing jobs handle cancellation signals

---

### 3. **Graph Data Loss on Monitor Reopen**
**Root Cause**: Lines 68-96 in ExperimentMonitor.tsx - state resets on component remount
```typescript
useEffect(() => {
  setCurrentIteration(0);
  setConvergenceData([]);  // ❌ Data cleared!
  setProgress(0);
  // ...
}, [experimentId]);
```

**Issue**: When user navigates away and returns, `experimentId` changes → state resets
**Fix**: Fetch historical convergence data from API on mount, don't clear existing data

---

### 4. **Graph Starting from Later Iterations**
**Cause**: WebSocket only sends NEW convergence points, not historical ones
**Fix**: On component mount, fetch all existing convergence data before listening to WebSocket

---

### 5. **No Multi-Experiment Job Queue**
**Current**: Can only monitor one experiment at a time
**Needed**:
- View all running experiments
- Click to monitor any running experiment
- Queue visualization
- Parallel experiment execution

---

## Architecture Issues

### Current Architecture Problems:

1. **State Management**:
   - No persistent state (Redux/Zustand)
   - State resets on navigation
   - No caching of experiment data
   - Convergence data not persisted

2. **WebSocket Handling**:
   - No reconnection logic
   - No message queuing during disconnection
   - State lost if WebSocket closes
   - Falls back to polling (inefficient)

3. **API Integration**:
   - Experiment and Job data in separate endpoints
   - No unified data fetching
   - No optimistic updates
   - Polling is inefficient

4. **UI/UX Issues**:
   - Not professional/scientific appearance
   - Poor responsive design
   - Inconsistent sizing
   - No keyboard shortcuts
   - Limited accessibility

---

## Proposed Solutions

### Phase 1: Critical Bug Fixes (Immediate)

#### 1.1 Fix Progress Bar
**Backend Change**:
```python
# api/routes/experiments.py
@router.get("/{experiment_id}")
async def get_experiment(experiment_id: str):
    experiment = ExperimentDB.get_by_id(experiment_id)
    job = JobDB.get_by_experiment_id(experiment_id)  # NEW

    return {
        **experiment,
        "job": {  # NEW: Include job data
            "progress": job.progress,
            "current_iteration": job.current_iteration,
            "status": job.status
        }
    }
```

**Frontend Change**:
```typescript
// Use job.progress instead of experiment.status-based progress
setProgress(experiment.job?.progress || 0);
setCurrentIteration(experiment.job?.current_iteration || 0);
```

#### 1.2 Fix Cancellation
**Backend Investigation**:
- Check `check_cancellation()` in experiment_service.py
- Verify signal propagation to VQE/SQD solvers
- Add cancellation polling in solver loops

**Frontend**: Add timeout for cancellation (show error if not cancelled after 10s)

#### 1.3 Fix Graph Data Persistence
**Remove State Reset**:
```typescript
useEffect(() => {
  if (!experimentId) return;

  // DON'T reset state if data exists
  // ONLY fetch if we don't have data yet
  if (convergenceData.length === 0) {
    fetchHistoricalConvergenceData(experimentId);
  }
}, [experimentId]);
```

**Add Historical Data Fetching**:
```typescript
const fetchHistoricalConvergenceData = async (expId: string) => {
  const exp = await api.getExperiment(expId);
  if (exp.results?.convergence_history) {
    setConvergenceData(exp.results.convergence_history);
  }
};
```

---

### Phase 2: Architecture Improvements

#### 2.1 Implement State Management (Zustand)
```typescript
// store/experimentStore.ts
import create from 'zustand';
import { persist } from 'zustand/middleware';

interface ExperimentStore {
  experiments: Map<string, ExperimentData>;
  convergenceData: Map<string, ConvergencePoint[]>;
  activeExperiments: string[];

  addConvergencePoint: (expId: string, point: ConvergencePoint) => void;
  updateExperiment: (expId: string, data: Partial<ExperimentData>) => void;
  // ...
}

export const useExperimentStore = create<ExperimentStore>()(
  persist(
    (set) => ({
      experiments: new Map(),
      convergenceData: new Map(),
      activeExperiments: [],

      addConvergencePoint: (expId, point) => set((state) => {
        const data = state.convergenceData.get(expId) || [];
        state.convergenceData.set(expId, [...data, point]);
        return { convergenceData: new Map(state.convergenceData) };
      }),
    }),
    { name: 'experiment-store' }
  )
);
```

#### 2.2 Robust WebSocket Manager
```typescript
//lib/websocket-manager.ts
class WebSocketManager {
  private connections: Map<string, WebSocket> = new Map();
  private reconnectTimers: Map<string, NodeJS.Timeout> = new Map();
  private messageQueues: Map<string, any[]> = new Map();

  connect(experimentId: string, onMessage: (msg: any) => void) {
    const ws = new WebSocket(`ws://localhost:8000/api/ws/experiments/${experimentId}`);

    ws.onopen = () => {
      console.log(`✅ WebSocket connected: ${experimentId}`);
      this.flushMessageQueue(experimentId);
    };

    ws.onmessage = (event) => {
      const message = JSON.parse(event.data);
      onMessage(message);
    };

    ws.onclose = () => {
      console.log(`🔌 WebSocket closed: ${experimentId}, attempting reconnect...`);
      this.scheduleReconnect(experimentId, onMessage);
    };

    this.connections.set(experimentId, ws);
  }

  private scheduleReconnect(experimentId: string, onMessage: (msg: any) => void) {
    const timer = setTimeout(() => {
      console.log(`🔄 Reconnecting WebSocket: ${experimentId}`);
      this.connect(experimentId, onMessage);
    }, 3000);

    this.reconnectTimers.set(experimentId, timer);
  }

  disconnect(experimentId: string) {
    const ws = this.connections.get(experimentId);
    if (ws) {
      ws.close();
      this.connections.delete(experimentId);
    }

    const timer = this.reconnectTimers.get(experimentId);
    if (timer) {
      clearTimeout(timer);
      this.reconnectTimers.delete(experimentId);
    }
  }
}

export const wsManager = new WebSocketManager();
```

#### 2.3 Unified Data Fetching with React Query
```typescript
// hooks/useExperiment.ts
import { useQuery, useMutation } from '@tanstack/react-query';

export const useExperiment = (experimentId: string) => {
  return useQuery({
    queryKey: ['experiment', experimentId],
    queryFn: () => api.getExperiment(experimentId),
    refetchInterval: (data) => {
      // Poll only if running
      return data?.status === 'running' ? 5000 : false;
    },
  });
};

export const useCancelExperiment = () => {
  return useMutation({
    mutationFn: (experimentId: string) => api.cancelExperiment(experimentId),
    onSuccess: () => {
      queryClient.invalidateQueries(['experiment']);
    },
  });
};
```

---

### Phase 3: UI/UX Redesign

#### 3.1 Professional Scientific Tool Design

**Design Principles**:
- Clean, data-dense interface
- Monospace fonts for data
- High-contrast colors
- Clear visual hierarchy
- Scientific color schemes (sequential/diverging)

**Layout Improvements**:
```typescript
// New responsive grid system
<div className="grid grid-cols-12 gap-6 h-screen">
  {/* Left Sidebar - Config */}
  <div className="col-span-3 border-r">
    <ConfigPanel />
  </div>

  {/* Center - Main Viz */}
  <div className="col-span-6 flex flex-col">
    <ConvergenceGraph className="flex-1 min-h-[400px]" />
    <MetricsPanel className="h-32" />
  </div>

  {/* Right Sidebar - Logs/Queue */}
  <div className="col-span-3 border-l flex flex-col">
    <JobQueue className="flex-1" />
    <ExecutionLogs className="h-48" />
  </div>
</div>
```

#### 3.2 Professional Components

**Progress Visualization**:
- Linear progress bar (current)
- + Circular progress indicator
- + Stage indicators (Init → Optimize → Analyze)
- + ETA display
- + Throughput metrics (iters/sec)

**Convergence Graph**:
- Interactive zoom/pan
- Export to SVG/PNG
- Multiple series (energy, gradient, parameters)
- Log scale option
- Error bars
- Comparison mode (multiple experiments)

**Metrics Dashboard**:
- Real-time energy display (large, prominent)
- Iteration counter
- Wall time / CPU time
- Convergence criteria status
- Optimization history (mini graph)

---

### Phase 4: Multi-Experiment Job Queue

#### 4.1 Job Queue Visualization
```typescript
// components/queue/JobQueuePanel.tsx
interface JobQueueItem {
  id: string;
  experimentId: string;
  status: 'queued' | 'running' | 'completed' | 'failed';
  progress: number;
  startedAt?: Date;
  estimatedCompletion?: Date;
}

export const JobQueuePanel = () => {
  const { jobs } = useJobQueue();
  const [selectedJob, setSelectedJob] = useState<string | null>(null);

  return (
    <div className="flex flex-col h-full">
      <h3 className="text-lg font-semibold p-4 border-b">
        Job Queue ({jobs.filter(j => j.status === 'running').length} running)
      </h3>

      <div className="flex-1 overflow-y-auto">
        {jobs.map(job => (
          <JobQueueItem
            key={job.id}
            job={job}
            selected={selectedJob === job.id}
            onClick={() => setSelectedJob(job.id)}
            onMonitor={() => openMonitor(job.experimentId)}
          />
        ))}
      </div>
    </div>
  );
};
```

#### 4.2 Multi-Monitor Support
```typescript
// components/monitor/MultiExperimentMonitor.tsx
export const MultiExperimentMonitor = () => {
  const [monitoredExperiments, setMonitoredExperiments] = useState<string[]>([]);

  return (
    <div className="grid grid-cols-2 gap-4">
      {monitoredExperiments.map(expId => (
        <ExperimentMonitorCard
          key={expId}
          experimentId={expId}
          onClose={() => removeMonitor(expId)}
          compact={true}
        />
      ))}
    </div>
  );
};
```

---

## Implementation Timeline

### Week 1: Critical Fixes
- [ ] Day 1-2: Fix progress bar (backend + frontend)
- [ ] Day 3: Fix cancellation
- [ ] Day 4: Fix graph data persistence
- [ ] Day 5: Testing and bug fixes

### Week 2: Architecture
- [ ] Day 1-2: Implement Zustand state management
- [ ] Day 3-4: Robust WebSocket manager
- [ ] Day 5: React Query integration

### Week 3: UI/UX
- [ ] Day 1-2: New design system (colors, fonts, spacing)
- [ ] Day 3-4: Professional components (graphs, metrics)
- [ ] Day 5: Responsive layout improvements

### Week 4: Job Queue
- [ ] Day 1-2: Job queue backend enhancements
- [ ] Day 3-4: Job queue UI and multi-monitor
- [ ] Day 5: Integration testing and polish

---

## Success Metrics

1. **Progress Bar**: Updates smoothly from 0% to 100%
2. **Cancellation**: Experiment stops within 2 seconds
3. **Graph Persistence**: Data retained on navigation
4. **Multi-Experiment**: Monitor 4+ experiments simultaneously
5. **Performance**: UI remains responsive with 100+ convergence points
6. **Professional Look**: Passes design review for scientific tool

---

## Next Steps

1. Approve this plan
2. Start with Phase 1 (Critical Fixes)
3. Iterate based on feedback
4. Deploy incrementally

