# Frontend Complete Audit & Improvement Plan
**Date**: 2025-10-28
**Scope**: Full analysis of `/home/mk/deeprealm/kanad/web/src`

---

## 🚨 Critical Issues Found

### 1. Graph Type Confusion (Line vs Bar Chart)

**Problem**: Excited States VQE shows LINE chart instead of BAR chart

**Root Cause**:
```typescript
// Line 815 in ExperimentMonitor.tsx
if (method === "SQD" || method === "EXCITED_STATES") {
    return <BarChart />  // Shows bar chart
}
// Line 872
return <LineChart />  // Shows line chart
```

When user selects "Excited States" with "VQE" method:
- Config has `method: "VQE"` (not "EXCITED_STATES")
- Frontend checks `if (method === "EXCITED_STATES")` → FALSE
- Falls through to LineChart → WRONG!

**Solution**: Check BOTH method AND excited states configuration
```typescript
const isExcitedStates = method === "EXCITED_STATES" ||
                        method === "VQE" && (config?.excitedMethod || config?.excited_method);

const isSQD = method === "SQD";

if (isExcitedStates || isSQD) {
    return <BarChart />
} else {
    return <LineChart />
}
```

---

### 2. Data Loss on Reopen (Missing Iteration 1)

**Problem**: Graph shows iterations 2-3 instead of 1-3

**Root Cause**: Look at line 118 in ExperimentMonitor.tsx:
```typescript
// Only restore convergence data if we don't have it yet
if (convergenceData.length === 0) {
    if (exp.results.convergence_history) {
        setConvergenceData(exp.results.convergence_history);
    }
}
```

**BUT**: `convergenceData.length === 0` check fails if component has ANY existing data!

**Also**: The component fetches experiment data which might not have ALL convergence points if some were only sent via WebSocket.

**Solution**:
1. ALWAYS merge historical data with current data
2. Backend should save ALL convergence points to results
3. Frontend should deduplicate by iteration number

```typescript
// Merge historical with current, deduplicate by iteration
if (exp.results.convergence_history) {
    setConvergenceData(prev => {
        const merged = [...prev, ...exp.results.convergence_history];
        const unique = Array.from(
            new Map(merged.map(item => [item.iteration, item])).values()
        );
        return unique.sort((a, b) => a.iteration - b.iteration);
    });
}
```

---

### 3. Iteration Count Mismatch

**Problem**:
- Progress shows "Iteration 3 / 100"
- Current Iteration shows "2"
- Graph shows 2-3 (missing 1)

**Root Cause**: Multiple state update paths not synchronized

**Lines Involved**:
- Line 109: `setCurrentIteration(exp.job.current_iteration)`
- Line 256: `setCurrentIteration((prev) => Math.max(prev, exp.job.current_iteration))`
- Line 280: `setCurrentIteration((prev) => Math.max(prev, exp.convergenceData.length))`

**Solution**: Use SINGLE source of truth
```typescript
// ALWAYS use job.current_iteration if available
if (exp.job?.current_iteration !== undefined) {
    setCurrentIteration(exp.job.current_iteration);
} else {
    // Fallback to convergence data length
    setCurrentIteration(convergenceData.length);
}
```

---

### 4. Hardcoded Values Throughout Codebase

#### ExperimentMonitor.tsx

**Line 655**: `{experimentConfig?.backendSettings?.method || "VQE"}` - Default "VQE"
**Line 661**: `{experimentConfig?.backendSettings?.backend || "Classical"}` - Default "Classical"
**Line 670**: `{experimentConfig?.backendSettings?.ansatz || "HEA"}` - Default "HEA"
**Line 681**: `{experimentConfig?.backendSettings?.subspaceDim || 10}` - Default 10
**Line 687**: `{experimentConfig?.backendSettings?.nStates || 3}` - Default 3
**Line 699**: `...|| "CIS"}` - Default "CIS"
**Line 705**: `...|| 5}` - Default 5 states
**Line 714**: `{experimentConfig?.molecule?.basis || "STO-3G"}` - Default "STO-3G"
**Line 741**: `maxIterations = experimentConfig?.backendSettings?.maxIterations || 100;` - Default 100
**Line 743**: `maxIterations = 6;` - Hardcoded SQD stages
**Line 748**: `maxIterations = 6;` - Hardcoded for bluequbit

**Problem**: If backend config is missing, shows misleading default values

**Solution**: Show "N/A" or "Loading..." instead of fake defaults
```typescript
{experimentConfig?.backendSettings?.method || <span className="text-muted">N/A</span>}
```

---

#### ConfigurationSelector.tsx

Need to audit this file for:
- Default basis sets
- Default method selections
- Default backend selections

---

#### SettingsModal.tsx

Need to audit for:
- Default VQE parameters
- Default SQD parameters
- Default excited states configuration

---

### 5. Chart Label Issues

**Line 836**: `label: { value: "Energy State"` - Should say "Excited State #" for EXCITED_STATES

**Line 858**: `labelFormatter={(label) => \`State ${label - stateOffset}\``
- This calculates state number by subtracting offset
- If iteration 1 is missing, State 0 will be shown (WRONG!)

**Solution**: Use actual state number from iteration
```typescript
labelFormatter={(label) => `S${label}`}  // S1, S2, S3 instead of calculation
```

---

### 6. Progress Bar Edge Cases

**Line 725**: Shows `{progress}%` - could be NaN or undefined

**Solution**: Safe rendering
```typescript
{progress !== undefined ? `${Math.round(progress)}%` : '0%'}
```

---

### 7. Live Metrics Display

**Lines 767-794**: Current Energy and Current Iteration

**Issues**:
- Uses `-` when no data (should use consistent placeholder)
- Energy formatting inconsistent (sometimes scientific, sometimes decimal)
- No units shown for large/small energies

**Solution**: Consistent formatting utilities
```typescript
const formatEnergy = (e: number | undefined) => {
    if (e === undefined) return 'N/A';
    if (Math.abs(e) > 1000 || Math.abs(e) < 0.001) {
        return e.toExponential(6) + ' Ha';
    }
    return e.toFixed(6) + ' Ha';
};
```

---

## 📊 Hardcoded Values - Complete List

### Configuration Defaults
| Location | Value | Issue |
|----------|-------|-------|
| ExperimentMonitor:655 | "VQE" | Method default |
| ExperimentMonitor:661 | "Classical" | Backend default |
| ExperimentMonitor:670 | "HEA" | Ansatz default |
| ExperimentMonitor:681 | 10 | Subspace dim default |
| ExperimentMonitor:687 | 3 | N states default |
| ExperimentMonitor:699 | "CIS" | ES method default |
| ExperimentMonitor:705 | 5 | ES states default |
| ExperimentMonitor:714 | "STO-3G" | Basis default |
| ExperimentMonitor:741 | 100 | Max iterations default |
| ExperimentMonitor:743 | 6 | SQD stages hardcoded |
| ExperimentMonitor:748 | 6 | BluQubit redirect hardcoded |

### UI/UX Constants
| Location | Value | Issue |
|----------|-------|-------|
| ExperimentMonitor:728-732 | Progress bar styling | Fixed height, color |
| ExperimentMonitor:862 | "#ea580c" | Bar color hardcoded |
| ExperimentMonitor:898 | "#ea580c" | Line color hardcoded |
| ExperimentMonitor:900 | `dot={false}` | No data points on line |
| ExperimentMonitor:840 | `fontSize: 11` | Small, hard to read |

---

## 🎨 UI/UX Issues

### Color Scheme
**Problem**: Uses generic orange (`#ea580c`) everywhere
**Solution**: Professional scientific color palette
```css
/* Energy levels (sequential) */
--energy-0: #0d47a1;  /* Deep blue - ground state */
--energy-1: #1976d2;  /* Blue */
--energy-2: #42a5f5;  /* Light blue */
--energy-3: #64b5f6;  /* Lighter blue */
/* ... more states */

/* Convergence (diverging) */
--converged: #2e7d32;    /* Green */
--converging: #f57c00;   /* Orange */
--not-converged: #c62828; /* Red */
```

### Typography
**Problem**: Font sizes too small (11px), hard to read
**Solution**: Scientific publication standards
```css
--font-size-label: 14px;
--font-size-axis: 12px;
--font-size-title: 16px;
--font-size-data: 18px;  /* For energy values */
```

### Layout
**Problem**: Fixed heights, doesn't respond to content
**Current**: `min-h-0 flex flex-col` with manual sizing
**Solution**: CSS Grid with proper fr units
```css
.monitor-layout {
    display: grid;
    grid-template-columns: 300px 1fr 350px;
    grid-template-rows: auto 1fr auto;
    gap: 1.5rem;
    height: 100vh;
}
```

### Spacing
**Problem**: Inconsistent padding/margins
**Current**: Mix of `p-4`, `mb-3`, `space-y-2`
**Solution**: Design system with consistent scale
```typescript
const spacing = {
    xs: '0.25rem',  // 4px
    sm: '0.5rem',   // 8px
    md: '1rem',     // 16px
    lg: '1.5rem',   // 24px
    xl: '2rem',     // 32px
};
```

---

## 🔧 Proposed Solutions

### Immediate Fixes (Week 1)

#### Fix 1: Graph Type Detection
```typescript
// ExperimentMonitor.tsx line ~810
const isExcitedStatesExperiment = () => {
    const method = experimentConfig?.backendSettings?.method;
    const excitedMethod = experimentConfig?.backendSettings?.excited_method ||
                         experimentConfig?.backendSettings?.excitedMethod;

    return method === "EXCITED_STATES" ||
           (method === "VQE" && excitedMethod);
};

// Use in chart selection
{isExcitedStatesExperiment() || method === "SQD" ? (
    <BarChart data={convergenceData}>
        {/* Bar chart rendering */}
    </BarChart>
) : (
    <LineChart data={convergenceData}>
        {/* Line chart rendering */}
    </LineChart>
)}
```

#### Fix 2: Data Persistence with Merge
```typescript
// ExperimentMonitor.tsx line ~118
if (exp.results) {
    setResults(exp.results);

    // MERGE historical data with current (don't replace!)
    if (exp.results.convergence_history) {
        setConvergenceData(prev => {
            const historical = exp.results.convergence_history;

            // Combine and deduplicate
            const allData = [...prev, ...historical];
            const uniqueData = Array.from(
                new Map(allData.map(item => [item.iteration, item])).values()
            );

            // Sort by iteration
            return uniqueData.sort((a, b) => a.iteration - b.iteration);
        });

        // Update iteration count to highest value
        const maxIter = Math.max(...exp.results.convergence_history.map(d => d.iteration));
        setCurrentIteration(prev => Math.max(prev, maxIter));
    }
}
```

#### Fix 3: Remove Hardcoded Defaults
```typescript
// Create utility component
const ConfigValue = ({ value, fallback = "N/A", unit = "" }) => (
    <span className={value ? "font-medium" : "text-muted-foreground italic"}>
        {value ? `${value}${unit}` : fallback}
    </span>
);

// Use throughout
<div className="flex justify-between">
    <span className="text-muted-foreground">Method:</span>
    <ConfigValue value={experimentConfig?.backendSettings?.method} />
</div>
```

#### Fix 4: Consistent Iteration Display
```typescript
// Single source of truth for current iteration
const getCurrentIteration = () => {
    // Priority 1: Job data
    if (experiment?.job?.current_iteration !== undefined) {
        return experiment.job.current_iteration;
    }

    // Priority 2: Convergence data length
    if (convergenceData.length > 0) {
        return Math.max(...convergenceData.map(d => d.iteration));
    }

    // Priority 3: State variable
    return currentIteration;
};

// Use everywhere
<div>Iteration {getCurrentIteration()} / {getMaxIterations()}</div>
```

---

### Medium-term Improvements (Week 2-3)

#### Professional Chart Component
```typescript
// components/charts/EnergyChart.tsx
interface EnergyChartProps {
    data: ConvergencePoint[];
    type: 'convergence' | 'spectrum';
    method: string;
    isLive?: boolean;
}

export const EnergyChart = ({ data, type, method, isLive }: EnergyChartProps) => {
    const Chart component = type === 'spectrum' ? BarChart : LineChart;

    const chartConfig = {
        convergence: {
            title: 'Energy Convergence',
            xLabel: 'Iteration',
            yLabel: 'Energy (Hartree)',
            color: COLORS.primary,
            showDots: true,
        },
        spectrum: {
            title: 'Energy Spectrum',
            xLabel: 'State Number',
            yLabel: 'Energy (Hartree)',
            color: COLORS.states,
            showLabels: true,
        }
    };

    const config = chartConfig[type];

    return (
        <ResponsiveContainer>
            <ChartComponent data={data}>
                <CartesianGrid {...GRID_CONFIG} />
                <XAxis label={config.xLabel} {...AXIS_CONFIG} />
                <YAxis label={config.yLabel} {...AXIS_CONFIG} />
                <Tooltip {...TOOLTIP_CONFIG} />
                {type === 'convergence' ? (
                    <Line {...config} />
                ) : (
                    <Bar {...config} />
                )}
            </ChartComponent>
        </ResponsiveContainer>
    );
};
```

#### Design System
```typescript
// lib/design-system.ts
export const COLORS = {
    // Energy states (sequential blue scale)
    states: [
        '#0d47a1', '#1565c0', '#1976d2', '#1e88e5',
        '#2196f3', '#42a5f5', '#64b5f6', '#90caf9'
    ],

    // Convergence
    converged: '#2e7d32',
    converging: '#f57c00',
    notConverged: '#c62828',

    // Methods
    vqe: '#7c4dff',
    hf: '#00acc1',
    sqd: '#f4511e',
    excited: '#ff6f00',
};

export const TYPOGRAPHY = {
    title: { fontSize: 18, fontWeight: 600 },
    subtitle: { fontSize: 16, fontWeight: 500 },
    label: { fontSize: 14, fontWeight: 400 },
    data: { fontSize: 20, fontWeight: 700, fontFamily: 'monospace' },
};

export const SPACING = {
    xs: 4, sm: 8, md: 16, lg: 24, xl: 32, xxl: 48
};
```

---

## 📁 Files Requiring Updates

### Critical Priority
1. ✅ `ExperimentMonitor.tsx` - Graph type, data persistence, iteration count
2. `SettingsModal.tsx` - Remove defaults, improve UX
3. `ConfigurationSelector.tsx` - Better method selection
4. `lib/api.ts` - Ensure all data properly fetched

### High Priority
5. `DashboardHome.tsx` - Consistent styling
6. `AnalysisResults.tsx` - Professional formatting
7. `ExperimentReport.tsx` - Better layout

### Medium Priority
8. `lib/utils.ts` - Add formatting utilities
9. `lib/types.ts` - Stronger type definitions
10. New: `components/charts/` - Professional chart components
11. New: `lib/design-system.ts` - Centralized design tokens

---

## 🚀 Implementation Plan

### Sprint 1: Critical Fixes (Days 1-3)
- [x] Fix graph type detection for excited states
- [ ] Fix data merge/deduplication on reopen
- [ ] Fix iteration count consistency
- [ ] Remove/replace hardcoded defaults

### Sprint 2: Professional UI (Days 4-7)
- [ ] Implement design system
- [ ] Create reusable chart components
- [ ] Improve typography and spacing
- [ ] Responsive layout improvements

### Sprint 3: Polish (Days 8-10)
- [ ] Animations and transitions
- [ ] Loading states and skeletons
- [ ] Error handling and user feedback
- [ ] Performance optimization

---

## ✅ Success Criteria

1. **Graph Type**: Excited states always show bar chart
2. **Data Persistence**: ALL convergence points visible on reopen (including iteration 1)
3. **Iteration Count**: Consistent across all displays
4. **No Hardcoded Defaults**: Show actual values or N/A
5. **Professional Appearance**: Matches scientific software standards
6. **Responsive**: Works on laptop (1920x1080), desktop (2560x1440), ultrawide
7. **Accessible**: Keyboard navigation, screen reader support
8. **Performance**: < 100ms render time with 1000 data points

