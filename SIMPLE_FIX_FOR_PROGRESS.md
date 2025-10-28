# Simple Fix for VQE Progress Broadcasting

## The Problem

VQE solver (in `kanad/` package) cannot import from `api/` package due to Python module structure. We've been fighting with sys.path manipulation which doesn't work reliably.

## The SIMPLE Solution

**DON'T try to broadcast from inside VQE solver!**

Instead:
1. VQE solver already has a CALLBACK mechanism for progress
2. The experiment service (in `api/`) creates the VQE solver
3. The experiment service CAN import broadcast functions
4. Pass a callback from experiment service to VQE that broadcasts!

## Implementation

### In `/home/mk/deeprealm/kanad/kanad/utils/vqe_solver.py` line ~959:

**REMOVE all the broadcasting code** - just keep the callback:

```python
# Call user callback if provided
if hasattr(self, '_callback') and self._callback is not None:
    self._callback(self.iteration_count, energy, parameters)
```

### In `/home/mk/deeprealm/kanad/api/services/experiment_service.py` where VQE is created:

**CREATE a progress callback** and pass it to VQE:

```python
# Create progress callback for VQE
def vqe_progress_callback(iteration, energy, parameters):
    """Broadcast VQE progress to frontend"""
    try:
        broadcast_convergence_sync(experiment_id, iteration, float(energy))

        # Update progress bar
        if max_iterations > 0:
            progress = 20.0 + (iteration / max_iterations) * 60.0
            progress = min(progress, 80.0)
            JobDB.update_progress(experiment_id, progress=progress, current_iteration=iteration)
    except Exception as e:
        logger.error(f"Progress broadcast failed: {e}")

# Create VQE solver WITH CALLBACK
solver = VQESolver(
    bond=bond,
    ansatz_type=ansatz_type,
    ...
    callback=vqe_progress_callback  # <-- ADD THIS
)
```

This is the CLEAN solution - no import hacks, no sys.path manipulation, just using the existing callback mechanism!
