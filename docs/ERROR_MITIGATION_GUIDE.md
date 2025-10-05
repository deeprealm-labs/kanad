# Error Mitigation Guide for IBM Quantum

## Overview

The IBM VQE solver includes comprehensive error mitigation techniques to improve accuracy on NISQ hardware.

## Error Mitigation Stack

### 1. **Resilience Levels** (Base Mitigation)

```python
resilience_level=0  # No mitigation
resilience_level=1  # Measurement error mitigation
resilience_level=2  # ZNE + Measurement calibration (BEST)
```

### 2. **Dynamical Decoupling** (Circuit-level Mitigation)

Inserts pulse sequences during idle times to suppress decoherence:

```python
enable_dynamical_decoupling=True
dynamical_decoupling_sequence='XY4'  # Options: 'XX', 'XY4', 'CPMG'
```

**Sequences:**
- `XX`: Simple 2-pulse sequence (fastest)
- `XY4`: 4-pulse sequence (better suppression)
- `CPMG`: Carr-Purcell-Meiboom-Gill (advanced)

### 3. **Gate Twirling** (Noise Averaging)

Randomizes gate implementations to average out coherent errors:

```python
enable_twirling=True
twirling_strategy='all'  # Options: 'active', 'passive', 'all'
```

**Strategies:**
- `active`: Twirl gate operations
- `passive`: Twirl measurement basis
- `all`: Both (RECOMMENDED)

### 4. **Zero-Noise Extrapolation (ZNE)**

Extrapolates to zero-noise limit by running at multiple noise levels:

```python
resilience_level=2  # Required for ZNE
zne_extrapolator='exponential'  # Options: 'linear', 'exponential', 'polynomial'
```

**Extrapolators:**
- `linear`: Simple linear fit (fast)
- `exponential`: Exponential fit (better for most cases)
- `polynomial`: Higher-order polynomial (most flexible)

## Usage Examples

### Minimal Mitigation (Fastest)
```python
from kanad.backends.ibm import IBMVQESolver

vqe = IBMVQESolver(
    hamiltonian, ansatz,
    backend_name='ibm_torino',
    resilience_level=0,
    enable_dynamical_decoupling=False,
    enable_twirling=False
)
```

### Moderate Mitigation (Balanced)
```python
vqe = IBMVQESolver(
    hamiltonian, ansatz,
    backend_name='ibm_torino',
    resilience_level=1,  # Measurement calibration
    enable_dynamical_decoupling=True,
    dynamical_decoupling_sequence='XX',
    enable_twirling=True,
    twirling_strategy='passive'
)
```

### Maximum Mitigation (Best Accuracy)
```python
vqe = IBMVQESolver(
    hamiltonian, ansatz,
    backend_name='ibm_torino',
    shots=8192,  # More shots
    optimization_level=3,  # Max transpiler optimization
    resilience_level=2,  # ZNE + measurement
    enable_dynamical_decoupling=True,
    dynamical_decoupling_sequence='XY4',  # Best DD
    enable_twirling=True,
    twirling_strategy='all',  # All twirling
    zne_extrapolator='exponential'  # ZNE enabled
)
```

## Expected Impact

| Technique | Error Reduction | Runtime Overhead | Recommended For |
|-----------|----------------|------------------|-----------------|
| Resilience 1 | 20-40% | 1.2x | Always use |
| Resilience 2 (ZNE) | 40-60% | 2-3x | Final production runs |
| Dynamical Decoupling | 10-30% | 1.5x | Medium+ circuits |
| Gate Twirling | 15-35% | 1.3x | Always use |
| Full Stack | 60-80% | 3-5x | Best accuracy needed |

## Performance Considerations

### Shot Count Recommendations
- **No mitigation**: 1024-2048 shots
- **Moderate mitigation**: 4096-8192 shots
- **Maximum mitigation**: 8192-16384 shots

### Runtime Estimates (per iteration)
- **No mitigation**: ~30 seconds
- **Moderate**: ~60 seconds
- **Maximum**: ~120-180 seconds

## Troubleshooting

### Error still high after mitigation?
1. Increase `max_iterations` (try 50-100)
2. Increase `shots` (try 16384)
3. Try different `optimizer` (COBYLA → SLSQP)
4. Use better initial parameters
5. Consider simpler ansatz

### Circuit too deep for DD?
- Use `optimization_level=3` to reduce depth first
- Try simpler DD sequence (`XX` instead of `XY4`)
- Disable DD if circuit depth > 1000

### Runtime too slow?
- Reduce `resilience_level` to 1
- Disable ZNE (`zne_extrapolator=None`)
- Use fewer shots (4096 instead of 8192)
- Use Session mode for batch jobs

## Advanced: Custom Error Mitigation

For custom mitigation strategies, modify the `IBMVQESolver` directly:

```python
# Add custom transpiler passes
from qiskit.transpiler import PassManager
from qiskit.transpiler.passes import YourCustomPass

# In solve() method, after ISA transpilation:
custom_pm = PassManager([YourCustomPass()])
isa_circuit = custom_pm.run(isa_circuit)
```

## References

- [IBM Quantum Error Mitigation](https://quantum.ibm.com/docs/guides/error-mitigation)
- [Qiskit Runtime Primitives](https://quantum.ibm.com/docs/guides/primitives)
- [Dynamical Decoupling](https://quantum.ibm.com/docs/guides/dynamical-decoupling)
- [Zero-Noise Extrapolation](https://quantum.ibm.com/docs/guides/zero-noise-extrapolation)
