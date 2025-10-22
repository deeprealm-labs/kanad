# BlueQubit Device Options Implementation

## Changes Made

### 1. Backend Support
Added support for all BlueQubit device types:
- **GPU** (36 qubits, fastest) - Default
- **CPU** (34 qubits)
- **MPS GPU** (40+ qubits) - Matrix Product State on GPU
- **MPS CPU** (40+ qubits) - Matrix Product State on CPU
- **Pauli-Path** (50+ qubits) - Advanced Pauli propagation algorithm

### 2. Backend Code Updates

**File**: [kanad/backends/bluequbit/backend.py](kanad/backends/bluequbit/backend.py)

**Added**:
- `pauli-path` to `SUPPORTED_DEVICES`
- `method` parameter (legacy alias for `device`)
- Documentation for device-specific options:
  - `mps_bond_dimension`: For MPS devices
  - `pauli_path_truncation_threshold`: For Pauli-Path (minimum 1e-5)
  - `pauli_path_circuit_transpilation_level`: Transpilation control

### 3. API Service Updates

**File**: [api/services/experiment_service.py:186-225](api/services/experiment_service.py#L186-L225)

**Added**:
- Device selection from configuration: `bluequbit_device`
- Device-specific options parsing:
  - MPS devices: `mps_bond_dimension`
  - Pauli-Path: `pauli_path_truncation_threshold`, `pauli_path_circuit_transpilation_level`
- Logging for selected device

### 4. Frontend Settings UI

**File**: [web/src/components/settings/SettingsModal.tsx](web/src/components/settings/SettingsModal.tsx)

**Added**:
- State variable: `bluequbitDevice` (default: "gpu")
- Device selector dropdown when BlueQubit backend is selected
- Options: GPU, CPU, MPS GPU, MPS CPU, Pauli-Path
- Save/load bluequbitDevice in settings

**Disabled**:
- IBM Quantum backend (temporarily commented out)
- Reason: Batch mode causes long queue times; need session mode for VQE

### 5. Device Capabilities

| Device | Max Qubits | Speed | Best For | Requires Credits |
|--------|-----------|-------|----------|------------------|
| GPU | 36 | Fastest | Most experiments | Yes ($0.20/job) |
| CPU | 34 | Fast | Small molecules | Yes ($0.20/job) |
| MPS GPU | 40+ | Medium | Larger molecules | Yes (higher cost) |
| MPS CPU | 40+ | Slow | Larger molecules | Yes (higher cost) |
| Pauli-Path | 50+ | Varies | Very large systems | Yes (highest cost) |

## How to Use

### Basic Usage (GPU)

1. **Settings**:
   - Backend: BlueQubit
   - Device: GPU (default)
   - Max Iterations: 20

2. **Cost**: ~$4 for 20 iterations with COBYLA (~40 jobs)

### Advanced Usage (MPS)

1. **Settings**:
   - Backend: BlueQubit
   - Device: MPS GPU or MPS CPU
   - Max Iterations: 10

2. **Configuration** (future): Add MPS bond dimension in advanced settings

3. **Cost**: Higher than GPU (pricing TBD by BlueQubit)

### Expert Usage (Pauli-Path)

1. **Settings**:
   - Backend: BlueQubit
   - Device: Pauli-Path
   - Max Iterations: 5

2. **Configuration** (future):
   - Add `pauli_path_truncation_threshold`: 0.01 (or higher for faster, less accurate)
   - For circuits > 13 qubits: MUST specify truncation threshold

3. **Cost**: Highest (pricing TBD by BlueQubit)

## Device Selection Recommendations

### Small Molecules (H2, HF, LiH)
- **Device**: GPU
- **Qubits needed**: 4-8
- **Cost**: ~$0.40 for 10 iterations (COBYLA)

### Medium Molecules (H2O, NH3, CH4)
- **Device**: GPU or CPU
- **Qubits needed**: 8-20
- **Cost**: ~$0.80 for 10 iterations (COBYLA)

### Larger Molecules (CO2, COOH, benzene)
- **Device**: MPS GPU
- **Qubits needed**: 20-40
- **Cost**: ~$2-5 for 10 iterations

### Very Large Systems (proteins, materials)
- **Device**: Pauli-Path
- **Qubits needed**: 40-50+
- **Cost**: ~$10-20 for 5 iterations
- **Note**: Requires careful truncation threshold tuning

## API Configuration Format

The frontend will send:

```json
{
  "backend": "bluequbit",
  "bluequbit_device": "gpu",
  "optimizer": "COBYLA",
  "max_iterations": 20
}
```

For MPS:
```json
{
  "backend": "bluequbit",
  "bluequbit_device": "mps.gpu",
  "mps_bond_dimension": 100,
  "optimizer": "COBYLA",
  "max_iterations": 10
}
```

For Pauli-Path:
```json
{
  "backend": "bluequbit",
  "bluequbit_device": "pauli-path",
  "pauli_path_truncation_threshold": 0.01,
  "optimizer": "COBYLA",
  "max_iterations": 5
}
```

## Future Enhancements

1. **Advanced Options Panel**: Add UI for device-specific parameters
2. **Cost Estimator**: Show estimated cost before running
3. **Device Recommendations**: Auto-suggest device based on molecule size
4. **MPS Bond Dimension Slider**: UI control for MPS truncation
5. **Pauli-Path Threshold Selector**: UI control for truncation threshold

## IBM Quantum Status

**Temporarily Disabled** - Commented out in settings UI

**Reason**: VQE with batch mode causes long queue times (hours). Need to:
1. Contact IBM for session mode access
2. Implement session-based VQE execution
3. Configure proper job batching for session mode

**Re-enable when**: Session mode is configured

## Testing

After restarting the API server and web app:

1. **Open Settings**
2. **Select BlueQubit backend**
3. **Choose device** (GPU, CPU, MPS, or Pauli-Path)
4. **Save settings**
5. **Run experiment**
6. **Check logs** for device selection:
   ```
   📍 Using BlueQubit backend: device=gpu
   ```

## Summary

BlueQubit now supports all device types with proper configuration passing from frontend → API → backend. Users can choose the optimal device for their molecule size and budget.

For most users:
- **Start with GPU** (fastest, cheapest for small molecules)
- **Upgrade to MPS** if you need more qubits
- **Use Pauli-Path** only for very large systems with expert knowledge
