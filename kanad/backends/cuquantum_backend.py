"""
cuQuantum backend for GPU-accelerated quantum simulation.

NVIDIA cuQuantum provides GPU-accelerated quantum circuit simulation
using CUDA, offering significant speedup for large-scale simulations.

Prerequisites:
- NVIDIA GPU with CUDA support
- cuQuantum SDK installed
- cuQuantum-Python bindings
"""

from typing import Optional, Dict, Any
import warnings
import logging

logger = logging.getLogger(__name__)


class CuQuantumBackend:
    """
    NVIDIA cuQuantum backend for GPU-accelerated simulation.

    Features:
    - GPU-accelerated statevector simulation
    - Tensor network simulation
    - Support for large qubit counts
    - Significant speedup over CPU simulation

    Supports two simulation methods:
    1. Statevector: Full statevector simulation on GPU
    2. Tensor Network: Memory-efficient for large circuits
    """

    def __init__(
        self,
        backend_name: str = 'cuquantum_statevector',
        simulation_method: str = 'statevector',
        gpu_device: int = 0,
        precision: str = 'single',
        **options
    ):
        """
        Initialize cuQuantum backend.

        Args:
            backend_name: Backend identifier
                - 'cuquantum_statevector': GPU statevector simulation
                - 'cuquantum_tensornet': Tensor network simulation
            simulation_method: Simulation method ('statevector' or 'tensornet')
            gpu_device: GPU device ID (default: 0)
            precision: Precision ('single' or 'double')
            **options: Additional cuQuantum-specific options
        """
        self.backend_name = backend_name
        self.simulation_method = simulation_method
        self.gpu_device = gpu_device
        self.precision = precision
        self.options = options

        # Initialize cuQuantum components
        self.backend = None
        self._estimator = None
        self._cuquantum_available = self._check_cuquantum()

        if self._cuquantum_available:
            self._initialize_backend()
        else:
            raise ImportError(
                "cuQuantum not available. Install with:\n"
                "  pip install cuquantum-python\n"
                "Requires: NVIDIA GPU + CUDA toolkit"
            )

    def _check_cuquantum(self) -> bool:
        """Check if cuQuantum is available."""
        try:
            import cuquantum
            import cupy as cp

            # Check GPU availability
            gpu_count = cp.cuda.runtime.getDeviceCount()
            if gpu_count == 0:
                warnings.warn("No CUDA GPUs detected. cuQuantum requires NVIDIA GPU.")
                return False

            # Check CUDA version
            cuda_version = cp.cuda.runtime.runtimeGetVersion()
            logger.info(f"CUDA runtime version: {cuda_version}")
            logger.info(f"Available GPUs: {gpu_count}")
            logger.info(f"cuQuantum version: {cuquantum.__version__}")

            return True

        except ImportError as e:
            warnings.warn(
                f"cuQuantum dependencies not found: {e}\n"
                "Install: pip install cuquantum-python cupy-cuda11x"
            )
            return False

    def _initialize_backend(self):
        """Initialize cuQuantum backend."""
        try:
            # Import cuQuantum components
            import cupy as cp
            from qiskit_aer import AerSimulator

            # Set GPU device
            cp.cuda.Device(self.gpu_device).use()
            logger.info(f"Using GPU device: {self.gpu_device}")

            # Create Aer simulator with GPU acceleration
            # Qiskit Aer 0.17+ automatically uses cuQuantum if available
            if self.simulation_method == 'statevector':
                self.backend = AerSimulator(
                    method='statevector',
                    device='GPU',
                    precision=self.precision,
                    **self.options
                )
                logger.info("Initialized cuQuantum statevector backend (GPU)")

            elif self.simulation_method == 'tensornet':
                self.backend = AerSimulator(
                    method='tensor_network',
                    device='GPU',
                    precision=self.precision,
                    **self.options
                )
                logger.info("Initialized cuQuantum tensor network backend (GPU)")

            else:
                raise ValueError(
                    f"Unknown simulation method: {self.simulation_method}. "
                    "Use 'statevector' or 'tensornet'"
                )

        except ImportError as e:
            raise ImportError(
                f"Failed to initialize cuQuantum backend: {e}\n"
                "Ensure qiskit-aer-gpu is installed: pip install qiskit-aer-gpu"
            )

    def get_estimator(self, **options):
        """
        Get Estimator primitive for cuQuantum backend.

        Uses Qiskit Aer's GPU-accelerated Estimator with cuStateVec.

        Args:
            **options: Estimator-specific options

        Returns:
            GPU-accelerated Estimator primitive
        """
        if self._estimator is not None:
            return self._estimator

        # Use Aer Estimator with GPU backend (compatible with Qiskit 2.2+)
        from qiskit_aer.primitives import Estimator

        # Pass GPU backend options to Estimator
        # In Qiskit Aer, Estimator uses backend_options to configure the internal AerSimulator
        backend_opts = {
            'method': 'statevector',
            'device': 'GPU',
            'cuStateVec_enable': True,
            'precision': 'double',
            'blocking_enable': True,
            'blocking_qubits': 20
        }
        backend_opts.update(options)

        self._estimator = Estimator(backend_options=backend_opts)

        logger.info(f"Created cuQuantum Estimator (GPU-accelerated on device {self.gpu_device})")

        return self._estimator

    def get_sampler(self, **options):
        """
        Get Sampler primitive for cuQuantum backend.

        Args:
            **options: Sampler-specific options

        Returns:
            GPU-accelerated Sampler primitive
        """
        from qiskit_aer.primitives import Sampler

        return Sampler()

    def benchmark_gpu_speedup(self, circuit, n_runs: int = 10):
        """
        Benchmark GPU speedup vs CPU.

        Args:
            circuit: Qiskit QuantumCircuit
            n_runs: Number of benchmark runs

        Returns:
            Dictionary with timing results
        """
        import time
        import numpy as np
        from qiskit_aer import AerSimulator

        # CPU simulation
        cpu_backend = AerSimulator(method='statevector', device='CPU')
        cpu_times = []

        print(f"\nBenchmarking {circuit.num_qubits} qubits, {circuit.depth()} depth...")
        print("Running CPU simulation...")

        for i in range(n_runs):
            start = time.time()
            cpu_backend.run(circuit, shots=1).result()
            cpu_times.append(time.time() - start)

        cpu_avg = np.mean(cpu_times)

        # GPU simulation
        gpu_times = []
        print("Running GPU simulation...")

        for i in range(n_runs):
            start = time.time()
            self.backend.run(circuit, shots=1).result()
            gpu_times.append(time.time() - start)

        gpu_avg = np.mean(gpu_times)
        speedup = cpu_avg / gpu_avg

        results = {
            'cpu_time': cpu_avg,
            'gpu_time': gpu_avg,
            'speedup': speedup,
            'qubits': circuit.num_qubits,
            'depth': circuit.depth()
        }

        print(f"\nBenchmark Results:")
        print(f"  CPU time: {cpu_avg*1000:.2f} ms")
        print(f"  GPU time: {gpu_avg*1000:.2f} ms")
        print(f"  Speedup:  {speedup:.2f}x")

        return results

    def get_backend_info(self) -> Dict[str, Any]:
        """
        Get cuQuantum backend information.

        Returns:
            Dictionary with backend properties
        """
        import cupy as cp

        device = cp.cuda.Device(self.gpu_device)
        device_props = device.attributes

        info = {
            'name': self.backend_name,
            'simulation_method': self.simulation_method,
            'gpu_device': self.gpu_device,
            'gpu_name': device.name,
            'gpu_memory_total': device_props.get('TotalMemory', 'Unknown'),
            'precision': self.precision,
        }

        # Add CUDA info
        try:
            cuda_version = cp.cuda.runtime.runtimeGetVersion()
            info['cuda_version'] = cuda_version
        except:
            pass

        # Add cuQuantum info
        try:
            import cuquantum
            info['cuquantum_version'] = cuquantum.__version__
        except:
            pass

        return info

    def estimate_max_qubits(self) -> int:
        """
        Estimate maximum number of qubits for current GPU.

        Based on available GPU memory.

        Returns:
            Estimated maximum qubits
        """
        import cupy as cp

        device = cp.cuda.Device(self.gpu_device)
        total_memory = device.mem_info[1]  # Total memory in bytes

        # Statevector size: 2^n * (8 bytes for single, 16 for double)
        bytes_per_amplitude = 8 if self.precision == 'single' else 16

        # Conservative estimate (use 80% of memory)
        usable_memory = total_memory * 0.8

        # 2^n * bytes_per_amplitude <= usable_memory
        # n <= log2(usable_memory / bytes_per_amplitude)
        import math
        max_qubits = int(math.log2(usable_memory / bytes_per_amplitude))

        print(f"\nGPU Memory Analysis:")
        print(f"  Total GPU memory: {total_memory / 1e9:.2f} GB")
        print(f"  Usable memory: {usable_memory / 1e9:.2f} GB")
        print(f"  Precision: {self.precision}")
        print(f"  Estimated max qubits: {max_qubits}")

        return max_qubits

    def __repr__(self) -> str:
        """String representation."""
        return (f"CuQuantumBackend(method='{self.simulation_method}', "
                f"device={self.gpu_device}, precision='{self.precision}')")


def check_cuquantum_available() -> bool:
    """
    Check if cuQuantum is available and properly configured.

    Returns:
        True if cuQuantum can be used, False otherwise
    """
    try:
        import cuquantum
        import cupy as cp

        gpu_count = cp.cuda.runtime.getDeviceCount()
        if gpu_count == 0:
            print("❌ No CUDA GPUs detected")
            return False

        print(f"✅ cuQuantum available")
        print(f"   Version: {cuquantum.__version__}")
        print(f"   GPUs: {gpu_count}")

        return True

    except ImportError as e:
        print(f"❌ cuQuantum not available: {e}")
        print("   Install: pip install cuquantum-python cupy-cuda11x")
        return False


def get_gpu_info():
    """
    Get information about available GPUs.

    Returns:
        List of GPU information dictionaries
    """
    try:
        import cupy as cp

        gpu_count = cp.cuda.runtime.getDeviceCount()
        gpu_info = []

        for i in range(gpu_count):
            device = cp.cuda.Device(i)
            mem_info = device.mem_info

            # Get GPU properties using runtime API
            props = cp.cuda.runtime.getDeviceProperties(i)

            info = {
                'id': i,
                'name': props['name'].decode('utf-8') if isinstance(props['name'], bytes) else str(props['name']),
                'compute_capability': f"{props['major']}.{props['minor']}",
                'total_memory_gb': mem_info[1] / 1e9,
                'free_memory_gb': mem_info[0] / 1e9,
            }
            gpu_info.append(info)

        return gpu_info

    except ImportError:
        return []
