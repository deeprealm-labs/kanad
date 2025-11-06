"""
Kanad Molecular Dynamics Module

🌟 WORLD'S FIRST: Quantum-Enhanced MD with Governance Protocols 🌟

This module provides molecular dynamics simulation capabilities with unique features:
- Quantum forces from VQE/SQD (correlation effects beyond HF)
- Governance-aware dynamics (bond-type-specific evolution)
- Environment integration (temperature, solvent, pH, pressure)
- Hybrid quantum-classical propagation

Components:
-----------
- md_simulator: Main MD simulation engine
- integrators: Time evolution algorithms (Velocity Verlet, Leapfrog, etc.)
- thermostats: Temperature control (Berendsen, Nose-Hoover, Langevin)
- barostats: Pressure control (NPT ensemble)
- trajectory: Trajectory storage and management
- initialization: Initial conditions (Maxwell-Boltzmann velocities)
- quantum_md: VQE/SQD-driven dynamics
- governance_md: Bond-aware MD constraints
- solvated_md: Implicit/explicit solvent MD

Example Usage:
-------------
```python
from kanad.bonds import BondFactory
from kanad.dynamics import MDSimulator

# Create molecule
bond = BondFactory.create_bond('H', 'H', distance=0.74)

# Run classical MD
md = MDSimulator(
    bond,
    temperature=300.0,  # K
    timestep=0.5,  # fs
    integrator='velocity_verlet',
    thermostat='berendsen'
)

trajectory = md.run(n_steps=1000)

# Or use quantum forces
md_quantum = MDSimulator(
    bond,
    temperature=300.0,
    force_method='vqe',  # Use quantum correlation!
    use_governance=True   # Bond-aware constraints
)

trajectory = md_quantum.run(n_steps=100)
```

References:
----------
- Velocity Verlet: Swope et al. (1982) J. Chem. Phys. 76, 637
- Berendsen thermostat: Berendsen et al. (1984) J. Chem. Phys. 81, 3684
- Nose-Hoover: Hoover (1985) Phys. Rev. A 31, 1695
- Ab initio MD: Car & Parrinello (1985) Phys. Rev. Lett. 55, 2471
"""

# Core MD components
from kanad.dynamics.integrators import (
    VelocityVerletIntegrator,
    LeapfrogIntegrator,
    RungeKuttaIntegrator
)

from kanad.dynamics.thermostats import (
    BerendsenThermostat,
    NoseHooverThermostat,
    LangevinThermostat,
    VelocityRescaling
)

from kanad.dynamics.trajectory import (
    Trajectory,
    TrajectoryFrame,
    TrajectoryWriter
)

from kanad.dynamics.initialization import (
    MaxwellBoltzmannInitializer,
    remove_com_motion,
    equilibrate_system
)

from kanad.dynamics.md_simulator import (
    MDSimulator,
    MDResult
)

__all__ = [
    # Main interface
    'MDSimulator',
    'MDResult',

    # Integrators
    'VelocityVerletIntegrator',
    'LeapfrogIntegrator',
    'RungeKuttaIntegrator',

    # Thermostats
    'BerendsenThermostat',
    'NoseHooverThermostat',
    'LangevinThermostat',
    'VelocityRescaling',

    # Trajectory
    'Trajectory',
    'TrajectoryFrame',
    'TrajectoryWriter',

    # Initialization
    'MaxwellBoltzmannInitializer',
    'remove_com_motion',
    'equilibrate_system',
]

__version__ = '1.0.0'
__author__ = 'Kanad Team'
__description__ = "Quantum-Enhanced Molecular Dynamics with Governance Protocols"
