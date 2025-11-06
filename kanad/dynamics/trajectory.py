"""
Trajectory Storage and Management for MD Simulations

Handles storage, retrieval, and analysis of molecular dynamics trajectories.

A trajectory is a sequence of molecular configurations (frames) sampled during
an MD simulation. Each frame contains:
- Atomic positions
- Velocities
- Forces
- Energies (kinetic, potential, total)
- Temperature
- Time

Storage Formats:
---------------
- **HDF5**: Efficient binary format for large trajectories (recommended)
- **XYZ**: Simple text format for visualization
- **NumPy**: In-memory arrays for analysis

Features:
--------
- Compressed HDF5 storage (saves disk space)
- Chunked I/O (memory-efficient for large trajectories)
- Frame extraction and slicing
- Property time series extraction
- Trajectory concatenation
- Memory-efficient streaming

Example Usage:
-------------
```python
from kanad.dynamics import Trajectory, TrajectoryWriter

# Create trajectory
traj = Trajectory()

# Add frames during MD
for step in range(n_steps):
    traj.add_frame(positions, velocities, forces, energy, temperature, time)

# Save to file
writer = TrajectoryWriter('trajectory.h5', format='hdf5')
writer.write(traj)

# Load and analyze
traj2 = writer.read('trajectory.h5')
energies = traj2.get_property('potential_energy')
```

References:
----------
- HDF5: https://www.hdfgroup.org/solutions/hdf5/
- MDAnalysis: Michaud-Agrawal et al. (2011) J. Comput. Chem. 32, 2319
- MDTraj: McGibbon et al. (2015) Biophys. J. 109, 1528
"""

import numpy as np
import h5py
import logging
from typing import List, Optional, Union, Dict, Any
from dataclasses import dataclass, field
from pathlib import Path

logger = logging.getLogger(__name__)


@dataclass
class TrajectoryFrame:
    """
    Single snapshot from MD trajectory.

    Attributes:
        positions: Atomic positions (N_atoms, 3) in Bohr
        velocities: Atomic velocities (N_atoms, 3) in Bohr/fs
        forces: Forces on atoms (N_atoms, 3) in Ha/Bohr
        kinetic_energy: Kinetic energy in Hartree
        potential_energy: Potential energy in Hartree
        total_energy: Total energy in Hartree
        temperature: Instantaneous temperature in K
        time: Simulation time in fs
        box: Simulation box dimensions (optional, for PBC)
        metadata: Additional properties (dict)
    """
    positions: np.ndarray
    velocities: np.ndarray
    forces: np.ndarray
    kinetic_energy: float
    potential_energy: float
    total_energy: float
    temperature: float
    time: float
    box: Optional[np.ndarray] = None
    metadata: Dict[str, Any] = field(default_factory=dict)


class Trajectory:
    """
    Container for MD trajectory data.

    Stores a sequence of TrajectoryFrame objects with efficient access methods.
    """

    def __init__(self):
        """Initialize empty trajectory."""
        self.frames: List[TrajectoryFrame] = []
        self.n_atoms: Optional[int] = None
        self.atom_symbols: Optional[List[str]] = None
        self.atom_masses: Optional[np.ndarray] = None

        logger.debug("Initialized empty Trajectory")

    def add_frame(
        self,
        positions: np.ndarray,
        velocities: np.ndarray,
        forces: np.ndarray,
        kinetic_energy: float,
        potential_energy: float,
        temperature: float,
        time: float,
        **kwargs
    ):
        """
        Add a new frame to trajectory.

        Args:
            positions: (N_atoms, 3) in Bohr
            velocities: (N_atoms, 3) in Bohr/fs
            forces: (N_atoms, 3) in Ha/Bohr
            kinetic_energy: KE in Hartree
            potential_energy: PE in Hartree
            temperature: T in K
            time: Time in fs
            **kwargs: Additional metadata
        """
        # Validate n_atoms consistency
        if self.n_atoms is None:
            self.n_atoms = len(positions)
        elif len(positions) != self.n_atoms:
            raise ValueError(f"Position array has {len(positions)} atoms, expected {self.n_atoms}")

        total_energy = kinetic_energy + potential_energy

        frame = TrajectoryFrame(
            positions=positions.copy(),
            velocities=velocities.copy(),
            forces=forces.copy(),
            kinetic_energy=kinetic_energy,
            potential_energy=potential_energy,
            total_energy=total_energy,
            temperature=temperature,
            time=time,
            metadata=kwargs
        )

        self.frames.append(frame)

    def __len__(self) -> int:
        """Return number of frames."""
        return len(self.frames)

    def __getitem__(self, idx: Union[int, slice]) -> Union[TrajectoryFrame, List[TrajectoryFrame]]:
        """Get frame(s) by index or slice."""
        return self.frames[idx]

    def get_positions(self, frame_indices: Optional[List[int]] = None) -> np.ndarray:
        """
        Extract positions for specified frames.

        Args:
            frame_indices: Frame indices to extract (None = all)

        Returns:
            Array of shape (n_frames, n_atoms, 3)
        """
        if frame_indices is None:
            frame_indices = range(len(self.frames))

        positions = np.array([self.frames[i].positions for i in frame_indices])
        return positions

    def get_velocities(self, frame_indices: Optional[List[int]] = None) -> np.ndarray:
        """Extract velocities (n_frames, n_atoms, 3)."""
        if frame_indices is None:
            frame_indices = range(len(self.frames))

        velocities = np.array([self.frames[i].velocities for i in frame_indices])
        return velocities

    def get_forces(self, frame_indices: Optional[List[int]] = None) -> np.ndarray:
        """Extract forces (n_frames, n_atoms, 3)."""
        if frame_indices is None:
            frame_indices = range(len(self.frames))

        forces = np.array([self.frames[i].forces for i in frame_indices])
        return forces

    def get_property(self, property_name: str) -> np.ndarray:
        """
        Extract time series of a property.

        Args:
            property_name: Property name ('kinetic_energy', 'potential_energy',
                          'total_energy', 'temperature', 'time')

        Returns:
            Array of property values (n_frames,)
        """
        valid_properties = [
            'kinetic_energy', 'potential_energy', 'total_energy',
            'temperature', 'time'
        ]

        if property_name not in valid_properties:
            raise ValueError(f"Unknown property '{property_name}'. Valid: {valid_properties}")

        values = np.array([getattr(frame, property_name) for frame in self.frames])
        return values

    def get_time_series(self) -> Dict[str, np.ndarray]:
        """
        Extract all time series as dictionary.

        Returns:
            Dictionary with keys: time, kinetic_energy, potential_energy,
            total_energy, temperature
        """
        return {
            'time': self.get_property('time'),
            'kinetic_energy': self.get_property('kinetic_energy'),
            'potential_energy': self.get_property('potential_energy'),
            'total_energy': self.get_property('total_energy'),
            'temperature': self.get_property('temperature'),
        }

    def slice(self, start: int = 0, stop: Optional[int] = None, step: int = 1) -> 'Trajectory':
        """
        Create new trajectory from frame slice.

        Args:
            start: Start frame
            stop: Stop frame (None = end)
            step: Step size

        Returns:
            New Trajectory with sliced frames
        """
        new_traj = Trajectory()
        new_traj.n_atoms = self.n_atoms
        new_traj.atom_symbols = self.atom_symbols
        new_traj.atom_masses = self.atom_masses
        new_traj.frames = self.frames[start:stop:step]

        logger.debug(f"Sliced trajectory: {len(self)} → {len(new_traj)} frames")
        return new_traj

    def concatenate(self, other: 'Trajectory'):
        """
        Concatenate another trajectory to this one.

        Args:
            other: Trajectory to append
        """
        if self.n_atoms != other.n_atoms:
            raise ValueError(f"Cannot concatenate: different n_atoms ({self.n_atoms} vs {other.n_atoms})")

        self.frames.extend(other.frames)
        logger.debug(f"Concatenated trajectories: {len(other)} frames added")

    def compute_statistics(self) -> Dict[str, Any]:
        """
        Compute statistical properties over trajectory.

        Returns:
            Dictionary with mean, std, min, max for energies and temperature
        """
        time_series = self.get_time_series()

        stats = {}
        for key in ['kinetic_energy', 'potential_energy', 'total_energy', 'temperature']:
            values = time_series[key]
            stats[key] = {
                'mean': np.mean(values),
                'std': np.std(values),
                'min': np.min(values),
                'max': np.max(values)
            }

        return stats

    # Convenient property accessors
    @property
    def times(self) -> np.ndarray:
        """Get time series."""
        return self.get_property('time')

    @property
    def kinetic_energies(self) -> np.ndarray:
        """Get kinetic energy series."""
        return self.get_property('kinetic_energy')

    @property
    def potential_energies(self) -> np.ndarray:
        """Get potential energy series."""
        return self.get_property('potential_energy')

    @property
    def total_energies(self) -> np.ndarray:
        """Get total energy series."""
        return self.get_property('total_energy')

    @property
    def temperatures(self) -> np.ndarray:
        """Get temperature series."""
        return self.get_property('temperature')


class TrajectoryWriter:
    """
    Write trajectories to file (HDF5 or XYZ format).

    HDF5 Format:
    -----------
    - Efficient binary storage
    - Compression (gzip, level 4)
    - Chunked storage (memory-efficient)
    - Stores all properties

    XYZ Format:
    ----------
    - Simple text format
    - Positions only (for visualization)
    - Compatible with VMD, PyMOL, etc.
    """

    def __init__(self, format: str = 'hdf5'):
        """
        Initialize writer.

        Args:
            format: Output format ('hdf5' or 'xyz')
        """
        self.format = format.lower()
        if self.format not in ['hdf5', 'xyz']:
            raise ValueError(f"Unknown format '{format}'. Use 'hdf5' or 'xyz'")

        logger.debug(f"Initialized TrajectoryWriter (format={self.format})")

    def write(self, trajectory: Trajectory, filename: Union[str, Path]):
        """
        Write trajectory to file.

        Args:
            trajectory: Trajectory object
            filename: Output filename
        """
        filename = Path(filename)

        if self.format == 'hdf5':
            self._write_hdf5(trajectory, filename)
        elif self.format == 'xyz':
            self._write_xyz(trajectory, filename)

        logger.info(f"Wrote {len(trajectory)} frames to {filename}")

    def _write_hdf5(self, trajectory: Trajectory, filename: Path):
        """Write to HDF5 format with compression."""
        n_frames = len(trajectory)
        n_atoms = trajectory.n_atoms

        with h5py.File(filename, 'w') as f:
            # Create datasets with compression
            # Chunk size: (1 frame, all atoms, 3 coords) for efficient frame access
            chunk_shape = (1, n_atoms, 3)

            positions = f.create_dataset(
                'positions',
                shape=(n_frames, n_atoms, 3),
                dtype='f8',
                compression='gzip',
                compression_opts=4,
                chunks=chunk_shape
            )

            velocities = f.create_dataset(
                'velocities',
                shape=(n_frames, n_atoms, 3),
                dtype='f8',
                compression='gzip',
                compression_opts=4,
                chunks=chunk_shape
            )

            forces = f.create_dataset(
                'forces',
                shape=(n_frames, n_atoms, 3),
                dtype='f8',
                compression='gzip',
                compression_opts=4,
                chunks=chunk_shape
            )

            # Scalar properties (no chunking needed)
            time = f.create_dataset('time', shape=(n_frames,), dtype='f8')
            kinetic = f.create_dataset('kinetic_energy', shape=(n_frames,), dtype='f8')
            potential = f.create_dataset('potential_energy', shape=(n_frames,), dtype='f8')
            total = f.create_dataset('total_energy', shape=(n_frames,), dtype='f8')
            temperature = f.create_dataset('temperature', shape=(n_frames,), dtype='f8')

            # Write data
            for i, frame in enumerate(trajectory.frames):
                positions[i] = frame.positions
                velocities[i] = frame.velocities
                forces[i] = frame.forces
                time[i] = frame.time
                kinetic[i] = frame.kinetic_energy
                potential[i] = frame.potential_energy
                total[i] = frame.total_energy
                temperature[i] = frame.temperature

            # Metadata
            f.attrs['n_frames'] = n_frames
            f.attrs['n_atoms'] = n_atoms
            f.attrs['format'] = 'kanad_md_trajectory'
            f.attrs['version'] = '1.0'

            if trajectory.atom_symbols is not None:
                f.attrs['atom_symbols'] = trajectory.atom_symbols
            if trajectory.atom_masses is not None:
                f.create_dataset('atom_masses', data=trajectory.atom_masses)

    def _write_xyz(self, trajectory: Trajectory, filename: Path):
        """Write to XYZ format (positions only, for visualization)."""
        n_atoms = trajectory.n_atoms

        # Convert Bohr to Angstrom for visualization
        BOHR_TO_ANGSTROM = 0.529177

        with open(filename, 'w') as f:
            for frame in trajectory.frames:
                # XYZ format:
                # Line 1: N_atoms
                # Line 2: Comment (time, energy)
                # Lines 3+: element x y z

                f.write(f"{n_atoms}\n")
                comment = f"time={frame.time:.2f}fs E={frame.total_energy:.6f}Ha T={frame.temperature:.1f}K"
                f.write(f"{comment}\n")

                # Atom lines
                positions_angstrom = frame.positions * BOHR_TO_ANGSTROM

                for i in range(n_atoms):
                    # Use atom symbols if available, otherwise use 'X'
                    if trajectory.atom_symbols is not None:
                        symbol = trajectory.atom_symbols[i]
                    else:
                        symbol = 'X'

                    x, y, z = positions_angstrom[i]
                    f.write(f"{symbol:2s} {x:12.6f} {y:12.6f} {z:12.6f}\n")

    def read(self, filename: Union[str, Path]) -> Trajectory:
        """
        Read trajectory from file.

        Args:
            filename: Input filename

        Returns:
            Trajectory object
        """
        filename = Path(filename)

        if filename.suffix == '.h5' or filename.suffix == '.hdf5':
            traj = self._read_hdf5(filename)
        elif filename.suffix == '.xyz':
            raise NotImplementedError("XYZ reading not yet implemented (write-only for visualization)")
        else:
            raise ValueError(f"Unknown file extension: {filename.suffix}")

        logger.info(f"Read {len(traj)} frames from {filename}")
        return traj

    def _read_hdf5(self, filename: Path) -> Trajectory:
        """Read from HDF5 format."""
        trajectory = Trajectory()

        with h5py.File(filename, 'r') as f:
            n_frames = f.attrs['n_frames']
            n_atoms = f.attrs['n_atoms']

            trajectory.n_atoms = n_atoms

            # Read atom metadata if present
            if 'atom_symbols' in f.attrs:
                trajectory.atom_symbols = list(f.attrs['atom_symbols'])
            if 'atom_masses' in f:
                trajectory.atom_masses = f['atom_masses'][:]

            # Read frames
            for i in range(n_frames):
                frame = TrajectoryFrame(
                    positions=f['positions'][i],
                    velocities=f['velocities'][i],
                    forces=f['forces'][i],
                    kinetic_energy=float(f['kinetic_energy'][i]),
                    potential_energy=float(f['potential_energy'][i]),
                    total_energy=float(f['total_energy'][i]),
                    temperature=float(f['temperature'][i]),
                    time=float(f['time'][i])
                )
                trajectory.frames.append(frame)

        return trajectory


class TrajectoryAnalyzer:
    """
    Analyze MD trajectories.

    Provides methods for computing structural and dynamical properties.
    """

    def __init__(self, trajectory: Trajectory):
        """
        Initialize analyzer.

        Args:
            trajectory: Trajectory to analyze
        """
        self.trajectory = trajectory
        logger.debug(f"Initialized TrajectoryAnalyzer ({len(trajectory)} frames)")

    def compute_rdf(
        self,
        atom_i: int,
        atom_j: int,
        r_max: float = 10.0,
        n_bins: int = 100
    ) -> tuple:
        """
        Compute radial distribution function g(r) between atoms i and j.

        RDF measures the probability of finding atom j at distance r from atom i.

        Args:
            atom_i: Index of first atom
            atom_j: Index of second atom
            r_max: Maximum distance in Bohr
            n_bins: Number of histogram bins

        Returns:
            (r, g_r): Distance bins and RDF values
        """
        # Distance bins
        r_bins = np.linspace(0, r_max, n_bins + 1)
        r_centers = 0.5 * (r_bins[1:] + r_bins[:-1])

        # Compute distances across all frames
        distances = []
        for frame in self.trajectory.frames:
            pos_i = frame.positions[atom_i]
            pos_j = frame.positions[atom_j]
            r = np.linalg.norm(pos_j - pos_i)
            distances.append(r)

        # Histogram
        hist, _ = np.histogram(distances, bins=r_bins)

        # Normalize to get g(r)
        # For single pair: g(r) = hist / (density * volume_shell * n_frames)
        # For diatomic: just normalized histogram
        n_frames = len(self.trajectory)
        g_r = hist / n_frames

        # Normalize to 1 at large r
        if np.max(g_r) > 0:
            g_r = g_r / np.mean(g_r[-10:]) if np.mean(g_r[-10:]) > 0 else g_r

        return r_centers, g_r

    def compute_msd(self, atom_indices: Optional[List[int]] = None) -> tuple:
        """
        Compute mean squared displacement (MSD) for diffusion.

        MSD(t) = ⟨|r(t) - r(0)|²⟩

        Args:
            atom_indices: Atoms to include (None = all)

        Returns:
            (times, msd): Time points and MSD values
        """
        if atom_indices is None:
            atom_indices = range(self.trajectory.n_atoms)

        n_frames = len(self.trajectory)
        times = np.array([frame.time for frame in self.trajectory.frames])

        # Reference positions (frame 0)
        r0 = self.trajectory.frames[0].positions[atom_indices]

        # MSD at each time
        msd = np.zeros(n_frames)
        for i in range(n_frames):
            rt = self.trajectory.frames[i].positions[atom_indices]
            displacement = rt - r0
            msd[i] = np.mean(np.sum(displacement**2, axis=1))

        return times, msd

    def compute_vacf(self, atom_indices: Optional[List[int]] = None) -> tuple:
        """
        Compute velocity autocorrelation function (VACF).

        VACF(t) = ⟨v(t) · v(0)⟩

        Args:
            atom_indices: Atoms to include (None = all)

        Returns:
            (times, vacf): Time points and VACF values
        """
        if atom_indices is None:
            atom_indices = range(self.trajectory.n_atoms)

        n_frames = len(self.trajectory)
        times = np.array([frame.time for frame in self.trajectory.frames])

        # Reference velocities (frame 0)
        v0 = self.trajectory.frames[0].velocities[atom_indices]

        # VACF at each time
        vacf = np.zeros(n_frames)
        for i in range(n_frames):
            vt = self.trajectory.frames[i].velocities[atom_indices]
            vacf[i] = np.mean(np.sum(vt * v0, axis=1))

        # Normalize to VACF(0) = 1
        if vacf[0] > 0:
            vacf = vacf / vacf[0]

        return times, vacf
