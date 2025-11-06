"""
Test XYZ Trajectory Reading/Writing

Validates that XYZ trajectory files can be written and read back correctly.

CRITICAL FIX: Issue #6 - XYZ trajectory reading was not implemented
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from kanad.dynamics.trajectory import Trajectory, TrajectoryWriter
import numpy as np
import tempfile
import os


def test_xyz_write_read_roundtrip():
    """Test that XYZ write/read gives back the same data."""
    print("=" * 80)
    print("TEST: XYZ Write/Read Roundtrip")
    print("=" * 80)

    # Create a simple trajectory
    traj_original = Trajectory()
    traj_original.atom_symbols = ['H', 'H']
    traj_original.n_atoms = 2

    # Add 3 frames
    for i in range(3):
        time = i * 0.5  # fs
        positions = np.array([[0.0, 0.0, 0.0], [1.4 + i*0.1, 0.0, 0.0]])  # Bohr
        velocities = np.array([[0.0, 0.0, 0.0], [0.01, 0.0, 0.0]])  # Bohr/fs
        forces = np.array([[0.0, 0.0, 0.0], [-0.001, 0.0, 0.0]])  # Ha/Bohr

        kinetic_energy = 0.001 * i
        potential_energy = -1.1 + i * 0.01
        temperature = 300.0 + i * 10.0

        traj_original.add_frame(
            positions=positions,
            velocities=velocities,
            forces=forces,
            kinetic_energy=kinetic_energy,
            potential_energy=potential_energy,
            temperature=temperature,
            time=time
        )

    print(f"\nOriginal trajectory: {len(traj_original)} frames")
    print(f"Atom symbols: {traj_original.atom_symbols}")
    print(f"N atoms: {traj_original.n_atoms}")

    # Write to XYZ
    with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz', delete=False) as f:
        xyz_file = f.name

    try:
        writer = TrajectoryWriter(format='xyz')
        writer.write(traj_original, xyz_file)
        print(f"\n✅ Wrote to: {xyz_file}")

        # Read back
        traj_read = writer.read(xyz_file)
        print(f"✅ Read back: {len(traj_read)} frames")

        # Validate
        print("\n" + "=" * 80)
        print("VALIDATION")
        print("=" * 80)

        # Check number of frames
        assert len(traj_read) == len(traj_original), \
            f"Frame count mismatch: {len(traj_read)} vs {len(traj_original)}"
        print(f"✅ Frame count: {len(traj_read)} frames")

        # Check atom symbols
        assert traj_read.atom_symbols == traj_original.atom_symbols, \
            f"Atom symbols mismatch: {traj_read.atom_symbols} vs {traj_original.atom_symbols}"
        print(f"✅ Atom symbols: {traj_read.atom_symbols}")

        # Check n_atoms
        assert traj_read.n_atoms == traj_original.n_atoms, \
            f"N atoms mismatch: {traj_read.n_atoms} vs {traj_original.n_atoms}"
        print(f"✅ N atoms: {traj_read.n_atoms}")

        # Check each frame
        all_passed = True
        for i in range(len(traj_read)):
            frame_orig = traj_original[i]
            frame_read = traj_read[i]

            # Positions should match (within numerical precision)
            pos_diff = np.max(np.abs(frame_orig.positions - frame_read.positions))
            if pos_diff > 1e-6:
                print(f"❌ Frame {i}: Position mismatch (max diff: {pos_diff:.2e} Bohr)")
                all_passed = False
            else:
                print(f"✅ Frame {i}: Positions match (max diff: {pos_diff:.2e} Bohr)")

            # Time should match
            time_diff = abs(frame_orig.time - frame_read.time)
            if time_diff > 1e-6:
                print(f"❌ Frame {i}: Time mismatch ({frame_read.time} vs {frame_orig.time})")
                all_passed = False

            # Total energy should match
            energy_diff = abs(frame_orig.total_energy - frame_read.total_energy)
            if energy_diff > 1e-6:
                print(f"❌ Frame {i}: Energy mismatch ({frame_read.total_energy} vs {frame_orig.total_energy})")
                all_passed = False

            # Temperature should match
            temp_diff = abs(frame_orig.temperature - frame_read.temperature)
            if temp_diff > 0.1:
                print(f"❌ Frame {i}: Temperature mismatch ({frame_read.temperature} vs {frame_orig.temperature})")
                all_passed = False

        # Note: XYZ doesn't store velocities or forces
        print(f"\n⚠️  Note: XYZ format doesn't store velocities or forces (set to zero on read)")

        if all_passed:
            print("\n" + "=" * 80)
            print("✅ ALL VALIDATION CHECKS PASSED")
            print("=" * 80)
            return True
        else:
            print("\n" + "=" * 80)
            print("❌ SOME VALIDATION CHECKS FAILED")
            print("=" * 80)
            return False

    finally:
        # Clean up
        if os.path.exists(xyz_file):
            os.remove(xyz_file)
            print(f"\n🧹 Cleaned up: {xyz_file}")


def test_xyz_read_external():
    """Test reading an externally-created XYZ file."""
    print("\n\n")
    print("=" * 80)
    print("TEST: Read External XYZ File")
    print("=" * 80)

    # Create a simple XYZ file manually
    xyz_content = """2
time=0.00fs E=-1.110000Ha T=300.0K
H        0.000000    0.000000    0.000000
H        1.400000    0.000000    0.000000
2
time=0.50fs E=-1.100000Ha T=310.0K
H        0.000000    0.000000    0.000000
H        1.500000    0.000000    0.000000
"""

    with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz', delete=False) as f:
        f.write(xyz_content)
        xyz_file = f.name

    try:
        writer = TrajectoryWriter(format='xyz')
        traj = writer.read(xyz_file)

        print(f"\n✅ Read {len(traj)} frames from external XYZ file")
        print(f"Atom symbols: {traj.atom_symbols}")
        print(f"N atoms: {traj.n_atoms}")

        # Validate
        assert len(traj) == 2, f"Expected 2 frames, got {len(traj)}"
        assert traj.atom_symbols == ['H', 'H'], f"Expected ['H', 'H'], got {traj.atom_symbols}"
        assert traj.n_atoms == 2, f"Expected 2 atoms, got {traj.n_atoms}"

        # Check frame 0
        frame0 = traj[0]
        assert abs(frame0.time - 0.0) < 1e-6, f"Frame 0 time: expected 0.0, got {frame0.time}"
        assert abs(frame0.total_energy - (-1.110000)) < 1e-6, f"Frame 0 energy: expected -1.11, got {frame0.total_energy}"
        assert abs(frame0.temperature - 300.0) < 0.1, f"Frame 0 temp: expected 300.0, got {frame0.temperature}"

        # Check frame 1
        frame1 = traj[1]
        assert abs(frame1.time - 0.5) < 1e-6, f"Frame 1 time: expected 0.5, got {frame1.time}"
        assert abs(frame1.total_energy - (-1.100000)) < 1e-6, f"Frame 1 energy: expected -1.10, got {frame1.total_energy}"
        assert abs(frame1.temperature - 310.0) < 0.1, f"Frame 1 temp: expected 310.0, got {frame1.temperature}"

        print("\n" + "=" * 80)
        print("✅ EXTERNAL XYZ FILE READ SUCCESSFULLY")
        print("=" * 80)
        return True

    except Exception as e:
        print(f"\n❌ FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False

    finally:
        if os.path.exists(xyz_file):
            os.remove(xyz_file)


if __name__ == '__main__':
    print("\n🔬 XYZ TRAJECTORY READING/WRITING TEST\n")
    print("Testing Issue #6: XYZ trajectory reading implementation")
    print("CRITICAL: Was NotImplementedError, now implemented!\n")

    passed_roundtrip = test_xyz_write_read_roundtrip()
    passed_external = test_xyz_read_external()

    print("\n\n")
    print("=" * 80)
    print("FINAL SUMMARY")
    print("=" * 80)
    print(f"Roundtrip Test: {'✅ PASS' if passed_roundtrip else '❌ FAIL'}")
    print(f"External Test:  {'✅ PASS' if passed_external else '❌ FAIL'}")

    if passed_roundtrip and passed_external:
        print("\n✅ XYZ TRAJECTORY READING IMPLEMENTED!")
        print("   ✓ Can write XYZ files")
        print("   ✓ Can read XYZ files")
        print("   ✓ Positions preserved")
        print("   ✓ Time/energy/temperature preserved")
        print("   ✓ Atom symbols preserved")
        print("\n🎉 Issue #6 RESOLVED!")
        sys.exit(0)
    else:
        print("\n❌ SOME TESTS FAILED")
        sys.exit(1)
