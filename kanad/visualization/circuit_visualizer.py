"""
Circuit visualization module for generating text-based circuit diagrams.

Provides ASCII art representation of quantum circuits for preview.
"""

from typing import List, Dict, Any, Optional
import numpy as np
import base64
import io
import tempfile
import os


class CircuitVisualizer:
    """Text-based quantum circuit visualizer."""

    def __init__(self, n_qubits: int):
        """
        Initialize visualizer.

        Args:
            n_qubits: Number of qubits in circuit
        """
        self.n_qubits = n_qubits
        self.circuit_lines = [[] for _ in range(n_qubits)]
        self.max_depth = 0

    def add_gate(self, gate_type: str, qubits: List[int], params: List[Any] = None):
        """
        Add gate to visualization.

        Args:
            gate_type: Type of gate (h, x, y, z, rx, ry, rz, cx, etc.)
            qubits: List of qubit indices
            params: Gate parameters (if any)
        """
        gate_type = gate_type.lower()

        if len(qubits) == 1:
            # Single-qubit gate
            self._add_single_qubit_gate(gate_type, qubits[0], params)
        elif len(qubits) == 2:
            # Two-qubit gate
            self._add_two_qubit_gate(gate_type, qubits[0], qubits[1])

        self.max_depth += 1

    def _add_single_qubit_gate(self, gate_type: str, qubit: int, params: List[Any] = None):
        """Add single-qubit gate."""
        # Map gate types to symbols
        gate_symbols = {
            'h': 'H',
            'x': 'X',
            'y': 'Y',
            'z': 'Z',
            'rx': 'Rx',
            'ry': 'Ry',
            'rz': 'Rz',
            's': 'S',
            't': 'T',
            'sdg': 'S†',
            'tdg': 'T†'
        }

        symbol = gate_symbols.get(gate_type, gate_type.upper())

        # Add parameter notation if present
        if params and len(params) > 0:
            param = params[0]
            if hasattr(param, 'name'):
                symbol = f"{symbol}(θ)"
            elif isinstance(param, (int, float)):
                symbol = f"{symbol}({param:.2f})"

        # Add gate to circuit line
        self.circuit_lines[qubit].append(f"[{symbol}]")

        # Add spacers to other qubits
        for i in range(self.n_qubits):
            if i != qubit:
                self.circuit_lines[i].append("─────")

    def _add_two_qubit_gate(self, gate_type: str, control: int, target: int):
        """Add two-qubit gate (CNOT, CZ, etc.)."""
        min_q = min(control, target)
        max_q = max(control, target)

        for i in range(self.n_qubits):
            if i == control:
                if gate_type in ['cx', 'cnot']:
                    self.circuit_lines[i].append("──●──")
                else:
                    self.circuit_lines[i].append("──□──")
            elif i == target:
                if gate_type in ['cx', 'cnot']:
                    self.circuit_lines[i].append("──⊕──")
                elif gate_type == 'cz':
                    self.circuit_lines[i].append("──●──")
                else:
                    self.circuit_lines[i].append("──■──")
            elif min_q < i < max_q:
                # Connection line between control and target
                self.circuit_lines[i].append("──│──")
            else:
                self.circuit_lines[i].append("─────")

    def to_string(self) -> str:
        """
        Convert circuit to ASCII string representation.

        Returns:
            ASCII art circuit diagram
        """
        lines = []
        for i in range(self.n_qubits):
            # Start with qubit label
            qubit_line = f"q{i}: "
            # Add all gates
            qubit_line += "".join(self.circuit_lines[i])
            lines.append(qubit_line)

        return "\n".join(lines)

    def get_stats(self) -> Dict[str, Any]:
        """
        Get circuit statistics.

        Returns:
            Dictionary with circuit stats
        """
        # Count gate types
        gate_counts = {}
        total_gates = 0

        for qubit_gates in self.circuit_lines:
            for gate_str in qubit_gates:
                if gate_str != "─────" and gate_str != "──│──":
                    total_gates += 1
                    # Extract gate type
                    gate_type = gate_str.strip('[').strip(']').split('(')[0]
                    gate_counts[gate_type] = gate_counts.get(gate_type, 0) + 1

        return {
            'n_qubits': self.n_qubits,
            'depth': self.max_depth,
            'total_gates': total_gates,
            'gate_counts': gate_counts
        }


def generate_matplotlib_circuit(circuit) -> Optional[str]:
    """
    Generate circuit diagram using Qiskit's native Matplotlib drawer.

    IMPORTANT: This function expects a REAL Qiskit QuantumCircuit object,
    not the custom internal circuit representation. Use this for generating
    production-quality circuit diagrams.

    Args:
        circuit: Qiskit QuantumCircuit instance (from ansatz.circuit.to_qiskit())

    Returns:
        Base64 encoded PNG image string, or None if generation fails
    """
    try:
        import matplotlib
        matplotlib.use('Agg')  # Non-interactive backend
        import matplotlib.pyplot as plt
        from qiskit.circuit import QuantumCircuit
    except ImportError as e:
        print(f"❌ Import failed: {e}")
        return None

    try:
        # Check if circuit is a Qiskit QuantumCircuit
        if not isinstance(circuit, QuantumCircuit):
            # Try to convert if it's our internal circuit representation
            if hasattr(circuit, 'to_qiskit'):
                circuit = circuit.to_qiskit()
            else:
                print(f"❌ Circuit is not a Qiskit QuantumCircuit: {type(circuit)}")
                return None

        # Generate circuit diagram using Qiskit's native drawer
        # Use 'mpl' style for high-quality Matplotlib output
        fig = circuit.draw('mpl', style='iqp', fold=-1, scale=0.8)

        # Convert to base64 PNG
        buf = io.BytesIO()
        fig.savefig(buf, format='png', dpi=150, bbox_inches='tight',
                    facecolor='white', edgecolor='none')
        plt.close(fig)

        buf.seek(0)
        img_bytes = buf.read()

        # Validate PNG signature
        if len(img_bytes) < 8:
            print("❌ PNG generation produced invalid data")
            return None

        png_signature = img_bytes[:8]
        expected_signature = b'\x89PNG\r\n\x1a\n'
        if png_signature != expected_signature:
            print("❌ PNG signature validation failed")
            return None

        img_base64 = base64.b64encode(img_bytes).decode('utf-8')
        print(f"✅ Qiskit circuit image generated: {len(img_base64)} bytes (base64)")

        return img_base64

    except Exception as e:
        print(f"❌ Qiskit circuit generation failed: {e}")
        import traceback
        traceback.print_exc()
        return None


def visualize_from_circuit(circuit) -> tuple[str, Dict[str, Any], List[Dict[str, Any]]]:
    """
    Create visualization from QuantumCircuit object.

    Args:
        circuit: QuantumCircuit instance

    Returns:
        Tuple of (circuit_diagram_string, stats_dict, structured_gates_list)
    """
    visualizer = CircuitVisualizer(circuit.n_qubits)

    # Collect structured gate data for frontend rendering
    structured_gates = []

    # Add all gates from circuit
    for gate in circuit.gates:
        # Store structured gate data
        structured_gates.append({
            'type': gate['type'],
            'qubits': gate['qubits'],
            'params': gate.get('params', [])
        })

        visualizer.add_gate(
            gate['type'],
            gate['qubits'],
            gate.get('params', [])
        )

    return visualizer.to_string(), visualizer.get_stats(), structured_gates


def create_preview_circuit(molecule_config: Dict[str, Any],
                           backend_config: Dict[str, Any]) -> Dict[str, Any]:
    """
    Create a preview circuit based on configuration.

    Args:
        molecule_config: Molecule configuration (atoms, basis, etc.)
        backend_config: Backend configuration (method, ansatz, mapper, etc.)

    Returns:
        Dictionary with circuit preview information
    """
    from kanad.core.molecule import Molecule
    from kanad.core.atom import Atom
    from kanad.ansatze import (
        HardwareEfficientAnsatz,
        RealAmplitudesAnsatz,
        EfficientSU2Ansatz,
        UCCAnsatz,
        TwoLocalAnsatz,
        IonicGovernanceAnsatz,
        CovalentGovernanceAnsatz,
        AdaptiveGovernanceAnsatz
    )

    # Create molecule - handle both SMILES and atoms
    if molecule_config.get('smiles'):
        # Create from SMILES using RDKit
        from rdkit import Chem
        from rdkit.Chem import AllChem

        smiles = molecule_config['smiles']
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES string: {smiles}")

        # Add hydrogens and generate 3D coordinates
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)

        # Extract atoms and positions
        conf = mol.GetConformer()
        atoms = []
        for i in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(i)
            pos = conf.GetAtomPosition(i)
            atoms.append(
                Atom(
                    atom.GetSymbol(),
                    position=np.array([pos.x, pos.y, pos.z])
                )
            )
    elif molecule_config.get('atoms'):
        # Create from atoms
        atoms = []
        for atom_data in molecule_config['atoms']:
            atom = Atom(
                symbol=atom_data['symbol'],
                position=np.array([atom_data['x'], atom_data['y'], atom_data['z']])
            )
            atoms.append(atom)
    else:
        raise ValueError("Must provide either 'smiles' or 'atoms' in molecule_config")

    # Convert multiplicity to spin (spin = 2S where multiplicity = 2S+1)
    multiplicity = molecule_config.get('multiplicity', 1)
    spin = multiplicity - 1
    charge = molecule_config.get('charge', 0)
    basis = molecule_config.get('basis', 'sto-3g')

    # Try to create molecule, auto-correct multiplicity for odd electrons if needed
    try:
        molecule = Molecule(
            atoms=atoms,
            charge=charge,
            spin=spin,
            basis=basis
        )
    except RuntimeError as e:
        if "Electron number" in str(e) and "spin" in str(e) and "not consistent" in str(e):
            # Auto-correct: odd electrons need multiplicity=2 (doublet)
            # Calculate total electrons
            total_electrons = sum(Atom.ATOMIC_NUMBERS.get(atom.symbol, 0) for atom in atoms) - charge
            if total_electrons % 2 == 1:  # Odd number of electrons
                # Use doublet state (multiplicity=2, spin=1)
                spin = 1
                molecule = Molecule(
                    atoms=atoms,
                    charge=charge,
                    spin=spin,
                    basis=basis
                )
            else:
                raise  # Re-raise if it's a different issue
        else:
            raise  # Re-raise if it's a different error

    # Calculate qubits needed
    n_qubits = molecule.n_orbitals * 2  # Spin orbitals

    # Get method and ansatz
    method = backend_config.get('method', 'VQE')

    if method != 'VQE':
        # Non-VQE methods (SQD, SKQD): Generate real Qiskit circuit for HF state preparation
        # These methods use quantum circuits for subspace generation
        try:
            from qiskit import QuantumCircuit

            # Create a simple Hartree-Fock state preparation circuit
            qiskit_circuit = QuantumCircuit(n_qubits)

            # Add X gates to prepare HF state (occupy the first n_electrons qubits)
            for i in range(molecule.n_electrons):
                qiskit_circuit.x(i)

            # Add measurement for visualization
            qiskit_circuit.barrier()

            # Generate real Qiskit circuit image
            circuit_image = generate_matplotlib_circuit(qiskit_circuit)

            # Also generate ASCII diagram for backwards compatibility
            visualizer = CircuitVisualizer(n_qubits)
            for i in range(molecule.n_electrons):
                visualizer.add_gate('x', [i])
            diagram = visualizer.to_string()
            stats = visualizer.get_stats()

            return {
                'has_circuit': True,
                'diagram': diagram,
                'stats': stats,
                'circuit_image': circuit_image,  # Real Qiskit circuit image
                'n_qubits': n_qubits,
                'n_electrons': molecule.n_electrons,
                'n_parameters': 0,
                'method': method,
                'note': f'{method} uses quantum subspace expansion with HF state preparation'
            }

        except Exception as e:
            print(f"❌ Error generating {method} circuit: {e}")
            import traceback
            traceback.print_exc()

            # Fallback to ASCII only
            visualizer = CircuitVisualizer(n_qubits)
            for i in range(molecule.n_electrons):
                visualizer.add_gate('x', [i])
            diagram = visualizer.to_string()
            stats = visualizer.get_stats()

            return {
                'has_circuit': True,
                'diagram': diagram,
                'stats': stats,
                'n_qubits': n_qubits,
                'n_electrons': molecule.n_electrons,
                'n_parameters': 0,
                'method': method,
                'note': f'{method} uses classical computation with quantum state preparation'
            }

    ansatz_type = backend_config.get('ansatz', 'hardware_efficient')

    # Create ansatz circuit for VQE - support all ansatz types
    if ansatz_type == 'ucc':
        ansatz = UCCAnsatz(n_qubits=n_qubits, n_electrons=molecule.n_electrons)
    elif ansatz_type == 'real_amplitudes':
        ansatz = RealAmplitudesAnsatz(
            n_qubits=n_qubits,
            n_electrons=molecule.n_electrons,
            n_layers=2
        )
    elif ansatz_type == 'efficient_su2':
        ansatz = EfficientSU2Ansatz(
            n_qubits=n_qubits,
            n_electrons=molecule.n_electrons,
            n_layers=2
        )
    elif ansatz_type == 'two_local':
        ansatz = TwoLocalAnsatz(
            n_qubits=n_qubits,
            n_electrons=molecule.n_electrons,
            n_layers=2,
            rotation_gates='ry',
            entanglement='linear'
        )
    elif ansatz_type == 'ionic_governance':
        ansatz = IonicGovernanceAnsatz(
            n_qubits=n_qubits,
            n_electrons=molecule.n_electrons,
            n_layers=2
        )
    elif ansatz_type == 'covalent_governance':
        ansatz = CovalentGovernanceAnsatz(
            n_qubits=n_qubits,
            n_electrons=molecule.n_electrons,
            n_layers=2
        )
    elif ansatz_type == 'adaptive_governance':
        ansatz = AdaptiveGovernanceAnsatz(
            n_qubits=n_qubits,
            n_electrons=molecule.n_electrons,
            n_layers=2
        )
    else:
        # Default: hardware_efficient
        ansatz = HardwareEfficientAnsatz(
            n_qubits=n_qubits,
            n_electrons=molecule.n_electrons,
            n_layers=2  # Preview with 2 layers
        )

    # Build circuit
    circuit = ansatz.build_circuit()

    # Visualize (ASCII and structured gates for backwards compatibility)
    diagram, stats, structured_gates = visualize_from_circuit(circuit)

    # Generate REAL Qiskit circuit image using Qiskit's native drawer
    # Convert internal circuit to Qiskit QuantumCircuit
    qiskit_circuit = None
    circuit_image = None

    try:
        # Get the real Qiskit circuit from ansatz
        if hasattr(circuit, 'to_qiskit'):
            qiskit_circuit = circuit.to_qiskit()
        elif hasattr(ansatz, 'to_qiskit'):
            qiskit_circuit = ansatz.to_qiskit()
        else:
            print("⚠️  Ansatz does not have to_qiskit() method")

        # Generate Matplotlib circuit image using Qiskit's native drawer
        if qiskit_circuit is not None:
            print(f"🎨 Generating Qiskit circuit diagram for {qiskit_circuit.num_qubits} qubits, {qiskit_circuit.depth()} depth...")
            circuit_image = generate_matplotlib_circuit(qiskit_circuit)

            if circuit_image:
                print(f"✅ Real Qiskit circuit image generated: {len(circuit_image)} bytes (base64)")
            else:
                print("❌ Qiskit circuit image generation failed!")
        else:
            print("❌ Could not obtain Qiskit circuit from ansatz")

    except Exception as e:
        print(f"❌ Error generating Qiskit circuit: {e}")
        import traceback
        traceback.print_exc()

    return {
        'has_circuit': True,
        'diagram': diagram,
        'stats': stats,
        'gates': structured_gates,  # Add structured gate data
        'circuit_image': circuit_image,  # Base64 PNG image from REAL Qiskit circuit
        'n_qubits': n_qubits,
        'n_electrons': molecule.n_electrons,
        'n_parameters': len(circuit.parameters),
        'ansatz_type': ansatz_type,
        'mapper_type': backend_config.get('mapper', 'jordan_wigner')
    }
