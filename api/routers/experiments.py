"""
Experiments API router.

Handles experiment creation, listing, retrieval, and deletion.
"""

import logging
from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy.orm import Session
from typing import List, Optional
from datetime import datetime

from api.database import get_db
from api.models.experiment import Experiment
from api.models.queue import QueueItem
from api.utils.validators import ExperimentCreate
from api.utils.exceptions import ExperimentNotFoundError
from api.services.job_queue import job_queue
from api.services.experiment_service import experiment_service

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/experiments", tags=["experiments"])


@router.post("/", status_code=201)
def create_experiment(
    experiment_data: ExperimentCreate,
    db: Session = Depends(get_db)
):
    """
    Create a new experiment.

    If execute_immediately=True, adds to job queue for immediate execution.
    Otherwise, creates experiment record for later queuing.
    """
    # Create experiment record
    experiment = Experiment(
        name=experiment_data.name or f"Experiment {datetime.now().strftime('%Y%m%d_%H%M%S')}",
        status="queued" if experiment_data.execute_immediately else "pending",
        smiles=experiment_data.molecule.smiles,
        molecule_data=experiment_data.molecule.dict(),
        configuration=experiment_data.configuration.dict()
    )

    db.add(experiment)
    db.commit()
    db.refresh(experiment)

    # If execute immediately, add to job queue
    if experiment_data.execute_immediately:
        # Create queue item
        queue_item = QueueItem(
            experiment_id=experiment.id,
            status="queued",
            priority=0
        )
        db.add(queue_item)
        db.commit()

        # Add to background queue
        job_queue.add_job(experiment.id, priority=0)

    return {
        "id": experiment.id,
        "name": experiment.name,
        "status": experiment.status,
        "message": "Experiment created and queued for execution" if experiment_data.execute_immediately else "Experiment created"
    }


@router.get("/")
def list_experiments(
    status: Optional[str] = Query(None, description="Filter by status"),
    limit: int = Query(100, ge=1, le=1000, description="Max number of results"),
    offset: int = Query(0, ge=0, description="Number of results to skip"),
    db: Session = Depends(get_db)
):
    """
    List experiments with optional filtering.

    Supports filtering by status and pagination.
    """
    query = db.query(Experiment)

    # Filter by status
    if status:
        query = query.filter(Experiment.status == status)

    # Order by creation date (newest first)
    query = query.order_by(Experiment.created_at.desc())

    # Pagination
    total = query.count()
    experiments = query.offset(offset).limit(limit).all()

    return {
        "total": total,
        "offset": offset,
        "limit": limit,
        "experiments": [exp.to_dict() for exp in experiments]
    }


@router.get("/{experiment_id}")
def get_experiment(
    experiment_id: int,
    db: Session = Depends(get_db)
):
    """
    Get experiment details by ID.

    Returns complete experiment data including results and convergence history.
    """
    experiment = db.query(Experiment).filter(Experiment.id == experiment_id).first()

    if not experiment:
        raise ExperimentNotFoundError(experiment_id)

    return experiment.to_dict()


@router.delete("/{experiment_id}")
def delete_experiment(
    experiment_id: int,
    db: Session = Depends(get_db)
):
    """
    Delete an experiment.

    Also removes associated queue items.
    """
    experiment = db.query(Experiment).filter(Experiment.id == experiment_id).first()

    if not experiment:
        raise ExperimentNotFoundError(experiment_id)

    # Delete associated queue items
    db.query(QueueItem).filter(QueueItem.experiment_id == experiment_id).delete()

    # Delete experiment
    db.delete(experiment)
    db.commit()

    return {
        "message": f"Experiment {experiment_id} deleted successfully"
    }


@router.get("/{experiment_id}/status")
def get_experiment_status(
    experiment_id: int,
    db: Session = Depends(get_db)
):
    """
    Get experiment status and progress.

    Useful for polling during execution.
    """
    experiment = db.query(Experiment).filter(Experiment.id == experiment_id).first()

    if not experiment:
        raise ExperimentNotFoundError(experiment_id)

    # Calculate progress percentage
    progress = 0
    if experiment.status == "completed":
        progress = 100
    elif experiment.status == "running" and experiment.convergence_data:
        # Estimate progress based on convergence data
        # Assume max 1000 iterations (from config)
        max_iter = experiment.configuration.get('max_iterations', 1000)
        current_iter = len(experiment.convergence_data)
        progress = min(100, int(100 * current_iter / max_iter))

    return {
        "id": experiment.id,
        "status": experiment.status,
        "progress": progress,
        "current_iteration": len(experiment.convergence_data) if experiment.convergence_data else 0,
        "started_at": experiment.started_at.isoformat() if experiment.started_at else None,
        "error_message": experiment.error_message
    }


@router.get("/{experiment_id}/convergence")
def get_convergence_data(
    experiment_id: int,
    db: Session = Depends(get_db)
):
    """
    Get real-time convergence data.

    Returns energy convergence history for plotting.
    """
    experiment = db.query(Experiment).filter(Experiment.id == experiment_id).first()

    if not experiment:
        raise ExperimentNotFoundError(experiment_id)

    return {
        "id": experiment.id,
        "status": experiment.status,
        "convergence_data": experiment.convergence_data or []
    }


@router.get("/{experiment_id}/circuit")
def get_circuit_visualization(
    experiment_id: int,
    format: str = Query("json", description="Output format: json, qasm, or ascii"),
    db: Session = Depends(get_db)
):
    """
    Get quantum circuit visualization data.

    Returns the parametrized quantum circuit used in the experiment.
    Supports multiple formats: JSON structure, QASM string, or ASCII art.
    """
    experiment = db.query(Experiment).filter(Experiment.id == experiment_id).first()

    if not experiment:
        raise ExperimentNotFoundError(experiment_id)

    # Check if experiment has completed or has circuit data stored
    if not experiment.results:
        raise HTTPException(
            status_code=400,
            detail="Experiment has no results yet. Circuit data available after execution."
        )

    # Extract circuit data from experiment configuration
    config = experiment.configuration
    method = config.get('method', 'VQE').upper()

    # Only VQE and SQD have circuits
    if method not in ['VQE', 'SQD']:
        raise HTTPException(
            status_code=400,
            detail=f"Method {method} does not use quantum circuits"
        )

    # Reconstruct circuit from configuration
    try:
        from api.services.experiment_service import experiment_service

        # Recreate molecule
        molecule = experiment_service.create_molecule(experiment.molecule_data)

        # Get circuit parameters from config
        ansatz_type = config.get('ansatz', 'ucc')
        n_qubits = 2 * molecule.n_orbitals
        n_electrons = molecule.n_electrons

        # Build circuit structure based on ansatz type
        circuit_data = {
            "n_qubits": n_qubits,
            "n_electrons": n_electrons,
            "ansatz_type": ansatz_type,
            "method": method,
            "metadata": {
                "depth": None,
                "gate_count": None,
                "parameter_count": None
            }
        }

        # Generate circuit representation based on format
        if format == "json":
            # Return structured circuit data
            circuit_data["gates"] = _generate_circuit_gates(ansatz_type, n_qubits, n_electrons)
            circuit_data["metadata"]["parameter_count"] = len([g for g in circuit_data["gates"] if g.get("parameterized")])
            circuit_data["metadata"]["gate_count"] = len(circuit_data["gates"])

        elif format == "ascii":
            # Generate ASCII art representation
            circuit_data["ascii"] = _generate_circuit_ascii(ansatz_type, n_qubits, n_electrons)

        elif format == "qasm":
            # Generate OpenQASM representation
            circuit_data["qasm"] = _generate_circuit_qasm(ansatz_type, n_qubits, n_electrons)

        else:
            raise HTTPException(status_code=400, detail=f"Invalid format: {format}. Use 'json', 'ascii', or 'qasm'")

        return circuit_data

    except Exception as e:
        logger.error(f"Failed to generate circuit visualization: {e}")
        raise HTTPException(status_code=500, detail=f"Failed to generate circuit: {str(e)}")


@router.get("/{experiment_id}/report")
def generate_experiment_report(
    experiment_id: int,
    format: str = Query("json", description="Report format: json or markdown"),
    db: Session = Depends(get_db)
):
    """
    Generate comprehensive experiment report.

    Returns detailed report including:
    - Molecular structure and properties
    - Calculation method and parameters
    - Energy results and convergence
    - Analysis (bond orders, dipole moments, etc.)
    - Visualizations (if available)
    """
    experiment = db.query(Experiment).filter(Experiment.id == experiment_id).first()

    if not experiment:
        raise ExperimentNotFoundError(experiment_id)

    if experiment.status != "completed":
        raise HTTPException(
            status_code=400,
            detail="Report can only be generated for completed experiments"
        )

    # Build comprehensive report
    report = {
        "experiment_id": experiment.id,
        "name": experiment.name,
        "created_at": experiment.created_at.isoformat(),
        "completed_at": experiment.completed_at.isoformat() if experiment.completed_at else None,
        "execution_time": (experiment.completed_at - experiment.started_at).total_seconds() if experiment.completed_at and experiment.started_at else None,

        "molecule": {
            "smiles": experiment.smiles,
            "formula": experiment.results.get('molecule_formula'),
            "n_electrons": experiment.results.get('n_electrons'),
            "n_orbitals": experiment.results.get('n_orbitals'),
            "basis": experiment.molecule_data.get('basis', 'sto-3g'),
            "charge": experiment.molecule_data.get('charge', 0),
            "multiplicity": experiment.molecule_data.get('multiplicity', 1)
        },

        "method": {
            "type": experiment.configuration.get('method', 'VQE'),
            "ansatz": experiment.configuration.get('ansatz'),
            "mapper": experiment.configuration.get('mapper'),
            "optimizer": experiment.configuration.get('optimizer'),
            "backend": experiment.configuration.get('backend'),
            "max_iterations": experiment.configuration.get('max_iterations'),
        },

        "results": {
            "energy": experiment.results.get('energy'),
            "hf_energy": experiment.results.get('hf_energy'),
            "correlation_energy": experiment.results.get('correlation_energy'),
            "converged": experiment.results.get('converged'),
            "iterations": experiment.results.get('iterations'),
        },

        "convergence": {
            "data": experiment.convergence_data or [],
            "final_iteration": len(experiment.convergence_data) if experiment.convergence_data else 0,
        },

        "analysis": experiment.results.get('analysis', {}),
        "properties": experiment.results.get('properties', {}),
        "validation": experiment.results.get('validation', {}),
    }

    if format == "json":
        return report
    elif format == "markdown":
        # Generate markdown report
        markdown = _generate_markdown_report(report)
        return {"format": "markdown", "content": markdown}
    else:
        raise HTTPException(status_code=400, detail=f"Invalid format: {format}. Use 'json' or 'markdown'")


# Helper functions for circuit generation
def _generate_circuit_gates(ansatz_type: str, n_qubits: int, n_electrons: int) -> list:
    """Generate structured gate list for circuit visualization."""
    gates = []

    # Initial state preparation
    n_alpha = n_electrons // 2
    for i in range(n_alpha):
        gates.append({"type": "x", "qubits": [i], "parameterized": False})
    for i in range(n_electrons - n_alpha):
        gates.append({"type": "x", "qubits": [n_qubits // 2 + i], "parameterized": False})

    if ansatz_type.lower() == 'ucc':
        # UCC ansatz gates (simplified representation)
        param_idx = 0
        # Singles
        for i in range(n_electrons // 2):
            for a in range(n_electrons // 2, n_qubits // 2):
                gates.append({"type": "ry", "qubits": [i], "parameter": f"theta_{param_idx}", "parameterized": True})
                gates.append({"type": "cx", "qubits": [i, a], "parameterized": False})
                param_idx += 1

        # Doubles (simplified)
        for i in range(n_electrons // 2):
            for j in range(i+1, n_electrons // 2):
                for a in range(n_electrons // 2, n_qubits // 2):
                    gates.append({"type": "ry", "qubits": [i], "parameter": f"theta_{param_idx}", "parameterized": True})
                    gates.append({"type": "cx", "qubits": [i, j], "parameterized": False})
                    gates.append({"type": "cx", "qubits": [j, a], "parameterized": False})
                    param_idx += 1

    elif ansatz_type.lower() == 'hardware_efficient':
        # Hardware-efficient ansatz (layered structure)
        n_layers = 3
        param_idx = 0
        for layer in range(n_layers):
            # Rotation layer
            for q in range(n_qubits):
                gates.append({"type": "ry", "qubits": [q], "parameter": f"theta_{param_idx}", "parameterized": True})
                param_idx += 1
                gates.append({"type": "rz", "qubits": [q], "parameter": f"theta_{param_idx}", "parameterized": True})
                param_idx += 1

            # Entanglement layer (linear)
            for q in range(n_qubits - 1):
                gates.append({"type": "cx", "qubits": [q, q+1], "parameterized": False})

    elif ansatz_type.lower() == 'governance':
        # Governance-aware ansatz (simplified)
        param_idx = 0
        for q in range(n_qubits):
            gates.append({"type": "ry", "qubits": [q], "parameter": f"theta_{param_idx}", "parameterized": True})
            param_idx += 1
        for q in range(0, n_qubits-1, 2):
            gates.append({"type": "cx", "qubits": [q, q+1], "parameterized": False})

    return gates


def _generate_circuit_ascii(ansatz_type: str, n_qubits: int, n_electrons: int) -> str:
    """Generate ASCII art representation of circuit."""
    ascii_lines = []
    ascii_lines.append(f"# {ansatz_type.upper()} Ansatz Circuit ({n_qubits} qubits)")
    ascii_lines.append("")

    for q in range(n_qubits):
        line = f"q{q}: "
        if q < n_electrons // 2 or (n_qubits // 2 <= q < n_qubits // 2 + (n_electrons - n_electrons // 2)):
            line += "|0>--[X]--"
        else:
            line += "|0>-------"

        if ansatz_type.lower() == 'ucc':
            line += "--[RY(θ)]--[CX]--"
        elif ansatz_type.lower() == 'hardware_efficient':
            line += "--[RY(θ)]--[RZ(θ)]--[CX]--"
        else:
            line += "--[Gates]--"

        line += "--|M|"
        ascii_lines.append(line)

    ascii_lines.append("")
    ascii_lines.append("Legend: X=Pauli-X, RY=Y-rotation, RZ=Z-rotation, CX=CNOT, M=Measurement")

    return "\n".join(ascii_lines)


def _generate_circuit_qasm(ansatz_type: str, n_qubits: int, n_electrons: int) -> str:
    """Generate OpenQASM representation of circuit."""
    qasm_lines = []
    qasm_lines.append("OPENQASM 2.0;")
    qasm_lines.append('include "qelib1.inc";')
    qasm_lines.append(f"qreg q[{n_qubits}];")
    qasm_lines.append(f"creg c[{n_qubits}];")
    qasm_lines.append("")
    qasm_lines.append(f"// {ansatz_type.upper()} Ansatz")
    qasm_lines.append("// State preparation")

    # State preparation
    n_alpha = n_electrons // 2
    for i in range(n_alpha):
        qasm_lines.append(f"x q[{i}];")
    for i in range(n_electrons - n_alpha):
        qasm_lines.append(f"x q[{n_qubits // 2 + i}];")

    qasm_lines.append("")
    qasm_lines.append("// Parametrized gates")

    if ansatz_type.lower() == 'ucc':
        qasm_lines.append("// UCC excitations (singles and doubles)")
        qasm_lines.append("ry(theta_0) q[0];")
        qasm_lines.append("cx q[0], q[2];")
    elif ansatz_type.lower() == 'hardware_efficient':
        qasm_lines.append("// Hardware-efficient layers")
        for q in range(n_qubits):
            qasm_lines.append(f"ry(theta_{q}) q[{q}];")
            qasm_lines.append(f"rz(theta_{q + n_qubits}) q[{q}];")
        for q in range(n_qubits - 1):
            qasm_lines.append(f"cx q[{q}], q[{q+1}];")

    qasm_lines.append("")
    qasm_lines.append("// Measurement")
    for q in range(n_qubits):
        qasm_lines.append(f"measure q[{q}] -> c[{q}];")

    return "\n".join(qasm_lines)


def _generate_markdown_report(report: dict) -> str:
    """Generate markdown formatted report."""
    md = []
    md.append(f"# Experiment Report: {report['name']}")
    md.append("")
    md.append(f"**Experiment ID:** {report['experiment_id']}")
    md.append(f"**Created:** {report['created_at']}")
    md.append(f"**Completed:** {report['completed_at']}")
    md.append(f"**Execution Time:** {report['execution_time']:.2f} seconds" if report['execution_time'] else "")
    md.append("")

    md.append("## Molecule")
    md.append("")
    mol = report['molecule']
    md.append(f"- **SMILES:** `{mol['smiles']}`")
    md.append(f"- **Formula:** {mol['formula']}")
    md.append(f"- **Electrons:** {mol['n_electrons']}")
    md.append(f"- **Orbitals:** {mol['n_orbitals']}")
    md.append(f"- **Basis Set:** {mol['basis']}")
    md.append(f"- **Charge:** {mol['charge']}")
    md.append(f"- **Multiplicity:** {mol['multiplicity']}")
    md.append("")

    md.append("## Calculation Method")
    md.append("")
    method = report['method']
    md.append(f"- **Method:** {method['type']}")
    if method.get('ansatz'):
        md.append(f"- **Ansatz:** {method['ansatz']}")
    if method.get('mapper'):
        md.append(f"- **Qubit Mapper:** {method['mapper']}")
    if method.get('optimizer'):
        md.append(f"- **Optimizer:** {method['optimizer']}")
    if method.get('backend'):
        md.append(f"- **Backend:** {method['backend']}")
    md.append("")

    md.append("## Results")
    md.append("")
    results = report['results']
    md.append(f"- **Final Energy:** {results['energy']:.8f} Hartree")
    if results.get('hf_energy'):
        md.append(f"- **HF Reference:** {results['hf_energy']:.8f} Hartree")
    if results.get('correlation_energy'):
        corr_kcal = results['correlation_energy'] * 627.509
        md.append(f"- **Correlation Energy:** {results['correlation_energy']:.8f} Hartree ({corr_kcal:.2f} kcal/mol)")
    md.append(f"- **Converged:** {'Yes' if results['converged'] else 'No'}")
    md.append(f"- **Iterations:** {results['iterations']}")
    md.append("")

    if report['convergence']['data']:
        md.append("## Convergence")
        md.append("")
        md.append(f"Final iteration: {report['convergence']['final_iteration']}")
        md.append("")

    if report.get('properties'):
        md.append("## Molecular Properties")
        md.append("")
        for key, value in report['properties'].items():
            md.append(f"- **{key}:** {value}")
        md.append("")

    return "\n".join(md)
