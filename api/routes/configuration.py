"""
Configuration endpoints - Returns ALL available options from Kanad framework
Based on actual framework inspection
"""

from fastapi import APIRouter

router = APIRouter()


@router.get("/options")
async def get_configuration_options():
    """
    Get all available configuration options from the Kanad framework.
    """
    return {
        # Methods from kanad/solvers/ (Excited states moved to Advanced Analysis)
        "methods": [
            {"value": "HF", "label": "Hartree-Fock", "description": "Classical mean-field (ground state)", "requires_ansatz": False, "requires_mapper": False, "implemented": True},
            {"value": "VQE", "label": "VQE", "description": "Variational Quantum Eigensolver (ground state)", "requires_ansatz": True, "requires_mapper": True, "implemented": True},
            {"value": "SQD", "label": "SQD", "description": "Subspace Quantum Diagonalization (ground state)", "requires_ansatz": False, "requires_mapper": False, "implemented": True},
        ],

        # Ansatze from kanad/ansatze/
        "ansatze": [
            {"value": "ucc", "label": "UCC", "description": "Unitary Coupled Cluster"},
            {"value": "uccsd", "label": "UCCSD", "description": "UCC Singles + Doubles"},
            {"value": "ucc_s", "label": "UCC Singles", "description": "Singles only"},
            {"value": "ucc_d", "label": "UCC Doubles", "description": "Doubles only"},
            {"value": "hardware_efficient", "label": "Hardware Efficient", "description": "Low-depth NISQ ansatz"},
            {"value": "real_amplitudes", "label": "Real Amplitudes", "description": "Real parameters only"},
            {"value": "efficient_su2", "label": "EfficientSU2", "description": "SU2 rotations"},
            {"value": "two_local", "label": "Two Local", "description": "General two-local"},
            {"value": "ionic_governance", "label": "Ionic Governance", "description": "For ionic bonds"},
            {"value": "covalent_governance", "label": "Covalent Governance", "description": "For covalent bonds"},
            {"value": "adaptive_governance", "label": "Adaptive Governance", "description": "Adaptive ansatz"},
        ],

        # Mappers from kanad/core/mappers/
        "mappers": [
            {"value": "jordan_wigner", "label": "Jordan-Wigner", "description": "Standard mapping (best for ionic)"},
            {"value": "bravyi_kitaev", "label": "Bravyi-Kitaev", "description": "Reduced Pauli weight"},
            {"value": "hybrid_orbital", "label": "Hybrid Orbital", "description": "Paired orbitals (best for covalent)"},
        ],

        # Optimizers - scipy.optimize.minimize methods
        "optimizers": [
            {"value": "SLSQP", "label": "SLSQP", "description": "Sequential Least Squares (default)"},
            {"value": "COBYLA", "label": "COBYLA", "description": "Constrained optimization"},
            {"value": "L-BFGS-B", "label": "L-BFGS-B", "description": "Limited-memory BFGS"},
            {"value": "Powell", "label": "Powell", "description": "Powell's method"},
            {"value": "Nelder-Mead", "label": "Nelder-Mead", "description": "Simplex method"},
            {"value": "CG", "label": "CG", "description": "Conjugate Gradient"},
            {"value": "BFGS", "label": "BFGS", "description": "BFGS"},
            {"value": "TNC", "label": "TNC", "description": "Truncated Newton"},
        ],

        # Basis sets from kanad/core/integrals/basis_registry.py
        "basis_sets": [
            {"value": "sto-3g", "label": "STO-3G", "category": "minimal"},
            {"value": "sto-6g", "label": "STO-6G", "category": "minimal"},
            {"value": "minao", "label": "MinAO", "category": "minimal"},
            {"value": "3-21g", "label": "3-21G", "category": "split_valence"},
            {"value": "6-31g", "label": "6-31G", "category": "split_valence"},
            {"value": "6-311g", "label": "6-311G", "category": "split_valence"},
            {"value": "6-31g*", "label": "6-31G*", "category": "polarization"},
            {"value": "6-31g**", "label": "6-31G**", "category": "polarization"},
            {"value": "6-311g**", "label": "6-311G**", "category": "polarization"},
            {"value": "6-31+g", "label": "6-31+G", "category": "diffuse"},
            {"value": "6-31++g", "label": "6-31++G", "category": "diffuse"},
            {"value": "6-31++g**", "label": "6-31++G**", "category": "diffuse"},
            {"value": "cc-pvdz", "label": "cc-pVDZ", "category": "correlation_consistent"},
            {"value": "cc-pvtz", "label": "cc-pVTZ", "category": "correlation_consistent"},
            {"value": "cc-pvqz", "label": "cc-pVQZ", "category": "correlation_consistent"},
            {"value": "aug-cc-pvdz", "label": "aug-cc-pVDZ", "category": "correlation_consistent"},
            {"value": "aug-cc-pvtz", "label": "aug-cc-pVTZ", "category": "correlation_consistent"},
            {"value": "def2-svp", "label": "def2-SVP", "category": "def2"},
            {"value": "def2-tzvp", "label": "def2-TZVP", "category": "def2"},
            {"value": "def2-tzvpp", "label": "def2-TZVPP", "category": "def2"},
            {"value": "def2-qzvp", "label": "def2-QZVP", "category": "def2"},
            {"value": "gth-dzvp", "label": "GTH-DZVP", "category": "periodic"},
            {"value": "gth-tzvp", "label": "GTH-TZVP", "category": "periodic"},
        ],

        # Backends from kanad/backends/
        "backends": [
            {"value": "classical", "label": "Classical (Statevector)", "requires_credentials": False, "type": "simulator"},
            {"value": "qasm", "label": "QASM Simulator", "requires_credentials": False, "type": "simulator"},
            {"value": "ibm_quantum", "label": "IBM Quantum", "requires_credentials": True, "type": "hardware"},
            {"value": "bluequbit", "label": "BlueQubit", "requires_credentials": True, "type": "gpu_simulator"},
        ],

        "ibm_backends": [
            {"value": "ibm_torino", "label": "IBM Torino", "qubits": 133},
            {"value": "ibm_brisbane", "label": "IBM Brisbane", "qubits": 127},
            {"value": "ibm_kyoto", "label": "IBM Kyoto", "qubits": 127},
            {"value": "ibm_osaka", "label": "IBM Osaka", "qubits": 127},
        ],

        "bond_types": [
            {"value": "auto", "label": "Auto-detect"},
            {"value": "ionic", "label": "Ionic", "best_mapper": "jordan_wigner", "best_ansatz": "ionic_governance"},
            {"value": "covalent", "label": "Covalent", "best_mapper": "hybrid_orbital", "best_ansatz": "covalent_governance"},
            {"value": "metallic", "label": "Metallic", "best_mapper": "jordan_wigner", "best_ansatz": "hardware_efficient"},
        ],

        "hamiltonians": [
            {"value": "molecular", "label": "Molecular Hamiltonian"},
            {"value": "ionic", "label": "Ionic Hamiltonian"},
            {"value": "covalent", "label": "Covalent Hamiltonian"},
            {"value": "metallic", "label": "Metallic Hamiltonian"},
            {"value": "periodic", "label": "Periodic Hamiltonian"},
        ],

        "analysis_options": [
            {"value": "energy_decomposition", "label": "Energy Decomposition"},
            {"value": "bond_analysis", "label": "Bond Analysis"},
            {"value": "dipole_moment", "label": "Dipole Moment"},
            {"value": "polarizability", "label": "Polarizability"},
            {"value": "thermochemistry", "label": "Thermochemistry"},
            {"value": "spectroscopy", "label": "Spectroscopy"},
            {"value": "vibrational", "label": "Vibrational Analysis"},
            {"value": "uncertainty", "label": "Uncertainty Analysis"},
            {"value": "bond_scan", "label": "Bond Scan"},
            {"value": "dos", "label": "Density of States"},
        ],

        "optimization_options": [
            {"value": "geometry", "label": "Geometry Optimization"},
            {"value": "orbitals", "label": "Orbital Optimization"},
            {"value": "circuit", "label": "Circuit Optimization"},
            {"value": "adaptive", "label": "Adaptive VQE"},
        ],
    }
