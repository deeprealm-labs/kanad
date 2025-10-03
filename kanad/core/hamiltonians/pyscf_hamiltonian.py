"""
PySCF-based Hamiltonian builder for Kanad framework.

Uses PySCF to handle molecular integrals and OpenFermion for fermionic-to-qubit
transformations. This ensures correctness by leveraging battle-tested libraries.

Requires:
    - pyscf
    - openfermion
"""

import numpy as np
from typing import TYPE_CHECKING, List

if TYPE_CHECKING:
    from kanad.core.atom import Atom

try:
    from pyscf import gto, scf
    from openfermion import FermionOperator, jordan_wigner, bravyi_kitaev
    from qiskit.quantum_info import SparsePauliOp
    PYSCF_AVAILABLE = True
except ImportError:
    PYSCF_AVAILABLE = False


class PySCFHamiltonian:
    """
    Build quantum Hamiltonians using PySCF + OpenFermion.

    This class provides a validated path from molecular geometry to Pauli operators
    by using PySCF's molecular integrals and OpenFermion's fermionic transformations.

    Usage:
        >>> from kanad.core.atom import Atom
        >>> import numpy as np
        >>>
        >>> h1 = Atom('H', np.array([0.0, 0.0, 0.0]))
        >>> h2 = Atom('H', np.array([0.0, 0.0, 0.74]))
        >>>
        >>> pauli_op = PySCFHamiltonian.from_atoms(
        >>>     [h1, h2],
        >>>     basis='sto-3g',
        >>>     mapper='jordan_wigner'
        >>> )
    """

    @staticmethod
    def from_atoms(
        atoms: List['Atom'],
        basis: str = 'sto-3g',
        mapper: str = 'jordan_wigner',
        charge: int = 0,
        spin: int = 0
    ) -> SparsePauliOp:
        """
        Build Pauli operator Hamiltonian from list of atoms.

        Args:
            atoms: List of Atom objects with symbols and positions
            basis: Basis set name (e.g., 'sto-3g', '6-31g')
            mapper: Transformation type ('jordan_wigner' or 'bravyi_kitaev')
            charge: Total molecular charge
            spin: Spin multiplicity - 1 (0 for singlet, 1 for doublet, etc.)

        Returns:
            SparsePauliOp representing the molecular Hamiltonian

        Raises:
            ImportError: If PySCF or OpenFermion not installed
        """
        if not PYSCF_AVAILABLE:
            raise ImportError(
                "PySCF and/or OpenFermion not installed.\n"
                "Install with: pip install pyscf openfermion"
            )

        # Build PySCF molecule
        mol = PySCFHamiltonian._build_pyscf_molecule(atoms, basis, charge, spin)

        # Run SCF to get molecular orbitals
        mf = scf.RHF(mol)
        mf.verbose = 0
        mf.kernel()

        # Build fermionic Hamiltonian in MO basis
        fermion_hamiltonian = PySCFHamiltonian._build_fermion_hamiltonian(mol, mf)

        # Apply fermionic-to-qubit transformation
        if mapper == 'jordan_wigner':
            qubit_hamiltonian = jordan_wigner(fermion_hamiltonian)
        elif mapper == 'bravyi_kitaev':
            qubit_hamiltonian = bravyi_kitaev(fermion_hamiltonian)
        else:
            raise ValueError(f"Unknown mapper: {mapper}. Use 'jordan_wigner' or 'bravyi_kitaev'")

        # Convert to Qiskit
        n_orbitals = mol.nao
        n_qubits = n_orbitals * 2  # Spin orbitals
        pauli_op = PySCFHamiltonian._openfermion_to_qiskit(qubit_hamiltonian, n_qubits)

        return pauli_op

    @staticmethod
    def _build_pyscf_molecule(
        atoms: List['Atom'],
        basis: str,
        charge: int,
        spin: int
    ):
        """
        Build PySCF Mole object from Kanad atoms.

        Args:
            atoms: List of Atom objects
            basis: Basis set name
            charge: Molecular charge
            spin: Spin (2S)

        Returns:
            PySCF Mole object
        """
        # Convert atoms to PySCF format: [(symbol, (x, y, z)), ...]
        geometry = []
        for atom in atoms:
            symbol = atom.symbol
            position = atom.position  # In Angstroms
            geometry.append((symbol, tuple(position)))

        mol = gto.M(
            atom=geometry,
            basis=basis,
            charge=charge,
            spin=spin,
            unit='angstrom'
        )

        return mol

    @staticmethod
    def _build_fermion_hamiltonian(mol, mf) -> FermionOperator:
        """
        Build OpenFermion FermionOperator from PySCF calculation.

        This uses molecular orbital (MO) basis integrals which are properly
        orthogonalized and easier to work with than atomic orbital (AO) integrals.

        Args:
            mol: PySCF Mole object
            mf: PySCF mean-field object (SCF result)

        Returns:
            FermionOperator representing the molecular Hamiltonian
        """
        from pyscf import ao2mo

        # Get number of spatial orbitals
        n_orbitals = mol.nao

        # Nuclear repulsion energy
        nuclear_repulsion = mol.energy_nuc()

        # One-electron integrals in MO basis
        h_core_ao = mf.get_hcore()
        mo_coeff = mf.mo_coeff
        h_mo = mo_coeff.T @ h_core_ao @ mo_coeff

        # Two-electron integrals in MO basis
        # ao2mo.kernel returns ERI in physicist notation (pq|rs)
        # in a compressed format for efficiency
        eri_mo = ao2mo.kernel(mol, mo_coeff)

        # Restore to 4D tensor: eri_mo is in chemist notation [pq,rs]
        # ao2mo.restore creates full 4D tensor [p,q,r,s] = (pr|qs) physicist notation
        eri_mo_4d = ao2mo.restore(1, eri_mo, n_orbitals)  # Chemist notation [i,j,k,l] = (ij|kl)

        # Build fermionic Hamiltonian
        ferm_ham = FermionOperator()

        # 1. Nuclear repulsion (constant)
        ferm_ham += FermionOperator('', nuclear_repulsion)

        # 2. One-electron terms: h_pq a†_p a_q for both spins
        # Spin orbital convention: 0↑, 0↓, 1↑, 1↓, ... (interleaved)
        for p in range(n_orbitals):
            for q in range(n_orbitals):
                if abs(h_mo[p, q]) > 1e-12:
                    # Spin-up: 2*p
                    ferm_ham += FermionOperator(f'{2*p}^ {2*q}', h_mo[p, q])
                    # Spin-down: 2*p+1
                    ferm_ham += FermionOperator(f'{2*p+1}^ {2*q+1}', h_mo[p, q])

        # 3. Two-electron terms: (1/2) Σ ⟨pq|rs⟩ a†_p a†_q a_s a_r
        # eri_mo_4d[p,q,r,s] = (pq|rs) in chemist notation
        for p in range(n_orbitals):
            for q in range(n_orbitals):
                for r in range(n_orbitals):
                    for s in range(n_orbitals):
                        v_pqrs = eri_mo_4d[p, q, r, s]

                        if abs(v_pqrs) > 1e-12:
                            # Hamiltonian: (1/2) Σ ⟨pq|rs⟩ a†_p a†_q a_s a_r
                            # Add for all spin combinations
                            # Both spin-up
                            ferm_ham += FermionOperator(
                                f'{2*p}^ {2*q}^ {2*s} {2*r}',
                                0.5 * v_pqrs
                            )
                            # p,r spin-up; q,s spin-down
                            ferm_ham += FermionOperator(
                                f'{2*p}^ {2*q+1}^ {2*s+1} {2*r}',
                                0.5 * v_pqrs
                            )
                            # p,r spin-down; q,s spin-up
                            ferm_ham += FermionOperator(
                                f'{2*p+1}^ {2*q}^ {2*s} {2*r+1}',
                                0.5 * v_pqrs
                            )
                            # Both spin-down
                            ferm_ham += FermionOperator(
                                f'{2*p+1}^ {2*q+1}^ {2*s+1} {2*r+1}',
                                0.5 * v_pqrs
                            )

        return ferm_ham

    @staticmethod
    def _openfermion_to_qiskit(qubit_op, n_qubits: int) -> SparsePauliOp:
        """
        Convert OpenFermion QubitOperator to Qiskit SparsePauliOp.

        Args:
            qubit_op: OpenFermion QubitOperator
            n_qubits: Number of qubits

        Returns:
            Qiskit SparsePauliOp
        """
        pauli_dict = {}

        for term, coeff in qubit_op.terms.items():
            if len(term) == 0:
                # Identity term
                pauli_str = 'I' * n_qubits
            else:
                # Build Pauli string
                pauli_list = ['I'] * n_qubits
                for qubit_idx, pauli_op in term:
                    pauli_list[qubit_idx] = pauli_op
                pauli_str = ''.join(pauli_list)

            # Accumulate coefficients
            if pauli_str in pauli_dict:
                pauli_dict[pauli_str] += coeff
            else:
                pauli_dict[pauli_str] = coeff

        # Filter negligible terms
        pauli_dict = {k: v for k, v in pauli_dict.items() if abs(v) > 1e-12}

        # Convert to Qiskit
        if len(pauli_dict) == 0:
            return SparsePauliOp(['I' * n_qubits], [0.0])

        pauli_strings = list(pauli_dict.keys())
        coefficients = list(pauli_dict.values())

        return SparsePauliOp(pauli_strings, coefficients)
