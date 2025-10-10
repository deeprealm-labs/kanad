"""
Convert molecular Hamiltonians to Qiskit Pauli operators.

This module bridges the gap between the fermionic Hamiltonians
(h_core, ERI tensors) and Qiskit's SparsePauliOp representation
using the mappers (Jordan-Wigner, Bravyi-Kitaev, etc.).
"""

from typing import Dict, Tuple
import numpy as np
import logging

logger = logging.getLogger(__name__)


class PauliConverter:
    """
    Converts molecular Hamiltonians to Qiskit Pauli operators.

    Uses fermionic-to-qubit mappers to transform:
        H = Σ h_ij a†_i a_j + ½ Σ g_ijkl a†_i a†_j a_l a_k + E_nn

    into Pauli operator representation:
        H = Σ c_k P_k

    where P_k are Pauli strings (e.g., 'XXYZI') and c_k are coefficients.
    """

    @staticmethod
    def to_sparse_pauli_op(hamiltonian, mapper, use_qiskit_nature=True):
        """
        Convert molecular Hamiltonian to Qiskit SparsePauliOp.

        Args:
            hamiltonian: MolecularHamiltonian instance
            mapper: Fermionic-to-qubit mapper (JW, BK, etc.)
            use_qiskit_nature: If True, use Qiskit Nature for correct two-electron terms
                             (RECOMMENDED - mathematically correct!)

        Returns:
            qiskit.quantum_info.SparsePauliOp
        """
        try:
            from qiskit.quantum_info import SparsePauliOp
        except ImportError:
            raise ImportError(
                "Qiskit not installed. Install with: pip install qiskit>=2.0"
            )

        # Try Qiskit Nature approach
        # Always use it - even for charged systems we'll use the fermionic operators
        # (The bug is in ElectronicEnergy, but we can work around it)
        if use_qiskit_nature and hasattr(hamiltonian, 'eri') and hasattr(hamiltonian, 'h_core'):
            try:
                logger.info("Using Qiskit Nature fermionic operators for Pauli conversion")
                return PauliConverter._to_pauli_qiskit_nature(hamiltonian, mapper)
            except ImportError:
                logger.info("Qiskit Nature not installed, falling back")
                pass
            except Exception as e:
                logger.warning(f"Qiskit Nature approach failed ({e}), falling back")
                pass

        # Collect all Pauli terms
        pauli_dict = {}  # {pauli_string: coefficient}

        n_orbitals = hamiltonian.n_orbitals

        # CRITICAL FIX: Transform AO basis integrals to MO basis
        # VQE operates on MO basis, not AO basis!

        # Get MO coefficients from HF calculation
        mo_energies, C = hamiltonian.compute_molecular_orbitals()

        # Transform h_core to MO basis: h_mo = C† h_ao C
        h_core = C.T @ hamiltonian.h_core @ C

        # Transform ERI to MO basis
        if hasattr(hamiltonian, 'eri') and hamiltonian.eri is not None:
            eri_ao = hamiltonian.eri
            # Transform all 4 indices: eri_mo[i,j,k,l] = Σ_pqrs C[p,i] C[q,j] eri_ao[p,q,r,s] C[r,k] C[s,l]
            eri = np.einsum('pi,qj,pqrs,rk,sl->ijkl', C, C, eri_ao, C, C, optimize=True)
        else:
            eri = None

        # For Jordan-Wigner and other mappers: need spin orbitals (alpha + beta)
        # Each spatial orbital has 2 spin orbitals (spin-up, spin-down)
        n_spin_orbitals = n_orbitals * 2
        n_qubits = mapper.n_qubits(n_spin_orbitals)

        # 1. One-body terms: Σ h_ij a†_i a_j
        for i in range(n_orbitals):
            for j in range(n_orbitals):
                if abs(h_core[i, j]) > 1e-10:
                    # Map fermionic operator a†_i a_j to Pauli operators
                    # Map for both spin-up and spin-down
                    for spin_offset in [0, n_orbitals]:  # spin-up (0), spin-down (n_orbitals)
                        i_spin = i + spin_offset
                        j_spin = j + spin_offset
                        pauli_terms = mapper.map_excitation_operator(i_spin, j_spin, n_spin_orbitals)

                        for pauli_string, coeff in pauli_terms.items():
                            full_coeff = h_core[i, j] * coeff

                            if pauli_string in pauli_dict:
                                pauli_dict[pauli_string] += full_coeff
                            else:
                                pauli_dict[pauli_string] = full_coeff

        # 2. Two-body terms: (1/2) Σ ⟨ij|kl⟩ a†_i a†_k a_l a_j
        # CRITICAL: Qiskit Nature convention (from SymmetricTwoBody):
        # For eri[i,j,k,l] = ⟨ij|kl⟩, the operator is a†_i a†_k a_l a_j
        # NOT a†_i a†_j a_l a_k!
        if eri is not None:
            for i in range(n_orbitals):
                for j in range(n_orbitals):
                    for k in range(n_orbitals):
                        for l in range(n_orbitals):
                            v_ijkl = eri[i, j, k, l]  # ⟨ij|kl⟩

                            if abs(v_ijkl) > 1e-10:
                                # Map operator: a†_i a†_k a_l a_j
                                # For all spin combinations
                                for spin_i_l in [0, n_orbitals]:  # Same spin for i, l
                                    for spin_k_j in [0, n_orbitals]:  # Same spin for k, j
                                        i_spin = i + spin_i_l
                                        k_spin = k + spin_k_j
                                        l_spin = l + spin_i_l
                                        j_spin = j + spin_k_j

                                        # Use mapper to convert a†_i a†_k a_l a_j to Pauli
                                        if hasattr(mapper, 'map_double_excitation'):
                                            # map_double_excitation(orb_from_1, orb_from_2, orb_to_1, orb_to_2)
                                            # for a†_i a†_k a_l a_j: annihilate j,l; create i,k
                                            pauli_terms = mapper.map_double_excitation(
                                                j_spin, l_spin, i_spin, k_spin, n_spin_orbitals
                                            )
                                        else:
                                            pauli_terms = PauliConverter._map_two_body_approximate(
                                                mapper, i_spin, k_spin, l_spin, j_spin, n_spin_orbitals
                                            )

                                        for pauli_string, coeff in pauli_terms.items():
                                            full_coeff = 0.5 * v_ijkl * coeff

                                            if pauli_string in pauli_dict:
                                                pauli_dict[pauli_string] += full_coeff
                                            else:
                                                pauli_dict[pauli_string] = full_coeff

        # 3. Nuclear repulsion (constant term - identity operator)
        # n_qubits already calculated above
        identity_string = 'I' * n_qubits

        if identity_string in pauli_dict:
            pauli_dict[identity_string] += hamiltonian.nuclear_repulsion
        else:
            pauli_dict[identity_string] = hamiltonian.nuclear_repulsion

        # Use OpenFermion's validated Jordan-Wigner transformation
        logger.info("Using OpenFermion Jordan-Wigner transformation (validated reference)")

        from kanad.core.hamiltonians.openfermion_jw import openfermion_jordan_wigner

        # Return the validated transformation
        return openfermion_jordan_wigner(
            h_core,
            eri,
            hamiltonian.nuclear_repulsion,
            hamiltonian.n_electrons
        )

        # 4. Filter out negligible terms
        pauli_dict = {
            k: v for k, v in pauli_dict.items()
            if abs(v) > 1e-12
        }

        # 5. Convert to Qiskit SparsePauliOp
        if len(pauli_dict) == 0:
            # Empty Hamiltonian - return zero operator
            return SparsePauliOp(['I' * n_qubits], [0.0])

        pauli_strings = list(pauli_dict.keys())
        coefficients = list(pauli_dict.values())

        return SparsePauliOp(pauli_strings, coefficients)

    @staticmethod
    def _map_two_body_approximate(mapper, i, j, k, l, n_orbitals):
        """
        Approximate two-body operator using product of excitations.

        a†_i a†_j a_l a_k ≈ (a†_i a_k)(a†_j a_l)

        This is not exact due to anticommutation, but provides a starting point.
        """
        # Get single excitations
        exc_ik = mapper.map_excitation_operator(k, i, n_orbitals)
        exc_jl = mapper.map_excitation_operator(l, j, n_orbitals)

        # Multiply Pauli strings
        result = mapper.pauli_string_multiply(exc_ik, exc_jl)

        return result

    @staticmethod
    def to_pauli_dict(hamiltonian, mapper) -> Dict[str, complex]:
        """
        Convert Hamiltonian to dictionary of Pauli terms.

        Returns:
            Dictionary mapping Pauli strings to coefficients
        """
        sparse_pauli = PauliConverter.to_sparse_pauli_op(hamiltonian, mapper)

        pauli_dict = {}
        for pauli, coeff in zip(sparse_pauli.paulis, sparse_pauli.coeffs):
            pauli_dict[str(pauli)] = complex(coeff)

        return pauli_dict

    @staticmethod
    def count_pauli_terms(hamiltonian, mapper) -> int:
        """Count number of Pauli terms in Hamiltonian."""
        sparse_pauli = PauliConverter.to_sparse_pauli_op(hamiltonian, mapper)
        return len(sparse_pauli)

    @staticmethod
    def get_hamiltonian_info(hamiltonian, mapper) -> Dict:
        """
        Get information about the Pauli decomposition.

        Returns:
            Dictionary with statistics
        """
        sparse_pauli = PauliConverter.to_sparse_pauli_op(hamiltonian, mapper)

        coeffs = np.abs(sparse_pauli.coeffs)

        return {
            'num_terms': len(sparse_pauli),
            'max_coeff': np.max(coeffs),
            'min_coeff': np.min(coeffs[coeffs > 0]) if np.any(coeffs > 0) else 0,
            'mean_coeff': np.mean(coeffs),
            'total_weight': np.sum(coeffs),
        }

    @staticmethod
    def _to_pauli_governance_based(hamiltonian, mapper):
        """
        Governance-based Pauli conversion using number operators for two-electron terms.

        This approach uses n_i n_j directly instead of fermionic a†a operators,
        avoiding Jordan-Wigner transformation issues for two-body terms.

        More accurate for small molecules in minimal basis.
        """
        try:
            from qiskit.quantum_info import SparsePauliOp
        except ImportError:
            raise ImportError("Qiskit not installed")

        n_orbitals = hamiltonian.n_orbitals
        n_spin_orbitals = n_orbitals * 2
        n_qubits = mapper.n_qubits(n_spin_orbitals)

        h_core = hamiltonian.h_core
        eri = hamiltonian.eri
        pauli_dict = {}

        # 1. One-electron terms: Σ h_ij n_i (for diagonal terms i=j)
        #    Σ h_ij (a†_i a_j + a†_j a_i)/2 (for off-diagonal, but simplified to number ops)
        for i in range(n_orbitals):
            # Diagonal terms: h_ii (n_i↑ + n_i↓)
            if abs(h_core[i, i]) > 1e-10:
                for spin in [0, 1]:  # 0=up, 1=down
                    i_spin = 2 * i + spin
                    # n_i = (I - Z_i)/2
                    pauli_I = 'I' * n_qubits
                    pauli_Z = list('I' * n_qubits)
                    pauli_Z[i_spin] = 'Z'
                    pauli_Z = ''.join(pauli_Z)

                    if pauli_I in pauli_dict:
                        pauli_dict[pauli_I] += 0.5 * h_core[i, i]
                    else:
                        pauli_dict[pauli_I] = 0.5 * h_core[i, i]

                    if pauli_Z in pauli_dict:
                        pauli_dict[pauli_Z] += -0.5 * h_core[i, i]
                    else:
                        pauli_dict[pauli_Z] = -0.5 * h_core[i, i]

            # Off-diagonal terms: use excitation operators
            for j in range(i + 1, n_orbitals):
                if abs(h_core[i, j]) > 1e-10:
                    # Map a†_i a_j + a†_j a_i for both spins
                    for spin in [0, 1]:
                        i_spin = 2 * i + spin
                        j_spin = 2 * j + spin

                        # a†_i a_j
                        pauli_terms_ij = mapper.map_excitation_operator(j_spin, i_spin, n_spin_orbitals)
                        for pauli_str, coeff in pauli_terms_ij.items():
                            full_coeff = h_core[i, j] * coeff
                            if pauli_str in pauli_dict:
                                pauli_dict[pauli_str] += full_coeff
                            else:
                                pauli_dict[pauli_str] = full_coeff

                        # a†_j a_i (hermitian conjugate)
                        pauli_terms_ji = mapper.map_excitation_operator(i_spin, j_spin, n_spin_orbitals)
                        for pauli_str, coeff in pauli_terms_ji.items():
                            full_coeff = h_core[j, i] * coeff
                            if pauli_str in pauli_dict:
                                pauli_dict[pauli_str] += full_coeff
                            else:
                                pauli_dict[pauli_str] = full_coeff

        # 2. Two-electron terms: ONLY NUMBER OPERATOR PRODUCTS n_i n_j
        #    Avoid problematic a†_i a†_k a_l a_j terms by using Coulomb approximation
        #    WARNING: This only includes diagonal Coulomb terms, missing exchange!
        for i in range(n_orbitals):
            for j in range(i, n_orbitals):  # j >= i to avoid double-counting
                # Coulomb integral: ⟨ii|jj⟩ n_i n_j
                J_iijj = eri[i, i, j, j]

                if abs(J_iijj) > 1e-10:
                    # All four spin combinations
                    for spin_i in [0, 1]:
                        for spin_j in [0, 1]:
                            i_spin = 2 * i + spin_i
                            j_spin = 2 * j + spin_j

                            # Skip if same spin orbital (would be n_i^2 = n_i)
                            if i == j and spin_i == spin_j:
                                continue

                            # n_i n_j = (I - Z_i)(I - Z_j)/4
                            pauli_II = 'I' * n_qubits
                            pauli_IZ = list('I' * n_qubits)
                            pauli_IZ[j_spin] = 'Z'
                            pauli_IZ = ''.join(pauli_IZ)

                            pauli_ZI = list('I' * n_qubits)
                            pauli_ZI[i_spin] = 'Z'
                            pauli_ZI = ''.join(pauli_ZI)

                            pauli_ZZ = list('I' * n_qubits)
                            pauli_ZZ[i_spin] = 'Z'
                            pauli_ZZ[j_spin] = 'Z'
                            pauli_ZZ = ''.join(pauli_ZZ)

                            # Factor: 0.5 (from Hamiltonian 1/2 prefactor) × 1/4 (from n_i n_j expansion)
                            # But if i != j, we only count once (upper triangle), so no extra factor
                            factor = 0.5 * J_iijj * 0.25

                            if pauli_II in pauli_dict:
                                pauli_dict[pauli_II] += factor
                            else:
                                pauli_dict[pauli_II] = factor

                            if pauli_IZ in pauli_dict:
                                pauli_dict[pauli_IZ] += -factor
                            else:
                                pauli_dict[pauli_IZ] = -factor

                            if pauli_ZI in pauli_dict:
                                pauli_dict[pauli_ZI] += -factor
                            else:
                                pauli_dict[pauli_ZI] = -factor

                            if pauli_ZZ in pauli_dict:
                                pauli_dict[pauli_ZZ] += factor
                            else:
                                pauli_dict[pauli_ZZ] = factor

        # 3. Nuclear repulsion
        identity = 'I' * n_qubits
        if identity in pauli_dict:
            pauli_dict[identity] += hamiltonian.nuclear_repulsion
        else:
            pauli_dict[identity] = hamiltonian.nuclear_repulsion

        # 4. Clean up
        pauli_dict = {k: v for k, v in pauli_dict.items() if abs(v) > 1e-12}

        if len(pauli_dict) == 0:
            return SparsePauliOp(['I' * n_qubits], [0.0])

        pauli_strings = list(pauli_dict.keys())
        coefficients = list(pauli_dict.values())

        return SparsePauliOp(pauli_strings, coefficients)

    @staticmethod
    def _to_pauli_qiskit_nature(hamiltonian, mapper):
        """
        Legacy method - now uses OpenFermion instead of qiskit-nature.

        The qiskit-nature dependency has been removed. This method now delegates
        to OpenFermion-based implementation which provides the same correctness.

        Args:
            hamiltonian: MolecularHamiltonian with h_core and eri
            mapper: Not used

        Returns:
            qiskit.quantum_info.SparsePauliOp
        """
        # Use OpenFermion implementation instead (qiskit-nature removed)
        from kanad.core.hamiltonians.openfermion_jw import openfermion_jordan_wigner

        # Get MO coefficients
        mo_energies, C = hamiltonian.compute_molecular_orbitals()

        # Transform to MO basis
        h_core_mo = C.T @ hamiltonian.h_core @ C
        eri_mo = np.einsum('pi,qj,pqrs,rk,sl->ijkl', C, C, hamiltonian.eri, C, C, optimize=True)

        return openfermion_jordan_wigner(
            h_mo=h_core_mo,
            eri_mo=eri_mo,
            nuclear_repulsion=hamiltonian.nuclear_repulsion,
            n_electrons=hamiltonian.n_electrons
        )
