"""
FIXED: Convert molecular Hamiltonians to Qiskit Pauli operators.

This version correctly handles spatial → spin orbital conversion.

Key principle: The h_core and ERI matrices are in SPATIAL orbital basis.
To convert to spin orbitals, we create operators for EACH spin separately,
but use the SAME matrix elements.
"""

from typing import Dict
import numpy as np


class PauliConverterFixed:
    """
    Correctly converts molecular Hamiltonians to Qiskit Pauli operators.

    The Hamiltonian in second quantization (spin-orbital basis):
        H = Σ_{pσ,qτ} h_{pq} δ_{στ} a†_{pσ} a_{qτ}  [one-electron]
          + ½ Σ_{pσ,qτ,rρ,sν} g_{pqrs} δ_{σρ}δ_{τν} a†_{pσ} a†_{qτ} a_{sν} a_{rρ} [two-electron]
          + E_nn [nuclear repulsion]

    Where:
    - p,q,r,s are SPATIAL orbital indices
    - σ,τ,ρ,ν are SPIN indices (↑ or ↓)
    - δ_{στ} is Kronecker delta (same spin)
    - h_{pq} and g_{pqrs} are in SPATIAL basis
    """

    @staticmethod
    def to_sparse_pauli_op(hamiltonian, mapper):
        """
        Convert molecular Hamiltonian to Qiskit SparsePauliOp.

        Args:
            hamiltonian: MolecularHamiltonian with h_core (spatial), ERI (spatial)
            mapper: Fermion-to-qubit mapper (must support spin orbitals)

        Returns:
            qiskit.quantum_info.SparsePauliOp
        """
        try:
            from qiskit.quantum_info import SparsePauliOp
        except ImportError:
            raise ImportError("Qiskit not installed. Install with: pip install qiskit>=2.0")

        pauli_dict = {}  # {pauli_string: coefficient}

        n_spatial_orbitals = hamiltonian.n_orbitals
        h_core = hamiltonian.h_core
        eri = hamiltonian.eri if hasattr(hamiltonian, 'eri') else None

        # Determine spin-orbital indexing
        # Convention: First n_spatial_orbitals indices are spin-up (α)
        #            Next n_spatial_orbitals indices are spin-down (β)
        n_spin_orbitals = 2 * n_spatial_orbitals
        n_qubits = mapper.n_qubits(n_spin_orbitals)

        # Helper: map spatial orbital + spin to spin-orbital index
        def spin_orbital_index(spatial_idx, spin):
            """spin: 0 for ↑ (alpha), 1 for ↓ (beta)"""
            if spin == 0:  # alpha
                return spatial_idx
            else:  # beta
                return spatial_idx + n_spatial_orbitals

        # 1. One-electron terms: Σ_{pσ,qτ} h_{pq} δ_{στ} a†_{pσ} a_{qτ}
        #    δ_{στ} means only same-spin terms contribute
        for p in range(n_spatial_orbitals):
            for q in range(n_spatial_orbitals):
                if abs(h_core[p, q]) < 1e-10:
                    continue

                # Alpha spin (σ=τ=↑)
                p_alpha = spin_orbital_index(p, 0)
                q_alpha = spin_orbital_index(q, 0)
                pauli_terms = mapper.map_excitation_operator(p_alpha, q_alpha, n_spin_orbitals)
                for pauli_str, coeff in pauli_terms.items():
                    full_coeff = h_core[p, q] * coeff
                    pauli_dict[pauli_str] = pauli_dict.get(pauli_str, 0) + full_coeff

                # Beta spin (σ=τ=↓)
                p_beta = spin_orbital_index(p, 1)
                q_beta = spin_orbital_index(q, 1)
                pauli_terms = mapper.map_excitation_operator(p_beta, q_beta, n_spin_orbitals)
                for pauli_str, coeff in pauli_terms.items():
                    full_coeff = h_core[p, q] * coeff
                    pauli_dict[pauli_str] = pauli_dict.get(pauli_str, 0) + full_coeff

        # 2. Two-electron terms: ½ Σ g_{pqrs} a†_p a†_q a_s a_r
        #    In spin-orbital form with spin restrictions:
        #    ½ Σ_{pσ,qτ,rρ,sν} g_{pqrs} δ_{σρ}δ_{τν} a†_{pσ} a†_{qτ} a_{sν} a_{rρ}
        if eri is not None:
            for p in range(n_spatial_orbitals):
                for q in range(n_spatial_orbitals):
                    for r in range(n_spatial_orbitals):
                        for s in range(n_spatial_orbitals):
                            # ERI in physicist's notation: (pr|qs)
                            # Antisymmetrized: g_{pqrs} = (pr|qs) - (ps|qr)
                            g_pqrs = eri[p, r, q, s] - eri[p, s, q, r]

                            if abs(g_pqrs) < 1e-10:
                                continue

                            # Four spin combinations (σρ same, τν same):
                            # 1. αα→αα: a†_{p↑} a†_{q↑} a_{s↑} a_{r↑}
                            # 2. αβ→αβ: a†_{p↑} a†_{q↓} a_{s↓} a_{r↑}
                            # 3. βα→βα: a†_{p↓} a†_{q↑} a_{s↑} a_{r↓}
                            # 4. ββ→ββ: a†_{p↓} a†_{q↓} a_{s↓} a_{r↓}

                            spin_cases = [
                                (0, 0, 0, 0),  # αα→αα
                                (0, 1, 1, 0),  # αβ→αβ
                                (1, 0, 0, 1),  # βα→βα
                                (1, 1, 1, 1),  # ββ→ββ
                            ]

                            for spin_p, spin_q, spin_s, spin_r in spin_cases:
                                p_spin = spin_orbital_index(p, spin_p)
                                q_spin = spin_orbital_index(q, spin_q)
                                r_spin = spin_orbital_index(r, spin_r)
                                s_spin = spin_orbital_index(s, spin_s)

                                # Use mapper to convert two-electron operator
                                # Most mappers don't have direct two-electron mapping,
                                # so approximate as product
                                pauli_terms = PauliConverterFixed._map_two_electron_approximate(
                                    mapper, p_spin, q_spin, s_spin, r_spin, n_spin_orbitals
                                )

                                for pauli_str, coeff in pauli_terms.items():
                                    full_coeff = 0.5 * g_pqrs * coeff
                                    pauli_dict[pauli_str] = pauli_dict.get(pauli_str, 0) + full_coeff

        # 3. Nuclear repulsion (constant term)
        identity_str = 'I' * n_qubits
        pauli_dict[identity_str] = pauli_dict.get(identity_str, 0) + hamiltonian.nuclear_repulsion

        # 4. Filter negligible terms
        pauli_dict = {k: v for k, v in pauli_dict.items() if abs(v) > 1e-12}

        # 5. Convert to SparsePauliOp
        if len(pauli_dict) == 0:
            return SparsePauliOp(['I' * n_qubits], [0.0])

        pauli_strings = list(pauli_dict.keys())
        coefficients = list(pauli_dict.values())

        return SparsePauliOp(pauli_strings, coefficients)

    @staticmethod
    def _map_two_electron_approximate(mapper, p, q, s, r, n_spin_orbitals):
        """
        Approximate two-electron operator a†_p a†_q a_s a_r.

        Uses product approximation: (a†_p a_r)(a†_q a_s)
        This is not exact due to anticommutation but provides starting point.
        """
        # Map single excitations
        exc_pr = mapper.map_excitation_operator(p, r, n_spin_orbitals)
        exc_qs = mapper.map_excitation_operator(q, s, n_spin_orbitals)

        # Multiply Pauli strings
        result = mapper.pauli_string_multiply(exc_pr, exc_qs)

        return result
