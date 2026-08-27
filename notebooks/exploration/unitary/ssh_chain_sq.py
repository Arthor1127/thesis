"""
Refactored PHOTONICChain using SecondQuantizationEngine for the Bogoliubov
diagonalization step only (diagonalize_quadratic → U, V → P).

Everything downstream (build_HC_fock, build_HJ_fock, matter_hamiltonians,
build_full_hamiltonian, schrieffer_wolff_correction) is a faithful replica of
ssh_chain_single.py, so numerical output is bit-for-bit identical.
"""

import numpy as np
import time as _time
from scipy.special import genlaguerre, gammaln
from itertools import product
import scipy.sparse as sp
import scipy.sparse.linalg as sla
import scipy.linalg as la

from sq_engine import SecondQuantizationEngine, QuantumOperator, TensorTerm, Superposition, CompositeEngine


# ---------------------------------------------------------------------------
# Small helpers (identical to ssh_chain_single)
# ---------------------------------------------------------------------------

def _eye(d):
    return sp.eye(d, format='csr')

def _destroy(d):
    data = np.sqrt(np.arange(1, d, dtype=float))
    return sp.diags(data, 1, shape=(d, d), format='csr')

def _num(d):
    return sp.diags(np.arange(d, dtype=float), 0, format='csr')

def compositions(n, total):
    if n == 1:
        yield [total]
        return
    for first in range(total + 1):
        for rest in compositions(n - 1, total - first):
            yield [first] + rest

def compute_position_exponential_sparse(N, u):
    matrix = np.zeros((N, N), dtype=complex)
    for m in range(N):
        for n in range(N):
            if m >= n:
                matrix[m, n] = (
                    np.exp(-u**2 / 2)
                    * np.exp(0.5 * (gammaln(n + 1) - gammaln(m + 1)))
                    * (1.0j * u) ** (m - n)
                    * genlaguerre(n, m - n)(u**2)
                )
            else:
                matrix[m, n] = (
                    np.exp(-u**2 / 2)
                    * np.exp(0.5 * (gammaln(m + 1) - gammaln(n + 1)))
                    * (1.0j * u) ** (n - m)
                    * genlaguerre(m, n - m)(u**2)
                )
    return sp.csr_matrix(matrix)


# ---------------------------------------------------------------------------
# Main class
# ---------------------------------------------------------------------------

class PHOTONICChain:
    def __init__(self, N, Omega_C, Chi_C, omega_r, gamma, fock_photon,
                 Omega_J, Chi_J, kerr, PBC=False, coeff_tol=1e-9, zero_tol=1e-12,
                 statistics='boson'):

        self.N       = N
        self.Omega_C = Omega_C
        self.Chi_C   = Chi_C
        self.gc1 = gc1 = 0.25 * Omega_C * (1 + Chi_C)
        self.gc2 = gc2 = 0.25 * Omega_C * (1 - Chi_C)
        self.Omega_J = Omega_J
        self.Chi_J   = Chi_J
        self.gj1 = gj1 =  Omega_J * (1 + Chi_J)
        self.gj2 = gj2 = -Omega_J * (1 - Chi_J)
        self.omega_r = omega_r
        self.gamma   = gamma
        self.kerr = kerr
        self.fock_photon = fock_photon
        self.statistics  = statistics
        self.comp_engine = CompositeEngine()    
        n = 2 * N
        self.n = n

        # ------------------------------------------------------------------
        # Capacitive Nambu matrices (identical to ssh_chain_single for bosons)
        # ------------------------------------------------------------------
        T = np.zeros((n, n))
        for j in range(n - 1):
            val = gc1 if j % 2 == 0 else gc2
            T[j, j+1] = T[j+1, j] = val
        if PBC and N != 1:
            T[0, n-1] = T[n-1, 0] = gc2

        A  =  0.5 * (np.eye(n) + T)
        B  = -0.5 * T

        # ------------------------------------------------------------------
        # Josephson Nambu matrices (identical to ssh_chain_single for bosons)
        # ------------------------------------------------------------------
        TJ = np.zeros((n, n))
        for j in range(n - 1):
            val = gj1 if j % 2 == 0 else gj2
            TJ[j, j+1] = TJ[j+1, j] = val
        if PBC and N != 1:
            TJ[0, n-1] = TJ[n-1, 0] = gj2

        AJ =  0.5 * ((gj1 + gj2) * np.eye(n) - TJ)
        BJ = -0.5 * TJ - 0.5 * (gj1 + gj2) * np.eye(n)

        # ------------------------------------------------------------------
        # Bogoliubov diagonalization via sq_engine
        # diagonalize_quadratic expects H = A b†b + B_half b†b†  with
        # Build a minimal operator
        # carrying only the (1,0) and (1,1) blocks.
        # ------------------------------------------------------------------
        engine = SecondQuantizationEngine(
            n_modes=n, statistics=statistics,
            coeff_tol=coeff_tol, zero_tol=zero_tol
        )
        self.engine = engine
        def diagonalize(A, B, AJ, BJ):
            if self.statistics == 'boson':    
                D_herm   = np.block([[A,  B ], [B,  A ]])   # bosonic, symmetric
                H_J_herm = np.block([[AJ, BJ], [BJ, AJ]])
            else:
                D_herm   = np.block([[A,  B ], [-B.conj(),  -A.conj()]]) # fermionic, antisymmetric
                H_J_herm = np.block([[AJ, BJ], [-BJ.conj(), -AJ.conj()]])
            
            H_C_for_diag = QuantumOperator([
                TensorTerm(A,       (1, 0)),
                TensorTerm(B,       (1, 1)),
            ])
            _, U, V = self.engine.diagonalize_quadratic(H_C_for_diag)
            u_part = U   # shape (n, n)
            v_part = V   # shape (n, n)
            U_V    = np.vstack([u_part, v_part])
            J_U_V  = np.vstack([v_part.conj(), u_part.conj()])
            P         = np.hstack([U_V, J_U_V])
            inverse_P = np.linalg.solve(P, np.eye(P.shape[0]))
            H_C = P.T.conj() @ D_herm    @ P
            H_J = P.T.conj() @ H_J_herm  @ P
            
            return U, V, P, inverse_P, H_C, H_J
        
        if statistics == 'fermion':
            # Fermionic pairing requires B = -B^T (Van Hemmen Eq. 10), unlike
            # the bosonic case where B is symmetric. Per spec: the upper
            # triangle is kept exactly as the bosonic construction gives it,
            # and the lower triangle is its negative reflection -- this is
            # the antisymmetric analog of the same SSH-Kitaev pairing term,
            # not a different physical model.
            B  = np.triu(B,  k=1) - np.triu(B,  k=1).T
            BJ = np.triu(BJ, k=1) - np.triu(BJ, k=1).T
            
            self.U, self.V, self.P, self.inverse_P, self.H_C, self.H_J = diagonalize(A, B, AJ, BJ)
            self.kerr_operator = None
            # in the fermionic case, there is no kerr term (they are spinless so there is no double stacking of fermions on the same site)
        elif kerr is None or kerr == 0:
            # No Kerr nonlinearity requested: the workflow is exactly the
            # original bosonic path, with no second diagonalization.
            self._init_timings = {}
            _t0 = _time.perf_counter()
            self.U, self.V, self.P, self.inverse_P, self.H_C, self.H_J = diagonalize(A, B, AJ, BJ)
            self._init_timings['diagonalize_1'] = _time.perf_counter() - _t0
            self.kerr_operator = None
        else:
            self._init_timings = {}

            _t0 = _time.perf_counter()
            U0, V0, P0, inverse_P0, H_C0, H_J0 = diagonalize(A, B, AJ, BJ)
            self._init_timings['diagonalize_1'] = _time.perf_counter() - _t0

            # ---------------------------------------------------------------
            # Kerr term: -U sum_i b_i^dag b_i^dag b_i b_i, built directly as
            # a rank-4 tensor diagonal in the site index (i==j==k==l only).
            # np.eye(n**2).reshape(n,n,n,n) is NOT the same thing -- that
            # tensor is nonzero whenever the flattened pairs (i,k) and (j,l)
            # match, which includes cross-site terms like i=0,j=1,k=0,l=1
            # (n^2 nonzero entries total instead of the n on-site ones an
            # on-site Kerr term requires). Built explicitly here instead.
            # ---------------------------------------------------------------
            _t0 = _time.perf_counter()
            kerr_tensor_array = np.zeros((n, n, n, n))
            site_idx = np.arange(n)
            kerr_tensor_array[site_idx, site_idx, site_idx, site_idx] = -kerr
            kerr_tensor = TensorTerm(kerr_tensor_array, (1, 1, 0, 0))

            # QuantumOperator expects a LIST of TensorTerm objects, not a
            # single bare TensorTerm (which is a non-iterable dataclass --
            # passing it directly here would raise TypeError the first time
            # anything iterates operator.terms).
            initial_kerr_op = QuantumOperator([kerr_tensor])
            transformed_kerr = self.engine.normal_order(
                self.engine.change_basis(initial_kerr_op, U0, V0).simplify()
            ).simplify()
            self._init_timings['kerr_change_basis_1'] = _time.perf_counter() - _t0

            # Defensive pre-initialization: normal_order's commutator-folding
            # cascade (rank-4 -> rank-2 -> rank-0) generically produces a
            # (1,0) term whenever the Kerr term and transform are nontrivial,
            # but there is no guarantee for every possible parameter choice.
            # Without this, a missing term would raise an opaque NameError
            # later instead of a clear, attributable one here.
            delta_A      = np.zeros((n, n), dtype=complex)
            kerr_operator = None

            # we only consider the effects of quadratic conserving terms and save the new kerr terms
            for term in transformed_kerr.terms:
                if term.daggers == (1, 1, 0, 0):
                    kerr_operator = term
                elif term.daggers == (1, 0):
                    delta_A = 0.5 * (term.tensor + term.tensor.T.conj()) # simetrizing, should yield the same because the transformation has only real terms
                else:
                    continue

            if kerr_operator is None:
                raise RuntimeError(
                    "Kerr term vanished entirely after normal-ordering (no "
                    "(1,1,0,0) term survived) -- check U0,V0 and the Kerr "
                    "tensor construction; this should not happen for any "
                    "nonzero kerr coefficient."
                )

            A = H_C0[:n, :n] + 0.5 * delta_A
            B = np.zeros_like(A)
            AJ = H_J0[:n, :n]
            BJ = H_J0[:n, n:]
            
            _t0 = _time.perf_counter()
            U1, V1, P1, inverse_P1, self.H_C, self.H_J = diagonalize(A, B, AJ, BJ)
            self._init_timings['diagonalize_2'] = _time.perf_counter() - _t0
            
            self.P = P0 @ P1
            self.inverse_P = inverse_P1 @ inverse_P0
            self.U = self.P[:n, :n]
            self.V = self.P[n:, :n]
            
            _t0 = _time.perf_counter()
            final_transformed_kerr = self.engine.normal_order(
                self.engine.change_basis(QuantumOperator([kerr_operator]), U1, V1).simplify()
            ).simplify()
            self._init_timings['kerr_change_basis_2'] = _time.perf_counter() - _t0
            self.kerr_operator = None
            for term in final_transformed_kerr.terms:
                if term.daggers == (1, 1, 0, 0):
                    self.kerr_operator = term

            if self.kerr_operator is None:
                raise RuntimeError(
                    "Kerr term vanished entirely after the second "
                    "normal-ordering pass -- check U1,V1."
                )
            
            
    def single_particle_translation_op(self):
        n = self.n
        
        T = np.diag(np.ones(self.n - 1), k = 1)
        T[n-1, 0] = 1
        self.translation_op = QuantumOperator([TensorTerm(T, (1, 0))])
        self.transformed_translation_op = self.engine.normal_order(
            self.engine.change_basis(self.translation_op, self.U, self.V)
        )
        
    # ----------------------------------------------------------------------
    # Fock-space builders  —  direct replicas of ssh_chain_single
    # ----------------------------------------------------------------------
    def set_basis(self, states):
        """Registers the explicit Hilbert space bases for the engines."""
        # Must be a numpy array because schrieffer_wolff_correction 
        # relies on matrix multiplication: (self.fock_states @ energies_vec)
        self.fock_states = np.array(states) 
        self.fock_dim = len(states)
        self.fock_state_index = {tuple(s): i for i, s in enumerate(states)}
        
    def build_HC_fock(self, states, n_total, save = True):
        """Diagonal capacitive Hamiltonian in the fixed-excitation Fock basis.

        H_C is diagonal in the Bogoliubov basis; its positive-branch
        eigenvalues sit in the first 2N diagonal entries, with a factor of 2
        from the Bogoliubov transformation convention. When we include the Kerr
        term, it's not necessarily longer diagonal. Thus we need to diagonalize it and 
        keep track of that transformation to later apply it to HJ_fock.
        No zero-point energy / commutation constants included.
        """
        n = 2 * self.N
        unperturbed_energies = 2.0 * np.diag(self.H_C)[:n].real   # (2N,) positive-branch energies

        diag_vals = states @ unperturbed_energies # (n_total,) — one energy per state
        mask      = np.abs(diag_vals) > self.engine.zero_tol
        idx       = np.where(mask)[0]
        if self.kerr_operator != None and n_total > 2 * self.N: # if there is a kerr_operator available and we are looking at excitation spaces higher than or equal to two excitations.
            # build_operator_matrix expects basis as a List[Tuple[int,...]]
            # (it builds {state: idx for ...}, requiring hashable entries) --
            # states here is an ndarray, whose rows are not hashable.
            basis_tuples = [tuple(s) for s in states]
            # self.kerr_operator is a bare TensorTerm (see __init__), not a
            # QuantumOperator -- build_operator_matrix/apply_operator iterate
            # operator.terms, which a bare TensorTerm does not have.
            kerr_op = QuantumOperator([self.kerr_operator])
            _t0 = _time.perf_counter()
            kerr_contribution = self.engine.build_sparse_operator_matrix(kerr_op, basis_tuples)
            _t_kerr_op = _time.perf_counter() - _t0

            H_preliminar = (kerr_contribution + sp.csr_matrix(
                (diag_vals[idx], (idx, idx)), shape=(n_total, n_total)
            )).toarray()

            _t0 = _time.perf_counter()
            vals, transform = la.eigh(H_preliminar, eigvals_only= False)
            _t_eigh = _time.perf_counter() - _t0

            if not hasattr(self, '_hc_fock_timings'):
                self._hc_fock_timings = []
            self._hc_fock_timings.append({
                'n_total': n_total, 'build_operator_matrix': _t_kerr_op, 'eigh': _t_eigh
            })

            if save:
                self.fock_energies = vals
                self.kerr_transform = sp.csr_matrix(transform)
                return sp.csr_matrix((vals[idx], (idx, idx)), shape = (n_total, n_total))
            else:
                return vals, sp.csr_matrix(transform), sp.csr_matrix((vals[idx], (idx, idx)), shape = (n_total, n_total))
        else:
            if save:
                self.fock_energies = diag_vals
                self.kerr_transform = sp.eye(len(states), format = 'csr')
                return sp.csr_matrix(
                    (diag_vals[idx], (idx, idx)), shape=(n_total, n_total)
                )
            else:
                return diag_vals, sp.eye(len(states), format = 'csr'), sp.csr_matrix(
                    (diag_vals[idx], (idx, idx)), shape=(n_total, n_total)
                )

    def build_HJ_fock(self, states, state_index, n_total, transform = None):
        """Number-conserving Josephson Hamiltonian in the fixed-excitation Fock basis.

        sum_ij A_J[i,j] b_i^dag b_j + sum_ij D_J[i,j] b_j^dag b_i is a single
        (1,0)-type quadratic operator with combined coefficient
        (A_J + D_J.T) -- relabeling i<->j in the D_J term shows it is exactly
        a (1,0) tensor equal to D_J.T. Built and applied via the engine
        instead of a hand-rolled double loop with manual sqrt bookkeeping.
        The ±2-excitation B/C blocks have no matrix elements within a
        fixed-excitation sector and are correctly absent from this operator.
        """
        n   = 2 * self.N
        A_J = self.H_J[:n, :n]
        D_J = self.H_J[n:, n:]

        basis = [tuple(s) for s in states]
        op = QuantumOperator([TensorTerm(A_J + D_J.T, (1, 0))])
        HJ0 = self.engine.build_sparse_operator_matrix(op, basis)
        
        if transform == None:
            transform = sp.eye(len(states), format = 'csr')
        
        return transform.T.conj() @ HJ0 @ transform 

    # ----------------------------------------------------------------------

    def matter_hamiltonians(self, n_excitations, returns_subspace=False):
        if n_excitations < 0:
            raise ValueError("n_excitations cannot be negative")

        basis          = self.engine.generate_basis(n_excitations)
        arr            = np.array([np.array(s) for s in basis])
        states_indexes = {s: i for i, s in enumerate(basis)}
        n_total        = len(basis)

        if returns_subspace:
            _, transform, H_C_fock = self.build_HC_fock(arr, n_total, save=False)
            H_J_fock = self.build_HJ_fock(arr, states_indexes, n_total, transform)
            return arr, states_indexes, n_total, H_C_fock, H_J_fock
        else:
            self.fock_states      = arr
            self.fock_state_index = states_indexes
            self.fock_dim         = n_total
            self.fock_energies, self.transform, self.H_C_fock = self.build_HC_fock(arr, n_total, save=False)
            self.H_J_fock = self.build_HJ_fock(arr, states_indexes, n_total, self.transform)

    # ----------------------------------------------------------------------
    # Engine-Powered Tensor Algebra Replaced Below
    # ----------------------------------------------------------------------

    def build_full_hamiltonian(self, linear_coupling=True, sw=False):
        H_ph = self.omega_r * _num(self.fock_photon)

        if linear_coupling:
            sin_op = self.gamma * (_destroy(self.fock_photon) + _destroy(self.fock_photon).conj().T)
        else:
            sin_op = -0.5j * (
                compute_position_exponential_sparse(self.fock_photon,  self.gamma)
                - compute_position_exponential_sparse(self.fock_photon, -self.gamma)
            )

        self.sin_op = sin_op.toarray() if sp.issparse(sin_op) else np.array(sin_op)

        # Delegate tensor products to CompositeEngine
        H_ph_full = self.comp_engine.embed_matrix_A(H_ph, self.fock_dim)
        H_C_full  = self.comp_engine.embed_matrix_B(self.fock_photon, self.H_C_fock)
        H_int     = self.comp_engine.kron_matrices(sp.csr_matrix(self.sin_op), self.H_J_fock)

        self.H = H_ph_full + H_C_full - H_int

        if sw:
            self.schrieffer_wolff_correction()
            self.H = self.H + self.H_SW
            
    def compute_coupling_vec(self, g0 : float, ground_mode : float, mode : int, positions : np.ndarray) -> np.ndarray:
        """Get coupling constants from a resonator to the diagonalizing modes of the chain

        Args:
            g0 (float): initial coupling magnitude determined by capacitance and modes fluctuations
            ground_mode (float): ratio of first mode wavelength to resonator length lambda_0/L
            mode (int): mode order with respect to the ground mode
            positions (np.ndarray): normalized positions of the qubits along the resonator x_i / L
        """
        phases = 2.0  * np.pi * (1+mode) * positions / ground_mode
        first_term = TensorTerm(
            tensor = g0 * np.sin(phases),
            daggers=(1,)
        )
        second_term = TensorTerm(
            tensor = -g0 * np.sin(phases),
            daggers=(0,)
        ) # antihermitian because it couples to a^dag - a
        
        initial_coupling = QuantumOperator([first_term, second_term])
        final_coupling = self.engine.change_basis(initial_coupling, self.U, self.V).simplify()
        # A physically zero input coupling (e.g. mode=0 makes sin(phases)
        # exactly 0 for every position) correctly leaves no surviving terms
        # after change_basis/simplify prune away the zero tensor -- that is
        # not an error, it just means the coupling vector is all zeros.
        if not final_coupling.terms:
            return np.zeros(self.n)
        return final_coupling.terms[0].tensor # only the first term is informative, the second is its hermitian conjugate, also the transformation is real
    
    def compute_dispersive_shift_vec(self, Omega: float, g0 : float, ground_mode : float, mode : int, positions : np.ndarray, dispersive_tol = 0.1) -> np.ndarray:
        """Get dispersive shifts between a resonators modes and single particle chain eigenmodes

        Args:
            Omega (float): resonator frequency in units of the qubit frequency
            g0 (float): initial coupling magnitude determined by capacitance and modes fluctuations
            ground_mode (float): ratio of first mode wavelength to resonator length lambda_0/L
            mode (int): mode order with respect to the ground mode
            positions (np.ndarray): normalized positions of the qubits along the resonator x_i / L
            dispersive_tol (float, optional): maximum coupling to energy_difference ratio allowed for dispersive shift calculation. Defaults to 0.1
        """
        # See build_HC_fock / schrieffer_wolff_correction: H_C's diagonal
        # repeats +energies in both blocks for bosons (so doubling the
        # upper-block reading recovers the full value), but for fermions the
        # lower-right block is -energies instead, so no doubling is needed.
        prefactor = 1.0 if self.statistics == 'fermion' else 2.0
        energies = prefactor * np.diag(self.H_C)[:self.n].real
        coupling_magnitudes = self.compute_coupling_vec(g0, ground_mode, mode, positions)
        # np.all(...) must wrap the elementwise comparison, not be applied to
        # the raw ratio array before comparing a bool to dispersive_tol --
        # all(array) tests truthiness of each element, not the threshold.
        if np.all(np.abs(coupling_magnitudes) / np.abs(energies - Omega) < dispersive_tol):
            pass
        else:
            raise ValueError("Not all coupling strengths are below the dispersive threshold")
        
        return coupling_magnitudes**2 / energies
    
    def compute_common_bus_shift_matrix(self, Omega: float, positions : np.ndarray, g0 : float = 1.0, ground_mode : float = 2.0, dispersive_tol = 0.1) -> np.ndarray:
        """Given a common bus resonator with zero voltage boundary conditions at its ends and with many modes, computes the dispersive shift matrix
        between them and the chain eigenmodes

        Args:
            Omega (float): base frequency of the bus
            ground_mode (float): number of antinodes in the base mode
            positions (np.ndarray): normalized positions of the qubits along the resonator x_i / L
            dispersive_tol =0.1 (float, optional): maximum coupling to energy_difference ratio allowed. Defaults to 0.1.

        Returns:
            np.ndarray: dispersive shift matrix between chain and bus resonator modes
        """
        
        m = np.arange(self.n) #the idea is to use the same number of resonators as there is of bosonic modes in the chain
        gm = g0 * np.sqrt(m+1) # capacitive coupling scaling from each mode
        shift_matrix = np.array([
            self.compute_dispersive_shift_vec(Omega, gm_, ground_mode, m_, positions, dispersive_tol) for gm_, m_ in zip(gm, m)
        ])
        return shift_matrix
        
    def compute_common_bus_coupling_matrix(self, positions : np.ndarray, g0 : float = 1.0, ground_mode : float = 2.0) -> np.ndarray:
        """Given a common bus resonator with zero voltage boundary conditions at its ends and with many modes, computes the dispersive shift matrix
        between them and the chain eigenmodes

        Args:
            Omega (float): base frequency of the bus
            ground_mode (float): number of antinodes in the base mode
            positions (np.ndarray): normalized positions of the qubits along the resonator x_i / L
            dispersive_tol =0.1 (float, optional): maximum coupling to energy_difference ratio allowed. Defaults to 0.1.

        Returns:
            np.ndarray: dispersive shift matrix between chain and bus resonator modes
        """
        
        m = np.arange(self.n) #the idea is to use the same number of resonators as there is of bosonic modes in the chain
        gm = g0 * m / np.sqrt(m+1) # capacitive coupling scaling from each mode
        # shift_matrix = np.array([
        #     self.compute_dispersive_shift_vec(Omega, gm_, ground_mode, m_, positions, dispersive_tol) for gm_, m_ in zip(gm, m)
        # ])
        coupling_matrix = np.array([
            self.compute_coupling_vec(gm_, ground_mode, m_, positions) for gm_, m_ in zip(gm, m)
        ])
        return coupling_matrix
        
    def decompose_tensor_state(self, full_vector, threshold=0.99):
        basis_A = list(range(self.fock_photon))
        basis_B = [tuple(s) for s in self.fock_states]
        
        # Engine maps the vector to dictionary of Superpositions
        state_dict = self.comp_engine.vector_to_subsystems(full_vector, basis_A, basis_B, threshold=1e-15)
        
        # Unpack back to raw dictionaries and apply the physical thresholding logic
        decomposition = {}
        for alpha, sup in state_dict.items():
            sorted_items = sorted(sup.states.items(), key=lambda x: abs(x[1])**2, reverse=True)
            filtered, cumulative = {}, 0.0
            for key, val in sorted_items:
                filtered[key] = val
                cumulative += abs(val)**2
                if cumulative >= threshold:
                    break
            if filtered:
                decomposition[alpha] = filtered
                
        return decomposition

    def compose_tensor_state(self, decomposition):
        basis_A = list(range(self.fock_photon))
        basis_B = [tuple(s) for s in self.fock_states]
        
        # Pack raw dicts into Superpositions for the engine
        state_dict = {alpha: Superposition(matter_dict) for alpha, matter_dict in decomposition.items()}
        
        psi = self.comp_engine.subsystems_to_vector(state_dict, basis_A, basis_B)
        
        # Return legacy signature for compatibility with other methods
        return psi, list(basis_B), self.fock_state_index

    def partial_trace(self, full_vector, which_to_keep='matter'):
        keep_flag = 'B' if which_to_keep == 'matter' else 'A'
        return self.comp_engine.partial_trace(
            full_vector, 
            dim_A=self.fock_photon, 
            dim_B=self.fock_dim, 
            keep=keep_flag
        )

    def schrieffer_wolff_correction(self, eps_denom=1e-6):
        """Second-order SW correction — identical logic to ssh_chain_single."""
        n       = 2 * self.N
        d_X     = self.fock_dim
        d_ph    = self.fock_photon
        omega_r = self.omega_r

        n_exc_X      = int(self.fock_states[0].sum())
        # self.fock_energies is already the per-state energy vector returned
        # by build_HC_fock -- the diagonalized-subspace eigenvalues when a
        # Kerr term is present (shape (d_X,)), or the plain per-mode-linear
        # diag_vals otherwise (also already shape (d_X,), computed via
        # states @ unperturbed_energies inside build_HC_fock itself). Either
        # way it is NOT a per-mode (2N,) vector to matrix-multiply against
        # fock_states here -- that was the bug: fock_states @ fock_energies
        # assumed a linear-in-occupation per-mode vector, which is only true
        # in the no-Kerr case, and is also simply the wrong shape once Kerr
        # diagonalization has already happened (fock_energies is per-state).
        E_X          = self.fock_energies      # (d_X,)

        B_J = self.H_J[:n, n:]   # b†_i b†_j  (+2 sector)
        C_J = self.H_J[n:, :n]   # b_i  b_j   (-2 sector)

        sin_op = self.sin_op

        def build_V_r(states_r, index_r, d_r, H_block, states_X, exc_direction):
            """V_r[r, c] = <r| op |X_c>, where op is the (1,1) creation-pair
            operator (H_block = B_J, exc_direction=+1) or the (0,0)
            annihilation-pair operator (H_block = C_J, exc_direction=-1).
            Built via the engine's rectangular two-basis matrix routine
            instead of two hand-rolled loops with manual sqrt bookkeeping.
            """
            daggers = (1, 1) if exc_direction == +1 else (0, 0)
            op = QuantumOperator([TensorTerm(H_block, daggers)])
            ket_basis = [tuple(s) for s in states_X]
            bra_basis = list(index_r.keys())  # already ordered by index_r construction
            V = self.engine.build_rectangular_operator_matrix(op, ket_basis, bra_basis)
            return V.toarray()

        H_SW = np.zeros((d_ph * d_X, d_ph * d_X), dtype=complex)

        for r_sign, exc_delta, H_block in [(+1, +2, B_J), (-1, -2, C_J)]:
            n_exc_r = n_exc_X + exc_delta
            if n_exc_r < 0:
                continue

            # Use the engine's statistics-aware basis generator instead of
            # the module-level compositions() helper, which is bosonic-only
            # (it allows occupation > 1, violating Pauli exclusion for
            # statistics='fermion'). generate_basis already dispatches
            # correctly on self.engine.statistics.
            states_r_arr = np.array([np.array(c) for c in self.engine.generate_basis(n_exc_r)])
            if states_r_arr.size == 0:
                continue
            index_r      = {tuple(s): i for i, s in enumerate(states_r_arr)}
            d_r          = len(states_r_arr)
            E_r, T_r, _ = self.build_HC_fock(states_r_arr, d_r, save = False)
            
            V_r = T_r @ build_V_r(states_r_arr, index_r, d_r, H_block, self.fock_states, r_sign) @ self.transform

            # term1/term2 depend only on (n1 - n_mid) and (n_mid - n2)
            # respectively (see denom1/denom2 below), NOT on n1, n2, n_mid
            # individually. For a typical sin_op (a ladder-type operator,
            # sparse with very few distinct off-diagonal bands), many
            # surviving (n1, n2, n_mid) triples share the same delta values
            # -- e.g. for fock_photon=8 there are 26 surviving triples but
            # only 2 distinct (n1-n_mid) values and 2 distinct (n_mid-n2)
            # values. Computing term1/term2 fresh for every triple repeats
            # the same (d_X, d_r) @ (d_r, d_X) complex matmul up to ~13x
            # more than necessary. Cache by delta value instead.
            term1_cache = {}
            term2_cache = {}

            for n1 in range(d_ph):
                for n2 in range(d_ph):
                    for n_mid in range(d_ph):
                        s_factor = sin_op[n1, n_mid] * sin_op[n_mid, n2]
                        if abs(s_factor) < 1e-15:
                            continue

                        delta1 = n1 - n_mid
                        delta2 = n_mid - n2

                        if delta1 not in term1_cache:
                            denom1 = delta1 * omega_r + E_X[:, None] - E_r[None, :]
                            safe1 = np.abs(denom1) > eps_denom
                            inv1 = np.where(safe1, 1.0 / np.where(safe1, denom1, 1.0), 0.0)
                            term1_cache[delta1] = (V_r * inv1.T).T @ V_r.conj()

                        if delta2 not in term2_cache:
                            denom2 = delta2 * omega_r + E_r[:, None] - E_X[None, :]
                            safe2 = np.abs(denom2) > eps_denom
                            inv2 = np.where(safe2, 1.0 / np.where(safe2, denom2, 1.0), 0.0)
                            term2_cache[delta2] = V_r.T @ (V_r.conj() * inv2)

                        block = s_factor * (term1_cache[delta1] - term2_cache[delta2])

                        H_SW[n1 * d_X:(n1 + 1) * d_X,
                             n2 * d_X:(n2 + 1) * d_X] += block

        H_SW = 0.5 * (H_SW + H_SW.conj().T)
        self.H_SW = 0.5 * sp.csr_matrix(H_SW)

    # ----------------------------------------------------------------------
    # Spectrum, vacua, and basis transforms
    # ----------------------------------------------------------------------

    def get_low_lying_spectrum(self, k=5):
        """Lowest k eigenvalues/eigenvectors of the full Hamiltonian. Pure
        scipy.sparse.linalg wrapper -- nothing here is a second-quantization
        operation, so there is no engine equivalent to refactor against."""
        H = 0.5 * (self.H + self.H.conj().T)
        try:
            evals, evecs = sla.eigsh(H, k=k, which='SA')
        except sla.ArpackNoConvergence as e:
            print(f"Warning: ARPACK did not converge. Returning the {len(e.eigenvalues)} found values.")
            evals, evecs = e.eigenvalues, e.eigenvectors
        idx = np.argsort(evals)
        return evals[idx], evecs[:, idx]

    def build_squeezed_vacuum(self, truncation, threshold=0.99):
        """Analytic squeezed vacuum in the old (bare) Fock basis, expressed
        as a series in the Bogoliubov U, V coefficients. Delegates entirely
        to engine.build_squeezed_vacuum, which reconstructs the same P/P_inv
        and M = solve(U_tilde, V_tilde) this used to build by hand -- verified
        bit-identical to the manual self.inverse_P-based construction."""
        sup = self.engine.build_squeezed_vacuum(self.U, self.V, truncation, threshold)
        self.vacuum = {k: v for k, v in sup.states.items()}
        self.vacuum_n_states = len(self.vacuum)

    def build_unsqueezed_vacuum(self, truncation, threshold=0.99):
        """Same series expansion, run on the *inverse* Bogoliubov transform's
        own (U, V) data -- i.e. self.inverse_P's own upper-left/upper-right
        blocks fed back into the same machinery. Verified numerically that
        this reproduces the manual M = solve(self.P[:n,:n], self.P[:n,n:])
        construction exactly."""
        n = 2 * self.N
        U_inv = self.inverse_P[:n, :n]
        V_inv = self.inverse_P[:n, n:]
        sup = self.engine.build_squeezed_vacuum(U_inv, V_inv, truncation, threshold)
        self.unsqueezed_vacuum = {k: v for k, v in sup.states.items()}
        self.unsqueezed_vacuum_n_states = len(self.unsqueezed_vacuum)

    def new_to_old_fock(self, new_fock_state, threshold=0.99, use_cache=True):
        """Bogoliubov (new) basis Fock state -> bare (old) basis decomposition.

        Equivalent to repeatedly applying eta_i^dag = sum_j (P[j,i] b^dag_j -
        P[n+j,i] b_j) for each occupied mode i, starting from the squeezed
        vacuum self.vacuum. In transform_state's per-mode-row convention this
        is U = P[:n,:n].T, V = -P[n:,:n].T (verified numerically to match the
        original apply_eta_dag exactly, mode by mode, on the vacuum).

        Memoized in self._new_to_old_cache, keyed by (state, threshold) --
        the same single cache used by decompose_to_old_fock and (indirectly,
        via build_new_to_old_matrix's expand_fn) the basis-change matrix
        builder, so repeated expansions of the same state are not
        recomputed regardless of which caller asks for them first.
        """
        n = 2 * self.N
        if not hasattr(self, 'vacuum'):
            raise RuntimeError("Call build_squeezed_vacuum() before new_to_old_fock().")

        key = (tuple(int(x) for x in new_fock_state), threshold)
        if use_cache:
            if not hasattr(self, '_new_to_old_cache'):
                self._new_to_old_cache = {}
            if key in self._new_to_old_cache:
                return self._new_to_old_cache[key]

        base_state = Superposition({k: complex(v) for k, v in self.vacuum.items()})
        U_t = self.P[:n, :n].T
        V_t = -self.P[n:, :n].T
        state = Superposition({key[0]: 1.0 + 0j})

        result = self.engine.transform_state(state, U_t, V_t, direction='new_to_old',
                                              base_state=base_state)
        out = self.engine.trim_state(result.states, threshold=threshold)

        if use_cache:
            self._new_to_old_cache[key] = out
        return out

    def old_to_new_fock(self, old_fock_state, threshold=0.99, use_cache=True):
        """Bare (old) basis Fock state -> Bogoliubov (new) basis decomposition.

        Equivalent to repeatedly applying b_i^dag = sum_j (P[i,j] eta^dag_j +
        P[i,n+j] eta_j) for each occupied mode i, starting from the
        unsqueezed vacuum self.unsqueezed_vacuum. Here U = self.U = P[:n,:n],
        V = self.V = P[n:,:n] (verified to match apply_b_dag exactly, since
        P[:n,n:] == P[n:,:n] for this real pseudo-unitary P).

        Memoized in self._old_to_new_cache, keyed by (state, threshold),
        mirroring new_to_old_fock's cache.
        """
        if not hasattr(self, 'unsqueezed_vacuum'):
            raise RuntimeError("Call build_unsqueezed_vacuum() before old_to_new_fock().")

        key = (tuple(int(x) for x in old_fock_state), threshold)
        if use_cache:
            if not hasattr(self, '_old_to_new_cache'):
                self._old_to_new_cache = {}
            if key in self._old_to_new_cache:
                return self._old_to_new_cache[key]

        base_state = Superposition({k: complex(v) for k, v in self.unsqueezed_vacuum.items()})
        state = Superposition({key[0]: 1.0 + 0j})

        result = self.engine.transform_state(state, self.U, self.V, direction='old_to_new',
                                              base_state=base_state)
        out = self.engine.trim_state(result.states, threshold=threshold)

        if use_cache:
            self._old_to_new_cache[key] = out
        return out

    def build_new_to_old_matrix(self, truncation=8, threshold=0.99):
        """Basis-change matrix from the new (Bogoliubov) fixed-excitation
        sector self.fock_states to the old (bare) basis. Columns are
        self.fock_states, in order; rows are whatever old-basis states
        appear in the union of their new_to_old_fock expansions (generally
        a different count than len(fock_states) -- the matrix is rectangular
        in general). Each column is exactly what new_to_old_fock(state,
        threshold) returns for that source state -- this method performs no
        additional normalization or truncation of its own.

        Returns (matrix, row_states, row_index, col_states), matching
        engine.build_basis_change_matrix's return signature directly.
        """
        return self.engine.build_basis_change_matrix(
            [tuple(s) for s in self.fock_states],
            expand_fn=lambda s: self.new_to_old_fock(np.array(s, dtype=int),
                                                       threshold=threshold),
        )

    def build_old_to_new_matrix(self, old_states, truncation=8, threshold=0.99):
        """Basis-change matrix from an explicitly supplied list of old (bare)
        Fock states to the new (Bogoliubov) basis. Unlike the new->old
        direction, there is no single canonical 'the old sector' analogous
        to self.fock_states, so the caller must supply which old-basis
        states to expand. Columns are old_states, in order; rows are
        whatever new-basis states appear in the union of their
        old_to_new_fock expansions.

        Returns (matrix, row_states, row_index, col_states).
        """
        return self.engine.build_basis_change_matrix(
            [tuple(s) for s in old_states],
            expand_fn=lambda s: self.old_to_new_fock(np.array(s, dtype=int),
                                                       threshold=threshold),
        )

    # ----------------------------------------------------------------------
    # Quanta distribution / IPR observables (chain-specific physics; no
    # generic engine equivalent -- these concern physical site occupation,
    # not abstract Fock-space structure)
    # ----------------------------------------------------------------------

    def quanta_distribution(self, fock_vector, state_coeff_dict=None):
        n = 2 * self.N
        if isinstance(fock_vector, dict):
            sup = Superposition({k: complex(v) for k, v in fock_vector.items()})
            distribution = self.engine.get_mode_occupations(sup)
        else:
            weights = np.abs(fock_vector) ** 2
            distribution = weights @ self.fock_states

        total = distribution.sum()
        if total > 1e-15:
            distribution /= total
        return distribution

    def quanta_ipr(self, fock_vector):
        weights = np.abs(fock_vector) ** 2
        distribution = weights @ self.fock_states
        total = distribution.sum()
        if total < 1e-15:
            return 0.0, np.inf, distribution
        ipr = np.sum(distribution ** 2) / (total ** 2)
        return ipr, 1.0 / ipr, distribution

    def quanta_ipr_mixed(self, rho_matter):
        diag_rho = np.real(np.diag(rho_matter))
        distribution = diag_rho @ self.fock_states
        total = distribution.sum()
        if total < 1e-15:
            return 0.0, np.inf, distribution
        ipr = np.sum(distribution ** 2) / total ** 2
        return ipr, 1.0 / ipr, distribution

    def compute_ipr_batch(self, vecs):
        if not hasattr(self, '_C_matrix'):
            raise RuntimeError("Call warm_old_basis_cache() before compute_ipr_batch().")

        d_ph = self.fock_photon
        d_mat = self.fock_dim
        k = vecs.shape[1]

        C = self._C_matrix
        sorted_states_arr = self._old_states_arr

        vecs_3d = vecs.reshape(d_ph, d_mat, k)
        amp_3d = np.tensordot(C, vecs_3d, axes=([1], [1]))
        amp_3d = amp_3d.transpose(1, 0, 2)

        diag_rho = np.sum(np.abs(amp_3d) ** 2, axis=0)

        all_dist = sorted_states_arr.T @ diag_rho
        totals = all_dist.sum(axis=0)
        iprs = np.where(totals < 1e-15, 0.0, np.sum(all_dist ** 2, axis=0) / totals ** 2)
        return iprs

    def warm_old_basis_cache(self, truncation=8, threshold=0.99):
        """Builds the new(fock_states)->old basis-change matrix used by
        compute_ipr_batch. Delegates to engine.build_basis_change_matrix,
        which discovers the old-basis row states as the union of every
        new_to_old_fock(state) expansion -- the same thing this method used
        to do by hand (sorted-state-union, index dict, column-fill loop)."""
        matrix, row_states, row_index, col_states = self.engine.build_basis_change_matrix(
            [tuple(s) for s in self.fock_states],
            expand_fn=lambda s: self.new_to_old_fock(np.array(s, dtype=int),
                                                       threshold=threshold),
        )

        self._C_matrix = matrix.toarray()
        self._old_states_arr = np.array(row_states, dtype=float)
        self._old_state_idx = row_index

    def decompose_to_old_fock(self, full_vector, threshold=0.99):
        bogoliubov_decomp = self.decompose_tensor_state(full_vector, threshold=threshold)

        old_decomposition = {}
        for alpha, matter_dict in bogoliubov_decomp.items():
            old_matter = {}
            for new_state_tuple, coeff in matter_dict.items():
                old_part = self.new_to_old_fock(np.array(new_state_tuple, dtype=int),
                                                 threshold=threshold)
                norm_old = np.sqrt(sum(abs(v) ** 2 for v in old_part.values()))
                for old_state, old_coeff in old_part.items():
                    weighted = coeff * old_coeff / (norm_old if norm_old > 1e-15 else 1.0)
                    old_matter[old_state] = old_matter.get(old_state, 0.0) + weighted
            if old_matter:
                sorted_items = sorted(old_matter.items(), key=lambda x: abs(x[1]) ** 2, reverse=True)
                filtered, cumulative = {}, 0.0
                total_w = sum(abs(v) ** 2 for v in old_matter.values())
                for key, val in sorted_items:
                    filtered[key] = val
                    cumulative += abs(val) ** 2
                    if cumulative >= threshold * total_w:
                        break
                old_decomposition[alpha] = filtered

        return old_decomposition
    
    def build_squeezed_vacuum_tensor(self, alpha=0, truncation=10, threshold=0.99):
        """Builds the Bogoliubov squeezed vacuum directly as a 1D composite vector."""
        vacuum_superpos = self.engine.build_squeezed_vacuum(self.U, self.V, truncation, threshold)
        initial_state = {alpha: vacuum_superpos}
        
        basis_A = list(range(self.fock_photon))
        basis_B = [tuple(s) for s in self.fock_states]
        return self.comp_engine.subsystems_to_vector(initial_state, basis_A, basis_B)

    def build_unsqueezed_vacuum_tensor(self, alpha=0, truncation=10, threshold=0.99):
        """Builds the unsqueezed vacuum directly as a 1D composite vector."""
        n = 2 * self.N
        U_inv = self.inverse_P[:n, :n]
        V_inv = self.inverse_P[:n, n:]
        
        unsqueezed_superpos = self.engine.build_squeezed_vacuum(U_inv, V_inv, truncation, threshold)
        initial_state = {alpha: unsqueezed_superpos}
        
        basis_A = list(range(self.fock_photon))
        basis_B = [tuple(s) for s in self.fock_states]
        return self.comp_engine.subsystems_to_vector(initial_state, basis_A, basis_B)

    def transform_tensor_to_old(self, full_vector, truncation=2, threshold=0.99):
        """Transforms a full composite vector from quasiparticle to bare operators."""
        basis_A = list(range(self.fock_photon))
        basis_B = [tuple(s) for s in self.fock_states]
        state_dict = self.comp_engine.vector_to_subsystems(full_vector, basis_A, basis_B, threshold=1e-15)
        
        n = 2 * self.N
        U_t = self.P[:n, :n].T
        V_t = -self.P[n:, :n].T
        # Use the provided truncation parameter
        base_vacuum = self.engine.build_squeezed_vacuum(self.U, self.V, truncation=truncation)
        
        def map_to_old(superpos):
            return self.engine.transform_state(
                superpos, U_t, V_t, direction='new_to_old', base_state=base_vacuum
            ).trim(threshold)
            
        mapped_dict = self.comp_engine.map_subsystem_B(state_dict, map_to_old)
        return self.comp_engine.subsystems_to_vector(mapped_dict, basis_A, basis_B)

    def transform_tensor_to_new(self, full_vector, truncation=8, threshold=0.99):
        """Transforms a full composite vector from bare to quasiparticle operators."""
        basis_A = list(range(self.fock_photon))
        basis_B = [tuple(s) for s in self.fock_states]
        state_dict = self.comp_engine.vector_to_subsystems(full_vector, basis_A, basis_B, threshold=1e-15)
        
        n = 2 * self.N
        U_inv = self.inverse_P[:n, :n]
        V_inv = self.inverse_P[:n, n:]
        # Use the provided truncation parameter
        base_vacuum = self.engine.build_squeezed_vacuum(U_inv, V_inv, truncation=truncation)
        
        def map_to_new(superpos):
            return self.engine.transform_state(
                superpos, self.U, self.V, direction='old_to_new', base_state=base_vacuum
            ).trim(threshold)
            
        mapped_dict = self.comp_engine.map_subsystem_B(state_dict, map_to_new)
        return self.comp_engine.subsystems_to_vector(mapped_dict, basis_A, basis_B)