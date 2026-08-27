"""
Clase de la cadena ssh acoplada a un modo de cavidad
La idea es que se diagonaliza el hamiltoniano en espacio real realizando una tranformacion de bogoliubov,
El hamiltoniano Josephson resultante contiene elementos que conservan el número de partículas y otros que lo cambian en +/-2
Se aplica una transformación canónica sobre todo el hamiltoniano para eliminar la parte anómala del hamiltoniano Josephson a orden cuadrático.

La idea es entonces qudarse con la parte del hamiltoniano truncada a un número fijo de excitaciones habiendo introducido la corrección inducida por la
transformación canónica a posteriori.
"""

import numpy as np
from scipy.special import genlaguerre, gammaln
from itertools import product
import scipy.sparse as sp
import scipy.sparse.linalg as sla
import scipy.linalg as la

def ssh_energies(band, omega_q, gc1, gc2, k):
    delta = np.sqrt(gc1**2 + gc2**2 + 2 * gc1 * gc2 * np.cos(k))
    if band == 0:
        return omega_q * np.sqrt(1 + 2 * delta / omega_q)
    else:
        return omega_q * np.sqrt(1 - 2 * delta / omega_q)
    
    
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
    
class PHOTONICChain:
    """Clase general para la cadena de transmons acoplada al resonador
    Se trabaja en la representación adimensionalizada y se puede elegir si se trabaja con OBC o PBC
    
    La idea es diagonalizar la matriz dinámica asociada a la parte capacitiva del hamiltoniano y encontrar la base que lo diagonaliza.
    
    """
    def __init__(self, N, Omega_C, Chi_C, omega_r, gamma, fock_photon, Omega_J, Chi_J, PBC = False):
        self.N       = N
        self.Omega_C = Omega_C
        self.Chi_C   = Chi_C
        self.gc1 = gc1 = 0.25 * Omega_C * (1 + Chi_C)
        self.gc2 = gc2 = 0.25 * Omega_C * (1 - Chi_C)
        self.Omega_J = Omega_J
        self.Chi_J   = Chi_J
        self.gj1 = gj1 =  Omega_J * (1 + Chi_J) # omega_J no tiene un factor 4 como lo tiene Omega_c
        self.gj2 = gj2 = -Omega_J * (1 - Chi_J)
        self.omega_r = omega_r
        self.gamma = gamma
        self.fock_photon = fock_photon

        n = 2 * N
        T = np.zeros((n, n))
        for j in range(n - 1):
            val = gc1 if j % 2 == 0 else gc2
            T[j, j+1] = T[j+1, j] = val
        if PBC and N != 1:
            T[0, n-1] = T[n-1, 0] = gc2

        A =  0.5 * (np.eye(n) + T)
        B = -0.5 * T

        tau_z = np.block([
            [ np.eye(n),         np.zeros((n, n))],
            [ np.zeros((n, n)), -np.eye(n)       ]
        ])

        # Van Hemmen: diagonalize the dynamical matrix I_hat @ D = tau_z @ D_herm
        D_herm = np.block([[A, B],
                   [B, A]])  # real symmetric, keep it real

        L = np.linalg.cholesky(D_herm)          # D_herm = L @ L.T
        M_sym = L.T @ tau_z @ L                 # real symmetric

        evals_sym, evecs_sym = la.eigh(M_sym)   # real evals and evecs, sorted

        evecs_real = la.solve_triangular(L.T, evecs_sym)  # back to x = (L.T)^{-1} y

        eigenvals = evals_sym  # already real

        # Keep positive branch
        pos_idx   = eigenvals > 1e-15
        pos_evals = eigenvals[pos_idx]
        pos_evecs = evecs_real[:, pos_idx]

        sort_idx  = np.argsort(pos_evals)
        pos_evals = pos_evals[sort_idx]
        pos_evecs = pos_evecs[:, sort_idx]

        # Symplectic Gram-Schmidt — stays real since inputs are real
        def symplectic_gram_schmidt(vecs, tau_z):
            n_cols = vecs.shape[1]
            result = np.zeros_like(vecs, dtype=float)   # <-- float, not complex
            for i in range(n_cols):
                v = vecs[:, i].copy()
                for j in range(i):
                    e_j = result[:, j]
                    proj = e_j @ tau_z @ v              # no .conj() needed
                    v    = v - proj * e_j
                norm = v @ tau_z @ v                    # should be positive
                result[:, i] = v / np.sqrt(norm)
            return result

        U_V   = symplectic_gram_schmidt(pos_evecs, tau_z)

        # Gauge fixing: LAPACK only guarantees each eigenvector up to an
        # overall sign, not a canonical choice. Force the largest-magnitude
        # component of each column to be positive so the result is
        # reproducible across calls/LAPACK routines, and matches the
        # convention used by SecondQuantizationEngine.diagonalize_quadratic.
        for col in range(U_V.shape[1]):
            v = U_V[:, col]
            idx = np.argmax(np.abs(v))
            if v[idx] < 0:
                U_V[:, col] = -v

        u_part = U_V[:n, :]
        v_part = U_V[n:, :]
        J_U_V  = np.vstack([v_part, u_part])           # J acts as swap on real vecs

        self.P = np.hstack([U_V, J_U_V])               # real pseudo-unitary matrix
        self.inverse_P = np.linalg.solve(self.P, np.eye(self.P.shape[0]))
        # self.inverse_P = np.linalg.inv(self.P)
        # Transform both Hermitian Hamiltonians into the Bogoliubov basis
        # H_C should come out diagonal; H_J generally will not
        self.H_C = self.P.T.conj() @ D_herm       @ self.P

        TJ = np.zeros((n, n))
        for j in range(n - 1):
            val = gj1 if j % 2 == 0 else gj2
            TJ[j, j+1] = TJ[j+1, j] = val
        if PBC and N != 1:
            TJ[0, n-1] = TJ[n-1, 0] = gj2

        AJ =  0.5 * ((gj1 + gj2) * np.eye(n) - TJ)
        BJ = -0.5 * TJ - 0.5 * (gj1 + gj2) * np.eye(n)

        H_J_herm  = np.block([[AJ, BJ],
                            [BJ, AJ]])
        self.H_J  = self.P.T.conj() @ H_J_herm @ self.P
    
    def build_HC_fock(self, states, n_total):
        n = 2 * self.N
        energies = np.diag(self.H_C)[:n]

        # states is a flat ndarray (n_total, 2N); H_C is diagonal so this is just a matmul
        diag_vals = states @ energies          # (n_total,)
        mask      = np.abs(diag_vals) > 1e-15
        idx       = np.where(mask)[0]

        return 2.0 * sp.csr_matrix(
            (diag_vals[idx], (idx, idx)), shape=(n_total, n_total)
        )
    
    def build_HJ_fock(self, states, state_index, n_total):
        n   = 2 * self.N
        A_J = self.H_J[:n, :n]   # b†_i b_j  (number-conserving)
        D_J = self.H_J[n:, n:]   # b_i b†_j  (number-conserving, delta dropped)
        # B and C (±2 excitation pairing) are dead with a single fixed sector.

        data_dict = {}

        def add(row, col, val):
            if abs(val) > 1e-9:
                if (row, col) in data_dict:
                    data_dict[(row, col)] += val
                else:
                    data_dict[(row, col)]  = val

        for s_idx, state in enumerate(states):
            global_col = s_idx

            # A_J[i,j]: b†_i b_j |state>
            for j in range(n):
                if state[j] == 0:
                    continue
                factor_j = np.sqrt(state[j])
                for i in range(n):
                    val = A_J[i, j]
                    if abs(val) < 1e-9:
                        continue
                    new_state    = state.copy()
                    new_state[j] -= 1
                    factor_i     = np.sqrt(new_state[i] + 1)
                    new_state[i] += 1
                    key = tuple(new_state)
                    if key in state_index:
                        add(state_index[key], global_col, val * factor_j * factor_i)

            # D_J[i,j]: b†_j b_i |state>  (normal-ordered, delta dropped)
            for i in range(n):
                if state[i] == 0:
                    continue
                factor_i = np.sqrt(state[i])
                for j in range(n):
                    val = D_J[i, j]
                    if abs(val) < 1e-9:
                        continue
                    new_state    = state.copy()
                    new_state[i] -= 1
                    factor_j     = np.sqrt(new_state[j] + 1)
                    new_state[j] += 1
                    key = tuple(new_state)
                    if key in state_index:
                        add(state_index[key], global_col, val * factor_i * factor_j)

        rows = [k[0] for k in data_dict]
        cols = [k[1] for k in data_dict]
        data = list(data_dict.values())

        return sp.csr_matrix((data, (rows, cols)), shape=(n_total, n_total))


    def matter_hamiltonians(self, n_excitations, returns_subspace=False):
        n = 2 * self.N

        if n_excitations < 0:
            raise ValueError("Excitations cannot be negative")

        arr           = np.array([np.array(combo) for combo in compositions(n, n_excitations)])
        states_indexes = {tuple(s): i for i, s in enumerate(arr)}
        n_total        = len(arr)

        if returns_subspace:
            H_C_fock = self.build_HC_fock(arr, n_total)
            H_J_fock = self.build_HJ_fock(arr, states_indexes, n_total)
            return arr, states_indexes, n_total, H_C_fock, H_J_fock
        else:
            self.fock_states      = arr
            self.fock_state_index = states_indexes
            self.fock_dim         = n_total
            self.H_C_fock = self.build_HC_fock(arr, n_total)
            self.H_J_fock = self.build_HJ_fock(arr, states_indexes, n_total)

    def build_full_hamiltonian(self, linear_coupling=True, sw=False):
        # --- photon sector ---
        H_ph = self.omega_r * _num(self.fock_photon)  # ω_r a†a

        # --- Interaction operator ---
        if linear_coupling:
            # Linearized coupling: γ(a + a†)
            sin_op = self.gamma * (_destroy(self.fock_photon) + _destroy(self.fock_photon).conj().T)
        else:
            # Exact non-linear coupling: sin(γ(a + a†))
            sin_op = -0.5j * (
                compute_position_exponential_sparse(self.fock_photon,  self.gamma)
                - compute_position_exponential_sparse(self.fock_photon, -self.gamma)
            )

        # Store sin_op as dense array for reuse in schrieffer_wolff_correction
        self.sin_op = sin_op.toarray() if sp.issparse(sin_op) else np.array(sin_op)

        eye_ph = _eye(self.fock_photon)
        eye_tr = sp.eye(self.fock_dim, format='csr')

        # --- chain sectors ---
        H0 = self.H_C_fock
        H1 = self.H_J_fock

        # --- combine with kron (photon is the leftmost/slowest index) ---
        self.H = (
            sp.kron(H_ph,   eye_tr,  format='csr')   # H_ph  ⊗ I_chain
            + sp.kron(eye_ph, H0,    format='csr')   # I_ph  ⊗ H_C
            - sp.kron(sin_op, H1,    format='csr')   # sin   ⊗ H_J
        )

        if sw:
            self.schrieffer_wolff_correction()
            self.H = self.H + self.H_SW
        
    def schrieffer_wolff_correction(self, eps_denom=1e-6):
        """
        Second-order Schrieffer-Wolff correction from the B/C (±2 excitation) blocks
        of H_J, mediated by sin_op.

        The perturbation is V = -sin_op ⊗ H_J_BC, connecting excitation sector X
        to virtual sectors X±2.  The correction lives entirely within sector X and
        has matrix elements:

          <n1,j1|H_SW|n2,j2> = 0.5 * sum_{r=±1, j_r, n} sin[n1,n]*sin[n,n2]
                                * V_r[j_r,j1] * V_r[j_r,j2].conj()
                                * (1/denom1 - 1/denom2)

        where denom1 = (n1-n)*omega_r + E_X[j1] - E_r[j_r]
              denom2 = (n -n2)*omega_r + E_r[j_r] - E_X[j2]

        Indexing convention: photon is the slow (left) index, matching sp.kron order
        throughout the class.  Full index = n_photon * d_X + j_matter.

        Stores result as self.H_SW (sparse, same shape as self.H).
        """
        n          = 2 * self.N
        d_X        = self.fock_dim
        d_ph       = self.fock_photon
        omega_r    = self.omega_r

        # n_excitations of the target sector X
        n_exc_X    = int(self.fock_states[0].sum())

        # energies_vec: factor of 2 from the Bogoliubov transformation
        energies_vec = 2.0 * np.diag(self.H_C)[:n]   # shape (2N,)

        # Fock-state energies in sector X: E_X[j] = fock_states[j] @ energies_vec
        E_X = self.fock_states @ energies_vec          # (d_X,)

        # B and C blocks of H_J (throwing away the ±2 parts that build_HJ_fock ignores)
        B_J = self.H_J[:n, n:]   # b†_i b†_j  raises excitation count by +2
        C_J = self.H_J[n:, :n]   # b_i  b_j   lowers excitation count by -2

        # sin_op dense matrix (fock_photon, fock_photon) — stored by build_full_hamiltonian
        sin_op = self.sin_op      # shape (d_ph, d_ph), already dense

        # -----------------------------------------------------------------------
        # Build V_r: matter matrix element between virtual sector X+2r and X.
        # V_r[j_r, j] = <j_r| H_J_BC |j>  (with the sign from the -sin⊗H_J coupling)
        # The sign (-1) from the coupling is folded in after: overall V = -sin ⊗ H_J_BC,
        # so the matrix element squared picks up (-1)^2 = 1 — no sign needed in V_r itself.
        # -----------------------------------------------------------------------
        def build_V_r(states_r, index_r, d_r, H_block, states_X, exc_direction):
            """
            Build the (d_r, d_X) matter matrix <j_r | H_block | j>.

            For exc_direction = +1 (B block, X -> X+2):
                H_block = B_J, action is b†_i b†_j |state_X>  -> state in X+2
            For exc_direction = -1 (C block, X -> X-2):
                H_block = C_J, action is b_i b_j |state_X>    -> state in X-2
            """
            rows, cols, data = [], [], []

            if exc_direction == +1:
                # B block: b†_i b†_j |state>  -- raises by 2
                for col_idx, state in enumerate(states_X):
                    for j in range(n):
                        for i in range(n):
                            val = H_block[i, j]
                            if abs(val) < 1e-9:
                                continue
                            new_state    = state.copy()
                            new_state[j] += 1
                            factor_j     = np.sqrt(new_state[j])
                            new_state[i] += 1
                            factor_i     = np.sqrt(new_state[i])
                            key = tuple(new_state)
                            if key in index_r:
                                rows.append(index_r[key])
                                cols.append(col_idx)
                                data.append(val * factor_j * factor_i)

            else:
                # C block: b_i b_j |state>  -- lowers by 2
                for col_idx, state in enumerate(states_X):
                    for j in range(n):
                        if state[j] == 0:
                            continue
                        factor_j  = np.sqrt(state[j])
                        tmp_state = state.copy()
                        tmp_state[j] -= 1
                        for i in range(n):
                            val = H_block[i, j]
                            if abs(val) < 1e-9:
                                continue
                            if tmp_state[i] == 0:
                                continue
                            new_state    = tmp_state.copy()
                            factor_i     = np.sqrt(new_state[i])
                            new_state[i] -= 1
                            key = tuple(new_state)
                            if key in index_r:
                                rows.append(index_r[key])
                                cols.append(col_idx)
                                data.append(val * factor_j * factor_i)

            if not data:
                return np.zeros((d_r, d_X))
            V = np.zeros((d_r, d_X), dtype=complex)
            for r, c, v in zip(rows, cols, data):
                V[r, c] += v
            return V

        # -----------------------------------------------------------------------
        # Main correction matrix in the X subspace
        # -----------------------------------------------------------------------
        H_SW = np.zeros((d_ph * d_X, d_ph * d_X), dtype=complex)

        for r_sign, exc_delta, H_block in [(+1, +2, B_J), (-1, -2, C_J)]:
            n_exc_r = n_exc_X + exc_delta
            if n_exc_r < 0:
                continue

            # Build virtual sector
            states_r_arr  = np.array([np.array(c) for c in compositions(n, n_exc_r)])
            index_r       = {tuple(s): i for i, s in enumerate(states_r_arr)}
            d_r           = len(states_r_arr)

            E_r = states_r_arr @ energies_vec   # (d_r,)

            # V_r: (d_r, d_X)  — matter matrix elements
            V_r = build_V_r(states_r_arr, index_r, d_r, H_block, self.fock_states, r_sign)

            # Triple loop over photon indices (at most fock_photon^3, typically ≤ 64)
            for n1 in range(d_ph):
                for n2 in range(d_ph):
                    for n_mid in range(d_ph):
                        s_factor = sin_op[n1, n_mid] * sin_op[n_mid, n2]
                        if abs(s_factor) < 1e-15:
                            continue

                        # denom1[j1, j_r] = (n1 - n_mid)*omega_r + E_X[j1] - E_r[j_r]
                        denom1 = (n1 - n_mid) * omega_r + E_X[:, None] - E_r[None, :]  # (d_X, d_r)
                        # denom2[j_r, j2] = (n_mid - n2)*omega_r + E_r[j_r] - E_X[j2]
                        denom2 = (n_mid - n2) * omega_r + E_r[:, None] - E_X[None, :]  # (d_r, d_X)

                        # Guard against (near-)degeneracy
                        safe1 = np.abs(denom1) > eps_denom
                        safe2 = np.abs(denom2) > eps_denom

                        inv1 = np.where(safe1, 1.0 / np.where(safe1, denom1, 1.0), 0.0)  # (d_X, d_r)
                        inv2 = np.where(safe2, 1.0 / np.where(safe2, denom2, 1.0), 0.0)  # (d_r, d_X)

                        # Contribution: sum_{j_r} V_r[j_r,j1]*V_r[j_r,j2].conj() * (inv1[j1,j_r] - inv2[j_r,j2])
                        # = (V_r.T * inv1).T @ V_r.conj()  -  V_r.T @ (V_r.conj() * inv2)
                        #   both are (d_X, d_X)
                        term1 = (V_r * inv1.T).T @ V_r.conj()          # (d_X, d_X)
                        term2 = V_r.T @ (V_r.conj() * inv2)          # (d_X, d_X)

                        block = s_factor * (term1 - term2)

                        # Place into full matrix at photon block (n1, n2)
                        H_SW[n1 * d_X:(n1 + 1) * d_X,
                             n2 * d_X:(n2 + 1) * d_X] += block

        # Symmetrize to enforce Hermiticity
        H_SW = 0.5 * (H_SW + H_SW.conj().T)

        self.H_SW = 0.5 * sp.csr_matrix(H_SW)

    def get_low_lying_spectrum(self, k=5):
        """
        Computes the lowest 'k' eigenvalues and eigenvectors of the full Hamiltonian.
        
        Parameters:
        k (int): Number of low-lying states to compute (ground state + k-1 excited states).
        linear_coupling (bool): Whether to use the linearized cavity coupling.
        
        Returns:
        evals (ndarray): The 'k' lowest eigenvalues.
        evecs (ndarray): The corresponding 'k' eigenvectors (as columns).
        """

        H = 0.5 * (self.H + self.H.conj().T)
        
        # 3. Solve for the lowest eigenvalues using the Lanczos method
        # which='SA' means "Smallest Algebraic" (the most negative/lowest energies)
        try:
            evals, evecs = sla.eigsh(H, k=k, which='SA')
        except sla.ArpackNoConvergence as e:
            print(f"Warning: ARPACK did not converge. Returning the {len(e.eigenvalues)} found values.")
            evals, evecs = e.eigenvalues, e.eigenvectors
            
        # ARPACK doesn't always guarantee perfectly sorted output, so we sort them here
        idx = np.argsort(evals)
        return evals[idx], evecs[:, idx]


    def build_squeezed_vacuum(self, truncation, threshold=0.99):
        n       = 2 * self.N
        U_tilde = self.inverse_P[:n, :n]
        V_tilde = self.inverse_P[:n, n:]
        M       = np.linalg.solve(U_tilde, V_tilde)

        vacuum                   = tuple(np.zeros(n, dtype=int))
        state_coeff_dict         = {vacuum: 1.0}
        current_terms            = {vacuum: 1.0}

        def add(d, key, val):
            if key in d: d[key] += val
            else:        d[key]  = val

        for k in range(1, truncation + 1):
            next_terms = {}
            for state_tuple, coeff in current_terms.items():
                state = np.array(state_tuple, dtype=int)
                for j in range(n):
                    for i in range(n):
                        val = M[i, j]
                        if abs(val) < 1e-9:
                            continue
                        new_state    = state.copy()
                        new_state[j] += 1
                        factor_j     = np.sqrt(new_state[j])
                        new_state[i] += 1
                        factor_i     = np.sqrt(new_state[i])
                        new_key      = tuple(new_state)
                        add(next_terms, new_key, coeff * (-0.5) * val * factor_j * factor_i / float(k))
            for key, val in next_terms.items():
                if abs(val) > 1e-9:
                    add(state_coeff_dict, key, val)
            current_terms = next_terms

        # --- check weight before normalizing ---
        raw_norm_sq = sum(abs(v)**2 for v in state_coeff_dict.values())
        if raw_norm_sq < threshold:
            import warnings
            warnings.warn(
                f"build_squeezed_vacuum: truncation={truncation} captured only "
                f"{raw_norm_sq:.4f} of the total weight (threshold={threshold}). "
                f"Consider increasing truncation.", RuntimeWarning
            )

        norm = np.sqrt(raw_norm_sq)
        state_coeff_dict = {k: v / norm for k, v in state_coeff_dict.items()}

        # --- threshold filter by descending weight ---
        sorted_items = sorted(state_coeff_dict.items(), key=lambda x: abs(x[1])**2, reverse=True)
        filtered   = {}
        cumulative = 0.0
        for key, val in sorted_items:
            filtered[key] = val
            cumulative   += abs(val)**2
            if cumulative >= threshold:
                break

        norm_filtered = np.sqrt(sum(abs(v)**2 for v in filtered.values()))
        self.vacuum          = {k: v / norm_filtered for k, v in filtered.items()}
        self.vacuum_n_states = len(filtered)


    def build_unsqueezed_vacuum(self, truncation, threshold=0.99):
        n       = 2 * self.N
        U_tilde = self.P[:n, :n]
        V_tilde = self.P[:n, n:]
        M       = np.linalg.solve(U_tilde, V_tilde)

        vacuum                   = tuple(np.zeros(n, dtype=int))
        state_coeff_dict         = {vacuum: 1.0}
        current_terms            = {vacuum: 1.0}

        def add(d, key, val):
            if key in d: d[key] += val
            else:        d[key]  = val

        for k in range(1, truncation + 1):
            next_terms = {}
            for state_tuple, coeff in current_terms.items():
                state = np.array(state_tuple, dtype=int)
                for j in range(n):
                    for i in range(n):
                        val = M[i, j]
                        if abs(val) < 1e-9:
                            continue
                        new_state    = state.copy()
                        new_state[j] += 1
                        factor_j     = np.sqrt(new_state[j])
                        new_state[i] += 1
                        factor_i     = np.sqrt(new_state[i])
                        new_key      = tuple(new_state)
                        add(next_terms, new_key, coeff * (-0.5) * val * factor_j * factor_i / float(k))
            for key, val in next_terms.items():
                if abs(val) > 1e-9:
                    add(state_coeff_dict, key, val)
            current_terms = next_terms

        # --- check weight before normalizing ---
        raw_norm_sq = sum(abs(v)**2 for v in state_coeff_dict.values())
        if raw_norm_sq < threshold:
            import warnings
            warnings.warn(
                f"build_unsqueezed_vacuum: truncation={truncation} captured only "
                f"{raw_norm_sq:.4f} of the total weight (threshold={threshold}). "
                f"Consider increasing truncation.", RuntimeWarning
            )

        norm = np.sqrt(raw_norm_sq)
        state_coeff_dict = {k: v / norm for k, v in state_coeff_dict.items()}

        sorted_items = sorted(state_coeff_dict.items(), key=lambda x: abs(x[1])**2, reverse=True)
        filtered   = {}
        cumulative = 0.0
        for key, val in sorted_items:
            filtered[key] = val
            cumulative   += abs(val)**2
            if cumulative >= threshold:
                break

        norm_filtered = np.sqrt(sum(abs(v)**2 for v in filtered.values()))
        self.unsqueezed_vacuum          = {k: v / norm_filtered for k, v in filtered.items()}
        self.unsqueezed_vacuum_n_states = len(filtered)
        
    def apply_eta_dag(self, i, state_coeff_dict, threshold=1e-9):
        n   = 2 * self.N
        U_i = self.P[:n, i]       # U[:, i], shape (n,)
        V_i = -self.P[n:, i]      # -V[:, i]

        # precompute nonzero indices once
        u_nz = np.where(np.abs(U_i) > threshold)[0]
        v_nz = np.where(np.abs(V_i) > threshold)[0]

        result = {}
        for state_tuple, coeff in state_coeff_dict.items():
            state = np.array(state_tuple, dtype=int)

            for j in u_nz:
                new_state    = state.copy()
                new_state[j] += 1
                factor        = np.sqrt(new_state[j])
                key           = tuple(new_state)
                val           = coeff * U_i[j] * factor
                result[key]   = result.get(key, 0.0) + val

            for j in v_nz:
                if state[j] == 0:
                    continue
                new_state    = state.copy()
                factor        = np.sqrt(state[j])
                new_state[j] -= 1
                key           = tuple(new_state)
                val           = coeff * V_i[j] * factor
                result[key]   = result.get(key, 0.0) + val

        return {k: v for k, v in result.items() if abs(v) > threshold}

    def apply_b_dag(self, i, state_coeff_dict, threshold=1e-9):
        n   = 2 * self.N
        U_i = self.P[i, :n]
        V_i = self.P[i, n:]     

        # precompute nonzero indices once
        u_nz = np.where(np.abs(U_i) > threshold)[0]
        v_nz = np.where(np.abs(V_i) > threshold)[0]

        result = {}
        for state_tuple, coeff in state_coeff_dict.items():
            state = np.array(state_tuple, dtype=int)

            for j in u_nz:
                new_state    = state.copy()
                new_state[j] += 1
                factor        = np.sqrt(new_state[j])
                key           = tuple(new_state)
                val           = coeff * U_i[j] * factor
                result[key]   = result.get(key, 0.0) + val

            for j in v_nz:
                if state[j] == 0:
                    continue
                new_state    = state.copy()
                factor        = np.sqrt(state[j])
                new_state[j] -= 1
                key           = tuple(new_state)
                val           = coeff * V_i[j] * factor
                result[key]   = result.get(key, 0.0) + val

        return {k: v for k, v in result.items() if abs(v) > threshold}
    
    def new_to_old_fock(self, new_fock_state, threshold=0.99):
        n           = 2 * self.N
        occupations = np.array(new_fock_state, dtype=int)

        state_coeff_dict = dict(self.vacuum)

        for i in range(n):
            r_i = occupations[i]
            if r_i == 0:
                continue
            for k in range(1, r_i + 1):
                state_coeff_dict = self.apply_eta_dag(i, state_coeff_dict)
                state_coeff_dict = {s: v / np.sqrt(float(k)) for s, v in state_coeff_dict.items()}

        norm = np.sqrt(sum(abs(v)**2 for v in state_coeff_dict.values()))
        if norm > 1e-14:
            state_coeff_dict = {k: v / norm for k, v in state_coeff_dict.items()}

        # --- threshold filter ---
        sorted_items = sorted(state_coeff_dict.items(), key=lambda x: abs(x[1])**2, reverse=True)
        filtered   = {}
        cumulative = 0.0
        for key, val in sorted_items:
            filtered[key] = val
            cumulative   += abs(val)**2
            if cumulative >= threshold:
                break

        norm_filtered = np.sqrt(sum(abs(v)**2 for v in filtered.values()))
        if norm_filtered > 1e-14:
            filtered = {k: v / norm_filtered for k, v in filtered.items()}

        return filtered


    def old_to_new_fock(self, old_fock_state, threshold=0.99):
        n           = 2 * self.N
        occupations = np.array(old_fock_state, dtype=int)

        state_coeff_dict = dict(self.unsqueezed_vacuum)

        for i in range(n):
            r_i = occupations[i]
            if r_i == 0:
                continue
            for k in range(1, r_i + 1):
                state_coeff_dict = self.apply_b_dag(i, state_coeff_dict)
                state_coeff_dict = {s: v / np.sqrt(float(k)) for s, v in state_coeff_dict.items()}

        norm = np.sqrt(sum(abs(v)**2 for v in state_coeff_dict.values()))
        if norm > 1e-14:
            state_coeff_dict = {k: v / norm for k, v in state_coeff_dict.items()}

        # --- threshold filter ---
        sorted_items = sorted(state_coeff_dict.items(), key=lambda x: abs(x[1])**2, reverse=True)
        filtered   = {}
        cumulative = 0.0
        for key, val in sorted_items:
            filtered[key] = val
            cumulative   += abs(val)**2
            if cumulative >= threshold:
                break

        norm_filtered = np.sqrt(sum(abs(v)**2 for v in filtered.values()))
        if norm_filtered > 1e-14:
            filtered = {k: v / norm_filtered for k, v in filtered.items()}

        return filtered

    
    def quanta_distribution(self, fock_vector, state_coeff_dict=None):
        """
        Computes the quanta probability distribution over old bosonic modes.
        
        Given a state |Ψ> = Σ_s c_s |s>, computes:
            n_i = Σ_s |c_s|² * s_i   (expected occupation per mode)
        then normalizes so that Σ_i n_i = 1.
        
        Parameters
        ----------
        fock_vector      : dict {tuple(fock_state) -> coefficient}
                        OR 1D array in the global Fock basis ordering
        
        Returns
        -------
        distribution : 1D array of length n, normalized quanta distribution
        """
        n = 2 * self.N

        if isinstance(fock_vector, dict):
            state_coeff_dict = fock_vector
            distribution = np.zeros(n)
            for state_tuple, coeff in state_coeff_dict.items():
                distribution += (np.abs(coeff)**2) * np.array(state_tuple)
        else:
            # fock_states is a flat ndarray (fock_dim, 2N); vectorized path
            weights      = np.abs(fock_vector)**2          # (fock_dim,)
            distribution = weights @ self.fock_states      # (2N,)

        total = distribution.sum()
        if total > 1e-15:
            distribution /= total

        return distribution


    def quanta_distribution_subspace(self, fock_vector, n_excitations, cumulative=False):
        """
        With fock_layers=0 there is exactly one excitation sector, so this is
        identical to quanta_distribution.  The n_excitations and cumulative
        arguments are kept for API compatibility but are ignored.
        """
        return self.quanta_distribution(fock_vector)


    def quanta_ipr(self, fock_vector):
        """
        Computes the quanta-weighted Inverse Participation Ratio (IPR):
        
            IPR = Σ_i <n_i>² / (Σ_i <n_i>)²
        
        where <n_i> = Σ_s |c_s|² * s_i is the expected occupation of mode i.
        
        IPR -> 1/n_modes : fully delocalized
        IPR -> 1         : fully localized on one mode
        
        Parameters
        ----------
        fock_vector : 1D array in the global Fock basis ordering
        
        Returns
        -------
        ipr          : float, the IPR value
        xi           : float, effective localization length (= 1/IPR in units of modes)
        distribution : 1D array, the unnormalized <n_i> used for the IPR
        """
        weights      = np.abs(fock_vector)**2          # (fock_dim,)
        distribution = weights @ self.fock_states      # (2N,)

        total = distribution.sum()
        if total < 1e-15:
            return 0.0, np.inf, distribution

        ipr = np.sum(distribution**2) / (total**2)
        xi  = 1.0 / ipr

        return ipr, xi, distribution
    
    def quanta_ipr_mixed(self, rho_matter):
        """
        Generalization of quanta_ipr to a mixed matter state.
        
        Computes IPR = Σ_i <n_i>² / (Σ_i <n_i>)²
        where <n_i> = Tr[rho_matter * n_i] = Σ_s rho_ss * s_i
        (diagonal of rho_matter replaces |c_s|² from the pure case).

        Parameters
        ----------
        rho_matter : 2D array (fock_dim, fock_dim), reduced density matrix
                    from partial_trace(..., which_to_keep='matter')

        Returns
        -------
        ipr          : float
        xi           : float, effective number of modes (= 1/IPR)
        distribution : 1D array of length 2*N, the <n_i> per mode
        """
        diag_rho     = np.real(np.diag(rho_matter))   # (fock_dim,)
        distribution = diag_rho @ self.fock_states     # (2N,)

        total = distribution.sum()
        if total < 1e-15:
            return 0.0, np.inf, distribution

        ipr = np.sum(distribution**2) / total**2
        xi  = 1.0 / ipr

        return ipr, xi, distribution
    
    def decompose_tensor_state(self, full_vector, threshold=0.99):
        d_ph  = self.fock_photon
        d_mat = self.fock_dim

        assert len(full_vector) == d_ph * d_mat, (
            f"Vector length {len(full_vector)} != fock_photon*fock_dim = {d_ph*d_mat}"
        )

        psi = full_vector.reshape(d_ph, d_mat)

        decomposition = {}
        for alpha in range(d_ph):
            matter_slice = psi[alpha, :]
            matter_dict  = {}
            for s_idx, state in enumerate(self.fock_states):
                coeff = matter_slice[s_idx]
                if abs(coeff) > 1e-15:
                    matter_dict[tuple(state)] = matter_dict.get(tuple(state), 0.0) + coeff

            if not matter_dict:
                continue

            # --- threshold filter per photon sector ---
            sorted_items = sorted(matter_dict.items(), key=lambda x: abs(x[1])**2, reverse=True)
            filtered   = {}
            cumulative = 0.0
            for key, val in sorted_items:
                filtered[key] = val
                cumulative   += abs(val)**2
                if cumulative >= threshold:
                    break

            # keep original coefficients (no renorm) to preserve global norm
            norm_filtered = np.sqrt(sum(abs(v)**2 for v in filtered.values()))
            if norm_filtered > 1e-15:
                decomposition[alpha] = filtered

        return decomposition
            

    def compute_ipr_batch(self, vecs):
        """
        Vectorized IPR computation over a batch of eigenvectors in the old Fock basis.

        Requires self._new_to_old_cache and self._old_basis_index to be pre-built
        via warm_old_basis_cache().

        Parameters
        ----------
        vecs : ndarray of shape (fock_photon * fock_dim, k)
               Columns are eigenvectors in the Bogoliubov-photon product basis.

        Returns
        -------
        iprs : ndarray of shape (k,)
        """
        if not hasattr(self, '_C_matrix'):
            raise RuntimeError("Call warm_old_basis_cache() before compute_ipr_batch().")

        d_ph  = self.fock_photon
        d_mat = self.fock_dim
        k     = vecs.shape[1]
        n     = 2 * self.N

        C                 = self._C_matrix          # (d_mat_old, d_mat)
        sorted_states_arr = self._old_states_arr    # (d_mat_old, n)

        # project all eigenvectors: (d_ph, d_mat, k) -> (d_ph, d_mat_old, k)
        vecs_3d = vecs.reshape(d_ph, d_mat, k)
        amp_3d  = np.tensordot(C, vecs_3d, axes=([1], [1]))   # (d_mat_old, d_ph, k)
        amp_3d  = amp_3d.transpose(1, 0, 2)                   # (d_ph, d_mat_old, k)

        # diagonal of reduced matter density matrix
        diag_rho = np.sum(np.abs(amp_3d)**2, axis=0)          # (d_mat_old, k)

        # quanta distribution and IPR
        all_dist = sorted_states_arr.T @ diag_rho             # (n, k)
        totals   = all_dist.sum(axis=0)                        # (k,)
        iprs     = np.where(
            totals < 1e-15,
            0.0,
            np.sum(all_dist**2, axis=0) / totals**2
        )
        return iprs

    def warm_old_basis_cache(self):
        """
        Pre-populate _new_to_old_cache for every Bogoliubov Fock state in
        self.fock_states, then build the dense change-of-basis matrix _C_matrix
        and old-state array _old_states_arr needed by compute_ipr_batch().

        Call this once after build_squeezed_vacuum() and before the IPR sweep.
        """
        if not hasattr(self, '_new_to_old_cache'):
            self._new_to_old_cache = {}

        # populate cache for every new Fock state
        for state in self.fock_states:
            key = tuple(state)
            if key not in self._new_to_old_cache:
                self._new_to_old_cache[key] = self.new_to_old_fock(
                    np.array(state, dtype=int)
                )

        # collect all old-basis states
        all_old_states = set()
        for old_dict in self._new_to_old_cache.values():
            all_old_states.update(old_dict.keys())

        sorted_states     = sorted(all_old_states, key=lambda s: (sum(s), s))
        state_to_idx      = {s: i for i, s in enumerate(sorted_states)}
        d_mat_old         = len(sorted_states)
        d_mat             = self.fock_dim

        # build C: (d_mat_old, d_mat)
        C = np.zeros((d_mat_old, d_mat), dtype=complex)
        for s_idx, state in enumerate(self.fock_states):
            new_key  = tuple(state)
            old_dict = self._new_to_old_cache[new_key]
            norm     = np.sqrt(sum(abs(v)**2 for v in old_dict.values()))
            if norm < 1e-15:
                continue
            for old_state, coeff in old_dict.items():
                C[state_to_idx[old_state], s_idx] += coeff / norm

        self._C_matrix       = C
        self._old_states_arr = np.array(sorted_states, dtype=float)   # (d_mat_old, 2N)
        self._old_state_idx  = state_to_idx

    def decompose_to_old_fock(self, full_vector, threshold=0.99):
        bogoliubov_decomp = self.decompose_tensor_state(full_vector, threshold=threshold)

        # build cache on first call, keyed by new Fock state tuple
        if not hasattr(self, '_new_to_old_cache'):
            self._new_to_old_cache = {}

        old_decomposition = {}
        for alpha, matter_dict in bogoliubov_decomp.items():
            old_matter = {}
            for new_state_tuple, coeff in matter_dict.items():
                if new_state_tuple not in self._new_to_old_cache:
                    self._new_to_old_cache[new_state_tuple] = self.new_to_old_fock(
                        np.array(new_state_tuple, dtype=int), threshold=threshold
                    )
                old_part = self._new_to_old_cache[new_state_tuple]
                norm_old = np.sqrt(sum(abs(v)**2 for v in old_part.values()))
                for old_state, old_coeff in old_part.items():
                    weighted = coeff * old_coeff / (norm_old if norm_old > 1e-15 else 1.0)
                    old_matter[old_state] = old_matter.get(old_state, 0.0) + weighted
            if old_matter:
                # threshold filter as before
                sorted_items = sorted(old_matter.items(), key=lambda x: abs(x[1])**2, reverse=True)
                filtered, cumulative, total_w = {}, 0.0, sum(abs(v)**2 for v in old_matter.values())
                for key, val in sorted_items:
                    filtered[key] = val
                    cumulative += abs(val)**2
                    if cumulative >= threshold * total_w:
                        break
                old_decomposition[alpha] = filtered

        return old_decomposition

    def compose_tensor_state(self, decomposition):
        """
        Inverse of decompose_tensor_state, generalized to handle matter Fock spaces
        that may differ from self.fock_dim (e.g. output of decompose_to_old_fock).

        The matter Hilbert space dimension and basis are inferred directly from
        the union of all state_tuples appearing across all photon sectors.

        Parameters
        ----------
        decomposition : dict { alpha (int) -> dict { tuple(matter_state) -> complex } }

        Returns
        -------
        full_vector : 1D array of length fock_photon * d_mat_inferred
        """
        d_ph = self.fock_photon

        # --- collect all unique matter states across all photon sectors ---
        all_states = set()
        for matter_dict in decomposition.values():
            all_states.update(matter_dict.keys())

        # --- build a canonical ordering for the inferred basis ---
        # sort by excitation number first, then lexicographically within each sector
        sorted_states = sorted(all_states, key=lambda s: (sum(s), s))
        state_to_idx  = {s: i for i, s in enumerate(sorted_states)}
        d_mat_inferred = len(sorted_states)

        # --- fill amplitude matrix ---
        psi = np.zeros((d_ph, d_mat_inferred), dtype=complex)

        for alpha, matter_dict in decomposition.items():
            if alpha >= d_ph:
                continue
            for state_tuple, coeff in matter_dict.items():
                m = state_to_idx[state_tuple]
                psi[alpha, m] += coeff

        return psi.reshape(d_ph * d_mat_inferred), sorted_states, state_to_idx

    def partial_trace(self, full_vector, which_to_keep='matter'):
        """
        Computes the reduced density matrix by tracing out one subsystem.

        The full Hilbert space is H_photon ⊗ H_matter, with global index
            idx = alpha * fock_dim + m
        (photon is the slow/left index from sp.kron).

        Parameters
        ----------
        full_vector   : 1D array of length fock_photon * fock_dim
        which_to_keep : 'matter' or 'photon'

        Returns
        -------
        rho : 2D dense array of shape (fock_dim, fock_dim)   if which_to_keep == 'matter'
                                    (fock_photon, fock_photon) if which_to_keep == 'photon'
        """
        d_ph  = self.fock_photon
        d_mat = self.fock_dim

        psi = full_vector.reshape(d_ph, d_mat)   # (alpha, m)

        if which_to_keep == 'matter':
            # rho_matter_{m,m'} = Σ_alpha psi[alpha,m]* psi[alpha,m'].conj()
            return psi.conj().T @ psi            # (fock_dim, fock_dim)
        elif which_to_keep == 'photon':
            # rho_photon_{alpha,alpha'} = Σ_m psi[alpha,m] * psi[alpha',m].conj()
            return psi @ psi.conj().T            # (fock_photon, fock_photon)
        else:
            raise ValueError(f"which_to_keep must be 'matter' or 'photon', got '{which_to_keep}'")