import numpy as np
from dataclasses import dataclass
from typing import Tuple, List, Dict
from collections import defaultdict
import itertools
import scipy.sparse as sp
import math
import warnings
import scipy.linalg as la
import scipy.sparse.linalg as sla

try:
    import cupy as cp
    _CUPY_AVAILABLE = True
except ImportError:
    _CUPY_AVAILABLE = False

@dataclass
class TensorTerm:
    tensor: np.ndarray
    daggers: Tuple[int, ...]

class Superposition:
    """An object-oriented wrapper for a linear combination of Fock states."""
    def __init__(self, state_dict=None):
        self.states = defaultdict(complex, state_dict if state_dict else {})

    def __add__(self, other, tol=1e-12):
        result = self.states.copy()
        for st, amp in other.states.items():
            result[st] += amp
        return Superposition(self._clean_dict(result, tol=tol))

    def __rmul__(self, scalar, tol=1e-12):
        result = {st: scalar * amp for st, amp in self.states.items()}
        return Superposition(self._clean_dict(result, tol=tol))

    def _clean_dict(self, d, tol=1e-12):
        return {st: amp for st, amp in d.items() if abs(amp) > tol}

    def normalize(self):
        norm = np.sqrt(sum(abs(amp)**2 for amp in self.states.values()))
        if norm > 0:
            self.states = {st: amp / norm for st, amp in self.states.items()}
        return self

    def trim(self, threshold=0.99):
        self.normalize()
        sorted_items = sorted(self.states.items(), key=lambda item: abs(item[1])**2, reverse=True)
        trimmed_state = {}
        accumulated_prob = 0.0
        
        for st, amp in sorted_items:
            trimmed_state[st] = amp
            accumulated_prob += abs(amp)**2
            if accumulated_prob >= threshold:
                break
                
        self.states = trimmed_state
        return self.normalize()

class QuantumOperator:
    """An algebraically aware container for tensor terms."""
    def __init__(self, terms=None):
        self.terms = terms if terms is not None else []

    def __add__(self, other):
        return QuantumOperator(self.terms + other.terms)

    def __rmul__(self, scalar):
        new_terms = [TensorTerm(term.tensor * scalar, term.daggers) for term in self.terms]
        return QuantumOperator(new_terms)

    def __mul__(self, other):
        """Computes the product of two operators via tensor outer product."""
        new_terms = []
        for tA in self.terms:
            for tB in other.terms:
                # np.tensordot with axes=0 computes the full outer product
                new_tensor = np.tensordot(tA.tensor, tB.tensor, axes=0)
                new_daggers = tA.daggers + tB.daggers
                new_terms.append(TensorTerm(new_tensor, new_daggers))
        return QuantumOperator(new_terms)

    def simplify(self, tol: float = 1e-12):
        """Optimizes the operator by merging tensors with identical dagger structures."""
        accumulated = {}
        for term in self.terms:
            if term.daggers not in accumulated:
                accumulated[term.daggers] = np.zeros_like(term.tensor)
            accumulated[term.daggers] += term.tensor
        self.terms = []
        for dags, tensor in accumulated.items():
            tensor = tensor.copy()
            tensor[np.abs(tensor) < tol] = 0.0   # elementwise cleanup
            if np.any(np.abs(tensor) > 0):
                self.terms.append(TensorTerm(tensor, dags))
        return self

class SecondQuantizationEngine:
    def __init__(self, n_modes: int, statistics: str = 'boson',
                 coeff_tol: float = 1e-9, zero_tol: float = 1e-12,
                 use_gpu: bool = False):
        if statistics not in ['boson', 'fermion']:
            raise ValueError("Statistics must be 'boson' or 'fermion'")
        self.n_modes = n_modes
        self.statistics = statistics
        self.sign = 1 if statistics == 'boson' else -1
        self.coeff_tol = coeff_tol
        self.zero_tol = zero_tol
        if use_gpu and not _CUPY_AVAILABLE:
            warnings.warn(
                "use_gpu=True requested but cupy is not installed. "
                "Falling back to CPU. Install cupy matching your CUDA version to enable GPU.",
                RuntimeWarning
            )
            self.use_gpu = False
        else:
            self.use_gpu = use_gpu

    def _get_parity_phase(self, mode: int, state: Tuple[int, ...]) -> int:
        """Calculates the Jordan-Wigner phase for fermions."""
        if self.statistics == 'boson':
            return 1
        # For fermions, count particles in modes prior to the target mode
        particles_before = sum(state[:mode])
        return -1 if particles_before % 2 == 1 else 1

    def _apply_annihilation(self, mode: int, state: Tuple[int, ...]) -> Tuple[Tuple[int, ...], complex]:
        n = state[mode]
        if n == 0:
            return None, 0.0 # Annihilating the vacuum yields 0
        
        phase = self._get_parity_phase(mode, state)
        new_state = list(state)
        new_state[mode] -= 1
        
        amplitude = phase * np.sqrt(n)
        return tuple(new_state), amplitude

    def _apply_creation(self, mode: int, state: Tuple[int, ...]) -> Tuple[Tuple[int, ...], complex]:
        n = state[mode]
        if self.statistics == 'fermion' and n == 1:
            return None, 0.0 # Pauli exclusion principle
            
        phase = self._get_parity_phase(mode, state)
        new_state = list(state)
        new_state[mode] += 1
        
        amplitude = phase * np.sqrt(n + 1)
        return tuple(new_state), amplitude

    def apply_operator(self, operator: QuantumOperator, superposition: Superposition) -> Superposition:
        """Updated to interface seamlessly with the new Superposition class."""
        result_state = defaultdict(complex)
        active_states = superposition._clean_dict(superposition.states, tol=self.coeff_tol)
        
        if not active_states:
            return Superposition()
            
        for term in operator.terms:
            non_zero_indices = np.argwhere(np.abs(term.tensor) > self.coeff_tol)
            for idx in non_zero_indices:
                idx_tuple = tuple(idx)
                coeff = term.tensor[idx_tuple]
                
                current_branches = {st: amp * coeff for st, amp in active_states.items()}
                
                for i in reversed(range(len(term.daggers))):
                    mode = idx_tuple[i]
                    is_dagger = term.daggers[i]
                    next_branches = defaultdict(complex)
                    
                    for curr_state, curr_amp in current_branches.items():
                        if is_dagger:
                            new_st, amp_factor = self._apply_creation(mode, curr_state)
                        else:
                            new_st, amp_factor = self._apply_annihilation(mode, curr_state)
                            
                        if new_st is not None:
                            next_branches[new_st] += curr_amp * amp_factor
                            
                    current_branches = next_branches
                    if not current_branches:
                        break
                        
                for final_st, final_amp in current_branches.items():
                    result_state[final_st] += final_amp
                    
        return Superposition(result_state)
    
    def normalize_state(self, superposition: Dict[Tuple[int, ...], complex]) -> Dict[Tuple[int, ...], complex]:
        """
        Normalizes a linear combination of Fock states so that the sum of |c|^2 is 1.
        """
        # Calculate the norm: sqrt(sum(|c|^2))
        norm_sq = sum(abs(amp)**2 for amp in superposition.values())
        
        if np.isclose(norm_sq, 0.0):
            return {} # Cannot normalize the vacuum/zero-vector
            
        norm = np.sqrt(norm_sq)
        return {st: amp / norm for st, amp in superposition.items()}

    def trim_state(self, superposition: Dict[Tuple[int, ...], complex], threshold: float = 0.99) -> Dict[Tuple[int, ...], complex]:
        """
        Normalizes, sorts by magnitude descending, and truncates to keep the most 
        significant terms that sum up to the probability threshold.
        """
        # Step 1: Normalize the incoming state
        normalized_state = self.normalize_state(superposition)
        
        # Step 2: Sort descending by absolute magnitude squared (|c|^2)
        sorted_items = sorted(normalized_state.items(), key=lambda item: abs(item[1])**2, reverse=True)
        
        # Step 3: Accumulate probabilities and trim
        trimmed_state = {}
        accumulated_prob = 0.0
        
        for st, amp in sorted_items:
            prob = abs(amp)**2
            trimmed_state[st] = amp
            accumulated_prob += prob
            
            # Stop once the accumulated probability hits or exceeds the threshold
            if accumulated_prob >= threshold:
                break
                
        # Step 4: Re-normalize the trimmed state so it is completely unitary again
        return self.normalize_state(trimmed_state)
    
    def normal_order(self, operator: QuantumOperator) -> QuantumOperator:
        """
        Reduces a QuantumOperator to normal order using Wick's theorem.
        Returns a newly optimized QuantumOperator object.
        """
        pending_terms = operator.terms.copy()
        accumulated_terms = {}

        while pending_terms:
            term = pending_terms.pop(0)
            daggers = term.daggers
            tensor = np.asarray(term.tensor)
            
            swap_idx = -1
            for i in range(len(daggers) - 1):
                if daggers[i] == 0 and daggers[i+1] == 1:
                    swap_idx = i
                    break
                    
            if swap_idx != -1:
                # Swapped Term
                new_daggers_swap = list(daggers)
                new_daggers_swap[swap_idx], new_daggers_swap[swap_idx+1] = 1, 0
                new_tensor_swap = np.swapaxes(tensor, swap_idx, swap_idx+1) * self.sign
                pending_terms.append(TensorTerm(new_tensor_swap, tuple(new_daggers_swap)))
                
                # Contracted Term
                new_daggers_contract = tuple(daggers[:swap_idx] + daggers[swap_idx+2:])
                new_tensor_contract = np.trace(tensor, axis1=swap_idx, axis2=swap_idx+1)
                pending_terms.append(TensorTerm(new_tensor_contract, new_daggers_contract))
            else:
                if daggers not in accumulated_terms:
                    accumulated_terms[daggers] = np.zeros_like(tensor)
                accumulated_terms[daggers] += tensor
                
        final_terms = [
            TensorTerm(tensor, dags) 
            for dags, tensor in accumulated_terms.items() 
            if not np.allclose(tensor, 0, atol=self.zero_tol)
        ]
        return QuantumOperator(final_terms)
    
    def change_basis(self, operator: QuantumOperator, U: np.ndarray, V: np.ndarray) -> QuantumOperator:
        """
        Applies a Bogoliubov transformation to the operator.
        b_i = U_ij b_j + V*_ij b^dag_j
        """
        new_terms = []
        
        for original_term in operator.terms:
            k = len(original_term.daggers)
            
            # Start the branching with the current tensor and daggers
            # We use lists for daggers temporarily so we can mutate them during branching
            branches = [(original_term.tensor, list(original_term.daggers))]
            
            for axis in range(k):
                next_branches = []
                
                for tensor, daggers in branches:
                    is_dagger = daggers[axis]
                    
                    if is_dagger == 0:
                        # Annihilation transforms: b_i = U_ij b_j + V*_ij b^dag_j
                        
                        # Branch 1: U (annihilation)
                        t_U = np.tensordot(tensor, U, axes=([axis], [0]))
                        t_U = np.moveaxis(t_U, -1, axis)
                        dag_U = daggers.copy()
                        dag_U[axis] = 0
                        next_branches.append((t_U, dag_U))
                        
                        # Branch 2: V^* (creation)
                        t_V_conj = np.tensordot(tensor, np.conj(V), axes=([axis], [0]))
                        t_V_conj = np.moveaxis(t_V_conj, -1, axis)
                        dag_V = daggers.copy()
                        dag_V[axis] = 1
                        next_branches.append((t_V_conj, dag_V))
                        
                    else:
                        # Creation transforms: b^dag_i = U*_ij b^dag_j + V_ij b_j
                        
                        # Branch 1: U^* (creation)
                        t_U_conj = np.tensordot(tensor, np.conj(U), axes=([axis], [0]))
                        t_U_conj = np.moveaxis(t_U_conj, -1, axis)
                        dag_U_conj = daggers.copy()
                        dag_U_conj[axis] = 1
                        next_branches.append((t_U_conj, dag_U_conj))
                        
                        # Branch 2: V (annihilation)
                        t_V = np.tensordot(tensor, V, axes=([axis], [0]))
                        t_V = np.moveaxis(t_V, -1, axis)
                        dag_V = daggers.copy()
                        dag_V[axis] = 0
                        next_branches.append((t_V, dag_V))
                        
                # Update branches for the next axis
                branches = next_branches
                
            # After iterating all axes, commit the 2^k branches to the new operator
            for tensor, dags in branches:
                if not np.allclose(tensor, 0, atol=self.zero_tol): # Filter out exact zeros to save memory
                    new_terms.append(TensorTerm(tensor, tuple(dags)))
                    
        return QuantumOperator(new_terms)

    def generate_basis(self, n_excitations: int) -> List[Tuple[int, ...]]:
        """
        Generates all valid Fock states for the system with exactly `n_excitations`.
        Automatically handles Bosonic vs Fermionic statistics.
        """
        if n_excitations < 0:
            return []
            
        if self.statistics == 'fermion':
            if n_excitations > self.n_modes:
                return []  # Pauli exclusion: cannot have more fermions than modes
                
            # For fermions, we simply choose `n_excitations` modes to occupy
            basis = []
            for occupied_indices in itertools.combinations(range(self.n_modes), n_excitations):
                state = [0] * self.n_modes
                for idx in occupied_indices:
                    state[idx] = 1
                basis.append(tuple(state))
            return basis
            
        else:
            # For bosons, we use the recursive stars-and-bars integer partitioning
            def compositions(n, total):
                if n == 1:
                    yield (total,)
                    return
                for first in range(total + 1):
                    for rest in compositions(n - 1, total - first):
                        yield (first,) + rest
                        
            return list(compositions(self.n_modes, n_excitations))

    def _precompute_term_nonzeros(self, operator: QuantumOperator):
        """
        Scans each term's tensor for entries above coeff_tol ONCE, returning
        a list of (daggers, [(idx_tuple, coeff), ...]) per term. Shared by
        build_operator_matrix and build_sparse_operator_matrix so that
        building a matrix over a basis of size dim does one tensor scan per
        term instead of dim scans of the same fixed tensor -- the dominant
        cost for any operator whose tensor is dense (e.g. a higher-rank term
        after a basis transformation that destroys its original sparsity).
        """
        term_nonzeros = []
        for term in operator.terms:
            non_zero_indices = np.argwhere(np.abs(term.tensor) > self.coeff_tol)
            entries = [(tuple(idx), term.tensor[tuple(idx)]) for idx in non_zero_indices]
            term_nonzeros.append((term.daggers, entries))
        return term_nonzeros

    def _apply_precomputed_to_state(self, term_nonzeros, ket_state):
        """
        Applies a precomputed (from _precompute_term_nonzeros) operator to a
        single basis ket, returning {result_state: amplitude}. Factored out
        so build_operator_matrix and build_sparse_operator_matrix share the
        identical per-column logic, differing only in how they write the
        result into a dense vs. sparse matrix.
        """
        current_total = defaultdict(complex)
        for daggers, entries in term_nonzeros:
            for idx_tuple, coeff in entries:
                current_branches = {ket_state: coeff}
                for i in reversed(range(len(daggers))):
                    mode = idx_tuple[i]
                    is_dagger = daggers[i]
                    next_branches = defaultdict(complex)
                    for curr_state, curr_amp in current_branches.items():
                        if is_dagger:
                            new_st, amp_factor = self._apply_creation(mode, curr_state)
                        else:
                            new_st, amp_factor = self._apply_annihilation(mode, curr_state)
                        if new_st is not None:
                            next_branches[new_st] += curr_amp * amp_factor
                    current_branches = next_branches
                    if not current_branches:
                        break
                for final_st, final_amp in current_branches.items():
                    current_total[final_st] += final_amp
        return current_total

    def _build_operator_matrix_vectorized_boson(self, operator: QuantumOperator,
                                                 basis: List[Tuple[int, ...]]) -> sp.csr_matrix:
        """
        Vectorized fast path for build_sparse_operator_matrix, valid ONLY for
        statistics='boson' (no Jordan-Wigner phase to track, so amplitude
        factors are plain sqrt(occupation) products with no sign bookkeeping
        -- see _get_parity_phase, which is unconditionally +1 for bosons).

        Operates on the whole basis as a single (dim, n_modes) occupation
        array per (term, nonzero-index) entry, instead of walking a Python
        dict per basis column. This is the fix for the case that actually
        dominates runtime: an operator whose tensor is dense after a basis
        transformation (e.g. Kerr post-Bogoliubov-rotation), applied across
        a basis with hundreds of states -- the per-column dict-walking loop
        in _apply_precomputed_to_state is what costs seconds there, not the
        tensor scan itself (which _precompute_term_nonzeros already
        amortizes to once per term).

        Verified bit-identical to _apply_precomputed_to_state's per-column
        result before being used anywhere -- see test_kerr_performance.py.
        """
        assert self.statistics == 'boson', "vectorized fast path is boson-only"

        dim = len(basis)
        n = self.n_modes
        basis_arr = np.array(basis, dtype=np.int64)   # (dim, n)
        state_to_index = {state: idx for idx, state in enumerate(basis)}

        term_nonzeros = self._precompute_term_nonzeros(operator)

        rows, cols, data = [], [], []

        for daggers, entries in term_nonzeros:
            k = len(daggers)
            for idx_tuple, coeff in entries:
                # Batch-apply this single tensor entry's operator string to
                # every basis column at once.
                occ = basis_arr.copy()                       # (dim, n)
                amp = np.full(dim, coeff, dtype=complex)      # (dim,)
                alive = np.ones(dim, dtype=bool)              # rows still nonzero

                for i in reversed(range(k)):
                    mode = idx_tuple[i]
                    is_dagger = daggers[i]
                    if is_dagger:
                        amp[alive] *= np.sqrt(occ[alive, mode] + 1)
                        occ[alive, mode] += 1
                    else:
                        zero_mask = alive & (occ[:, mode] == 0)
                        alive[zero_mask] = False
                        amp[zero_mask] = 0.0
                        amp[alive] *= np.sqrt(occ[alive, mode])
                        occ[alive, mode] -= 1
                    if not np.any(alive):
                        break

                if not np.any(alive):
                    continue

                # Map surviving result states back to row indices; states
                # that fall outside the given basis contribute nothing
                # (matching the dict-based path's "if result_state in
                # state_to_index" guard).
                live_idx = np.where(alive)[0]
                for col_idx in live_idx:
                    result_state = tuple(occ[col_idx])
                    row_idx = state_to_index.get(result_state)
                    if row_idx is not None and abs(amp[col_idx]) > self.coeff_tol:
                        rows.append(row_idx)
                        cols.append(col_idx)
                        data.append(amp[col_idx])

        matrix = sp.coo_matrix((data, (rows, cols)), shape=(dim, dim), dtype=complex)
        return matrix.tocsr()

    def build_operator_matrix(self, operator: QuantumOperator, basis: List[Tuple[int, ...]]) -> np.ndarray:
        """
        Computes the dense matrix representation of an operator within a given basis.
        Rows correspond to the 'bra' (output state), Columns correspond to the 'ket' (input state).

        See _precompute_term_nonzeros: each term's tensor is scanned once,
        not once per basis column.
        """
        dim = len(basis)
        matrix = np.zeros((dim, dim), dtype=complex)
        state_to_index = {state: idx for idx, state in enumerate(basis)}
        term_nonzeros = self._precompute_term_nonzeros(operator)

        for col_idx, ket_state in enumerate(basis):
            current_total = self._apply_precomputed_to_state(term_nonzeros, ket_state)
            for result_state, amplitude in current_total.items():
                if result_state in state_to_index:
                    row_idx = state_to_index[result_state]
                    matrix[row_idx, col_idx] += amplitude

        return matrix

    def build_sparse_operator_matrix(self, operator: QuantumOperator, basis: List[Tuple[int, ...]]) -> sp.csr_matrix:
        """
        Computes the sparse matrix representation of an operator within a given basis.
        Returns a CSR format sparse matrix optimized for fast arithmetic and diagonalization.

        For statistics='boson', delegates to a vectorized batch-apply fast
        path (no Jordan-Wigner phase to track), which operates on the whole
        basis as numpy arrays instead of walking a Python dict per basis
        column -- the dominant cost for any operator whose tensor is dense
        after a basis transformation (e.g. a Kerr term post-Bogoliubov-
        rotation). Fermions fall back to the per-column path, which has to
        track the Jordan-Wigner sign per state.
        """
        if self.statistics == 'boson':
            return self._build_operator_matrix_vectorized_boson(operator, basis)

        dim = len(basis)
        
        # DOK (Dictionary of Keys) format is the fastest for constructing sparse matrices 
        # when you are adding elements one by one.
        matrix = sp.dok_matrix((dim, dim), dtype=complex)
        
        # Fast lookup dictionary
        state_to_index = {state: idx for idx, state in enumerate(basis)}
        term_nonzeros = self._precompute_term_nonzeros(operator)

        # Iterate over the basis (Columns/Kets)
        for col_idx, ket_state in enumerate(basis):
            current_total = self._apply_precomputed_to_state(term_nonzeros, ket_state)

            for result_state, amplitude in current_total.items():
                if abs(amplitude) <= self.coeff_tol:
                    continue
                if result_state in state_to_index:
                    row_idx = state_to_index[result_state]
                    matrix[row_idx, col_idx] += amplitude
                    
        # Convert to CSR format before returning. 
        # CSR is required for fast Scipy linear algebra operations (like eigsh).
        return matrix.tocsr()

    def build_rectangular_operator_matrix(self, operator: QuantumOperator,
                                           ket_basis: List[Tuple[int, ...]],
                                           bra_basis: List[Tuple[int, ...]]) -> sp.csr_matrix:
        """
        Computes the matrix of `operator` mapping ket_basis -> bra_basis when
        the two live in different excitation sectors (e.g. a creation-pair
        operator taking an n-excitation state to an (n+2)-excitation one).
        build_sparse_operator_matrix assumes a single shared basis for both
        rows and columns, which doesn't apply here -- this is the rectangular
        generalization, built the same way (apply_operator per ket column).
        """
        n_cols = len(ket_basis)
        n_rows = len(bra_basis)
        matrix = sp.dok_matrix((n_rows, n_cols), dtype=complex)
        state_to_row = {state: idx for idx, state in enumerate(bra_basis)}

        for col_idx, ket_state in enumerate(ket_basis):
            initial_superposition = Superposition({ket_state: 1.0 + 0j})
            result_superposition = self.apply_operator(operator, initial_superposition)
            for result_state, amplitude in result_superposition.states.items():
                if abs(amplitude) <= self.coeff_tol:
                    continue
                if result_state in state_to_row:
                    matrix[state_to_row[result_state], col_idx] += amplitude

        return matrix.tocsr()

    def build_fast_sparse_matrix(self, operator: QuantumOperator, basis: list) -> sp.csr_matrix:
        """Massively optimized sparse matrix builder using COO format."""
        dim = len(basis)
        state_to_index = {state: idx for idx, state in enumerate(basis)}
        
        operator.simplify() # Ensure redundant tensors are merged before loop
        
        rows, cols, data = [], [], []
        
        for col_idx, ket_state in enumerate(basis):
            initial_super = Superposition({ket_state: 1.0 + 0j})
            result_super = self.apply_operator(operator, initial_super)
            
            for result_state, amplitude in result_super.states.items():
                if abs(amplitude) <= self.coeff_tol:
                    continue
                if result_state in state_to_index:
                    rows.append(state_to_index[result_state])
                    cols.append(col_idx)
                    data.append(amplitude)
                    
        # Construct the COO matrix in a single vectorized burst, then convert to CSR
        coo = sp.coo_matrix((data, (rows, cols)), shape=(dim, dim), dtype=complex)
        return coo.tocsr()
    
    def exp_operator(self, operator: QuantumOperator, order: int) -> QuantumOperator:
        """
        Computes the Taylor expansion exp(A) = I + A + A^2/2! + ... + A^n/n!
        and returns the fully normal-ordered QuantumOperator.
        """
        # Identity operator: 0-order tensor, scalar 1.0
        identity_term = TensorTerm(np.array(1.0 + 0j), ())
        result_op = QuantumOperator([identity_term])
        
        if order == 0:
            return result_op
            
        current_power = operator
        result_op = result_op + (1.0 * current_power)
        
        for n in range(2, order + 1):
            # Multiply operator by itself to get next power
            current_power = current_power * operator
            
            # Crucial optimization: simplify before the next multiplication 
            # to prevent the tensor array sizes from bottlenecking RAM
            current_power.simplify() 
            
            term_to_add = (1.0 / math.factorial(n)) * current_power
            result_op = result_op + term_to_add
            
        result_op.simplify()
        
        # Pass the final summation through the Wick's theorem engine
        return self.normal_order(result_op)

    def _extract_quadratic_blocks(self, operator: QuantumOperator) -> Tuple[np.ndarray, np.ndarray]:
        """Extracts the A (b^dag b) and B (b^dag b^dag) matrices from a quadratic operator."""
        operator.simplify()
        is_real = all(np.isrealobj(term.tensor) for term in operator.terms)
        dtype = np.float64 if is_real else complex
        A = np.zeros((self.n_modes, self.n_modes), dtype=dtype)
        B = np.zeros((self.n_modes, self.n_modes), dtype=dtype)
        
        for term in operator.terms:
            if len(term.daggers) != 2:
                raise ValueError("Operator contains non-quadratic terms.")
                
            if term.daggers == (1, 0):    # b^dag_i b_j
                A += term.tensor
            elif term.daggers == (1, 1):  # b^dag_i b^dag_j
                B += term.tensor
                
        return A, B

    def diagonalize_quadratic(self, operator: QuantumOperator):
        A, B = self._extract_quadratic_blocks(operator)
        n = self.n_modes

        is_real = np.isrealobj(A) and np.isrealobj(B)
        dtype = np.float64 if is_real else complex

        if self.statistics == 'fermion':
            return self._diagonalize_quadratic_fermion(A, B, n, dtype)
        return self._diagonalize_quadratic_boson(A, B, n, dtype)

    def _diagonalize_quadratic_boson(self, A, B, n, dtype):
        # Build the Hermitian dynamical matrix. Conjugation is a no-op for
        # real input but kept dtype-aware: forcing complex unconditionally
        # here would route real problems through a different LAPACK routine
        # (zheevd vs dsyevd) than the original hand-rolled implementation,
        # which has no guaranteed agreement on eigenvector sign for
        # degenerate/near-degenerate clusters.
        D_herm = np.block([
            [A, B],
            [B.conj().T, A.conj()]
        ]).astype(dtype)
        
        tau_z = np.block([
            [np.eye(n), np.zeros((n, n))],
            [np.zeros((n, n)), -np.eye(n)]
        ]).astype(dtype)
        
        # Van Hemmen Cholesky decomposition
        L = np.linalg.cholesky(D_herm)
        M_sym = L.T @ tau_z @ L
        
        evals_sym, evecs_sym = la.eigh(M_sym)
        evecs_real = la.solve_triangular(L.T, evecs_sym)
        
        # Keep positive branch
        pos_idx = evals_sym > self.zero_tol
        pos_evals = evals_sym[pos_idx]
        pos_evecs = evecs_real[:, pos_idx]
        
        sort_idx = np.argsort(pos_evals)
        pos_evals = pos_evals[sort_idx]
        pos_evecs = pos_evecs[:, sort_idx]
        
        # Symplectic Gram-Schmidt
        def symplectic_gram_schmidt(vecs):
            n_cols = vecs.shape[1]
            result = np.zeros_like(vecs, dtype=vecs.dtype)
            for i in range(n_cols):
                v = vecs[:, i].copy()
                for j in range(i):
                    e_j = result[:, j]
                    proj = e_j.conj().T @ tau_z @ v
                    v = v - proj * e_j
                norm = np.real(v.conj().T @ tau_z @ v)
                result[:, i] = v / np.sqrt(norm)
            return result
            
        U_V = symplectic_gram_schmidt(pos_evecs)
        
        # Canonicalize sign: force the largest-magnitude component of each
        # column to be positive (real part, for the real-valued case). LAPACK
        # only guarantees eigenspaces up to sign/phase, not a specific choice,
        # so different call sites (or different dtypes/routines) can
        # otherwise disagree on sign for quantities that are odd in it
        # (e.g. the B_J/C_J pairing blocks used downstream).
        for col in range(U_V.shape[1]):
            v = U_V[:, col]
            idx = np.argmax(np.abs(v))
            ref = v[idx]
            if np.real(ref) < 0 or (np.real(ref) == 0 and np.imag(ref) < 0):
                U_V[:, col] = -v
        
        U = U_V[:n, :]
        V = U_V[n:, :]
        
        return pos_evals, U, V

    def _diagonalize_quadratic_fermion(self, A, B, n, dtype):
        """
        Van Hemmen's fermion procedure (Sect. 3 of the 1980 paper): because
        the anticommutator algebra forces D = -A^t, B = -B^t, C = -C^t (Eq. 10)
        and Hermiticity forces A^+ = A, B^+ = -C (Eq. 11), the dynamical
        matrix D = [[A, B], [-B_conj, -A_conj]] (Eq. 12) is directly Hermitian
        -- no symplectic/Cholesky machinery is needed, unlike bosons. A
        canonical T is simply a unitary matrix diagonalizing D (Eq. 19:
        T^+ D T = Omega), built from an orthonormal eigenbasis. Unlike the
        boson case, fermion eigenvalues need not be positive -- there is no
        analogous stability requirement (Sect. 5). Eigenvalues come in
        +/-omega pairs via the J operator (Eq. 53), and one representative
        per pair is selected by block dominance (see below), recovering a
        complete set of n quasiparticle modes.
        """
        D_herm = np.block([
            [A, B],
            [-B.conj(), -A.conj()]
        ]).astype(dtype)

        evals, evecs = la.eigh(D_herm)

        # Eigenvalues of D_herm come in +/-omega pairs (Eq. 53: D x = wx =>
        # D(Jx) = -w(Jx), with J(u,v)=(v_bar,u_bar)), but eigenvalue sign is
        # NOT a reliable way to pick one representative per pair -- e.g. for
        # B=0, D is block-diagonal as [[A,0],[0,-A]], so a negative
        # eigenvalue of A produces a genuinely "upper-block" eigenvector with
        # a negative eigenvalue, indistinguishable by sign alone from its
        # lower-block J-partner. The reliable invariant is which of the two
        # n-blocks the eigenvector actually lives in: select the n
        # eigenvectors whose upper-block (creation/particle) component
        # dominates the lower-block (annihilation/hole) component.
        upper_norm_sq = np.sum(np.abs(evecs[:n, :]) ** 2, axis=0)
        lower_norm_sq = np.sum(np.abs(evecs[n:, :]) ** 2, axis=0)
        keep = np.where(upper_norm_sq > lower_norm_sq)[0]

        if len(keep) != n:
            # Degenerate/borderline cases (upper_norm == lower_norm, e.g. an
            # exactly-zero mode) can misassign; fall back to whichever n
            # indices have the larger margin.
            margin = upper_norm_sq - lower_norm_sq
            order = np.argsort(margin)[::-1]
            keep = np.sort(order[:n])

        pos_evals = evals[keep]
        T_cols = evecs[:, keep]

        sort_idx = np.argsort(pos_evals)
        pos_evals = pos_evals[sort_idx]
        T_cols = T_cols[:, sort_idx]

        # T is already unitary (eigh returns an orthonormal eigenbasis), so
        # no Gram-Schmidt/normalization step is needed -- this is the direct
        # simplification the fermion case affords over the boson one.
        for col in range(T_cols.shape[1]):
            v = T_cols[:, col]
            idx = np.argmax(np.abs(v))
            ref = v[idx]
            if np.real(ref) < 0 or (np.real(ref) == 0 and np.imag(ref) < 0):
                T_cols[:, col] = -v

        U = T_cols[:n, :]
        V = T_cols[n:, :]

        return pos_evals, U, V
    

    def build_squeezed_vacuum(self, U: np.ndarray, V: np.ndarray, truncation: int, threshold: float = 0.99) -> Superposition:
        """
        Generates the exact analytical squeezed vacuum state for a quadratic Hamiltonian.

        For statistics='boson', delegates to a vectorized fast path (see
        _build_squeezed_vacuum_vectorized_boson): the per-level term count
        grows combinatorially (e.g. 1 -> 78 -> 1365 -> 12376 -> 75582 for a
        12-mode system by truncation level 4), and the original dict-based
        loop rescans M's full n x n entries for every single term at every
        level -- by level 4 that is tens of millions of Python-level dict
        operations. Fermions fall back to the original path (not addressed
        here; the M matrix's structure/meaning differs for fermions and
        this fast path has only been verified for the boson case).
        """
        if self.statistics == 'boson':
            return self._build_squeezed_vacuum_vectorized_boson(U, V, truncation, threshold)

        n = self.n_modes
        
        # RECONSTRUCT U_V HERE:
        U_V = np.vstack([U, V])
        
        # Build the full pseudo-unitary transformation matrix to find the inverse
        tau_x = np.block([[np.zeros((n, n)), np.eye(n)], [np.eye(n), np.zeros((n, n))]])
        J_U_V = tau_x @ U_V.conj() # Pseudo-time reversal for the negative branch
        P = np.hstack([U_V, J_U_V])
        P_inv = np.linalg.inv(P)
        
        U_tilde = P_inv[:n, :n]
        V_tilde = P_inv[:n, n:]
        M = np.linalg.solve(U_tilde, V_tilde)
        
        vacuum = tuple([0] * n)
        current_terms = {vacuum: 1.0 + 0j}

        final_state = defaultdict(complex)
        final_state[vacuum] = 1.0 + 0j
        
        for k in range(1, truncation + 1):
            next_terms = defaultdict(complex)
            for state_tuple, coeff in current_terms.items():
                state = list(state_tuple)
                for j in range(n):
                    for i in range(n):
                        val = M[i, j]
                        if abs(val) < self.coeff_tol:
                            continue
                            
                        new_state = state.copy()
                        new_state[j] += 1
                        factor_j = np.sqrt(new_state[j])
                        new_state[i] += 1
                        factor_i = np.sqrt(new_state[i])
                        
                        next_terms[tuple(new_state)] += coeff * (-0.5) * val * factor_j * factor_i / float(k)
                        
            for st, val in next_terms.items():
                if abs(val) > self.coeff_tol:
                    final_state[st] += val
            current_terms = next_terms
            
        return Superposition(final_state).trim(threshold)

    def _build_squeezed_vacuum_vectorized_boson(self, U: np.ndarray, V: np.ndarray,
                                                  truncation: int, threshold: float = 0.99) -> Superposition:
        """
        Vectorized fast path for build_squeezed_vacuum, boson-only (pure
        sqrt(occupation) factors, no Jordan-Wigner phase). At each
        truncation level, represents the current term set as a
        (num_states, n) occupation array + (num_states,) coefficient array,
        applies every nonzero (i,j) pair of M to the WHOLE array at once via
        broadcasting, then deduplicates resulting states (different (state,
        i, j) combinations can land on the same output state) by summing
        coefficients per unique row -- this replaces the dict's automatic
        accumulation, which is the one piece that doesn't trivially
        vectorize away.

        Verified bit-identical (state-for-state) to the original dict-based
        path -- see test_vectorized_squeezed_vacuum.py.
        """
        n = self.n_modes

        U_V = np.vstack([U, V])
        tau_x = np.block([[np.zeros((n, n)), np.eye(n)], [np.eye(n), np.zeros((n, n))]])
        J_U_V = tau_x @ U_V.conj()
        P = np.hstack([U_V, J_U_V])
        P_inv = np.linalg.inv(P)
        U_tilde = P_inv[:n, :n]
        V_tilde = P_inv[:n, n:]
        M = np.linalg.solve(U_tilde, V_tilde)

        # Precompute M's nonzero (i, j, val) triples once, not once per state.
        nz_i, nz_j = np.where(np.abs(M) > self.coeff_tol)
        nz_val = M[nz_i, nz_j]

        vacuum = tuple([0] * n)
        final_state = defaultdict(complex)
        final_state[vacuum] = 1.0 + 0j

        current_states = np.zeros((1, n), dtype=np.int64)
        current_coeffs = np.array([1.0 + 0j])

        # Hoist loop-invariant arrays to device once if GPU is enabled.
        # nz_i/nz_j/nz_val don't change across truncation levels.
        if self.use_gpu:
            d_nz_i   = cp.asarray(nz_i)
            d_nz_j   = cp.asarray(nz_j)
            d_nz_val = cp.asarray(nz_val)

        for k in range(1, truncation + 1):
            num_states = current_states.shape[0]
            num_pairs = len(nz_val)

            if num_states == 0 or num_pairs == 0:
                current_states = np.zeros((0, n), dtype=np.int64)
                current_coeffs = np.zeros((0,), dtype=complex)
                continue

            if self.use_gpu:
                current_states, current_coeffs = self._gpu_level_step(
                    current_states, current_coeffs,
                    d_nz_i, d_nz_j, d_nz_val,
                    k, n
                )
            else:
                current_states, current_coeffs = self._cpu_level_step(
                    current_states, current_coeffs,
                    nz_i, nz_j, nz_val,
                    k, n
                )

            for st, val in zip(current_states, current_coeffs):
                final_state[tuple(int(x) for x in st)] += val

        return Superposition(final_state).trim(threshold)

    def _cpu_level_step(self, current_states, current_coeffs, nz_i, nz_j, nz_val,
                         k, n):
        """One truncation level on CPU. Extracted from the original loop body
        so the GPU path can share the same return contract without code duplication."""
        num_states = current_states.shape[0]
        num_pairs  = len(nz_val)

        new_states = np.broadcast_to(current_states, (num_pairs, num_states, n)).copy()
        rows = np.arange(num_pairs)

        new_states[rows, :, nz_j] += 1
        factor_j = np.sqrt(new_states[rows, :, nz_j].astype(float))
        new_states[rows, :, nz_i] += 1
        factor_i = np.sqrt(new_states[rows, :, nz_i].astype(float))

        new_coeffs = (current_coeffs[None, :] * (-0.5) * nz_val[:, None]
                      * factor_j * factor_i / float(k))

        flat_states = new_states.reshape(num_pairs * num_states, n)
        flat_coeffs = new_coeffs.reshape(num_pairs * num_states)

        max_occ = int(flat_states.max()) if flat_states.size > 0 else 0
        base    = max_occ + 2

        if n * np.log2(base) >= 62:
            # Overflow guard: fall back to row-wise unique
            unique_states, inverse = np.unique(flat_states, axis=0, return_inverse=True)
            inverse = inverse.reshape(-1)
            summed_coeffs = np.zeros(unique_states.shape[0], dtype=complex)
            np.add.at(summed_coeffs, inverse, flat_coeffs)
            keep_mask   = np.abs(summed_coeffs) > self.coeff_tol
            next_states = unique_states[keep_mask]
            next_coeffs = summed_coeffs[keep_mask]
        else:
            powers = base ** np.arange(n - 1, -1, -1, dtype=np.int64)
            packed = flat_states @ powers

            unique_packed, inverse = np.unique(packed, return_inverse=True)
            inverse = inverse.reshape(-1)
            summed_coeffs = np.zeros(unique_packed.shape[0], dtype=complex)
            np.add.at(summed_coeffs, inverse, flat_coeffs)

            keep_mask   = np.abs(summed_coeffs) > self.coeff_tol
            kept_packed = unique_packed[keep_mask]
            next_coeffs = summed_coeffs[keep_mask]

            next_states = np.zeros((kept_packed.shape[0], n), dtype=np.int64)
            remaining   = kept_packed.copy()
            for mode in reversed(range(n)):
                next_states[:, mode] = remaining % base
                remaining //= base

        return next_states, next_coeffs

    def _gpu_level_step(self, current_states, current_coeffs,
                         d_nz_i, d_nz_j, d_nz_val, k, n):
        """One truncation level on GPU. Data arrives as numpy; nz_i/nz_j/nz_val
        are already on-device (hoisted by the caller). Single host->device
        transfer for states/coeffs at the start, single device->host at the end.

        cp.add.at does not support complex128 on current cupy versions. Split
        into float64 real/imag parts for accumulation (mathematically exact
        since addition never mixes real and imaginary parts), then recombine.
        """
        num_states = current_states.shape[0]
        num_pairs  = len(d_nz_val)

        # Transfer current level's states/coeffs to device
        cs = cp.asarray(current_states)
        cc = cp.asarray(current_coeffs)

        new_states = cp.broadcast_to(cs, (num_pairs, num_states, n)).copy()
        rows = cp.arange(num_pairs)

        new_states[rows, :, d_nz_j] += 1
        factor_j = cp.sqrt(new_states[rows, :, d_nz_j].astype(cp.float64))
        new_states[rows, :, d_nz_i] += 1
        factor_i = cp.sqrt(new_states[rows, :, d_nz_i].astype(cp.float64))

        new_coeffs = (cc[None, :] * (-0.5) * d_nz_val[:, None]
                      * factor_j * factor_i / float(k))

        flat_states = new_states.reshape(num_pairs * num_states, n)
        flat_coeffs = new_coeffs.reshape(num_pairs * num_states)

        max_occ = int(flat_states.max())
        base    = max_occ + 2

        powers = base ** cp.arange(n - 1, -1, -1, dtype=cp.int64)
        # cp.matmul on int64 hits CUDA_ERROR_INVALID_VALUE for large arrays;
        # elementwise multiply + sum is equivalent and avoids the integer matmul kernel.
        packed = (flat_states * powers[None, :]).sum(axis=1)

        unique_packed, inverse = cp.unique(packed, return_inverse=True)
        inverse = inverse.reshape(-1)

        # Split complex into real + imag for cp.add.at (no complex128 support)
        real_part = cp.zeros(unique_packed.shape[0], dtype=cp.float64)
        imag_part = cp.zeros(unique_packed.shape[0], dtype=cp.float64)
        cp.add.at(real_part, inverse, flat_coeffs.real)
        cp.add.at(imag_part, inverse, flat_coeffs.imag)
        summed_coeffs = real_part + 1j * imag_part

        keep_mask   = cp.abs(summed_coeffs) > self.coeff_tol
        kept_packed = unique_packed[keep_mask]
        next_coeffs = summed_coeffs[keep_mask]

        next_states = cp.zeros((kept_packed.shape[0], n), dtype=cp.int64)
        remaining   = kept_packed.copy()
        for mode in reversed(range(n)):
            next_states[:, mode] = remaining % base
            remaining //= base

        # Single transfer back to host for this level
        return cp.asnumpy(next_states), cp.asnumpy(next_coeffs)


    def transform_state(self, state: Superposition, U: np.ndarray, V: np.ndarray,
                         direction: str = 'old_to_new',
                         base_state: 'Superposition' = None) -> Superposition:
        """
        Transforms a state vector between the bare Fock basis and the Bogoliubov
        eigenbasis, by repeated application of a single quasiparticle creation
        operator (built from rows of U, V) to a base/reference state.

        direction='old_to_new': U, V are used as-given (rows = per-mode
        coefficients), and the base state defaults to the analytic squeezed
        vacuum built from U, V.

        direction='new_to_old' (or anything else): U, V are used as-given too
        -- the caller is responsible for passing whatever U, V already
        correctly encode the inverse transform (e.g. P[:n,:n].T, -P[n:,:n].T
        for a real pseudo-unitary P), since re-deriving P^-1 from U, V here
        assumes a specific complex Bogoliubov convention (P = [[U,V*],[V,U*]])
        that does not match every caller's P construction. The base state
        defaults to the bare vacuum.

        base_state: optional explicit starting Superposition (e.g. a
        precomputed squeezed vacuum) to use instead of either default. This
        is what lets callers reuse a cached vacuum (self.vacuum,
        self.unsqueezed_vacuum) instead of recomputing it via
        build_squeezed_vacuum on every call.

        Returns the untrimmed result (only the 1e-12 numerical-noise cleanup
        from Superposition arithmetic is applied). Callers that want a
        threshold-based truncation should call trim_state()/trim() explicitly
        with their own threshold -- this method no longer trims internally,
        since a hardcoded internal threshold silently overrides whatever
        threshold the caller actually wants, and double-trimming (once here,
        once by the caller) corrupts the result.
        """
        n = self.n_modes
        U_mat, V_mat = U, V
        if base_state is not None:
            base_vacuum = base_state
        elif direction == 'old_to_new':
            base_vacuum = self.build_squeezed_vacuum(U, V, truncation=10) # Approximated base
        else:
            base_vacuum = Superposition({tuple([0] * n): 1.0 + 0j})
            
        result = Superposition()
        
        for st_tuple, amp in state.states.items():
            occupations = np.array(st_tuple)
            current_branch = base_vacuum
            
            for mode_idx in range(n):
                r_i = occupations[mode_idx]
                if r_i == 0:
                    continue
                    
                for k in range(1, r_i + 1):
                    # Apply the Bogoliubov transformed creation operator
                    next_branch = defaultdict(complex)
                    for br_st, br_amp in current_branch.states.items():
                        # U component (Creation)
                        for j in range(n):
                            if abs(U_mat[mode_idx, j]) > self.coeff_tol:
                                new_st, factor = self._apply_creation(j, br_st)
                                if new_st:
                                    next_branch[new_st] += br_amp * U_mat[mode_idx, j] * factor
                        # V component (Annihilation)
                        for j in range(n):
                            if abs(V_mat[mode_idx, j]) > self.coeff_tol:
                                new_st, factor = self._apply_annihilation(j, br_st)
                                if new_st:
                                    next_branch[new_st] += br_amp * V_mat[mode_idx, j] * factor
                                    
                    # Normalize factorial scaling
                    current_branch = Superposition({s: v / np.sqrt(k) for s, v in next_branch.items()})
                    
            result = result + (amp * current_branch)
            
        return result
    
    def build_quadratic_sparse_matrix(self, operator: QuantumOperator, basis: list) -> sp.csr_matrix:
        """
        Highly optimized sparse matrix builder specifically for quadratic Hamiltonians.
        Bypasses generalized tensor contractions for massive speedup.
        """
        A, B = self._extract_quadratic_blocks(operator)
        n = self.n_modes
        dim = len(basis)
        state_to_index = {state: idx for idx, state in enumerate(basis)}
        
        rows, cols, data = [], [], []
        
        def add_element(r, c, val):
            if abs(val) > self.coeff_tol:
                rows.append(r)
                cols.append(c)
                data.append(val)
                
        for col_idx, ket in enumerate(basis):
            ket_arr = np.array(ket)
            
            # Number-conserving terms: A (b^dag_i b_j)
            for j in range(n):
                if ket_arr[j] == 0:
                    continue
                for i in range(n):
                    if abs(A[i, j]) < self.coeff_tol:
                        continue
                    new_st = ket_arr.copy()
                    new_st[j] -= 1
                    factor_j = np.sqrt(ket_arr[j])
                    factor_i = np.sqrt(new_st[i] + 1)
                    new_st[i] += 1
                    
                    st_tuple = tuple(new_st)
                    if st_tuple in state_to_index:
                        add_element(state_to_index[st_tuple], col_idx, A[i, j] * factor_i * factor_j)
                        
            # Pairing terms: B (b^dag_i b^dag_j)
            for j in range(n):
                for i in range(n):
                    if abs(B[i, j]) < self.coeff_tol:
                        continue
                    new_st = ket_arr.copy()
                    new_st[j] += 1
                    factor_j = np.sqrt(new_st[j])
                    new_st[i] += 1
                    factor_i = np.sqrt(new_st[i])
                    
                    st_tuple = tuple(new_st)
                    if st_tuple in state_to_index:
                        add_element(state_to_index[st_tuple], col_idx, B[i, j] * factor_i * factor_j)
                        
            # Conjugate Pairing terms: B^* (b_i b_j)
            B_conj = B.conj().T
            for j in range(n):
                if ket_arr[j] == 0:
                    continue
                for i in range(n):
                    if abs(B_conj[i, j]) < self.coeff_tol:
                        continue
                    new_st = ket_arr.copy()
                    factor_j = np.sqrt(new_st[j])
                    new_st[j] -= 1
                    if new_st[i] == 0:
                        continue
                    factor_i = np.sqrt(new_st[i])
                    new_st[i] -= 1
                    
                    st_tuple = tuple(new_st)
                    if st_tuple in state_to_index:
                        add_element(state_to_index[st_tuple], col_idx, B_conj[i, j] * factor_i * factor_j)

        return sp.coo_matrix((data, (rows, cols)), shape=(dim, dim), dtype=complex).tocsr()
    
    def partial_trace(self, state: Superposition, keep_modes: List[int]) -> np.ndarray:
        """
        Computes the reduced density matrix by tracing out all modes NOT in `keep_modes`.
        """
        keep_modes = sorted(keep_modes)
        trace_modes = [i for i in range(self.n_modes) if i not in keep_modes]
        
        # Determine the dimension of the kept subspace dynamically based on the state
        max_excitations = {m: 0 for m in keep_modes}
        for st in state.states.keys():
            for m in keep_modes:
                max_excitations[m] = max(max_excitations[m], st[m])
                
        # Create an index mapping for the kept basis
        import itertools
        ranges = [range(max_excitations[m] + 1) for m in keep_modes]
        kept_basis = list(itertools.product(*ranges))
        kept_idx = {st: i for i, st in enumerate(kept_basis)}
        dim = len(kept_basis)
        
        rho_reduced = np.zeros((dim, dim), dtype=complex)
        
        # Group states by their traced-out mode configuration
        trace_groups = defaultdict(list)
        for st_tuple, amp in state.states.items():
            traced_config = tuple(st_tuple[m] for m in trace_modes)
            kept_config = tuple(st_tuple[m] for m in keep_modes)
            trace_groups[traced_config].append((kept_config, amp))
            
        # Perform the partial trace: sum over outer products of matching traced configurations
        for traced_config, kept_components in trace_groups.items():
            for (kept_1, amp_1) in kept_components:
                for (kept_2, amp_2) in kept_components:
                    i = kept_idx[kept_1]
                    j = kept_idx[kept_2]
                    rho_reduced[i, j] += amp_1 * np.conj(amp_2)
                    
        return rho_reduced
    
    def get_mode_occupations(self, state: Superposition) -> np.ndarray:
        """Calculates the expected photon/excitation number <n_i> for each mode."""
        occupations = np.zeros(self.n_modes, dtype=float)
        for st_tuple, amp in state.states.items():
            prob = abs(amp)**2
            for i in range(self.n_modes):
                occupations[i] += st_tuple[i] * prob
        return occupations

    def get_total_excitations(self, state: Superposition) -> float:
        """Returns the total expected number of excitations in the system."""
        return np.sum(self.get_mode_occupations(state))

    def build_basis_change_matrix(self, source_states: List[Tuple[int, ...]],
                                   expand_fn, threshold: float = 0.99):
        """
        Generic, direction-agnostic basis-change matrix builder. The engine
        does not store or know about any specific transformation (Bogoliubov,
        old<->new, or otherwise) -- it only assembles a matrix out of
        whatever single-state expansion function the caller supplies.

        source_states: the basis whose states become the matrix's COLUMNS,
        in the given order.

        expand_fn: a callable, state_tuple -> {target_state: amplitude},
        e.g. a single Fock state's decomposition in some other basis. This
        is supplied by the caller (the chain-specific logic that knows how
        to do the expansion); the engine has no opinion about its meaning.

        The target/row basis is NOT given in advance -- it is discovered as
        the union of every state appearing in any column's expansion, since
        expand_fn can return a different number of nonzero target states for
        each source state (some columns may be sparser/denser than others).
        Rows are ordered deterministically by (total excitation, state
        tuple) for reproducibility.

        Returns (matrix, row_states, row_index, col_states):
          matrix:     sparse CSR matrix, shape (len(row_states), len(source_states))
          row_states: list of discovered target-basis Fock states, in row order
          row_index:  dict mapping each row_states entry to its row index
          col_states: source_states, returned unchanged for convenience so
                      callers always have both bases tracked alongside the
                      matrix itself.

        Each column is whatever expand_fn returns for that source state,
        un-modified (e.g. already normalized/truncated by expand_fn itself,
        if that's what the caller's expand_fn does) -- this method does not
        re-normalize or re-truncate on its own.
        """
        col_states = list(source_states)
        column_expansions = [expand_fn(state) for state in col_states]

        all_target_states = set()
        for expansion in column_expansions:
            all_target_states.update(expansion.keys())

        row_states = sorted(all_target_states, key=lambda s: (sum(s), s))
        row_index = {state: idx for idx, state in enumerate(row_states)}

        rows, cols, data = [], [], []
        for col_idx, expansion in enumerate(column_expansions):
            for target_state, amplitude in expansion.items():
                if abs(amplitude) <= self.coeff_tol:
                    continue
                rows.append(row_index[target_state])
                cols.append(col_idx)
                data.append(amplitude)

        matrix = sp.coo_matrix(
            (data, (rows, cols)),
            shape=(len(row_states), len(col_states)),
            dtype=complex
        ).tocsr()

        return matrix, row_states, row_index, col_states
    

class CompositeEngine:
    """A general-purpose engine for bipartite tensor-product Hilbert spaces H_A ⊗ H_B."""
    
    def vector_to_subsystems(self, full_vector: np.ndarray, basis_A: list, basis_B: list, threshold: float = 1e-12) -> dict:
        """
        Maps a 1D state vector of H_A ⊗ H_B into a dictionary of {st_A: Superposition_B}.
        """
        dim_A = len(basis_A)
        dim_B = len(basis_B)
        
        if full_vector.shape[0] != dim_A * dim_B:
            raise ValueError(f"Vector size {full_vector.shape[0]} does not match basis dimensions {dim_A}x{dim_B}")
            
        reshaped = full_vector.reshape((dim_A, dim_B))
        state_dict = {}
        
        for i, st_A in enumerate(basis_A):
            matter_dict = {}
            for j, st_B in enumerate(basis_B):
                amp = reshaped[i, j]
                if abs(amp) > threshold:
                    matter_dict[st_B] = amp
            if matter_dict:
                state_dict[st_A] = Superposition(matter_dict)
                
        return state_dict

    def subsystems_to_vector(self, state_dict: dict, basis_A: list, basis_B: list) -> np.ndarray:
        """
        Flattens a dictionary of {st_A: Superposition_B} back into a 1D complex numpy array.
        """
        dim_A = len(basis_A)
        dim_B = len(basis_B)
        
        # Fast reverse lookups
        idx_A = {st: i for i, st in enumerate(basis_A)}
        idx_B = {st: i for i, st in enumerate(basis_B)}
        
        psi = np.zeros((dim_A, dim_B), dtype=complex)
        
        for st_A, superpos in state_dict.items():
            if st_A not in idx_A: 
                continue
            i = idx_A[st_A]
            for st_B, amp in superpos.states.items():
                if st_B in idx_B:
                    j = idx_B[st_B]
                    psi[i, j] += amp
                    
        return psi.reshape(dim_A * dim_B)

    def map_subsystem_B(self, state_dict: dict, map_func) -> dict:
        """
        Applies a mapping function (e.g., Bogoliubov transforms) to the Subsystem B Superposition.
        """
        new_dict = {}
        for st_A, superpos_B in state_dict.items():
            mapped_superpos = map_func(superpos_B)
            if mapped_superpos.states: # Only store if not empty
                new_dict[st_A] = mapped_superpos
        return new_dict

    def partial_trace(self, full_vector: np.ndarray, dim_A: int, dim_B: int, keep: str = 'B') -> np.ndarray:
        """
        Computes the reduced density matrix by tracing out one subsystem.
        """
        psi = full_vector.reshape(dim_A, dim_B)
        if keep == 'B':
            return psi.conj().T @ psi
        elif keep == 'A':
            return psi @ psi.conj().T
        else:
            raise ValueError("keep must be 'A' or 'B'")

    def embed_matrix_A(self, matrix_A: sp.csr_matrix, dim_B: int) -> sp.csr_matrix:
        I_B = sp.eye(dim_B, format='csr')
        return sp.kron(matrix_A, I_B, format='csr')

    def embed_matrix_B(self, dim_A: int, matrix_B: sp.csr_matrix) -> sp.csr_matrix:
        I_A = sp.eye(dim_A, format='csr')
        return sp.kron(I_A, matrix_B, format='csr')

    def kron_matrices(self, matrix_A: sp.csr_matrix, matrix_B: sp.csr_matrix) -> sp.csr_matrix:
        return sp.kron(matrix_A, matrix_B, format='csr')