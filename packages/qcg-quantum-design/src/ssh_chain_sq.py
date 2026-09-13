"""ssh_chain_sq -- SSH photonic chain with a second-quantization engine.

Single-file merge of the group's ``sq_engine.py`` (SecondQuantizationEngine +
Bogoliubov machinery) and ``ssh_chain_sq.py`` (the PHOTONICChain that uses it),
from ``notebooks/exploration/unitary/``. Identical code, one import::

    from ssh_chain_sq import PHOTONICChain, SecondQuantizationEngine

CuPy is used automatically when importable; otherwise everything runs on NumPy.

Public names: TensorTerm, Superposition, QuantumOperator,
SecondQuantizationEngine, CompositeEngine, PHOTONICChain,
compute_position_exponential_sparse, compositions.
"""
import math
import time as _time
import warnings
import itertools
from itertools import product
from collections import defaultdict
from dataclasses import dataclass
from typing import Tuple, List, Dict

import numpy as np
import scipy.sparse as sp
import scipy.linalg as la
import scipy.sparse.linalg as sla
from scipy.special import genlaguerre, gammaln

try:
    import cupy as cp
    _CUPY_AVAILABLE = True
except ImportError:
    _CUPY_AVAILABLE = False


# ===========================================================================
# Part 1 -- second-quantization engine        (was sq_engine.py)
# ===========================================================================
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
            # rho_B[b, b'] = sum_a psi[a, b] conj(psi[a, b'])   (NOT psi^dag psi,
            # which is the transpose/conjugate of that).
            return psi.T @ psi.conj()
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


# ===========================================================================
# Part 2 -- SSH photonic chain               (was ssh_chain_sq.py)
# ===========================================================================

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
        # +/- Nambu convention: factor 2 for bosons, 1 for fermions (the lower
        # block holds -energies there, not +energies). See _nambu_prefactor.
        unperturbed_energies = self._nambu_prefactor() * np.diag(self.H_C)[:n].real   # (2N,) positive-branch energies

        diag_vals = states @ unperturbed_energies # (n_total,) — one energy per state
        mask      = np.abs(diag_vals) > self.engine.zero_tol
        idx       = np.where(mask)[0]
        if self.kerr_operator is not None and n_total > 2 * self.N: # if there is a kerr_operator available and we are looking at excitation spaces higher than or equal to two excitations.
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
        
        if transform is None:
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
            
    def compute_coupling_vec_from_sites(self, g_sites : np.ndarray) -> np.ndarray:
        """Transform an arbitrary per-site coupling vector into the chain's
        diagonalizing (eta) basis.

        Given a resonator coupled as (a^dag - a) sum_j g_sites[j] (b^dag_j - b_j)
        in the bare site basis, returns g'_j such that after the chain's
        Bogoliubov diagonalization the coupling reads
        (a^dag - a) sum_j g'_j (eta^dag_j - eta_j).

        Args:
            g_sites (np.ndarray): length-(2N) per-site coupling coefficients.
        """
        g_sites = np.asarray(g_sites, dtype=float)
        if g_sites.shape != (self.n,):
            raise ValueError(f"g_sites must have shape ({self.n},), got {g_sites.shape}")

        initial_coupling = QuantumOperator([
            TensorTerm(tensor=g_sites.copy(),  daggers=(1,)),
            TensorTerm(tensor=-g_sites.copy(), daggers=(0,)),   # anti-hermitian: couples to a^dag - a
        ])
        final_coupling = self.engine.change_basis(initial_coupling, self.U, self.V).simplify()
        # A physically zero input coupling correctly leaves no surviving terms
        # after change_basis/simplify prune away the zero tensor -- that is not
        # an error, it just means the coupling vector is all zeros.
        if not final_coupling.terms:
            return np.zeros(self.n)
        # first term is (1,); the (0,) term is its h.c. The Bogoliubov
        # transform is effectively real (the Kerr path just routes it through
        # the complex LAPACK branch), so g' is real up to numerical noise.
        return np.asarray(final_coupling.terms[0].tensor).real.astype(float)

    def compute_coupling_vec(self, g0 : float, ground_mode : float, mode : int, positions : np.ndarray) -> np.ndarray:
        """Get coupling constants from a resonator to the diagonalizing modes of the chain

        Standing-wave parametrization: g_sites[j] = g0 * sin(2 pi (1+mode) x_j / (L * ground_mode)).
        For a fully general per-site vector use compute_coupling_vec_from_sites.

        Args:
            g0 (float): initial coupling magnitude determined by capacitance and modes fluctuations
            ground_mode (float): ratio of first mode wavelength to resonator length lambda_0/L
            mode (int): mode order with respect to the ground mode
            positions (np.ndarray): normalized positions of the qubits along the resonator x_i / L
        """
        phases = 2.0 * np.pi * (1 + mode) * np.asarray(positions) / ground_mode
        return self.compute_coupling_vec_from_sites(g0 * np.sin(phases))

    # ------------------------------------------------------------------
    # Coupling of the whole (already-diagonalized) chain to a NEW resonator
    # ------------------------------------------------------------------
    #
    # The constructor has already Bogoliubov-diagonalized the chain's
    # quadratic part, so in the eta basis the chain reads
    #
    #     H_chain = sum_j omega_j eta^dag_j eta_j
    #               + sum_{ijkl} V_{ijkl} eta^dag_i eta^dag_j eta_k eta_l
    #
    # with omega_j = prefactor * diag(self.H_C)[:n] (self.H_C is diagonal
    # here) and V = self.kerr_operator.tensor, daggers (1,1,0,0). We now add a
    # new resonator alpha of frequency omega_r_new, coupled as
    #
    #     (alpha^dag - alpha) sum_j g'_j (eta^dag_j - eta_j),   g'_j from
    #     compute_coupling_vec.
    #
    # The methods below, assuming the chain is NOT coupled to the original
    # cavity (so self.U/self.V/self.H_C describe the whole quadratic part):
    #   1. compute_resonator_mode_transformation -- diagonalize the (n+1)-mode
    #      quadratic part into modes beta_k and pick beta_0, the one whose
    #      frequency is closest to omega_r_new (unambiguous in the dispersive
    #      regime),
    #   2. transform_kerr_to_resonator_basis -- push V through that second
    #      Bogoliubov transform -> V'_{ijkl} in the beta basis,
    #   3. _resonator_number_kerr_matrix -- keep only the beta_0^dag beta_0 -
    #      weighted quadratic chain piece,
    #         (sum_ij V''_ij beta^dag_i beta_j) beta_0^dag beta_0,
    #      dropping every other term (we work in the chain's single-excitation
    #      subspace),
    #   4. _diagonalize_degenerate_blocks -- diagonalize V'' inside each
    #      degenerate chain subspace so step 5 never divides by zero,
    #   5. compute_resonator_dispersive_shift_vec -- second-order shift
    #         chi_i = V''_ii + sum_{j != i} |V''_ij|^2 / (E_i - E_j).

    def _nambu_prefactor(self):
        """diag(self.H_C)[:n] carries +energy in both Nambu blocks for bosons
        (double it to get the physical frequency) but +/-energy for fermions
        (no doubling). Use this everywhere a physical mode frequency is read
        back off H_C -- build_HC_fock, compute_dispersive_shift_vec, and the
        new-resonator machinery."""
        return 1.0 if self.statistics == 'fermion' else 2.0

    def _extended_engine(self):
        """A SecondQuantizationEngine sized for the chain plus one extra mode."""
        ext = getattr(self, "_ext_engine", None)
        if ext is None or ext.n_modes != self.n + 1:
            ext = SecondQuantizationEngine(
                n_modes=self.n + 1,
                statistics=self.statistics,
                coeff_tol=self.engine.coeff_tol,
                zero_tol=self.engine.zero_tol,
            )
            self._ext_engine = ext
        return ext

    def _diagonalize_extended_quadratic(self, A, B):
        """Bogoliubov-diagonalize a quadratic form on n+1 modes from its (1,0)
        block A and (1,1) block B, following the exact convention of the
        constructor's internal ``diagonalize`` closure: Nambu matrix
        D_herm = [[A, B], [B, A]] (boson) or [[A, B], [-B*, -A*]] (fermion),
        H in the new basis read back as P^dag D_herm P.

        Returns (U, V, P, inverse_P, H_new), U/V shape (n+1, n+1).
        """
        ext = self._extended_engine()
        if self.statistics == 'boson':
            D_herm = np.block([[A, B], [B, A]])
        else:
            D_herm = np.block([[A, B], [-B.conj(), -A.conj()]])

        op = QuantumOperator([TensorTerm(A, (1, 0)), TensorTerm(B, (1, 1))])
        _, U, V = ext.diagonalize_quadratic(op)
        U_V = np.vstack([U, V])
        J_U_V = np.vstack([V.conj(), U.conj()])
        P = np.hstack([U_V, J_U_V])
        inverse_P = np.linalg.solve(P, np.eye(P.shape[0]))
        H_new = P.conj().T @ D_herm @ P
        return U, V, P, inverse_P, H_new

    def _resonator_mode_transformation_from_gprime(self, omega_r_new, g_prime):
        """Core of compute_resonator_mode_transformation, keyed directly on the
        already-eta-transformed coupling vector g_prime (length n).

        Returns (omega_beta, U2, V2, P2, b0) -- see
        compute_resonator_mode_transformation for the meaning of each.
        """
        n = self.n
        pref = self._nambu_prefactor()
        g_prime = np.asarray(g_prime, dtype=float)
        if g_prime.shape != (n,):
            raise ValueError(f"g_prime must have shape ({n},), got {g_prime.shape}")

        A = np.zeros((n + 1, n + 1))
        B = np.zeros((n + 1, n + 1))

        # chain block: exactly the already-diagonalized quadratic form
        A[:n, :n] = np.asarray(self.H_C[:n, :n]).real
        B[:n, :n] = np.asarray(self.H_C[:n, n:]).real

        # self.H_C's diagonal is (physical frequency / pref); everything added
        # for the new resonator is scaled the same way so the whole matrix
        # stays on one consistent scale, then omega_beta multiplies pref back.
        A[n, n] = omega_r_new / pref

        # (alpha^dag - alpha) sum_j g'_j (eta^dag_j - eta_j)
        #   = sum_j g'_j [ alpha^dag eta^dag_j - alpha^dag eta_j
        #                  - alpha eta^dag_j + alpha eta_j ]
        A[n, :n] = -g_prime / pref
        A[:n, n] = -g_prime / pref
        B[n, :n] = g_prime / pref
        B[:n, n] = g_prime / pref

        U2, V2, P2, _, H_beta = self._diagonalize_extended_quadratic(A, B)
        omega_beta = pref * np.diag(H_beta)[: n + 1].real
        b0 = int(np.argmin(np.abs(omega_beta - omega_r_new)))
        return omega_beta, U2, V2, P2, b0

    def compute_resonator_mode_transformation(self, omega_r_new, g0, ground_mode,
                                              mode, positions):
        """Diagonalize the quadratic part of chain + new resonator.

        Assembles the (n+1)-mode quadratic Hamiltonian

            sum_j omega_j eta^dag_j eta_j + omega_r_new alpha^dag alpha
            + (alpha^dag - alpha) sum_j g'_j (eta^dag_j - eta_j)

        (the chain block taken straight from self.H_C, the new resonator at
        index n) and Bogoliubov-diagonalizes it into modes beta_k. The
        standing-wave coupling parametrization (g0, ground_mode, mode,
        positions) goes through compute_coupling_vec; for an arbitrary per-site
        coupling call _resonator_mode_transformation_from_gprime directly with
        compute_coupling_vec_from_sites(g_sites).

        Returns
        -------
        omega_beta : (n+1,) ndarray
            Physical beta-mode frequencies, ascending.
        U2, V2 : (n+1, n+1) ndarray
            Bogoliubov coefficients eta+alpha -> beta, in the convention
            change_basis expects (b_i = U_ij b_j + V*_ij b^dag_j).
        P2 : (2(n+1), 2(n+1)) ndarray
            The full pseudo-unitary.
        b0 : int
            Index of the beta mode whose frequency is closest to omega_r_new.
        g_prime : (n,) ndarray
            Resonator-to-eta-mode couplings (compute_coupling_vec output).
        """
        g_prime = self.compute_coupling_vec(g0, ground_mode, mode, positions)
        omega_beta, U2, V2, P2, b0 = self._resonator_mode_transformation_from_gprime(
            omega_r_new, g_prime
        )
        return omega_beta, U2, V2, P2, b0, np.asarray(g_prime, dtype=float)

    def transform_kerr_to_resonator_basis(self, U2, V2):
        """Push the eta-basis Kerr tensor V_{ijkl} through the second
        Bogoliubov transform (eta+alpha -> beta) and return the normal-ordered
        (1,1,0,0) tensor V'_{ijkl}, shape (n+1,)*4. The new resonator carries
        no Kerr, so its slices of the input tensor are zero.
        """
        if self.kerr_operator is None:
            raise RuntimeError(
                "No Kerr tensor available (kerr was None/0, or "
                "statistics='fermion'). The resonator dispersive shift is "
                "driven entirely by the Kerr nonlinearity and is identically "
                "zero here."
            )
        n = self.n
        ext = self._extended_engine()

        V_ext = np.zeros((n + 1, n + 1, n + 1, n + 1), dtype=complex)
        V_ext[:n, :n, :n, :n] = self.kerr_operator.tensor

        op = QuantumOperator([TensorTerm(V_ext, (1, 1, 0, 0))])
        transformed = ext.normal_order(
            ext.change_basis(op, U2, V2).simplify()
        ).simplify()

        for term in transformed.terms:
            if term.daggers == (1, 1, 0, 0):
                return term.tensor
        return np.zeros((n + 1, n + 1, n + 1, n + 1), dtype=complex)

    @staticmethod
    def _resonator_number_kerr_matrix(V_prime, b0):
        """From sum V'_{pqrs} beta^dag_p beta^dag_q beta_r beta_s isolate the
        terms with exactly one creation and one annihilation on the resonator
        mode b0, i.e.

            (sum_{ij} V''_ij beta^dag_i beta_j) beta_0^dag beta_0 .

        For one index of {p,q} and one of {r,s} equal to b0 and the other two
        chain modes i, j (which commute with beta_0, beta_0^dag), normal
        ordering always yields beta_0^dag beta_0 beta^dag_i beta_j, so V''_ij
        sums the four placements. Terms with 0 or 2 resonator excitations, or
        that change the resonator number, are dropped (single-excitation
        subspace / dispersive regime).

        Returns (V_dprime, chain_idx): V_dprime shape (n, n), chain_idx the
        beta indices != b0 in ascending order.
        """
        m1 = V_prime.shape[0]
        chain_idx = [k for k in range(m1) if k != b0]
        m = len(chain_idx)
        V = np.zeros((m, m), dtype=complex)
        for a, i in enumerate(chain_idx):
            for b, j in enumerate(chain_idx):
                V[a, b] = (V_prime[b0, i, b0, j] + V_prime[b0, i, j, b0]
                           + V_prime[i, b0, b0, j] + V_prime[i, b0, j, b0])
        V = 0.5 * (V + V.conj().T)   # Hermitian by construction; kill noise
        return V, chain_idx

    @staticmethod
    def _diagonalize_degenerate_blocks(M, energies, tol):
        """Block-diagonalize the Hermitian M within each set of near-degenerate
        ``energies`` (connected under |E_i - E_j| <= tol). Returns (M_rot, W)
        with M_rot = W^dag M W and W a block-diagonal unitary, so M_rot has no
        off-diagonal element linking two equal-energy states.
        """
        m = len(energies)
        W = np.eye(m, dtype=complex)

        remaining = set(range(m))
        groups = []
        while remaining:
            group = [remaining.pop()]
            stack = list(group)
            while stack:
                a = stack.pop()
                for b in list(remaining):
                    if abs(energies[a] - energies[b]) <= tol:
                        remaining.discard(b)
                        group.append(b)
                        stack.append(b)
            groups.append(sorted(group))

        for group in groups:
            if len(group) == 1:
                continue
            sub = M[np.ix_(group, group)]
            _, vecs = la.eigh(0.5 * (sub + sub.conj().T))
            W[np.ix_(group, group)] = vecs

        M_rot = W.conj().T @ M @ W
        return M_rot, W

    def compute_resonator_dispersive_shift_from_gprime(self, omega_r_new, g_prime,
                                                       degeneracy_tol=1e-6,
                                                       dispersive_tol=0.1):
        """Dispersive shift imprinted on a new resonator by each single
        excitation of the chain, given the eta-basis coupling vector g_prime
        directly (length n).

        Diagonalize chain + resonator -> pick beta_0 (resonator-like) ->
        transform the Kerr tensor into the beta basis -> isolate the
        beta_0^dag beta_0 - weighted chain quadratic V'' -> diagonalize V'' on
        degenerate chain subspaces -> second-order shift

            chi_i = V''_ii + sum_{j != i} |V''_ij|^2 / (E_i - E_j).

        This is the shared core of compute_resonator_dispersive_shift_vec
        (standing-wave parametrization) and
        compute_resonator_dispersive_shift_vec_from_sites (arbitrary per-site
        coupling).

        Parameters
        ----------
        omega_r_new : float
            New-resonator frequency, in qubit-frequency units (like omega_j).
        g_prime : (n,) array_like
            Resonator-to-eta-mode couplings, e.g. from
            compute_coupling_vec / compute_coupling_vec_from_sites.
        degeneracy_tol : float
            Beta modes within this frequency window are treated as degenerate;
            V'' is diagonalized inside each such block and those pairs are
            skipped in the second-order sum.
        dispersive_tol : float
            Raises ValueError if any |g'_j| / |omega_j - omega_r_new| reaches
            this, matching compute_dispersive_shift_vec.

        Returns
        -------
        chi : (n,) ndarray
            Real dispersive shift per chain mode, ordered by ascending beta
            frequency (resonator removed) after the degenerate-subspace
            rotation.
        """
        n = self.n
        g_prime = np.asarray(g_prime, dtype=float)
        if g_prime.shape != (n,):
            raise ValueError(f"g_prime must have shape ({n},), got {g_prime.shape}")

        omega_chain = self._nambu_prefactor() * np.diag(self.H_C)[:n].real
        coupled = np.abs(g_prime) > 1e-14   # modes the resonator actually sees
        with np.errstate(divide='ignore', invalid='ignore'):
            ratio = np.where(coupled,
                             np.abs(g_prime) / np.abs(omega_chain - omega_r_new),
                             0.0)
        if coupled.any() and np.max(ratio) >= dispersive_tol:
            raise ValueError(
                "Not all resonator-chain couplings are below the dispersive "
                f"threshold (max |g'|/|Delta| = {np.max(ratio):.3g} "
                f">= {dispersive_tol})."
            )

        omega_beta, U2, V2, _P2, b0 = self._resonator_mode_transformation_from_gprime(
            omega_r_new, g_prime
        )
        V_prime = self.transform_kerr_to_resonator_basis(U2, V2)
        V_dprime, chain_idx = self._resonator_number_kerr_matrix(V_prime, b0)
        E = omega_beta[chain_idx].real

        V_rot, _ = self._diagonalize_degenerate_blocks(V_dprime, E, degeneracy_tol)

        chi = np.zeros(n)
        for i in range(n):
            shift = V_rot[i, i].real
            for j in range(n):
                if j == i:
                    continue
                denom = E[i] - E[j]
                if abs(denom) <= degeneracy_tol:
                    continue
                shift += np.abs(V_rot[i, j]) ** 2 / denom
            chi[i] = shift
        return chi

    def compute_resonator_dispersive_shift_vec(self, omega_r_new, g0, ground_mode,
                                               mode, positions,
                                               degeneracy_tol=1e-6,
                                               dispersive_tol=0.1):
        """Dispersive shift imprinted on a new resonator by each single
        excitation of the chain, standing-wave coupling parametrization.

        g_sites[j] = g0 * sin(2 pi (1+mode) x_j / (L * ground_mode)); the rest
        of the pipeline is compute_resonator_dispersive_shift_from_gprime.
        For an arbitrary per-site coupling use
        compute_resonator_dispersive_shift_vec_from_sites.

        Returns the length-n vector chi (see
        compute_resonator_dispersive_shift_from_gprime).
        """
        g_prime = self.compute_coupling_vec(g0, ground_mode, mode, positions)
        return self.compute_resonator_dispersive_shift_from_gprime(
            omega_r_new, g_prime, degeneracy_tol, dispersive_tol
        )

    def compute_resonator_dispersive_shift_vec_from_sites(self, omega_r_new, g_sites,
                                                          degeneracy_tol=1e-6,
                                                          dispersive_tol=0.1):
        """Dispersive shift imprinted on a new resonator by each single
        excitation of the chain, for an arbitrary per-site coupling vector.

        The resonator is coupled as
        (alpha^dag - alpha) sum_j g_sites[j] (b^dag_j - b_j) in the bare site
        basis; g_sites is transformed to the eta basis via
        compute_coupling_vec_from_sites and the rest of the pipeline is
        compute_resonator_dispersive_shift_from_gprime.

        Parameters
        ----------
        omega_r_new : float
            New-resonator frequency, in qubit-frequency units.
        g_sites : (2N,) array_like
            Per-site coupling coefficients (e.g. a constant vector for uniform
            coupling across all qubits).
        degeneracy_tol, dispersive_tol : float
            As in compute_resonator_dispersive_shift_from_gprime.

        Returns
        -------
        chi : (n,) ndarray
        """
        g_prime = self.compute_coupling_vec_from_sites(g_sites)
        return self.compute_resonator_dispersive_shift_from_gprime(
            omega_r_new, g_prime, degeneracy_tol, dispersive_tol
        )

    def probe_freq_in_gap(self, gap_index, fraction=0.5, margin=0.4):
        """A probe frequency placed relative to the chain's eta-mode spectrum.

        The chain has n = 2N eta modes with (sorted) frequencies omega. This
        returns one omega_r_new positioned by ``gap_index``:

          gap_index = -1          -> below every mode
          gap_index = k, 0..n-2   -> inside the gap between sorted modes k, k+1
          gap_index = n-1         -> above every mode

        Interior gaps land at omega[k] + fraction*(omega[k+1] - omega[k]). The
        below/above positions sit ``margin`` times the full mode bandwidth
        (omega[-1] - omega[0]) past the extreme mode -- a bandwidth fraction
        rather than a level-spacing count, so "just above" stays comfortably
        detuned even when the probe couples mostly to a band-edge mode (a
        uniform g_j couples dominantly to one edge, so "above all" would
        otherwise be near-resonant with it while "below all" is a full
        bandwidth away). Because the eta frequencies move with (Omega_C,
        Chi_C), calling this per parameter point keeps "the probe is in gap k"
        true across a whole sweep.
        """
        omega = np.sort(self._nambu_prefactor() * np.diag(self.H_C)[:self.n].real)
        bandwidth = float(omega[-1] - omega[0])
        if gap_index == -1:
            return float(omega[0] - margin * bandwidth)
        if gap_index == self.n - 1:
            return float(omega[-1] + margin * bandwidth)
        if not 0 <= gap_index <= self.n - 2:
            raise ValueError(
                f"gap_index must be -1..{self.n - 1}, got {gap_index}"
            )
        lo, hi = omega[gap_index], omega[gap_index + 1]
        return float(lo + fraction * (hi - lo))

    def mode_discriminability(self, omega_r_new, g_sites, normalize=1.0,
                              degeneracy_tol=1e-6, dispersive_tol=0.1):
        """How well a single dispersive measurement resolves the chain modes.

        Returns min_{i<j} |chi_i - chi_j| / normalize, the smallest pairwise
        gap between the per-mode dispersive shifts (the bottleneck for reading
        off which mode is excited). Pass normalize = g**2 for a coupling-
        independent figure of merit, since chi ~ O(g^2 * kerr).

        Propagates the dispersive-regime ValueError from the underlying shift
        calculation (callers sweeping a grid should catch it -> NaN).
        """
        chi = self.compute_resonator_dispersive_shift_vec_from_sites(
            omega_r_new, g_sites,
            degeneracy_tol=degeneracy_tol, dispersive_tol=dispersive_tol,
        )
        return min_pairwise_gap(chi) / normalize

    def min_abs_shift(self, omega_r_new, g_sites, normalize=1.0,
                      degeneracy_tol=1e-6, dispersive_tol=0.1):
        """The smallest dispersive-shift magnitude over the chain modes.

        Returns min_i |chi_i| / normalize -- the *weakest* imprint any single
        chain excitation leaves on the probe, i.e. a floor on what is
        resolvable at all. (mode_discriminability, by contrast, is about
        telling two excited modes apart: min_{i<j}|chi_i - chi_j|.) Pass
        normalize = g**2 for a coupling-independent figure of merit, since
        chi ~ O(g^2 * kerr).

        Propagates the dispersive-regime ValueError from the underlying shift
        calculation (callers sweeping a grid should catch it -> NaN).
        """
        chi = self.compute_resonator_dispersive_shift_vec_from_sites(
            omega_r_new, g_sites,
            degeneracy_tol=degeneracy_tol, dispersive_tol=dispersive_tol,
        )
        return float(np.min(np.abs(chi))) / normalize

    def compute_dispersive_shift_vec(self, Omega: float, g0 : float, ground_mode : float, mode : int, positions : np.ndarray, dispersive_tol = 0.1) -> np.ndarray:
        """Get dispersive shifts between a resonator mode and the single-particle chain eigenmodes.

        The shift of chain mode j is g'_j^2 / (omega_j - Omega), the standard
        second-order (dispersive) result, with omega_j the chain eigenmode
        frequency and Omega the resonator frequency (flip the overall sign for
        the resonator-minus-qubit convention).

        Args:
            Omega (float): resonator frequency in units of the qubit frequency
            g0 (float): initial coupling magnitude determined by capacitance and modes fluctuations
            ground_mode (float): ratio of first mode wavelength to resonator length lambda_0/L
            mode (int): mode order with respect to the ground mode
            positions (np.ndarray): normalized positions of the qubits along the resonator x_i / L
            dispersive_tol (float, optional): maximum coupling to energy_difference ratio allowed for dispersive shift calculation. Defaults to 0.1
        """
        energies = self._nambu_prefactor() * np.diag(self.H_C)[:self.n].real
        detuning = energies - Omega
        coupling_magnitudes = self.compute_coupling_vec(g0, ground_mode, mode, positions)
        # np.all(...) must wrap the elementwise comparison, not be applied to
        # the raw ratio array before comparing a bool to dispersive_tol --
        # all(array) tests truthiness of each element, not the threshold.
        if not np.all(np.abs(coupling_magnitudes) / np.abs(detuning) < dispersive_tol):
            raise ValueError("Not all coupling strengths are below the dispersive threshold")

        return coupling_magnitudes**2 / detuning
    
    def compute_common_bus_shift_matrix(self, Omega: float, positions : np.ndarray, g0 : float = 1.0, ground_mode : float = 2.0, dispersive_tol = 0.1) -> np.ndarray:
        """Given a common bus resonator with zero voltage boundary conditions at its ends and with many modes, computes the dispersive shift matrix
        between them and the chain eigenmodes.

        Bus mode m (m = 0 .. n-1) sits at frequency Omega * (m + 1) -- a
        harmonic ladder above the base frequency Omega -- and its capacitive
        coupling scales as sqrt(m + 1) (the same scaling used by
        compute_common_bus_coupling_matrix).

        Args:
            Omega (float): base (fundamental) frequency of the bus
            ground_mode (float): number of antinodes in the base mode
            positions (np.ndarray): normalized positions of the qubits along the resonator x_i / L
            dispersive_tol =0.1 (float, optional): maximum coupling to energy_difference ratio allowed. Defaults to 0.1.

        Returns:
            np.ndarray: dispersive shift matrix between chain and bus resonator modes
        """

        m = np.arange(self.n) #the idea is to use the same number of resonators as there is of bosonic modes in the chain
        omega_m = Omega * (m + 1)  # bus mode m: harmonic (m+1) of the base frequency
        gm = g0 * np.sqrt(m + 1)   # capacitive coupling ~ sqrt(omega_m)
        shift_matrix = np.array([
            self.compute_dispersive_shift_vec(om_, gm_, ground_mode, m_, positions, dispersive_tol)
            for om_, gm_, m_ in zip(omega_m, gm, m)
        ])
        return shift_matrix

    def compute_common_bus_coupling_matrix(self, positions : np.ndarray, g0 : float = 1.0, ground_mode : float = 2.0) -> np.ndarray:
        """Given a common bus resonator with zero voltage boundary conditions at its ends and with many modes, computes the coupling matrix
        between the bus modes and the chain eigenmodes.

        Bus mode m (m = 0 .. n-1) couples with capacitive strength
        g0 * sqrt(m + 1) (same scaling as compute_common_bus_shift_matrix).

        Args:
            ground_mode (float): number of antinodes in the base mode
            positions (np.ndarray): normalized positions of the qubits along the resonator x_i / L

        Returns:
            np.ndarray: coupling matrix between chain and bus resonator modes
        """

        m = np.arange(self.n) #the idea is to use the same number of resonators as there is of bosonic modes in the chain
        gm = g0 * np.sqrt(m + 1) # capacitive coupling scaling from each mode
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


# ===========================================================================
#  Momentum-space analysis of the periodic SSH photonic chain
#  (developed in notebooks/unitary/zak_phase/k_space_diag.ipynb)
#
#  Small numpy utilities live at module level; the physics-result methods
#  live on PeriodicPHOTONICChain.  OpenPHOTONICChain is the open-boundary
#  counterpart (no momentum-space methods -- they assume translational
#  invariance).
# ===========================================================================


def k_grid(num_cells, shift=True):
    """Canonical Brillouin-zone momenta on the FFT grid: k_m = 2*pi*fftfreq(num_cells).
    shift=True returns them fftshift-ed onto the symmetric interval (-pi, pi]."""
    momenta = 2 * np.pi * np.fft.fftfreq(num_cells)
    return np.fft.fftshift(momenta) if shift else momenta


def nambu_metric(half_dim):
    """tau_3 = diag(+1,...,+1, -1,...,-1)  (half_dim entries of each sign)."""
    return np.diag(np.concatenate([np.ones(half_dim), -np.ones(half_dim)]))


def real_space_nambu(num_cells, Omega_C, Chi_C, Omega_J, Chi_J, PBC=True):
    """Real-space Nambu blocks (cap_normal, cap_pairing, joseph_normal, joseph_pairing),
    byte-for-byte the construction in PHOTONICChain.__init__ (bosonic branch).
    num_sites = 2 * num_cells; site l = 2*cell + sublattice."""
    num_sites = 2 * num_cells
    g_cap_intra = 0.25 * Omega_C * (1 + Chi_C)
    g_cap_inter = 0.25 * Omega_C * (1 - Chi_C)
    g_joseph_intra =  Omega_J * (1 + Chi_J)
    g_joseph_inter = -Omega_J * (1 - Chi_J)

    cap_hopping = np.zeros((num_sites, num_sites))
    for site in range(num_sites - 1):
        cap_hopping[site, site + 1] = cap_hopping[site + 1, site] = (
            g_cap_intra if site % 2 == 0 else g_cap_inter)
    if PBC and num_cells != 1:
        cap_hopping[0, num_sites - 1] = cap_hopping[num_sites - 1, 0] = g_cap_inter
    cap_normal  = 0.5 * (np.eye(num_sites) + cap_hopping)
    cap_pairing = -0.5 * cap_hopping

    joseph_hopping = np.zeros((num_sites, num_sites))
    for site in range(num_sites - 1):
        joseph_hopping[site, site + 1] = joseph_hopping[site + 1, site] = (
            g_joseph_intra if site % 2 == 0 else g_joseph_inter)
    if PBC and num_cells != 1:
        joseph_hopping[0, num_sites - 1] = joseph_hopping[num_sites - 1, 0] = g_joseph_inter
    joseph_normal  =  0.5 * ((g_joseph_intra + g_joseph_inter) * np.eye(num_sites) - joseph_hopping)
    joseph_pairing = -0.5 * joseph_hopping - 0.5 * (g_joseph_intra + g_joseph_inter) * np.eye(num_sites)
    return cap_normal, cap_pairing, joseph_normal, joseph_pairing


def negated_k_index(momenta):
    """For each index m, the index m' with momenta[m'] == -momenta[m]  (mod 2*pi)."""
    num = len(momenta)
    partner = np.empty(num, dtype=int)
    for m in range(num):
        partner[m] = int(np.argmin(np.abs(np.angle(np.exp(1j * (momenta + momenta[m]))))))
    return partner


def fourier_operator_matrix(num_cells, shift=True):
    """Unitary fourier_matrix (num_sites x num_sites) with b_{j,s} =
    (1/sqrt(num_cells)) sum_m e^{i k_m j} c_{m,s}; block-diagonal in sublattice."""
    momenta = k_grid(num_cells, shift=shift)
    cell_indices = np.arange(num_cells)
    per_sublattice = np.exp(1j * np.outer(cell_indices, momenta)) / np.sqrt(num_cells)
    num_sites = 2 * num_cells
    fourier_matrix = np.zeros((num_sites, num_sites), dtype=complex)
    fourier_matrix[0::2, 0::2] = per_sublattice
    fourier_matrix[1::2, 1::2] = per_sublattice
    return fourier_matrix


def to_k_space_nambu(cap_normal, cap_pairing, joseph_normal, joseph_pairing, num_cells, shift=True):
    """Full (2*num_sites x 2*num_sites) Nambu matrices in k space:
    H_k = nambu_fourier^dag @ H @ nambu_fourier,  nambu_fourier = blkdiag(F, conj(F))."""
    fourier_matrix = fourier_operator_matrix(num_cells, shift=shift)
    num_sites = 2 * num_cells
    zeros = np.zeros((num_sites, num_sites))
    nambu_fourier = np.block([[fourier_matrix, zeros], [zeros, fourier_matrix.conj()]])
    H_cap    = np.block([[cap_normal, cap_pairing], [cap_pairing, cap_normal]])
    H_joseph = np.block([[joseph_normal, joseph_pairing], [joseph_pairing, joseph_normal]])
    return (nambu_fourier.conj().T @ H_cap    @ nambu_fourier,
            nambu_fourier.conj().T @ H_joseph @ nambu_fourier)


def sector_permutation(num_cells, shift=True):
    """Permutation of the 2*num_sites Nambu indices grouping them by unordered pair {k, -k}.
    Returns (permutation, sector_blocks) with sector_blocks = [(start, size, k_label), ...];
    applying it block-diagonalizes a k-space Nambu matrix (8x8 generic, 4x4 at a TRIM)."""
    momenta = k_grid(num_cells, shift=shift)
    partner = negated_k_index(momenta)
    num_sites = 2 * num_cells
    permutation, sector_blocks, placed = [], [], set()
    for m in range(num_cells):
        if m in placed:
            continue
        m_partner = int(partner[m])
        if m_partner == m:
            placed.add(m)
            indices = [2 * m, 2 * m + 1, num_sites + 2 * m, num_sites + 2 * m + 1]
        else:
            placed.update((m, m_partner))
            indices = [2 * m, 2 * m + 1, 2 * m_partner, 2 * m_partner + 1,
                       num_sites + 2 * m, num_sites + 2 * m + 1,
                       num_sites + 2 * m_partner, num_sites + 2 * m_partner + 1]
        sector_blocks.append((len(permutation), len(indices), float(momenta[m])))
        permutation.extend(indices)
    return np.array(permutation), sector_blocks


def off_block_norm(matrix, sector_blocks):
    """Frobenius norm of everything outside the declared diagonal blocks
    (0 <=> matrix is exactly block diagonal for that layout)."""
    outside = np.ones(matrix.shape, dtype=bool)
    for start, size, _ in sector_blocks:
        outside[start:start + size, start:start + size] = False
    return np.linalg.norm(matrix[outside])


def paraunitary_defect(transform, half_dim):
    """0  <=>  transform is paraunitary (transform^dag tau_3 transform = tau_3)."""
    tau3 = nambu_metric(half_dim)
    return np.abs(transform.conj().T @ tau3 @ transform - tau3).max()


def bogoliubov_frequencies(nambu_hamiltonian, half_dim):
    """physical mode frequencies = positive eigenvalues of tau_3 @ nambu_hamiltonian."""
    eigenvalues = np.linalg.eigvals(nambu_metric(half_dim) @ nambu_hamiltonian).real
    return np.sort(eigenvalues[eigenvalues > 1e-9])


def bosonic_annihilation_ops(fock_cutoffs):
    """annihilation operators for a list of bosonic modes with the given Fock dimensions."""
    operators = []
    for target_mode, target_cutoff in enumerate(fock_cutoffs):
        single_mode = np.diag(np.sqrt(np.arange(1, target_cutoff)), 1)
        full = np.array([[1.0]])
        for mode, cutoff in enumerate(fock_cutoffs):
            full = np.kron(full, single_mode if mode == target_mode else np.eye(cutoff))
        operators.append(full.astype(complex))
    return operators


def passive_nambu_transform(mode_mixing):
    """new = mode_mixing @ old  (mode_mixing unitary)  ->  Nambu transform Psi = transform Psi'."""
    inverse = np.linalg.inv(mode_mixing)
    zeros = np.zeros_like(inverse)
    return np.block([[inverse, zeros], [zeros, inverse.conj()]])


def phased_sublattice_rotation(bloch_phase, half_dim):
    """rows = (alpha, beta); alpha_k = (e^{-i phase/2} b_kA + e^{+i phase/2} b_kB)/sqrt(2).
    phase(-k) = -phase(k).  half_dim in {2 (TRIM), 4 (generic k)}."""
    at_plus_k = np.array([[np.exp(-1j * bloch_phase / 2),  np.exp(1j * bloch_phase / 2)],
                          [np.exp(-1j * bloch_phase / 2), -np.exp(1j * bloch_phase / 2)]]) / np.sqrt(2)
    if half_dim == 2:
        return at_plus_k
    at_minus_k = np.array([[np.exp(1j * bloch_phase / 2),  np.exp(-1j * bloch_phase / 2)],
                           [np.exp(1j * bloch_phase / 2), -np.exp(-1j * bloch_phase / 2)]]) / np.sqrt(2)
    zeros = np.zeros((2, 2))
    return np.block([[at_plus_k, zeros], [zeros, at_minus_k]])


def squeeze_nambu_transform(half_dim, squeeze_specs):
    """squeeze_specs = [(mode_i, mode_j, squeeze), ...];
    gamma_i -> cosh(squeeze) gamma_i + sinh(squeeze) gamma_j^dag.
    mode_i == mode_j  => single-mode squeeze (TRIM)."""
    diag = np.eye(half_dim, dtype=complex)
    anti_diag = np.zeros((half_dim, half_dim), dtype=complex)
    for mode_i, mode_j, squeeze in squeeze_specs:
        diag[mode_i, mode_i] = diag[mode_j, mode_j] = np.cosh(squeeze)
        anti_diag[mode_i, mode_j] = anti_diag[mode_j, mode_i] = np.sinh(squeeze)
    return np.block([[diag, anti_diag], [anti_diag.conj(), diag.conj()]])


def compute_zak_phase(bloch_spinor_at, num_k_points, sublattice_positions=(-0.25, 0.25)):
    """Discrete Berry / Zak phase of a band, Wilson loop in the orbital gauge:
    gamma_Zak = -Im ln [ prod_j <u_j|u_{j+1}> ],  u_k = e^{-i k r} psi_k.
    bloch_spinor_at(k) returns psi_k: a 2-vector (bare band) or a 4-vector Bogoliubov
    spinor (u_A, u_B, v_A, v_B); the hole block picks up the opposite orbital phase."""
    momenta = np.linspace(0.0, 2 * np.pi, num_k_points, endpoint=False)
    spinors = [bloch_spinor_at(k) for k in momenta]
    positions = np.asarray(sublattice_positions)
    bz_boundary_phase = (np.exp(-1j * 2 * np.pi * positions) if len(spinors[0]) == 2
                         else np.exp(-1j * 2 * np.pi * np.r_[positions, -positions]))
    wilson_product = np.prod([np.vdot(spinors[j], spinors[j + 1]) for j in range(num_k_points - 1)])
    wilson_product *= np.vdot(spinors[-1], bz_boundary_phase * spinors[0])
    return -np.angle(wilson_product)


class OpenPHOTONICChain(PHOTONICChain):
    """PHOTONICChain with open boundary conditions (PBC=False forced).

    Real-space edge-mode analysis would live here; the momentum-space methods
    are only defined on PeriodicPHOTONICChain, since they assume translational
    invariance.  The Schrieffer-Wolff cavity elimination is on OpenChainSW.
    """

    def __init__(self, N, Omega_C, Chi_C, omega_r, gamma, fock_photon,
                 Omega_J, Chi_J, kerr, coeff_tol=1e-9, zero_tol=1e-12, statistics='boson'):
        super().__init__(N=N, Omega_C=Omega_C, Chi_C=Chi_C, omega_r=omega_r, gamma=gamma,
                         fock_photon=fock_photon, Omega_J=Omega_J, Chi_J=Chi_J, kerr=kerr,
                         PBC=False, coeff_tol=coeff_tol, zero_tol=zero_tol, statistics=statistics)


class PeriodicPHOTONICChain(PHOTONICChain):
    """PHOTONICChain with periodic boundary conditions (PBC=True forced), plus the
    momentum-space analysis developed in k_space_diag.ipynb:

      * k_space_blocks / sector_nambu     -- per-{k,-k} sector Nambu Hamiltonians
      * diagonalize_cap_sector            -- phased sublattice rotation + two-mode squeeze
      * joseph_in_qp_gauge                -- H_J in the H_C-diagonalising alpha/beta gauge
      * qp_band_frequencies               -- the two bands  sqrt(1 +/- 2|delta_C_k|)
      * ssh_bloch_spinors / *_spinor      -- Bloch spinors for the Wilson loop
      * zak_phase                         -- bare-SSH or BdG-beta Zak phase
      * interband_resonant_momentum       -- k* of the cavity interband anticrossing

    The cavity Schrieffer-Wolff correction, the dressed effective Hamiltonian and
    the two-level reduction are on the subclass PeriodicChainSW.

    Every analysis method resolves its couplings from the instance
    (Omega_C, Chi_C, Omega_J, Chi_J, omega_r, gamma), but each accepts keyword
    overrides so a parameter sweep never rebuilds the chain.
    """

    def __init__(self, N, Omega_C, Chi_C, omega_r, gamma, fock_photon,
                 Omega_J, Chi_J, kerr, coeff_tol=1e-9, zero_tol=1e-12, statistics='boson'):
        super().__init__(N=N, Omega_C=Omega_C, Chi_C=Chi_C, omega_r=omega_r, gamma=gamma,
                         fock_photon=fock_photon, Omega_J=Omega_J, Chi_J=Chi_J, kerr=kerr,
                         PBC=True, coeff_tol=coeff_tol, zero_tol=zero_tol, statistics=statistics)

    # ------------------------------------------------------------------
    # parameter resolution
    # ------------------------------------------------------------------
    def _resolve_couplings(self, Omega_C=None, Chi_C=None, Omega_J=None, Chi_J=None):
        Omega_C = self.Omega_C if Omega_C is None else Omega_C
        Chi_C   = self.Chi_C   if Chi_C   is None else Chi_C
        Omega_J = self.Omega_J if Omega_J is None else Omega_J
        Chi_J   = self.Chi_J   if Chi_J   is None else Chi_J
        return (0.25 * Omega_C * (1 + Chi_C), 0.25 * Omega_C * (1 - Chi_C),
                Omega_J * (1 + Chi_J), -Omega_J * (1 - Chi_J))

    def momentum_grid(self, shift=True):
        """The discrete Brillouin-zone momenta (2*N of them, FFT grid, fftshift-ed)."""
        return k_grid(self.N, shift=shift)

    # ------------------------------------------------------------------
    # per-{k,-k} sector Nambu Hamiltonians
    # ------------------------------------------------------------------
    def k_space_blocks(self, k, **couplings):
        """Reduced per-sector k-space blocks (cap_normal, cap_pairing, joseph_normal,
        joseph_pairing) in the factor-2 / all-annihilation convention:
        4x4 for generic k, 2x2 at a TRIM (sin(k) == 0)."""
        g_cap_intra, g_cap_inter, g_joseph_intra, g_joseph_inter = self._resolve_couplings(**couplings)
        g_joseph_sum = g_joseph_intra + g_joseph_inter
        bloch_cap    = g_cap_intra    + g_cap_inter    * np.exp(-1.0j * k)
        bloch_joseph = g_joseph_intra + g_joseph_inter * np.exp(-1.0j * k)
        bloch_cap_conj, bloch_joseph_conj = np.conj(bloch_cap), np.conj(bloch_joseph)

        if np.isclose(np.sin(k), 0.0):
            cap_normal  = np.array([[1.0, bloch_cap], [bloch_cap_conj, 1.0]])
            cap_pairing = np.array([[0.0, -bloch_cap], [-bloch_cap_conj, 0.0]])
            joseph_normal  = np.array([[g_joseph_sum, -bloch_joseph],
                                       [-bloch_joseph_conj, g_joseph_sum]])
            joseph_pairing = np.array([[-g_joseph_sum, -bloch_joseph],
                                       [-bloch_joseph_conj, -g_joseph_sum]])
            return cap_normal, cap_pairing, joseph_normal, joseph_pairing

        cap_normal = np.array([
            [1.0,            bloch_cap,      0.0,            0.0],
            [bloch_cap_conj, 1.0,            0.0,            0.0],
            [0.0,            0.0,            1.0,            bloch_cap_conj],
            [0.0,            0.0,            bloch_cap,      1.0]])
        cap_pairing = np.array([
            [0.0,             0.0,             0.0,            -bloch_cap],
            [0.0,             0.0,            -bloch_cap_conj,  0.0],
            [0.0,            -bloch_cap_conj,  0.0,             0.0],
            [-bloch_cap,      0.0,             0.0,             0.0]])
        joseph_normal = np.array([
            [g_joseph_sum,       -bloch_joseph,       0.0,                 0.0],
            [-bloch_joseph_conj,  g_joseph_sum,       0.0,                 0.0],
            [0.0,                 0.0,                 g_joseph_sum,       -bloch_joseph_conj],
            [0.0,                 0.0,                -bloch_joseph,        g_joseph_sum]])
        joseph_pairing = np.array([
            [0.0,                 0.0,                -g_joseph_sum,       -bloch_joseph],
            [0.0,                 0.0,                -bloch_joseph_conj,  -g_joseph_sum],
            [-g_joseph_sum,      -bloch_joseph_conj,   0.0,                 0.0],
            [-bloch_joseph,      -g_joseph_sum,        0.0,                 0.0]])
        return cap_normal, cap_pairing, joseph_normal, joseph_pairing

    def sector_nambu(self, k, **couplings):
        """8x8 capacitive + Josephson Nambu blocks for the {k,-k} sector,
        basis (b_kA, b_kB, b_-kA, b_-kB, + daggers), honest 1/2 Psi^dag H Psi convention."""
        g_cap_intra, g_cap_inter, g_joseph_intra, g_joseph_inter = self._resolve_couplings(**couplings)
        g_joseph_sum = g_joseph_intra + g_joseph_inter
        bloch_cap    = g_cap_intra    + g_cap_inter    * np.exp(-1j * k)
        bloch_joseph = g_joseph_intra + g_joseph_inter * np.exp(-1j * k)
        zero_2x2 = np.zeros((2, 2))

        def sector_from_bloch(on_site, bloch_off_diagonal):
            at_plus_k  = np.array([[on_site, bloch_off_diagonal],
                                   [np.conj(bloch_off_diagonal), on_site]], dtype=complex)
            at_minus_k = np.array([[on_site, np.conj(bloch_off_diagonal)],
                                   [bloch_off_diagonal, on_site]], dtype=complex)
            return np.block([[at_plus_k, zero_2x2], [zero_2x2, at_minus_k]])

        cap_normal  = 0.5 * sector_from_bloch(1.0, bloch_cap)
        cap_pairing = np.zeros((4, 4), dtype=complex)
        cap_pairing[0:2, 2:4] = -0.5 * np.array([[0, bloch_cap], [np.conj(bloch_cap), 0]])
        cap_pairing[2:4, 0:2] = -0.5 * np.array([[0, np.conj(bloch_cap)], [bloch_cap, 0]])

        joseph_normal = 0.5 * sector_from_bloch(g_joseph_sum, -bloch_joseph)
        joseph_pairing = np.zeros((4, 4), dtype=complex)
        joseph_pairing[0:2, 2:4] = (-0.5 * np.array([[0, bloch_joseph], [np.conj(bloch_joseph), 0]])
                                    - 0.5 * g_joseph_sum * np.eye(2))
        joseph_pairing[2:4, 0:2] = (-0.5 * np.array([[0, np.conj(bloch_joseph)], [bloch_joseph, 0]])
                                    - 0.5 * g_joseph_sum * np.eye(2))
        return (np.block([[cap_normal, cap_pairing], [cap_pairing.conj(), cap_normal.conj()]]),
                np.block([[joseph_normal, joseph_pairing], [joseph_pairing.conj(), joseph_normal.conj()]]))

    def diagonalize_cap_sector(self, k, **couplings):
        """(transform, H_cap_diag) for the {k,-k} sector: a paraunitary transform
        (phased sublattice rotation then two-mode squeeze) that makes the capacitive
        Nambu block numerically diagonal."""
        g_cap_intra, g_cap_inter, _, _ = self._resolve_couplings(**couplings)
        H_cap, _ = self.sector_nambu(k, **couplings)
        half_dim = H_cap.shape[0] // 2
        bloch_cap = g_cap_intra + g_cap_inter * np.exp(-1j * k)
        bloch_magnitude, bloch_phase = abs(bloch_cap), np.angle(bloch_cap)
        squeeze_alpha = 0.25 * np.log(1 + 2 * bloch_magnitude)
        squeeze_beta  = 0.25 * np.log(1 - 2 * bloch_magnitude)

        rotation = passive_nambu_transform(phased_sublattice_rotation(bloch_phase, half_dim))
        if half_dim == 2:
            squeeze = squeeze_nambu_transform(2, [(0, 0, squeeze_alpha), (1, 1, squeeze_beta)])
        else:
            squeeze = squeeze_nambu_transform(4, [(0, 2, squeeze_alpha), (1, 3, squeeze_beta)])
        transform = rotation @ squeeze
        return transform, transform.conj().T @ H_cap @ transform

    def joseph_in_qp_gauge(self, k, **couplings):
        """2 * transform^dag H_joseph^(k) transform  (physical convention): the Josephson
        Nambu block in the H_cap-diagonalising alpha/beta gauge.  Number block A_J is
        momentum-block-diagonal; pairing block B_J is momentum-block-anti-diagonal."""
        _, H_joseph = self.sector_nambu(k, **couplings)
        transform, _ = self.diagonalize_cap_sector(k, **couplings)
        return 2.0 * (transform.conj().T @ H_joseph @ transform)

    def qp_band_frequencies(self, k, **couplings):
        """(omega_alpha, omega_beta) = sqrt(1 +/- 2|delta_C_k|), the two physical bands."""
        _, H_cap_diag = self.diagonalize_cap_sector(k, **couplings)
        return tuple(2.0 * np.diag(H_cap_diag[:2, :2]).real)

    # ------------------------------------------------------------------
    # Bloch spinors + Zak phase
    # ------------------------------------------------------------------
    def ssh_bloch_spinors(self, k, **couplings):
        """[u_alpha (+|delta|) | u_beta (-|delta|)] eigenvectors of [[0, delta],[delta*, 0]]."""
        g_cap_intra, g_cap_inter, _, _ = self._resolve_couplings(**couplings)
        bloch_cap = g_cap_intra + g_cap_inter * np.exp(-1j * k)
        _, eigenvectors = np.linalg.eigh(np.array([[0.0 + 0j, bloch_cap], [np.conj(bloch_cap), 0.0]]))
        return eigenvectors[:, [1, 0]]

    def ssh_lower_band_spinor(self, k, **couplings):
        """lower (beta) band Bloch spinor of the bare SSH 2x2."""
        return self.ssh_bloch_spinors(k, **couplings)[:, 1]

    def beta_bogoliubov_spinor(self, k, **couplings):
        """positive beta-band Bogoliubov spinor (u_A, u_B, v_A, v_B), normalised
        Phi^dag tau_3 Phi = +1."""
        g_cap_intra, g_cap_inter, _, _ = self._resolve_couplings(**couplings)
        bloch_cap = g_cap_intra + g_cap_inter * np.exp(-1j * k)
        normal_block = np.array([[1.0 + 0j, bloch_cap], [np.conj(bloch_cap), 1.0]])
        bdg = np.block([[normal_block, np.eye(2) - normal_block],
                        [np.eye(2) - normal_block, normal_block]])
        tau3 = nambu_metric(2)
        eigenvalues, eigenvectors = np.linalg.eig(tau3 @ bdg)
        beta_freq = np.sqrt(1 - 2 * abs(bloch_cap))
        spinor = eigenvectors[:, np.argmin(np.abs(eigenvalues - beta_freq))]
        return spinor / np.sqrt(abs(np.real(spinor.conj() @ tau3 @ spinor)))

    def zak_phase(self, paired=False, num_k_points=400, sublattice_positions=(-0.25, 0.25),
                  **couplings):
        """Zak phase of the lower band. paired=False: bare SSH lower band.
        paired=True: BdG beta-band Bogoliubov spinor (the pairing leaves it unchanged)."""
        spinor_at = ((lambda k: self.beta_bogoliubov_spinor(k, **couplings)) if paired
                     else (lambda k: self.ssh_lower_band_spinor(k, **couplings)))
        return compute_zak_phase(spinor_at, num_k_points, sublattice_positions)

    # ------------------------------------------------------------------
    # interband resonance
    # ------------------------------------------------------------------
    def interband_resonant_momentum(self, omega_r=None, **couplings):
        """k* where omega_alpha(k) - omega_beta(k) == omega_r  (the cavity-mediated
        interband anticrossing), plus |delta_C_{k*}|.  Returns (k_star, bloch_magnitude)."""
        g_cap_intra, g_cap_inter, _, _ = self._resolve_couplings(**couplings)
        omega_r = self.omega_r if omega_r is None else omega_r
        bloch_magnitude = 0.5 * omega_r * np.sqrt(1 - omega_r ** 2 / 4)
        k_star = np.arccos(np.clip(
            (bloch_magnitude ** 2 - g_cap_intra ** 2 - g_cap_inter ** 2)
            / (2 * g_cap_intra * g_cap_inter), -1, 1))
        return k_star, bloch_magnitude


class OpenChainSW(OpenPHOTONICChain):
    """Open-boundary PHOTONICChain with the 2nd-order Schrieffer-Wolff elimination of
    the cavity coupling.

    The elimination itself is the inherited ``PHOTONICChain.schrieffer_wolff_correction``
    (real space, all matter modes at once): call ``build_full_hamiltonian(sw=True)``
    then read ``self.H``/``self.H_SW``.  Unlike the periodic case there is no analytic
    alpha/beta diagonalisation gauge for open boundaries -- the SSH open-chain modes are
    not plane waves -- so the correction is expressed in whatever gauge the numerical
    Bogoliubov diagonalisation produces.
    """


class PeriodicChainSW(PeriodicPHOTONICChain):
    """Periodic PHOTONICChain with the 2nd-order Schrieffer-Wolff correction of the
    cavity coupling, done per {k,-k} sector in the analytic H_C-diagonalising
    (alpha_k / beta_k) gauge -- the same physics as
    ``PHOTONICChain.schrieffer_wolff_correction`` but in a fixed, meaningful gauge
    rather than the arbitrary one LAPACK imposes.  Adds:

      * sw_correction_qp        -- eliminate the resonator; k-preserving vs k-changing pieces
      * sw_correction_composite -- eliminate the chain pairing into (cavity Fock) x QP
      * dressed_hamiltonian     -- cavity-dressed H_eff per momentum (+/- the SW fold)
      * dressed_sector          -- same on the full {k,-k} sector (to check +k<->-k = 0)
      * dressed_state / dressed_zak_phase -- Wilson loop of a cavity-dressed band
      * tls_downfold / tls_three_gaps     -- two-level reduction of the {|alpha_k,0>,|beta_k,1>} doublet
    """

    # ------------------------------------------------------------------
    # Schrieffer-Wolff correction of the cavity coupling
    # ------------------------------------------------------------------
    def sw_correction_qp(self, k, coupling_mode="full", qp_fock_cutoff=3, photon_fock_cutoff=2,
                         gamma=None, omega_r=None, denom_floor=1e-9, **couplings):
        """2nd-order SW elimination of the resonator in the single-excitation subspace,
        {k,-k} sector.  Returns (qp_freqs, sw_number_correction, sw_pairing_correction):
        4x4 over (alpha_k, beta_k, alpha_-k, beta_-k); sw_number is the k-preserving
        (number-conserving) piece, sw_pairing the k-changing (induced-pairing) piece.
        coupling_mode='full' eliminates the whole coupling; 'pairing' the B_J part only."""
        gamma   = self.gamma   if gamma   is None else gamma
        omega_r = self.omega_r if omega_r is None else omega_r
        H_cap, H_joseph = self.sector_nambu(k, **couplings)
        transform, H_cap_diag = self.diagonalize_cap_sector(k, **couplings)
        H_joseph_qp = transform.conj().T @ H_joseph @ transform
        num_qp_modes = 4
        qp_freqs = 2.0 * np.diag(H_cap_diag[:num_qp_modes, :num_qp_modes]).real
        joseph_number_qp  = 2.0 * H_joseph_qp[:num_qp_modes, :num_qp_modes]
        joseph_pairing_qp = 2.0 * H_joseph_qp[:num_qp_modes, num_qp_modes:]

        fock_cutoffs = [photon_fock_cutoff] + [qp_fock_cutoff] * num_qp_modes
        photon_annihilation, *qp_annihilation = bosonic_annihilation_ops(fock_cutoffs)
        photon_number_op = photon_annihilation.conj().T @ photon_annihilation
        H_unperturbed = omega_r * photon_number_op + sum(
            qp_freqs[m] * (qp_annihilation[m].conj().T @ qp_annihilation[m])
            for m in range(num_qp_modes))

        matter_coupling = np.zeros_like(H_unperturbed)
        for mode_i in range(num_qp_modes):
            for mode_j in range(num_qp_modes):
                if coupling_mode == "full":
                    matter_coupling += joseph_number_qp[mode_i, mode_j] * (
                        qp_annihilation[mode_i].conj().T @ qp_annihilation[mode_j])
                matter_coupling += joseph_pairing_qp[mode_i, mode_j] * (
                    qp_annihilation[mode_i].conj().T @ qp_annihilation[mode_j].conj().T)
                matter_coupling += np.conj(joseph_pairing_qp[mode_i, mode_j]) * (
                    qp_annihilation[mode_j] @ qp_annihilation[mode_i])
        coupling_operator = -gamma * ((photon_annihilation + photon_annihilation.conj().T)
                                      @ matter_coupling)

        photon_occupation = np.round(np.real(np.diag(photon_number_op))).astype(int)
        in_P, in_Q = (photon_occupation == 0), (photon_occupation != 0)
        energies = np.real(np.diag(H_unperturbed))
        energies_P, energies_Q = energies[in_P], energies[in_Q]
        coupling_PQ = coupling_operator[np.ix_(in_P, in_Q)]
        denominator = energies_P[:, None] - energies_Q[None, :]
        generator = np.where(np.abs(denominator) > denom_floor,
                             coupling_PQ / np.where(np.abs(denominator) > denom_floor, denominator, 1.0),
                             0.0)
        H_sw = 0.5 * (generator @ coupling_PQ.conj().T + coupling_PQ @ generator.conj().T)

        P_states = [state for state, keep in zip(list(np.ndindex(*fock_cutoffs)), in_P) if keep]
        vacuum_index = next(i for i, s in enumerate(P_states) if sum(s[1:]) == 0)
        one_qp_index = {s[1:].index(1): i for i, s in enumerate(P_states) if sum(s[1:]) == 1}
        two_qp_index = {tuple(sorted(m for m in range(num_qp_modes) for _ in range(s[1 + m]))): i
                        for i, s in enumerate(P_states) if sum(s[1:]) == 2}
        sw_number_correction = np.array(
            [[H_sw[one_qp_index[i], one_qp_index[j]] for j in range(num_qp_modes)]
             for i in range(num_qp_modes)])
        sw_pairing_correction = np.array(
            [[H_sw[two_qp_index[tuple(sorted((i, j)))], vacuum_index]
              if tuple(sorted((i, j))) in two_qp_index else 0.0
              for j in range(num_qp_modes)] for i in range(num_qp_modes)])
        return qp_freqs, sw_number_correction, sw_pairing_correction

    def sw_correction_composite(self, k, photon_cutoff=6, qp_fock_cutoff=4, denom_floor=1e-6,
                                gamma=None, omega_r=None, **couplings):
        """2nd-order SW elimination of the chain pairing ~B_J (single-QP -> 3-QP virtual),
        photon ladder = gamma(a + a^dag), per {k,-k} sector.  Returns H_SW on
        (cavity Fock photon_cutoff) x (alpha_k, beta_k, alpha_-k, beta_-k), index n*4 + s.
        Denominators use the soft regulator  denom / (denom^2 + denom_floor^2):
        denom_floor=1e-6 reproduces PHOTONICChain.schrieffer_wolff_correction; a physical
        denom_floor (~5% of omega_r) keeps H_SW finite when a 3-QP virtual state crosses
        into resonance (where 2nd-order perturbation theory is invalid anyway)."""
        gamma   = self.gamma   if gamma   is None else gamma
        omega_r = self.omega_r if omega_r is None else omega_r
        _, H_joseph = self.sector_nambu(k, **couplings)
        transform, H_cap_diag = self.diagonalize_cap_sector(k, **couplings)
        H_joseph_qp = transform.conj().T @ H_joseph @ transform
        num_qp_modes = 4
        qp_freqs = 2.0 * np.diag(H_cap_diag[:num_qp_modes, :num_qp_modes]).real
        joseph_pairing_qp = 2.0 * H_joseph_qp[:num_qp_modes, num_qp_modes:]

        single_mode_lower = np.diag(np.sqrt(np.arange(1, qp_fock_cutoff)), 1)
        def qp_creation_op(target_mode):
            full = np.array([[1.0]])
            for mode in range(num_qp_modes):
                full = np.kron(full, single_mode_lower.conj().T if mode == target_mode
                               else np.eye(qp_fock_cutoff))
            return full
        qp_creation = [qp_creation_op(m) for m in range(num_qp_modes)]
        qp_fock_states = list(np.ndindex(*([qp_fock_cutoff] * num_qp_modes)))
        total_excitations = np.array([sum(s) for s in qp_fock_states])
        pairing_operator = sum(joseph_pairing_qp[i, j] * (qp_creation[i] @ qp_creation[j])
                               for i in range(num_qp_modes) for j in range(num_qp_modes))
        idx_single_qp = np.where(total_excitations == 1)[0]
        idx_virtual   = np.where(total_excitations == 3)[0]
        fock_energies = np.array(qp_fock_states) @ qp_freqs
        energy_single_qp, energy_virtual = fock_energies[idx_single_qp], fock_energies[idx_virtual]
        pairing_matrix_elements = pairing_operator[np.ix_(idx_virtual, idx_single_qp)]
        order = np.argsort([qp_fock_states[x].index(1) for x in idx_single_qp])
        energy_single_qp = energy_single_qp[order]
        pairing_matrix_elements = pairing_matrix_elements[:, order]
        dim_single_qp = 4

        photon_ladder = np.zeros((photon_cutoff, photon_cutoff))
        for n in range(photon_cutoff - 1):
            photon_ladder[n, n + 1] = photon_ladder[n + 1, n] = gamma * np.sqrt(n + 1)

        H_SW = np.zeros((photon_cutoff * dim_single_qp, photon_cutoff * dim_single_qp), dtype=complex)
        for n_bra in range(photon_cutoff):
            for n_ket in range(photon_cutoff):
                for n_mid in range(photon_cutoff):
                    ladder_factor = photon_ladder[n_bra, n_mid] * photon_ladder[n_mid, n_ket]
                    if abs(ladder_factor) < 1e-15:
                        continue
                    denom_first  = ((n_bra - n_mid) * omega_r
                                    + energy_single_qp[:, None] - energy_virtual[None, :])
                    denom_second = ((n_mid - n_ket) * omega_r
                                    + energy_virtual[:, None] - energy_single_qp[None, :])
                    resolvent_first  = denom_first  / (denom_first  ** 2 + denom_floor ** 2)
                    resolvent_second = denom_second / (denom_second ** 2 + denom_floor ** 2)
                    block = ladder_factor * (
                        (pairing_matrix_elements.T * resolvent_first) @ pairing_matrix_elements.conj()
                        - pairing_matrix_elements.T @ (pairing_matrix_elements.conj() * resolvent_second))
                    H_SW[n_bra * dim_single_qp:(n_bra + 1) * dim_single_qp,
                         n_ket * dim_single_qp:(n_ket + 1) * dim_single_qp] += block
        return 0.25 * (H_SW + H_SW.conj().T)

    # ------------------------------------------------------------------
    # cavity-dressed effective Hamiltonian
    # ------------------------------------------------------------------
    def _sector_freqs_and_coupling(self, k, num_qp_modes, **couplings):
        _, H_joseph = self.sector_nambu(k, **couplings)
        transform, H_cap_diag = self.diagonalize_cap_sector(k, **couplings)
        H_joseph_qp = transform.conj().T @ H_joseph @ transform
        qp_freqs = 2.0 * np.diag(H_cap_diag[:num_qp_modes, :num_qp_modes]).real
        joseph_number_qp = H_joseph_qp[:num_qp_modes, :num_qp_modes]
        return qp_freqs, joseph_number_qp + joseph_number_qp.conj().T

    def dressed_hamiltonian(self, k, photon_cutoff=6, include_sw=True, gamma=None, omega_r=None,
                            qp_fock_cutoff=4, denom_floor=1e-6, **couplings):
        """Cavity-dressed effective Hamiltonian on {alpha_k, beta_k} x Fock(photon_cutoff),
        index n*2 + s (s: 0=alpha, 1=beta).  With the chain pairing dropped
        (include_sw=False) or folded in via 2nd-order SW (include_sw=True).
        This is exactly the +k block -- +k and -k never couple (the coupling is
        momentum-block-diagonal and so is H_SW)."""
        gamma   = self.gamma   if gamma   is None else gamma
        omega_r = self.omega_r if omega_r is None else omega_r
        qp_freqs, matter_coupling = self._sector_freqs_and_coupling(k, 2, **couplings)
        photon_lower = np.diag(np.sqrt(np.arange(1, photon_cutoff)), 1)
        H_eff = np.diag(np.array([omega_r * n + qp_freqs[s]
                                  for n in range(photon_cutoff) for s in (0, 1)], dtype=complex))
        H_eff += -gamma * np.kron(photon_lower + photon_lower.conj().T, matter_coupling)
        if include_sw:
            plus_k = [n * 4 + s for n in range(photon_cutoff) for s in (0, 1)]
            H_eff += self.sw_correction_composite(k, photon_cutoff, qp_fock_cutoff, denom_floor,
                                                  gamma=gamma, omega_r=omega_r,
                                                  **couplings)[np.ix_(plus_k, plus_k)]
        return 0.5 * (H_eff + H_eff.conj().T)

    def dressed_sector(self, k, photon_cutoff=6, include_sw=True, gamma=None, omega_r=None,
                       qp_fock_cutoff=4, denom_floor=1e-6, **couplings):
        """Same as dressed_hamiltonian but on the full {alpha_k, beta_k, alpha_-k, beta_-k}
        x Fock(photon_cutoff) sector (index n*4 + s), so the (identically zero) +k <-> -k
        coupling can be measured."""
        gamma   = self.gamma   if gamma   is None else gamma
        omega_r = self.omega_r if omega_r is None else omega_r
        qp_freqs, matter_coupling = self._sector_freqs_and_coupling(k, 4, **couplings)
        photon_lower = np.diag(np.sqrt(np.arange(1, photon_cutoff)), 1)
        H_eff = np.diag(np.array([omega_r * n + qp_freqs[s]
                                  for n in range(photon_cutoff) for s in range(4)], dtype=complex))
        H_eff += -gamma * np.kron(photon_lower + photon_lower.conj().T, matter_coupling)
        if include_sw:
            H_eff += self.sw_correction_composite(k, photon_cutoff, qp_fock_cutoff, denom_floor,
                                                  gamma=gamma, omega_r=omega_r, **couplings)
        return 0.5 * (H_eff + H_eff.conj().T)

    def dressed_state(self, k, band, photon_cutoff=6, gamma=None, omega_r=None, denom_floor=1e-6,
                      **couplings):
        """Eigenstate(s) of dressed_hamiltonian, lifted into the fixed
        (sublattice x cavity Fock) basis via the SSH Bloch spinors -- ready for a
        Wilson loop.  `band` may be an int or an iterable of ints."""
        H_eff = self.dressed_hamiltonian(k, photon_cutoff, include_sw=True, gamma=gamma,
                                         omega_r=omega_r, denom_floor=denom_floor, **couplings)
        bloch_spinors = self.ssh_bloch_spinors(k, **couplings)
        _, eigenvectors = np.linalg.eigh(H_eff)
        bands = band if np.iterable(band) else [band]
        lifted = np.stack([(bloch_spinors @ eigenvectors[:, b].reshape(photon_cutoff, 2).T).reshape(-1)
                           for b in bands], axis=1)
        return lifted if np.iterable(band) else lifted[:, 0]

    def dressed_zak_phase(self, band, num_k_points=300, sublattice_positions=(-0.25, 0.25),
                          photon_cutoff=6, gamma=None, omega_r=None, denom_floor=1e-6, **couplings):
        """Zak phase of a cavity-dressed band (int) or a near-degenerate set of them
        (iterable -> det Wilson loop)."""
        momenta = np.linspace(0.0, 2 * np.pi, num_k_points, endpoint=False)
        loop = [self.dressed_state(k, band, photon_cutoff, gamma, omega_r, denom_floor, **couplings)
                for k in momenta]
        bz_boundary_phase = np.repeat(
            np.exp(-1j * 2 * np.pi * np.asarray(sublattice_positions)), photon_cutoff)
        if np.iterable(band):
            wilson = np.eye(len(band), dtype=complex)
            for j in range(num_k_points):
                nxt = ((bz_boundary_phase[:, None] * loop[0]) if j == num_k_points - 1
                       else loop[(j + 1) % num_k_points])
                wilson = wilson @ (loop[j].conj().T @ nxt)
            return -np.angle(np.linalg.det(wilson))
        wilson_product = (np.prod([np.vdot(loop[j], loop[j + 1]) for j in range(num_k_points - 1)])
                          * np.vdot(loop[-1], bz_boundary_phase * loop[0]))
        return -np.angle(wilson_product)

    # ------------------------------------------------------------------
    # two-level (doublet) reduction
    # ------------------------------------------------------------------
    @staticmethod
    def tls_downfold(hamiltonian, target_indices, denom_floor=1e-9):
        """2x2 effective doublet from a Hermitian block: diagonalise the target (P) and
        spectator (Q) sub-blocks, rotate the P-Q coupling into those eigenbases, apply
        the diagonal-denominator 2nd-order Schrieffer-Wolff shift, rotate back."""
        spectator_indices = [i for i in range(hamiltonian.shape[0]) if i not in target_indices]
        energies_P, basis_P = np.linalg.eigh(hamiltonian[np.ix_(target_indices, target_indices)])
        energies_Q, basis_Q = np.linalg.eigh(hamiltonian[np.ix_(spectator_indices, spectator_indices)])
        coupling = (basis_P.conj().T
                    @ hamiltonian[np.ix_(target_indices, spectator_indices)]
                    @ basis_Q)
        denom = energies_P[:, None] - energies_Q[None, :]
        resolvent = np.where(np.abs(denom) > denom_floor,
                             1.0 / np.where(np.abs(denom) > denom_floor, denom, 1.0), 0.0)
        second_order = (coupling * resolvent) @ coupling.conj().T
        return basis_P @ (np.diag(energies_P) + 0.5 * (second_order + second_order.conj().T)) @ basis_P.conj().T

    def tls_three_gaps(self, k, photon_cutoff=6, qp_fock_cutoff=4, denom_floor=0.05,
                       gamma=None, omega_r=None, **couplings):
        """The gap of the {|alpha_k,0>, |beta_k,1>} doublet, three ways:
        (exact) the two eigenvalues of the full dressed block with the most weight on
        the doublet; (bare) the raw 2x2 projection; (SW) the tls_downfold 2x2.
        Returns (gap_exact, gap_bare, gap_sw); gap_exact is nan where the doublet is
        not spectrally isolated."""
        H_eff = self.dressed_hamiltonian(k, photon_cutoff, include_sw=True, gamma=gamma,
                                         omega_r=omega_r, qp_fock_cutoff=qp_fock_cutoff,
                                         denom_floor=denom_floor, **couplings)
        triples = [(n, 1, 0) if s == 0 else (n, 0, 1)
                   for n in range(photon_cutoff) for s in (0, 1)]
        doublet = [triples.index((0, 1, 0)), triples.index((1, 0, 1))]
        eigenvalues, eigenvectors = np.linalg.eigh(H_eff)
        weight = (np.abs(eigenvectors[doublet, :]) ** 2).sum(0)
        top_two = np.argsort(weight)[-2:]
        gap_exact = np.ptp(eigenvalues[top_two]) if weight[top_two].min() > 0.7 else np.nan
        return (gap_exact,
                np.ptp(np.linalg.eigvalsh(H_eff[np.ix_(doublet, doublet)])),
                np.ptp(np.linalg.eigvalsh(self.tls_downfold(H_eff, doublet))))


# ---------------------------------------------------------------------------
# Dispersive-shift mode-discriminability sweeps
# ---------------------------------------------------------------------------

def min_pairwise_gap(values):
    """Smallest |v_i - v_j| over all pairs of a 1-D array.

    After sorting, the closest pair is always adjacent, so this is just
    min(diff(sort(values))). Returns +inf for fewer than two values.
    """
    s = np.sort(np.asarray(values, dtype=float).ravel())
    return float(np.min(np.diff(s))) if s.size > 1 else np.inf


def sweep_discriminability(make_chain, omega_C_values, chi_C_values, g_sites, *,
                           gap_index=None, omega_r_new=None, normalize=1.0,
                           degeneracy_tol=1e-6, dispersive_tol=0.3):
    """Mode-discriminability heat map over an (Omega_C, Chi_C) grid.

    For every grid point: build the chain, place the probe, and record
    min_{i<j}|chi_i - chi_j| / normalize (see PHOTONICChain.mode_discriminability).
    Grid points that fall outside the dispersive regime are left as NaN.

    Parameters
    ----------
    make_chain : callable
        ``make_chain(Omega_C, Chi_C) -> PHOTONICChain``. The caller binds N,
        kerr, PBC, statistics, ... so this function stays model-agnostic.
    omega_C_values, chi_C_values : 1-D arrays
        Grid axes (columns and rows of every returned 2-D array).
    g_sites : (2N,) array_like
        Per-site probe coupling (e.g. a constant vector for uniform coupling).
    gap_index : int or None
        If given, the probe frequency is adaptive per grid point,
        ``chain.probe_freq_in_gap(gap_index)`` -- so "the probe sits in gap k"
        stays true even as the eta spectrum moves. If None, ``omega_r_new`` is
        used as a fixed frequency everywhere.
    omega_r_new : float or None
        Fixed probe frequency; required when ``gap_index`` is None.
    normalize : float
        Divisor for the gap (pass ``g_const**2`` for a coupling-independent
        figure of merit).

    Returns
    -------
    dict with
        gap_grid          : (len(chi_C_values), len(omega_C_values))
        chi_cube          : (len(chi_C_values), len(omega_C_values), 2N)
        omega_r_new_grid  : (len(chi_C_values), len(omega_C_values))  -- the
                            probe frequency actually used at each point
        omega_C_values, chi_C_values : the axes, as arrays
    """
    if gap_index is None and omega_r_new is None:
        raise ValueError("pass either gap_index (adaptive) or omega_r_new (fixed)")

    omega_C_values = np.asarray(omega_C_values, dtype=float)
    chi_C_values = np.asarray(chi_C_values, dtype=float)
    shape = (len(chi_C_values), len(omega_C_values))

    gap_grid = np.full(shape, np.nan)
    wr_grid = np.full(shape, np.nan)
    chi_cube = None

    for a, cc in enumerate(chi_C_values):
        for b, oc in enumerate(omega_C_values):
            chain = make_chain(oc, cc)
            if chi_cube is None:
                chi_cube = np.full(shape + (chain.n,), np.nan)
            wr = (chain.probe_freq_in_gap(gap_index)
                  if gap_index is not None else omega_r_new)
            wr_grid[a, b] = wr
            try:
                chi = chain.compute_resonator_dispersive_shift_vec_from_sites(
                    wr, g_sites,
                    degeneracy_tol=degeneracy_tol, dispersive_tol=dispersive_tol,
                )
            except ValueError:
                continue
            chi_cube[a, b] = chi
            gap_grid[a, b] = min_pairwise_gap(chi) / normalize

    return dict(gap_grid=gap_grid, chi_cube=chi_cube, omega_r_new_grid=wr_grid,
                omega_C_values=omega_C_values, chi_C_values=chi_C_values)
