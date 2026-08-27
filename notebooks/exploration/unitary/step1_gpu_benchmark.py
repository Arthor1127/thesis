"""
Step 1 benchmark (per the GPU handoff plan) for
SecondQuantizationEngine._build_squeezed_vacuum_vectorized_boson.

Run this on a machine with an actual GPU + cupy installed. It does NOT
modify sq_engine.py or ssh_chain_sq.py -- it is a standalone diagnostic
that imports the real engine purely to extract realistic
(current_states, current_coeffs, nz_i, nz_j, nz_val) inputs at each
truncation level, matching what the function actually processes for a
real PHOTONICChain at N=8 (16 modes), the size used in the original
profiling table.

Usage:
    pip install cupy-cudaXXX   # matching your CUDA version
    python step1_gpu_benchmark.py

Decision rule (from the handoff): if GPU (including transfer) isn't at
least ~3-5x faster than CPU at these real sizes, STOP -- do not proceed
to Step 2. This script prints a clear verdict at the end.
"""

import time
import numpy as np
import sys

try:
    import cupy as cp
    CUPY_AVAILABLE = True
except ImportError:
    CUPY_AVAILABLE = False

sys.path.insert(0, '.')
from sq_engine import SecondQuantizationEngine
from ssh_chain_sq import PHOTONICChain


# --------------------------------------------------------------------- #
# Step A: extract real (M, nz_i, nz_j, nz_val) from an actual chain
# --------------------------------------------------------------------- #

def get_real_M_nonzeros(N=8, coeff_tol=1e-9):
    """
    Reproduces the exact M-matrix construction from
    _build_squeezed_vacuum_vectorized_boson, using a real chain's U, V.
    """
    chain = PHOTONICChain(
        N=N, Omega_C=0.3, Chi_C=0.5, omega_r=0.13, gamma=0.02,
        fock_photon=4, Omega_J=1.0, Chi_J=0.0, kerr=None,
        PBC=False, statistics='boson'
    )
    U, V = chain.U, chain.V
    n = chain.n

    U_V = np.vstack([U, V])
    tau_x = np.block([[np.zeros((n, n)), np.eye(n)], [np.eye(n), np.zeros((n, n))]])
    J_U_V = tau_x @ U_V.conj()
    P = np.hstack([U_V, J_U_V])
    P_inv = np.linalg.inv(P)
    U_tilde = P_inv[:n, :n]
    V_tilde = P_inv[:n, n:]
    M = np.linalg.solve(U_tilde, V_tilde)

    nz_i, nz_j = np.where(np.abs(M) > coeff_tol)
    nz_val = M[nz_i, nz_j]

    print(f"N={N}  n_modes={n}  |M| nonzeros={len(nz_val)} / {n*n} "
          f"({100*len(nz_val)/(n*n):.1f}% dense)")
    return n, nz_i, nz_j, nz_val


# --------------------------------------------------------------------- #
# Step B: run the real per-level recursion (CPU) up to truncation, and
# snapshot (current_states, current_coeffs) entering each level -- these
# are the actual realistic inputs to the broadcast step at each k.
# --------------------------------------------------------------------- #

def generate_level_snapshots(n, nz_i, nz_j, nz_val, max_truncation, coeff_tol=1e-9):
    """
    Re-runs the CPU recursion from _build_squeezed_vacuum_vectorized_boson
    (the dedup/pack logic only, in numpy) just far enough to capture the
    (current_states, current_coeffs) array ENTERING each truncation level
    -- i.e. the actual sizes/values the broadcast step would process in
    the real function, at the real N=8 problem size.
    """
    snapshots = {}
    current_states = np.zeros((1, n), dtype=np.int64)
    current_coeffs = np.array([1.0 + 0j])

    for k in range(1, max_truncation + 1):
        snapshots[k] = (current_states.copy(), current_coeffs.copy())

        num_states = current_states.shape[0]
        num_pairs = len(nz_val)
        if num_states == 0 or num_pairs == 0:
            current_states = np.zeros((0, n), dtype=np.int64)
            current_coeffs = np.zeros((0,), dtype=complex)
            continue

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
        base = max_occ + 2
        powers = base ** np.arange(n - 1, -1, -1, dtype=np.int64)
        packed = flat_states @ powers
        unique_packed, inverse = np.unique(packed, return_inverse=True)
        inverse = inverse.reshape(-1)
        summed_coeffs = np.zeros(unique_packed.shape[0], dtype=complex)
        np.add.at(summed_coeffs, inverse, flat_coeffs)
        keep_mask = np.abs(summed_coeffs) > coeff_tol
        kept_packed = unique_packed[keep_mask]
        next_coeffs = summed_coeffs[keep_mask]

        next_states = np.zeros((kept_packed.shape[0], n), dtype=np.int64)
        remaining = kept_packed.copy()
        for mode in reversed(range(n)):
            next_states[:, mode] = remaining % base
            remaining //= base

        current_states, current_coeffs = next_states, next_coeffs

    return snapshots


# --------------------------------------------------------------------- #
# CPU implementations of the three sub-steps (mirrors the production code)
# --------------------------------------------------------------------- #

def cpu_broadcast_step(current_states, current_coeffs, nz_i, nz_j, nz_val, k, n):
    num_states = current_states.shape[0]
    num_pairs = len(nz_val)
    new_states = np.broadcast_to(current_states, (num_pairs, num_states, n)).copy()
    rows = np.arange(num_pairs)
    new_states[rows, :, nz_j] += 1
    factor_j = np.sqrt(new_states[rows, :, nz_j].astype(float))
    new_states[rows, :, nz_i] += 1
    factor_i = np.sqrt(new_states[rows, :, nz_i].astype(float))
    new_coeffs = current_coeffs[None, :] * (-0.5) * nz_val[:, None] * factor_j * factor_i / float(k)
    return new_states, new_coeffs


def cpu_pack_step(flat_states, base, n):
    powers = base ** np.arange(n - 1, -1, -1, dtype=np.int64)
    return flat_states @ powers


def cpu_dedup_step(packed, flat_coeffs):
    unique_packed, inverse = np.unique(packed, return_inverse=True)
    inverse = inverse.reshape(-1)
    summed_coeffs = np.zeros(unique_packed.shape[0], dtype=complex)
    np.add.at(summed_coeffs, inverse, flat_coeffs)
    return unique_packed, summed_coeffs


def cpu_fused_level_step(current_states, current_coeffs, nz_i, nz_j, nz_val,
                          k, n, coeff_tol=1e-9):
    """CPU equivalent of gpu_fused_level_step (defined below, only if cupy
    is available), for a fair comparison at the LEVEL granularity (not the
    isolated-sub-step granularity used by the per-step benchmarks above).
    This is pure CPU/numpy logic and must always be importable regardless
    of whether cupy is installed."""
    num_states = current_states.shape[0]
    num_pairs = len(nz_val)

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
    base = max_occ + 2
    powers = base ** np.arange(n - 1, -1, -1, dtype=np.int64)
    packed = flat_states @ powers

    unique_packed, inverse = np.unique(packed, return_inverse=True)
    inverse = inverse.reshape(-1)
    summed_coeffs = np.zeros(unique_packed.shape[0], dtype=complex)
    np.add.at(summed_coeffs, inverse, flat_coeffs)

    keep_mask = np.abs(summed_coeffs) > coeff_tol
    kept_packed = unique_packed[keep_mask]
    next_coeffs = summed_coeffs[keep_mask]

    next_states = np.zeros((kept_packed.shape[0], n), dtype=np.int64)
    remaining = kept_packed.copy()
    for mode in reversed(range(n)):
        next_states[:, mode] = remaining % base
        remaining //= base

    return next_states, next_coeffs


# --------------------------------------------------------------------- #
# GPU implementations (only defined if cupy is available)
# --------------------------------------------------------------------- #

if CUPY_AVAILABLE:
    def gpu_broadcast_step(current_states, current_coeffs, nz_i, nz_j, nz_val, k, n):
        # Transfer in
        cs = cp.asarray(current_states)
        cc = cp.asarray(current_coeffs)
        ci = cp.asarray(nz_i)
        cj = cp.asarray(nz_j)
        cv = cp.asarray(nz_val)

        num_states = cs.shape[0]
        num_pairs = len(cv)
        new_states = cp.broadcast_to(cs, (num_pairs, num_states, n)).copy()
        rows = cp.arange(num_pairs)
        new_states[rows, :, cj] += 1
        factor_j = cp.sqrt(new_states[rows, :, cj].astype(float))
        new_states[rows, :, ci] += 1
        factor_i = cp.sqrt(new_states[rows, :, ci].astype(float))
        new_coeffs = cc[None, :] * (-0.5) * cv[:, None] * factor_j * factor_i / float(k)

        # Transfer out
        return cp.asnumpy(new_states), cp.asnumpy(new_coeffs)

    def gpu_pack_step(flat_states, base, n):
        fs = cp.asarray(flat_states)
        powers = base ** cp.arange(n - 1, -1, -1, dtype=cp.int64)
        packed = fs @ powers
        return cp.asnumpy(packed)

    def gpu_dedup_step(packed, flat_coeffs):
        p = cp.asarray(packed)
        fc = cp.asarray(flat_coeffs)
        unique_packed, inverse = cp.unique(p, return_inverse=True)
        inverse = inverse.reshape(-1)

        # cp.add.at does not support complex128 (confirmed: raises
        # "cupy.add.at only supports int32, int64, float16, float32,
        # float64, uint32, uint64 as data type"). cupyx.scatter_add's
        # documented supported-dtype list also excludes float64/complex
        # (only float32/int32/uint32/uint64/ulonglong via CUDA atomicAdd),
        # so it is not a usable fallback here either. Split into real and
        # imaginary float64 parts (both ARE supported), accumulate each
        # separately, and recombine -- mathematically identical to direct
        # complex accumulation since addition is linear and real/imaginary
        # parts never mix under +. Verified against a numpy reference
        # before relying on this.
        real_part = cp.zeros(unique_packed.shape[0], dtype=cp.float64)
        imag_part = cp.zeros(unique_packed.shape[0], dtype=cp.float64)
        cp.add.at(real_part, inverse, fc.real)
        cp.add.at(imag_part, inverse, fc.imag)
        summed_coeffs = real_part + 1j * imag_part

        return cp.asnumpy(unique_packed), cp.asnumpy(summed_coeffs)

    def gpu_fused_level_step(current_states, current_coeffs, nz_i, nz_j, nz_val,
                              k, n, coeff_tol=1e-9):
        """
        Single device round-trip covering broadcast + pack + dedup for one
        truncation level -- the actual shape Step 2's implementation would
        take (per the handoff: move data on-device once before the k-loop,
        only come back to host once after). The per-sub-step benchmark
        above pays transfer overhead three separate times per level, which
        is NOT representative of this -- it isolates each operation to
        measure it individually, at the cost of realism. This function is
        the realistic comparison point.

        Returns (next_states, next_coeffs) as numpy arrays, the same
        return contract as one iteration of the production CPU loop body.
        """
        cs = cp.asarray(current_states)
        cc = cp.asarray(current_coeffs)
        ci = cp.asarray(nz_i)
        cj = cp.asarray(nz_j)
        cv = cp.asarray(nz_val)

        num_states = cs.shape[0]
        num_pairs = len(cv)

        new_states = cp.broadcast_to(cs, (num_pairs, num_states, n)).copy()
        rows = cp.arange(num_pairs)
        new_states[rows, :, cj] += 1
        factor_j = cp.sqrt(new_states[rows, :, cj].astype(float))
        new_states[rows, :, ci] += 1
        factor_i = cp.sqrt(new_states[rows, :, ci].astype(float))
        new_coeffs = cc[None, :] * (-0.5) * cv[:, None] * factor_j * factor_i / float(k)

        flat_states = new_states.reshape(num_pairs * num_states, n)
        flat_coeffs = new_coeffs.reshape(num_pairs * num_states)

        max_occ = int(flat_states.max()) if flat_states.size > 0 else 0
        base = max_occ + 2
        powers = base ** cp.arange(n - 1, -1, -1, dtype=cp.int64)
        packed = flat_states @ powers

        unique_packed, inverse = cp.unique(packed, return_inverse=True)
        inverse = inverse.reshape(-1)
        real_part = cp.zeros(unique_packed.shape[0], dtype=cp.float64)
        imag_part = cp.zeros(unique_packed.shape[0], dtype=cp.float64)
        cp.add.at(real_part, inverse, flat_coeffs.real)
        cp.add.at(imag_part, inverse, flat_coeffs.imag)
        summed_coeffs = real_part + 1j * imag_part

        keep_mask = cp.abs(summed_coeffs) > coeff_tol
        kept_packed = unique_packed[keep_mask]
        next_coeffs = summed_coeffs[keep_mask]

        next_states = cp.zeros((kept_packed.shape[0], n), dtype=cp.int64)
        remaining = kept_packed.copy()
        for mode in reversed(range(n)):
            next_states[:, mode] = remaining % base
            remaining //= base

        # Single transfer out for the whole level
        return cp.asnumpy(next_states), cp.asnumpy(next_coeffs)


# --------------------------------------------------------------------- #
# Timing harness
# --------------------------------------------------------------------- #

def time_it(fn, *args, n_repeats=5):
    # warm-up (especially important for GPU: first call pays kernel
    # compilation / context init cost that subsequent calls don't)
    fn(*args)
    times = []
    for _ in range(n_repeats):
        t0 = time.perf_counter()
        result = fn(*args)
        times.append(time.perf_counter() - t0)
    return result, min(times), np.mean(times)


def run_benchmark():
    print("=" * 70)
    print("Step 1 GPU benchmark: _build_squeezed_vacuum_vectorized_boson")
    print("=" * 70)

    if not CUPY_AVAILABLE:
        print("\ncupy is NOT installed in this environment.")
        print("Install with: pip install cupy-cudaXXX (matching your CUDA version)")
        print("Cannot proceed with GPU timing. Exiting.")
        return

    print(f"\ncupy version: {cp.__version__}")
    try:
        device = cp.cuda.Device()
        print(f"GPU device: {device.id}  "
              f"compute capability: {device.compute_capability}")
    except Exception as e:
        print(f"WARNING: could not query GPU device info: {e}")

    # verify cp.unique(return_inverse=True) and cp.add.at exist before
    # relying on them anywhere below -- per the handoff's explicit caveat
    print("\n--- API availability checks ---")
    try:
        _, inv = cp.unique(cp.array([3, 1, 2, 1]), return_inverse=True)
        print(f"  cp.unique(return_inverse=True): OK  (sample inverse={cp.asnumpy(inv)})")
    except Exception as e:
        print(f"  cp.unique(return_inverse=True): FAILED -- {e}")
        print("  Cannot proceed: dedup step has no verified GPU path.")
        return

    try:
        # Test the ACTUAL approach gpu_dedup_step uses: split complex into
        # real/imaginary float64 parts, accumulate each separately, recombine.
        # A naive direct cp.add.at on a complex128 array is confirmed to fail
        # ("cupy.add.at only supports int32, int64, float16, float32, float64,
        # uint32, uint64 as data type") -- complex128 is NOT in that list.
        # cupyx.scatter_add was considered as an alternative but its own
        # documented supported-dtype list (float32/int32/uint32/uint64/
        # ulonglong via CUDA atomicAdd) also excludes float64/complex, so it
        # is not a usable fallback either -- ruled out by reading its source,
        # not assumed.
        test_real = cp.zeros(3, dtype=cp.float64)
        test_imag = cp.zeros(3, dtype=cp.float64)
        test_vals = cp.array([1.0+2.0j, 2.0+1.0j, 3.0-1.0j])
        test_idx  = cp.array([0, 1, 0])
        cp.add.at(test_real, test_idx, test_vals.real)
        cp.add.at(test_imag, test_idx, test_vals.imag)
        test_arr = test_real + 1j * test_imag
        expected = np.array([4.0+1.0j, 2.0+1.0j, 0.0+0.0j])
        ok = np.allclose(cp.asnumpy(test_arr), expected)
        print(f"  cp.add.at (float64, split real/imag for complex): "
              f"{'OK' if ok else 'WRONG RESULT'}  "
              f"(got {cp.asnumpy(test_arr)}, expected {expected})")
        if not ok:
            print("  Cannot proceed: split-complex add.at does not match expected semantics.")
            return
    except Exception as e:
        print(f"  cp.add.at (float64, split real/imag): FAILED -- {e}")
        print("  Cannot proceed: scatter-accumulate has no verified GPU path "
              "even via the float64 split workaround.")
        return

    # ------------------------------------------------------------------ #
    # Extract real inputs from an actual chain at N=8 (handoff's target)
    # ------------------------------------------------------------------ #
    print("\n--- Extracting real M-matrix from PHOTONICChain(N=8) ---")
    n, nz_i, nz_j, nz_val = get_real_M_nonzeros(N=8)

    max_truncation = 4
    print(f"\n--- Generating real level-{1}..{max_truncation} snapshots ---")
    snapshots = generate_level_snapshots(n, nz_i, nz_j, nz_val, max_truncation)
    for k, (states, coeffs) in snapshots.items():
        print(f"  level {k}: current_states.shape={states.shape}  "
              f"current_coeffs.shape={coeffs.shape}")

    # ------------------------------------------------------------------ #
    # Benchmark each sub-step at k=3 and k=4 (the expensive levels)
    # ------------------------------------------------------------------ #
    results = {}

    for k in [3, 4]:
        if k not in snapshots:
            print(f"\nSkipping k={k}: truncation only generated up to "
                  f"{max_truncation}")
            continue

        current_states, current_coeffs = snapshots[k]
        if current_states.shape[0] == 0:
            print(f"\nSkipping k={k}: empty state set at this level "
                  f"(truncation may have terminated early)")
            continue

        print(f"\n{'='*70}")
        print(f"k={k}  num_states={current_states.shape[0]}  "
              f"num_pairs={len(nz_val)}  "
              f"broadcast_shape=({len(nz_val)}, {current_states.shape[0]}, {n})")
        print(f"{'='*70}")

        # --- sub-step 1: broadcast + increment + sqrt ---
        print("\n[1] Broadcast + increment + sqrt step")
        (cpu_new_states, cpu_new_coeffs), cpu_min, cpu_mean = time_it(
            cpu_broadcast_step, current_states, current_coeffs,
            nz_i, nz_j, nz_val, k, n
        )
        print(f"    CPU: min={cpu_min*1000:.2f}ms  mean={cpu_mean*1000:.2f}ms")

        (gpu_new_states, gpu_new_coeffs), gpu_min, gpu_mean = time_it(
            gpu_broadcast_step, current_states, current_coeffs,
            nz_i, nz_j, nz_val, k, n
        )
        print(f"    GPU: min={gpu_min*1000:.2f}ms  mean={gpu_mean*1000:.2f}ms  "
              f"(includes host<->device transfer)")

        correctness = (np.allclose(cpu_new_states, gpu_new_states) and
                       np.allclose(cpu_new_coeffs, gpu_new_coeffs, atol=1e-8))
        speedup = cpu_min / gpu_min if gpu_min > 0 else float('inf')
        print(f"    Correctness: {'OK' if correctness else 'MISMATCH -- DO NOT TRUST THIS PATH'}")
        print(f"    Speedup (CPU_min / GPU_min): {speedup:.2f}x")
        results[(k, 'broadcast')] = (speedup, correctness)

        # --- sub-step 2: pack ---
        flat_states = cpu_new_states.reshape(-1, n)
        flat_coeffs = cpu_new_coeffs.reshape(-1)
        max_occ = int(flat_states.max()) if flat_states.size > 0 else 0
        base = max_occ + 2

        print("\n[2] Integer packing step (flat_states @ powers)")
        cpu_packed, cpu_min, cpu_mean = time_it(cpu_pack_step, flat_states, base, n)
        print(f"    CPU: min={cpu_min*1000:.2f}ms  mean={cpu_mean*1000:.2f}ms")

        gpu_packed, gpu_min, gpu_mean = time_it(gpu_pack_step, flat_states, base, n)
        print(f"    GPU: min={gpu_min*1000:.2f}ms  mean={gpu_mean*1000:.2f}ms")

        correctness = np.array_equal(cpu_packed, gpu_packed)
        speedup = cpu_min / gpu_min if gpu_min > 0 else float('inf')
        print(f"    Correctness: {'OK' if correctness else 'MISMATCH -- DO NOT TRUST THIS PATH'}")
        print(f"    Speedup: {speedup:.2f}x")
        results[(k, 'pack')] = (speedup, correctness)

        # --- sub-step 3: unique/dedup ---
        print("\n[3] np.unique + add.at dedup step")
        (cpu_unique, cpu_summed), cpu_min, cpu_mean = time_it(
            cpu_dedup_step, cpu_packed, flat_coeffs
        )
        print(f"    CPU: min={cpu_min*1000:.2f}ms  mean={cpu_mean*1000:.2f}ms")

        (gpu_unique, gpu_summed), gpu_min, gpu_mean = time_it(
            gpu_dedup_step, cpu_packed, flat_coeffs
        )
        print(f"    GPU: min={gpu_min*1000:.2f}ms  mean={gpu_mean*1000:.2f}ms")

        correctness = (np.array_equal(cpu_unique, gpu_unique) and
                       np.allclose(cpu_summed, gpu_summed, atol=1e-8))
        speedup = cpu_min / gpu_min if gpu_min > 0 else float('inf')
        print(f"    Correctness: {'OK' if correctness else 'MISMATCH -- DO NOT TRUST THIS PATH'}")
        print(f"    Speedup: {speedup:.2f}x")
        results[(k, 'dedup')] = (speedup, correctness)

        # --- FUSED level step: single device round-trip, matching Step 2's
        # actual intended implementation (move data on-device once before
        # the k-loop, come back to host once after). The three isolated
        # sub-step measurements above each pay transfer overhead
        # independently, which double/triple-counts a fixed cost that
        # Step 2's real implementation would only pay ONCE per level --
        # this is the realistic comparison point, not the isolated steps. ---
        print("\n[FUSED] Full level step (broadcast+pack+dedup, ONE transfer in/out)")
        (cpu_fused_states, cpu_fused_coeffs), cpu_min, cpu_mean = time_it(
            cpu_fused_level_step, current_states, current_coeffs,
            nz_i, nz_j, nz_val, k, n
        )
        print(f"    CPU: min={cpu_min*1000:.2f}ms  mean={cpu_mean*1000:.2f}ms")

        (gpu_fused_states, gpu_fused_coeffs), gpu_min, gpu_mean = time_it(
            gpu_fused_level_step, current_states, current_coeffs,
            nz_i, nz_j, nz_val, k, n
        )
        print(f"    GPU: min={gpu_min*1000:.2f}ms  mean={gpu_mean*1000:.2f}ms  "
              f"(single host<->device round trip for the whole level)")

        # states/coeffs may come out in different orders (np.unique on CPU
        # vs GPU need not sort identically in edge cases) -- compare as sets
        # of (state, coeff) pairs rather than assuming row order matches.
        def to_state_dict(states, coeffs):
            return {tuple(int(x) for x in s): c for s, c in zip(states, coeffs)}

        cpu_dict = to_state_dict(cpu_fused_states, cpu_fused_coeffs)
        gpu_dict = to_state_dict(gpu_fused_states, gpu_fused_coeffs)
        same_keys = set(cpu_dict.keys()) == set(gpu_dict.keys())
        same_vals = same_keys and all(
            np.allclose(cpu_dict[key], gpu_dict[key], atol=1e-8) for key in cpu_dict
        )
        correctness = same_keys and same_vals
        speedup = cpu_min / gpu_min if gpu_min > 0 else float('inf')
        print(f"    Correctness: {'OK' if correctness else 'MISMATCH -- DO NOT TRUST THIS PATH'}")
        print(f"    Speedup (CPU_min / GPU_min): {speedup:.2f}x")
        results[(k, 'FUSED_LEVEL')] = (speedup, correctness)

    # ------------------------------------------------------------------ #
    # Verdict
    # ------------------------------------------------------------------ #
    print(f"\n{'='*70}")
    print("VERDICT")
    print(f"{'='*70}")

    if not results:
        print("No levels were benchmarked (truncation may have terminated "
              "before reaching k=3/4 -- check level snapshot shapes above).")
        return

    all_correct = all(c for _, c in results.values())

    print(f"\nPer-step speedups (CPU_min / GPU_min):")
    for (k, step), (speedup, correct) in sorted(results.items()):
        flag = "" if correct else "  ** CORRECTNESS FAILED, IGNORE THIS NUMBER **"
        print(f"  k={k}  {step:12s}  {speedup:6.2f}x{flag}")

    if not all_correct:
        print("\n*** At least one result's GPU output did not match CPU. ***")
        print("*** DO NOT PROCEED to Step 2 until this is resolved. ***")
        return

    # The isolated sub-step numbers (broadcast/pack/dedup measured
    # separately, each paying its own transfer cost) are diagnostic --
    # they show WHERE the cost lives, but they are NOT representative of
    # Step 2's actual implementation, which keeps data on-device across
    # the whole per-level computation and only transfers once per level.
    # The FUSED_LEVEL numbers are the realistic comparison and are what
    # the proceed/stop decision should actually be based on.
    isolated_speedups = [s for (k, step), (s, c) in results.items()
                        if step != 'FUSED_LEVEL' and c]
    fused_speedups = [s for (k, step), (s, c) in results.items()
                      if step == 'FUSED_LEVEL' and c]

    if isolated_speedups:
        print(f"\nIsolated sub-step speedups (diagnostic only, NOT the decision basis):")
        print(f"  min={min(isolated_speedups):.2f}x  mean={np.mean(isolated_speedups):.2f}x")
        print(f"  (pack is consistently <1x in isolation -- this is expected: it's a")
        print(f"   tiny matmul whose isolated cost is almost entirely fixed transfer")
        print(f"   overhead, not compute. This does NOT mean pack should stay on CPU")
        print(f"   in a fused pipeline -- see FUSED_LEVEL below for the real answer.)")

    if not fused_speedups:
        print("\nNo FUSED_LEVEL results available -- cannot make a real-pipeline "
              "decision. Check the [FUSED] section output above for errors.")
        return

    min_fused = min(fused_speedups)
    mean_fused = np.mean(fused_speedups)
    print(f"\nFUSED_LEVEL speedups (the realistic, decision-relevant numbers):")
    print(f"  min={min_fused:.2f}x  mean={mean_fused:.2f}x")

    if min_fused >= 3.0:
        print("\n--> GPU clears the 3-5x bar on the FUSED level computation. Proceed to Step 2,")
        print("    implementing _build_squeezed_vacuum_vectorized_boson's per-level loop body")
        print("    as a single on-device computation (broadcast+pack+dedup), matching")
        print("    gpu_fused_level_step above, NOT as three separately-dispatched GPU calls.")
    elif mean_fused >= 3.0:
        print("\n--> GPU clears the bar on AVERAGE across truncation levels but not on every")
        print("    level measured. Consider whether the levels that do clear the bar (likely")
        print("    the higher-k, larger-state-count ones) are worth a conditional dispatch")
        print("    (e.g. only use GPU once num_states exceeds some threshold, falling back")
        print("    to CPU for small early levels where transfer overhead dominates).")
    else:
        print("\n--> GPU does NOT clear the 3-5x bar even with a fused, single-transfer")
        print("    pipeline. STOP HERE. Per the handoff plan, implementing Step 2 would add")
        print("    a dependency and complexity for no real benefit. Do not proceed.")


if __name__ == "__main__":
    run_benchmark()