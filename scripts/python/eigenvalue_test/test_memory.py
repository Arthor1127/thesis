import subprocess
import time
import numpy as np
from scipy.optimize import curve_fit

def get_peak_memory_mb(cmd):
    """Runs a command and returns the peak Resident Set Size (RSS) in MB."""
    # We use GNU 'time -v' or a simple polling loop if time -v is available,
    # but a robust cross-platform way in Linux is wrapping with time -f "%M"
    time_cmd = ["/usr/bin/time", "-f", "%M"] + cmd
    
    start_time = time.time()
    result = subprocess.run(time_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
    if result.returncode != 0:
        print(f"Error running command: {result.stderr}")
        return None
    
    # %M in GNU time outputs max resident set size in Kilobytes
    lines = result.stderr.strip().split("\n")
    # Usually the last line or the line containing the metric holds it, 
    # let's find the numeric value in stderr
    for line in reversed(lines):
        try:
            peak_kb = float(line.strip())
            return peak_kb / 1024.0 # Convert to MB
        except ValueError:
            continue
    return None

def memory_model(P, c1, c2):
    """
    Fitting function based on the scaling law: M(P) = c1 + c2 * P^(1/3)
    (Assuming N is fixed for a given geometry/mesh, so N factors are absorbed into c1 and c2)
    """
    return c1 + c2 * (P ** (1.0 / 3.0))

def run_benchmark():
    core_counts = [1, 2, 4]
    measured_mem_mb = []
    
    base_command_template = [
        "python", "build_and_eigenmode.py",
        "--ground-pos", "3.95",
        "--palace-bin", "/home/ruiz/repo/spack/opt/spack/linux-zen3/palace-develop-blpbtxe43ne2or3yo7gou6ycmdwdxuzl/bin/palace-x86_64.bin",
        "--starting-freq", "13e9",
        "--target-freq-ghz", "15",
        "--number-of-freqs", "5",
        "--amr-max-its", "0",
        "--bg-min-size-um", "8", "--bg-max-size-um", "100",
        "--resonator-bg-min-size-um", "15", "--resonator-bg-max-size-um", "150",
        "--strap-taper-dist-max-um", "40",
        "--outdir", "/home/ruiz/Documents/thesis/scripts/python/eigenvalue_test/output_bench"
    ]

    print("=== Starting Memory Scaling Benchmark ===")
    for p in core_counts:
        cmd = base_command_template + ["--num-cpus", str(p)]
        print(f"Running simulation with P = {p} cores...")
        
        # Note: We wrap mpirun or let python handle it depending on how build_and_eigenmode calls it
        mem_mb = get_peak_memory_mb(cmd)
        if mem_mb is None:
            print(f"Failed to capture memory for P = {p}")
            return
        
        print(f"  -> Peak Memory: {mem_mb:.2f} MB")
        measured_mem_mb.append(mem_mb)

    P_arr = np.array(core_counts, dtype=float)
    M_arr = np.array(measured_mem_mb, dtype=float)

    # Fit curve to M(P) = c1 + c2 * P^(1/3)
    popt, _ = curve_fit(memory_model, P_arr, M_arr)
    c1, c2 = popt

    print("\n=== Fit Results ===")
    print(f"Fitted model: M(P) = {c1:.2f} + {c2:.2f} * P^(1/3)")

    # Extrapolate for larger core counts (e.g., 8, 16, 32)
    target_cores = [8, 16, 32, 64]
    print("\n=== Extrapolated Memory Predictions ===")
    for target_p in target_cores:
        pred_mem = memory_model(target_p, c1, c2)
        print(f"  P = {target_p:2d} cores: ~{pred_mem:.2f} MB (~{pred_mem/1024:.2f} GB total)")

if __name__ == "__main__":
    run_benchmark()