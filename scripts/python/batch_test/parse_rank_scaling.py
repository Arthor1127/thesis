"""
Summarise the rank-scaling experiment.

Usage:
    python parse_rank_scaling.py rank08.o* rank16.o* rank32.o*

Reads Palace's Elapsed Time Report and the mesh/rank info out of each
log, then reports speedup, parallel efficiency, and the Amdahl parallel
fraction implied by each consecutive pair.

The number that matters is Preconditioner - it is ~73% of runtime on
this problem, so it sets the ceiling for everything else.
"""
import re
import sys

# Lines of the Elapsed Time Report worth extracting. Palace prints this
# as a FLAT list that sums to Total (the indentation is cosmetic, not
# nesting - verified against the reported totals).
WANT = ["Preconditioner", "Coarse Solve", "Linear Solve",
        "Div.-Free Projection", "Total"]


def parse(path):
    with open(path, errors="replace") as f:
        txt = f.read()

    out = {"path": path}

    m = re.search(r"Running with (\d+) MPI process", txt)
    if m:
        out["ranks"] = int(m.group(1))
    else:
        m = re.search(r"NSLOTS=(\d+)", txt)
        out["ranks"] = int(m.group(1)) if m else None

    m = re.search(r"Host:\s*(\S+)", txt)
    out["host"] = m.group(1) if m else "?"

    m = re.search(r"global unknowns\s*=\s*([\d,]+)", txt)
    out["unknowns"] = int(m.group(1).replace(",", "")) if m else None

    # "  Preconditioner   328.410   328.868   328.771" -> take Avg (last)
    for key in WANT:
        pat = re.escape(key) + r"\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)"
        m = re.search(pat, txt)
        if m:
            out[key] = float(m.group(3))
        else:
            # Total is printed with a single value repeated
            m2 = re.search(re.escape(key) + r"\s+([\d.]+)", txt)
            out[key] = float(m2.group(1)) if m2 else None

    m = re.search(r"All modes found \(GHz\): \[([^\]]*)\]", txt)
    if m:
        out["f1"] = float(m.group(1).split(",")[0])

    return out


def amdahl_fraction(t1, t2, n1, n2):
    """Parallel fraction p implied by two timings, from
    T(n) = T1*((1-p) + p/n). Returns None if the pair is degenerate."""
    if not (t1 and t2) or n1 == n2:
        return None
    # T2/T1 = ((1-p) + p/n2) / ((1-p) + p/n1)
    r = t2 / t1
    denom = (1 / n2 - 1) - r * (1 / n1 - 1)
    if abs(denom) < 1e-12:
        return None
    return (r - 1) / denom


def main():
    paths = sys.argv[1:]
    if not paths:
        print(__doc__)
        sys.exit(1)

    runs = [parse(p) for p in paths]
    runs = [r for r in runs if r.get("ranks")]
    runs.sort(key=lambda r: r["ranks"])

    if not runs:
        print("No parseable runs found.")
        sys.exit(1)

    hosts = {r["host"] for r in runs}
    unks = {r["unknowns"] for r in runs if r["unknowns"]}
    print("=== Rank scaling ===")
    if len(hosts) > 1:
        print(f"!! WARNING: runs span multiple hosts {hosts} - this compares")
        print("!! HARDWARE as much as rank count. Re-run pinned to one host.")
    if len(unks) > 1:
        print(f"!! WARNING: unknown counts differ {unks} - the mesh changed")
        print("!! between runs, so this is NOT a clean rank comparison.")
    if len(hosts) == 1 and len(unks) <= 1:
        print(f"host={hosts.pop()}  unknowns={unks.pop() if unks else '?'}  "
              f"(clean comparison)")
    print()

    hdr = f"{'ranks':>6} {'Precond':>9} {'Total':>9} {'f1(GHz)':>9}"
    print(hdr)
    print("-" * len(hdr))
    for r in runs:
        pc = f"{r['Preconditioner']:.1f}" if r.get("Preconditioner") else "?"
        tt = f"{r['Total']:.1f}" if r.get("Total") else "?"
        f1 = f"{r['f1']:.4f}" if r.get("f1") else "?"
        print(f"{r['ranks']:>6} {pc:>9} {tt:>9} {f1:>9}")

    base = runs[0]
    print(f"\n=== Speedup vs {base['ranks']} ranks ===")
    print(f"{'ranks':>6} {'ideal':>7} {'Precond':>9} {'eff':>6} "
          f"{'Total':>9} {'eff':>6}")
    print("-" * 48)
    for r in runs:
        ideal = r["ranks"] / base["ranks"]
        row = [f"{r['ranks']:>6}", f"{ideal:>6.2f}x"]
        for key in ("Preconditioner", "Total"):
            if base.get(key) and r.get(key):
                sp = base[key] / r[key]
                row.append(f"{sp:>8.2f}x")
                row.append(f"{sp/ideal:>5.0%}")
            else:
                row += ["?".rjust(9), "?".rjust(6)]
        print(" ".join(row))

    print("\n=== Implied Amdahl parallel fraction (consecutive pairs) ===")
    for a, b in zip(runs, runs[1:]):
        for key in ("Preconditioner", "Total"):
            p = amdahl_fraction(a.get(key), b.get(key), a["ranks"], b["ranks"])
            if p is not None:
                cap = 1 / (1 - p) if p < 1 else float("inf")
                print(f"  {a['ranks']:>2}->{b['ranks']:<2} {key:<16} "
                      f"p={p:6.1%}  max speedup at infinite ranks = {cap:.1f}x")

    print("\nReading it:")
    print("  eff >~80% to 32 ranks  -> scaling is fine, use more ranks.")
    print("  eff falls off sharply  -> subdomain overhead (ghost cells +")
    print("                            multigrid coarsening fragmentation)")
    print("                            dominates; that is the case FOR a")
    print("                            hybrid MPI+OpenMP build, which cuts")
    print("                            subdomain count at equal core count.")
    print("  f1 must match across all runs - if it does not, something")
    print("  other than rank count changed and the comparison is void.")


if __name__ == "__main__":
    main()
