"""Aggregate sweep_run/*/result.json into a per-region table."""
import glob
import json
import os
import sys

SWEEP_DIR = sys.argv[1] if len(sys.argv) > 1 else "sweep_run"

rows = []
for path in glob.glob(os.path.join(SWEEP_DIR, "*", "result.json")):
    with open(path) as f:
        rows.append(json.load(f))

by_region = {}
for r in rows:
    by_region.setdefault(r["region"], []).append(r)

for region, items in sorted(by_region.items()):
    items.sort(key=lambda r: -r["min_size_um"])
    print(f"\n=== {region} ===")
    print(f"{'min_um':>8} {'status':>10} {'wall_s':>8} {'elements':>10}")
    for r in items:
        wall = r.get("wall_s")
        wall_s = f"{wall:.1f}" if wall is not None else "-"
        elements = r.get("elements", "-")
        print(f"{r['min_size_um']:>8} {r['status']:>10} {wall_s:>8} {str(elements):>10}")
