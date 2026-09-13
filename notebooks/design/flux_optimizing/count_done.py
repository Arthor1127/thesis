#!/usr/bin/env python
"""Print '<n_done> <n_error>' across sweep/shard_*.csv."""
import csv
import glob
import sys

sweep_dir = sys.argv[1] if len(sys.argv) > 1 else 'sweep'
n = err = 0
for f in glob.glob(f'{sweep_dir}/shard_*.csv'):
    try:
        for r in csv.DictReader(open(f)):
            if (r.get('flux_per_amp') or '').strip():
                n += 1
            elif (r.get('error') or '').strip():
                err += 1
    except Exception:
        pass
print(f"{n} {err}")
