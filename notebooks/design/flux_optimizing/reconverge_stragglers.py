#!/usr/bin/env python
"""Re-solve the sweep points that hit SOLVER_MAXITS without converging.

All 46 stragglers found (2026-09-04) hit exactly n_ksp=1201 -- they ran BEFORE
SOLVER_MAXITS was bumped 1200->2000 mid-sweep, not a genuinely-unreachable
tolerance. This script re-solves each with a further-boosted iteration cap
(default 6000, since it's only 46 points and correctness matters more than
wall time here) and APPENDS the corrected row to the row's *original* shard
CSV. It never overwrites -- aggregate_sweep.py's
``drop_duplicates(subset=AXES, keep='last')`` picks the newer (appended) row
automatically.

One SGE array task = one straggler shard, so concurrent tasks never write the
same file (avoids any append race). Discovery re-scans the CSVs at every
launch, so this is safe to resubmit if a handful still don't converge.

    python reconverge_stragglers.py --list                      # just show them
    python reconverge_stragglers.py --task-id 3 --n-tasks 25 \
        --outdir sweep --palace-bin /path/to/palace --num-cpus 12 \
        --solver-maxits 6000 --gmres-maxsize 300
"""
import argparse
import csv
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sweep_flux as sf  # noqa: E402


def find_stragglers(outdir):
    """{shard_path: [param_dict, ...]} for every row with converged == 'False'."""
    by_shard = {}
    for f in sorted(glob.glob(os.path.join(outdir, 'shard_*.csv'))):
        rows = []
        for r in csv.DictReader(open(f)):
            if (r.get('flux_per_amp') or '').strip() and r.get('converged') == 'False':
                rows.append({a: float(r[a]) for a in sf.AXIS_ORDER})
        if rows:
            by_shard[f] = rows
    return by_shard


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--outdir', default='sweep')
    ap.add_argument('--list', action='store_true', help='print stragglers and exit')
    ap.add_argument('--task-id', type=int)
    ap.add_argument('--n-tasks', type=int)
    ap.add_argument('--palace-bin')
    ap.add_argument('--num-cpus', type=int, default=12)
    ap.add_argument('--solver-maxits', type=int, default=6000)
    ap.add_argument('--gmres-maxsize', type=int, default=300)
    a = ap.parse_args()

    by_shard = find_stragglers(a.outdir)
    total = sum(len(v) for v in by_shard.values())

    if a.list or not a.task_id:
        print(f"{total} not-converged rows across {len(by_shard)} shards:")
        for shard, rows in by_shard.items():
            for p in rows:
                tag = ' '.join(f"{k}={p[k]:g}" for k in sf.AXIS_ORDER)
                print(f"  {shard}: {tag}")
        if not a.task_id:
            return

    shards = sorted(by_shard)  # deterministic -- same scan, same order every task
    my_shards = shards[a.task_id - 1::a.n_tasks]
    if not my_shards:
        print(f"task {a.task_id}/{a.n_tasks}: no shards assigned, nothing to do")
        return

    sf._patch_gmsh_multipoly()
    for shard in my_shards:
        pts = by_shard[shard]
        workdir = shard[:-4] + '_reconv_work'
        os.makedirs(workdir, exist_ok=True)
        fh = open(shard, 'a', newline='')
        w = csv.DictWriter(fh, fieldnames=sf.FIELDS)
        for i, params in enumerate(pts):
            tag = ' '.join(f"{k}={v:g}" for k, v in params.items())
            print(f"[{shard}] {i+1}/{len(pts)}  {tag}", flush=True)
            row = sf.run_one(params, a.palace_bin, a.num_cpus, workdir,
                             solver_maxits=a.solver_maxits,
                             gmres_maxsize=a.gmres_maxsize)
            w.writerow(row)
            fh.flush()
            status = row.get('converged')
            print(f"  -> Phi0/mA={row.get('Phi0_per_mA')}  converged={status}  "
                 f"n_ksp={row.get('n_ksp')}  {row.get('seconds')}s", flush=True)
        fh.close()
        import shutil
        shutil.rmtree(workdir, ignore_errors=True)
    print(f"task {a.task_id}/{a.n_tasks} done: {len(my_shards)} shards, "
         f"{sum(len(by_shard[s]) for s in my_shards)} points")


if __name__ == '__main__':
    main()
