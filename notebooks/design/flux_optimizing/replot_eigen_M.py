#!/usr/bin/env python
"""Redraw an eigen_M_cyc*.npz bar chart in the thesis house style (no need to
touch the cluster -- everything needed was already saved to disk).

    python replot_eigen_M.py eigen/small/eigen_M_cyc1.npz
    python replot_eigen_M.py eigen/best/eigen_M_cyc1.npz
"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import thesis_style as ts

ts.use()

path = sys.argv[1]
d = np.load(path)
n_loops = len(d['loop'])
f_r = float(d['f_r_Hz'])
M_total = float(d['M_total_H'])
M_magneto = float(d['M_magneto_H']) if 'M_magneto_H' in d and np.isfinite(d['M_magneto_H']) else None

fig, ax = plt.subplots(2, 1, figsize=ts.figsize(1.5, ratio=0.6), sharex=True)
ax[0].bar(d['loop'], d['Phi_signed'] * 1e15, color=ts.COLORS['blue'])
ax[0].set_title(f"$\\Phi_k$ (fWb)   z={float(d['z_offset_um']):.1f} um above plane")
ax[1].bar(d['loop'], d['M_k'] * 1e12, color=ts.COLORS['teal'])
ax[1].set_title('$M_k = \\Phi_k / I_{short}$  (pH)')
ax[1].set_xlabel('loop index')

ttl = f"{path}\ncycle {int(d['cycle'])}   $f_r$={f_r/1e9:.3f} GHz   $M_{{total}}$={M_total*1e12:.2f} pH"
if M_magneto is not None:
    ttl += f"   $M_{{magnetostatic}}$={M_magneto*1e12:.2f} pH   ratio={M_total/M_magneto:.3f}"
fig.suptitle(ttl, fontsize=9)
ts.finalize(ax)
fig.tight_layout()
out = path.rsplit('.', 1)[0] + '_thesis.png'
fig.savefig(out, dpi=160, facecolor='white')
print('wrote', out)
