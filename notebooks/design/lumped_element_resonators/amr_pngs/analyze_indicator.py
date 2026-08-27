import numpy as np
import pyvista as pv

# Component bounding boxes (mm, matching the mesh's coordinate units), from
# build_design.py's qgeometry tables, padded slightly.
CAP_BOX = dict(xmin=-0.06, xmax=0.06, ymin=-0.16, ymax=0.015)
IND_BOX = dict(xmin=-0.01, xmax=1.12, ymin=-0.36, ymax=0.21)


def classify(x, y):
    in_cap = (CAP_BOX['xmin'] <= x <= CAP_BOX['xmax']) and (CAP_BOX['ymin'] <= y <= CAP_BOX['ymax'])
    in_ind = (IND_BOX['xmin'] <= x <= IND_BOX['xmax']) and (IND_BOX['ymin'] <= y <= IND_BOX['ymax'])
    if in_cap:
        return 'cap'
    if in_ind:
        return 'inductor'
    return 'other'


for name in ['indicator_iter1', 'indicator_iter2', 'indicator_final']:
    m = pv.read(f'{name}/data.pvtu')
    pts = m.points
    ind = np.asarray(m.point_data['Indicator'])

    cats = np.array([classify(x, y) for x, y, z in pts])

    print(f'\n=== {name} ({m.n_points} points) ===')
    for c in ['cap', 'inductor', 'other']:
        mask = cats == c
        n = mask.sum()
        if n == 0:
            print(f'  {c:10s}: 0 points')
            continue
        vals = ind[mask]
        print(f'  {c:10s}: n={n:8d}  sum={vals.sum():.4e}  mean={vals.mean():.4e}  max={vals.max():.4e}')
