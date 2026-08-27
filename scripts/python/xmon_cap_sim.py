"""
Xmon capacitance simulation (SQDMetal + Palace), headless/cluster version of
notebooks/design/tests/test.ipynb.

Builds the same qiskit-metal planar design (launchpads, transmission line,
transmon cross, readout resonator), meshes it with GMSH, runs the Palace
capacitance simulation, and writes the resulting capacitance matrix to a
.dat file. No Qt/GUI is created so this runs on a display-less cluster node.

Usage:
    python xmon_cap_sim.py --num-cpus 8 --palace-dir /path/to/palace \
        --sim-dir xmon_cap_sim --outdir .
"""
import os

os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
os.environ["PMIX_MCA_gds"] = "hash"
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import argparse

import numpy as np

from qiskit_metal import designs, Dict
from qiskit_metal.qlibrary.qubits.transmon_cross import TransmonCross
from qiskit_metal.qlibrary.tlines.meandered import RouteMeander
from qiskit_metal.qlibrary.tlines.pathfinder import RoutePathfinder
from qiskit_metal.qlibrary.terminations.launchpad_wb import LaunchpadWirebond
from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround

from SQDMetal.PALACE.Capacitance_Simulation import PALACE_Capacitance_Simulation


def build_design():
    design = designs.DesignPlanar({}, overwrite_enabled=True)

    design.chips.main.size.size_x = '4.6mm'
    design.chips.main.size.size_y = '2.4mm'
    design.chips.main.size.size_z = '-280um'
    design.chips.main.size.center_x = '0mm'
    design.chips.main.size.center_y = '-1mm'

    design.variables['cpw_width'] = '10 um'  # S
    design.variables['cpw_gap'] = '6 um'     # W

    launch_options1 = dict(chip='main', pos_x='-2mm', pos_y='0mm', orientation='360',
                            lead_length='30um', pad_height='103um', pad_width='103um', pad_gap='60um')
    LaunchpadWirebond(design, 'LP1', options=launch_options1)

    launch_options2 = dict(chip='main', pos_x='2mm', pos_y='0mm', orientation='180',
                            lead_length='30um', pad_height='103um', pad_width='103um', pad_gap='60um')
    LaunchpadWirebond(design, 'LP2', options=launch_options2)

    RoutePathfinder(design, 'TL', options=dict(
        chip='main', trace_width='10um', trace_gap='6um', fillet='90um',
        hfss_wire_bonds=True, lead=dict(end_straight='0.1mm'),
        pin_inputs=Dict(
            start_pin=Dict(component='LP1', pin='tie'),
            end_pin=Dict(component='LP2', pin='tie'))))

    TransmonCross(design, 'Q1', options=dict(
        pos_x='0.6075mm', pos_y='-1.464',
        connection_pads=dict(
            bus_01=dict(connector_location='180', claw_length='95um'),
            readout=dict(connector_location='0')),
        fl_options=dict()))

    OpenToGround(design, 'otg1', options=dict(chip='main', pos_x='-0.2mm', pos_y='-40um', orientation=180))

    RouteMeander(design, 'resonator1', Dict(
        trace_width='10um', trace_gap='6um', total_length='3.7mm',
        hfss_wire_bonds=False, fillet='99.9 um',
        lead=dict(start_straight='300um'),
        pin_inputs=Dict(
            start_pin=Dict(component='otg1', pin='open'),
            end_pin=Dict(component='Q1', pin='readout'))))

    design.rebuild()
    return design


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--num-cpus', type=int, default=int(os.environ.get('NSLOTS', 8)),
                         help='Number of MPI ranks Palace uses (defaults to $NSLOTS from SGE)')
    parser.add_argument('--palace-dir', default=os.environ.get(
        'PALACE_DIR', '/home/ruiz/repo/spack/opt/spack/linux-zen3/palace-develop-blpbtxe43ne2or3yo7gou6ycmdwdxuzl/bin/palace'),
        help='Path to the palace binary')
    parser.add_argument('--sim-dir', default='xmon_cap_sim',
                         help='Simulation name; SQDMetal creates a folder with this name for mesh/config/output')
    parser.add_argument('--sim-parent-dir', default='',
                         help='Parent directory in which the simulation folder is created')
    parser.add_argument('--outdir', default='.', help='Directory to write the .dat result files into')
    return parser.parse_args()


def main():
    args = parse_args()

    design = build_design()

    user_defined_options = {
        "mesh_refinement": 0,
        "dielectric_material": "silicon",
        "solver_order": 2,
        "solns_to_save": 3,
        "solver_tol": 1.0e-8,
        "solver_maxits": 200,
        "fillet_resolution": 12,
        "palace_dir": args.palace_dir,
        "num_cpus": args.num_cpus,
    }

    cap_sim = PALACE_Capacitance_Simulation(
        name=args.sim_dir,
        metal_design=design,
        sim_parent_directory=args.sim_parent_dir,
        mode='simPC',
        meshing='GMSH',
        user_options=user_defined_options,
        create_files=True)

    cap_sim.add_metallic(1)
    cap_sim.add_ground_plane()
    cap_sim.fine_mesh_components(['Q1'], min_size=12e-6, max_size=100e-6, taper_dist_min=10e-6, metals_only=False)
    cap_sim.fine_mesh_components(['resonator1', 'TL', 'LP1', 'LP2'], min_size=12e-6, max_size=120e-6, taper_dist_min=10e-6)

    cap_sim.prepare_simulation()
    cap_mat = cap_sim.run()

    os.makedirs(args.outdir, exist_ok=True)
    out_F = os.path.join(args.outdir, f'{args.sim_dir}_capacitance_F.dat')
    out_fF = os.path.join(args.outdir, f'{args.sim_dir}_capacitance_fF.dat')
    np.savetxt(out_F, cap_mat, header='Capacitance matrix (F)')
    np.savetxt(out_fF, cap_mat * 1e15, header='Capacitance matrix (fF)')

    print(f'Saved capacitance matrix to {out_F} and {out_fF}')


if __name__ == '__main__':
    main()
