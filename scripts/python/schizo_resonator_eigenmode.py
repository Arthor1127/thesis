"""
Eigenmode simulation (SQDMetal + Palace) of the periodically-grounded
lambda/2 resonator ("schizo resonator"), headless/cluster version of
notebooks/design/schizo_resonators/schizo_resonators.ipynb.

Builds the qiskit-metal planar design (open-open RouteMeanderGrounded
resonator with grounding straps along its length), meshes it with GMSH,
runs the Palace eigenmode simulation, and prints/saves the resulting
eigenfrequencies + Q (from eig.csv). No Qt/GUI is created so this runs on
a display-less cluster node.

Usage:
    python schizo_resonator_eigenmode.py --num-cpus 32 --palace-dir /path/to/palace \
        --sim-dir res_eigen_test --outdir .

Requires RouteMeanderGrounded to be installed in this environment's
qiskit_metal (qlibrary/tlines/meandered_grounded.py) - copy it over if
it's not already there.
"""
import os

os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
os.environ["PMIX_MCA_gds"] = "hash"
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import argparse
import shutil

import pandas as pd

from qiskit_metal import designs, Dict
from qiskit_metal.qlibrary.tlines.meandered_grounded import RouteMeanderGrounded
from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround

from SQDMetal.PALACE.Eigenmode_Simulation import PALACE_Eigenmode_Simulation
from SQDMetal.Utilities.Materials import MaterialInterface


def build_design(total_length, ground_strap_positions, ground_strap_width):
    design = designs.DesignPlanar({}, overwrite_enabled=True)

    design.chips.main.size.size_x = '4.8mm'
    design.chips.main.size.size_y = '2.4mm'
    design.chips.main.size.size_z = '500um'
    design.chips.main.size.center_x = '0mm'
    design.chips.main.size.center_y = '-1mm'

    design.variables['cpw_width'] = '10 um'  # S
    design.variables['cpw_gap'] = '6 um'     # W

    OpenToGround(design, 'otg1', options=dict(chip='main', pos_x='-0.2mm', pos_y='-40um', orientation=180))
    OpenToGround(design, 'otg2', options=dict(chip='main', pos_x='0mm', pos_y='-1.35mm', orientation=-90))

    RouteMeanderGrounded(design, 'resonator1', Dict(
        trace_width='10um',
        trace_gap='6um',
        total_length=total_length,
        hfss_wire_bonds=False,
        fillet='99.9 um',
        lead=dict(start_straight='300um'),
        ground_straps=dict(positions=ground_strap_positions, width=ground_strap_width),
        pin_inputs=Dict(
            start_pin=Dict(component='otg1', pin='open'),
            end_pin=Dict(component='otg2', pin='open'))))

    design.rebuild()
    return design


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--num-cpus', type=int, default=int(os.environ.get('NSLOTS', 8)),
                         help='Number of MPI ranks Palace uses (defaults to $NSLOTS from SGE)')
    parser.add_argument('--palace-dir', default=os.environ.get('PALACE_DIR', 'palace'),
                         help='Path to the palace binary')
    parser.add_argument('--sim-dir', default='res_eigen_test',
                         help='Simulation name; SQDMetal creates a folder with this name for mesh/config/output')
    parser.add_argument('--sim-parent-dir', default='',
                         help='Parent directory in which the simulation folder is created')
    parser.add_argument('--outdir', default='.', help='Directory to copy the eig.csv/field-plot results into')
    parser.add_argument('--total-length', default='3.6mm')
    parser.add_argument('--ground-strap-positions', nargs='+', default=['1.2mm', '2.4mm'])
    parser.add_argument('--ground-strap-width', default='10um')
    parser.add_argument('--starting-freq', type=float, default=10e9)
    parser.add_argument('--number-of-freqs', type=int, default=5)
    parser.add_argument('--solver-tol', type=float, default=1.0e-7)
    return parser.parse_args()


def main():
    args = parse_args()

    design = build_design(args.total_length, args.ground_strap_positions, args.ground_strap_width)

    user_defined_options = {
        "dielectric_material": "silicon",
        "starting_freq": args.starting_freq,
        "number_of_freqs": args.number_of_freqs,
        "solns_to_save": args.number_of_freqs,
        "solver_order": 2,
        "solver_tol": args.solver_tol,
        "solver_maxits": 200,
        "fillet_resolution": 12,
        "palace_dir": args.palace_dir,
        "num_cpus": args.num_cpus,
    }

    eigen_sim = PALACE_Eigenmode_Simulation(
        name=args.sim_dir,
        metal_design=design,
        sim_parent_directory=args.sim_parent_dir,
        mode='simPC',
        meshing='GMSH',
        user_options=user_defined_options,
        create_files=True)

    eigen_sim.add_metallic(1)
    eigen_sim.add_ground_plane()
    eigen_sim.fine_mesh_components(['resonator1'], min_size=16e-6, max_size=100e-6, taper_dist_min=10e-6, metals_only=False)
    eigen_sim.setup_EPR_interfaces(
        metal_air=MaterialInterface('Aluminium-Vacuum'),
        substrate_air=MaterialInterface('Silicon-Vacuum'),
        substrate_metal=MaterialInterface('Silicon-Aluminium'))

    eigen_sim.prepare_simulation()
    eigen_sim.run()

    eig_df = pd.read_csv(eigen_sim._output_data_dir + '/eig.csv')
    eig_df.columns = eig_df.columns.str.strip()
    print(eig_df)

    os.makedirs(args.outdir, exist_ok=True)
    shutil.copy(eigen_sim._output_data_dir + '/eig.csv', os.path.join(args.outdir, f'{args.sim_dir}_eig.csv'))
    print(f'Saved eigenfrequencies to {os.path.join(args.outdir, f"{args.sim_dir}_eig.csv")}')


if __name__ == '__main__':
    main()
