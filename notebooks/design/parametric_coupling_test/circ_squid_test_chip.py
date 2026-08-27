import qiskit_metal as metal
from qiskit_metal import designs, draw, MetalGUI, Dict
import numpy as np
import pint
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science', 'no-latex', 'notebook'])
from qiskit_metal.qlibrary.qubits.circle_transmon_squid import CircTransmonSQUID

design = designs.DesignPlanar()
gui = MetalGUI(design)

design.overwrite_enabled = True
design.chips.main.size.size_x = '5mm'
design.chips.main.size.size_y = '2.5mm'

options = Dict(
    pos_x='0um', pos_y='0um', orientation='0', chip='main',
    cpw_width='10um', cpw_gap='6um',

    pad_radius='90um',
    pad_gap='10um',

    jj_options=Dict(
        jj_width='5um',
        jj_angle='0',      # single JJ, straight down, away from the SQUID
        L_j='10nH',
        C_j='2fF',
        export_mask=False,
        jj_sim_gap='3um',
    ),

    squid_options=Dict(
        theta_2='180',        # SQUID sits at the top of the pad
        delta_angle='5.0',    # wide enough to see the two prongs clearly
        g1='4um',

        l1='60um',
        w1='6um',
        l2='60um',

        jj_a_options=Dict(
            jj_width='4um', L_j='15nH', C_j='1.5fF',
            export_mask=False, jj_sim_gap='3um',
        ),
        jj_b_options=Dict(
            jj_width='6um', L_j='8nH', C_j='1.8fF',
            export_mask=False, jj_sim_gap='3um',
        ),
    ),
)
try:
    qubit_1.delete()
except NameError : pass
qubit_1 = CircTransmonSQUID(design, 'Q1_SQUID', options=options)

gui.rebuild()
gui.autoscale()

a_gds = design.renderers.gds
a_gds.options['path_filename'] = '../resources/Fake_Junctions.GDS'  # only matters if export touches junction cells
a_gds.export_to_gds('parametric_coupling_test.gds', highlight_qcomponents=['Q1_SQUID'])