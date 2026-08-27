"""Tests for components.py / chip.py against toy capacitance matrices.

The headline check: compute_coupling's 0.5*Q^T Cinv Q result is verified
against a fully independent manual calculation using scipy.constants (not
just checked for "doesn't crash" -- see test_coupling_matches_manual_calc).
"""

import numpy as np
from scipy.constants import e, h, hbar

from design_helper import DesignHelper
from chip import Chip


def make_toy_chip(dh):
    cap_matrix_F = np.array([
        [100e-15, -5e-15],
        [-5e-15, 100e-15],
    ])
    node_labels = ['qubit_pad', 'res_pad']
    specs = [
        {'type': 'qubit', 'name': 'Q1', 'nodes': ['qubit_pad'],
         'EJ': dh._as_quantity(15, 'GHz'), 'EC': dh._as_quantity(0.25, 'GHz')},
        {'type': 'resonator', 'name': 'R1', 'nodes': ['res_pad'],
         'capacitance': dh._as_quantity(100, 'fF'),
         'omega_r': (2 * np.pi * dh._as_quantity(6, 'GHz')).to('1/s')},
    ]
    return Chip.from_simulation(dh, cap_matrix_F, node_labels, specs)


def test_coupling_matches_manual_calc():
    dh = DesignHelper(light_speed='3e8 m/s')
    chip = make_toy_chip(dh)
    g = chip.compute_coupling('Q1', 'R1').to('Hz').magnitude

    # fully independent manual calc, no DesignHelper/Chip code reused
    C = np.array([[100e-15, -5e-15], [-5e-15, 100e-15]])
    cinv_12 = np.linalg.inv(C)[0, 1]
    EJ_J, EC_J = 15e9 * h, 0.25e9 * h
    n_zpf = 0.5 * (EJ_J / (2 * EC_J)) ** 0.25
    q_zpf_qubit = 2 * e * n_zpf
    omega_res = 2 * np.pi * 6e9
    Z_res = 1 / (omega_res * 100e-15)
    q_zpf_res = np.sqrt(hbar / (2 * Z_res))
    g_manual_Hz = (0.5 * cinv_12 * q_zpf_qubit * q_zpf_res) / h

    assert abs(g - g_manual_Hz) / abs(g_manual_Hz) < 1e-9, (g, g_manual_Hz)
    print(f"PASS: g = {g/1e6:.4f} MHz matches manual calc {g_manual_Hz/1e6:.4f} MHz")


def test_duplicate_component_raises():
    dh = DesignHelper(light_speed='3e8 m/s')
    chip = make_toy_chip(dh)
    try:
        chip.add_component({'type': 'ground', 'name': 'Q1', 'nodes': ['qubit_pad']})
        assert False, "should have raised"
    except ValueError:
        print("PASS: duplicate component name raises ValueError")


def test_unknown_node_label_raises():
    dh = DesignHelper(light_speed='3e8 m/s')
    chip = make_toy_chip(dh)
    try:
        chip.add_component({'type': 'ground', 'name': 'GND', 'nodes': ['nope']})
        assert False, "should have raised"
    except ValueError:
        print("PASS: unknown node label raises ValueError")


def test_update_component_changes_live_coupling():
    dh = DesignHelper(light_speed='3e8 m/s')
    chip = make_toy_chip(dh)
    g_before = chip.compute_coupling('Q1', 'R1').to('MHz').magnitude
    chip.update_component('Q1', EJ=dh._as_quantity(16, 'GHz'))
    g_after = chip.compute_coupling('Q1', 'R1').to('MHz').magnitude
    assert g_before != g_after
    print(f"PASS: update_component changed g from {g_before:.4f} to {g_after:.4f} MHz")


def test_remove_component():
    dh = DesignHelper(light_speed='3e8 m/s')
    chip = make_toy_chip(dh)
    chip.remove_component('R1')
    assert 'R1' not in chip.components
    try:
        chip.compute_coupling('Q1', 'R1')
        assert False, "should have raised"
    except KeyError:
        print("PASS: removed component absent + compute_coupling raises KeyError")


def test_matrix_size_mismatch_raises():
    dh = DesignHelper(light_speed='3e8 m/s')
    chip = Chip(dh, ['a', 'b', 'c'])
    try:
        chip.set_capacitance_matrix(np.array([[1e-13, -5e-15], [-5e-15, 1e-13]]))
        assert False, "should have raised"
    except ValueError:
        print("PASS: matrix/node_labels size mismatch raises ValueError")


def test_floating_transmon_coupling_matches_manual_calc():
    """Multi-node (floating transmon) coupling: independent manual calc check.

    Convention confirmed against a real SQDMetal/Palace TransmonPocket sim:
    two-pad floating transmon has differential charge operator Q=Q_pad0-Q_pad1,
    so the effective Cinv element for coupling to a single-node resonator is
        Cinv_eff = Cinv[pad0, res] - Cinv[pad1, res]

    Toy 3-node matrix: [float_pad0, float_pad1, res_pad]
    Asymmetric mutual caps (-8 fF vs -3 fF) so the two terms don't cancel.
    """
    dh = DesignHelper(light_speed='3e8 m/s')

    # Asymmetric: pad0-res = -8 fF, pad1-res = -3 fF
    C = np.array([
        [100e-15, -40e-15, -8e-15],
        [-40e-15, 100e-15, -3e-15],
        [-8e-15,  -3e-15, 100e-15],
    ])
    chip = Chip.from_simulation(dh, C, ['float_pad0', 'float_pad1', 'res_pad'], [
        {'type': 'qubit', 'name': 'Qfloat', 'nodes': ['float_pad0', 'float_pad1'],
         'EJ': dh._as_quantity(15, 'GHz'), 'EC': dh._as_quantity(0.25, 'GHz')},
        {'type': 'resonator', 'name': 'R1', 'nodes': ['res_pad'],
         'capacitance': dh._as_quantity(100, 'fF'),
         'omega_r': (2 * np.pi * dh._as_quantity(6, 'GHz')).to('1/s')},
    ])

    g = chip.compute_coupling('Qfloat', 'R1').to('Hz').magnitude

    # Independent manual calc -- no DesignHelper/Chip code reused
    from scipy.constants import e, h, hbar
    Cinv = np.linalg.inv(C)
    cinv_eff = Cinv[0, 2] - Cinv[1, 2]   # differential-mode convention
    EJ_J, EC_J = 15e9 * h, 0.25e9 * h
    n_zpf = 0.5 * (EJ_J / (2 * EC_J)) ** 0.25
    q_zpf_qubit = 2 * e * n_zpf
    omega_res = 2 * np.pi * 6e9
    Z_res = 1 / (omega_res * 100e-15)
    q_zpf_res = np.sqrt(hbar / (2 * Z_res))
    g_manual_Hz = (0.5 * cinv_eff * q_zpf_qubit * q_zpf_res) / h

    assert abs(g - g_manual_Hz) / abs(g_manual_Hz) < 1e-9, (g, g_manual_Hz)
    print(f"PASS: floating transmon g = {g/1e6:.4f} MHz matches manual calc {g_manual_Hz/1e6:.4f} MHz")


def test_floating_transmon_coupling_symmetric_cancels():
    """When both pads couple equally to a neighbor, differential mode gives zero.

    Physical sanity check: a perfectly symmetric floating transmon (C[pad0,res]
    == C[pad1,res]) has no net coupling to a neighbor via the differential charge
    operator -- both contributions cancel. g should be exactly 0.
    """
    dh = DesignHelper(light_speed='3e8 m/s')

    # Perfectly symmetric: both pads couple equally to res_pad
    C = np.array([
        [100e-15, -40e-15, -5e-15],
        [-40e-15, 100e-15, -5e-15],
        [-5e-15,  -5e-15, 100e-15],
    ])
    chip = Chip.from_simulation(dh, C, ['float_pad0', 'float_pad1', 'res_pad'], [
        {'type': 'qubit', 'name': 'Qfloat', 'nodes': ['float_pad0', 'float_pad1'],
         'EJ': dh._as_quantity(15, 'GHz'), 'EC': dh._as_quantity(0.25, 'GHz')},
        {'type': 'resonator', 'name': 'R1', 'nodes': ['res_pad'],
         'capacitance': dh._as_quantity(100, 'fF'),
         'omega_r': (2 * np.pi * dh._as_quantity(6, 'GHz')).to('1/s')},
    ])

    g = chip.compute_coupling('Qfloat', 'R1').to('Hz').magnitude
    # Cinv[pad0, res] == Cinv[pad1, res] by symmetry -> cinv_eff = 0 -> g = 0
    assert abs(g) < 1e-6, f"expected ~0 Hz, got {g:.3e} Hz"
    print(f"PASS: symmetric floating transmon g = {g:.3e} Hz (correctly ~0)")


def test_coupling_matrix_excludes_ground():
    dh = DesignHelper(light_speed='3e8 m/s')
    n = 4
    C = np.eye(n) * 100e-15
    C[0, 1] = C[1, 0] = -5e-15
    C[0, 2] = C[2, 0] = -2e-15
    C[1, 2] = C[2, 1] = -1e-15
    C[3, 3] = 500e-15
    node_labels = ['qubit_pad', 'res1_pad', 'res2_pad', 'ground']
    specs = [
        {'type': 'qubit', 'name': 'Q1', 'nodes': ['qubit_pad'],
         'EJ': dh._as_quantity(15, 'GHz'), 'EC': dh._as_quantity(0.25, 'GHz')},
        {'type': 'resonator', 'name': 'R1', 'nodes': ['res1_pad'],
         'capacitance': dh._as_quantity(100, 'fF'), 'omega_r': (2 * np.pi * dh._as_quantity(6, 'GHz')).to('1/s')},
        {'type': 'resonator', 'name': 'R2', 'nodes': ['res2_pad'],
         'capacitance': dh._as_quantity(100, 'fF'), 'omega_r': (2 * np.pi * dh._as_quantity(7, 'GHz')).to('1/s')},
        {'type': 'ground', 'name': 'GND', 'nodes': ['ground']},
    ]
    chip = Chip.from_simulation(dh, C, node_labels, specs)
    cm = chip.coupling_matrix()
    assert len(cm) == 3
    assert all('GND' not in pair for pair in cm)
    print(f"PASS: coupling_matrix has {len(cm)} pairs, ground excluded")



# ===========================================================================
# New tests: dispersive shift, dressed frequencies, EPR bridge, generate_output

# ===========================================================================

def test_dispersive_shift_matches_formula():
    """chi = -g^2 * alpha / (Delta*(Delta+alpha)), independent manual calc."""
    from scipy.constants import h
    dh = DesignHelper(light_speed='3e8 m/s')
    chip = make_toy_chip(dh)

    g = chip.compute_coupling('Q1', 'R1').to('Hz').magnitude
    # resonator linear freq from omega_r
    omega_r = (2 * np.pi * dh._as_quantity(6, 'GHz')).to('1/s').magnitude
    f_res = omega_r / (2 * np.pi)
    # qubit freq: sqrt(8*EJ*EC) - EC
    EJ_Hz = 15e9; EC_Hz = 0.25e9
    f_qubit = (np.sqrt(8 * EJ_Hz * EC_Hz) - EC_Hz)
    alpha = -EC_Hz  # leading-order anharmonicity

    Delta = f_res - f_qubit
    chi_manual = -(g**2 * alpha) / (Delta * (Delta + alpha))

    chi = chip.compute_dispersive_shift('Q1', 'R1').to('Hz').magnitude
    assert abs(chi - chi_manual) / abs(chi_manual) < 1e-6, (chi, chi_manual)
    print(f"PASS: chi = {chi/1e6:.4f} MHz matches manual {chi_manual/1e6:.4f} MHz")


def test_dressed_frequencies_keys_and_signs():
    """dressed_frequencies returns Qubit↔Resonator pairs only, with correct sign."""
    dh = DesignHelper(light_speed='3e8 m/s')
    chip = make_toy_chip(dh)
    df = chip.dressed_frequencies()

    assert len(df) == 1
    assert ('Q1', 'R1') in df
    entry = df[('Q1', 'R1')]
    assert 'chi' in entry and 'qubit' in entry and 'resonator' in entry

    chi = entry['chi'].to('Hz').magnitude
    f_q_bare = chip.components['Q1'].frequency().to('Hz').magnitude
    f_r_bare = (chip.components['R1'].omega_r / (2 * np.pi)).to('Hz').magnitude

    # qubit dressed = f_q - chi, resonator dressed = f_r + chi
    assert abs(entry['qubit'].to('Hz').magnitude - (f_q_bare - chi)) < 1, \
        "dressed qubit frequency wrong"
    assert abs(entry['resonator'].to('Hz').magnitude - (f_r_bare + chi)) < 1, \
        "dressed resonator frequency wrong"
    print(f"PASS: dressed_frequencies: chi={chi/1e6:.4f} MHz, "
          f"f_q_dressed={entry['qubit'].to('GHz').magnitude:.4f} GHz, "
          f"f_r_dressed={entry['resonator'].to('GHz').magnitude:.4f} GHz")


def test_epr_bridge_coupling_consistent():
    """from_epr_params: compute_coupling on the bridged Chip reproduces g_epr."""
    dh = DesignHelper(light_speed='3e8 m/s')
    from scipy.constants import hbar as hbar_si, e
    # phi_zpf_sq for a 6 GHz resonator at 50 Ohm: Z=50, phi_zpf^2=hbar*Z/2, in phi0^2
    phi0 = hbar_si / (2 * e)
    Z0 = 50.0
    phi_zpf_sq_dimless = (hbar_si * Z0 / 2) / phi0**2

    epr_params = {
        "Ej": dh._as_quantity(15, 'GHz'),
        "Ec": dh._as_quantity(0.25, 'GHz'),
        "g": dh._as_quantity(63.234, 'MHz'),
        "cavity_frequency_linear": dh._as_quantity(6, 'GHz'),
        "phi_zpf_sq": phi_zpf_sq_dimless,
    }
    chip_epr = Chip.from_epr_params(dh, epr_params, qubit_name='Q1', resonator_name='R1')
    g_roundtrip = chip_epr.compute_coupling('Q1', 'R1').to('MHz').magnitude
    assert abs(g_roundtrip - 63.234) < 1e-3, \
        f"EPR bridge g roundtrip failed: {g_roundtrip:.4f} MHz vs 63.234 MHz"
    print(f"PASS: from_epr_params roundtrip g = {g_roundtrip:.4f} MHz")


def test_generate_output_structure():
    """generate_output returns expected keys and non-empty component table."""
    dh = DesignHelper(light_speed='3e8 m/s')
    chip = make_toy_chip(dh)
    out = chip.generate_output(title="Test Chip")

    assert out['title'] == "Test Chip"
    assert 'components' in out and 'couplings' in out and 'dispersive' in out

    # check component rows include Q1 and R1
    try:
        import pandas as pd
        names = list(out['components']['name'])
    except ImportError:
        names = [r['name'] for r in out['components']]
    assert 'Q1' in names and 'R1' in names

    # couplings has at least one entry
    try:
        import pandas as pd
        assert not out['couplings'].empty
    except ImportError:
        assert len(out['couplings']) > 0

    print("PASS: generate_output structure correct")
    chip.print_output(title="Test Chip")



# ===========================================================================
# New tests: dispersive guardrail, Res-Res hybridization, Purcell, higher-order
# ===========================================================================

def test_dispersive_guardrail_valid():
    """is_dispersive_valid=True when g/|Delta| << 0.1."""
    dh = DesignHelper(light_speed='3e8 m/s')
    chip = make_toy_chip(dh)
    df = chip.dressed_frequencies()
    entry = df[('Q1', 'R1')]
    assert 'is_dispersive_valid' in entry
    assert entry['is_dispersive_valid'] == True, \
        f"Expected valid dispersive regime, ratio={entry['dispersive_ratio']:.4f}"
    print(f"PASS: dispersive guardrail valid, ratio={entry['dispersive_ratio']:.4f}")


def test_dispersive_guardrail_invalid():
    """is_dispersive_valid=False when g/|Delta| > 0.1 (near-resonant case)."""
    dh = DesignHelper(light_speed='3e8 m/s')
    # Put qubit very close to resonator so g/Delta >> 0.1
    cap_matrix_F = np.array([[100e-15, -30e-15], [-30e-15, 100e-15]])
    node_labels = ['qubit_pad', 'res_pad']
    specs = [
        {'type': 'qubit', 'name': 'Q1', 'nodes': ['qubit_pad'],
         'EJ': dh._as_quantity(15, 'GHz'), 'EC': dh._as_quantity(0.25, 'GHz')},
        {'type': 'resonator', 'name': 'R1', 'nodes': ['res_pad'],
         'capacitance': dh._as_quantity(100, 'fF'),
         # Resonator very close to qubit (~5.23 GHz) -> small Delta, large g
         'omega_r': (2 * np.pi * dh._as_quantity(5.3, 'GHz')).to('1/s')},
    ]
    chip = Chip.from_simulation(dh, cap_matrix_F, node_labels, specs)
    df = chip.dressed_frequencies()
    entry = df[('Q1', 'R1')]
    assert entry['is_dispersive_valid'] == False, \
        f"Expected invalid dispersive regime, ratio={entry['dispersive_ratio']:.4f}"
    print(f"PASS: dispersive guardrail triggered, ratio={entry['dispersive_ratio']:.4f}")


def test_resonator_resonator_hybridization():
    """Res↔Res: exact normal-mode frequencies match analytic formula."""
    dh = DesignHelper(light_speed='3e8 m/s')
    C = np.array([
        [100e-15, -5e-15, -5e-15],
        [-5e-15, 100e-15, -5e-15],
        [-5e-15, -5e-15, 100e-15],
    ])
    node_labels = ['qubit_pad', 'res1_pad', 'res2_pad']
    specs = [
        {'type': 'qubit', 'name': 'Q1', 'nodes': ['qubit_pad'],
         'EJ': dh._as_quantity(15, 'GHz'), 'EC': dh._as_quantity(0.25, 'GHz')},
        {'type': 'resonator', 'name': 'R1', 'nodes': ['res1_pad'],
         'capacitance': dh._as_quantity(100, 'fF'),
         'omega_r': (2 * np.pi * dh._as_quantity(6.0, 'GHz')).to('1/s')},
        {'type': 'resonator', 'name': 'R2', 'nodes': ['res2_pad'],
         'capacitance': dh._as_quantity(100, 'fF'),
         'omega_r': (2 * np.pi * dh._as_quantity(7.0, 'GHz')).to('1/s')},
    ]
    chip = Chip.from_simulation(dh, C, node_labels, specs)
    df = chip.dressed_frequencies()
    assert ('R1', 'R2') in df, f"R1↔R2 pair missing, got keys: {list(df.keys())}"
    entry = df[('R1', 'R2')]
    assert 'omega_plus' in entry and 'omega_minus' in entry

    # Manual check: omega_pm = (6+7)/2 +/- sqrt((1/2)^2 + 4*g^2)/2
    g = chip.compute_coupling('R1', 'R2').to('GHz').magnitude
    f1, f2 = 6.0, 7.0
    disc = np.sqrt((f2 - f1)**2 + 4*g**2) / 2
    omega_plus_manual = (f1 + f2)/2 + disc
    omega_minus_manual = (f1 + f2)/2 - disc

    op = entry['omega_plus'].to('GHz').magnitude
    om = entry['omega_minus'].to('GHz').magnitude
    assert abs(op - omega_plus_manual) < 1e-6, f"omega_plus mismatch: {op} vs {omega_plus_manual}"
    assert abs(om - omega_minus_manual) < 1e-6, f"omega_minus mismatch: {om} vs {omega_minus_manual}"
    print(f"PASS: R↔R hybridization: omega_+ = {op:.4f} GHz, omega_- = {om:.4f} GHz, "
          f"mixing_angle = {np.degrees(entry['mixing_angle']):.2f}°")


def test_purcell_decay():
    """Purcell: gamma = kappa*(g/Delta)^2, None when kappa unset."""
    dh = DesignHelper(light_speed='3e8 m/s')
    # Without kappa -- should return None
    chip = make_toy_chip(dh)
    assert chip.compute_purcell('Q1', 'R1') is None
    print("PASS: compute_purcell returns None when kappa not set")

    # Now set kappa and check formula
    chip.update_component('R1', kappa=dh._as_quantity(1.0, 'MHz'))
    result = chip.compute_purcell('Q1', 'R1')
    assert result is not None

    # Manual: g/Delta
    g_Hz = chip.compute_coupling('Q1', 'R1').to('Hz').magnitude
    comp_r = chip.components['R1']
    f_res = (comp_r.omega_r / (2 * np.pi)).to('Hz').magnitude
    f_q = chip.components['Q1'].frequency().to('Hz').magnitude
    Delta = f_res - f_q
    gamma_manual = 1e6 * (g_Hz / Delta)**2  # kappa=1MHz * ratio^2, in Hz
    gamma_computed = result['gamma_purcell'].to('Hz').magnitude
    assert abs(gamma_computed - gamma_manual) / abs(gamma_manual) < 1e-6, \
        f"Purcell gamma mismatch: {gamma_computed:.3e} vs {gamma_manual:.3e}"
    print(f"PASS: Purcell gamma = {gamma_computed/1e3:.4f} kHz, "
          f"T1_purcell = {result['T1_purcell_us']:.2f} us")


def test_higher_order_dispersive():
    """chi_prime and vacuum shifts have correct signs and scaling."""
    dh = DesignHelper(light_speed='3e8 m/s')
    chip = make_toy_chip(dh)
    ho = chip.compute_higher_order_dispersive('Q1', 'R1')

    assert 'chi' in ho and 'chi_prime' in ho
    assert 'vacuum_shift_res' in ho and 'vacuum_shift_q' in ho
    assert 'is_dispersive_valid' in ho

    # vacuum_shift_res should be negative (resonator pulled down), qubit positive
    assert ho['vacuum_shift_res'].to('Hz').magnitude < 0, "vac shift res should be negative"
    assert ho['vacuum_shift_q'].to('Hz').magnitude > 0, "vac shift q should be positive"

    # chi_prime: for transmon (alpha < 0, Delta > 0): chi < 0, chi' should be positive
    # chi' = chi * 2*alpha / (Delta + alpha); alpha=-EC<0, Delta>0 => 2*alpha<0, Delta+alpha>0
    # => chi' = chi * (negative/positive) = chi * negative; chi<0 => chi'>0
    chi_prime_Hz = ho['chi_prime'].to('Hz').magnitude
    chi_Hz = ho['chi'].to('Hz').magnitude
    assert chi_prime_Hz * chi_Hz < 0, \
        f"chi_prime sign unexpected: chi={chi_Hz:.3e}, chi'={chi_prime_Hz:.3e}"
    print(f"PASS: higher-order dispersive: chi={chi_Hz/1e6:.4f} MHz, "
          f"chi'={chi_prime_Hz/1e3:.4f} kHz, "
          f"vac_shift_res={ho['vacuum_shift_res'].to('Hz').magnitude/1e3:.4f} kHz")


def test_inductance_matrix_inverse_matches_manual_calc():
    """linv is the matrix inverse of inductance_matrix, independent check."""
    dh = DesignHelper(light_speed='3e8 m/s')
    chip = make_toy_chip(dh)  # capacitance matrix + components, unrelated

    L = np.array([
        [12.0, 0.8, 0.3],
        [0.8, 15.0, 0.5],
        [0.3, 0.5, 9.0],
    ])  # nH
    chip.set_inductance_matrix(L, node_labels=['loop_a', 'loop_b', 'loop_c'])

    linv = chip.linv.to('1/nH').magnitude
    linv_manual = np.linalg.inv(L)  # fully independent, no Chip code reused

    assert np.allclose(linv, linv_manual, rtol=1e-9), (linv, linv_manual)

    m_ab = chip.inductance_element('loop_a', 'loop_b').to('nH').magnitude
    assert abs(m_ab - 0.8) < 1e-9
    print(f"PASS: inductance matrix inverse matches manual np.linalg.inv, "
          f"M[a,b]={m_ab:.4f} nH")


def test_inductance_node_labels_independent_of_capacitance_node_labels():
    """Inductance-matrix node labels are a separate index space from the
    capacitance matrix's node_labels -- setting one must not disturb the
    other, and they're allowed to differ in count/names/order entirely.
    """
    dh = DesignHelper(light_speed='3e8 m/s')
    chip = make_toy_chip(dh)  # node_labels = ['qubit_pad', 'res_pad']
    cap_node_labels_before = list(chip.node_labels)

    # Inductance matrix over a totally different, larger set of loops
    L = np.eye(4) * 10.0
    L[0, 1] = L[1, 0] = 1.0
    chip.set_inductance_matrix(L, node_labels=['jj_arm_1', 'jj_arm_2', 'loop_x', 'loop_y'])

    assert chip.node_labels == cap_node_labels_before, \
        "capacitance node_labels were disturbed by set_inductance_matrix"
    assert chip.inductance_node_labels == ['jj_arm_1', 'jj_arm_2', 'loop_x', 'loop_y']
    assert chip.capacitance_matrix is not None, \
        "capacitance matrix was wiped by set_inductance_matrix"
    print("PASS: inductance/capacitance node-label spaces are independent")


def test_inductance_matrix_size_mismatch_raises():
    dh = DesignHelper(light_speed='3e8 m/s')
    chip = make_toy_chip(dh)
    try:
        chip.set_inductance_matrix(np.eye(3), node_labels=['a', 'b'])
        assert False, "should have raised"
    except ValueError:
        print("PASS: inductance matrix/node_labels size mismatch raises ValueError")


def test_inductance_unknown_label_raises():
    dh = DesignHelper(light_speed='3e8 m/s')
    chip = make_toy_chip(dh)
    chip.set_inductance_matrix(np.eye(2) * 10.0, node_labels=['a', 'b'])
    try:
        chip.inductance_element('a', 'nope')
        assert False, "should have raised"
    except ValueError:
        print("PASS: unknown inductance node label raises ValueError")


def test_inductance_matrix_not_loaded_raises():
    dh = DesignHelper(light_speed='3e8 m/s')
    chip = make_toy_chip(dh)
    try:
        _ = chip.linv
        assert False, "should have raised"
    except RuntimeError:
        print("PASS: accessing linv before set_inductance_matrix raises RuntimeError")


def test_from_inductance_simulation_standalone():
    """from_inductance_simulation builds a Chip with no capacitance matrix
    or components -- inductance-matrix-only ingestion, as requested.
    """
    dh = DesignHelper(light_speed='3e8 m/s')
    L = np.array([[20.0, 2.0], [2.0, 25.0]])
    chip = Chip.from_inductance_simulation(dh, L, node_labels=['squid_a', 'squid_b'])

    assert chip.capacitance_matrix is None
    assert len(chip.components) == 0
    m = chip.inductance_element('squid_a', 'squid_b').to('nH').magnitude
    assert abs(m - 2.0) < 1e-9
    print(f"PASS: from_inductance_simulation standalone Chip, M[a,b]={m:.4f} nH")


def test_floating_squid_EJ_matches_compute_squid_EJ():
    """FloatingSQUID.EJ tracks DesignHelper.compute_squid_EJ exactly as
    flux changes -- independent check, not reusing FloatingSQUID's own
    property in the assertion.
    """
    from components import FloatingSQUID
    dh = DesignHelper(light_speed='3e8 m/s')
    sq = FloatingSQUID('Qsq', [0, 1], dh, EJ1=dh._as_quantity(10, 'GHz'),
                        EJ2=dh._as_quantity(8, 'GHz'), EC=dh._as_quantity(0.25, 'GHz'),
                        flux=0.0)

    for flux in [0.0, 0.25, 0.4]:
        sq.flux = flux
        EJ_via_property = sq.EJ.to('GHz').magnitude
        EJ_independent = dh.compute_squid_EJ(
            dh._as_quantity(10, 'GHz'), dh._as_quantity(8, 'GHz'), flux
        ).to('GHz').magnitude
        assert abs(EJ_via_property - EJ_independent) < 1e-9, (flux, EJ_via_property, EJ_independent)
    print(f"PASS: FloatingSQUID.EJ tracks compute_squid_EJ across flux sweep "
          f"(EJ(0)={sq.EJ.to('GHz').magnitude:.4f} GHz at flux={sq.flux})")


def test_floating_squid_frequency_changes_with_flux():
    """Inherited Qubit.frequency()/.anharmonicity() pick up the flux-tuned
    EJ with zero overrides -- this is the actual point of making EJ a
    property instead of duplicating Qubit's methods.
    """
    from components import FloatingSQUID
    dh = DesignHelper(light_speed='3e8 m/s')
    sq = FloatingSQUID('Qsq', [0, 1], dh, EJ1=dh._as_quantity(10, 'GHz'),
                        EJ2=dh._as_quantity(8, 'GHz'), EC=dh._as_quantity(0.25, 'GHz'),
                        flux=0.0)
    f_at_zero_flux = sq.frequency().to('GHz').magnitude
    sq.flux = 0.4
    f_at_flux = sq.frequency().to('GHz').magnitude
    assert f_at_zero_flux != f_at_flux, "frequency should change with flux"

    # independent cross-check via compute_flux_tunable_frequency directly
    f_manual = dh.compute_flux_tunable_frequency(
        dh._as_quantity(10, 'GHz'), dh._as_quantity(8, 'GHz'),
        dh._as_quantity(0.25, 'GHz'), 0.4
    ).to('GHz').magnitude
    assert abs(f_at_flux - f_manual) < 1e-9, (f_at_flux, f_manual)
    print(f"PASS: FloatingSQUID frequency at flux=0 -> {f_at_zero_flux:.4f} GHz, "
          f"at flux=0.4 -> {f_at_flux:.4f} GHz (matches compute_flux_tunable_frequency)")


def test_floating_squid_coupling_matches_manual_calc():
    """FloatingSQUID combined with the differential-mode Chip coupling
    convention: independent manual calc, same pattern as
    test_floating_transmon_coupling_matches_manual_calc but with
    flux-tunable EJ on top.
    """
    dh = DesignHelper(light_speed='3e8 m/s')
    C = np.array([
        [100e-15, -40e-15, -8e-15],
        [-40e-15, 100e-15, -3e-15],
        [-8e-15,  -3e-15, 100e-15],
    ])
    chip = Chip(dh, ['sq_pad0', 'sq_pad1', 'res_pad'])
    chip.set_capacitance_matrix(C)
    chip.add_component({
        'type': 'floating_squid', 'name': 'Qsq', 'nodes': ['sq_pad0', 'sq_pad1'],
        'EJ1': dh._as_quantity(10, 'GHz'), 'EJ2': dh._as_quantity(8, 'GHz'),
        'EC': dh._as_quantity(0.25, 'GHz'), 'flux': 0.3,
    })
    chip.add_component({
        'type': 'resonator', 'name': 'R1', 'nodes': ['res_pad'],
        'capacitance': dh._as_quantity(100, 'fF'),
        'omega_r': (2 * np.pi * dh._as_quantity(6, 'GHz')).to('1/s'),
    })

    g = chip.compute_coupling('Qsq', 'R1').to('Hz').magnitude

    from scipy.constants import e, h, hbar
    Cinv = np.linalg.inv(C)
    cinv_eff = Cinv[0, 2] - Cinv[1, 2]
    EJ_flux_GHz = dh.compute_squid_EJ(
        dh._as_quantity(10, 'GHz'), dh._as_quantity(8, 'GHz'), 0.3
    ).to('GHz').magnitude
    EJ_J, EC_J = EJ_flux_GHz * 1e9 * h, 0.25e9 * h
    n_zpf = 0.5 * (EJ_J / (2 * EC_J)) ** 0.25
    q_zpf_qubit = 2 * e * n_zpf
    omega_res = 2 * np.pi * 6e9
    Z_res = 1 / (omega_res * 100e-15)
    q_zpf_res = np.sqrt(hbar / (2 * Z_res))
    g_manual_Hz = (0.5 * cinv_eff * q_zpf_qubit * q_zpf_res) / h

    assert abs(g - g_manual_Hz) / abs(g_manual_Hz) < 1e-9, (g, g_manual_Hz)
    print(f"PASS: FloatingSQUID coupling g = {g/1e6:.4f} MHz matches manual calc "
          f"{g_manual_Hz/1e6:.4f} MHz (flux=0.3, EJ(flux)={EJ_flux_GHz:.4f} GHz)")


def test_floating_squid_requires_two_nodes():
    from components import FloatingSQUID
    dh = DesignHelper(light_speed='3e8 m/s')
    try:
        FloatingSQUID('bad', [0], dh, EJ1=dh._as_quantity(10, 'GHz'),
                       EJ2=dh._as_quantity(8, 'GHz'), EC=dh._as_quantity(0.25, 'GHz'))
        assert False, "should have raised"
    except ValueError:
        print("PASS: FloatingSQUID with 1 node raises ValueError")
    try:
        FloatingSQUID('bad', [0, 1, 2], dh, EJ1=dh._as_quantity(10, 'GHz'),
                       EJ2=dh._as_quantity(8, 'GHz'), EC=dh._as_quantity(0.25, 'GHz'))
        assert False, "should have raised"
    except ValueError:
        print("PASS: FloatingSQUID with 3 nodes raises ValueError")


def test_floating_squid_EJ_cannot_be_set_directly():
    from components import FloatingSQUID
    dh = DesignHelper(light_speed='3e8 m/s')
    sq = FloatingSQUID('Qsq', [0, 1], dh, EJ1=dh._as_quantity(10, 'GHz'),
                        EJ2=dh._as_quantity(8, 'GHz'), EC=dh._as_quantity(0.25, 'GHz'))
    try:
        sq.EJ = dh._as_quantity(15, 'GHz')
        assert False, "should have raised"
    except AttributeError:
        print("PASS: FloatingSQUID.EJ direct assignment raises AttributeError")


def test_floating_squid_update_component_flux():
    """Chip.update_component('Qsq', flux=...) is the intended way to
    flux-tune a FloatingSQUID already registered on a Chip.
    """
    dh = DesignHelper(light_speed='3e8 m/s')
    chip = Chip(dh, ['sq_pad0', 'sq_pad1', 'res_pad'])
    chip.set_capacitance_matrix(np.array([
        [100e-15, -40e-15, -8e-15],
        [-40e-15, 100e-15, -3e-15],
        [-8e-15,  -3e-15, 100e-15],
    ]))
    chip.add_component({
        'type': 'floating_squid', 'name': 'Qsq', 'nodes': ['sq_pad0', 'sq_pad1'],
        'EJ1': dh._as_quantity(10, 'GHz'), 'EJ2': dh._as_quantity(8, 'GHz'),
        'EC': dh._as_quantity(0.25, 'GHz'), 'flux': 0.0,
    })
    chip.add_component({
        'type': 'resonator', 'name': 'R1', 'nodes': ['res_pad'],
        'capacitance': dh._as_quantity(100, 'fF'),
        'omega_r': (2 * np.pi * dh._as_quantity(6, 'GHz')).to('1/s'),
    })
    g_before = chip.compute_coupling('Qsq', 'R1').to('MHz').magnitude
    chip.update_component('Qsq', flux=0.4)
    g_after = chip.compute_coupling('Qsq', 'R1').to('MHz').magnitude
    assert g_before != g_after
    print(f"PASS: update_component(flux=...) changed g from {g_before:.4f} to {g_after:.4f} MHz")


if __name__ == "__main__":
    test_coupling_matches_manual_calc()
    test_duplicate_component_raises()
    test_unknown_node_label_raises()
    test_update_component_changes_live_coupling()
    test_remove_component()
    test_matrix_size_mismatch_raises()
    test_floating_transmon_coupling_matches_manual_calc()
    test_floating_transmon_coupling_symmetric_cancels()
    test_coupling_matrix_excludes_ground()
    test_dispersive_shift_matches_formula()
    test_dressed_frequencies_keys_and_signs()
    test_epr_bridge_coupling_consistent()
    test_generate_output_structure()
    test_dispersive_guardrail_valid()
    test_dispersive_guardrail_invalid()
    test_resonator_resonator_hybridization()
    test_purcell_decay()
    test_higher_order_dispersive()
    test_inductance_matrix_inverse_matches_manual_calc()
    test_inductance_node_labels_independent_of_capacitance_node_labels()
    test_inductance_matrix_size_mismatch_raises()
    test_inductance_unknown_label_raises()
    test_inductance_matrix_not_loaded_raises()
    test_from_inductance_simulation_standalone()
    test_floating_squid_EJ_matches_compute_squid_EJ()
    test_floating_squid_frequency_changes_with_flux()
    test_floating_squid_coupling_matches_manual_calc()
    test_floating_squid_requires_two_nodes()
    test_floating_squid_EJ_cannot_be_set_directly()
    test_floating_squid_update_component_flux()
    print("\nAll tests passed.")