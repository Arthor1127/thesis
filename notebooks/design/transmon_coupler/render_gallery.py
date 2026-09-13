"""Parameter-sweep gallery for the two-transmon interdigital coupler.

Builds the CircTransmon + CapNInterdigitalRound design over several sweeps of the
coupling geometry (qubit spacing, tangent-flare radius, moat fillets, finger
count / length, ...) and writes a 2-panel image (full design + junction zoom)
for every variant into ``gallery/``, plus a contact sheet per sweep and a
gallery/README.md.

Run:  conda run -n qcg-quantum-design python render_gallery.py
"""
import os
import sys
from dataclasses import dataclass, replace, asdict

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from qiskit_metal import designs, Dict
from qiskit_metal.qlibrary.qubits.circle_transmon import CircTransmon
from qiskit_metal.qlibrary.couplers.cap_n_interdigital_round import CapNInterdigitalRound

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from render_geometry import render, _design_bounds_mm  # noqa: E402

OUT = os.path.join(HERE, "gallery")

# CHIP_MARGIN_FACTOR analogue, kept small so the full panel stays readable
CHIP_MARGIN = 2.4
CHIP_Z_UM = 500.0
ABUT_OVERLAP_UM = 4.0


@dataclass(frozen=True)
class Cfg:
    # transmon pad
    pad_radius_um: float = 90.0
    pad_gap_um: float = 10.0
    # coupler fork
    finger_count: int = 8
    finger_width_um: float = 8.0
    finger_length_um: float = 80.0
    finger_gap_um: float = 4.0
    finger_tip_gap_um: float = 6.0
    spine_width_um: float = 12.0
    fillet_radius_um: float = 6.0
    round_resolution: int = 32
    # coupler ground moat
    cap_gap_ground_um: float = 6.0
    etch_corner_radius_um: float = 16.0
    etch_blend_radius_um: float = 18.0
    # neck / tangent flare into the pad
    lead_length_um: float = 40.0
    lead_width_um: float = 12.0
    lead_flare_width_um: float = 44.0
    lead_blend_radius_um: float = 40.0

    # --- derived ---
    @property
    def half_h_um(self):
        return 0.5 * (self.finger_length_um + self.finger_tip_gap_um)

    @property
    def island_center_offset_um(self):
        span = self.half_h_um + self.spine_width_um + self.lead_length_um
        return self.pad_radius_um + span - ABUT_OVERLAP_UM


def build(cfg: Cfg):
    d = designs.DesignPlanar({}, overwrite_enabled=True)
    d.variables["cpw_width"] = "20um"
    d.variables["cpw_gap"] = "12.25um"

    D = cfg.island_center_offset_um

    def transmon(pos_x_um):
        return Dict(
            pos_x=f"{pos_x_um}um", pos_y="0um", orientation="0", chip="main",
            cpw_width="cpw_width", cpw_gap="cpw_gap",
            pad_radius=f"{cfg.pad_radius_um}um", pad_gap=f"{cfg.pad_gap_um}um",
            jj_options=Dict(jj_width="5um", jj_angle="270", L_j="10nH",
                            C_j="2fF", export_mask=False, jj_sim_gap="3um"),
        )

    CircTransmon(d, "Q1", options=transmon(-D))
    CircTransmon(d, "Q2", options=transmon(+D))
    CapNInterdigitalRound(d, "coupler", options=Dict(
        pos_x="0um", pos_y="0um", orientation="90", chip="main",
        finger_count=f"{cfg.finger_count}",
        finger_width=f"{cfg.finger_width_um}um",
        finger_length=f"{cfg.finger_length_um}um",
        finger_gap=f"{cfg.finger_gap_um}um",
        finger_tip_gap=f"{cfg.finger_tip_gap_um}um",
        spine_width=f"{cfg.spine_width_um}um",
        fillet_radius=f"{cfg.fillet_radius_um}um",
        cap_gap_ground=f"{cfg.cap_gap_ground_um}um",
        etch_corner_radius=f"{cfg.etch_corner_radius_um}um",
        etch_blend_radius=f"{cfg.etch_blend_radius_um}um",
        lead_length=f"{cfg.lead_length_um}um",
        lead_width=f"{cfg.lead_width_um}um",
        lead_flare_width=f"{cfg.lead_flare_width_um}um",
        lead_blend_radius=f"{cfg.lead_blend_radius_um}um",
        round_resolution=f"{cfg.round_resolution}",
        island_radius=f"{cfg.pad_radius_um}um",
        island_gap=f"{cfg.pad_gap_um}um",
        island_center_offset=f"{cfg.island_center_offset_um}um",
    ))

    d.rebuild()

    minx, miny, maxx, maxy = _design_bounds_mm(d)
    cx, cy = (minx + maxx) / 2, (miny + maxy) / 2
    d.chips.main.size.size_x = f"{(maxx - minx) * CHIP_MARGIN}mm"
    d.chips.main.size.size_y = f"{(maxy - miny) * CHIP_MARGIN}mm"
    d.chips.main.size.size_z = f"{CHIP_Z_UM}um"
    d.chips.main.size.center_x = f"{cx}mm"
    d.chips.main.size.center_y = f"{cy}mm"
    return d


def _sanity(d):
    """(#metal nets, min gap um, #ground pieces) -- cheap geometry health check."""
    from shapely.ops import unary_union
    from shapely.geometry import box as _box

    p = d.qgeometry.tables["poly"]
    sub = p["subtract"].astype(bool)
    metal = unary_union(list(p[~sub]["geometry"]))
    big = sorted((g for g in getattr(metal, "geoms", [metal]) if g.area * 1e6 > 500),
                 key=lambda g: -g.area)
    gap = min((a.distance(b) for i, a in enumerate(big) for b in big[i + 1:]),
              default=float("nan"))
    sx = float(d.chips.main.size.size_x[:-2]); sy = float(d.chips.main.size.size_y[:-2])
    cx = float(d.chips.main.size.center_x[:-2]); cy = float(d.chips.main.size.center_y[:-2])
    ground = _box(cx - sx / 2, cy - sy / 2, cx + sx / 2, cy + sy / 2).difference(
        unary_union(list(p[sub]["geometry"])))
    gp = [g for g in getattr(ground, "geoms", [ground]) if g.area > 1e-7]
    return len(big), gap * 1e3, len(gp)


def figure(cfg: Cfg, title: str, path: str):
    d = build(cfg)
    nets, gap_um, gpieces = _sanity(d)
    minx, miny, maxx, maxy = _design_bounds_mm(d)
    D = cfg.island_center_offset_um * 1e-3
    pad_edge = -(D - cfg.pad_radius_um * 1e-3)

    fig, ax = plt.subplots(1, 2, figsize=(15, 6), gridspec_kw={"width_ratios": [1.55, 1]})
    render(d, ax[0])
    ax[0].set_xlim(minx - 0.03, maxx + 0.03)
    ax[0].set_ylim(miny - 0.03, maxy + 0.03)
    ax[0].set_title("full design")

    render(d, ax[1])
    zx = 0.070
    ax[1].set_xlim(pad_edge - 0.030, pad_edge + zx)
    ax[1].set_ylim(-0.5 * (zx + 0.030), 0.5 * (zx + 0.030))
    ax[1].set_title("Q1 junction")
    ax[1].grid(alpha=0.15)

    flag = "" if (nets == 2 and gpieces == 1 and gap_um > 1.5) else "   [!] "
    fig.suptitle(f"{flag}{title}\n"
                 f"qubit spacing 2D = {2 * cfg.island_center_offset_um:.0f} um   |   "
                 f"metal nets {nets}   min gap {gap_um:.1f} um   ground pieces {gpieces}",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return dict(nets=nets, gap_um=round(gap_um, 1), gpieces=gpieces,
                spacing_um=round(2 * cfg.island_center_offset_um))


BASE = Cfg()

# (dir name, human title, knob attr, "{value} label" fmt, [values])
SWEEPS = [
    ("01_qubit_spacing", "Qubit spacing (lead_length: spine-to-pad gap)",
     "lead_length_um", "lead_length = {v:g} um", [14, 24, 40, 60, 90]),
    ("02_tangent_flare_radius", "Tangent-flare radius where the neck leaves the pad",
     "lead_blend_radius_um", "lead_blend_radius = {v:g} um", [10, 20, 40, 60]),
    ("03_neck_width", "Neck width at the pad edge",
     "lead_width_um", "lead_width = {v:g} um", [6, 12, 22, 36]),
    ("04_moat_corner_radius", "Ground-moat rectangle corner radius",
     "etch_corner_radius_um", "etch_corner_radius = {v:g} um", [3, 12, 24, 44]),
    ("05_moat_pocket_cusp_fillet", "Fillet where the moat meets each transmon pocket",
     "etch_blend_radius_um", "etch_blend_radius = {v:g} um", [2, 10, 20, 35]),
    ("06_moat_width", "Ground-moat width (dielectric gap to ground)",
     "cap_gap_ground_um", "cap_gap_ground = {v:g} um", [3, 6, 12, 20]),
    ("07_finger_root_fillet", "Finger-root / spine-corner fillet radius",
     "fillet_radius_um", "fillet_radius = {v:g} um", [0, 3, 6, 10]),
    ("08_finger_count", "Number of interdigital fingers (coupling strength)",
     "finger_count", "finger_count = {v}", [4, 6, 8, 12, 16]),
    ("09_finger_length", "Finger length / overlap (coupling strength)",
     "finger_length_um", "finger_length = {v:g} um", [40, 60, 80, 110, 150]),
    ("10_finger_gap", "Gap between interleaved fingers (coupling strength)",
     "finger_gap_um", "finger_gap = {v:g} um", [3, 4, 6, 10]),
    ("11_pad_radius", "Transmon pad radius",
     "pad_radius_um", "pad_radius = {v:g} um", [60, 90, 120]),
]


def main():
    os.makedirs(OUT, exist_ok=True)

    figure(BASE, "Baseline configuration", os.path.join(OUT, "00_baseline.png"))

    lines = ["# Interdigital transmon-coupler geometry gallery", "",
             "Regenerate with `conda run -n qcg-quantum-design python render_gallery.py`.",
             "",
             "Every panel: **left** = full two-transmon design, **right** = zoom on the "
             "Q1 junction. Blue = metal, red hatch = ground-plane cutout. A `[!]` in a "
             "title flags a variant whose quick geometry check is off "
             "(want: 2 metal nets, 1 connected ground piece, comb-to-comb gap > "
             "`finger_gap`).", "",
             f"**Baseline** ([00_baseline.png](00_baseline.png)): pad_radius "
             f"{BASE.pad_radius_um:g} um, pad_gap {BASE.pad_gap_um:g} um, "
             f"finger_count {BASE.finger_count}, finger_length {BASE.finger_length_um:g} um, "
             f"finger_gap {BASE.finger_gap_um:g} um, spine_width {BASE.spine_width_um:g} um, "
             f"lead_length {BASE.lead_length_um:g} um, lead_blend_radius "
             f"{BASE.lead_blend_radius_um:g} um.", "",
             "Each sweep below also has a stacked contact sheet "
             "(`<sweep>.png`).", ""]

    for dname, title, attr, labelfmt, values in SWEEPS:
        sweep_dir = os.path.join(OUT, dname)
        os.makedirs(sweep_dir, exist_ok=True)
        print(f"\n=== {title} ===")
        lines += [f"## {title}", "",
                  f"`{attr}` swept over {values}. Contact sheet: [{dname}.png]({dname}.png)",
                  ""]
        panels = []
        for v in values:
            cfg = replace(BASE, **{attr: v})
            label = labelfmt.format(v=v)
            fname = f"{attr}_{str(v).replace('.', 'p')}.png"
            path = os.path.join(sweep_dir, fname)
            info = figure(cfg, f"{title}  --  {label}", path)
            print(f"  {label:34s} spacing={info['spacing_um']}um  nets={info['nets']} "
                  f"gap={info['gap_um']}um  ground={info['gpieces']}")
            panels.append((label, path, info))
            lines.append(f"- `{label}` -> spacing {info['spacing_um']} um, "
                         f"nets {info['nets']}, min gap {info['gap_um']} um, "
                         f"ground {info['gpieces']}  "
                         f"([{fname}]({dname}/{fname}))")
        lines.append("")

        # contact sheet for the sweep
        n = len(panels)
        cs, axes = plt.subplots(n, 1, figsize=(13, 4.6 * n))
        axes = np.atleast_1d(axes)
        for ax, (label, path, _info) in zip(axes, panels):
            ax.imshow(plt.imread(path))
            ax.set_axis_off()
            ax.set_title(label, fontsize=12)
        cs.suptitle(title, fontsize=14, y=1.0)
        cs.tight_layout()
        cs.savefig(os.path.join(OUT, f"{dname}.png"), dpi=90, bbox_inches="tight")
        plt.close(cs)

    with open(os.path.join(OUT, "README.md"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"\nWrote gallery to {OUT}/  ({sum(len(v) for *_, v in SWEEPS)} panels + "
          f"{len(SWEEPS)} contact sheets + README.md)")


if __name__ == "__main__":
    main()
