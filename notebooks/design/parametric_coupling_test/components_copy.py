"""Component classes for the superconducting design workflow.

Per the locked-in architecture decisions:
- Components are dumb data holders for *bare* physical parameters (EJ/EC,
  omega_r/L/C). They hold a shared DesignHelper instance and call into it
  for derived quantities (ZPFs, frequencies) -- they never reimplement
  physics formulas themselves.
- Components know which capacitance-matrix node index/indices are "theirs",
  but never see the capacitance matrix itself -- that's Chip's job.
- Resonator.zpfs() and Qubit.zpfs() both return (phi_zpf, q_zpf) with
  q_zpf in COULOMBS, even though the transmon's natural ZPF (n_zpf) is
  dimensionless. This is a deliberate contract: Chip.compute_coupling can
  then treat any two components identically without caring whether they're
  a resonator or a qubit. See Qubit.zpfs() for the 2e conversion.
"""

import pint


class Component:
    """Base class. Holds bare physical parameters + node bookkeeping only."""

    def __init__(self, name: str, node_indices, design_helper):
        self.name = name
        # Supports multi-node components (e.g. a SQUID-loop transmon with
        # two junction nodes) from the start, even though single-node is
        # the only case confirmed against real SQDMetal output so far.
        self.node_indices = (
            [node_indices] if isinstance(node_indices, int) else list(node_indices)
        )
        self.dh = design_helper  # shared, not owned -- Component never constructs its own

    def zpfs(self):
        """Return (phi_zpf, q_zpf), q_zpf always in Coulombs. Subclasses
        must override; base class has no physical content to derive this
        from.
        """
        raise NotImplementedError(
            f"{type(self).__name__} does not implement zpfs()"
        )

    def impedance(self):
        """Characteristic impedance of this component's mode (resistance
        units). Subclasses must override; base class has no physical
        content to derive this from. See Resonator.impedance()/
        Qubit.impedance() for the concrete formulas -- this is a genuinely
        different quantity per component type (Z=sqrt(L/C) for a resonator,
        Z_q=sqrt(L_J/C_Sigma) for a qubit), unlike zpfs() which has a
        uniform Coulombs-valued contract across types.
        """
        raise NotImplementedError(
            f"{type(self).__name__} does not implement impedance()"
        )

    def __repr__(self):
        return f"{type(self).__name__}(name={self.name!r}, nodes={self.node_indices})"


class Resonator(Component):
    """An LC (or CPW) resonator mode. Store whichever two of
    (omega_r, inductance, capacitance) are known; compute_resonator_impedance
    resolves the third internally each time zpfs() is called.
    """

    def __init__(self, name, node_indices, design_helper,
                 omega_r=None, inductance=None, capacitance=None, kappa=None):
        super().__init__(name, node_indices, design_helper)
        self.omega_r = omega_r
        self.inductance = inductance
        self.capacitance = capacitance
        self.kappa = kappa  # resonator linewidth (FWHM), in design_helper.freq_units.
        # None means unmeasured/unknown -- Purcell calculation is silently
        # skipped rather than crashing. Set via update_component('R1', kappa=...)
        # once you have the measured or simulated linewidth.

    def zpfs(self):
        phi_zpf, q_zpf = self.dh.compute_resonator_zpfs(
            omega_r=self.omega_r,
            inductance=self.inductance,
            capacitance=self.capacitance,
        )
        return phi_zpf, q_zpf  # already Coulombs -- no conversion needed

    def impedance(self):
        """Characteristic impedance Z = sqrt(L/C) (or the omega_r-derived
        equivalent), via DesignHelper.compute_resonator_impedance -- the
        same resolution-from-any-two-of-three logic zpfs() already uses
        internally, just exposed directly rather than only as an
        intermediate step toward phi_zpf/q_zpf.
        """
        return self.dh.compute_resonator_impedance(
            omega_r=self.omega_r,
            inductance=self.inductance,
            capacitance=self.capacitance,
        )


class Qubit(Component):
    """A transmon-like qubit mode, parametrized by EJ, EC."""

    def __init__(self, name, node_indices, design_helper, EJ=None, EC=None):
        super().__init__(name, node_indices, design_helper)
        self.EJ = EJ
        self.EC = EC

    def zpfs(self):
        phi_zpf, n_zpf = self.dh.compute_transmon_zpfs(self.EJ, self.EC)
        # n_zpf is dimensionless (charge number in units of 2e). Convert to
        # an actual charge ZPF in Coulombs so Chip.compute_coupling can
        # combine it with Cinv (units 1/F) the same way it does for a
        # Resonator's q_zpf, with no special-casing by component type.
        # DO NOT skip this conversion -- a dimensionless number multiplying
        # into Cinv won't raise a pint DimensionalityError, it'll just
        # silently give a wrong magnitude.
        q_zpf = (2.0 * self.dh.electron_charge * n_zpf).to("C")
        return phi_zpf, q_zpf

    def frequency(self, exact=False):
        return self.dh.compute_transmon_frequency(self.EJ, self.EC, exact=exact)

    def anharmonicity(self, exact=False):
        return self.dh.compute_anharmonicity(
            self.EC, exact_EJ=self.EJ if exact else None
        )


class GroundPlane(Component):
    """A node (or set of nodes) tied to ground. No EJ/EC/omega_r -- exists
    so Chip's component-list bookkeeping doesn't need a special case to
    skip ground nodes. zpfs() is intentionally undefined: ground has no
    mode, so asking for its ZPF is a caller error, not a 0/None to paper
    over silently.
    """

    def zpfs(self):
        raise NotImplementedError(
            f"GroundPlane '{self.name}' has no mode -- zpfs() is undefined. "
            "Exclude ground components before calling coupling_matrix()."
        )


class FloatingSQUID(Qubit):
    """A floating (two-island) transmon whose junction is itself a SQUID
    loop -- i.e. flux-tunable EJ, combined with the differential-charge
    multi-node coupling convention. Two node_indices required, same as a
    plain floating transmon (see Chip._cinv_element's differential-mode
    path); NOT the same case as CircTransmonSQUID, which is single-node
    (grounded, SQUID loop replaces one junction on ONE island -- see
    chip.py module docstring). This class is for an island PAIR joined
    by a SQUID loop instead of a single JJ.

    Deliberately does NOT store a fixed self.EJ like the base Qubit class.
    Instead, EJ1/EJ2/flux are stored and self.EJ is a computed property
    via DesignHelper.compute_squid_EJ, re-evaluated on every access (so
    changing self.flux -- e.g. via Chip.update_component(..., flux=...) --
    immediately changes EJ, frequency, anharmonicity, and zpfs() with no
    stale-cache risk). Qubit.zpfs()/.anharmonicity() are inherited
    UNCHANGED and read self.EJ like any other Qubit -- this is the whole
    point of making EJ a property rather than overriding those methods:
    the two previously-separate, separately-validated mechanisms (SQUID
    flux-tunability, differential-mode coupling) compose for free with
    zero duplicated physics.
    """

    def __init__(self, name, node_indices, design_helper, EJ1=None, EJ2=None,
                 EC=None, flux=0.0, flux_units="Phi0"):
        node_indices = (
            [node_indices] if isinstance(node_indices, int) else list(node_indices)
        )
        if len(node_indices) != 2:
            raise ValueError(
                f"FloatingSQUID '{name}' requires exactly 2 node_indices "
                f"(the two islands joined by the SQUID loop), got "
                f"{len(node_indices)}: {node_indices}. For a single-island "
                "SQUID-loop transmon (e.g. CircTransmonSQUID), use Qubit "
                "with EJ=compute_squid_EJ(...) instead -- that case is "
                "single-node, not FloatingSQUID."
            )
        # Bypass Qubit.__init__ (it expects a fixed EJ); call Component's
        # init directly and set EJ1/EJ2/EC/flux ourselves.
        Component.__init__(self, name, node_indices, design_helper)
        self.EJ1 = EJ1
        self.EJ2 = EJ2
        self.EC = EC
        self.flux = flux
        self.flux_units = flux_units

    @property
    def EJ(self):
        """Flux-tunable EJ, recomputed on every access from EJ1/EJ2/flux.
        Inherited Qubit.zpfs()/.frequency()/.anharmonicity() read this
        exactly as they would a plain Qubit's fixed self.EJ -- no
        override of those methods needed.
        """
        return self.dh.compute_squid_EJ(
            self.EJ1, self.EJ2, self.flux, flux_units=self.flux_units
        )

    @EJ.setter
    def EJ(self, value):
        # Qubit.__init__ would normally do self.EJ = EJ; since we bypass
        # that, this setter only exists so nothing breaks if anything
        # generic (e.g. update_component's setattr) ever tries to assign
        # to .EJ directly. Assigning EJ directly on a flux-tunable SQUID
        # is almost certainly a mistake -- raise rather than silently
        # creating a stale, disconnected attribute that EJ's getter would
        # then immediately shadow and ignore.
        raise AttributeError(
            f"FloatingSQUID '{self.name}'.EJ is computed from EJ1/EJ2/flux "
            "and cannot be set directly -- set self.flux (or EJ1/EJ2) "
            "instead, e.g. chip.update_component(name, flux=0.3)."
        )

    def phase_offset(self):
        """Flux-dependent phase offset of the SQUID's potential, via
        DesignHelper.compute_squid_phase_offset. Not used by zpfs()/
        frequency()/anharmonicity() (those only need the flux-tunable
        EJ magnitude) -- provided for completeness/diagnostics.
        """
        return self.dh.compute_squid_phase_offset(
            self.EJ1, self.EJ2, self.flux, flux_units=self.flux_units
        )