"""Chip: owns the full capacitance matrix (+ inverse), node labels, and the
component list. Cross-component physics (coupling) lives here because it
needs the full Cinv matrix, which no single Component has access to.

Confirmed against real SQDMetal/Palace output (SQuADDS Tutorial 7):
- Palace's capacitance-sim output is a Maxwell capacitance matrix in
  FARADS, not fF. set_capacitance_matrix defaults bare arrays to "F".
- Node labels are NEVER in the simulator output -- they're assigned by the
  user, conventionally after checking cap_sim.display_conductor_indices()
  to see which conductor index is which named piece of geometry. Treat
  node_labels as user-supplied; there is nothing to auto-infer from.
- This covers the capacitance-matrix coupling route only (0.5 * Q^T Cinv Q).
  The separate EPR/eigenmode route (participation ratios + eigenfrequencies)
  is implemented as Chip.from_eigenmode_epr() -- see handoff_plan.md Path B.
- Multi-node (floating transmon) coupling convention confirmed against a
  real SQDMetal/Palace TransmonPocket simulation: the two pads appear as
  two separate conductor indices (Cond2=bottom pad, Cond3=top pad). The
  physical charge operator is the differential mode Q=Q_pad0-Q_pad1, so the
  effective Cinv coupling element is Cinv[pad0,k] - Cinv[pad1,k] for any
  neighbor node k. CircTransmonSQUID is single-node (one island pad, two JJ
  leads to ground) and uses the n=1 path -- no differential mode needed.
- Inductance-matrix ingestion (set_inductance_matrix/linv) mirrors the
  capacitance-matrix machinery but is deliberately BARE: matrix + cached
  inverse + element lookup only, no coupling/EJ/Hamiltonian-parameter
  derivation built on top yet. Palace's magnetostatic sim
  (PALACE_Inductance_Simulation) returns one mutual inductance M=Phi/I per
  current-drive run, in inductance units directly (not inverted, unlike
  the capacitance route) -- an NxN matrix needs N separate Palace runs.
  Uses its OWN node_labels (inductance_node_labels), independent of the
  capacitance matrix's node_labels -- see set_inductance_matrix docstring.
"""

import numpy as np
import pint

from components import Component, GroundPlane


class Chip:
    def __init__(self, design_helper, node_labels):
        self.dh = design_helper
        self.node_labels = list(node_labels)
        self.capacitance_matrix = None  # pint.Quantity, NxN, set via set_capacitance_matrix
        self._cinv = None
        self.components = {}

        # Inductance-matrix machinery (Palace magnetostatic / mutual-
        # inductance sims). Deliberately SEPARATE node-label list from
        # node_labels above -- the conductors driven in an inductance sim
        # (current-carrying loops/leads) are not necessarily the same set
        # as the capacitance sim's conductors (islands/pads), so forcing
        # them to share indexing would be an unjustified assumption.
        # Bare matrix + inverse only for now -- intentionally NOT wired
        # into compute_coupling/EJ/any Hamiltonian-parameter derivation.
        self.inductance_node_labels = None
        self.inductance_matrix = None  # pint.Quantity, NxN, set via set_inductance_matrix
        self._linv = None

    # ------------------------------------------------------------------
    # Construction / mutation (Decision 4)
    # ------------------------------------------------------------------
    @classmethod
    def from_simulation(cls, design_helper, cap_matrix, node_labels, component_specs):
        """Build a Chip from a raw simulation output.

        cap_matrix : array-like or pint.Quantity, shape (N, N)
            Maxwell capacitance matrix. Bare arrays are interpreted as
            Farads (Palace's native units) -- NOT fF.
        node_labels : list[str], length N
            User-supplied node names matching cap_matrix's axis order.
            Recommended workflow: run cap_sim.display_conductor_indices()
            first to confirm which conductor index is which named piece of
            geometry, then build this list in that order.
        component_specs : list[dict], each one of:
            {"type": "qubit",     "name": str, "nodes": [str, ...], "EJ": ..., "EC": ...}
            {"type": "resonator", "name": str, "nodes": [str, ...],
             "omega_r": ... (or "inductance"/"capacitance")}
            {"type": "ground",    "name": str, "nodes": [str, ...]}
            "nodes" entries must match node_labels entries exactly.
        """
        chip = cls(design_helper, node_labels)
        chip.set_capacitance_matrix(cap_matrix)
        for spec in component_specs:
            chip.add_component(spec)
        return chip

    def set_capacitance_matrix(self, cap_matrix, units="F"):
        """(Re)load the capacitance matrix; invalidates the cached inverse.

        units: default "F" -- matches Palace's native capacitance-sim
        output (confirmed, see module docstring). Pass units="fF" explicitly
        if you've already rescaled the matrix yourself.
        """
        if isinstance(cap_matrix, self.dh.ureg.Quantity):
            cq = cap_matrix
        else:
            arr = np.asarray(cap_matrix, dtype=float)
            cq = arr * self.dh.ureg(units)

        n = cq.magnitude.shape[0]
        if cq.magnitude.shape != (n, n):
            raise ValueError(
                f"capacitance matrix must be square, got shape {cq.magnitude.shape}"
            )
        if n != len(self.node_labels):
            raise ValueError(
                f"capacitance matrix is {n}x{n} but {len(self.node_labels)} "
                "node_labels were given -- these must match."
            )

        self.capacitance_matrix = cq
        self._cinv = None  # invalidate cache

    def set_inductance_matrix(self, ind_matrix, node_labels, units="nH"):
        """Load a mutual-inductance matrix (e.g. from N separate Palace
        PALACE_Inductance_Simulation runs, one current-drive per row/col).

        ind_matrix : array-like or pint.Quantity, shape (N, N)
            Mutual inductance matrix. Palace's magnetostatic sim returns
            flux-per-ampere (M = Phi/I) directly in inductance units, NOT
            an already-inverted quantity -- unlike the capacitance route,
            there's no Farads-vs-inverse-Farads subtlety here. Bare arrays
            are interpreted in nH.
        node_labels : list[str], length N
            User-supplied labels for the current-carrying elements/loops
            this matrix indexes. Stored separately from self.node_labels
            (the capacitance matrix's node labels) -- these two index sets
            are not assumed to coincide; see __init__ docstring note.
        units : default "nH", matching DesignHelper.inductance_units.

        This only stores the matrix and invalidates the cached inverse.
        No coupling/EJ/Hamiltonian-parameter derivation consumes this yet
        -- bare matrix ingestion only, by design, for now.
        """
        if isinstance(ind_matrix, self.dh.ureg.Quantity):
            lq = ind_matrix
        else:
            arr = np.asarray(ind_matrix, dtype=float)
            lq = arr * self.dh.ureg(units)

        n = lq.magnitude.shape[0]
        if lq.magnitude.shape != (n, n):
            raise ValueError(
                f"inductance matrix must be square, got shape {lq.magnitude.shape}"
            )
        if n != len(node_labels):
            raise ValueError(
                f"inductance matrix is {n}x{n} but {len(node_labels)} "
                "node_labels were given -- these must match."
            )

        self.inductance_node_labels = list(node_labels)
        self.inductance_matrix = lq
        self._linv = None  # invalidate cache

    @classmethod
    def from_inductance_simulation(cls, design_helper, ind_matrix, node_labels):
        """Build a standalone inductance-matrix-only Chip (no capacitance
        matrix, no components) -- a convenience constructor for the case
        where you only have inductance-sim output so far and want to
        inspect/validate it on its own before a capacitance matrix or
        component list exists.

        For a Chip that already has a capacitance matrix/components and
        you're adding an inductance matrix to it, call
        chip.set_inductance_matrix(...) directly on the existing instance
        instead of this classmethod.
        """
        chip = cls(design_helper, node_labels=[])
        chip.set_inductance_matrix(ind_matrix, node_labels)
        return chip

    @property
    def linv(self):
        """Inverse mutual-inductance matrix, cached. Mirrors the cinv
        property exactly. Bare matrix only -- nothing currently consumes
        this; it's provided for inspection/future use.
        """
        if self.inductance_matrix is None:
            raise RuntimeError(
                "no inductance matrix loaded -- call set_inductance_matrix() "
                "or build via Chip.from_inductance_simulation() first."
            )
        if self._linv is None:
            mag = self.inductance_matrix.magnitude
            inv_mag = np.linalg.inv(mag)
            inv_units = 1.0 / self.inductance_matrix.units
            self._linv = inv_mag * inv_units
        return self._linv

    def inductance_element(self, label_a: str, label_b: str):
        """Direct L[i,j] lookup by inductance-matrix node label (NOT
        Component name -- this matrix has its own, separate label set,
        see set_inductance_matrix). Plain matrix access, no physics.
        """
        if self.inductance_matrix is None:
            raise RuntimeError("no inductance matrix loaded")
        i = self._resolve_inductance_node_index(label_a)
        j = self._resolve_inductance_node_index(label_b)
        return self.inductance_matrix[i, j]

    def _resolve_inductance_node_index(self, label):
        if label not in self.inductance_node_labels:
            raise ValueError(
                f"node label '{label}' not found in this Chip's "
                f"inductance_node_labels ({self.inductance_node_labels})."
            )
        return self.inductance_node_labels.index(label)

    def add_component(self, spec: dict):
        """Construct the right Component subclass from a spec dict, resolve
        its 'nodes' labels to indices via self.node_labels, and register it.
        """
        from components import Resonator, Qubit, FloatingSQUID  # local import: keeps
        # components.py ignorant of Chip (no circular import), Chip is the
        # only side that needs to know both module's classes.

        name = spec["name"]
        if name in self.components:
            raise ValueError(
                f"component '{name}' already exists -- use update_component "
                "to modify it, or remove_component first."
            )

        node_indices = self._resolve_node_indices(spec["nodes"])
        ctype = spec["type"]

        if ctype == "qubit":
            comp = Qubit(name, node_indices, self.dh, EJ=spec.get("EJ"), EC=spec.get("EC"))
        elif ctype == "floating_squid":
            comp = FloatingSQUID(
                name, node_indices, self.dh,
                EJ1=spec.get("EJ1"), EJ2=spec.get("EJ2"), EC=spec.get("EC"),
                flux=spec.get("flux", 0.0), flux_units=spec.get("flux_units", "Phi0"),
            )
        elif ctype == "resonator":
            comp = Resonator(
                name, node_indices, self.dh,
                omega_r=spec.get("omega_r"),
                inductance=spec.get("inductance"),
                capacitance=spec.get("capacitance"),
                kappa=spec.get("kappa"),
            )
        elif ctype == "ground":
            comp = GroundPlane(name, node_indices, self.dh)
        else:
            raise ValueError(
                f"unknown component type '{ctype}' for '{name}' -- "
                "expected 'qubit', 'floating_squid', 'resonator', or 'ground'."
            )

        self.components[name] = comp
        return comp

    def update_component(self, name: str, **kwargs):
        """Update an existing component's bare parameters in place (e.g.
        EJ/EC refined after a fit). Raises if `name` isn't found -- no
        silent no-op. Does NOT touch the capacitance matrix or node wiring;
        use remove_component + add_component if nodes need to change.
        """
        if name not in self.components:
            raise KeyError(f"no component named '{name}' on this Chip")
        comp = self.components[name]
        for key, value in kwargs.items():
            if not hasattr(comp, key):
                raise AttributeError(
                    f"{type(comp).__name__} '{name}' has no attribute '{key}'"
                )
            setattr(comp, key, value)
        return comp

    def remove_component(self, name: str):
        if name not in self.components:
            raise KeyError(f"no component named '{name}' on this Chip")
        del self.components[name]

    def _resolve_node_indices(self, node_names):
        indices = []
        for node_name in node_names:
            if node_name not in self.node_labels:
                raise ValueError(
                    f"node label '{node_name}' not found in this Chip's "
                    f"node_labels ({self.node_labels}) -- check spelling "
                    "against display_conductor_indices() output."
                )
            indices.append(self.node_labels.index(node_name))
        return indices

    # ------------------------------------------------------------------
    # Capacitance matrix access
    # ------------------------------------------------------------------
    @property
    def cinv(self):
        if self.capacitance_matrix is None:
            raise RuntimeError(
                "no capacitance matrix loaded -- call set_capacitance_matrix() "
                "or build via Chip.from_simulation() first."
            )
        if self._cinv is None:
            mag = self.capacitance_matrix.magnitude
            inv_mag = np.linalg.inv(mag)
            inv_units = 1.0 / self.capacitance_matrix.units
            self._cinv = inv_mag * inv_units
        return self._cinv

    def _cinv_element(self, comp_a: Component, comp_b: Component):
        """Effective Cinv coupling element between two components.

        Implements the projection of the capacitive interaction
        H_int = 0.5 * Q^T Cinv Q onto each component's charge operator.

        For a SINGLE-node component (e.g. grounded transmon, resonator pad),
        the charge operator is simply Q_i, so the effective Cinv entry is
        just Cinv[i, k].

        For a MULTI-node component (e.g. floating transmon with two
        distinguishable pads), the physical charge operator is the
        DIFFERENTIAL mode across the junction: Q = Q_pad0 - Q_pad1.

        CONFIRMED against a real SQDMetal/Palace capacitance simulation of a
        TransmonPocket (floating transmon). Palace conductor indexing showed:
          - Cond1: ground plane
          - Cond2: bottom pad  (node_indices[0])
          - Cond3: top pad     (node_indices[1])
        The two pads couple with OPPOSITE signs to any neighbor node k:
            Cinv_eff(A->k) = Cinv[pad0, k] - Cinv[pad1, k]

        NOT YET CONFIRMED for floating-SQUID transmons (a floating transmon
        where both pads are bridged by a SQUID loop instead of a single JJ).
        That geometry would still be two-node and the differential convention
        is physically expected to apply, but no cap sim has been run on it yet.
        For CircTransmonSQUID specifically: that component has a single circular
        island pad with two JJ leads running to ground -- it is single-node and
        uses the n=1 path here, no differential mode needed.

        Convention generalizes to N nodes: each component's nodes contribute
        with sign (+1) for node_indices[0] and (-1) for node_indices[1:].
        This is consistent with how compute_transmon_zpfs defines n_zpf (a
        single ZPF for the junction degree of freedom, not per-pad). A future
        3+ node component would need its own sign convention documented here
        before use.

        Single-node components are recovered as the n=1 special case (sign
        is just +1), so all existing single-node behavior is unchanged.

        Node ordering convention: pass nodes=['island_pad', 'ground_side_pad']
        in the component spec so that the island (the pad physically closer to
        the readout resonator) gets the +1 sign.
        """
        signs_a = [1 if idx == 0 else -1
                   for idx in range(len(comp_a.node_indices))]
        signs_b = [1 if idx == 0 else -1
                   for idx in range(len(comp_b.node_indices))]

        total = None
        for sign_p, node_p in zip(signs_a, comp_a.node_indices):
            for sign_q, node_q in zip(signs_b, comp_b.node_indices):
                term = (sign_p * sign_q) * self.cinv[node_p, node_q]
                total = term if total is None else total + term

        return total.to("1/F")

    # ------------------------------------------------------------------
    # Coupling: 0.5 * Q^T Cinv Q projected onto ZPFs
    # ------------------------------------------------------------------
    def compute_coupling(self, name_a: str, name_b: str):
        """g_ab from H_int = 0.5 * Q^T Cinv Q, projected onto each
        component's charge ZPF.

        Both Resonator.zpfs() and Qubit.zpfs() return q_zpf already in
        Coulombs (see components.py), so this method treats every component
        type identically -- no per-type branching here.

        Returns a frequency-units pint.Quantity.
        """
        if name_a == name_b:
            raise ValueError("compute_coupling requires two distinct components")
        comp_a = self.components[name_a]
        comp_b = self.components[name_b]

        _, q_zpf_a = comp_a.zpfs()
        _, q_zpf_b = comp_b.zpfs()
        cinv_ab = self._cinv_element(comp_a, comp_b).to("1/F")

        g_energy = (0.5 * cinv_ab * q_zpf_a * q_zpf_b).to("J")
        return self.dh._energy_as_frequency(g_energy)

    def coupling_matrix(self):
        """All pairwise couplings among non-ground components.

        Returns dict {(name_a, name_b): g} with each unordered pair
        appearing once. This is a generic, architecture-agnostic
        representation -- it doesn't assume any particular number of
        components, connectivity pattern, or downstream consumer.
        """
        names = [
            name for name, comp in self.components.items()
            if not isinstance(comp, GroundPlane)
        ]
        result = {}
        for idx_a, name_a in enumerate(names):
            for name_b in names[idx_a + 1:]:
                result[(name_a, name_b)] = self.compute_coupling(name_a, name_b)
        return result

    # ------------------------------------------------------------------
    # Dispersive regime
    # ------------------------------------------------------------------
    def compute_dispersive_shift(self, qubit_name: str, resonator_name: str):
        """Chi for a qubit-resonator pair, via DesignHelper.compute_dispersive_shift.

        Requires the qubit component to have a .frequency() and .anharmonicity()
        method (i.e. must be a Qubit), and the resonator to have omega_r set
        so its linear frequency can be recovered. Both must be registered on
        this Chip and the cap matrix must be loaded (for compute_coupling).

        Returns chi in freq_units (pint Quantity).
        """
        from components import Qubit, Resonator
        comp_q = self.components[qubit_name]
        comp_r = self.components[resonator_name]

        if not isinstance(comp_q, Qubit):
            raise TypeError(
                f"'{qubit_name}' is not a Qubit -- compute_dispersive_shift "
                "requires a Qubit component with .frequency() and .anharmonicity()."
            )
        if not isinstance(comp_r, Resonator):
            raise TypeError(
                f"'{resonator_name}' is not a Resonator."
            )
        if comp_r.omega_r is None:
            raise ValueError(
                f"Resonator '{resonator_name}' has no omega_r set -- needed "
                "to recover the linear resonator frequency for chi."
            )

        g = self.compute_coupling(qubit_name, resonator_name)
        f_qubit = comp_q.frequency()
        alpha = comp_q.anharmonicity()
        # omega_r is angular; convert to linear frequency
        omega_r = self.dh._as_quantity(comp_r.omega_r, "1/s")
        f_res = (omega_r / (2.0 * np.pi)).to(self.dh.freq_units)

        return self.dh.compute_dispersive_shift(g, f_qubit, f_res, alpha)

    def dressed_frequencies(self):
        """Dressed frequencies for all coupled pairs on this Chip.

        Handles two cases:

        QUBIT↔RESONATOR: dispersive approximation (Koch two-pole chi).
            Also sets 'is_dispersive_valid' and 'dispersive_ratio' in the
            sub-dict. If g/|Delta| > 0.1 the dispersive approximation is
            failing -- chi is still computed and returned but
            is_dispersive_valid=False flags it explicitly. Skipped if
            the resonator has no omega_r set.

        RESONATOR↔RESONATOR: exact 2x2 diagonalization.
            Normal-mode frequencies omega_plus/omega_minus + mixing angle.
            No dispersive approximation is made for linear-mode hybridization.

        QUBIT↔QUBIT: skipped (chi between two anharmonic modes requires
            a full Hamiltonian diagonalization beyond the scope here).

        Ground components are always excluded.
        Returns dict keyed by (name_a, name_b) in the order they appear
        in self.components.
        """
        from components import Qubit, Resonator
        result = {}
        names = [n for n, c in self.components.items() if not isinstance(c, GroundPlane)]
        for idx_a, name_a in enumerate(names):
            for name_b in names[idx_a + 1:]:
                comp_a = self.components[name_a]
                comp_b = self.components[name_b]

                # Qubit↔Resonator
                if (isinstance(comp_a, Qubit) and isinstance(comp_b, Resonator)):
                    q_name, r_name = name_a, name_b
                elif (isinstance(comp_b, Qubit) and isinstance(comp_a, Resonator)):
                    q_name, r_name = name_b, name_a
                else:
                    q_name, r_name = None, None

                if q_name is not None:
                    try:
                        comp_r = self.components[r_name]
                        comp_q = self.components[q_name]
                        if comp_r.omega_r is None:
                            continue
                        g = self.compute_coupling(q_name, r_name)
                        f_qubit = comp_q.frequency()
                        alpha = comp_q.anharmonicity()
                        omega_r = self.dh._as_quantity(comp_r.omega_r, "1/s")
                        f_res = (omega_r / (2.0 * np.pi)).to(self.dh.freq_units)
                        dressed = self.dh.compute_dressed_frequencies(
                            g, f_qubit, f_res, anharmonicity=alpha
                        )
                        dressed["g"] = g
                        result[(q_name, r_name)] = dressed
                    except (ValueError, TypeError):
                        continue

                # Resonator↔Resonator
                elif isinstance(comp_a, Resonator) and isinstance(comp_b, Resonator):
                    try:
                        if comp_a.omega_r is None or comp_b.omega_r is None:
                            continue
                        g = self.compute_coupling(name_a, name_b)
                        omega_a = self.dh._as_quantity(comp_a.omega_r, "1/s")
                        omega_b = self.dh._as_quantity(comp_b.omega_r, "1/s")
                        f_a = (omega_a / (2.0 * np.pi)).to(self.dh.freq_units)
                        f_b = (omega_b / (2.0 * np.pi)).to(self.dh.freq_units)
                        dressed = self.dh.compute_dressed_frequencies(
                            g, f_a, None, anharmonicity=None, mode_b_freq=f_b
                        )
                        dressed["g"] = g
                        result[(name_a, name_b)] = dressed
                    except (ValueError, TypeError):
                        continue
        return result

    def compute_purcell(self, qubit_name: str, resonator_name: str):
        """Purcell-limited T1 for a qubit coupled to a resonator with known kappa.

        Returns None if the resonator has no kappa set -- this is not an error,
        just means the linewidth hasn't been measured/provided yet.

        Returns dict from DesignHelper.compute_purcell_decay, with keys:
          'gamma_purcell', 'T1_purcell_s', 'T1_purcell_us', 'dispersive_ratio'
        """
        from components import Qubit, Resonator
        comp_q = self.components[qubit_name]
        comp_r = self.components[resonator_name]

        if not isinstance(comp_q, Qubit):
            raise TypeError(f"'{qubit_name}' is not a Qubit")
        if not isinstance(comp_r, Resonator):
            raise TypeError(f"'{resonator_name}' is not a Resonator")
        if comp_r.kappa is None:
            return None

        g = self.compute_coupling(qubit_name, resonator_name)
        f_qubit = comp_q.frequency()
        omega_r = self.dh._as_quantity(comp_r.omega_r, "1/s")
        f_res = (omega_r / (2.0 * np.pi)).to(self.dh.freq_units)
        kappa = self.dh._as_quantity(comp_r.kappa, self.dh.freq_units)

        return self.dh.compute_purcell_decay(g, f_qubit, f_res, kappa)

    def compute_higher_order_dispersive(self, qubit_name: str, resonator_name: str):
        """Higher-order (4th-order) dispersive corrections for a qubit-resonator pair.

        Returns dict with chi, chi_prime, vacuum_shift_res, vacuum_shift_q,
        is_dispersive_valid, dispersive_ratio.
        """
        from components import Qubit, Resonator
        comp_q = self.components[qubit_name]
        comp_r = self.components[resonator_name]

        if not isinstance(comp_q, Qubit):
            raise TypeError(f"'{qubit_name}' is not a Qubit")
        if not isinstance(comp_r, Resonator):
            raise TypeError(f"'{resonator_name}' is not a Resonator")
        if comp_r.omega_r is None:
            raise ValueError(f"Resonator '{resonator_name}' has no omega_r set")

        g = self.compute_coupling(qubit_name, resonator_name)
        f_qubit = comp_q.frequency()
        alpha = comp_q.anharmonicity()
        omega_r = self.dh._as_quantity(comp_r.omega_r, "1/s")
        f_res = (omega_r / (2.0 * np.pi)).to(self.dh.freq_units)

        return self.dh.compute_higher_order_dispersive(g, f_qubit, f_res, alpha)


    @classmethod
    def from_epr_params(cls, design_helper, epr_params,
                        qubit_name="Q1", resonator_name="R1",
                        qubit_nodes=None, resonator_nodes=None):
        """Build a minimal Chip from EPR-derived Hamiltonian parameters.

        Takes the dict returned by Chip.from_eigenmode_epr() and constructs
        a Chip with one Qubit and one Resonator populated from it -- no cap
        matrix, but compute_coupling/compute_dispersive_shift/dressed_frequencies
        all work because the components carry their own ZPF-generating params
        (EJ, EC, omega_r) and the Chip's Cinv is constructed from the EPR
        phi_zpf and g directly rather than from a simulated cap matrix.

        Specifically: we back-solve a synthetic 2×2 cap matrix whose Cinv
        element [0,1] reproduces the EPR g exactly when combined with the
        ZPFs from EJ/EC and omega_r. This is the minimum information needed
        to make compute_coupling return a value consistent with the EPR result,
        without requiring a separate cap sim.

        Parameters
        ----------
        epr_params : dict
            Output of Chip.from_eigenmode_epr(). Must contain at minimum:
            'Ej', 'Ec', 'g', 'cavity_frequency_linear'.
        qubit_name, resonator_name : str
            Names to register the components under.
        qubit_nodes, resonator_nodes : list[str] or None
            Node label lists for the synthetic components. Defaults to
            ['qubit_pad'] and ['res_pad'] respectively -- only used for
            bookkeeping, not for a real cap matrix lookup.
        """
        from components import Qubit, Resonator

        if qubit_nodes is None:
            qubit_nodes = ["qubit_pad"]
        if resonator_nodes is None:
            resonator_nodes = ["res_pad"]

        EJ = epr_params["Ej"]
        EC = epr_params["Ec"]
        f_res = epr_params["cavity_frequency_linear"]
        g_epr = epr_params["g"]

        omega_r = (2.0 * np.pi * f_res).to("1/s")

        # Build chip with synthetic node labels (no real cap sim)
        node_labels = qubit_nodes + resonator_nodes
        chip = cls(design_helper, node_labels)

        # Qubit component
        q_idx = list(range(len(qubit_nodes)))
        q_comp = Qubit(qubit_name, q_idx, design_helper, EJ=EJ, EC=EC)
        chip.components[qubit_name] = q_comp

        # Resonator's node indices into the combined node_labels list
        # (qubit_nodes + resonator_nodes), offset past the qubit's nodes.
        r_idx = list(range(len(qubit_nodes), len(qubit_nodes) + len(resonator_nodes)))
        # needs L or C too. If phi_zpf_sq is in epr_params (it always is from
        # from_eigenmode_epr), back-solve an effective capacitance from it:
        # phi_zpf = sqrt(hbar*Z/2) and Z = 1/(omega_r * C_eff)
        # => C_eff = hbar / (2 * phi_zpf^2 * omega_r)
        # This makes r_comp.zpfs() return values consistent with the EPR result.
        if "phi_zpf_sq" in epr_params and epr_params["phi_zpf_sq"] is not None:
            phi_zpf_sq = float(epr_params["phi_zpf_sq"])  # dimensionless (rad^2)
            # phi_zpf here is the reduced flux ZPF (dimensionless, in units of phi0=hbar/2e)
            # Z = 2 * phi_zpf^2 * phi0^2 / hbar = 2 * phi_zpf^2 * hbar / (2e)^2 / hbar
            # simpler: Z = hbar/(2 * q_zpf^2) and q_zpf from the resonator side
            # Just compute C_eff from omega_r and Z=hbar/(2*q_zpf^2):
            # phi_zpf^2 (in Wb^2) = hbar*Z/2; phi_zpf (Wb) = phi0 * phi_zpf_dimensionless
            phi0_J = design_helper.hbar.to("J*s").magnitude / (
                2 * design_helper.electron_charge.to("C").magnitude
            )  # hbar/(2e) in Wb
            phi_zpf_Wb_sq = phi_zpf_sq * phi0_J ** 2
            hbar_J = design_helper.hbar.to("J*s").magnitude
            Z_eff = 2 * phi_zpf_Wb_sq / hbar_J  # Ohms
            omega_r_val = omega_r.to("1/s").magnitude
            C_eff = 1.0 / (omega_r_val * Z_eff)  # Farads
            C_eff_q = design_helper._as_quantity(C_eff * 1e15, "fF")
        else:
            # No phi_zpf_sq available: use a default 50Ω CPW impedance
            # Z = 1/(omega_r * C) => C = 1/(omega_r * Z0)
            Z0 = 50.0
            C_eff_q = design_helper._as_quantity(
                1e15 / (omega_r.to("1/s").magnitude * Z0), "fF"
            )

        r_comp = Resonator(resonator_name, r_idx, design_helper,
                           omega_r=omega_r, capacitance=C_eff_q)
        chip.components[resonator_name] = r_comp

        # Back-solve synthetic Cinv[0,1] from g = 0.5 * Cinv_01 * q_zpf_q * q_zpf_r
        # so that compute_coupling returns g_epr exactly.
        _, q_zpf_q = q_comp.zpfs()
        _, q_zpf_r = r_comp.zpfs()
        g_energy = design_helper._frequency_as_energy(g_epr).to("J")
        cinv_01 = (2.0 * g_energy / (q_zpf_q * q_zpf_r)).to("1/F").magnitude

        # Construct minimal synthetic cap matrix whose inverse has the
        # right off-diagonal element(s). cinv_01 is in 1/F (typically
        # ~1e11-1e12 for realistic qubit-resonator couplings), so the
        # natural CAPACITANCE scale is 1/cinv_01 (~pF), NOT cinv_01 itself
        # -- using cinv_01's own (inverse-Farad) magnitude as a capacitance
        # scale makes the self-cap astronomically large (~1e13 "F") and the
        # exact quadratic below numerically degenerate (verified: this was
        # the actual bug -- roundtrip g came back as ~0 instead of g_epr).
        # The self-capacitance values don't enter compute_coupling, only
        # Cinv[i,j] with i≠j does -- so the scale choice is otherwise
        # arbitrary, just needs to be well-conditioned.
        n = len(node_labels)
        C_mag = np.zeros((n, n))
        c_scale = 100.0 / abs(cinv_01)  # Farads; well-conditioned self-cap
        for i in range(n):
            C_mag[i, i] = c_scale
        # Back-solve: we need inv(C)[i0, j0] == cinv_01 for the qubit/
        # resonator node pair, with both self-caps = c_scale.
        # For the 2x2 (or reduced 2x2 block) case:
        # inv([[c,m],[m,c]]) = [[c,-m],[-m,c]] / (c^2-m^2)
        # We want -m/(c^2-m^2) = cinv_01 => cinv_01*m^2 - m - cinv_01*c^2 = 0
        # Solved exactly (not a first-order approximation) for any n, since
        # only the qubit/resonator node pair carries a nonzero mutual term.
        i0, j0 = q_idx[0], r_idx[0]
        c = c_scale
        coeffs = [cinv_01, -1.0, -cinv_01 * c**2]
        roots = np.roots(coeffs)
        m = min(roots, key=abs)  # smaller-magnitude root -> well-conditioned matrix
        C_mag[i0, j0] = m
        C_mag[j0, i0] = m

        chip.set_capacitance_matrix(C_mag * design_helper.ureg("F"))
        return chip

    # ------------------------------------------------------------------
    # Output summary
    # ------------------------------------------------------------------
    def generate_output(self, title="Chip Summary"):
        """Generate a structured summary of all chip parameters as a
        pandas DataFrame (or plain dict-of-dicts if pandas isn't available).

        Columns: component name, type, nodes, bare parameters (EJ/EC or
        omega_r), frequency, anharmonicity (qubits only), plus all pairwise
        g and chi values as additional rows.

        Returns
        -------
        dict with keys:
          'components'  : DataFrame (or list of dicts) of per-component params
          'couplings'   : DataFrame (or list of dicts) of pairwise g values
          'dispersive'  : DataFrame (or list of dicts) of chi + dressed freqs
                          for all Qubit↔Resonator pairs
          'title'       : str
        """
        from components import Qubit, Resonator

        comp_rows = []
        for name, comp in self.components.items():
            if isinstance(comp, GroundPlane):
                ctype = "ground"
                row = {"name": name, "type": ctype, "nodes": comp.node_indices,
                       "EJ": None, "EC": None, "omega_r": None,
                       "frequency": None, "anharmonicity": None}
            elif isinstance(comp, Qubit):
                ctype = "qubit"
                freq = comp.frequency().to(self.dh.freq_units) if comp.EJ and comp.EC else None
                anharm = comp.anharmonicity().to(self.dh.freq_units) if comp.EC else None
                row = {
                    "name": name, "type": ctype, "nodes": comp.node_indices,
                    "EJ": str(comp.EJ.to(self.dh.freq_units)) if comp.EJ else None,
                    "EC": str(comp.EC.to(self.dh.freq_units)) if comp.EC else None,
                    "omega_r": None,
                    "frequency": str(freq) if freq else None,
                    "anharmonicity": str(anharm) if anharm else None,
                }
            elif isinstance(comp, Resonator):
                ctype = "resonator"
                if comp.omega_r is not None:
                    omega_r = self.dh._as_quantity(comp.omega_r, "1/s")
                    f_res_str = str((omega_r / (2.0 * np.pi)).to(self.dh.freq_units))
                else:
                    f_res_str = None
                row = {
                    "name": name, "type": ctype, "nodes": comp.node_indices,
                    "EJ": None, "EC": None,
                    "omega_r": str(comp.omega_r) if comp.omega_r else None,
                    "frequency": f_res_str,
                    "anharmonicity": None,
                }
            else:
                continue
            comp_rows.append(row)

        coupling_rows = []
        if self.capacitance_matrix is not None:
            try:
                cm = self.coupling_matrix()
                for (na, nb), g in cm.items():
                    coupling_rows.append({
                        "component_a": na,
                        "component_b": nb,
                        "g": str(g.to(self.dh.freq_units)),
                        "g_MHz": round(g.to("MHz").magnitude, 4),
                    })
            except Exception:
                pass

        dispersive_rows = []
        purcell_rows = []
        higher_order_rows = []
        try:
            df = self.dressed_frequencies()
            from components import Qubit, Resonator as Res
            for (name_a, name_b), vals in df.items():
                comp_a = self.components[name_a]
                comp_b = self.components[name_b]
                g_val = vals["g"]
                is_qr = (isinstance(comp_a, Qubit) or isinstance(comp_b, Qubit))

                if is_qr:
                    # Qubit-Resonator row
                    row = {
                        "pair": f"{name_a}↔{name_b}",
                        "type": "Q↔R",
                        "g_MHz": round(g_val.to("MHz").magnitude, 4),
                        "chi_MHz": round(vals["chi"].to("MHz").magnitude, 4),
                        "dressed_qubit_GHz": round(vals["qubit"].to("GHz").magnitude, 6),
                        "dressed_resonator_GHz": round(vals["resonator"].to("GHz").magnitude, 6),
                        "dispersive_ratio": round(vals.get("dispersive_ratio", float("nan")), 4),
                        "is_dispersive_valid": vals.get("is_dispersive_valid", True),
                    }
                    dispersive_rows.append(row)

                    # Purcell
                    q_name = name_a if isinstance(comp_a, Qubit) else name_b
                    r_name = name_b if isinstance(comp_a, Qubit) else name_a
                    try:
                        p = self.compute_purcell(q_name, r_name)
                        if p is not None:
                            purcell_rows.append({
                                "qubit": q_name,
                                "resonator": r_name,
                                "gamma_purcell_MHz": round(p["gamma_purcell"].to("MHz").magnitude, 6),
                                "T1_purcell_us": round(p["T1_purcell_us"], 4),
                            })
                    except Exception:
                        pass

                    # Higher-order
                    try:
                        ho = self.compute_higher_order_dispersive(q_name, r_name)
                        higher_order_rows.append({
                            "qubit": q_name,
                            "resonator": r_name,
                            "chi_MHz": round(ho["chi"].to("MHz").magnitude, 4),
                            "chi_prime_MHz": round(ho["chi_prime"].to("MHz").magnitude, 6),
                            "vac_shift_res_MHz": round(ho["vacuum_shift_res"].to("MHz").magnitude, 6),
                            "vac_shift_q_MHz": round(ho["vacuum_shift_q"].to("MHz").magnitude, 6),
                        })
                    except Exception:
                        pass
                else:
                    # Resonator-Resonator row
                    row = {
                        "pair": f"{name_a}↔{name_b}",
                        "type": "R↔R",
                        "g_MHz": round(g_val.to("MHz").magnitude, 4),
                        "omega_plus_GHz": round(vals["omega_plus"].to("GHz").magnitude, 6),
                        "omega_minus_GHz": round(vals["omega_minus"].to("GHz").magnitude, 6),
                        "mixing_angle_deg": round(np.degrees(vals["mixing_angle"]), 4),
                    }
                    dispersive_rows.append(row)
        except Exception:
            pass

        try:
            import pandas as pd
            out = {
                "title": title,
                "components": pd.DataFrame(comp_rows),
                "couplings": pd.DataFrame(coupling_rows) if coupling_rows else pd.DataFrame(),
                "dispersive": pd.DataFrame(dispersive_rows) if dispersive_rows else pd.DataFrame(),
                "purcell": pd.DataFrame(purcell_rows) if purcell_rows else pd.DataFrame(),
                "higher_order": pd.DataFrame(higher_order_rows) if higher_order_rows else pd.DataFrame(),
            }
        except ImportError:
            out = {
                "title": title,
                "components": comp_rows,
                "couplings": coupling_rows,
                "dispersive": dispersive_rows,
                "purcell": purcell_rows,
                "higher_order": higher_order_rows,
            }
        return out

    def print_output(self, title="Chip Summary"):
        """Print a human-readable summary of the chip to stdout."""
        out = self.generate_output(title=title)
        print(f"\n{'='*60}")
        print(f"  {out['title']}")
        print(f"{'='*60}")

        def _print_table(label, key):
            section = out.get(key)
            if section is None:
                return
            try:
                import pandas as pd
                if isinstance(section, pd.DataFrame) and not section.empty:
                    print(f"\n--- {label} ---")
                    print(section.to_string(index=False))
            except ImportError:
                if section:
                    print(f"\n--- {label} ---")
                    for row in section:
                        print("  " + ", ".join(
                            f"{k}={v}" for k, v in row.items() if v is not None))

        _print_table("Components", "components")
        _print_table("Pairwise Couplings", "couplings")
        _print_table("Dressed Frequencies & Dispersive", "dispersive")
        _print_table("Purcell Limits", "purcell")
        _print_table("Higher-Order Dispersive", "higher_order")
        print(f"{'='*60}\n")

    # ------------------------------------------------------------------
    # Path B: eigenmode/EPR ingestion (Track C)
    # ------------------------------------------------------------------
    # Deliberately NOT unified with the Path A (capacitance-matrix) route
    # above -- different raw simulation output (participation ratios +
    # complex eigenfrequencies, not a capacitance matrix), different
    # g-derivation (EPR/chi-based, not Q^T Cinv Q). See handoff_plan.md
    # Path B and multi_agent_split_plan.md Track C.
    @classmethod
    def from_eigenmode_epr(cls, design_helper, participation_ratios, eigenfrequencies,
                            Lj, qubit_idx, cavity_idx):
        """Derive qubit/cavity Hamiltonian parameters from a PALACE eigenmode
        + EPR simulation, via the SQuADDS Tutorial 7 formula chain.

        This does NOT build a Chip with a capacitance matrix/component list
        -- the EPR route never produces or needs a Cinv. It returns a plain
        dict of derived parameters (mirroring the tutorial's
        hamiltonian_params table), which the caller can then use however
        they like (e.g. to populate a Qubit/Resonator's bare EJ/EC by hand).

        Parameters
        ----------
        design_helper : DesignHelper
            Used only for the Lj -> EJ conversion (via Lj_to_EJ) and the
            hbar/electron_charge constants it already holds -- no physics
            formula here is rederived inline that DesignHelper already owns.
        participation_ratios : array-like, shape (n_modes, n_ports)
            `mat_mode_port` from eigen_sim.retrieve_mode_port_EPR_from_file().
            Only the port-0 column is used here (single-junction qubit) --
            this is the 2-mode case; see note below on N-mode.
        eigenfrequencies : array-like, length n_modes
            Linear (not angular) mode frequencies in Hz. Bare floats/arrays
            are interpreted as Hz, matching eig.csv's 'Re{f} (GHz)' column
            after the usual *1e9 conversion (consistent with the tutorial,
            which works in linear_qubit_freq/linear_res_freq Hz throughout).
        Lj : str, float, or pint.Quantity
            Junction inductance (the *linear* L_J used as a lumped element
            in the eigenmode sim, e.g. PALACE_Eigenmode_Simulation's
            create_port_JosephsonJunction(L_J=...)). Bare floats are
            interpreted in nH. Converted to E_J via design_helper.Lj_to_EJ.
        qubit_idx, cavity_idx : int
            Row indices into participation_ratios / eigenfrequencies
            identifying which simulated mode is the qubit-like mode and
            which is the cavity-like mode (determined by the caller, e.g.
            by inspecting mode frequencies/field plots as in the tutorial --
            this is NOT auto-detected, there is nothing in the raw EPR
            output that labels a mode "qubit" vs "cavity").

        Returns
        -------
        dict with keys: Ej, Ec, Lj, Ej_over_Ec, qubit_frequency,
        qubit_anharmonicity, cavity_frequency_linear, cavity_anharmonicity,
        chi, g, g_rwa, detuning, P_qubit, P_qubit_normalized, P_cav,
        P_cav_normalized, phi_zpf_sq. All frequency-like values are
        pint.Quantity in design_helper.freq_units; phi_zpf_sq and the
        participation ratios are dimensionless.

        N-mode note: this implementation is intentionally the 2-mode
        (one qubit mode, one cavity mode) case from the tutorial. The
        underlying chi/g formulas are pairwise by construction (a
        detuning-and-anharmonicity-denominator expression between exactly
        two modes), so generalizing to N participating modes isn't a
        matter of broadcasting this same formula over more rows -- it would
        need a real decision about how a multi-mode chi matrix collapses to
        pairwise g's, which the Tutorial 7 source doesn't address. Per
        multi_agent_split_plan.md Track C item 2, that generalization is
        left as a follow-up rather than guessed at here.
        """
        ureg = design_helper.ureg
        P = np.asarray(participation_ratios, dtype=float)
        freqs_hz = np.asarray(eigenfrequencies, dtype=float)

        P_qubit = abs(P[qubit_idx, 0])
        P_cav = abs(P[cavity_idx, 0])
        P_qubit_normalized = P_qubit / (P_qubit + P_cav)
        P_cav_normalized = P_cav / (P_qubit + P_cav)

        linear_qubit_freq = freqs_hz[qubit_idx]
        linear_res_freq = freqs_hz[cavity_idx]
        omega_qubit = 2.0 * np.pi * linear_qubit_freq
        omega_res = 2.0 * np.pi * linear_res_freq

        EJ = design_helper._frequency_as_energy(design_helper.Lj_to_EJ(Lj))
        hbar = design_helper.hbar

        # All of these are plain floats in SI (rad/s, J) -- omega_qubit/
        # omega_res are bare Python floats (not pint), so the only pint
        # quantity in this chain is EJ; keep it as a bare-magnitude float
        # in Joules for the arithmetic below and reattach units only at
        # the end, the same spirit as the tutorial's scipy.constants-based
        # calculation (and avoids carrying a fractional-hbar-power pint
        # quantity through several chained multiplications/divisions --
        # see the pint footgun note in handoff_plan.md / multi_agent_split_plan.md).
        EJ_J = EJ.to("J").magnitude
        hbar_J = hbar.to("J*s").magnitude

        phi_zpf_sq = P_qubit * hbar_J * 2.0 * omega_qubit / (2.0 * EJ_J)
        anharm_qubit = P_qubit ** 2 * hbar_J * omega_qubit ** 2 / (8.0 * EJ_J)
        anharm_res = P_cav ** 2 * hbar_J * omega_res ** 2 / (8.0 * EJ_J)
        Ec_J = anharm_qubit * hbar_J
        cross_kerr = (P_qubit * P_cav * hbar_J * omega_qubit * omega_res) / (4.0 * EJ_J)

        lamb_shift_qubit = anharm_qubit - cross_kerr / 2.0
        lamb_shift_res = anharm_res - cross_kerr / 2.0

        Delta = omega_res - omega_qubit
        Sigma = omega_res + omega_qubit
        alpha = -anharm_qubit
        chi = -cross_kerr
        denom = (alpha / (Delta * (Delta - alpha))) - (alpha / (Sigma * (Sigma + alpha)))
        g_rad = np.sqrt(chi / (2.0 * denom))

        disp_shift_qubit = cross_kerr
        delta_lin = (linear_res_freq - anharm_res / (2.0 * np.pi)) - \
                    (linear_qubit_freq - anharm_qubit / (2.0 * np.pi))
        g_rwa_rad = np.sqrt(disp_shift_qubit * delta_lin * (1.0 + delta_lin / anharm_qubit))

        qubit_freq_corrected = linear_qubit_freq - lamb_shift_qubit / (2.0 * np.pi)
        cavity_freq_corrected = linear_res_freq - lamb_shift_res / (2.0 * np.pi)

        f = design_helper.freq_units
        return {
            "Ej": design_helper._as_quantity(EJ_J / hbar_J / (2.0 * np.pi), "Hz").to(f),
            "Ec": design_helper._as_quantity(Ec_J / hbar_J / (2.0 * np.pi), "Hz").to(f),
            "Ej_over_Ec": EJ_J / Ec_J,
            "Lj": design_helper._as_quantity(Lj, "nH"),
            "qubit_frequency": design_helper._as_quantity(qubit_freq_corrected, "Hz").to(f),
            "qubit_frequency_linear": design_helper._as_quantity(linear_qubit_freq, "Hz").to(f),
            "qubit_anharmonicity": design_helper._as_quantity(anharm_qubit / (2.0 * np.pi), "Hz").to(f),
            "cavity_frequency": design_helper._as_quantity(cavity_freq_corrected, "Hz").to(f),
            "cavity_frequency_linear": design_helper._as_quantity(linear_res_freq, "Hz").to(f),
            "cavity_anharmonicity": design_helper._as_quantity(anharm_res / (2.0 * np.pi), "Hz").to(f),
            "chi": design_helper._as_quantity(cross_kerr / (2.0 * np.pi), "Hz").to(f),
            "g": design_helper._as_quantity(g_rad / (2.0 * np.pi), "Hz").to(f),
            "g_rwa": design_helper._as_quantity(g_rwa_rad / (2.0 * np.pi), "Hz").to(f),
            "detuning": design_helper._as_quantity(delta_lin, "Hz").to(f),
            "P_qubit": P_qubit,
            "P_qubit_normalized": P_qubit_normalized,
            "P_cav": P_cav,
            "P_cav_normalized": P_cav_normalized,
            "phi_zpf_sq": phi_zpf_sq,
        }

    def __repr__(self):
        return (
            f"Chip(nodes={self.node_labels}, "
            f"components={list(self.components.keys())})"
        )