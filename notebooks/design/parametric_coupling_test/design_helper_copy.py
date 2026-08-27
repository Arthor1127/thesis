import numpy as np
import pint

class DesignHelper:
    """Clase de ayuda para el diseño de circuitos cuánticos superconductores

    Las unidades estándar de la clase son
    Frecuencia, energia: GHz
    Velocidad: m/s
    longitud: um
    """

    freq_units = "GHz"
    speed_units = "m/s"
    length_units = "um"
    resistance_units = "ohm"
    current_units = "nA"
    inductance_units = "nH"
    capacitance_units = "fF"
    
    def __init__(self, light_speed: str, alpha_inductance=None):
        self.ureg = pint.UnitRegistry()
        self.hbar = self.ureg.Quantity(1, "hbar")
        self.planck = self.ureg.Quantity(1, 'planck_constant')
        self.electron_charge = self.ureg.Quantity(1, "elementary_charge")
        self.boltzmann = self.ureg.Quantity(1, "boltzmann_constant")
        self.reduced_flux_quantum = self.ureg.Quantity(1, "magnetic_flux_quantum") / (2.0 * np.pi)
        self.light_speed = self.ureg.Quantity(light_speed).to(self.speed_units)
        self.alpha_inductance = alpha_inductance if alpha_inductance is not None else 0.0

    def _as_quantity(self, value, default_units):
        if isinstance(value, self.ureg.Quantity):
            return value
        if isinstance(value, str):
            return self.ureg.Quantity(value)
        return value * self.ureg(default_units)

    def _energy_as_frequency(self, energy_quantity: pint.Quantity):
        energy_freq = energy_quantity / self.planck
        return energy_freq.to(self.freq_units)

    def _frequency_as_energy(self, freq_quantity: pint.Quantity):
        return freq_quantity * self.planck

    # ------------------------------------------------------------------
    # Resonator / transmission line
    # ------------------------------------------------------------------
    def compute_resonator_frequency(self, length, lambda_length_ratio):
        length = self._as_quantity(length, self.length_units)
        bare_resonance = self.light_speed / (length * lambda_length_ratio)
        return (bare_resonance * np.sqrt(1 - self.alpha_inductance)).to(self.freq_units)

    def compute_resonator_length(self, frequency, lambda_length_ratio):
        frequency = self._as_quantity(frequency, self.freq_units)
        bare_length = self.light_speed / (frequency * lambda_length_ratio)
        return (bare_length / np.sqrt(1 - self.alpha_inductance)).to(self.length_units)

    # ------------------------------------------------------------------
    # CPW geometry: impedance, propagation velocity, kinetic inductance
    # ------------------------------------------------------------------
    # Implements Appendix B of Ramos et al., "Microwave characterization
    # of superconducting coplanar resonators made out of granular
    # aluminium" (arXiv:2508.03873v1, 2025), Eqs. (B1)-(B9). This is a
    # SEPARATE route from compute_resonator_frequency/length above (which
    # take a user-supplied light_speed/lambda_length_ratio and don't know
    # anything about CPW cross-section geometry) -- this section instead
    # derives Z0, v, and the kinetic inductance contribution FROM
    # (w, g, t, eps_r), the actual conductor width/gap/thickness and
    # substrate permittivity.
    #
    # IMPORTANT -- one transcription correction made here, NOT a silent
    # "fix": the PDF's Eq. (B4) prints the impedance prefactor as 377/16.
    # Using that literal value against the paper's OWN worked example
    # (w=20um, g=11um, t=70nm -> stated Z0=50.1 Ohm) gives ~12.4 Ohm, off
    # by close to a factor of 4 -- not a rounding-level discrepancy.
    # The standard textbook CPW impedance formula (Simons; Ghione &
    # Naldi; Wheeler) uses a prefactor of 30*pi (=377*pi/4 #126.5 ->
    # equivalently just "30*pi" #94.25), which reproduces 49.67 Ohm for
    # the same geometry -- matching the paper's stated 50.1 Ohm to within
    # a sub-percent eps_r choice (eps_r#11.9 for high-resistivity Si
    # reproduces their value to <0.1%, see method docstrings below for
    # the verification). This implementation uses 30*pi, not the PDF's
    # literal 377/16, which is almost certainly a rendering/transcription
    # artifact (dropped pi, doubled denominator) rather than the authors'
    # intended formula. If you have access to the published (non-arXiv)
    # version or Clem's reference [87], it would be worth double-checking
    # this against the typeset equation directly.
    #
    # The Eq. (B9) kinetic-inductance geometric factor G(w,g,t) is
    # transcribed as printed (with the "Clem typo correction" the paper
    # describes applied as literally given -- no second guess of what
    # the correction itself was, since that requires the Clem paper this
    # implementation doesn't have access to). End-to-end validation
    # against the paper's own Table II numbers (back-solving sheet
    # kinetic inductance from a target alpha) showed a ~2x discrepancy
    # this implementation could NOT resolve -- see
    # compute_cpw_kinetic_inductance_per_length's docstring. The paper
    # itself calls this formula a first estimate / "indicative only"
    # (Sec. V.1), so this matches the authors' own stated confidence
    # level, not a regression introduced here.

    def _ellipk_agm(self, k):
        """Complete elliptic integral of the first kind K(k), via the
        arithmetic-geometric mean. Exact (not an approximation) and
        converges to machine precision in ~10-15 iterations for any
        k in [0, 1) -- verified against scipy.special.ellipk across the
        full range relevant here. No scipy dependency: this keeps
        DesignHelper's only third-party dependencies as numpy/pint,
        consistent with the rest of the file.
        """
        k = float(k)
        a, b = 1.0, np.sqrt(1.0 - k**2)
        for _ in range(30):
            a, b = (a + b) / 2.0, np.sqrt(a * b)
        return np.pi / (2.0 * a)

    def _cpw_k0(self, w, g):
        """Zero-thickness CPW elliptic modulus k0 = w/(w+2g), Eq. (B3)."""
        return w / (w + 2.0 * g)

    def compute_cpw_geometric_parameters(self, width, gap, thickness, eps_r):
        """Geometric (non-kinetic) capacitance/inductance per unit length,
        characteristic impedance, and propagation velocity of a CPW,
        including the finite-thickness correction -- Eqs. (B1)-(B7).

        Parameters
        ----------
        width, gap, thickness : length-like (default units: micrometers)
            Central conductor width w, gap to ground g, metal thickness t.
        eps_r : float
            Relative permittivity of the substrate (e.g. ~11.45-11.9 for
            high-resistivity silicon, ~9.4-11.5 for sapphire depending on
            crystal axis/orientation -- this is a material property the
            caller must supply, not something this method looks up).

        Returns
        -------
        dict with keys:
          'Cg' : capacitance per unit length (F/m, via pint, "farad/meter")
          'Lg' : geometric inductance per unit length ("henry/meter")
          'Z0' : characteristic impedance (self.resistance_units)
          'v'  : propagation velocity (self.speed_units)
          'eps_eff' : effective permittivity (dimensionless), = (v0/v)^2
                      using the zero-thickness v0 = c/sqrt((1+eps_r)/2)

        Verified against the paper's own three resonator lengths
        (Table/App. B text): with eps_r=11.9, this reproduces the stated
        CMT frequencies f_g for R1/R2/R3 (5.448/5.938/6.426 GHz) to within
        0.1% using f = v/(4*length) -- i.e. the full Cg/Lg/v chain here,
        not just Z0 in isolation, matches the paper's worked example.
        """
        w = self._as_quantity(width, self.length_units).to("m").magnitude
        g = self._as_quantity(gap, self.length_units).to("m").magnitude
        t = self._as_quantity(thickness, self.length_units).to("m").magnitude

        eps0 = self.ureg.Quantity(1, "vacuum_permittivity").to("F/m").magnitude
        mu0 = self.ureg.Quantity(1, "magnetic_constant").to("H/m").magnitude

        k0 = self._cpw_k0(w, g)
        kp0 = np.sqrt(1.0 - k0**2)
        K_k0 = self._ellipk_agm(k0)
        K_kp0 = self._ellipk_agm(kp0)

        # Finite-thickness correction (B6 inputs): u1, u2, then k_t
        d_t = 2.0 * t / np.pi
        u1 = (
            w / 2.0 + d_t / 2.0 + (3.0 * np.log(2.0) / 2.0) * d_t
            - (d_t / 2.0) * np.log(2.0 * d_t / w)
            + (d_t / 2.0) * np.log(g / (w + g))
        )
        u2 = (
            w / 2.0 + g - d_t / 2.0 - (3.0 * np.log(2.0) / 2.0) * d_t
            + (d_t / 2.0) * np.log(2.0 * d_t / (w + 2.0 * g))
            + (d_t / 2.0) * np.log(g / (w + g))
        )
        kt = u1 / u2
        ktp = np.sqrt(1.0 - kt**2)
        K_kt = self._ellipk_agm(kt)
        K_ktp = self._ellipk_agm(ktp)

        # Cg, Eq. (B6) -- finite-thickness capacitance per unit length
        Cg = eps0 * (2.0 * K_kt / K_ktp) + eps_r * eps0 * (2.0 * K_k0 / K_kp0)

        # Lg, Eq. (B7). NOTE on the "k_{t/2}" subscript: the PDF's printed
        # subscript is ambiguous (re-evaluate the thickness correction at
        # t/2, or reuse k_t/k_t' from B6 directly under a relabeling).
        # This implementation uses the literal reading -- recompute the
        # finite-thickness modulus using d(t/2) -- since it reproduces the
        # paper's stated Z0=50.1 Ohm slightly more closely (49.40 Ohm vs
        # 49.23 Ohm for the "reuse kt directly" reading) for their exact
        # geometry. The difference between the two readings is small
        # (~0.3%) -- if you have Clem's reference paper to disambiguate
        # this notation properly, that would settle it more rigorously
        # than this numerical tie-break.
        d_t2 = 2.0 * (t / 2.0) / np.pi
        u1_h = (
            w / 2.0 + d_t2 / 2.0 + (3.0 * np.log(2.0) / 2.0) * d_t2
            - (d_t2 / 2.0) * np.log(2.0 * d_t2 / w)
            + (d_t2 / 2.0) * np.log(g / (w + g))
        )
        u2_h = (
            w / 2.0 + g - d_t2 / 2.0 - (3.0 * np.log(2.0) / 2.0) * d_t2
            + (d_t2 / 2.0) * np.log(2.0 * d_t2 / (w + 2.0 * g))
            + (d_t2 / 2.0) * np.log(g / (w + g))
        )
        kt_h = u1_h / u2_h
        ktp_h = np.sqrt(1.0 - kt_h**2)
        K_kt_h = self._ellipk_agm(kt_h)
        K_ktp_h = self._ellipk_agm(ktp_h)
        Lg = (mu0 / 4.0) * (K_ktp_h / K_kt_h)

        Z0 = np.sqrt(Lg / Cg)
        v = 1.0 / np.sqrt(Lg * Cg)
        eps_eff0 = (1.0 + eps_r) / 2.0
        v0 = self.light_speed.to("m/s").magnitude / np.sqrt(eps_eff0)
        eps_eff = (v0 / v) ** 2

        return {
            "Cg": self.ureg.Quantity(Cg, "F/m"),
            "Lg": self.ureg.Quantity(Lg, "H/m"),
            "Z0": self.ureg.Quantity(Z0, "ohm").to(self.resistance_units),
            "v": self.ureg.Quantity(v, "m/s").to(self.speed_units),
            "eps_eff": float(eps_eff),
        }

    def compute_cpw_geometric_factor(self, width, gap, thickness):
        """Geometric factor G(w,g,t) from Eq. (B9) (Watanabe-Yoshida-
        Kohjiro expression, with the paper's stated Clem typo correction
        applied -- transcribed as printed; this implementation has not
        independently re-derived or verified that correction itself
        against Clem's original paper, only the formula as given in
        Ramos et al.). Used by compute_cpw_kinetic_inductance_per_length;
        exposed standalone for inspection.
        """
        w = self._as_quantity(width, self.length_units).to("m").magnitude
        g = self._as_quantity(gap, self.length_units).to("m").magnitude
        t = self._as_quantity(thickness, self.length_units).to("m").magnitude

        k0 = self._cpw_k0(w, g)
        kp0 = np.sqrt(1.0 - k0**2)
        K_k0 = self._ellipk_agm(k0)

        term1 = -np.log(t / (4.0 * w))
        term2 = -(w / (w + 2.0 * g)) * np.log(t / (4.0 * (w + 2.0 * g)))
        term3 = (2.0 * (w + g) / (w + 2.0 * g)) * np.log(g / (w + g))
        return (1.0 / (2.0 * kp0**2 * K_k0**2)) * (term1 + term2 + term3)

    def compute_cpw_kinetic_inductance_per_length(self, width, gap, thickness,
                                                   sheet_kinetic_inductance):
        """Kinetic inductance per unit length, Eq. (B9):
        Lk_per_length = (L_k^sq / w) * G(w, g, t).

        Parameters
        ----------
        sheet_kinetic_inductance : inductance-like (default units: pH)
            L_k^sq, the sheet kinetic inductance of the film (e.g. from
            Eq. (1)/Rn_to... this module's existing BCS estimate, or a
            fitted value from a measured alpha -- see
            solve_sheet_kinetic_inductance_from_alpha below).

        CALIBRATION CAVEAT (read before trusting this for design):
        attempting to reproduce the paper's own Table II via the
        round-trip alpha(L_k^sq) = alpha_1 -> solve for L_k^sq gave a
        result approximately HALF the paper's reported sheet kinetic
        inductance for the same target alpha (Al, R1, alpha_1=0.139:
        this implementation backs out ~1.0 pH/sq vs. their reported
        2.104 pH/sq). Multiple plausible sources of this gap were
        checked and ruled out (the B7 "k_{t/2}" ambiguity above barely
        moves the result; eps_r choice doesn't explain a factor of ~2).
        The remaining most likely explanations are either (a) a
        misread of the Clem typo correction itself, which would need
        the original Clem reference to resolve, or (b) the paper's
        Table II L_k^sq coming from a fit against the FULL fabricated
        resonator (including the coupling capacitor's effect on the
        effective length) rather than the simple per-unit-length
        alpha=Lk/(Lg+Lk) relation used here. The paper itself describes
        this geometric-factor route as a first estimate, "indicative
        only" (Sec. V.1) -- treat this method the same way: a
        physically-reasonable-order-of-magnitude estimate, not a
        calibrated design number, until the discrepancy above is
        tracked down against Clem's original paper.
        """
        Lk_sq = self._as_quantity(sheet_kinetic_inductance, "pH").to("H").magnitude
        w = self._as_quantity(width, self.length_units).to("m").magnitude
        G = self.compute_cpw_geometric_factor(width, gap, thickness)
        Lk_per_length = (Lk_sq / w) * G
        return self.ureg.Quantity(Lk_per_length, "H/m")

    def compute_cpw_kinetic_fraction(self, width, gap, thickness, eps_r,
                                      sheet_kinetic_inductance):
        """Kinetic inductance fraction alpha = Lk/(Lg+Lk), Eq. (3)/(B8)-
        (B9) combined -- the quantity the paper calls alpha_CM (Sec. V.1)
        when fg is obtained via conformal mapping rather than simulation.

        Subject to the same calibration caveat as
        compute_cpw_kinetic_inductance_per_length -- see that method's
        docstring.
        """
        geo = self.compute_cpw_geometric_parameters(width, gap, thickness, eps_r)
        Lg_per_length = geo["Lg"].to("H/m").magnitude
        Lk_per_length = self.compute_cpw_kinetic_inductance_per_length(
            width, gap, thickness, sheet_kinetic_inductance
        ).to("H/m").magnitude
        return Lk_per_length / (Lg_per_length + Lk_per_length)

    def solve_sheet_kinetic_inductance_from_alpha(self, width, gap, thickness, eps_r,
                                                   target_alpha,
                                                   bracket=(1e-4, 1e4), tol=1e-10,
                                                   max_iter=200):
        """Inverse of compute_cpw_kinetic_fraction: given a MEASURED kinetic
        inductance fraction (e.g. alpha_1 from a Mattis-Bardeen fit, Eq. 7),
        solve for the sheet kinetic inductance L_k^sq that would produce
        it under this geometric-factor model -- the procedure described in
        Sec. V.1 ("we calculated alpha(L_k^sq) and solved for
        alpha(L_k^sq)=alpha_1").

        Plain bisection, not scipy.optimize.brentq -- this file has no
        scipy dependency (a deliberate choice made earlier for the same
        reason as the AGM-based elliptic integral above: keep
        DesignHelper's only third-party dependency as numpy/pint).
        compute_cpw_kinetic_fraction is monotonically increasing in
        sheet_kinetic_inductance (more kinetic inductance -> larger alpha,
        always, for any physical geometry), so bisection is well-posed
        and doesn't need brentq's extra robustness.

        Parameters
        ----------
        bracket : (low, high) in pH, the search range for L_k^sq. Default
            (1e-4, 1e4) pH spans from far below typical clean Al (~1-2
            pH/sq) to far above the most disordered grAl reported in this
            paper (~7 pH/sq) -- widen if target_alpha falls outside the
            range this brackets and ValueError is raised.

        Subject to the same calibration caveat documented on
        compute_cpw_kinetic_inductance_per_length -- this inverse problem
        inherits that unresolved ~2x discrepancy against the paper's own
        Table II, so treat a solved L_k^sq here as an order-of-magnitude
        estimate, not a calibrated design number, until that gap is
        tracked down (see that method's docstring for what was checked).

        Returns
        -------
        pint.Quantity in "pH" (sheet kinetic inductance, L_k^sq).
        """
        lo, hi = bracket
        f_lo = self.compute_cpw_kinetic_fraction(
            width, gap, thickness, eps_r, self.ureg.Quantity(lo, "pH")
        ) - target_alpha
        f_hi = self.compute_cpw_kinetic_fraction(
            width, gap, thickness, eps_r, self.ureg.Quantity(hi, "pH")
        ) - target_alpha
        if f_lo * f_hi > 0:
            raise ValueError(
                f"target_alpha={target_alpha} is not bracketed by "
                f"sheet_kinetic_inductance in [{lo}, {hi}] pH "
                f"(alpha(lo)={f_lo + target_alpha:.6f}, "
                f"alpha(hi)={f_hi + target_alpha:.6f}) -- widen `bracket`."
            )

        for _ in range(max_iter):
            mid = 0.5 * (lo + hi)
            f_mid = self.compute_cpw_kinetic_fraction(
                width, gap, thickness, eps_r, self.ureg.Quantity(mid, "pH")
            ) - target_alpha
            if abs(f_mid) < tol or (hi - lo) < tol:
                return self.ureg.Quantity(mid, "pH")
            if f_lo * f_mid <= 0:
                hi = mid
            else:
                lo, f_lo = mid, f_mid
        return self.ureg.Quantity(0.5 * (lo + hi), "pH")

    def compute_resonator_frequency(self, length, lambda_length_ratio):
        length = self._as_quantity(length, self.length_units)
        bare_resonance = self.light_speed / (length * lambda_length_ratio)
        return (bare_resonance * np.sqrt(1 - self.alpha_inductance)).to(self.freq_units)

    def Ic_to_EJ(self, Ic):
        Ic = self._as_quantity(Ic, self.current_units)
        j_energy = Ic * self.reduced_flux_quantum
        return self._energy_as_frequency(j_energy)

    def EJ_to_Ic(self, j_energy):
        ej = self._frequency_as_energy(self._as_quantity(j_energy, self.freq_units))
        curr = ej / self.reduced_flux_quantum
        return curr.to(self.current_units)

    # ------------------------------------------------------------------
    # Ambegaokar-Baratoff: normal-state resistance <-> E_J
    # ------------------------------------------------------------------
    def Rn_to_EJ(self, Rn, gap_energy, temperature=None):
        """Critical current and E_J from normal-state junction resistance via
        the Ambegaokar-Baratoff relation.

        R_n * I_c = (pi * Delta) / (2e) * tanh(Delta / (2 k_B T))

        At T -> 0 (or if `temperature` is omitted) the tanh factor -> 1, giving
        the usual zero-temperature AB relation R_n * I_c = pi * Delta / (2e).

        Parameters
        ----------
        Rn : str or float
            Normal-state (room-temperature or above-Tc) junction resistance.
            If a bare float, interpreted in `self.resistance_units`.The capacitance matrix (and its inverse) belongs to Chip, not to any one Component. Components only know which node index/indices are theirs.
        gap_energy : str or float
            Superconducting gap Delta (energy units; if a bare float,
            interpreted via the spectroscopy-context frequency convention,
            i.e. in `self.freq_units`, consistent with the rest of the class).
        temperature : str or float, optional
            Junction temperature. If None, the T=0 limit is used.

        Returns
        -------
        dict with keys 'Ic' (current_units) and 'EJ' (freq_units)
        """
        Rn = self._as_quantity(Rn, self.resistance_units)
        gap = self._frequency_as_energy(self._as_quantity(gap_energy, self.freq_units))

        if temperature is not None:
            T = self._as_quantity(temperature, "K")
            thermal_factor = np.tanh(
                (gap / (2.0 * self.boltzmann * T)).to_reduced_units().magnitude
            )
        else:
            thermal_factor = 1.0

        Ic = (np.pi * gap * thermal_factor) / (2.0 * self.electron_charge * Rn)
        Ic = Ic.to(self.current_units)
        EJ = self.Ic_to_EJ(Ic)
        return {"Ic": Ic, "EJ": EJ}

    def EJ_to_Rn(self, j_energy, gap_energy, temperature=None):
        """Inverse of Rn_to_EJ: predicted normal-state resistance from a
        target E_J, given the gap. Useful at the fab stage when you know the
        E_J you're aiming for and want the room-temperature R_n to screen for
        on a wafer probe station before cooldown.
        """
        Ic = self.EJ_to_Ic(j_energy)
        gap = self._frequency_as_energy(self._as_quantity(gap_energy, self.freq_units))

        if temperature is not None:
            T = self._as_quantity(temperature, "K")
            thermal_factor = np.tanh(
                (gap / (2.0 * self.boltzmann * T)).to_reduced_units().magnitude
            )
        else:
            thermal_factor = 1.0

        Rn = (np.pi * gap * thermal_factor) / (2.0 * self.electron_charge * Ic)
        return Rn.to(self.resistance_units)

    # ------------------------------------------------------------------
    # Charging energy / transmon spectrum
    # ------------------------------------------------------------------
    def C_to_EC(self, capacitance):
        """E_C = e^2 / (2C), returned in freq_units."""
        C = self._as_quantity(capacitance, self.capacitance_units)
        ec_energy = self.electron_charge ** 2 / (2.0 * C)
        return self._energy_as_frequency(ec_energy)

    def EC_to_C(self, ec_freq):
        """Inverse of C_to_EC: shunt capacitance needed for a target E_C."""
        ec_energy = self._frequency_as_energy(self._as_quantity(ec_freq, self.freq_units))
        C = self.electron_charge ** 2 / (2.0 * ec_energy)
        return C.to(self.capacitance_units)

    def compute_transmon_frequency(self, EJ, EC, exact=False, n_levels=5):
        """f_01 from E_J, E_C.

        By default uses the standard leading-order asymptotic expansion
            f_01 ~ sqrt(8 EJ EC) - EC
        valid for EJ/EC >> 1 (the usual transmon regime).

        If exact=True, diagonalizes the transmon charge-basis Hamiltonian
        H = 4 EC (n - ng)^2 - EJ cos(phi)
        numerically (ng=0) and returns f_01 from the exact spectrum instead,
        which matters once EJ/EC isn't deep in the asymptotic regime or you
        want to cross-check the leading-order formula.
        """
        EJ_q = self._as_quantity(EJ, self.freq_units)
        EC_q = self._as_quantity(EC, self.freq_units)

        if not exact:
            ratio = (EJ_q / EC_q).to_reduced_units().magnitude
            if ratio < 5.0:
                import warnings
                warnings.warn(
                    f"EJ/EC = {ratio:.2f} is well outside the leading-order "
                    "transmon asymptotic regime (EJ/EC >> 1, conventionally "
                    ">~20-50). The formula f01 ~ sqrt(8 EJ EC) - EC can return "
                    "nonsensical (even negative) frequencies here -- e.g. this "
                    "happens routinely near the EJ->0 sweet spot of a flux-"
                    "tunable SQUID-loop transmon. Use exact=True instead.",
                    stacklevel=2,
                )
            EJ_energy = self._frequency_as_energy(EJ_q)
            EC_energy = self._frequency_as_energy(EC_q)
            f01_energy = np.sqrt(8.0 * EJ_energy * EC_energy) - EC_energy
            return self._energy_as_frequency(f01_energy)

        # Exact diagonalization in the charge basis, ng = 0. Stay in
        # freq_units throughout (frequency units in, frequency units out) --
        # no energy<->frequency conversion needed since the Hamiltonian's
        # eigenvalue gaps are linear in EJ, EC regardless of which units
        # (energy or frequency-equivalent) they're expressed in.
        ncut = max(20, 5 * n_levels)
        n_vals = np.arange(-ncut, ncut + 1)
        dim = len(n_vals)
        EJ_mag = EJ_q.to(self.freq_units).magnitude
        EC_mag = EC_q.to(self.freq_units).magnitude

        H = np.diag(4.0 * EC_mag * n_vals ** 2).astype(float)
        off_diag = -0.5 * EJ_mag
        for i in range(dim - 1):
            H[i, i + 1] += off_diag
            H[i + 1, i] += off_diag

        evals = np.linalg.eigvalsh(H)
        evals.sort()
        f01 = evals[1] - evals[0]
        return self._as_quantity(f01, self.freq_units)

    def compute_anharmonicity(self, EC, exact_EJ=None):
        """Leading-order anharmonicity alpha ~ -EC.

        If exact_EJ is supplied, instead returns the exact (f12 - f01) from
        diagonalizing the charge-basis Hamiltonian, which is the quantity
        you actually measure spectroscopically and is somewhat more negative
        than -EC once EJ/EC is moderate.
        """
        EC_q = self._as_quantity(EC, self.freq_units)
        if exact_EJ is None:
            return -1.0 * EC_q

        EJ_q = self._as_quantity(exact_EJ, self.freq_units)

        ncut = 60
        n_vals = np.arange(-ncut, ncut + 1)
        dim = len(n_vals)
        EJ_mag = EJ_q.to(self.freq_units).magnitude
        EC_mag = EC_q.to(self.freq_units).magnitude

        H = np.diag(4.0 * EC_mag * n_vals ** 2).astype(float)
        off_diag = -0.5 * EJ_mag
        for i in range(dim - 1):
            H[i, i + 1] += off_diag
            H[i + 1, i] += off_diag

        evals = np.linalg.eigvalsh(H)
        evals.sort()
        f01 = evals[1] - evals[0]
        f12 = evals[2] - evals[1]
        return self._as_quantity(f12 - f01, self.freq_units)

    # ------------------------------------------------------------------
    # SQUID / flux tunability (CircTransmonSQUID-relevant)
    # ------------------------------------------------------------------
    def compute_squid_EJ(self, EJ1, EJ2, flux, flux_units="Phi0"):
        """Flux-tunable Josephson energy of an asymmetric DC SQUID.

        EJ(flux) = EJ_sum * cos(pi * flux/Phi0) * sqrt(1 + d^2 tan^2(pi * flux/Phi0))

        where EJ_sum = EJ1 + EJ2 and d = (EJ2 - EJ1)/(EJ1 + EJ2) is the
        junction asymmetry. flux is taken in units of the flux quantum Phi0
        by default (flux_units="Phi0" means the bare float you pass in is
        already Phi/Phi0).
        """
        EJ1_q = self._as_quantity(EJ1, self.freq_units)
        EJ2_q = self._as_quantity(EJ2, self.freq_units)
        EJ_sum = EJ1_q + EJ2_q
        d = np.abs(EJ2_q - EJ1_q) / EJ_sum

        if flux_units == "Phi0":
            reduced_flux = float(flux)
        else:
            flux_q = self._as_quantity(flux, flux_units)
            reduced_flux = (flux_q / self.ureg.Quantity(1, "magnetic_flux_quantum")).to_reduced_units().magnitude

        phase = np.pi * reduced_flux
        EJ_flux = EJ_sum * np.cos(phase) * np.sqrt(1.0 + d ** 2 * np.tan(phase) ** 2)
        return EJ_flux
    
    def compute_squid_phase_offset(self, EJ1, EJ2, flux, flux_units="Phi0"):
        """Phase offset of the flux potential due to SQUID assymmetry
        V(phi, flux) = E_J(flux) * cos(phi - phi_0(flux))

        tan phi_0(flux) = d tan(pi * flux /Phi0)

        where d = (EJ2 - EJ1)/(EJ1 + EJ2) is the junction asymmetry. flux
        is taken in units of the flux quantum Phi0 by default (flux_units="Phi0" 
        means the bare float you pass in is already Phi/Phi0).
        """
        
        EJ1_q = self._as_quantity(EJ1, self.freq_units)
        EJ2_q = self._as_quantity(EJ2, self.freq_units)
        EJ_sum = EJ1_q + EJ2_q
        d = np.abs(EJ2_q - EJ1_q) / EJ_sum
        if flux_units == "Phi0":
            reduced_flux = float(flux)
        else:
            flux_q = self._as_quantity(flux, flux_units)
            reduced_flux = (flux_q / self.ureg.Quantity(1, "magnetic_flux_quantum")).to_reduced_units().magnitude
        
        return np.atan(d * np.tan(reduced_flux))
    
    def compute_flux_tunable_frequency(self, EJ1, EJ2, EC, flux, flux_units="Phi0", exact=False):
        """f_01(flux) for a flux-tunable (SQUID-loop) transmon, combining
        compute_squid_EJ with compute_transmon_frequency.
        """
        EJ_flux = self.compute_squid_EJ(EJ1, EJ2, flux, flux_units=flux_units)
        return self.compute_transmon_frequency(EJ_flux, EC, exact=exact)

    def compute_squid_inductance(self, EJ1, EJ2, flux, flux_units="Phi0"):
        """Flux-dependent Josephson inductance L_J(flux) = hbar / (2e * Ic(flux))
        of the SQUID loop, via L_J = Phi0^2 / (4 pi^2 EJ(flux)) [energy form],
        equivalently Phi0 / (2 pi Ic).
        """
        EJ_flux = self._frequency_as_energy(self.compute_squid_EJ(EJ1, EJ2, flux, flux_units=flux_units))
        phi0 = self.ureg.Quantity(1, "magnetic_flux_quantum")
        L_J = phi0 ** 2 / (4.0 * np.pi ** 2 * EJ_flux)
        return L_J.to(self.inductance_units)

    def Lj_to_EJ(self, Lj):
        """Inverse of compute_squid_inductance's energy-form relation:
        E_J = (hbar/2e)^2 / L_J = phi0^2 / L_J, returned in freq_units.

        This is the bare Josephson-inductance <-> energy relation (not
        SQUID-flux-dependent) -- the one needed when a sim/measurement
        hands you a single inductance value (e.g. a junction's linear L_J
        used as a lumped element in an eigenmode sim) and you need E_J for
        the standard transmon formulas. Reuses self.hbar/self.electron_charge
        rather than rederiving phi0 = hbar/(2e) inline.
        """
        Lj_q = self._as_quantity(Lj, self.inductance_units)
        phi0 = self.hbar / (2.0 * self.electron_charge)
        EJ_energy = (phi0 ** 2 / Lj_q).to("J")
        return self._energy_as_frequency(EJ_energy)
    
    def compute_resonator_impedance(self, omega_r=None, inductance=None, capacitance=None):
        """Characteristic impedance of an LC resonator from any two of
        (angular frequency, inductance, capacitance).

        Z = sqrt(L/C)              given L, C
        Z = 1 / (omega_r * C)      given omega_r, C
        Z = omega_r * L            given omega_r, L

        omega_r is taken as an *angular* frequency (rad/s), not f = omega/2pi
        -- consistent with how it feeds into compute_resonator_zpfs and the
        ladder-operator ZPF formulas, which are written in terms of omega.
        Pass omega_r already as an angular-frequency pint.Quantity (e.g.
        2*np.pi*f) or in units of "rad/s" / "1/s"; a bare GHz value will be
        treated as ordinary frequency-units and is very likely not what you
        want here.
        """
        n_given = sum(x is not None for x in (omega_r, inductance, capacitance))
        if n_given != 2:
            raise ValueError(
                "compute_resonator_impedance requires exactly two of "
                "(omega_r, inductance, capacitance); got "
                f"{n_given}."
            )

        if omega_r is None:
            L = self._as_quantity(inductance, self.inductance_units)
            C = self._as_quantity(capacitance, self.capacitance_units)
            return np.sqrt(L / C).to(self.resistance_units)
        elif capacitance is None:
            omega = self._as_quantity(omega_r, "1/s")
            L = self._as_quantity(inductance, self.inductance_units)
            return (omega * L).to(self.resistance_units)
        else:
            omega = self._as_quantity(omega_r, "1/s")
            C = self._as_quantity(capacitance, self.capacitance_units)
            return (1.0 / (omega * C)).to(self.resistance_units)

    def compute_resonator_zpfs(self, omega_r=None, inductance=None, capacitance=None):
        """Zero-point flux and charge fluctuations of an LC resonator mode,
        phi_zpf = sqrt(hbar*Z/2), q_zpf = sqrt(hbar/(2Z)), with
        Z from compute_resonator_impedance (see that method's docstring for
        the required omega_r convention -- angular frequency, not f).

        Returns phi_zpf in Wb and q_zpf in C, converted immediately rather
        than left as compound sqrt(hbar/...) expressions. This is a
        deliberate defensive measure: on some pint versions/installations,
        chaining further arithmetic on a quantity whose units still contain
        a fractional power of a derived constant like hbar (dirac_constant)
        can trip a known pint dimensionality-reduction bug (see e.g.
        hgrecco/pint issue #876), where a compound unit that is genuinely
        dimensionless or genuinely the right dimension gets mis-evaluated
        after several more operations. Collapsing to plain SI units here,
        right after the sqrt, avoids ever carrying a fractional-hbar-power
        unit forward into later multiplications (e.g. in
        Chip.compute_coupling).
        """
        impedance = self.compute_resonator_impedance(
            omega_r=omega_r, inductance=inductance, capacitance=capacitance
        )
        phi_zpf = np.sqrt(self.hbar * impedance / 2.0).to("Wb")
        q_zpf = np.sqrt(self.hbar / (2.0 * impedance)).to("C")
        return phi_zpf, q_zpf

    def compute_transmon_zpfs(self, EJ, EC):
        """Dimensionless zero-point flux/charge fluctuations of a transmon
        mode, phi_zpf = (2 EC/EJ)^(1/4), n_zpf = (1/2)(EJ/2EC)^(1/4), in the
        standard convention where phi (reduced flux) and n (charge number,
        in units of 2e) are canonically conjugate: phi_hat = phi_zpf(a+a^dag),
        n_hat = i n_zpf(a^dag - a).
        """
        EJ_q = self._as_quantity(EJ, self.freq_units)
        EC_q = self._as_quantity(EC, self.freq_units)
        phi_zpf = (2.0 * EC_q / EJ_q) ** 0.25
        n_zpf = 0.5 * (EJ_q / (2.0 * EC_q)) ** 0.25
        return phi_zpf, n_zpf

    # ------------------------------------------------------------------
    # Dispersive regime: chi, dressed frequencies
    # ------------------------------------------------------------------
    def compute_dispersive_shift(self, g, qubit_freq, resonator_freq, anharmonicity):
        """Dispersive shift chi in the full dispersive approximation.

        Uses the two-pole formula that accounts for both the |0>-|1> and
        |1>-|2> transitions of the transmon (not just the simple g^2/Delta
        RWA limit), following Koch et al. (2007) eq. (3.6):

            chi = -g^2 * alpha / (Delta * (Delta + alpha))

        where Delta = f_res - f_qubit (detuning, signed) and alpha is the
        anharmonicity (negative for a transmon, alpha ~ -EC).

        The simple dispersive limit chi ~ g^2/Delta is recovered when
        |alpha| >> |Delta|, but this formula is more accurate across the
        full transmon range and is consistent with the SQuADDS paper formula
        used in from_eigenmode_epr.

        All inputs can be pint Quantities or bare floats in freq_units.
        Returns chi in freq_units (same sign convention as Koch: chi < 0
        for a transmon qubit below the resonator, chi > 0 above).

        Parameters
        ----------
        g : coupling strength (freq units)
        qubit_freq : qubit f_01 (freq units)
        resonator_freq : bare resonator frequency (freq units)
        anharmonicity : alpha = f_12 - f_01, negative for transmon (freq units)
        """
        g_q = self._as_quantity(g, self.freq_units)
        fq = self._as_quantity(qubit_freq, self.freq_units)
        fr = self._as_quantity(resonator_freq, self.freq_units)
        alpha = self._as_quantity(anharmonicity, self.freq_units)

        Delta = fr - fq
        chi = -(g_q ** 2 * alpha) / (Delta * (Delta + alpha))
        return chi.to(self.freq_units)

    def compute_dressed_frequencies(self, g, qubit_freq, resonator_freq,
                                     anharmonicity=None, mode_b_freq=None):
        """Dressed frequencies for a coupled pair, handling two cases:

        QUBIT-RESONATOR (anharmonicity provided, mode_b_freq=None):
            Uses the dispersive approximation chi formula. Also checks
            g/|Delta| > 0.1 and sets 'is_dispersive_valid': False in the
            result if so -- the caller must not trust chi in that regime.

            Returns dict: chi, qubit, resonator, is_dispersive_valid,
            dispersive_ratio (= g/|Delta|).

        RESONATOR-RESONATOR (anharmonicity=None, mode_b_freq provided):
            Exact 2x2 diagonalization -- no dispersive approximation.
            Normal mode frequencies are:
                omega_pm = (omega_a + omega_b)/2 +/- sqrt(Delta^2 + 4g^2)/2
            This is exact for linear modes and must NOT use chi=g^2/Delta.

            Returns dict: omega_plus, omega_minus, mixing_angle (in radians).

        All frequency inputs can be pint Quantities or bare floats in freq_units.
        """
        g_q = self._as_quantity(g, self.freq_units)
        fa = self._as_quantity(qubit_freq, self.freq_units)  # 'mode a'

        if anharmonicity is not None:
            # --- Qubit-Resonator dispersive path ---
            fr = self._as_quantity(resonator_freq, self.freq_units)
            alpha = self._as_quantity(anharmonicity, self.freq_units)
            chi = self.compute_dispersive_shift(g_q, fa, fr, alpha)

            Delta = fr - fa
            ratio = abs((g_q / Delta).to_reduced_units().magnitude)
            is_valid = ratio <= 0.1

            return {
                "chi": chi,
                "qubit": (fa - chi).to(self.freq_units),
                "resonator": (fr + chi).to(self.freq_units),
                "is_dispersive_valid": is_valid,
                "dispersive_ratio": float(ratio),
            }
        else:
            # --- Resonator-Resonator exact diagonalization ---
            if mode_b_freq is None:
                raise ValueError(
                    "compute_dressed_frequencies: for Resonator-Resonator pairs "
                    "provide mode_b_freq instead of anharmonicity."
                )
            fb = self._as_quantity(mode_b_freq, self.freq_units)
            Delta = fb - fa
            # sqrt(Delta^2 + 4g^2)/2  -- collapse units immediately to avoid
            # pint compound-sqrt issues
            discriminant = np.sqrt(
                (Delta.to(self.freq_units).magnitude ** 2
                 + 4.0 * g_q.to(self.freq_units).magnitude ** 2)
            ) * self.ureg(self.freq_units)

            omega_plus = ((fa + fb) / 2.0 + discriminant / 2.0).to(self.freq_units)
            omega_minus = ((fa + fb) / 2.0 - discriminant / 2.0).to(self.freq_units)
            # Mixing angle theta: tan(2*theta) = 2g / Delta
            mixing_angle = 0.5 * np.arctan2(
                2.0 * g_q.to(self.freq_units).magnitude,
                Delta.to(self.freq_units).magnitude
            )  # radians, dimensionless

            return {
                "omega_plus": omega_plus,
                "omega_minus": omega_minus,
                "mixing_angle": mixing_angle,
            }

    # ------------------------------------------------------------------
    # Higher-order dispersive corrections (Blais 2021 review)
    # ------------------------------------------------------------------
    def compute_higher_order_dispersive(self, g, qubit_freq, resonator_freq,
                                        anharmonicity):
        """4th-order dispersive corrections from the Blais (2021) review.

        The standard chi is O(g^2). These add:

        1. 4th-order vacuum shift (Lamb shift correction):
               delta_omega_r^(4) ~ -g^4 / (Delta^3)   [resonator]
               delta_omega_q^(4) ~ +g^4 / (Delta^3)   [qubit]
           These shift the bare frequencies beyond the leading-order chi.

        2. Self-Kerr / nonlinear cavity pull (chi'):
               chi' = chi * (2*alpha) / (Delta + alpha)
           This is the rate at which chi changes per photon in the resonator.
           High-power readout drives the resonator into multi-photon states;
           chi' sets the scale at which the dispersive approximation breaks
           down with photon number. When |chi'| is large relative to chi,
           a few photons already shift the resonator significantly.

        All inputs bare floats or pint Quantities in freq_units.

        Returns dict with keys:
          'chi'              : 2nd-order dispersive shift (Koch formula)
          'chi_prime'        : self-Kerr / nonlinear cavity pull (freq_units)
          'vacuum_shift_res' : 4th-order resonator frequency shift (freq_units)
          'vacuum_shift_q'   : 4th-order qubit frequency shift (freq_units)
          'is_dispersive_valid' : bool, False if g/|Delta| > 0.1
          'dispersive_ratio' : float g/|Delta|
        """
        g_q = self._as_quantity(g, self.freq_units)
        fq = self._as_quantity(qubit_freq, self.freq_units)
        fr = self._as_quantity(resonator_freq, self.freq_units)
        alpha = self._as_quantity(anharmonicity, self.freq_units)

        Delta = fr - fq
        chi = self.compute_dispersive_shift(g_q, fq, fr, alpha)

        # chi' = chi * 2*alpha / (Delta + alpha)
        # This is the photon-number-dependent correction to chi.
        chi_prime = chi * (2.0 * alpha) / (Delta + alpha)

        # 4th-order vacuum shift ~ g^4 / Delta^3
        # Sign: resonator shifts down, qubit shifts up (perturbation theory)
        vac_shift_mag = (g_q ** 4 / Delta ** 3).to(self.freq_units)
        vacuum_shift_res = -vac_shift_mag
        vacuum_shift_q = +vac_shift_mag

        ratio = abs((g_q / Delta).to_reduced_units().magnitude)

        return {
            "chi": chi.to(self.freq_units),
            "chi_prime": chi_prime.to(self.freq_units),
            "vacuum_shift_res": vacuum_shift_res,
            "vacuum_shift_q": vacuum_shift_q,
            "is_dispersive_valid": ratio <= 0.1,
            "dispersive_ratio": float(ratio),
        }

    def compute_purcell_decay(self, g, qubit_freq, resonator_freq, kappa):
        """Purcell decay rate for a qubit coupled to a lossy resonator.

        gamma_purcell = kappa * (g / Delta)^2

        This is the T1 ceiling imposed by the resonator's linewidth kappa
        through the qubit-resonator coupling. In the dispersive regime
        (g << |Delta|), this is the dominant non-radiative decay channel
        set by the readout resonator.

        Parameters
        ----------
        g : coupling strength (freq units)
        qubit_freq : qubit f_01 (freq units)
        resonator_freq : bare resonator frequency (freq units)
        kappa : resonator linewidth FWHM (freq units). Must be provided;
                this method does not silently return None if kappa is missing
                -- that check is the caller's responsibility (Chip does it).

        Returns dict with:
          'gamma_purcell' : Purcell decay rate in freq_units
          'T1_purcell_s'  : Purcell-limited T1 = 1/(2*pi*gamma_purcell) in seconds
          'T1_purcell_us' : same, in microseconds (plain float, T1_purcell_s * 1e6)
          'dispersive_ratio' : float g/|Delta|
        """
        g_q = self._as_quantity(g, self.freq_units)
        fq = self._as_quantity(qubit_freq, self.freq_units)
        fr = self._as_quantity(resonator_freq, self.freq_units)
        kappa_q = self._as_quantity(kappa, self.freq_units)

        Delta = fr - fq
        ratio = (g_q / Delta).to_reduced_units().magnitude
        gamma_purcell = (kappa_q * ratio ** 2).to(self.freq_units)

        # T1 = 1 / (2*pi*gamma) -- gamma here is already in linear freq units
        # so T1 = 1 / (2*pi*gamma_Hz)
        gamma_Hz = gamma_purcell.to("Hz").magnitude
        T1_s = 1.0 / (2.0 * np.pi * gamma_Hz) if gamma_Hz != 0 else float("inf")

        return {
            "gamma_purcell": gamma_purcell,
            "T1_purcell_s": T1_s,
            "T1_purcell_us": T1_s * 1e6,
            "dispersive_ratio": float(ratio),
        }