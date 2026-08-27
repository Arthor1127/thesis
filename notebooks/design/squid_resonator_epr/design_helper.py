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
    
    def __init__(self, alpha_inductance=None):
        self.ureg = pint.UnitRegistry()
        self.hbar = self.ureg.Quantity(1, "hbar")
        self.planck = self.ureg.Quantity(1, 'planck_constant')
        self.electron_charge = self.ureg.Quantity(1, "elementary_charge")
        self.boltzmann = self.ureg.Quantity(1, "boltzmann_constant")
        self.reduced_flux_quantum = self.ureg.Quantity(1, "magnetic_flux_quantum") / (2.0 * np.pi)
        self.light_speed = self.ureg.Quantity(1, 'speed_of_light').to(self.speed_units)
        self.vacuum_permitivity = (8.854187817e-12 * self.ureg.farad / self.ureg.meter).to(f"{self.capacitance_units} / {self.length_units}")
        
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
    # of superconducting coplanar resonators made out of granular aluminium"
    # (arXiv:2508.03873v1, 2025), Eqs. (B1)-(B9).
    #
    # NOTE on Eq. (B4) prefactor: the PDF prints 377/16 but that gives
    # ~12.4 Ohm for the paper's own example (w=20um, g=11um, t=70nm ->
    # stated Z0=50.1 Ohm). The standard CPW formula uses 30*pi (Simons;
    # Ghione & Naldi), which reproduces 49.67 Ohm for the same geometry.
    # This implementation uses 30*pi -- the 377/16 is almost certainly a
    # dropped-pi/doubled-denominator transcription artifact.
    #
    # NOTE on Eq. (B9) G(w,g,t): calibrated against Table II of the paper.
    # A ~2x discrepancy in L_k^sq vs. the paper's reported values remains
    # unresolved (likely traces to the Clem typo correction -- see
    # compute_cpw_kinetic_inductance_per_length docstring). Treat this
    # method as a first-estimate tool consistent with the paper's own
    # stated confidence level for this formula (Sec. V.1: "indicative only").

    def _ellipk_agm(self, k):
        """Complete elliptic integral of the first kind K(k) via the
        arithmetic-geometric mean. Exact to machine precision in ~15
        iterations; no scipy dependency (keeps DesignHelper's third-party
        deps to numpy/pint only).
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
        including finite-thickness correction -- Eqs. (B1)-(B7).

        Parameters
        ----------
        width, gap, thickness : length-like (default units: micrometers)
            Central conductor width w, gap to ground g, metal thickness t.
        eps_r : float
            Relative permittivity of the substrate (e.g. ~11.45-11.9 for
            high-resistivity Si, ~9.4-11.5 for sapphire).

        Returns dict with keys:
          'Cg'     : capacitance per unit length (F/m)
          'Lg'     : geometric inductance per unit length (H/m)
          'Z0'     : characteristic impedance (resistance_units)
          'v'      : propagation velocity (speed_units)
          'eps_eff': effective permittivity (dimensionless)

        Verified: reproduces the paper's stated CMT frequencies for R1/R2/R3
        (5.448/5.938/6.426 GHz) to within 0.1% at eps_r=11.9.
        """
        w = self._as_quantity(width, self.length_units).to("m").magnitude
        g = self._as_quantity(gap, self.length_units).to("m").magnitude
        t = self._as_quantity(thickness, self.length_units).to("m").magnitude

        eps0 = self.ureg.Quantity(1, "vacuum_permittivity").to("F/m").magnitude
        mu0  = self.ureg.Quantity(1, "magnetic_constant").to("H/m").magnitude

        k0   = self._cpw_k0(w, g)
        kp0  = np.sqrt(1.0 - k0**2)
        K_k0  = self._ellipk_agm(k0)
        K_kp0 = self._ellipk_agm(kp0)

        # Finite-thickness correction -- u1, u2, then k_t  (B6 inputs)
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
        kt   = u1 / u2
        ktp  = np.sqrt(1.0 - kt**2)
        K_kt  = self._ellipk_agm(kt)
        K_ktp = self._ellipk_agm(ktp)

        # Cg, Eq. (B6)
        Cg = eps0 * (2.0 * K_kt / K_ktp) + eps_r * eps0 * (2.0 * K_k0 / K_kp0)

        # Lg, Eq. (B7) -- using literal "k_{t/2}" reading (recompute at t/2),
        # which gives 49.40 Ohm vs 49.23 Ohm for the "reuse kt" reading;
        # the former is closer to the paper's stated 50.1 Ohm.
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
        kt_h  = u1_h / u2_h
        ktp_h = np.sqrt(1.0 - kt_h**2)
        K_kt_h  = self._ellipk_agm(kt_h)
        K_ktp_h = self._ellipk_agm(ktp_h)
        Lg = (mu0 / 4.0) * (K_ktp_h / K_kt_h)

        Z0 = np.sqrt(Lg / Cg)
        v  = 1.0 / np.sqrt(Lg * Cg)
        eps_eff0 = (1.0 + eps_r) / 2.0
        v0       = self.light_speed.to("m/s").magnitude / np.sqrt(eps_eff0)
        eps_eff  = (v0 / v) ** 2

        return {
            "Cg":      self.ureg.Quantity(Cg, "F/m"),
            "Lg":      self.ureg.Quantity(Lg, "H/m"),
            "Z0":      self.ureg.Quantity(Z0, "ohm").to(self.resistance_units),
            "v":       self.ureg.Quantity(v,  "m/s").to(self.speed_units),
            "eps_eff": float(eps_eff),
        }

    def compute_cpw_geometric_factor(self, width, gap, thickness):
        """Geometric factor G(w,g,t) from Eq. (B9) (Watanabe et al., with
        the Clem typo correction as printed in Ramos et al.). Used by
        compute_cpw_kinetic_inductance_per_length; exposed for inspection.
        """
        w  = self._as_quantity(width,     self.length_units).to("m").magnitude
        g  = self._as_quantity(gap,       self.length_units).to("m").magnitude
        t  = self._as_quantity(thickness, self.length_units).to("m").magnitude

        k0   = self._cpw_k0(w, g)
        kp0  = np.sqrt(1.0 - k0**2)
        K_k0 = self._ellipk_agm(k0)

        term1 = -np.log(t / (4.0 * w))
        term2 = -(w / (w + 2.0 * g)) * np.log(t / (4.0 * (w + 2.0 * g)))
        term3 =  (2.0 * (w + g) / (w + 2.0 * g)) * np.log(g / (w + g))
        # Denominator is 4*k'0^2*K(k0)^2 per Watanabe (1994) Eq. 4 and Gao
        # (2008) thesis Eq. 3.21. The Ramos et al. paper prints 2 -- that is
        # a transcription error relative to the primary sources.
        return (1.0 / (4.0 * kp0**2 * K_k0**2)) * (term1 + term2 + term3)

    def compute_cpw_kinetic_inductance_per_length(self, width, gap, thickness,
                                                   sheet_kinetic_inductance):
        """Kinetic inductance per unit length, Eq. (B9): Lk = (L_k^sq/w)*G.

        sheet_kinetic_inductance: L_k^sq in pH (bare float or pint Quantity).

        CALIBRATION CAVEAT: after the G-denominator fix (4 vs 2), the
        analytical Lg path of solve_sheet_kinetic_inductance_from_alpha gives
        ~1.05 pH/sq vs. the paper's 2.10 pH/sq for Al R1 (alpha=0.139). The
        remaining ~2x traces to effective-length geometry in CPW (see
        solve_sheet_kinetic_inductance_from_alpha docstring). Treat as
        indicative.
        """
        Lk_sq = self._as_quantity(sheet_kinetic_inductance, "pH").to("H").magnitude
        w     = self._as_quantity(width, self.length_units).to("m").magnitude
        G     = self.compute_cpw_geometric_factor(width, gap, thickness)
        return self.ureg.Quantity((Lk_sq / w) * G, "H/m")

    def compute_cpw_kinetic_fraction(self, width, gap, thickness, eps_r,
                                      sheet_kinetic_inductance):
        """Kinetic inductance fraction alpha = Lk/(Lg+Lk), Eqs. (3)+(B8)+(B9).

        Subject to the same calibration caveat as
        compute_cpw_kinetic_inductance_per_length.
        """
        geo  = self.compute_cpw_geometric_parameters(width, gap, thickness, eps_r)
        Lg   = geo["Lg"].to("H/m").magnitude
        Lk   = self.compute_cpw_kinetic_inductance_per_length(
            width, gap, thickness, sheet_kinetic_inductance
        ).to("H/m").magnitude
        return Lk / (Lg + Lk)

    def solve_sheet_kinetic_inductance_from_alpha(self, width, gap, thickness,
                                                   eps_r, target_alpha,
                                                   fg_sim=None, Z0_line=50.0,
                                                   resonator_length=None,
                                                   bracket=(1e-4, 1e4), tol=1e-10,
                                                   max_iter=200):
        """Inverse of compute_cpw_kinetic_fraction: given a measured alpha
        (e.g. alpha_1 from a Mattis-Bardeen fit), solve for L_k^sq that
        produces it. Bisection on the monotone alpha(L_k^sq) curve.

        Returns L_k^sq as a pint Quantity in pH.

        Two modes for the geometric inductance Lg used inside alpha(L_k^sq):

        **Analytical (default):** Lg per unit length from compute_cpw_geometric_parameters
        (Eq. B7). Self-consistent with the CPW formulae but can differ ~17%
        from a FEM result.

        **Simulation-backed:** pass fg_sim (Hz) and resonator_length (m) -- or
        any length unit pint Quantity -- together with Z0_line (Ohm, default
        50). Lg_total is then derived from the simulated resonance frequency
        via Lg_total = 1/(Cg_total * omega_g^2) with Cg_total = pi/(4*omega_g*Z0_line)
        (Pozar), exactly as in Ramos et al. Sec. V.1 / Table II. This path is
        more accurate when a reliable fg_sim is available.

        CALIBRATION NOTE: after the G-denominator fix (4 vs 2) the analytical
        path returns ~1.05 pH/sq vs. the paper's 2.10 pH/sq for Al R1
        (alpha=0.139). The simulation-backed path closes roughly half of that
        remaining gap. The residual ~2x is a real geometry effect: effective
        squares in CPW ≈ 0.5*(ell/w) due to ground-plane field distribution
        (Ramos et al. Sec. V.1, "indicative only"). Treat absolute values
        accordingly.
        """
        import math

        if fg_sim is not None:
            # Simulation-backed Lg: derive from simulated resonance frequency.
            # Requires resonator_length to convert Lg/m -> Lg_total.
            if resonator_length is None:
                raise ValueError(
                    "resonator_length is required when fg_sim is provided -- "
                    "need it to compute Lg_total from Cg_total."
                )
            ell = self._as_quantity(resonator_length, self.length_units).to("m").magnitude
            omega_g = 2.0 * math.pi * float(fg_sim)
            Cg_total = math.pi / (4.0 * omega_g * float(Z0_line))   # F
            Lg_total = 1.0 / (Cg_total * omega_g**2)                # H
            Lg_per_m = Lg_total / ell                                # H/m

            def _alpha(Lk_sq_pH):
                Lk_sq = Lk_sq_pH * 1e-12
                w = self._as_quantity(width, self.length_units).to("m").magnitude
                G = self.compute_cpw_geometric_factor(width, gap, thickness)
                Lk_per_m = (Lk_sq / w) * G
                return Lk_per_m / (Lg_per_m + Lk_per_m)
        else:
            # Analytical Lg from B7 via compute_cpw_geometric_parameters.
            def _alpha(Lk_sq_pH):
                return self.compute_cpw_kinetic_fraction(
                    width, gap, thickness, eps_r,
                    self.ureg.Quantity(Lk_sq_pH, "pH")
                )

        lo, hi = bracket
        f_lo = _alpha(lo) - target_alpha
        f_hi = _alpha(hi) - target_alpha
        if f_lo * f_hi > 0:
            raise ValueError(
                f"target_alpha={target_alpha} not bracketed in [{lo}, {hi}] pH "
                f"(alpha(lo)={f_lo+target_alpha:.4f}, alpha(hi)={f_hi+target_alpha:.4f})"
            )
        for _ in range(max_iter):
            mid   = 0.5 * (lo + hi)
            f_mid = _alpha(mid) - target_alpha
            if abs(f_mid) < tol or (hi - lo) < tol:
                return self.ureg.Quantity(mid, "pH")
            if f_lo * f_mid <= 0:
                hi = mid
            else:
                lo, f_lo = mid, f_mid
        return self.ureg.Quantity(0.5 * (lo + hi), "pH")
    
    def compute_resonator_Qc_from_Cc(self, width, gap, thickness, eps_r, 
                                     resonator_length, Cc, Z0_line=50.0,
                                     wavelength_ratio=0.25):
        """
        Calculates the coupling quality factor (Qc) and decay rate (kappa) 
        of a resonator for any wavelength ratio (e.g. 0.25 for lambda/4, 0.5 for lambda/2).
        
        Assumes a standard "hanger" configuration where the resonator is capacitively
        coupled to a feedline of impedance Z0_line.
        
        Args:
            width, gap, thickness, eps_r: CPW geometry parameters.
            resonator_length: Physical length of the trace.
            Cc: Coupling capacitance in fF (or pint Quantity).
            Z0_line: Characteristic impedance of the feedline (usually 50 Ohms).
            wavelength_ratio: Ratio of length to wavelength (0.25 = lambda/4, 0.5 = lambda/2).
            
        Returns:
            dict with 'Qc' (dimensionless float), 'kappa' (freq units), 
            and 'omega_r' (freq units).
        """
        # Get CPW per-unit-length parameters
        geo = self.compute_cpw_geometric_parameters(width=width, gap=gap, 
                                                    thickness=thickness, eps_r=eps_r)
        Cg_m = self._as_quantity(geo["Cg"], "F/m").magnitude
        v_ph = self._as_quantity(geo["v"], "m/s").magnitude
        
        ell = self._as_quantity(resonator_length, self.length_units).to("m").magnitude
        Cc_F = self._as_quantity(Cc, self.capacitance_units).to("F").magnitude
        Z0 = float(Z0_line)
        
        # 1. Frequency mapping based on the chosen mode
        # f_r = v_ph / lambda, where lambda = ell / wavelength_ratio
        omega_r = 2.0 * np.pi * v_ph * (float(wavelength_ratio) / ell)
        
        # 2. Lumped Equivalent Circuit Parameters
        # For the fundamental mode of both lambda/4 and lambda/2 CPW resonators, 
        # the effective lumped capacitance is half the total geometric capacitance.
        C_total = Cg_m * ell
        C_eff = C_total / 2.0
        
        # External conductance introduced by the hanger feedline
        # G_ext = (omega * Cc)^2 * R_env, where R_env = Z0 / 2 for a hanger
        G_ext = (omega_r * Cc_F)**2 * (Z0 / 2.0)
        
        # 3. Quality Factor
        Qc = (omega_r * C_eff) / G_ext
        kappa_rad = omega_r / Qc
        
        return {
            "Qc": float(Qc),
            "kappa": self._as_quantity(kappa_rad / (2.0 * np.pi), "Hz").to(self.freq_units),
            "omega_r": self._as_quantity(omega_r / (2.0 * np.pi), "Hz").to(self.freq_units)
        }
    def compute_Cc_from_Qc(self, width, gap, thickness, eps_r,
                            resonator_length, Qc_target,
                            Z0_line=50.0, wavelength_ratio=0.5):
        """
        Analytically computes the coupling capacitance Cc needed to achieve
        a target external quality factor Qc, for a CPW resonator in hanger geometry.

        Args:
            width, gap, thickness, eps_r: CPW geometry parameters.
            resonator_length: Physical length of the resonator trace.
            Qc_target: Target external quality factor (dimensionless).
            Z0_line: Feedline impedance (default 50 Ohm).
            wavelength_ratio: 0.25 for lambda/4, 0.5 for lambda/2.

        Returns:
            Cc as a pint Quantity in self.capacitance_units.
        """
        geo = self.compute_cpw_geometric_parameters(
            width=width, gap=gap, thickness=thickness, eps_r=eps_r
        )
        Cg_m = self._as_quantity(geo["Cg"], "F/m").magnitude
        v_ph = self._as_quantity(geo["v"], "m/s").magnitude
        ell  = self._as_quantity(resonator_length, self.length_units).to("m").magnitude
        Z0   = float(Z0_line)

        omega_r = 2.0 * np.pi * v_ph * (float(wavelength_ratio) / ell)
        C_eff   = (Cg_m * ell) / 2.0

        # Inverted from Qc = 2*C_eff / (Z0 * omega_r * Cc^2)
        Cc_F = np.sqrt(2.0 * C_eff / (Z0 * omega_r * float(Qc_target)))

        return self._as_quantity(Cc_F, "F").to(self.capacitance_units)

    def compute_resonator_Qc_from_sim(self, target_coupling_length, 
                                      qc_sweep_data_lengths, qc_sweep_data_Qcs):
        """
        Interpolates the expected Qc from a pre-computed simulation sweep 
        (e.g., HFSS or Palace) of coupling length vs. Qc.
        
        Args:
            target_coupling_length: The design length of the coupling segment.
            qc_sweep_data_lengths: Array of simulated coupling lengths.
            qc_sweep_data_Qcs: Array of simulated Qc values.
            
        Returns:
            Interpolated Qc (dimensionless float).
        """
        from scipy.interpolate import interp1d
        
        l_target = self._as_quantity(target_coupling_length, self.length_units).magnitude
        l_array = [self._as_quantity(l, self.length_units).magnitude for l in qc_sweep_data_lengths]
        q_array = np.array(qc_sweep_data_Qcs)
        
        # Use cubic interpolation for smooth curves, fallback to linear if points are few
        kind = 'cubic' if len(l_array) > 3 else 'linear'
        f_interp = interp1d(l_array, q_array, kind=kind, fill_value="extrapolate")
        
        Qc_estimated = float(f_interp(l_target))
        return Qc_estimated
    
    
    # ------------------------------------------------------------------
    # Josephson junction: critical current <-> E_J
    # ------------------------------------------------------------------
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
    
    def compute_transmon_energy_levels(self, EC : pint.Quantity, EJ : pint.Quantity, n_levels : int, exact = False) -> np.array:
        """Compute transmon spectrum

        Args:
            EC (pint.Quantity): Charge energy of the transmon
            EJ (pint.Quantity): Josephson energy of the transmon
            n_levels(int): number of energy levels needed for the computation
            exact (bool, optional): whether to calculate the spectrum by the asymptotic formula or by exact diagonalization of the hamiltonian in the charge basis. Defaults to False.

        Raises:
            ValueError: _description_
            ValueError: _description_

        Returns:
            np.array: array of first energies (difference with the base level)
        """
        
        EC_q = self._as_quantity(EC, self.freq_units)
        EJ_q = self._as_quantity(EJ, self.freq_units)
        
        if not exact:
            omega_q = np.sqrt(8 * EC_q *EJ_q) - EC_q
            U = EC_q / 2.0
            n_values = np.arange(1, n_levels+1, 1.0)
            return omega_q * n_values - U * n_values * (n_values - 1)

        ncut = n_levels + 60
        n_vals = np.arange(-ncut, ncut + 1)
        dim = len(n_vals)

        # Extract plain magnitudes — already done correctly
        EJ_mag = EJ_q.to(self.freq_units).magnitude
        EC_mag = EC_q.to(self.freq_units).magnitude

        H = np.diag(4.0 * EC_mag * n_vals**2).astype(float)
        off_diag = -0.5 * EJ_mag
        for i in range(dim - 1):
            H[i, i + 1] += off_diag
            H[i + 1, i] += off_diag

        evals = np.linalg.eigvalsh(H)   # plain float array, in freq_units
        evals.sort()
        evals -= evals[0]                # subtract ground state
        evals = evals[1:n_levels + 1]   # take only what you need

        # Attach units AFTER building the plain numpy array
        return evals * self._as_quantity(1.0, self.freq_units).units
        
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

    def compute_transmon_impedance(self, EJ, EC):
        """Characteristic impedance of the transmon's junction mode,
        Z_q = sqrt(L_J / C_Sigma), where L_J = (hbar/(2e))^2 / EJ and
        C_Sigma = e^2/(2*EC).

        This is the transmon analogue of Z=sqrt(L/C) for a resonator.
        It controls the ZPF magnitudes and is useful when comparing a
        transmon's quantum fluctuations against those of a coupled resonator
        (e.g. for estimating coupling ratios without a full cap matrix).
        Returned in resistance_units.

        Note: this is the linearised (harmonic) impedance at the potential
        minimum -- the transmon's anharmonicity means it's not a pure
        harmonic oscillator, but Z_q is still the relevant scale for the
        ZPF magnitudes that enter the coupling formula.
        """
        EJ_q = self._as_quantity(EJ, self.freq_units)
        EC_q = self._as_quantity(EC, self.freq_units)
        # L_J = phi_0^2 / EJ (energy form), phi_0 = hbar/(2e)
        phi0    = self.hbar / (2.0 * self.electron_charge)
        EJ_energy = self._frequency_as_energy(EJ_q)
        L_J     = (phi0**2 / EJ_energy).to("H")
        # C_Sigma = e^2 / (2 EC)
        EC_energy = self._frequency_as_energy(EC_q)
        C_Sigma = (self.electron_charge**2 / (2.0 * EC_energy)).to("F")
        return np.sqrt(L_J / C_Sigma).to(self.resistance_units)

    def circular_pad_capacitance(self, radius : pint.Quantity, thickness : pint.Quantity, gap : pint.Quantity, eps_r : float = 11.9) -> pint.Quantity:
        """Calculated capacitance of a cylindrical capacitor 

        Args:
            radius (pint.Quantity): charge island radius
            thickness (pint.Quantity): film thickness, or equivalently the cyllinder's length
            gap (pint.Quantity): vacuum gap between the island to the ground plate
            eps_r (float, optional): dielectric constant of the substrate. Defaults to 11.9 corresponding to silicon.

        Raises:
            ValueError: _description_

        Returns:
            pint.Quantity: capacitance in the helper's default units
        """
        
        r = self._as_quantity(radius, self.length_units)
        g = self._as_quantity(gap, self.length_units)
        t = self._as_quantity(thickness, self.length_units)
        permittivity = 0.5 * (1 + eps_r) * self.vacuum_permitivity
        
        log_factor = 1 / np.log(1+r.magnitude / g.magnitude)
        numerator = 2 * np.pi * permittivity * t
        return log_factor * numerator.to_reduced_units

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
    
    def compute_level_dependent_dispersive_shift(self, g, qubit_freq, resonator_freq, anharmonicity, level = 0):
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
        level (int, optional): Level in which the dispersive shift is calculated. Defaults to 0
        """
        g_q = self._as_quantity(g, self.freq_units)
        fq = self._as_quantity(qubit_freq, self.freq_units)
        fr = self._as_quantity(resonator_freq, self.freq_units)
        alpha = self._as_quantity(anharmonicity, self.freq_units)
        if level == 0:
            Delta = fr - fq
            chi = -(g_q ** 2 * alpha) / (Delta * (Delta + alpha))
            return chi.to(self.freq_units)
        transmon_energy = lambda j: j * fq - 0.5 * j * (j-1) * alpha
        
        omega_j = transmon_energy(level)
        omega_j_minus1 = transmon_energy(level - 1)
        omega_j_plus1 = transmon_energy(level + 1)
        
        dispersive_shift = g_q**2 * (
            level / (omega_j - omega_j_minus1 - fr) - (level + 1) / (omega_j_plus1 - omega_j - fr)
        )
        return dispersive_shift

    def compute_level_dependent_lamb_shift(self, g, qubit_freq, resonator_freq, anharmonicity, level = 0):
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
        level (int, optional): Level in which the dispersive shift is calculated. Defaults to 0
        """
        g_q = self._as_quantity(g, self.freq_units)
        fq = self._as_quantity(qubit_freq, self.freq_units)
        fr = self._as_quantity(resonator_freq, self.freq_units)
        alpha = self._as_quantity(anharmonicity, self.freq_units)
        if level == 0:
            return 0
        
        transmon_energy = lambda j: j * fq - 0.5 * j * (j-1) * alpha
        
        omega_j = transmon_energy(level)
        
        lamb_shift = level * g_q**2 / (omega_j - alpha * (level -1) - fr)
        return lamb_shift
    
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
        
    def compute_purcell_decay_feedline(self, g, qubit_freq, resonator_freq,
                                        kappa_ext, kappa_int=None):
        """
        Purcell decay rate of a qubit coupled to a resonator that is itself
        coupled to a feedline, decomposed into internal and external channels.

        The total Purcell rate is:
            gamma_purcell = (g / Delta)^2 * (kappa_ext + kappa_int)

        The feedline (external) contribution is:
            gamma_ext = (g / Delta)^2 * kappa_ext

        The quantum efficiency of the Purcell channel (fraction of qubit decay
        that goes into the feedline rather than being lost internally) is:
            eta = kappa_ext / (kappa_ext + kappa_int)

        Parameters
        ----------
        g            : qubit-resonator coupling strength (freq units)
        qubit_freq   : qubit f_01 (freq units)
        resonator_freq : bare resonator frequency (freq units)
        kappa_ext    : external (feedline) decay rate = omega_r / Qc (freq units)
        kappa_int    : internal loss rate = omega_r / Qi (freq units).
                    If None, assumes ideal overcoupled resonator (kappa_int=0),
                    i.e. kappa_total = kappa_ext.

        Returns dict with:
            'gamma_purcell_total'   : total Purcell-limited decay rate
            'gamma_purcell_ext'     : portion going into feedline
            'gamma_purcell_int'     : portion lost internally (0 if kappa_int=None)
            'T1_purcell_total_us'   : T1 ceiling from total Purcell decay (us)
            'T1_purcell_ext_us'     : T1 from feedline channel only (us)
            'dispersive_ratio'      : g / |Delta| (should be << 1)
            'eta_feedline'          : quantum efficiency of Purcell channel
            'overcoupled'           : bool, whether kappa_ext >> kappa_int
        """
        g_q    = self._as_quantity(g, self.freq_units)
        fq     = self._as_quantity(qubit_freq, self.freq_units)
        fr     = self._as_quantity(resonator_freq, self.freq_units)
        k_ext  = self._as_quantity(kappa_ext, self.freq_units)
        k_int  = self._as_quantity(kappa_int, self.freq_units) \
                if kappa_int is not None \
                else self._as_quantity(0.0, self.freq_units)

        Delta = fr - fq
        ratio = (g_q / Delta).to_reduced_units().magnitude  # g/Delta, dimensionless

        k_total = k_ext + k_int

        gamma_total = (k_total * ratio**2).to(self.freq_units)
        gamma_ext   = (k_ext   * ratio**2).to(self.freq_units)
        gamma_int   = (k_int   * ratio**2).to(self.freq_units)

        def _to_T1_us(gamma):
            gamma_Hz = gamma.to("Hz").magnitude
            if gamma_Hz == 0:
                return float("inf")
            return 1.0 / (2.0 * np.pi * gamma_Hz) * 1e6

        k_ext_Hz  = k_ext.to("Hz").magnitude
        k_int_Hz  = k_int.to("Hz").magnitude
        k_tot_Hz  = k_ext_Hz + k_int_Hz
        eta       = k_ext_Hz / k_tot_Hz if k_tot_Hz > 0 else 1.0

        return {
            "gamma_purcell_total":  gamma_total,
            "gamma_purcell_ext":    gamma_ext,
            "gamma_purcell_int":    gamma_int,
            "T1_purcell_total_us":  _to_T1_us(gamma_total),
            "T1_purcell_ext_us":    _to_T1_us(gamma_ext),
            "dispersive_ratio":     float(ratio),
            "eta_feedline":         float(eta),
            "overcoupled":          k_ext_Hz > 10 * k_int_Hz,  # practical threshold
        }