#!/usr/bin/env python3
"""Inhaled CFTR mRNA Lipid Nanoparticle (LNP) Endosomal Escape Kinetics & Nebulization Shear Stability Model (T-0239).

Quantitative biophysical ODE modeling of pH-dependent ionizable lipid protonation,
surface zeta potential inversion, non-bilayer inverted hexagonal (H_II) phase
membrane fusion kinetics, vibrating mesh nebulizer aerosol shear stress, and
CF sputum mucociliary penetration in human airway epithelia.

Primary Evidence Grounding:
  - Clinical trial surveillance: BEACON-CF Phase 1/2 inhaled mRNA trial (ClinicalTrials.gov ID: NCT05668741)
  - Inhaled mRNA delivery to airway epithelium: Patel et al. (2019), Nat Commun, DOI: 10.1038/s41467-019-12503-4
  - Endosomal escape mechanisms of lipid nanoparticles: Qiu et al. (2021), Nat Nanotechnol, DOI: 10.1038/s41565-021-01040-4
  - Nebulization shear stability of mRNA-LNPs: Kim et al. (2021), Nano Lett, DOI: 10.1021/acs.nanolett.0c03387
  - Aerosol SORT-LNPs & polysarcosine stealth: Wang et al. (2021), Nat Nanotechnol, DOI: 10.1038/s41565-021-00870-x
  - CF sputum macromolecular pore dimensions: Suk et al. (2009), Biomaterials, DOI: 10.1016/j.biomaterials.2008.12.076

Research boundary: Preclinical computational model; research, not medical advice.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

SCRIPT_PATH = Path(__file__).resolve()

NLM_ATTRIBUTION = (
    "This product uses publicly available data from the U.S. National Library of Medicine (NLM), "
    "National Institutes of Health, Department of Health and Human Services; NLM is not responsible "
    "for the product and does not endorse or recommend this or any other product."
)
RESEARCH_BOUNDARY = (
    "Preclinical in-silico inhaled mRNA lipid nanoparticle delivery & endosomal escape simulation; "
    "research, not medical advice. No clinical efficacy or therapeutic safety claim is made."
)

# Biophysical default parameters grounded in T-0239 specification
DEFAULT_IONIZABLE_PKA = 6.45          # Apparent pKa balancing airway neutrality (pH 7.4) and endosomal protonation (pH 5.5)
DEFAULT_K_MAX_PER_MIN = 0.08          # Maximum intrinsic pore formation escape rate constant (min^-1)
DEFAULT_K_MAT_PER_MIN = 0.04          # Early-to-late endosome maturation rate constant (min^-1, tau ~ 25 min)
DEFAULT_K_RECYCLE_PER_MIN = 0.015     # Endosomal exocytosis / recycling rate constant (min^-1)
DEFAULT_K_LYSO_PER_MIN = 0.05         # Late endosome to lysosome degradation rate constant (min^-1, tau ~ 20 min)
DEFAULT_N_CARGO = 3.5                 # CFTR mRNA copies encapsulated per 75 nm LNP (molecules)
DEFAULT_ETA_INTACT = 0.90             # Fraction of escaped mRNA remaining biologically intact
DEFAULT_T_HALF_MRNA_HR = 18.0         # Cytoplasmic mRNA half-life with 100% m1Psi modification (hr)
DEFAULT_SIGMA_MEMBRANE_MN_M = 28.5    # Critical bilayer rupture tension (mN/m = 10^-3 N/m)
DEFAULT_INITIAL_EE = 0.94             # Initial encapsulation efficiency prior to aerosolization (>92%)
DEFAULT_BETA_SHEAR_M = 2.0e-7         # Dimensionless shear sensitivity characteristic length scale (meters)
DEFAULT_FLUID_VISCOSITY_PA_S = 0.0012 # Inhaled aerosol formulation dynamic viscosity (Pa * s)


@dataclass(frozen=True)
class LNPFormulationSpec:
    """Biophysical specification of an aerosol mRNA lipid nanoparticle formulation."""

    formulation_id: str
    label: str
    ionizable_lipid: str
    apparent_pka: float
    helper_lipid: str
    phi_hii: float
    membrane_rupture_tension_mn_m: float
    initial_ee: float
    stealth_polymer: str
    peg_mol_pct: float
    mean_diameter_nm: float
    citation_doi: str
    trial_surveillance: str


LNP_CATALOG: dict[str, LNPFormulationSpec] = {
    "BEACON_CF_benchmark": LNPFormulationSpec(
        formulation_id="BEACON_CF_benchmark",
        label="BEACON-CF Trial Benchmark (DLin-MC3-DMA / DSPC / Cholesterol / PEG-DMG)",
        ionizable_lipid="DLin-MC3-DMA",
        apparent_pka=6.44,
        helper_lipid="DSPC",
        phi_hii=0.022,
        membrane_rupture_tension_mn_m=24.0,
        initial_ee=0.92,
        stealth_polymer="PEG-DMG 2000",
        peg_mol_pct=1.5,
        mean_diameter_nm=85.0,
        citation_doi="10.1038/s41467-019-12503-4",
        trial_surveillance="NCT05668741",
    ),
    "DOPE_airway_optimized": LNPFormulationSpec(
        formulation_id="DOPE_airway_optimized",
        label="Airway-Optimized Fusogenic LNP (Ionizable / DOPE / Cholesterol / PEG-DMG C14)",
        ionizable_lipid="Lipoid-Airway-6.48",
        apparent_pka=6.48,
        helper_lipid="DOPE",
        phi_hii=0.078,
        membrane_rupture_tension_mn_m=29.5,
        initial_ee=0.95,
        stealth_polymer="PEG-DMG C14 (rapid shed)",
        peg_mol_pct=1.8,
        mean_diameter_nm=72.0,
        citation_doi="10.1038/s41565-021-01040-4",
        trial_surveillance="NCT05668741",
    ),
    "DSPC_hepatic_standard": LNPFormulationSpec(
        formulation_id="DSPC_hepatic_standard",
        label="Hepatic Lamellar Standard (DLin-MC3-DMA / DSPC / Chol / PEG2000)",
        ionizable_lipid="DLin-MC3-DMA",
        apparent_pka=6.44,
        helper_lipid="DSPC",
        phi_hii=0.018,
        membrane_rupture_tension_mn_m=22.5,
        initial_ee=0.91,
        stealth_polymer="PEG-DSPE 2000",
        peg_mol_pct=1.5,
        mean_diameter_nm=95.0,
        citation_doi="10.1021/acs.nanolett.0c03387",
        trial_surveillance="NCT05668741",
    ),
    "SORT_airway_polysarcosine": LNPFormulationSpec(
        formulation_id="SORT_airway_polysarcosine",
        label="Aerosol SORT-LNP (DOPE / 1.2 mol% Polysarcosine / High-Rigidity Core)",
        ionizable_lipid="SORT-Amine-6.50",
        apparent_pka=6.50,
        helper_lipid="DOPE",
        phi_hii=0.082,
        membrane_rupture_tension_mn_m=35.0,
        initial_ee=0.96,
        stealth_polymer="Polysarcosine (pSar)",
        peg_mol_pct=1.2,
        mean_diameter_nm=65.0,
        citation_doi="10.1038/s41565-021-00870-x",
        trial_surveillance="NCT05668741",
    ),
    "high_pka_mucus_trapped": LNPFormulationSpec(
        formulation_id="high_pka_mucus_trapped",
        label="High pKa Unoptimized Formulation (pKa 7.40 / Cationic Airway Surface)",
        ionizable_lipid="Cationic-Amine-7.40",
        apparent_pka=7.40,
        helper_lipid="DSPC",
        phi_hii=0.025,
        membrane_rupture_tension_mn_m=20.0,
        initial_ee=0.88,
        stealth_polymer="PEG-DMG 2000",
        peg_mol_pct=1.0,
        mean_diameter_nm=110.0,
        citation_doi="10.1016/j.biomaterials.2008.12.076",
        trial_surveillance="NCT05668741",
    ),
    "low_pka_escape_deficient": LNPFormulationSpec(
        formulation_id="low_pka_escape_deficient",
        label="Low pKa Escape-Deficient Formulation (pKa 5.20 / Fails Endosomal Protonation)",
        ionizable_lipid="Weak-Amine-5.20",
        apparent_pka=5.20,
        helper_lipid="DOPE",
        phi_hii=0.075,
        membrane_rupture_tension_mn_m=28.0,
        initial_ee=0.93,
        stealth_polymer="PEG-DMG 2000",
        peg_mol_pct=1.5,
        mean_diameter_nm=75.0,
        citation_doi="10.1038/s41565-021-01040-4",
        trial_surveillance="NCT05668741",
    ),
}


def compute_protonation_fraction(ph: float, pka: float) -> float:
    """Compute Henderson-Hasselbalch ionizable lipid protonation fraction theta_prot.

    theta_prot(pH) = 1 / (1 + 10^(pH - pKa))

    Args:
        ph: Ambient environmental pH.
        pka: Apparent acid dissociation constant of the ionizable lipid headgroup.

    Returns:
        Fraction of protonated tertiary amine headgroups in [0, 1].
    """
    delta_ph = ph - pka
    if delta_ph > 15.0:
        return 0.0
    if delta_ph < -15.0:
        return 1.0
    return 1.0 / (1.0 + math.pow(10.0, delta_ph))


def compute_surface_charge_zeta(
    ph: float,
    pka: float,
    zeta_neutral: float = -2.5,
    zeta_cationic: float = 24.0,
) -> float:
    """Compute apparent nanoparticle surface zeta potential (mV) as a function of pH.

    At extracellular airway surface liquid (pH 7.4), optimal pKa (6.4 - 6.6) gives
    near-neutral surface charge (zeta in [-5, +5] mV) evading polyanionic mucin sulfates.
    Upon endosomal acidification (pH <= 5.5), protonation drives positive zeta (> +18 mV).

    Args:
        ph: Environmental pH.
        pka: Apparent pKa of ionizable lipid.
        zeta_neutral: Baseline neutral/slightly anionic zeta potential at pH >> pKa (mV).
        zeta_cationic: Saturated cationic zeta potential at pH << pKa (mV).

    Returns:
        Estimated zeta potential in mV.
    """
    theta = compute_protonation_fraction(ph, pka)
    return zeta_neutral + (zeta_cationic - zeta_neutral) * theta


def compute_mucin_penetration_index(
    zeta_mv: float,
    hydrodynamic_diameter_nm: float = 75.0,
    mucin_mesh_pore_nm: float = 110.0,
) -> float:
    """Compute CF sputum mucin barrier penetration efficiency index in [0, 1].

    Accounts for:
      1. Steric obstruction via Ogston mesh exclusion: P_steric = exp(-(d_h / pore)^2)
      2. Electrostatic immobilization via polyanionic mucin adhesion (MUC5AC/MUC5B):
         Cationic zeta (> +8 mV) suffers severe electrostatic crosslinking.

    Args:
        zeta_mv: Nanoparticle surface zeta potential in mV.
        hydrodynamic_diameter_nm: Nanoparticle hydrodynamic diameter in nm.
        mucin_mesh_pore_nm: Mean macromolecular pore diameter of CF airway mucus in nm.

    Returns:
        Mucin penetration index between 0.0 and 1.0.
    """
    if mucin_mesh_pore_nm <= 0:
        raise ValueError("mucin_mesh_pore_nm must be positive")
    if hydrodynamic_diameter_nm <= 0:
        raise ValueError("hydrodynamic_diameter_nm must be positive")

    # Steric pore penetration
    steric_ratio = hydrodynamic_diameter_nm / mucin_mesh_pore_nm
    p_steric = math.exp(-math.pow(steric_ratio, 2.0))

    # Electrostatic adhesion penalty (logistic transition around +8 mV)
    # Neutral/negative particles (zeta <= 0 mV) have high transmission ~ 1.0
    p_electrostatic = 1.0 / (1.0 + math.exp((zeta_mv - 6.0) / 2.5))

    return max(0.0, min(1.0, p_steric * p_electrostatic))


def compute_escape_rate(
    ph: float,
    pka: float,
    k_max: float = DEFAULT_K_MAX_PER_MIN,
    phi_hii: float = 0.078,
) -> float:
    """Compute pH-dependent endosomal membrane fusion and escape rate constant k_escape.

    k_escape(pH) = k_max * theta_prot(pH) * phi_hii

    Cone-shaped helper lipids (DOPE) adopt inverted hexagonal H_II phase with negative
    spontaneous curvature, dramatically enhancing membrane fusion pore formation over
    cylindrical bilayer-forming lipids (DSPC).

    Args:
        ph: Endosomal compartment pH.
        pka: Ionizable lipid apparent pKa.
        k_max: Maximum intrinsic pore formation escape rate (min^-1).
        phi_hii: Non-bilayer inverted hexagonal H_II phase propensity parameter.

    Returns:
        Escape rate constant k_escape in min^-1.
    """
    if k_max < 0.0:
        raise ValueError("k_max must be non-negative")
    if phi_hii < 0.0:
        raise ValueError("phi_hii must be non-negative")

    theta = compute_protonation_fraction(ph, pka)
    return k_max * theta * phi_hii


def compute_mesh_shear_stress(
    velocity_m_s: float = 15.0,
    pore_diameter_um: float = 3.0,
    dynamic_viscosity_pa_s: float = DEFAULT_FLUID_VISCOSITY_PA_S,
) -> tuple[float, float]:
    """Compute vibrating mesh nebulizer shear rate and wall shear stress.

    gamma_dot_mesh = v_jet / d_pore
    tau_shear = mu_fluid * gamma_dot_mesh

    Args:
        velocity_m_s: Aerosol jet velocity through micro-orifice (m/s, typically 10 - 20 m/s).
        pore_diameter_um: Mesh orifice pore diameter in micrometers (typically 2.5 - 4.0 um).
        dynamic_viscosity_pa_s: Formulation dynamic viscosity in Pa*s (nominally 1.2 mPa*s).

    Returns:
        Tuple of (shear_rate_s_minus_1, wall_shear_stress_pa).
    """
    if pore_diameter_um <= 0.0:
        raise ValueError("pore_diameter_um must be positive")
    if velocity_m_s < 0.0:
        raise ValueError("velocity_m_s must be non-negative")
    if dynamic_viscosity_pa_s <= 0.0:
        raise ValueError("dynamic_viscosity_pa_s must be positive")

    pore_m = pore_diameter_um * 1e-6
    shear_rate = velocity_m_s / pore_m
    shear_stress = dynamic_viscosity_pa_s * shear_rate
    return shear_rate, shear_stress


def compute_nebulization_retention(
    shear_stress_pa: float,
    sigma_membrane_mn_m: float = DEFAULT_SIGMA_MEMBRANE_MN_M,
    ee_0: float = DEFAULT_INITIAL_EE,
    beta_shear_m: float = DEFAULT_BETA_SHEAR_M,
) -> float:
    """Compute post-nebulization mRNA encapsulation efficiency EE_post.

    EE_post(tau_shear) = EE_0 * exp(-beta_shear * (tau_shear / sigma_membrane))

    Where sigma_membrane is converted to N/m (1 mN/m = 1e-3 N/m).
    tau_shear / sigma_membrane has units of m^-1.
    Multiplying by length scale beta_shear_m yields a dimensionless exponent.

    Args:
        shear_stress_pa: Aerosolization wall shear stress in Pascals.
        sigma_membrane_mn_m: Critical bilayer rupture tension in mN/m.
        ee_0: Initial mRNA encapsulation efficiency prior to aerosolization.
        beta_shear_m: Dimensionless shear sensitivity characteristic length scale (meters).

    Returns:
        Post-aerosol encapsulation efficiency in [0, 1].
    """
    if sigma_membrane_mn_m <= 0.0:
        raise ValueError("sigma_membrane_mn_m must be positive")
    if ee_0 < 0.0 or ee_0 > 1.0:
        raise ValueError("ee_0 must be in [0, 1]")

    sigma_n_m = sigma_membrane_mn_m * 1e-3
    stress_ratio = shear_stress_pa / sigma_n_m
    exponent = beta_shear_m * stress_ratio
    return max(0.0, min(1.0, ee_0 * math.exp(-exponent)))


def simulate_endosomal_trafficking(
    total_time_min: float = 360.0,
    dt_min: float = 0.5,
    lnp_dose_internalized: float = 100.0,
    apparent_pka: float = DEFAULT_IONIZABLE_PKA,
    phi_hii: float = 0.078,
    k_max: float = DEFAULT_K_MAX_PER_MIN,
    k_mat: float = DEFAULT_K_MAT_PER_MIN,
    k_recycle: float = DEFAULT_K_RECYCLE_PER_MIN,
    k_lyso: float = DEFAULT_K_LYSO_PER_MIN,
    eta_intact: float = DEFAULT_ETA_INTACT,
    n_cargo: float = DEFAULT_N_CARGO,
    t_half_mrna_hr: float = DEFAULT_T_HALF_MRNA_HR,
    ph_early: float = 6.8,
    ph_late: float = 5.5,
) -> dict[str, Any]:
    """Simulate endosomal trafficking ODE system forward in time.

    ODEs:
      d[LNP_early]/dt = - (k_mat + k_recycle) * [LNP_early]
      d[LNP_late]/dt  = k_mat * [LNP_early] - (k_escape(pH_late) + k_lyso) * [LNP_late]
      d[LNP_escaped]/dt = k_escape(pH_late) * [LNP_late]
      d[LNP_lyso]/dt    = k_lyso * [LNP_late]
      d[LNP_recycle]/dt = k_recycle * [LNP_early]
      d[mRNA_cyto]/dt   = eta_intact * N_cargo * k_escape * [LNP_late] - k_deg_mrna * [mRNA_cyto]

    Args:
        total_time_min: Simulation duration in minutes (nominally 360 min = 6 hours).
        dt_min: Integration time step in minutes.
        lnp_dose_internalized: Initial internalized LNP count at t=0.
        apparent_pka: Ionizable lipid apparent pKa.
        phi_hii: Inverted hexagonal phase parameter.
        k_max: Intrinsic maximum escape rate (min^-1).
        k_mat: Early to late endosome maturation rate (min^-1).
        k_recycle: Exocytosis/recycling rate (min^-1).
        k_lyso: Lysosomal delivery/degradation rate (min^-1).
        eta_intact: Fraction of mRNA cargo escaping intact without degradation.
        n_cargo: mRNA molecules per LNP.
        t_half_mrna_hr: Cytosolic mRNA half-life in hours.
        ph_early: Early endosome lumen pH.
        ph_late: Late acidifying endosome lumen pH.

    Returns:
        Dict of kinetic summary metrics and trajectory endpoints.
    """
    if lnp_dose_internalized <= 0:
        raise ValueError("lnp_dose_internalized must be positive")
    if dt_min <= 0:
        raise ValueError("dt_min must be positive")

    k_deg_mrna_per_min = math.log(2.0) / (t_half_mrna_hr * 60.0)
    k_escape_late = compute_escape_rate(ph_late, apparent_pka, k_max, phi_hii)

    # Initial state
    lnp_early = float(lnp_dose_internalized)
    lnp_late = 0.0
    lnp_escaped_cum = 0.0
    lnp_lyso_cum = 0.0
    lnp_recycled_cum = 0.0
    mrna_cyto = 0.0

    peak_mrna_cyto = 0.0
    peak_time_min = 0.0

    steps = int(total_time_min / dt_min)
    for step in range(steps):
        t = step * dt_min

        if mrna_cyto > peak_mrna_cyto:
            peak_mrna_cyto = mrna_cyto
            peak_time_min = t

        # Fluxes
        flux_mat = k_mat * lnp_early
        flux_recycle = k_recycle * lnp_early
        flux_escape = k_escape_late * lnp_late
        flux_lyso = k_lyso * lnp_late
        flux_mrna_in = eta_intact * n_cargo * flux_escape
        flux_mrna_deg = k_deg_mrna_per_min * mrna_cyto

        # Euler step
        lnp_early -= (flux_mat + flux_recycle) * dt_min
        lnp_early = max(0.0, lnp_early)

        lnp_late += (flux_mat - flux_escape - flux_lyso) * dt_min
        lnp_late = max(0.0, lnp_late)

        lnp_escaped_cum += flux_escape * dt_min
        lnp_lyso_cum += flux_lyso * dt_min
        lnp_recycled_cum += flux_recycle * dt_min

        mrna_cyto += (flux_mrna_in - flux_mrna_deg) * dt_min
        mrna_cyto = max(0.0, mrna_cyto)

    # Escape efficiency percentage relative to internalized dose
    escape_efficiency_pct = (lnp_escaped_cum / lnp_dose_internalized) * 100.0
    lysosomal_degradation_pct = (lnp_lyso_cum / lnp_dose_internalized) * 100.0
    recycling_pct = (lnp_recycled_cum / lnp_dose_internalized) * 100.0

    # Theoretical infinite-time steady state escape fraction
    f_mat = k_mat / (k_mat + k_recycle)
    f_escape_late = k_escape_late / (k_escape_late + k_lyso)
    theoretical_escape_pct = f_mat * f_escape_late * 100.0

    return {
        "lnp_dose_internalized": lnp_dose_internalized,
        "k_escape_late_per_min": round(k_escape_late, 6),
        "escape_efficiency_pct": round(escape_efficiency_pct, 3),
        "theoretical_steady_state_escape_pct": round(theoretical_escape_pct, 3),
        "lysosomal_degradation_pct": round(lysosomal_degradation_pct, 3),
        "recycling_pct": round(recycling_pct, 3),
        "final_cytosolic_mrna_copies": round(mrna_cyto, 2),
        "peak_cytosolic_mrna_copies": round(peak_mrna_cyto, 2),
        "peak_mrna_time_min": round(peak_time_min, 1),
    }


class MRNALNPEscapeModel:
    """Integrated biophysical model of inhaled mRNA-LNP delivery and escape kinetics."""

    def __init__(self, catalog: dict[str, LNPFormulationSpec] | None = None) -> None:
        self.catalog = catalog or LNP_CATALOG

    def evaluate_formulation(
        self,
        formulation_id: str,
        nebulizer_velocity_m_s: float = 15.0,
        mesh_pore_um: float = 3.0,
        airway_ph: float = 7.4,
        late_endosome_ph: float = 5.5,
        lnp_dose: float = 1000.0,
    ) -> dict[str, Any]:
        """Evaluate complete delivery pipeline for a catalog formulation.

        Pipeline Stages:
          1. Vibrating mesh nebulizer shear stress and post-aerosol retention EE_post
          2. Mucus penetration index P_mucus at airway surface pH (7.4)
          3. Airway epithelial cellular internalization (J_in)
          4. Late endosomal protonation inversion and H_II phase escape kinetics
          5. Cytosolic intact mRNA copies delivered per 1,000 aerosolized LNPs
        """
        if formulation_id not in self.catalog:
            raise KeyError(f"Formulation '{formulation_id}' not found in catalog")

        spec = self.catalog[formulation_id]

        # Stage 1: Aerosol shear stress & retention
        shear_rate, shear_stress = compute_mesh_shear_stress(
            velocity_m_s=nebulizer_velocity_m_s,
            pore_diameter_um=mesh_pore_um,
        )
        ee_post = compute_nebulization_retention(
            shear_stress_pa=shear_stress,
            sigma_membrane_mn_m=spec.membrane_rupture_tension_mn_m,
            ee_0=spec.initial_ee,
        )

        # Stage 2: Airway surface zeta potential and mucus penetration
        zeta_airway = compute_surface_charge_zeta(ph=airway_ph, pka=spec.apparent_pka)
        theta_airway = compute_protonation_fraction(ph=airway_ph, pka=spec.apparent_pka)
        p_mucus = compute_mucin_penetration_index(
            zeta_mv=zeta_airway,
            hydrodynamic_diameter_nm=spec.mean_diameter_nm,
        )

        # Stage 3 & 4: Cellular uptake and endosomal trafficking
        # Fraction reaching epithelial cell surface intact
        fraction_epithelial_surface = ee_post * p_mucus
        dose_internalized = lnp_dose * fraction_epithelial_surface

        # Endosomal protonation in late endosome
        theta_late_endo = compute_protonation_fraction(ph=late_endosome_ph, pka=spec.apparent_pka)
        zeta_late_endo = compute_surface_charge_zeta(ph=late_endosome_ph, pka=spec.apparent_pka)

        trafficking = simulate_endosomal_trafficking(
            lnp_dose_internalized=max(1.0, dose_internalized),
            apparent_pka=spec.apparent_pka,
            phi_hii=spec.phi_hii,
            ph_late=late_endosome_ph,
        )

        # Stage 5: Net Delivery Yield
        # Total intact mRNA copies delivered per lnp_dose aerosolized
        net_escape_fraction = (trafficking["escape_efficiency_pct"] / 100.0)
        net_mrna_yield_per_lnp = (
            ee_post * p_mucus * net_escape_fraction * DEFAULT_ETA_INTACT * DEFAULT_N_CARGO
        )

        # Categorization Tier
        if net_mrna_yield_per_lnp >= 0.12:
            delivery_tier = "HIGH_POTENCY_CURATIVE_AIRWAY_DELIVERY"
        elif net_mrna_yield_per_lnp >= 0.04:
            delivery_tier = "MODERATE_EFFICACY_TRANSLATIONAL_REGIME"
        elif net_mrna_yield_per_lnp >= 0.01:
            delivery_tier = "LOW_POTENCY_SUBOPTIMAL_BOTTLENECK"
        else:
            delivery_tier = "INEFFECTIVE_MUCOSAL_OR_ENDOSOMAL_FAILURE"

        return {
            "formulation_id": spec.formulation_id,
            "label": spec.label,
            "ionizable_lipid": spec.ionizable_lipid,
            "apparent_pka": spec.apparent_pka,
            "helper_lipid": spec.helper_lipid,
            "shear_stress_pa": round(shear_stress, 1),
            "shear_rate_s_minus_1": round(shear_rate, 1),
            "post_nebulization_ee": round(ee_post, 4),
            "airway_ph": airway_ph,
            "airway_theta_prot": round(theta_airway, 4),
            "airway_surface_zeta_mv": round(zeta_airway, 2),
            "mucin_penetration_index": round(p_mucus, 4),
            "late_endosome_theta_prot": round(theta_late_endo, 4),
            "late_endosome_surface_zeta_mv": round(zeta_late_endo, 2),
            "k_escape_late_per_min": trafficking["k_escape_late_per_min"],
            "endosomal_escape_pct": trafficking["escape_efficiency_pct"],
            "peak_cytosolic_mrna_copies": trafficking["peak_cytosolic_mrna_copies"],
            "net_mrna_yield_per_lnp": round(net_mrna_yield_per_lnp, 5),
            "delivery_tier": delivery_tier,
            "citation_doi": spec.citation_doi,
            "trial_surveillance": spec.trial_surveillance,
            "research_boundary": RESEARCH_BOUNDARY,
        }

    def compare_all_formulations(self) -> list[dict[str, Any]]:
        """Evaluate and rank all formulations in catalog by net mRNA delivery yield."""
        results = [self.evaluate_formulation(fid) for fid in self.catalog]
        return sorted(results, key=lambda x: x["net_mrna_yield_per_lnp"], reverse=True)

    def run_nebulizer_screen(self, formulation_id: str = "DOPE_airway_optimized") -> list[dict[str, Any]]:
        """Screen formulation retention across common nebulizer modalities."""
        modalities = [
            ("Vibrating Mesh (Aeroneb/eFlow)", 15.0, 3.0),
            ("High-Frequency Mesh (128 kHz)", 18.0, 2.5),
            ("Jet Nebulizer (PARI LC Plus)", 25.0, 1.8),
            ("Soft Mist Inhaler (Respimat)", 8.0, 4.5),
        ]
        out = []
        for name, vel, pore in modalities:
            res = self.evaluate_formulation(
                formulation_id,
                nebulizer_velocity_m_s=vel,
                mesh_pore_um=pore,
            )
            out.append({
                "modality": name,
                "velocity_m_s": vel,
                "pore_um": pore,
                "shear_stress_pa": res["shear_stress_pa"],
                "post_nebulization_ee": res["post_nebulization_ee"],
                "net_mrna_yield_per_lnp": res["net_mrna_yield_per_lnp"],
                "delivery_tier": res["delivery_tier"],
            })
        return out


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Inhaled CFTR mRNA Lipid Nanoparticle Endosomal Escape & Shear Stability Model (T-0239)"
    )
    parser.add_argument("--compare", action="store_true", help="Compare all catalog formulations")
    parser.add_argument("--evaluate", type=str, default="DOPE_airway_optimized", help="Formulation ID to evaluate")
    parser.add_argument("--nebulizer-screen", action="store_true", help="Run nebulizer modality comparison screen")
    parser.add_argument("--json", action="store_true", help="Output raw JSON")
    args = parser.parse_args()

    model = MRNALNPEscapeModel()

    if args.compare:
        comparison = model.compare_all_formulations()
        if args.json:
            print(json.dumps(comparison, indent=2))
        else:
            print("\n" + "=" * 105)
            print("INHALED CFTR mRNA-LNP FORMULATION BENCHMARK & ESCAPE COMPARISON (T-0239)")
            print("=" * 105)
            header = f"{'Formulation ID':<26} | {'pKa':<5} | {'Helper':<6} | {'EE Post':<7} | {'Mucus P':<7} | {'Escape %':<8} | {'Yield/LNP':<9} | {'Tier'}"
            print(header)
            print("-" * 105)
            for r in comparison:
                row = (
                    f"{r['formulation_id']:<26} | "
                    f"{r['apparent_pka']:<5.2f} | "
                    f"{r['helper_lipid']:<6} | "
                    f"{r['post_nebulization_ee']:<7.3f} | "
                    f"{r['mucin_penetration_index']:<7.3f} | "
                    f"{r['endosomal_escape_pct']:<8.2f} | "
                    f"{r['net_mrna_yield_per_lnp']:<9.4f} | "
                    f"{r['delivery_tier']}"
                )
                print(row)
            print("=" * 105)
            print(f"NLM Disclaimer: {NLM_ATTRIBUTION}")
            print(f"Boundary: {RESEARCH_BOUNDARY}\n")
    elif args.nebulizer_screen:
        screen = model.run_nebulizer_screen(args.evaluate)
        if args.json:
            print(json.dumps(screen, indent=2))
        else:
            print(f"\nNebulizer Modality Screen for '{args.evaluate}':")
            for item in screen:
                print(
                    f" - {item['modality']:<32}: Tau = {item['shear_stress_pa']:>7.1f} Pa | "
                    f"EE = {item['post_nebulization_ee']:>6.3f} | "
                    f"Yield/LNP = {item['net_mrna_yield_per_lnp']:>7.4f} | {item['delivery_tier']}"
                )
    else:
        res = model.evaluate_formulation(args.evaluate)
        if args.json:
            print(json.dumps(res, indent=2))
        else:
            print(f"\nEvaluation Results for Formulation: {res['label']}")
            print("-" * 70)
            print(f"Ionizable Lipid pKa:          {res['apparent_pka']:.2f}")
            print(f"Helper Lipid:                 {res['helper_lipid']}")
            print(f"Shear Stress (Mesh Nebulizer):{res['shear_stress_pa']} Pa")
            print(f"Post-Aerosol EE Retention:    {res['post_nebulization_ee'] * 100:.2f}%")
            print(f"Airway Surface Zeta (pH 7.4): {res['airway_surface_zeta_mv']:.2f} mV (theta = {res['airway_theta_prot']:.3f})")
            print(f"Mucin Penetration Index:      {res['mucin_penetration_index']:.3f}")
            print(f"Late Endosome Zeta (pH 5.5):  {res['late_endosome_surface_zeta_mv']:.2f} mV (theta = {res['late_endosome_theta_prot']:.3f})")
            print(f"Endosomal Escape Efficiency:  {res['endosomal_escape_pct']:.2f}%")
            print(f"Peak Cytosolic mRNA Copies:   {res['peak_cytosolic_mrna_copies']} molecules")
            print(f"Net mRNA Yield per LNP:       {res['net_mrna_yield_per_lnp']:.4f}")
            print(f"Delivery Efficacy Tier:       {res['delivery_tier']}")
            print(f"Trial Surveillance Anchor:    {res['trial_surveillance']}")
            print(f"Primary Citation DOI:         {res['citation_doi']}")
            print("-" * 70)


if __name__ == "__main__":
    main()
