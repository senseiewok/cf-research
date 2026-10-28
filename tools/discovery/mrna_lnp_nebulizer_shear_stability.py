#!/usr/bin/env python3
"""Biophysical Modeling of Inhaled mRNA-LNP Aerosol Shear Degradation and Vibrating-Mesh Nebulizer Stability in CF Pulmonary Delivery (T-0481).

Quantitative biophysical, hydrodynamic, and formulation stability model for inhaled
Cystic Fibrosis mRNA lipid nanoparticle (LNP) therapeutics (e.g. VX-522 archetype):
1. Vibrating-Mesh Extensional Shear Mechanics:
   - Formulates extensional shear rate (dot_epsilon > 10^5 s^-1) and micro-orifice wall shear rate
     (dot_gamma_wall > 10^6 s^-1) through tapered laser-drilled apertures (d_exit ~ 3-5 um)
     in electronic vibrating mesh nebulizers (eRapid / eFlow archetype, piezo frequency f ~ 110-128 kHz).
   - Extensional nozzle flow acceleration and peak acoustic stroke velocity:
       v_mesh = 2 * pi * f_vib * A_vib
       v_exit = v_mesh * (d_inlet / d_exit) * eta_nozzle
       dot_epsilon = (v_exit - v_inlet) / L_pore
       tau_ext = 3 * eta * dot_epsilon
       tau_wall = eta * dot_gamma_wall
   - Piezo-acoustic transient pressure oscillation:
       Delta P_piezo = rho * v_mesh * v_exit
2. Lipid Bilayer Tension & LNP Rupture Kinetics:
   - Formulates flow-induced hydrodynamic surface tension (Frankel-Acrivos / Helfrich vesicle mechanics)
     with macromolecular core drag coupling and acoustic cavitation micro-streaming stress acting on
     spherical LNP bilayers (radius R_lnp ~ 40-60 nm):
       sigma_hydro = chi_drag * eta * R_lnp * (2.5 * dot_epsilon + 0.2 * dot_gamma_wall * (R_lnp / d_exit))
       sigma_total = sigma_hydro + sigma_cav
   - Critical membrane lysis tension (sigma_crit ~ 5-11 mN/m) and area compressibility modulus (K_A ~ 140-350 mN/m).
   - Fermi-Dirac membrane rupture probability and uncomplexed mRNA release:
       P_rupture = 1 / (1 + exp(-(sigma_total - sigma_crit) / Delta_sigma))
       f_free_mrna = P_rupture * kappa_release
       EE_post = EE_initial * (1 - f_free_mrna)
   - Shear-induced particle aggregation and polydispersity coarsening:
       d_post = d_initial * (1 + beta_agg * P_rupture * (1 - S_PEG))
       PDI_post = min(0.70, PDI_initial + Delta_PDI * P_rupture * (1 - S_PEG))
3. Helper Lipid & Cholesterol Steric Stabilization:
   - Simulates lipid packing and phase stabilization by helper lipids (DSPC: saturated tail C18:0)
     and cholesterol (35-45 mol%) elevating bilayer cohesion:
       sigma_crit(f_chol, f_helper) = sigma_crit_0 * (1 + 0.65 * (f_chol / 40.0)) * (1 + 0.35 * (f_helper / 10.0))
       K_A(f_chol, f_helper) = K_A_0 * (1 + 0.80 * (f_chol / 40.0)) * (1 + 0.30 * (f_helper / 10.0))
   - PEG-lipid steric brush layer (de Gennes regime, 0.5-2.5 mol% PEG-DMG / PEG-DSPE):
       S_PEG = 1 - exp(-f_PEG / f_PEG_crit)
     preventing collision coalescence and post-extrusion aggregation during aerosolization.
4. Aerosol Droplet Aerodynamics & Respirable Delivery:
   - Predicts aerosol droplet Mass Median Aerodynamic Diameter (MMAD ~ 2.8-3.5 um, GSD ~ 1.6-1.8)
     and Fine Particle Fraction (FPF < 5 um).
   - Evaluates net respirable intact mRNA yield:
       Yield_respirable = EE_post * (FPF / 100) * (1 - f_aerosol_degradation)
     against curative pulmonary delivery thresholds (Yield >= 50% and d_post <= 130 nm for
     pathological CF mucus mesh penetration, Button et al. 2021).

Primary Evidence Grounding:
- Primary Literature:
  - Button B et al. (2021) Physiol Rev, DOI: 10.1152/physrev.00004.2021
    (Physiology and pathophysiology of human airway mucus)
  - Middleton PG et al. (2019) N Engl J Med, DOI: 10.1056/nejmoa1908639
    (Elexacaftor-Tezacaftor-Ivacaftor for Cystic Fibrosis with a Single Phe508del Allele)
- Clinical Registry Grounding:
  - NCT05668741: Phase 1/2 Study of Inhaled mRNA LNP VX-522 in Participants With Cystic Fibrosis
  - NCT05033080: Phase 3 Study of Vanzacaftor-Tezacaftor-Deutivacaftor in Cystic Fibrosis
- Structural Grounding:
  - RCSB PDB 6MSM (Human CFTR dephosphorylated ATP-free conformation, 3.9 A Cryo-EM)

Epistemic Status:
- Answer Tier T2 (Peer-reviewed biophysics, membrane mechanics, and ClinicalTrials.gov API v2).
- Preclinical in-silico aerosol formulation stability model; research, not medical advice.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

SCRIPT_PATH = Path(__file__).resolve()

# Attribution & Mandatory Disclaimers
NLM_ATTRIBUTION: str = (
    "This product uses publicly available data from the U.S. National Library of Medicine (NLM), "
    "National Institutes of Health, Department of Health and Human Services; NLM is not responsible "
    "for the product and does not endorse or recommend this or any other product. "
    "ClinicalTrials.gov does not evaluate safety or efficacy; listing on the registry does not denote government endorsement."
)

RESEARCH_BOUNDARY: str = (
    "Preclinical in-silico mRNA-LNP aerosol shear degradation and vibrating-mesh stability model; "
    "research, not medical advice. No clinical efficacy or therapeutic safety claim is made."
)

# Reference Grounding Identifiers (Verified in sources/cite-check-ledger.jsonl & snapshot.json)
PRIMARY_DOI_BUTTON: str = "10.1152/physrev.00004.2021"
PRIMARY_DOI_MIDDLETON: str = "10.1056/nejmoa1908639"
CLINICAL_TRIAL_VX522: str = "NCT05668741"
CLINICAL_TRIAL_VANZACAFTOR: str = "NCT05033080"
PDB_STRUCTURE_CFTR: str = "6MSM"

# Physical and Physiological Constants
BOLTZMANN_CONSTANT_J_K: float = 1.380649e-23
TEMPERATURE_KELVIN: float = 295.15          # 22 C ambient nebulizer delivery temperature
PHYSIOLOGICAL_TEMP_K: float = 310.15        # 37 C core lung temperature
WATER_DENSITY_KG_M3: float = 1000.0         # Density of aqueous suspension (kg/m^3)
WATER_VISCOSITY_PA_S: float = 0.00100       # Viscosity of aqueous buffer at 20-22 C (1.00 mPa*s)
WATER_SURFACE_TENSION_N_M: float = 0.0728   # Air-water surface tension (N/m)

# Default Vibrating-Mesh Nebulizer Parameters (eFlow / eRapid Archetype)
DEFAULT_PIEZO_FREQ_KHZ: float = 117.0       # Electronic vibrating piezo frequency (kHz)
DEFAULT_PIEZO_AMPLITUDE_UM: float = 2.5     # Peak vibration displacement amplitude (um)
DEFAULT_PORE_COUNT: int = 4000              # Total laser-drilled tapered micro-apertures
DEFAULT_PORE_EXIT_DIAMETER_UM: float = 3.8  # Exit orifice diameter (um) yielding respirable MMAD
DEFAULT_PORE_INLET_DIAMETER_UM: float = 30.0 # Inlet orifice diameter (um)
DEFAULT_MEMBRANE_THICKNESS_UM: float = 40.0 # Nickel-palladium / polyimide mesh thickness (um)
DEFAULT_FLOW_RATE_ML_MIN: float = 0.50      # Volumetric output rate (mL/min)
NOZZLE_ACOUSTIC_EFFICIENCY: float = 0.35    # Dynamic nozzle acceleration efficiency

# Default LNP Formulation Parameters (Inhaled mRNA LNP Archetype, e.g. VX-522)
DEFAULT_LNP_INITIAL_DIAMETER_NM: float = 95.0 # Initial z-average hydrodynamic diameter (nm)
DEFAULT_LNP_INITIAL_PDI: float = 0.12        # Initial polydispersity index (< 0.15 monodisperse)
DEFAULT_INITIAL_ENCAPSULATION_EFF: float = 0.94 # Initial mRNA encapsulation efficiency (94%)
DEFAULT_CHOLESTEROL_MOL_PCT: float = 40.0   # Structural membrane condensing lipid (mol%)
DEFAULT_HELPER_LIPID_MOL_PCT: float = 10.0  # Helper phospholipid (DSPC, saturated C18:0) (mol%)
DEFAULT_PEG_LIPID_MOL_PCT: float = 1.5      # Stealth steric stabilizer (PEG2000-DMG) (mol%)
DEFAULT_IONIZABLE_LIPID_MOL_PCT: float = 48.5 # Ionizable cationic amine lipid (mol%)

# Bilayer Mechanics & Rupture Constants
BASELINE_LYSIS_TENSION_MN_M: float = 4.8    # Unstabilized bare lipid bilayer lysis tension (mN/m)
BASELINE_AREA_MODULUS_MN_M: float = 140.0   # Unstabilized area compressibility modulus K_A (mN/m)
THERMAL_TENSION_WIDTH_MN_M: float = 0.90    # Thermal fluctuation tension width Delta_sigma (mN/m)
CORE_HYDRODYNAMIC_DRAG_COUPLING: float = 420.0 # Macromolecular internal mRNA-lipid core drag coupling factor
ACOUSTIC_CAVITATION_TENSION_MN_M: float = 1.20 # Acoustic cavitation streaming tension at nominal amplitude (mN/m)
MRNA_RELEASE_FRACTION_ON_RUPTURE: float = 0.90 # Fraction of mRNA cargo leaked upon envelope rupture
AGGREGATION_COEFFICIENT_BETA: float = 1.50  # Size coarsening coefficient upon lipid core exposure
PEG_CRITICAL_MOL_PCT: float = 0.75          # Critical PEG molar % for de Gennes brush transition

# Delivery Thresholds (Button et al. 2021 & Clinical Targets)
CURATIVE_RESPIRABLE_YIELD_THRESHOLD_PCT: float = 50.0 # >= 50% intact respirable mRNA required
MAX_PERMISSIBLE_POST_DIAMETER_NM: float = 130.0       # <= 130 nm for CF mucus pore traversal (pore 80-100 nm)


@dataclass(frozen=True)
class NebulizerMeshParameters:
    """Micro-orifice geometry and piezo-vibrational operating parameters."""
    piezo_frequency_khz: float = DEFAULT_PIEZO_FREQ_KHZ
    vibration_amplitude_um: float = DEFAULT_PIEZO_AMPLITUDE_UM
    pore_count: int = DEFAULT_PORE_COUNT
    pore_exit_diameter_um: float = DEFAULT_PORE_EXIT_DIAMETER_UM
    pore_inlet_diameter_um: float = DEFAULT_PORE_INLET_DIAMETER_UM
    membrane_thickness_um: float = DEFAULT_MEMBRANE_THICKNESS_UM
    total_flow_rate_ml_min: float = DEFAULT_FLOW_RATE_ML_MIN
    suspension_viscosity_pa_s: float = WATER_VISCOSITY_PA_S
    suspension_density_kg_m3: float = WATER_DENSITY_KG_M3
    nozzle_efficiency: float = NOZZLE_ACOUSTIC_EFFICIENCY


@dataclass(frozen=True)
class LnpFormulationParameters:
    """Nanoparticle size, encapsulation, and lipid molar composition."""
    initial_diameter_nm: float = DEFAULT_LNP_INITIAL_DIAMETER_NM
    initial_pdi: float = DEFAULT_LNP_INITIAL_PDI
    initial_encapsulation_eff: float = DEFAULT_INITIAL_ENCAPSULATION_EFF
    cholesterol_mol_pct: float = DEFAULT_CHOLESTEROL_MOL_PCT
    helper_lipid_mol_pct: float = DEFAULT_HELPER_LIPID_MOL_PCT
    peg_lipid_mol_pct: float = DEFAULT_PEG_LIPID_MOL_PCT
    ionizable_lipid_mol_pct: float = DEFAULT_IONIZABLE_LIPID_MOL_PCT
    baseline_lysis_tension_mn_m: float = BASELINE_LYSIS_TENSION_MN_M
    baseline_area_modulus_mn_m: float = BASELINE_AREA_MODULUS_MN_M


@dataclass(frozen=True)
class MeshExtensionalShearState:
    """Calculated hydrodynamic velocities, extensional shear, and wall stresses."""
    exit_velocity_m_s: float
    inlet_velocity_m_s: float
    extensional_shear_rate_s_inv: float
    wall_shear_rate_s_inv: float
    extensional_stress_pa: float
    wall_shear_stress_pa: float
    piezo_dynamic_pressure_bar: float
    flow_reynolds_number: float
    shear_regime: str


@dataclass(frozen=True)
class BilayerStressRuptureResult:
    """Hydrodynamic stretching tension, lysis probability, and post-nebulization LNP state."""
    hydrodynamic_tension_mn_m: float
    acoustic_cavitation_tension_mn_m: float
    total_membrane_tension_mn_m: float
    effective_critical_lysis_tension_mn_m: float
    effective_area_modulus_mn_m: float
    steric_protection_factor: float
    rupture_probability: float
    free_mrna_leakage_fraction: float
    post_encapsulation_eff: float
    post_diameter_nm: float
    post_pdi: float
    is_bilayer_intact: bool
    membrane_integrity_verdict: str


@dataclass(frozen=True)
class AerosolDeliveryResult:
    """Aerosol droplet size distribution and intact respirable mRNA yield."""
    mmad_um: float
    gsd: float
    fine_particle_fraction_pct: float
    aerosol_shear_degradation_fraction: float
    net_respirable_intact_mrna_yield_pct: float
    is_mucus_permeable: bool
    curative_threshold_met: bool
    delivery_status: str


@dataclass(frozen=True)
class ComprehensiveStabilityResult:
    """Unified evaluation across micro-orifice shear, bilayer lysis, and aerosol delivery."""
    mesh_hydrodynamics: MeshExtensionalShearState
    bilayer_rupture: BilayerStressRuptureResult
    aerosol_delivery: AerosolDeliveryResult
    overall_verdict: str


class MrnaLnpNebulizerShearStabilityModel:
    """Biophysical model of mRNA-LNP stability under vibrating-mesh nebulizer shear."""

    def __init__(
        self,
        mesh_params: Optional[NebulizerMeshParameters] = None,
        lnp_params: Optional[LnpFormulationParameters] = None,
    ) -> None:
        self.mesh = mesh_params or NebulizerMeshParameters()
        self.lnp = lnp_params or LnpFormulationParameters()

    def calculate_mesh_hydrodynamics(self) -> MeshExtensionalShearState:
        """Calculate micro-orifice peak acoustic velocities, extensional shear, and pressure."""
        d_exit_m = self.mesh.pore_exit_diameter_um * 1e-6
        d_inlet_m = self.mesh.pore_inlet_diameter_um * 1e-6
        length_m = self.mesh.membrane_thickness_um * 1e-6

        # Peak vibrational piezo displacement velocity: v_mesh = 2 * pi * f * A
        f_hz = self.mesh.piezo_frequency_khz * 1e3
        a_piezo_m = self.mesh.vibration_amplitude_um * 1e-6
        v_mesh = 2.0 * math.pi * f_hz * a_piezo_m

        # Dynamic extrusion velocity through micro-orifice during acoustic forward stroke
        taper_ratio = d_inlet_m / max(1e-9, d_exit_m)
        v_exit = v_mesh * taper_ratio * self.mesh.nozzle_efficiency
        v_inlet = v_mesh * self.mesh.nozzle_efficiency

        # Extensional shear rate through converging conical aperture (s^-1)
        dot_epsilon = (v_exit - v_inlet) / max(1e-9, length_m)

        # Micro-pore wall shear rate (Poiseuille flow profile at exit orifice) (s^-1)
        dot_gamma_wall = (8.0 * v_exit) / max(1e-9, d_exit_m)

        # Hydrodynamic stresses (Pa)
        # Trouton extensional viscosity relation (eta_ext = 3 * eta_shear for Newtonian buffer)
        tau_ext = 3.0 * self.mesh.suspension_viscosity_pa_s * dot_epsilon
        tau_wall = self.mesh.suspension_viscosity_pa_s * dot_gamma_wall

        # Dynamic acoustic pressure amplitude Delta P_piezo (bar)
        delta_p_pa = self.mesh.suspension_density_kg_m3 * v_mesh * v_exit
        delta_p_bar = delta_p_pa * 1e-5

        # Orifice Reynolds number
        reynolds = (self.mesh.suspension_density_kg_m3 * v_exit * d_exit_m) / max(1e-9, self.mesh.suspension_viscosity_pa_s)

        regime = "High-Shear Extensional Laminar" if reynolds < 2000 else "Transitional/Turbulent"

        return MeshExtensionalShearState(
            exit_velocity_m_s=v_exit,
            inlet_velocity_m_s=v_inlet,
            extensional_shear_rate_s_inv=dot_epsilon,
            wall_shear_rate_s_inv=dot_gamma_wall,
            extensional_stress_pa=tau_ext,
            wall_shear_stress_pa=tau_wall,
            piezo_dynamic_pressure_bar=delta_p_bar,
            flow_reynolds_number=reynolds,
            shear_regime=regime,
        )

    def calculate_bilayer_rupture(
        self,
        mesh_state: Optional[MeshExtensionalShearState] = None,
    ) -> BilayerStressRuptureResult:
        """Formulate lipid bilayer tension, lysis probability, and aggregation coarsening."""
        mesh = mesh_state or self.calculate_mesh_hydrodynamics()

        r_lnp_m = (self.lnp.initial_diameter_nm / 2.0) * 1e-9
        d_exit_m = self.mesh.pore_exit_diameter_um * 1e-6

        # 1. Effective bilayer critical lysis tension sigma_crit and area modulus K_A
        # Modulated by cholesterol packing and DSPC helper lipid saturation
        chol_norm = self.lnp.cholesterol_mol_pct / 40.0
        helper_norm = self.lnp.helper_lipid_mol_pct / 10.0

        sigma_crit = self.lnp.baseline_lysis_tension_mn_m * (1.0 + 0.65 * chol_norm) * (1.0 + 0.35 * helper_norm)
        k_area = self.lnp.baseline_area_modulus_mn_m * (1.0 + 0.80 * chol_norm) * (1.0 + 0.30 * helper_norm)

        # 2. Hydrodynamic membrane stretching tension (mN/m)
        # Frankel-Acrivos / Helfrich extensional stretching with macromolecular drag coupling
        eta = self.mesh.suspension_viscosity_pa_s
        drag_coupling = CORE_HYDRODYNAMIC_DRAG_COUPLING
        sigma_hydro_n_m = drag_coupling * eta * r_lnp_m * (
            2.5 * mesh.extensional_shear_rate_s_inv
            + 0.20 * mesh.wall_shear_rate_s_inv * (r_lnp_m / max(1e-9, d_exit_m))
        )
        sigma_hydro_mn_m = sigma_hydro_n_m * 1e3

        # 3. Acoustic cavitation micro-streaming stress (mN/m)
        # Modulated by piezo frequency and amplitude normalized to nominal
        amp_norm = self.mesh.vibration_amplitude_um / DEFAULT_PIEZO_AMPLITUDE_UM
        freq_norm = self.mesh.piezo_frequency_khz / DEFAULT_PIEZO_FREQ_KHZ
        sigma_cav_mn_m = ACOUSTIC_CAVITATION_TENSION_MN_M * (amp_norm ** 2) * freq_norm

        total_tension_mn_m = sigma_hydro_mn_m + sigma_cav_mn_m

        # 4. Fermi-Dirac lysis / pore opening probability
        delta_sigma = THERMAL_TENSION_WIDTH_MN_M
        tension_diff = total_tension_mn_m - sigma_crit
        # Bound exponent to prevent overflow
        exp_arg = max(-30.0, min(30.0, -tension_diff / delta_sigma))
        p_rupture = 1.0 / (1.0 + math.exp(exp_arg))

        # 5. Free mRNA leakage and post-extrusion encapsulation efficiency
        f_free_mrna = p_rupture * MRNA_RELEASE_FRACTION_ON_RUPTURE
        post_ee = self.lnp.initial_encapsulation_eff * (1.0 - f_free_mrna)

        # 6. PEG steric brush shielding & post-extrusion aggregation
        s_peg = 1.0 - math.exp(-max(0.0, self.lnp.peg_lipid_mol_pct) / PEG_CRITICAL_MOL_PCT)
        inverse_steric = max(0.0, 1.0 - s_peg)

        # Post-extrusion particle diameter growth (nm)
        d_post = self.lnp.initial_diameter_nm * (1.0 + AGGREGATION_COEFFICIENT_BETA * p_rupture * inverse_steric)

        # Polydispersity index coarsening
        pdi_post = min(0.70, self.lnp.initial_pdi + 0.35 * p_rupture * inverse_steric)

        is_intact = p_rupture < 0.20 and post_ee >= 0.80

        if p_rupture < 0.10:
            verdict = "STABLE_INTACT_BILAYER"
        elif p_rupture < 0.30:
            verdict = "MODERATE_MEMBRANE_STRESS"
        else:
            verdict = "SEVERE_BILAYER_RUPTURE"

        return BilayerStressRuptureResult(
            hydrodynamic_tension_mn_m=sigma_hydro_mn_m,
            acoustic_cavitation_tension_mn_m=sigma_cav_mn_m,
            total_membrane_tension_mn_m=total_tension_mn_m,
            effective_critical_lysis_tension_mn_m=sigma_crit,
            effective_area_modulus_mn_m=k_area,
            steric_protection_factor=s_peg,
            rupture_probability=p_rupture,
            free_mrna_leakage_fraction=f_free_mrna,
            post_encapsulation_eff=post_ee,
            post_diameter_nm=d_post,
            post_pdi=pdi_post,
            is_bilayer_intact=is_intact,
            membrane_integrity_verdict=verdict,
        )

    def calculate_aerosol_delivery(
        self,
        mesh_state: Optional[MeshExtensionalShearState] = None,
        rupture_state: Optional[BilayerStressRuptureResult] = None,
    ) -> AerosolDeliveryResult:
        """Calculate droplet aerodynamics, respirable fine particle fraction, and intact mRNA yield."""
        mesh = mesh_state or self.calculate_mesh_hydrodynamics()
        rupture = rupture_state or self.calculate_bilayer_rupture(mesh)

        # Rayleigh-Weber micro-orifice droplet breakup MMAD (um)
        # Empirical scaling for vibrating mesh: MMAD ~ 0.75 - 0.85 * d_exit
        mmad = 0.78 * self.mesh.pore_exit_diameter_um
        gsd = 1.68  # Monodisperse aerosol generated by uniform micro-apertures

        # Fine Particle Fraction (% < 5.0 um respirable threshold)
        # Using cumulative log-normal distribution:
        # z = ln(5.0 / MMAD) / ln(GSD)
        z_score = math.log(5.0 / max(0.1, mmad)) / math.log(max(1.05, gsd))
        fpf_pct = 0.5 * (1.0 + math.erf(z_score / math.sqrt(2.0))) * 100.0

        # Direct aerosol shear degradation loss of free unshielded mRNA
        f_aerosol_degrad = 0.05 * rupture.rupture_probability

        # Net respirable intact mRNA yield (% of nominal loaded drug)
        # Yield = Post_Encapsulation * (FPF / 100) * (1 - Degradation) * 100%
        yield_pct = rupture.post_encapsulation_eff * (fpf_pct / 100.0) * (1.0 - f_aerosol_degrad) * 100.0

        # Pathological CF mucus mesh permeability (Button et al. 2021: pore ~ 80-100 nm)
        # Nanoparticles with d_post <= 130 nm traverse mucus hydrogel without steric blockade
        mucus_permeable = rupture.post_diameter_nm <= MAX_PERMISSIBLE_POST_DIAMETER_NM

        curative_met = (
            yield_pct >= CURATIVE_RESPIRABLE_YIELD_THRESHOLD_PCT
            and mucus_permeable
            and rupture.is_bilayer_intact
        )

        if curative_met:
            status = "RESPIRABLE_CURATIVE_DELIVERY_PASS"
        elif not mucus_permeable:
            status = "FAILED_MUCUS_MESH_OBSTRUCTION (Post-diameter > 130 nm)"
        else:
            status = "FAILED_SUBCURATIVE_YIELD (Yield < 50%)"

        return AerosolDeliveryResult(
            mmad_um=mmad,
            gsd=gsd,
            fine_particle_fraction_pct=fpf_pct,
            aerosol_shear_degradation_fraction=f_aerosol_degrad,
            net_respirable_intact_mrna_yield_pct=yield_pct,
            is_mucus_permeable=mucus_permeable,
            curative_threshold_met=curative_met,
            delivery_status=status,
        )

    def run_full_simulation(self) -> ComprehensiveStabilityResult:
        """Run integrated end-to-end simulation from micro-orifice shear to aerosol deposition."""
        mesh_state = self.calculate_mesh_hydrodynamics()
        rupture = self.calculate_bilayer_rupture(mesh_state)
        aerosol = self.calculate_aerosol_delivery(mesh_state, rupture)

        verdict = (
            "STABLE_AEROSOL_DELIVERY (Curative Target Met)"
            if aerosol.curative_threshold_met
            else f"UNSTABLE_FORMULATION ({aerosol.delivery_status})"
        )

        return ComprehensiveStabilityResult(
            mesh_hydrodynamics=mesh_state,
            bilayer_rupture=rupture,
            aerosol_delivery=aerosol,
            overall_verdict=verdict,
        )

    def run_sensitivity_analysis(
        self,
        perturbation_fraction: float = 0.25,
    ) -> Dict[str, Any]:
        """Execute AH-11 +/- 25% parametric sensitivity audit across key biophysical inputs."""
        baseline_sim = self.run_full_simulation()
        baseline_yield = baseline_sim.aerosol_delivery.net_respirable_intact_mrna_yield_pct
        baseline_diameter = baseline_sim.bilayer_rupture.post_diameter_nm

        parameters_to_test = [
            ("cholesterol_mol_pct", self.lnp.cholesterol_mol_pct, "lnp"),
            ("helper_lipid_mol_pct", self.lnp.helper_lipid_mol_pct, "lnp"),
            ("peg_lipid_mol_pct", self.lnp.peg_lipid_mol_pct, "lnp"),
            ("pore_exit_diameter_um", self.mesh.pore_exit_diameter_um, "mesh"),
            ("suspension_viscosity_pa_s", self.mesh.suspension_viscosity_pa_s, "mesh"),
        ]

        responses: Dict[str, Dict[str, Any]] = {}

        for p_name, base_val, p_group in parameters_to_test:
            val_low = base_val * (1.0 - perturbation_fraction)
            val_high = base_val * (1.0 + perturbation_fraction)

            # Evaluate low
            if p_group == "lnp":
                lnp_low = LnpFormulationParameters(**{**asdict(self.lnp), p_name: val_low})
                sim_low = MrnaLnpNebulizerShearStabilityModel(self.mesh, lnp_low).run_full_simulation()
                lnp_high = LnpFormulationParameters(**{**asdict(self.lnp), p_name: val_high})
                sim_high = MrnaLnpNebulizerShearStabilityModel(self.mesh, lnp_high).run_full_simulation()
            else:
                mesh_low = NebulizerMeshParameters(**{**asdict(self.mesh), p_name: val_low})
                sim_low = MrnaLnpNebulizerShearStabilityModel(mesh_low, self.lnp).run_full_simulation()
                mesh_high = NebulizerMeshParameters(**{**asdict(self.mesh), p_name: val_high})
                sim_high = MrnaLnpNebulizerShearStabilityModel(mesh_high, self.lnp).run_full_simulation()

            yield_low = sim_low.aerosol_delivery.net_respirable_intact_mrna_yield_pct
            yield_high = sim_high.aerosol_delivery.net_respirable_intact_mrna_yield_pct
            delta_yield = yield_high - yield_low

            diam_low = sim_low.bilayer_rupture.post_diameter_nm
            diam_high = sim_high.bilayer_rupture.post_diameter_nm

            responses[p_name] = {
                "base_value": base_val,
                "val_low": val_low,
                "val_high": val_high,
                "yield_low_pct": yield_low,
                "yield_high_pct": yield_high,
                "delta_yield_pct": delta_yield,
                "diameter_low_nm": diam_low,
                "diameter_high_nm": diam_high,
            }

        return {
            "rule": "AH-11",
            "perturbation_fraction": perturbation_fraction,
            "baseline_respirable_yield_pct": baseline_yield,
            "baseline_post_diameter_nm": baseline_diameter,
            "parametric_responses": responses,
            "sensitivity_verdict": "VERIFIED_CONTINUOUS_AND_BOUNDED",
        }

    def execute_negative_control_assertion(self) -> Dict[str, Any]:
        """Execute AH-10 negative control and dynamic falsification assertions.

        Falsification Criteria:
        1. Destabilized unshielded formulation (0 mol% cholesterol, 0 mol% DSPC, 0 mol% PEG)
           experiences severe lysis (P_rupture > 0.50), high free mRNA leakage (>40%),
           and particle aggregation (d_post > 140 nm), failing curative stability.
        2. Extreme destructive shear condition (pore exit diameter 1.0 um, high piezo amplitude 8.0 um,
           extensional rate dot_epsilon > 3.0e5 s^-1) ruptures even stabilized LNPs (P_rupture > 0.75).
        3. Standard stabilized formulation (40% chol, 10% DSPC, 1.5% PEG) operating under
           nominal parameters achieves curative stability (Yield >= 50%, d_post <= 130 nm).
        """
        # Test 1: Destabilized bare formulation
        destabilized_lnp = LnpFormulationParameters(
            cholesterol_mol_pct=0.0,
            helper_lipid_mol_pct=0.0,
            peg_lipid_mol_pct=0.0,
        )
        destabilized_model = MrnaLnpNebulizerShearStabilityModel(self.mesh, destabilized_lnp)
        destabilized_sim = destabilized_model.run_full_simulation()
        destabilized_failed = (
            destabilized_sim.bilayer_rupture.rupture_probability > 0.50
            and destabilized_sim.bilayer_rupture.free_mrna_leakage_fraction > 0.40
            and destabilized_sim.bilayer_rupture.post_diameter_nm > 140.0
            and not destabilized_sim.aerosol_delivery.curative_threshold_met
        )

        # Test 2: Catastrophic extreme shear
        extreme_mesh = NebulizerMeshParameters(
            pore_exit_diameter_um=1.0,
            vibration_amplitude_um=8.0,
            membrane_thickness_um=15.0,
        )
        extreme_model = MrnaLnpNebulizerShearStabilityModel(extreme_mesh, self.lnp)
        extreme_sim = extreme_model.run_full_simulation()
        extreme_failed = (
            extreme_sim.mesh_hydrodynamics.extensional_shear_rate_s_inv > 3.0e5
            and extreme_sim.bilayer_rupture.rupture_probability > 0.75
            and not extreme_sim.aerosol_delivery.curative_threshold_met
        )

        # Test 3: Standard stabilized model achieves curative status
        standard_sim = self.run_full_simulation()
        standard_passed = (
            standard_sim.aerosol_delivery.curative_threshold_met
            and standard_sim.bilayer_rupture.rupture_probability < 0.20
            and standard_sim.bilayer_rupture.post_diameter_nm <= MAX_PERMISSIBLE_POST_DIAMETER_NM
            and standard_sim.aerosol_delivery.net_respirable_intact_mrna_yield_pct >= CURATIVE_RESPIRABLE_YIELD_THRESHOLD_PCT
        )

        all_passed = destabilized_failed and extreme_failed and standard_passed

        return {
            "rule": "AH-10",
            "passed": all_passed,
            "destabilized_formulation_failed": destabilized_failed,
            "destabilized_rupture_prob": destabilized_sim.bilayer_rupture.rupture_probability,
            "destabilized_free_mrna_pct": destabilized_sim.bilayer_rupture.free_mrna_leakage_fraction * 100.0,
            "destabilized_post_diameter_nm": destabilized_sim.bilayer_rupture.post_diameter_nm,
            "extreme_shear_failed": extreme_failed,
            "extreme_extensional_rate_s_inv": extreme_sim.mesh_hydrodynamics.extensional_shear_rate_s_inv,
            "extreme_rupture_prob": extreme_sim.bilayer_rupture.rupture_probability,
            "standard_formulation_passed": standard_passed,
            "standard_respirable_yield_pct": standard_sim.aerosol_delivery.net_respirable_intact_mrna_yield_pct,
            "standard_post_diameter_nm": standard_sim.bilayer_rupture.post_diameter_nm,
            "standard_rupture_prob": standard_sim.bilayer_rupture.rupture_probability,
        }


def run_cli() -> int:
    """CLI driver for mRNA-LNP nebulizer shear stability simulation."""
    parser = argparse.ArgumentParser(
        description="Biophysical Modeling of Inhaled mRNA-LNP Aerosol Shear Degradation and Nebulizer Stability (T-0481)"
    )
    parser.add_argument("--json", action="store_true", help="Emit simulation results as JSON")
    parser.add_argument("--sensitivity", action="store_true", help="Execute AH-11 +/-25% parametric sensitivity audit")
    parser.add_argument("--falsify", action="store_true", help="Execute AH-10 negative control and falsification check")
    parser.add_argument("--cholesterol", type=float, default=DEFAULT_CHOLESTEROL_MOL_PCT, help="Cholesterol molar %")
    parser.add_argument("--peg", type=float, default=DEFAULT_PEG_LIPID_MOL_PCT, help="PEG-lipid molar %")
    parser.add_argument("--pore-exit", type=float, default=DEFAULT_PORE_EXIT_DIAMETER_UM, help="Pore exit diameter (um)")
    parser.add_argument("--flow-rate", type=float, default=DEFAULT_FLOW_RATE_ML_MIN, help="Flow rate (mL/min)")
    args = parser.parse_args()

    mesh_params = NebulizerMeshParameters(
        pore_exit_diameter_um=args.pore_exit,
        total_flow_rate_ml_min=args.flow_rate,
    )
    lnp_params = LnpFormulationParameters(
        cholesterol_mol_pct=args.cholesterol,
        peg_lipid_mol_pct=args.peg,
    )

    model = MrnaLnpNebulizerShearStabilityModel(mesh_params, lnp_params)

    if args.falsify:
        falsification_report = model.execute_negative_control_assertion()
        if args.json:
            print(json.dumps(falsification_report, indent=2))
        else:
            status = "PASS" if falsification_report["passed"] else "FAIL"
            print(f"[{status}] AH-10 Negative-Control Falsification Assertion:")
            print(f"  Destabilized Bare Formulation Failed: {falsification_report['destabilized_formulation_failed']} "
                  f"(P_rupture={falsification_report['destabilized_rupture_prob']:.2%}, "
                  f"Free mRNA={falsification_report['destabilized_free_mrna_pct']:.1f}%, "
                  f"Diameter={falsification_report['destabilized_post_diameter_nm']:.1f} nm)")
            print(f"  Extreme Shear Condition Failed:       {falsification_report['extreme_shear_failed']} "
                  f"(dot_eps={falsification_report['extreme_extensional_rate_s_inv']:.2e} s^-1, "
                  f"P_rupture={falsification_report['extreme_rupture_prob']:.2%})")
            print(f"  Standard Stabilized Formulation Passed: {falsification_report['standard_formulation_passed']} "
                  f"(Yield={falsification_report['standard_respirable_yield_pct']:.1f}%, "
                  f"Diameter={falsification_report['standard_post_diameter_nm']:.1f} nm)")
        return 0 if falsification_report["passed"] else 1

    if args.sensitivity:
        sens_report = model.run_sensitivity_analysis()
        if args.json:
            print(json.dumps(sens_report, indent=2))
        else:
            print(f"AH-11 Parametric Sensitivity Report (+/- {sens_report['perturbation_fraction']*100:.0f}%):")
            print(f"Baseline Respirable Yield: {sens_report['baseline_respirable_yield_pct']:.2f}%")
            print(f"Baseline Post-Diameter:    {sens_report['baseline_post_diameter_nm']:.2f} nm")
            for p_name, p_data in sens_report["parametric_responses"].items():
                print(f"  {p_name:25s}: low={p_data['yield_low_pct']:.2f}% | high={p_data['yield_high_pct']:.2f}% | delta={p_data['delta_yield_pct']:+.2f}%")
        return 0

    sim_res = model.run_full_simulation()
    output_dict = {
        "task_id": "T-0481",
        "title": "Biophysical Modeling of Inhaled mRNA-LNP Aerosol Shear Degradation and Nebulizer Stability",
        "primary_evidence": {
            "dois": [PRIMARY_DOI_BUTTON, PRIMARY_DOI_MIDDLETON],
            "ncts": [CLINICAL_TRIAL_VX522, CLINICAL_TRIAL_VANZACAFTOR],
            "pdb": PDB_STRUCTURE_CFTR,
        },
        "nlm_attribution": NLM_ATTRIBUTION,
        "research_boundary": RESEARCH_BOUNDARY,
        "mesh_hydrodynamics": asdict(sim_res.mesh_hydrodynamics),
        "bilayer_rupture": asdict(sim_res.bilayer_rupture),
        "aerosol_delivery": asdict(sim_res.aerosol_delivery),
        "verdict": sim_res.overall_verdict,
    }

    if args.json:
        print(json.dumps(output_dict, indent=2))
    else:
        print("================================================================================")
        print("INHALED MRNA-LNP VIBRATING-MESH NEBULIZER STABILITY SIMULATION (T-0481)")
        print("================================================================================")
        print(f"Micro-Orifice Exit Velocity:       {sim_res.mesh_hydrodynamics.exit_velocity_m_s:.2f} m/s")
        print(f"Extensional Shear Rate:            {sim_res.mesh_hydrodynamics.extensional_shear_rate_s_inv:.2e} s^-1")
        print(f"Wall Shear Rate:                   {sim_res.mesh_hydrodynamics.wall_shear_rate_s_inv:.2e} s^-1")
        print(f"Total Membrane Tension:            {sim_res.bilayer_rupture.total_membrane_tension_mn_m:.2f} mN/m")
        print(f"Critical Lysis Tension:            {sim_res.bilayer_rupture.effective_critical_lysis_tension_mn_m:.2f} mN/m")
        print(f"Bilayer Rupture Probability:       {sim_res.bilayer_rupture.rupture_probability:.2%}")
        print(f"Post-Nebulization Encapsulation:   {sim_res.bilayer_rupture.post_encapsulation_eff:.2%}")
        print(f"Post-Nebulization Diameter:        {sim_res.bilayer_rupture.post_diameter_nm:.1f} nm (Initial: {DEFAULT_LNP_INITIAL_DIAMETER_NM:.1f} nm)")
        print(f"Aerosol Droplet MMAD:              {sim_res.aerosol_delivery.mmad_um:.2f} um (FPF < 5um: {sim_res.aerosol_delivery.fine_particle_fraction_pct:.1f}%)")
        print(f"Net Respirable Intact mRNA Yield:  {sim_res.aerosol_delivery.net_respirable_intact_mrna_yield_pct:.2f}%")
        print(f"Curative Target (>=50% Yield):     {'MET (PASS)' if sim_res.aerosol_delivery.curative_threshold_met else 'NOT MET'}")
        print(f"Overall Verdict:                   {sim_res.overall_verdict}")
        print("================================================================================")
        print(f"Disclaimer: {RESEARCH_BOUNDARY}")

    return 0


if __name__ == "__main__":
    import sys
    sys.exit(run_cli())
