#!/usr/bin/env python3
"""Inhaled TMEM16A Alternative Chloride Conductance & Intracellular Calcium Coupling Model (T-0392).

Quantitative biophysical model of TMEM16A (ANO1) calcium-activated chloride channel (CaCC)
conductance and intracellular calcium coupling in airway epithelial cells:
  1. Dual Calcium Cooperativity & Voltage-Dependent Gating:
     Homodimeric TMEM16A gating parameterized by Hill open probability (P_o) with Hill
     coefficient n_H = 2.0 (dual Ca2+ binding per subunit in Cryo-EM PDB 6QP6 / 5OYB) and
     exponential voltage-dependent apparent affinity shift (z_eff = 0.40).
  2. Allosteric Potentiator Left-Shifting:
     Models small-molecule potentiators (e.g., ETD002, ETD001) that allosterically shift
     apparent calcium affinity (K_0.5) from >1500 nM down into resting physiological basal
     range (100-200 nM), permitting therapeutic chloride efflux at basal cytosolic [Ca2+]_i
     without triggering smooth muscle contraction or MUC5AC hypersecretion.
  3. Transepithelial Secretion & Airway Surface Liquid (ASL) Hydration:
     Calculates apical chloride current density (I_Cl,TMEM16A) and integrates osmotic fluid
     replenishment over time to restore periciliary liquid layer height (h_ASL >= 7.0 um).
  4. Genotype-Agnostic Rescue:
     Evaluates conductance restoration across Class I nonsense, Class II refractory, and
     Class III gating mutations unresponsive or refractory to direct CFTR modulators.

Primary Evidence Grounding:
  - Clinical trial surveillance:
    - NCT04488705: Safety, Tolerability and Pharmacokinetics of Inhaled ETD002 in Healthy Volunteers and People With Cystic Fibrosis
    - NCT06478706: Study to Evaluate Efficacy, Safety, Tolerability and PK of Inhaled ETD001 in People With Cystic Fibrosis
  - Cryo-EM structures:
    - PDB 6QP6: Mouse TMEM16A in complex with 1PBC inhibitor at 3.3 A (open/inhibited state)
    - PDB 5OYB: TMEM16A calcium-bound homodimeric architecture at 3.2 A
  - Primary literature:
    - Cabrita I et al. (2019) J Clin Invest, DOI: 10.1172/JCI128564
    - Danahay H et al. (2020) Eur Respir J, DOI: 10.1183/13993003.congress-2020.4449

Epistemic status:
  - Answer Tier T2 (Open-access primary literature, Cryo-EM coordinates, ClinicalTrials.gov API v2).
  - Preclinical computational model; research, not medical advice.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

SCRIPT_PATH = Path(__file__).resolve()

NLM_ATTRIBUTION = (
    "This product uses publicly available data from the U.S. National Library of Medicine (NLM), "
    "National Institutes of Health, Department of Health and Human Services; NLM is not responsible "
    "for the product and does not endorse or recommend this or any other product. "
    "ClinicalTrials.gov does not evaluate safety or efficacy; listing on the registry does not denote government endorsement."
)

RESEARCH_BOUNDARY = (
    "Preclinical in-silico TMEM16A calcium conductance and airway hydration biophysical model; "
    "research, not medical advice. No clinical efficacy or therapeutic safety claim is made."
)

# Physical constants
FARADAY_CONSTANT_C_PER_MOL = 96485.3321
GAS_CONSTANT_J_PER_MOL_K = 8.3144626
BODY_TEMPERATURE_K = 310.15  # 37.0 C
THERMAL_VOLTAGE_MV = (GAS_CONSTANT_J_PER_MOL_K * BODY_TEMPERATURE_K / FARADAY_CONSTANT_C_PER_MOL) * 1000.0  # ~26.726 mV

# Nominal biophysical parameters for TMEM16A gating & epithelial transport (T-0392)
DEFAULT_BASAL_CA_NM = 120.0                # Resting basal epithelial cytosolic [Ca2+]_i (80 - 150 nM)
DEFAULT_K05_ZERO_BASAL_NM = 900.0          # Baseline unpotentiated K_0.5 at 0 mV (nM)
DEFAULT_HILL_COEFFICIENT_NH = 2.0          # Hill cooperativity coefficient (dual Ca2+ binding per subunit)
DEFAULT_EFFECTIVE_VALENCE_ZEFF = 0.40      # Effective gating charge valence for voltage-dependent Ca2+ sensitivity
DEFAULT_POTENTIATOR_KD_ACTIVE_NM = 15.0    # Potentiator (e.g. ETD002) binding affinity for active state K_d,active (nM)
DEFAULT_POTENTIATOR_KD_INACTIVE_NM = 300.0 # Potentiator binding affinity for inactive state K_d,inactive (nM)
DEFAULT_ALLOSTERIC_ALPHA = 0.05            # Allosteric left-shift ratio alpha = K_d,act / K_d,inact (~20-fold max shift)
DEFAULT_APICAL_VM_MV = -40.0               # Resting airway apical membrane electrical potential V_m (mV)
DEFAULT_CL_NERNST_MV = -65.0               # Apical chloride equilibrium reversal potential E_Cl (mV)
DEFAULT_UNITARY_GAMMA_PS = 2.8             # TMEM16A single-channel chloride conductance gamma_Cl (pS)
DEFAULT_CHANNEL_DENSITY_PER_CM2 = 1.8e9    # Apical channel density (channels/cm^2)
DEFAULT_WT_CFTR_ICL_UA_CM2 = 45.0          # Wild-type CFTR secretory chloride current benchmark (uA/cm^2)
DEFAULT_CURATIVE_ICL_PCT_WT = 20.0         # Threshold chloride current (% WT CFTR) to restore hydration
DEFAULT_CURATIVE_ASL_THRESHOLD_UM = 7.0    # Periciliary liquid layer height threshold for ciliary beat clearance (um)
DEFAULT_CF_BASAL_ASL_UM = 3.2              # Dehydrated cystic fibrosis baseline ASL height (um)
DEFAULT_MAX_ASL_UM = 9.8                   # Maximum physiological ASL ceiling under saturating secretion (um)
DEFAULT_BRONCHOCONSTRICTION_CA_NM = 350.0  # Threshold cytosolic Ca2+ triggering smooth muscle spasm / MUC5AC release


@dataclass(frozen=True)
class CFTRGenotypeCohort:
    """Genotype class specification for evaluating CFTR-independent TMEM16A bypass."""

    cohort_id: str
    class_label: str
    representative_genotypes: str
    residual_cftr_pct_wt: float
    enac_hyperabsorption_factor: float
    description: str


COHORT_REGISTRY: Dict[str, CFTRGenotypeCohort] = {
    "Class_I_Nonsense_Null": CFTRGenotypeCohort(
        cohort_id="Class_I_Nonsense_Null",
        class_label="Class I (Premature Termination Codon / Splice Null)",
        representative_genotypes="W1282X/W1282X, G542X/G542X, c.1585-1G>A",
        residual_cftr_pct_wt=0.0,
        enac_hyperabsorption_factor=1.0,
        description="Complete absence of CFTR channel synthesis; relies entirely on alternative chloride conductance.",
    ),
    "Class_II_Refractory_Folding": CFTRGenotypeCohort(
        cohort_id="Class_II_Refractory_Folding",
        class_label="Class II (Modulator-Refractory Misfolding)",
        representative_genotypes="N1303K/N1303K, I507del/N1303K",
        residual_cftr_pct_wt=1.5,
        enac_hyperabsorption_factor=0.98,
        description="Severe trafficking defect with negligible response to standard tri-combination modulators.",
    ),
    "Class_III_Gating_Defect": CFTRGenotypeCohort(
        cohort_id="Class_III_Gating_Defect",
        class_label="Class III (Gating / Modulator-Intolerant)",
        representative_genotypes="G551D/F508del (intolerant to potentiators)",
        residual_cftr_pct_wt=3.0,
        enac_hyperabsorption_factor=0.95,
        description="Impaired ATP-dependent channel opening requiring orthogonal non-CFTR secretory pathway.",
    ),
}


@dataclass(frozen=True)
class TMEM16ACalciumSimulationResult:
    """Simulation output metrics for TMEM16A calcium conductance and ASL hydration."""

    cohort_id: str
    regimen_label: str
    potentiator_conc_nm: float
    cytosolic_ca_nm: float
    membrane_voltage_mv: float
    k05_apparent_nm: float
    fold_affinity_sensitization: float
    open_probability: float
    tmem16a_current_ua_cm2: float
    total_apical_current_ua_cm2: float
    pct_wt_cftr_current: float
    steady_state_asl_um: float
    ciliary_clearance_restored: bool
    curative_conductance_achieved: bool
    bronchoconstriction_risk: bool
    cryoem_pdb_structures: Tuple[str, str]
    registry_trials: Tuple[str, str]
    primary_dois: Tuple[str, str]


def compute_voltage_shifted_k05_nm(
    vm_mv: float = DEFAULT_APICAL_VM_MV,
    k05_zero_nm: float = DEFAULT_K05_ZERO_BASAL_NM,
    z_eff: float = DEFAULT_EFFECTIVE_VALENCE_ZEFF,
) -> float:
    """Calculate voltage-dependent apparent calcium sensitivity K_0.5(V_m) at zero potentiator.

    K_0.5(V_m) = K_0.5(0) * exp(-z_eff * F * V_m / (R * T))
    Depolarization (less negative V_m) enhances apparent calcium affinity (lowers K_0.5).
    """
    if k05_zero_nm <= 0:
        raise ValueError(f"k05_zero_nm must be positive, got {k05_zero_nm}")
    if z_eff < 0:
        raise ValueError(f"z_eff must be non-negative, got {z_eff}")

    exponent = -(z_eff * vm_mv) / THERMAL_VOLTAGE_MV
    exponent = max(-30.0, min(30.0, exponent))
    return k05_zero_nm * math.exp(exponent)


def compute_potentiated_k05_nm(
    potentiator_conc_nm: float,
    vm_mv: float = DEFAULT_APICAL_VM_MV,
    k05_zero_basal_nm: float = DEFAULT_K05_ZERO_BASAL_NM,
    z_eff: float = DEFAULT_EFFECTIVE_VALENCE_ZEFF,
    kd_active_nm: float = DEFAULT_POTENTIATOR_KD_ACTIVE_NM,
    kd_inactive_nm: float = DEFAULT_POTENTIATOR_KD_INACTIVE_NM,
) -> float:
    """Calculate allosterically shifted apparent calcium affinity K_0.5([Pot], V_m).

    K_0.5([Pot], V_m) = K_0.5(V_m, basal) * (1 + [Pot] / K_d,inactive) / (1 + [Pot] / K_d,active)
    where K_d,active < K_d,inactive gives left-shift toward higher calcium sensitivity (lower K_0.5).
    """
    if potentiator_conc_nm < 0:
        raise ValueError(f"potentiator_conc_nm must be non-negative, got {potentiator_conc_nm}")
    if kd_active_nm <= 0 or kd_inactive_nm <= 0:
        raise ValueError("Dissociation constants must be positive")

    k05_vm_basal = compute_voltage_shifted_k05_nm(
        vm_mv=vm_mv,
        k05_zero_nm=k05_zero_basal_nm,
        z_eff=z_eff,
    )

    ratio_active = potentiator_conc_nm / kd_active_nm
    ratio_inactive = potentiator_conc_nm / kd_inactive_nm

    # Potentiator binds preferentially to active state, shifting K_0.5 lower
    inverted_shift = (1.0 + ratio_inactive) / (1.0 + ratio_active)
    return k05_vm_basal * inverted_shift


def compute_tmem16a_open_probability(
    ca_conc_nm: float = DEFAULT_BASAL_CA_NM,
    potentiator_conc_nm: float = 0.0,
    vm_mv: float = DEFAULT_APICAL_VM_MV,
    k05_zero_basal_nm: float = DEFAULT_K05_ZERO_BASAL_NM,
    n_hill: float = DEFAULT_HILL_COEFFICIENT_NH,
    z_eff: float = DEFAULT_EFFECTIVE_VALENCE_ZEFF,
    kd_active_nm: float = DEFAULT_POTENTIATOR_KD_ACTIVE_NM,
    kd_inactive_nm: float = DEFAULT_POTENTIATOR_KD_INACTIVE_NM,
) -> float:
    """Calculate cooperative calcium- and voltage-dependent open probability P_o.

    P_o([Ca2+]_i, V_m, [Pot]) = 1 / (1 + (K_0.5([Pot], V_m) / [Ca2+]_i)^n_H)
    """
    if ca_conc_nm <= 0:
        raise ValueError(f"ca_conc_nm must be positive, got {ca_conc_nm}")
    if n_hill <= 0:
        raise ValueError(f"n_hill must be positive, got {n_hill}")

    k05_app = compute_potentiated_k05_nm(
        potentiator_conc_nm=potentiator_conc_nm,
        vm_mv=vm_mv,
        k05_zero_basal_nm=k05_zero_basal_nm,
        z_eff=z_eff,
        kd_active_nm=kd_active_nm,
        kd_inactive_nm=kd_inactive_nm,
    )

    ratio = k05_app / ca_conc_nm
    if ratio > 1e4:
        return 0.0
    term = ratio ** n_hill
    return 1.0 / (1.0 + term)


def compute_apical_chloride_flux(
    open_probability: float,
    vm_mv: float = DEFAULT_APICAL_VM_MV,
    e_cl_mv: float = DEFAULT_CL_NERNST_MV,
    gamma_ps: float = DEFAULT_UNITARY_GAMMA_PS,
    channel_density_per_cm2: float = DEFAULT_CHANNEL_DENSITY_PER_CM2,
    residual_cftr_pct_wt: float = 0.0,
) -> Dict[str, float]:
    """Calculate macroscopic apical chloride current density (uA/cm^2) and % WT CFTR equivalence.

    I_Cl,TMEM16A = N_channel * gamma_single * P_o * |V_m - E_Cl|
    """
    if not (0.0 <= open_probability <= 1.0):
        raise ValueError(f"open_probability must be in [0, 1], got {open_probability}")
    if gamma_ps <= 0 or channel_density_per_cm2 <= 0:
        raise ValueError("gamma_ps and channel_density_per_cm2 must be positive")

    driving_force_v = abs(vm_mv - e_cl_mv) * 1e-3  # V
    gamma_siemens = gamma_ps * 1e-12                # Siemens (S)
    i_tmem16a_a_cm2 = channel_density_per_cm2 * open_probability * gamma_siemens * driving_force_v
    i_tmem16a_ua_cm2 = i_tmem16a_a_cm2 * 1e6

    i_cftr_residual_ua = (max(0.0, residual_cftr_pct_wt) / 100.0) * DEFAULT_WT_CFTR_ICL_UA_CM2
    total_icl_ua_cm2 = i_tmem16a_ua_cm2 + i_cftr_residual_ua
    pct_wt_cftr = (total_icl_ua_cm2 / DEFAULT_WT_CFTR_ICL_UA_CM2) * 100.0

    return {
        "tmem16a_current_ua_cm2": round(i_tmem16a_ua_cm2, 4),
        "total_apical_current_ua_cm2": round(total_icl_ua_cm2, 4),
        "pct_wt_cftr_current": round(pct_wt_cftr, 4),
    }


def integrate_asl_hydration_ode(
    total_icl_ua_cm2: float,
    enac_hyperabsorption_factor: float = 1.0,
    t_max_hr: float = 8.0,
    dt_hr: float = 0.1,
    h_init_um: float = DEFAULT_CF_BASAL_ASL_UM,
) -> Dict[str, Any]:
    """Integrate osmotic airway surface liquid (ASL) height recovery over time.

    dh_ASL/dt = k_hydr * (h_target(I_Cl, J_ENaC) - h_ASL)
    Target steady-state ASL height reaches >= 7.0 um when I_Cl >= 20% WT CFTR (9.0 uA/cm^2).
    """
    if total_icl_ua_cm2 < 0:
        raise ValueError(f"total_icl_ua_cm2 must be non-negative, got {total_icl_ua_cm2}")
    if t_max_hr <= 0 or dt_hr <= 0:
        raise ValueError("Time parameters must be positive")

    k_half_icl = 6.4 * max(0.5, enac_hyperabsorption_factor)
    osmotic_gain_um = (DEFAULT_MAX_ASL_UM - DEFAULT_CF_BASAL_ASL_UM) * (
        total_icl_ua_cm2 / (k_half_icl + total_icl_ua_cm2)
    )
    h_target_um = DEFAULT_CF_BASAL_ASL_UM + osmotic_gain_um

    k_hydr_rate_hr = 0.90
    current_h_um = h_init_um
    steps = int(round(t_max_hr / dt_hr))
    time_series: List[Dict[str, float]] = []

    for step in range(steps + 1):
        t_hr = step * dt_hr
        dh_dt = k_hydr_rate_hr * (h_target_um - current_h_um)
        time_series.append(
            {
                "time_hr": round(t_hr, 2),
                "asl_height_um": round(current_h_um, 4),
                "rate_um_per_hr": round(dh_dt, 4),
            }
        )
        current_h_um = max(2.0, current_h_um + dh_dt * dt_hr)

    return {
        "steady_state_asl_um": round(h_target_um, 4),
        "final_asl_um": round(time_series[-1]["asl_height_um"], 4),
        "trajectory": time_series,
    }


class TMEM16ACalciumConductanceModel:
    """Biophysical model for TMEM16A allosteric conductance and intracellular calcium coupling (T-0392)."""

    def __init__(self, cohorts: Optional[Dict[str, CFTRGenotypeCohort]] = None) -> None:
        self.cohorts = cohorts or COHORT_REGISTRY

    def simulate(
        self,
        cohort_id: str = "Class_I_Nonsense_Null",
        potentiator_conc_nm: float = 200.0,
        cytosolic_ca_nm: float = DEFAULT_BASAL_CA_NM,
        membrane_voltage_mv: float = DEFAULT_APICAL_VM_MV,
    ) -> TMEM16ACalciumSimulationResult:
        if cohort_id not in self.cohorts:
            raise KeyError(f"Unknown cohort_id '{cohort_id}'. Valid keys: {list(self.cohorts.keys())}")

        cohort = self.cohorts[cohort_id]

        k05_basal_vm = compute_voltage_shifted_k05_nm(
            vm_mv=membrane_voltage_mv,
            k05_zero_nm=DEFAULT_K05_ZERO_BASAL_NM,
            z_eff=DEFAULT_EFFECTIVE_VALENCE_ZEFF,
        )
        k05_app = compute_potentiated_k05_nm(
            potentiator_conc_nm=potentiator_conc_nm,
            vm_mv=membrane_voltage_mv,
            k05_zero_basal_nm=DEFAULT_K05_ZERO_BASAL_NM,
            z_eff=DEFAULT_EFFECTIVE_VALENCE_ZEFF,
            kd_active_nm=DEFAULT_POTENTIATOR_KD_ACTIVE_NM,
            kd_inactive_nm=DEFAULT_POTENTIATOR_KD_INACTIVE_NM,
        )
        sensitization_fold = k05_basal_vm / k05_app if k05_app > 0 else 1.0

        p_open = compute_tmem16a_open_probability(
            ca_conc_nm=cytosolic_ca_nm,
            potentiator_conc_nm=potentiator_conc_nm,
            vm_mv=membrane_voltage_mv,
            k05_zero_basal_nm=DEFAULT_K05_ZERO_BASAL_NM,
            n_hill=DEFAULT_HILL_COEFFICIENT_NH,
            z_eff=DEFAULT_EFFECTIVE_VALENCE_ZEFF,
            kd_active_nm=DEFAULT_POTENTIATOR_KD_ACTIVE_NM,
            kd_inactive_nm=DEFAULT_POTENTIATOR_KD_INACTIVE_NM,
        )

        flux_dict = compute_apical_chloride_flux(
            open_probability=p_open,
            vm_mv=membrane_voltage_mv,
            residual_cftr_pct_wt=cohort.residual_cftr_pct_wt,
        )

        hydration_dict = integrate_asl_hydration_ode(
            total_icl_ua_cm2=flux_dict["total_apical_current_ua_cm2"],
            enac_hyperabsorption_factor=cohort.enac_hyperabsorption_factor,
        )
        ss_asl_um = float(hydration_dict["steady_state_asl_um"])

        if potentiator_conc_nm == 0.0 and cytosolic_ca_nm <= 150.0:
            regimen = "Resting Basal State (Unpotentiated)"
        elif potentiator_conc_nm == 0.0 and cytosolic_ca_nm > 350.0:
            regimen = f"Non-Selective Calcium Elevation ({cytosolic_ca_nm:.0f} nM)"
        else:
            regimen = f"Allosteric Potentiator ETD002-like ({potentiator_conc_nm:.0f} nM, {cytosolic_ca_nm:.0f} nM Ca2+)"

        return TMEM16ACalciumSimulationResult(
            cohort_id=cohort_id,
            regimen_label=regimen,
            potentiator_conc_nm=potentiator_conc_nm,
            cytosolic_ca_nm=cytosolic_ca_nm,
            membrane_voltage_mv=membrane_voltage_mv,
            k05_apparent_nm=round(k05_app, 2),
            fold_affinity_sensitization=round(sensitization_fold, 2),
            open_probability=round(p_open, 5),
            tmem16a_current_ua_cm2=flux_dict["tmem16a_current_ua_cm2"],
            total_apical_current_ua_cm2=flux_dict["total_apical_current_ua_cm2"],
            pct_wt_cftr_current=flux_dict["pct_wt_cftr_current"],
            steady_state_asl_um=round(ss_asl_um, 2),
            ciliary_clearance_restored=ss_asl_um >= DEFAULT_CURATIVE_ASL_THRESHOLD_UM,
            curative_conductance_achieved=flux_dict["pct_wt_cftr_current"] >= DEFAULT_CURATIVE_ICL_PCT_WT,
            bronchoconstriction_risk=cytosolic_ca_nm >= DEFAULT_BRONCHOCONSTRICTION_CA_NM,
            cryoem_pdb_structures=("6QP6", "5OYB"),
            registry_trials=("NCT04488705", "NCT06478706"),
            primary_dois=("10.1172/JCI128564", "10.1183/13993003.congress-2020.4449"),
        )

    def run_dose_response_sweep(
        self,
        cohort_ids: Optional[List[str]] = None,
        doses_nm: Optional[List[float]] = None,
    ) -> List[TMEM16ACalciumSimulationResult]:
        target_cohorts = cohort_ids or list(self.cohorts.keys())
        target_doses = doses_nm or [0.0, 25.0, 50.0, 100.0, 200.0, 400.0]

        results: List[TMEM16ACalciumSimulationResult] = []
        for cid in target_cohorts:
            for dose in target_doses:
                results.append(self.simulate(cohort_id=cid, potentiator_conc_nm=dose))
        return results


def format_markdown_table(results: List[TMEM16ACalciumSimulationResult]) -> str:
    """Format simulation results into a grounded Markdown report."""
    lines = [
        "# TMEM16A Alternative Chloride Conductance & Calcium Coupling Report (T-0392)",
        "",
        "**Clinical Trial Surveillance**: `NCT04488705` (Inhaled ETD002), `NCT06478706` (Inhaled ETD001)",
        "**Cryo-EM Structural Grounding**: PDB `6QP6` (mouse TMEM16A open/inhibited at 3.3 A), PDB `5OYB` (calcium-bound dimer at 3.2 A)",
        "**Primary Literature DOIs**: `10.1172/JCI128564` (Cabrita 2019), `10.1183/13993003.congress-2020.4449` (Danahay 2020)",
        f"**Curative Target Threshold**: ASL Height >= {DEFAULT_CURATIVE_ASL_THRESHOLD_UM:.1f} um at I_Cl >= {DEFAULT_CURATIVE_ICL_PCT_WT:.1f}% WT CFTR",
        "",
        "| Cohort | [Pot] (nM) | [Ca2+]_i (nM) | K_0.5 App. (nM) | Sensitization | P_o | I_Cl (uA/cm2) | % WT CFTR | ASL (um) | Hydration State |",
        "| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |",
    ]

    for r in results:
        status_badge = "**RESTORED (>=7.0 um)**" if r.ciliary_clearance_restored else "Dehydrated (<7.0 um)"
        lines.append(
            f"| {r.cohort_id} | {r.potentiator_conc_nm:.0f} | {r.cytosolic_ca_nm:.0f} | "
            f"{r.k05_apparent_nm:.1f} | {r.fold_affinity_sensitization:.1f}x | {r.open_probability * 100:.2f}% | "
            f"{r.total_apical_current_ua_cm2:.2f} | {r.pct_wt_cftr_current:.1f}% | **{r.steady_state_asl_um:.2f} um** | {status_badge} |"
        )

    lines.extend(
        [
            "",
            "## Safety & Epistemic Verification Boundaries",
            f"- {NLM_ATTRIBUTION}",
            f"- {RESEARCH_BOUNDARY}",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Inhaled TMEM16A Alternative Chloride Conductance & Calcium Coupling Model (T-0392)"
    )
    parser.add_argument(
        "--cohort",
        choices=list(COHORT_REGISTRY.keys()),
        default="Class_I_Nonsense_Null",
        help="Target CFTR genotype cohort (default: Class_I_Nonsense_Null)",
    )
    parser.add_argument(
        "--potentiator-nm",
        type=float,
        default=200.0,
        help="Allosteric potentiator concentration in nM (default: 200.0 nM)",
    )
    parser.add_argument(
        "--ca-nm",
        type=float,
        default=DEFAULT_BASAL_CA_NM,
        help="Basal epithelial cytosolic calcium [Ca2+]_i in nM (default: 120.0 nM)",
    )
    parser.add_argument(
        "--vm-mv",
        type=float,
        default=DEFAULT_APICAL_VM_MV,
        help="Apical membrane electrical potential V_m in mV (default: -40.0 mV)",
    )
    parser.add_argument(
        "--matrix",
        action="store_true",
        help="Run dose-response matrix across all cohorts",
    )
    parser.add_argument(
        "--report",
        action="store_true",
        help="Generate Markdown summary report",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Output raw JSON format",
    )

    args = parser.parse_args()
    model = TMEM16ACalciumConductanceModel()

    if args.matrix or args.report:
        matrix_res = model.run_dose_response_sweep()
        if args.json:
            print(json.dumps([asdict(item) for item in matrix_res], indent=2))
        else:
            print(format_markdown_table(matrix_res))
    else:
        res = model.simulate(
            cohort_id=args.cohort,
            potentiator_conc_nm=args.potentiator_nm,
            cytosolic_ca_nm=args.ca_nm,
            membrane_voltage_mv=args.vm_mv,
        )
        if args.json:
            print(json.dumps(asdict(res), indent=2))
        else:
            status_text = "RESTORED (>=7.0 um)" if res.ciliary_clearance_restored else "DEPLETED (<7.0 um)"
            print("=" * 80)
            print("TMEM16A ALTERNATIVE CHLORIDE CONDUCTANCE & CALCIUM COUPLING MODEL (T-0392)")
            print("=" * 80)
            print(f"Cohort:                       {res.cohort_id}")
            print(f"Regimen:                      {res.regimen_label}")
            print(f"Potentiator Conc [Pot]:       {res.potentiator_conc_nm:.1f} nM")
            print(f"Basal Cytosolic [Ca2+]_i:     {res.cytosolic_ca_nm:.1f} nM")
            print(f"Apparent Sensitivity K_0.5:   {res.k05_apparent_nm:.1f} nM ({res.fold_affinity_sensitization:.2f}x fold shift)")
            print(f"Open Probability P_o:         {res.open_probability * 100:.2f}% (Hill n_H = {DEFAULT_HILL_COEFFICIENT_NH})")
            print(f"Apical Chloride Current:      {res.total_apical_current_ua_cm2:.2f} uA/cm^2 ({res.pct_wt_cftr_current:.1f}% WT CFTR)")
            print(f"Steady-State ASL Depth:       {res.steady_state_asl_um:.2f} um ({status_text})")
            print(f"Bronchospasm Risk Flag:       {'HIGH RISK' if res.bronchoconstriction_risk else 'NORMAL (Basal Ca2+ Safe)'}")
            print("-" * 80)
            print(f"Cryo-EM Coordinates:          PDB {res.cryoem_pdb_structures[0]}, PDB {res.cryoem_pdb_structures[1]}")
            print(f"Clinical Trial Anchors:       {res.registry_trials[0]} (ETD002), {res.registry_trials[1]} (ETD001)")
            print(f"Primary Evidence DOIs:        {res.primary_dois[0]}, {res.primary_dois[1]}")
            print("Scope: Preclinical research model, not medical advice.")
            print("=" * 80)


if __name__ == "__main__":
    main()
