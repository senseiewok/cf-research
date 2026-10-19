#!/usr/bin/env python3
"""Extracellular Vesicle (EV) Basal Stem Cell Tropism & Trans-Tight-Junction Permeation Model (T-0234).

Models trans-epithelial paracellular EV permeation across intact and transiently
permeabilized airway tight junctions (ZO-1, Claudin-1/4), basolateral surface
receptor binding (alpha6beta4 integrin, CD9/CD81 tetraspanins), and durable
chimeric stem cell correction accumulation in human bronchial basal cells (TP63+, KRT5+).

Primary Evidence Grounding:
  - Ivleva et al. (2026), Front Bioeng Biotechnol, DOI: 10.3389/fbioe.2026.1855907 (PMID: 42529376)
  - Vijaykumar et al. (2026), JCI Insight, DOI: 10.1172/jci.insight.206375 (PMID: 42611556)
  - ClinicalTrials.gov ID: NCT05248230 (4D-710 aerosolized gene delivery surveillance)
  - Discovery Dossier: cf-research/landscape/frontier-discovery-2026-10-18.md

Research boundary: Preclinical computational transport model; research, not medical advice.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

SCRIPT_PATH = Path(__file__).resolve()

NLM_ATTRIBUTION = (
    "This product uses publicly available data from the U.S. National Library of Medicine (NLM), "
    "National Institutes of Health, Department of Health and Human Services; NLM is not responsible "
    "for the product and does not endorse or recommend this or any other product."
)
RESEARCH_BOUNDARY = (
    "Preclinical in-silico EV basal cell transport and chimerism model; research, not medical advice. "
    "No clinical efficacy or safety claim is made."
)

# Biophysical and physiological default constants
P_BASELINE_CM_PER_S = 1.2e-7        # Baseline tight junction EV permeability (~1000 Ohm*cm^2 TEER)
H_APICAL_CM = 0.05                  # Apical ASL / mucus liquid film depth (cm)
H_BASOLATERAL_CM = 0.02             # Basolateral interstitial diffusion layer thickness (cm)
KM_ENTRY_NM = 28.0                  # Basolateral binding/entry half-saturation constant K_D (nM)
VMAX_ENTRY_PER_HR = 0.45            # Maximum receptor-mediated basal cell uptake rate (hr^-1)
ETA_PAYLOAD_EFFICIENCY = 0.72       # Intracellular delivery and correction efficiency per EV uptake event
MU_CLEARANCE_PER_DAY = 0.0045       # Basal stem cell mitotic senescence/replacement rate (day^-1)
CURATIVE_CHIMERISM_THRESHOLD = 0.08  # 8% basal stem cell correction restores functional mucociliary transport (Vijaykumar 2026)


@dataclass(frozen=True)
class EVPseudotypeProfile:
    """Biophysical characteristics of an extracellular vesicle surface formulation."""

    profile_id: str
    label: str
    mean_diameter_nm: float
    surface_ligand_factor: float      # Phi_pseudotype: binding/uptake multiplier
    transcytosis_enhancement: float   # Transcellular vesicular transport boost
    citation_doi: str


PSEUDOTYPE_CATALOG: dict[str, EVPseudotypeProfile] = {
    "unmodified_exosome": EVPseudotypeProfile(
        profile_id="unmodified_exosome",
        label="Native Naive Exosome (CD63+ wild-type)",
        mean_diameter_nm=110.0,
        surface_ligand_factor=1.00,
        transcytosis_enhancement=1.00,
        citation_doi="10.3389/fbioe.2026.1855907",
    ),
    "tetraspanin_cd9_hi": EVPseudotypeProfile(
        profile_id="tetraspanin_cd9_hi",
        label="CD9/CD81-Clustered Fusogenic EV",
        mean_diameter_nm=95.0,
        surface_ligand_factor=2.35,
        transcytosis_enhancement=1.45,
        citation_doi="10.3389/fbioe.2026.1855907",
    ),
    "integrin_a6b4_targeted": EVPseudotypeProfile(
        profile_id="integrin_a6b4_targeted",
        label="alpha6beta4-Integrin Targeted Engineered EV",
        mean_diameter_nm=88.0,
        surface_ligand_factor=4.60,
        transcytosis_enhancement=2.10,
        citation_doi="10.3389/fbioe.2026.1855907",
    ),
    "vsvg_pseudotyped_ev": EVPseudotypeProfile(
        profile_id="vsvg_pseudotyped_ev",
        label="VSV-G Glycoprotein Pseudotyped EV",
        mean_diameter_nm=125.0,
        surface_ligand_factor=6.20,
        transcytosis_enhancement=2.80,
        citation_doi="10.3389/fbioe.2026.1855907",
    ),
}


def compute_tight_junction_permeability(
    t_hours: float,
    beta_loosening: float = 0.0,
    tau_recovery_hours: float = 2.5,
    p_baseline: float = P_BASELINE_CM_PER_S,
) -> float:
    """Compute time-dependent tight junction permeability P_TJ(t) in cm/s.

    P_TJ(t) = P_baseline * (1 + beta * exp(-t / tau_recovery))
    """
    if t_hours < 0.0:
        raise ValueError(f"t_hours must be non-negative, got {t_hours}")
    if beta_loosening < 0.0:
        raise ValueError(f"beta_loosening must be non-negative, got {beta_loosening}")
    if tau_recovery_hours <= 0.0:
        raise ValueError(f"tau_recovery_hours must be positive, got {tau_recovery_hours}")

    decay = math.exp(-t_hours / tau_recovery_hours)
    return p_baseline * (1.0 + beta_loosening * decay)


def compute_basal_uptake_velocity(
    c_basolateral_nm: float,
    phi_pseudotype: float = 1.0,
    vmax: float = VMAX_ENTRY_PER_HR,
    kd_nm: float = KM_ENTRY_NM,
) -> float:
    """Compute saturating receptor-mediated basal stem cell entry velocity v_entry (hr^-1).

    v_entry = V_max * (C_basolateral / (K_D + C_basolateral)) * Phi_pseudotype
    """
    if c_basolateral_nm < 0.0:
        raise ValueError(f"c_basolateral_nm must be non-negative, got {c_basolateral_nm}")
    if phi_pseudotype < 0.0:
        raise ValueError(f"phi_pseudotype must be non-negative, got {phi_pseudotype}")
    if c_basolateral_nm == 0.0:
        return 0.0

    saturating_flux = c_basolateral_nm / (kd_nm + c_basolateral_nm)
    return vmax * saturating_flux * phi_pseudotype


class EVBasalTropismModel:
    """Computational simulator for EV paracellular permeation and basal stem cell chimerism (T-0234)."""

    def __init__(self) -> None:
        self.profiles = PSEUDOTYPE_CATALOG

    def simulate_delivery(
        self,
        profile_id: str,
        initial_apical_dose_nm: float = 100.0,
        beta_loosening: float = 18.0,
        tau_recovery_hours: float = 2.5,
        duration_hours: float = 48.0,
        dt_hours: float = 0.05,
    ) -> dict[str, Any]:
        """Simulate trans-epithelial transport, basal cell uptake, and durable chimerism."""
        if profile_id not in self.profiles:
            raise ValueError(
                f"Unknown profile_id '{profile_id}'. Valid profiles: {sorted(self.profiles)}"
            )
        if initial_apical_dose_nm <= 0.0:
            raise ValueError(f"initial_apical_dose_nm must be positive, got {initial_apical_dose_nm}")
        if duration_hours <= 0.0:
            raise ValueError(f"duration_hours must be positive, got {duration_hours}")
        if dt_hours <= 0.0 or dt_hours > 0.5:
            raise ValueError(f"dt_hours must be in (0, 0.5], got {dt_hours}")

        profile = self.profiles[profile_id]

        # Compartment state variables:
        c_apical = float(initial_apical_dose_nm)
        c_basolateral = 0.0
        f_basal_chimerism = 0.0
        peak_c_basolateral = 0.0

        mu_clearance_per_hr = MU_CLEARANCE_PER_DAY / 24.0
        total_steps = int(round(duration_hours / dt_hours))

        # Conversion: cm/s to cm/hr for kinetic rate matching
        sec_per_hour = 3600.0

        time_snapshots: dict[str, float] = {}
        target_check_hours = {1.0, 2.0, 4.0, 8.0, 12.0, 24.0, 48.0}

        for step in range(1, total_steps + 1):
            t_curr = (step - 1) * dt_hours

            # Effective permeability including paracellular opening and transcellular vesicular flux
            p_tj_cm_s = compute_tight_junction_permeability(
                t_hours=t_curr,
                beta_loosening=beta_loosening,
                tau_recovery_hours=tau_recovery_hours,
            )
            # Transcellular vesicular transport contribution
            p_effective_cm_s = p_tj_cm_s * profile.transcytosis_enhancement
            p_effective_cm_hr = p_effective_cm_s * sec_per_hour

            # Driving concentration gradient
            conc_diff = max(0.0, c_apical - c_basolateral)
            j_para = p_effective_cm_hr * conc_diff  # nmol / (cm^2 * hr)

            v_entry = compute_basal_uptake_velocity(
                c_basolateral_nm=c_basolateral,
                phi_pseudotype=profile.surface_ligand_factor,
            )

            # Differential equations:
            # dC_apical/dt = -J_para / h_apical
            # dC_basolateral/dt = J_para / h_basolateral - v_entry * C_basolateral
            # dF_basal/dt = eta_payload * v_entry * (1 - F_basal) - mu_clearance * F_basal
            dc_apical = -j_para / H_APICAL_CM
            dc_basolateral = (j_para / H_BASOLATERAL_CM) - (v_entry * c_basolateral)
            df_basal = (
                ETA_PAYLOAD_EFFICIENCY * v_entry * (1.0 - f_basal_chimerism)
                - mu_clearance_per_hr * f_basal_chimerism
            )

            c_apical = max(0.0, c_apical + dc_apical * dt_hours)
            c_basolateral = max(0.0, c_basolateral + dc_basolateral * dt_hours)
            f_basal_chimerism = min(1.0, max(0.0, f_basal_chimerism + df_basal * dt_hours))

            if c_basolateral > peak_c_basolateral:
                peak_c_basolateral = c_basolateral

            elapsed = round(step * dt_hours, 4)
            for ch in target_check_hours:
                if abs(elapsed - ch) < (dt_hours / 2.0):
                    time_snapshots[f"chimerism_hr_{int(ch)}"] = round(f_basal_chimerism, 4)

        meets_curative = bool(f_basal_chimerism >= CURATIVE_CHIMERISM_THRESHOLD)

        return {
            "profile_id": profile.profile_id,
            "label": profile.label,
            "mean_diameter_nm": profile.mean_diameter_nm,
            "surface_ligand_factor": profile.surface_ligand_factor,
            "transcytosis_enhancement": profile.transcytosis_enhancement,
            "initial_apical_dose_nm": round(initial_apical_dose_nm, 2),
            "beta_loosening": round(beta_loosening, 2),
            "tau_recovery_hours": round(tau_recovery_hours, 2),
            "peak_basolateral_conc_nm": round(peak_c_basolateral, 3),
            "final_apical_retention_pct": round((c_apical / initial_apical_dose_nm) * 100.0, 2),
            "final_basal_chimerism_pct": round(f_basal_chimerism * 100.0, 2),
            "curative_chimerism_threshold_pct": round(CURATIVE_CHIMERISM_THRESHOLD * 100.0, 2),
            "meets_curative_chimerism_target": meets_curative,
            "time_snapshots": time_snapshots,
            "citation_doi": profile.citation_doi,
            "clinical_surveillance_nct": "NCT05248230",
            "nlm_attribution": NLM_ATTRIBUTION,
            "research_boundary": RESEARCH_BOUNDARY,
        }

    def compare_all_profiles(
        self,
        initial_apical_dose_nm: float = 100.0,
        beta_loosening: float = 18.0,
    ) -> dict[str, Any]:
        """Compare all engineered EV pseudotype formulations under equivalent transient loosening."""
        profiles_res = {
            pid: self.simulate_delivery(
                pid,
                initial_apical_dose_nm=initial_apical_dose_nm,
                beta_loosening=beta_loosening,
            )
            for pid in self.profiles
        }
        naive_chimerism = profiles_res["unmodified_exosome"]["final_basal_chimerism_pct"]
        targeted_chimerism = profiles_res["integrin_a6b4_targeted"]["final_basal_chimerism_pct"]
        fold_increase = targeted_chimerism / max(0.01, naive_chimerism)

        return {
            "initial_apical_dose_nm": initial_apical_dose_nm,
            "beta_loosening": beta_loosening,
            "profiles": profiles_res,
            "integrin_targeted_vs_naive_fold_increase": round(fold_increase, 2),
            "nlm_attribution": NLM_ATTRIBUTION,
            "research_boundary": RESEARCH_BOUNDARY,
        }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Simulate Extracellular Vesicle (EV) Basal Stem Cell Tropism & Trans-Tight-Junction Permeation (T-0234)."
    )
    parser.add_argument(
        "--profile",
        default="integrin_a6b4_targeted",
        choices=sorted(PSEUDOTYPE_CATALOG.keys()),
        help="EV pseudotype profile identifier",
    )
    parser.add_argument(
        "--dose-nm",
        type=float,
        default=100.0,
        help="Initial apical EV dose (nM)",
    )
    parser.add_argument(
        "--beta",
        type=float,
        default=18.0,
        help="Tight junction transient loosening factor beta (0 = intact, 18 = acoustic/osmotic opening)",
    )
    parser.add_argument(
        "--compare-all",
        action="store_true",
        help="Compare all EV pseudotype profiles across permeation and basal chimerism",
    )
    parser.add_argument("--json", action="store_true", help="Emit JSON output")
    args = parser.parse_args(argv)

    model = EVBasalTropismModel()
    result = (
        model.compare_all_profiles(initial_apical_dose_nm=args.dose_nm, beta_loosening=args.beta)
        if args.compare_all
        else model.simulate_delivery(args.profile, initial_apical_dose_nm=args.dose_nm, beta_loosening=args.beta)
    )
    if args.json:
        print(json.dumps(result, indent=2))
    else:
        print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
