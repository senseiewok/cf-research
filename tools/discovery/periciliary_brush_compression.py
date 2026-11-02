#!/usr/bin/env python3
"""Biophysical Modeling of Periciliary Brush Osmotic Compression vs Airway Surface Liquid Height (T-0562).

Quantitative biophysical model of Flory-Huggins grafted periciliary brush osmotic pressure,
mucus layer compressive osmotic pressure, and compressive equilibrium governing ciliary
collapse and transport arrest in CFTR-deficient airway microenvironments.

Mathematical & Biophysical Formulation:
1. Flory-Huggins Grafted Periciliary Brush Osmotic Pressure:
   Pi_PCL(h) = Pi_0 * (h_0 / h)^(9/4)
   where:
     Pi_0 = Resting PCL osmotic pressure (~500 Pa, Button et al. Science 2012)
     h_0  = Uncompressed cilia height (~7.0 um)
     h    = Instantaneous periciliary liquid (PCL) layer height (um)

2. Mucus Layer Compressive Osmotic Pressure:
   Pi_mucus(phi_solids) = A * (phi_solids / phi_0)^B
   where:
     phi_solids = Mucus solids fraction (weight %, typically 1.5 - 2.5% healthy, >= 8.0% CF)
     phi_0      = Baseline hydrated solids fraction (~2.0%)
     A          = Characteristic compressive scaling modulus (~200 Pa)
     B          = Osmotic power-law exponent (~2.5)

3. Compressive Equilibrium and Critical Ciliary Motility Collapse:
   Delta_Pi_net = Pi_mucus(phi_solids) - Pi_PCL(h)
   At mechanical equilibrium:
     h_eq = h_0 * (Pi_0 / Pi_mucus(phi_solids))^(4/9)   [capped at h_0 when Pi_mucus <= Pi_0]
   When phi_solids >= 8.0%, mucus osmotic pressure exceeds 6,400 Pa >> 500 Pa, compressing the
   PCL layer below critical cilia height h_crit (~4.0 um). Cilia bend flat against the apical
   epithelial surface, inducing complete mucociliary transport arrest.

Primary Evidence Grounding:
- DOI: 10.1126/science.1223012 (A Periciliary Brush Promotes the Lung Health by Separating the Mucus Layer from Airway Epithelia)
- DOI: 10.1152/physrev.00004.2021 (Physiology and pathophysiology of human airway mucus)
- DOI: 10.1016/j.biomaterials.2008.12.076 (The penetration of fresh undiluted sputum expectorated by cystic fibrosis patients)

Clinical Trials:
- NCT05668741 (A Phase 1/2 Study of VX-522 in Participants With Cystic Fibrosis)

Mandatory Attribution & Disclaimers:
ClinicalTrials.gov data are cited here with attribution to the U.S. National Library of Medicine (NLM):
"ClinicalTrials.gov does not evaluate safety or efficacy; listing on the registry does not denote government endorsement."
Research boundary: Preclinical computational model; research, not medical advice. No clinical efficacy or therapeutic safety claim is made.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
import sys
from typing import Any, Dict, List, Optional, Tuple

NLM_REGISTRY_DISCLAIMER = (
    "ClinicalTrials.gov does not evaluate safety or efficacy; "
    "listing on the registry does not denote government endorsement."
)
RESEARCH_BOUNDARY = (
    "Preclinical in-silico biophysical modeling of periciliary brush osmotic compression "
    "and airway surface liquid height; research, not medical advice. "
    "No clinical efficacy or therapeutic safety claim is made."
)

CLINICAL_TRIAL_VX522 = "NCT05668741"
PRIMARY_DOI_PERICILIARY_BRUSH = "10.1126/science.1223012"
PRIMARY_DOI_MUCUS_PHYSIOLOGY = "10.1152/physrev.00004.2021"
PRIMARY_DOI_SPUTUM_PENETRATION = "10.1016/j.biomaterials.2008.12.076"


@dataclass(frozen=True)
class PericiliaryBrushParameters:
    """Biophysical periciliary brush and mucus compressive rheology parameters."""
    profile_name: str = "Healthy_Hydrated_Epithelium"
    resting_pcl_osmotic_pressure_pa: float = 500.0  # Pi_0: Resting grafted PCL osmotic pressure (Pa)
    uncompressed_cilia_height_um: float = 7.0       # h_0: Uncompressed cilium height (um)
    critical_cilia_height_um: float = 4.0           # h_crit: Threshold height for motility collapse (um)
    ciliary_arrest_height_um: float = 2.0           # h_arrest: Total mechanical immobilization height (um)
    mucus_solids_fraction_pct: float = 2.0          # phi_solids: Mucus weight solids fraction (%)
    baseline_solids_pct: float = 2.0                # phi_0: Baseline healthy solids fraction (%)
    mucus_pressure_coeff_a_pa: float = 200.0        # A: Mucus compressive prefactor (Pa)
    mucus_pressure_exponent_b: float = 2.5          # B: Mucus compressive scaling exponent
    cbf_baseline_hz: float = 12.0                   # f_0: Baseline uncompressed ciliary beat frequency (Hz)


@dataclass
class PericiliaryBrushProfile:
    """Quantitative periciliary compression and ciliary clearance output."""
    profile_name: str
    mucus_solids_pct: float
    mucus_osmotic_pressure_pa: float
    pcl_osmotic_pressure_pa: float
    equilibrium_height_um: float
    compression_ratio: float
    effective_cbf_hz: float
    net_osmotic_gradient_pa: float
    ciliary_motility_status: str


def compute_mucus_osmotic_pressure(
    solids_pct: float,
    baseline_pct: float = 2.0,
    coeff_a: float = 200.0,
    exp_b: float = 2.5,
) -> float:
    """Compute mucus layer compressive osmotic pressure: Pi_mucus = A * (phi / phi_0)^B in Pa."""
    phi = max(solids_pct, 0.1)
    phi_0 = max(baseline_pct, 0.1)
    return coeff_a * ((phi / phi_0) ** exp_b)


def compute_pcl_osmotic_pressure(
    height_um: float,
    resting_pi_0: float = 500.0,
    uncompressed_h_0: float = 7.0,
) -> float:
    """Compute Flory-Huggins grafted periciliary brush osmotic pressure: Pi_PCL = Pi_0 * (h_0 / h)^(9/4) in Pa."""
    h = max(height_um, 0.5)
    h_0 = max(uncompressed_h_0, 1.0)
    return resting_pi_0 * ((h_0 / h) ** 2.25)


def calculate_equilibrium_pcl_height(params: PericiliaryBrushParameters) -> float:
    """Calculate mechanical equilibrium PCL height where Pi_PCL(h_eq) = Pi_mucus(phi).

    Formula:
      h_eq = h_0 * (Pi_0 / Pi_mucus)^(4/9)
    Capped at h_0 when Pi_mucus <= Pi_0 (brush is uncompressed).
    """
    pi_mucus = compute_mucus_osmotic_pressure(
        params.mucus_solids_fraction_pct,
        params.baseline_solids_pct,
        params.mucus_pressure_coeff_a_pa,
        params.mucus_pressure_exponent_b,
    )
    pi_0 = params.resting_pcl_osmotic_pressure_pa
    h_0 = params.uncompressed_cilia_height_um

    if pi_mucus <= pi_0:
        return h_0

    # Compressive regime
    ratio = pi_0 / pi_mucus
    exponent = 4.0 / 9.0  # (9/4)^(-1)
    h_eq = h_0 * (ratio ** exponent)
    return max(h_eq, 0.5)


def simulate_periciliary_brush(params: PericiliaryBrushParameters) -> PericiliaryBrushProfile:
    """Simulate periciliary brush osmotic compression and ciliary motility state."""
    pi_mucus = compute_mucus_osmotic_pressure(
        params.mucus_solids_fraction_pct,
        params.baseline_solids_pct,
        params.mucus_pressure_coeff_a_pa,
        params.mucus_pressure_exponent_b,
    )
    h_eq = calculate_equilibrium_pcl_height(params)
    pi_pcl = compute_pcl_osmotic_pressure(
        h_eq,
        params.resting_pcl_osmotic_pressure_pa,
        params.uncompressed_cilia_height_um,
    )
    net_gradient = pi_mucus - pi_pcl

    compression_ratio = min(max(h_eq / params.uncompressed_cilia_height_um, 0.0), 1.0)

    # Ciliary Beat Frequency (CBF) response to PCL height
    h_crit = params.critical_cilia_height_um
    h_arrest = params.ciliary_arrest_height_um
    h_0 = params.uncompressed_cilia_height_um

    if h_eq >= h_0:
        cbf_effective = params.cbf_baseline_hz
        status = "ACTIVE_UNCOMPRESSED_MOTILITY"
    elif h_eq >= h_crit:
        # Partial compression: linear degradation of beat frequency
        fraction = (h_eq - h_crit) / (h_0 - h_crit)
        cbf_effective = params.cbf_baseline_hz * (0.60 + 0.40 * fraction)
        status = "COMPROMISED_PARTIAL_COMPRESSION"
    elif h_eq > h_arrest:
        # Severe compression below critical threshold
        fraction = (h_eq - h_arrest) / (h_crit - h_arrest)
        cbf_effective = params.cbf_baseline_hz * 0.60 * fraction
        status = "ARRESTED_CILIARY_COLLAPSE"
    else:
        cbf_effective = 0.0
        status = "ARRESTED_CILIARY_COLLAPSE"

    return PericiliaryBrushProfile(
        profile_name=params.profile_name,
        mucus_solids_pct=round(params.mucus_solids_fraction_pct, 1),
        mucus_osmotic_pressure_pa=round(pi_mucus, 2),
        pcl_osmotic_pressure_pa=round(pi_pcl, 2),
        equilibrium_height_um=round(h_eq, 2),
        compression_ratio=round(compression_ratio, 3),
        effective_cbf_hz=round(cbf_effective, 2),
        net_osmotic_gradient_pa=round(net_gradient, 2),
        ciliary_motility_status=status,
    )


def evaluate_periciliary_brush_compression() -> Dict[str, Any]:
    """Evaluate periciliary brush compression across healthy, CF, and rescued phenotypes."""
    profiles = {
        "Healthy_Hydrated": PericiliaryBrushParameters(
            profile_name="Healthy_Hydrated",
            mucus_solids_fraction_pct=2.0,
        ),
        "CF_Untreated_Hyperconcentrated": PericiliaryBrushParameters(
            profile_name="CF_Untreated_Hyperconcentrated",
            mucus_solids_fraction_pct=8.5,
        ),
        "CF_Partially_Rescued": PericiliaryBrushParameters(
            profile_name="CF_Partially_Rescued",
            mucus_solids_fraction_pct=4.5,
        ),
        "Severe_CF_Dehydrated": PericiliaryBrushParameters(
            profile_name="Severe_CF_Dehydrated",
            mucus_solids_fraction_pct=10.0,
        ),
    }

    results: Dict[str, Dict[str, Any]] = {}
    for name, p in profiles.items():
        res = simulate_periciliary_brush(p)
        results[name] = asdict(res)

    healthy = results["Healthy_Hydrated"]
    cf_untreated = results["CF_Untreated_Hyperconcentrated"]
    cf_rescued = results["CF_Partially_Rescued"]

    # Biophysical invariants verification
    healthy_uncompressed = (
        healthy["equilibrium_height_um"] >= 7.0
        and healthy["ciliary_motility_status"] == "ACTIVE_UNCOMPRESSED_MOTILITY"
    )
    cf_arrested = (
        cf_untreated["equilibrium_height_um"] < 4.0
        and cf_untreated["ciliary_motility_status"] == "ARRESTED_CILIARY_COLLAPSE"
    )
    rescued_intermediate = (
        4.0 <= cf_rescued["equilibrium_height_um"] < 7.0
        and cf_rescued["ciliary_motility_status"] == "COMPROMISED_PARTIAL_COMPRESSION"
    )

    verified = healthy_uncompressed and cf_arrested and rescued_intermediate

    return {
        "task_id": "T-0562",
        "primary_evidence": {
            "dois": [
                PRIMARY_DOI_PERICILIARY_BRUSH,
                PRIMARY_DOI_MUCUS_PHYSIOLOGY,
                PRIMARY_DOI_SPUTUM_PENETRATION,
            ],
            "clinical_trials": [CLINICAL_TRIAL_VX522],
        },
        "regulatory_disclaimer": NLM_REGISTRY_DISCLAIMER,
        "research_boundary": RESEARCH_BOUNDARY,
        "profiles": results,
        "healthy_equilibrium_height_um": healthy["equilibrium_height_um"],
        "cf_untreated_equilibrium_height_um": cf_untreated["equilibrium_height_um"],
        "periciliary_brush_verified": verified,
    }


def falsify_periciliary_brush_compression() -> Dict[str, Any]:
    """AH-10 dynamic negative-control falsification check."""
    # Negative control: extreme hyperconcentration (phi_solids = 12.0%) must cause severe collapse below arrest height
    neg_params = PericiliaryBrushParameters(
        profile_name="Poisoned_Negative_Control",
        mucus_solids_fraction_pct=12.0,
    )
    neg_profile = simulate_periciliary_brush(neg_params)
    passed = (
        neg_profile.equilibrium_height_um < 2.0
        and neg_profile.effective_cbf_hz == 0.0
        and neg_profile.ciliary_motility_status == "ARRESTED_CILIARY_COLLAPSE"
    )

    return {
        "rule": "AH-10",
        "task_id": "T-0562",
        "test_name": "falsify_periciliary_brush_compression",
        "negative_control_passed": passed,
        "collapsed_equilibrium_height_um": neg_profile.equilibrium_height_um,
        "effective_cbf_hz": neg_profile.effective_cbf_hz,
        "verdict": "PASS" if passed else "FAIL",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Periciliary brush osmotic compression simulation (T-0562)")
    parser.add_argument("--json", action="store_true", help="Emit raw JSON")
    parser.add_argument("--falsify", action="store_true", help="Run AH-10 negative-control falsification")
    args = parser.parse_args()

    if args.falsify:
        res = falsify_periciliary_brush_compression()
        print(json.dumps(res, indent=2))
        return 0 if res["negative_control_passed"] else 1

    report = evaluate_periciliary_brush_compression()
    if args.json:
        print(json.dumps(report, indent=2))
    else:
        print(f"Periciliary Brush Osmotic Compression (T-0562): {'PASS' if report['periciliary_brush_verified'] else 'FAIL'}")
        print(f"  Healthy Equilibrium Height: {report['healthy_equilibrium_height_um']} um")
        print(f"  CF Untreated Height: {report['cf_untreated_equilibrium_height_um']} um")
    return 0 if report["periciliary_brush_verified"] else 1


if __name__ == "__main__":
    sys.exit(main())
