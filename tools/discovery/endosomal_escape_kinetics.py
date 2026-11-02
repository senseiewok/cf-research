#!/usr/bin/env python3
"""Biophysical Modeling of Endosomal Escape Kinetics & Proton Sponge Buffering Capacity (T-0557).

Quantitative biophysical simulation coupling vacuolar H+-ATPase proton pumping,
ionizable lipid proton sponge buffering capacity, counterion chloride osmotic swelling,
endosomal membrane rupture probability, and cytosolic mRNA escape yield in inhaled
mRNA-LNP delivery to human airway basal stem cells (e.g. VX-522 CFTR mRNA).

Mathematical & Biophysical Formulation:
1. Endosomal Proton Influx and Buffering:
   d[H+]/dt = (J_vATPase([H+]) - J_leak([H+])) / (1 + beta_eff([H+]))
   where:
     J_vATPase = J_max / (1 + 10^(pK_pump - pH_endo))    [Active proton pumping flux]
     J_leak    = P_leak * max([H+]_endo - [H+]_cyto, 0)   [Passive proton leak flux]
     beta_eff  = (C_lipid * K_a) / (K_a + [H+])^2        [Dimensionless buffering scaling]
     beta_VanSlyke = ln(10) * C_lipid_mM * (K_a * [H+]) / (K_a + [H+])^2 [mM / pH unit]
     with K_a  = 10^(-pKa_lipid)

2. Osmotic Swelling and Membrane Rupture Probability:
   P_rupture(t) = 1 - exp( - ( Delta_Pi(t) / sigma_lyso )^m )
   where:
     Delta_Pi(t) = Delta_C_counterion(t) * R * T         [Osmotic swelling pressure in kPa]
     sigma_lyso  = Critical endosomal membrane lytic tension threshold (kPa)
     m           = Weibull rupture shape exponent (m >= 2.0)
     Counterion Cl- influx via ClC-7 maintains luminal electroneutrality.

3. Cytosolic mRNA Escape Yield:
   Y_escape = integral_0^t_sort P_rupture(t) * (1 - P_degrade(t)) * k_release dt
   where:
     t_sort    = Endosomal sorting / maturation time window (s, typically 1800 - 3600 s)
     P_degrade = Probability of lysosomal acid hydrolase degradation (pH_endo < 4.8)

Primary Literature Grounding:
- DOI: 10.1038/s41565-020-0669-6 (Selective organ targeting SORT lipid nanoparticles for tissue-specific mRNA delivery)
- DOI: 10.1038/s41587-023-01861-1 (Inhaled delivery of lipid nanoparticles to respiratory epithelia)

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
    "Preclinical in-silico biophysical modeling of endosomal escape kinetics and "
    "proton sponge buffering capacity; research, not medical advice. "
    "No clinical efficacy or therapeutic safety claim is made."
)

CLINICAL_TRIAL_VX522 = "NCT05668741"
PRIMARY_DOI_SORT_LNP = "10.1038/s41565-020-0669-6"
PRIMARY_DOI_INHALED_LNP = "10.1038/s41587-023-01861-1"

# Universal physical and chemical constants
GAS_CONSTANT_R = 8.314462618              # J / (mol * K)
BODY_TEMPERATURE_K = 310.15               # 37 C (physiological body temp)
CYTOSOLIC_PH = 7.2                        # Physiological cytoplasmic pH
CYTOSOLIC_H_MOLAR = 10.0 ** (-CYTOSOLIC_PH)
RT_CONVERSION_KPA_PER_MOLAR = GAS_CONSTANT_R * BODY_TEMPERATURE_K  # 1 mol/L * RT = 2578.5 kPa


@dataclass(frozen=True)
class EndosomalEscapeParameters:
    """Biophysical and physicochemical parameters governing endosomal acidification and rupture."""
    formulation_name: str = "SORT_Airway_LNP"
    pka_ionizable_lipid: float = 6.45          # Apparent pKa of ionizable SORT lipid (6.2 - 6.8)
    lipid_concentration_mm: float = 35.0       # Luminal lipid concentration (mM)
    j_vatpase_max_molar_per_s: float = 1.2e-5  # Maximum v-ATPase proton flux (mol / (L * s))
    pk_pump: float = 5.8                       # v-ATPase thermodynamic backpressure inflection pK
    p_leak_per_s: float = 0.02                 # Endosomal proton membrane permeability coefficient (s^-1)
    sigma_lyso_kpa: float = 22.0               # Endosomal critical rupture osmotic tension (kPa)
    weibull_exponent_m: float = 2.8            # Weibull membrane rupture steepness parameter
    t_sort_seconds: float = 3600.0             # Endocytic maturation sorting window (60 min)
    dt_seconds: float = 1.0                    # Numerical integration time step (s)
    degradation_ph_threshold: float = 4.8      # Critical pH activating lysosomal nucleases/cathepsins
    degradation_rate_per_s: float = 0.0015     # First-order degradation rate constant when pH < threshold
    release_rate_scaling: float = 0.08         # Characteristic cytosolic transit fraction


@dataclass
class EndosomalKineticTrajectory:
    """Time-course output of endosomal acidification, osmotic swelling, and escape yield."""
    formulation_name: str
    pka: float
    lipid_concentration_mm: float
    final_ph: float
    peak_osmotic_pressure_kpa: float
    peak_rupture_probability: float
    cytosolic_escape_yield_pct: float
    lysosomal_degradation_pct: float
    rupture_time_seconds: Optional[float]
    buffering_capacity_peak: float
    escape_verdict: str


def compute_buffering_capacity(h_conc: float, pka: float, c_lipid_mm: float) -> float:
    """Compute Van Slyke differential proton sponge buffering capacity beta_VanSlyke([H+]) in mM/pH unit.

    beta_VanSlyke = ln(10) * C_lipid_mM * (K_a * [H+]) / (K_a + [H+])^2
    """
    if c_lipid_mm <= 0.0 or pka <= 0.0:
        return 0.0
    k_a = 10.0 ** (-pka)
    denominator = (k_a + h_conc) ** 2
    if denominator <= 0.0:
        return 0.0
    beta = math.log(10.0) * c_lipid_mm * (k_a * h_conc) / denominator
    return beta


def simulate_endosomal_escape(params: EndosomalEscapeParameters) -> EndosomalKineticTrajectory:
    """Simulate endosomal acidification, osmotic swelling, and mRNA cytosolic escape."""
    c_lipid_molar = params.lipid_concentration_mm * 1e-3
    k_a = 10.0 ** (-params.pka_ionizable_lipid) if params.pka_ionizable_lipid > 0.0 else 1e-14
    t_max = params.t_sort_seconds
    dt = params.dt_seconds
    steps = int(t_max / dt)

    h_endo = CYTOSOLIC_H_MOLAR
    accumulated_counterions = 0.0
    intact_fraction = 1.0
    escaped_cumulative = 0.0
    degraded_cumulative = 0.0
    peak_delta_pi = 0.0
    peak_p_rupture = 0.0
    peak_buffering = 0.0
    rupture_t: Optional[float] = None

    for step in range(steps):
        t = step * dt
        ph_endo = -math.log10(max(h_endo, 1e-14))

        # 1. Proton pump influx vs leak
        j_pump = params.j_vatpase_max_molar_per_s / (1.0 + 10.0 ** (params.pk_pump - ph_endo))
        j_leak = params.p_leak_per_s * max(h_endo - CYTOSOLIC_H_MOLAR, 0.0)
        net_flux = max(j_pump - j_leak, 0.0)

        # 2. Buffering capacity:
        # beta_eff is dimensionless d[B]/d[H+] for proton ODE
        # van_slyke is in mM/pH unit for reporting
        if params.lipid_concentration_mm > 0.0 and params.pka_ionizable_lipid > 0.0:
            beta_eff = (c_lipid_molar * k_a) / ((k_a + h_endo) ** 2)
            van_slyke = compute_buffering_capacity(h_endo, params.pka_ionizable_lipid, params.lipid_concentration_mm)
        else:
            beta_eff = 0.0
            van_slyke = 0.0

        if van_slyke > peak_buffering:
            peak_buffering = van_slyke

        dh_dt = net_flux / (1.0 + beta_eff)
        h_endo = max(h_endo + dh_dt * dt, 1e-14)

        # 3. Counterion (Cl-) accumulation balancing proton pumping
        accumulated_counterions += net_flux * dt

        # 4. Osmotic swelling pressure Delta_Pi in kPa
        delta_pi_kpa = accumulated_counterions * RT_CONVERSION_KPA_PER_MOLAR
        if delta_pi_kpa > peak_delta_pi:
            peak_delta_pi = delta_pi_kpa

        # 5. Weibull membrane rupture probability
        ratio = max(delta_pi_kpa / params.sigma_lyso_kpa, 0.0)
        p_rupture = 1.0 - math.exp(-(ratio ** params.weibull_exponent_m))
        if p_rupture > peak_p_rupture:
            peak_p_rupture = p_rupture

        if rupture_t is None and p_rupture >= 0.50:
            rupture_t = t

        # 6. Lysosomal degradation hazard if endosome acidifies below threshold
        if ph_endo < params.degradation_ph_threshold:
            p_deg_rate = params.degradation_rate_per_s * (params.degradation_ph_threshold - ph_endo)
        else:
            p_deg_rate = 0.0

        # Fractional transfer
        rate_escape = (params.release_rate_scaling / 1200.0) * p_rupture * max(1.0 - p_deg_rate * 100.0, 0.0) * intact_fraction
        rate_degrade = p_deg_rate * intact_fraction

        escaped_step = rate_escape * dt
        degraded_step = rate_degrade * dt

        # Bound transitions
        total_loss = escaped_step + degraded_step
        if total_loss > intact_fraction:
            factor = intact_fraction / max(total_loss, 1e-12)
            escaped_step *= factor
            degraded_step *= factor

        escaped_cumulative += escaped_step
        degraded_cumulative += degraded_step
        intact_fraction = max(intact_fraction - (escaped_step + degraded_step), 0.0)

    escape_pct = escaped_cumulative * 100.0
    degrade_pct = degraded_cumulative * 100.0
    final_ph = -math.log10(max(h_endo, 1e-14))

    if escape_pct >= 8.0:
        verdict = "EFFICIENT_ENDOSOMAL_ESCAPE"
    elif escape_pct >= 2.0:
        verdict = "MODERATE_ESCAPE_SUB_OPTIMAL"
    else:
        verdict = "ENTRAPPED_LYSOSOMAL_DEGRADATION"

    return EndosomalKineticTrajectory(
        formulation_name=params.formulation_name,
        pka=params.pka_ionizable_lipid,
        lipid_concentration_mm=params.lipid_concentration_mm,
        final_ph=round(final_ph, 2),
        peak_osmotic_pressure_kpa=round(peak_delta_pi, 2),
        peak_rupture_probability=round(peak_p_rupture, 4),
        cytosolic_escape_yield_pct=round(escape_pct, 2),
        lysosomal_degradation_pct=round(degrade_pct, 2),
        rupture_time_seconds=round(rupture_t, 1) if rupture_t is not None else None,
        buffering_capacity_peak=round(peak_buffering, 2),
        escape_verdict=verdict,
    )


def evaluate_endosomal_escape_kinetics() -> Dict[str, Any]:
    """Evaluate endosomal escape kinetics across distinct LNP lipid formulations."""
    formulations = {
        "SORT_Airway_LNP": EndosomalEscapeParameters(
            formulation_name="SORT_Airway_LNP",
            pka_ionizable_lipid=6.45,
            lipid_concentration_mm=35.0,
            sigma_lyso_kpa=22.0,
            release_rate_scaling=0.08,
        ),
        "Standard_MC3_LNP": EndosomalEscapeParameters(
            formulation_name="Standard_MC3_LNP",
            pka_ionizable_lipid=6.44,
            lipid_concentration_mm=20.0,
            sigma_lyso_kpa=28.0,
            release_rate_scaling=0.05,
        ),
        "Unbuffered_Cationic_LNP": EndosomalEscapeParameters(
            formulation_name="Unbuffered_Cationic_LNP",
            pka_ionizable_lipid=0.0,
            lipid_concentration_mm=0.0,
            sigma_lyso_kpa=35.0,
            release_rate_scaling=0.08,
        ),
        "Neutral_Liposome": EndosomalEscapeParameters(
            formulation_name="Neutral_Liposome",
            pka_ionizable_lipid=0.0,
            lipid_concentration_mm=0.0,
            sigma_lyso_kpa=45.0,
            release_rate_scaling=0.08,
        ),
    }

    trajectories: Dict[str, Dict[str, Any]] = {}
    for name, params in formulations.items():
        res = simulate_endosomal_escape(params)
        trajectories[name] = asdict(res)

    sort_res = trajectories["SORT_Airway_LNP"]
    cationic_res = trajectories["Unbuffered_Cationic_LNP"]

    # Falsification assertions
    escape_superiority = sort_res["cytosolic_escape_yield_pct"] > (cationic_res["cytosolic_escape_yield_pct"] * 3.0)
    osmotic_swelling_verified = sort_res["peak_osmotic_pressure_kpa"] > 15.0
    proton_sponge_buffering_verified = sort_res["buffering_capacity_peak"] > 10.0

    return {
        "task_id": "T-0557",
        "primary_evidence": {
            "dois": [PRIMARY_DOI_SORT_LNP, PRIMARY_DOI_INHALED_LNP],
            "clinical_trials": [CLINICAL_TRIAL_VX522],
        },
        "regulatory_disclaimer": NLM_REGISTRY_DISCLAIMER,
        "research_boundary": RESEARCH_BOUNDARY,
        "formulations": trajectories,
        "sort_airway_lnp_escape_yield_pct": sort_res["cytosolic_escape_yield_pct"],
        "unbuffered_cationic_escape_yield_pct": cationic_res["cytosolic_escape_yield_pct"],
        "endosomal_escape_verified": escape_superiority and osmotic_swelling_verified and proton_sponge_buffering_verified,
    }


def falsify_endosomal_escape_kinetics() -> Dict[str, Any]:
    """AH-10 dynamic negative-control falsification check."""
    # Negative control: formulation with zero lipid buffering and high membrane resistance
    neg_params = EndosomalEscapeParameters(
        formulation_name="Poisoned_ZeroBuffering_Control",
        pka_ionizable_lipid=0.0,
        lipid_concentration_mm=0.0,
        sigma_lyso_kpa=100.0,
    )
    neg_res = simulate_endosomal_escape(neg_params)
    passed = neg_res.cytosolic_escape_yield_pct < 1.0 and neg_res.escape_verdict == "ENTRAPPED_LYSOSOMAL_DEGRADATION"

    return {
        "rule": "AH-10",
        "task_id": "T-0557",
        "test_name": "falsify_endosomal_escape_kinetics",
        "negative_control_passed": passed,
        "unbuffered_escape_yield_pct": neg_res.cytosolic_escape_yield_pct,
        "verdict": "PASS" if passed else "FAIL",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Endosomal escape kinetics simulation (T-0557)")
    parser.add_argument("--json", action="store_true", help="Emit raw JSON")
    parser.add_argument("--falsify", action="store_true", help="Run AH-10 negative-control falsification")
    args = parser.parse_args()

    if args.falsify:
        res = falsify_endosomal_escape_kinetics()
        print(json.dumps(res, indent=2))
        return 0 if res["negative_control_passed"] else 1

    report = evaluate_endosomal_escape_kinetics()
    if args.json:
        print(json.dumps(report, indent=2))
    else:
        print(f"Endosomal Escape Kinetics (T-0557): {'PASS' if report['endosomal_escape_verified'] else 'FAIL'}")
        print(f"  SORT Airway LNP Yield: {report['sort_airway_lnp_escape_yield_pct']}%")
        print(f"  Unbuffered Cationic Yield: {report['unbuffered_cationic_escape_yield_pct']}%")
    return 0 if report["endosomal_escape_verified"] else 1


if __name__ == "__main__":
    sys.exit(main())
