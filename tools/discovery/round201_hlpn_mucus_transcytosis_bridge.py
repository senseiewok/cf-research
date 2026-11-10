#!/usr/bin/env python3
"""
Round 201: Hybrid Lipid-Polymer Nanoparticle (hLPN) Airway Epithelial Transcytosis & Mucus Penetration Invariant Bridge (T-0892).

Mathematical & Physical Pillar (Column A - Deductive Theorem):
  1. Steric Hydrodynamic Mucus Penetration of hLPNs:
     Hybrid lipid-polymer nanoparticles (hLPNs) encapsulate therapeutics (such as siRNA or peptide nucleic acid,
     PNA) within a biodegradable polymeric core (e.g., poly(lactic-co-glycolic acid), PLGA) enveloped by a lipid/PEG
     corona (e.g., DPPC and DSPE-PEG2000). In cystic fibrosis airway mucus, a concentrated gel network of mucins
     (MUC5AC and MUC5B) presents both a steric mesh filtration barrier (pore diameter xi_m, typically 25-80 nm in
     CF sputum) and an adhesive hydrophobic trap (Conte et al., ACS Appl Mater Interfaces 2022,
     DOI: 10.1021/acsami.1c14975; Comegna et al., Sci Rep 2021, DOI: 10.1038/s41598-021-85549-z).
     Dense grafting of polyethylene glycol (PEG) creates a hydrophilic, neutrally charged brush that shields the
     hydrophobic core:
       Gamma_peg >= 0 (chains/nm^2)
       d_p > 0 (nm, nanoparticle core diameter)
       S_peg = Gamma_peg^2 / (Gamma_peg^2 + K_peg^2)
     where K_peg = 0.45 chains/nm^2.
     The effective diffusion coefficient in CF mucus relative to aqueous buffer D_eff / D_0 satisfies an Ogston-like
     steric mesh obstruction coupled to adhesive interaction attenuation:
       D_eff / D_0 = exp(- (pi / 4.0) * ((d_p + d_f) / (xi_m + d_f))^2 * (1.0 - sigma_peg * S_peg)) *
                     (1.0 / (1.0 + kappa_adh * (1.0 - S_peg) * (c_m / c_m0)))
     where d_f = 6.0 nm is the mucin fiber diameter, xi_m is the CF mucus pore size (nm), c_m is the mucin solids
     concentration (mg/mL with baseline c_m0 = 30.0 mg/mL), sigma_peg = 0.82 is the steric shielding efficiency,
     and kappa_adh = 8.5 is the mucin hydrophobic adhesion coefficient.
     The steady-state mucus penetration probability P_mucus_pen across the mucus layer to the periciliary liquid is:
       P_mucus_pen = (D_eff / D_0) * (xi_m^2 / (xi_m^2 + (d_p / (1.0 + 0.5 * Gamma_peg))^2)).

  2. Apical Epithelial Endocytosis & Cytosolic Bioavailability:
     Upon reaching the apical plasma membrane of airway epithelial cells, hLPNs undergo cellular uptake
     (predominantly clathrin-mediated endocytosis). While dense PEG enhances mucus penetration, excessive PEGylation
     creates steric hindrance against cell membrane association (the "PEG dilemma"):
       k_uptake = k_up0 * (C_apical / (K_m_up + C_apical)) * (1.0 / (1.0 + omega_peg * Gamma_peg))
     where k_up0 = 1.85 h^-1, K_m_up = 25.0 ug/mL, and omega_peg = 0.28 (chains/nm^2)^-1.
     Following endosomal internalization, the lipid shell facilitates membrane fusion and endosomal escape:
       eta_escape = eta_max * (phi_lipid / (phi_lipid + K_lipid)) * (1.0 / (1.0 + beta_core * (d_p / 100.0)))
     where eta_max = 0.58, phi_lipid is the helper/fusogenic lipid molar fraction, K_lipid = 0.12, and
     beta_core = 0.22.
     The net intracellular functional oligonucleotide delivery index Phi_net satisfies:
       Phi_net = P_mucus_pen * (k_uptake / k_up0) * eta_escape * f_payload_intact
     where f_payload_intact = 0.92 is the protected encapsulated payload fraction.

  3. Kahn-Kalai Curative Threshold for hLPN Mucus Transcytosis & Delivery:
     - Effective mucus diffusivity ratio D_eff / D_0 >= 0.15 (>= 15.0% aqueous mobility in CF mucus).
     - Mucus layer penetration probability P_mucus_pen >= 0.20 (>= 20.0% penetration).
     - Endosomal escape efficiency eta_escape >= 0.25 (>= 25.0% escape).
     - Net intracellular functional delivery index Phi_net >= 0.040 (>= 4.0% applied aerosol dose delivered).

  4. Invariants & Deterministic Bounds:
     - Range conservation: 0.0 <= D_eff/D_0, P_mucus_pen, eta_escape, Phi_net <= 1.0.
     - Directional derivatives: d(D_eff/D_0)/d(Gamma_peg) > 0; d(D_eff/D_0)/d(c_m) < 0; d(P_mucus_pen)/d(d_p) < 0.
     - Complex-step differentiation: d(Phi_net)/d(Gamma_peg) verified at step h = 1e-30 against central finite
       differences with relative error <= 1e-6.
     - Poisoned negative controls: negative PEG density, negative nanoparticle diameter, negative pore size,
       or negative lipid fraction fail deterministically with ValueError.

CF Biophysical Domain Mapping (Column B - Popperian Falsifiable Hypothesis):
  - In Dense PEG-Shielded hLPN Curative Regimen (Gamma_peg=1.35, d_p=75.0 nm, xi_m=65.0 nm, c_m=30.0 mg/mL):
    D_eff/D_0 = 0.4135 >= 0.15, P_mucus_pen = 0.2805 >= 0.20, eta_escape = 0.3830 >= 0.25,
    Phi_net = 0.0493 >= 0.040 (Pass).
  - In Ultra-Small Core hLPN High-Penetration Curative Regimen (Gamma_peg=1.65, d_p=60.0 nm, xi_m=65.0 nm, c_m=28.0 mg/mL):
    D_eff/D_0 = 0.5497, P_mucus_pen = 0.4377, eta_escape = 0.4045, Phi_net = 0.0805 (Pass).
  - In Unshielded Hydrophobic PLGA Nanoparticle Regimen (Gamma_peg=0.02, d_p=120.0 nm, xi_m=32.0 nm, c_m=65.0 mg/mL):
    D_eff/D_0 = 0.0000 < 0.15, P_mucus_pen = 0.0000 < 0.20, trapped by mucin hydrophobic adhesion (Fail).
  - In Oversized Aggregated Nanoparticle Regimen (Gamma_peg=0.60, d_p=220.0 nm, xi_m=25.0 nm, c_m=70.0 mg/mL):
    P_mucus_pen = 0.0000 < 0.20, excluded by mucin steric mesh (Fail).
  - In Lipid-Deficient Endosomal Entrapment Regimen (Gamma_peg=1.10, d_p=90.0 nm, phi_lipid=0.02):
    eta_escape = 0.0692 < 0.25, Phi_net = 0.0026 < 0.040, failing cytosolic delivery (Fail).

Research simulation only; not medical advice.
"""
from __future__ import annotations

import cmath
import math
from dataclasses import dataclass
from typing import Any, Dict, Tuple

PRIMARY_CITATIONS: Tuple[str, ...] = (
    "DOI: 10.1021/acsami.1c14975",
    "DOI: 10.1038/s41598-021-85549-z",
)

# Biophysical Constants for hLPN & Airway Mucus Transport
DEFAULT_MUCIN_FIBER_DIAMETER_NM: float = 6.0
DEFAULT_BASELINE_MUCIN_CONC_MG_ML: float = 30.0
DEFAULT_K_PEG_SHIELD: float = 0.45
DEFAULT_SIGMA_PEG: float = 0.82
DEFAULT_KAPPA_ADHESION: float = 8.5
DEFAULT_K_UPTAKE_0_PER_H: float = 1.85
DEFAULT_KM_UPTAKE_UG_ML: float = 25.0
DEFAULT_OMEGA_PEG_INHIBITION: float = 0.28
DEFAULT_ETA_MAX_ESCAPE: float = 0.58
DEFAULT_K_LIPID_ESCAPE: float = 0.12
DEFAULT_BETA_CORE_ESCAPE: float = 0.22
DEFAULT_PAYLOAD_INTACT_FRACTION: float = 0.92

PARAMETER_PROVENANCE: Dict[str, Dict[str, object]] = {
    "peg_steric_shielding_half_saturation": {
        "value": DEFAULT_K_PEG_SHIELD,
        "tag": "phenomenological_heuristic",
        "source": (
            "Analytical half-saturation PEG grafting density K_peg parameterizing transition to dense steric brush "
            "regime (exact 0.45 chains/nm^2 scalar not stated in Conte et al., ACS Appl Mater Interfaces 2022, "
            "DOI: 10.1021/acsami.1c14975)"
        ),
        "bounds": (0.20, 0.90),
    },
    "steric_shielding_transmission_factor": {
        "value": DEFAULT_SIGMA_PEG,
        "tag": "phenomenological_heuristic",
        "source": (
            "Analytical steric mesh shielding transmission efficiency sigma_peg preventing mucin fiber entanglement "
            "(exact 0.82 dimensionless factor not stated in Conte et al., ACS Appl Mater Interfaces 2022, "
            "DOI: 10.1021/acsami.1c14975)"
        ),
        "bounds": (0.50, 0.95),
    },
    "mucin_hydrophobic_adhesion_coefficient": {
        "value": DEFAULT_KAPPA_ADHESION,
        "tag": "phenomenological_heuristic",
        "source": (
            "Analytical hydrophobic adhesion coefficient kappa_adh for unshielded PLGA cores in CF sputum "
            "(exact 8.5 dimensionless scalar not stated in Conte et al., ACS Appl Mater Interfaces 2022, "
            "DOI: 10.1021/acsami.1c14975)"
        ),
        "bounds": (3.0, 18.0),
    },
    "peg_endocytosis_hindrance_coefficient": {
        "value": DEFAULT_OMEGA_PEG_INHIBITION,
        "tag": "phenomenological_heuristic",
        "source": (
            "Analytical PEG dilemma membrane association attenuation coefficient omega_peg "
            "(exact 0.28 (chains/nm^2)^-1 scalar not stated in Comegna et al., Sci Rep 2021, "
            "DOI: 10.1038/s41598-021-85549-z)"
        ),
        "bounds": (0.10, 0.60),
    },
    "lipid_endosomal_escape_half_saturation": {
        "value": DEFAULT_K_LIPID_ESCAPE,
        "tag": "phenomenological_heuristic",
        "source": (
            "Analytical lipid fraction half-saturation K_lipid for endosomal membrane destabilization "
            "(exact 0.12 molar fraction not stated in Comegna et al., Sci Rep 2021, "
            "DOI: 10.1038/s41598-021-85549-z)"
        ),
        "bounds": (0.05, 0.25),
    },
}

TWO_COLUMN_EPISTEMIC_HYGIENE: Dict[str, Dict[str, str]] = {
    "deductive_column_a": {
        "formal_theorem": (
            "Coupling dense PEG steric brush shielding Gamma_peg to Ogston mesh obstruction and endosomal membrane "
            "fusion yields a smooth, non-negative delivery operator Phi_net(Gamma_peg, d_p, xi_m, phi_lipid) on [0, 1] "
            "with d(D_eff/D_0)/d(Gamma_peg) > 0 and d(P_mucus_pen)/d(d_p) < 0, resolving the trade-off between mucus "
            "penetration and cellular endocytosis."
        ),
        "invariants": (
            "0 <= D_eff/D_0, P_mucus_pen, eta_escape, Phi_net <= 1; d(D_eff/D_0)/d(Gamma_peg) > 0; "
            "d(D_eff/D_0)/d(c_m) < 0; d(P_mucus_pen)/d(d_p) < 0; complex-step relative error <= 1e-6."
        ),
    },
    "popperian_column_b": {
        "falsifiable_hypothesis": (
            "Aerosol delivery of hybrid lipid-polymer nanoparticles achieves >= 4.0% net cytosolic delivery in CF airway "
            "epithelia only when PEG grafting density >= 0.8 chains/nm^2 shields hydrophobic adhesion, nanoparticle diameter "
            "<= 100 nm permits mesh penetration, and fusogenic lipid fraction >= 20% mediates endosomal escape."
        ),
        "empirical_falsification_threshold": (
            "Falsified if unshielded hydrophobic nanoparticles (Gamma_peg < 0.1) achieve effective mucus diffusion "
            "D_eff/D_0 >= 0.15 in human CF sputum, or if dense PEGylation (Gamma_peg >= 1.2) completely abolishes epithelial "
            "internalization without delivering functional oligonucleotides."
        ),
    },
}

FIVE_LAYER_LADDER: Dict[str, str] = {
    "cool_human_fact": (
        "In cystic fibrosis lungs, mucus is up to 60 times thicker than normal, acting like an impenetrable microscopic "
        "spiderweb that catches inhaled medicine before it ever reaches the airway cells."
    ),
    "everyday_analogy": (
        "It's like wrapping a sticky tennis ball in a slick, water-loving silicone coat: without the coat, the ball gets "
        "glued to the first Velcro fiber it touches; with the slick coat, it slips effortlessly right through the netting."
    ),
    "engineering_lens": (
        "hLPNs balance a biophysical optimization trade-off known as the PEG dilemma: dense PEG brush grafting minimizes "
        "mucin hydrophobic adhesion and steric collision, while helper lipid shells ensure endosomal membrane disruption."
    ),
    "formal_theorem": (
        "The composite transfer operator Phi_net = P_pen(Gamma_peg, d_p, xi_m) * (k_up/k_0) * eta_esc(phi_lipid, d_p) * f_p "
        "is real-analytic on (0, inf)^n and possesses well-defined complex-step derivatives enabling sensitivity optimization."
    ),
    "falsifiable_code_and_toy": (
        "evaluate_round201_hlpn_mucus_transcytosis(regime_key) evaluates steric shielding, mucus penetration, "
        "and cytosolic delivery across 5 clinical nanoparticle formulations."
    ),
}

REGIME_PRESETS: Dict[str, Dict[str, float]] = {
    "hlpn_dense_peg_optimized_curative": {
        "gamma_peg": 1.35,
        "d_p_nm": 75.0,
        "xi_m_nm": 65.0,
        "c_m_mg_ml": 30.0,
        "c_apical_ug_ml": 55.0,
        "phi_lipid": 0.40,
    },
    "hlpn_ultra_small_cleavable_curative": {
        "gamma_peg": 1.65,
        "d_p_nm": 60.0,
        "xi_m_nm": 65.0,
        "c_m_mg_ml": 28.0,
        "c_apical_ug_ml": 65.0,
        "phi_lipid": 0.45,
    },
    "unshielded_plga_hydrophobic_adhesion_cf": {
        "gamma_peg": 0.02,
        "d_p_nm": 120.0,
        "xi_m_nm": 32.0,
        "c_m_mg_ml": 65.0,
        "c_apical_ug_ml": 20.0,
        "phi_lipid": 0.10,
    },
    "oversized_mesh_excluded_aggregate_cf": {
        "gamma_peg": 0.60,
        "d_p_nm": 220.0,
        "xi_m_nm": 25.0,
        "c_m_mg_ml": 70.0,
        "c_apical_ug_ml": 30.0,
        "phi_lipid": 0.25,
    },
    "lipid_deficient_endosomal_entrapment": {
        "gamma_peg": 1.10,
        "d_p_nm": 90.0,
        "xi_m_nm": 50.0,
        "c_m_mg_ml": 38.0,
        "c_apical_ug_ml": 45.0,
        "phi_lipid": 0.02,
    },
}


@dataclass(frozen=True)
class Round201HLPNDeliveryResult:
    regime_key: str
    gamma_peg: float
    d_p_nm: float
    xi_m_nm: float
    c_m_mg_ml: float
    c_apical_ug_ml: float
    phi_lipid: float
    s_peg_shielding: float
    d_eff_over_d0: float
    p_mucus_pen: float
    k_uptake_per_h: float
    eta_escape: float
    phi_net_delivery: float
    crosses_kahn_kalai_threshold: bool
    falsification_invariant_satisfied: bool
    complex_step_relative_error: float


def compute_peg_shielding_factor(gamma_peg: float | complex) -> float | complex:
    """Computes the dense PEG brush steric shielding factor S_peg on [0, 1]."""
    if isinstance(gamma_peg, (int, float)) and gamma_peg < 0.0:
        raise ValueError(f"gamma_peg must be non-negative, got {gamma_peg}")
    k2 = DEFAULT_K_PEG_SHIELD**2
    g2 = gamma_peg**2
    return g2 / (g2 + k2)


def compute_mucus_diffusivity_and_penetration(
    gamma_peg: float | complex,
    d_p_nm: float | complex,
    xi_m_nm: float | complex,
    c_m_mg_ml: float | complex,
) -> Tuple[float | complex, float | complex, float | complex]:
    """Computes S_peg, D_eff / D_0, and mucus penetration probability P_mucus_pen."""
    if isinstance(gamma_peg, (int, float)) and gamma_peg < 0.0:
        raise ValueError(f"gamma_peg must be non-negative, got {gamma_peg}")
    if isinstance(d_p_nm, (int, float)) and d_p_nm <= 0.0:
        raise ValueError(f"d_p_nm must be positive, got {d_p_nm}")
    if isinstance(xi_m_nm, (int, float)) and xi_m_nm <= 0.0:
        raise ValueError(f"xi_m_nm must be positive, got {xi_m_nm}")
    if isinstance(c_m_mg_ml, (int, float)) and c_m_mg_ml < 0.0:
        raise ValueError(f"c_m_mg_ml must be non-negative, got {c_m_mg_ml}")

    s_peg = compute_peg_shielding_factor(gamma_peg)
    df = DEFAULT_MUCIN_FIBER_DIAMETER_NM
    mesh_ratio = (d_p_nm + df) / (xi_m_nm + df)
    steric_exponent = -(math.pi / 4.0) * (mesh_ratio**2) * (1.0 - DEFAULT_SIGMA_PEG * s_peg)

    # Exponentiation supporting complex-step arithmetic
    if isinstance(steric_exponent, complex):
        steric_trans = cmath.exp(steric_exponent)
    else:
        steric_trans = math.exp(steric_exponent)

    adh_denom = 1.0 + DEFAULT_KAPPA_ADHESION * (1.0 - s_peg) * (c_m_mg_ml / DEFAULT_BASELINE_MUCIN_CONC_MG_ML)
    d_eff_ratio = steric_trans / adh_denom

    eff_dp = d_p_nm / (1.0 + 0.5 * gamma_peg)
    p_pen = d_eff_ratio * (xi_m_nm**2 / (xi_m_nm**2 + eff_dp**2))
    return s_peg, d_eff_ratio, p_pen


def compute_uptake_and_escape(
    gamma_peg: float | complex,
    d_p_nm: float | complex,
    c_apical_ug_ml: float | complex,
    phi_lipid: float | complex,
) -> Tuple[float | complex, float | complex]:
    """Computes apical endocytosis rate k_uptake and endosomal escape efficiency eta_escape."""
    if isinstance(c_apical_ug_ml, (int, float)) and c_apical_ug_ml < 0.0:
        raise ValueError(f"c_apical_ug_ml must be non-negative, got {c_apical_ug_ml}")
    if isinstance(phi_lipid, (int, float)) and (phi_lipid < 0.0 or phi_lipid > 1.0):
        raise ValueError(f"phi_lipid must be in [0, 1], got {phi_lipid}")

    c_sat = c_apical_ug_ml / (DEFAULT_KM_UPTAKE_UG_ML + c_apical_ug_ml)
    peg_inhibition = 1.0 / (1.0 + DEFAULT_OMEGA_PEG_INHIBITION * gamma_peg)
    k_up = DEFAULT_K_UPTAKE_0_PER_H * c_sat * peg_inhibition

    lipid_sat = phi_lipid / (phi_lipid + DEFAULT_K_LIPID_ESCAPE)
    core_penalty = 1.0 / (1.0 + DEFAULT_BETA_CORE_ESCAPE * (d_p_nm / 100.0))
    eta_esc = DEFAULT_ETA_MAX_ESCAPE * lipid_sat * core_penalty

    return k_up, eta_esc


def compute_net_delivery_index(
    p_pen: float | complex,
    k_up: float | complex,
    eta_esc: float | complex,
) -> float | complex:
    """Computes net cytosolic oligonucleotide delivery fraction Phi_net."""
    norm_up = k_up / DEFAULT_K_UPTAKE_0_PER_H
    return p_pen * norm_up * eta_esc * DEFAULT_PAYLOAD_INTACT_FRACTION


def complex_step_gamma_delivery_derivative(
    gamma_peg: float,
    d_p_nm: float,
    xi_m_nm: float,
    c_m_mg_ml: float,
    c_apical_ug_ml: float,
    phi_lipid: float,
    h: float = 1e-30,
) -> Tuple[float, float, float]:
    """Verifies d(Phi_net)/d(Gamma_peg) via complex-step perturbation at step h=1e-30."""
    gamma_c = complex(gamma_peg, h)
    _, _, p_pen_c = compute_mucus_diffusivity_and_penetration(gamma_c, d_p_nm, xi_m_nm, c_m_mg_ml)
    k_up_c, eta_esc_c = compute_uptake_and_escape(gamma_c, d_p_nm, c_apical_ug_ml, phi_lipid)
    phi_c = compute_net_delivery_index(p_pen_c, k_up_c, eta_esc_c)
    cs_deriv = phi_c.imag / h

    # Central finite differences at delta = 1e-5
    delta = 1e-5
    _, _, p_plus = compute_mucus_diffusivity_and_penetration(gamma_peg + delta, d_p_nm, xi_m_nm, c_m_mg_ml)
    k_plus, eta_plus = compute_uptake_and_escape(gamma_peg + delta, d_p_nm, c_apical_ug_ml, phi_lipid)
    phi_plus = float(compute_net_delivery_index(p_plus, k_plus, eta_plus))

    _, _, p_minus = compute_mucus_diffusivity_and_penetration(gamma_peg - delta, d_p_nm, xi_m_nm, c_m_mg_ml)
    k_minus, eta_minus = compute_uptake_and_escape(gamma_peg - delta, d_p_nm, c_apical_ug_ml, phi_lipid)
    phi_minus = float(compute_net_delivery_index(p_minus, k_minus, eta_minus))

    fd_deriv = (phi_plus - phi_minus) / (2.0 * delta)
    denom = max(abs(cs_deriv), abs(fd_deriv), 1e-12)
    rel_err = abs(cs_deriv - fd_deriv) / denom

    return cs_deriv, fd_deriv, rel_err


def evaluate_round201_hlpn_mucus_transcytosis(regime_key: str) -> Round201HLPNDeliveryResult:
    """Evaluates hLPN mucus penetration and epithelial transcytosis for a named regimen."""
    if regime_key not in REGIME_PRESETS:
        raise KeyError(f"Unknown regime '{regime_key}'. Valid regimes: {sorted(REGIME_PRESETS)}")

    p = REGIME_PRESETS[regime_key]
    gamma_peg = p["gamma_peg"]
    d_p = p["d_p_nm"]
    xi_m = p["xi_m_nm"]
    c_m = p["c_m_mg_ml"]
    c_apical = p["c_apical_ug_ml"]
    phi_lipid = p["phi_lipid"]

    s_peg_val, d_eff_val, p_pen_val = compute_mucus_diffusivity_and_penetration(gamma_peg, d_p, xi_m, c_m)
    k_up_val, eta_esc_val = compute_uptake_and_escape(gamma_peg, d_p, c_apical, phi_lipid)
    phi_net_val = compute_net_delivery_index(p_pen_val, k_up_val, eta_esc_val)

    s_peg = float(s_peg_val)
    d_eff = float(d_eff_val)
    p_pen = float(p_pen_val)
    k_up = float(k_up_val)
    eta_esc = float(eta_esc_val)
    phi_net = float(phi_net_val)

    eval_gamma = max(gamma_peg, 0.40)
    _, _, rel_err = complex_step_gamma_delivery_derivative(eval_gamma, d_p, xi_m, c_m, c_apical, phi_lipid)

    inv_satisfied = (
        0.0 <= s_peg <= 1.0
        and 0.0 <= d_eff <= 1.0
        and 0.0 <= p_pen <= 1.0
        and 0.0 <= eta_esc <= 1.0
        and 0.0 <= phi_net <= 1.0
        and rel_err <= 1e-6
    )

    crosses_threshold = (
        d_eff >= 0.15
        and p_pen >= 0.20
        and eta_esc >= 0.25
        and phi_net >= 0.040
    )

    return Round201HLPNDeliveryResult(
        regime_key=regime_key,
        gamma_peg=gamma_peg,
        d_p_nm=d_p,
        xi_m_nm=xi_m,
        c_m_mg_ml=c_m,
        c_apical_ug_ml=c_apical,
        phi_lipid=phi_lipid,
        s_peg_shielding=s_peg,
        d_eff_over_d0=d_eff,
        p_mucus_pen=p_pen,
        k_uptake_per_h=k_up,
        eta_escape=eta_esc,
        phi_net_delivery=phi_net,
        crosses_kahn_kalai_threshold=crosses_threshold,
        falsification_invariant_satisfied=inv_satisfied,
        complex_step_relative_error=rel_err,
    )


def evaluate_falsification_negative_control(
    poison_negative_gamma: bool = False,
    poison_negative_diameter: bool = False,
    poison_negative_pore: bool = False,
    poison_invalid_lipid: bool = False,
) -> Dict[str, Any]:
    """Tests poisoned negative controls that must deterministically raise ValueError."""
    reasons: list[str] = []
    if poison_negative_gamma:
        try:
            compute_peg_shielding_factor(-0.5)
            return {"passed": False, "reason": "Negative PEG density accepted unexpectedly"}
        except ValueError:
            reasons.append("Deterministically rejected negative PEG density")

    if poison_negative_diameter:
        try:
            compute_mucus_diffusivity_and_penetration(1.0, -80.0, 50.0, 30.0)
            return {"passed": False, "reason": "Negative nanoparticle diameter accepted unexpectedly"}
        except ValueError:
            reasons.append("Deterministically rejected negative nanoparticle diameter")

    if poison_negative_pore:
        try:
            compute_mucus_diffusivity_and_penetration(1.0, 80.0, -20.0, 30.0)
            return {"passed": False, "reason": "Negative mucus pore size accepted unexpectedly"}
        except ValueError:
            reasons.append("Deterministically rejected negative mucus pore size")

    if poison_invalid_lipid:
        try:
            compute_uptake_and_escape(1.0, 80.0, 40.0, 1.45)
            return {"passed": False, "reason": "Lipid fraction > 1.0 accepted unexpectedly"}
        except ValueError:
            reasons.append("Deterministically rejected lipid fraction > 1.0")

    return {"passed": True, "reason": "; ".join(reasons) if reasons else "Nominal baseline"}
