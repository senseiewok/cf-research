#!/usr/bin/env python3
"""
Round 200: Bicellular Tight Junction Claudin-2/Claudin-4 Paracellular Ion Permselectivity Invariant Bridge (T-0889).

Mathematical & Physical Pillar (Column A - Deductive Theorem):
  1. Claudin-2 / Claudin-4 Paracellular Permselectivity & Nernst-Planck Shunt Resistance:
     Human conducting airway tight junctions regulate paracellular ion and water movement via opposing claudin
     strands: Claudin-2 forms cation-selective, low-resistance paracellular pores (increasing Na+ permeability
     P_Na and P_Na/P_HCO3 selectivity), whereas Claudin-4 tightens the bicellular barrier and restricts paracellular
     Na+ leak, raising tight-junction resistance R_TJ (Schlingmann et al., Semin Cell Dev Biol 2015,
     DOI: 10.1016/j.semcdb.2015.04.009). Proinflammatory cytokines (TNF-alpha / IFN-gamma) disrupt tight-junction
     strands and open non-selective paracellular leak pathways I_cyt >= 0 (Coyne et al., Mol Biol Cell 2002,
     DOI: 10.1091/mbc.e02-03-0134).
     Given normalized Claudin-2 expression x_cldn2 >= 0, Claudin-4 expression x_cldn4 >= 0, and cytokine strand
     disruption I_cyt >= 0, paracellular permeabilities (in 10^-6 cm/s) are:
       P_Na   = (P_Na0 + alpha_cldn2 * x_cldn2) / (1.0 + beta_cldn4 * x_cldn4) + gamma_na * I_cyt
       P_HCO3 = P_HCO3_0 / (1.0 + 0.28 * x_cldn4 + 0.30 * x_cldn2) + gamma_hco3 * I_cyt
       P_Cl   = P_Cl0 / (1.0 + 0.32 * x_cldn4) + gamma_cl * I_cyt
     By the Nernst-Planck conductance relation G_i = (F^2 / RT) * P_i * C_bar at T = 310.15 K
     (where F^2 / RT = 0.0036102 mS*cm / (10^-6 cm/s * mM)):
       G_TJ = 0.0036102 * (P_Na * 122.5 + P_HCO3 * 18.0 + P_Cl * 115.0)   [mS/cm^2]
       R_TJ = 1000.0 / G_TJ                                               [Ohm*cm^2].

  2. Goldman-Hodgkin-Katz (GHK) Paracellular Electrodiffusion & Steady-State ASL pH Tempering:
     Across transepithelial potential V_te (mV, lumen relative to basolateral bath) with thermal voltage
     V_T = RT / F = 26.726 mV at 310.15 K, let u_i = z_i * V_te / V_T (z_Na = +1, z_HCO3 = -1) and define the
     complex-analytic Bernoulli factor B(u) = u / (exp(u) - 1) (with Taylor expansion 1 - u/2 + u^2/12 - u^4/720
     for |u| < 1e-5). Basolateral-to-apical paracellular flux (nmol * cm^-2 * h^-1) satisfies:
       J_i_para = 3.6 * P_i * B(u_i) * (C_i_baso - C_i_asl * exp(u_i)).
     Balancing active apical CFTR bicarbonate transport J_hco3_trans against unopposed ATP12A proton load
     J_proton_atp12a = 24.0 nmol * cm^-2 * h^-1 and paracellular GHK exchange yields the steady-state luminal
     bicarbonate concentration C_HCO3_asl (mM) and Henderson-Hasselbalch pH_ASL (Thornell et al., J Physiol 2020,
     DOI: 10.1113/jp280120):
       raw_num    = C_HCO3_baso + (J_hco3_trans - J_proton_atp12a) / (3.6 * P_HCO3 * B(u_HCO3))
       C_HCO3_asl = 0.5 * (raw_num + sqrt(raw_num^2 + 1.0)) / exp(u_HCO3)
       pH_ASL     = 6.10 + log10(C_HCO3_asl / 1.20)
     where C_HCO3_baso = 24.0 mM and dissolved CO2 buffer is 1.20 mM (40 mmHg * 0.03 mM/mmHg).
     In untreated hyperpolarized CF airway (J_hco3_trans = 2.0 < J_proton_atp12a = 24.0, V_te = -24.0 mV),
     paracellular HCO3- influx (J_HCO3_para = +21.95 nmol * cm^-2 * h^-1 > 0) flows from submucosa into ASL to
     temper acidification around pH_ASL = 6.874 rather than unbuffered collapse (Thornell et al., J Physiol 2020).
     Steady-state ASL hydration depth h_ASL (um) reflects active anion transport retained by intact R_TJ:
       h_ASL = 2.4 + 6.5 * (R_TJ / (R_TJ + 220.0)) * (J_hco3_trans / (J_hco3_trans + 18.0)) / (1.0 + (J_Na_para / 650.0)^2).

  3. Kahn-Kalai Curative Threshold for Bicellular Tight-Junction Barrier & ASL Alkalinization:
     - Paracellular tight-junction resistance R_TJ >= 500.0 Ohm*cm^2.
     - Steady-state luminal pH_ASL >= 7.25 (C_HCO3_asl >= 16.0 mM).
     - Paracellular sodium leak magnitude |J_Na_para| <= 450.0 nmol * cm^-2 * h^-1.
     - Steady-state airway surface liquid height h_ASL >= 5.20 um.

  4. Invariants & Deterministic Bounds:
     - Selectivity & barrier derivatives: d(P_Na)/d(x_cldn2) > 0; d(P_Na / P_HCO3)/d(x_cldn2) > 0;
       d(R_TJ)/d(x_cldn4) > 0; d(R_TJ)/d(I_cyt) < 0.
     - Removable GHK singularity at V_te -> 0 mV smoothly recovers Fick's first law without NaN/ZeroDivisionError.
     - Complex-step differentiation: d(pH_ASL)/d(J_hco3_trans) verified at step h = 1e-30 against central finite
       differences with relative error <= 1e-6.
     - Poisoned negative controls: negative claudin expression, negative cytokine index, or negative transcellular
       flux fail deterministically with ValueError.

CF Biophysical Domain Mapping (Column B - Popperian Falsifiable Hypothesis):
  - In Claudin-4 Intact Barrier + CFTR Rescued Curative Regimen (x_cldn2=0.25, x_cldn4=1.35, I_cyt=0.0, J_trans=56.0, V_te=-8.5 mV):
    P_Na = 1.293, P_HCO3 = 0.998, R_TJ = 863.2 Ohm*cm^2 >= 500, C_HCO3_asl = 25.10 mM, pH_ASL = 7.420 >= 7.25,
    J_Na_para = 345.60 nmol/(cm^2*h) <= 450, h_ASL = 5.45 um >= 5.20 um (Pass).
  - In Claudin-4 Enhanced Barrier Stabilized Curative Regimen (x_cldn2=0.15, x_cldn4=1.80, I_cyt=0.0, J_trans=62.0, V_te=-7.5 mV):
    P_Na = 0.919, P_HCO3 = 0.936, R_TJ = 1062.4 Ohm*cm^2, C_HCO3_asl = 27.96 mM, pH_ASL = 7.467,
    J_Na_para = 230.18 nmol/(cm^2*h), h_ASL = 6.11 um (Pass).
  - In Untreated Hyperpolarized CF Airway Paracellular Tempering Regimen (x_cldn2=0.45, x_cldn4=0.90, I_cyt=0.15, J_trans=2.0, V_te=-24.0 mV):
    R_TJ = 474.4 Ohm*cm^2 < 500, J_HCO3_para = +21.95 nmol/(cm^2*h) (basolateral-to-apical tempering influx),
    C_HCO3_asl = 7.14 mM, pH_ASL = 6.874 < 7.25, h_ASL = 2.48 um < 5.20 um (Fail).
  - In Cytokine-Disrupted Leaky Tight Junction Regimen (x_cldn2=0.80, x_cldn4=0.20, I_cyt=2.20, J_trans=50.0, V_te=-14.0 mV):
    R_TJ = 93.1 Ohm*cm^2 < 500, C_HCO3_asl = 14.91 mM, pH_ASL = 7.194 < 7.25, J_Na_para = 5063.24 nmol/(cm^2*h),
    h_ASL = 2.42 um < 5.20 um, failing barrier and hydration thresholds despite active CFTR transport.
  - In Claudin-2 Overexpression Cation-Selective Shunt Leak Regimen (x_cldn2=2.60, x_cldn4=0.15, I_cyt=0.25, J_trans=42.0, V_te=-10.0 mV):
    P_Na / P_HCO3 = 8.45, R_TJ = 139.8 Ohm*cm^2 < 500, J_Na_para = 3944.01 nmol/(cm^2*h) > 450,
    h_ASL = 2.44 um < 5.20 um, failing barrier resistance and hydration retention.

Research simulation only; not medical advice.
"""
from __future__ import annotations

import cmath
import math
from dataclasses import dataclass
from typing import Any, Dict, Tuple

PRIMARY_CITATIONS: Tuple[str, ...] = (
    "DOI: 10.1113/jp280120",
    "DOI: 10.1016/j.semcdb.2015.04.009",
    "DOI: 10.1091/mbc.e02-03-0134",
)

# Physical & Paracellular Transport Constants
THERMAL_VOLTAGE_MV_310K: float = 26.726       # RT / F at 310.15 K (mV)
NERNST_PLANCK_G_SCALE_MS: float = 0.0036102   # F^2 / RT conversion to mS/cm^2 per (10^-6 cm/s * mM)
C_NA_BASO_MM: float = 140.0
C_NA_APICAL_MM: float = 105.0
C_HCO3_BASO_MM: float = 24.0
C_CO2_DISSOLVED_MM: float = 1.20              # 40 mmHg * 0.03 mM/mmHg Henry's law at 37 C
PKA_CARBONIC_37C: float = 6.10
DEFAULT_J_PROTON_ATP12A: float = 24.0         # Apical non-gastric H+/K+-ATPase acid load (nmol * cm^-2 * h^-1)

DEFAULT_P_NA_BASE: float = 2.10
DEFAULT_ALPHA_CLDN2: float = 4.80
DEFAULT_BETA_CLDN4: float = 1.15
DEFAULT_GAMMA_NA_CYT: float = 4.20
DEFAULT_P_HCO3_BASE: float = 1.45
DEFAULT_GAMMA_HCO3_CYT: float = 3.20
DEFAULT_P_CL_BASE: float = 1.80
DEFAULT_GAMMA_CL_CYT: float = 3.60

PARAMETER_PROVENANCE: Dict[str, Dict[str, object]] = {
    "thermal_voltage_rt_over_f_310k": {
        "value": THERMAL_VOLTAGE_MV_310K,
        "tag": "derived",
        "source": (
            "Thermodynamic thermal voltage RT/F at physiological airway temperature T = 310.15 K "
            "(CODATA R = 8.31446 J/(mol*K), F = 96485.3 C/mol; exact 26.726 mV rounding not stated in "
            "Thornell et al., J Physiol 2020, DOI: 10.1113/jp280120)"
        ),
        "bounds": (26.5, 27.0),
    },
    "claudin2_cation_pore_gain": {
        "value": DEFAULT_ALPHA_CLDN2,
        "tag": "phenomenological_heuristic",
        "source": (
            "Analytical cation-selective paracellular Na+ permeability gain per unit Claudin-2 expression "
            "(exact 4.80 x 10^-6 cm/s scalar not stated in Schlingmann et al., Semin Cell Dev Biol 2015, "
            "DOI: 10.1016/j.semcdb.2015.04.009)"
        ),
        "bounds": (2.0, 8.0),
    },
    "claudin4_barrier_sealing_coefficient": {
        "value": DEFAULT_BETA_CLDN4,
        "tag": "phenomenological_heuristic",
        "source": (
            "Analytical bicellular cation-sealing attenuation factor per unit Claudin-4 expression "
            "(exact 1.15 dimensionless scalar not stated in Schlingmann et al., Semin Cell Dev Biol 2015, "
            "DOI: 10.1016/j.semcdb.2015.04.009)"
        ),
        "bounds": (0.50, 2.50),
    },
    "cytokine_strand_leak_conductance_gain": {
        "value": DEFAULT_GAMMA_HCO3_CYT,
        "tag": "phenomenological_heuristic",
        "source": (
            "Analytical non-selective paracellular permeability increment under TNF-alpha/IFN-gamma tight-junction "
            "remodeling (exact 3.20 x 10^-6 cm/s scalar not stated in Coyne et al., Mol Biol Cell 2002, "
            "DOI: 10.1091/mbc.e02-03-0134)"
        ),
        "bounds": (1.0, 6.0),
    },
    "baseline_paracellular_bicarbonate_permeability": {
        "value": DEFAULT_P_HCO3_BASE,
        "tag": "phenomenological_heuristic",
        "source": (
            "Analytical baseline airway epithelial paracellular bicarbonate permeability constant P_HCO3_0 "
            "(exact 1.45 x 10^-6 cm/s scalar not stated in Thornell et al., J Physiol 2020, DOI: 10.1113/jp280120)"
        ),
        "bounds": (0.50, 3.00),
    },
}

TWO_COLUMN_EPISTEMIC_HYGIENE: Dict[str, Dict[str, str]] = {
    "deductive_column_a": {
        "formal_theorem": (
            "Under Goldman-Hodgkin-Katz electrodiffusion across a bicellular tight junction with ion permeabilities "
            "P_Na(x_cldn2, x_cldn4, I_cyt) and P_HCO3(x_cldn2, x_cldn4, I_cyt), the Bernoulli-regularized flux "
            "J_i = 3.6 * P_i * B(z_i V_te / V_T) * (C_i,baso - C_i,asl * exp(z_i V_te / V_T)) is entire in V_te "
            "(including V_te = 0), strictly recovers Nernst-Donnan equilibrium when active transport vanishes, and "
            "satisfies d(P_Na / P_HCO3)/d(x_cldn2) > 0 and d(R_TJ)/d(x_cldn4) > 0."
        ),
        "invariants": (
            "P_Na > 0; P_HCO3 > 0; R_TJ > 0; C_HCO3_asl > 0; d(P_Na)/d(x_cldn2) > 0; "
            "d(P_Na / P_HCO3)/d(x_cldn2) > 0; d(R_TJ)/d(x_cldn4) > 0; smooth limit as V_te -> 0."
        ),
    },
    "popperian_column_b": {
        "falsifiable_hypothesis": (
            "In hyperpolarized CF airway epithelia lacking apical CFTR bicarbonate transport, paracellular HCO3- "
            "influx from the 24 mM basolateral bath tempers unopposed ATP12A acidification near pH 6.8–6.9, whereas "
            "durable restoration of pH_ASL >= 7.25 and ASL hydration >= 5.2 um requires both restored transcellular "
            "anion transport and an intact Claudin-4-sealed tight junction (R_TJ >= 500 Ohm*cm^2)."
        ),
        "empirical_falsification_threshold": (
            "Falsified if cytokine-disrupted airway epithelia with R_TJ < 120 Ohm*cm^2 sustain physiological "
            "ASL hydration height >= 5.2 um against hypertonic paracellular Na+ shunt leak, or if removing "
            "basolateral bicarbonate fails to further acidify untreated CF airway surface liquid."
        ),
    },
}

FIVE_LAYER_LADDER: Dict[str, str] = {
    "cool_human_fact": (
        "Airway cells are zipped together at their top edges by molecular Velcro proteins called claudins: "
        "Claudin-4 acts like a waterproof gasket, while Claudin-2 opens tiny cation-selective pores—and in CF, "
        "small amounts of bicarbonate sneaking between cells actually keep airway liquid from turning dangerously acidic."
    ),
    "everyday_analogy": (
        "Imagine trying to refill and balance the pH of a rooftop swimming pool: turning on the clean water pump "
        "(CFTR) only works if the tile grout between panels (Claudin-4) is sealed; if inflammation cracks the grout, "
        "water and buffers leak right back through the seams."
    ),
    "engineering_lens": (
        "Epithelial ion homeostasis is a parallel circuit of active transcellular pumps (CFTR, ATP12A) shunted by "
        "a permselective paracellular resistor R_TJ: Claudin-2 increases cation-over-anion selectivity P_Na/P_HCO3, "
        "while Claudin-4 maximizes shunt impedance to preserve transepithelial osmotic gradients."
    ),
    "formal_theorem": (
        "Let B(u) = u / (e^u - 1) be the exponential Bernoulli generating kernel for dimensionless potential "
        "u_i = z_i F V_te / (RT). Stationary apical bicarbonate balance admits an exact closed-form root "
        "C_HCO3,asl(J_trans, P_HCO3, V_te) analytic across u = 0 with exact complex-step sensitivity at h = 1e-30."
    ),
    "falsifiable_code_and_toy": (
        "evaluate_round200_tight_junction_claudin_flux(regime_key) computes ion-specific paracellular permeabilities, "
        "R_TJ, GHK Na+ and HCO3- fluxes, steady-state pH_ASL, and hydration depth across 5 tight-junction states."
    ),
}

REGIME_PRESETS: Dict[str, Dict[str, float]] = {
    "claudin4_intact_barrier_cftr_rescued_curative": {
        "x_cldn2": 0.25,
        "x_cldn4": 1.35,
        "i_cytokine": 0.00,
        "j_hco3_transcellular": 56.0,
        "v_te_mv": -8.5,
    },
    "claudin4_enhanced_barrier_stabilized_curative": {
        "x_cldn2": 0.15,
        "x_cldn4": 1.80,
        "i_cytokine": 0.00,
        "j_hco3_transcellular": 62.0,
        "v_te_mv": -7.5,
    },
    "untreated_cf_hyperpolarized_paracellular_tempering": {
        "x_cldn2": 0.45,
        "x_cldn4": 0.90,
        "i_cytokine": 0.15,
        "j_hco3_transcellular": 2.0,
        "v_te_mv": -24.0,
    },
    "cytokine_disrupted_leaky_tight_junction": {
        "x_cldn2": 0.80,
        "x_cldn4": 0.20,
        "i_cytokine": 2.20,
        "j_hco3_transcellular": 50.0,
        "v_te_mv": -14.0,
    },
    "claudin2_overexpression_cation_shunt_leak": {
        "x_cldn2": 2.60,
        "x_cldn4": 0.15,
        "i_cytokine": 0.25,
        "j_hco3_transcellular": 42.0,
        "v_te_mv": -10.0,
    },
}


@dataclass(frozen=True)
class Round200TightJunctionResult:
    regime_key: str
    x_cldn2: float
    x_cldn4: float
    i_cytokine: float
    j_hco3_transcellular: float
    v_te_mv: float
    p_na_1e6_cm_s: float
    p_hco3_1e6_cm_s: float
    cation_anion_selectivity: float
    r_tj_ohm_cm2: float
    c_hco3_asl_mm: float
    j_hco3_para_nmol_cm2_h: float
    j_na_para_nmol_cm2_h: float
    ph_asl: float
    h_asl_um: float
    crosses_kahn_kalai_threshold: bool
    falsification_invariant_satisfied: bool
    complex_step_relative_error: float
    research_boundary: str = "Research simulation only; not medical advice."


def bernoulli_ghk_kernel(u: float | complex) -> float | complex:
    """Evaluates B(u) = u / (exp(u) - 1) smoothly across u = 0 for both real and complex arguments."""
    if abs(u) < 1e-5:
        return 1.0 - 0.5 * u + (u ** 2) / 12.0 - (u ** 4) / 720.0
    exp_fn = cmath.exp if isinstance(u, complex) else math.exp
    return u / (exp_fn(u) - 1.0)


def compute_claudin_permeabilities_and_resistance(
    x_cldn2: float | complex,
    x_cldn4: float | complex,
    i_cytokine: float | complex = 0.0,
) -> Tuple[float | complex, float | complex, float | complex]:
    """Computes (P_Na, P_HCO3, R_TJ) in (10^-6 cm/s, 10^-6 cm/s, Ohm*cm^2)."""
    for name, val in (("x_cldn2", x_cldn2), ("x_cldn4", x_cldn4), ("i_cytokine", i_cytokine)):
        if isinstance(val, (int, float)) and val < 0.0:
            raise ValueError(f"{name} must be non-negative, got {val}")

    p_na = (DEFAULT_P_NA_BASE + DEFAULT_ALPHA_CLDN2 * x_cldn2) / (
        1.0 + DEFAULT_BETA_CLDN4 * x_cldn4
    ) + DEFAULT_GAMMA_NA_CYT * i_cytokine
    p_hco3 = DEFAULT_P_HCO3_BASE / (
        1.0 + 0.28 * x_cldn4 + 0.30 * x_cldn2
    ) + DEFAULT_GAMMA_HCO3_CYT * i_cytokine
    p_cl = DEFAULT_P_CL_BASE / (1.0 + 0.32 * x_cldn4) + DEFAULT_GAMMA_CL_CYT * i_cytokine

    g_tj_ms = NERNST_PLANCK_G_SCALE_MS * (p_na * 122.5 + p_hco3 * 18.0 + p_cl * 115.0)
    r_tj_ohm_cm2 = 1000.0 / g_tj_ms
    return p_na, p_hco3, r_tj_ohm_cm2


def compute_ghk_paracellular_fluxes_and_asl_state(
    p_na: float | complex,
    p_hco3: float | complex,
    r_tj_ohm_cm2: float | complex,
    j_hco3_transcellular: float | complex,
    v_te_mv: float | complex,
) -> Tuple[float | complex, float | complex, float | complex, float | complex, float | complex]:
    """Computes steady-state (C_HCO3_asl, J_HCO3_para, J_Na_para, pH_ASL, h_ASL)."""
    for name, val in (
        ("p_na", p_na),
        ("p_hco3", p_hco3),
        ("r_tj_ohm_cm2", r_tj_ohm_cm2),
        ("j_hco3_transcellular", j_hco3_transcellular),
    ):
        if isinstance(val, (int, float)) and val < 0.0:
            raise ValueError(f"{name} must be non-negative, got {val}")

    is_complex = any(isinstance(v, complex) for v in (p_na, p_hco3, r_tj_ohm_cm2, j_hco3_transcellular, v_te_mv))
    exp_fn = cmath.exp if is_complex else math.exp
    log10_fn = cmath.log10 if is_complex else math.log10

    u_hco3 = (-1.0) * v_te_mv / THERMAL_VOLTAGE_MV_310K
    u_na = (+1.0) * v_te_mv / THERMAL_VOLTAGE_MV_310K
    b_hco3 = bernoulli_ghk_kernel(u_hco3)
    b_na = bernoulli_ghk_kernel(u_na)

    raw_num = C_HCO3_BASO_MM + (j_hco3_transcellular - DEFAULT_J_PROTON_ATP12A) / (3.6 * p_hco3 * b_hco3)
    c_hco3_asl = 0.5 * (raw_num + (raw_num ** 2 + 1.0) ** 0.5) / exp_fn(u_hco3)

    j_hco3_para = 3.6 * p_hco3 * b_hco3 * (C_HCO3_BASO_MM - c_hco3_asl * exp_fn(u_hco3))
    j_na_para = 3.6 * p_na * b_na * (C_NA_BASO_MM - C_NA_APICAL_MM * exp_fn(u_na))

    ph_asl = PKA_CARBONIC_37C + log10_fn(c_hco3_asl / C_CO2_DISSOLVED_MM)
    shunt_retention = r_tj_ohm_cm2 / (r_tj_ohm_cm2 + 220.0)
    h_asl_um = (
        2.4
        + 6.5
        * shunt_retention
        * (j_hco3_transcellular / (j_hco3_transcellular + 18.0))
        / (1.0 + (j_na_para / 650.0) ** 2)
    )
    return c_hco3_asl, j_hco3_para, j_na_para, ph_asl, h_asl_um


def complex_step_ph_transcellular_derivative(
    x_cldn2: float,
    x_cldn4: float,
    i_cytokine: float,
    j_hco3_transcellular: float,
    v_te_mv: float,
    step_h: float = 1e-30,
) -> Tuple[float, float, float]:
    """Verifies d(pH_ASL)/d(J_hco3_transcellular) at h = 1e-30 against central finite differences."""
    p_na, p_hco3, r_tj = compute_claudin_permeabilities_and_resistance(x_cldn2, x_cldn4, i_cytokine)

    def _eval_ph(j_val: float | complex) -> float | complex:
        _, _, _, ph_val, _ = compute_ghk_paracellular_fluxes_and_asl_state(p_na, p_hco3, r_tj, j_val, v_te_mv)
        return ph_val

    cs_val = _eval_ph(complex(j_hco3_transcellular, step_h))
    cs_deriv = cs_val.imag / step_h

    fd_h = 1e-5
    fd_deriv = float(_eval_ph(j_hco3_transcellular + fd_h) - _eval_ph(j_hco3_transcellular - fd_h)) / (2.0 * fd_h)
    rel_err = abs(cs_deriv - fd_deriv) / max(abs(cs_deriv), 1e-12)
    return cs_deriv, fd_deriv, rel_err


def evaluate_round200_tight_junction_claudin_flux(regime_key: str) -> Round200TightJunctionResult:
    """Evaluates the Round 200 Claudin-2/Claudin-4 paracellular electrodiffusion and ASL pH bridge."""
    if regime_key not in REGIME_PRESETS:
        raise KeyError(f"Unknown regime '{regime_key}'. Valid regimes: {list(REGIME_PRESETS.keys())}")

    cfg = REGIME_PRESETS[regime_key]
    x_cldn2 = cfg["x_cldn2"]
    x_cldn4 = cfg["x_cldn4"]
    i_cyt = cfg["i_cytokine"]
    j_trans = cfg["j_hco3_transcellular"]
    v_te = cfg["v_te_mv"]

    p_na_r, p_hco3_r, r_tj_r = compute_claudin_permeabilities_and_resistance(x_cldn2, x_cldn4, i_cyt)
    p_na = float(p_na_r)
    p_hco3 = float(p_hco3_r)
    r_tj = float(r_tj_r)
    selectivity = p_na / p_hco3

    chco3_r, jhpara_r, jnapara_r, ph_r, hasl_r = compute_ghk_paracellular_fluxes_and_asl_state(
        p_na, p_hco3, r_tj, j_trans, v_te
    )
    c_hco3 = float(chco3_r)
    j_hco3_para = float(jhpara_r)
    j_na_para = float(jnapara_r)
    ph_asl = float(ph_r)
    h_asl = float(hasl_r)

    _, _, rel_err = complex_step_ph_transcellular_derivative(x_cldn2, x_cldn4, i_cyt, max(j_trans, 5.0), v_te)

    inv_satisfied = (
        p_na > 0.0
        and p_hco3 > 0.0
        and selectivity > 0.0
        and r_tj > 0.0
        and c_hco3 > 0.0
        and 5.5 <= ph_asl <= 8.5
        and 2.0 <= h_asl <= 10.0
        and rel_err <= 1e-6
    )

    crosses_threshold = (
        r_tj >= 500.0
        and ph_asl >= 7.25
        and c_hco3 >= 16.0
        and abs(j_na_para) <= 450.0
        and h_asl >= 5.20
    )

    return Round200TightJunctionResult(
        regime_key=regime_key,
        x_cldn2=x_cldn2,
        x_cldn4=x_cldn4,
        i_cytokine=i_cyt,
        j_hco3_transcellular=j_trans,
        v_te_mv=v_te,
        p_na_1e6_cm_s=p_na,
        p_hco3_1e6_cm_s=p_hco3,
        cation_anion_selectivity=selectivity,
        r_tj_ohm_cm2=r_tj,
        c_hco3_asl_mm=c_hco3,
        j_hco3_para_nmol_cm2_h=j_hco3_para,
        j_na_para_nmol_cm2_h=j_na_para,
        ph_asl=ph_asl,
        h_asl_um=h_asl,
        crosses_kahn_kalai_threshold=crosses_threshold,
        falsification_invariant_satisfied=inv_satisfied,
        complex_step_relative_error=rel_err,
    )


def evaluate_falsification_negative_control(
    poison_negative_claudin: bool = False,
    poison_negative_cytokine: bool = False,
    poison_negative_transcellular_flux: bool = False,
) -> Dict[str, Any]:
    """Tests poisoned negative controls that must deterministically raise ValueError."""
    reasons: list[str] = []
    if poison_negative_claudin:
        try:
            compute_claudin_permeabilities_and_resistance(-0.5, 1.0, 0.0)
            return {"passed": False, "reason": "Negative Claudin-2 expression accepted unexpectedly"}
        except ValueError:
            reasons.append("Deterministically rejected negative Claudin-2 expression")

    if poison_negative_cytokine:
        try:
            compute_claudin_permeabilities_and_resistance(0.2, 1.0, -1.0)
            return {"passed": False, "reason": "Negative cytokine index accepted unexpectedly"}
        except ValueError:
            reasons.append("Deterministically rejected negative cytokine index")

    if poison_negative_transcellular_flux:
        try:
            compute_ghk_paracellular_fluxes_and_asl_state(1.0, 1.0, 600.0, -10.0, -8.0)
            return {"passed": False, "reason": "Negative transcellular HCO3 flux accepted unexpectedly"}
        except ValueError:
            reasons.append("Deterministically rejected negative transcellular flux")

    return {"passed": True, "reason": "; ".join(reasons) if reasons else "Nominal baseline"}
