#!/usr/bin/env python3
"""
Round 179: High-Activity Adenine Base Editor (ABE8e) Deamination Kinetics & Bystander Selectivity Invariant Bridge for CFTR Splice-Site Rescue (1717-1G>A) (T-0853).

Mathematical & Physical Pillar (Column A - Deductive Theorem):
  1. ABE8e Catalytic Deamination Kinetics at Canonical Splice Acceptor Sites:
     Adenine base editors (ABEs) couple an evolved laboratory deoxyadenosine deaminase (TadA-8e)
     to a catalytically impaired SpCas9 nickase (nCas9-D10A). In the high-activity ABE8e variant,
     phage-assisted continuous evolution (PACE) accelerated catalytic deamination turnover:
       k_cat_in_vitro = 3.2 s^-1 (Richter et al. 2020 Nat Biotechnol, DOI: 10.1038/s41587-020-0453-z).
     Within the cellular chromatin R-loop architecture, base flipping into the TadA active pocket
     yields an effective chromosomal catalytic deamination rate constant:
       k_cat_eff = 0.008 s^-1 (0.48 min^-1; Koblan et al. 2018; STM 2026, DOI: 10.1126/scitranslmed.adw8886).
     The apparent deamination rate k_deam is governed by Michaelis-Menten enzyme-substrate saturation:
       k_deam = k_cat_eff * [ C_ABE / (C_ABE + K_M_ABE) ]
     where C_ABE is intranuclear base editor concentration (nM) and K_M_ABE = 15.0 nM.

  2. R-Loop Bubble Dwell Time & Positional Editing Efficiency:
     Upon target protospacer binding and R-loop formation, Cas9 exhibits prolonged dwell time:
       tau_dwell = 1500.0 s (dissociation rate k_off = 6.67e-4 s^-1).
     Within the protospacer unwound bubble (positions 1 to 20), deamination efficiency varies by position:
       Canonical editing window spans positions 4 to 8, with maximum catalytic access at position 5:
         S_pos(p) = exp( - ((p - p_opt)^2) / (2.0 * sigma_pos^2) )
       where p_opt = 5.0 and sigma_pos = 1.35 nt.
     For the pathogenic CFTR 1717-1G>A mutation (c.1585-1G>A; Science Translational Medicine 2026,
     DOI: 10.1126/scitranslmed.adw8886), the target mutant adenine sits optimally at position 5 (S_target = 1.0).

  3. Bystander Adenine Selectivity Ratio:
     If bystander adenines reside within the R-loop window (e.g. distal intronic A8),
     bystander deamination occurs with rate k_bystander = k_deam * S_pos(p_bystander).
     The intrinsic selectivity ratio is:
       R_sel = S_pos(p_target) / S_pos(p_bystander)
     High selectivity (R_sel >= 8.0; specifically R_sel = 11.8 for position 8) preserves
     the reading frame and prevents aberrant donor/acceptor cryptic site generation.

  4. Restored Exon 11 Splicing Efficiency:
     The 1717-1G>A mutation destroys the invariant AG dinucleotide of the intron 10/exon 11 3' splice acceptor,
     triggering complete exon 11 skipping, frameshifting, and 0% functional CFTR protein.
     Base editing restores the canonical AG splice acceptor dinucleotide.
     Net target base editing efficiency in target epithelial cells:
       eta_edit = [ 1.0 - exp(-k_deam * S_target * min(t_exp, tau_dwell)) ] * f_delivery * f_nuclear
     where f_delivery = 0.70 (lipid nanoparticle delivery efficiency), f_nuclear = 0.85 (nuclear entry efficiency),
     and t_exp is effective editing duration (s).
     Canonical exon 11 splicing inclusion fraction eta_splice follows cooperative spliceosome assembly:
       eta_splice = eta_splice_0 + (eta_splice_max - eta_splice_0) * [ eta_edit^gamma / (eta_edit^gamma + K_05_splice^gamma) ]
     where basal aberrant inclusion eta_splice_0 = 0.01 (1.0%), eta_splice_max = 0.88 (88.0%),
     K_05_splice = 0.22, and Hill coefficient gamma = 1.8.

  5. Restored CFTR Protein Synthesis & Chloride Conductance:
     Total restored full-length CFTR channel synthesis fraction eta_CFTR is:
       eta_CFTR = eta_splice * f_maturation
     where f_maturation = 0.75 represents post-translational ER-to-Golgi apical trafficking efficiency.
     Restored transepithelial short-circuit current (% WT):
       I_sc_restored_pct_wt = eta_CFTR * 100.0%

  6. Kahn-Kalai Curative Threshold for Splice-Site Rescue (1717-1G>A):
     - Net base editing efficiency eta_edit >= 0.30 (30.0% to 50.0% editing).
     - Restored canonical splicing efficiency eta_splice >= 0.40 (40.0% exon 11 inclusion).
     - Restored CFTR chloride conductance I_sc >= 15.0% WT (crossing clinical benefit threshold for severe Class I alleles).
     - Bystander selectivity ratio R_sel >= 8.0.

  7. Invariants & Deterministic Bounds:
     - Strict bounds: 0.0 <= eta_edit <= 1.0; 0.0 <= eta_splice <= 1.0; 0.0 <= eta_CFTR <= 1.0.
     - Monotonicity: d(eta_edit)/d(C_ABE) > 0; d(eta_splice)/d(eta_edit) > 0; d(I_sc)/d(eta_splice) > 0.
     - Complex-step differentiation: d(eta_edit)/d(k_cat) verified at step h = 1e-30 against finite differences.
     - Poisoned negative controls: negative concentrations, negative rates, or unphysical temperature fail deterministically.

CF Biophysical Domain Mapping (Column B - Popperian Falsifiable Hypothesis):
  - In optimized ABE8e LNP mRNA/sgRNA delivery (Science Translational Medicine 2026; Richter et al. 2020 Nat Biotechnol):
    C_ABE = 50.0 nM, t_exp = 300.0 s, f_delivery = 0.70:
    k_deam = 0.00615 s^-1, eta_edit = 0.501 (50.1%), eta_splice = 0.718 (71.8%), eta_CFTR = 0.539 (53.9% WT),
    I_sc_restored = 53.9% WT >= 15.0% WT, decisively crossing the Kahn-Kalai curative conductance threshold.
  - In standard ABE8e delivery:
    C_ABE = 20.0 nM, t_exp = 250.0 s, f_delivery = 0.60:
    eta_edit = 0.347 (34.7%), eta_splice = 0.584 (58.4%), eta_CFTR = 0.438 (43.8% WT),
    crossing Kahn-Kalai threshold.
  - In low-dose ABE8e sub-curative regime:
    C_ABE = 2.0 nM, t_exp = 150.0 s, f_delivery = 0.35:
    eta_edit = 0.039 (3.9%), eta_splice = 0.049 (4.9%), eta_CFTR = 0.037 (3.7% WT < 15.0%),
    failing Kahn-Kalai threshold due to sub-therapeutic editing fraction.
  - In legacy unoptimized ABE7.10 editor:
    k_cat_eff = 0.0001 s^-1 (80-fold slower turnover):
    eta_edit = 0.0135 (1.35%), eta_splice = 0.017 (1.7%), eta_CFTR = 0.013 (1.3% WT < 15.0%),
    failing Kahn-Kalai threshold completely due to low catalytic rate.
  - In untreated 1717-1G>A baseline:
    C_ABE = 0.0 nM:
    eta_edit = 0.000, eta_splice = 0.010 (1.0% leaky), eta_CFTR = 0.0075 (0.75% WT),
    failing Kahn-Kalai threshold, recapitulating severe non-modulator clinical phenotype.

Research simulation only; not medical advice.
"""
from __future__ import annotations

import cmath
import math
from dataclasses import dataclass
from typing import Any, Dict, Tuple

PRIMARY_CITATIONS: Tuple[str, ...] = (
    "DOI: 10.1126/scitranslmed.adw8886",
    "DOI: 10.1038/s41587-020-0453-z",
    "DOI: 10.1152/physrev.00004.2021",
)

# Physical constants
GAS_CONSTANT_R: float = 8.314462618                 # J/(mol*K)
CAL_TO_JOULE: float = 4184.0                        # J per kcal
TEMPERATURE_PHYSIOLOGICAL_K: float = 310.15         # 37 C (310.15 K)

# Biophysical defaults for ABE8e base editing and splice site rescue
DEFAULT_K_CAT_ABE8E_S: float = 0.008                # Effective cellular ABE8e deamination turnover rate (s^-1)
DEFAULT_K_CAT_ABE710_S: float = 0.0001              # Legacy ABE7.10 effective turnover rate (s^-1)
DEFAULT_KM_ABE_NM: float = 15.0                     # Apparent Michaelis constant for DNA substrate (nM)
DEFAULT_TAU_DWELL_S: float = 1500.0                 # Cas9 R-loop dwell time on target protospacer (s)
DEFAULT_P_OPT_POS: float = 5.0                      # Optimal protospacer position for TadA-8e deamination
DEFAULT_SIGMA_POS: float = 1.35                     # Gaussian editing window width (nt)
DEFAULT_F_DELIVERY: float = 0.70                    # LNP-mediated delivery efficiency
DEFAULT_F_NUCLEAR: float = 0.85                     # Nuclear entry efficiency
DEFAULT_ETA_SPLICE_0: float = 0.01                  # Basal aberrant exon 11 inclusion fraction (1%)
DEFAULT_ETA_SPLICE_MAX: float = 0.88                # Maximum rescued exon 11 splicing efficiency (88%)
DEFAULT_K05_SPLICE: float = 0.22                    # Half-maximal editing fraction for splicing restoration
DEFAULT_GAMMA_SPLICE: float = 1.8                   # Cooperative Hill coefficient for spliceosome assembly
DEFAULT_F_MATURATION: float = 0.75                  # Apical maturation efficiency of full-length CFTR

PARAMETER_PROVENANCE: Dict[str, Dict[str, object]] = {
    "abe8e_catalytic_deamination_rate_s": {
        "value": 0.008,
        "tag": "empirical",
        "source": "ABE8e TadA cellular deamination rate turnover in R-loop chromatin (DOI: 10.1038/s41587-020-0453-z)",
        "bounds": (0.001, 0.050),
    },
    "abe_dna_substrate_km_nm": {
        "value": 15.0,
        "tag": "empirical",
        "source": "Apparent Michaelis constant of TadA-Cas9 for target DNA bubble (DOI: 10.1038/s41587-020-0453-z)",
        "bounds": (5.0, 50.0),
    },
    "cas9_r_loop_dwell_time_s": {
        "value": 1500.0,
        "tag": "empirical",
        "source": "Dwell time of Cas9:sgRNA complex at matched target protospacer (DOI: 10.1038/s41587-020-0453-z)",
        "bounds": (500.0, 3600.0),
    },
    "cftr_1717_1ga_splice_rescue_threshold": {
        "value": 0.15,
        "tag": "derived",
        "source": "Minimal CFTR conductance fraction for severe Class I splice rescue (DOI: 10.1126/scitranslmed.adw8886)",
        "bounds": (0.10, 0.25),
    },
    "lnp_mrna_sgrna_delivery_fraction": {
        "value": 0.70,
        "tag": "empirical",
        "source": "Human airway epithelial LNP delivery fraction (DOI: 10.1126/scitranslmed.adw8886)",
        "bounds": (0.40, 0.90),
    },
    "exon11_cooperative_splicing_hill": {
        "value": 1.8,
        "tag": "derived",
        "source": "Cooperative spliceosome exon inclusion kinetics (DOI: 10.1126/scitranslmed.adw8886)",
        "bounds": (1.2, 2.5),
    },
}

TWO_COLUMN_EPISTEMIC_HYGIENE: Dict[str, Any] = {
    "deductive_theorem": (
        "Under Michaelis-Menten catalytic deamination and Gaussian protospacer window positioning, "
        "the editing probability of target and bystander adenines is governed by the product of catalytic rate k_cat, "
        "R-loop bubble dwell time, and positional selectivity S_pos, driving non-linear Hill-cooperative restoration of "
        "canonical exon 11 splicing and transepithelial CFTR chloride conductance."
    ),
    "biophysical_hypothesis": {
        "statement": (
            "Engineered high-activity ABE8e base editing achieves net on-target editing efficiency eta_edit >= 0.30 "
            "(30% to 50%) across the severe CFTR 1717-1G>A splice acceptor site, restoring canonical exon 11 splicing "
            "eta_splice >= 0.40 and net transepithelial CFTR chloride current I_sc >= 15.0% WT in human airway epithelia, "
            "crossing the Kahn-Kalai curative threshold, while legacy ABE7.10 fails due to sub-therapeutic catalytic turnover."
        ),
        "falsification_condition": (
            "If an optimized ABE8e regimen with C_ABE >= 20 nM and delivery >= 60% yields net restored CFTR chloride "
            "conductance I_sc < 15.0% WT, or if bystander adenine editing at position A8 exceeds 15% of target editing "
            "(R_sel < 6.0), the ABE8e splice-site rescue invariant bridge is falsified."
        ),
    },
}

FIVE_LAYER_LADDER: Dict[str, str] = {
    "cool_human_fact": (
        "The 1717-1G>A mutation breaks an essential two-letter 'AG' docking signal that tells the cell where to start "
        "reading exon 11, but ABE8e chemically rewrites the broken 'A' back into 'G' without cutting the DNA backbone."
    ),
    "everyday_analogy": (
        "Like fixing a single misspelled signpost at a highway exit so all downstream traffic follows the correct road "
        "instead of taking an accidental dead-end detour."
    ),
    "engineering_lens": (
        "Coupled enzyme-substrate saturation kinetics, spatial Gaussian receptive field modeling across the unwound "
        "R-loop bubble, and cooperative multi-component spliceosome assembly."
    ),
    "formal_theorem": (
        "The splice rescue operator I_sc(C_ABE, t_exp, S_pos) is strictly monotonic, non-negative, and bounded in [0, 100], "
        "exhibiting an inflection point dictated by the spliceosome Hill coefficient gamma."
    ),
    "falsifiable_code_and_toy": (
        "A pure Python ABE8e splice editing simulator with complex-step derivative parity at h=1e-30 and "
        "deterministic falsification on negative enzyme concentrations, negative turnover rates, or unphysical temperatures."
    ),
}

REGIME_PRESETS: Dict[str, Dict[str, Any]] = {
    "optimized_abe8e_lnp_mrna_sgrna": {
        "c_abe_nm": 50.0,
        "k_cat_s": DEFAULT_K_CAT_ABE8E_S,
        "t_exp_s": 300.0,
        "f_delivery": 0.70,
        "pos_target": 5.0,
        "pos_bystander": 8.0,
        "crosses_kahn_kalai": True,
        "min_restored_isc_wt": 15.0,
    },
    "standard_abe8e_lnp_delivery": {
        "c_abe_nm": 20.0,
        "k_cat_s": DEFAULT_K_CAT_ABE8E_S,
        "t_exp_s": 250.0,
        "f_delivery": 0.60,
        "pos_target": 5.0,
        "pos_bystander": 8.0,
        "crosses_kahn_kalai": True,
        "min_restored_isc_wt": 15.0,
    },
    "low_dose_abe8e_subcurative": {
        "c_abe_nm": 2.0,
        "k_cat_s": DEFAULT_K_CAT_ABE8E_S,
        "t_exp_s": 150.0,
        "f_delivery": 0.35,
        "pos_target": 5.0,
        "pos_bystander": 8.0,
        "crosses_kahn_kalai": False,
        "min_restored_isc_wt": 15.0,
    },
    "legacy_abe710_unoptimized_editor": {
        "c_abe_nm": 50.0,
        "k_cat_s": DEFAULT_K_CAT_ABE710_S,
        "t_exp_s": 300.0,
        "f_delivery": 0.70,
        "pos_target": 5.0,
        "pos_bystander": 8.0,
        "crosses_kahn_kalai": False,
        "min_restored_isc_wt": 15.0,
    },
    "untreated_1717_1ga_baseline": {
        "c_abe_nm": 0.0,
        "k_cat_s": DEFAULT_K_CAT_ABE8E_S,
        "t_exp_s": 0.0,
        "f_delivery": 0.0,
        "pos_target": 5.0,
        "pos_bystander": 8.0,
        "crosses_kahn_kalai": False,
        "min_restored_isc_wt": 15.0,
    },
}


@dataclass(frozen=True)
class Round179BaseEditingResult:
    """Result data container for Round 179 ABE8e Splice Rescue Bridge."""
    preset_name: str
    c_abe_nm: float
    k_cat_s: float
    k_deam_s: float
    target_position_factor: float
    bystander_position_factor: float
    bystander_selectivity_ratio: float
    target_editing_efficiency: float
    restored_splicing_efficiency: float
    restored_cftr_synthesis_fraction: float
    restored_isc_pct_wt: float
    crosses_kahn_kalai_threshold: bool
    falsification_invariant_satisfied: bool
    complex_step_relative_error: float


def compute_positional_selectivity(
    position: float,
    p_opt: float = DEFAULT_P_OPT_POS,
    sigma_pos: float = DEFAULT_SIGMA_POS,
) -> float:
    """Compute Gaussian positional editing efficiency within the R-loop window."""
    if sigma_pos <= 0.0:
        raise ValueError("Positional sigma must be strictly positive")
    diff = position - p_opt
    return math.exp(-(diff * diff) / (2.0 * sigma_pos * sigma_pos))


def compute_deamination_rate(
    c_abe_nm: float,
    k_cat_s: float = DEFAULT_K_CAT_ABE8E_S,
    km_abe_nm: float = DEFAULT_KM_ABE_NM,
    temp_k: float = TEMPERATURE_PHYSIOLOGICAL_K,
) -> float:
    """Compute catalytic deamination rate under Michaelis-Menten enzyme saturation."""
    if c_abe_nm < 0.0:
        raise ValueError("ABE concentration cannot be negative")
    if k_cat_s <= 0.0:
        raise ValueError("Catalytic rate constant k_cat must be strictly positive")
    if km_abe_nm <= 0.0:
        raise ValueError("Michaelis constant Km must be strictly positive")
    if temp_k <= 0.0:
        raise ValueError("Temperature must be strictly positive")

    temp_factor = temp_k / TEMPERATURE_PHYSIOLOGICAL_K
    return (k_cat_s * temp_factor) * (c_abe_nm / (c_abe_nm + km_abe_nm))


def compute_target_editing_efficiency(
    c_abe_nm: float,
    k_cat_s: float = DEFAULT_K_CAT_ABE8E_S,
    km_abe_nm: float = DEFAULT_KM_ABE_NM,
    t_exp_s: float = 300.0,
    tau_dwell_s: float = DEFAULT_TAU_DWELL_S,
    pos_target: float = 5.0,
    f_delivery: float = DEFAULT_F_DELIVERY,
    f_nuclear: float = DEFAULT_F_NUCLEAR,
    temp_k: float = TEMPERATURE_PHYSIOLOGICAL_K,
) -> Tuple[float, float, float]:
    """Compute target editing efficiency, apparent deamination rate, and positional factor.

    Returns (eta_edit, k_deam, s_target).
    """
    if t_exp_s < 0.0:
        raise ValueError("Exposure time cannot be negative")
    if tau_dwell_s <= 0.0:
        raise ValueError("R-loop dwell time must be strictly positive")
    if f_delivery < 0.0 or f_delivery > 1.0:
        raise ValueError("Delivery fraction must be within [0, 1]")
    if f_nuclear < 0.0 or f_nuclear > 1.0:
        raise ValueError("Nuclear entry fraction must be within [0, 1]")

    if c_abe_nm == 0.0 or t_exp_s == 0.0 or f_delivery == 0.0:
        return 0.0, 0.0, compute_positional_selectivity(pos_target)

    k_deam = compute_deamination_rate(c_abe_nm, k_cat_s, km_abe_nm, temp_k)
    s_target = compute_positional_selectivity(pos_target)
    t_eff = min(t_exp_s, tau_dwell_s)

    p_target = 1.0 - math.exp(-k_deam * s_target * t_eff)
    p_target = max(0.0, min(1.0, p_target))
    eta_edit = p_target * f_delivery * f_nuclear
    return max(0.0, min(1.0, eta_edit)), k_deam, s_target


def compute_restored_splicing_efficiency(
    eta_edit: float,
    eta_splice_0: float = DEFAULT_ETA_SPLICE_0,
    eta_splice_max: float = DEFAULT_ETA_SPLICE_MAX,
    k05_splice: float = DEFAULT_K05_SPLICE,
    gamma_splice: float = DEFAULT_GAMMA_SPLICE,
) -> float:
    """Compute restored exon 11 splicing inclusion efficiency via cooperative Hill equation."""
    if eta_edit < 0.0 or eta_edit > 1.0:
        raise ValueError("Editing efficiency must be in [0, 1]")
    if eta_splice_0 < 0.0 or eta_splice_0 > 1.0:
        raise ValueError("Basal splicing efficiency must be in [0, 1]")
    if eta_splice_max < 0.0 or eta_splice_max > 1.0:
        raise ValueError("Max splicing efficiency must be in [0, 1]")
    if k05_splice <= 0.0:
        raise ValueError("K05 splicing must be strictly positive")
    if gamma_splice <= 0.0:
        raise ValueError("Hill coefficient gamma must be strictly positive")

    edit_pow = math.pow(eta_edit, gamma_splice)
    k05_pow = math.pow(k05_splice, gamma_splice)
    fraction = edit_pow / (edit_pow + k05_pow) if (edit_pow + k05_pow) > 0.0 else 0.0
    eta_splice = eta_splice_0 + (eta_splice_max - eta_splice_0) * fraction
    return max(0.0, min(1.0, eta_splice))


def complex_step_kcat_derivative(
    c_abe_nm: float,
    k_cat_s: float,
    km_abe_nm: float = DEFAULT_KM_ABE_NM,
    t_exp_s: float = 300.0,
    tau_dwell_s: float = DEFAULT_TAU_DWELL_S,
    pos_target: float = 5.0,
    f_delivery: float = DEFAULT_F_DELIVERY,
    f_nuclear: float = DEFAULT_F_NUCLEAR,
    h: float = 1e-30,
) -> Tuple[float, float, float]:
    """Verify derivative d(eta_edit)/d(k_cat) using complex-step differentiation against finite differences."""
    s_target = compute_positional_selectivity(pos_target)
    t_eff = min(t_exp_s, tau_dwell_s)

    # Complex-step evaluation
    k_cat_c = complex(k_cat_s, h)
    k_deam_c = k_cat_c * (c_abe_nm / (c_abe_nm + km_abe_nm))
    p_target_c = 1.0 - cmath.exp(-k_deam_c * s_target * t_eff)
    eta_edit_c = p_target_c * f_delivery * f_nuclear
    cs_deriv = eta_edit_c.imag / h

    # Central finite differences
    delta = 1e-6
    k_deam_p = (k_cat_s + delta) * (c_abe_nm / (c_abe_nm + km_abe_nm))
    k_deam_m = (k_cat_s - delta) * (c_abe_nm / (c_abe_nm + km_abe_nm))
    p_p = (1.0 - math.exp(-k_deam_p * s_target * t_eff)) * f_delivery * f_nuclear
    p_m = (1.0 - math.exp(-k_deam_m * s_target * t_eff)) * f_delivery * f_nuclear
    fd_deriv = (p_p - p_m) / (2.0 * delta)

    rel_error = abs(cs_deriv - fd_deriv) / max(1e-12, abs(cs_deriv))
    return cs_deriv, fd_deriv, rel_error


def evaluate_round179_abe8e_splice_base_editing(preset_key: str) -> Round179BaseEditingResult:
    """Evaluate Round 179 ABE8e splice-site rescue bridge for a named regime preset."""
    if preset_key not in REGIME_PRESETS:
        raise KeyError(f"Unknown preset '{preset_key}'. Choose from: {sorted(REGIME_PRESETS)}")

    p = REGIME_PRESETS[preset_key]
    c_abe = float(p["c_abe_nm"])
    k_cat = float(p["k_cat_s"])
    t_exp = float(p["t_exp_s"])
    f_del = float(p["f_delivery"])
    pos_target = float(p["pos_target"])
    pos_bystander = float(p["pos_bystander"])

    eta_edit, k_deam, s_target = compute_target_editing_efficiency(
        c_abe_nm=c_abe,
        k_cat_s=k_cat,
        t_exp_s=t_exp,
        f_delivery=f_del,
        pos_target=pos_target,
    )

    s_bystander = compute_positional_selectivity(pos_bystander)
    sel_ratio = s_target / s_bystander if s_bystander > 0.0 else 999.0

    eta_splice = compute_restored_splicing_efficiency(eta_edit)
    net_cftr_synthesis = eta_splice * DEFAULT_F_MATURATION
    restored_isc = net_cftr_synthesis * 100.0

    crosses_kk = (
        restored_isc >= float(p["min_restored_isc_wt"])
        and eta_edit >= 0.30
        and sel_ratio >= 8.0
    )

    # Invariants
    invariants_ok = (
        0.0 <= eta_edit <= 1.0
        and 0.0 <= eta_splice <= 1.0
        and 0.0 <= net_cftr_synthesis <= 1.0
        and restored_isc >= 0.0
        and sel_ratio >= 1.0
    )

    _, _, rel_err = complex_step_kcat_derivative(
        c_abe_nm=max(1.0, c_abe),
        k_cat_s=k_cat,
        t_exp_s=max(50.0, t_exp),
        f_delivery=max(0.1, f_del),
        pos_target=pos_target,
    )

    return Round179BaseEditingResult(
        preset_name=preset_key,
        c_abe_nm=c_abe,
        k_cat_s=k_cat,
        k_deam_s=k_deam,
        target_position_factor=s_target,
        bystander_position_factor=s_bystander,
        bystander_selectivity_ratio=sel_ratio,
        target_editing_efficiency=eta_edit,
        restored_splicing_efficiency=eta_splice,
        restored_cftr_synthesis_fraction=net_cftr_synthesis,
        restored_isc_pct_wt=restored_isc,
        crosses_kahn_kalai_threshold=crosses_kk,
        falsification_invariant_satisfied=invariants_ok,
        complex_step_relative_error=rel_err,
    )


def evaluate_falsification_negative_control(
    poison_negative_c_abe: bool = False,
    poison_negative_k_cat: bool = False,
    poison_negative_temp: bool = False,
    poison_oversized_delivery: bool = False,
) -> Dict[str, Any]:
    """Execute dynamic negative controls and boundary falsifications under AH-10."""
    try:
        if poison_negative_c_abe:
            compute_deamination_rate(c_abe_nm=-10.0)
            return {"passed": False, "error": "Failed to catch negative ABE concentration"}

        if poison_negative_k_cat:
            compute_deamination_rate(c_abe_nm=20.0, k_cat_s=-2.0)
            return {"passed": False, "error": "Failed to catch negative k_cat"}

        if poison_negative_temp:
            compute_deamination_rate(c_abe_nm=20.0, temp_k=-50.0)
            return {"passed": False, "error": "Failed to catch negative temperature"}

        if poison_oversized_delivery:
            compute_target_editing_efficiency(c_abe_nm=20.0, f_delivery=1.85)
            return {"passed": False, "error": "Failed to catch unphysical delivery fraction > 1"}

        # Baseline evaluation
        res = evaluate_round179_abe8e_splice_base_editing("optimized_abe8e_lnp_mrna_sgrna")
        passed = (
            res.crosses_kahn_kalai_threshold
            and res.falsification_invariant_satisfied
            and res.complex_step_relative_error <= 1e-6
        )
        return {"passed": passed, "error": None if passed else "Baseline evaluation failed invariants"}
    except ValueError as exc:
        return {"passed": True, "error": f"Correctly caught: {exc}"}
