#!/usr/bin/env python3
"""
Round 203: Site-Specific m1A Chemically Modified Suppressor tRNA & Microfluidic Ionizable LNP Nonsense Rescue Bridge (T-0896).

Mathematical & Physical Pillar (Column A - Deductive Theorem):
  1. Endosomal Nuclease Persistence & TLR7/8 Innate Cytokine Evasion via Site-Specific m1A Modification:
     Nonviral delivery of engineered suppressor tRNAs (sup-tRNAs) to rescue CFTR premature termination codons (PTCs;
     W1282X in exon 23 and G542X in exon 12, addressing part of the CFFPR 2024 8.0% genotype-ineligible cohort) is
     rate-limited by rapid intracellular/endosomal nuclease decay and Toll-like receptor 7/8 (TLR7/8) innate sensing
     (Chen, Cui, Li et al., Science 2026, DOI: 10.1126/science.aeb0054).
     Site-specific N1-methyladenosine (m1A) chemical modification in the T-Psi-C loop (normalized T-loop modification
     index f_mod in [0, 1]) stabilizes the tertiary L-shaped tRNA elbow and attenuates nuclease cleavage:
       k_rnase(f_mod) = k_rnase0 * ((1.0 - sigma_rnase) + sigma_rnase * exp(-alpha_stab * f_mod))
       t_half_h(f_mod) = ln(2) / k_rnase(f_mod)
     where k_rnase0 = 0.69315 h^-1 (unmodified half-life ~1.0 h), sigma_rnase = 0.90, and alpha_stab = 5.5.
     Encapsulation inside microfluidic-engineered ionizable lipid nanoparticles (LNPs; Hanna, Shepherd, Datto et al.,
     ACS Nano 2025, DOI: 10.1021/acsnano.5c15613) with encapsulation fraction eta_encap in [0, 1] and endosomal
     escape efficiency eta_esc in [0, 1] yields well-mixed first-order cytosolic intact sup-tRNA bioavailability B_trna:
       B_trna = eta_encap * eta_esc * (t_half_h / (t_half_h + ln(2) * tau_transit_h))
     where tau_transit_h = 2.5 h (since 1 / (1 + k_rnase * tau_transit_h) = t_half_h / (t_half_h + ln(2) * tau_transit_h)).
     Endosomal TLR7/8 innate cytokine activation index I_tlr78 in [0, 1] is suppressed by site-specific m1A chemical
     modification minimizing endosomal innate immune recognition:
       I_tlr78 = I_tlr0 * exp(-beta_tlr * f_mod) * (0.65 + 0.35 * (1.0 - eta_esc))
     where I_tlr0 = 0.92 and beta_tlr = 5.2.

  2. Aminoacyl-tRNA Synthetase Charging Fidelity, Evo 2 Structural Fitness, & PTC Decoding:
     Chemical modifications encroaching on the anticodon triplet (positions 34-36) or discriminator/CCA 3'-acceptor
     stem (positions 73-76), parameterized by core-blocking fraction f_core_block in [0, 1], sterically abolish
     aminoacylation and ribosomal A-site decoding (Chen et al., Science 2026, DOI: 10.1126/science.aeb0054), while
     L-shaped cloverleaf coaxial stacking fitness S_evo2 in [0, 1] (contextualized by Evo 2 pan-domain sequence
     likelihood modeling; Brixi, Durrant, Ku et al., Nature 2026, DOI: 10.1038/s41586-026-10176-5) bounds productive
     synthetase charging:
       eta_charge = S_evo2 * exp(-kappa_core * f_core_block) / (1.0 + (f_mod / K_overmod)^4)
     where kappa_core = 7.5 and K_overmod = 0.78.
     Productive ribosomal PTC readthrough fraction alpha_suptrna competing against eRF1/eRF3 release factors satisfies:
       alpha_suptrna = alpha_max * (B_trna * eta_charge) / (B_trna * eta_charge + K_erf1)
     where alpha_max = 0.38 and K_erf1 = 0.085.

  3. Cognate vs. Near-Cognate NBD2 Thermodynamic Folding Penalty & Restored Anion Efflux:
     For CFTR W1282X (NBD2 beta-strand 1) and G542X (NBD1), inserting cognate amino acid (Trp at 1282; Gly at 542)
     preserves native domain folding (ddg_insert_kcal = 0.0 kcal/mol), whereas near-cognate cavity-forming insertion
     at W1282X imposes a thermodynamic destabilization ddg_insert_kcal > 0 mitigated by allosteric corrector binding
     ddg_corr_kcal >= 0:
       P_fold = 1.0 / (1.0 + exp((ddg_insert_kcal - ddg_corr_kcal - dg0_margin_kcal) / RT_kcal))
     where dg0_margin_kcal = 1.10 kcal/mol and RT_kcal = 0.616 kcal/mol (310.15 K).
     Active ribosomal readthrough shields PTC mRNA from UPF1-mediated nonsense-mediated decay (NMD):
       f_mrna_eff = f_mrna0 + (1.0 - f_mrna0) * (alpha_suptrna / (alpha_suptrna + K_nmd_shield))
     where K_nmd_shield = 0.12. Functional apical CFTR density P_cftr = f_mrna_eff * alpha_suptrna * P_fold drives
     equivalent transepithelial anion efflux current (% WT):
       I_sc_equiv = 100.0 * (P_cftr / (P_cftr + K_05_cftr * (1.0 - P_cftr)))
     where K_05_cftr = 0.18.

  4. Kahn-Kalai Curative Threshold & Deterministic Invariants:
     - Endosomal sup-tRNA half-life t_half_h >= 4.0 h.
     - Endosomal TLR7/8 innate activation index I_tlr78 <= 0.20.
     - Synthetase charging & decoding competence eta_charge >= 0.65.
     - Restored equivalent anion efflux current I_sc_equiv >= 15.0% WT.
     - Complex-step derivative d(I_sc_equiv)/d(f_mod) at h = 1e-30 matches central finite differences (rel_err <= 1e-6).
     - Poisoned negative controls: out-of-range modification fractions, negative encapsulation, or negative energetic
       inputs raise ValueError deterministically; unmodified naked sup-tRNA and anticodon-blocked sup-tRNA fail
       the curative threshold deterministically (AH-10).

CF Biophysical Domain Mapping (Column B - Popperian Falsifiable Hypothesis):
  - W1282X Chemically Modified Cognate sup-tRNA + Microfluidic LNP Curative Regimen (f_mod=0.42, f_core_block=0.0):
    t_half = 5.28 h >= 4.0 h, I_tlr78 = 0.0898 <= 0.20, eta_charge = 0.8763 >= 0.65, I_sc_equiv = 59.98% WT >= 15.0% (Pass).
  - G542X Chemically Modified Cognate sup-tRNA + Microfluidic LNP Curative Regimen (f_mod=0.45, f_core_block=0.0):
    t_half = 5.69 h >= 4.0 h, I_tlr78 = 0.0762 <= 0.20, eta_charge = 0.8643 >= 0.65, I_sc_equiv = 56.09% WT >= 15.0% (Pass).
  - Unmodified Naked sup-tRNA RNase Degradation & TLR7/8 Spike Control (f_mod=0.0, eta_encap=0.08, eta_esc=0.05):
    t_half = 1.00 h < 4.0 h, I_tlr78 = 0.9039 > 0.20, I_sc_equiv = 0.75% WT < 15.0% (Fail).
  - Anticodon-Blocked Overmodified Synthetase Inactivation Control (f_mod=0.55, f_core_block=0.70):
    eta_charge = 0.0037 < 0.65, alpha_suptrna = 0.0047, I_sc_equiv = 0.85% WT < 15.0% (Fail).
  - W1282X Near-Cognate NBD2 Thermodynamic Misfolding Unrescued Control (ddg_insert=2.85 kcal/mol, ddg_corr=0.0):
    P_fold = 0.0552, P_cftr = 0.0117, I_sc_equiv = 6.19% WT < 15.0% (Fail).

Research simulation only; not medical advice.
"""
from __future__ import annotations

import cmath
import math
from dataclasses import dataclass
from typing import Any, Dict, Tuple

RESEARCH_BOUNDARY: str = "Research simulation only; not medical advice."

PRIMARY_CITATIONS: Tuple[str, ...] = (
    "DOI: 10.1126/science.aeb0054",
    "DOI: 10.1021/acsnano.5c15613",
    "DOI: 10.1038/s41586-026-10176-5",
)

# Biophysical Constants for Chemically Modified sup-tRNA & Microfluidic LNP Delivery
DEFAULT_K_RNASE0_PER_H: float = 0.6931471805599453  # Unmodified tRNA endosomal decay (~1.0 h half-life)
DEFAULT_SIGMA_RNASE_MAX: float = 0.90
DEFAULT_ALPHA_STAB_MOD: float = 5.5
DEFAULT_TAU_TRANSIT_H: float = 2.5
DEFAULT_I_TLR0: float = 0.92
DEFAULT_BETA_TLR_SUPPRESS: float = 5.2
DEFAULT_KAPPA_CORE_BLOCK: float = 7.5
DEFAULT_K_OVERMOD: float = 0.78
DEFAULT_ALPHA_MAX_READTHROUGH: float = 0.38
DEFAULT_K_ERF1_COMPETITION: float = 0.085
DEFAULT_DG0_FOLD_MARGIN_KCAL: float = 1.10
DEFAULT_RT_310K_KCAL: float = 0.616
DEFAULT_K_NMD_SHIELD: float = 0.12
DEFAULT_K05_CFTR_ANION_EFFLUX: float = 0.18

PARAMETER_PROVENANCE: Dict[str, Dict[str, object]] = {
    "rnase_stabilization_attenuation_slope": {
        "value": DEFAULT_ALPHA_STAB_MOD,
        "tag": "phenomenological_heuristic",
        "source": (
            "Analytical exponential slope alpha_stab relating T-loop site-specific m1A modification index to nuclease "
            "persistence (exact 5.5 dimensionless scalar not stated in Chen et al., Science 2026, "
            "DOI: 10.1126/science.aeb0054)"
        ),
        "bounds": (2.5, 9.0),
    },
    "tlr78_evasion_exponential_slope": {
        "value": DEFAULT_BETA_TLR_SUPPRESS,
        "tag": "phenomenological_heuristic",
        "source": (
            "Analytical innate immune suppression slope beta_tlr for site-specific m1A-modified sup-tRNA attenuation "
            "of TLR7/8 sensing (exact 5.2 dimensionless scalar not stated in Chen et al., Science 2026, "
            "DOI: 10.1126/science.aeb0054)"
        ),
        "bounds": (2.5, 8.5),
    },
    "anticodon_cca_core_inactivation_slope": {
        "value": DEFAULT_KAPPA_CORE_BLOCK,
        "tag": "phenomenological_heuristic",
        "source": (
            "Analytical synthetase/decoding steric inactivation penalty kappa_core when chemical modifications "
            "encroach on anticodon or 3'-CCA core positions (exact 7.5 scalar not stated in Chen et al., Science 2026, "
            "DOI: 10.1126/science.aeb0054; structural constraint context in Brixi et al., Nature 2026, "
            "DOI: 10.1038/s41586-026-10176-5)"
        ),
        "bounds": (4.0, 12.0),
    },
    "lnp_endosomal_transit_timescale_h": {
        "value": DEFAULT_TAU_TRANSIT_H,
        "tag": "phenomenological_heuristic",
        "source": (
            "Characteristic endosomal-to-cytosolic trafficking timescale tau_transit_h for microfluidic ionizable LNPs "
            "(exact 2.5 h scalar not stated in Hanna et al., ACS Nano 2025, DOI: 10.1021/acsnano.5c15613)"
        ),
        "bounds": (1.0, 6.0),
    },
    "erf1_ribosomal_competition_half_saturation": {
        "value": DEFAULT_K_ERF1_COMPETITION,
        "tag": "phenomenological_heuristic",
        "source": (
            "Analytical half-saturation constant K_erf1 for charged sup-tRNA competition against eRF1/eRF3 at CFTR "
            "premature stop codons (exact 0.085 dimensionless scalar not stated in Chen et al., Science 2026, "
            "DOI: 10.1126/science.aeb0054)"
        ),
        "bounds": (0.03, 0.20),
    },
    "cftr_anion_efflux_spare_capacity_k05": {
        "value": DEFAULT_K05_CFTR_ANION_EFFLUX,
        "tag": "phenomenological_heuristic",
        "source": (
            "Nonlinear epithelial spare-capacity constant K_05_cftr mapping rescued apical CFTR fraction to short-circuit "
            "anion efflux current % WT (exact 0.18 scalar not stated in Chen et al., Science 2026, "
            "DOI: 10.1126/science.aeb0054)"
        ),
        "bounds": (0.08, 0.30),
    },
}

TWO_COLUMN_EPISTEMIC_HYGIENE: Dict[str, Dict[str, str]] = {
    "deductive_column_a": {
        "formal_theorem": (
            "On the interior domain (f_mod, eta_encap, eta_esc) in (0, 1)^3 with f_core_block = 0, site-specific "
            "T-loop m1A modification simultaneously increases sup-tRNA persistence half-life t_half_h(f_mod) "
            "and decreases TLR7/8 activation I_tlr78(f_mod) while preserving synthetase charging up to K_overmod, "
            "yielding a real-analytic anion efflux rescue operator I_sc_equiv with exact complex-step sensitivity."
        ),
        "invariants": (
            "t_half_h >= ln(2)/k_rnase0; 0 <= B_trna, I_tlr78, eta_charge, alpha_suptrna, P_fold, P_cftr <= 1; "
            "0 <= I_sc_equiv <= 100%; d(t_half_h)/d(f_mod) > 0; d(I_tlr78)/d(f_mod) < 0; "
            "d(eta_charge)/d(f_core_block) < 0; complex-step relative error <= 1e-6."
        ),
    },
    "popperian_column_b": {
        "falsifiable_hypothesis": (
            "Cargo-tailored microfluidic ionizable LNP delivery of site-specifically m1A-modified cognate "
            "sup-tRNA rescues >= 15.0% WT transepithelial anion efflux in W1282X and G542X CF airway epithelia "
            "(addressing part of the CFFPR 2024 8.0% modulator-ineligible cohort) while keeping endosomal TLR7/8 "
            "activation <= 0.20, whereas unmodified naked tRNA or anticodon-blocked tRNA fails to cross threshold."
        ),
        "empirical_falsification_threshold": (
            "Falsified if heavy chemical modification of anticodon positions 34-36 or 3'-CCA positions 73-76 retains "
            ">= 65% aminoacylation efficiency without synthetase engineering, or if unmodified unencapsulated tRNA "
            "achieves t_half_h >= 4.0 h and I_tlr78 <= 0.20 in primary CF bronchial epithelial endosomes."
        ),
    },
}

FIVE_LAYER_LADDER: Dict[str, str] = {
    "cool_human_fact": (
        "In the CFFPR 2024 annual report, 8.0% of registered people with cystic fibrosis remain genotype-ineligible "
        "for small-molecule CFTR modulators, often because nonsense stop mutations cut the protein short before it finishes."
    ),
    "everyday_analogy": (
        "An m1A-modified suppressor tRNA inside a lipid nanoparticle acts like a courier wearing a rainproof "
        "stealth jacket who replaces a misprinted 'STOP' sign on the assembly line with the missing amino acid brick—so "
        "long as the modification never covers the courier's hands."
    ),
    "engineering_lens": (
        "Site-specific T-loop m1A modification outside the synthetase identity and anticodon elements enhances "
        "intracellular persistence, aminoacylation competence, and TLR7/8 evasion without steric inactivation, while "
        "cognate residue insertion avoids the W1282X NBD2 cavity folding penalty."
    ),
    "formal_theorem": (
        "The coupled bioavailability-charging-folding operator P_cftr(f_mod, f_core_block, ddg_insert) is C^infinity "
        "on (0, 1)^2 x R_+, exhibiting a non-monotonic interior optimum in f_mod bounded above by K_overmod^4."
    ),
    "falsifiable_code_and_toy": (
        "evaluate_round203_chem_mod_suptrna_lnp(regime_key) evaluates RNase half-life, TLR7/8 innate sensing, "
        "synthetase charging, NBD2 folding, and restored anion efflux across 5 curative and negative-control regimens."
    ),
}

REGIME_PRESETS: Dict[str, Dict[str, float]] = {
    "w1282x_chem_mod_suptrna_lnp_curative": {
        "f_mod": 0.42,
        "f_core_block": 0.0,
        "eta_encap": 0.94,
        "eta_esc": 0.38,
        "s_evo2_fitness": 0.95,
        "f_mrna0": 0.32,
        "ddg_insert_kcal": 0.0,
        "ddg_corr_kcal": 0.80,
    },
    "g542x_chem_mod_suptrna_lnp_curative": {
        "f_mod": 0.45,
        "f_core_block": 0.0,
        "eta_encap": 0.95,
        "eta_esc": 0.40,
        "s_evo2_fitness": 0.96,
        "f_mrna0": 0.22,
        "ddg_insert_kcal": 0.0,
        "ddg_corr_kcal": 0.0,
    },
    "unmodified_naked_suptrna_rnase_tlr78_fail": {
        "f_mod": 0.0,
        "f_core_block": 0.0,
        "eta_encap": 0.08,
        "eta_esc": 0.05,
        "s_evo2_fitness": 0.95,
        "f_mrna0": 0.22,
        "ddg_insert_kcal": 0.0,
        "ddg_corr_kcal": 0.0,
    },
    "anticodon_blocked_overmodified_inactivated_fail": {
        "f_mod": 0.55,
        "f_core_block": 0.70,
        "eta_encap": 0.94,
        "eta_esc": 0.38,
        "s_evo2_fitness": 0.88,
        "f_mrna0": 0.32,
        "ddg_insert_kcal": 0.0,
        "ddg_corr_kcal": 0.80,
    },
    "w1282x_nearcognate_nbd2_misfolding_unrescued_fail": {
        "f_mod": 0.40,
        "f_core_block": 0.0,
        "eta_encap": 0.92,
        "eta_esc": 0.36,
        "s_evo2_fitness": 0.94,
        "f_mrna0": 0.28,
        "ddg_insert_kcal": 2.85,
        "ddg_corr_kcal": 0.0,
    },
}


@dataclass(frozen=True)
class Round203ChemModSupTRNAResult:
    regime_key: str
    f_mod: float
    f_core_block: float
    eta_encap: float
    eta_esc: float
    s_evo2_fitness: float
    f_mrna0: float
    ddg_insert_kcal: float
    ddg_corr_kcal: float
    t_half_trna_h: float
    b_trna_bioavail: float
    i_tlr78_index: float
    eta_charge: float
    alpha_suptrna_readthrough: float
    p_fold_nbd: float
    f_mrna_rescued: float
    p_cftr_apical: float
    i_sc_equiv_pct_wt: float
    crosses_kahn_kalai_threshold: bool
    falsification_invariant_satisfied: bool
    complex_step_relative_error: float
    research_boundary: str = RESEARCH_BOUNDARY


def _exp_num(x: float | complex) -> float | complex:
    return cmath.exp(x) if isinstance(x, complex) else math.exp(x)


def _logistic_neg(z: float | complex) -> float | complex:
    """Overflow-safe holomorphic evaluation of 1 / (1 + exp(z)) on C."""
    z_real = z.real if isinstance(z, complex) else float(z)
    if z_real >= 0.0:
        e_neg = _exp_num(-z)
        return e_neg / (1.0 + e_neg)
    return 1.0 / (1.0 + _exp_num(z))


def compute_endosomal_stability_and_tlr78(
    f_mod: float | complex,
    eta_encap: float | complex,
    eta_esc: float | complex,
) -> Tuple[float | complex, float | complex, float | complex]:
    """Computes endosomal tRNA half-life (h), cytosolic bioavailability B_trna, and TLR7/8 index I_tlr78."""
    if isinstance(f_mod, (int, float)) and (f_mod < 0.0 or f_mod > 1.0):
        raise ValueError(f"f_mod must be in [0, 1], got {f_mod}")
    if isinstance(eta_encap, (int, float)) and (eta_encap < 0.0 or eta_encap > 1.0):
        raise ValueError(f"eta_encap must be in [0, 1], got {eta_encap}")
    if isinstance(eta_esc, (int, float)) and (eta_esc < 0.0 or eta_esc > 1.0):
        raise ValueError(f"eta_esc must be in [0, 1], got {eta_esc}")

    k_rnase = DEFAULT_K_RNASE0_PER_H * (
        (1.0 - DEFAULT_SIGMA_RNASE_MAX) + DEFAULT_SIGMA_RNASE_MAX * _exp_num(-DEFAULT_ALPHA_STAB_MOD * f_mod)
    )
    ln2 = math.log(2.0)
    t_half_h = ln2 / k_rnase
    b_trna = eta_encap * eta_esc * (t_half_h / (t_half_h + ln2 * DEFAULT_TAU_TRANSIT_H))
    i_tlr78 = DEFAULT_I_TLR0 * _exp_num(-DEFAULT_BETA_TLR_SUPPRESS * f_mod) * (0.65 + 0.35 * (1.0 - eta_esc))
    return t_half_h, b_trna, i_tlr78


def compute_synthetase_charging_and_decoding(
    f_mod: float | complex,
    f_core_block: float | complex,
    s_evo2_fitness: float | complex,
    b_trna: float | complex,
) -> Tuple[float | complex, float | complex]:
    """Computes aminoacylation efficiency eta_charge and ribosomal PTC readthrough alpha_suptrna."""
    if isinstance(f_mod, (int, float)) and (f_mod < 0.0 or f_mod > 1.0):
        raise ValueError(f"f_mod must be in [0, 1], got {f_mod}")
    if isinstance(f_core_block, (int, float)) and (f_core_block < 0.0 or f_core_block > 1.0):
        raise ValueError(f"f_core_block must be in [0, 1], got {f_core_block}")
    if isinstance(s_evo2_fitness, (int, float)) and (s_evo2_fitness < 0.0 or s_evo2_fitness > 1.0):
        raise ValueError(f"s_evo2_fitness must be in [0, 1], got {s_evo2_fitness}")
    if isinstance(b_trna, (int, float)) and (b_trna < 0.0 or b_trna > 1.0):
        raise ValueError(f"b_trna must be in [0, 1], got {b_trna}")

    overmod_penalty = 1.0 + (f_mod / DEFAULT_K_OVERMOD) ** 4
    eta_charge = s_evo2_fitness * _exp_num(-DEFAULT_KAPPA_CORE_BLOCK * f_core_block) / overmod_penalty
    active_pool = b_trna * eta_charge
    alpha_suptrna = DEFAULT_ALPHA_MAX_READTHROUGH * (active_pool / (active_pool + DEFAULT_K_ERF1_COMPETITION))
    return eta_charge, alpha_suptrna


def compute_nbd_folding_and_anion_efflux(
    alpha_suptrna: float | complex,
    f_mrna0: float | complex,
    ddg_insert_kcal: float | complex,
    ddg_corr_kcal: float | complex,
) -> Tuple[float | complex, float | complex, float | complex, float | complex]:
    """Computes NBD thermodynamic folding P_fold, NMD-shielded mRNA, apical CFTR density, and % WT anion efflux."""
    if isinstance(alpha_suptrna, (int, float)) and (alpha_suptrna < 0.0 or alpha_suptrna > 1.0):
        raise ValueError(f"alpha_suptrna must be in [0, 1], got {alpha_suptrna}")
    if isinstance(f_mrna0, (int, float)) and (f_mrna0 <= 0.0 or f_mrna0 > 1.0):
        raise ValueError(f"f_mrna0 must be in (0, 1], got {f_mrna0}")
    if isinstance(ddg_insert_kcal, (int, float)) and ddg_insert_kcal < 0.0:
        raise ValueError(f"ddg_insert_kcal must be non-negative, got {ddg_insert_kcal}")
    if isinstance(ddg_corr_kcal, (int, float)) and ddg_corr_kcal < 0.0:
        raise ValueError(f"ddg_corr_kcal must be non-negative, got {ddg_corr_kcal}")

    boltz_arg = (ddg_insert_kcal - ddg_corr_kcal - DEFAULT_DG0_FOLD_MARGIN_KCAL) / DEFAULT_RT_310K_KCAL
    p_fold = _logistic_neg(boltz_arg)
    f_mrna_eff = f_mrna0 + (1.0 - f_mrna0) * (alpha_suptrna / (alpha_suptrna + DEFAULT_K_NMD_SHIELD))
    p_cftr = f_mrna_eff * alpha_suptrna * p_fold
    denom = p_cftr + DEFAULT_K05_CFTR_ANION_EFFLUX * (1.0 - p_cftr)
    i_sc_equiv = 100.0 * (p_cftr / denom)
    return p_fold, f_mrna_eff, p_cftr, i_sc_equiv


def _eval_isc_from_f_mod(
    f_mod: float | complex,
    f_core_block: float,
    eta_encap: float,
    eta_esc: float,
    s_evo2_fitness: float,
    f_mrna0: float,
    ddg_insert_kcal: float,
    ddg_corr_kcal: float,
) -> float | complex:
    _, b_trna, _ = compute_endosomal_stability_and_tlr78(f_mod, eta_encap, eta_esc)
    _, alpha_sup = compute_synthetase_charging_and_decoding(f_mod, f_core_block, s_evo2_fitness, b_trna)
    _, _, _, i_sc = compute_nbd_folding_and_anion_efflux(alpha_sup, f_mrna0, ddg_insert_kcal, ddg_corr_kcal)
    return i_sc


def complex_step_f_mod_derivative(
    f_mod: float,
    f_core_block: float,
    eta_encap: float,
    eta_esc: float,
    s_evo2_fitness: float,
    f_mrna0: float,
    ddg_insert_kcal: float,
    ddg_corr_kcal: float,
    h: float = 1e-30,
) -> Tuple[float, float, float]:
    """Verifies d(I_sc_equiv)/d(f_mod) via complex-step perturbation f(x + 1j*1e-30).imag / 1e-30."""
    i_sc_c = _eval_isc_from_f_mod(
        complex(f_mod, h),
        f_core_block,
        eta_encap,
        eta_esc,
        s_evo2_fitness,
        f_mrna0,
        ddg_insert_kcal,
        ddg_corr_kcal,
    )
    cs_deriv = complex(i_sc_c).imag / h

    delta = 1e-6
    if f_mod - delta >= 0.0 and f_mod + delta <= 1.0:
        i_plus = float(
            _eval_isc_from_f_mod(
                f_mod + delta,
                f_core_block,
                eta_encap,
                eta_esc,
                s_evo2_fitness,
                f_mrna0,
                ddg_insert_kcal,
                ddg_corr_kcal,
            )
        )
        i_minus = float(
            _eval_isc_from_f_mod(
                f_mod - delta,
                f_core_block,
                eta_encap,
                eta_esc,
                s_evo2_fitness,
                f_mrna0,
                ddg_insert_kcal,
                ddg_corr_kcal,
            )
        )
        fd_deriv = (i_plus - i_minus) / (2.0 * delta)
    elif f_mod - delta < 0.0:
        # Second-order forward finite difference at f_mod near 0 (e.g. f_mod = 0.0)
        i_0 = float(
            _eval_isc_from_f_mod(
                f_mod, f_core_block, eta_encap, eta_esc, s_evo2_fitness, f_mrna0, ddg_insert_kcal, ddg_corr_kcal
            )
        )
        i_1 = float(
            _eval_isc_from_f_mod(
                f_mod + delta, f_core_block, eta_encap, eta_esc, s_evo2_fitness, f_mrna0, ddg_insert_kcal, ddg_corr_kcal
            )
        )
        i_2 = float(
            _eval_isc_from_f_mod(
                f_mod + 2.0 * delta, f_core_block, eta_encap, eta_esc, s_evo2_fitness, f_mrna0, ddg_insert_kcal, ddg_corr_kcal
            )
        )
        fd_deriv = (-i_2 + 4.0 * i_1 - 3.0 * i_0) / (2.0 * delta)
    else:
        # Second-order backward finite difference at f_mod near 1 (e.g. f_mod = 1.0)
        i_0 = float(
            _eval_isc_from_f_mod(
                f_mod, f_core_block, eta_encap, eta_esc, s_evo2_fitness, f_mrna0, ddg_insert_kcal, ddg_corr_kcal
            )
        )
        i_m1 = float(
            _eval_isc_from_f_mod(
                f_mod - delta, f_core_block, eta_encap, eta_esc, s_evo2_fitness, f_mrna0, ddg_insert_kcal, ddg_corr_kcal
            )
        )
        i_m2 = float(
            _eval_isc_from_f_mod(
                f_mod - 2.0 * delta, f_core_block, eta_encap, eta_esc, s_evo2_fitness, f_mrna0, ddg_insert_kcal, ddg_corr_kcal
            )
        )
        fd_deriv = (3.0 * i_0 - 4.0 * i_m1 + i_m2) / (2.0 * delta)
    denom = max(abs(cs_deriv), abs(fd_deriv), 1e-12)
    rel_err = abs(cs_deriv - fd_deriv) / denom
    return cs_deriv, fd_deriv, rel_err


def evaluate_round203_chem_mod_suptrna_lnp(regime_key: str) -> Round203ChemModSupTRNAResult:
    """Evaluates chemically modified sup-tRNA microfluidic LNP delivery and PTC rescue for a named regimen."""
    if regime_key not in REGIME_PRESETS:
        raise KeyError(f"Unknown regime '{regime_key}'. Valid regimes: {sorted(REGIME_PRESETS)}")

    p = REGIME_PRESETS[regime_key]
    f_mod = p["f_mod"]
    f_core = p["f_core_block"]
    eta_encap = p["eta_encap"]
    eta_esc = p["eta_esc"]
    s_evo2 = p["s_evo2_fitness"]
    f_mrna0 = p["f_mrna0"]
    ddg_ins = p["ddg_insert_kcal"]
    ddg_corr = p["ddg_corr_kcal"]

    t_half_v, b_trna_v, i_tlr_v = compute_endosomal_stability_and_tlr78(f_mod, eta_encap, eta_esc)
    eta_chg_v, alpha_sup_v = compute_synthetase_charging_and_decoding(f_mod, f_core, s_evo2, b_trna_v)
    p_fold_v, f_mrna_v, p_cftr_v, i_sc_v = compute_nbd_folding_and_anion_efflux(
        alpha_sup_v, f_mrna0, ddg_ins, ddg_corr
    )

    t_half = float(t_half_v)
    b_trna = float(b_trna_v)
    i_tlr = float(i_tlr_v)
    eta_chg = float(eta_chg_v)
    alpha_sup = float(alpha_sup_v)
    p_fold = float(p_fold_v)
    f_mrna = float(f_mrna_v)
    p_cftr = float(p_cftr_v)
    i_sc = float(i_sc_v)

    _, _, rel_err = complex_step_f_mod_derivative(
        f_mod, f_core, eta_encap, eta_esc, s_evo2, f_mrna0, ddg_ins, ddg_corr
    )

    inv_satisfied = (
        t_half >= 0.99
        and 0.0 <= b_trna <= 1.0
        and 0.0 <= i_tlr <= 1.0
        and 0.0 <= eta_chg <= 1.0
        and 0.0 <= alpha_sup <= 1.0
        and 0.0 <= p_fold <= 1.0
        and 0.0 <= f_mrna <= 1.0
        and 0.0 <= p_cftr <= 1.0
        and 0.0 <= i_sc <= 100.0
        and rel_err <= 1e-6
    )

    crosses_threshold = (
        t_half >= 4.0
        and i_tlr <= 0.20
        and eta_chg >= 0.65
        and i_sc >= 15.0
    )

    return Round203ChemModSupTRNAResult(
        regime_key=regime_key,
        f_mod=f_mod,
        f_core_block=f_core,
        eta_encap=eta_encap,
        eta_esc=eta_esc,
        s_evo2_fitness=s_evo2,
        f_mrna0=f_mrna0,
        ddg_insert_kcal=ddg_ins,
        ddg_corr_kcal=ddg_corr,
        t_half_trna_h=t_half,
        b_trna_bioavail=b_trna,
        i_tlr78_index=i_tlr,
        eta_charge=eta_chg,
        alpha_suptrna_readthrough=alpha_sup,
        p_fold_nbd=p_fold,
        f_mrna_rescued=f_mrna,
        p_cftr_apical=p_cftr,
        i_sc_equiv_pct_wt=i_sc,
        crosses_kahn_kalai_threshold=crosses_threshold,
        falsification_invariant_satisfied=inv_satisfied,
        complex_step_relative_error=rel_err,
    )


def evaluate_falsification_negative_control(
    poison_negative_f_mod: bool = False,
    poison_exceeded_core_block: bool = False,
    poison_invalid_encapsulation: bool = False,
    poison_negative_ddg: bool = False,
    verify_naked_and_anticodon_controls: bool = False,
) -> Dict[str, Any]:
    """Tests poisoned negative controls (ValueError + biophysical AH-10 failure modes)."""
    reasons: list[str] = []
    if poison_negative_f_mod:
        try:
            compute_endosomal_stability_and_tlr78(-0.15, 0.9, 0.35)
            return {"passed": False, "reason": "Negative f_mod accepted unexpectedly"}
        except ValueError:
            reasons.append("Deterministically rejected negative f_mod")

    if poison_exceeded_core_block:
        try:
            compute_synthetase_charging_and_decoding(0.4, 1.25, 0.9, 0.25)
            return {"passed": False, "reason": "f_core_block > 1.0 accepted unexpectedly"}
        except ValueError:
            reasons.append("Deterministically rejected f_core_block > 1.0")

    if poison_invalid_encapsulation:
        try:
            compute_endosomal_stability_and_tlr78(0.4, 1.30, 0.35)
            return {"passed": False, "reason": "eta_encap > 1.0 accepted unexpectedly"}
        except ValueError:
            reasons.append("Deterministically rejected eta_encap > 1.0")

    if poison_negative_ddg:
        try:
            compute_nbd_folding_and_anion_efflux(0.10, 0.30, -1.2, 0.0)
            return {"passed": False, "reason": "Negative ddg_insert_kcal accepted unexpectedly"}
        except ValueError:
            reasons.append("Deterministically rejected negative ddg_insert_kcal")

    if verify_naked_and_anticodon_controls:
        naked = evaluate_round203_chem_mod_suptrna_lnp("unmodified_naked_suptrna_rnase_tlr78_fail")
        blocked = evaluate_round203_chem_mod_suptrna_lnp("anticodon_blocked_overmodified_inactivated_fail")
        if naked.crosses_kahn_kalai_threshold or naked.i_tlr78_index <= 0.60 or naked.t_half_trna_h >= 2.0:
            return {"passed": False, "reason": "Unmodified naked sup-tRNA failed RNase/TLR7/8 negative control"}
        if blocked.crosses_kahn_kalai_threshold or blocked.eta_charge >= 0.05:
            return {"passed": False, "reason": "Anticodon-blocked sup-tRNA failed synthetase inactivation control"}
        reasons.append("Verified naked tRNA degradation/TLR7/8 spike and anticodon-block inactivation")

    return {"passed": True, "reason": "; ".join(reasons) if reasons else "Nominal baseline"}
