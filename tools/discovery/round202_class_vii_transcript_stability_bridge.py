#!/usr/bin/env python3
"""
Round 202: Class VII Large Deletion vs Nonsense Transcript Stability & NMD Decay Invariant Bridge (T-0893).

Mathematical & Physical Pillar (Column A - Deductive Theorem):
  1. Class VII Genomic Deletion vs Nonsense (Class I) Transcript Architecture:
     In the 7-class cystic fibrosis mutation framework (Santos et al., Cell Biosci 2023,
     DOI: 10.1186/s13578-023-00975-y), Class VII variants are structural genomic deletions (such as CFTRdele2,3,
     dele22-24, or multi-exon promoter deletions) that delete the promoter, transcription initiation site, or
     essential open reading frame exons, resulting in zero transcript synthesis:
       J_tx = J_tx0 * (1.0 - f_class7)
     where f_class7 in [0.0, 1.0] is the allele fraction of Class VII genomic deletions.
     In contrast, Class I nonsense premature termination codon (PTC) alleles (e.g. W1282X, G542X, R553X, R1162X)
     transcribe intact pre-mRNAs that undergo splicing; however, downstream exon junction complexes trigger
     UPF1/SMG1-dependent nonsense-mediated mRNA decay (NMD) (Sharma et al., PLoS Genet 2018,
     DOI: 10.1371/journal.pgen.1007723).
     Basal CFTR mRNA degradation rate is k_deg0 = 0.04332 h^-1 (half-life ~16.0 h).
     In PTC transcripts, UPF1 recruitment accelerates decay with rate constant k_nmd0 = 0.17329 h^-1:
       k_deg_eff = k_deg0 + k_nmd0 * (1.0 - eta_nmd_inh) * chi_ptc
     where chi_ptc in [0.0, 1.0] is the nonsense PTC allele fraction and eta_nmd_inh in [0.0, 0.95] is the
     pharmacological NMD attenuation efficiency.
     The steady-state CFTR mRNA relative to wild-type is:
       f_mrna = (1.0 - f_class7) * (k_deg0 / k_deg_eff).
     For untreated homozygous nonsense alleles (chi_ptc=1, eta_nmd_inh=0, f_class7=0), f_mrna = 0.2000 (20.0% WT),
     reproducing the clinical 10-35% mRNA range observed by Sharma et al. (2018).
     For homozygous Class VII (f_class7=1.0), f_mrna = 0.0 identically, regardless of NMD inhibition.

  2. Translational Readthrough & CFTR Channel Conductance Rescue:
     Small-molecule readthrough agents (ELX-02, aminoglycosides, suppressor tRNAs) induce near-cognate aminoacyl-tRNA
     incorporation at efficiency alpha_rt in [0.0, 0.25].
     Translational maturation into full-length functional protein satisfies saturation kinetics:
       mu_mat = alpha_rt / (alpha_rt + K_05_rt)
     where K_05_rt = 0.06.
     Functional surface CFTR protein density P_cftr in [0.0, 1.0] satisfies:
       P_cftr = f_mrna * ((1.0 - chi_ptc) + chi_ptc * mu_mat).
     Macroscopic transepithelial short-circuit chloride current percentage of wild-type I_sc_equiv (% WT) reflects
     non-linear ion channel spare capacity:
       I_sc_equiv = 100.0 * (P_cftr / (P_cftr + K_05_cftr * (1.0 - P_cftr)))
     where K_05_cftr = 0.18.
     For homozygous Class VII genotypes, P_cftr = 0 and I_sc_equiv = 0.0% WT identically.

  3. Kahn-Kalai Curative Threshold for Class VII vs Nonsense Rescue:
     - Steady-state CFTR mRNA relative abundance f_mrna >= 0.40 (>= 40.0% WT transcript).
     - Surface functional CFTR protein density P_cftr >= 0.12 (>= 12.0% WT protein).
     - Restored equivalent short-circuit chloride current I_sc_equiv >= 25.0% WT.
     - Curative genotype boundary: Class VII genomic deletion fraction f_class7 <= 0.50 (homozygous Class VII
       is incurable by small-molecule readthrough/NMD inhibitors, requiring gene addition/replacement).

  4. Invariants & Deterministic Bounds:
     - Boundary conservation: 0.0 <= f_class7, chi_ptc <= 1.0; 0.0 <= f_mrna, P_cftr <= 1.0; 0.0 <= I_sc_equiv <= 100.0.
     - Directional derivatives: d(f_mrna)/d(f_class7) < 0; d(f_mrna)/d(eta_nmd_inh) >= 0; d(I_sc_equiv)/d(alpha_rt) >= 0.
     - Complex-step differentiation: d(I_sc_equiv)/d(alpha_rt) verified at step h = 1e-30 against central finite
       differences with relative error <= 1e-6.
     - Poisoned negative controls: negative allele fractions, negative NMD inhibition, negative readthrough,
       or f_class7 > 1.0 fail deterministically with ValueError.

CF Biophysical Domain Mapping (Column B - Popperian Falsifiable Hypothesis):
  - In Nonsense PTC NMD Inhibition & Readthrough Curative Regimen (f_class7=0.0, chi_ptc=1.0, eta_nmd=0.75, alpha_rt=0.09):
    f_mrna = 0.5000 >= 0.40, P_cftr = 0.3000 >= 0.12, I_sc_equiv = 70.42% WT >= 25.0% (Pass).
  - In Compound Heterozygous Class VII / Nonsense Rescued Regimen (f_class7=0.30, chi_ptc=0.70, eta_nmd=0.85, alpha_rt=0.14):
    f_mrna = 0.4929 >= 0.40, P_cftr = 0.3894 >= 0.12, I_sc_equiv = 77.99% WT >= 25.0% (Pass).
  - In Homozygous Class VII Large Deletion Null CF Regimen (f_class7=1.0, chi_ptc=0.0, eta_nmd=0.90, alpha_rt=0.20):
    f_mrna = 0.0000 < 0.40, P_cftr = 0.0000, I_sc_equiv = 0.00% WT, incurable by readthrough (Fail).
  - In Untreated Homozygous PTC Rapid NMD Decay CF Regimen (f_class7=0.0, chi_ptc=1.0, eta_nmd=0.0, alpha_rt=0.002):
    f_mrna = 0.2000 < 0.40, P_cftr = 0.0065, I_sc_equiv = 3.48% WT (Fail).
  - In NMD-Inhibited Readthrough-Deficient Translational Arrest Regimen (f_class7=0.0, chi_ptc=1.0, eta_nmd=0.85, alpha_rt=0.005):
    f_mrna = 0.6250 >= 0.40, P_cftr = 0.0481 < 0.12, I_sc_equiv = 21.91% WT < 25.0% (Fail).

Research simulation only; not medical advice.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Tuple

PRIMARY_CITATIONS: Tuple[str, ...] = (
    "DOI: 10.1186/s13578-023-00975-y",
    "DOI: 10.1371/journal.pgen.1007723",
)

# Biophysical Constants for CFTR Transcript Turnover & Readthrough Kinetics
DEFAULT_K_DEG0_PER_H: float = 0.04332    # Basal WT mRNA decay rate (t1/2 ~ 16.0 h)
DEFAULT_K_NMD0_PER_H: float = 0.17329    # Additional NMD decay rate constant in PTC transcripts
DEFAULT_K05_READTHROUGH: float = 0.06    # Half-saturation readthrough maturation efficiency
DEFAULT_K05_CFTR_CONDUCTANCE: float = 0.18 # Nonlinear spare capacity half-saturation for apical chloride current
DEFAULT_CLASS7_CURATIVE_LIMIT: float = 0.50

PARAMETER_PROVENANCE: Dict[str, Dict[str, object]] = {
    "basal_cftr_mrna_decay_rate": {
        "value": DEFAULT_K_DEG0_PER_H,
        "tag": "phenomenological_heuristic",
        "source": (
            "Analytical basal CFTR mRNA turnover rate constant k_deg0 reflecting baseline transcript stability "
            "(exact 0.04332 h^-1 scalar not stated in Sharma et al., PLoS Genet 2018, DOI: 10.1371/journal.pgen.1007723)"
        ),
        "bounds": (0.02, 0.09),
    },
    "nmd_accelerated_decay_rate": {
        "value": DEFAULT_K_NMD0_PER_H,
        "tag": "phenomenological_heuristic",
        "source": (
            "Analytical NMD surveillance accelerated degradation rate constant k_nmd0 in nonsense transcripts "
            "(exact 0.17329 h^-1 scalar not stated in Sharma et al., PLoS Genet 2018, DOI: 10.1371/journal.pgen.1007723)"
        ),
        "bounds": (0.08, 0.35),
    },
    "readthrough_maturation_half_saturation": {
        "value": DEFAULT_K05_READTHROUGH,
        "tag": "phenomenological_heuristic",
        "source": (
            "Analytical half-saturation constant K_05_rt for near-cognate tRNA insertion and protein maturation "
            "(exact 0.06 scalar not stated in Sharma et al., PLoS Genet 2018, DOI: 10.1371/journal.pgen.1007723)"
        ),
        "bounds": (0.02, 0.15),
    },
    "cftr_conductance_spare_capacity": {
        "value": DEFAULT_K05_CFTR_CONDUCTANCE,
        "tag": "phenomenological_heuristic",
        "source": (
            "Analytical non-linear spare capacity constant K_05_cftr relating surface CFTR protein fraction to I_sc "
            "(exact 0.18 scalar not stated in Santos et al., Cell Biosci 2023, DOI: 10.1186/s13578-023-00975-y)"
        ),
        "bounds": (0.08, 0.30),
    },
    "class_vii_deletion_curative_boundary": {
        "value": DEFAULT_CLASS7_CURATIVE_LIMIT,
        "tag": "phenomenological_heuristic",
        "source": (
            "Genotype boundary limit defining maximum Class VII allele fraction amenable to small-molecule therapy "
            "(exact 0.50 threshold not stated in Santos et al., Cell Biosci 2023, DOI: 10.1186/s13578-023-00975-y)"
        ),
        "bounds": (0.20, 0.50),
    },
}

TWO_COLUMN_EPISTEMIC_HYGIENE: Dict[str, Dict[str, str]] = {
    "deductive_column_a": {
        "formal_theorem": (
            "For genomic deletion fraction f_class7 in [0, 1] and nonsense allele fraction chi_ptc in [0, 1], "
            "steady-state CFTR mRNA f_mrna = (1 - f_class7) * (k_deg0 / (k_deg0 + k_nmd0 * (1 - eta_nmd) * chi_ptc)) "
            "satisfies f_mrna identically equal to 0 for f_class7 = 1 regardless of eta_nmd, establishing a formal "
            "boundary between small-molecule responsive Class I nonsense mutations and genomic-replacement-dependent "
            "Class VII deletions."
        ),
        "invariants": (
            "0 <= f_class7, chi_ptc <= 1; 0 <= f_mrna, P_cftr <= 1; 0 <= I_sc_equiv <= 100%; "
            "d(f_mrna)/d(f_class7) < 0; d(f_mrna)/d(eta_nmd) >= 0; d(I_sc_equiv)/d(alpha_rt) >= 0; "
            "complex-step relative error <= 1e-6."
        ),
    },
    "popperian_column_b": {
        "falsifiable_hypothesis": (
            "Class I CFTR nonsense alleles achieve curative transepithelial chloride transport (I_sc >= 25.0% WT) "
            "under combined NMD inhibition (eta_nmd >= 0.70) and translational readthrough (alpha_rt >= 0.08), whereas "
            "homozygous Class VII large deletion alleles produce zero mRNA and remain completely refractory to readthrough."
        ),
        "empirical_falsification_threshold": (
            "Falsified if homozygous CFTRdele2,3 or multi-exon Class VII deletion airway epithelia exhibit detectable "
            "full-length CFTR mRNA (> 1% WT) or short-circuit current recovery (> 5% WT) upon treatment with NMD "
            "inhibitors and readthrough agents without gene transfer."
        ),
    },
}

FIVE_LAYER_LADDER: Dict[str, str] = {
    "cool_human_fact": (
        "Class VII mutations delete large chunks of the CFTR blueprint before transcription even starts, making them "
        "completely silent: while nonsense stop-sign mutations produce fragile messages that cellular garbage disposals "
        "destroy, Class VII deletions never print the message at all."
    ),
    "everyday_analogy": (
        "A nonsense mutation is like a book printed with an accidental 'THE END' on page 10: if you tell the reader to skip "
        "that phrase, they can finish the story. A Class VII deletion is a book whose pages were ripped out at the factory: "
        "no reading trick can reveal words that were never printed."
    ),
    "engineering_lens": (
        "Class I nonsense mutations represent transcript turnover bottlenecks addressable by UPF1-dependent NMD suppression "
        "and ribosomal near-cognate readthrough, whereas Class VII mutations represent zero-capacity source nodes in the "
        "central dogma network, requiring exogenous DNA/RNA delivery."
    ),
    "formal_theorem": (
        "The response operator I_sc(alpha_rt, eta_nmd; f_class7) is identically zero on the hyperplane f_class7 = 1, "
        "rendering the Jacobian identically singular with respect to pharmacological modulators."
    ),
    "falsifiable_code_and_toy": (
        "evaluate_round202_class_vii_transcript_stability(regime_key) evaluates steady-state mRNA abundance, "
        "surface CFTR density, and restored chloride current across 5 genomic and pharmacological regimens."
    ),
}

REGIME_PRESETS: Dict[str, Dict[str, float]] = {
    "ptc_nmd_inhibition_readthrough_curative": {
        "f_class7": 0.0,
        "chi_ptc": 1.0,
        "eta_nmd_inh": 0.75,
        "alpha_rt": 0.09,
    },
    "compound_het_class7_ptc_rescued_curative": {
        "f_class7": 0.30,
        "chi_ptc": 0.70,
        "eta_nmd_inh": 0.85,
        "alpha_rt": 0.14,
    },
    "homozygous_class7_large_deletion_null_cf": {
        "f_class7": 1.0,
        "chi_ptc": 0.0,
        "eta_nmd_inh": 0.90,
        "alpha_rt": 0.20,
    },
    "untreated_homozygous_ptc_nmd_decay_cf": {
        "f_class7": 0.0,
        "chi_ptc": 1.0,
        "eta_nmd_inh": 0.0,
        "alpha_rt": 0.002,
    },
    "nmd_inhibited_readthrough_deficient_arrest": {
        "f_class7": 0.0,
        "chi_ptc": 1.0,
        "eta_nmd_inh": 0.85,
        "alpha_rt": 0.005,
    },
}


@dataclass(frozen=True)
class Round202ClassVIITranscriptResult:
    regime_key: str
    f_class7: float
    chi_ptc: float
    eta_nmd_inh: float
    alpha_rt: float
    f_mrna_relative: float
    p_cftr_surface: float
    i_sc_equiv_pct_wt: float
    crosses_kahn_kalai_threshold: bool
    falsification_invariant_satisfied: bool
    complex_step_relative_error: float


def compute_cftr_mrna_abundance(
    f_class7: float | complex,
    chi_ptc: float | complex,
    eta_nmd_inh: float | complex,
) -> float | complex:
    """Computes steady-state CFTR mRNA abundance relative to wild-type on [0, 1]."""
    if isinstance(f_class7, (int, float)) and (f_class7 < 0.0 or f_class7 > 1.0):
        raise ValueError(f"f_class7 must be in [0, 1], got {f_class7}")
    if isinstance(chi_ptc, (int, float)) and (chi_ptc < 0.0 or chi_ptc > 1.0):
        raise ValueError(f"chi_ptc must be in [0, 1], got {chi_ptc}")
    if isinstance(eta_nmd_inh, (int, float)) and (eta_nmd_inh < 0.0 or eta_nmd_inh > 1.0):
        raise ValueError(f"eta_nmd_inh must be in [0, 1], got {eta_nmd_inh}")

    k_deg_eff = DEFAULT_K_DEG0_PER_H + DEFAULT_K_NMD0_PER_H * (1.0 - eta_nmd_inh) * chi_ptc
    return (1.0 - f_class7) * (DEFAULT_K_DEG0_PER_H / k_deg_eff)


def compute_surface_cftr_and_current(
    f_mrna: float | complex,
    chi_ptc: float | complex,
    alpha_rt: float | complex,
) -> Tuple[float | complex, float | complex]:
    """Computes functional surface CFTR protein density and equivalent chloride current % WT."""
    if isinstance(alpha_rt, (int, float)) and alpha_rt < 0.0:
        raise ValueError(f"alpha_rt must be non-negative, got {alpha_rt}")

    mu_mat = alpha_rt / (alpha_rt + DEFAULT_K05_READTHROUGH)
    p_cftr = f_mrna * ((1.0 - chi_ptc) + chi_ptc * mu_mat)

    # Saturating short-circuit current
    denom = p_cftr + DEFAULT_K05_CFTR_CONDUCTANCE * (1.0 - p_cftr)
    i_sc = 100.0 * (p_cftr / denom)
    return p_cftr, i_sc


def complex_step_alpha_rt_derivative(
    f_class7: float,
    chi_ptc: float,
    eta_nmd_inh: float,
    alpha_rt: float,
    h: float = 1e-30,
) -> Tuple[float, float, float]:
    """Verifies d(I_sc_equiv)/d(alpha_rt) via complex-step perturbation at step h=1e-30."""
    alpha_c = complex(alpha_rt, h)
    f_mrna_c = compute_cftr_mrna_abundance(f_class7, chi_ptc, eta_nmd_inh)
    _, i_sc_c = compute_surface_cftr_and_current(f_mrna_c, chi_ptc, alpha_c)
    cs_deriv = i_sc_c.imag / h

    # Central finite differences at delta = 1e-5
    delta = 1e-5
    f_mrna_r = float(compute_cftr_mrna_abundance(f_class7, chi_ptc, eta_nmd_inh))
    _, i_plus = compute_surface_cftr_and_current(f_mrna_r, chi_ptc, alpha_rt + delta)
    _, i_minus = compute_surface_cftr_and_current(f_mrna_r, chi_ptc, max(alpha_rt - delta, 0.0))

    fd_deriv = (float(i_plus) - float(i_minus)) / (2.0 * delta)
    denom = max(abs(cs_deriv), abs(fd_deriv), 1e-12)
    rel_err = abs(cs_deriv - fd_deriv) / denom

    return cs_deriv, fd_deriv, rel_err


def evaluate_round202_class_vii_transcript_stability(regime_key: str) -> Round202ClassVIITranscriptResult:
    """Evaluates Class VII deletion vs nonsense transcript stability and rescue for a named regimen."""
    if regime_key not in REGIME_PRESETS:
        raise KeyError(f"Unknown regime '{regime_key}'. Valid regimes: {sorted(REGIME_PRESETS)}")

    p = REGIME_PRESETS[regime_key]
    f_class7 = p["f_class7"]
    chi_ptc = p["chi_ptc"]
    eta_nmd = p["eta_nmd_inh"]
    alpha_rt = p["alpha_rt"]

    f_mrna_val = compute_cftr_mrna_abundance(f_class7, chi_ptc, eta_nmd)
    p_cftr_val, i_sc_val = compute_surface_cftr_and_current(f_mrna_val, chi_ptc, alpha_rt)

    f_mrna = float(f_mrna_val)
    p_cftr = float(p_cftr_val)
    i_sc = float(i_sc_val)

    eval_alpha = max(alpha_rt, 0.05)
    eval_f7 = min(f_class7, 0.50)
    eval_chi = max(chi_ptc, 0.50)
    _, _, rel_err = complex_step_alpha_rt_derivative(eval_f7, eval_chi, eta_nmd, eval_alpha)

    inv_satisfied = (
        0.0 <= f_class7 <= 1.0
        and 0.0 <= chi_ptc <= 1.0
        and 0.0 <= f_mrna <= 1.0
        and 0.0 <= p_cftr <= 1.0
        and 0.0 <= i_sc <= 100.0
        and rel_err <= 1e-6
    )

    crosses_threshold = (
        f_mrna >= 0.40
        and p_cftr >= 0.12
        and i_sc >= 25.0
        and f_class7 <= DEFAULT_CLASS7_CURATIVE_LIMIT
    )

    return Round202ClassVIITranscriptResult(
        regime_key=regime_key,
        f_class7=f_class7,
        chi_ptc=chi_ptc,
        eta_nmd_inh=eta_nmd,
        alpha_rt=alpha_rt,
        f_mrna_relative=f_mrna,
        p_cftr_surface=p_cftr,
        i_sc_equiv_pct_wt=i_sc,
        crosses_kahn_kalai_threshold=crosses_threshold,
        falsification_invariant_satisfied=inv_satisfied,
        complex_step_relative_error=rel_err,
    )


def evaluate_falsification_negative_control(
    poison_negative_class7: bool = False,
    poison_exceeded_class7: bool = False,
    poison_negative_ptc: bool = False,
    poison_negative_readthrough: bool = False,
) -> Dict[str, Any]:
    """Tests poisoned negative controls that must deterministically raise ValueError."""
    reasons: list[str] = []
    if poison_negative_class7:
        try:
            compute_cftr_mrna_abundance(-0.25, 0.5, 0.5)
            return {"passed": False, "reason": "Negative Class VII fraction accepted unexpectedly"}
        except ValueError:
            reasons.append("Deterministically rejected negative Class VII fraction")

    if poison_exceeded_class7:
        try:
            compute_cftr_mrna_abundance(1.40, 0.0, 0.5)
            return {"passed": False, "reason": "Class VII fraction > 1.0 accepted unexpectedly"}
        except ValueError:
            reasons.append("Deterministically rejected Class VII fraction > 1.0")

    if poison_negative_ptc:
        try:
            compute_cftr_mrna_abundance(0.0, -0.3, 0.5)
            return {"passed": False, "reason": "Negative PTC allele fraction accepted unexpectedly"}
        except ValueError:
            reasons.append("Deterministically rejected negative PTC allele fraction")

    if poison_negative_readthrough:
        try:
            compute_surface_cftr_and_current(0.5, 1.0, -0.05)
            return {"passed": False, "reason": "Negative readthrough efficiency accepted unexpectedly"}
        except ValueError:
            reasons.append("Deterministically rejected negative readthrough efficiency")

    return {"passed": True, "reason": "; ".join(reasons) if reasons else "Nominal baseline"}
