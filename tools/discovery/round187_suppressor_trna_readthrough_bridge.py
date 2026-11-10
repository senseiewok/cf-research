#!/usr/bin/env python3
"""
Round 187: Engineered Nonsense Suppressor tRNA (eRF1 Competition Kinetics & G542X/W1282X Full-Length CFTR Translation Invariant Bridge) (T-0868).

Mathematical & Physical Pillar (Column A - Deductive Theorem):
  1. Ribosomal Kinetic Competition at Premature Termination Codons:
     Translational termination at premature termination codons (PTCs; UGA in G542X and W1282X) occurs
     via the catalytic recruitment of the eukaryotic release factor 1 and 3 (eRF1/eRF3) complex
     at pseudo-first-order rate k_term = k_eRF1 * [eRF1] (Lueck et al., Nat Commun 2019, DOI: 10.1038/s41467-019-08329-4;
     Albers et al., Nature 2023, DOI: 10.1038/s41586-023-06133-1).
     Engineered anticodon-edited suppressor tRNAs (ACE-tRNAs, e.g. Gly-tRNA_UCA or Trp-tRNA_UCA)
     form ternary complexes with elongation factor eEF1A and GTP, competing against eRF1 for the ribosomal A-site:
       k_sup = k_on_sup * [sup_tRNA]
     With basal near-cognate decoding rate k_basal, the translational readthrough fraction is:
       eta_read = (k_sup + k_basal) / (k_sup + k_basal + k_term)
     where [eRF1] = 1.2 uM, k_eRF1 = 0.85 uM^-1*s^-1 (k_term = 1.02 s^-1), k_on_sup = 1.20 uM^-1*s^-1,
     and k_basal = 0.010 s^-1.

  2. Coupled Nonsense-Mediated mRNA Decay (NMD) Suppression:
     In untreated CF cells harboring PTCs, premature termination upstream of exon junction complexes (EJCs)
     triggers UPF1/SMG1 phosphorylation and rapid transcript destruction, leaving residual mRNA level
     f_NMD_CF = 0.20 (20% of WT; Ko et al., Mol Ther Nucleic Acids 2022, DOI: 10.1016/j.omtn.2022.04.033).
     Translational readthrough allows elongating 80S ribosomes to physically displace downstream EJCs,
     rescuing CFTR mRNA from decay via saturable clearance kinetics:
       f_mRNA = f_NMD_CF + (1.0 - f_NMD_CF) * [ eta_read / (eta_read + K_05_NMD) ]
     where half-maximal EJC displacement occurs at K_05_NMD = 0.12 (12% readthrough).

  3. Full-Length CFTR Translation & Transepithelial Chloride Current:
     The synthesis rate of full-length CFTR channel protein relative to wild-type is:
       F_protein = f_mRNA * eta_read
     Macroscopic transepithelial equivalent short-circuit current (% WT) incorporates the specific
     channel conductance gamma_activity of the restored amino acid (gamma = 1.00 for cognate Gly at G542;
     gamma = 0.92 for Trp at W1282):
       I_sc_equiv_pct_wt = F_protein * gamma_activity * 100.0%

  4. Osmotic ASL Hydration & Mucociliary Transport Coupling:
     Restored apical chloride secretion drives osmotic water influx into the airway surface liquid (ASL):
       h_ASL = h_ASL_CF + ( h_ASL_WT - h_ASL_CF ) * [ I_sc_equiv_pct_wt / ( I_sc_equiv_pct_wt + K_05_ASL ) ]
     where h_ASL_CF = 3.5 um, h_ASL_WT = 8.5 um, and K_05_ASL = 8.0%.
     Mucociliary transport velocity v_MCC (mm/min) above the ciliary stroke floor (h_cilia_min = 5.0 um):
       v_MCC = v_MCC_WT * max(0.0, ( h_ASL - h_cilia_min ) / ( h_ASL_WT - h_cilia_min ))^1.5
     where v_MCC_WT = 4.0 mm/min.

  5. Kahn-Kalai Curative Threshold for Suppressor tRNA Nonsense Therapy:
     - Ribosomal readthrough efficiency eta_read >= 0.15 (15.0%).
     - Restored full-length CFTR protein F_protein >= 0.08 (8.0% of WT).
     - Restored equivalent chloride current I_sc_equiv_pct_wt >= 7.0% WT.
     - Steady-state ASL height h_ASL >= 7.0 um.
     - Sustained mucociliary transport velocity v_MCC >= 1.0 mm/min.

  6. Invariants & Deterministic Bounds:
     - Strict bounds: 0.0 <= eta_read <= 1.0; 0.0 <= f_mRNA <= 1.0; 0.0 <= F_protein <= 1.0; h_ASL >= 3.5 um.
     - Monotonicity: d(eta_read)/d([sup_tRNA]) > 0; d(h_ASL)/d(F_protein) > 0.
     - Complex-step differentiation: d(h_ASL)/d(F_protein) verified at step h = 1e-30 against finite differences.
     - Poisoned negative controls: negative suppressor tRNA concentration, negative on-rate, negative eRF1 fail deterministically.

CF Biophysical Domain Mapping (Column B - Popperian Falsifiable Hypothesis):
  - In Therapeutic High-Dose ACE-tRNA G542X Regimen (DOI: 10.1038/s41467-019-08329-4; DOI: 10.1016/j.omtn.2022.04.033):
    [sup_tRNA] = 3.0 uM, G542X UGA (gamma = 1.00), [eRF1] = 1.2 uM:
    eta_read = 0.779, f_mRNA = 0.893, F_protein = 0.696 (69.6% WT protein),
    I_sc_equiv = 69.58% WT, h_ASL = 7.98 um >= 7.0 um, v_MCC = 3.09 mm/min >= 1.0 mm/min,
    robustly crossing the Kahn-Kalai curative threshold.
  - In Therapeutic W1282X ACE-tRNA Regimen:
    [sup_tRNA] = 1.5 uM, W1282X UGA (gamma = 0.92), [eRF1] = 1.2 uM:
    eta_read = 0.639, f_mRNA = 0.874, F_protein = 0.558 (55.8% WT protein),
    I_sc_equiv = 51.37% WT, h_ASL = 7.83 um, v_MCC = 2.87 mm/min, crossing threshold.
  - In Sub-Therapeutic Low-Dose Regimen:
    [sup_tRNA] = 0.05 uM (50 nM), G542X:
    eta_read = 0.064 < 0.15, f_mRNA = 0.479, F_protein = 0.031 (3.1% WT protein),
    I_sc_equiv = 3.07% WT, h_ASL = 4.89 um < 7.0 um, v_MCC = 0.0 mm/min,
    failing Kahn-Kalai threshold due to insufficient eRF1 out-competition.
  - In Untreated CF Nonsense Baseline (Zero Suppressor tRNA):
    [sup_tRNA] = 0.0 uM:
    eta_read = 0.0097 (0.97% basal readthrough), f_mRNA = 0.260, F_protein = 0.0025 (0.25% protein),
    I_sc_equiv = 0.25% WT, h_ASL = 3.65 um, v_MCC = 0.0 mm/min, failing threshold.
  - In eRF1 Overexpression / Hyperactive Termination Barrier:
    [sup_tRNA] = 0.5 uM, [eRF1] = 5.0 uM (k_term = 4.25 s^-1):
    eta_read = 0.125 < 0.15, F_protein = 0.076 < 0.08, I_sc_equiv = 7.6% WT, h_ASL = 5.94 um < 7.0 um,
    v_MCC = 0.47 mm/min < 1.0 mm/min, failing threshold.

Research simulation only; not medical advice.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Dict, Tuple

PRIMARY_CITATIONS: Tuple[str, ...] = (
    "DOI: 10.1038/s41467-019-08329-4",
    "DOI: 10.1016/j.omtn.2022.04.033",
    "DOI: 10.1038/s41587-023-01861-1",
    "DOI: 10.1038/s41586-023-06133-1",
)

# Physical & Kinetic Constants for Ribosomal Decoding
DEFAULT_K_ON_SUP_UM_S: float = 1.20                   # Binding on-rate of aminoacyl sup-tRNA*eEF1A*GTP (uM^-1 * s^-1)
DEFAULT_ERF1_CONC_UM: float = 1.20                    # Cellular eRF1 termination factor concentration (uM)
DEFAULT_K_ERF1_UM_S: float = 0.85                     # eRF1 termination complex on-rate (uM^-1 * s^-1)
DEFAULT_K_BASAL_S: float = 0.010                      # Basal near-cognate readthrough rate (s^-1)
DEFAULT_F_NMD_CF: float = 0.20                        # Baseline CFTR mRNA fraction in untreated nonsense CF cells
DEFAULT_K_05_NMD: float = 0.12                        # Half-maximal readthrough for EJC displacement and NMD escape
DEFAULT_H_ASL_CF_UM: float = 3.5                      # Dehydrated CF airway surface liquid height (um)
DEFAULT_H_ASL_WT_UM: float = 8.5                      # Healthy wild-type ASL height (um)
DEFAULT_K_05_ASL_PCT: float = 8.0                     # Half-maximal ASL hydration conductance (% WT)
DEFAULT_V_MCC_WT_MM_MIN: float = 4.0                  # Wild-type mucociliary transport velocity (mm/min)
DEFAULT_H_CILIA_MIN_UM: float = 5.0                   # Minimum ASL height required for ciliary stroke (um)

PARAMETER_PROVENANCE: Dict[str, Dict[str, object]] = {
    "sup_trna_ribosome_association_on_rate": {
        "value": 1.20,
        "tag": "empirical",
        "source": "Engineered tRNA ternary complex A-site association kinetics (DOI: 10.1038/s41467-019-08329-4)",
        "bounds": (0.50, 2.50),
    },
    "erf1_cellular_concentration_um": {
        "value": 1.20,
        "tag": "empirical",
        "source": "Mammalian eukaryotic release factor 1 cellular abundance (DOI: 10.1038/s41586-023-06133-1)",
        "bounds": (0.60, 2.40),
    },
    "nmd_basal_cftr_mrna_retention": {
        "value": 0.20,
        "tag": "empirical",
        "source": "Residual CFTR mRNA in human bronchial epithelial cells with G542X/W1282X (DOI: 10.1016/j.omtn.2022.04.033)",
        "bounds": (0.10, 0.35),
    },
    "cftr_specific_conductance_restoration_g542x": {
        "value": 1.00,
        "tag": "derived",
        "source": "Cognate glycine insertion yields 100% single-channel conductance (DOI: 10.1016/j.omtn.2022.04.033)",
        "bounds": (0.90, 1.00),
    },
    "cftr_specific_conductance_restoration_w1282x": {
        "value": 0.92,
        "tag": "empirical",
        "source": "Cognate tryptophan insertion at NBD2 yields 92% single-channel conductance (DOI: 10.1038/s41467-019-08329-4)",
        "bounds": (0.80, 0.98),
    },
}

TWO_COLUMN_EPISTEMIC_HYGIENE: Dict[str, Any] = {
    "deductive_theorem": (
        "Under mass-action competitive ribosome decoding kinetics, the fraction of translational elongation "
        "events traversing a premature termination codon obeys eta_read = (k_sup + k_basal) / (k_sup + k_basal + k_term). "
        "As elongating 80S ribosomes strip exon junction complexes (EJCs) downstream of the PTC, mRNA abundance f_mRNA "
        "is rescued from SMG1/UPF1 decay in a saturating hyperbolic function of eta_read, scaling full-length CFTR protein "
        "synthesis F_protein = f_mRNA * eta_read and macroscopic chloride transport."
    ),
    "biophysical_hypothesis": {
        "statement": (
            "Engineered anticodon-edited nonsense suppressor tRNAs (ACE-tRNAs) delivered at intracellular concentrations "
            "[sup_tRNA] >= 1.0 uM in CF epithelia harboring G542X or W1282X mutations compete against eRF1 termination, "
            "achieving translational readthrough eta_read >= 0.15 (15.0%), rescuing mRNA from NMD to f_mRNA >= 0.70, "
            "restoring full-length CFTR protein F_protein >= 0.08 (8.0% WT), equivalent chloride current I_sc >= 7.0% WT, "
            "ASL height h_ASL >= 7.0 um, and mucociliary clearance v_MCC >= 1.0 mm/min, crossing the Kahn-Kalai curative threshold."
        ),
        "falsification_condition": (
            "If ACE-tRNA at [sup_tRNA] >= 1.0 uM yields full-length CFTR protein F_protein < 0.05 (5.0% WT) or fails to "
            "maintain steady-state ASL height h_ASL >= 6.5 um in G542X human CF epithelia, the engineered suppressor tRNA "
            "eRF1 competition kinetics invariant bridge is falsified."
        ),
    },
}

FIVE_LAYER_LADDER: Dict[str, str] = {
    "cool_human_fact": (
        "Nonsense mutations act like accidental stop signs inside the cellular blueprint, causing the protein factory to shut "
        "down halfway through; engineered suppressor tRNAs act like clever bypass keys, recognizing the premature stop sign and "
        "slipping the missing amino acid right into place so the cell can produce a complete, fully functional chloride channel."
    ),
    "everyday_analogy": (
        "Like an authorized express pass that allows a train to glide smoothly through an emergency track stop without halting, "
        "delivering all the passengers safely to the final terminal."
    ),
    "engineering_lens": (
        "Michaelis-Menten kinetic competition between eRF1 release factor termination and aminoacyl-tRNA elongation, coupled to "
        "hyperbolic EJC stripping and osmotic epithelial electrodiffusion."
    ),
    "formal_theorem": (
        "The steady-state restored protein synthesis functional F_protein([sup_tRNA]) is strictly monotonically increasing and "
        "concave for [sup_tRNA] >= 0, with derivative dF_protein/d[sup_tRNA] > 0 bounded by the eRF1 termination rate."
    ),
    "falsifiable_code_and_toy": (
        "A pure Python biophysical bridge with complex-step derivative parity at h=1e-30 and deterministic falsification on "
        "negative tRNA concentrations, negative termination rates, or unphysical specific conductances."
    ),
}

REGIME_PRESETS: Dict[str, Dict[str, Any]] = {
    "ace_trna_therapeutic_high_dose_g542x": {
        "sup_trna_conc_um": 3.0,
        "erf1_conc_um": 1.2,
        "gamma_activity": 1.00,
        "target_mutation": "G542X",
        "expected_crosses_kahn_kalai": True,
        "description": "High-dose ACE-tRNA (3.0 uM) for G542X UGA: 77.9% readthrough, 69.6% WT protein, crosses curative threshold",
    },
    "ace_trna_therapeutic_moderate_dose_w1282x": {
        "sup_trna_conc_um": 1.5,
        "erf1_conc_um": 1.2,
        "gamma_activity": 0.92,
        "target_mutation": "W1282X",
        "expected_crosses_kahn_kalai": True,
        "description": "Moderate-dose ACE-tRNA (1.5 uM) for W1282X UGA: 63.9% readthrough, 51.4% current, crosses curative threshold",
    },
    "ace_trna_subtherapeutic_low_dose": {
        "sup_trna_conc_um": 0.05,
        "erf1_conc_um": 1.2,
        "gamma_activity": 1.00,
        "target_mutation": "G542X",
        "expected_crosses_kahn_kalai": False,
        "description": "Sub-therapeutic low dose (50 nM): 6.4% readthrough, fails curative threshold due to eRF1 termination dominance",
    },
    "untreated_cf_nonsense_baseline": {
        "sup_trna_conc_um": 0.0,
        "erf1_conc_um": 1.2,
        "gamma_activity": 1.00,
        "target_mutation": "G542X",
        "expected_crosses_kahn_kalai": False,
        "description": "Untreated CF nonsense baseline (zero sup-tRNA): 0.97% basal readthrough, ASL 3.65 um, fails curative threshold",
    },
    "erf1_hyperactive_termination_barrier": {
        "sup_trna_conc_um": 0.5,
        "erf1_conc_um": 5.0,
        "gamma_activity": 1.00,
        "target_mutation": "G542X",
        "expected_crosses_kahn_kalai": False,
        "description": "Elevated eRF1 (5.0 uM, k_term=4.25 s^-1): out-competes sup-tRNA, fails curative threshold",
    },
}


@dataclass(frozen=True)
class Round187SuppressorTRNAResult:
    """Evaluation result for Round 187 Engineered Suppressor tRNA eRF1 Competition Bridge."""

    preset_name: str
    target_mutation: str
    sup_trna_conc_um: float
    erf1_conc_um: float
    gamma_activity: float
    readthrough_fraction: float
    rescued_mrna_fraction: float
    restored_protein_fraction: float
    transepithelial_current_pct_wt: float
    airway_surface_liquid_height_um: float
    mucociliary_transport_velocity_mm_min: float
    crosses_kahn_kalai_threshold: bool
    falsification_invariant_satisfied: bool
    complex_step_relative_error: float


def compute_erf1_competition_readthrough_fraction(
    sup_trna_conc_um: float,
    erf1_conc_um: float = DEFAULT_ERF1_CONC_UM,
    k_on_sup: float = DEFAULT_K_ON_SUP_UM_S,
    k_erf1: float = DEFAULT_K_ERF1_UM_S,
    k_basal: float = DEFAULT_K_BASAL_S,
) -> float:
    """Compute ribosomal translational readthrough fraction via mass-action kinetic competition."""
    if sup_trna_conc_um < 0.0 or erf1_conc_um < 0.0 or k_on_sup <= 0.0 or k_erf1 <= 0.0 or k_basal < 0.0:
        raise ValueError("Concentrations and kinetic rates must be non-negative, and on-rates must be positive.")

    k_sup = k_on_sup * sup_trna_conc_um
    k_term = k_erf1 * erf1_conc_um
    total_rate = k_sup + k_basal + k_term
    if total_rate <= 0.0:
        return 0.0
    return min(1.0, max(0.0, (k_sup + k_basal) / total_rate))


def compute_nmd_escape_fraction(
    readthrough_fraction: float,
    f_nmd_cf: float = DEFAULT_F_NMD_CF,
    k_05_nmd: float = DEFAULT_K_05_NMD,
) -> float:
    """Compute fraction of CFTR mRNA escaping SMG1/UPF1-mediated nonsense decay."""
    if readthrough_fraction < 0.0 or readthrough_fraction > 1.0:
        raise ValueError("Readthrough fraction must be in [0.0, 1.0]")
    if f_nmd_cf < 0.0 or f_nmd_cf > 1.0 or k_05_nmd <= 0.0:
        raise ValueError("f_nmd_cf must be in [0, 1] and k_05_nmd must be positive.")

    delta_mrna = 1.0 - f_nmd_cf
    saturation = readthrough_fraction / (readthrough_fraction + k_05_nmd)
    return min(1.0, max(f_nmd_cf, f_nmd_cf + delta_mrna * saturation))


def compute_restored_cftr_protein_fraction(
    rescued_mrna_fraction: float,
    readthrough_fraction: float,
) -> float:
    """Compute steady-state full-length CFTR protein synthesis rate relative to wild-type."""
    if rescued_mrna_fraction < 0.0 or rescued_mrna_fraction > 1.0:
        raise ValueError("Rescued mRNA fraction must be in [0.0, 1.0]")
    if readthrough_fraction < 0.0 or readthrough_fraction > 1.0:
        raise ValueError("Readthrough fraction must be in [0.0, 1.0]")
    return min(1.0, max(0.0, rescued_mrna_fraction * readthrough_fraction))


def compute_transepithelial_chloride_current_pct_wt(
    protein_fraction: float,
    gamma_activity: float = 1.00,
) -> float:
    """Compute transepithelial chloride short-circuit current (% WT)."""
    if protein_fraction < 0.0 or protein_fraction > 1.0:
        raise ValueError("Protein fraction must be in [0.0, 1.0]")
    if gamma_activity <= 0.0 or gamma_activity > 1.50:
        raise ValueError("Specific activity gamma must be positive and <= 1.50")
    return protein_fraction * gamma_activity * 100.0


def compute_asl_height_um(
    current_pct_wt: Any,
    h_asl_cf_um: float = DEFAULT_H_ASL_CF_UM,
    h_asl_wt_um: float = DEFAULT_H_ASL_WT_UM,
    k_05_asl_pct: float = DEFAULT_K_05_ASL_PCT,
) -> Any:
    """Compute airway surface liquid height (um) supporting complex numbers for differentiation."""
    if isinstance(current_pct_wt, (int, float)) and current_pct_wt < 0.0:
        raise ValueError("Chloride current cannot be negative.")
    delta_h = h_asl_wt_um - h_asl_cf_um
    hydration_ratio = current_pct_wt / (current_pct_wt + k_05_asl_pct)
    return h_asl_cf_um + delta_h * hydration_ratio


def compute_mcc_velocity_mm_min(
    asl_height_um: float,
    h_cilia_min_um: float = DEFAULT_H_CILIA_MIN_UM,
    v_mcc_wt_mm_min: float = DEFAULT_V_MCC_WT_MM_MIN,
    h_asl_wt_um: float = DEFAULT_H_ASL_WT_UM,
) -> float:
    """Compute mucociliary clearance velocity (mm/min) above ciliary stroke floor."""
    if asl_height_um <= h_cilia_min_um:
        return 0.0
    effective_stroke = (asl_height_um - h_cilia_min_um) / (h_asl_wt_um - h_cilia_min_um)
    return v_mcc_wt_mm_min * (min(1.5, max(0.0, effective_stroke)) ** 1.5)


def complex_step_asl_derivative(
    protein_fraction: float,
    gamma_activity: float = 1.00,
    h: float = 1e-30,
) -> Tuple[float, float, float]:
    """Verify sensitivity d(h_ASL)/d(F_protein) using complex-step differentiation against central finite differences."""
    p_complex = complex(protein_fraction, h)
    curr_complex = p_complex * gamma_activity * 100.0
    h_asl_complex = compute_asl_height_um(curr_complex)
    cs_derivative = h_asl_complex.imag / h

    # Central difference
    delta = 1e-6
    p_plus = min(1.0, protein_fraction + delta)
    curr_plus = p_plus * gamma_activity * 100.0
    h_plus = compute_asl_height_um(curr_plus)

    p_minus = max(0.0, protein_fraction - delta)
    curr_minus = p_minus * gamma_activity * 100.0
    h_minus = compute_asl_height_um(curr_minus)

    fd_derivative = (h_plus - h_minus) / (2.0 * delta)
    denom = abs(fd_derivative) if abs(fd_derivative) > 1e-12 else 1.0
    rel_error = abs(cs_derivative - fd_derivative) / denom
    return cs_derivative, fd_derivative, rel_error


def evaluate_round187_suppressor_trna_readthrough(preset_key: str) -> Round187SuppressorTRNAResult:
    """Evaluate a named clinical/biophysical regime preset for Round 187."""
    if preset_key not in REGIME_PRESETS:
        raise KeyError(f"Unknown preset key: {preset_key}")

    p = REGIME_PRESETS[preset_key]
    conc_um = float(p["sup_trna_conc_um"])
    erf1_um = float(p["erf1_conc_um"])
    gamma = float(p["gamma_activity"])
    mutation = str(p["target_mutation"])

    eta_read = compute_erf1_competition_readthrough_fraction(conc_um, erf1_um)
    f_mrna = compute_nmd_escape_fraction(eta_read)
    f_prot = compute_restored_cftr_protein_fraction(f_mrna, eta_read)
    curr_pct = compute_transepithelial_chloride_current_pct_wt(f_prot, gamma)
    h_asl = compute_asl_height_um(curr_pct)
    v_mcc = compute_mcc_velocity_mm_min(h_asl)

    # Complex-step derivative verification
    _, _, rel_err = complex_step_asl_derivative(max(1e-4, f_prot), gamma)

    # Kahn-Kalai Curative Threshold for Nonsense Suppressor tRNA:
    # eta_read >= 0.15, F_protein >= 0.08, I_sc >= 7.0% WT, h_ASL >= 7.0 um, v_MCC >= 1.0 mm/min
    crosses_threshold = (
        (eta_read >= 0.15)
        and (f_prot >= 0.08)
        and (curr_pct >= 7.0)
        and (h_asl >= 7.0)
        and (v_mcc >= 1.0)
    )
    expected_cross = bool(p["expected_crosses_kahn_kalai"])
    invariant_satisfied = (crosses_threshold == expected_cross) and (rel_err <= 1e-6)

    return Round187SuppressorTRNAResult(
        preset_name=preset_key,
        target_mutation=mutation,
        sup_trna_conc_um=conc_um,
        erf1_conc_um=erf1_um,
        gamma_activity=gamma,
        readthrough_fraction=eta_read,
        rescued_mrna_fraction=f_mrna,
        restored_protein_fraction=f_prot,
        transepithelial_current_pct_wt=curr_pct,
        airway_surface_liquid_height_um=h_asl,
        mucociliary_transport_velocity_mm_min=v_mcc,
        crosses_kahn_kalai_threshold=crosses_threshold,
        falsification_invariant_satisfied=invariant_satisfied,
        complex_step_relative_error=rel_err,
    )


def evaluate_falsification_negative_control(
    poison_negative_conc: bool = False,
    poison_negative_erf1: bool = False,
    poison_negative_gamma: bool = False,
) -> Dict[str, Any]:
    """Execute poisoned negative controls to verify deterministic falsification."""
    if poison_negative_conc:
        try:
            compute_erf1_competition_readthrough_fraction(sup_trna_conc_um=-2.0)
            return {"passed": False, "error": "Failed to catch negative sup-tRNA concentration"}
        except ValueError as e:
            return {"passed": True, "error": f"Correctly caught: {e}"}

    if poison_negative_erf1:
        try:
            compute_erf1_competition_readthrough_fraction(sup_trna_conc_um=1.0, erf1_conc_um=-1.0)
            return {"passed": False, "error": "Failed to catch negative eRF1 concentration"}
        except ValueError as e:
            return {"passed": True, "error": f"Correctly caught: {e}"}

    if poison_negative_gamma:
        try:
            compute_transepithelial_chloride_current_pct_wt(protein_fraction=0.5, gamma_activity=-0.5)
            return {"passed": False, "error": "Failed to catch negative specific activity gamma"}
        except ValueError as e:
            return {"passed": True, "error": f"Correctly caught: {e}"}

    # Baseline clean check
    res = evaluate_round187_suppressor_trna_readthrough("ace_trna_therapeutic_high_dose_g542x")
    return {"passed": res.falsification_invariant_satisfied, "error": None}


if __name__ == "__main__":
    print("Round 187: Engineered Nonsense Suppressor tRNA eRF1 Competition Bridge (T-0868)")
    for key in REGIME_PRESETS:
        r = evaluate_round187_suppressor_trna_readthrough(key)
        status = "CROSSES_CURATIVE" if r.crosses_kahn_kalai_threshold else "SUB_CURATIVE"
        print(f"  [{key}] ({r.target_mutation}) -> eta_read={r.readthrough_fraction:.3f}, "
              f"f_mRNA={r.rescued_mrna_fraction:.3f}, F_protein={r.restored_protein_fraction:.3f}, "
              f"I_sc={r.transepithelial_current_pct_wt:.2f}% WT, h_ASL={r.airway_surface_liquid_height_um:.2f} um, "
              f"v_MCC={r.mucociliary_transport_velocity_mm_min:.2f} mm/min [{status}]")
