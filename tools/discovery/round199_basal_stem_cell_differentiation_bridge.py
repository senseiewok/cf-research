#!/usr/bin/env python3
"""
Round 199: Basal Stem Cell Notch/Wnt Signaling & Multiciliated Cell Lineage Differentiation Invariant Bridge (T-0888).

Mathematical & Physical Pillar (Column A - Deductive Theorem):
  1. Airway Basal Stem Cell Plasticity & Lineage Bifurcation:
     Human airway basal stem cells (BCs, marked by p63 and cytokeratin 5/14, KRT5/KRT14) reside on the basal lamina of the
     pseudostratified respiratory epithelium, constituting the self-renewing progenitor pool that regenerates secretory,
     intermediate, and multiciliated cell lineages (Rock et al., PNAS 2009; Mou et al., Cell Stem Cell 2016;
     Purkayastha et al., Cell Stem Cell 2020, DOI: 10.1016/j.stem.2020.11.010; Bulcaen et al., Nat Biotechnol 2020,
     DOI: 10.1038/s41587-020-0565-5).
     Lineage commitment between secretory (goblet / club) and multiciliated cell fates is orchestrated by an antagonistic
     signaling switch governed by Notch and canonical Wnt / beta-catenin pathways:
       - Active Notch signaling (Notch_ON) triggers HES1/HEY1 transcriptional cascades that commit progenitors toward secretory
         fates while repressing the master multiciliated transcription factor multicilin (MCIDAS) and downstream FOXJ1.
       - Attenuation or pharmacological inhibition of Notch signaling (Notch_OFF, e.g. via gamma-secretase inhibitors DAPT/DBZ)
         paired with canonical Wnt activation (e.g. via GSK3b inhibition using CHIR99021 or Wnt-3a ligands) releases MCIDAS
         inhibition, driving massive centriole amplification and multiciliated cell differentiation.

  2. Hill-Type Cooperative Lineage Commitment Probability:
     Multiciliated lineage commitment probability P_ciliated is formulated as a coupled cooperative dual-input Hill function:
       P_ciliated = ( W^n_w / (K_w^n_w + W^n_w) ) * ( K_n^n_n / (K_n^n_n + N^n_n) )
     where:
       W in [0, 1] is the normalized canonical Wnt / beta-catenin transcriptional activation level,
       N in [0, 1] is the normalized Notch pathway receptor cleavage and nuclear signaling intensity,
       n_w = 2.0 and K_w = 0.35 represent Wnt cooperativity and half-activation threshold,
       n_n = 2.0 and K_n = 0.30 represent Notch repressive cooperativity and half-inhibition threshold.

  3. Surface Cellular Coverage Fractions:
     Given competent basal progenitor fraction f_basal_competent in [0, 1] and terminal differentiation efficiency eta_diff in [0, 1]:
       f_ciliated = f_basal_competent * P_ciliated * eta_diff
       f_secretory = f_basal_competent * ( N^n_n / (K_n^n_n + N^n_n) ) * eta_diff
       f_basal_reserve = 1.0 - (f_ciliated + f_secretory)
     with physical conservation law:
       f_ciliated + f_secretory + f_basal_reserve = 1.0.

  4. Coordinated Mucociliary Transport Velocity & Metachronal Wave Emergence:
     Multiciliated cell clearance requires a critical surface coverage percolation threshold (f_thresh approx 0.20) below which
     cilia beat discoordinately without collective hydrodynamic coupling. Above f_thresh, coordinated metachronal wave
     propagation emerges, restoring directional mucus transport velocity V_mcc (mm/min):
       V_mcc = V_max * ( max(0.0, f_ciliated - f_thresh)^gamma / ( K_mcc^gamma + max(0.0, f_ciliated - f_thresh)^gamma ) )
     where V_max = 10.0 mm/min (normal human tracheobronchial clearance ceiling), f_thresh = 0.20, K_mcc = 0.35, and gamma = 2.0.

  5. Kahn-Kalai Curative Threshold for Airway Regeneration:
     - Multiciliated lineage commitment P_ciliated >= 0.60 (>= 60.0% of committing progenitors).
     - Multiciliated surface coverage fraction f_ciliated >= 0.45 (>= 45.0% surface ciliation).
     - Coordinated mucociliary transport velocity V_mcc >= 4.0 mm/min (>= 40% of physiological maximum).
     - Secretory-to-ciliated metaplasia ratio f_secretory / f_ciliated <= 0.80.

  6. Invariants & Deterministic Bounds:
     - Strict probability and fraction bounds: 0.0 <= P_ciliated, f_ciliated, f_secretory, f_basal_reserve <= 1.0.
     - Transport velocity bounds: 0.0 <= V_mcc <= V_max (10.0 mm/min).
     - Monotonicity: d(P_ciliated)/d(W) >= 0; d(P_ciliated)/d(N) <= 0; d(V_mcc)/d(f_ciliated) >= 0.
     - Complex-step differentiation: d(V_mcc)/d(W) verified at step h = 1e-30 against central finite differences
       with relative error <= 1e-3.
     - Poisoned negative controls: negative Wnt, negative Notch, or negative basal competence fail deterministically with ValueError.

CF Biophysical Domain Mapping (Column B - Popperian Falsifiable Hypothesis):
  - In High-Wnt / Notch-Inhibited Curative Regeneration Regimen:
    W = 0.85, N = 0.10, f_basal_competent = 0.90, eta_diff = 0.85:
    P_ciliated = 0.772 >= 0.60, f_ciliated = 0.591 >= 0.45, f_secretory = 0.076,
    ratio = 0.129 <= 0.80, V_mcc = 5.54 mm/min >= 4.0 mm/min,
    robustly crossing the Kahn-Kalai curative threshold for functional mucociliary restoration.
  - In Maximum-Wnt / Notch-Ablated Super-Curative Regimen:
    W = 0.95, N = 0.05, f_basal_competent = 0.95, eta_diff = 0.90:
    P_ciliated = 0.856, f_ciliated = 0.732, f_secretory = 0.023, V_mcc = 6.98 mm/min,
    achieving near-complete physiological airway surface rematuration.
  - In Untreated CF Inflammatory Hyper-Notch Goblet Metaplasia Regimen:
    W = 0.20, N = 0.85, f_basal_competent = 0.80, eta_diff = 0.80:
    P_ciliated = 0.027, f_ciliated = 0.017 < 0.20, f_secretory = 0.569,
    ratio = 33.5 >> 0.80, V_mcc = 0.0 mm/min, failing clearance threshold completely.
  - In Sub-Curative Partial Differentiation Regimen:
    W = 0.45, N = 0.35, f_basal_competent = 0.85, eta_diff = 0.80:
    P_ciliated = 0.262, f_ciliated = 0.178 < 0.45, V_mcc = 0.0 mm/min (below percolation threshold),
    failing coordinated transport threshold.
  - In Wnt-Deficient Quiescent Reserve Regimen:
    W = 0.05, N = 0.15, f_basal_competent = 0.75, eta_diff = 0.75:
    P_ciliated = 0.016, f_ciliated = 0.009, f_basal_reserve = 0.887, V_mcc = 0.0 mm/min,
    basal cells remain dormant without productive multiciliated differentiation.

Research simulation only; not medical advice.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Dict, Tuple

PRIMARY_CITATIONS: Tuple[str, ...] = (
    "DOI: 10.1016/j.stem.2020.11.010",
    "DOI: 10.1038/s41587-020-0565-5",
)

# Biophysical Differentiation & Clearance Constants
DEFAULT_K_W: float = 0.35                # Half-activation Wnt signaling threshold
DEFAULT_N_W: float = 2.0                 # Wnt pathway Hill cooperativity coefficient
DEFAULT_K_N: float = 0.30                # Half-repression Notch signaling threshold
DEFAULT_N_N: float = 2.0                 # Notch pathway Hill cooperativity coefficient
DEFAULT_BASAL_COMPETENCE: float = 0.90   # Fraction of basal progenitor pool competent for differentiation
DEFAULT_ETA_DIFF: float = 0.85           # Terminal differentiation yield efficiency
DEFAULT_V_MAX_MCC_MM_MIN: float = 10.0   # Maximum physiological tracheobronchial clearance velocity (mm/min)
DEFAULT_F_THRESH_CILIATED: float = 0.20  # Minimum ciliated coverage required for metachronal wave propagation
DEFAULT_K_MCC_FRACTION: float = 0.35     # Half-saturation ciliated coverage fraction for MCC velocity
DEFAULT_GAMMA_MCC: float = 2.0           # Metachronal hydrodynamic coupling exponent

PARAMETER_PROVENANCE: Dict[str, Dict[str, object]] = {
    "wnt_half_activation_threshold": {
        "value": 0.35,
        "tag": "phenomenological_heuristic",
        "source": "Analytical canonical Wnt/beta-catenin transcriptional half-activation constant K_w (exact 0.35 scalar not stated in Purkayastha et al., Cell Stem Cell 2020, DOI: 10.1016/j.stem.2020.11.010)",
        "bounds": (0.15, 0.60),
    },
    "notch_half_repression_threshold": {
        "value": 0.30,
        "tag": "phenomenological_heuristic",
        "source": "Analytical Notch signaling repression threshold K_n governing multicilin/FOXJ1 silencing (exact 0.30 scalar not stated in Purkayastha et al., Cell Stem Cell 2020, DOI: 10.1016/j.stem.2020.11.010)",
        "bounds": (0.10, 0.50),
    },
    "basal_stem_competence_fraction": {
        "value": 0.90,
        "tag": "phenomenological_heuristic",
        "source": "Analytical viable basal progenitor differentiation competence fraction in human airway cultures (exact 0.90 constant not stated in Bulcaen et al., Nat Biotechnol 2020, DOI: 10.1038/s41587-020-0565-5)",
        "bounds": (0.60, 0.98),
    },
    "mcc_metachronal_threshold": {
        "value": 0.20,
        "tag": "phenomenological_heuristic",
        "source": "Analytical critical percolation threshold f_thresh of surface ciliation for coordinated mucociliary wave clearance (exact 0.20 threshold not stated in Purkayastha et al., Cell Stem Cell 2020, DOI: 10.1016/j.stem.2020.11.010)",
        "bounds": (0.10, 0.30),
    },
    "mcc_max_velocity_mm_min": {
        "value": 10.0,
        "tag": "empirical",
        "source": "Physiological human tracheobronchial mucociliary clearance velocity ceiling (4.0-10.0 mm/min in healthy human airways; Purkayastha et al., Cell Stem Cell 2020, DOI: 10.1016/j.stem.2020.11.010)",
        "bounds": (4.0, 15.0),
    },
}

TWO_COLUMN_EPISTEMIC_HYGIENE: Dict[str, Dict[str, str]] = {
    "deductive_column_a": {
        "formal_theorem": (
            "Given dual antagonistic Wnt-activation and Notch-inhibition inputs on airway basal stem cells, the joint "
            "differentiation probability P_ciliated is the product of cooperative Hill activations and repressions: "
            "P_ciliated = (W^n_w / (K_w^n_w + W^n_w)) * (K_n^n_n / (K_n^n_n + N^n_n)). Surface ciliated density "
            "scales proportionally, and mucociliary transport velocity V_mcc emerges as a non-linear percolation "
            "order parameter above the critical ciliary density threshold f_thresh."
        ),
        "invariants": (
            "0 <= P_ciliated <= 1; 0 <= f_ciliated, f_secretory, f_basal_reserve <= 1; "
            "f_ciliated + f_secretory + f_basal_reserve == 1.0; 0 <= V_mcc <= V_max; "
            "d(P_ciliated)/d(W) >= 0; d(P_ciliated)/d(N) <= 0; d(V_mcc)/d(f_ciliated) >= 0."
        ),
    },
    "popperian_column_b": {
        "falsifiable_hypothesis": (
            "Pharmacological Notch inhibition (N <= 0.15) combined with canonical Wnt activation (W >= 0.70) in human "
            "CF airway basal stem cell cultures restores multiciliated cell surface fraction to >= 45% and rescues "
            "coordinated mucociliary clearance transport velocity to >= 4.0 mm/min."
        ),
        "empirical_falsification_threshold": (
            "The hypothesis is falsified if Notch inhibition (N = 0.10) and Wnt activation (W = 0.85) yield < 45% "
            "ciliated surface coverage, produce a secretory/ciliated ratio > 0.80, or fail to achieve V_mcc >= 4.0 mm/min."
        ),
    },
}

FIVE_LAYER_LADDER: Dict[str, str] = {
    "cool_human_fact": (
        "Human airways are lined by millions of microscopic hair-like cilia beating 15 times a second in continuous "
        "waves that sweep mucus and trapped bacteria up and out of the lungs. In cystic fibrosis, chronic inflammation "
        "tricks airway stem cells into making sticky mucus-secreting goblet cells instead of multiciliated cells, "
        "paralyzing the natural lung cleaning escalators."
    ),
    "everyday_analogy": (
        "Think of a stadium crowd performing 'the wave.' If only 1 out of 10 people stand up, no collective wave can travel. "
        "Once more than half the people join in, the wave sweeps effortlessly around the entire stadium. Airway cilia need a "
        "critical team density (at least 45%) to coordinate and propel mucus out of the airways."
    ),
    "engineering_lens": (
        "From a developmental dynamical systems viewpoint, the basal stem cell fate decision operates as a bistable genetic "
        "toggle switch. Notch signaling acts as a repressor of the multicilin-FOXJ1 ciliary program while driving secretory "
        "fate. Suppressing Notch while driving Wnt switches the basin of attraction entirely toward multiciliated morphogenesis."
    ),
    "formal_theorem": (
        "Let stem cell fate vector X = [f_ciliated, f_secretory, f_reserve]^T satisfy |X|_1 = 1 on the unit simplex Delta^2. "
        "Multiciliated commitment probability P_c(W, N) = H_W(W) * H_N(N) satisfies dP_c/dW > 0 and dP_c/dN < 0 everywhere "
        "on (0, 1)^2. Macroscopic transport velocity V_mcc(f_c) is a continuous piecewise-smooth Hill order parameter with "
        "zero transport for f_c <= f_thresh and positive monotonic curvature d^2V_mcc/df_c^2 > 0 near percolation onset."
    ),
    "falsifiable_code_and_toy": (
        "evaluate_round199_basal_stem_cell_differentiation(regime_key) computes exact commitment probabilities, surface "
        "cellular coverage fractions, and mucociliary transport velocities across 5 defined regimes with complex-step verification."
    ),
}

REGIME_PRESETS: Dict[str, Dict[str, float]] = {
    "wnt_high_notch_low_curative_regeneration": {
        "wnt_activity": 0.85,
        "notch_activity": 0.10,
        "f_basal_competent": 0.90,
        "eta_diff": 0.85,
    },
    "wnt_max_notch_inhibited_super_curative": {
        "wnt_activity": 0.95,
        "notch_activity": 0.05,
        "f_basal_competent": 0.95,
        "eta_diff": 0.90,
    },
    "cf_inflammatory_notch_hyperactive_metaplasia": {
        "wnt_activity": 0.20,
        "notch_activity": 0.85,
        "f_basal_competent": 0.80,
        "eta_diff": 0.80,
    },
    "sub_curative_partial_differentiation": {
        "wnt_activity": 0.45,
        "notch_activity": 0.35,
        "f_basal_competent": 0.85,
        "eta_diff": 0.80,
    },
    "wnt_deficient_dormant_basal_state": {
        "wnt_activity": 0.05,
        "notch_activity": 0.15,
        "f_basal_competent": 0.75,
        "eta_diff": 0.75,
    },
}


@dataclass(frozen=True)
class Round199BasalStemCellResult:
    regime_key: str
    wnt_activity: float
    notch_activity: float
    p_ciliated: float
    f_ciliated: float
    f_secretory: float
    f_basal_reserve: float
    secretory_to_ciliated_ratio: float
    v_mcc_mm_min: float
    crosses_kahn_kalai_threshold: bool
    falsification_invariant_satisfied: bool
    complex_step_relative_error: float


def compute_multiciliated_commitment_probability(
    wnt_activity: float | complex,
    notch_activity: float | complex,
    k_w: float = DEFAULT_K_W,
    n_w: float = DEFAULT_N_W,
    k_n: float = DEFAULT_K_N,
    n_n: float = DEFAULT_N_N,
) -> float | complex:
    """Computes cooperative multiciliated lineage commitment probability from Wnt and Notch inputs."""
    if isinstance(wnt_activity, (int, float)) and (wnt_activity < 0.0 or wnt_activity > 1.0):
        raise ValueError(f"Wnt activity must be in [0, 1], got {wnt_activity}")
    if isinstance(notch_activity, (int, float)) and (notch_activity < 0.0 or notch_activity > 1.0):
        raise ValueError(f"Notch activity must be in [0, 1], got {notch_activity}")
    if k_w <= 0.0 or k_n <= 0.0:
        raise ValueError(f"Half-saturation thresholds must be positive (k_w={k_w}, k_n={k_n})")

    wnt_term = (wnt_activity ** n_w) / ((k_w ** n_w) + (wnt_activity ** n_w))
    notch_repression = (k_n ** n_n) / ((k_n ** n_n) + (notch_activity ** n_n))
    return wnt_term * notch_repression


def compute_lineage_coverage_fractions(
    p_ciliated: float | complex,
    notch_activity: float | complex,
    f_basal_competent: float = DEFAULT_BASAL_COMPETENCE,
    eta_diff: float = DEFAULT_ETA_DIFF,
    k_n: float = DEFAULT_K_N,
    n_n: float = DEFAULT_N_N,
) -> Tuple[float | complex, float | complex, float | complex]:
    """Computes surface ciliated, secretory, and reserve basal cellular coverage fractions."""
    if isinstance(p_ciliated, (int, float)) and (p_ciliated < 0.0 or p_ciliated > 1.0):
        raise ValueError(f"Commitment probability must be in [0, 1], got {p_ciliated}")
    if f_basal_competent < 0.0 or f_basal_competent > 1.0:
        raise ValueError(f"Basal competence must be in [0, 1], got {f_basal_competent}")
    if eta_diff < 0.0 or eta_diff > 1.0:
        raise ValueError(f"Differentiation efficiency must be in [0, 1], got {eta_diff}")

    f_ciliated = f_basal_competent * p_ciliated * eta_diff
    notch_activation = (notch_activity ** n_n) / ((k_n ** n_n) + (notch_activity ** n_n))
    f_secretory = f_basal_competent * notch_activation * eta_diff
    f_basal_reserve = 1.0 - (f_ciliated + f_secretory)
    return f_ciliated, f_secretory, f_basal_reserve


def compute_mucociliary_transport_velocity(
    f_ciliated: float | complex,
    v_max: float = DEFAULT_V_MAX_MCC_MM_MIN,
    f_thresh: float = DEFAULT_F_THRESH_CILIATED,
    k_mcc: float = DEFAULT_K_MCC_FRACTION,
    gamma: float = DEFAULT_GAMMA_MCC,
) -> float | complex:
    """Computes macroscopic coordinated mucociliary transport velocity (mm/min)."""
    if isinstance(f_ciliated, (int, float)) and f_ciliated < 0.0:
        raise ValueError(f"Ciliated coverage fraction must be non-negative, got {f_ciliated}")
    if v_max <= 0.0:
        raise ValueError(f"Maximum velocity must be positive, got {v_max}")

    if isinstance(f_ciliated, complex):
        # Smooth activation for complex-step autodiff
        excess = f_ciliated - f_thresh
        smooth_excess = 0.5 * (excess + (excess ** 2 + 1e-16) ** 0.5)
        num = smooth_excess ** gamma
        den = (k_mcc ** gamma) + num
        return v_max * (num / den)
    else:
        if f_ciliated <= f_thresh:
            return 0.0
        excess_real = f_ciliated - f_thresh
        num_real = excess_real ** gamma
        den_real = (k_mcc ** gamma) + num_real
        return v_max * (num_real / den_real)


def complex_step_mcc_derivative(
    wnt_activity: float,
    notch_activity: float,
    f_basal_competent: float = DEFAULT_BASAL_COMPETENCE,
    eta_diff: float = DEFAULT_ETA_DIFF,
    step_h: float = 1e-30,
) -> Tuple[float, float, float]:
    """Computes d(V_mcc)/d(Wnt) via complex-step differentiation vs central finite difference."""
    # Complex step evaluation
    z_wnt = complex(wnt_activity, step_h)
    p_c_complex = compute_multiciliated_commitment_probability(z_wnt, notch_activity)
    f_c_complex, _, _ = compute_lineage_coverage_fractions(
        p_c_complex, notch_activity, f_basal_competent, eta_diff
    )
    v_complex = compute_mucociliary_transport_velocity(f_c_complex)
    cs_deriv = v_complex.imag / step_h

    # Finite difference evaluation
    fd_step = 1e-4
    def _eval_v(w: float) -> float:
        p_c = float(compute_multiciliated_commitment_probability(w, notch_activity))
        f_c, _, _ = compute_lineage_coverage_fractions(
            p_c, notch_activity, f_basal_competent, eta_diff
        )
        return float(compute_mucociliary_transport_velocity(f_c))

    w_hi = min(1.0, wnt_activity + fd_step)
    w_lo = max(0.0, wnt_activity - fd_step)
    fd_deriv = (_eval_v(w_hi) - _eval_v(w_lo)) / (w_hi - w_lo)
    rel_err = abs(cs_deriv - fd_deriv) / max(abs(cs_deriv), 1e-9)
    return cs_deriv, fd_deriv, rel_err


def evaluate_round199_basal_stem_cell_differentiation(regime_key: str) -> Round199BasalStemCellResult:
    """Evaluates the Round 199 basal stem cell Notch/Wnt differentiation model for a preset regime."""
    if regime_key not in REGIME_PRESETS:
        raise KeyError(f"Unknown regime '{regime_key}'. Valid regimes: {list(REGIME_PRESETS.keys())}")

    cfg = REGIME_PRESETS[regime_key]
    wnt = cfg["wnt_activity"]
    notch = cfg["notch_activity"]
    f_competent = cfg.get("f_basal_competent", DEFAULT_BASAL_COMPETENCE)
    eta = cfg.get("eta_diff", DEFAULT_ETA_DIFF)

    p_ciliated = float(compute_multiciliated_commitment_probability(wnt, notch))
    f_ciliated, f_secretory, f_reserve = compute_lineage_coverage_fractions(
        p_ciliated, notch, f_competent, eta
    )
    f_ciliated = float(f_ciliated)
    f_secretory = float(f_secretory)
    f_reserve = float(f_reserve)

    sec_to_cil_ratio = f_secretory / f_ciliated if f_ciliated > 0.0 else float("inf")
    v_mcc = float(compute_mucociliary_transport_velocity(f_ciliated))

    # Complex step derivative check (evaluated at operative or midpoint Wnt)
    eval_wnt = max(0.20, min(0.85, wnt))
    _, _, rel_err = complex_step_mcc_derivative(eval_wnt, notch, f_competent, eta)

    # Invariants verification
    inv_satisfied = (
        0.0 <= p_ciliated <= 1.0
        and 0.0 <= f_ciliated <= 1.0
        and 0.0 <= f_secretory <= 1.0
        and -1e-9 <= f_reserve <= 1.0 + 1e-9
        and abs((f_ciliated + f_secretory + f_reserve) - 1.0) <= 1e-6
        and 0.0 <= v_mcc <= DEFAULT_V_MAX_MCC_MM_MIN
        and rel_err <= 1e-3
    )

    # Kahn-Kalai curative threshold:
    # 1. P_ciliated >= 0.60
    # 2. f_ciliated >= 0.45
    # 3. V_mcc >= 4.0 mm/min
    # 4. sec_to_cil_ratio <= 0.80
    crosses_threshold = (
        p_ciliated >= 0.60
        and f_ciliated >= 0.45
        and v_mcc >= 4.0
        and sec_to_cil_ratio <= 0.80
    )

    return Round199BasalStemCellResult(
        regime_key=regime_key,
        wnt_activity=wnt,
        notch_activity=notch,
        p_ciliated=p_ciliated,
        f_ciliated=f_ciliated,
        f_secretory=f_secretory,
        f_basal_reserve=f_reserve,
        secretory_to_ciliated_ratio=sec_to_cil_ratio,
        v_mcc_mm_min=v_mcc,
        crosses_kahn_kalai_threshold=crosses_threshold,
        falsification_invariant_satisfied=inv_satisfied,
        complex_step_relative_error=rel_err,
    )


def evaluate_falsification_negative_control(
    poison_negative_wnt: bool = False,
    poison_negative_notch: bool = False,
    poison_negative_competence: bool = False,
) -> Dict[str, Any]:
    """Tests poisoned negative controls that must deterministically raise ValueError."""
    if poison_negative_wnt:
        try:
            compute_multiciliated_commitment_probability(-0.5, 0.2)
            return {"passed": False, "reason": "Negative Wnt activity accepted unexpectedly"}
        except ValueError:
            return {"passed": True, "reason": "Deterministically rejected negative Wnt activity"}

    if poison_negative_notch:
        try:
            compute_multiciliated_commitment_probability(0.5, -0.2)
            return {"passed": False, "reason": "Negative Notch activity accepted unexpectedly"}
        except ValueError:
            return {"passed": True, "reason": "Deterministically rejected negative Notch activity"}

    if poison_negative_competence:
        try:
            compute_lineage_coverage_fractions(0.5, 0.2, f_basal_competent=-0.1)
            return {"passed": False, "reason": "Negative basal competence accepted unexpectedly"}
        except ValueError:
            return {"passed": True, "reason": "Deterministically rejected negative basal competence"}

    return {"passed": True, "reason": "Baseline control pass"}
