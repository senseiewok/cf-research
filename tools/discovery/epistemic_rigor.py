#!/usr/bin/env python3
"""Epistemic Rigor, Parameter Provenance (AH-9), Dimensional Bounds, and Parametric Sensitivity (AH-11).

Implements the Wave C Scientific Grounding & Epistemic Questioning upgrades from the
Whole-Workspace 5-Chair Council Review:
  - AH-9: BiophysicalParameter dataclass distinguishing 'empirical', 'derived', and
    'phenomenological_heuristic' parameters with explicit source location and uncertainty bounds.
  - AH-11: Parametric sensitivity & rank-order stability harness (+/-25% perturbation sweep
    on phenomenological heuristics to prove qualitative ranking invariance).
  - Dimensional & physical boundary validator (probabilities in [0, 1], non-negative
    conductances/diffusivities, Kelvin temperatures > 0, single-AAV <= 4,700 bp).
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Any, Callable, Dict, List, Literal, Mapping, Sequence, Tuple

ParameterCategory = Literal["empirical", "derived", "phenomenological_heuristic"]
VALID_CATEGORIES: Tuple[ParameterCategory, ...] = (
    "empirical",
    "derived",
    "phenomenological_heuristic",
)


@dataclass(frozen=True)
class BiophysicalParameter:
    """Explicitly tagged biophysical parameter enforcing AH-9 provenance rules."""

    name: str
    value: float
    units: str
    category: ParameterCategory
    doi_or_pdb: str
    location_in_source: str
    uncertainty_range: Tuple[float, float] | None = None
    notes: str = ""

    def __post_init__(self) -> None:
        if not self.name or not self.name.strip():
            raise ValueError("BiophysicalParameter.name must be non-empty")
        if not math.isfinite(self.value):
            raise ValueError(f"Parameter '{self.name}' value must be finite, got {self.value}")
        if not self.units or not self.units.strip():
            raise ValueError(f"Parameter '{self.name}' must declare physical units (use 'dimensionless' if unitless)")
        if self.category not in VALID_CATEGORIES:
            raise ValueError(
                f"Parameter '{self.name}' has invalid category '{self.category}'; must be one of {VALID_CATEGORIES}"
            )
        if self.category in ("empirical", "derived"):
            if not self.doi_or_pdb or not self.doi_or_pdb.strip():
                raise ValueError(
                    f"AH-9 violation: '{self.name}' ({self.category}) must cite a primary DOI, PDB, NCT, or physical constant source"
                )
            if not self.location_in_source or not self.location_in_source.strip():
                raise ValueError(
                    f"AH-9 violation: '{self.name}' ({self.category}) must specify location_in_source (e.g. Table/Figure/Equation)"
                )
        if self.uncertainty_range is not None:
            low, high = self.uncertainty_range
            if not (math.isfinite(low) and math.isfinite(high)) or low > high:
                raise ValueError(
                    f"Parameter '{self.name}' has invalid uncertainty_range ({low}, {high}): expected finite low <= high"
                )
            if not (low - 1e-12 <= self.value <= high + 1e-12):
                raise ValueError(
                    f"Parameter '{self.name}' value {self.value} falls outside uncertainty_range [{low}, {high}]"
                )

    def perturbed(self, pct: float) -> Tuple[float, float]:
        """Return (low_value, high_value) under +/- pct fractional perturbation."""
        if self.uncertainty_range is not None and pct <= 0.0:
            return self.uncertainty_range
        delta = abs(self.value) * abs(pct)
        if delta == 0.0:
            delta = abs(pct)
        return (self.value - delta, self.value + delta)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def validate_physical_bounds(
    name: str,
    value: float,
    *,
    kind: Literal[
        "probability",
        "percentage",
        "non_negative",
        "positive_kelvin",
        "aav_single_capsid_bp",
        "free_energy_kcal_mol",
    ],
) -> float:
    """Validate physical and thermodynamic invariants on a scalar quantity."""
    if not math.isfinite(value):
        raise ValueError(f"Physical bound violation: '{name}' is non-finite ({value})")
    if kind == "probability" and not (0.0 - 1e-9 <= value <= 1.0 + 1e-9):
        raise ValueError(f"Physical bound violation: probability '{name}' = {value} outside [0, 1]")
    elif kind == "percentage" and not (0.0 - 1e-9 <= value <= 100.0 + 1e-9):
        raise ValueError(f"Physical bound violation: percentage '{name}' = {value}% outside [0, 100]")
    elif kind == "non_negative" and value < -1e-12:
        raise ValueError(f"Physical bound violation: non-negative quantity '{name}' = {value} < 0")
    elif kind == "positive_kelvin" and value <= 0.0:
        raise ValueError(f"Physical bound violation: Kelvin temperature '{name}' = {value} K <= 0")
    elif kind == "aav_single_capsid_bp" and (value <= 0 or value > 4700):
        raise ValueError(
            f"Physical bound violation: single-AAV genome length '{name}' = {value} bp exceeds 4,700 bp packaging ceiling"
        )
    elif kind == "free_energy_kcal_mol" and not (-250.0 <= value <= 250.0):
        raise ValueError(
            f"Physical bound violation: physiological free energy '{name}' = {value} kcal/mol outside [-250, +250]"
        )
    return value


CORE_BIOPHYSICAL_PARAMETERS: Dict[str, BiophysicalParameter] = {
    "PHYSIOLOGICAL_TEMP_K": BiophysicalParameter(
        name="PHYSIOLOGICAL_TEMP_K",
        value=310.15,
        units="K",
        category="empirical",
        doi_or_pdb="CODATA-2022 / Standard Airway Physiology",
        location_in_source="37.0 C = 310.15 K",
        uncertainty_range=(308.15, 311.15),
    ),
    "RT_37C_KCAL_MOL": BiophysicalParameter(
        name="RT_37C_KCAL_MOL",
        value=0.001987204 * 310.15,
        units="kcal/mol",
        category="derived",
        doi_or_pdb="CODATA-2022",
        location_in_source="R * T = 1.987204e-3 kcal/(mol*K) * 310.15 K = 0.61633 kcal/mol",
        uncertainty_range=(0.612, 0.619),
    ),
    "WT_CFTR_SINGLE_CHANNEL_CONDUCTANCE_PS": BiophysicalParameter(
        name="WT_CFTR_SINGLE_CHANNEL_CONDUCTANCE_PS",
        value=8.4,
        units="pS",
        category="empirical",
        doi_or_pdb="10.1126/science.ade2216",
        location_in_source="Electrophysiology patch-clamp single-channel Cl- conductance (8.4-8.5 pS)",
        uncertainty_range=(7.8, 9.2),
    ),
    "NBD1_WT_TM_CELSIUS": BiophysicalParameter(
        name="NBD1_WT_TM_CELSIUS",
        value=49.5,
        units="degC",
        category="empirical",
        doi_or_pdb="PDB:5UAK | 10.1126/science.ade2216",
        location_in_source="Differential scanning calorimetry (DSC) WT human NBD1 melting midpoint",
        uncertainty_range=(48.0, 51.0),
    ),
    "NBD1_F508DEL_DELTA_TM_CELSIUS": BiophysicalParameter(
        name="NBD1_F508DEL_DELTA_TM_CELSIUS",
        value=-7.5,
        units="degC",
        category="empirical",
        doi_or_pdb="PDB:5UAK | 10.1126/science.ade2216",
        location_in_source="F508del thermal destabilization relative to WT NBD1 (-6.5 to -8.5 C)",
        uncertainty_range=(-8.5, -6.5),
    ),
    "CF_SPUTUM_MUCIN_MESH_PORE_NM": BiophysicalParameter(
        name="CF_SPUTUM_MUCIN_MESH_PORE_NM",
        value=110.0,
        units="nm",
        category="empirical",
        doi_or_pdb="10.1038/s41565-020-0669-6",
        location_in_source="CF purulent sputum mucin pore size distribution (90-140 nm; severe 62 nm at 7.5% solids)",
        uncertainty_range=(62.0, 140.0),
    ),
    "AAV_SINGLE_CAPSID_PACKAGING_LIMIT_BP": BiophysicalParameter(
        name="AAV_SINGLE_CAPSID_PACKAGING_LIMIT_BP",
        value=4700.0,
        units="bp",
        category="empirical",
        doi_or_pdb="10.1073/pnas.261714599 | NCT05248230",
        location_in_source="Wild-type AAV genome packaging ceiling including 2x145 bp ITRs",
        uncertainty_range=(4650.0, 4750.0),
    ),
    "CURATIVE_ISC_WT_THRESHOLD_PCT": BiophysicalParameter(
        name="CURATIVE_ISC_WT_THRESHOLD_PCT",
        value=10.0,
        units="%WT",
        category="empirical",
        doi_or_pdb="10.1016/j.jcf.2021.08.014 | 10.1038/s41586-023-06133-1",
        location_in_source="HBE Ussing chamber transepithelial Cl- current benchmark associated with mild/pancreatic-sufficient CF",
        uncertainty_range=(8.0, 15.0),
    ),
    "MUCUS_ELECTROSTATIC_PENALTY_COEFF": BiophysicalParameter(
        name="MUCUS_ELECTROSTATIC_PENALTY_COEFF",
        value=0.045,
        units="mV^-2",
        category="phenomenological_heuristic",
        doi_or_pdb="10.1038/s41565-020-0669-6",
        location_in_source="Phenomenological Gaussian electrostatic binding penalty fit to mucin fibre tracking",
        uncertainty_range=(0.03375, 0.05625),
        notes="Subject to AH-11 +/-25% parametric sensitivity verification.",
    ),
    "ERAD_COMPETITIVE_INHIBITION_ALPHA": BiophysicalParameter(
        name="ERAD_COMPETITIVE_INHIBITION_ALPHA",
        value=0.18,
        units="dimensionless",
        category="phenomenological_heuristic",
        doi_or_pdb="10.1021/acs.biochem.5c00590",
        location_in_source="Effective competitive chaperone sequestration coupling constant for W1282X CFTR_1281",
        uncertainty_range=(0.135, 0.225),
        notes="Subject to AH-11 +/-25% parametric sensitivity verification.",
    ),
    "COMPOSITE_SAFETY_EFFICACY_WEIGHT_RATIO": BiophysicalParameter(
        name="COMPOSITE_SAFETY_EFFICACY_WEIGHT_RATIO",
        value=1.5,
        units="dimensionless",
        category="phenomenological_heuristic",
        doi_or_pdb="cf-research/tools/discovery/alpha_evolve_pareto_regimen.py",
        location_in_source="Multi-objective scalarization weight balancing %WT Cl- restoration vs DILI/TLR9 risk",
        uncertainty_range=(1.125, 1.875),
        notes="Subject to AH-11 +/-25% parametric sensitivity verification.",
    ),
}


def run_parametric_sensitivity(
    eval_candidates_fn: Callable[[Mapping[str, float]], Mapping[str, float]],
    parameters: Sequence[BiophysicalParameter],
    *,
    perturbation_pct: float = 0.25,
    only_heuristics: bool = True,
) -> Dict[str, Any]:
    """Execute AH-11 +/-25% parametric sensitivity & rank-order stability verification.

    Args:
        eval_candidates_fn: Function mapping {param_name: float} -> {candidate_id: score_float}
            (higher score = better candidate).
        parameters: Sequence of BiophysicalParameter instances to perturb.
        perturbation_pct: Fractional perturbation (default 0.25 for +/-25%).
        only_heuristics: If True, only perturbs parameters with category == 'phenomenological_heuristic'.

    Returns:
        Dictionary with baseline ranking, per-parameter tornado sensitivities, max relative drift %,
        and whether the top candidate and full qualitative rank order are preserved.
    """
    if perturbation_pct <= 0.0 or perturbation_pct >= 1.0:
        raise ValueError(f"perturbation_pct must be in (0, 1), got {perturbation_pct}")

    base_map: Dict[str, float] = {p.name: p.value for p in parameters}
    base_scores = dict(eval_candidates_fn(base_map))
    if not base_scores:
        raise ValueError("eval_candidates_fn returned empty candidate score map")

    base_ranking = [
        k for k, _ in sorted(base_scores.items(), key=lambda kv: (-kv[1], kv[0]))
    ]
    base_top = base_ranking[0]

    target_params = [
        p for p in parameters
        if (not only_heuristics or p.category == "phenomenological_heuristic")
    ]

    tornado: List[Dict[str, Any]] = []
    rank_flips: List[Dict[str, Any]] = []
    top_flips: List[Dict[str, Any]] = []
    max_rel_drift_pct = 0.0

    for p in target_params:
        low_val, high_val = p.perturbed(perturbation_pct)
        param_max_drift = 0.0
        for label, trial_val in (("minus_pct", low_val), ("plus_pct", high_val)):
            trial_map = dict(base_map)
            trial_map[p.name] = trial_val
            trial_scores = dict(eval_candidates_fn(trial_map))
            trial_ranking = [
                k for k, _ in sorted(trial_scores.items(), key=lambda kv: (-kv[1], kv[0]))
            ]
            for cand_id, b_val in base_scores.items():
                t_val = trial_scores[cand_id]
                denom = max(abs(b_val), 1e-9)
                rel_pct = abs(t_val - b_val) / denom * 100.0
                if rel_pct > param_max_drift:
                    param_max_drift = rel_pct
                if rel_pct > max_rel_drift_pct:
                    max_rel_drift_pct = rel_pct

            if trial_ranking != base_ranking:
                rank_flips.append({
                    "parameter": p.name,
                    "direction": label,
                    "trial_value": round(trial_val, 6),
                    "baseline_ranking": base_ranking,
                    "perturbed_ranking": trial_ranking,
                })
            if trial_ranking[0] != base_top:
                top_flips.append({
                    "parameter": p.name,
                    "direction": label,
                    "trial_value": round(trial_val, 6),
                    "baseline_top": base_top,
                    "perturbed_top": trial_ranking[0],
                })

        tornado.append({
            "parameter": p.name,
            "category": p.category,
            "baseline_value": p.value,
            "low_value": round(low_val, 6),
            "high_value": round(high_val, 6),
            "max_relative_score_drift_pct": round(param_max_drift, 3),
        })

    tornado.sort(key=lambda item: -item["max_relative_score_drift_pct"])

    return {
        "rule": "AH-11",
        "perturbation_pct": perturbation_pct * 100.0,
        "parameters_tested": len(target_params),
        "baseline_ranking": base_ranking,
        "baseline_scores": {k: round(v, 4) for k, v in base_scores.items()},
        "top_candidate_preserved": len(top_flips) == 0,
        "full_ranking_preserved": len(rank_flips) == 0,
        "max_relative_drift_pct": round(max_rel_drift_pct, 3),
        "tornado_sensitivities": tornado,
        "rank_flips": rank_flips,
        "top_flips": top_flips,
        "status": "PASS" if len(top_flips) == 0 else "WARN_RANK_SENSITIVE",
    }
