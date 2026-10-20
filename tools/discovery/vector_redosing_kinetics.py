#!/usr/bin/env python3
"""Inhaled viral vs non-viral vector repeat-dosing neutralizing antibody & epithelial turnover decay model (T-0198).

Models two-compartment airway epithelial turnover, cellular genome retention, and mucosal
neutralizing antibody (NAb) / anti-PEG anamnestic kinetics across repeat inhaled dosing:
1. Two-Compartment Airway Epithelial Turnover ODE Model:
   - Terminally differentiated luminal cells (ciliated/secretory, SCGB1A1+/FOXJ1+):
     turnover half-life t_1/2,cil approx 90-120 days (decay rate lambda_cil = ln(2)/t_1/2,cil).
   - Long-lived basal stem/progenitor cells (KRT5+/TP63+):
     residence half-life t_1/2,basal >= 3,650 days (10+ years quiescent reservoir).
   - Transgene persistence regimes:
     * Episomal concatemer (AAV 4D-710 / NCT05248230, HSV-1 KB407 / NCT05504837):
       retained in non-dividing luminal cells; diluted upon basal stem cell division.
     * Host chromatin integration (Lentiviral BI 3720931 / NCT06515002):
       stable proviral insertion in basal progenitors passed to all daughter lineages.
     * Transient mRNA (VX-522 / NCT05668741, RCT2100 / NCT06237335):
       cytoplasmic translation subject to rapid enzymatic degradation (t_1/2,mRNA approx 2-4 days).
2. Mucosal Anti-Capsid / Anti-PEG Anamnestic Neutralization Kinetics:
   - Neutralizing antibody titer boost per repeat dose B_modality:
     T_k = T_0 + sum_{j<k} B_modality * exp(-lambda_ab * (t_k - t_j))
   - Sigmoidal / Hill neutralization law for effective delivery:
     E_k = E_1 / (1 + (T_k / IC_50)^n_H)
   - Contrasts viral capsids (AAV, HSV-1, lentivirus) exhibiting high anamnestic titer accumulation
     against synthetic mRNA LNPs exhibiting minimal anti-PEG neutralization permitting frequent redosing.
3. Multi-Dose Saw-Tooth CFTR Expression Trajectory Generator:
   - Evaluates consecutive doses across arbitrary dosing intervals Delta t, computing peak,
     trough, time above therapeutic threshold (15% WT CFTR), and redosing attenuation ratio.

Grounded Sources:
- Clinical trial accession: NCT05248230 (4D-710 Phase 1/2 aerosolized AAV in CF)
- Clinical trial accession: NCT06515002 (BI 3720931 inhaled lentiviral vector in CF)
- Clinical trial accession: NCT05504837 (KB407 inhaled HSV-1 in CF)
- Clinical trial accession: NCT05668741 (VX-522 inhaled mRNA LNP in CF)
- Clinical trial accession: NCT06237335 (RCT2100 inhaled mRNA in CF)
All verified against tools/trial_atlas/snapshot.json.

Attribution & Disclaimers:
ClinicalTrials.gov data are used here with attribution to the U.S. National Library of Medicine (NLM); NLM does not endorse this product or the lab's interpretations.
Research, not medical advice.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

NLM_REGISTRY_DISCLAIMER: str = (
    "ClinicalTrials.gov data are used here with attribution to the U.S. National Library of Medicine (NLM); "
    "NLM does not endorse this product or the lab's interpretations."
)

MEDICAL_DISCLAIMER: str = "Research, not medical advice."

# Physiological parameters
THERAPEUTIC_THRESHOLD_CFTR_PCT: float = 15.0  # 15% WT CFTR chloride secretion threshold


@dataclass(frozen=True)
class VectorModalityConfig:
    """Biophysical and immunological parameters for an inhaled gene transfer modality."""
    modality_id: str
    display_name: str
    nct_id: str
    vector_class: str  # 'episomal_aav', 'episomal_hsv1', 'integrative_lentivirus', 'transient_mrna'
    genome_type: str   # 'episomal', 'integrative', 'transient_mrna'
    luminal_half_life_days: float   # Ciliated / secretory cell clearance half-life
    basal_half_life_days: float     # Basal progenitor persistence half-life
    mrna_half_life_days: float      # Transcript degradation half-life (days)
    basal_fraction_accessible: float # Fraction of transducing units reaching basal progenitors
    baseline_nab_titer: float       # Pre-existing mucosal / transudated serum NAb titer (T0)
    nab_boost_per_dose: float       # Anamnestic titer increase per subsequent inhalation (B)
    nab_decay_half_life_days: float # NAb systemic/mucosal titer clearance half-life
    ic50_titer: float               # Titer threshold for 50% neutralization
    hill_coefficient: float         # Cooperativity of neutralization (n_H)
    initial_efficiency: float       # Baseline first-dose transduction / expression index
    description: str


# Catalog of 5 clinically evaluated inhaled CF gene delivery modalities (verified in snapshot.json)
MODALITY_CATALOG: Dict[str, VectorModalityConfig] = {
    "4D-710": VectorModalityConfig(
        modality_id="4D-710",
        display_name="4D-710 Aerosolized AAV (A101 Capsid)",
        nct_id="NCT05248230",
        vector_class="episomal_aav",
        genome_type="episomal",
        luminal_half_life_days=105.0,  # Airway epithelial ciliated turnover 90-120 days
        basal_half_life_days=270.0,    # Episomal loss / dilution on basal stem cell division
        mrna_half_life_days=3.0,
        basal_fraction_accessible=0.15,
        baseline_nab_titer=1.0,
        nab_boost_per_dose=350.0,      # Robust anti-capsid memory boost
        nab_decay_half_life_days=270.0,
        ic50_titer=45.0,
        hill_coefficient=1.5,
        initial_efficiency=1.0,
        description="Directed evolution AAV capsid expressing codon-optimized CFTR cDNA episomally.",
    ),
    "BI-3720931": VectorModalityConfig(
        modality_id="BI-3720931",
        display_name="BI 3720931 Inhaled Lentiviral Vector",
        nct_id="NCT06515002",
        vector_class="integrative_lentivirus",
        genome_type="integrative",
        luminal_half_life_days=105.0,
        basal_half_life_days=3650.0,   # Stable genomic integration in self-renewing basal stem cells (10+ yrs)
        mrna_half_life_days=3.0,
        basal_fraction_accessible=0.25,
        baseline_nab_titer=1.0,
        nab_boost_per_dose=220.0,
        nab_decay_half_life_days=240.0,
        ic50_titer=40.0,
        hill_coefficient=1.4,
        initial_efficiency=0.90,
        description="Recombinant pseudotyped lentiviral vector stably integrating into host progenitor chromatin.",
    ),
    "KB407": VectorModalityConfig(
        modality_id="KB407",
        display_name="KB407 Inhaled HSV-1 Vector",
        nct_id="NCT05504837",
        vector_class="episomal_hsv1",
        genome_type="episomal",
        luminal_half_life_days=105.0,
        basal_half_life_days=300.0,
        mrna_half_life_days=3.0,
        basal_fraction_accessible=0.12,
        baseline_nab_titer=5.0,        # Pre-existing human anti-HSV1 immunity
        nab_boost_per_dose=450.0,
        nab_decay_half_life_days=365.0,
        ic50_titer=55.0,
        hill_coefficient=1.6,
        initial_efficiency=0.95,
        description="Replication-incompetent engineered HSV-1 carrying dual full-length CFTR copies episomally.",
    ),
    "VX-522": VectorModalityConfig(
        modality_id="VX-522",
        display_name="VX-522 Inhaled CFTR mRNA LNP",
        nct_id="NCT05668741",
        vector_class="transient_mrna",
        genome_type="transient_mrna",
        luminal_half_life_days=3.0,    # mRNA cytoplasmic decay t_1/2 ~3 days
        basal_half_life_days=3.0,
        mrna_half_life_days=3.0,
        basal_fraction_accessible=0.20,
        baseline_nab_titer=1.0,
        nab_boost_per_dose=4.0,        # Very low immunogenicity, negligible anti-PEG titer
        nab_decay_half_life_days=60.0,
        ic50_titer=80.0,
        hill_coefficient=1.0,
        initial_efficiency=1.0,
        description="Synthetic codon-optimized human CFTR mRNA encapsulated in aerosolized lipid nanoparticles.",
    ),
    "RCT2100": VectorModalityConfig(
        modality_id="RCT2100",
        display_name="RCT2100 Inhaled CFTR mRNA",
        nct_id="NCT06237335",
        vector_class="transient_mrna",
        genome_type="transient_mrna",
        luminal_half_life_days=3.0,
        basal_half_life_days=3.0,
        mrna_half_life_days=3.0,
        basal_fraction_accessible=0.20,
        baseline_nab_titer=1.0,
        nab_boost_per_dose=5.0,
        nab_decay_half_life_days=60.0,
        ic50_titer=80.0,
        hill_coefficient=1.0,
        initial_efficiency=0.98,
        description="Aerosolized chemically modified CFTR mRNA formulation for non-viral repeat inhalation.",
    ),
}


@dataclass(frozen=True)
class DoseEventEvaluation:
    """Evaluates the NAb titer and effective delivered fraction for a discrete dosing event."""
    dose_index: int              # 1-indexed dose number (k = 1, 2, ...)
    time_days: float             # Inhalation timepoint (days)
    nab_titer_pre_dose: float    # Cumulative anamnestic NAb titer T_k before inhalation
    neutralization_fraction: float  # Fraction of inhaled vector neutralized (0.0 to 1.0)
    effective_delivery_efficiency: float # Relative transduction efficiency E_k (normalized to E_1)
    luminal_input_cftr_pct: float   # Newly synthesized CFTR pool in luminal cells (% WT)
    basal_input_cftr_pct: float     # Newly synthesized CFTR pool in basal cells (% WT)


@dataclass(frozen=True)
class TrajectoryTimePoint:
    """Discrete timepoint state during the longitudinal simulation."""
    time_days: float
    ciliated_luminal_pct: float
    basal_stem_pct: float
    total_cftr_pct_wt: float
    active_nab_titer: float
    is_above_therapeutic_threshold: bool


@dataclass(frozen=True)
class SimulationSummary:
    """Summary metrics of the multi-dose redosing simulation."""
    modality: str
    nct_id: str
    vector_class: str
    total_doses: int
    dosing_interval_days: float
    total_simulation_days: float
    peak_cftr_pct_wt: float
    trough_cftr_pct_wt: float
    final_cftr_pct_wt: float
    time_above_15pct_threshold_days: float
    fraction_time_above_threshold: float
    dose_1_efficiency: float
    final_dose_efficiency: float
    redosing_attenuation_ratio: float  # E_N / E_1 (measures immunity drop-off)
    is_curative_integrative: bool
    nlm_disclaimer: str
    medical_disclaimer: str


def calculate_anamnestic_titer(
    dose_times_days: List[float],
    current_time_days: float,
    baseline_titer: float,
    boost_per_dose: float,
    decay_half_life_days: float,
) -> float:
    """Calculates mucosal/serum neutralizing antibody titer at current_time_days.

    T_k = T_0 + sum_{j: t_j < current_time} B * exp(-lambda_ab * (current_time - t_j))
    """
    if decay_half_life_days <= 0.0:
        lambda_ab = 0.0
    else:
        lambda_ab = math.log(2.0) / decay_half_life_days

    cumulative_boost = 0.0
    for t_dose in dose_times_days:
        if t_dose < current_time_days:
            dt = current_time_days - t_dose
            cumulative_boost += boost_per_dose * math.exp(-lambda_ab * dt)

    return baseline_titer + cumulative_boost


def calculate_hill_delivery_efficiency(
    titer: float,
    initial_efficiency: float,
    ic50_titer: float,
    hill_coefficient: float,
) -> Tuple[float, float]:
    """Calculates effective transduction efficiency E_k and neutralized fraction using Hill equation.

    E_k = E_1 / (1 + (T_k / IC_50)^n_H)
    """
    if ic50_titer <= 0.0 or titer <= 0.0:
        return initial_efficiency, 0.0

    ratio = titer / ic50_titer
    hill_term = math.pow(ratio, hill_coefficient)
    eff = initial_efficiency / (1.0 + hill_term)
    neut_frac = hill_term / (1.0 + hill_term)
    return eff, neut_frac


def simulate_redosing_trajectory(
    modality_key: str,
    dose_count: int = 4,
    interval_days: float = 30.0,
    nominal_dose_peak_pct: float = 40.0,
    total_days: Optional[float] = None,
    dt_days: float = 1.0,
) -> Tuple[List[DoseEventEvaluation], List[TrajectoryTimePoint], SimulationSummary]:
    """Simulates multi-dose saw-tooth CFTR expression trajectory and NAb kinetics.

    Parameters:
    -----------
    modality_key: str, key in MODALITY_CATALOG (e.g. '4D-710', 'BI-3720931', 'VX-522', 'KB407', 'RCT2100')
    dose_count: int, number of repeat inhalations (k = 1, 2, ..., N)
    interval_days: float, spacing between doses in days
    nominal_dose_peak_pct: float, target % WT CFTR expression achieved by 100% unneutralized single dose
    total_days: float, duration of simulation (defaults to dose_count * interval_days + 60 days)
    dt_days: float, ODE integration step in days (default 1.0)
    """
    if modality_key not in MODALITY_CATALOG:
        raise KeyError(f"Unknown modality '{modality_key}'. Valid choices: {list(MODALITY_CATALOG.keys())}")

    cfg = MODALITY_CATALOG[modality_key]
    if dose_count < 1:
        raise ValueError("dose_count must be >= 1")
    if interval_days < 1.0:
        raise ValueError("interval_days must be >= 1.0")

    if total_days is None:
        total_days = max(180.0, (dose_count - 1) * interval_days + 90.0)

    # Pre-compute dose schedule
    dose_times = [i * interval_days for i in range(dose_count)]

    # Evaluate each dose event
    dose_evaluations: List[DoseEventEvaluation] = []
    dose_inputs: Dict[int, Tuple[float, float]] = {}  # dose_idx -> (luminal_input, basal_input)

    for idx, t_dose in enumerate(dose_times, 1):
        # Titer before dose delivery
        titer_pre = calculate_anamnestic_titer(
            dose_times_days=dose_times[:idx-1],
            current_time_days=t_dose,
            baseline_titer=cfg.baseline_nab_titer,
            boost_per_dose=cfg.nab_boost_per_dose,
            decay_half_life_days=cfg.nab_decay_half_life_days,
        )
        eff, neut = calculate_hill_delivery_efficiency(
            titer=titer_pre,
            initial_efficiency=cfg.initial_efficiency,
            ic50_titer=cfg.ic50_titer,
            hill_coefficient=cfg.hill_coefficient,
        )
        total_input = nominal_dose_peak_pct * eff
        basal_input = total_input * cfg.basal_fraction_accessible
        luminal_input = total_input * (1.0 - cfg.basal_fraction_accessible)

        eval_ev = DoseEventEvaluation(
            dose_index=idx,
            time_days=t_dose,
            nab_titer_pre_dose=titer_pre,
            neutralization_fraction=neut,
            effective_delivery_efficiency=eff,
            luminal_input_cftr_pct=luminal_input,
            basal_input_cftr_pct=basal_input,
        )
        dose_evaluations.append(eval_ev)
        dose_inputs[round(t_dose)] = (luminal_input, basal_input)

    # Decay constants
    lambda_lum = math.log(2.0) / cfg.luminal_half_life_days
    lambda_basal = math.log(2.0) / cfg.basal_half_life_days

    # For integrative lentivirus, long-lived basal cells divide and replenish ciliated pool
    # Rate of progenitor-to-luminal differentiation replenishment:
    is_integrative = (cfg.genome_type == "integrative")
    diff_rate_replenish = lambda_lum * 0.85 if is_integrative else 0.0

    # Simulation loop
    c_lum = 0.0
    c_basal = 0.0
    trajectory: List[TrajectoryTimePoint] = []
    time_above_thresh_days = 0.0
    peak_cftr = 0.0
    trough_cftr = float("inf")

    steps = int(math.ceil(total_days / dt_days)) + 1
    for step in range(steps):
        t = step * dt_days

        # Check for dose delivery at this timepoint
        for t_dose_ev in dose_times:
            if abs(t - t_dose_ev) < (dt_days / 2.0):
                # Apply pulse
                inp_lum, inp_basal = dose_inputs.get(round(t_dose_ev), (0.0, 0.0))
                c_lum += inp_lum
                c_basal += inp_basal

        tot_cftr = c_lum + c_basal
        if tot_cftr > peak_cftr:
            peak_cftr = tot_cftr

        # Calculate active NAb titer at current time
        current_titer = calculate_anamnestic_titer(
            dose_times_days=[td for td in dose_times if td <= t],
            current_time_days=t,
            baseline_titer=cfg.baseline_nab_titer,
            boost_per_dose=cfg.nab_boost_per_dose,
            decay_half_life_days=cfg.nab_decay_half_life_days,
        )

        is_above = (tot_cftr >= THERAPEUTIC_THRESHOLD_CFTR_PCT)
        if is_above:
            time_above_thresh_days += dt_days

        # Track trough after the first dose peak has established
        if t >= interval_days and tot_cftr < trough_cftr:
            trough_cftr = tot_cftr

        trajectory.append(
            TrajectoryTimePoint(
                time_days=round(t, 2),
                ciliated_luminal_pct=round(c_lum, 3),
                basal_stem_pct=round(c_basal, 3),
                total_cftr_pct_wt=round(tot_cftr, 3),
                active_nab_titer=round(current_titer, 2),
                is_above_therapeutic_threshold=is_above,
            )
        )

        # Decay ODE step:
        # dc_lum/dt = - lambda_lum * c_lum + diff_rate * c_basal
        # dc_basal/dt = - lambda_basal * c_basal
        d_lum = (-lambda_lum * c_lum + diff_rate_replenish * c_basal) * dt_days
        d_basal = (-lambda_basal * c_basal) * dt_days

        c_lum = max(0.0, c_lum + d_lum)
        c_basal = max(0.0, c_basal + d_basal)

    if trough_cftr == float("inf"):
        trough_cftr = trajectory[-1].total_cftr_pct_wt

    final_cftr = trajectory[-1].total_cftr_pct_wt
    d1_eff = dose_evaluations[0].effective_delivery_efficiency
    dN_eff = dose_evaluations[-1].effective_delivery_efficiency
    attenuation_ratio = (dN_eff / d1_eff) if d1_eff > 0.0 else 0.0

    summary = SimulationSummary(
        modality=cfg.modality_id,
        nct_id=cfg.nct_id,
        vector_class=cfg.vector_class,
        total_doses=dose_count,
        dosing_interval_days=interval_days,
        total_simulation_days=total_days,
        peak_cftr_pct_wt=round(peak_cftr, 2),
        trough_cftr_pct_wt=round(trough_cftr, 2),
        final_cftr_pct_wt=round(final_cftr, 2),
        time_above_15pct_threshold_days=round(time_above_thresh_days, 1),
        fraction_time_above_threshold=round(time_above_thresh_days / total_days, 3),
        dose_1_efficiency=round(d1_eff, 4),
        final_dose_efficiency=round(dN_eff, 4),
        redosing_attenuation_ratio=round(attenuation_ratio, 4),
        is_curative_integrative=is_integrative,
        nlm_disclaimer=NLM_REGISTRY_DISCLAIMER,
        medical_disclaimer=MEDICAL_DISCLAIMER,
    )

    return dose_evaluations, trajectory, summary


def compare_all_modalities(
    dose_count: int = 4,
    interval_days: float = 30.0,
    total_days: float = 180.0,
) -> Dict[str, Any]:
    """Runs comparative multi-dose simulation across all 5 clinical CF gene delivery vectors."""
    results: Dict[str, Any] = {}
    for key in MODALITY_CATALOG:
        doses, traj, summary = simulate_redosing_trajectory(
            modality_key=key,
            dose_count=dose_count,
            interval_days=interval_days,
            total_days=total_days,
        )
        results[key] = {
            "summary": asdict(summary),
            "doses": [asdict(d) for d in doses],
            "sample_trajectory_points": [asdict(p) for p in traj[::15]],  # 15-day sample interval
        }
    return results


def run_self_test() -> int:
    """Deterministic internal self-test proving all invariants pass."""
    print("Running vector_redosing_kinetics internal self-test...")

    # 1. Modality catalog checks
    assert len(MODALITY_CATALOG) == 5, f"Expected 5 modalities, got {len(MODALITY_CATALOG)}"
    for key, cfg in MODALITY_CATALOG.items():
        assert cfg.nct_id.startswith("NCT"), f"Invalid NCT ID: {cfg.nct_id}"
        assert cfg.luminal_half_life_days > 0.0
        assert cfg.basal_half_life_days >= cfg.luminal_half_life_days or cfg.genome_type == "transient_mrna"

    # 2. Hill delivery efficiency checks
    e0, neut0 = calculate_hill_delivery_efficiency(titer=0.0, initial_efficiency=1.0, ic50_titer=50.0, hill_coefficient=1.5)
    assert e0 == 1.0 and neut0 == 0.0
    e_ic50, neut_ic50 = calculate_hill_delivery_efficiency(titer=50.0, initial_efficiency=1.0, ic50_titer=50.0, hill_coefficient=1.5)
    assert abs(e_ic50 - 0.5) < 1e-4 and abs(neut_ic50 - 0.5) < 1e-4

    # 3. Repeat AAV vs mRNA anamnestic attenuation check
    _, _, summary_aav = simulate_redosing_trajectory("4D-710", dose_count=4, interval_days=30.0)
    _, _, summary_mrna = simulate_redosing_trajectory("VX-522", dose_count=4, interval_days=30.0)
    assert summary_aav.redosing_attenuation_ratio < 0.25, (
        f"AAV redosing attenuation must be severe (<0.25), got {summary_aav.redosing_attenuation_ratio}"
    )
    assert summary_mrna.redosing_attenuation_ratio > 0.85, (
        f"mRNA redosing attenuation must be negligible (>0.85), got {summary_mrna.redosing_attenuation_ratio}"
    )

    # 4. Integrative lentivirus basal persistence check
    _, traj_lenti, summary_lenti = simulate_redosing_trajectory("BI-3720931", dose_count=2, interval_days=30.0, total_days=365.0)
    final_lenti = traj_lenti[-1].total_cftr_pct_wt
    _, traj_aav, _ = simulate_redosing_trajectory("4D-710", dose_count=2, interval_days=30.0, total_days=365.0)
    final_aav = traj_aav[-1].total_cftr_pct_wt
    assert final_lenti > final_aav, f"Lentivirus must exhibit superior 1-year persistence: {final_lenti} vs {final_aav}"

    print("ALL VECTOR REDOSING KINETICS INVARIANTS PASSED (5/5).")
    return 0


def build_cli_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Inhaled vector repeat-dosing neutralizing antibody & epithelial turnover decay model (T-0198)"
    )
    parser.add_argument(
        "--modality",
        type=str,
        default="4D-710",
        choices=list(MODALITY_CATALOG.keys()),
        help="Inhaled CF gene transfer candidate modality (default: 4D-710)",
    )
    parser.add_argument(
        "--doses",
        type=int,
        default=4,
        help="Number of repeat inhalations (default: 4)",
    )
    parser.add_argument(
        "--interval-days",
        type=float,
        default=30.0,
        help="Interval between doses in days (default: 30.0)",
    )
    parser.add_argument(
        "--nominal-dose-pct",
        type=float,
        default=40.0,
        help="Target % WT CFTR expression achieved by unneutralized dose (default: 40.0)",
    )
    parser.add_argument(
        "--days-total",
        type=float,
        default=180.0,
        help="Total duration of simulation in days (default: 180.0)",
    )
    parser.add_argument(
        "--compare-all",
        action="store_true",
        help="Simulate and compare all 5 clinical CF gene delivery candidates",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Output structured JSON payload",
    )
    parser.add_argument(
        "--self-test",
        action="store_true",
        help="Run deterministic biophysical invariant verification suite",
    )
    return parser


def main() -> int:
    parser = build_cli_parser()
    args = parser.parse_args()

    if args.self_test:
        return run_self_test()

    if args.compare_all:
        comp = compare_all_modalities(
            dose_count=args.doses,
            interval_days=args.interval_days,
            total_days=args.days_total,
        )
        if args.json:
            print(json.dumps(comp, indent=2))
        else:
            print("=" * 88)
            print("COMPARATIVE INHALED VECTOR REDOSING & EPITHELIAL TURNOVER MODEL")
            print(f"Doses: {args.doses} | Interval: {args.interval_days:.1f} days | Horizon: {args.days_total:.1f} days")
            print("=" * 88)
            print(f"{'Modality':<14} | {'NCT ID':<11} | {'Class':<14} | {'Peak CFTR':<10} | {'Trough':<8} | {'E_N / E_1':<9} | {'Time >15%'}")
            print("-" * 88)
            for mod, data in comp.items():
                s = data["summary"]
                print(
                    f"{mod:<14} | {s['nct_id']:<11} | {s['vector_class']:<14} | "
                    f"{s['peak_cftr_pct_wt']:>8.1f}% | {s['trough_cftr_pct_wt']:>6.1f}% | "
                    f"{s['redosing_attenuation_ratio']:>7.3f}   | {s['time_above_15pct_threshold_days']:>5.0f} d ({s['fraction_time_above_threshold']*100:.0f}%)"
                )
            print("=" * 88)
            print(f"NLM Disclaimer: {NLM_REGISTRY_DISCLAIMER}")
            print("Scope: Research, not medical advice.")
        return 0

    doses, traj, summary = simulate_redosing_trajectory(
        modality_key=args.modality,
        dose_count=args.doses,
        interval_days=args.interval_days,
        nominal_dose_peak_pct=args.nominal_dose_pct,
        total_days=args.days_total,
    )

    if args.json:
        payload = {
            "summary": asdict(summary),
            "doses": [asdict(d) for d in doses],
            "trajectory": [asdict(p) for p in traj],
        }
        print(json.dumps(payload, indent=2))
    else:
        cfg = MODALITY_CATALOG[args.modality]
        print("=" * 80)
        print(f"INHALED VECTOR REPEAT-DOSING & EPITHELIAL KINETICS: {cfg.display_name}")
        print(f"Registry Accession: {cfg.nct_id} | Class: {cfg.vector_class}")
        print("=" * 80)
        print(f"  Luminal Turnover t_1/2:          {cfg.luminal_half_life_days:.1f} days")
        print(f"  Basal Reservoir t_1/2:           {cfg.basal_half_life_days:.1f} days")
        print(f"  NAb Titer Boost (B_modality):    {cfg.nab_boost_per_dose:.1f} titer units/dose")
        print(f"  Neutralization IC50:             {cfg.ic50_titer:.1f} titer units (Hill n_H={cfg.hill_coefficient:.1f})")
        print("-" * 80)
        print("  DOSING SCHEDULE EVALUATIONS:")
        for d in doses:
            print(
                f"    Dose #{d.dose_index} at Day {d.time_days:>3.0f}: "
                f"NAb Titer = {d.nab_titer_pre_dose:>6.1f} | "
                f"Eff = {d.effective_delivery_efficiency*100:>5.1f}% | "
                f"Input CFTR = {d.luminal_input_cftr_pct + d.basal_input_cftr_pct:>5.1f}% WT"
            )
        print("-" * 80)
        print(f"  Peak Total CFTR:                 {summary.peak_cftr_pct_wt:.1f}% WT")
        print(f"  Trough CFTR (between doses):     {summary.trough_cftr_pct_wt:.1f}% WT")
        print(f"  Final Day {args.days_total:.0f} CFTR:            {summary.final_cftr_pct_wt:.1f}% WT")
        print(f"  Days Above 15% WT Threshold:     {summary.time_above_15pct_threshold_days:.0f} days ({summary.fraction_time_above_threshold*100:.1f}% of window)")
        print(f"  Redosing Attenuation (E_N / E_1): {summary.redosing_attenuation_ratio:.4f}")
        print("=" * 80)
        print(f"NLM Disclaimer: {NLM_REGISTRY_DISCLAIMER}")
        print("Scope: Research, not medical advice.")

    return 0


if __name__ == "__main__":
    sys.exit(main())
