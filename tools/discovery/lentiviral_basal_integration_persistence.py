#!/usr/bin/env python3
"""Inhaled Lentiviral Vector Epigenetic Silencing & Airway Basal Progenitor Integration Persistence Model (T-0407).

Quantitative biophysical model of lentiviral proviral integration in long-lived airway
basal stem cells (KRT5+/TP63+) and long-term promoter epigenetic silencing kinetics:
  1. Basal Progenitor Stem-Cell Retention & Clonal Transmission:
     Models proviral genomic integration retention across self-renewing basal stem cells
     (t_1/2 >= 3,650 days / ~10 years at homeostatic self-renewal) and clonal transmission
     to transit-amplifying (TA) progenitors and differentiated luminal daughter cells
     (ciliated FOXJ1+, secretory MUC5AC+/SCGB1A1+, and FOXI1+ pulmonary ionocytes; t_1/2 in [90, 120] days).
  2. CpG Methylation & Epigenetic Silencing Kinetics:
     Simulates DNA methyltransferase (DNMT)-mediated CpG island methylation and histone deacetylation
     attenuation of transgene expression over time across promoter designs (EF1a, MND, synthetic hCFTR).
  3. Proviral Integration Copy-Number Safety & Mutagenesis Balance:
     Evaluates vector copy number (VCN in [1, 4] per cell) balancing therapeutic apical CFTR
     channel density against insertional mutagenesis safety thresholds.
  4. Durable Transepithelial Chloride Secretion (% WT Delta I_sc):
     Couples persistent basal CFTR Band C (170 kDa complex-glycosylated, PDB 6MSM) maturation
     to steady-state apical chloride current, predicting multi-year duration above the curative
     threshold (>= 10% WT).
  5. Multi-Modality Comparative Longevity:
     Contrasts integrating lentiviral vectors (BI 3720931) against non-integrating episomal AAV
     (4D-710) and episodic aerosolized mRNA LNPs (VX-522).

Primary Evidence Grounding:
  - Discovery Brief: cf-research/landscape/mission-grounding-2026-10-24.md (Curative Discovery Council Track A Wave 18).
  - Structural Model: RCSB PDB 6MSM (Human CFTR phosphorylated open active state at 3.2 A).
  - Primary Literature:
    - Alton EW et al. (2015) Lancet Respir Med, DOI: 10.1016/S2213-2600(15)00245-3
    - Cooney AL et al. (2018) Mol Ther, DOI: 10.1016/j.ymthe.2018.04.014
    - Vijaykumar et al., JCI Insight 2026, DOI: 10.1172/jci.insight.206375
  - Clinical Registry Grounding:
    - NCT06515002 (Phase 1/2 trial of BI 3720931 inhaled lentiviral vector in CF)
    - NCT05668741 (Phase 1/2 aerosolized mRNA LNP VX-522 in CF)
    - NCT05248230 (Phase 1/2 aerosolized AAV capsid 4D-710 in CF)

Epistemic Status:
  - Answer Tier T2 (Open-access primary literature, RCSB PDB structures, ClinicalTrials.gov API v2).
  - Preclinical computational model; research, not medical advice. No clinical efficacy or therapeutic safety claim is made.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

SCRIPT_PATH = Path(__file__).resolve()

NLM_ATTRIBUTION = (
    "This product uses publicly available data from the U.S. National Library of Medicine (NLM), "
    "National Institutes of Health, Department of Health and Human Services; NLM is not responsible "
    "for the product and does not endorse or recommend this or any other product. "
    "ClinicalTrials.gov does not evaluate safety or efficacy; listing on the registry does not denote government endorsement."
)

RESEARCH_BOUNDARY = (
    "Preclinical in-silico lentiviral proviral integration & epigenetic persistence model; "
    "research, not medical advice. No clinical efficacy or therapeutic safety claim is made."
)

# Reference Grounding Identifiers
PRIMARY_DOI_ALTON = "10.1016/S2213-2600(15)00245-3"
PRIMARY_DOI_COONEY = "10.1016/j.ymthe.2018.04.014"
PRIMARY_DOI_VIJAYKUMAR = "10.1172/jci.insight.206375"
CLINICAL_TRIAL_BI3720931 = "NCT06515002"
CLINICAL_TRIAL_VX522 = "NCT05668741"
CLINICAL_TRIAL_4D710 = "NCT05248230"
STRUCTURAL_MODEL_PDB = "6MSM"

# Physiological and Kinetic Constants for Human Bronchial Epithelium
DEFAULT_BASAL_HALF_LIFE_DAYS = 3650.0          # Basal stem cell renewal half-life t_1/2 >= 3650 days (10 yr)
DEFAULT_LUMINAL_HALF_LIFE_DAYS = 100.0         # Differentiated luminal cell turnover t_1/2 in [90, 120] days
DEFAULT_TA_DIVISION_RATE_PER_DAY = 0.025       # Transit-amplifying progenitor division rate
DEFAULT_TA_MATURATION_RATE_PER_DAY = 0.020     # TA to differentiated luminal maturation rate
DEFAULT_FRACTION_BASAL = 0.30                  # Basal stem cell proportion in pseudostratified epithelium
DEFAULT_FRACTION_TA = 0.20                     # Transit-amplifying cell proportion
DEFAULT_FRACTION_DIFFERENTIATED = 0.50         # Differentiated luminal cell proportion

# Functional CFTR Threshold
CURATIVE_THRESHOLD_DELTA_ISC_PCT = 10.0        # >= 10% WT CFTR functional current threshold
CFTR_BAND_C_MATURATION_EFFICIENCY = 0.72       # Band B (140 kDa) -> Band C (170 kDa mature, PDB 6MSM)
SINGLE_PROVIRUS_ONCOGENIC_HIT_RATE = 1.5e-6    # Estimated insertional mutagenesis probability per integration
REGULATORY_MAX_INSERTIONAL_RISK = 1.0e-5       # Standard regulatory threshold (10 in a million)


@dataclass(frozen=True)
class PromoterParameters:
    """Parameters governing CpG island methylation and epigenetic silencing."""
    promoter_id: str
    display_name: str
    cpg_island_sites: int
    silencing_rate_per_day: float
    silencing_half_life_days: float
    expression_floor_fraction: float
    relative_transcription_strength: float
    insulator_protection: bool


PRESET_PROMOTERS: Dict[str, PromoterParameters] = {
    "EF1a": PromoterParameters(
        promoter_id="EF1a",
        display_name="Human Elongation Factor 1-alpha (EF1a)",
        cpg_island_sites=48,
        silencing_rate_per_day=0.00042,
        silencing_half_life_days=1650.0,
        expression_floor_fraction=0.15,
        relative_transcription_strength=1.00,
        insulator_protection=False,
    ),
    "MND": PromoterParameters(
        promoter_id="MND",
        display_name="Myeloproliferative Sarcoma Virus U3-deleted LTR (MND)",
        cpg_island_sites=32,
        silencing_rate_per_day=0.00075,
        silencing_half_life_days=924.0,
        expression_floor_fraction=0.08,
        relative_transcription_strength=1.35,
        insulator_protection=False,
    ),
    "synthetic_hCFTR": PromoterParameters(
        promoter_id="synthetic_hCFTR",
        display_name="Synthetic CpG-Free hCFTR Promoter with UCOE Insulator",
        cpg_island_sites=6,
        silencing_rate_per_day=0.00011,
        silencing_half_life_days=6300.0,
        expression_floor_fraction=0.45,
        relative_transcription_strength=0.85,
        insulator_protection=True,
    ),
}


@dataclass(frozen=True)
class EpithelialHierarchyState:
    """Cellular and transgene state in the pseudostratified airway epithelium."""
    day: float
    years: float
    basal_transduced_fraction: float
    ta_transduced_fraction: float
    luminal_transduced_fraction: float
    overall_corrected_fraction: float
    promoter_activity_fraction: float
    cftr_band_c_density_relative: float
    apical_chloride_current_pct_wt: float
    above_curative_threshold: bool


@dataclass(frozen=True)
class SimulationSummary:
    """Summary metrics of lentiviral persistence and safety simulation."""
    promoter_id: str
    vcn_per_cell: float
    initial_basal_transduction_pct: float
    peak_delta_isc_pct_wt: float
    final_delta_isc_pct_wt: float
    curative_duration_days: float
    curative_duration_years: float
    reaches_curative_threshold: bool
    insertional_risk_per_cell: float
    insertional_safety_verdict: str
    therapeutic_safety_index: float
    yearly_trajectory: List[EpithelialHierarchyState]


def compute_promoter_activity(
    day: float,
    promoter: PromoterParameters,
) -> float:
    """Compute promoter epigenetic expression factor at given time t (days).

    Simulates DNMT-driven progressive promoter CpG methylation and HDAC silencing:
      eta(t) = (1 - f_floor) * exp(-k_silence * t) + f_floor
    """
    if day <= 0.0:
        return 1.0
    floor = promoter.expression_floor_fraction
    k = promoter.silencing_rate_per_day
    attenuation = math.exp(-k * day)
    return (1.0 - floor) * attenuation + floor


def evaluate_insertional_mutagenesis_risk(
    vcn: float,
    hit_rate: float = SINGLE_PROVIRUS_ONCOGENIC_HIT_RATE,
) -> Tuple[float, str]:
    """Calculate cumulative insertional oncogenesis risk per transduced basal stem cell.

    P_hit = 1 - (1 - hit_rate)^VCN
    Assesses against regulatory safety threshold (1.0e-5).
    """
    if vcn <= 0.0:
        return 0.0, "PASS_ZERO_RISK"
    p_hit = 1.0 - math.pow(1.0 - hit_rate, vcn)
    verdict = "PASS_WITHIN_SAFETY_MARGIN" if p_hit <= REGULATORY_MAX_INSERTIONAL_RISK else "EXCEEDS_RECOMMENDED_RISK"
    return p_hit, verdict


def simulate_lentiviral_persistence(
    promoter_id: str = "synthetic_hCFTR",
    vcn: float = 2.0,
    initial_basal_transduction: float = 0.20,
    total_days: float = 3650.0,
    step_days: float = 10.0,
    basal_half_life_days: float = DEFAULT_BASAL_HALF_LIFE_DAYS,
    luminal_half_life_days: float = DEFAULT_LUMINAL_HALF_LIFE_DAYS,
) -> SimulationSummary:
    """Simulate multi-year lentiviral vector persistence in airway basal stem cells.

    Solves coupled lineage transmission and promoter silencing ODEs:
      dB_trans/dt = -lambda_B * B_trans (self-renewal turnover)
      dT_trans/dt = lambda_T * B_trans - gamma_T * T_trans
      dD_trans/dt = gamma_T * T_trans - mu_D * D_trans
    Coupled to CFTR Band C maturation and transepithelial chloride current (% WT Delta I_sc).
    """
    if promoter_id not in PRESET_PROMOTERS:
        raise ValueError(f"Unknown promoter '{promoter_id}'. Available: {list(PRESET_PROMOTERS.keys())}")
    if vcn < 0.5 or vcn > 8.0:
        raise ValueError(f"VCN {vcn} outside realistic range [0.5, 8.0]")

    promoter = PRESET_PROMOTERS[promoter_id]

    lambda_b = math.log(2.0) / max(100.0, basal_half_life_days)
    mu_d = math.log(2.0) / max(10.0, luminal_half_life_days)
    lambda_t = DEFAULT_TA_DIVISION_RATE_PER_DAY
    gamma_t = DEFAULT_TA_MATURATION_RATE_PER_DAY

    # Initial state
    b_trans = initial_basal_transduction
    t_trans = initial_basal_transduction * 0.50
    d_trans = initial_basal_transduction * 0.10

    trajectory: List[EpithelialHierarchyState] = []
    current_day = 0.0
    curative_duration = 0.0
    peak_current = 0.0

    steps = int(math.ceil(total_days / step_days))

    for step_idx in range(steps + 1):
        day = step_idx * step_days

        promoter_act = compute_promoter_activity(day, promoter)

        # Epithelial-wide weighted corrected fraction
        overall_fraction = (
            DEFAULT_FRACTION_BASAL * b_trans +
            DEFAULT_FRACTION_TA * t_trans +
            DEFAULT_FRACTION_DIFFERENTIATED * d_trans
        )

        # CFTR Band C mature channel density relative to wild-type
        # Band C density scales with VCN, promoter strength, promoter epigenetic activity, and corrected cells
        band_c_density = (
            vcn *
            promoter.relative_transcription_strength *
            promoter_act *
            CFTR_BAND_C_MATURATION_EFFICIENCY *
            overall_fraction
        )

        # Apical chloride current rescue (% WT Delta I_sc)
        # 1.0 WT Band C density delivers 100% WT Delta I_sc
        delta_isc_pct = band_c_density * 100.0 / CFTR_BAND_C_MATURATION_EFFICIENCY

        is_curative = delta_isc_pct >= CURATIVE_THRESHOLD_DELTA_ISC_PCT
        if is_curative:
            curative_duration = day
        if delta_isc_pct > peak_current:
            peak_current = delta_isc_pct

        # Record at annual intervals or specific steps
        if step_idx % int(max(1, 365.0 / step_days)) == 0 or step_idx == steps:
            trajectory.append(
                EpithelialHierarchyState(
                    day=round(day, 2),
                    years=round(day / 365.0, 2),
                    basal_transduced_fraction=round(b_trans, 4),
                    ta_transduced_fraction=round(t_trans, 4),
                    luminal_transduced_fraction=round(d_trans, 4),
                    overall_corrected_fraction=round(overall_fraction, 4),
                    promoter_activity_fraction=round(promoter_act, 4),
                    cftr_band_c_density_relative=round(band_c_density, 4),
                    apical_chloride_current_pct_wt=round(delta_isc_pct, 2),
                    above_curative_threshold=is_curative,
                )
            )

        # RK4 numerical integration of compartment dynamics
        def derivatives(b: float, t: float, d: float) -> Tuple[float, float, float]:
            db = -lambda_b * b
            dt = lambda_t * b - gamma_t * t
            dd = gamma_t * t - mu_d * d
            return db, dt, dd

        dt = step_days
        db1, dt1, dd1 = derivatives(b_trans, t_trans, d_trans)
        db2, dt2, dd2 = derivatives(b_trans + 0.5 * dt * db1, t_trans + 0.5 * dt * dt1, d_trans + 0.5 * dt * dd1)
        db3, dt3, dd3 = derivatives(b_trans + 0.5 * dt * db2, t_trans + 0.5 * dt * dt2, d_trans + 0.5 * dt * dd2)
        db4, dt4, dd4 = derivatives(b_trans + dt * db3, t_trans + dt * dt3, d_trans + dt * dd3)

        b_trans = max(0.0, b_trans + (dt / 6.0) * (db1 + 2 * db2 + 2 * db3 + db4))
        t_trans = max(0.0, t_trans + (dt / 6.0) * (dt1 + 2 * dt2 + 2 * dt3 + dt4))
        d_trans = max(0.0, d_trans + (dt / 6.0) * (dd1 + 2 * dd2 + 2 * dd3 + dd4))

    final_current = trajectory[-1].apical_chloride_current_pct_wt if trajectory else 0.0
    risk, verdict = evaluate_insertional_mutagenesis_risk(vcn)
    tsi = (peak_current / (risk * 1.0e5)) if risk > 0.0 else float("inf")

    return SimulationSummary(
        promoter_id=promoter_id,
        vcn_per_cell=vcn,
        initial_basal_transduction_pct=round(initial_basal_transduction * 100.0, 2),
        peak_delta_isc_pct_wt=round(peak_current, 2),
        final_delta_isc_pct_wt=round(final_current, 2),
        curative_duration_days=round(curative_duration, 1),
        curative_duration_years=round(curative_duration / 365.0, 2),
        reaches_curative_threshold=peak_current >= CURATIVE_THRESHOLD_DELTA_ISC_PCT,
        insertional_risk_per_cell=round(risk, 8),
        insertional_safety_verdict=verdict,
        therapeutic_safety_index=round(tsi, 2),
        yearly_trajectory=trajectory,
    )


def compare_gene_therapy_modalities(
    total_days: float = 3650.0,
) -> Dict[str, Any]:
    """Compare long-term durability of integrating lentivirus vs episomal AAV vs mRNA LNP.

    Benchmarks:
      - BI 3720931 Lentiviral Vector (integrating stem cell retention, t_1/2 >= 3,650 days)
      - 4D-710 Episomal AAV (diluted upon division, t_1/2 ~ 100 days)
      - VX-522 Aerosolized mRNA LNP (transient cytoplasmic decay, t_1/2 ~ 10 days)
    """
    days = [0, 30, 90, 180, 365, 730, 1460, 2190, 2920, 3650]
    days = [d for d in days if d <= total_days]

    # Inhaled Lentiviral model (synthetic hCFTR, VCN=2.0, 20% basal transduction)
    lenti_sim = simulate_lentiviral_persistence(
        promoter_id="synthetic_hCFTR",
        vcn=2.0,
        initial_basal_transduction=0.20,
        total_days=total_days,
    )
    lenti_map = {int(round(pt.day)): pt.apical_chloride_current_pct_wt for pt in lenti_sim.yearly_trajectory}

    # Episomal AAV model (4D-710, NCT05248230): initial expression 30% WT, mitotic dilution t_1/2 = 100 days
    aav_lambda = math.log(2.0) / 100.0
    aav_initial = 30.0

    # Inhaled mRNA LNP model (VX-522, NCT05668741): initial expression 40% WT, cytoplasmic decay t_1/2 = 10 days
    mrna_lambda = math.log(2.0) / 10.0
    mrna_initial = 40.0

    comparison: List[Dict[str, Any]] = []
    for d in days:
        lenti_val = lenti_map.get(d)
        if lenti_val is None:
            # Interpolate or evaluate directly
            sim_pt = simulate_lentiviral_persistence(
                promoter_id="synthetic_hCFTR",
                vcn=2.0,
                initial_basal_transduction=0.20,
                total_days=float(d),
            )
            lenti_val = sim_pt.final_delta_isc_pct_wt

        aav_val = max(0.0, aav_initial * math.exp(-aav_lambda * d))
        mrna_val = max(0.0, mrna_initial * math.exp(-mrna_lambda * d))

        comparison.append({
            "day": d,
            "years": round(d / 365.0, 2),
            "lentiviral_bi3720931_pct_wt": round(lenti_val, 2),
            "episomal_aav_4d710_pct_wt": round(aav_val, 2),
            "mrna_lnp_vx522_pct_wt": round(mrna_val, 4),
            "curative_threshold_pct_wt": CURATIVE_THRESHOLD_DELTA_ISC_PCT,
        })

    return {
        "title": "Airway Epithelial Gene Addition Modality Longevity Benchmark",
        "answer_tier": "Tier T2",
        "primary_trials": [CLINICAL_TRIAL_BI3720931, CLINICAL_TRIAL_VX522, CLINICAL_TRIAL_4D710],
        "primary_dois": [PRIMARY_DOI_ALTON, PRIMARY_DOI_COONEY, PRIMARY_DOI_VIJAYKUMAR],
        "comparison": comparison,
    }


def main() -> None:
    """CLI execution for inhaled lentiviral epigenetic persistence model."""
    parser = argparse.ArgumentParser(
        description="Inhaled Lentiviral Epigenetic Silencing & Basal Progenitor Integration Model (T-0407)",
    )
    parser.add_argument(
        "--promoter",
        choices=list(PRESET_PROMOTERS.keys()),
        default="synthetic_hCFTR",
        help="Lentiviral internal promoter design (default: synthetic_hCFTR)",
    )
    parser.add_argument(
        "--vcn",
        type=float,
        default=2.0,
        help="Vector copy number per cell (default: 2.0, range [1.0, 4.0])",
    )
    parser.add_argument(
        "--basal-transduction",
        type=float,
        default=0.20,
        help="Initial transduced fraction of basal stem cells (default: 0.20)",
    )
    parser.add_argument(
        "--days",
        type=float,
        default=3650.0,
        help="Total simulation duration in days (default: 3650.0 / 10 years)",
    )
    parser.add_argument(
        "--compare-modalities",
        action="store_true",
        help="Compare integrating lentivirus against episomal AAV and mRNA LNP",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Output results in machine-readable JSON format",
    )
    parser.add_argument(
        "--summary",
        action="store_true",
        help="Print concise clinical and biophysical summary",
    )
    args = parser.parse_args()

    if args.compare_modalities:
        comparison_res = compare_gene_therapy_modalities(total_days=args.days)
        if args.json:
            print(json.dumps(comparison_res, indent=2))
        else:
            print("=" * 72)
            print("AIRWAY GENE THERAPY MODALITY LONGEVITY BENCHMARK")
            print("=" * 72)
            for row in comparison_res["comparison"]:
                print(
                    f"Day {row['day']:4d} ({row['years']:4.1f} yr) | "
                    f"Lentivirus: {row['lentiviral_bi3720931_pct_wt']:5.1f}% WT | "
                    f"AAV: {row['episomal_aav_4d710_pct_wt']:5.1f}% WT | "
                    f"mRNA: {row['mrna_lnp_vx522_pct_wt']:6.3f}% WT"
                )
        return

    sim = simulate_lentiviral_persistence(
        promoter_id=args.promoter,
        vcn=args.vcn,
        initial_basal_transduction=args.basal_transduction,
        total_days=args.days,
    )

    if args.json:
        payload = asdict(sim)
        payload["nlm_attribution"] = NLM_ATTRIBUTION
        payload["research_boundary"] = RESEARCH_BOUNDARY
        payload["primary_trials"] = [CLINICAL_TRIAL_BI3720931]
        payload["primary_dois"] = [PRIMARY_DOI_ALTON, PRIMARY_DOI_COONEY, PRIMARY_DOI_VIJAYKUMAR]
        payload["structural_model_pdb"] = STRUCTURAL_MODEL_PDB
        print(json.dumps(payload, indent=2))
    elif args.summary:
        print("=" * 72)
        print("LENTIVIRAL BASAL PROGENITOR INTEGRATION & PERSISTENCE SUMMARY (T-0407)")
        print("=" * 72)
        print(f"Promoter:                 {sim.promoter_id} ({PRESET_PROMOTERS[sim.promoter_id].display_name})")
        print(f"Vector Copy Number (VCN): {sim.vcn_per_cell:.2f} copies/cell")
        print(f"Basal Transduction:       {sim.initial_basal_transduction_pct:.1f}%")
        print(f"Peak Apical Current:      {sim.peak_delta_isc_pct_wt:.2f}% WT Delta I_sc")
        print(f"Final Apical Current:     {sim.final_delta_isc_pct_wt:.2f}% WT Delta I_sc (at {args.days / 365.0:.1f} yr)")
        print(f"Curative Threshold Met:   {sim.reaches_curative_threshold} (Threshold >= {CURATIVE_THRESHOLD_DELTA_ISC_PCT}% WT)")
        print(f"Curative Duration:        {sim.curative_duration_years:.2f} years ({sim.curative_duration_days:.0f} days)")
        print(f"Insertional Risk:         {sim.insertional_risk_per_cell:.2e} ({sim.insertional_safety_verdict})")
        print(f"Therapeutic Safety Index: {sim.therapeutic_safety_index:.2f}")
        print("=" * 72)
    else:
        print(f"Lentiviral Persistence Simulation: {sim.promoter_id}, VCN={sim.vcn_per_cell}")
        print(f"Peak: {sim.peak_delta_isc_pct_wt}% WT | Curative Duration: {sim.curative_duration_years} years")


if __name__ == "__main__":
    main()
