#!/usr/bin/env python3
"""Unit tests for inhaled lentiviral vector epigenetic silencing & basal integration persistence model (T-0407).

Tests biophysical invariants, stem cell lineage kinetics, CpG island methylation attenuation,
proviral copy-number mutagenesis safety bounds, transepithelial chloride current rescue,
and multi-modality longevity benchmarks.
"""

from __future__ import annotations

import math
import pytest

from lentiviral_basal_integration_persistence import (
    CFTR_BAND_C_MATURATION_EFFICIENCY,
    CLINICAL_TRIAL_4D710,
    CLINICAL_TRIAL_BI3720931,
    CLINICAL_TRIAL_VX522,
    CURATIVE_THRESHOLD_DELTA_ISC_PCT,
    DEFAULT_BASAL_HALF_LIFE_DAYS,
    DEFAULT_LUMINAL_HALF_LIFE_DAYS,
    NLM_ATTRIBUTION,
    PRESET_PROMOTERS,
    PRIMARY_DOI_ALTON,
    PRIMARY_DOI_COONEY,
    PRIMARY_DOI_VIJAYKUMAR,
    REGULATORY_MAX_INSERTIONAL_RISK,
    RESEARCH_BOUNDARY,
    STRUCTURAL_MODEL_PDB,
    compare_gene_therapy_modalities,
    compute_promoter_activity,
    evaluate_insertional_mutagenesis_risk,
    simulate_lentiviral_persistence,
)


def test_grounding_identifiers_and_disclaimers():
    """Verify primary DOIs, clinical trial IDs, structural model, and mandatory disclaimers."""
    assert PRIMARY_DOI_ALTON == "10.1016/S2213-2600(15)00245-3"
    assert PRIMARY_DOI_COONEY == "10.1016/j.ymthe.2018.04.014"
    assert PRIMARY_DOI_VIJAYKUMAR == "10.1172/jci.insight.206375"
    assert CLINICAL_TRIAL_BI3720931 == "NCT06515002"
    assert CLINICAL_TRIAL_VX522 == "NCT05668741"
    assert CLINICAL_TRIAL_4D710 == "NCT05248230"
    assert STRUCTURAL_MODEL_PDB == "6MSM"

    assert "ClinicalTrials.gov does not evaluate safety or efficacy" in NLM_ATTRIBUTION
    assert "research, not medical advice" in RESEARCH_BOUNDARY
    assert CURATIVE_THRESHOLD_DELTA_ISC_PCT == 10.0
    assert CFTR_BAND_C_MATURATION_EFFICIENCY == 0.72


def test_promoter_presets_and_cpg_attenuation():
    """Verify promoter presets, CpG island profiles, and epigenetic attenuation kinetics."""
    assert "EF1a" in PRESET_PROMOTERS
    assert "MND" in PRESET_PROMOTERS
    assert "synthetic_hCFTR" in PRESET_PROMOTERS

    ef1a = PRESET_PROMOTERS["EF1a"]
    mnd = PRESET_PROMOTERS["MND"]
    syn = PRESET_PROMOTERS["synthetic_hCFTR"]

    # synthetic_hCFTR has lowest CpG sites, lowest silencing rate, and highest floor
    assert syn.cpg_island_sites < mnd.cpg_island_sites < ef1a.cpg_island_sites
    assert syn.silencing_rate_per_day < ef1a.silencing_rate_per_day < mnd.silencing_rate_per_day
    assert syn.expression_floor_fraction > ef1a.expression_floor_fraction > mnd.expression_floor_fraction
    assert syn.insulator_protection is True
    assert ef1a.insulator_protection is False

    # Day 0 activity must be 1.0
    assert compute_promoter_activity(0.0, syn) == 1.0
    assert compute_promoter_activity(0.0, ef1a) == 1.0

    # Over time, activity must decay monotonically towards expression floor
    act_1yr = compute_promoter_activity(365.0, syn)
    act_5yr = compute_promoter_activity(1825.0, syn)
    act_10yr = compute_promoter_activity(3650.0, syn)

    assert 1.0 > act_1yr > act_5yr > act_10yr >= syn.expression_floor_fraction


def test_insertional_mutagenesis_risk():
    """Verify insertional oncogenesis risk evaluation and safety thresholds."""
    zero_risk, zero_verdict = evaluate_insertional_mutagenesis_risk(0.0)
    assert zero_risk == 0.0
    assert zero_verdict == "PASS_ZERO_RISK"

    # Within clinical VCN window [1.0, 4.0]
    for vcn in [1.0, 2.0, 3.0, 4.0]:
        risk, verdict = evaluate_insertional_mutagenesis_risk(vcn)
        assert 0.0 < risk <= REGULATORY_MAX_INSERTIONAL_RISK
        assert verdict == "PASS_WITHIN_SAFETY_MARGIN"

    # Monotonic scaling with copy number
    risk_1, _ = evaluate_insertional_mutagenesis_risk(1.0)
    risk_2, _ = evaluate_insertional_mutagenesis_risk(2.0)
    risk_4, _ = evaluate_insertional_mutagenesis_risk(4.0)
    assert risk_1 < risk_2 < risk_4


def test_simulation_dynamics_and_basal_retention():
    """Verify coupled lineage ODE simulation, basal stem cell retention, and curative duration."""
    sim = simulate_lentiviral_persistence(
        promoter_id="synthetic_hCFTR",
        vcn=2.0,
        initial_basal_transduction=0.20,
        total_days=3650.0,
    )

    assert sim.promoter_id == "synthetic_hCFTR"
    assert sim.vcn_per_cell == 2.0
    assert sim.initial_basal_transduction_pct == 20.0
    assert sim.reaches_curative_threshold is True
    assert sim.peak_delta_isc_pct_wt >= CURATIVE_THRESHOLD_DELTA_ISC_PCT
    assert sim.curative_duration_years >= 5.0
    assert sim.insertional_safety_verdict == "PASS_WITHIN_SAFETY_MARGIN"
    assert sim.therapeutic_safety_index > 0.0

    # Inspect annual trajectory
    traj = sim.yearly_trajectory
    assert len(traj) >= 10
    year0 = traj[0]
    year10 = traj[-1]

    assert year0.years == 0.0
    assert year10.years == 10.0
    # Basal cells retain provirus over multi-year timescale (half-life 3650 days = ~50% at 10 years)
    assert year10.basal_transduced_fraction >= 0.08
    # Luminal differentiated progeny express CFTR via clonal transmission
    assert year10.luminal_transduced_fraction > 0.05
    assert year10.apical_chloride_current_pct_wt >= CURATIVE_THRESHOLD_DELTA_ISC_PCT


def test_promoter_longevity_comparison():
    """Verify synthetic_hCFTR promoter exhibits superior long-term retention vs EF1a and MND."""
    sim_syn = simulate_lentiviral_persistence("synthetic_hCFTR", vcn=2.0, total_days=3650.0)
    sim_ef1a = simulate_lentiviral_persistence("EF1a", vcn=2.0, total_days=3650.0)
    sim_mnd = simulate_lentiviral_persistence("MND", vcn=2.0, total_days=3650.0)

    # At year 10 (day 3650), synthetic_hCFTR should have higher residual current due to CpG depletion
    assert sim_syn.final_delta_isc_pct_wt > sim_ef1a.final_delta_isc_pct_wt
    assert sim_syn.final_delta_isc_pct_wt > sim_mnd.final_delta_isc_pct_wt
    assert sim_syn.curative_duration_years >= sim_ef1a.curative_duration_years


def test_vcn_parameter_boundaries():
    """Verify input validation and boundary checks on VCN and promoter identifiers."""
    with pytest.raises(ValueError, match="Unknown promoter"):
        simulate_lentiviral_persistence(promoter_id="nonexistent_promoter")

    with pytest.raises(ValueError, match="outside realistic range"):
        simulate_lentiviral_persistence(vcn=0.1)

    with pytest.raises(ValueError, match="outside realistic range"):
        simulate_lentiviral_persistence(vcn=10.0)


def test_gene_therapy_modality_comparison():
    """Verify multi-modality longevity benchmark contrasts integrating lentivirus with AAV and mRNA."""
    benchmark = compare_gene_therapy_modalities(total_days=3650.0)

    assert benchmark["answer_tier"] == "Tier T2"
    assert CLINICAL_TRIAL_BI3720931 in benchmark["primary_trials"]
    assert CLINICAL_TRIAL_4D710 in benchmark["primary_trials"]
    assert CLINICAL_TRIAL_VX522 in benchmark["primary_trials"]

    rows = {r["day"]: r for r in benchmark["comparison"]}

    # Day 0: all active
    assert rows[0]["lentiviral_bi3720931_pct_wt"] > 0.0
    assert rows[0]["episomal_aav_4d710_pct_wt"] == 30.0
    assert rows[0]["mrna_lnp_vx522_pct_wt"] == 40.0

    # Day 180: mRNA decayed, AAV and Lentivirus still active
    assert rows[180]["mrna_lnp_vx522_pct_wt"] < 0.01
    assert rows[180]["episomal_aav_4d710_pct_wt"] > 5.0
    assert rows[180]["lentiviral_bi3720931_pct_wt"] >= CURATIVE_THRESHOLD_DELTA_ISC_PCT

    # Day 3650 (10 years): Lentivirus remains curative (>= 10% WT), while AAV and mRNA are near zero
    assert rows[3650]["lentiviral_bi3720931_pct_wt"] >= CURATIVE_THRESHOLD_DELTA_ISC_PCT
    assert rows[3650]["episomal_aav_4d710_pct_wt"] < 0.001
    assert rows[3650]["mrna_lnp_vx522_pct_wt"] == 0.0
