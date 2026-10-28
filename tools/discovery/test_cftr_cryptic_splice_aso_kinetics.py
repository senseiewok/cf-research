#!/usr/bin/env python3
"""Unit tests for Thermodynamic and Kinetic Modeling of Intron-Targeted ASO Splice Correction (T-0470).

Validates:
1. Grounding identifiers, DOIs, NCT registry IDs, PDB structural anchors, and disclaimers.
2. Boltzmann spliceosome partitioning on deep-intronic variant 3849+10kbC->T (c.3718-2477C>T).
3. 2'-MOE vs PMO chemistry thermodynamics (affinity, steric displacement, half-life).
4. Aberrant pseudoexon NMD escape vs canonical transcript stability and total mRNA dynamics.
5. Compound heterozygote (3849+10kbC->T / W1282X) biallelic rescue with suppressor tRNA.
6. ER chaperone competition and ERAD titration by truncated W1282X polypeptides.
7. Modulator (ETI / Trikafta and VTD / Vanzacaftor) folding synergy and Band C maturation.
8. Functional apical chloride transport (I_sc % of WT) and clinical threshold crossings.
9. Extended pharmacokinetics, half-life washout, and therapeutic duration.
10. AH-10 deterministic negative control assertion.
11. CLI execution and JSON schema serialization.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys
import pytest

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import cftr_cryptic_splice_aso_kinetics as ck


@pytest.fixture
def engine() -> ck.CftrCrypticSpliceAsoEngine:
    return ck.CftrCrypticSpliceAsoEngine()


def test_grounding_identifiers_and_disclaimers():
    """Verify primary DOIs, clinical registry NCTs, PDB codes, and mandatory disclaimers."""
    assert ck.CLINICAL_TRIAL_VX522 == "NCT05668741"
    assert ck.CLINICAL_TRIAL_VANZACAFTOR == "NCT05033080"
    assert ck.CLINICAL_TRIAL_SPL84 == "NCT06429176"
    assert ck.PDB_DEPHOSPHO_CFTR == "6MSM"
    assert ck.PDB_TRIKAFTA_CFTR == "8EIQ"

    assert ck.DOI_TRNA_SUPPRESSION == "10.1038/s41586-023-06133-1"
    assert ck.DOI_W1282X_BARRIERS == "10.1021/acs.biochem.5c00590"
    assert ck.DOI_ETI_TRIAL == "10.1056/nejmoa1908639"
    assert ck.DOI_SPL84_PHASE1 == "10.1016/j.jcf.2024.10.004"
    assert ck.DOI_SPL84_DEVELOPMENT == "10.1016/j.jcf.2021.06.003"
    assert ck.DOI_ASO_MUTATION_CORRECTION == "10.1002/humu.22931"
    assert ck.DOI_3849_MOLECULAR_BASIS == "10.1006/geno.1998.5517"

    assert "ClinicalTrials.gov does not evaluate safety or efficacy" in ck.NLM_ATTRIBUTION
    assert "research, not medical advice" in ck.RESEARCH_BOUNDARY
    assert ck.SYMPTOM_RELIEF_THRESHOLD_PCT_WT == 10.0
    assert ck.CLINICAL_RESPONDER_THRESHOLD_PCT_WT == 15.0


def test_boltzmann_splice_partitioning_dose_response(engine: ck.CftrCrypticSpliceAsoEngine):
    """Verify monotonic increase in canonical splicing as a function of 2'-MOE ASO concentration."""
    doses = [0.0, 10.0, 38.0, 100.0, 300.0]
    results = [engine.compute_boltzmann_splice_partition(d, "2MOE") for d in doses]

    # Baseline: without ASO, cryptic donor is energetically favored (Delta G = -9.30 vs -8.20 kcal/mol)
    assert results[0].target_occupancy == 0.0
    assert results[0].canonical_splice_fraction < 0.20
    assert results[0].aberrant_splice_fraction > 0.80

    # At Kd (38 nM), target occupancy ~ 50%
    kd_res = results[2]
    assert 0.45 <= kd_res.target_occupancy <= 0.55

    # Monotonicity check
    for i in range(len(results) - 1):
        assert results[i + 1].canonical_splice_fraction > results[i].canonical_splice_fraction
        assert results[i + 1].restored_isc_pct_wt > results[i].restored_isc_pct_wt

    # At high dose (300 nM), canonical splicing should dominate (>80%)
    high_res = results[-1]
    assert high_res.canonical_splice_fraction > 0.80
    assert high_res.clinical_relief_threshold_met is True
    assert high_res.robust_responder_threshold_met is True


def test_chemistry_comparison_2moe_vs_pmo(engine: ck.CftrCrypticSpliceAsoEngine):
    """Verify 2'-MOE exhibits higher steric displacement potency than PMO at equimolar doses."""
    test_conc = 50.0  # nM
    res_2moe = engine.compute_boltzmann_splice_partition(test_conc, "2MOE")
    res_pmo = engine.compute_boltzmann_splice_partition(test_conc, "PMO")

    # 2'-MOE has lower Kd (38 vs 64 nM) and higher steric penalty (3.50 vs 2.85 kcal/mol)
    assert res_2moe.target_occupancy > res_pmo.target_occupancy
    assert res_2moe.canonical_splice_fraction > res_pmo.canonical_splice_fraction
    assert res_2moe.restored_isc_pct_wt > res_pmo.restored_isc_pct_wt


def test_nmd_escape_and_steady_state_mrna(engine: ck.CftrCrypticSpliceAsoEngine):
    """Verify NMD degradation of aberrant pseudoexon isoform and mRNA pool expansion upon ASO."""
    baseline = engine.compute_boltzmann_splice_partition(0.0, "2MOE")
    treated = engine.compute_boltzmann_splice_partition(100.0, "2MOE")

    # In untreated state, majority of pre-mRNA includes pseudoexon and is degraded by NMD
    assert baseline.steady_state_mrna_pct_wt < 40.0

    # ASO treatment switches splicing to canonical transcript, escaping NMD
    assert treated.steady_state_mrna_pct_wt > 75.0
    assert treated.steady_state_mrna_pct_wt > baseline.steady_state_mrna_pct_wt


def test_compound_het_rescue_and_erad_chaperone_competition(engine: ck.CftrCrypticSpliceAsoEngine):
    """Verify biallelic rescue on 3849+10kbC->T / W1282X and ERAD competition mitigation by modulators."""
    # Untreated baseline
    untreated = engine.simulate_compound_het_rescue(
        aso_concentration_nm=0.0,
        trna_readthrough_efficiency_pct=0.0,
        nmd_inhibition_pct=0.0,
        modulator_regimen="NONE",
    )
    assert untreated.clinical_relief_threshold_met is False
    assert untreated.robust_responder_threshold_met is False
    assert untreated.apical_chloride_current_isc_pct_wt < 10.0

    # Monotherapy ASO alone (100 nM ASO, 0% tRNA, NONE modulator)
    aso_only = engine.simulate_compound_het_rescue(
        aso_concentration_nm=100.0,
        trna_readthrough_efficiency_pct=0.0,
        nmd_inhibition_pct=0.0,
        modulator_regimen="NONE",
    )
    assert aso_only.allele1_canonical_fraction > 0.65
    assert aso_only.allele2_fl_synthesis_rate == 0.0

    # Dual ASO + tRNA (100 nM ASO, 12% tRNA readthrough converting stop codon to full length)
    dual_trna = engine.simulate_compound_het_rescue(
        aso_concentration_nm=100.0,
        trna_readthrough_efficiency_pct=12.0,
        nmd_inhibition_pct=0.0,
        modulator_regimen="NONE",
    )
    assert dual_trna.allele2_fl_synthesis_rate > 0.0
    assert dual_trna.total_mature_band_c_pct_wt > aso_only.total_mature_band_c_pct_wt

    # Manfredi et al. (2025) ERAD competition phenomenon:
    # NMD inhibition without corrector expands truncated W1282X polypeptide, increasing chaperone titration
    dual_nmd = engine.simulate_compound_het_rescue(
        aso_concentration_nm=100.0,
        trna_readthrough_efficiency_pct=12.0,
        nmd_inhibition_pct=25.0,
        modulator_regimen="NONE",
    )
    assert dual_nmd.allele2_truncated_rate > dual_trna.allele2_truncated_rate
    assert dual_nmd.erad_chaperone_competition_index > dual_trna.erad_chaperone_competition_index

    # Triple ASO + tRNA + ETI corrector/potentiator relieves ERAD chaperone competition
    triple = engine.simulate_compound_het_rescue(
        aso_concentration_nm=100.0,
        trna_readthrough_efficiency_pct=12.0,
        nmd_inhibition_pct=25.0,
        modulator_regimen="ETI",
    )
    assert triple.net_band_b_to_band_c_efficiency > dual_nmd.net_band_b_to_band_c_efficiency
    assert triple.apical_chloride_current_isc_pct_wt > dual_nmd.apical_chloride_current_isc_pct_wt
    assert triple.clinical_relief_threshold_met is True
    assert triple.robust_responder_threshold_met is True
    assert triple.synergy_factor > 1.0


def test_vanzacaftor_next_gen_regimen(engine: ck.CftrCrypticSpliceAsoEngine):
    """Verify next-generation VTD regimen provides incremental maturation boost over ETI."""
    res_eti = engine.simulate_compound_het_rescue(
        aso_concentration_nm=100.0,
        trna_readthrough_efficiency_pct=12.0,
        nmd_inhibition_pct=25.0,
        modulator_regimen="ETI",
    )
    res_vtd = engine.simulate_compound_het_rescue(
        aso_concentration_nm=100.0,
        trna_readthrough_efficiency_pct=12.0,
        nmd_inhibition_pct=25.0,
        modulator_regimen="VTD",
    )
    assert res_vtd.net_band_b_to_band_c_efficiency > res_eti.net_band_b_to_band_c_efficiency
    assert res_vtd.apical_chloride_current_isc_pct_wt >= res_eti.apical_chloride_current_isc_pct_wt


def test_pharmacokinetics_washout_dynamics(engine: ck.CftrCrypticSpliceAsoEngine):
    """Verify mono-exponential decay and therapeutic duration of 2'-MOE ASO."""
    pts_2moe = engine.simulate_aso_pharmacokinetics_washout(
        c0_nm=100.0,
        total_days=42.0,
        step_days=7.0,
        chemistry="2MOE",
        trna_readthrough_efficiency_pct=10.0,
        modulator_regimen="ETI",
    )
    # Day 0 concentration is 100 nM
    assert pts_2moe[0].intracellular_aso_nm == 100.0
    # Day 21 (half-life of 2'-MOE) concentration should be ~50 nM
    day21_pt = [p for p in pts_2moe if abs(p.day - 21.0) < 0.1][0]
    assert 48.0 <= day21_pt.intracellular_aso_nm <= 52.0

    # At day 21, therapeutic coverage is maintained above 15% responder threshold
    assert day21_pt.restored_isc_pct_wt >= ck.CLINICAL_RESPONDER_THRESHOLD_PCT_WT


def test_curative_equity_matrix(engine: ck.CftrCrypticSpliceAsoEngine):
    """Verify curative equity scenarios evaluate properly."""
    matrix = engine.evaluate_curative_equity_matrix()
    assert matrix["genotype"] == "3849+10kbC->T / W1282X"
    assert len(matrix["scenarios"]) == 7
    # Verify triple combo is clinical responder
    triple_row = [r for r in matrix["scenarios"] if "Triple ASO + tRNA + ETI" in r["regimen_name"]][0]
    assert triple_row["clinical_responder"] is True
    assert triple_row["restored_isc_pct_wt"] >= ck.CLINICAL_RESPONDER_THRESHOLD_PCT_WT


def test_ah10_deterministic_negative_control(engine: ck.CftrCrypticSpliceAsoEngine):
    """AH-10 negative control: 0 nM ASO + 0% tRNA on 3849/W1282X must fail clinical thresholds."""
    null_res = engine.simulate_compound_het_rescue(
        aso_concentration_nm=0.0,
        trna_readthrough_efficiency_pct=0.0,
        nmd_inhibition_pct=0.0,
        modulator_regimen="NONE",
    )
    assert null_res.clinical_relief_threshold_met is False
    assert null_res.robust_responder_threshold_met is False
    assert null_res.apical_chloride_current_isc_pct_wt < 10.0


def test_parameter_provenance_and_validation():
    """Verify AH-9 parameter provenance dictionary entries."""
    if ck.BiophysicalParameter is None:
        pytest.skip("epistemic_rigor module not available")

    assert len(ck.CORE_BIOPHYSICAL_PARAMETERS) >= 10
    param = ck.CORE_BIOPHYSICAL_PARAMETERS["U1_CRYPTIC_3849_DELTA_G"]
    assert param.value == -9.30
    assert param.units == "kcal/mol"
    assert param.category == "empirical"
    assert "10.1016/j.jcf.2021.06.003" in param.doi_or_pdb


def test_cli_execution_and_json_output(capsys: pytest.CaptureFixture):
    """Verify CLI argument handling and JSON schema."""
    ret = ck.main(["--aso-nm", "75.0", "--trna-readthrough-pct", "10.0", "--json"])
    assert ret == 0
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert payload["task"] == "T-0470"
    assert payload["discipline"] == "rna_gene_editing"
    assert "single_allele_3849_partition" in payload
    assert "compound_heterozygote_rescue" in payload
    assert "curative_equity_scenarios" in payload
    assert "washout_pharmacokinetics" in payload
