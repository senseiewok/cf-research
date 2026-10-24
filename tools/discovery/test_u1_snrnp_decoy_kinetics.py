#!/usr/bin/env python3
"""Unit test suite for CFTR Pre-mRNA U1 snRNP Decoy Kinetics & Cryptic Exon Exclusion Model (T-0292).

Acceptance Criteria:
  pytest tools/discovery/test_u1_snrnp_decoy_kinetics.py -v; python tools/claims/check_anti_hallucination.py
"""
import math
import sys
from pathlib import Path
import pytest

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from u1_snrnp_decoy_kinetics import (
    MUTATION_SPLICE_CATALOG,
    NLM_ATTRIBUTION,
    RESEARCH_BOUNDARY,
    SSO_CHEMISTRY_CATALOG,
    MutationSpliceProfile,
    SSOChemistryProfile,
    calculate_splicing_restoration_dose_response,
    compute_competitive_occupancy,
    compute_u1_dissociation_constant,
    evaluate_discovery_panel,
    main,
)


class TestCatalogAndGrounding:
    """Validate deep-intronic mutation profiles, SSO chemistry profiles, and regulatory disclaimers."""

    def test_mutation_catalog_entries(self):
        assert "3849_plus_10kb_C_to_T" in MUTATION_SPLICE_CATALOG
        assert "1811_plus_1_6kb_A_to_G" in MUTATION_SPLICE_CATALOG

    def test_sso_chemistry_catalog_entries(self):
        assert "SPL84_2MOE_donor_decoy" in SSO_CHEMISTRY_CATALOG
        assert "PMO_pseudoexon_blocker" in SSO_CHEMISTRY_CATALOG
        assert "MOE_ISS_enhancer_blocker" in SSO_CHEMISTRY_CATALOG

    def test_regulatory_and_boundary_statements(self):
        assert "ClinicalTrials.gov does not evaluate safety or efficacy" in NLM_ATTRIBUTION
        assert "Preclinical in-silico biophysical modeling" in RESEARCH_BOUNDARY
        prohibited = ["proven to cure", "guarantees clinical benefit", "cures cystic fibrosis"]
        for p in prohibited:
            assert p not in NLM_ATTRIBUTION.lower()
            assert p not in RESEARCH_BOUNDARY.lower()

    def test_clinical_trial_grounding(self):
        for name, mut in MUTATION_SPLICE_CATALOG.items():
            assert mut.clinical_trial_nct == "NCT06429176"
            assert mut.evidence_doi in ("10.1016/j.jcf.2024.10.004", "10.1056/NEJM199410133311503")


class TestBiophysicalMechanisms:
    """Validate base-pairing thermodynamics, dissociation constants, and competitive occupancy."""

    def test_cryptic_donor_thermodynamic_contrast(self):
        mut = MUTATION_SPLICE_CATALOG["3849_plus_10kb_C_to_T"]
        # Mutant donor has consensus GTAAGT sequence with strong negative Delta G
        assert mut.mut_delta_g_kcal_mol < -9.0
        # WT donor GCAAGT is weak non-consensus
        assert mut.wt_delta_g_kcal_mol > -4.5
        assert mut.mut_delta_g_kcal_mol < mut.wt_delta_g_kcal_mol

    def test_kd_calculation_from_delta_g(self):
        kd_wt = compute_u1_dissociation_constant(-3.2)
        kd_mut = compute_u1_dissociation_constant(-10.4)
        # Higher affinity (more negative Delta G) means much lower Kd
        assert kd_mut < kd_wt
        assert kd_mut < 100.0  # Nanomolar range for consensus donor
        assert kd_wt > 1000.0  # Micromolar range for non-consensus

    def test_competitive_mass_action_occupancy(self):
        # At zero SSO, U1 occupies majority of available donor sites
        t_sso_0, t_u1_0 = compute_competitive_occupancy(sso_conc_nm=0.0, kd_sso_nm=0.8, u1_conc_nm=30.0, kd_u1_nm=25.0)
        assert math.isclose(t_sso_0, 0.0, abs_tol=1e-5)
        assert t_u1_0 > 0.40

        # At high SSO (100 nM), SSO displaces U1 from donor sites
        t_sso_hi, t_u1_hi = compute_competitive_occupancy(sso_conc_nm=100.0, kd_sso_nm=0.8, u1_conc_nm=30.0, kd_u1_nm=25.0)
        assert t_sso_hi > 0.90
        assert t_u1_hi < 0.05
        assert t_sso_hi + t_u1_hi <= 1.0


class TestDoseTitrationAndPanel:
    """Validate dose-response curves, chloride current rescue, and discovery panel."""

    def test_dose_dependent_splicing_restoration(self):
        mut = MUTATION_SPLICE_CATALOG["3849_plus_10kb_C_to_T"]
        sso = SSO_CHEMISTRY_CATALOG["SPL84_2MOE_donor_decoy"]

        titration = calculate_splicing_restoration_dose_response(mutation=mut, sso=sso)
        assert len(titration) == 10

        t_zero = titration[0]
        t_high = titration[-1]

        # Baseline reflects untreated defective splicing
        assert math.isclose(t_zero.full_length_mrna_fraction, mut.baseline_full_length_fraction, rel_tol=1e-3)
        assert not t_zero.therapeutic_threshold_met

        # High dose restores >80% full-length mRNA and surpasses therapeutic threshold
        assert t_high.full_length_mrna_fraction > 0.80
        assert t_high.restoration_relative_to_wt >= 0.25
        assert t_high.therapeutic_threshold_met

        # Check monotonic increase in full-length mRNA
        for i in range(len(titration) - 1):
            assert titration[i + 1].full_length_mrna_fraction >= titration[i].full_length_mrna_fraction
            assert titration[i + 1].theta_u1_occupancy <= titration[i].theta_u1_occupancy

    def test_discovery_panel_coverage(self):
        panel = evaluate_discovery_panel()
        assert panel["metadata"]["task"] == "T-0292"
        assert "3849+10kb C->T" in panel["mutations"]
        assert "1811+1.6kb A->G" in panel["mutations"]

        prof_3849 = panel["mutations"]["3849+10kb C->T"]
        assert prof_3849["clinical_trial_nct"] == "NCT06429176"
        assert prof_3849["restored_full_length_100nm_pct"] > 75.0
        assert prof_3849["therapeutic_threshold_met"] is True

    def test_cli_execution(self, capsys):
        ret = main(["--mutation", "3849_plus_10kb_C_to_T", "--sso", "SPL84_2MOE_donor_decoy"])
        assert ret == 0
        captured = capsys.readouterr().out
        assert "CFTR U1 snRNP DECOY KINETICS" in captured
        assert "NCT06429176" in captured
        assert "SSO DOSE TITRATION" in captured

    def test_cli_json_mode(self, capsys):
        ret = main(["--json"])
        assert ret == 0
        captured = capsys.readouterr().out
        import json
        data = json.loads(captured)
        assert data["metadata"]["task"] == "T-0292"
        assert "3849+10kb C->T" in data["mutations"]
