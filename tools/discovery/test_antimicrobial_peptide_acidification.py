#!/usr/bin/env python3
"""Unit tests for Biophysical Modeling of Airway Antimicrobial Peptide (LL-37 / Beta-Defensin-2) Inactivation by Acidified Airway Surface Liquid in Cystic Fibrosis (T-0461)."""

from __future__ import annotations

import json
from pathlib import Path
import sys
import pytest

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
for p in (str(HERE), str(REPO_ROOT)):
    if p not in sys.path:
        sys.path.insert(0, p)

from antimicrobial_peptide_acidification import (
    CF_ACIDIFIED_ASL_PH,
    CF_MUCUS_SOLIDS_WT_PCT,
    CURATIVE_LOG10_CFU_REDUCTION_THRESHOLD,
    DOI_AIRWAY_MUCUS,
    DOI_ASL_REGULATION,
    DOI_SUPPRESSOR_TRNA,
    HEALTHY_ASL_PH,
    HEALTHY_MUCUS_SOLIDS_WT_PCT,
    NCT_4D710_AAV_TRIAL,
    NCT_VX522_MRNA_TRIAL,
    NLM_ATTRIBUTION,
    PDB_ACTIVE_CFTR,
    PDB_DEPHOSPHO_CFTR,
    RESEARCH_BOUNDARY,
    AntimicrobialPeptideAcidificationEngine,
    main,
)


@pytest.fixture
def engine() -> AntimicrobialPeptideAcidificationEngine:
    return AntimicrobialPeptideAcidificationEngine()


def test_grounding_identifiers_and_disclaimers() -> None:
    """Verify DOIs, NCT registry IDs, PDB structures, and mandatory NLM/research disclaimers."""
    assert NCT_VX522_MRNA_TRIAL == "NCT05668741"
    assert NCT_4D710_AAV_TRIAL == "NCT05248230"
    assert PDB_DEPHOSPHO_CFTR == "6MSM"
    assert PDB_ACTIVE_CFTR == "5UAK"
    assert DOI_AIRWAY_MUCUS == "10.1152/physrev.00004.2021"
    assert DOI_ASL_REGULATION == "10.1164/rccm.202102-0461oc"
    assert DOI_SUPPRESSOR_TRNA == "10.1038/s41586-023-06133-1"
    assert CURATIVE_LOG10_CFU_REDUCTION_THRESHOLD == 3.0
    assert "ClinicalTrials.gov does not evaluate safety or efficacy" in NLM_ATTRIBUTION
    assert "research, not medical advice" in RESEARCH_BOUNDARY


def test_ph_dependent_electrostatics_and_helicity(engine: AntimicrobialPeptideAcidificationEngine) -> None:
    """Domain 1: Verify amphipathic folding and membrane affinity reduction of LL-37 and hBD-2 at pH 6.5 vs 7.35."""
    for pep in ("LL-37", "hBD-2"):
        cf_state = engine.evaluate_peptide_electrostatics(pep, asl_ph=CF_ACIDIFIED_ASL_PH, bicarbonate_mm=3.5)
        healthy_state = engine.evaluate_peptide_electrostatics(pep, asl_ph=HEALTHY_ASL_PH, bicarbonate_mm=24.0)

        assert healthy_state.amphipathic_helix_fraction > cf_state.amphipathic_helix_fraction * 1.8
        assert healthy_state.intrinsic_membrane_affinity_ka > cf_state.intrinsic_membrane_affinity_ka * 2.0


def test_mucin_polyanionic_sequestration(engine: AntimicrobialPeptideAcidificationEngine) -> None:
    """Domain 2: Verify competitive sequestration by hyper-concentrated 8.5 wt% CF mucins vs 2.5 wt% healthy mucins."""
    cf_trap = engine.calculate_mucin_sequestration(
        peptide_id="LL-37",
        asl_ph=CF_ACIDIFIED_ASL_PH,
        mucus_solids_wt_pct=CF_MUCUS_SOLIDS_WT_PCT,
        bicarbonate_mm=3.5,
    )
    healthy_trap = engine.calculate_mucin_sequestration(
        peptide_id="LL-37",
        asl_ph=HEALTHY_ASL_PH,
        mucus_solids_wt_pct=HEALTHY_MUCUS_SOLIDS_WT_PCT,
        bicarbonate_mm=24.0,
    )
    assert cf_trap.mucin_sequestration_fraction > 0.70
    assert healthy_trap.mucin_sequestration_fraction < 0.35
    assert healthy_trap.free_bioavailable_peptide_ug_ml > cf_trap.free_bioavailable_peptide_ug_ml * 2.2


def test_bactericidal_killing_and_bicarbonate_rescue(engine: AntimicrobialPeptideAcidificationEngine) -> None:
    """Domain 3 & 4: Verify bactericidal failure in acidified CF ASL and restoration via apical bicarbonate transport."""
    rescue = engine.simulate_bicarbonate_rescue_regimens(
        tmem16a_potentiation_pct=75.0,
        slc26a4_pendrin_activation_pct=70.0,
    )
    assert rescue["rescued_asl_ph"] >= 7.25
    assert rescue["rescued_bicarbonate_mm"] >= 20.0

    pa_u = rescue["pseudomonas_ll37_untreated"]
    pa_r = rescue["pseudomonas_ll37_rescued"]
    sa_u = rescue["staph_hbd2_untreated"]
    sa_r = rescue["staph_hbd2_rescued"]

    assert pa_u["bactericidal_clearance_achieved"] is False
    assert pa_u["log10_cfu_reduction_60min"] < 1.5
    assert pa_r["bactericidal_clearance_achieved"] is True
    assert pa_r["log10_cfu_reduction_60min"] >= CURATIVE_LOG10_CFU_REDUCTION_THRESHOLD

    assert sa_u["bactericidal_clearance_achieved"] is False
    assert sa_r["bactericidal_clearance_achieved"] is True
    assert sa_r["effective_mbc99_ug_ml"] < sa_u["effective_mbc99_ug_ml"] * 0.25


def test_validation_and_cli_execution(capsys: pytest.CaptureFixture) -> None:
    """Verify input validation errors and CLI execution flags."""
    with pytest.raises(ValueError, match="Default peptide concentration must be positive"):
        AntimicrobialPeptideAcidificationEngine(default_peptide_dose_ug_ml=0.0)

    eng = AntimicrobialPeptideAcidificationEngine()
    with pytest.raises(ValueError, match="Peptide ID must be"):
        eng.evaluate_peptide_electrostatics("UNKNOWN", 7.2)
    with pytest.raises(ValueError, match="ASL pH must be within"):
        eng.evaluate_peptide_electrostatics("LL-37", 3.0)

    assert main(["--check"]) == 0
    assert "Self-check PASS" in capsys.readouterr().out

    assert main(["--json"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["grounding"]["trials"][0] == "NCT05668741"
    assert payload["bicarbonate_rescue"]["pseudomonas_ll37_rescued"]["bactericidal_clearance_achieved"] is True

    assert main([]) == 0
    assert "Sensei Ewok CF Lab - Airway AMP Acidification & Bicarbonate Rescue (T-0461)" in capsys.readouterr().out
