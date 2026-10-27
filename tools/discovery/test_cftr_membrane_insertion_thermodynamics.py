#!/usr/bin/env python3
"""Unit tests for Thermodynamic Partitioning of Hydrophobic Surface Residues in CFTR-Membrane Insertion under Trikafta Allosteric Remodeling (T-0453)."""

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

from cftr_membrane_insertion_thermodynamics import (
    APICAL_RAFT_THICKNESS_A,
    CURATIVE_MEMBRANE_STABILIZATION_THRESHOLD_KCAL,
    DOI_AIRWAY_MUCUS,
    DOI_CFTR_HOTSPOT,
    DOI_ETI_NEJM,
    ER_BILAYER_THICKNESS_A,
    NCT_4D710_AAV_TRIAL,
    NCT_VX522_MRNA_TRIAL,
    NLM_ATTRIBUTION,
    PDB_ACTIVE_CFTR,
    PDB_DEPHOSPHO_CFTR,
    RESEARCH_BOUNDARY,
    CFTRMembraneInsertionThermodynamicsEngine,
    main,
)


@pytest.fixture
def engine() -> CFTRMembraneInsertionThermodynamicsEngine:
    return CFTRMembraneInsertionThermodynamicsEngine()


def test_grounding_identifiers_and_disclaimers() -> None:
    """Verify DOIs, NCT registry IDs, PDB structures, and mandatory NLM/research disclaimers."""
    assert NCT_VX522_MRNA_TRIAL == "NCT05668741"
    assert NCT_4D710_AAV_TRIAL == "NCT05248230"
    assert PDB_DEPHOSPHO_CFTR == "6MSM"
    assert PDB_ACTIVE_CFTR == "5UAK"
    assert DOI_CFTR_HOTSPOT == "10.1126/science.aaw7611"
    assert DOI_ETI_NEJM == "10.1056/NEJMoa1908639"
    assert DOI_AIRWAY_MUCUS == "10.1152/physrev.00004.2021"
    assert CURATIVE_MEMBRANE_STABILIZATION_THRESHOLD_KCAL == -3.5
    assert "ClinicalTrials.gov does not evaluate safety or efficacy" in NLM_ATTRIBUTION
    assert "research, not medical advice" in RESEARCH_BOUNDARY


def test_wimley_white_tm1_to_tm12_partitioning(engine: CFTRMembraneInsertionThermodynamicsEngine) -> None:
    """Domain 1 & 2: Verify TM1-12 transfer free energy in 6MSM vs 5UAK and ETI lipid-shell stabilization <= -3.5 kcal/mol."""
    part = engine.evaluate_wimley_white_helix_partitioning(
        tezacaftor_uM=3.0,
        elexacaftor_uM=3.0,
        ivacaftor_uM=1.0,
    )
    assert len(part["helices"]) == 12
    # Phosphorylated active state (5UAK) is thermodynamically more favorably partitioned than dephosphorylated 6MSM
    assert part["total_dg_5uak_kcal"] < part["total_dg_6msm_kcal"]
    assert part["tm3_tm4_tm8_tm11_shell_ddg_kcal"] <= CURATIVE_MEMBRANE_STABILIZATION_THRESHOLD_KCAL
    assert part["meets_shell_remodeling_threshold"] is True

    # Verify TM8 has the least negative (most unstable) intrinsic transfer free energy due to its unwound polar groove
    helices_by_id = {h["helix_id"]: h for h in part["helices"]}
    assert helices_by_id["TM8"]["dg_transfer_6msm_kcal"] > -3.0
    assert helices_by_id["TM8"]["corrector_stabilization_ddg_kcal"] < -1.0


def test_bilayer_hydrophobic_mismatch_and_curvature(engine: CFTRMembraneInsertionThermodynamicsEngine) -> None:
    """Domain 3: Verify hydrophobic mismatch between inward-facing 6MSM and active 5UAK in apical rafts."""
    pm_6msm = engine.calculate_bilayer_mismatch_energetics(
        bilayer_thickness_a=APICAL_RAFT_THICKNESS_A,
        conformation_pdb="6MSM",
        compartment="Apical_PM_Raft",
    )
    pm_5uak = engine.calculate_bilayer_mismatch_energetics(
        bilayer_thickness_a=APICAL_RAFT_THICKNESS_A,
        conformation_pdb="5UAK",
        compartment="Apical_PM_Raft",
    )
    er_6msm = engine.calculate_bilayer_mismatch_energetics(
        bilayer_thickness_a=ER_BILAYER_THICKNESS_A,
        conformation_pdb="6MSM",
        compartment="ER_Membrane",
    )
    # Active 5UAK matches the 31.8 A apical raft far better than tilted 6MSM
    assert abs(pm_5uak.hydrophobic_mismatch_a) < abs(pm_6msm.hydrophobic_mismatch_a)
    assert pm_5uak.total_deformation_penalty_kcal_mol < pm_6msm.total_deformation_penalty_kcal_mol
    assert er_6msm.mismatch_penalty_kcal_mol < pm_6msm.mismatch_penalty_kcal_mol


def test_rare_transmembrane_variant_rescue_l206w_and_v232d(engine: CFTRMembraneInsertionThermodynamicsEngine) -> None:
    """Domain 4: Verify Delta Tm and folded maturation rescue across L206W (TM3) and V232D (TM4)."""
    for variant in ("L206W", "V232D"):
        untreated = engine.evaluate_variant_thermostability(variant, "Untreated")
        mono = engine.evaluate_variant_thermostability(variant, "Tezacaftor_Mono")
        eti = engine.evaluate_variant_thermostability(variant, "ETI_Trikafta")

        assert untreated.folded_maturation_fraction < 0.45
        assert untreated.curative_stabilization_achieved is False
        assert mono.folded_maturation_fraction > untreated.folded_maturation_fraction
        assert eti.folded_maturation_fraction > mono.folded_maturation_fraction
        assert eti.delta_tm_vs_untreated_c > 8.0
        assert eti.curative_stabilization_achieved is True


def test_validation_and_cli_execution(capsys: pytest.CaptureFixture) -> None:
    """Verify input validation errors and CLI execution modes."""
    with pytest.raises(ValueError, match="Temperature in Kelvin must be positive"):
        CFTRMembraneInsertionThermodynamicsEngine(temperature_k=-10.0)

    eng = CFTRMembraneInsertionThermodynamicsEngine()
    with pytest.raises(ValueError, match="Modulator concentrations cannot be negative"):
        eng.evaluate_wimley_white_helix_partitioning(tezacaftor_uM=-1.0)
    with pytest.raises(ValueError, match="Unsupported genotype"):
        eng.evaluate_variant_thermostability("INVALID_VAR")

    assert main(["--check"]) == 0
    assert "Self-check PASS" in capsys.readouterr().out

    assert main(["--json", "--genotype", "L206W"]) == 0
    data = json.loads(capsys.readouterr().out)
    assert data["variant_rescue"]["eti_trikafta"]["genotype"] == "L206W"
    assert data["grounding"]["trials"][0] == "NCT05668741"

    assert main(["--genotype", "V232D"]) == 0
    assert "Sensei Ewok CF Lab - CFTR Membrane Insertion Thermodynamics (T-0453)" in capsys.readouterr().out
