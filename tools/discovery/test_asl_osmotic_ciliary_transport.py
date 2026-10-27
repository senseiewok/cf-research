#!/usr/bin/env python3
"""Unit tests for Biophysical Modeling of Mucociliary Transport Dynamics and Periciliary Liquid Layer Restoration Under ENaC Blockade and Inhaled Osmotic Therapy (T-0459)."""

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

from asl_osmotic_ciliary_transport import (
    CF_BASELINE_CBF_HZ,
    CF_BASELINE_MCT_VELOCITY_MM_MIN,
    CF_BASELINE_PCL_HEIGHT_UM,
    CRITICAL_PCL_STASIS_HEIGHT_UM,
    CURATIVE_CBF_THRESHOLD_HZ,
    CURATIVE_MCT_THRESHOLD_MM_MIN,
    DOI_AIRWAY_MUCUS,
    DOI_ASL_REGULATION,
    DOI_SUPPRESSOR_TRNA,
    HEALTHY_CBF_HZ,
    HEALTHY_MCT_VELOCITY_MM_MIN,
    HEALTHY_PCL_HEIGHT_UM,
    NCT_4D710_AAV_TRIAL,
    NCT_VX522_MRNA_TRIAL,
    NLM_ATTRIBUTION,
    PDB_ACTIVE_CFTR,
    PDB_DEPHOSPHO_CFTR,
    RESEARCH_BOUNDARY,
    SYNERGY_INDEX_THRESHOLD,
    ASLOsmoticCiliaryTransportEngine,
    main,
)


@pytest.fixture
def engine() -> ASLOsmoticCiliaryTransportEngine:
    return ASLOsmoticCiliaryTransportEngine()


def test_grounding_identifiers_and_disclaimers():
    """Verify primary evidence identifiers, registries, and mandatory disclaimers."""
    assert NCT_VX522_MRNA_TRIAL == "NCT05668741"
    assert NCT_4D710_AAV_TRIAL == "NCT05248230"
    assert PDB_DEPHOSPHO_CFTR == "6MSM"
    assert PDB_ACTIVE_CFTR == "5UAK"
    assert DOI_AIRWAY_MUCUS == "10.1152/physrev.00004.2021"
    assert DOI_ASL_REGULATION == "10.1164/rccm.202102-0461oc"
    assert DOI_SUPPRESSOR_TRNA == "10.1038/s41586-023-06133-1"

    assert "ClinicalTrials.gov does not evaluate safety or efficacy" in NLM_ATTRIBUTION
    assert "research, not medical advice" in RESEARCH_BOUNDARY
    assert SYNERGY_INDEX_THRESHOLD == 1.25
    assert CURATIVE_CBF_THRESHOLD_HZ == 8.0
    assert CURATIVE_MCT_THRESHOLD_MM_MIN == 3.0


def test_pcl_hydrodynamics_and_osmotic_influx(engine: ASLOsmoticCiliaryTransportEngine):
    """Domain 1: Assert trans-epithelial osmotic water flux and ENaC absorption kinetics."""
    # Osmotic gradient decays exponentially from peak
    grad_0 = engine.calculate_osmotic_gradient(inhaled_saline_pct=7.0, time_elapsed_min=0.0)
    grad_30 = engine.calculate_osmotic_gradient(inhaled_saline_pct=7.0, time_elapsed_min=30.0)
    assert grad_0 > 200.0  # ~213.5 mOsm/kg
    assert grad_30 == pytest.approx(grad_0 * 0.5, rel=0.05)

    # Osmotic water efflux scales linearly with gradient
    flux_0 = engine.calculate_osmotic_water_efflux(grad_0)
    flux_half = engine.calculate_osmotic_water_efflux(grad_0 * 0.5)
    assert flux_0 > 0.25  # >0.25 um/min initial influx
    assert flux_half == pytest.approx(flux_0 * 0.5, rel=0.01)

    # ENaC absorption decreases under blockade
    absorb_unblocked = engine.calculate_enac_absorption(enac_block_pct=0.0, time_elapsed_min=0.0)
    absorb_blocked_80 = engine.calculate_enac_absorption(enac_block_pct=80.0, time_elapsed_min=0.0)
    assert absorb_unblocked > 0.15
    assert absorb_blocked_80 < absorb_unblocked * 0.30


def test_ciliary_mechanics_and_viscous_drag(engine: ASLOsmoticCiliaryTransportEngine):
    """Domain 2: Assert ciliary beat frequency recovery from collapsed PCL to hydrated extension."""
    # Compressed baseline CF PCL (2.5 um) with hyper-concentrated mucus (8.5 wt%)
    cbf_cf = engine.calculate_cbf(pcl_height_um=2.5, mucus_solids_wt_pct=8.5)
    assert cbf_cf < CURATIVE_CBF_THRESHOLD_HZ
    assert cbf_cf < 4.5

    # Critical stasis threshold (4.0 um)
    cbf_stasis = engine.calculate_cbf(pcl_height_um=4.0, mucus_solids_wt_pct=6.0)
    assert cbf_stasis > cbf_cf

    # Physiological hydrated PCL (7.0 um) with diluted mucus (2.8 wt%)
    cbf_healthy = engine.calculate_cbf(pcl_height_um=7.0, mucus_solids_wt_pct=2.8)
    assert cbf_healthy >= CURATIVE_CBF_THRESHOLD_HZ
    assert cbf_healthy >= 10.0


def test_mucociliary_clearance_velocity_and_stasis(engine: ASLOsmoticCiliaryTransportEngine):
    """Domain 3: Assert MCT transport velocity stasis vs physiological clearance."""
    # Untreated CF stasis (< 0.5 mm/min)
    mct_cf = engine.calculate_mucociliary_velocity(pcl_height_um=2.5, cbf_hz=3.5, mucus_solids_wt_pct=8.5)
    assert mct_cf < 0.5
    assert mct_cf <= CF_BASELINE_MCT_VELOCITY_MM_MIN * 1.5

    # Restored physiological clearance (>= 3.0 mm/min)
    mct_restored = engine.calculate_mucociliary_velocity(pcl_height_um=7.0, cbf_hz=11.0, mucus_solids_wt_pct=2.8)
    assert mct_restored >= CURATIVE_MCT_THRESHOLD_MM_MIN
    assert mct_restored >= 3.5


def test_therapeutic_synergy_dual_vs_monotherapy(engine: ASLOsmoticCiliaryTransportEngine):
    """Domain 4: Assert dual therapy achieves marked synergy over individual monotherapies."""
    synergy_report = engine.evaluate_synergy(enac_block_pct=80.0, inhaled_saline_pct=7.0, duration_min=240.0)

    # Assert synergy index exceeds 1.25 threshold
    assert synergy_report.synergy_index >= SYNERGY_INDEX_THRESHOLD
    assert synergy_report.synergy_verdict == "SYNERGISTIC"
    assert synergy_report.recommendation == "DUAL_THERAPY_STRONGLY_RECOMMENDED"

    # Dual therapy maintains sustained normalization longer than HS alone or ENaC alone
    assert synergy_report.sustained_hours_dual >= 2.0
    assert synergy_report.sustained_hours_dual > synergy_report.sustained_hours_hs_alone * 1.5
    assert synergy_report.dual_therapy_auc > synergy_report.expected_additive_auc


def test_dosing_sweep(engine: ASLOsmoticCiliaryTransportEngine):
    """Assert parametric dosing sweep demonstrates progressive restoration with ENaC blockade."""
    sweep = engine.sweep_enac_concentrations(saline_pct=7.0, enac_steps=4, duration_min=180.0)
    assert len(sweep) == 5

    # Sustained duration increases monotonically with ENaC blockade level
    durations = [row["sustained_duration_min"] for row in sweep]
    for i in range(len(durations) - 1):
        assert durations[i] <= durations[i + 1]


def test_physical_boundary_and_validation_errors():
    """Verify input validation and physical boundaries."""
    with pytest.raises(ValueError, match="Baseline PCL height must be positive"):
        ASLOsmoticCiliaryTransportEngine(baseline_pcl_height_um=-1.0)

    with pytest.raises(ValueError, match="Initial mucus solids must be positive"):
        ASLOsmoticCiliaryTransportEngine(initial_mucus_solids_wt_pct=0.0)

    with pytest.raises(ValueError, match="Cilia length must be positive"):
        ASLOsmoticCiliaryTransportEngine(cilia_length_um=0.0)

    eng = ASLOsmoticCiliaryTransportEngine()
    with pytest.raises(ValueError, match="Inhaled saline percentage cannot be negative"):
        eng.calculate_osmotic_gradient(-1.0, 10.0)

    with pytest.raises(ValueError, match="ENaC blockade percentage must be between 0 and 100"):
        eng.calculate_enac_absorption(120.0, 10.0)


def test_cli_execution(capsys: pytest.CaptureFixture):
    """Verify CLI interface with standard, --check, --json, and --sweep modes."""
    # 1. Standard run
    ret_default = main([])
    assert ret_default == 0
    captured = capsys.readouterr()
    assert "Sensei Ewok CF Lab - ASL Osmotic Ciliary Transport Simulation" in captured.out
    assert "NCT05668741" in captured.out
    assert "Peak PCL Height" in captured.out

    # 2. Self-check run
    ret_check = main(["--check"])
    assert ret_check == 0
    captured_check = capsys.readouterr()
    assert "Self-check PASS" in captured_check.out

    # 3. JSON run
    ret_json = main(["--json", "--enac", "75.0", "--saline", "7.0"])
    assert ret_json == 0
    captured_json = capsys.readouterr()
    data = json.loads(captured_json.out)
    assert "simulation" in data
    assert "synergy_report" in data
    assert data["grounding"]["trials"][0] == "NCT05668741"

    # 4. Sweep run
    ret_sweep = main(["--sweep", "--saline", "7.0"])
    assert ret_sweep == 0
    captured_sweep = capsys.readouterr()
    assert "ENaC Blockade & Osmotic Therapy Dosing Sweep" in captured_sweep.out
