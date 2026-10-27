#!/usr/bin/env python3
"""Unit tests for Biophysical Model of Lipid Nanoparticle Corona Remodeling by Mucus Surfactant Proteins (T-0456)."""

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

from lnp_surfactant_corona_remodeling import (
    CF_MUCUS_MEAN_PORE_NM,
    CURATIVE_EPITHELIAL_DELIVERY_THRESHOLD,
    DOI_AIRWAY_MUCUS,
    DOI_PEGYLATION,
    DOI_SORT_LNP,
    DOI_SUPPRESSOR_TRNA,
    NCT_4D710_AAV_TRIAL,
    NCT_VX522_MRNA_TRIAL,
    NLM_ATTRIBUTION,
    OPTIMAL_SHIELDING_MAX_MOL_PCT,
    OPTIMAL_SHIELDING_MIN_MOL_PCT,
    PDB_ACTIVE_CFTR,
    PDB_DEPHOSPHO_CFTR,
    RESEARCH_BOUNDARY,
    SP_A_MAX_CORONA_THICKNESS_NM,
    SP_D_MAX_CORONA_THICKNESS_NM,
    LNPSurfactantCoronaRemodelingEngine,
    main,
)


@pytest.fixture
def engine() -> LNPSurfactantCoronaRemodelingEngine:
    return LNPSurfactantCoronaRemodelingEngine(core_diameter_nm=65.0, mesh_pore_nm=90.0)


def test_grounding_identifiers_and_disclaimers():
    """Verify primary evidence identifiers, registries, and mandatory disclaimers."""
    assert NCT_VX522_MRNA_TRIAL == "NCT05668741"
    assert NCT_4D710_AAV_TRIAL == "NCT05248230"
    assert PDB_DEPHOSPHO_CFTR == "6MSM"
    assert PDB_ACTIVE_CFTR == "5UAK"
    assert DOI_AIRWAY_MUCUS == "10.1152/physrev.00004.2021"
    assert DOI_SORT_LNP == "10.1038/s41565-020-0669-6"
    assert DOI_PEGYLATION == "10.1016/j.addr.2015.09.012"
    assert DOI_SUPPRESSOR_TRNA == "10.1038/s41586-023-06133-1"

    assert "ClinicalTrials.gov does not evaluate safety or efficacy" in NLM_ATTRIBUTION
    assert "research, not medical advice" in RESEARCH_BOUNDARY
    assert CURATIVE_EPITHELIAL_DELIVERY_THRESHOLD == 0.80
    assert OPTIMAL_SHIELDING_MIN_MOL_PCT == 0.5
    assert OPTIMAL_SHIELDING_MAX_MOL_PCT == 2.5


def test_collectin_opsonization_unshielded_vs_shielded(engine: LNPSurfactantCoronaRemodelingEngine):
    """Domain 1: Assert competitive collectin adsorption and +15 to +35 nm hydrodynamic diameter increase."""
    # Unshielded LNP (0.0 mol% stealth shield)
    unshielded = engine.calculate_collectin_adsorption(shielding_mol_pct=0.0)
    assert unshielded.steric_brush_factor == 1.0
    assert unshielded.theta_sp_a > 0.40
    assert unshielded.theta_sp_d > 0.10
    assert unshielded.theta_total > 0.55
    # Hydrodynamic diameter expansion: +15 to +35 nm
    assert 15.0 <= unshielded.delta_d_h_nm <= 35.0
    assert unshielded.hydrodynamic_diameter_nm == pytest.approx(65.0 + unshielded.delta_d_h_nm, abs=0.1)

    # Stealth-shielded LNP (1.8 mol% stealth shield)
    shielded = engine.calculate_collectin_adsorption(shielding_mol_pct=1.8)
    assert shielded.steric_brush_factor < 0.20
    assert shielded.theta_total < 0.15
    # Corona expansion suppressed to < 5 nm
    assert shielded.delta_d_h_nm < 5.0
    assert shielded.hydrodynamic_diameter_nm < 70.0


def test_mucus_mesh_entrapment_and_diffusivity(engine: LNPSurfactantCoronaRemodelingEngine):
    """Domain 2: Assert Stokes-Einstein diffusivity and steric filtering through CF mucus hydrogels."""
    unshielded_adsorp = engine.calculate_collectin_adsorption(0.0)
    unshielded_diff = engine.calculate_mucus_diffusivity(unshielded_adsorp)

    # Bulk diffusivity is finite and physical
    assert unshielded_diff.bulk_diffusivity_d0_um2_s > 0.0

    # Unshielded LNP suffers severe steric obstruction and adhesive drag
    assert unshielded_diff.diffusivity_ratio < 0.15
    assert unshielded_diff.entrapment_verdict in ("ENTRAPPED", "RETARDED")

    # Shielded LNP preserves mucus mesh penetration
    shielded_adsorp = engine.calculate_collectin_adsorption(2.0)
    shielded_diff = engine.calculate_mucus_diffusivity(shielded_adsorp)
    assert shielded_diff.diffusivity_ratio >= 0.35
    assert shielded_diff.entrapment_verdict == "PENETRATING"
    assert shielded_diff.effective_diffusivity_um2_s > unshielded_diff.effective_diffusivity_um2_s * 2.0


def test_macrophage_clearance_vs_epithelial_uptake(engine: LNPSurfactantCoronaRemodelingEngine):
    """Domain 3: Assert opsonin-driven macrophage phagocytosis vs productive epithelial endocytosis."""
    # Unshielded LNP
    unshielded_eval = engine.evaluate_formulation("Unshielded-LNP", shielding_mol_pct=0.0)
    fate_unshielded = unshielded_eval.cellular_fate
    # Macrophage phagocytosis dominates
    assert fate_unshielded.fraction_macrophage > 0.45
    assert fate_unshielded.fraction_epithelial < 0.45
    assert not fate_unshielded.curative_delivery_goal_met
    assert "SUB_THERAPEUTIC" in unshielded_eval.recommended_disposition

    # Shielded LNP (1.8 mol%)
    shielded_eval = engine.evaluate_formulation("Shielded-1.8molPct", shielding_mol_pct=1.8)
    fate_shielded = shielded_eval.cellular_fate
    # Productive epithelial delivery exceeds >80% threshold
    assert fate_shielded.fraction_epithelial >= CURATIVE_EPITHELIAL_DELIVERY_THRESHOLD
    assert fate_shielded.fraction_macrophage < 0.15
    assert fate_shielded.curative_delivery_goal_met
    assert shielded_eval.recommended_disposition == "OPTIMAL_STEALTH_CURATIVE"


def test_shielding_optimization_sweep(engine: LNPSurfactantCoronaRemodelingEngine):
    """Domain 4: Assert parameter sweep identifies optimal shielding window in [0.5, 2.5] mol%."""
    sweep = engine.sweep_shielding_optimization(mol_pct_start=0.0, mol_pct_end=3.0, steps=10)
    assert len(sweep) == 10

    # Monotonic reduction in corona diameter expansion as shielding increases
    delta_dh_values = [e.collectin_adsorption.delta_d_h_nm for e in sweep]
    for i in range(len(delta_dh_values) - 1):
        assert delta_dh_values[i] >= delta_dh_values[i + 1]

    # Check that formulations in [1.5, 2.5] mol% achieve curative goal
    curative_candidates = [
        e for e in sweep
        if OPTIMAL_SHIELDING_MIN_MOL_PCT <= e.shielding_mol_pct <= OPTIMAL_SHIELDING_MAX_MOL_PCT
        and e.cellular_fate.curative_delivery_goal_met
    ]
    assert len(curative_candidates) > 0


def test_physical_boundary_and_validation_errors():
    """Verify input validation and physical boundaries."""
    with pytest.raises(ValueError, match="Core diameter must be positive"):
        LNPSurfactantCoronaRemodelingEngine(core_diameter_nm=0.0)

    with pytest.raises(ValueError, match="Mesh pore size must be positive"):
        LNPSurfactantCoronaRemodelingEngine(mesh_pore_nm=-10.0)

    eng = LNPSurfactantCoronaRemodelingEngine()
    with pytest.raises(ValueError, match="Shielding molar percent cannot be negative"):
        eng.calculate_collectin_adsorption(-0.5)


def test_cli_execution(capsys: pytest.CaptureFixture):
    """Verify CLI interface with standard and JSON outputs."""
    ret_default = main([])
    assert ret_default == 0
    captured = capsys.readouterr()
    assert "Sensei Ewok CF Lab - LNP Collectin Corona Remodeling Simulation" in captured.out
    assert "NCT05668741" in captured.out

    ret_json = main(["--json", "--shielding", "2.0"])
    assert ret_json == 0
    captured_json = capsys.readouterr()
    data = json.loads(captured_json.out)
    assert data["core_diameter_nm"] == 65.0
    assert data["shielding_mol_pct"] == 2.0
    assert "collectin_adsorption" in data

    ret_sweep = main(["--sweep", "--core-diameter", "70.0"])
    assert ret_sweep == 0
    captured_sweep = capsys.readouterr()
    assert "Optimization Sweep" in captured_sweep.out
