#!/usr/bin/env python3
"""Unit tests for CFTR Truncation Rescuability and Cryo-EM Interdomain Packing (T-0465).

Validates:
1. Structural coordinate parsing and residue mapping against PDB 6MSM and 8EIQ.
2. Truncation catalog invariants and C-terminal breakpoint properties.
3. Thermodynamic stabilization conferred by VX-661, VX-445, and next-gen VX-121.
4. Allosteric cooperativity and pocket disruption specificity (e.g. Y1092X).
5. Compound heterozygote theratyping and suppressor tRNA synergy.
6. Epistemic rigor and deterministic negative controls (AH-10).
"""
from __future__ import annotations

import json
from pathlib import Path
import pytest
import sys

# Ensure discovery package is on sys.path
HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import cftr_truncation_corrector_rescuability as ctcr


def test_truncation_catalog_invariants():
    """Verify all cataloged variants contain physical parameters and valid bounds."""
    catalog = ctcr.TRUNCATION_CATALOG
    assert "W1282X" in catalog
    assert "S1255X" in catalog
    assert "R1162X" in catalog
    assert "Y1092X" in catalog
    assert "G542X" in catalog
    assert "F508del" in catalog
    assert "WT" in catalog

    for var, bp in catalog.items():
        assert 500 <= bp.synthesized_length <= 1480
        assert bp.delta_delta_g_destab_kcal >= 0.0
        assert 0.0 <= bp.cryo_em_bsa_retention_ratio <= 1.0
        assert isinstance(bp.retains_tmd1_pocket, bool)
        assert isinstance(bp.retains_tmd2_pocket, bool)


def test_corrector_regimens_invariants():
    """Verify corrector profiles specify negative free-energy stabilization values."""
    regimens = ctcr.CORRECTOR_REGIMENS
    assert "NONE" in regimens
    assert "VX-661" in regimens
    assert "VX-445" in regimens
    assert "VX-121" in regimens
    assert "DUAL_ET" in regimens
    assert "TRIPLE_NEXTGEN" in regimens

    for name, prof in regimens.items():
        assert prof.effective_total_delta_g_kcal <= 0.0
        assert prof.delta_delta_g_tmd1_kcal <= 0.0
        assert prof.delta_delta_g_tmd2_kcal <= 0.0


def test_coordinate_loading_and_packing():
    """Verify Cryo-EM coordinate extraction from cftr-data.js for W1282X and Y1092X."""
    data = ctcr.load_cftr_coordinates()
    assert "structures" in data
    assert any(s["id"] == "6MSM" for s in data["structures"])
    assert any(s["id"] == "8EIQ" for s in data["structures"])

    pack_w1282x = ctcr.compute_cryoem_domain_packing("W1282X")
    assert pack_w1282x["variant"] == "W1282X"
    assert pack_w1282x["synthesized_length"] == 1281
    assert pack_w1282x["retains_tmd1_pocket"] is True
    assert pack_w1282x["retains_tmd2_pocket"] is True
    assert pack_w1282x["pdb_6msm_residues_resolved"] > 1000

    pack_y1092x = ctcr.compute_cryoem_domain_packing("Y1092X")
    assert pack_y1092x["retains_tmd2_pocket"] is False
    assert pack_y1092x["bsa_retention_ratio"] < pack_w1282x["bsa_retention_ratio"]


def test_w1282x_dual_corrector_rescuability():
    """Verify thermodynamic rescue of W1282X under DUAL_ET vs uncorrected baseline."""
    base = ctcr.compute_truncation_rescuability("W1282X", regimen="NONE")
    assert base.net_folding_delta_g_kcal > 0.0
    assert base.erad_escape_probability < 0.20
    assert base.band_c_fraction_pct < 5.0
    assert not base.therapeutic_threshold_met

    rescued = ctcr.compute_truncation_rescuability("W1282X", regimen="DUAL_ET")
    assert rescued.corrector_delta_g_stabilization_kcal == -5.00
    assert rescued.net_folding_delta_g_kcal < base.net_folding_delta_g_kcal
    assert rescued.erad_escape_probability > base.erad_escape_probability
    assert rescued.band_c_fraction_pct > base.band_c_fraction_pct


def test_nextgen_vanzacaftor_superiority():
    """Verify next-generation VX-121 confers greater stabilization than first-generation VX-661."""
    res_661 = ctcr.compute_truncation_rescuability("W1282X", regimen="VX-661")
    res_121 = ctcr.compute_truncation_rescuability("W1282X", regimen="VX-121")
    assert res_121.corrector_delta_g_stabilization_kcal < res_661.corrector_delta_g_stabilization_kcal
    assert res_121.band_c_fraction_pct > res_661.band_c_fraction_pct

    res_dual = ctcr.compute_truncation_rescuability("W1282X", regimen="DUAL_ET")
    res_triple = ctcr.compute_truncation_rescuability("W1282X", regimen="TRIPLE_NEXTGEN")
    assert res_triple.corrector_delta_g_stabilization_kcal < res_dual.corrector_delta_g_stabilization_kcal
    assert res_triple.band_c_fraction_pct > res_dual.band_c_fraction_pct


def test_truncation_series_hierarchy():
    """Verify rescue hierarchy aligns with domain integrity: WT > F508del > W1282X > S1255X > R1162X > Y1092X > G542X."""
    series = ctcr.simulate_truncation_series(regimen="DUAL_ET")
    assert series["WT"].band_c_fraction_pct > series["F508del"].band_c_fraction_pct
    assert series["F508del"].band_c_fraction_pct > series["W1282X"].band_c_fraction_pct
    assert series["W1282X"].band_c_fraction_pct > series["S1255X"].band_c_fraction_pct
    assert series["S1255X"].band_c_fraction_pct > series["R1162X"].band_c_fraction_pct
    assert series["R1162X"].band_c_fraction_pct >= series["Y1092X"].band_c_fraction_pct
    assert series["Y1092X"].band_c_fraction_pct >= series["G542X"].band_c_fraction_pct
    assert series["R1162X"].net_folding_delta_g_kcal < series["Y1092X"].net_folding_delta_g_kcal < series["G542X"].net_folding_delta_g_kcal


def test_pocket_disruption_specificity():
    """Verify Y1092X (lacking intact TMD2 pocket) receives zero TMD2 stabilization from VX-445."""
    res_y1092x = ctcr.compute_truncation_rescuability("Y1092X", regimen="VX-445")
    assert res_y1092x.corrector_delta_g_stabilization_kcal == 0.0
    assert res_y1092x.allosteric_coupling_efficiency == 0.0

    res_w1282x = ctcr.compute_truncation_rescuability("W1282X", regimen="VX-445")
    assert res_w1282x.corrector_delta_g_stabilization_kcal == -2.40
    assert res_w1282x.allosteric_coupling_efficiency == 1.0


def test_suppressor_trna_synergy():
    """Verify combination of dual corrector with suppressor tRNA readthrough crosses curative threshold on W1282X."""
    no_trna = ctcr.compute_truncation_rescuability("W1282X", regimen="DUAL_ET", suppressor_trna_efficiency_pct=0.0)
    assert not no_trna.therapeutic_threshold_met

    with_trna = ctcr.compute_truncation_rescuability("W1282X", regimen="DUAL_ET", suppressor_trna_efficiency_pct=8.0)
    assert with_trna.therapeutic_threshold_met
    assert with_trna.apical_isc_conductance_pct_wt >= ctcr.CURATIVE_EFFICACY_THRESHOLD_PCT


def test_compound_heterozygote_modeling():
    """Verify compound heterozygote modeling (W1282X / F508del vs W1282X / W1282X)."""
    het = ctcr.model_compound_heterozygote("W1282X", "F508del", regimen="DUAL_ET")
    hom = ctcr.model_compound_heterozygote("W1282X", "W1282X", regimen="DUAL_ET")

    assert het["composite_apical_isc_pct_wt"] > hom["composite_apical_isc_pct_wt"]
    assert het["curative_threshold_met"] is True
    assert hom["curative_threshold_met"] is False


def test_negative_controls_falsification_suite():
    """Verify all 4 epistemic negative controls pass without failure."""
    res = ctcr.run_falsification_negative_controls()
    assert res["status"] == "PASS"
    assert res["negative_controls_passed"] == 4


def test_parameter_validation_and_errors():
    """Verify strict validation of temperature and readthrough efficiency."""
    with pytest.raises(ValueError, match="Unphysical temperature"):
        ctcr.compute_truncation_rescuability("W1282X", temperature_k=-5.0)

    with pytest.raises(ValueError, match="Unphysical temperature"):
        ctcr.compute_truncation_rescuability("W1282X", temperature_k=600.0)

    with pytest.raises(ValueError, match="Invalid suppressor tRNA efficiency"):
        ctcr.compute_truncation_rescuability("W1282X", suppressor_trna_efficiency_pct=-1.0)

    with pytest.raises(KeyError, match="not found in TRUNCATION_CATALOG"):
        ctcr.get_truncation_breakpoint("NON_EXISTENT_MUTATION")

    with pytest.raises(KeyError, match="not found in CORRECTOR_REGIMENS"):
        ctcr.get_corrector_profile("UNKNOWN_DRUG")


def test_cli_execution():
    """Verify CLI flags execute and return exit code 0."""
    assert ctcr.main(["--check-controls"]) == 0
    assert ctcr.main(["--variant", "W1282X", "--regimen", "DUAL_ET"]) == 0
    assert ctcr.main(["--variant", "W1282X", "--json"]) == 0
    assert ctcr.main(["--series", "--regimen", "TRIPLE_NEXTGEN"]) == 0
    assert ctcr.main(["--heterozygote", "W1282X", "F508del", "--regimen", "DUAL_ET"]) == 0
