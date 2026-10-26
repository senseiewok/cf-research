#!/usr/bin/env python3
"""Comprehensive Unit Tests for Airway Mucus-Penetrating Particle (MPP) Obstruction Scaling (T-0244).

Verifies:
  1. Stokes-Einstein hydrodynamic diffusivity (D_0) in water at 37 °C.
  2. Flory polymer radius (R_F) and mushroom-to-brush transition density (Gamma*).
  3. Mucin & eDNA hydrogel mesh pore scaling (xi_mesh) in healthy vs CF sputum.
  4. Enzymatic rhDNase I (dornase alfa) mesh pore expansion.
  5. Ogston / Brinkman hydrodynamic steric exclusion factor (F_steric).
  6. Dense PEG brush shielding and polyvalent adhesive obstruction factor (F_adhesive).
  7. Nanoparticle mucus transit kinetics and mucociliary clearance escape boundaries.
  8. Conventional Particle (CPP) vs Mucus-Penetrating Particle (MPP) comparative ratio.
  9. Deterministic bit-for-bit reproducibility.
  10. CLI execution, JSON payload schema, and sweep consistency.
  11. NLM attribution and anti-hallucination guardrail compliance.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
import subprocess
import sys

import pytest

DISCOVERY_DIR = Path(__file__).resolve().parent
if str(DISCOVERY_DIR) not in sys.path:
    sys.path.insert(0, str(DISCOVERY_DIR))

import mpp_obstruction_scaling as mpp


# ---------------------------------------------------------------------------
# 1. Stokes-Einstein Diffusivity Tests
# ---------------------------------------------------------------------------
def test_stokes_einstein_d0_inverse_scaling() -> None:
    """Verify D_0 scales inversely with hydrodynamic diameter."""
    d0_50 = mpp.calculate_stokes_einstein_d0(50.0)
    d0_100 = mpp.calculate_stokes_einstein_d0(100.0)
    d0_200 = mpp.calculate_stokes_einstein_d0(200.0)

    assert d0_50 > d0_100 > d0_200
    # Exactly inverse linear scaling: d0_50 / d0_100 ~= 2.0
    assert math.isclose(d0_50 / d0_100, 2.0, rel_tol=1e-3)
    assert math.isclose(d0_100 / d0_200, 2.0, rel_tol=1e-3)


def test_stokes_einstein_nominal_65nm() -> None:
    """Verify nominal 65 nm nanoparticle diffusivity matches physical expectation (~9.7 um^2/s)."""
    d0_65 = mpp.calculate_stokes_einstein_d0(65.0)
    # At 37 °C, mu = 6.91e-4 Pa*s, r = 32.5 nm -> D0 = 9.68 um^2/s
    assert 9.0 <= d0_65 <= 10.5


def test_stokes_einstein_invalid_input() -> None:
    """Ensure non-positive diameter raises ValueError."""
    with pytest.raises(ValueError, match="diameter must be positive"):
        mpp.calculate_stokes_einstein_d0(0.0)
    with pytest.raises(ValueError, match="diameter must be positive"):
        mpp.calculate_stokes_einstein_d0(-20.0)


# ---------------------------------------------------------------------------
# 2. Flory Radius & Brush Transition Tests
# ---------------------------------------------------------------------------
def test_flory_radius_mw_scaling() -> None:
    """Verify Flory radius R_F scales with M_w^0.6."""
    rf_2k, gamma_star_2k = mpp.calculate_flory_radius_and_brush_transition(2.0)
    rf_5k, gamma_star_5k = mpp.calculate_flory_radius_and_brush_transition(5.0)
    rf_10k, gamma_star_10k = mpp.calculate_flory_radius_and_brush_transition(10.0)

    assert rf_2k < rf_5k < rf_10k
    assert gamma_star_2k > gamma_star_5k > gamma_star_10k

    # For 5 kDa PEG: N ~ 113.5, R_F ~ 0.35 * (113.5)^0.6 ~= 6.0 nm
    assert 5.0 <= rf_5k <= 7.0
    assert 0.005 <= gamma_star_5k <= 0.015


# ---------------------------------------------------------------------------
# 3. Mucin & eDNA Mesh Pore Scaling Tests
# ---------------------------------------------------------------------------
def test_mucus_mesh_pore_cf_vs_healthy() -> None:
    """Verify CF mucus mesh pore size is significantly tighter than healthy mucus."""
    healthy_mesh = mpp.calculate_mucus_mesh_spacing(
        mucin_volume_fraction=0.02,
        edna_conc_mg_ml=0.0,
    )
    cf_mesh = mpp.calculate_mucus_mesh_spacing(
        mucin_volume_fraction=0.065,
        edna_conc_mg_ml=2.8,
    )

    # Healthy is ~220 nm baseline, CF is ~70-100 nm
    assert healthy_mesh >= 200.0
    assert 55.0 <= cf_mesh <= 110.0
    assert cf_mesh < (healthy_mesh * 0.5)


def test_rhdnase_mesh_expansion() -> None:
    """Verify enzymatic rhDNase I digestion expands effective mesh pores."""
    cf_baseline = mpp.calculate_mucus_mesh_spacing(
        mucin_volume_fraction=0.065,
        edna_conc_mg_ml=2.8,
        rhdnase_treated=False,
    )
    cf_digested = mpp.calculate_mucus_mesh_spacing(
        mucin_volume_fraction=0.065,
        edna_conc_mg_ml=2.8,
        rhdnase_treated=True,
        rhdnase_cleavage_efficiency=0.85,
    )

    assert cf_digested > cf_baseline
    # Digestion expands mesh by >1.3x
    assert (cf_digested / cf_baseline) >= 1.25


# ---------------------------------------------------------------------------
# 4. Steric Exclusion Tests
# ---------------------------------------------------------------------------
def test_steric_exclusion_size_dependent_filtering() -> None:
    """Verify small particles pass steric mesh while large particles are filtered."""
    mesh_nm = 80.0
    f_steric_30 = mpp.calculate_steric_exclusion_factor(30.0, mesh_nm)
    f_steric_65 = mpp.calculate_steric_exclusion_factor(65.0, mesh_nm)
    f_steric_120 = mpp.calculate_steric_exclusion_factor(120.0, mesh_nm)
    f_steric_200 = mpp.calculate_steric_exclusion_factor(200.0, mesh_nm)

    assert f_steric_30 > f_steric_65 > f_steric_120 >= f_steric_200
    assert f_steric_30 >= 0.40
    # 200 nm particle severely filtered in 80 nm mesh
    assert f_steric_200 < 0.01


# ---------------------------------------------------------------------------
# 5. PEG Brush Shielding & Adhesive Trapping Tests
# ---------------------------------------------------------------------------
def test_peg_brush_shielding_cpp_vs_mpp() -> None:
    """Verify unshielded CPP is adhesively trapped whereas dense PEG brush rescues F_adhesive."""
    _, regime_cpp, exposed_cpp, f_adh_cpp = mpp.calculate_peg_brush_shielding(0.0)
    _, regime_mpp, exposed_mpp, f_adh_mpp = mpp.calculate_peg_brush_shielding(0.35)

    assert regime_cpp == "uncoated_cpp"
    assert regime_mpp == "dense_brush_mpp"

    # CPP: exposed area ~1.0, F_adhesive severely retarded (< 0.005)
    assert exposed_cpp == 1.0
    assert f_adh_cpp < 0.005

    # MPP: exposed core area screened (< 0.01), F_adhesive dramatically rescued (> 0.50)
    assert exposed_mpp < 0.01
    assert f_adh_mpp > 0.50
    assert (f_adh_mpp / f_adh_cpp) > 100.0


# ---------------------------------------------------------------------------
# 6. Full MPP Transport Simulation Tests
# ---------------------------------------------------------------------------
def test_simulate_mpp_transport_nominal_65nm() -> None:
    """Verify complete transport evaluation for optimized 65 nm dense brush MPP."""
    particle = mpp.NanoparticleConfig(
        diameter_nm=65.0,
        peg_mw_kda=5.0,
        grafting_density_chains_nm2=0.35,
    )
    mucus = mpp.MucusHydrogelConfig(
        mucin_volume_fraction=0.065,
        edna_conc_mg_ml=2.8,
        mucus_thickness_um=30.0,
    )

    result = mpp.simulate_mpp_transport(particle, mucus)

    assert result.particle_diameter_nm == 65.0
    assert result.brush_regime == "dense_brush_mpp"
    assert result.effective_diffusivity_ratio >= mpp.TARGET_MINIMUM_DIFFUSION_RATIO
    assert result.mucus_transit_time_minutes <= mpp.MUCOCILIARY_CLEARANCE_WINDOW_MIN
    assert result.escapes_mucociliary_clearance is True
    assert "rapid_mpp_penetrator" in result.transport_classification


def test_simulate_cpp_fails_clearance() -> None:
    """Verify conventional unshielded particle fails to escape mucociliary clearance."""
    particle = mpp.NanoparticleConfig(
        diameter_nm=65.0,
        grafting_density_chains_nm2=0.0,
    )
    mucus = mpp.MucusHydrogelConfig()
    result = mpp.simulate_mpp_transport(particle, mucus)

    assert result.escapes_mucociliary_clearance is False
    assert result.effective_diffusivity_ratio < 0.01
    assert result.mucus_transit_time_minutes > 100.0
    assert result.transport_classification == "adhesively_immobilized_cpp"


# ---------------------------------------------------------------------------
# 7. Comparative Scaling Summary Tests
# ---------------------------------------------------------------------------
def test_compare_cpp_vs_mpp_scaling() -> None:
    """Verify comparative scaling returns valid summary metrics and citations."""
    summary = mpp.compare_cpp_vs_mpp_scaling(diameter_nm=65.0)

    assert summary.particle_diameter_nm == 65.0
    assert summary.fold_diffusivity_enhancement >= 50.0
    assert summary.cpp_escapes_clearance is False
    assert summary.mpp_escapes_clearance is True
    assert summary.tier == "Tier T2 Model"

    # Clinical grounding contains required NCTs
    assert "NCT00843817" in summary.clinical_grounding
    assert "NCT05248230" in summary.clinical_grounding

    # Primary literature DOIs
    assert "10.1016/j.addr.2008.11.002" in summary.primary_dois
    assert "10.1016/j.addr.2008.09.009" in summary.primary_dois


# ---------------------------------------------------------------------------
# 8. Deterministic Reproducibility Tests
# ---------------------------------------------------------------------------
def test_simulation_determinism() -> None:
    """Verify two independent runs yield bit-for-bit identical results."""
    particle = mpp.NanoparticleConfig(diameter_nm=70.0, grafting_density_chains_nm2=0.30)
    mucus = mpp.MucusHydrogelConfig(mucin_volume_fraction=0.06, edna_conc_mg_ml=3.0)

    run_1 = mpp.simulate_mpp_transport(particle, mucus)
    run_2 = mpp.simulate_mpp_transport(particle, mucus)

    assert run_1 == run_2


# ---------------------------------------------------------------------------
# 9. CLI Execution & Parameter Sweep Tests
# ---------------------------------------------------------------------------
def test_run_size_sweep() -> None:
    """Verify diameter sweep produces ordered results across size ranges."""
    sweep = mpp.run_size_sweep([30.0, 65.0, 150.0])
    assert len(sweep) == 3
    assert sweep[0]["diameter_nm"] == 30.0
    assert sweep[1]["diameter_nm"] == 65.0
    assert sweep[2]["diameter_nm"] == 150.0

    # 30 nm and 65 nm MPP should have higher diffusion ratio than 150 nm
    assert sweep[0]["mpp_deff_ratio"] > sweep[1]["mpp_deff_ratio"] > sweep[2]["mpp_deff_ratio"]


def test_cli_execution_json() -> None:
    """Verify running the script via CLI returns valid JSON with required fields."""
    cmd = [sys.executable, str(DISCOVERY_DIR / "mpp_obstruction_scaling.py"), "--diameter", "65.0", "--json"]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=True)

    data = json.loads(proc.stdout)
    assert "simulation" in data
    assert "comparison" in data
    assert "attribution" in data
    assert data["simulation"]["particle_diameter_nm"] == 65.0


def test_cli_execution_sweep() -> None:
    """Verify running the CLI sweep mode exits with code 0."""
    cmd = [sys.executable, str(DISCOVERY_DIR / "mpp_obstruction_scaling.py"), "--run-sweep"]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
    assert proc.returncode == 0
    assert "MPP DIAMETER SWEEP" in proc.stdout


# ---------------------------------------------------------------------------
# 10. Anti-Hallucination & Provenance Guardrail Tests
# ---------------------------------------------------------------------------
def test_anti_hallucination_provenance() -> None:
    """Verify NLM disclaimer and absence of prohibited clinical overclaim language."""
    code_text = (DISCOVERY_DIR / "mpp_obstruction_scaling.py").read_text(encoding="utf-8")

    # Required NLM disclaimer verbatim
    assert "ClinicalTrials.gov does not evaluate safety or efficacy;" in code_text
    assert "listing on the registry does not denote government endorsement." in code_text

    # No prohibited clinical promises
    prohibited_phrases = [
        "proven to cure",
        "safe and effective for patients",
        "guarantees clinical benefit",
        "eliminates cystic fibrosis",
    ]
    for phrase in prohibited_phrases:
        assert phrase not in code_text.lower()
