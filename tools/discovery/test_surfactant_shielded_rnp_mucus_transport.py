#!/usr/bin/env python3
"""Unit tests for Surfactant-Shielded RNP Mucus Transport and Delivery Model (T-0482)."""

from __future__ import annotations

import json
import sys
from pathlib import Path
import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
DISCOVERY_DIR = Path(__file__).resolve().parent
if str(DISCOVERY_DIR) not in sys.path:
    sys.path.insert(0, str(DISCOVERY_DIR))

try:
    from tools.discovery.surfactant_shielded_rnp_mucus_transport import (
        CLINICAL_TRIAL_VANZACAFTOR,
        CLINICAL_TRIAL_VX522,
        CURATIVE_THRESHOLD_CFTR_PCT,
        DEFAULT_BARE_RNP_RADIUS_NM,
        DEFAULT_BARE_ZETA_POTENTIAL_MV,
        DEFAULT_CF_MUCUS_PORE_SIZE_NM,
        DEFAULT_CORONA_ZETA_POTENTIAL_MV,
        DEFAULT_SURFACTANT_CONC_MG_ML,
        NLM_ATTRIBUTION,
        PDB_STRUCTURE_CFTR,
        PRIMARY_DOI_BUTTON,
        PRIMARY_DOI_COONEY,
        PRIMARY_DOI_MIDDLETON,
        RESEARCH_BOUNDARY,
        IntracellularKineticsParameters,
        MucusMeshParameters,
        SurfactantCoronaParameters,
        SurfactantShieldedRnpMucusTransportModel,
    )
except ImportError:
    from surfactant_shielded_rnp_mucus_transport import (
        CLINICAL_TRIAL_VANZACAFTOR,
        CLINICAL_TRIAL_VX522,
        CURATIVE_THRESHOLD_CFTR_PCT,
        DEFAULT_BARE_RNP_RADIUS_NM,
        DEFAULT_BARE_ZETA_POTENTIAL_MV,
        DEFAULT_CF_MUCUS_PORE_SIZE_NM,
        DEFAULT_CORONA_ZETA_POTENTIAL_MV,
        DEFAULT_SURFACTANT_CONC_MG_ML,
        NLM_ATTRIBUTION,
        PDB_STRUCTURE_CFTR,
        PRIMARY_DOI_BUTTON,
        PRIMARY_DOI_COONEY,
        PRIMARY_DOI_MIDDLETON,
        RESEARCH_BOUNDARY,
        IntracellularKineticsParameters,
        MucusMeshParameters,
        SurfactantCoronaParameters,
        SurfactantShieldedRnpMucusTransportModel,
    )


class TestSurfactantShieldedRnpMucusTransport:
    """Test suite for T-0482 surfactant-shielded base-editor RNP delivery."""

    def test_grounding_constants_and_attributions(self) -> None:
        """Verify peer-reviewed literature DOIs, clinical trial NCTs, and structural grounding."""
        assert PRIMARY_DOI_BUTTON == "10.1152/physrev.00004.2021"
        assert PRIMARY_DOI_COONEY == "10.1016/j.ymthe.2018.04.014"
        assert PRIMARY_DOI_MIDDLETON == "10.1056/nejmoa1908639"
        assert CLINICAL_TRIAL_VX522 == "NCT05668741"
        assert CLINICAL_TRIAL_VANZACAFTOR == "NCT05033080"
        assert PDB_STRUCTURE_CFTR == "6MSM"
        assert "ClinicalTrials.gov does not evaluate safety or efficacy" in NLM_ATTRIBUTION
        assert "research, not medical advice" in RESEARCH_BOUNDARY
        assert CURATIVE_THRESHOLD_CFTR_PCT == 10.0

    def test_surfactant_corona_adsorption_equilibrium(self) -> None:
        """Verify Hill cooperative adsorption isotherm and surface charge remodeling."""
        model = SurfactantShieldedRnpMucusTransportModel()

        # Zero surfactant concentration -> zero coverage, bare cationic charge
        corona_zero = model.calculate_surfactant_corona(conc_mg_ml=0.0)
        assert corona_zero.fractional_coverage == 0.0
        assert corona_zero.effective_zeta_potential_mv == DEFAULT_BARE_ZETA_POTENTIAL_MV
        assert not corona_zero.is_electrostatically_shielded
        assert corona_zero.effective_hydrodynamic_radius_nm == DEFAULT_BARE_RNP_RADIUS_NM

        # Standard physiological surfactant concentration (2.0 mg/mL)
        corona_std = model.calculate_surfactant_corona(conc_mg_ml=2.0)
        assert corona_std.fractional_coverage > 0.85
        assert corona_std.is_electrostatically_shielded
        assert corona_std.effective_zeta_potential_mv < 0.0  # Shielded to near-neutral / negative
        assert corona_std.effective_hydrodynamic_radius_nm > DEFAULT_BARE_RNP_RADIUS_NM
        assert corona_std.adsorption_free_energy_kj_mol < 0.0  # Spontaneous thermodynamic adsorption

        # High saturation concentration (10.0 mg/mL) -> asymptotic coverage
        corona_sat = model.calculate_surfactant_corona(conc_mg_ml=10.0)
        assert corona_sat.fractional_coverage >= corona_std.fractional_coverage
        assert corona_sat.effective_zeta_potential_mv <= DEFAULT_CORONA_ZETA_POTENTIAL_MV + 1.0

    def test_cf_mucus_diffusion_and_sieving(self) -> None:
        """Verify Stokes-Einstein diffusivity, Ogston steric sieving, and electrostatic penalty."""
        model = SurfactantShieldedRnpMucusTransportModel()

        # Test shielded carrier diffusion
        corona_shielded = model.calculate_surfactant_corona(conc_mg_ml=2.0)
        mucus_shielded = model.calculate_mucus_transport(corona_shielded)

        assert mucus_shielded.diffusivity_water_cm2_s > 1e-8
        assert 0.0 < mucus_shielded.steric_sieving_factor <= 1.0
        assert mucus_shielded.electrostatic_transmission_factor == 1.0  # Non-cationic suffers no penalty
        assert mucus_shielded.diffusivity_ratio > 0.05
        assert mucus_shielded.mucus_penetration_probability > 0.30

        # Test unshielded cationic carrier (+26 mV) suffers severe electrostatic immobilization
        corona_bare = model.calculate_surfactant_corona(conc_mg_ml=0.0)
        mucus_bare = model.calculate_mucus_transport(corona_bare)

        assert mucus_bare.electrostatic_transmission_factor < 0.05
        assert mucus_bare.diffusivity_ratio < 0.01
        assert mucus_bare.mucus_penetration_probability < 0.01
        assert mucus_bare.transport_regime == "ENTRAPPED_STERIC_ELECTROSTATIC"

    def test_intracellular_uncoating_and_basal_editing(self) -> None:
        """Verify pH-triggered endosomal uncoating and allelic editing in basal stem cells."""
        model = SurfactantShieldedRnpMucusTransportModel()
        sim_res = model.run_full_simulation()

        editing = sim_res.intracellular_editing
        assert editing.uncoating_fraction > 0.70  # pH 5.5 triggers uncoating
        assert editing.effective_endosomal_escape > 0.08
        assert editing.nuclear_translocation_fraction > 0.70
        assert editing.net_allelic_editing_efficiency > 0.10
        assert editing.restored_cftr_pct_wt >= CURATIVE_THRESHOLD_CFTR_PCT
        assert editing.curative_threshold_met is True
        assert sim_res.verdict == "CURATIVE_RESCUE_ACHIEVED"

    def test_ah10_negative_control_assertion(self) -> None:
        """Assert AH-10 negative control falsification passes cleanly."""
        model = SurfactantShieldedRnpMucusTransportModel()
        neg_ctrl = model.execute_negative_control_assertion()

        assert neg_ctrl["rule"] == "AH-10"
        assert neg_ctrl["passed"] is True
        assert neg_ctrl["unshielded_cationic_entrapped"] is True
        assert neg_ctrl["unshielded_diffusivity_ratio"] < 0.02
        assert neg_ctrl["unshielded_restored_cftr_pct"] < 2.0
        assert neg_ctrl["oversized_carrier_excluded"] is True
        assert neg_ctrl["oversized_steric_factor"] < 0.05
        assert neg_ctrl["shielded_curative_restored_cftr_pct"] >= CURATIVE_THRESHOLD_CFTR_PCT

    def test_ah11_parametric_sensitivity(self) -> None:
        """Verify +/-25% parametric sensitivity perturbation maintains monotonic physical bounds."""
        model = SurfactantShieldedRnpMucusTransportModel()
        sens = model.run_sensitivity_analysis(perturbation_fraction=0.25)

        assert sens["rule"] == "AH-11"
        assert sens["parametric_stability_verified"] is True
        assert sens["baseline_restored_cftr_pct"] >= CURATIVE_THRESHOLD_CFTR_PCT

        responses = sens["parametric_responses"]
        assert "surfactant_conc_mg_ml" in responses
        assert "pore_size_nm" in responses
        assert "mucus_thickness_um" in responses

        # Increasing pore size must increase delivery
        assert responses["pore_size_nm"]["delta_response"] > 0
        # Increasing mucus thickness must decrease delivery
        assert responses["mucus_thickness_um"]["delta_response"] < 0

    def test_cli_execution_flags(self) -> None:
        """Verify CLI driver support for --json, --falsify, and --sensitivity."""
        import subprocess

        script = SCRIPT_PATH = Path(__file__).resolve().parent / "surfactant_shielded_rnp_mucus_transport.py"

        # Test --json flag
        res_json = subprocess.run(
            [sys.executable, str(script), "--json"],
            capture_output=True,
            text=True,
            check=True,
        )
        data = json.loads(res_json.stdout)
        assert data["task_id"] == "T-0482"
        assert data["verdict"] == "CURATIVE_RESCUE_ACHIEVED"

        # Test --falsify flag
        res_falsify = subprocess.run(
            [sys.executable, str(script), "--falsify", "--json"],
            capture_output=True,
            text=True,
            check=True,
        )
        falsify_data = json.loads(res_falsify.stdout)
        assert falsify_data["passed"] is True

        # Test --sensitivity flag
        res_sens = subprocess.run(
            [sys.executable, str(script), "--sensitivity", "--json"],
            capture_output=True,
            text=True,
            check=True,
        )
        sens_data = json.loads(res_sens.stdout)
        assert sens_data["parametric_stability_verified"] is True
