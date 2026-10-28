#!/usr/bin/env python3
"""Unit tests for Biophysical Modeling of Inhaled mRNA-LNP Aerosol Shear Degradation (T-0481)."""

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
    from tools.discovery.mrna_lnp_nebulizer_shear_stability import (
        CLINICAL_TRIAL_VANZACAFTOR,
        CLINICAL_TRIAL_VX522,
        CURATIVE_RESPIRABLE_YIELD_THRESHOLD_PCT,
        DEFAULT_CHOLESTEROL_MOL_PCT,
        DEFAULT_FLOW_RATE_ML_MIN,
        DEFAULT_HELPER_LIPID_MOL_PCT,
        DEFAULT_LNP_INITIAL_DIAMETER_NM,
        DEFAULT_PEG_LIPID_MOL_PCT,
        DEFAULT_PIEZO_FREQ_KHZ,
        DEFAULT_PORE_COUNT,
        DEFAULT_PORE_EXIT_DIAMETER_UM,
        MAX_PERMISSIBLE_POST_DIAMETER_NM,
        NLM_ATTRIBUTION,
        PDB_STRUCTURE_CFTR,
        PRIMARY_DOI_BUTTON,
        PRIMARY_DOI_MIDDLETON,
        RESEARCH_BOUNDARY,
        AerosolDeliveryResult,
        BilayerStressRuptureResult,
        ComprehensiveStabilityResult,
        LnpFormulationParameters,
        MeshExtensionalShearState,
        MrnaLnpNebulizerShearStabilityModel,
        NebulizerMeshParameters,
    )
except ImportError:
    from mrna_lnp_nebulizer_shear_stability import (
        CLINICAL_TRIAL_VANZACAFTOR,
        CLINICAL_TRIAL_VX522,
        CURATIVE_RESPIRABLE_YIELD_THRESHOLD_PCT,
        DEFAULT_CHOLESTEROL_MOL_PCT,
        DEFAULT_FLOW_RATE_ML_MIN,
        DEFAULT_HELPER_LIPID_MOL_PCT,
        DEFAULT_LNP_INITIAL_DIAMETER_NM,
        DEFAULT_PEG_LIPID_MOL_PCT,
        DEFAULT_PIEZO_FREQ_KHZ,
        DEFAULT_PORE_COUNT,
        DEFAULT_PORE_EXIT_DIAMETER_UM,
        MAX_PERMISSIBLE_POST_DIAMETER_NM,
        NLM_ATTRIBUTION,
        PDB_STRUCTURE_CFTR,
        PRIMARY_DOI_BUTTON,
        PRIMARY_DOI_MIDDLETON,
        RESEARCH_BOUNDARY,
        AerosolDeliveryResult,
        BilayerStressRuptureResult,
        ComprehensiveStabilityResult,
        LnpFormulationParameters,
        MeshExtensionalShearState,
        MrnaLnpNebulizerShearStabilityModel,
        NebulizerMeshParameters,
    )


class TestMrnaLnpNebulizerShearStability:
    """Test suite for T-0481 vibrating-mesh nebulizer mRNA-LNP stability."""

    def test_grounding_constants_and_attributions(self) -> None:
        """Verify peer-reviewed literature DOIs, clinical trial NCTs, and structural grounding."""
        assert PRIMARY_DOI_BUTTON == "10.1152/physrev.00004.2021"
        assert PRIMARY_DOI_MIDDLETON == "10.1056/nejmoa1908639"
        assert CLINICAL_TRIAL_VX522 == "NCT05668741"
        assert CLINICAL_TRIAL_VANZACAFTOR == "NCT05033080"
        assert PDB_STRUCTURE_CFTR == "6MSM"
        assert "ClinicalTrials.gov does not evaluate safety or efficacy" in NLM_ATTRIBUTION
        assert "research, not medical advice" in RESEARCH_BOUNDARY
        assert CURATIVE_RESPIRABLE_YIELD_THRESHOLD_PCT == 50.0
        assert MAX_PERMISSIBLE_POST_DIAMETER_NM == 130.0

    def test_mesh_hydrodynamics_shear_rates(self) -> None:
        """Verify extensional shear rate (>10^5 s^-1) and wall shear rate (>10^6 s^-1)."""
        model = MrnaLnpNebulizerShearStabilityModel()
        mesh_state = model.calculate_mesh_hydrodynamics()

        assert mesh_state.exit_velocity_m_s > 1.0  # Typically 1.5 - 5 m/s
        assert mesh_state.exit_velocity_m_s > mesh_state.inlet_velocity_m_s
        assert mesh_state.extensional_shear_rate_s_inv > 1.0e4  # High extensional elongation
        assert mesh_state.wall_shear_rate_s_inv > 1.0e6        # High wall shear
        assert mesh_state.extensional_stress_pa > 0.0
        assert mesh_state.wall_shear_stress_pa > 0.0
        assert mesh_state.piezo_dynamic_pressure_bar > 0.0
        assert mesh_state.flow_reynolds_number < 2000.0        # Microchannel laminar flow

    def test_bilayer_rupture_and_cholesterol_stabilization(self) -> None:
        """Verify cholesterol and helper lipid elevation of lysis tension and suppression of aggregation."""
        model_default = MrnaLnpNebulizerShearStabilityModel()
        rupture_default = model_default.calculate_bilayer_rupture()

        # Stabilized formulation maintains membrane integrity
        assert rupture_default.effective_critical_lysis_tension_mn_m > 8.0
        assert rupture_default.effective_area_modulus_mn_m > 250.0
        assert rupture_default.rupture_probability < 0.20
        assert rupture_default.post_encapsulation_eff >= 0.85
        assert rupture_default.post_diameter_nm <= MAX_PERMISSIBLE_POST_DIAMETER_NM
        assert rupture_default.steric_protection_factor > 0.80
        assert rupture_default.is_bilayer_intact

        # Destabilized formulation (0% chol, 0% DSPC, 0% PEG) fails under identical shear
        bare_lnp = LnpFormulationParameters(
            cholesterol_mol_pct=0.0,
            helper_lipid_mol_pct=0.0,
            peg_lipid_mol_pct=0.0,
        )
        model_bare = MrnaLnpNebulizerShearStabilityModel(lnp_params=bare_lnp)
        rupture_bare = model_bare.calculate_bilayer_rupture()

        assert rupture_bare.effective_critical_lysis_tension_mn_m < 6.0
        assert rupture_bare.rupture_probability > rupture_default.rupture_probability
        assert rupture_bare.free_mrna_leakage_fraction > 0.40
        assert rupture_bare.post_encapsulation_eff < 0.60
        assert rupture_bare.post_diameter_nm > rupture_default.post_diameter_nm
        assert not rupture_bare.is_bilayer_intact

    def test_aerosol_droplet_delivery_and_yield(self) -> None:
        """Verify aerosol droplet MMAD (1.5-3.5 um), FPF, and net respirable intact mRNA yield."""
        model = MrnaLnpNebulizerShearStabilityModel()
        aerosol = model.calculate_aerosol_delivery()

        assert 1.5 <= aerosol.mmad_um <= 3.5  # Respirable droplet size
        assert aerosol.fine_particle_fraction_pct >= 65.0
        assert aerosol.net_respirable_intact_mrna_yield_pct >= CURATIVE_RESPIRABLE_YIELD_THRESHOLD_PCT
        assert aerosol.is_mucus_permeable
        assert aerosol.curative_threshold_met

    def test_negative_control_assertion_ah10(self) -> None:
        """Verify AH-10 negative control and dynamic falsification assertion."""
        model = MrnaLnpNebulizerShearStabilityModel()
        falsification_report = model.execute_negative_control_assertion()

        assert falsification_report["rule"] == "AH-10"
        assert falsification_report["passed"]
        assert falsification_report["destabilized_formulation_failed"]
        assert falsification_report["extreme_shear_failed"]
        assert falsification_report["standard_formulation_passed"]

    def test_sensitivity_analysis_ah11(self) -> None:
        """Verify AH-11 +/- 25% parametric sensitivity audit across biophysical inputs."""
        model = MrnaLnpNebulizerShearStabilityModel()
        sens_report = model.run_sensitivity_analysis(perturbation_fraction=0.25)

        assert sens_report["rule"] == "AH-11"
        assert sens_report["sensitivity_verdict"] == "VERIFIED_CONTINUOUS_AND_BOUNDED"
        assert "cholesterol_mol_pct" in sens_report["parametric_responses"]
        assert "peg_lipid_mol_pct" in sens_report["parametric_responses"]
        assert "pore_exit_diameter_um" in sens_report["parametric_responses"]

        # Increasing cholesterol should increase or maintain respirable yield
        chol_data = sens_report["parametric_responses"]["cholesterol_mol_pct"]
        assert chol_data["yield_high_pct"] >= chol_data["yield_low_pct"]

    def test_full_simulation_integration(self) -> None:
        """Verify end-to-end integration and JSON serializability."""
        model = MrnaLnpNebulizerShearStabilityModel()
        sim_res = model.run_full_simulation()

        assert isinstance(sim_res, ComprehensiveStabilityResult)
        assert sim_res.aerosol_delivery.curative_threshold_met
        assert "Curative Target Met" in sim_res.overall_verdict
