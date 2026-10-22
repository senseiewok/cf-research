#!/usr/bin/env python3
"""Unit test suite for Airway Ciliated Epithelial Glycocalyx Permeabilization Model (T-0291).

Acceptance Criteria:
  pytest tools/discovery/test_glycocalyx_permeabilization_model.py -v; python tools/claims/check_anti_hallucination.py
"""
import math
import sys
from pathlib import Path
import pytest

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from glycocalyx_permeabilization_model import (
    ENZYME_CATALOG,
    NANOPARTICLE_CATALOG,
    NLM_ATTRIBUTION,
    RESEARCH_BOUNDARY,
    EnzymeTreatmentProfile,
    GlycocalyxBrushParameters,
    NanoparticleProfile,
    calculate_brush_thickness_collapse,
    calculate_de_gennes_mesh_dilation,
    calculate_electrostatic_partition,
    calculate_ogston_sieving,
    calculate_stokes_einstein_diffusivity,
    evaluate_discovery_panel,
    main,
    simulate_glycocalyx_permeabilization,
)


class TestCatalogAndGrounding:
    """Validate nanoparticle profiles, enzyme profiles, and regulatory disclaimers."""

    def test_nanoparticle_catalog_entries(self):
        assert "LNP_mRNA_BEACON_benchmark" in NANOPARTICLE_CATALOG
        assert "AAV2_capsid_reference" in NANOPARTICLE_CATALOG
        assert "AAV5_capsid_reference" in NANOPARTICLE_CATALOG
        assert "Polymeric_NP_PEGylated" in NANOPARTICLE_CATALOG

    def test_regulatory_and_boundary_statements(self):
        assert "ClinicalTrials.gov does not evaluate safety or efficacy" in NLM_ATTRIBUTION
        assert "Preclinical in-silico biophysical modeling" in RESEARCH_BOUNDARY
        # Ensure no prohibited overclaim phrases exist in attribution
        prohibited = ["proven to cure", "guarantees clinical benefit", "cures cystic fibrosis"]
        for p in prohibited:
            assert p not in NLM_ATTRIBUTION.lower()
            assert p not in RESEARCH_BOUNDARY.lower()

    def test_lnp_beacon_grounding(self):
        lnp = NANOPARTICLE_CATALOG["LNP_mRNA_BEACON_benchmark"]
        assert lnp.clinical_trial_nct == "NCT05668741"
        assert lnp.evidence_doi == "10.1038/nbt.3298"
        assert 80.0 <= lnp.diameter_nm <= 120.0
        assert lnp.zeta_potential_mv < 0.0

    def test_enzyme_catalog_entries(self):
        assert "Sialidase_therapeutic" in ENZYME_CATALOG
        assert "Heparinase_therapeutic" in ENZYME_CATALOG
        assert "Dual_Sialidase_Heparinase" in ENZYME_CATALOG
        for name, enz in ENZYME_CATALOG.items():
            assert enz.concentration_u_ml > 0.0
            assert enz.kcat_per_s > 0.0
            assert enz.km_um > 0.0


class TestBiophysicalMechanisms:
    """Validate core mathematical formulas for sieving, electrostatics, and polymer scaling."""

    def test_stokes_einstein_diffusivity_scaling(self):
        d_small = calculate_stokes_einstein_diffusivity(25.0)
        d_large = calculate_stokes_einstein_diffusivity(100.0)
        assert d_small > d_large
        # 4x diameter difference should yield 4x inverse diffusivity
        ratio = d_small / d_large
        assert math.isclose(ratio, 4.0, rel_tol=1e-3)

    def test_ogston_sieving_monotonicity(self):
        mesh_nm = 30.0
        fiber_nm = 3.5
        # AAV (25 nm) vs LNP (95 nm)
        k_aav = calculate_ogston_sieving(25.0, mesh_nm, fiber_nm)
        k_lnp = calculate_ogston_sieving(95.0, mesh_nm, fiber_nm)
        assert 0.0 <= k_lnp < k_aav <= 1.0
        # When mesh expands, sieving partition increases
        k_lnp_dilated = calculate_ogston_sieving(95.0, 70.0, fiber_nm)
        assert k_lnp_dilated > k_lnp

    def test_electrostatic_partition_repulsion(self):
        # Like negative charges repel (K_elec < 1.0)
        k_repel = calculate_electrostatic_partition(
            particle_zeta_mv=-20.0,
            glycocalyx_potential_mv=-30.0,
        )
        # Neutral particle experiences no electrostatic penalty
        k_neutral = calculate_electrostatic_partition(
            particle_zeta_mv=0.0,
            glycocalyx_potential_mv=-30.0,
        )
        assert k_repel < k_neutral
        assert math.isclose(k_neutral, 1.0, rel_tol=1e-5)

    def test_de_gennes_mesh_dilation(self):
        xi_0 = 28.0
        xi_undigested = calculate_de_gennes_mesh_dilation(xi_0, 0.0, 0.0)
        assert math.isclose(xi_undigested, xi_0, rel_tol=1e-5)

        xi_digested = calculate_de_gennes_mesh_dilation(xi_0, 0.8, 0.9)
        assert xi_digested > xi_0 * 1.5

    def test_brush_thickness_collapse(self):
        h_0 = 350.0
        h_collapsed = calculate_brush_thickness_collapse(h_0, 0.9, 0.8)
        assert h_collapsed < h_0
        assert h_collapsed >= h_0 * 0.20


class TestSimulationAndPanel:
    """Validate dynamic trajectory simulation and multi-vehicle comparison panel."""

    def test_lnp_permeabilization_enhancement(self):
        lnp = NANOPARTICLE_CATALOG["LNP_mRNA_BEACON_benchmark"]
        sim = simulate_glycocalyx_permeabilization(
            nanoparticle=lnp,
            enzyme_regimen=ENZYME_CATALOG["Sialidase_therapeutic"],
            total_time_minutes=45.0,
        )
        assert len(sim) > 1
        t_init = sim[0]
        t_final = sim[-1]

        # Mesh expands and thickness collapses
        assert t_final.mesh_spacing_nm > t_init.mesh_spacing_nm
        assert t_final.brush_thickness_nm < t_init.brush_thickness_nm
        assert t_final.sialic_acid_cleavage_fraction > 0.50

        # Collision frequency increases
        assert t_final.fold_enhancement_vs_baseline > 1.0
        assert t_final.apical_collision_frequency_m2_s >= t_init.apical_collision_frequency_m2_s

    def test_dual_enzyme_synergy(self):
        lnp = NANOPARTICLE_CATALOG["LNP_mRNA_BEACON_benchmark"]
        sim_sial = simulate_glycocalyx_permeabilization(
            nanoparticle=lnp,
            enzyme_regimen=ENZYME_CATALOG["Sialidase_therapeutic"],
            total_time_minutes=30.0,
        )
        sim_dual = simulate_glycocalyx_permeabilization(
            nanoparticle=lnp,
            enzyme_regimen=ENZYME_CATALOG["Dual_Sialidase_Heparinase"],
            total_time_minutes=30.0,
        )
        # Dual enzyme should produce greater mesh dilation due to GAG core depolymerization
        assert sim_dual[-1].mesh_spacing_nm > sim_sial[-1].mesh_spacing_nm

    def test_discovery_panel_generation(self):
        panel = evaluate_discovery_panel(duration_min=30.0)
        assert panel["metadata"]["task"] == "T-0291"
        assert "BEACON-CF Inhaled CFTR mRNA Lipid Nanoparticle" in panel["profiles"]
        assert "Recombinant Adeno-Associated Virus Serotype 2 (AAV2)" in panel["profiles"]

        lnp_prof = panel["profiles"]["BEACON-CF Inhaled CFTR mRNA Lipid Nanoparticle"]
        assert lnp_prof["fold_enhancement"] > 1.0
        assert lnp_prof["clinical_trial_nct"] == "NCT05668741"

    def test_cli_execution(self, capsys):
        ret = main(["--duration-min", "10", "--particle", "AAV2_capsid_reference"])
        assert ret == 0
        captured = capsys.readouterr().out
        assert "AIRWAY GLYCOCALYX PERMEABILIZATION" in captured
        assert "xi_final" in captured

    def test_cli_json_mode(self, capsys):
        ret = main(["--json", "--duration-min", "15"])
        assert ret == 0
        captured = capsys.readouterr().out
        import json
        data = json.loads(captured)
        assert data["metadata"]["task"] == "T-0291"
