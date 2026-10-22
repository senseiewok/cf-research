#!/usr/bin/env python3
"""Unit test suite for Inhaled LPN Airway Macrophage Evasion & Scavenger Receptor Kinetics (T-0271).

Acceptance Criteria:
  pytest tools/discovery/test_macrophage_evasion_kinetics.py -v; python tools/claims/check_anti_hallucination.py
"""
import math
import sys
from pathlib import Path
import pytest

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from macrophage_evasion_kinetics import (
    NLM_ATTRIBUTION,
    RESEARCH_BOUNDARY,
    SURFACE_CATALOG,
    AirwayMicroenvironment,
    NanoparticleSurfaceProfile,
    calculate_epithelial_penetration_rate,
    calculate_macrophage_clearance_rate,
    calculate_mucociliary_clearance_rate,
    calculate_opsonin_adsorption,
    calculate_sirp_alpha_inhibition,
    evaluate_discovery_panel,
    main,
    simulate_macrophage_evasion_kinetics,
)


class TestSurfaceCatalog:
    """Validate nanoparticle surface parameters and trial grounding."""

    def test_catalog_entries_present(self):
        assert "uncoated_core_np" in SURFACE_CATALOG
        assert "standard_peg_lnp" in SURFACE_CATALOG
        assert "dense_brush_lpn" in SURFACE_CATALOG
        assert "stealth_cd47_lpn_hybrid" in SURFACE_CATALOG

    def test_biophysical_invariants(self):
        for pid, p in SURFACE_CATALOG.items():
            assert p.peg_density_chains_nm2 >= 0.0
            assert p.cd47_peptide_density_peptides_um2 >= 0.0
            assert p.mean_diameter_nm > 50.0 and p.mean_diameter_nm < 200.0
            assert p.trial_surveillance == "NCT05668741"
            assert p.citation_doi.startswith("10.")


class TestOpsoninAdsorptionAndSIRPInhibition:
    """Validate de Gennes steric brush exclusion and CD47-SIRP-alpha Hill inhibition."""

    def test_peg_brush_reduces_opsonin_adsorption(self):
        # 0 PEG vs dense PEG (0.34 chains/nm^2)
        theta_bare = calculate_opsonin_adsorption(peg_density_chains_nm2=0.0, sp_conc_ug_ml=18.0, c3b_conc_ug_ml=40.0)
        theta_dense = calculate_opsonin_adsorption(peg_density_chains_nm2=0.34, sp_conc_ug_ml=18.0, c3b_conc_ug_ml=40.0)
        assert theta_dense < theta_bare
        assert theta_dense < 0.10  # > 90% reduction in protein corona

    def test_cd47_sirp_alpha_inhibition_hill_kinetics(self):
        # 0 density -> 0 inhibition
        inhib_zero = calculate_sirp_alpha_inhibition(0.0)
        assert inhib_zero == 0.0

        # At IC50 (150 peptides/um^2) -> exactly 50% inhibition
        inhib_ic50 = calculate_sirp_alpha_inhibition(150.0, ic50_density=150.0)
        assert pytest.approx(inhib_ic50, abs=1e-3) == 0.50

        # At high density (320 peptides/um^2) -> >75% inhibition
        inhib_high = calculate_sirp_alpha_inhibition(320.0, ic50_density=150.0)
        assert inhib_high > 0.75


class TestClearanceRateConstants:
    """Validate competing kinetic rate constants for macrophage, MCC, and penetration."""

    def test_macrophage_clearance_charge_and_cd47_dependence(self):
        # High opsonin & negative charge without CD47 -> fast clearance
        k_fast = calculate_macrophage_clearance_rate(
            opsonin_coverage=0.80,
            sirp_inhibition=0.0,
            zeta_potential_mv=-30.0,
            macrophage_density_cells_ml=2.4e6,
        )
        # Low opsonin & neutral charge with CD47 -> slow clearance
        k_slow = calculate_macrophage_clearance_rate(
            opsonin_coverage=0.05,
            sirp_inhibition=0.80,
            zeta_potential_mv=-2.0,
            macrophage_density_cells_ml=2.4e6,
        )
        assert k_fast > k_slow
        assert k_fast / k_slow > 10.0  # > 10-fold reduction in phagocytic velocity

    def test_mucus_penetration_enhancement(self):
        # Dense PEG vs bare
        k_pen_dense = calculate_epithelial_penetration_rate(
            peg_density_chains_nm2=0.32,
            mean_diameter_nm=95.0,
            mucus_thickness_um=32.0,
            mucus_mesh_pore_nm=140.0,
        )
        k_pen_bare = calculate_epithelial_penetration_rate(
            peg_density_chains_nm2=0.0,
            mean_diameter_nm=110.0,
            mucus_thickness_um=32.0,
            mucus_mesh_pore_nm=140.0,
        )
        assert k_pen_dense > k_pen_bare


class TestFullSimulationOutcomes:
    """Validate full model integration across surface catalog."""

    def test_stealth_hybrid_outperforms_uncoated(self):
        res_bare = simulate_macrophage_evasion_kinetics(SURFACE_CATALOG["uncoated_core_np"])
        res_stealth = simulate_macrophage_evasion_kinetics(SURFACE_CATALOG["stealth_cd47_lpn_hybrid"])

        assert res_stealth.evasion_half_life_hr > res_bare.evasion_half_life_hr
        assert res_stealth.epithelial_retained_fraction > res_bare.epithelial_retained_fraction
        assert res_stealth.macrophage_cleared_fraction < res_bare.macrophage_cleared_fraction
        assert "SUPERIOR_STEALTH" in res_stealth.biophysical_verdict
        assert "RAPID_CLEARANCE" in res_bare.biophysical_verdict

    def test_conservation_of_fates(self):
        for pid, p in SURFACE_CATALOG.items():
            res = simulate_macrophage_evasion_kinetics(p)
            total_fraction = (
                res.macrophage_cleared_fraction
                + res.mucociliary_cleared_fraction
                + res.epithelial_retained_fraction
            )
            assert pytest.approx(total_fraction, abs=1e-3) == 1.0


class TestDiscoveryPanelAndCLI:
    """Validate discovery panel synthesis and CLI verification."""

    def test_evaluate_discovery_panel(self):
        panel = evaluate_discovery_panel()
        assert panel["status"] == "PASS"
        assert panel["panel_count"] == 4
        assert NLM_ATTRIBUTION in panel["nlm_attribution"]
        assert RESEARCH_BOUNDARY in panel["research_boundary"]

    def test_main_cli_verify(self, monkeypatch, capsys):
        monkeypatch.setattr(sys, "argv", ["macrophage_evasion_kinetics.py", "--verify"])
        code = main()
        assert code == 0
        captured = capsys.readouterr()
        assert "PASS" in captured.out
        assert "stealth_cd47_lpn_hybrid" in captured.out
