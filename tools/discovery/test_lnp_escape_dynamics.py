#!/usr/bin/env python3
"""Unit test suite for Airway Basal Stem Cell LNP Endosomal Escape Dynamics Model (T-0285).

Acceptance Criteria:
  pytest tools/discovery/test_lnp_escape_dynamics.py -v; python tools/claims/check_anti_hallucination.py
"""
import math
import sys
from pathlib import Path
import pytest

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from lnp_escape_dynamics import (
    FORMULATION_CATALOG,
    NLM_ATTRIBUTION,
    RESEARCH_BOUNDARY,
    BasalStemCellParameters,
    IonizableLipidProfile,
    LNPFormulation,
    calculate_endosomal_ph,
    calculate_hii_escape_rate,
    calculate_lipid_protonation,
    calculate_osmotic_swelling_kpa,
    evaluate_discovery_panel,
    main,
    simulate_lnp_escape_dynamics,
)


class TestFormulationCatalog:
    """Validate formulation parameters, biophysical invariants, and trial grounding."""

    def test_catalog_entries_present(self):
        assert "BEACON_CF_benchmark" in FORMULATION_CATALOG
        assert "DOPE_fusogenic_stem_optimized" in FORMULATION_CATALOG
        assert "ALC0315_standard" in FORMULATION_CATALOG
        assert "SM102_standard" in FORMULATION_CATALOG

    def test_biophysical_invariants(self):
        for fid, form in FORMULATION_CATALOG.items():
            assert 6.2 <= form.ionizable_lipid.apparent_pka <= 6.8, f"{fid} pKa outside physiological range"
            assert form.mean_diameter_nm > 50.0 and form.mean_diameter_nm < 150.0
            assert form.trial_surveillance == "NCT05668741"
            assert form.citation_doi.startswith("10.")
            assert form.h_ii_propensity_factor > 0.0


class TestProtonationAndAcidification:
    """Validate Henderson-Hasselbalch ionization and V-ATPase acidification kinetics."""

    def test_lipid_protonation_at_neutral_and_acidic_ph(self):
        pka = 6.48
        # At pH == pKa, exactly 50% protonated
        alpha_half = calculate_lipid_protonation(ph=6.48, apparent_pka=pka)
        assert pytest.approx(alpha_half, rel=1e-3) == 0.50

        # At neutral extracellular pH 7.4, minimal protonation (< 15%)
        alpha_neutral = calculate_lipid_protonation(ph=7.40, apparent_pka=pka)
        assert alpha_neutral < 0.15

        # At late endosomal pH 5.0, high protonation (> 90%)
        alpha_acidic = calculate_lipid_protonation(ph=5.00, apparent_pka=pka)
        assert alpha_acidic > 0.90

    def test_endosomal_acidification_trajectory(self):
        # At t=0, pH is extracellular 7.4
        ph_0 = calculate_endosomal_ph(0.0)
        assert pytest.approx(ph_0, abs=1e-2) == 7.40

        # Over time, pH drops toward lysosomal pH 4.8
        ph_30 = calculate_endosomal_ph(30.0)
        ph_120 = calculate_endosomal_ph(120.0)
        assert ph_30 < ph_0
        assert ph_120 < ph_30
        assert ph_120 >= 4.80


class TestOsmoticPressureAndEscapeRates:
    """Validate proton-sponge osmotic swelling and H_II phase transition escape rates."""

    def test_osmotic_swelling_monotonicity(self):
        p_low = calculate_osmotic_swelling_kpa(alpha_ion=0.10, charge_ratio_n_p=6.0, particle_radius_nm=40.0)
        p_high = calculate_osmotic_swelling_kpa(alpha_ion=0.90, charge_ratio_n_p=6.0, particle_radius_nm=40.0)
        assert p_high > p_low
        assert p_high > 15.0  # Significant osmotic pressure in kPa

    def test_dope_vs_dspc_escape_rates(self):
        # DOPE has high H_II propensity, DSPC has low H_II propensity
        k_dope = calculate_hii_escape_rate(alpha_ion=0.85, h_ii_propensity=2.45, osmotic_swelling_kpa=25.0)
        k_dspc = calculate_hii_escape_rate(alpha_ion=0.85, h_ii_propensity=0.65, osmotic_swelling_kpa=25.0)
        assert k_dope > k_dspc


class TestSimulationDynamics:
    """Validate full time-course simulation of endosomal escape in basal stem cells."""

    def test_bioavailability_within_literature_range(self):
        # Wittrup et al. 2015 and Gilleron et al. 2013: 1% to 5% escape fraction
        for fid, form in FORMULATION_CATALOG.items():
            res = simulate_lnp_escape_dynamics(form)
            assert 0.010 <= res.bioavailability_fraction_phi <= 0.055, f"{fid} escape fraction outside 1-5% range"
            assert 1.0 <= res.translation_onset_half_life_hr <= 4.5, f"{fid} translation onset outside 1-4.5 hr range"
            assert res.max_gal3_puncta_per_cell > 0.5, "Expected detectable Gal3 puncta formation"
            assert res.trial_surveillance == "NCT05668741"

    def test_fusogenic_outperforms_benchmark(self):
        res_beacon = simulate_lnp_escape_dynamics(FORMULATION_CATALOG["BEACON_CF_benchmark"])
        res_dope = simulate_lnp_escape_dynamics(FORMULATION_CATALOG["DOPE_fusogenic_stem_optimized"])
        assert res_dope.bioavailability_fraction_phi > res_beacon.bioavailability_fraction_phi
        assert "OPTIMAL_FUSOGENIC" in res_dope.biophysical_verdict

    def test_trajectory_conservation(self):
        form = FORMULATION_CATALOG["DOPE_fusogenic_stem_optimized"]
        res = simulate_lnp_escape_dynamics(form)
        assert len(res.trajectory) > 50
        # Check non-negativity across all states
        for pt in res.trajectory:
            assert pt.early_endosome_fraction >= 0.0
            assert pt.late_endosome_fraction >= 0.0
            assert pt.lysosome_fraction >= 0.0
            assert pt.recycled_fraction >= 0.0
            assert pt.cytosolic_intact_fraction >= 0.0


class TestDiscoveryPanelAndCLI:
    """Validate discovery panel synthesis and CLI verification."""

    def test_evaluate_discovery_panel(self):
        panel = evaluate_discovery_panel()
        assert panel["status"] == "PASS"
        assert panel["panel_count"] == 4
        assert NLM_ATTRIBUTION in panel["nlm_attribution"]
        assert RESEARCH_BOUNDARY in panel["research_boundary"]

    def test_main_cli_verify(self, monkeypatch, capsys):
        monkeypatch.setattr(sys, "argv", ["lnp_escape_dynamics.py", "--verify"])
        code = main()
        assert code == 0
        captured = capsys.readouterr()
        assert "PASS" in captured.out
        assert "BEACON_CF_benchmark" in captured.out
