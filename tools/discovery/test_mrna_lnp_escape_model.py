#!/usr/bin/env python3
"""Unit tests for Inhaled CFTR mRNA Lipid Nanoparticle (LNP) Endosomal Escape Kinetics Model (T-0239).

Verifies pH-dependent ionizable lipid protonation, surface zeta potential inversion,
mucin penetration mechanics, non-bilayer H_II phase escape kinetics, vibrating mesh
nebulizer shear stress retention, and comparison against BEACON-CF (NCT05668741).
"""

from __future__ import annotations

import sys
from pathlib import Path
import pytest

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO_ROOT))

from mrna_lnp_escape_model import (
    MRNALNPEscapeModel,
    compute_protonation_fraction,
    compute_surface_charge_zeta,
    compute_mucin_penetration_index,
    compute_escape_rate,
    compute_mesh_shear_stress,
    compute_nebulization_retention,
    simulate_endosomal_trafficking,
    LNP_CATALOG,
    NLM_ATTRIBUTION,
    RESEARCH_BOUNDARY,
)


class TestMRNALNPEscapeModel:
    @pytest.fixture
    def model(self) -> MRNALNPEscapeModel:
        return MRNALNPEscapeModel()

    def test_protonation_fraction_bounds_and_symmetry(self) -> None:
        """Verify Henderson-Hasselbalch protonation fraction theta_prot."""
        # Exact equivalence at pH == pKa
        assert pytest.approx(compute_protonation_fraction(6.45, 6.45), 1e-5) == 0.50

        # Acidic extreme (pH << pKa) -> 1.0
        assert compute_protonation_fraction(2.0, 6.45) > 0.9999

        # Basic extreme (pH >> pKa) -> 0.0
        assert compute_protonation_fraction(11.0, 6.45) < 0.0001

        # Extracellular airway surface (pH 7.4, pKa 6.45) -> 0.08 - 0.12
        theta_airway = compute_protonation_fraction(7.4, 6.45)
        assert 0.08 <= theta_airway <= 0.12

        # Late acidifying endosome (pH 5.5, pKa 6.45) -> 0.85 - 0.95
        theta_endo = compute_protonation_fraction(5.5, 6.45)
        assert 0.85 <= theta_endo <= 0.95

    def test_surface_charge_zeta_inversion(self) -> None:
        """Surface charge transitions from near-neutral in airway to cationic in endosome."""
        zeta_airway = compute_surface_charge_zeta(ph=7.4, pka=6.45)
        zeta_endo = compute_surface_charge_zeta(ph=5.5, pka=6.45)

        # Extracellular airway surface charge should be near-neutral: [-5.0, +5.0] mV
        assert -5.0 <= zeta_airway <= 5.0

        # Late endosome surface charge should be strongly cationic: > +18.0 mV
        assert zeta_endo >= 18.0

        # Strong inversion delta > 18 mV
        assert (zeta_endo - zeta_airway) >= 18.0

    def test_mucin_penetration_index(self) -> None:
        """Near-neutral nanoparticles penetrate CF mucus; cationic particles adhere."""
        # Near-neutral 72 nm particle (airway-optimized)
        p_opt = compute_mucin_penetration_index(zeta_mv=-0.5, hydrodynamic_diameter_nm=72.0)
        assert 0.40 <= p_opt <= 0.80

        # Cationic particle (zeta = +20.0 mV) has severe electrostatic entrapment
        p_cationic = compute_mucin_penetration_index(zeta_mv=20.0, hydrodynamic_diameter_nm=72.0)
        assert p_cationic < 0.10
        assert p_opt > p_cationic * 5.0

        # Exception on invalid dimensions
        with pytest.raises(ValueError, match="positive"):
            compute_mucin_penetration_index(0.0, hydrodynamic_diameter_nm=-5.0)
        with pytest.raises(ValueError, match="positive"):
            compute_mucin_penetration_index(0.0, mucin_mesh_pore_nm=0.0)

    def test_escape_rate_and_helper_lipid_phase_preference(self) -> None:
        """Cone-shaped DOPE (H_II phase) drives faster escape than cylindrical DSPC."""
        k_dope = compute_escape_rate(ph=5.5, pka=6.45, k_max=0.08, phi_hii=0.078)
        k_dspc = compute_escape_rate(ph=5.5, pka=6.45, k_max=0.08, phi_hii=0.022)

        assert k_dope > k_dspc * 3.0
        assert k_dope > 0.004  # Significant pore formation velocity

        # At pH 7.4 escape rate is very low (< 0.001 min^-1)
        k_extracellular = compute_escape_rate(ph=7.4, pka=6.45, k_max=0.08, phi_hii=0.078)
        assert k_extracellular < 0.001

        with pytest.raises(ValueError, match="non-negative"):
            compute_escape_rate(5.5, 6.45, k_max=-0.01)

    def test_nebulizer_shear_stress_and_retention(self) -> None:
        """Vibrating mesh nebulizer shear dynamics and critical tension protection."""
        shear_rate, shear_stress = compute_mesh_shear_stress(velocity_m_s=15.0, pore_diameter_um=3.0)

        # gamma_dot = 15 / (3.0e-6) = 5.0e6 s^-1
        assert pytest.approx(shear_rate, rel=1e-4) == 5.0e6
        # tau = 0.0012 * 5.0e6 = 6000 Pa
        assert pytest.approx(shear_stress, rel=1e-4) == 6000.0

        # Retention with nominal membrane tension
        ee_standard = compute_nebulization_retention(shear_stress_pa=6000.0, sigma_membrane_mn_m=28.5)
        # Retention with reinforced membrane (polysarcosine / high cholesterol)
        ee_reinforced = compute_nebulization_retention(shear_stress_pa=6000.0, sigma_membrane_mn_m=35.0)

        assert 0.85 <= ee_standard <= 0.95
        assert ee_reinforced > ee_standard

        with pytest.raises(ValueError, match="positive"):
            compute_mesh_shear_stress(velocity_m_s=15.0, pore_diameter_um=0.0)

    def test_endosomal_trafficking_conservation_and_trajectories(self) -> None:
        """Internalized LNP mass is conserved across escape, lysosomal degradation, and recycling."""
        dose = 100.0
        res = simulate_endosomal_trafficking(
            total_time_min=360.0,
            lnp_dose_internalized=dose,
            apparent_pka=6.48,
            phi_hii=0.078,
        )

        total_tracked_pct = (
            res["escape_efficiency_pct"]
            + res["lysosomal_degradation_pct"]
            + res["recycling_pct"]
        )
        # Over 6 hours, all early and late endosomes have resolved (> 99%)
        assert pytest.approx(total_tracked_pct, abs=1.5) == 100.0

        # Escape percentage is within biophysical bounds (5 - 12% for optimized DOPE)
        assert 4.0 <= res["escape_efficiency_pct"] <= 12.0
        assert res["peak_cytosolic_mrna_copies"] > 10.0
        assert 10.0 <= res["peak_mrna_time_min"] <= 120.0

    def test_beacon_cf_vs_airway_optimized_sort(self, model: MRNALNPEscapeModel) -> None:
        """Airway-optimized SORT-LNP outperforms BEACON-CF benchmark by >2.5x in net mRNA yield."""
        beacon_res = model.evaluate_formulation("BEACON_CF_benchmark")
        airway_res = model.evaluate_formulation("DOPE_airway_optimized")
        sort_res = model.evaluate_formulation("SORT_airway_polysarcosine")
        mucus_trapped = model.evaluate_formulation("high_pka_mucus_trapped")

        # BEACON-CF baseline escape is in the classic 1 - 3% bottleneck regime
        assert 1.0 <= beacon_res["endosomal_escape_pct"] <= 3.5

        # DOPE and SORT airway formulations exhibit enhanced escape
        assert airway_res["endosomal_escape_pct"] >= 5.0
        assert sort_res["endosomal_escape_pct"] >= 5.0

        # Net yield comparison (combining nebulizer EE * mucus penetration * escape * cargo)
        assert airway_res["net_mrna_yield_per_lnp"] > beacon_res["net_mrna_yield_per_lnp"] * 2.5
        assert sort_res["net_mrna_yield_per_lnp"] > beacon_res["net_mrna_yield_per_lnp"] * 2.5

        # High pKa formulation fails mucus penetration and ranks as ineffective
        assert mucus_trapped["mucin_penetration_index"] < 0.10
        assert mucus_trapped["delivery_tier"] in (
            "INEFFECTIVE_MUCOSAL_OR_ENDOSOMAL_FAILURE",
            "LOW_POTENCY_SUBOPTIMAL_BOTTLENECK",
        )

    def test_compare_all_and_nebulizer_screen(self, model: MRNALNPEscapeModel) -> None:
        """Catalog comparison and nebulizer modality screens execute deterministically."""
        ranked = model.compare_all_formulations()
        assert len(ranked) == len(LNP_CATALOG)
        # Top formulation should be airway optimized or SORT
        assert ranked[0]["formulation_id"] in ("SORT_airway_polysarcosine", "DOPE_airway_optimized")

        screen = model.run_nebulizer_screen("DOPE_airway_optimized")
        assert len(screen) == 4
        # Respimat soft mist should have lowest shear stress
        respimat = next(m for m in screen if "Respimat" in m["modality"])
        mesh = next(m for m in screen if "Aeroneb" in m["modality"])
        assert respimat["shear_stress_pa"] < mesh["shear_stress_pa"]

    def test_provenance_and_epistemic_guardrails(self) -> None:
        """Ensure mandatory NLM attribution and research boundary disclaimers are present."""
        assert "National Library of Medicine" in NLM_ATTRIBUTION
        assert "research, not medical advice" in RESEARCH_BOUNDARY
        for spec in LNP_CATALOG.values():
            assert spec.trial_surveillance.startswith("NCT")
            assert len(spec.citation_doi) > 0
