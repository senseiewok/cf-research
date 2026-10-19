#!/usr/bin/env python3
"""Unit tests for Extracellular Vesicle (EV) Basal Stem Cell Tropism & Trans-Tight-Junction Permeation Model (T-0234).

Verifies paracellular tight junction permeability dynamics, saturating basal stem cell
endocytic/fusogenic uptake, pseudotyping kinetic acceleration, and durable chimerism attainment.
"""
from __future__ import annotations

import sys
from pathlib import Path
import pytest

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO_ROOT))

from ev_basal_tropism import (
    EVBasalTropismModel,
    compute_tight_junction_permeability,
    compute_basal_uptake_velocity,
    PSEUDOTYPE_CATALOG,
    P_BASELINE_CM_PER_S,
    CURATIVE_CHIMERISM_THRESHOLD,
)


class TestEVBasalTropismModel:
    @pytest.fixture
    def model(self) -> EVBasalTropismModel:
        return EVBasalTropismModel()

    def test_tight_junction_permeability_dynamics(self) -> None:
        """P_TJ(t) begins at P_baseline * (1 + beta) and decays exponentially to P_baseline."""
        p_init = compute_tight_junction_permeability(t_hours=0.0, beta_loosening=10.0, tau_recovery_hours=2.0)
        assert pytest.approx(p_init, rel=1e-5) == P_BASELINE_CM_PER_S * 11.0

        p_recovered = compute_tight_junction_permeability(t_hours=20.0, beta_loosening=10.0, tau_recovery_hours=2.0)
        assert pytest.approx(p_recovered, rel=1e-3) == P_BASELINE_CM_PER_S

        # When beta = 0 (intact barrier), permeability remains strictly constant at baseline
        p_intact_0 = compute_tight_junction_permeability(t_hours=0.0, beta_loosening=0.0)
        p_intact_10 = compute_tight_junction_permeability(t_hours=10.0, beta_loosening=0.0)
        assert p_intact_0 == P_BASELINE_CM_PER_S
        assert p_intact_10 == P_BASELINE_CM_PER_S

        with pytest.raises(ValueError, match="non-negative"):
            compute_tight_junction_permeability(-1.0)
        with pytest.raises(ValueError, match="positive"):
            compute_tight_junction_permeability(0.0, tau_recovery_hours=0.0)

    def test_basal_uptake_velocity_and_saturation(self) -> None:
        """v_entry scales with basolateral concentration and pseudotype enhancement factor."""
        v_zero = compute_basal_uptake_velocity(c_basolateral_nm=0.0)
        assert v_zero == 0.0

        v_naive_10nm = compute_basal_uptake_velocity(c_basolateral_nm=10.0, phi_pseudotype=1.0)
        v_pseudotype_10nm = compute_basal_uptake_velocity(c_basolateral_nm=10.0, phi_pseudotype=4.6)
        assert v_pseudotype_10nm == pytest.approx(v_naive_10nm * 4.6, rel=1e-5)

        # Saturation check: v(1000 nM) should approach Vmax * phi
        v_high = compute_basal_uptake_velocity(c_basolateral_nm=2000.0, phi_pseudotype=1.0)
        assert 0.40 < v_high <= 0.45

    def test_paracellular_loosening_enhances_basal_delivery(self, model: EVBasalTropismModel) -> None:
        """Transient tight junction loosening accelerates peak basolateral concentration and early chimerism."""
        res_intact = model.simulate_delivery(
            "unmodified_exosome",
            initial_apical_dose_nm=10.0,
            beta_loosening=0.0,
            duration_hours=24.0,
        )
        res_loosened = model.simulate_delivery(
            "unmodified_exosome",
            initial_apical_dose_nm=10.0,
            beta_loosening=18.0,
            duration_hours=24.0,
        )

        assert res_loosened["peak_basolateral_conc_nm"] > res_intact["peak_basolateral_conc_nm"]
        # Early chimerism at 1 hr should be substantially higher under loosening
        chim_1hr_intact = res_intact["time_snapshots"]["chimerism_hr_1"]
        chim_1hr_loosened = res_loosened["time_snapshots"]["chimerism_hr_1"]
        assert chim_1hr_loosened > chim_1hr_intact * 10.0

    def test_pseudotyping_accelerates_stem_cell_chimerism(self, model: EVBasalTropismModel) -> None:
        """Engineered integrin-targeted EVs achieve significantly faster curative chimerism than naive EVs."""
        res_naive = model.simulate_delivery(
            "unmodified_exosome",
            initial_apical_dose_nm=25.0,
            beta_loosening=10.0,
        )
        res_targeted = model.simulate_delivery(
            "integrin_a6b4_targeted",
            initial_apical_dose_nm=25.0,
            beta_loosening=10.0,
        )

        # At hour 1, targeted EV chimerism exceeds naive EV chimerism
        assert (
            res_targeted["time_snapshots"]["chimerism_hr_1"]
            > res_naive["time_snapshots"]["chimerism_hr_1"]
        )
        assert res_targeted["meets_curative_chimerism_target"] is True
        assert res_targeted["final_basal_chimerism_pct"] >= CURATIVE_CHIMERISM_THRESHOLD * 100.0

    def test_compare_all_profiles_catalog(self, model: EVBasalTropismModel) -> None:
        """compare_all_profiles tests all 4 cataloged EV formulations."""
        comp = model.compare_all_profiles(initial_apical_dose_nm=50.0, beta_loosening=12.0)
        assert len(comp["profiles"]) == len(PSEUDOTYPE_CATALOG)
        for pid in PSEUDOTYPE_CATALOG:
            assert pid in comp["profiles"]
            assert comp["profiles"][pid]["meets_curative_chimerism_target"] is True

    def test_evidence_grounding_and_boundary_tags(self, model: EVBasalTropismModel) -> None:
        """Every simulation result retains NLM attribution, surveillance NCT, and clinical boundary."""
        res = model.simulate_delivery("integrin_a6b4_targeted")
        assert "NLM" in res["nlm_attribution"]
        assert res["clinical_surveillance_nct"] == "NCT05248230"
        assert "research, not medical advice" in res["research_boundary"]
        assert res["citation_doi"] == "10.3389/fbioe.2026.1855907"
