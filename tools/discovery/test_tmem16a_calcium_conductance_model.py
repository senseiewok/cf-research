#!/usr/bin/env python3
"""Unit test suite for TMEM16A Alternative Chloride Conductance & Calcium Coupling Model (T-0392).

Verifies:
  - Voltage-dependent gating and calcium sensitivity shift (z_eff = 0.40, Hill n_H = 2.0)
  - Allosteric potentiator left-shift into resting basal range (100-200 nM)
  - Resting cytosolic [Ca2+]_i (120 nM) activation without elevating global calcium or triggering bronchoconstriction
  - Apical chloride current I_Cl >= 20% of WT CFTR rescuing ASL height >= 7.0 um
  - Genotype-agnostic rescue across Class I nonsense, Class II misfolding, and Class III gating cohorts
  - Strict grounding in Cryo-EM structures (PDB 6QP6, 5OYB), ClinicalTrials.gov (NCT04488705, NCT06478706),
    and primary literature DOIs (10.3389/fphar.2019.00003, 10.1164/rccm.202004-0887le).
"""

from __future__ import annotations

import sys
from pathlib import Path
import pytest

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tmem16a_calcium_conductance_model import (
    TMEM16ACalciumConductanceModel,
    COHORT_REGISTRY,
    compute_voltage_shifted_k05_nm,
    compute_potentiated_k05_nm,
    compute_tmem16a_open_probability,
    compute_apical_chloride_flux,
    integrate_asl_hydration_ode,
    format_markdown_table,
    DEFAULT_BASAL_CA_NM,
    DEFAULT_K05_ZERO_BASAL_NM,
    DEFAULT_HILL_COEFFICIENT_NH,
    DEFAULT_EFFECTIVE_VALENCE_ZEFF,
    DEFAULT_POTENTIATOR_KD_ACTIVE_NM,
    DEFAULT_POTENTIATOR_KD_INACTIVE_NM,
    DEFAULT_CURATIVE_ASL_THRESHOLD_UM,
    DEFAULT_CURATIVE_ICL_PCT_WT,
    DEFAULT_BRONCHOCONSTRICTION_CA_NM,
    NLM_ATTRIBUTION,
    RESEARCH_BOUNDARY,
)


class TestTMEM16ACalciumConductanceModel:
    @pytest.fixture
    def model(self) -> TMEM16ACalciumConductanceModel:
        return TMEM16ACalciumConductanceModel()

    def test_voltage_dependent_calcium_sensitivity(self) -> None:
        """Membrane potential shifts K_0.5 exponentially with effective gating charge z_eff = 0.40."""
        # At 0 mV, K_0.5 equals the baseline reference value
        k05_zero = compute_voltage_shifted_k05_nm(vm_mv=0.0)
        assert pytest.approx(k05_zero, rel=1e-4) == DEFAULT_K05_ZERO_BASAL_NM

        # Hyperpolarization (-60 mV vs -40 mV vs -20 mV) increases K_0.5 (decreases apparent calcium affinity)
        k05_hyper60 = compute_voltage_shifted_k05_nm(vm_mv=-60.0)
        k05_rest40 = compute_voltage_shifted_k05_nm(vm_mv=-40.0)
        k05_depol20 = compute_voltage_shifted_k05_nm(vm_mv=-20.0)

        assert k05_hyper60 > k05_rest40 > k05_depol20 > k05_zero

        # Check invalid inputs
        with pytest.raises(ValueError, match="k05_zero_nm must be positive"):
            compute_voltage_shifted_k05_nm(k05_zero_nm=0.0)
        with pytest.raises(ValueError, match="z_eff must be non-negative"):
            compute_voltage_shifted_k05_nm(z_eff=-0.1)

    def test_allosteric_potentiator_left_shift(self) -> None:
        """Potentiator shifts apparent calcium K_0.5 from basal >1500 nM down into 100-200 nM range."""
        # Unpotentiated baseline at resting V_m = -40 mV
        k05_unpot = compute_potentiated_k05_nm(potentiator_conc_nm=0.0, vm_mv=-40.0)
        assert k05_unpot > 1500.0  # ~1638 nM

        # At intermediate dose (50 nM), K_0.5 shifts lower
        k05_50nm = compute_potentiated_k05_nm(potentiator_conc_nm=50.0, vm_mv=-40.0)
        assert k05_50nm < k05_unpot

        # At therapeutic dose (200-400 nM), K_0.5 reaches resting physiological range (100-200 nM)
        k05_200nm = compute_potentiated_k05_nm(potentiator_conc_nm=200.0, vm_mv=-40.0)
        k05_400nm = compute_potentiated_k05_nm(potentiator_conc_nm=400.0, vm_mv=-40.0)
        assert 100.0 <= k05_400nm <= 200.0
        assert 100.0 <= k05_200nm <= 250.0
        assert k05_unpot / k05_400nm > 8.0  # Robust left-shift

        # Invariant: monotonic decrease with increasing potentiator
        assert k05_unpot > k05_50nm > k05_200nm > k05_400nm

        with pytest.raises(ValueError, match="non-negative"):
            compute_potentiated_k05_nm(potentiator_conc_nm=-10.0)

    def test_dual_calcium_cooperativity_gating(self) -> None:
        """Open probability follows Hill equation with n_H = 2.0."""
        # At resting Ca2+ (120 nM) without potentiator, channel is closed (P_o < 0.01)
        po_basal = compute_tmem16a_open_probability(ca_conc_nm=DEFAULT_BASAL_CA_NM, potentiator_conc_nm=0.0)
        assert po_basal < 0.01

        # With 200 nM potentiator at basal Ca2+ (120 nM), P_o increases substantially (> 0.15)
        po_pot = compute_tmem16a_open_probability(ca_conc_nm=DEFAULT_BASAL_CA_NM, potentiator_conc_nm=200.0)
        assert po_pot > 0.15
        assert po_pot / po_basal > 20.0

        # High calcium (10 uM) saturates channel (P_o > 0.95)
        po_high_ca = compute_tmem16a_open_probability(ca_conc_nm=10000.0, potentiator_conc_nm=200.0)
        assert po_high_ca > 0.95

        # Error checks
        with pytest.raises(ValueError, match="ca_conc_nm must be positive"):
            compute_tmem16a_open_probability(ca_conc_nm=0.0)
        with pytest.raises(ValueError, match="n_hill must be positive"):
            compute_tmem16a_open_probability(n_hill=0.0)

    def test_apical_chloride_flux_and_cftr_equivalence(self) -> None:
        """Apical current density scales with open probability and meets >= 20% WT CFTR threshold."""
        # At P_o = 0, current is 0
        flux_zero = compute_apical_chloride_flux(open_probability=0.0)
        assert flux_zero["tmem16a_current_ua_cm2"] == 0.0
        assert flux_zero["total_apical_current_ua_cm2"] == 0.0

        # At P_o = 0.08, current is >= 9.0 uA/cm^2 (20% of 45.0 uA/cm^2)
        flux_pot = compute_apical_chloride_flux(open_probability=0.08)
        assert flux_pot["tmem16a_current_ua_cm2"] >= 9.0
        assert flux_pot["pct_wt_cftr_current"] >= DEFAULT_CURATIVE_ICL_PCT_WT

        with pytest.raises(ValueError, match="open_probability must be in"):
            compute_apical_chloride_flux(open_probability=1.2)

    def test_asl_hydration_ode_recovery(self) -> None:
        """ASL height recovers from 3.2 um CF baseline to >= 7.0 um when chloride secretion >= 9.0 uA/cm^2."""
        # Dehydrated CF baseline without current (0.2 uA/cm^2)
        asl_unpot = integrate_asl_hydration_ode(total_icl_ua_cm2=0.2)
        assert asl_unpot["steady_state_asl_um"] < 4.0

        # Potentiated current (11.0 uA/cm^2 ~ 24% WT CFTR) restores ASL >= 7.0 um
        asl_pot = integrate_asl_hydration_ode(total_icl_ua_cm2=11.0)
        assert asl_pot["steady_state_asl_um"] >= DEFAULT_CURATIVE_ASL_THRESHOLD_UM
        assert 7.0 <= asl_pot["steady_state_asl_um"] <= 8.5

        # Monotonic time series trajectory
        trajectory = asl_pot["trajectory"]
        assert len(trajectory) > 5
        assert trajectory[0]["asl_height_um"] == 3.2
        assert trajectory[-1]["asl_height_um"] >= 7.0
        for i in range(len(trajectory) - 1):
            assert trajectory[i + 1]["asl_height_um"] >= trajectory[i]["asl_height_um"]

        with pytest.raises(ValueError, match="total_icl_ua_cm2 must be non-negative"):
            integrate_asl_hydration_ode(total_icl_ua_cm2=-1.0)

    def test_genotype_agnostic_rescue_across_cohorts(self, model: TMEM16ACalciumConductanceModel) -> None:
        """TMEM16A potentiation rescues ASL height >= 7.0 um across Class I, II, and III cohorts."""
        # Untreated Class I null: collapsed ASL, sub-curative
        res_c1_untreated = model.simulate("Class_I_Nonsense_Null", potentiator_conc_nm=0.0)
        assert res_c1_untreated.steady_state_asl_um < 4.0
        assert not res_c1_untreated.ciliary_clearance_restored
        assert not res_c1_untreated.curative_conductance_achieved
        assert not res_c1_untreated.bronchoconstriction_risk

        # Potentiated (200 nM) at basal Ca2+ (120 nM): rescues ASL >= 7.0 um without raising Ca2+
        res_c1_pot = model.simulate("Class_I_Nonsense_Null", potentiator_conc_nm=200.0)
        assert res_c1_pot.steady_state_asl_um >= DEFAULT_CURATIVE_ASL_THRESHOLD_UM
        assert res_c1_pot.ciliary_clearance_restored
        assert res_c1_pot.curative_conductance_achieved
        assert not res_c1_pot.bronchoconstriction_risk

        # Rescues Class II and Class III cohorts as well
        res_c2_pot = model.simulate("Class_II_Refractory_Folding", potentiator_conc_nm=200.0)
        res_c3_pot = model.simulate("Class_III_Gating_Defect", potentiator_conc_nm=200.0)
        assert res_c2_pot.ciliary_clearance_restored
        assert res_c3_pot.ciliary_clearance_restored
        assert res_c3_pot.steady_state_asl_um >= res_c2_pot.steady_state_asl_um >= res_c1_pot.steady_state_asl_um

    def test_bronchoconstriction_safety_boundary(self, model: TMEM16ACalciumConductanceModel) -> None:
        """Global cytosolic calcium elevation >= 350 nM triggers bronchoconstriction risk flag."""
        # Basal Ca2+ (120 nM) with potentiator: safe
        safe_res = model.simulate("Class_I_Nonsense_Null", potentiator_conc_nm=200.0, cytosolic_ca_nm=120.0)
        assert not safe_res.bronchoconstriction_risk

        # Ca2+ spiked to 500 nM (e.g. non-selective calcium ionophore or P2Y2 hyperstimulation): flag raised
        unsafe_res = model.simulate("Class_I_Nonsense_Null", potentiator_conc_nm=0.0, cytosolic_ca_nm=500.0)
        assert unsafe_res.bronchoconstriction_risk

    def test_structural_and_literature_grounding(self, model: TMEM16ACalciumConductanceModel) -> None:
        """Verifies PDB coordinates, NCT IDs, DOIs, and disclaimers are present in results and reports."""
        res = model.simulate("Class_I_Nonsense_Null", potentiator_conc_nm=200.0)
        assert res.cryoem_pdb_structures == ("6QP6", "5OYB")
        assert res.registry_trials == ("NCT04488705", "NCT06478706")
        assert res.primary_dois == ("10.3389/fphar.2019.00003", "10.1164/rccm.202004-0887le")

        sweep = model.run_dose_response_sweep()
        assert len(sweep) == len(COHORT_REGISTRY) * 6

        report = format_markdown_table(sweep)
        assert "NCT04488705" in report
        assert "NCT06478706" in report
        assert "6QP6" in report
        assert "5OYB" in report
        assert "10.3389/fphar.2019.00003" in report
        assert "10.1164/rccm.202004-0887le" in report
        assert NLM_ATTRIBUTION in report
        assert RESEARCH_BOUNDARY in report
