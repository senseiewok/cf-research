#!/usr/bin/env python3
"""Unit tests for TMEM16A/ANO1 Allosteric Potentiation & CaCC Conductance Cooperativity Model (T-0240).

Verifies:
  - Allosteric EC50 left-shift kinetics (EC50_basal = 1850 nM, K_D = 42.0 nM, alpha = 0.075 -> 13.3x sensitization)
  - Cooperative calcium and voltage-dependent gating (Hill n_H = 2.1, PDB: 5OYB and 6QP6)
  - Resting cytosolic calcium ([Ca2+]_i = 95 nM) activation without crossing the 350 nM bronchoconstriction threshold
  - Transepithelial chloride current I_Cl,TMEM16A >= 20% WT CFTR restoring ASL height >= 7.0 um (nominal ~7.2 um)
  - Genotype-agnostic rescue across Class I pan-null, Class II refractory, and Class III cohorts
  - Grounding in NCT06478706, PDB 5OYB/6QP6, and DOIs 10.1038/nature24652, 10.1038/s41598-020-73905-2, 10.1016/j.jcf.2019.06.004
"""

from __future__ import annotations

import sys
from pathlib import Path
import pytest

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO_ROOT))

from tmem16a_conductance_model import (
    TMEM16AConductanceModel,
    GENOTYPE_COHORT_CATALOG,
    compute_apparent_ec50_nm,
    compute_open_probability,
    compute_apical_chloride_current,
    simulate_asl_hydration_kinetics,
    format_markdown_report,
    DEFAULT_BASAL_CA_NM,
    DEFAULT_EC50_BASAL_NM,
    DEFAULT_POTENTIATOR_KD_NM,
    DEFAULT_ALLOSTERIC_ALPHA,
    DEFAULT_HILL_COEFF_NH,
    DEFAULT_UNITARY_GAMMA_PS,
    DEFAULT_CURATIVE_ASL_THRESHOLD_UM,
    DEFAULT_CURATIVE_ICL_PCT_WT,
    NLM_ATTRIBUTION,
    RESEARCH_BOUNDARY,
)


class TestTMEM16AConductanceModel:
    @pytest.fixture
    def model(self) -> TMEM16AConductanceModel:
        return TMEM16AConductanceModel()

    def test_allosteric_ec50_left_shift(self) -> None:
        """Allosteric modulator shifts EC50 from 1850 nM down toward alpha * EC50_basal (~138.75 nM)."""
        # Unpotentiated baseline EC50 = 1850 nM (1.85 uM)
        ec50_0 = compute_apparent_ec50_nm(0.0)
        assert pytest.approx(ec50_0, 1e-4) == DEFAULT_EC50_BASAL_NM

        # At K_D (42.0 nM), EC50 is midway between EC50_basal and alpha * EC50_basal
        ec50_kd = compute_apparent_ec50_nm(DEFAULT_POTENTIATOR_KD_NM)
        expected_mid = DEFAULT_EC50_BASAL_NM * (1.0 + DEFAULT_ALLOSTERIC_ALPHA) / 2.0
        assert pytest.approx(ec50_kd, 1e-4) == expected_mid

        # At saturating ligand concentration, fold-sensitization approaches 1 / alpha = 13.33x (~138.75 nM)
        ec50_sat = compute_apparent_ec50_nm(100000.0)
        assert pytest.approx(ec50_sat, 1e-2) == DEFAULT_EC50_BASAL_NM * DEFAULT_ALLOSTERIC_ALPHA
        assert pytest.approx(DEFAULT_EC50_BASAL_NM / ec50_sat, 1e-2) == 1.0 / DEFAULT_ALLOSTERIC_ALPHA

        # Monotonically decreasing with ligand concentration
        assert compute_apparent_ec50_nm(10.0) > compute_apparent_ec50_nm(50.0) > compute_apparent_ec50_nm(250.0)

        with pytest.raises(ValueError, match="non-negative"):
            compute_apparent_ec50_nm(-5.0)

    def test_cooperative_calcium_and_voltage_gating(self) -> None:
        """Open probability increases cooperatively with [Ca2+]_i, [L], and membrane depolarization."""
        # At resting Ca2+ (95 nM) without potentiator, channel is nearly closed (P_open < 0.5%)
        p_basal = compute_open_probability(ca_conc_nm=DEFAULT_BASAL_CA_NM, ligand_conc_nm=0.0, vm_mv=-40.0)
        assert p_basal < 0.005

        # With 250 nM allosteric potentiator at resting Ca2+ (95 nM), P_open increases >20-fold
        p_pot = compute_open_probability(ca_conc_nm=DEFAULT_BASAL_CA_NM, ligand_conc_nm=250.0, vm_mv=-40.0)
        assert p_pot > 0.08
        assert (p_pot / p_basal) > 20.0

        # Depolarizing apical voltage (-20 mV vs -60 mV) increases P_open
        p_hyper = compute_open_probability(ca_conc_nm=DEFAULT_BASAL_CA_NM, ligand_conc_nm=250.0, vm_mv=-60.0)
        p_depol = compute_open_probability(ca_conc_nm=DEFAULT_BASAL_CA_NM, ligand_conc_nm=250.0, vm_mv=-20.0)
        assert p_depol > p_pot > p_hyper

        # High calcium (10 uM) saturates channel near 1.0
        p_sat = compute_open_probability(ca_conc_nm=10000.0, ligand_conc_nm=250.0, vm_mv=0.0)
        assert p_sat > 0.98

    def test_apical_chloride_current_and_wt_cftr_equivalence(self) -> None:
        """Macroscopic I_Cl scales with P_open and unitary conductance gamma_Cl = 2.8 pS."""
        curr_zero = compute_apical_chloride_current(p_open=0.0)
        assert curr_zero["tmem16a_icl_ua_cm2"] == 0.0
        assert curr_zero["icl_pct_of_wt_cftr"] == 0.0

        # At P_open = 0.0894 (~250 nM ETD001), I_Cl exceeds 20% of WT CFTR (9.0 uA/cm^2)
        curr_pot = compute_apical_chloride_current(p_open=0.0894)
        assert curr_pot["tmem16a_icl_ua_cm2"] >= 9.0
        assert curr_pot["icl_pct_of_wt_cftr"] >= DEFAULT_CURATIVE_ICL_PCT_WT

        with pytest.raises(ValueError, match="p_open"):
            compute_apical_chloride_current(p_open=1.5)

    def test_asl_hydration_kinetics_and_7um_threshold(self) -> None:
        """ASL height recovers from ~3.2 um CF baseline to >= 7.0 um when I_Cl >= 20% WT CFTR."""
        # Untreated Class I CF (0.33 uA/cm^2): remains collapsed around 3.2 - 3.6 um
        asl_untreated = simulate_asl_hydration_kinetics(total_icl_ua_cm2=0.33)
        assert 3.2 <= float(asl_untreated["steady_state_asl_um"]) <= 3.6

        # Threshold 20% WT CFTR (9.0 uA/cm^2) reaches >= 6.2 um, and 25% WT CFTR (11.27 uA/cm^2 at 250 nM) exceeds 7.2 um
        asl_potentiated = simulate_asl_hydration_kinetics(total_icl_ua_cm2=11.27)
        assert float(asl_potentiated["steady_state_asl_um"]) >= DEFAULT_CURATIVE_ASL_THRESHOLD_UM
        assert 7.1 <= float(asl_potentiated["steady_state_asl_um"]) <= 7.5

        # Trajectory monotonically increases toward steady state
        traj = asl_potentiated["trajectory"]
        assert isinstance(traj, list)
        assert len(traj) > 10
        assert traj[0]["asl_height_um"] == 3.2
        assert traj[-1]["asl_height_um"] > 6.8
        for i in range(len(traj) - 1):
            assert traj[i + 1]["asl_height_um"] >= traj[i]["asl_height_um"]

    def test_genotype_agnostic_rescue_and_bronchoconstriction_safety(
        self, model: TMEM16AConductanceModel
    ) -> None:
        """250 nM allosteric potentiator rescues Class I, II, and III cohorts at resting Ca2+ without bronchoconstriction."""
        # 1. Untreated Class I pan-null: collapsed ASL (< 4.0 um), sub-curative
        untreated = model.simulate_condition("Class_I_Pan_Null", potentiator_conc_nm=0.0)
        assert untreated.steady_state_asl_um < 4.0
        assert not untreated.ciliary_clearance_restored
        assert not untreated.curative_conductance_met
        assert not untreated.bronchoconstriction_risk

        # 2. Allosteric potentiation (250 nM ETD001-like) at resting 95 nM Ca2+:
        # Restores ASL >= 7.0 um (nominal ~7.22 um) and I_Cl >= 20% WT without elevating cytosolic Ca2+
        rescued_c1 = model.simulate_condition("Class_I_Pan_Null", potentiator_conc_nm=250.0)
        assert rescued_c1.icl_pct_of_wt_cftr >= DEFAULT_CURATIVE_ICL_PCT_WT
        assert rescued_c1.steady_state_asl_um >= DEFAULT_CURATIVE_ASL_THRESHOLD_UM
        assert rescued_c1.ciliary_clearance_restored
        assert rescued_c1.curative_conductance_met
        assert not rescued_c1.bronchoconstriction_risk

        # 3. Genotype-agnostic rescue across Class II and Class III cohorts
        rescued_c2 = model.simulate_condition("Class_II_Refractory", potentiator_conc_nm=250.0)
        rescued_c3 = model.simulate_condition("Class_III_Gating_Suboptimal", potentiator_conc_nm=250.0)
        assert rescued_c2.ciliary_clearance_restored
        assert rescued_c3.ciliary_clearance_restored
        assert rescued_c3.steady_state_asl_um >= rescued_c2.steady_state_asl_um >= rescued_c1.steady_state_asl_um

        # 4. Contrast with direct Ca2+ elevation agonist (e.g. P2Y2 agonist spiking Ca2+ to 600 nM without potentiator):
        # Triggers bronchoconstriction / goblet cell exocytosis risk flag
        ca_spike = model.simulate_condition("Class_I_Pan_Null", potentiator_conc_nm=0.0, cytosolic_ca_nm=600.0)
        assert ca_spike.bronchoconstriction_risk

    def test_matrix_and_markdown_report_grounding(self, model: TMEM16AConductanceModel) -> None:
        """Dose-response matrix and Markdown report include required Cryo-EM PDBs, NCT ID, DOIs, and disclaimers."""
        matrix = model.generate_dose_response_matrix()
        assert len(matrix) == len(GENOTYPE_COHORT_CATALOG) * 6

        report = format_markdown_report(matrix)
        assert "NCT06478706" in report
        assert "5OYB" in report
        assert "6QP6" in report
        assert "10.1038/nature24652" in report
        assert "10.1038/s41598-020-73905-2" in report
        assert "10.1016/j.jcf.2019.06.004" in report
        assert NLM_ATTRIBUTION in report
        assert RESEARCH_BOUNDARY in report
