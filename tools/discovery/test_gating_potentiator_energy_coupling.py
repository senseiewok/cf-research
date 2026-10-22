#!/usr/bin/env python3
"""Unit test suite for CFTR Gating Potentiator Allosteric Cooperativity & Energy Coupling Model (T-0270).

Acceptance Criteria:
  pytest tools/discovery/test_gating_potentiator_energy_coupling.py -v; python tools/claims/check_anti_hallucination.py
"""
import math
import pytest
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from gating_potentiator_energy_coupling import (
    R_GAS_KCAL,
    TEMP_PHYSIOLOGICAL_K,
    ALLELE_REGISTRY,
    POTENTIATOR_REGISTRY,
    AlleleProfile,
    PotentiatorProfile,
    calculate_coupling_energy,
    calculate_macroscopic_current_pa,
    calculate_po_dual,
    calculate_po_single,
    compute_mutant_cycle_analysis,
    evaluate_discovery_panel,
    main,
    simulate_dose_response,
)


class TestPotentiatorRegistry:
    """Validate thermodynamic registry integrity and biophysical invariants."""

    def test_potentiator_entries_present(self):
        assert "ivacaftor" in POTENTIATOR_REGISTRY
        assert "deutivacaftor" in POTENTIATOR_REGISTRY
        assert "glpg1837" in POTENTIATOR_REGISTRY

    def test_potentiator_thermodynamic_invariants(self):
        for name, p in POTENTIATOR_REGISTRY.items():
            assert p.kd_open_nm > 0
            assert p.kd_closed_nm > p.kd_open_nm, f"{name}: Kd closed must exceed Kd open"
            assert p.allosteric_ratio_alpha > 1.0
            assert p.coupling_free_energy_kcal < 0.0, f"{name}: DeltaDelta G must be negative to stabilize open pore"
            assert len(p.pdb_id) == 4
            assert p.evidence_doi.startswith("10.")


class TestAlleleRegistry:
    """Validate CFTR gating allele baseline parameters."""

    def test_allele_entries_present(self):
        assert "WT" in ALLELE_REGISTRY
        assert "G551D" in ALLELE_REGISTRY
        assert "G1244E" in ALLELE_REGISTRY
        assert "S549N" in ALLELE_REGISTRY
        assert "F508del_corrected" in ALLELE_REGISTRY

    def test_allele_basal_gating_invariants(self):
        for name, a in ALLELE_REGISTRY.items():
            assert 0.0 < a.basal_po < 1.0
            assert a.single_channel_conductance_ps > 0.0
            assert a.basal_k0 > 0.0
            # Severe gating mutants have severely lower basal Po than WT
            if name in {"G551D", "G1244E", "S549N"}:
                assert a.basal_po < ALLELE_REGISTRY["WT"].basal_po / 10.0


class TestCouplingEnergyCalculation:
    """Test state-dependent free energy coupling calculation."""

    def test_coupling_energy_exact(self):
        # When Kd_closed / Kd_open = e, DeltaDelta G = -R * T * 1 = -1.9872e-3 * 310.15 = -0.61633 kcal/mol
        kd_open = 10.0
        kd_closed = 10.0 * math.e
        ddg = calculate_coupling_energy(kd_open, kd_closed, temp_k=310.15)
        assert pytest.approx(ddg, rel=1e-4) == -1.98720425864083e-3 * 310.15

    def test_coupling_energy_invalid_inputs(self):
        with pytest.raises(ValueError):
            calculate_coupling_energy(-10.0, 100.0)
        with pytest.raises(ValueError):
            calculate_coupling_energy(10.0, 0.0)


class TestSinglePotentiatorPo:
    """Test channel open probability calculation under single potentiator."""

    def test_po_zero_concentration_equals_basal(self):
        basal = 0.0025
        po = calculate_po_single(basal_po=basal, conc_nm=0.0, kd_open_nm=24.0, kd_closed_nm=380.0)
        assert pytest.approx(po, rel=1e-6) == basal

    def test_po_monotonically_increases_with_concentration(self):
        basal = 0.0025
        concs = [0.0, 1.0, 10.0, 100.0, 1000.0, 10000.0]
        pos = [calculate_po_single(basal, c, 24.0, 380.0) for c in concs]
        for i in range(len(pos) - 1):
            assert pos[i] < pos[i + 1]

    def test_po_saturating_limit(self):
        basal = 0.0025
        kd_o = 24.0
        kd_c = 380.0
        # Saturating Po formula: 1 / (1 + ((1 - basal) / basal) * (kd_o / kd_c))
        expected_sat = 1.0 / (1.0 + ((1.0 - basal) / basal) * (kd_o / kd_c))
        po_high = calculate_po_single(basal, 1e7, kd_o, kd_c)
        assert pytest.approx(po_high, rel=1e-3) == expected_sat

    def test_po_invalid_inputs(self):
        with pytest.raises(ValueError):
            calculate_po_single(0.0, 10.0, 20.0, 200.0)
        with pytest.raises(ValueError):
            calculate_po_single(1.0, 10.0, 20.0, 200.0)
        with pytest.raises(ValueError):
            calculate_po_single(0.1, -5.0, 20.0, 200.0)


class TestDualPotentiatorPo:
    """Test dual-potentiator cooperativity model."""

    def test_dual_zero_concentration_equals_basal(self):
        basal = 0.012
        po = calculate_po_dual(
            basal_po=basal,
            conc1_nm=0.0, kd1_open_nm=24.0, kd1_closed_nm=380.0,
            conc2_nm=0.0, kd2_open_nm=9.5, kd2_closed_nm=260.0,
            alpha_coop=1.0,
        )
        assert pytest.approx(po, rel=1e-6) == basal

    def test_dual_cooperativity_effects(self):
        basal = 0.0025
        po_additive = calculate_po_dual(
            basal, 100.0, 24.0, 380.0, 50.0, 9.5, 260.0, alpha_coop=1.0
        )
        po_synergistic = calculate_po_dual(
            basal, 100.0, 24.0, 380.0, 50.0, 9.5, 260.0, alpha_coop=0.7
        )
        po_antagonistic = calculate_po_dual(
            basal, 100.0, 24.0, 380.0, 50.0, 9.5, 260.0, alpha_coop=1.4
        )
        assert po_synergistic > po_additive
        assert po_antagonistic < po_additive


class TestMacroscopicCurrent:
    """Test macroscopic chloride current calculation."""

    def test_current_scaling(self):
        # 1000 channels, 8.4 pS, Vm = -60 mV, ECl = -30 mV (driving force -30 mV)
        # single channel current = 8.4 * (-30) * 1e-3 = -0.252 pA
        # at Po = 0.5: 1000 * -0.252 * 0.5 = -126 pA
        curr = calculate_macroscopic_current_pa(po=0.5, n_channels=1000, gamma_conductance_ps=8.4)
        assert pytest.approx(curr, rel=1e-3) == -126.0


class TestMutantCycleAnalysis:
    """Test double-mutant cycle energetic coupling."""

    def test_mutant_cycle_coupling_detection(self):
        # Case with strong allosteric coupling
        res = compute_mutant_cycle_analysis(
            po_wt=0.38,
            po_mut_a=0.12,
            po_mut_b=0.0025,
            po_double=0.0002,
        )
        assert res.is_energetically_coupled is True
        assert abs(res.delta_delta_g_int_kcal) >= 0.5
        assert "coupling" in res.coupling_mechanism.lower()

    def test_mutant_cycle_additive_non_coupled(self):
        # Simulate exact additive case where ddg_int = 0
        wt_po = 0.38
        dg_wt = -R_GAS_KCAL * TEMP_PHYSIOLOGICAL_K * math.log(wt_po / (1.0 - wt_po))
        ddg_a = 1.0  # +1.0 kcal
        ddg_b = 2.0  # +2.0 kcal
        ddg_ab = 3.0 # +3.0 kcal (perfectly additive)

        def _dg_to_po(dg):
            k = math.exp(-dg / (R_GAS_KCAL * TEMP_PHYSIOLOGICAL_K))
            return k / (1.0 + k)

        po_a = _dg_to_po(dg_wt + ddg_a)
        po_b = _dg_to_po(dg_wt + ddg_b)
        po_ab = _dg_to_po(dg_wt + ddg_ab)

        res = compute_mutant_cycle_analysis(wt_po, po_a, po_b, po_ab)
        assert pytest.approx(res.delta_delta_g_int_kcal, abs=1e-3) == 0.0
        assert res.is_energetically_coupled is False
        assert "independent" in res.coupling_mechanism.lower()


class TestSimulationAndPanel:
    """Test dose-response simulation and multi-allele discovery panel."""

    def test_simulate_dose_response_g551d_ivacaftor(self):
        res = simulate_dose_response("G551D", "ivacaftor")
        assert res["allele"]["id"] == "G551D"
        assert res["potentiator"]["name"].startswith("Ivacaftor")
        assert res["coupling_energy_kcal"] < -1.5
        assert 12.0 < res["max_fold_increase"] < 20.0  # Fold rescue bounded by allosteric alpha ratio (~15.8x)
        assert len(res["dose_response_curve"]) == 13

    def test_simulate_dose_response_invalid_keys(self):
        with pytest.raises(KeyError):
            simulate_dose_response("UNKNOWN_ALLELE", "ivacaftor")
        with pytest.raises(KeyError):
            simulate_dose_response("G551D", "UNKNOWN_POTENTIATOR")

    def test_evaluate_discovery_panel(self):
        panel = evaluate_discovery_panel()
        assert "potentiators" in panel
        assert "alleles" in panel
        assert "single_potentiator_responses" in panel
        assert "dual_potentiator_synergy" in panel
        assert "mutant_cycle_coupling" in panel
        assert "G551D" in panel["dual_potentiator_synergy"]
        dual_g551d = panel["dual_potentiator_synergy"]["G551D"]
        assert dual_g551d["dual_synergistic_po"] > dual_g551d["dual_additive_po"]


class TestCliInterface:
    """Test CLI execution and formatting."""

    def test_cli_default(self, capsys):
        code = main(["--allele", "G551D", "--potentiator", "ivacaftor"])
        assert code == 0
        captured = capsys.readouterr()
        assert "CFTR Gating Coupling Simulation" in captured.out
        assert "G551D" in captured.out

    def test_cli_panel_json(self, capsys):
        code = main(["--panel", "--json"])
        assert code == 0
        captured = capsys.readouterr()
        assert '"dual_potentiator_synergy"' in captured.out
