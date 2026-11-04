#!/usr/bin/env python3
"""Unit tests for the standalone CFTR Biophysics Bridge (`formal_proof_thermodynamic_landscape.py`).

Pure-mathematics and Lean 4 formal proof derivation tests live in `senseiewok/math`, while
this suite verifies that `cf-research` remains 100% standalone and independent via its pinned
`tools/discovery/data/canonical_math_bridge_snapshot.json` bridge.
"""
from __future__ import annotations

from pathlib import Path
import sys
import unittest

DISCOVERY_DIR = Path(__file__).resolve().parent
if str(DISCOVERY_DIR) not in sys.path:
    sys.path.insert(0, str(DISCOVERY_DIR))

from formal_proof_thermodynamic_landscape import (  # noqa: E402
    CHALLENGE_STAGE_MAP,
    COUNTERFACTUAL_KNOCKOUT_CATALOG,
    PROOF_THERMODYNAMIC_PARAMETERS,
    SNAPSHOT_PATH,
    evaluate_challenge_landscape,
    evaluate_clausius_duhem_continuum_mucus_rheology,
    evaluate_falsification_negative_control,
    evaluate_harada_sasa_fluctuation_dissipation_gating,
    evaluate_hatano_sasa_excess_entropy_gating,
    evaluate_information_thermodynamics_landauer_allosteric_feedback,
    evaluate_jarzynski_crooks_fluctuation_power_stroke,
    evaluate_kramers_escape_and_brinkman_permeability,
    evaluate_onsager_electroosmotic_transport,
    evaluate_onsager_machlup_transition_paths,
    evaluate_prigogine_minimum_entropy_production,
    evaluate_glansdorff_prigogine_evolution,
    evaluate_schnakenberg_cycle_kinetics,
    evaluate_onsager_machlup_action,
    evaluate_freidlin_wentzell_quasipotential,
    evaluate_graham_haken_potential,
    evaluate_maes_netocny_response,
    evaluate_knst_nonequilibrium_work,
    evaluate_geometric_gauge_holonomy,
    evaluate_sasa_tasaki_equality,
    evaluate_seifert_trajectory_entropy_production,
    evaluate_prost_joanny_parrondo_fdt,
    evaluate_speck_seifert_fdt,
    evaluate_kpv_relative_entropy,
    evaluate_dechant_sasa_bounds,
    evaluate_current_fluctuation_tur,
    evaluate_pietzonka_seifert_bounds,
    evaluate_hartich_barato_seifert_bipartite,
    evaluate_parrondo_horowitz_sagawa_information,
    evaluate_sagawa_ueda_fluctuation_efficiency,
    evaluate_barato_hartich_seifert_dissipation,
    evaluate_horowitz_parrondo_feedback_control,
    evaluate_dechant_sasa_speed_limits,
    evaluate_hatano_sasa_steady_state_thermodynamics,
    evaluate_koyuk_seifert_cyclic_speed_dissipation,
    evaluate_falasco_esposito_dissipation_time,
    evaluate_dechant_sasa_fluctuation_response,
    evaluate_lebowitz_spohn_action,
    evaluate_speck_seifert_housekeeping,
    evaluate_polettini_esposito_cycle_decomposition,
    evaluate_rao_esposito_conservation_laws,
    evaluate_manikandan_gupta_krishnamurthy_bounds,
    evaluate_wachtel_rao_esposito_coarse_graining,
    evaluate_stochastic_thermodynamics_ness_ion_transport,
    evaluate_thermodynamic_uncertainty_relations_cftr_precision,
    evaluate_two_way_bridge_isomorphism,
    simulate_counterfactual_axiom_knockout,
)


class TestCFTRBiophysicsMathBridge(unittest.TestCase):
    def test_snapshot_exists_and_is_self_contained(self) -> None:
        self.assertTrue(SNAPSHOT_PATH.is_file())
        self.assertGreaterEqual(len(PROOF_THERMODYNAMIC_PARAMETERS), 15)

    def test_counterfactual_knockout_and_allosteric_bridge(self) -> None:
        intact = simulate_counterfactual_axiom_knockout("wild_type_intact_proof")
        p2 = simulate_counterfactual_axiom_knockout("pythagorean_exponent_p2_mutation")
        sorry = simulate_counterfactual_axiom_knockout("unverified_sorry_axiom_injection")
        self.assertTrue(intact["qed_comparator_verified"])
        self.assertFalse(p2["qed_comparator_verified"])
        self.assertFalse(sorry["qed_comparator_verified"])
        self.assertEqual(len(COUNTERFACTUAL_KNOCKOUT_CATALOG), 5)

    def test_challenge_landscapes_and_sensitivity_invariants(self) -> None:
        self.assertEqual(len(CHALLENGE_STAGE_MAP), 4)
        for c_id in CHALLENGE_STAGE_MAP:
            bridged = evaluate_challenge_landscape(c_id, apply_allosteric_bridge=True)
            unbridged = evaluate_challenge_landscape(c_id, apply_allosteric_bridge=False)
            self.assertLessEqual(
                bridged["rate_limiting_peak_barrier_kcal_mol"],
                unbridged["unbridged_peak_barrier_kcal_mol"],
            )

    def test_cftr_biophysical_bridges_and_ah10_negative_control(self) -> None:
        report = evaluate_two_way_bridge_isomorphism()
        self.assertTrue(report["all_challenges_rank_order_preserved_25pct"])
        self.assertTrue(evaluate_kramers_escape_and_brinkman_permeability()["kramers_brinkman_bridge_verified"])
        self.assertTrue(evaluate_jarzynski_crooks_fluctuation_power_stroke()["jarzynski_crooks_bridge_verified"])
        self.assertTrue(evaluate_onsager_machlup_transition_paths()["onsager_machlup_bridge_verified"])
        self.assertTrue(evaluate_stochastic_thermodynamics_ness_ion_transport()["stochastic_thermodynamics_ness_bridge_verified"])
        self.assertTrue(evaluate_information_thermodynamics_landauer_allosteric_feedback()["information_thermodynamics_landauer_bridge_verified"])
        self.assertTrue(evaluate_thermodynamic_uncertainty_relations_cftr_precision()["tur_precision_dissipation_bridge_verified"])
        self.assertTrue(evaluate_harada_sasa_fluctuation_dissipation_gating()["harada_sasa_bridge_verified"])
        self.assertTrue(evaluate_hatano_sasa_excess_entropy_gating()["hatano_sasa_bridge_verified"])
        self.assertTrue(evaluate_clausius_duhem_continuum_mucus_rheology()["clausius_duhem_bridge_verified"])
        self.assertTrue(evaluate_onsager_electroosmotic_transport()["onsager_transport_bridge_verified"])
        self.assertTrue(evaluate_prigogine_minimum_entropy_production()["prigogine_entropy_bridge_verified"])
        self.assertTrue(evaluate_glansdorff_prigogine_evolution()["glansdorff_prigogine_bridge_verified"])
        self.assertTrue(evaluate_schnakenberg_cycle_kinetics()["schnakenberg_bridge_verified"])
        self.assertTrue(evaluate_onsager_machlup_action()["onsager_machlup_action_bridge_verified"])
        self.assertTrue(evaluate_freidlin_wentzell_quasipotential()["freidlin_wentzell_bridge_verified"])
        self.assertTrue(evaluate_graham_haken_potential()["graham_haken_bridge_verified"])
        self.assertTrue(evaluate_maes_netocny_response()["maes_netocny_bridge_verified"])
        self.assertTrue(evaluate_knst_nonequilibrium_work()["knst_work_bridge_verified"])
        self.assertTrue(evaluate_geometric_gauge_holonomy()["geometric_gauge_bridge_verified"])
        self.assertTrue(evaluate_sasa_tasaki_equality()["sasa_tasaki_bridge_verified"])
        self.assertTrue(evaluate_seifert_trajectory_entropy_production()["seifert_trajectory_bridge_verified"])
        self.assertTrue(evaluate_prost_joanny_parrondo_fdt()["prost_joanny_parrondo_bridge_verified"])
        self.assertTrue(evaluate_speck_seifert_fdt()["speck_seifert_bridge_verified"])
        self.assertTrue(evaluate_kpv_relative_entropy()["kpv_bridge_verified"])
        self.assertTrue(evaluate_dechant_sasa_bounds()["dechant_sasa_bridge_verified"])
        self.assertTrue(evaluate_current_fluctuation_tur()["current_fluctuation_tur_bridge_verified"])
        self.assertTrue(evaluate_pietzonka_seifert_bounds()["pietzonka_seifert_bounds_bridge_verified"])
        self.assertTrue(evaluate_hartich_barato_seifert_bipartite()["hartich_barato_seifert_bipartite_bridge_verified"])
        self.assertTrue(evaluate_parrondo_horowitz_sagawa_information()["parrondo_horowitz_sagawa_information_bridge_verified"])
        self.assertTrue(evaluate_sagawa_ueda_fluctuation_efficiency()["sagawa_ueda_fluctuation_efficiency_bridge_verified"])
        self.assertTrue(evaluate_barato_hartich_seifert_dissipation()["barato_hartich_seifert_dissipation_bridge_verified"])
        self.assertTrue(evaluate_horowitz_parrondo_feedback_control()["horowitz_parrondo_feedback_control_bridge_verified"])
        self.assertTrue(evaluate_dechant_sasa_speed_limits()["dechant_sasa_speed_limits_bridge_verified"])
        self.assertTrue(evaluate_hatano_sasa_steady_state_thermodynamics()["hatano_sasa_steady_state_bridge_verified"])
        self.assertTrue(evaluate_koyuk_seifert_cyclic_speed_dissipation()["koyuk_seifert_cyclic_speed_bridge_verified"])
        self.assertTrue(evaluate_falasco_esposito_dissipation_time()["falasco_esposito_dissipation_time_bridge_verified"])
        self.assertTrue(evaluate_dechant_sasa_fluctuation_response()["dechant_sasa_fluctuation_response_bridge_verified"])
        self.assertTrue(evaluate_lebowitz_spohn_action()["lebowitz_spohn_action_bridge_verified"])
        self.assertTrue(evaluate_speck_seifert_housekeeping()["speck_seifert_housekeeping_bridge_verified"])
        self.assertTrue(evaluate_polettini_esposito_cycle_decomposition()["polettini_esposito_cycle_bridge_verified"])
        self.assertTrue(evaluate_rao_esposito_conservation_laws()["rao_esposito_conservation_bridge_verified"])
        self.assertTrue(evaluate_manikandan_gupta_krishnamurthy_bounds()["manikandan_gupta_krishnamurthy_bridge_verified"])
        self.assertTrue(evaluate_wachtel_rao_esposito_coarse_graining()["wachtel_rao_esposito_coarse_graining_bridge_verified"])

        neg = evaluate_falsification_negative_control()
        self.assertTrue(neg["passed"])
        self.assertEqual(neg["rule"], "AH-10")


if __name__ == "__main__":
    unittest.main()
