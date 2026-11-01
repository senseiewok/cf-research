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

        neg = evaluate_falsification_negative_control()
        self.assertTrue(neg["passed"])
        self.assertEqual(neg["rule"], "AH-10")


if __name__ == "__main__":
    unittest.main()
