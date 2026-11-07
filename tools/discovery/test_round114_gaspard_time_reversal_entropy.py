"""Unit tests for Round 114: Gaspard Time-Reversed Dynamical Entropy in CFTR Single-Channel Kinetics."""

import unittest
from tools.discovery.round114_gaspard_time_reversal_entropy import (
    CITATIONS,
    PARAMETER_PROVENANCE,
    PROFILES,
    ROUND_ID,
    TASK_ID,
    GaspardEntropyProfile,
    run_falsification,
    run_falsification_checks,
    simulate_gaspard_time_reversal_entropy,
)


class TestRound114GaspardTimeReversalEntropy(unittest.TestCase):
    def test_01_metadata_and_ah10_contract(self) -> None:
        self.assertEqual(ROUND_ID, "round114")
        self.assertEqual(TASK_ID, "T-0748")
        self.assertEqual(len(CITATIONS), 3)
        ah10 = run_falsification()
        self.assertTrue(ah10["passed"])

    def test_02_parameter_provenance_and_profiles(self) -> None:
        self.assertEqual(len(PROFILES), 5)
        for meta in PARAMETER_PROVENANCE.values():
            self.assertLess(meta["bounds"][0], meta["bounds"][1])

    def test_03_all_falsification_checks_pass(self) -> None:
        res = run_falsification_checks()
        self.assertTrue(res["passed"])
        for k, v in res["checks"].items():
            self.assertTrue(v, msg=k)

    def test_04_poisoned_negative_controls_rejected(self) -> None:
        r1 = simulate_gaspard_time_reversal_entropy(
            GaspardEntropyProfile("p1", "bad kl", 10.0, 1.0, 10.0, "negative_kl_rate")
        )
        self.assertFalse(r1.is_kl_non_negative)
        r2 = simulate_gaspard_time_reversal_entropy(
            GaspardEntropyProfile("p2", "bad id", 10.0, 1.0, 10.0, "broken_gaspard_identity")
        )
        self.assertFalse(r2.is_gaspard_identity_satisfied)
        r3 = simulate_gaspard_time_reversal_entropy(
            GaspardEntropyProfile("p3", "bad ks", 10.0, 1.0, 10.0, "negative_ks_entropy")
        )
        self.assertFalse(r3.is_forward_ks_positive)


if __name__ == "__main__":
    unittest.main()
