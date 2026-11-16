#!/usr/bin/env python3
"""Unit tests for Two-Stage AGILE + Mucus Barrier LNP Candidate Ranker."""
from __future__ import annotations

from pathlib import Path
import sys
import unittest

TOOLS_DIR = Path(__file__).resolve().parent
if str(TOOLS_DIR) not in sys.path:
    sys.path.insert(0, str(TOOLS_DIR))

import agile_mucus_lnp_ranker as aml  # noqa: E402


class TestAgileMucusLnpRanker(unittest.TestCase):
    def test_ranking_gates_cationic_and_oversized_lnps(self) -> None:
        rep = aml.rank_lnp_candidates()
        self.assertEqual(rep["feasible_count"], 2)
        self.assertEqual(rep["rejected_count"], 2)
        self.assertEqual(rep["top_candidate"], "HLPN-PEG-80nm-m1A-suptRNA")

    def test_falsification_controls_pass(self) -> None:
        falsify = aml.evaluate_falsification_controls()
        self.assertEqual(falsify["status"], "PASS")
        self.assertEqual(falsify["controls_passed"], 3)

    def test_cli_main_passes(self) -> None:
        self.assertEqual(aml.main(["--falsify"]), 0)


if __name__ == "__main__":
    unittest.main()
