#!/usr/bin/env python3
"""Unit tests for Coupled CFTR Variant -> Modality -> Delivery Deduction Engine."""
from __future__ import annotations

from pathlib import Path
import sys
import unittest

TOOLS_DIR = Path(__file__).resolve().parent
if str(TOOLS_DIR) not in sys.path:
    sys.path.insert(0, str(TOOLS_DIR))

import cftr_variant_to_delivery_deducer as cvd  # noqa: E402


class TestCftrVariantToDeliveryDeducer(unittest.TestCase):
    def test_w1282x_routes_to_suptrna_and_passes_coupled_rescue(self) -> None:
        res = cvd.deduce_end_to_end_rescue("CFTR-W1282X", pore_diameter_nm=140.0, f_mod_m1a=0.42)
        self.assertEqual(res["routed_modality"], "M1A_SUPTRNA_LNP_OR_MRNA_GENE")
        self.assertEqual(res["top_lnp_candidate"], "HLPN-PEG-80nm-m1A-suptRNA")
        self.assertTrue(res["coupled_rescue_threshold_met"])
        self.assertGreaterEqual(res["coupled_i_sc_equiv_pct"], 15.0)

    def test_splice_1585_routes_to_aso_or_mrna(self) -> None:
        res = cvd.deduce_end_to_end_rescue("CFTR-1585-1G>A")
        self.assertEqual(res["routed_modality"], "ASO_SPLICE_MODULATION_OR_MRNA")
        self.assertFalse(res["alphamissense_covered"])

    def test_falsification_controls_pass(self) -> None:
        falsify = cvd.evaluate_falsification_controls()
        self.assertEqual(falsify["status"], "PASS")
        self.assertEqual(falsify["controls_passed"], 3)

    def test_cli_main_passes(self) -> None:
        self.assertEqual(cvd.main(["--falsify"]), 0)


if __name__ == "__main__":
    unittest.main()
