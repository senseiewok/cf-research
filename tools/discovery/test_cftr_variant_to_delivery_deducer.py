#!/usr/bin/env python3
"""Unit tests for Coupled CFTR Variant -> Modality -> Delivery Deduction Engine."""
from __future__ import annotations

import io
import json
from contextlib import redirect_stdout
from pathlib import Path
import sys
import unittest

TOOLS_DIR = Path(__file__).resolve().parent
if str(TOOLS_DIR) not in sys.path:
    sys.path.insert(0, str(TOOLS_DIR))

import cftr_variant_to_delivery_deducer as cvd  # noqa: E402
import cftr_theratype_benchmark as ctb  # noqa: E402


class TestCftrVariantToDeliveryDeducer(unittest.TestCase):
    def test_w1282x_and_g542x_produce_distinct_kinetics(self) -> None:
        w1282x = cvd.deduce_end_to_end_rescue("CFTR-W1282X", pore_diameter_nm=140.0, f_mod_m1a=0.42)
        g542x = cvd.deduce_end_to_end_rescue("CFTR-G542X", pore_diameter_nm=140.0, f_mod_m1a=0.42)
        self.assertEqual(w1282x["routed_modality"], "M1A_SUPTRNA_LNP_OR_MRNA_GENE")
        self.assertEqual(w1282x["top_lnp_candidate"], "HLPN-PEG-80nm-m1A-suptRNA")
        self.assertTrue(w1282x["coupled_rescue_threshold_met"])
        self.assertGreaterEqual(w1282x["coupled_i_sc_equiv_pct"], 15.0)
        self.assertNotEqual(w1282x["f_mrna0_used"], g542x["f_mrna0_used"])
        self.assertNotEqual(w1282x["coupled_i_sc_equiv_pct"], g542x["coupled_i_sc_equiv_pct"])
        self.assertIsNotNone(g542x["grch38_locus_context"])
        self.assertIn("recruiting_trial_count", w1282x["trial_coverage"])

    def test_route_variant_modality_all_branches(self) -> None:
        by_id = {c.variant_id: c for c in ctb.OPEN_CFTR_THERATYPE_BENCHMARK}
        self.assertEqual(cvd.route_variant_modality(by_id["CFTR-F508del"]), "SMALL_MOLECULE_MODULATOR_ELIGIBLE")
        self.assertEqual(cvd.route_variant_modality(by_id["CFTR-S549N"]), "SMALL_MOLECULE_MODULATOR_ELIGIBLE")
        self.assertEqual(cvd.route_variant_modality(by_id["CFTR-3849+10kbC>T"]), "ASO_SPLICE_MODULATION_OR_MODULATOR")
        self.assertEqual(cvd.route_variant_modality(by_id["CFTR-1585-1G>A"]), "ASO_SPLICE_MODULATION_OR_MRNA")
        self.assertEqual(cvd.route_variant_modality(by_id["CFTR-W1282X"]), "M1A_SUPTRNA_LNP_OR_MRNA_GENE")
        self.assertEqual(cvd.route_variant_modality(by_id["CFTR-2051delAA"]), "GENOTYPE_AGNOSTIC_MRNA_OR_GENE_THERAPY")

    def test_dornase_counterfactual_reports_strict_delta(self) -> None:
        cf = cvd.deduce_dornase_counterfactual("CFTR-W1282X")
        expected_delta = cf["dornase_220nm_feasible_lnps"] - cf["untreated_140nm_feasible_lnps"]
        self.assertEqual(cf["carrier_window_expansion_delta"], expected_delta)
        self.assertEqual(cf["dornase_expands_carrier_window"], expected_delta > 0)

    def test_collapsed_pore_and_invalid_inputs(self) -> None:
        collapsed = cvd.deduce_end_to_end_rescue("CFTR-W1282X", pore_diameter_nm=50.0)
        self.assertEqual(collapsed["feasible_lnp_count"], 0)
        self.assertIsNone(collapsed["top_lnp_candidate"])
        self.assertIsNone(collapsed["coupled_i_sc_equiv_pct"])

        with self.assertRaises(ValueError):
            cvd.deduce_end_to_end_rescue("CFTR-UNKNOWN")
        with self.assertRaises(ValueError):
            cvd.deduce_end_to_end_rescue("CFTR-W1282X", pore_diameter_nm=0.0)
        with self.assertRaises(ValueError):
            cvd.deduce_end_to_end_rescue("CFTR-W1282X", f_mod_m1a=1.5)

    def test_falsification_controls_and_cli_json(self) -> None:
        falsify = cvd.evaluate_falsification_controls()
        self.assertEqual(falsify["status"], "PASS")
        self.assertEqual(falsify["controls_passed"], 4)

        buf = io.StringIO()
        with redirect_stdout(buf):
            rc = cvd.main(["--variant", "CFTR-G542X", "--falsify", "--json"])
        self.assertEqual(rc, 0)
        parsed = json.loads(buf.getvalue())
        self.assertEqual(parsed["primary_deduction"]["variant_id"], "CFTR-G542X")
        self.assertEqual(parsed["falsification"]["status"], "PASS")


if __name__ == "__main__":
    unittest.main()
