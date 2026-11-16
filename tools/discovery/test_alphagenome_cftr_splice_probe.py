#!/usr/bin/env python3
"""Unit tests for AlphaGenome CFTR Splicing & Regulatory Variant Probe."""
from __future__ import annotations

from pathlib import Path
import sys
import unittest

TOOLS_DIR = Path(__file__).resolve().parent
if str(TOOLS_DIR) not in sys.path:
    sys.path.insert(0, str(TOOLS_DIR))

import alphagenome_cftr_splice_probe as agp  # noqa: E402


class TestAlphagenomeCftrSpliceProbe(unittest.TestCase):
    def test_panel_coordinates_valid_grch38(self) -> None:
        panel = agp.evaluate_panel()
        self.assertEqual(panel["panel_size"], 5)

    def test_avi_json_parser_accepts_molecular_summary(self) -> None:
        parsed = agp.parse_alphagenome_avi_json(
            {
                "variant": "chr7:117642553:C>T",
                "avi_phred": 28.4,
                "top_modality": "splicing",
                "molecular_summary": "Predicted cryptic splice donor activation in CFTR intron 22.",
            }
        )
        self.assertEqual(parsed["avi_phred"], 28.4)
        self.assertEqual(parsed["top_modality"], "splicing")

    def test_falsification_controls_pass(self) -> None:
        falsify = agp.evaluate_falsification_controls()
        self.assertEqual(falsify["status"], "PASS")
        self.assertEqual(falsify["controls_passed"], 3)

    def test_cli_main_passes(self) -> None:
        self.assertEqual(agp.main(["--falsify"]), 0)


if __name__ == "__main__":
    unittest.main()
