#!/usr/bin/env python3
"""Unit and safety tests for tools/trial_atlas/whole_person_visit_prep.py."""
from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import whole_person_visit_prep as wpvp  # noqa: E402


def _load_safe_output_checker():
    candidates = (
        REPO_ROOT.parent / "cf-skills" / ".claude" / "skills" / "cf-ai-safe-use" / "scripts" / "check_safe_output.py",
        REPO_ROOT.parents[2] / "cf-skills" / ".claude" / "skills" / "cf-ai-safe-use" / "scripts" / "check_safe_output.py",
    )
    for script in candidates:
        if script.is_file():
            spec = importlib.util.spec_from_file_location("check_safe_output", script)
            if spec and spec.loader:
                mod = importlib.util.module_from_spec(spec)
                sys.modules["check_safe_output"] = mod
                spec.loader.exec_module(mod)
                return mod
    return None


class TestWholePersonVisitPrep(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.safe_mod = _load_safe_output_checker()

    def test_1_english_guide_covers_five_domains_and_passes_safe_output(self) -> None:
        pkt = wpvp.build_whole_person_packet(domain="all", lang="en")
        self.assertEqual(len(pkt["sections"]), 5)
        md_en = wpvp.render_whole_person_markdown(pkt)
        self.assertIn("(911 in the US)", md_en)
        self.assertIn("Please check this with your CF care team before changing anything.", md_en)
        self.assertIn("not stated", md_en)
        self.assertIn("does not evaluate safety or efficacy", md_en)
        self.assertIn("13.0% Hispanic", md_en)
        self.assertIn("48%", md_en)

        if self.safe_mod is not None:
            flags = self.safe_mod.check(md_en)
            self.assertEqual(flags, [], f"check_safe_output.py flagged English guide: {flags}")

    def test_2_spanish_es_mx_guide_structure_and_citations(self) -> None:
        pkt = wpvp.build_whole_person_packet(domain="all", lang="es-MX")
        self.assertEqual(pkt["lang"], "es-MX")
        self.assertEqual(len(pkt["sections"]), 5)
        md_es = wpvp.render_whole_person_markdown(pkt)
        self.assertIn("(911 in the US)", md_es)
        self.assertIn("Esto no es un diagnóstico.", md_es)
        self.assertIn("NCT07312734", md_es)
        self.assertIn("intérprete médico calificado", md_es)

    def test_3_falsification_suite_and_phi_refusal(self) -> None:
        rep = wpvp.run_falsification_suite()
        self.assertTrue(rep["falsify_passed"], msg=str(rep))
        self.assertEqual(len(rep["negative_controls"]), 3)
        for c in rep["negative_controls"]:
            self.assertTrue(c["passed"], msg=str(c))
        with self.assertRaises(SystemExit) as cm:
            wpvp.main(["--domain", "invalid_domain"])
        self.assertEqual(cm.exception.code, 2)


if __name__ == "__main__":
    unittest.main()
