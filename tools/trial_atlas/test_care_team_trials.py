#!/usr/bin/env python3
"""Unit tests for tools/trial_atlas/care_team_trials.py, including cf-ai-safe-use verification."""

from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import care_team_trials as ctt  # noqa: E402

SAFE_OUTPUT_SCRIPT = (
    REPO_ROOT.parent / "cf-skills" / ".claude" / "skills" / "cf-ai-safe-use" / "scripts" / "check_safe_output.py"
)


def _load_safe_output_checker():
    if not SAFE_OUTPUT_SCRIPT.exists():
        return None
    spec = importlib.util.spec_from_file_location("check_safe_output", SAFE_OUTPUT_SCRIPT)
    if spec is None or spec.loader is None:
        return None
    mod = importlib.util.module_from_spec(spec)
    sys.modules["check_safe_output"] = mod
    spec.loader.exec_module(mod)
    return mod


class TestCareTeamTrials(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.meta, cls.included = ctt.load_prisma_included_studies()
        cls.safe_mod = _load_safe_output_checker()

    def test_prisma_included_excludes_non_cf_bronchiectasis(self) -> None:
        # NCT01515007 is Non-CF Bronchiectasis (ORBIT-3) and must be excluded by PRISMA X5
        ncts = {s["nct_id"] for s in self.included}
        self.assertNotIn("NCT01515007", ncts)
        self.assertEqual(len(self.included), 940)

    def test_filter_by_status_and_modality(self) -> None:
        mrna_trials = ctt.filter_registered_trials(self.included, status="ALL", modality="mRNA", limit=50)
        self.assertGreaterEqual(len(mrna_trials), 1)
        for t in mrna_trials:
            self.assertTrue(any("mrna" in m.lower() for m in t["modalities"]))
            self.assertTrue(t["nct_id"].startswith("NCT"))

    def test_refuses_personal_identifiers(self) -> None:
        with self.assertRaises(ctt.RefusedInputError):
            ctt.filter_registered_trials(self.included, query="born 05/12/2014")
        with self.assertRaises(ctt.RefusedInputError):
            ctt.filter_registered_trials(self.included, query="patient@example.com")
        with self.assertRaises(ctt.RefusedInputError):
            ctt.filter_registered_trials(self.included, query="MRN 12345678")

    def test_english_sheet_passes_cf_ai_safe_use_checker(self) -> None:
        trials = ctt.filter_registered_trials(self.included, status="RECRUITING", limit=8)
        sheet_en = ctt.render_care_team_sheet(
            trials,
            {"status": "RECRUITING", "modality": "All modalities", "phase": "All phases", "query": "None"},
            self.meta.get("data_timestamp", "2026-10-14"),
            lang="en",
        )
        self.assertIn("Please check this with your CF care team before changing anything.", sheet_en)
        self.assertIn("(911 in the US)", sheet_en)
        self.assertNotIn("may be eligible", sheet_en.lower())
        self.assertNotIn("could qualify", sheet_en.lower())

        if self.safe_mod is not None:
            flags = self.safe_mod.check_text(sheet_en)
            self.assertEqual(flags, [], f"check_safe_output.py flagged English sheet: {flags}")

    def test_spanish_es_mx_sheet_contains_mandatory_safety_elements(self) -> None:
        trials = ctt.filter_registered_trials(self.included, status="RECRUITING", modality="CFTR Modulator", limit=5)
        sheet_es = ctt.render_care_team_sheet(
            trials,
            {"status": "RECRUITING", "modality": "CFTR Modulator", "phase": "Todas las fases", "query": "Ninguno"},
            self.meta.get("data_timestamp", "2026-10-14"),
            lang="es-MX",
        )
        self.assertIn("(911 in the US)", sheet_es)
        self.assertIn("equipo de atención de FQ", sheet_es)
        self.assertIn("intérprete médico calificado", sheet_es)
        self.assertIn("Esto no es un diagnóstico.", sheet_es)


if __name__ == "__main__":
    unittest.main()
