#!/usr/bin/env python3
"""Unit and planted-defect tests for benchmarks/cftr_dms_expression/ (Council Item #1)."""
from __future__ import annotations

import copy
import json
import math
import subprocess
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from cftr_dms_benchmark import (  # noqa: E402
    MANIFEST_FILE_NAME,
    RESULTS_FILE_NAME,
    VariantRecord,
    fit_partial_penalized_ridge,
    parse_protein_distinction,
    parse_xlsx_sheet,
    residue_clustered_bootstrap_delta_ci,
    resolve_tedman_data_dir,
    verify_checked_in_artifacts,
)
from independent_audit import run_independent_audit  # noqa: E402

BENCHMARK_SCRIPT = HERE / "cftr_dms_benchmark.py"
AUDIT_SCRIPT = HERE / "independent_audit.py"
REPO_ROOT = HERE.parents[1]


def _write_synthetic_sparse_xlsx(dest: Path) -> None:
    """Create a minimal .xlsx exercising sharedStrings, inlineStr, sci-notation, blank, absent <c>, merged."""
    content_types = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
        '<Default Extension="xml" ContentType="application/xml"/>'
        '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
        '<Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
        '<Override PartName="/xl/sharedStrings.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sharedStrings+xml"/>'
        '</Types>'
    )
    rels = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>'
        '</Relationships>'
    )
    wb = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        '<sheets><sheet name="TestSheet" sheetId="1" r:id="rId1"/></sheets>'
        '</workbook>'
    )
    wb_rels = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/>'
        '</Relationships>'
    )
    shared = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<sst xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" count="2" uniqueCount="2">'
        '<si><t>Lasso</t></si>'
        '<si><t>R3W</t></si>'
        '</sst>'
    )
    # Row 1: merged block A1:A2 = 'Lasso', shared string C1 = 'R3W', sci notation F1 = 1.25E-2
    # Row 2: A2 absent (filled by merge!), C2 uses inlineStr 'E60X', D2 absent (<c> omitted!), F2 = blank <c r="F2"/>, R2 = 0.75
    sheet1 = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        '<sheetData>'
        '<row r="1">'
        '<c r="A1" t="s"><v>0</v></c>'
        '<c r="C1" t="s"><v>1</v></c>'
        '<c r="F1"><v>1.25E-2</v></c>'
        '</row>'
        '<row r="2">'
        '<c r="C2" t="inlineStr"><is><t>E60X</t></is></c>'
        '<c r="F2"/>'
        '<c r="R2"><v>0.75</v></c>'
        '</row>'
        '</sheetData>'
        '<mergeCells count="1"><mergeCell ref="A1:A2"/></mergeCells>'
        '</worksheet>'
    )
    with zipfile.ZipFile(dest, "w") as zf:
        zf.writestr("[Content_Types].xml", content_types)
        zf.writestr("_rels/.rels", rels)
        zf.writestr("xl/workbook.xml", wb)
        zf.writestr("xl/_rels/workbook.xml.rels", wb_rels)
        zf.writestr("xl/sharedStrings.xml", shared)
        zf.writestr("xl/worksheets/sheet1.xml", sheet1)


class TestSparseXlsxAndVariantClassifier(unittest.TestCase):
    def test_sparse_xlsx_contract(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            xp = Path(tmp) / "sparse.xlsx"
            _write_synthetic_sparse_xlsx(xp)
            cells, max_row = parse_xlsx_sheet(xp, "TestSheet", expand_merged=True)
            self.assertEqual(max_row, 2)
            self.assertEqual(cells[(1, "A")], "Lasso")
            self.assertEqual(cells[(2, "A")], "Lasso")
            self.assertEqual(cells[(1, "C")], "R3W")
            self.assertEqual(cells[(2, "C")], "E60X")
            self.assertAlmostEqual(float(cells[(1, "F")]), 0.0125, places=9)
            self.assertNotIn((2, "F"), cells)
            self.assertNotIn((2, "D"), cells)
            self.assertAlmostEqual(float(cells[(2, "R")]), 0.75, places=9)

    def test_protein_distinction_classifier(self) -> None:
        self.assertEqual(parse_protein_distinction("R3W"), (3, "missense"))
        self.assertEqual(parse_protein_distinction("E60X"), (60, "nonsense"))
        self.assertEqual(parse_protein_distinction("W1282*"), (1282, "nonsense"))
        self.assertEqual(parse_protein_distinction("F508del"), (508, "inframe_indel"))
        self.assertEqual(parse_protein_distinction("F312DEL"), (312, "inframe_indel"))
        self.assertEqual(parse_protein_distinction("L138ins"), (138, "inframe_indel"))
        self.assertEqual(parse_protein_distinction("Q359K/T360K"), (359, "other"))
        self.assertEqual(parse_protein_distinction("Q966="), (966, "other"))
        self.assertEqual(parse_protein_distinction("T854T"), (854, "other"))


class TestModelAndLeakageGuards(unittest.TestCase):
    def test_partial_penalized_ridge_does_not_shrink_dmso_slope(self) -> None:
        # True noiseless relation: y = 1.5 + 0.9 * x_dmso + domain/class shift; alpha=1e4 shouldn't shrink b1.
        rows = []
        domains = ("Lasso", "TMD1", "NBD1", "R Domain", "TMD2", "NBD2")
        for i in range(30):
            xd = -1.5 + 0.1 * i
            sd = domains[i % len(domains)]
            yc = 1.5 + 0.9 * xd
            rows.append(
                VariantRecord(
                    row_num=i + 3,
                    sub_domain=sd,
                    nucleic_acid_mutation=f"{i}A>G",
                    protein_distinction=f"R{i+1}W",
                    codon_position=i + 1,
                    variant_class="missense",
                    dmso_rel_mut_wt=2.0 ** xd,
                    combo_rel_mut_wt=2.0 ** yc,
                    x_dmso_log2=xd,
                    y_combo_log2=yc,
                )
            )
        mod = fit_partial_penalized_ridge(rows, alpha=10000.0)
        self.assertAlmostEqual(mod.b0_intercept, 1.5, places=6)
        self.assertAlmostEqual(mod.b1_dmso_slope, 0.9, places=6)

    def test_cluster_bootstrap_counts_unique_codons(self) -> None:
        rows = [
            VariantRecord(3, "TMD1", "1A>G", "R100W", 100, "missense", 0.5, 1.0, -1.0, 0.0, "holdout"),
            VariantRecord(4, "TMD1", "2A>T", "R100Q", 100, "missense", 0.5, 1.0, -1.0, 0.0, "holdout"),
            VariantRecord(5, "NBD1", "3G>A", "G500D", 500, "missense", 0.25, 0.5, -2.0, -1.0, "holdout"),
        ]
        boot = residue_clustered_bootstrap_delta_ci(rows, [0.5, 0.5, -0.5], [0.0, 0.0, -1.0], resamples=50)
        self.assertEqual(boot["n_holdout_clusters"], 2)
        self.assertTrue(boot["ci95_low_gt_zero"])

    def test_planted_split_leakage_and_tamper_rejected(self) -> None:
        manifest = json.loads((HERE / MANIFEST_FILE_NAME).read_text(encoding="utf-8"))
        results = json.loads((HERE / RESULTS_FILE_NAME).read_text(encoding="utf-8"))
        done_txt = (HERE / "done.md").read_text(encoding="utf-8")
        th_txt = (REPO_ROOT / "benchmarks/theratype/theratype_open_annotations_2026-11-21.json").read_text(encoding="utf-8")

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bdir = root / "benchmarks" / "cftr_dms_expression"
            tdir = root / "benchmarks" / "theratype"
            bdir.mkdir(parents=True)
            tdir.mkdir(parents=True)
            (tdir / "theratype_open_annotations_2026-11-21.json").write_text(th_txt, encoding="utf-8")
            (bdir / "done.md").write_text(done_txt, encoding="utf-8")
            (bdir / RESULTS_FILE_NAME).write_text(json.dumps(results), encoding="utf-8")

            # Defect 1: put an exposed PRED-0010 codon (508) into holdout.
            bad_m = copy.deepcopy(manifest)
            bad_m["assignments"][0]["split"] = "holdout"
            bad_m["assignments"][0]["codon_position"] = 508
            (bdir / MANIFEST_FILE_NAME).write_text(json.dumps(bad_m), encoding="utf-8")
            with self.assertRaises(AssertionError):
                verify_checked_in_artifacts(repo_root=root)


class TestCheckedInArtifactsAndAudit(unittest.TestCase):
    def test_verify_checked_in_and_independent_audit_cli(self) -> None:
        summary = verify_checked_in_artifacts(repo_root=REPO_ROOT)
        self.assertEqual(summary["status"], "PASS")
        self.assertGreaterEqual(summary["holdout_unique_codons"], 40)

        audit = run_independent_audit(repo_root=REPO_ROOT)
        self.assertEqual(audit["status"], "PASS")
        if resolve_tedman_data_dir(repo_root=REPO_ROOT) is not None:
            self.assertTrue(summary["raw_local_data_parity_checked"])
            self.assertTrue(audit["raw_local_data_audited"])

        for script in (BENCHMARK_SCRIPT, AUDIT_SCRIPT):
            res = subprocess.run(
                [sys.executable, str(script), "--check", "--json"],
                capture_output=True,
                text=True,
            )
            self.assertEqual(res.returncode, 0, msg=res.stderr)
            payload = json.loads(res.stdout)
            self.assertEqual(payload["status"], "PASS")


if __name__ == "__main__":
    unittest.main()
