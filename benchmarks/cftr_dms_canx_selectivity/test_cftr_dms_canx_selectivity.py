#!/usr/bin/env python3
"""Unit, negative-control, and artifact parity tests for benchmarks/cftr_dms_canx_selectivity/."""
from __future__ import annotations

import json
import math
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
import unittest
import zipfile

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import cftr_dms_canx_selectivity as bench  # noqa: E402
import independent_audit  # noqa: E402


def _make_dual_sheet_xlsx(path: Path, ko_row3_prot: str = "R3W") -> None:
    content_types = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
        '<Default Extension="xml" ContentType="application/xml"/>'
        '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
        '<Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
        '<Override PartName="/xl/worksheets/sheet2.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
        "</Types>"
    )
    rels = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>'
        "</Relationships>"
    )
    wb = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        "<sheets>"
        '<sheet name="All_variants" sheetId="1" r:id="rId1"/>'
        '<sheet name="All_variants CANX KO" sheetId="2" r:id="rId2"/>'
        "</sheets></workbook>"
    )
    wb_rels = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/>'
        '<Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet2.xml"/>'
        "</Relationships>"
    )

    def _sheet_xml(prot_name: str, include_j: bool) -> str:
        j_cell = '<c r="J3"><v>1.0</v></c>' if include_j else ""
        return (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
            "<sheetData>"
            '<row r="3">'
            '<c r="A3" t="inlineStr"><is><t>Lasso</t></is></c>'
            '<c r="B3" t="inlineStr"><is><t>7A&gt;T</t></is></c>'
            f'<c r="C3" t="inlineStr"><is><t>{prot_name}</t></is></c>'
            '<c r="F3"><v>0.5</v></c>'
            f"{j_cell}"
            '<c r="N3"><v>2.0</v></c>'
            '<c r="R3"><v>4.0</v></c>'
            "</row>"
            "</sheetData></worksheet>"
        )

    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr("[Content_Types].xml", content_types)
        zf.writestr("_rels/.rels", rels)
        zf.writestr("xl/workbook.xml", wb)
        zf.writestr("xl/_rels/workbook.xml.rels", wb_rels)
        zf.writestr("xl/worksheets/sheet1.xml", _sheet_xml("R3W", include_j=False))
        zf.writestr("xl/worksheets/sheet2.xml", _sheet_xml(ko_row3_prot, include_j=True))


class TestCanxSelectivityUnits(unittest.TestCase):
    def test_dampened_selectivity_bounded_near_zero(self) -> None:
        s_eps = bench.compute_selectivity_s_eps(0.05, 1e-6, eps=0.05)
        self.assertLess(s_eps, 1.01)
        self.assertAlmostEqual(bench.compute_selectivity_s_eps(1.0, 1.0, eps=0.05), 0.0, places=9)
        with self.assertRaises(ValueError):
            bench.compute_selectivity_s_eps(0.0, 1.0)

    def test_dual_sheet_mismatch_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "bad_dual.xlsx"
            _make_dual_sheet_xlsx(p, ko_row3_prot="P5L")
            fake_manifest = {
                3: {
                    "row_num": 3,
                    "protein_distinction": "R3W",
                    "sub_domain": "Lasso",
                    "codon_position": 3,
                    "variant_class": "missense",
                    "split": "dev",
                    "exposed_pred0010_forced": False,
                }
            }
            with self.assertRaises(ValueError):
                bench.load_dual_sheet_tedman_variants(p, fake_manifest)

    def test_bonferroni_975_bootstrap_uses_unique_clusters_and_planted_edge_protection(self) -> None:
        codons = [10, 10, 10, 20, 30, 40, 50, 60]
        targets = [1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0]
        base_p = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        cand_p = [0.8, 0.8, 0.8, 1.2, 1.2, 1.2, 1.2, 1.2]
        res = bench.residue_clustered_bootstrap_delta_cis(codons, targets, base_p, cand_p, resamples=400, seed=7)
        self.assertEqual(res["n_holdout_clusters"], 6)
        self.assertTrue(res["delta_mae_ci975_bonferroni"][0] <= res["delta_mae_ci95"][0])
        self.assertTrue(res["delta_mae_ci975_bonferroni"][1] >= res["delta_mae_ci95"][1])
        with self.assertRaises(ValueError):
            bench._safe_rel_mae_reduction(0.0, 0.1)

    def test_shuffled_negative_controls_fail_to_reach_gate(self) -> None:
        rng = random.Random(42)
        rows = []
        domains = ("TMD1", "NBD1", "TMD2", "NBD2")
        for idx in range(60):
            codon = 10 + idx
            sd = domains[idx % len(domains)]
            p_dmso = 0.25 + 0.02 * (idx % 15) + rng.uniform(0.0, 0.03)
            k_dmso = 0.20 + 0.75 * p_dmso + rng.uniform(0.0, 0.02)
            p_445 = (2.8 if sd in ("TMD1", "NBD1") else 1.8) * p_dmso
            p_combo = 1.3 * p_445
            k_661 = 1.4 * k_dmso
            k_445 = (2.9 if sd in ("TMD1", "NBD1") else 1.7) * k_dmso
            k_combo = 1.25 * (k_661 + k_445)
            rows.append(
                bench.DualSheetVariant(
                    row_num=3 + idx,
                    sub_domain=sd,
                    nucleic_acid_mutation=f"{idx}A>G",
                    protein_distinction=f"A{codon}V",
                    codon_position=codon,
                    variant_class="missense",
                    parental_dmso_rel=p_dmso,
                    parental_vx445_rel=p_445,
                    parental_combo_rel=p_combo,
                    canx_ko_dmso_rel=k_dmso,
                    canx_ko_vx661_rel=k_661,
                    canx_ko_vx445_rel=k_445,
                    canx_ko_combo_rel=k_combo,
                    x_dmso_log2=math.log2(p_dmso),
                    x_canx_ko_dmso_log2=math.log2(k_dmso),
                    y_combo_log2=math.log2(p_combo),
                    s_eps_canx_ko=bench.compute_selectivity_s_eps(k_445, k_661),
                    s_raw_canx_ko=math.log2(k_445 / k_661),
                    split="dev",
                    exposed_pred0010_forced=False,
                )
            )
        ctrl = bench.run_pre_holdout_controls(rows, len(rows), 0)
        self.assertTrue(ctrl["dev_shuffled_negative_controls"]["passed"])
        self.assertLess(ctrl["dev_shuffled_negative_controls"]["arm_a_shuffled_canx_dev_cv_rel_reduction"], 0.10)
        self.assertLess(ctrl["dev_shuffled_negative_controls"]["arm_b_shuffled_features_dev_cv_rel_reduction"], 0.10)

    def test_planted_tamper_rejected(self) -> None:
        results_path = HERE / "results_2026-11-23.json"
        if not results_path.is_file():
            self.skipTest("results_2026-11-23.json not yet written")
        with tempfile.TemporaryDirectory() as tmp:
            troot = Path(tmp) / "repo"
            shutil.copytree(HERE, troot / "benchmarks" / "cftr_dms_canx_selectivity")
            (troot / "benchmarks" / "cftr_dms_expression").mkdir(parents=True)
            (troot / "benchmarks" / "theratype").mkdir(parents=True)
            shutil.copy2(
                REPO_ROOT / "benchmarks" / "cftr_dms_expression" / "split_manifest.json",
                troot / "benchmarks" / "cftr_dms_expression" / "split_manifest.json",
            )
            shutil.copy2(
                REPO_ROOT / "benchmarks" / "theratype" / "theratype_open_annotations_2026-11-21.json",
                troot / "benchmarks" / "theratype" / "theratype_open_annotations_2026-11-21.json",
            )
            rfile = troot / "benchmarks" / "cftr_dms_canx_selectivity" / "results_2026-11-23.json"
            payload = json.loads(rfile.read_text(encoding="utf-8"))
            payload["done_sha256"] = "0" * 64
            rfile.write_text(json.dumps(payload, indent=2), encoding="utf-8")
            with self.assertRaises(AssertionError):
                independent_audit.run_independent_audit(repo_root=troot, data_dir=Path(tmp) / "missing")

    def test_checked_in_and_audit_cli(self) -> None:
        if not (HERE / "results_2026-11-23.json").is_file():
            self.skipTest("results_2026-11-23.json not yet written")
        p1 = subprocess.run(
            [sys.executable, str(HERE / "cftr_dms_canx_selectivity.py"), "--check", "--json"],
            capture_output=True,
            text=True,
            check=True,
        )
        d1 = json.loads(p1.stdout)
        self.assertEqual(d1["status"], "PASS")
        p2 = subprocess.run(
            [sys.executable, str(HERE / "independent_audit.py"), "--check", "--json"],
            capture_output=True,
            text=True,
            check=True,
        )
        d2 = json.loads(p2.stdout)
        self.assertEqual(d2["status"], "PASS")


if __name__ == "__main__":
    unittest.main()
