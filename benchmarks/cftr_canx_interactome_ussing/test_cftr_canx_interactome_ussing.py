#!/usr/bin/env python3
"""Unit, statistical calibration, and planted-defect tests for benchmarks/cftr_canx_interactome_ussing/."""
from __future__ import annotations

import itertools
import json
import math
from pathlib import Path
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

import cftr_canx_interactome_ussing as bench  # noqa: E402
import independent_audit  # noqa: E402


def _write_minimal_xlsx(path: Path, sheets: dict[str, list[tuple[int, dict[str, tuple[str, str]]]]]) -> None:
    """Create a minimal valid OOXML `.xlsx` file with sparse cell coordinates."""
    content_types = [
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">',
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-officedocument.relationships+xml"/>',
        '<Default Extension="xml" ContentType="application/xml"/>',
        '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>',
    ]
    wb_sheets = []
    wb_rels = []
    for idx, sname in enumerate(sheets, start=1):
        content_types.append(
            f'<Override PartName="/xl/worksheets/sheet{idx}.xml" '
            'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
        )
        wb_sheets.append(f'<sheet name="{sname}" sheetId="{idx}" r:id="rId{idx}"/>')
        wb_rels.append(
            f'<Relationship Id="rId{idx}" '
            'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" '
            f'Target="worksheets/sheet{idx}.xml"/>'
        )
    content_types.append("</Types>")

    rels_root = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>'
        "</Relationships>"
    )
    wb_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        f"<sheets>{''.join(wb_sheets)}</sheets></workbook>"
    )
    wb_rels_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        f"{''.join(wb_rels)}</Relationships>"
    )

    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr("[Content_Types].xml", "".join(content_types))
        zf.writestr("_rels/.rels", rels_root)
        zf.writestr("xl/workbook.xml", wb_xml)
        zf.writestr("xl/_rels/workbook.xml.rels", wb_rels_xml)
        for idx, (_, rows) in enumerate(sheets.items(), start=1):
            row_xml = []
            for rnum, cols in rows:
                c_xml = []
                for col_let, (ctype, val) in cols.items():
                    coord = f"{col_let}{rnum}"
                    if ctype == "s":
                        c_xml.append(f'<c r="{coord}" t="inlineStr"><is><t>{val}</t></is></c>')
                    elif ctype == "e":
                        c_xml.append(f'<c r="{coord}" t="e"><v>{val}</v></c>')
                    else:
                        c_xml.append(f'<c r="{coord}"><v>{val}</v></c>')
                row_xml.append(f'<row r="{rnum}">{"".join(c_xml)}</row>')
            s_xml = (
                '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
                f"<sheetData>{''.join(row_xml)}</sheetData></worksheet>"
            )
            zf.writestr(f"xl/worksheets/sheet{idx}.xml", s_xml)


def _make_synthetic_ussing_xlsx(path: Path, *, dup_in_combo: bool = True, partial_col: bool = False) -> None:
    rows = [
        (1, {"A": ("s", "Sample ID"), "B": ("s", "Acute drugs"), "C": ("s", "Replicates"), "M": ("s", "AVG")}),
        (2, {col: ("n", str(i)) for i, col in enumerate("CDEFGHIJKL", start=1)}),
    ]
    for idx, gname in enumerate(bench.EXPECTED_USSING_GROUPS):
        r0 = 3 + 3 * idx
        if gname == "siCANX, 661+445" and dup_in_combo:
            fsk_vals = {"C": "22.5", "D": "22.5", "E": "12.8", "F": "20.8"}
            vx_vals = {"C": "8.4", "D": "8.4", "E": "6.78", "F": "1.4"}
            inh_vals = {"C": "-29.0", "D": "-29.0", "E": "-19.0", "F": "-21.5"}
        else:
            fsk_vals = {"C": "10.0", "D": "12.0", "E": "11.0", "F": "13.0"}
            vx_vals = {"C": "1.0", "D": "1.5", "E": "1.2", "F": "1.1"}
            inh_vals = {"C": "-10.5", "D": "-13.0", "E": "-11.8", "F": "-13.5"}
        if gname == "NT" and partial_col:
            del inh_vals["D"]
        rows.append((r0, {"A": ("s", gname), "B": ("s", "Fsk"), **{c: ("n", v) for c, v in fsk_vals.items()}}))
        rows.append((r0 + 1, {"B": ("s", "VX-770"), **{c: ("n", v) for c, v in vx_vals.items()}}))
        rows.append((r0 + 2, {"B": ("s", "Inh172"), **{c: ("n", v) for c, v in inh_vals.items()}}))
    _write_minimal_xlsx(path, {"Raw-data-analyzed": rows})


def _make_synthetic_proteomics_xlsx(path: Path, *, div0_count: int = 1, cftr_val: str = "0") -> None:
    rows = [
        (2, {"A": ("s", "Gene Symbol"), "B": ("s", "Pathway"), "C": ("s", "ΔF\xa0(-CX)"), "D": ("s", "DelF (-CX) +445"), "E": ("s", "DelF (-CX) +661"), "F": ("s", "DelF (-CX) + combo")}),
        (3, {"A": ("s", "CFTR"), "B": ("s", "Ion Transport"), "C": ("n", cftr_val), "D": ("n", "0"), "E": ("n", "0"), "F": ("n", "0")}),
    ]
    for i in range(1, 9):
        rnum = 3 + i
        c_cell = ("e", "#DIV/0!") if i <= div0_count else ("n", str(0.4 + 0.05 * i))
        rows.append(
            (
                rnum,
                {
                    "A": ("s", f"GENE{i}"),
                    "B": ("s", "Folding" if i <= 4 else "Metabolism"),
                    "C": c_cell,
                    "D": ("n", str(0.1 * i)),
                    "E": ("n", str(-0.1 * i)),
                    "F": ("n", str(0.2 + 0.02 * i)),
                },
            )
        )
    _write_minimal_xlsx(path, {"Sheet1": rows})


class TestCanxInteractomeUssingUnits(unittest.TestCase):
    def test_sparse_column_and_duplicate_dedup_handling(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "ussing.xlsx"
            _make_synthetic_ussing_xlsx(p, dup_in_combo=True, partial_col=False)
            groups, dups = bench.parse_ussing_workbook(p, "SYNTH", "M")
            self.assertEqual(len(dups), 1)
            self.assertEqual(dups[0]["group"], "siCANX, 661+445")
            self.assertEqual((dups[0]["kept_col"], dups[0]["dropped_col"]), ("C", "D"))
            self.assertEqual(len(groups["siCANX, 661+445"].as_deposited), 4)
            self.assertEqual(len(groups["siCANX, 661+445"].deduplicated), 3)

            bad = Path(tmp) / "partial.xlsx"
            _make_synthetic_ussing_xlsx(bad, dup_in_combo=False, partial_col=True)
            with self.assertRaises(ValueError):
                bench.parse_ussing_workbook(bad, "BAD", "M")

    def test_exact_permutation_and_welch_t_calibration(self) -> None:
        # Known exact 3 vs 3 extreme split: 2 / C(6,3) = 2 / 20 = 0.1
        p_ext, ge, tot = bench.exact_two_sample_permutation_p([1.0, 2.0, 3.0], [10.0, 11.0, 12.0])
        self.assertEqual((ge, tot), (2, 20))
        self.assertAlmostEqual(p_ext, 0.1, places=12)
        self.assertAlmostEqual(
            p_ext,
            independent_audit._bitmask_perm_two_sided([1.0, 2.0, 3.0], [10.0, 11.0, 12.0]),
            places=12,
        )

        # Stratified DiD exact brute-force verification against fast bisect & independent two-pointer
        fc, fs = [20.0, 22.0, 18.0], [34.0, 36.0, 38.0]
        wc, ws = [60.0, 65.0, 70.0], [85.0, 90.0, 95.0]
        p_1s, p_2s, n_splits = bench.exact_stratified_did_permutation_p(fc, fs, wc, ws)
        self.assertEqual(n_splits, 400)
        p_ind_1s = independent_audit._two_pointer_did_perm(fc, fs, wc, ws)
        self.assertAlmostEqual(p_1s, p_ind_1s, places=12)

        # Brute force check
        obs_did = math.log((sum(fc) / 3) / (sum(fs) / 3)) - math.log((sum(wc) / 3) / (sum(ws) / 3))
        f_pool = fc + fs
        w_pool = wc + ws
        bf_1s = 0
        bf_2s = 0
        for cf in itertools.combinations(range(6), 3):
            lf = math.log((sum(f_pool[i] for i in cf) / 3) / ((sum(f_pool) - sum(f_pool[i] for i in cf)) / 3))
            for cw in itertools.combinations(range(6), 3):
                lw = math.log((sum(w_pool[i] for i in cw) / 3) / ((sum(w_pool) - sum(w_pool[i] for i in cw)) / 3))
                d = lf - lw
                if d <= obs_did + 1e-12:
                    bf_1s += 1
                if abs(d) >= abs(obs_did) - 1e-12:
                    bf_2s += 1
        self.assertAlmostEqual(p_1s, bf_1s / 400, places=12)
        self.assertAlmostEqual(p_2s, bf_2s / 400, places=12)

        # Welch t vs independent Simpson integration & textbook df=5 critical value (2.570581835636315 -> p=0.05)
        p_crit = bench.regularized_incomplete_beta(2.5, 0.5, 5.0 / (5.0 + 2.570581835636315 ** 2))
        self.assertAlmostEqual(p_crit, 0.05, places=9)
        wt = bench.welch_t_test(fc, fs)
        _, wt_ind_p = independent_audit._welch_p_simpson(fc, fs)
        self.assertAlmostEqual(wt["p_two_sided"], wt_ind_p, places=8)

    def test_wilcoxon_mwu_and_bh_fdr_parity(self) -> None:
        diffs = [-3.0, -2.0, -1.5, -1.5, 0.0, 0.5, -4.0]
        w_res = bench.wilcoxon_signed_rank_two_sided(diffs)
        w_ind = independent_audit._wilcoxon_counter_dp(diffs)
        self.assertEqual(w_res["n_nonzero"], 6)
        self.assertEqual(w_res["method"], "exact_dp")
        self.assertAlmostEqual(w_res["p_two_sided"], w_ind, places=12)

        grp = [-0.5, -0.4, -0.3, -0.3, 0.1]
        bg = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, -0.1]
        mwu = bench.mann_whitney_u_two_sided(grp, bg)
        mwu_ind = independent_audit._pairwise_mwu_p(grp, bg)
        self.assertAlmostEqual(mwu["p_two_sided"], mwu_ind, places=12)

        pvals = [0.01, 0.04, 0.03, 0.50]
        qs = bench.benjamini_hochberg(pvals)
        qs_ind = independent_audit._bh_q_decimal(pvals)
        for q1, q2, exp in zip(qs, qs_ind, [0.04, 0.053333333333, 0.053333333333, 0.50]):
            self.assertAlmostEqual(q1, exp, places=9)
            self.assertAlmostEqual(q2, exp, places=9)

        # Zero-variance MWU (e.g. constant |LFC| = 1.0) returns z_signed=0.0 and p=1.0 cleanly
        mwu_zv = bench.mann_whitney_u_two_sided([1.0, 1.0, 1.0], [1.0, 1.0, 1.0, 1.0])
        self.assertEqual(mwu_zv["z_signed"], 0.0)
        self.assertEqual(mwu_zv["p_two_sided"], 1.0)
        self.assertEqual(independent_audit._pairwise_mwu_p([1.0, 1.0, 1.0], [1.0, 1.0, 1.0, 1.0]), 1.0)

    def test_div0_nbsp_and_fail_closed_proteomics_guards(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            ok_p = Path(tmp) / "drug_ok.xlsx"
            _make_synthetic_proteomics_xlsx(ok_p, div0_count=2, cftr_val="0")
            rows, anchor, non_num = bench.parse_proteomics_workbook(ok_p, "DRUG", bench.DRUG_COLS)
            self.assertEqual(len(rows), 8)
            self.assertEqual(anchor.gene_symbol, "CFTR")
            self.assertEqual(len(non_num), 2)
            self.assertIsNone(rows[0].values["ΔF (-CX)"])

            # NBSP inside Pathway cell normalizes cleanly
            nbsp_p = Path(tmp) / "nbsp_pw.xlsx"
            _write_minimal_xlsx(
                nbsp_p,
                {
                    "Sheet1": [
                        (2, {"A": ("s", "Gene Symbol"), "B": ("s", "Pathway"), "C": ("s", "ΔF\xa0(-CX)"), "D": ("s", "DelF (-CX) +445"), "E": ("s", "DelF (-CX) +661"), "F": ("s", "DelF (-CX) + combo")}),
                        (3, {"A": ("s", "CFTR"), "B": ("s", "Ion Transport"), "C": ("n", "0"), "D": ("n", "0"), "E": ("n", "0"), "F": ("n", "0")}),
                        (4, {"A": ("s", "DNAJA1"), "B": ("s", "Protein\xa0folding"), "C": ("n", "0.5"), "D": ("n", "0.2"), "E": ("n", "0.3"), "F": ("n", "0.1")}),
                    ]
                },
            )
            nbsp_rows, _, _ = bench.parse_proteomics_workbook(nbsp_p, "DRUG", bench.DRUG_COLS)
            self.assertEqual(nbsp_rows[0].pathway_norm, "protein folding")

            # Permutation null handles a module with 0 finite observations in one variant column without ZeroDivisionError
            v_anchor = bench.ProteomicsRow(3, "CFTR", "Ion", "ion", {c: 0.0 for c in bench.VARIANT_COLS})
            v_rows = []
            for idx in range(12):
                vals = {c: 0.2 * (idx - 5) for c in bench.VARIANT_COLS}
                if idx == 0:
                    vals["I1366N (-CX)"] = None
                pw = "folding" if idx == 0 else "metabolism"
                v_rows.append(bench.ProteomicsRow(4 + idx, f"G{idx}", pw, pw, vals))
            d_anchor = bench.ProteomicsRow(3, "CFTR", "Ion", "ion", {c: 0.0 for c in bench.DRUG_COLS})
            d_rows = [
                bench.ProteomicsRow(4 + idx, f"D{idx}", "folding" if idx < 2 else "metabolism", "folding" if idx < 2 else "metabolism", {c: 0.2 * (idx + 1) for c in bench.DRUG_COLS})
                for idx in range(6)
            ]
            tmt_ok = {"runs": 6, "channels_per_run": 18, "total_channel_entries": 108, "channels": ["126"] * 18}
            orig_draws = bench.PERM_NULL_DRAWS
            try:
                bench.PERM_NULL_DRAWS = 50
                null_res = bench.run_arm2_controls(v_rows, v_anchor, d_rows, d_anchor, tmt_ok)
                self.assertTrue(null_res["pathway_permutation_null"]["passed"])
            finally:
                bench.PERM_NULL_DRAWS = orig_draws

            # > 5 non-numeric cells fails closed
            too_many = Path(tmp) / "drug_bad_cnt.xlsx"
            _make_synthetic_proteomics_xlsx(too_many, div0_count=6, cftr_val="0")
            with self.assertRaises(RuntimeError):
                bench.parse_proteomics_workbook(too_many, "DRUG", bench.DRUG_COLS)

            # Non-zero CFTR anchor fails dispersion/anchor control
            bad_anchor = Path(tmp) / "drug_bad_cftr.xlsx"
            _make_synthetic_proteomics_xlsx(bad_anchor, div0_count=0, cftr_val="0.5")
            r2, a2, _ = bench.parse_proteomics_workbook(bad_anchor, "DRUG", bench.DRUG_COLS)
            self.assertNotEqual(a2.values["ΔF (-CX)"], 0.0)

    def test_planted_tamper_rejected_and_offline_ci_check(self) -> None:
        results_path = HERE / bench.RESULTS_FILE_NAME
        if not results_path.is_file():
            self.skipTest(f"{bench.RESULTS_FILE_NAME} not yet generated")
        with tempfile.TemporaryDirectory() as tmp:
            troot = Path(tmp) / "repo"
            shutil.copytree(HERE, troot / "benchmarks" / "cftr_canx_interactome_ussing")
            rfile = troot / "benchmarks" / "cftr_canx_interactome_ussing" / bench.RESULTS_FILE_NAME
            lfile = troot / "benchmarks" / "cftr_canx_interactome_ussing" / "ledger.json"
            orig_payload = json.loads(rfile.read_text(encoding="utf-8"))
            orig_ledger = json.loads(lfile.read_text(encoding="utf-8"))

            # Offline CI verification (no data_dir and no local tedman-2025 under temp repo) passes first
            ci_rep = bench.verify_checked_in_artifacts(repo_root=troot, data_dir=None)
            self.assertEqual(ci_rep["status"], "PASS")
            self.assertFalse(ci_rep["raw_local_data_parity_checked"])
            audit_ci = independent_audit.run_independent_audit(repo_root=troot, data_dir=None)
            self.assertEqual(audit_ci["status"], "PASS")
            self.assertFalse(audit_ci["raw_local_data_audited"])

            # Explicit nonexistent --data-dir fails closed
            with self.assertRaises(FileNotFoundError):
                bench.verify_checked_in_artifacts(repo_root=troot, data_dir=Path(tmp) / "missing")
            with self.assertRaises(FileNotFoundError):
                independent_audit.run_independent_audit(repo_root=troot, data_dir=Path(tmp) / "missing")

            # Tamper done_sha256 -> both scripts fail closed
            p_hash = json.loads(json.dumps(orig_payload))
            p_hash["done_sha256"] = "0" * 64
            rfile.write_text(json.dumps(p_hash, indent=2), encoding="utf-8")
            with self.assertRaises(AssertionError):
                bench.verify_checked_in_artifacts(repo_root=troot, data_dir=None)
            with self.assertRaises(AssertionError):
                independent_audit.run_independent_audit(repo_root=troot, data_dir=None)

            # Tamper Arm 1/2 verdict + fake bh_fdr_q while keeping raw p intact -> both fail closed
            p_fake = json.loads(json.dumps(orig_payload))
            p_fake["arm_1_frt_ussing_sicanx_asymmetry"]["primary_strict_dedup"]["overall_verdict"] = "supported"
            p_fake["arm_1_frt_ussing_sicanx_asymmetry"]["primary_strict_dedup"]["ledger_standing"] = "open"
            p_fake["arm_2_canx_ko_proteostasis_remodeling"]["primary_evaluation"]["family_2a_drug_attenuation"]["tests"][0]["bh_fdr_q"] = 0.01
            p_fake["arm_2_canx_ko_proteostasis_remodeling"]["primary_evaluation"]["overall_verdict"] = "supported"
            rfile.write_text(json.dumps(p_fake, indent=2), encoding="utf-8")
            l_fake = json.loads(json.dumps(orig_ledger))
            for cl in l_fake["claims"]:
                if cl["id"] in ("c4-arm-1-ussing-sicanx-asymmetry", "c5-arm-2-proteostasis-remodeling"):
                    cl["standing"] = "open"
            lfile.write_text(json.dumps(l_fake, indent=2), encoding="utf-8")
            with self.assertRaises(AssertionError):
                bench.verify_checked_in_artifacts(repo_root=troot, data_dir=None)
            with self.assertRaises(AssertionError):
                independent_audit.run_independent_audit(repo_root=troot, data_dir=None)

    def test_checked_in_and_audit_cli(self) -> None:
        if not (HERE / bench.RESULTS_FILE_NAME).is_file():
            self.skipTest(f"{bench.RESULTS_FILE_NAME} not yet generated")
        p1 = subprocess.run(
            [sys.executable, str(HERE / "cftr_canx_interactome_ussing.py"), "--check", "--json"],
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
