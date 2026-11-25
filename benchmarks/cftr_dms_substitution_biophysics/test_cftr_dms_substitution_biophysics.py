#!/usr/bin/env python3
"""Unit and planted-defect tests for benchmarks/cftr_dms_substitution_biophysics/."""
from __future__ import annotations

import copy
import json
import math
import random
import subprocess
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from cftr_dms_substitution_biophysics import (  # noqa: E402
    BLOSUM62,
    BLOSUM62_ORDER,
    FEATURE_NAMES,
    KYTE_DOOLITTLE,
    RESIDUE_CHARGE,
    RESULTS_FILE_NAME,
    VariantSubstRecord,
    extract_7_biophysics_features,
    fit_partial_penalized_subst_ridge,
    load_p13569_context,
    parse_protein_distinction,
    parse_variant_components,
    parse_xlsx_sheet,
    residue_clustered_bootstrap_delta_ci,
    resolve_tedman_data_dir,
    verify_checked_in_artifacts,
    verify_p13569_wt_parity,
)
from independent_audit import (  # noqa: E402
    _AA_AXIS,
    _B62,
    _KD,
    _fit_partial_ridge_cholesky,
    cholesky_solve_spd,
    run_independent_audit,
)

BENCHMARK_SCRIPT = HERE / "cftr_dms_substitution_biophysics.py"
AUDIT_SCRIPT = HERE / "independent_audit.py"
REPO_ROOT = HERE.parents[1]


def _write_synthetic_sparse_xlsx(dest: Path) -> None:
    content_types = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
        '<Default Extension="xml" ContentType="application/xml"/>'
        '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
        '<Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
        '<Override PartName="/xl/sharedStrings.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sharedStrings+xml"/>'
        "</Types>"
    )
    rels = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>'
        "</Relationships>"
    )
    wb = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
        'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        '<sheets><sheet name="All_variants" sheetId="1" r:id="rId1"/></sheets>'
        "</workbook>"
    )
    wb_rels = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/>'
        "</Relationships>"
    )
    shared = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<sst xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" count="2" uniqueCount="2">'
        "<si><t>TMD1</t></si><si><t>R117H</t></si></sst>"
    )
    sheet1 = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        "<sheetData>"
        '<row r="1"><c r="A1" t="s"><v>0</v></c><c r="C1" t="s"><v>1</v></c><c r="F1"><v>2.5E-1</v></c></row>'
        '<row r="2"><c r="C2" t="inlineStr"><is><t>G551D</t></is></c><c r="F2"/><c r="R2"><v>1.5</v></c></row>'
        "</sheetData>"
        '<mergeCells count="1"><mergeCell ref="A1:A2"/></mergeCells>'
        "</worksheet>"
    )
    with zipfile.ZipFile(dest, "w") as zf:
        zf.writestr("[Content_Types].xml", content_types)
        zf.writestr("_rels/.rels", rels)
        zf.writestr("xl/workbook.xml", wb)
        zf.writestr("xl/_rels/workbook.xml.rels", wb_rels)
        zf.writestr("xl/sharedStrings.xml", shared)
        zf.writestr("xl/worksheets/sheet1.xml", sheet1)


class TestBiophysicsScalesAndVariantParser(unittest.TestCase):
    def test_blosum62_symmetry_and_published_anchors(self) -> None:
        self.assertEqual(len(BLOSUM62_ORDER), 20)
        self.assertEqual(set(BLOSUM62_ORDER), set(_AA_AXIS))
        for a in BLOSUM62_ORDER:
            for b in BLOSUM62_ORDER:
                self.assertEqual(BLOSUM62[a][b], BLOSUM62[b][a], msg=f"Asymmetry at ({a}, {b})")
                self.assertEqual(BLOSUM62[a][b], _B62[(a, b)], msg=f"Audit mismatch at ({a}, {b})")

        # Published Henikoff & Henikoff (1992) PNAS 89:10915 anchors
        expected_diag = {
            "A": 4, "R": 5, "N": 6, "D": 6, "C": 9,
            "Q": 5, "E": 5, "G": 6, "H": 8, "I": 4,
            "L": 4, "K": 5, "M": 5, "F": 6, "P": 7,
            "S": 4, "T": 5, "W": 11, "Y": 7, "V": 4,
        }
        for aa, val in expected_diag.items():
            self.assertEqual(BLOSUM62[aa][aa], val)
        self.assertEqual(BLOSUM62["H"]["G"], -2)
        self.assertEqual(BLOSUM62["D"]["E"], 2)
        self.assertEqual(BLOSUM62["K"]["R"], 2)
        self.assertEqual(BLOSUM62["F"]["Y"], 3)
        self.assertEqual(BLOSUM62["I"]["V"], 3)
        self.assertEqual(BLOSUM62["G"]["I"], -4)
        self.assertEqual(BLOSUM62["F"]["P"], -4)
        self.assertEqual(BLOSUM62["M"]["V"], 1)

    def test_kyte_doolittle_and_charge_anchors(self) -> None:
        self.assertEqual(set(KYTE_DOOLITTLE.keys()), set(BLOSUM62_ORDER))
        for aa in BLOSUM62_ORDER:
            self.assertAlmostEqual(KYTE_DOOLITTLE[aa], _KD[aa], places=12)
        self.assertAlmostEqual(KYTE_DOOLITTLE["I"], 4.5, places=12)
        self.assertAlmostEqual(KYTE_DOOLITTLE["V"], 4.2, places=12)
        self.assertAlmostEqual(KYTE_DOOLITTLE["L"], 3.8, places=12)
        self.assertAlmostEqual(KYTE_DOOLITTLE["F"], 2.8, places=12)
        self.assertAlmostEqual(KYTE_DOOLITTLE["M"], 1.9, places=12)
        self.assertAlmostEqual(KYTE_DOOLITTLE["A"], 1.8, places=12)
        self.assertAlmostEqual(KYTE_DOOLITTLE["G"], -0.4, places=12)
        self.assertAlmostEqual(KYTE_DOOLITTLE["P"], -1.6, places=12)
        self.assertAlmostEqual(KYTE_DOOLITTLE["H"], -3.2, places=12)
        self.assertAlmostEqual(KYTE_DOOLITTLE["K"], -3.9, places=12)
        self.assertAlmostEqual(KYTE_DOOLITTLE["R"], -4.5, places=12)

        self.assertEqual(RESIDUE_CHARGE["K"], 1.0)
        self.assertEqual(RESIDUE_CHARGE["R"], 1.0)
        self.assertEqual(RESIDUE_CHARGE["D"], -1.0)
        self.assertEqual(RESIDUE_CHARGE["E"], -1.0)
        self.assertEqual(RESIDUE_CHARGE["H"], 0.0)

    def test_compound_missense_q359k_t360k_and_logical_or_topology(self) -> None:
        p13569_ctx = load_p13569_context(REPO_ROOT)
        comps = parse_variant_components("Q359K/T360K")
        self.assertEqual(len(comps), 2)
        self.assertEqual((comps[0].wt_aa, comps[0].codon, comps[0].mut_aa), ("Q", 359, "K"))
        self.assertEqual((comps[1].wt_aa, comps[1].codon, comps[1].mut_aa), ("T", 360, "K"))
        self.assertEqual(parse_protein_distinction("Q359K/T360K"), (359, "other", "Q"))

        feats = extract_7_biophysics_features("Q359K/T360K", p13569_ctx)
        # Q->K: B62(Q,K)-B62(Q,Q)=1-5=-4, KD=-3.9-(-3.5)=-0.4, |dChg|=1
        # T->K: B62(T,K)-B62(T,T)=-1-5=-6, KD=-3.9-(-0.7)=-3.2, |dChg|=1
        self.assertAlmostEqual(feats[0], -10.0, places=12)
        self.assertAlmostEqual(feats[1], -3.6, places=12)
        self.assertAlmostEqual(feats[2], 2.0, places=12)
        self.assertEqual(feats[3], 0.0)  # neither Q,T,K in {P,G}
        self.assertEqual(feats[4], 0.0)  # neither 359 nor 360 in Transmembrane
        self.assertEqual(feats[5], 0.0)
        self.assertEqual(feats[6], 0.0)

        # Synthetic compound missense testing logical-OR across separate component positions:
        # position 85 is in Transmembrane & Helix (G85E: G->E gives helix_breaker=1 & tm_shock=1)
        # combined with position 3 in Lasso (R3W: neither in Helix nor TM).
        syn_ctx = {
            "tm_positions": frozenset({85}),
            "helix_or_tm_positions": frozenset({85}),
        }
        feats_or = extract_7_biophysics_features("R3W/G85E", syn_ctx)
        # R3W: dB62 = -3 - 5 = -8; dKD = -0.9 - (-4.5) = +3.6; |dChg| = |0 - 1| = 1
        # G85E: dB62 = -2 - 6 = -8; dKD = -3.5 - (-0.4) = -3.1; |dChg| = |-1 - 0| = 1
        self.assertAlmostEqual(feats_or[0], -16.0, places=12)
        self.assertAlmostEqual(feats_or[1], 0.5, places=12)
        self.assertAlmostEqual(feats_or[2], 2.0, places=12)
        self.assertEqual(feats_or[3], 1.0)
        self.assertEqual(feats_or[4], 1.0)
        self.assertEqual(feats_or[5], 0.0)
        self.assertEqual(feats_or[6], 0.0)

    def test_m470v_experimental_wt_delta_and_class_behaviors(self) -> None:
        p13569_ctx = load_p13569_context(REPO_ROOT)
        # Row 89 M470V: must use experimental WT M -> V, NOT canonical P13569 V470
        m470v = extract_7_biophysics_features("M470V", p13569_ctx)
        self.assertAlmostEqual(m470v[0], float(BLOSUM62["M"]["V"] - BLOSUM62["M"]["M"]), places=12)
        self.assertAlmostEqual(m470v[0], -4.0, places=12)
        self.assertAlmostEqual(m470v[1], 4.2 - 1.9, places=12)  # +2.3
        self.assertAlmostEqual(m470v[2], 0.0, places=12)
        self.assertEqual(m470v[3], 0.0)
        self.assertEqual(m470v[4], 0.0)
        self.assertEqual(m470v[5], 0.0)
        self.assertEqual(m470v[6], 0.0)

        # Synonymous variants -> all 7 features == 0.0
        for syn in ("E528E", "T854T", "Q966=", "L1156="):
            self.assertEqual(extract_7_biophysics_features(syn, p13569_ctx), (0.0,) * 7)

        # In-frame indels -> continuous substitution scales None, is_inframe_indel=1.0
        for indel in ("L138ins", "F312DEL", "I507del", "F508del"):
            self.assertEqual(
                extract_7_biophysics_features(indel, p13569_ctx),
                (None, None, None, 0.0, 0.0, 0.0, 1.0),
            )

        # Nonsense -> continuous substitution scales None, is_nonsense=1.0
        for stop in ("E60X", "W1282X", "Q1313X"):
            self.assertEqual(
                extract_7_biophysics_features(stop, p13569_ctx),
                (None, None, None, 0.0, 0.0, 1.0, 0.0),
            )


class TestRidgeFoldImputationAndCholesky(unittest.TestCase):
    def test_partial_penalized_ridge_recovery_and_cholesky_parity(self) -> None:
        rng = random.Random(123)
        xs_dmso = [-2.0 + 0.08 * i for i in range(50)]
        raw_feats = [
            (
                float(rng.randint(-10, -1)),
                rng.uniform(-6.0, 6.0),
                float(rng.choice((0, 1, 2))),
                1.0 if i % 5 == 0 else 0.0,
                1.0 if i % 6 == 0 else 0.0,
                1.0 if i % 11 == 0 else 0.0,
                1.0 if i % 13 == 0 else 0.0,
            )
            for i in range(50)
        ]

        # Case A: noiseless linear in x_dmso only; high alpha=1e5 must NOT shrink b0 or b1.
        ys_line = [1.25 + 0.85 * x for x in xs_dmso]
        mod_hi = fit_partial_penalized_subst_ridge(xs_dmso, raw_feats, ys_line, alpha=100000.0)
        self.assertAlmostEqual(mod_hi.b0_intercept, 1.25, places=6)
        self.assertAlmostEqual(mod_hi.b1_dmso_slope, 0.85, places=6)

        # Case B: planted feature weights recovered at near-zero regularization + Cholesky parity check.
        means = [sum(r[j] for r in raw_feats) / len(raw_feats) for j in range(7)]
        stds = [math.sqrt(sum((r[j] - means[j]) ** 2 for r in raw_feats) / len(raw_feats)) for j in range(7)]
        true_w = (0.35, -0.22, 0.14, -0.18, 0.27, -0.45, 0.12)
        ys_planted = [
            1.1 + 0.8 * xs_dmso[i] + sum(true_w[j] * ((raw_feats[i][j] - means[j]) / stds[j]) for j in range(7))
            for i in range(50)
        ]
        mod_lo = fit_partial_penalized_subst_ridge(xs_dmso, raw_feats, ys_planted, alpha=1e-10)
        self.assertAlmostEqual(mod_lo.b0_intercept, 1.1, places=6)
        self.assertAlmostEqual(mod_lo.b1_dmso_slope, 0.8, places=6)
        for w_hat, w_true in zip(mod_lo.weights, true_w):
            self.assertAlmostEqual(w_hat, w_true, places=6)

        beta_chol, imp_c, fm_c, fs_c = _fit_partial_ridge_cholesky(xs_dmso, raw_feats, ys_planted, alpha=10.0)
        mod_10 = fit_partial_penalized_subst_ridge(xs_dmso, raw_feats, ys_planted, alpha=10.0)
        self.assertAlmostEqual(beta_chol[0], mod_10.b0_intercept, places=9)
        self.assertAlmostEqual(beta_chol[1], mod_10.b1_dmso_slope, places=9)
        for w_c, w_g in zip(beta_chol[2:], mod_10.weights):
            self.assertAlmostEqual(w_c, w_g, places=9)
        for j in range(7):
            self.assertAlmostEqual(imp_c[j], mod_10.impute_means[j], places=9)
            self.assertAlmostEqual(fm_c[j], mod_10.feat_means[j], places=9)
            self.assertAlmostEqual(fs_c[j], mod_10.feat_stds[j], places=9)

        sol = cholesky_solve_spd([[4.0, 1.0], [1.0, 3.0]], [9.0, 5.0])
        self.assertAlmostEqual(sol[0], 2.0, places=9)
        self.assertAlmostEqual(sol[1], 1.0, places=9)

    def test_fold_imputation_zero_leakage(self) -> None:
        train_x = [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5]
        train_feats = [
            (-4.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0),
            (-6.0, -2.0, 0.0, 1.0, 0.0, 0.0, 0.0),
            (-8.0, 0.0, 2.0, 0.0, 1.0, 0.0, 0.0),
            (None, None, None, 0.0, 0.0, 1.0, 0.0),
            (None, None, None, 0.0, 0.0, 0.0, 1.0),
            (-2.0, 4.0, 1.0, 0.0, 0.0, 0.0, 0.0),
        ]
        train_y = [0.5 + 0.8 * x for x in train_x]
        mod = fit_partial_penalized_subst_ridge(train_x, train_feats, train_y, alpha=10.0)

        self.assertAlmostEqual(mod.impute_means[0], (-4.0 - 6.0 - 8.0 - 2.0) / 4.0, places=9)
        self.assertAlmostEqual(mod.impute_means[1], (2.0 - 2.0 + 0.0 + 4.0) / 4.0, places=9)
        self.assertAlmostEqual(mod.impute_means[2], (1.0 + 0.0 + 2.0 + 1.0) / 4.0, places=9)

        # Held-out nonsense row maps continuous features j in 0..2 to z_j = 0.0 without leaking fold statistics.
        z_miss = mod.transform_features((None, None, None, 0.0, 0.0, 1.0, 0.0))
        for j in range(3):
            self.assertAlmostEqual(z_miss[j], 0.0, places=12)

        val_pred_1 = mod.predict_one(0.25, (None, None, None, 0.0, 0.0, 1.0, 0.0))
        mod_2 = fit_partial_penalized_subst_ridge(train_x, train_feats, train_y, alpha=10.0)
        val_pred_2 = mod_2.predict_one(0.25, (None, None, None, 0.0, 0.0, 1.0, 0.0))
        self.assertEqual(val_pred_1, val_pred_2)


class TestSparseXlsxAndP13569Parity(unittest.TestCase):
    def test_sparse_xlsx_and_p13569_wt_parity_control(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            xp = Path(tmp) / "sparse.xlsx"
            _write_synthetic_sparse_xlsx(xp)
            cells, max_row = parse_xlsx_sheet(xp, "All_variants", expand_merged=True)
            self.assertEqual(max_row, 2)
            self.assertEqual(cells[(1, "A")], "TMD1")
            self.assertEqual(cells[(2, "A")], "TMD1")
            self.assertEqual(cells[(2, "C")], "G551D")
            self.assertAlmostEqual(float(cells[(1, "F")]), 0.25, places=9)
            self.assertNotIn((2, "F"), cells)

        p13569_ctx = load_p13569_context(REPO_ROOT)
        manifest = json.loads((REPO_ROOT / "benchmarks/cftr_dms_expression/split_manifest.json").read_text(encoding="utf-8"))
        rep = verify_p13569_wt_parity(manifest["assignments"], p13569_ctx)
        self.assertTrue(rep["passed"])
        self.assertEqual(rep["variants_checked"], 234)
        self.assertEqual(rep["components_checked"], 235)  # 233 single + 2 in Q359K/T360K
        self.assertEqual(rep["unexpected_mismatches"], 0)
        self.assertEqual(rep["holdout_mismatches"], 0)

        # Planted WT mismatch in holdout must fail closed
        bad_assigns = copy.deepcopy(manifest["assignments"])
        for row in bad_assigns:
            if row["split"] == "holdout":
                row["protein_distinction"] = "K" + str(row["codon_position"]) + "A"
                if p13569_ctx["sequence"][int(row["codon_position"]) - 1] != "K":
                    break
        with self.assertRaises(AssertionError):
            verify_p13569_wt_parity(bad_assigns, p13569_ctx)

    def test_cluster_bootstrap_counts_unique_codons(self) -> None:
        rows = [
            VariantSubstRecord(3, "TMD1", "1A>G", "R100W", 100, "missense", "holdout", False, 0.5, 1.0, -1.0, 0.0, (0.0,) * 7),
            VariantSubstRecord(4, "TMD1", "2A>T", "R100Q", 100, "missense", "holdout", False, 0.5, 1.0, -1.0, 0.0, (0.0,) * 7),
            VariantSubstRecord(5, "NBD1", "3G>A", "G500D", 500, "missense", "holdout", False, 0.25, 0.5, -2.0, -1.0, (0.0,) * 7),
        ]
        boot = residue_clustered_bootstrap_delta_ci(rows, [0.5, 0.5, -0.5], [0.0, 0.0, -1.0], resamples=50)
        self.assertEqual(boot["n_holdout_clusters"], 2)
        self.assertTrue(boot["ci95_low_gt_zero"])


class TestCheckedInArtifactsAndTamperGuards(unittest.TestCase):
    def test_verify_checked_in_and_independent_audit_cli(self) -> None:
        self.assertEqual(len(FEATURE_NAMES), 7)
        summary = verify_checked_in_artifacts(repo_root=REPO_ROOT)
        self.assertEqual(summary["status"], "PASS")
        self.assertEqual(summary["p13569_sequence_length"], 1480)
        self.assertEqual(summary["wt_components_checked"], 235)
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

    def test_tamper_done_and_results_fail_closed(self) -> None:
        results = json.loads((HERE / RESULTS_FILE_NAME).read_text(encoding="utf-8"))
        done_txt = (HERE / "done.md").read_text(encoding="utf-8")
        manifest_txt = (REPO_ROOT / "benchmarks/cftr_dms_expression/split_manifest.json").read_text(encoding="utf-8")
        th_txt = (REPO_ROOT / "benchmarks/theratype/theratype_open_annotations_2026-11-21.json").read_text(encoding="utf-8")
        attr_txt = (REPO_ROOT / "tools/cftr-model/data/ATTRIBUTION.md").read_text(encoding="utf-8")
        p135_bytes = (REPO_ROOT / "tools/cftr-model/data/p13569.json").read_bytes()

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bdir = root / "benchmarks" / "cftr_dms_substitution_biophysics"
            edir = root / "benchmarks" / "cftr_dms_expression"
            tdir = root / "benchmarks" / "theratype"
            mdir = root / "tools" / "cftr-model" / "data"
            for d in (bdir, edir, tdir, mdir):
                d.mkdir(parents=True)
            (edir / "split_manifest.json").write_text(manifest_txt, encoding="utf-8")
            (tdir / "theratype_open_annotations_2026-11-21.json").write_text(th_txt, encoding="utf-8")
            (mdir / "ATTRIBUTION.md").write_text(attr_txt, encoding="utf-8")
            (mdir / "p13569.json").write_bytes(p135_bytes)

            # Defect 1: tamper done.md -> both verifiers reject
            (bdir / "done.md").write_text(done_txt + "\n# tampered\n", encoding="utf-8")
            (bdir / RESULTS_FILE_NAME).write_text(json.dumps(results), encoding="utf-8")
            with self.assertRaises(AssertionError):
                verify_checked_in_artifacts(repo_root=root)
            with self.assertRaises(AssertionError):
                run_independent_audit(repo_root=root)

            # Defect 2: restore done.md, tamper p13569_sha256 in results -> both verifiers reject
            (bdir / "done.md").write_text(done_txt, encoding="utf-8")
            bad_res = copy.deepcopy(results)
            bad_res["p13569_sha256"] = "0" * 64
            (bdir / RESULTS_FILE_NAME).write_text(json.dumps(bad_res), encoding="utf-8")
            with self.assertRaises(AssertionError):
                verify_checked_in_artifacts(repo_root=root)
            with self.assertRaises(AssertionError):
                run_independent_audit(repo_root=root)

            # Defect 3: tamper relative_mae_reduction_vs_stronger -> audit rejects
            bad_res2 = copy.deepcopy(results)
            bad_res2["holdout_evaluation"]["relative_mae_reduction_vs_stronger"] = 0.25
            (bdir / RESULTS_FILE_NAME).write_text(json.dumps(bad_res2), encoding="utf-8")
            with self.assertRaises(AssertionError):
                run_independent_audit(repo_root=root)

            # Defect 4: restore results, plant PRED-0010 exposed codon (508) in holdout -> both reject
            (bdir / RESULTS_FILE_NAME).write_text(json.dumps(results), encoding="utf-8")
            bad_man = json.loads(manifest_txt)
            bad_man["assignments"][0]["split"] = "holdout"
            bad_man["assignments"][0]["codon_position"] = 508
            (edir / "split_manifest.json").write_text(json.dumps(bad_man), encoding="utf-8")
            with self.assertRaises(AssertionError):
                verify_checked_in_artifacts(repo_root=root)
            with self.assertRaises(AssertionError):
                run_independent_audit(repo_root=root)

    def test_ci_mode_without_xlsx_downloads_passes_offline_checks(self) -> None:
        with tempfile.TemporaryDirectory() as empty_dir:
            summary = verify_checked_in_artifacts(repo_root=REPO_ROOT, data_dir=Path(empty_dir))
            self.assertEqual(summary["status"], "PASS")
            self.assertFalse(summary["raw_local_data_parity_checked"])

            audit = run_independent_audit(repo_root=REPO_ROOT, data_dir=Path(empty_dir))
            self.assertEqual(audit["status"], "PASS")
            self.assertFalse(audit["raw_local_data_audited"])


if __name__ == "__main__":
    unittest.main()
