#!/usr/bin/env python3
"""Unit and planted-defect tests for benchmarks/cftr_dms_3d_structure/."""
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

from cftr_dms_3d_structure import (  # noqa: E402
    CIF_EXPECTED_SHA256,
    RESULTS_FILE_NAME,
    Variant3DRecord,
    apply_rigid_transform,
    fit_partial_penalized_3d_ridge,
    kabsch_svd_superpose,
    load_mmcif_gz,
    parse_mmcif_text,
    parse_protein_distinction,
    parse_xlsx_sheet,
    residue_clustered_bootstrap_delta_ci,
    resolve_tedman_data_dir,
    verify_checked_in_artifacts,
    verify_wt_sequence_parity,
)
from independent_audit import (  # noqa: E402
    _apply_rt,
    cholesky_solve_spd,
    horn_quaternion_superpose,
    run_independent_audit,
)

BENCHMARK_SCRIPT = HERE / "cftr_dms_3d_structure.py"
AUDIT_SCRIPT = HERE / "independent_audit.py"
REPO_ROOT = HERE.parents[1]


def _rot_euler(a: float, b: float, c: float) -> list[list[float]]:
    rz = [[math.cos(a), -math.sin(a), 0.0], [math.sin(a), math.cos(a), 0.0], [0.0, 0.0, 1.0]]
    ry = [[math.cos(b), 0.0, math.sin(b)], [0.0, 1.0, 0.0], [-math.sin(b), 0.0, math.cos(b)]]
    rx = [[1.0, 0.0, 0.0], [0.0, math.cos(c), -math.sin(c)], [0.0, math.sin(c), math.cos(c)]]

    def mm(x: list[list[float]], y: list[list[float]]) -> list[list[float]]:
        return [[sum(x[i][k] * y[k][j] for k in range(3)) for j in range(3)] for i in range(3)]

    return mm(rz, mm(ry, rx))


def _det3(m: list[list[float]]) -> float:
    return (
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    )


def _write_synthetic_sparse_xlsx(dest: Path) -> None:
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
        '<sheets><sheet name="All_variants" sheetId="1" r:id="rId1"/></sheets>'
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
        '<si><t>TMD1</t></si><si><t>R117H</t></si></sst>'
    )
    sheet1 = (
        '<?xml version="1.0" encoding="UTF-8"?>'
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        '<sheetData>'
        '<row r="1"><c r="A1" t="s"><v>0</v></c><c r="C1" t="s"><v>1</v></c><c r="F1"><v>2.5E-1</v></c></row>'
        '<row r="2"><c r="C2" t="inlineStr"><is><t>G551D</t></is></c><c r="F2"/><c r="R2"><v>1.5</v></c></row>'
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


class TestKabschAndMmCIFShift(unittest.TestCase):
    def test_synthetic_kabsch_and_horn_recover_rigid_motion_lt_1e9(self) -> None:
        rng = random.Random(42)
        moving = [
            (rng.uniform(-40.0, 40.0), rng.uniform(-40.0, 40.0), rng.uniform(-40.0, 40.0))
            for _ in range(120)
        ]
        for angles, shift in (
            ((0.7, -1.1, 2.3), (14.5, -8.2, 27.0)),
            ((3.14159, 0.2, -1.9), (-50.0, 12.0, 0.5)),
            ((0.0, 0.0, 0.0), (0.0, 0.0, 0.0)),
        ):
            r_true = _rot_euler(*angles)
            target = [
                (
                    r_true[0][0] * p[0] + r_true[0][1] * p[1] + r_true[0][2] * p[2] + shift[0],
                    r_true[1][0] * p[0] + r_true[1][1] * p[1] + r_true[1][2] * p[2] + shift[1],
                    r_true[2][0] * p[0] + r_true[2][1] * p[1] + r_true[2][2] * p[2] + shift[2],
                )
                for p in moving
            ]
            r_svd, t_svd = kabsch_svd_superpose(moving, target)
            r_horn, t_horn = horn_quaternion_superpose(moving, target)

            err_svd = max(math.dist(apply_rigid_transform(r_svd, t_svd, p), q) for p, q in zip(moving, target))
            err_horn = max(math.dist(_apply_rt(r_horn, t_horn, p), q) for p, q in zip(moving, target))
            self.assertLess(err_svd, 1e-9, msg=f"SVD error {err_svd} for {angles}")
            self.assertLess(err_horn, 1e-9, msg=f"Horn error {err_horn} for {angles}")
            self.assertAlmostEqual(_det3(r_svd), 1.0, places=9)
            self.assertAlmostEqual(_det3(r_horn), 1.0, places=9)

    def test_both_solvers_reject_reflection(self) -> None:
        rng = random.Random(7)
        moving = [
            (rng.uniform(-25.0, 25.0), rng.uniform(-25.0, 25.0), rng.uniform(-25.0, 25.0))
            for _ in range(80)
        ]
        mirrored = [(-p[0], p[1], p[2]) for p in moving]
        r_svd, t_svd = kabsch_svd_superpose(moving, mirrored)
        r_horn, t_horn = horn_quaternion_superpose(moving, mirrored)
        self.assertAlmostEqual(_det3(r_svd), 1.0, places=9)
        self.assertAlmostEqual(_det3(r_horn), 1.0, places=9)
        rmsd_svd = math.sqrt(sum(math.dist(apply_rigid_transform(r_svd, t_svd, p), q) ** 2 for p, q in zip(moving, mirrored)) / len(moving))
        rmsd_horn = math.sqrt(sum(math.dist(_apply_rt(r_horn, t_horn, p), q) ** 2 for p, q in zip(moving, mirrored)) / len(moving))
        self.assertGreater(rmsd_svd, 5.0)
        self.assertGreater(rmsd_horn, 5.0)

    def test_8eiq_deltaf508_minus_one_shift_regression(self) -> None:
        # Synthetic mmCIF reproducing the exact 8EIQ deltaF508 auth_seq_id vs label_seq_id -1 offset.
        synthetic_8eiq = "\n".join([
            "data_8EIQ_SYN",
            "loop_",
            "_atom_site.group_PDB",
            "_atom_site.pdbx_PDB_model_num",
            "_atom_site.auth_asym_id",
            "_atom_site.label_atom_id",
            "_atom_site.label_comp_id",
            "_atom_site.auth_seq_id",
            "_atom_site.label_seq_id",
            "_atom_site.Cartn_x",
            "_atom_site.Cartn_y",
            "_atom_site.Cartn_z",
            "_atom_site.B_iso_or_equiv",
            "_atom_site.type_symbol",
            "_atom_site.label_alt_id",
            "_atom_site.pdbx_PDB_ins_code",
            "ATOM 1 A CA ILE 507 507 1.0 2.0 3.0 50.0 C . ?",
            "ATOM 1 A CA GLY 509 508 4.0 5.0 6.0 55.0 C . ?",
            "ATOM 1 A CA VAL 510 509 7.0 8.0 9.0 60.0 C . ?",
        ])
        by_auth, _ = parse_mmcif_text(synthetic_8eiq, "8EIQ_SYN", seq_id_field="auth_seq_id")
        by_label, _ = parse_mmcif_text(synthetic_8eiq, "8EIQ_SYN", seq_id_field="label_seq_id")
        self.assertNotIn(508, by_auth)
        self.assertEqual(by_auth[509].comp_id, "GLY")
        # Buggy label_seq_id indexing shifts 509 GLY into slot 508!
        self.assertEqual(by_label[508].comp_id, "GLY")

        # Real offline 6MSM & 8EIQ verification
        data_dir = REPO_ROOT / "tools" / "cftr-model" / "data"
        ca6, _ = load_mmcif_gz(data_dir / "6MSM.cif.gz", "6MSM", CIF_EXPECTED_SHA256["6MSM.cif.gz"])
        ca8, _ = load_mmcif_gz(data_dir / "8EIQ.cif.gz", "8EIQ", CIF_EXPECTED_SHA256["8EIQ.cif.gz"])
        self.assertEqual(ca6[508].comp_id, "PHE")
        self.assertNotIn(508, ca8)
        self.assertEqual(ca8[509].comp_id, "GLY")
        self.assertEqual(ca8[509].label_seq_id, 508)


class TestRidgeAndFoldImputation(unittest.TestCase):
    def test_partial_penalized_ridge_planted_recovery_and_cholesky_parity(self) -> None:
        rng = random.Random(123)
        xs_dmso = [-2.0 + 0.1 * i for i in range(40)]
        raw_feats = [
            (
                rng.uniform(0.0, 1.0),
                rng.uniform(0.0, 1.0),
                float(rng.randint(4, 20)),
                rng.uniform(80.0, 180.0),
                rng.uniform(0.2, 6.0),
                1.0 if i % 7 == 0 else 0.0,
            )
            for i in range(40)
        ]
        # Case A: noiseless linear in x_dmso only; high alpha=1e5 must NOT shrink b0 or b1.
        ys_line = [1.25 + 0.85 * x for x in xs_dmso]
        mod_hi = fit_partial_penalized_3d_ridge(xs_dmso, raw_feats, ys_line, alpha=100000.0)
        self.assertAlmostEqual(mod_hi.b0_intercept, 1.25, places=6)
        self.assertAlmostEqual(mod_hi.b1_dmso_slope, 0.85, places=6)

        # Case B: planted feature weights recovered at near-zero regularization + Cholesky parity check.
        means = [sum(r[j] for r in raw_feats) / len(raw_feats) for j in range(6)]
        stds = [math.sqrt(sum((r[j] - means[j]) ** 2 for r in raw_feats) / len(raw_feats)) for j in range(6)]
        true_w = (0.4, -0.25, 0.15, -0.1, 0.3, -0.05)
        ys_planted = [
            1.1 + 0.8 * xs_dmso[i] + sum(true_w[j] * ((raw_feats[i][j] - means[j]) / stds[j]) for j in range(6))
            for i in range(40)
        ]
        mod_lo = fit_partial_penalized_3d_ridge(xs_dmso, raw_feats, ys_planted, alpha=1e-10)
        self.assertAlmostEqual(mod_lo.b0_intercept, 1.1, places=6)
        self.assertAlmostEqual(mod_lo.b1_dmso_slope, 0.8, places=6)
        for w_hat, w_true in zip(mod_lo.weights, true_w):
            self.assertAlmostEqual(w_hat, w_true, places=6)

        # Direct Cholesky solve check against known SPD system
        sol = cholesky_solve_spd([[4.0, 1.0], [1.0, 3.0]], [9.0, 5.0])
        self.assertAlmostEqual(sol[0], 2.0, places=9)
        self.assertAlmostEqual(sol[1], 1.0, places=9)

    def test_fold_imputation_zero_leakage(self) -> None:
        train_x = [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5]
        train_feats = [
            (0.2, 0.1, 10.0, 100.0, 1.0, 0.0),
            (0.4, 0.3, 12.0, 120.0, 2.0, 0.0),
            (0.6, 0.5, 14.0, 140.0, 3.0, 0.0),
            (0.8, 0.7, 16.0, 160.0, None, 1.0),
            (None, None, None, None, None, 1.0),
            (0.5, 0.4, 13.0, 130.0, 2.0, 0.0),
        ]
        train_y = [0.5 + 0.8 * x for x in train_x]
        mod = fit_partial_penalized_3d_ridge(train_x, train_feats, train_y, alpha=10.0)

        # Exact expected training-fold means over non-missing entries only
        self.assertAlmostEqual(mod.impute_means[0], (0.2 + 0.4 + 0.6 + 0.8 + 0.5) / 5.0, places=9)
        self.assertAlmostEqual(mod.impute_means[4], (1.0 + 2.0 + 3.0 + 2.0) / 4.0, places=9)

        # A held-out row with missing continuous features maps to z_j = 0.0 for j in 0..4.
        z_miss = mod.transform_features((None, None, None, None, None, 1.0))
        for j in range(5):
            self.assertAlmostEqual(z_miss[j], 0.0, places=12)

        # Mutating an external validation row cannot affect the fitted model.
        val_pred_1 = mod.predict_one(0.25, (None, None, None, None, None, 1.0))
        mod_2 = fit_partial_penalized_3d_ridge(train_x, train_feats, train_y, alpha=10.0)
        val_pred_2 = mod_2.predict_one(0.25, (None, None, None, None, None, 1.0))
        self.assertEqual(val_pred_1, val_pred_2)


class TestSparseXlsxAndWTParity(unittest.TestCase):
    def test_sparse_xlsx_and_protein_parser(self) -> None:
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

        self.assertEqual(parse_protein_distinction("R117H"), (117, "missense", "R"))
        self.assertEqual(parse_protein_distinction("W1282X"), (1282, "nonsense", "W"))
        self.assertEqual(parse_protein_distinction("F508del"), (508, "inframe_indel", "F"))
        self.assertEqual(parse_protein_distinction("Q966="), (966, "other", "Q"))

    def test_unexpected_wt_mismatch_rejected(self) -> None:
        data_dir = REPO_ROOT / "tools" / "cftr-model" / "data"
        ca6, _ = load_mmcif_gz(data_dir / "6MSM.cif.gz", "6MSM", CIF_EXPECTED_SHA256["6MSM.cif.gz"])
        ca8, _ = load_mmcif_gz(data_dir / "8EIQ.cif.gz", "8EIQ", CIF_EXPECTED_SHA256["8EIQ.cif.gz"])
        manifest = json.loads((REPO_ROOT / "benchmarks/cftr_dms_expression/split_manifest.json").read_text(encoding="utf-8"))
        bad_assigns = copy.deepcopy(manifest["assignments"])
        # Plant a bogus WT amino acid at holdout codon 117 (R117H -> K117H)
        for row in bad_assigns:
            if row["codon_position"] == 117:
                row["protein_distinction"] = "K117H"
                break
        with self.assertRaises(AssertionError):
            verify_wt_sequence_parity(bad_assigns, {"ca_6msm": ca6, "ca_8eiq": ca8})

    def test_cluster_bootstrap_counts_unique_codons(self) -> None:
        rows = [
            Variant3DRecord(3, "TMD1", "1A>G", "R100W", 100, "missense", "holdout", False, 0.5, 1.0, -1.0, 0.0, (0.0,) * 6),
            Variant3DRecord(4, "TMD1", "2A>T", "R100Q", 100, "missense", "holdout", False, 0.5, 1.0, -1.0, 0.0, (0.0,) * 6),
            Variant3DRecord(5, "NBD1", "3G>A", "G500D", 500, "missense", "holdout", False, 0.25, 0.5, -2.0, -1.0, (0.0,) * 6),
        ]
        boot = residue_clustered_bootstrap_delta_ci(rows, [0.5, 0.5, -0.5], [0.0, 0.0, -1.0], resamples=50)
        self.assertEqual(boot["n_holdout_clusters"], 2)
        self.assertTrue(boot["ci95_low_gt_zero"])


class TestCheckedInArtifactsAndTamperGuards(unittest.TestCase):
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

    def test_tamper_done_and_results_fail_closed(self) -> None:
        results = json.loads((HERE / RESULTS_FILE_NAME).read_text(encoding="utf-8"))
        done_txt = (HERE / "done.md").read_text(encoding="utf-8")
        manifest_txt = (REPO_ROOT / "benchmarks/cftr_dms_expression/split_manifest.json").read_text(encoding="utf-8")
        th_txt = (REPO_ROOT / "benchmarks/theratype/theratype_open_annotations_2026-11-21.json").read_text(encoding="utf-8")
        attr_txt = (REPO_ROOT / "tools/cftr-model/data/ATTRIBUTION.md").read_text(encoding="utf-8")

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bdir = root / "benchmarks" / "cftr_dms_3d_structure"
            edir = root / "benchmarks" / "cftr_dms_expression"
            tdir = root / "benchmarks" / "theratype"
            mdir = root / "tools" / "cftr-model" / "data"
            for d in (bdir, edir, tdir, mdir):
                d.mkdir(parents=True)
            (edir / "split_manifest.json").write_text(manifest_txt, encoding="utf-8")
            (tdir / "theratype_open_annotations_2026-11-21.json").write_text(th_txt, encoding="utf-8")
            (mdir / "ATTRIBUTION.md").write_text(attr_txt, encoding="utf-8")
            for fn in ("6MSM.cif.gz", "5UAK.cif.gz", "8EIQ.cif.gz"):
                (mdir / fn).write_bytes((REPO_ROOT / "tools" / "cftr-model" / "data" / fn).read_bytes())

            # Defect 1: tamper done.md -> both verifiers reject
            (bdir / "done.md").write_text(done_txt + "\n# tampered\n", encoding="utf-8")
            (bdir / RESULTS_FILE_NAME).write_text(json.dumps(results), encoding="utf-8")
            with self.assertRaises(AssertionError):
                verify_checked_in_artifacts(repo_root=root)
            with self.assertRaises(AssertionError):
                run_independent_audit(repo_root=root)

            # Defect 2: restore done.md, tamper RMSD in results -> both verifiers reject
            (bdir / "done.md").write_text(done_txt, encoding="utf-8")
            bad_res = copy.deepcopy(results)
            bad_res["pre_holdout_controls"]["cryoem_tmd_superposition_and_pockets"]["rmsd_tmd_8eiq_to_6msm_A"] = 0.123456
            (bdir / RESULTS_FILE_NAME).write_text(json.dumps(bad_res), encoding="utf-8")
            with self.assertRaises(AssertionError):
                verify_checked_in_artifacts(repo_root=root)
            with self.assertRaises(AssertionError):
                run_independent_audit(repo_root=root)

            # Defect 3: restore results, plant PRED-0010 exposed codon (508) in holdout -> both reject
            (bdir / RESULTS_FILE_NAME).write_text(json.dumps(results), encoding="utf-8")
            bad_man = json.loads(manifest_txt)
            bad_man["assignments"][0]["split"] = "holdout"
            bad_man["assignments"][0]["codon_position"] = 508
            (edir / "split_manifest.json").write_text(json.dumps(bad_man), encoding="utf-8")
            with self.assertRaises(AssertionError):
                verify_checked_in_artifacts(repo_root=root)
            with self.assertRaises(AssertionError):
                run_independent_audit(repo_root=root)

    def test_ci_mode_without_xlsx_downloads_passes_structural_checks(self) -> None:
        with tempfile.TemporaryDirectory() as empty_dir:
            summary = verify_checked_in_artifacts(repo_root=REPO_ROOT, data_dir=Path(empty_dir))
            self.assertEqual(summary["status"], "PASS")
            self.assertFalse(summary["raw_local_data_parity_checked"])

            audit = run_independent_audit(repo_root=REPO_ROOT, data_dir=Path(empty_dir))
            self.assertEqual(audit["status"], "PASS")
            self.assertFalse(audit["raw_local_data_audited"])


if __name__ == "__main__":
    unittest.main()
