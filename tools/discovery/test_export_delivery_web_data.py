#!/usr/bin/env python3
"""Unit tests for export_delivery_web_data.py (T-0180)."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import export_delivery_web_data as exporter
from transcytosis_model import TranscytosisModel
from epithelial_plasticity import AirwayEpithelialPlasticitySimulator
from cross_domain_miner import CrossDomainInnovationMiner


class TestExportDeliveryWebData(unittest.TestCase):
    def test_models_determinism_and_content(self):
        # 1. TranscytosisModel
        tm = TranscytosisModel()
        flux1 = tm.simulate_flux()
        flux2 = tm.simulate_flux()
        self.assertEqual(flux1, flux2)
        self.assertGreaterEqual(len(flux1), 3 * 8)
        receptors = {r["receptor_id"] for r in flux1}
        self.assertEqual(receptors, {"fcrn", "tfr", "pigr"})
        for row in flux1:
            self.assertIn("concentration_nm", row)
            self.assertIn("surface_binding_pct", row)
            self.assertIn("internalized_pct", row)
            self.assertIn("basolateral_exocytosis_pct", row)
            self.assertIn("net_flux_pct_per_h", row)
            self.assertIn("verdict", row)

        # 2. AirwayEpithelialPlasticitySimulator
        sim = AirwayEpithelialPlasticitySimulator()
        traj1 = sim.simulate_180d_trajectory()
        traj2 = sim.simulate_180d_trajectory()
        self.assertEqual(traj1, traj2)
        # Days 0 to 180 inclusive -> 181 points
        self.assertEqual(len(traj1), 181)
        self.assertEqual(traj1[0]["day"], 0)
        self.assertEqual(traj1[-1]["day"], 180)
        for pt in traj1:
            self.assertIn("uncorrected_secretory", pt)
            self.assertIn("corrected_secretory", pt)
            self.assertIn("transient_intermediate", pt)
            self.assertIn("corrected_basal_stem_cells", pt)
            self.assertIn("ciliated", pt)
            self.assertIn("cftr_plus_stem_cell_pct", pt)

        # Day 0 should start with 0% CFTR+ basal stem cells
        self.assertEqual(traj1[0]["cftr_plus_stem_cell_pct"], 0.0)
        # Transient intermediate should peak during pulse and resolve back to 0.0
        intermediate_values = [pt["transient_intermediate"] for pt in traj1]
        self.assertGreater(max(intermediate_values), 0.0)
        self.assertEqual(traj1[-1]["transient_intermediate"], 0.0)
        # Day 180 should retain sustained CFTR+ basal stem cell chimerism
        self.assertGreaterEqual(traj1[-1]["cftr_plus_stem_cell_pct"], 15.0)

        # 3. CrossDomainInnovationMiner
        miner = CrossDomainInnovationMiner()
        rank1 = miner.rank_candidates()
        rank2 = miner.rank_candidates()
        self.assertEqual(rank1, rank2)
        self.assertGreaterEqual(len(rank1), 5)
        for item in rank1:
            self.assertIn("transfer_domain", item)
            self.assertIn("bfi_score", item)
            self.assertIn("physical_barriers", item)
            self.assertIn("mechanism", item)
            self.assertTrue(item["primary_doi"])

    def test_json_ld_schema(self):
        json_ld = exporter.JSON_LD
        self.assertEqual(json_ld["@context"], "https://schema.org")
        self.assertEqual(json_ld["@type"], "Dataset")
        self.assertTrue(json_ld["name"])
        self.assertTrue(json_ld["description"])
        self.assertTrue(json_ld["license"].startswith("http"))
        self.assertIsInstance(json_ld["citation"], list)
        self.assertGreater(len(json_ld["citation"]), 0)
        for cite in json_ld["citation"]:
            self.assertTrue(cite.startswith("https://doi.org/10."))

    def test_metadata_tier_and_schema(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            snapshot, _, _, _ = exporter.build_snapshot_and_tsv(tmp_path)
            meta = snapshot["metadata"]
            self.assertEqual(meta["tier"], "Tier T0/T1 Exploratory In Silico Hypothesis (Unverified Parametric Constants)")
            self.assertEqual(meta["schema_version"], "1.0.0")
            self.assertEqual(
                meta["terms_of_use"],
                "Exploratory in silico hypothesis / unverified parametric constants only. Not primary empirical evidence and not medical advice."
            )
            self.assertEqual(len(meta["source_code_hash"]), 64)
            self.assertIsInstance(meta["doi_citations"], list)
            self.assertGreaterEqual(len(meta["doi_citations"]), 10)
            for doi in meta["doi_citations"]:
                self.assertTrue(doi.startswith("10."))

    def test_sha256_checksums(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = Path(tmp_dir)
            snapshot, bfi_tsv, transcytosis_tsv, trajectory_tsv = exporter.build_snapshot_and_tsv(tmp_path)
            checksums = snapshot["checksums"]

            # Validate TSV sha256
            self.assertEqual(checksums["bfi_rankings_tsv"], hashlib.sha256(bfi_tsv.encode("utf-8")).hexdigest())
            self.assertEqual(checksums["transcytosis_kinetics_tsv"], hashlib.sha256(transcytosis_tsv.encode("utf-8")).hexdigest())
            self.assertEqual(checksums["epithelial_trajectory_tsv"], hashlib.sha256(trajectory_tsv.encode("utf-8")).hexdigest())

            # Validate datasets canonical data_block sha256
            data_block_json = json.dumps(snapshot["datasets"], sort_keys=True, separators=(",", ":"))
            expected_block_sha = hashlib.sha256(data_block_json.encode("utf-8")).hexdigest()
            self.assertEqual(checksums["data_block"], expected_block_sha)
            self.assertEqual(checksums["datasets_block"], expected_block_sha)

    def test_cli_export_and_check(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            script_path = HERE / "export_delivery_web_data.py"

            # 1. --check on empty directory should fail (exit code 1)
            res_missing = subprocess.run(
                [sys.executable, str(script_path), "--out-dir", tmp_dir, "--check"],
                capture_output=True,
                text=True,
            )
            self.assertEqual(res_missing.returncode, 1)
            self.assertIn("Missing file", res_missing.stdout)

            # 2. Write artifacts
            res_write = subprocess.run(
                [sys.executable, str(script_path), "--out-dir", tmp_dir],
                capture_output=True,
                text=True,
            )
            self.assertEqual(res_write.returncode, 0)
            self.assertIn("Export completed successfully", res_write.stdout)

            # Confirm all 4 files exist
            snap_file = Path(tmp_dir) / "delivery_discovery_snapshot.json"
            bfi_file = Path(tmp_dir) / "bfi_rankings.tsv"
            trans_file = Path(tmp_dir) / "transcytosis_kinetics.tsv"
            traj_file = Path(tmp_dir) / "epithelial_trajectory.tsv"
            for f in (snap_file, bfi_file, trans_file, traj_file):
                self.assertTrue(f.exists(), f"File {f.name} was not created")

            # 3. --check should now pass (exit code 0)
            res_check = subprocess.run(
                [sys.executable, str(script_path), "--out-dir", tmp_dir, "--check"],
                capture_output=True,
                text=True,
            )
            self.assertEqual(res_check.returncode, 0)
            self.assertIn("up to date and verified", res_check.stdout)

            # 4. Modify one TSV to simulate drift
            with open(bfi_file, "a", encoding="utf-8") as f:
                f.write("DRIFT_LINE\tcorrupted\n")

            res_drift = subprocess.run(
                [sys.executable, str(script_path), "--out-dir", tmp_dir, "--check"],
                capture_output=True,
                text=True,
            )
            self.assertEqual(res_drift.returncode, 1)
            self.assertIn("Drift detected in bfi_rankings.tsv", res_drift.stdout)

            # 5. Overwrite with valid export again
            res_restore = subprocess.run(
                [sys.executable, str(script_path), "--out-dir", tmp_dir],
                capture_output=True,
                text=True,
            )
            self.assertEqual(res_restore.returncode, 0)

            # 6. Delete a file to simulate deletion
            traj_file.unlink()
            res_del = subprocess.run(
                [sys.executable, str(script_path), "--out-dir", tmp_dir, "--check"],
                capture_output=True,
                text=True,
            )
            self.assertEqual(res_del.returncode, 1)
            self.assertIn("Missing file: epithelial_trajectory.tsv", res_del.stdout)


if __name__ == "__main__":
    unittest.main()
