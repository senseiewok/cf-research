"""Unit tests for tools/trial_atlas/consensus_audit.py."""
from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
import sys
sys.path.insert(0, str(HERE))
import consensus_audit as ca
import snapshot as snap


class TestConsensusAudit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ca.ensure_artifacts()

    def test_sampling_stratification(self):
        data = json.loads(ca.DEFAULT_SNAPSHOT_JSON.read_text(encoding="utf-8"))
        sampled = ca.sample_stratified_studies(data, seed=ca.DEFAULT_SEED, count=ca.DEFAULT_CENSUS_SIZE)
        self.assertEqual(len(sampled), 50)

        # Check phase representation
        phases = {x["phase"] for x in sampled}
        for p in ("PHASE1", "PHASE2", "PHASE3", "PHASE4"):
            self.assertIn(p, phases)

        # Check all 8 canonical classes
        classes = {c for x in sampled for c in x["classes"]}
        for c in ca.CANONICAL_8_CLASSES:
            self.assertIn(c, classes)

    def test_benchmark_studies_included(self):
        data = json.loads(ca.DEFAULT_SNAPSHOT_JSON.read_text(encoding="utf-8"))
        sampled = ca.sample_stratified_studies(data, seed=ca.DEFAULT_SEED, count=ca.DEFAULT_CENSUS_SIZE)
        sampled_ncts = {x["study"]["nct_id"] for x in sampled}
        # Check presence of key benchmark trials in snapshot
        self.assertIn("NCT03525444", sampled_ncts)
        self.assertIn("NCT05033080", sampled_ncts)
        self.assertIn("NCT05668741", sampled_ncts)

    def test_quote_fidelity_invariant(self):
        data = json.loads(ca.DEFAULT_SNAPSHOT_JSON.read_text(encoding="utf-8"))
        sampled = ca.sample_stratified_studies(data, seed=ca.DEFAULT_SEED, count=ca.DEFAULT_CENSUS_SIZE)
        summary, rows = ca.evaluate_dual_agent_census(sampled)

        self.assertEqual(summary["quote_fidelity"], 1.00)
        for r in rows:
            for t in r["consensus_tags"]:
                field = t["field"]
                raw_text = r[field]
                self.assertEqual(raw_text[t["start"] : t["end"]], t["quote"])

    def test_precision_recall_thresholds(self):
        data = json.loads(ca.DEFAULT_SNAPSHOT_JSON.read_text(encoding="utf-8"))
        sampled = ca.sample_stratified_studies(data, seed=ca.DEFAULT_SEED, count=ca.DEFAULT_CENSUS_SIZE)
        summary, _ = ca.evaluate_dual_agent_census(sampled)

        self.assertGreaterEqual(summary["precision"], ca.MIN_PRECISION)
        self.assertGreaterEqual(summary["recall"], ca.MIN_RECALL)
        self.assertGreaterEqual(summary["f1"], 0.85)

    def test_snapshot_loadable_by_snapshot_py(self):
        s = snap.load(ca.DEFAULT_SNAPSHOT_DIR)
        self.assertEqual(len(s.studies("condition")), 1053)
        self.assertTrue(s.manifest.get("config_verified"))

    def test_key_json_integrity(self):
        key_data = json.loads(ca.DEFAULT_KEY_PATH.read_text(encoding="utf-8"))
        self.assertEqual(key_data.get("design"), "consensus-v1")
        self.assertEqual(key_data.get("study_count"), 50)
        self.assertIn("key_sha256", key_data)

    def test_run_audit_success_and_failure(self):
        s = snap.load(ca.DEFAULT_SNAPSHOT_DIR)
        tags_data = json.loads(ca.DEFAULT_TAGS_PATH.read_text(encoding="utf-8"))
        code = ca.run_audit(s, tags_data, ca.DEFAULT_KEY_PATH)
        self.assertEqual(code, 0)

        with tempfile.TemporaryDirectory() as td:
            bad_key = Path(td) / "missing_key.json"
            code_fail = ca.run_audit(s, tags_data, bad_key)
            self.assertEqual(code_fail, 1)

    def test_run_audit_fails_on_corrupted_quote(self):
        s = snap.load(ca.DEFAULT_SNAPSHOT_DIR)
        tags_data = json.loads(ca.DEFAULT_TAGS_PATH.read_text(encoding="utf-8"))
        key_data = json.loads(ca.DEFAULT_KEY_PATH.read_text(encoding="utf-8"))

        # Corrupt a verbatim quote in consensus_tags for the first row that has tags
        modified = False
        for r in key_data["rows"]:
            if r.get("consensus_tags"):
                r["consensus_tags"][0]["quote"] = "NON_VERBATIM_CORRUPTED_TEXT"
                modified = True
                break
        self.assertTrue(modified)

        # Recompute valid key_sha256 so the failure is specifically due to quote fidelity, not key tampering
        body = {k: v for k, v in key_data.items() if k != "key_sha256"}
        key_data["key_sha256"] = hashlib.sha256(
            (json.dumps(body, indent=1, sort_keys=True, ensure_ascii=False) + "\n").encode("utf-8")
        ).hexdigest()

        with tempfile.TemporaryDirectory() as td:
            corrupt_key = Path(td) / "corrupt_quote_key.json"
            corrupt_key.write_text(json.dumps(key_data, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
            code = ca.run_audit(s, tags_data, corrupt_key)
            self.assertEqual(code, 1)

    def test_run_audit_fails_on_corrupted_classes(self):
        s = snap.load(ca.DEFAULT_SNAPSHOT_DIR)
        tags_data = json.loads(ca.DEFAULT_TAGS_PATH.read_text(encoding="utf-8"))
        key_data = json.loads(ca.DEFAULT_KEY_PATH.read_text(encoding="utf-8"))

        # Corrupt classes to an invalid class name
        key_data["rows"][0]["ground_truth_classes"] = ["invalid_unrecognized_class"]

        # Recompute key_sha256 so the failure is specifically due to invalid class
        body = {k: v for k, v in key_data.items() if k != "key_sha256"}
        key_data["key_sha256"] = hashlib.sha256(
            (json.dumps(body, indent=1, sort_keys=True, ensure_ascii=False) + "\n").encode("utf-8")
        ).hexdigest()

        with tempfile.TemporaryDirectory() as td:
            corrupt_key = Path(td) / "corrupt_classes_key.json"
            corrupt_key.write_text(json.dumps(key_data, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
            code = ca.run_audit(s, tags_data, corrupt_key)
            self.assertEqual(code, 1)

    def test_run_audit_fails_on_tampered_key_sha256(self):
        s = snap.load(ca.DEFAULT_SNAPSHOT_DIR)
        tags_data = json.loads(ca.DEFAULT_TAGS_PATH.read_text(encoding="utf-8"))
        key_data = json.loads(ca.DEFAULT_KEY_PATH.read_text(encoding="utf-8"))

        # Tamper with summary or rows without recomputing sha256
        key_data["summary"]["design"] = "tampered-design-v2"

        with tempfile.TemporaryDirectory() as td:
            tampered_key = Path(td) / "tampered_key.json"
            tampered_key.write_text(json.dumps(key_data, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
            code = ca.run_audit(s, tags_data, tampered_key)
            self.assertEqual(code, 1)


if __name__ == "__main__":
    unittest.main()
