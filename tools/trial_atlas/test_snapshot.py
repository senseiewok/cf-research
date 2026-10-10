"""Tests for snapshot.py: the manifest hash covers what is read. A page under another name, a path outside the snapshot, an unlisted
file of any kind, an edited manifest field and an invalid NCT id are all refused. SYNTHETIC fixtures, no network.
Usage: python -m unittest -v (inside tools/trial_atlas)."""
import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import snapshot as snap  # noqa: E402
import synthetic_fixtures as sf  # noqa: E402


def rewrite_manifest(root: Path, mutate):
    """Edit the manifest and recompute its digest the way the code under test does, as someone rewriting the manifest would."""
    mp = root / "manifest.json"
    m = json.loads(mp.read_text(encoding="utf-8"))
    mutate(m)
    if hasattr(snap, "manifest_digest"):
        m["snapshot_sha256"] = snap.manifest_digest(m)
    else:
        m["snapshot_sha256"] = snap.snapshot_digest(m["files"])
    mp.write_text(json.dumps(m), encoding="utf-8")


class SnapshotLoadTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory()
        cls.base = Path(cls._tmp.name) / "base"
        sf.build_cf_snapshot(cls.base)

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def setUp(self):
        sf.block_network(self)
        self.work = Path(tempfile.mkdtemp(dir=self._tmp.name))
        self.root = self.work / "snap"
        shutil.copytree(self.base, self.root)

    def test_the_untouched_snapshot_loads(self):
        self.assertEqual(len(snap.load(self.root).studies("condition")), 16)

    def test_a_page_under_another_name_is_refused(self):
        src = self.root / "condition" / "page-0001.json"
        evil = json.loads(src.read_text(encoding="utf-8"))
        evil["studies"] = evil["studies"][:1]
        (self.root / "condition" / "page-0001.txt").write_text(json.dumps(evil), encoding="utf-8")
        rewrite_manifest(self.root, lambda m: m["routes"][0]["pages"][0].update(file="condition/page-0001.txt"))
        with self.assertRaises(snap.SnapshotError):
            snap.load(self.root)

    def test_a_route_page_outside_the_snapshot_is_refused(self):
        shutil.copy(self.root / "condition" / "page-0001.json", self.work / "outside.json")
        rewrite_manifest(self.root, lambda m: m["routes"][0]["pages"][0].update(file="../outside.json"))
        with self.assertRaises(snap.SnapshotError):
            snap.load(self.root)

    def test_an_absolute_route_page_is_refused(self):
        target = (self.root / "condition" / "page-0001.json").resolve().as_posix()
        rewrite_manifest(self.root, lambda m: m["routes"][0]["pages"][0].update(file=target))
        with self.assertRaises(snap.SnapshotError):
            snap.load(self.root)

    def test_any_unlisted_file_is_refused(self):
        (self.root / "notes.txt").write_text("SYNTHETIC stray file", encoding="utf-8")
        with self.assertRaises(snap.SnapshotError) as cm:
            snap.load(self.root)
        self.assertIn("notes.txt", str(cm.exception))

    def test_the_digest_covers_the_whole_manifest(self):
        mp = self.root / "manifest.json"
        m = json.loads(mp.read_text(encoding="utf-8"))
        m["data_timestamp"] = "2001-01-01T00:00:00"          # edited without recomputing anything
        m["routes"][0]["total_count"] = 99
        mp.write_text(json.dumps(m), encoding="utf-8")
        with self.assertRaises(snap.SnapshotError):
            snap.load(self.root)

    def test_an_invalid_nct_id_is_refused_on_load(self):
        page = self.root / "condition" / "page-0001.json"
        data = json.loads(page.read_text(encoding="utf-8"))
        data["studies"][0]["protocolSection"]["identificationModule"]["nctId"] = "NCT1\x1b[31m"
        page.write_text(json.dumps(data), encoding="utf-8")
        rewrite_manifest(self.root, lambda m: m["files"].update({"condition/page-0001.json": snap.sha256_file(page)}))
        with self.assertRaises(snap.SnapshotError):
            snap.load(self.root)

    def test_study_record_refuses_an_invalid_nct_id(self):
        with self.assertRaises(snap.SnapshotError):
            snap.study_record({"protocolSection": {"identificationModule": {"nctId": "not-an-id"}}})


if __name__ == "__main__":
    unittest.main()
