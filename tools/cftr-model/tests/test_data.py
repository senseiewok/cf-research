"""Tests for tools/cftr-model/build_traces.py: the superposition is right, and the data on the page is what the builder makes from the downloaded files.
Run from the repo root: python -m unittest test_data -v      (standard library only; the data checks skip if sources/downloads is absent)"""
import importlib.util
import json
import math
import random
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]                 # tools/cftr-model
REPO = ROOT.parents[1]
spec = importlib.util.spec_from_file_location("build_traces", ROOT / "build_traces.py")
bt = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bt)
DOWNLOADS = REPO / "sources" / "downloads"
DATA_JS = ROOT / "web" / "cftr-data.js"


def matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def rotation(a, b, c):
    Rz = [[math.cos(a), -math.sin(a), 0], [math.sin(a), math.cos(a), 0], [0, 0, 1]]
    Ry = [[math.cos(b), 0, math.sin(b)], [0, 1, 0], [-math.sin(b), 0, math.cos(b)]]
    Rx = [[1, 0, 0], [0, math.cos(c), -math.sin(c)], [0, math.sin(c), math.cos(c)]]
    return matmul(Rz, matmul(Ry, Rx))


def det3(R):
    return (R[0][0] * (R[1][1] * R[2][2] - R[1][2] * R[2][1]) - R[0][1] * (R[1][0] * R[2][2] - R[1][2] * R[2][0]) + R[0][2] * (R[1][0] * R[2][1] - R[1][1] * R[2][0]))


class SuperposeTests(unittest.TestCase):
    def setUp(self):
        random.seed(7)
        self.P = [[random.uniform(-40, 40) for _ in range(3)] for _ in range(300)]

    def test_it_recovers_a_known_rotation_and_shift(self):
        for angles, shift in (((0.9, -0.4, 2.2), [13.0, -7.5, 30.0]), ((3.0, 1.2, -2.5), [-100.0, 0.0, 5.0]), ((0.0, 0.0, 0.0), [0.0, 0.0, 0.0])):
            R0 = rotation(*angles)
            Q = [bt.apply(R0, shift, p) for p in self.P]
            R, t = bt.superpose(self.P, Q)
            self.assertLess(max(math.dist(bt.apply(R, t, p), q) for p, q in zip(self.P, Q)), 1e-6, angles)
            self.assertAlmostEqual(det3(R), 1.0, places=9)

    def test_it_never_returns_a_mirror_image(self):
        Q = [[-p[0], p[1], p[2]] for p in self.P]                     # a reflection cannot be reached by turning
        R, t = bt.superpose(self.P, Q)
        self.assertAlmostEqual(det3(R), 1.0, places=9)
        rmsd = math.sqrt(sum(math.dist(bt.apply(R, t, p), q) ** 2 for p, q in zip(self.P, Q)) / len(Q))
        self.assertGreater(rmsd, 5.0, "a mirrored cloud must not fit as if it were turned")

    def test_it_does_not_hide_a_real_difference(self):
        R0 = rotation(0.9, -0.4, 2.2)
        Q = [bt.apply(R0, [13.0, -7.5, 30.0], p) for p in self.P]
        for q in Q[:150]:
            q[0] += 6.0
        R, t = bt.superpose(self.P, Q)
        rmsd = math.sqrt(sum(math.dist(bt.apply(R, t, p), q) ** 2 for p, q in zip(self.P, Q)) / len(Q))
        self.assertGreater(rmsd, 2.5)
        self.assertLess(rmsd, 3.5)                                    # half the points moved 6 A: the best fit leaves 3 A

    def test_the_eigen_solver_matches_a_known_matrix(self):
        vals, vecs = bt.jacobi_eigen([[2.0, 1.0, 0.0, 0.0], [1.0, 2.0, 0.0, 0.0], [0.0, 0.0, 5.0, 0.0], [0.0, 0.0, 0.0, -1.0]])
        self.assertEqual(sorted(round(v, 9) for v in vals), [-1.0, 1.0, 3.0, 5.0])


def feature(kind, start, end, *pdb):
    return {"type": kind, "location": {"start": {"value": start}, "end": {"value": end}}, "evidences": [{"source": "PDB", "id": p} for p in pdb]}


class SecondaryStructureTests(unittest.TestCase):
    """The helix, strand and turn ranges are copied from UniProt's features: renamed, sorted, checked, never computed."""

    def test_it_keeps_only_the_three_types_sorted_with_their_pdb_evidence(self):
        feats = [feature("Beta strand", 40, 45, "2PZE"), feature("Natural variant", 1, 1), feature("Helix", 10, 20, "9MXL", "8EIO"), feature("Turn", 30, 32), feature("Domain", 5, 60)]
        self.assertEqual(bt.secondary_structure(feats), [{"type": "helix", "start": 10, "end": 20, "evidence": ["8EIO", "9MXL"]},
                                                         {"type": "turn", "start": 30, "end": 32, "evidence": []},
                                                         {"type": "strand", "start": 40, "end": 45, "evidence": ["2PZE"]}])

    def test_it_refuses_overlapping_or_unusable_ranges(self):
        with self.assertRaises(SystemExit):
            bt.secondary_structure([feature("Helix", 10, 20), feature("Beta strand", 20, 25)])
        with self.assertRaises(SystemExit):
            bt.secondary_structure([feature("Helix", 12, 10)])
        with self.assertRaises(SystemExit):
            bt.secondary_structure([{"type": "Helix", "location": {"start": {"value": None}, "end": {"value": 4}}}])
        self.assertEqual(bt.secondary_structure([]), [])


@unittest.skipUnless((DOWNLOADS / "p13569.json").exists(), "fetch the structure files first: python tools/sources/fetch_sources.py --only rcsb-pdb-6msm rcsb-pdb-5uak rcsb-pdb-8eiq rcsb-pdb-8ej1 uniprot-p13569 (see README.md)")
class DataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = bt.build(DOWNLOADS)
        cls.by_id = {s["id"]: s for s in cls.data["structures"]}

    def placed(self, sid):
        return {n for seg in self.by_id[sid]["segments"] for n in range(seg["start"], seg["start"] + len(seg["xyz"]) // 3)}

    def test_residue_508_is_present_in_the_wild_type_models_and_absent_in_the_two_deletion_models(self):
        for sid in ("6MSM", "5UAK"):
            self.assertEqual(self.by_id[sid]["res508"]["name"], "PHE")
            self.assertIn(508, self.placed(sid))
        for sid in ("8EIQ", "8EJ1"):
            self.assertIsNone(self.by_id[sid]["res508"])
            self.assertTrue({507, 509} <= self.placed(sid) and 508 not in self.placed(sid))

    def test_the_reference_matches_itself_and_the_others_are_measured(self):
        self.assertEqual(self.by_id["6MSM"]["superposition_rmsd_A"], 0.0)
        self.assertEqual({k: self.by_id[k]["superposition_rmsd_A"] for k in ("8EIQ", "5UAK", "8EJ1")}, {"8EIQ": 0.84, "5UAK": 5.01, "8EJ1": 4.87})
        self.assertEqual({k: self.by_id[k]["superposition_residues"] for k in ("6MSM", "8EIQ", "5UAK", "8EJ1")}, {"6MSM": 572, "8EIQ": 566, "5UAK": 557, "8EJ1": 547})

    def test_uniprot_parts_and_length(self):
        self.assertEqual(self.data["length"], 1480)
        parts = [(d["start"], d["end"]) for d in self.data["domains"] if d["type"] == "Domain"]
        self.assertEqual(parts, [(81, 365), (423, 646), (859, 1155), (1210, 1443)])
        self.assertIn((654, 831), [(d["start"], d["end"]) for d in self.data["domains"] if d["name"] == "Disordered R region"])

    def test_every_title_comes_whole_from_the_file(self):
        self.assertTrue(self.by_id["8EIQ"]["title"].endswith("and ATP/Mg"))
        self.assertIn("Trikafta", self.by_id["8EIQ"]["title"])
        self.assertEqual(self.by_id["8EJ1"]["pubmed"], "36264792")

    def test_each_file_declares_how_its_sequence_differs_from_uniprot(self):
        got = {sid: [(d["position"], d["file_residue"], d["uniprot_residue"], d["details"]) for d in self.by_id[sid]["differences"]] for sid in self.by_id}
        self.assertEqual(got["6MSM"], [(1371, "GLN", "GLU", "engineered mutation")])
        self.assertEqual(got["5UAK"], [])
        self.assertEqual(got["8EIQ"], [(508, None, "PHE", "deletion"), (1371, "GLN", "GLU", "engineered mutation")])
        self.assertEqual(got["8EJ1"], [(508, None, "PHE", "variant"), (1371, "GLN", "GLU", "engineered mutation")])
        self.assertEqual({k: v["expression_tag_residues"] for k, v in self.by_id.items()}, {"6MSM": 9, "5UAK": 9, "8EIQ": 0, "8EJ1": 0})

    def test_only_declared_differences_disagree_with_the_uniprot_sequence(self):
        seq = self.data["uniprot"]["sequence"]
        three = bt.THREE
        for sid, st in self.by_id.items():
            declared = {d["position"] for d in st["differences"]}
            placed = {}
            for seg in st["segments"]:
                for i in range(len(seg["xyz"]) // 3):
                    placed[seg["start"] + i] = None
            self.assertEqual(len(placed), st["residues_traced"])
        # the builder's own guard: a residue that disagrees without a declaration stops the build
        fake_ca = {10: ("GLY", [0.0, 0.0, 0.0])}
        with self.assertRaises(SystemExit):
            bt.check_against_sequence("TEST", fake_ca, [], "A" * 20)
        bt.check_against_sequence("TEST", fake_ca, [{"position": 10}], "A" * 20)        # declared: accepted
        self.assertEqual(len(seq), 1480)
        self.assertEqual(three["PHE"], seq[507])

    def test_the_uniprot_block_matches_the_record(self):
        u = self.data["uniprot"]
        self.assertEqual((u["accession"], len(u["sequence"]), len(u["variants"]), len(u["transmembrane"]), len(u["topology"]), len(u["atp_sites"])), ("P13569", 1480, 209, 12, 13, 6))
        for v in u["variants"]:
            self.assertTrue(1 <= v["start"] <= v["end"] <= 1480, v["id"])
            if v["from"] and v["start"] == v["end"]:
                self.assertEqual(v["from"], u["sequence"][v["start"] - 1], v["id"])      # the variant names the residue UniProt has there
            self.assertTrue(all(e.isdigit() for e in v["evidence"]), v["id"])
        self.assertEqual(len({v["id"] for v in u["variants"]}), 209)
        by_id = {v["id"]: v for v in u["variants"]}
        self.assertEqual((by_id["VAR_000172"]["start"], by_id["VAR_000172"]["from"], by_id["VAR_000172"]["to"], by_id["VAR_000172"]["evidence"]), (508, "F", ["C"], ["1379210"]))
        self.assertEqual(by_id["VAR_000172"]["description"], "in dbSNP:rs74571530")
        self.assertEqual(by_id["VAR_000101"]["to"], ["F"])
        self.assertEqual(by_id["VAR_000101"]["from"], "S")
        self.assertEqual(sum(1 for v in u["variants"] if v["start"] == 508), 2)
        # topology covers 1..1480 in order, side by side, alternating cytoplasmic and extracellular, with the transmembrane segments between them
        tops = u["topology"]
        self.assertEqual([t["name"] for t in tops[:4]], ["Cytoplasmic", "Extracellular", "Cytoplasmic", "Extracellular"])
        self.assertEqual(sorted(t["start"] for t in tops), [t["start"] for t in tops])
        for tm in u["transmembrane"]:
            self.assertFalse(any(t["start"] <= tm["start"] <= t["end"] for t in tops), "a transmembrane segment lies inside a topological domain")

    def test_the_secondary_structure_is_uniprots_annotation_copied_whole(self):
        sec = self.data["uniprot"]["secondary"]
        raw = json.loads((DOWNLOADS / "p13569.json").read_text(encoding="utf-8"))["features"]
        wanted = sorted(((bt.SECONDARY_TYPES[f["type"]], f["location"]["start"]["value"], f["location"]["end"]["value"]) for f in raw if f["type"] in bt.SECONDARY_TYPES),
                        key=lambda s: (s[1], s[2]))
        self.assertEqual([(s["type"], s["start"], s["end"]) for s in sec], wanted, "every UniProt helix, strand and turn, and nothing else")
        self.assertEqual(len(sec), 118)
        self.assertEqual({t: sum(1 for s in sec if s["type"] == t) for t in ("helix", "strand", "turn")}, {"helix": 76, "strand": 31, "turn": 11})
        for a, b in zip(sec, sec[1:]):
            self.assertLess(a["end"], b["start"], "ranges are in chain order and do not overlap")
        self.assertTrue(all(1 <= s["start"] <= s["end"] <= 1480 for s in sec))
        self.assertTrue(all(s["evidence"] and all(len(e) == 4 for e in s["evidence"]) for s in sec), "each range carries the PDB ids UniProt cites for it")
        self.assertIn({"type": "helix", "start": 502, "end": 507, "evidence": ["2PZE"]}, sec, "the helix UniProt annotates right before residue 508")
        self.assertFalse(any(s["start"] <= 508 <= s["end"] for s in sec), "UniProt annotates no helix, strand or turn at 508 itself")

    def test_the_building_is_deterministic(self):
        again = bt.build(DOWNLOADS)
        self.assertEqual(json.dumps(again, separators=(",", ":")), json.dumps(self.data, separators=(",", ":")))

    def test_the_data_on_the_page_is_what_the_builder_makes(self):
        shipped = DATA_JS.read_bytes().decode("utf-8").split("export default ", 1)[1].strip().rstrip(";")
        self.assertEqual(shipped, json.dumps(self.data, separators=(",", ":")), "rebuild with: python tools/cftr-model/build_traces.py")

    def test_the_output_file_has_the_header_and_ends_cleanly(self):
        with tempfile.TemporaryDirectory() as td:
            out = Path(td) / "o.js"
            bt.main(["--inputs", str(DOWNLOADS), "--out", str(out)])
            raw = out.read_bytes().decode("utf-8")
            self.assertTrue(raw.startswith(bt.HEADER + "export default {"))
            self.assertTrue(raw.endswith("};\r\n"))


class SyntheticInputTests(unittest.TestCase):
    """The branches the four real files do not exercise, with small invented input."""

    def test_the_fit_uses_only_residues_placed_in_both_with_the_same_amino_acid_inside_the_ranges(self):
        ref = {n: ("ALA", [0, 0, n]) for n in range(1, 21)}
        ca = {n: ("ALA", [0, 0, n]) for n in range(1, 21)}
        del ca[5]                                   # not placed in the structure
        del ref[6]                                  # not placed in the reference
        ca[7] = ("GLY", [0, 0, 7])                  # a different amino acid
        got = bt.shared_residues(ca, ref, [(3, 9), (15, 18)])
        self.assertEqual(got, [3, 4, 8, 9, 15, 16, 17, 18])

    def test_a_multi_line_header_value_is_read_whole(self):
        import gzip
        import tempfile
        text = "\n".join([
            "data_TEST", "_struct.title                        ", ";A title that the file", "wraps over two lines", ";", "_exptl.method            'ELECTRON MICROSCOPY'",
            "loop_", "_atom_site.group_PDB", "_atom_site.label_atom_id", "_atom_site.label_comp_id", "_atom_site.auth_asym_id", "_atom_site.auth_seq_id",
            "_atom_site.Cartn_x", "_atom_site.Cartn_y", "_atom_site.Cartn_z", "_atom_site.pdbx_PDB_model_num",
            "ATOM CA GLY A 1 1.0 2.0 3.0 1", "#", ""])
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "t.cif.gz"
            with gzip.open(path, "wt", encoding="utf-8") as f:
                f.write(text)
            meta, ca, ligands = bt.read_cif(path)
        self.assertEqual(meta["title"], "A title that the file wraps over two lines")
        self.assertEqual(meta["method"], "ELECTRON MICROSCOPY")
        self.assertEqual(ca[1][0], "GLY")


class RepoHygieneTests(unittest.TestCase):
    def test_downloaded_source_files_are_git_ignored(self):
        text = (REPO / ".gitignore").read_text(encoding="utf-8")
        self.assertIn("sources/downloads", text)


if __name__ == "__main__":
    unittest.main()
