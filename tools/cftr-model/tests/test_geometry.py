"""Tests of the viewer's pure 3D geometry, called inside a real browser with synthetic and real input (the functions are exported as `internals`).

The page promises three things about the picture, and a picture that looks plausible can still break them: the tube passes through every measured point (between two points
its path is a drawing choice, but at them it is not), a helix or a strand is drawn wider than a coil, and the mesh stays inside its budget. Each is checked here against a value
computed independently in Python.

Run from this folder:   python -m unittest test_geometry -v        (needs Playwright for Python and a Chromium)"""
import math
import random
import threading
import unittest

import support as V
from support import sync_playwright

PAGE_URL = "cftr-viewer.js"
LOAD = "() => import('./cftr-viewer.js').then(m => m.internals)"


def helix_run(n=14):
    """A synthetic run of consecutive C-alpha-like points: a gentle helix, 3.8 apart or so."""
    pts = []
    for i in range(n):
        a = i * 1.75
        pts += [2.3 * math.cos(a), 2.3 * math.sin(a), 1.5 * i]
    return pts


def max_turn_degrees(pos):
    """The largest change of direction between two neighbouring pieces of a curve given as a flat list of ring centres."""
    dirs = []
    for r in range(len(pos) // 3 - 1):
        d = [pos[3 * r + 3 + k] - pos[3 * r + k] for k in range(3)]
        l = math.hypot(*d)
        dirs.append([x / l for x in d])
    return max(math.degrees(math.acos(max(-1, min(1, sum(a * b for a, b in zip(dirs[r], dirs[r + 1])))))) for r in range(len(dirs) - 1))


@unittest.skipIf(sync_playwright is None, "Playwright is not installed in this Python")
class GeometryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        chrome = V.find_chrome()
        if not chrome:
            raise unittest.SkipTest("no Chromium found")
        cls.httpd = V.model_server.make_server(V.WEB, 0)
        threading.Thread(target=cls.httpd.serve_forever, daemon=True).start()
        cls.base = f"http://127.0.0.1:{cls.httpd.server_address[1]}/"
        cls.pw = sync_playwright().start()
        cls.browser = cls.pw.chromium.launch(executable_path=chrome, headless=True, args=V.ARGS)
        cls.ctx = cls.browser.new_context()
        cls.page = cls.ctx.new_page()
        cls.page.goto(cls.base + V.PAGE)
        cls.scripts = "scripts/" if V.PAGE and (V.WEB / "scripts" / "cftr-viewer.js").exists() else ""
        cls.load = LOAD.replace("./cftr-viewer.js", "/" + cls.scripts + "cftr-viewer.js")

    @classmethod
    def tearDownClass(cls):
        cls.ctx.close()
        cls.browser.close()
        cls.pw.stop()
        cls.httpd.shutdown()
        cls.httpd.server_close()

    def call(self, body, arg=None):
        """Run `body` (a JavaScript arrow function taking (internals, arg)) in the page and return its JSON-able result."""
        return self.page.evaluate("async ([load, body, arg]) => { const m = await import(load); return (0, eval)(body)(m.internals, arg); }", [f"/{self.scripts}cftr-viewer.js", body, arg]) \
            if False else self.page.evaluate(f"async (arg) => {{ const I = await import('/{self.scripts}cftr-viewer.js').then(m => m.internals); return ({body})(I, arg); }}", arg)

    # ---------------------------------------------------------------- the curve
    def test_the_curve_passes_through_every_measured_point(self):
        rng = random.Random(4)
        for _ in range(5):
            n = rng.randint(3, 40)
            pts = [c for _ in range(n) for c in (rng.uniform(-40, 40), rng.uniform(-40, 40), rng.uniform(-40, 40))]
            out = self.call("(I, P) => { const pos = I.splinePositions(P); return { rings: pos.length / 3, sub: I.SUBDIV, at: Array.from({length: P.length / 3}, (_, i) => Array.from(pos.slice(3 * i * I.SUBDIV, 3 * i * I.SUBDIV + 3))) }; }", pts)
            self.assertEqual(out["rings"], (n - 1) * out["sub"] + 1)
            for i in range(n):
                for k in range(3):
                    self.assertAlmostEqual(out["at"][i][k], pts[3 * i + k], delta=2e-3, msg=(n, i, k))

    def test_between_the_points_the_curve_stays_close_to_them_and_does_not_loop(self):
        pts = helix_run(20)
        pos = self.call("(I, P) => Array.from(I.splinePositions(P))", pts)
        sub = self.call("(I) => I.SUBDIV")
        steps = [math.dist(pos[3 * r:3 * r + 3], pos[3 * r + 3:3 * r + 6]) for r in range(len(pos) // 3 - 1)]
        chords = [math.dist(pts[3 * i:3 * i + 3], pts[3 * i + 3:3 * i + 6]) for i in range(len(pts) // 3 - 1)]
        for i, chord in enumerate(chords):
            piece = steps[i * sub:(i + 1) * sub]
            self.assertAlmostEqual(sum(piece), chord, delta=chord * 0.25, msg="the path between two points is about as long as the straight line, a little more for a bend")
            self.assertGreater(min(piece), 0)
            self.assertLess(max(piece), chord, "no single piece is longer than the whole span: no overshoot")
        # no corner or reversal: an interpolating curve through points that turn 100 degrees a residue (an alpha helix) bends most AT the points, so the direction change between
        # neighbouring pieces peaks there (55 degrees at worst over all four real structures; 17 degrees on average). A real kink or a reversal is far beyond 75.
        self.assertLess(max_turn_degrees(pos), 75.0, "the curve has no corner")

    def test_no_run_of_any_real_structure_has_a_corner_or_a_reversal(self):
        worst = 0.0
        for st in V.DATA["structures"]:
            for seg in st["segments"]:
                if len(seg["xyz"]) < 9:
                    continue
                pos = self.call("(I, P) => Array.from(I.splinePositions(P))", seg["xyz"])
                worst = max(worst, max_turn_degrees(pos))
        self.assertLess(worst, 75.0)
        self.assertGreater(worst, 20.0, "a real helix does bend between pieces; a value this small would mean the curve is not following the data")

    def test_every_stretch_a_model_does_not_place_gets_a_dotted_straight_connector(self):
        stats = self.page.evaluate("JSON.parse(document.getElementById('cftr-canvas').dataset.mesh)")
        for st, got in zip(V.DATA["structures"], stats["dots"]):
            want = 0
            for a, b in zip(st["segments"], st["segments"][1:]):
                end, start = a["xyz"][-3:], b["xyz"][:3]
                want += max(1, math.floor(math.dist(end, start) / 3.0))
            self.assertEqual(got, want, st["id"])
            self.assertEqual(len(st["segments"]) - 1 > 0, got > 0)

    # ---------------------------------------------------------------- the mesh
    def test_the_mesh_counts_add_up_and_every_index_and_normal_is_valid(self):
        run1, run2 = helix_run(9), [x + 30 for x in helix_run(6)]
        st = {"segments": [{"start": 10, "xyz": run1}, {"start": 40, "xyz": run2}]}
        out = self.call("""(I, st) => { const sec = new Uint8Array(100); const ao = new Float32Array(100);
            const m = I.tubeMesh(st, sec, ao); const v = m.verts, idx = m.idx;
            let maxIdx = 0, minIdx = 1e9; for (const i of idx) { if (i > maxIdx) maxIdx = i; if (i < minIdx) minIdx = i; }
            let bad = 0, nanCount = 0; for (let k = 0; k < m.vertices; k++) { const nx = v[8 * k + 3], ny = v[8 * k + 4], nz = v[8 * k + 5]; const l = Math.hypot(nx, ny, nz);
              if (!(Math.abs(l - 1) < 1e-3)) bad++; for (let j = 0; j < 8; j++) if (Number.isNaN(v[8 * k + j])) nanCount++; }
            return { triangles: m.triangles, vertices: m.vertices, indexCount: idx.length, maxIdx, minIdx, badNormals: bad, nanCount, sides: I.SIDES, sub: I.SUBDIV }; }""", st)
        sides, sub = out["sides"], out["sub"]
        rings = (9 - 1) * sub + 1 + (6 - 1) * sub + 1
        self.assertEqual(out["vertices"], rings * sides + 2 * 2 * (sides + 1))
        self.assertEqual(out["triangles"], (rings - 2) * sides * 2 + 2 * 2 * sides)
        self.assertEqual(out["indexCount"], out["triangles"] * 3)
        self.assertLess(out["maxIdx"], out["vertices"])
        self.assertEqual(out["minIdx"], 0)
        self.assertEqual((out["badNormals"], out["nanCount"]), (0, 0))

    def test_the_tube_has_enough_sides_to_look_round_and_few_enough_for_a_phone(self):
        sides = self.call("(I) => I.SIDES")
        self.assertGreaterEqual(sides, 6)
        self.assertLessEqual(sides, 12)

    def test_a_helix_and_a_strand_are_drawn_wider_than_a_coil_and_a_turn_is_thin(self):
        pts = [3.8 * i if k == 0 else 0.0 for i in range(10) for k in range(3)]       # a straight run along x
        st = {"segments": [{"start": 1, "xyz": pts}]}
        widths = {}
        for name, code in (("coil", 0), ("helix", 1), ("strand", 2), ("turn", 3)):
            widths[name] = self.call("""(I, a) => { const sec = new Uint8Array(40).fill(a.code); const m = I.tubeMesh(a.st, sec, new Float32Array(40)); let r = 0;
                for (let k = 0; k < m.vertices; k++) { const y = m.verts[8 * k + 1], z = m.verts[8 * k + 2]; r = Math.max(r, Math.hypot(y, z)); } return r; }""", {"st": st, "code": code})
        shapes = self.call("(I) => I.SHAPES")
        for name, code in (("coil", 0), ("helix", 1), ("strand", 2), ("turn", 3)):
            self.assertAlmostEqual(widths[name], max(shapes[code]), delta=0.03, msg=name)
        self.assertGreater(widths["helix"], widths["strand"])
        self.assertGreater(widths["strand"], 2 * widths["coil"])
        self.assertLess(widths["coil"], 0.6)
        self.assertAlmostEqual(widths["turn"], widths["coil"], delta=0.03)
        self.assertGreater(widths["helix"], 1.0)

    def test_the_secondary_structure_codes_are_exactly_uniprots_ranges(self):
        data = V.DATA
        codes = self.call("(I, a) => Array.from(I.secondaryByResidue(a.data))", {"data": data})
        want = [0] * (data["length"] + 2)
        for s in data["uniprot"]["secondary"]:
            for n in range(s["start"], s["end"] + 1):
                if n < len(want):
                    want[n] = {"helix": 1, "strand": 2, "turn": 3}[s["type"]]
        self.assertEqual(codes, want)
        self.assertGreater(sum(1 for c in codes if c == 1), 200, "helix residues exist")
        self.assertGreater(sum(1 for c in codes if c == 2), 50, "strand residues exist")

    def test_the_real_meshes_are_inside_the_budget_and_the_crowding_is_a_fraction(self):
        out = self.call("""(I, a) => a.data.structures.map(st => { const L = a.data.length; const ao = I.crowding(st, L); const sec = I.secondaryByResidue(a.data);
              const m = I.tubeMesh(st, sec, ao); let lo = 1e9, hi = -1e9; for (const x of ao) { lo = Math.min(lo, x); hi = Math.max(hi, x); }
              return { id: st.id, triangles: m.triangles, vertices: m.vertices, lo, hi }; })""", {"data": V.DATA})
        for r in out:
            self.assertLess(r["triangles"], 150000, r["id"])
            self.assertGreater(r["triangles"], 50000, r["id"])
            self.assertLess(r["vertices"], 65536, r["id"])             # 16-bit indices on every structure
            self.assertGreaterEqual(r["lo"], 0)
            self.assertLessEqual(r["hi"], 1)
            self.assertGreater(r["hi"], 0.2, "a crowded core exists in every structure")
        # what the page publishes about its own mesh agrees with what the function makes
        stats = self.page.evaluate("JSON.parse(document.getElementById('cftr-canvas').dataset.mesh)")
        self.assertEqual(stats["triangles"], [r["triangles"] for r in out])


if __name__ == "__main__":
    unittest.main()
