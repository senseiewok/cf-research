"""Browser tests for the research tools: residue look-up, distance, UniProt layers, shareable address, exports, provenance.

Every answer a tool gives is recomputed here from the data module (test_data.py checks the data itself against the downloaded files), so a tool that says something
the data does not say fails.

Run from this folder:   python -m unittest test_research -v"""
import csv
import io
import json
import math
import random
import re
import struct
import threading
import unittest
from pathlib import Path

import support as V
from support import sync_playwright

DATA, BY_ID, norm = V.DATA, V.BY_ID, V.norm
U = DATA["uniprot"]
SEQ = U["sequence"]
NAMES = {"A": "alanine", "R": "arginine", "N": "asparagine", "D": "aspartate", "C": "cysteine", "Q": "glutamine", "E": "glutamate", "G": "glycine", "H": "histidine", "I": "isoleucine",
         "L": "leucine", "K": "lysine", "M": "methionine", "F": "phenylalanine", "P": "proline", "S": "serine", "T": "threonine", "W": "tryptophan", "Y": "tyrosine", "V": "valine"}
PLACED = {}
for st in DATA["structures"]:
    at = {}
    for seg in st["segments"]:
        for i in range(len(seg["xyz"]) // 3):
            at[seg["start"] + i] = seg["xyz"][3 * i:3 * i + 3]
    PLACED[st["id"]] = at
SINGLE = [v for v in U["variants"] if v["start"] == v["end"]]
SINGLE_POSITIONS = sorted({v["start"] for v in SINGLE})


def topology_kind(n):
    for t in U["transmembrane"]:
        if t["start"] <= n <= t["end"]:
            return "transmembrane"
    for t in U["topology"]:
        if t["start"] <= n <= t["end"]:
            return t["name"].lower()
    return "none"


class ResearchPage(V.Page):
    def wait_tools(self):
        self.page.wait_for_selector("#tools:not([hidden])")


@unittest.skipIf(sync_playwright is None, "Playwright is not installed in this Python")
class ResearchTests(unittest.TestCase):
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

    @classmethod
    def tearDownClass(cls):
        cls.browser.close()
        cls.pw.stop()
        cls.httpd.shutdown()
        cls.httpd.server_close()

    def open(self, path=V.PAGE, **ctx):
        ctx.setdefault("viewport", {"width": 1100, "height": 800})
        ctx.setdefault("accept_downloads", True)
        p = ResearchPage(self.browser, self.base, path=path, **ctx)
        self.addCleanup(p.close)
        return p

    # ---------------------------------------------------------------- loading
    def test_the_tools_appear_with_the_script_and_the_page_stays_clean(self):
        p = self.open()
        p.wait_drawn()
        p.wait_tools()
        self.assertEqual(p.problems, [])
        for sel in ("#cftr-lookup", "#cftr-measure", "#cftr-layer-topology", "#cftr-layer-variants", "#cftr-export-csv", "#cftr-export-png", "#cftr-export-json", "#cftr-copy-link"):
            self.assertTrue(p.page.is_visible(sel), sel)
        self.assertFalse(p.page.is_checked("#cftr-layer-topology"))
        self.assertFalse(p.page.is_checked("#cftr-layer-variants"))

    def test_without_javascript_the_tools_are_hidden_but_the_provenance_and_the_limits_are_not(self):
        p = self.open(java_script_enabled=False)
        self.assertTrue(p.page.is_hidden("#tools"))
        text = norm(p.page.inner_text("main"))
        self.assertIn("Where each number comes from", text)
        self.assertIn("What the model leaves out", text)
        self.assertIn("Only the alpha-carbon atom of each amino acid is drawn", text)
        self.assertIn("Cite this", text)

    # ---------------------------------------------------------------- the residue look-up
    def test_the_lookup_says_what_uniprot_and_the_files_say_for_a_sample_of_residues(self):
        p = self.open()
        p.wait_drawn()
        p.wait_tools()
        random.seed(11)
        sample = [1, 508, 1371, 1480, 654, 831, 81, 365, 400, 401, 462, 493, 1219, 1247, 1250] + random.sample(range(1, 1481), 21)      # includes residues inside the annotated ATP-binding sites
        for n in sample:
            p.page.fill("#cftr-lookup-n", str(n))
            p.page.click("#cftr-lookup button[type=submit]")
            text = norm(p.page.inner_text("#cftr-dossier"))
            letter = SEQ[n - 1]
            self.assertIn(f"Residue {n} \u00b7 {NAMES[letter]} ({letter})", text, n)
            tk = topology_kind(n)
            if tk == "transmembrane":
                self.assertIn("Topology: transmembrane segment", text, n)
            elif tk == "none":
                self.assertIn("no topology entry covers this residue", text, n)
            else:
                self.assertIn(f"Topology: {tk} side", text, n)
            sites = [s for s in U["atp_sites"] if s["start"] <= n <= s["end"]]
            self.assertEqual("ATP-binding site annotation: yes" in text, bool(sites), n)
            nv = len([v for v in U["variants"] if v["start"] <= n <= v["end"]])
            self.assertIn(f"Annotated natural variants here: {nv}", text, n)
            for st in DATA["structures"]:
                line = re.search(re.escape(st["id"]) + r": ([^\n]*?)(?= (?:6MSM|5UAK|8EIQ|8EJ1): |Mark residue|$)", text)
                self.assertIsNotNone(line, (n, st["id"]))
                said = line.group(1)
                if n in PLACED[st["id"]]:
                    self.assertTrue(said.startswith("placed in this model"), (n, st["id"], said))
                else:
                    self.assertFalse(said.startswith("placed"), (n, st["id"], said))
                    below = max([k for k in PLACED[st["id"]] if k < n], default=None)
                    above = min([k for k in PLACED[st["id"]] if k > n], default=None)
                    deleted = any(d["position"] == n and not d["file_residue"] for d in st["differences"])
                    if deleted:
                        self.assertIn("the file records it as", said, (n, st["id"]))
                    else:
                        self.assertIn(f"nearest placed residues are {below if below else 'none before'} and {above if above else 'none after'}", said, (n, st["id"], said))

    def test_the_508_dossier_is_the_story_of_the_deletion(self):
        p = self.open()
        p.wait_drawn()
        p.wait_tools()
        p.page.click("#cftr-lookup button[type=submit]")                # the field starts at 508
        text = norm(p.page.locator("#cftr-dossier").text_content())            # the group labels are upper-cased by the stylesheet, so read the text itself
        self.assertIn("Residue 508 \u00b7 phenylalanine (F)", text)
        self.assertIn("nucleotide-binding part 1", text)
        self.assertIn("6MSM: placed in this model", text)
        self.assertIn("5UAK: placed in this model", text)
        self.assertRegex(text, r"8EIQ: not in the file: the file records it as \"deletion\"")
        self.assertRegex(text, r"8EJ1: not in the file: the file records it as \"variant\"")
        self.assertIn("UniProt annotation (P13569, CC BY 4.0). This is annotation, not a clinical statement.", text)
        # the long UniProt wording is behind a disclosure, verbatim, and the short one shows first
        full = [v for v in U["variants"] if v["start"] == 508 and len(v["description"]) > 200][0]["description"]
        self.assertEqual(norm(p.page.locator("#cftr-dossier details p").first.text_content()), norm(full))      # closed disclosure: read the text, not what is shown
        self.assertIn("in CF and CBAVD", text)
        self.assertIn("VAR_080302 spanning residues 220 to 1480", text)             # an entry that covers many residues is described as spanning, not as sitting at 508

    def test_the_lookup_rejects_what_is_not_a_residue_number(self):
        p = self.open()
        p.wait_drawn()
        p.wait_tools()
        for bad in ("0", "1481", "-3", "2.5", ""):
            p.page.fill("#cftr-lookup-n", bad)
            p.page.click("#cftr-lookup button[type=submit]")
            self.assertEqual(norm(p.page.inner_text("#cftr-dossier")), "Please enter a whole number from 1 to 1480.", repr(bad))

    def test_marking_a_residue_from_the_lookup_marks_it_on_the_model(self):
        p = self.open()
        p.wait_drawn()
        p.wait_tools()
        p.page.fill("#cftr-lookup-n", "300")
        p.page.click("#cftr-lookup button[type=submit]")
        p.page.click("#cftr-dossier-mark")
        self.assertEqual(p.state()["walk"], 300)
        self.assertIn("Residue 300", p.page.inner_text("#cftr-walk-out"))

    def test_the_lookup_form_works_from_the_keyboard(self):
        p = self.open()
        p.wait_drawn()
        p.wait_tools()
        p.page.focus("#cftr-lookup-n")
        p.page.keyboard.press("Control+A")
        p.page.keyboard.type("1371")
        p.page.keyboard.press("Enter")
        text = norm(p.page.inner_text("#cftr-dossier"))
        self.assertIn("Residue 1371 \u00b7 glutamate (E)", text)
        self.assertIn("glutamine here where UniProt has glutamate", text)       # the engineered change the files declare
        self.assertIn("5UAK: placed in this model", text)

    # ---------------------------------------------------------------- distance
    def test_the_distance_equals_the_distance_recomputed_from_the_coordinates(self):
        p = self.open()
        p.wait_drawn()
        p.wait_tools()
        random.seed(5)
        for sid in ("6MSM", "5UAK", "8EIQ", "8EJ1"):
            p.page.click(f"#cftr-picks button:nth-child({list(BY_ID).index(sid) + 1})")
            placed = sorted(PLACED[sid])
            pairs = [(random.choice(placed), random.choice(placed)) for _ in range(6)]
            for a, b in pairs:
                p.page.fill("#cftr-measure-a", str(a))
                p.page.fill("#cftr-measure-b", str(b))
                p.page.click("#cftr-measure button[type=submit]")
                out = norm(p.page.inner_text("#cftr-measure-out"))
                pa, pb = PLACED[sid][a], PLACED[sid][b]
                want = math.dist(pa, pb)
                got = float(re.search(r"is ([\d.]+) \u00c5", out).group(1))
                self.assertAlmostEqual(got, want, delta=0.051, msg=(sid, a, b, out))
                self.assertIn("alpha carbon to alpha carbon", out)
                self.assertIn("not a contact or a bond", out)

    def test_a_distance_to_a_residue_that_is_not_placed_is_refused_and_names_it(self):
        p = self.open()
        p.wait_drawn()
        p.wait_tools()
        p.page.click("#cftr-picks button:nth-child(3)")                  # 8EIQ lacks 508
        p.page.fill("#cftr-measure-a", "507")
        p.page.fill("#cftr-measure-b", "508")
        p.page.click("#cftr-measure button[type=submit]")
        self.assertEqual(norm(p.page.inner_text("#cftr-measure-out")), "Not measurable in 8EIQ: residue 508 is not placed in this model.")
        p.page.fill("#cftr-measure-a", "700")
        p.page.fill("#cftr-measure-b", "701")
        p.page.click("#cftr-measure button[type=submit]")
        self.assertIn("residue 700 and 701 are not placed", norm(p.page.inner_text("#cftr-measure-out")))
        p.page.fill("#cftr-measure-a", "0")
        p.page.click("#cftr-measure button[type=submit]")
        self.assertIn("Please enter two whole numbers from 1 to 1480.", norm(p.page.inner_text("#cftr-measure-out")))

    def test_a_measured_distance_draws_a_line_on_the_model_and_clearing_it_removes_it(self):
        p = self.open()
        p.wait_drawn()
        p.wait_tools()
        p.hold()
        p.page.fill("#cftr-measure-a", "100")
        p.page.fill("#cftr-measure-b", "300")
        before = p.look()
        p.page.click("#cftr-measure button[type=submit]")
        p.settle()
        self.assertIn("measure", p.state()["layers"])
        self.assertNotEqual(p.look(), before, "the measured line and its end points are drawn")
        p.page.fill("#cftr-measure-a", "0")
        p.page.click("#cftr-measure button[type=submit]")
        p.settle()
        self.assertNotIn("measure", p.state()["layers"])
        self.assertEqual(p.look(), before, "an invalid request clears the line")

    # ---------------------------------------------------------------- UniProt layers
    def test_the_topology_layer_recolours_the_chain_and_says_where_the_sides_come_from(self):
        p = self.open()
        p.wait_drawn()
        p.wait_tools()
        p.hold()
        p.settle()
        plain = p.look()
        self.assertTrue(p.page.is_hidden("#cftr-topology-legend"))
        p.page.check("#cftr-layer-topology")
        p.settle()
        self.assertTrue(p.page.is_visible("#cftr-topology-legend"))
        self.assertIn("topology", p.state()["layers"])
        self.assertNotEqual(p.look(), plain)
        self.assertGreater(p.look().changed_pixels(plain), 3000)
        note = norm(p.page.locator("#tools").inner_text())
        self.assertIn(f"{len(U['topology'])} topological domains and {len(U['transmembrane'])} transmembrane segments", note)
        self.assertIn("not from these coordinates", note)
        self.assertIn("This is not a membrane drawn from the data.", note)
        # the ruler follows the colour function: its coloured runs split where the topology changes
        runs_topology = p.page.locator("#cftr-ruler-svg rect[height='16']").count()
        p.page.uncheck("#cftr-layer-topology")
        p.settle()
        runs_parts = p.page.locator("#cftr-ruler-svg rect[height='16']").count()
        self.assertGreater(runs_topology, runs_parts + 5)
        self.assertEqual(p.look(), plain, "switching the layer off gives the original colours back")
        self.assertTrue(p.page.is_hidden("#cftr-topology-legend"))

    def test_the_variant_layer_marks_exactly_the_single_residue_positions_this_structure_places(self):
        p = self.open()
        p.wait_drawn()
        p.wait_tools()
        p.hold()
        p.settle()
        plain = p.look()
        spanning = len(U["variants"]) - len(SINGLE)
        for i, sid in enumerate(BY_ID, start=1):
            p.page.click(f"#cftr-picks button:nth-child({i})")
            p.page.check("#cftr-layer-variants") if not p.page.is_checked("#cftr-layer-variants") else None
            p.settle()
            here = len([n for n in SINGLE_POSITIONS if n in PLACED[sid]])
            summary = norm(p.page.inner_text("#cftr-variant-summary"))
            self.assertEqual(summary, f"{len(U['variants'])} UniProt entries: {len(SINGLE)} at a single residue ({len(SINGLE_POSITIONS)} different positions) and {spanning} spanning several residues, "
                                      f"which are not drawn (look the residues up to see them). {sid} places {here} of the {len(SINGLE_POSITIONS)} single-residue positions; the other "
                                      f"{len(SINGLE_POSITIONS) - here} are not placed in it, so nothing is drawn for them.")
        self.assertIn("variants", p.state()["layers"])
        p.page.click("#cftr-picks button:nth-child(1)")
        p.settle()
        with_beads = p.look()
        self.assertNotEqual(with_beads, plain, "beads are drawn")
        p.page.uncheck("#cftr-layer-variants")
        p.settle()
        self.assertEqual(p.look(), plain)

    def test_the_variant_layer_carries_uniprots_own_caveat_and_no_clinical_wording(self):
        p = self.open()
        p.wait_drawn()
        p.wait_tools()
        text = norm(p.page.locator("#tools").inner_text())
        self.assertIn("This is annotation, not a clinical, frequency or severity statement, and it says nothing about any person.", text)
        self.assertIn("Any medical or genetic information is provided for research, educational and informational purposes only.", text)
        self.assertIn("UniProt, CC BY 4.0", text)
        for banned in ("pathogenic", "benign", "severe", "mild", "carrier", "prevalence", "frequency of", "your risk"):
            self.assertNotIn(banned, text.lower().replace("not a clinical, frequency or severity statement", ""), banned)

    # ---------------------------------------------------------------- the address
    def test_the_address_holds_the_view_and_a_reload_restores_it(self):
        p = self.open()
        p.wait_drawn()
        p.wait_tools()
        p.page.click("#cftr-picks button:nth-child(3)")
        p.page.click("#cftr-overlay")
        p.page.fill("#cftr-lookup-n", "508")
        p.page.click("#cftr-lookup button[type=submit]")
        p.page.click("#cftr-dossier-mark")
        p.page.check("#cftr-layer-variants")
        p.page.focus("#cftr-canvas")
        p.page.keyboard.press("ArrowRight")
        p.page.wait_for_timeout(600)
        h = p.page.evaluate("location.hash")
        q = dict(kv.split("=") for kv in h[1:].split("&"))
        self.assertEqual((q["s"], q["r"], q["o"], q["l"]), ("8EIQ", "508", "1", "variants"))
        self.assertTrue({"y", "p", "z"} <= set(q), "a view the visitor moved is in the address")
        q2 = self.open(path=V.PAGE + h)
        q2.wait_drawn()
        q2.wait_tools()
        s = q2.state()
        self.assertEqual((s["structure"], s["walk"], s["overlay"], s["layers"]), ("8EIQ", 508, True, ["variants"]))
        self.assertFalse(s["playing"], "a shared view stands still")
        self.assertTrue(q2.page.is_checked("#cftr-layer-variants"))
        self.assertEqual(q2.page.evaluate("location.hash"), h)

    def test_an_address_that_is_wrong_or_hostile_is_ignored_safely(self):
        for bad in ("#s=NOPE&r=99999&y=abc&p=5&z=99&l=%3Cscript%3E,topology,variants;DROP", "#r=1.5&o=1&s=6MSM", "#s=&r=&y=&p=&z=", "#l=" + "a" * 500, "#%E0%A4%A", "#s=8EIQ&r=0"):
            p = self.open(path=V.PAGE + bad)
            p.wait_drawn()
            self.assertEqual([m for m in p.problems if "pageerror" in m or "error" in m], [], bad)
            s = p.state()
            self.assertIn(s["structure"], BY_ID, bad)
            self.assertTrue(s["walk"] is None or 1 <= s["walk"] <= 1480, bad)
            for name in s["layers"]:
                self.assertIn(name, ("topology", "variants"), bad)

    def test_the_copy_link_button_copies_this_view_to_the_clipboard(self):
        p = self.open()
        p.ctx.grant_permissions(["clipboard-read", "clipboard-write"], origin=self.base.rstrip("/"))
        p.wait_drawn()
        p.wait_tools()
        p.page.click("#cftr-picks button:nth-child(2)")
        p.page.wait_for_timeout(400)
        p.page.click("#cftr-copy-link")
        p.wait_until(lambda: "copied" in p.page.inner_text("#cftr-export-status"), "the confirmation")
        copied = p.page.evaluate("navigator.clipboard.readText()")
        self.assertEqual(copied, p.page.evaluate("location.href"))
        self.assertIn("#s=5UAK", copied)

    # ---------------------------------------------------------------- exports
    def test_the_csv_is_the_placed_residues_with_their_coordinates(self):
        p = self.open()
        p.wait_drawn()
        p.wait_tools()
        for sid in ("6MSM", "8EIQ"):
            p.page.click(f"#cftr-picks button:nth-child({list(BY_ID).index(sid) + 1})")
            with p.page.expect_download() as dl:
                p.page.click("#cftr-export-csv")
            d = dl.value
            self.assertEqual(d.suggested_filename, f"{sid}-alpha-carbons.csv")
            rows = list(csv.DictReader(io.StringIO(Path(d.path()).read_text(encoding="utf-8"))))
            self.assertEqual(len(rows), BY_ID[sid]["residues_traced"])
            self.assertEqual(list(rows[0]), ["structure", "residue", "amino_acid", "x_angstrom", "y_angstrom", "z_angstrom", "uniprot_part", "uniprot_topology"])
            self.assertEqual({r["structure"] for r in rows}, {sid})
            self.assertEqual([int(r["residue"]) for r in rows], sorted(PLACED[sid]))
            for r in random.Random(3).sample(rows, 40):
                n = int(r["residue"])
                self.assertEqual(r["amino_acid"], SEQ[n - 1] if n != 1371 else SEQ[n - 1])
                self.assertEqual([float(r["x_angstrom"]), float(r["y_angstrom"]), float(r["z_angstrom"])], PLACED[sid][n])
                self.assertEqual(r["uniprot_topology"], topology_kind(n))
            self.assertEqual(("508" in {r["residue"] for r in rows}), sid == "6MSM")
            self.assertIn("status", p.page.inner_text("#cftr-export-status") or "status") if False else None
            self.assertIn(f"Saved {len(rows)} placed residues of {sid}.", p.page.inner_text("#cftr-export-status"))

    def test_the_json_view_round_trips_with_the_address(self):
        p = self.open()
        p.wait_drawn()
        p.wait_tools()
        p.page.click("#cftr-picks button:nth-child(4)")
        p.page.check("#cftr-layer-topology")
        with p.page.expect_download() as dl:
            p.page.click("#cftr-export-json")
        obj = json.loads(Path(dl.value.path()).read_text(encoding="utf-8"))
        self.assertEqual(dl.value.suggested_filename, "cftr-view.json")
        self.assertEqual((obj["structure"], obj["layers"]), ("8EJ1", ["topology"]))
        self.assertIn("RCSB PDB (CC0)", obj["data"])
        self.assertIn("UniProt P13569 (CC BY 4.0)", obj["data"])

    def test_the_png_is_a_picture_with_its_label_baked_in(self):
        p = self.open()
        p.wait_drawn()
        p.wait_tools()
        p.hold()
        p.page.click("#cftr-picks button:nth-child(3)")
        p.settle()
        box = p.page.locator("#cftr-canvas").bounding_box()
        with p.page.expect_download() as dl:
            p.page.click("#cftr-export-png")
        self.assertEqual(dl.value.suggested_filename, "cftr-8EIQ.png")
        data = Path(dl.value.path()).read_bytes()
        self.assertEqual(data[:8], b"\x89PNG\r\n\x1a\n")
        w, h = struct.unpack(">II", data[16:24])
        self.assertEqual(h - 64, round(box["height"]))
        self.assertEqual(w, round(box["width"]))
        import pngpixels
        ww, hh, px = pngpixels.decode_png(data)
        label = px[(hh - 40) * ww: (hh - 40) * ww + ww]
        self.assertTrue(any(max(c) > 150 for c in label), "text is drawn in the label bar")
        top = pngpixels.hue_counts(data)
        self.assertGreater(sum(top.values()), 500, "the protein is in the picture")

    # ---------------------------------------------------------------- provenance
    def test_the_provenance_table_matches_the_data(self):
        p = self.open(java_script_enabled=False)
        rows = {r.query_selector("th").inner_text(): norm(r.inner_text()) for r in p.page.query_selector_all(".cftr-prov tbody tr")}
        sts = DATA["structures"]
        placed_total = sum(s["residues_traced"] for s in sts)
        declared = sum(1 for s in sts for d in s["differences"] if d["file_residue"])
        self.assertIn("; ".join(f"{s['id']} {s['residues_traced']}" for s in sts), rows["Residues placed"])
        self.assertIn(f"{placed_total:,} placed residues checked; {declared} differ from UniProt, all declared by their files", rows["Residue numbering"])
        for s in sts[1:]:
            self.assertIn(f"{s['id']} {s['superposition_rmsd_A']:.2f} \u00c5 over {s['superposition_residues']} residues", rows["Distance after fitting"])
        spanning = len(U["variants"]) - len(SINGLE)
        self.assertIn(f"{len(U['variants'])} entries: {len(SINGLE)} at a single residue ({len(SINGLE_POSITIONS)} different positions) and {spanning} spanning several residues", rows["Variant positions"])
        self.assertIn(f"{len(U['transmembrane'])} transmembrane segments; {len(U['topology'])} topological domains", rows["Topology and helices"])
        self.assertIn(f"entry version {U['entry_version']}, annotation updated {U['last_annotation_update']}", rows["Length and part boundaries"])
        for what in ("Residues placed", "Method, resolution, date", "Paper", "Stretches not placed", "Residue numbering", "Distance after fitting", "Positions in the view", "Variant positions"):
            self.assertIn(what, rows)

    def test_the_table_of_structures_gives_each_paper_and_each_entry_link_correctly(self):
        p = self.open(java_script_enabled=False)
        rows = {r.query_selector("th").inner_text(): r for r in p.page.query_selector_all(".cftr-table:not(.cftr-prov) tbody tr")}
        for st in DATA["structures"]:
            row = rows[st["id"]]
            text = norm(row.inner_text())
            self.assertIn(f"{st['journal']} {st['paper_year']}, PubMed {st['pubmed']}", text, st["id"])
            hrefs = [a.get_attribute("href") for a in row.query_selector_all("a")]
            self.assertIn(f"https://www.rcsb.org/structure/{st['id']}", hrefs)
            self.assertIn(f"https://pubmed.ncbi.nlm.nih.gov/{st['pubmed']}/", hrefs)

    def test_the_list_of_what_the_model_leaves_out_is_complete(self):
        p = self.open(java_script_enabled=False)
        items = [norm(li.inner_text()) for li in p.page.query_selector_all("ul.cftr-leaves li")]
        self.assertEqual(len(items), 7)
        wanted = ["Only the alpha-carbon atom of each amino acid is drawn", "Only chain A of model 1 of each entry is used", "dotted grey lines across it are straight connectors",
                  "the three medicines) are not drawn", "solved from different samples and under different conditions", "author numbers, which equal UniProt", "no motion, no function and no mechanism"]
        for w in wanted:
            self.assertTrue(any(w in i for i in items), w)

    def test_the_cite_block_names_every_entry_and_every_paper_once(self):
        p = self.open(java_script_enabled=False)
        text = p.page.inner_text("pre.cftr-cite")
        for st in DATA["structures"]:
            self.assertEqual(len(re.findall(rf"^{st['id']}: ", text, re.M)), 1, st["id"])
            self.assertIn(f"PubMed {st['pubmed']}", text)
        self.assertEqual(len(re.findall(r"PubMed \d+", text)), 4)
        self.assertIn(f"UniProt P13569, entry version {U['entry_version']}", text)
        self.assertIn("CC0 1.0", text)
        self.assertIn("CC BY 4.0", text)

    # ---------------------------------------------------------------- phones and keyboards
    def test_a_phone_needs_no_sideways_scroll_and_the_tools_are_big_enough_to_use(self):
        p = self.open(viewport={"width": 360, "height": 740}, has_touch=True, is_mobile=True, device_scale_factor=2)
        p.wait_drawn()
        p.wait_tools()
        self.assertEqual(p.problems, [])
        self.assertLessEqual(p.page.evaluate("document.documentElement.scrollWidth"), 360)
        small = []
        for el in p.page.query_selector_all("#tools button, #tools input:not([type=checkbox]), #tools summary"):
            b = el.bounding_box()
            if b and (b["height"] < 43.5):
                small.append((el.get_attribute("id") or el.inner_text(), b["height"]))
        self.assertEqual(small, [])
        for sel in (".cftr-prov", ".cftr-table"):
            wrap = p.page.locator(f"{sel}").first.locator("xpath=ancestor::div[contains(@class,'cftr-tablewrap')]")
            self.assertEqual(wrap.get_attribute("tabindex"), "0")

    def test_every_control_in_the_tools_has_a_name(self):
        p = self.open()
        p.wait_drawn()
        p.wait_tools()
        for el in p.page.query_selector_all("#tools button, #tools input"):
            name = p.page.evaluate("""el => (el.getAttribute('aria-label') || el.innerText || (el.labels && el.labels[0] && el.labels[0].innerText) || '').trim()""", el)
            self.assertTrue(name, f"unnamed control {el.get_attribute('id')}")


if __name__ == "__main__":
    unittest.main()
