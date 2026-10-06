"""Browser tests for the CFTR structure viewer, run over HTTP with a strict Content-Security-Policy.

What they prove, in a real Chromium with a software GL: the picture is drawn and coloured by part, residue 508 is marked and its absence is shown, every control works
with a real key or pointer, reduced motion draws one still frame and then stops, a lost WebGL context is rebuilt, the page needs no sideways scroll on a phone, it reads
without JavaScript, and the numbers written in the captions are the numbers in the data.

Run from this folder:   python -m unittest test_viewer -v        (needs Playwright for Python and a Chromium; see support.py for running them against another page)"""
import hashlib
import json
import re
import threading
import unittest

from support import (ARGS, BY_ID, DATA, PAGE, SRC, WEB, Look, Page, find_chrome, model_server, norm, pngpixels, sync_playwright)  # noqa: F401


@unittest.skipIf(sync_playwright is None, "Playwright is not installed in this Python")
class ViewerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        chrome = find_chrome()
        if not chrome:
            raise unittest.SkipTest("no Chromium found")
        cls.httpd = model_server.make_server(WEB, 0)
        threading.Thread(target=cls.httpd.serve_forever, daemon=True).start()
        cls.base = f"http://127.0.0.1:{cls.httpd.server_address[1]}/"
        cls.pw = sync_playwright().start()
        cls.browser = cls.pw.chromium.launch(executable_path=chrome, headless=True, args=ARGS)

    @classmethod
    def tearDownClass(cls):
        cls.browser.close()
        cls.pw.stop()
        cls.httpd.shutdown()
        cls.httpd.server_close()

    def open(self, **ctx):
        ctx.setdefault("viewport", {"width": 1100, "height": 800})
        p = Page(self.browser, self.base, **ctx)
        self.addCleanup(p.close)
        return p

    # ---------------------------------------------------------------- loading and drawing
    def test_loads_clean_and_stays_on_this_computer(self):
        p = self.open()
        p.wait_drawn()
        self.assertEqual(p.problems, [])
        foreign = [u for u in p.requests if not u.startswith(self.base)]
        self.assertEqual(foreign, [])
        self.assertFalse(any("fetch" in u or u.endswith(".json") for u in p.requests), "the data is a module, not a request")
        self.assertTrue(p.page.is_visible("#cftr-viewer"))
        self.assertTrue(p.page.is_hidden("#cftr-fallback"))

    def test_the_first_picture_is_coloured_by_part_and_marks_residue_508(self):
        p = self.open()
        p.wait_drawn()
        p.hold()                                                      # hold still so the pixels can be read
        png = p.shot()
        self.assertGreater(pngpixels.nonblank(png), 0.01)
        hues = pngpixels.hue_counts(png)
        for family in ("cyan", "violet", "amber", "blue"):
            self.assertGreater(hues[family], 30, f"no {family} pixels: {dict(hues)}")
        self.assertGreater(hues["rose"], 20, f"no rose marker pixels: {dict(hues)}")

    def test_the_absence_of_508_is_drawn_differently_from_its_presence(self):
        p = self.open()
        p.wait_drawn()
        p.page.click("#cftr-play")
        p.page.click("#cftr-picks button:nth-child(1)")
        p.settle()
        present = pngpixels.hue_counts(p.shot())["rose"]
        p.page.click("#cftr-picks button:nth-child(3)")
        p.settle()
        absent = pngpixels.hue_counts(p.shot())["rose"]
        self.assertGreater(absent, 100, "an empty ring should still be drawn where 508 would be")
        self.assertLess(absent, present * 0.8, f"the ring alone ({absent} rose pixels) should have fewer than the ring plus the filled dot ({present})")
        self.assertGreater(present - absent, 100, "the filled dot should add a visible amount")

    def test_the_tube_is_shaded_and_the_mesh_stays_within_its_budget(self):
        """The chain is a lit 3D tube, not a flat line: the cyan part shows a real spread of brightness between lit and shadowed faces. The mesh that makes it
        stays within the budget a phone can draw, and the canvas publishes the numbers so anyone can check them."""
        import colorsys
        p = self.open()
        p.wait_drawn()
        p.hold()
        mesh = json.loads(p.page.get_attribute("#cftr-canvas", "data-mesh"))
        self.assertEqual(len(mesh["triangles"]), len(DATA["structures"]))
        for n in mesh["triangles"]:
            self.assertTrue(20000 < n <= 150000, f"{n} triangles in one structure")
        self.assertLess(mesh["buildMs"], 1500, "built once at load, even on this software-rendered test machine (the laptop budget is 250 ms)")
        w, h, px = pngpixels.pixels(p.shot())
        cyan = sorted(hsv[2] for hsv in (colorsys.rgb_to_hsv(r / 255, g / 255, b / 255) for r, g, b in px) if 165 <= hsv[0] * 360 < 197 and hsv[1] >= 0.45 and hsv[2] >= 0.12)
        self.assertGreater(len(cyan), 500)
        dim, bright = cyan[len(cyan) // 10], cyan[9 * len(cyan) // 10]
        self.assertGreater(bright - dim, 0.15, f"a shaded tube has lit and shadowed cyan (10th to 90th percentile brightness: {dim:.2f} to {bright:.2f}; measured 0.71 to 0.93 here)")

    def test_the_page_says_what_the_smoothing_and_the_bands_are(self):
        """The curve between measured points and the helix/strand bands are drawing choices and UniProt's annotation; the page must say so where the colours are explained."""
        p = self.open()
        p.wait_drawn()
        self.assertTrue(p.page.is_visible("#cftr-viewer .cftr-shape-note"))
        note = norm(p.page.inner_text("#cftr-viewer .cftr-shape-note"))
        for needle in ("one atom per amino acid", "drawing choice, not data", "UniProt annotates as a helix (wider) or a strand (narrower)", "not from a computation on these coordinates", "says nothing about the protein"):
            self.assertIn(needle, note, needle)
        legend = norm(p.page.inner_text("#cftr-legend"))
        self.assertIn("Flat band: helix or strand, as UniProt annotates", legend)
        self.assertIn("Dotted grey: not placed in this structure", legend)
        self.assertEqual(len(DATA["uniprot"]["secondary"]), 118, "the bands come from the data file, which carries UniProt's ranges")
        self.assertTrue(all(s["type"] in ("helix", "strand", "turn") for s in DATA["uniprot"]["secondary"]))

    # ---------------------------------------------------------------- controls
    def test_choosing_a_structure_changes_what_is_drawn_and_what_is_written(self):
        p = self.open()
        p.wait_drawn()
        p.hold()
        first = p.look()
        for i, sid in enumerate(BY_ID, start=1):
            p.page.click(f"#cftr-picks button:nth-child({i})")
            p.settle()
            self.assertEqual(p.state()["structure"], sid)
            self.assertEqual(p.page.get_attribute(f"#cftr-picks button:nth-child({i})", "aria-pressed"), "true")
            self.assertEqual(sum(1 for b in p.page.query_selector_all("#cftr-picks button") if b.get_attribute("aria-pressed") == "true"), 1)
            text = norm(p.page.inner_text("#cftr-readout"))
            self.assertIn(sid, text)
            self.assertIn(f"{BY_ID[sid]['residues_traced']} of {DATA['length']} residues", text)
            self.assertIn("absent" if BY_ID[sid]["res508"] is None else "phenylalanine", text.split("Residue 508:")[1])
        self.assertNotEqual(p.look(), first)

    def test_the_walk_says_placed_not_placed_or_absent_and_draws_a_marker(self):
        p = self.open()
        p.wait_drawn()
        p.page.click("#cftr-play")
        p.page.click("#cftr-walk-508")
        p.settle()
        self.assertIn("phenylalanine", p.page.inner_text("#cftr-walk-out"))
        self.assertEqual(p.page.get_attribute("#cftr-walk", "value") or p.page.input_value("#cftr-walk"), "508")
        p.page.click("#cftr-walk-clear")
        p.settle()
        self.assertIn("No residue is marked", p.page.inner_text("#cftr-walk-out"))
        cleared = p.look()
        p.page.focus("#cftr-walk")
        p.page.keyboard.press("Home")                                  # residue 1 is placed in 6MSM
        p.settle()
        self.assertIn("Residue 1:", p.page.inner_text("#cftr-walk-out"))
        self.assertIn("placed in this model", p.page.inner_text("#cftr-walk-out"))
        self.assertNotEqual(p.look(), cleared, "a placed residue should add a marker")
        p.page.keyboard.press("End")                                   # residue 1480 is not placed in any of the four
        p.settle()
        out = p.page.inner_text("#cftr-walk-out")
        self.assertIn("Residue 1480", out)
        self.assertIn("not placed in this model", out)
        self.assertIn("Residue 1480", p.page.get_attribute("#cftr-walk", "aria-valuetext"))
        p.page.click("#cftr-picks button:nth-child(3)")
        p.page.click("#cftr-walk-508")
        p.settle()
        self.assertIn("absent from this structure", p.page.inner_text("#cftr-walk-out"))
        self.assertIn("507 straight to 509", p.page.inner_text("#cftr-walk-out"))

    def test_every_unplaced_walk_position_names_a_real_stretch(self):
        """For each structure, a residue inside each unplaced stretch must be reported with that stretch's own limits."""
        p = self.open()
        p.wait_drawn()
        for i, (sid, st) in enumerate(BY_ID.items(), start=1):
            p.page.click(f"#cftr-picks button:nth-child({i})")
            placed = set()
            for seg in st["segments"]:
                placed.update(range(seg["start"], seg["start"] + len(seg["xyz"]) // 3))
            n = 700
            self.assertNotIn(n, placed)
            a = n
            while a - 1 not in placed and a > 1:
                a -= 1
            b = n
            while b + 1 not in placed and b < DATA["length"]:
                b += 1
            p.page.evaluate("n => { const r = document.getElementById('cftr-walk'); r.value = n; r.dispatchEvent(new Event('input', {bubbles: true})); }", n)
            self.assertIn(f"the stretch {a} to {b} is not modelled", p.page.inner_text("#cftr-walk-out"), sid)

    def test_overlay_is_only_for_the_others_and_reports_the_computed_distance(self):
        p = self.open()
        p.wait_drawn()
        self.assertTrue(p.page.is_disabled("#cftr-overlay"), "6MSM cannot be laid over itself")
        for i, sid in ((2, "5UAK"), (3, "8EIQ"), (4, "8EJ1")):
            p.page.click(f"#cftr-picks button:nth-child({i})")
            self.assertFalse(p.page.is_disabled("#cftr-overlay"))
            p.page.click("#cftr-overlay") if p.page.get_attribute("#cftr-overlay", "aria-pressed") == "false" else None
            text = norm(p.page.inner_text("#cftr-readout"))
            st = BY_ID[sid]
            self.assertIn(f"{st['superposition_residues']} membrane-spanning residues, {st['superposition_rmsd_A']:.2f} \u00c5", text)
            self.assertEqual(p.page.get_attribute("#cftr-overlay", "aria-pressed"), "true")
        p.page.click("#cftr-picks button:nth-child(1)")
        self.assertEqual(p.page.get_attribute("#cftr-overlay", "aria-pressed"), "false")

    def test_the_keyboard_works_on_the_picture(self):
        p = self.open()
        p.wait_drawn()
        p.hold()
        p.page.focus("#cftr-canvas")
        turn0, tilt0 = p.page.input_value("#cftr-turn"), p.page.input_value("#cftr-tilt")
        p.page.keyboard.press("ArrowRight")
        p.page.keyboard.press("ArrowDown")
        self.assertNotEqual(p.page.input_value("#cftr-turn"), turn0)
        self.assertNotEqual(p.page.input_value("#cftr-tilt"), tilt0)
        p.page.keyboard.press("3")
        self.assertEqual(p.state()["structure"], "8EIQ")
        p.page.keyboard.press("r")
        self.assertEqual(p.page.input_value("#cftr-turn"), "34", "R returns to the starting view (0.6 radians is 34 degrees)")
        self.assertEqual(p.page.input_value("#cftr-tilt"), "1.35")
        p.page.keyboard.press("9")                                     # there is no ninth structure: nothing may change
        self.assertEqual(p.state()["structure"], "8EIQ")

    def test_a_pointer_drag_turns_it_and_the_sliders_do_the_same_without_dragging(self):
        p = self.open()
        p.wait_drawn()
        p.hold()
        box = p.page.locator("#cftr-canvas").bounding_box()
        x, y = box["x"] + box["width"] / 2, box["y"] + box["height"] / 2
        before = p.look()
        turn0 = p.page.input_value("#cftr-turn")
        p.page.mouse.move(x, y)
        p.page.mouse.down()
        p.page.mouse.move(x + 120, y + 40, steps=6)
        p.page.mouse.up()
        p.settle()
        self.assertNotEqual(p.page.input_value("#cftr-turn"), turn0)
        self.assertNotEqual(p.look(), before)
        p.page.focus("#cftr-turn")
        mid = p.look()
        p.page.keyboard.press("ArrowRight")
        p.page.keyboard.press("ArrowRight")
        p.page.keyboard.press("ArrowRight")
        p.settle()
        self.assertNotEqual(p.look(), mid, "the Turn slider must move the picture")
        mid = p.look()
        p.page.focus("#cftr-tilt")
        for _ in range(8):
            p.page.keyboard.press("ArrowRight")
        p.settle()
        self.assertNotEqual(p.look(), mid, "the Tilt slider must move the picture")

    def test_zoom_buttons_change_the_picture_and_reset_restores_it(self):
        p = self.open()
        p.wait_drawn()
        p.hold()
        p.page.click("#cftr-reset")
        p.settle()
        home = p.look()
        p.page.click("#cftr-zoom-in")
        p.page.click("#cftr-zoom-in")
        p.settle()
        zoomed = p.look()
        self.assertNotEqual(zoomed, home)
        p.page.click("#cftr-reset")
        p.settle()
        self.assertEqual(p.look(), home)

    def test_an_absent_residue_is_a_dashed_ring_and_a_present_one_is_a_solid_ring(self):
        import math
        p = self.open()
        p.wait_drawn()
        p.hold()
        def ring_transitions(png):
            """How many times the marker's ring switches between painted and unpainted around its circumference: a solid ring none, a dashed ring one per dash edge."""
            pts = pngpixels.rose_points(png)
            self.assertGreater(len(pts), 150, "a rose marker is on the canvas")
            cx = sum(x for x, _ in pts) / len(pts)
            cy = sum(y for _, y in pts) / len(pts)
            radii = sorted(math.hypot(x - cx, y - cy) for x, y in pts)
            outer = radii[int(len(radii) * 0.9)]                          # the ring's radius: most rose pixels lie on or inside it
            have = set(pts)
            covered = []
            for k in range(360):
                a = 2 * math.pi * k / 360
                covered.append(any((round(cx + r * math.cos(a)), round(cy + r * math.sin(a))) in have for r in (outer - 3, outer - 2.5, outer - 2, outer - 1.5, outer - 1, outer - 0.5, outer)))
            return sum(1 for k in range(360) if covered[k] != covered[(k + 1) % 360]), sum(covered) / 360
        p.page.click("#cftr-picks button:nth-child(1)")
        p.settle()
        present_edges, present_cover = ring_transitions(p.shot())
        p.page.click("#cftr-picks button:nth-child(3)")
        p.settle()
        absent_edges, absent_cover = ring_transitions(p.shot())
        self.assertLessEqual(present_edges, 6, f"a residue the structure has is marked with a closed ring ({present_edges} gaps)")
        self.assertGreater(present_cover, 0.9)
        self.assertGreaterEqual(absent_edges, 12, f"a residue the structure lacks is marked with a broken ring, so it cannot be mistaken for a present one ({absent_edges} edges)")
        self.assertEqual(json.loads(p.page.get_attribute("#cftr-canvas", "data-draws"))["absentRingMode"], 3)

    def test_the_walk_marker_is_a_white_dot_on_the_residue_and_leaves_with_the_marker(self):
        p = self.open()
        p.wait_drawn()
        p.hold()
        p.page.click("#cftr-walk-clear")
        p.settle(250)
        before = pngpixels.white_count(p.shot())
        p.page.fill("#cftr-lookup-n", "300") if p.page.query_selector("#cftr-lookup-n") else None
        p.page.evaluate("() => { const r = document.getElementById('cftr-walk'); r.value = 100; r.dispatchEvent(new Event('input', {bubbles: true})); }")
        p.settle(250)
        with_marker = pngpixels.white_count(p.shot())
        self.assertGreater(with_marker - before, 12, f"a white dot appears on the chosen residue ({before} -> {with_marker})")
        p.page.click("#cftr-walk-clear")
        p.settle(250)
        self.assertLess(pngpixels.white_count(p.shot()) - before, 6, "and it goes when the marker is cleared")

    def test_every_frame_draws_the_connectors_for_unplaced_stretches_and_depth_tests_the_structure(self):
        p = self.open()
        p.wait_drawn()
        p.hold()
        d0 = json.loads(p.page.get_attribute("#cftr-canvas", "data-draws"))
        self.assertTrue(d0["depthOn"], "near parts of a structure hide far parts of it")
        self.assertGreaterEqual(d0["dots"], 1, "the dotted connectors across stretches the model does not place are drawn")
        self.assertGreaterEqual(d0["tube"], d0["dots"])
        for i in range(2, 5):
            p.page.click(f"#cftr-picks button:nth-child({i})")
            p.page.wait_for_timeout(700)
            d = json.loads(p.page.get_attribute("#cftr-canvas", "data-draws"))
            self.assertGreater(d["dots"], d0["dots"], "each structure has stretches it does not place, so each draws its connectors")
            d0 = d

    def test_switching_structure_says_it_is_a_cross_fade_and_the_note_goes_away(self):
        p = self.open()
        p.wait_drawn()
        p.hold()
        if p.page.query_selector("#cftr-note") is None:
            self.skipTest("this page has no note element")
        p.page.click("#cftr-picks button:nth-child(3)")
        self.assertEqual(norm(p.page.inner_text("#cftr-note")), "Cross-fade between two separate structures, 6MSM and 8EIQ, not a movement.")
        p.wait_until(lambda: norm(p.page.inner_text("#cftr-note")) == "", "the note to go", 5)
        quiet = self.open(reduced_motion="reduce")
        quiet.wait_drawn()
        quiet.page.click("#cftr-picks button:nth-child(2)")
        self.assertEqual(norm(quiet.page.inner_text("#cftr-note")), "", "with reduced motion there is no cross-fade and so no note")

    def test_the_wheel_only_zooms_when_the_picture_has_focus(self):
        p = self.open()
        p.wait_drawn()
        p.hold()
        before = p.look()
        box = p.page.locator("#cftr-canvas").bounding_box()
        p.page.mouse.move(box["x"] + 50, box["y"] + 50)
        scroll0 = p.page.evaluate("window.scrollY")
        p.page.mouse.wheel(0, 300)
        p.page.wait_for_timeout(800)
        scrolled_to = p.page.evaluate("window.scrollY")                # read before look(), which scrolls back to its fixed place
        self.assertNotEqual(scrolled_to, scroll0, "an unfocused picture must let the page scroll")
        self.assertEqual(p.look(), before, "an unfocused picture must not capture the wheel")
        # focused: the wheel zooms the picture and the page stays put
        p.page.focus("#cftr-canvas")
        p.look()                                                       # puts the picture at its fixed place
        box = p.page.locator("#cftr-canvas").bounding_box()
        p.page.mouse.move(box["x"] + 50, box["y"] + 50)
        scroll1 = p.page.evaluate("window.scrollY")
        p.page.mouse.wheel(0, -300)
        p.page.wait_for_timeout(500)
        self.assertEqual(p.page.evaluate("window.scrollY"), scroll1, "a focused picture keeps the wheel for zooming")
        self.assertNotEqual(p.look(), before, "a focused picture zooms with the wheel")

    # ---------------------------------------------------------------- guided steps and captions
    def test_each_guided_step_sets_the_picture(self):
        p = self.open()
        p.wait_drawn()
        p.hold()                                                       # a step sets the view; the slow turn would keep redrawing under the check
        steps = p.page.query_selector_all("[data-step]")
        self.assertEqual(len(steps), 5)
        for li in steps:
            self.assertTrue(li.query_selector("button").is_visible(), "step buttons show once the script runs")
        expect = [("6MSM", False, None), ("6MSM", False, None), ("6MSM", False, 508), ("8EIQ", False, 508), ("8EIQ", True, None)]
        for li, (sid, overlay, walk) in zip(steps, expect):
            li.query_selector("button").click()
            p.settle(100)
            s = p.state()
            self.assertEqual((s["structure"], s["overlay"], s["walk"]), (sid, overlay, walk))
            self.assertEqual(li.get_attribute("aria-current"), "step")
            self.assertEqual(sum(1 for o in steps if o.get_attribute("aria-current") == "step"), 1)

    def test_the_numbers_in_the_captions_are_the_numbers_in_the_data(self):
        text = norm(self.open().page.inner_text("main"))
        t2 = norm(self.open().page.inner_text("[data-step]:nth-child(2) p"))
        self.assertIn(f"{DATA['length']:,} amino acids", text)
        for d in DATA["domains"]:
            if d["type"] == "Domain":
                self.assertIn(f"{d['start']} to {d['end']}", t2)
            if d["type"] == "Region" and d["name"] == "Disordered R region":
                self.assertIn(f"{d['start']} to {d['end']}", text)
        for sid in ("8EIQ", "5UAK", "8EJ1"):
            st = BY_ID[sid]
            self.assertIn(f"{st['superposition_rmsd_A']:.2f} \u00c5", text)
            self.assertIn(str(st["superposition_residues"]), text)
        counts = [s["residues_traced"] for s in DATA["structures"]]
        self.assertIn(f"between {min(counts):,} and {max(counts):,} of the {DATA['length']:,} positions in the UniProt sequence", text)
        # the lower-resolution remark must be about the structure that really has the lowest resolution
        worst = max(DATA["structures"], key=lambda s: float(s["resolution_A"]))
        self.assertEqual(worst["id"], "8EJ1")
        self.assertIn(f"8EJ1 was solved at lower resolution ({BY_ID['8EJ1']['resolution_A']} \u00c5)", text)
        # "none of that stretch is placed in any of these four models" must be true of the data
        r = next(d for d in DATA["domains"] if d["name"] == "Disordered R region")
        for st in DATA["structures"]:
            for seg in st["segments"]:
                span = range(seg["start"], seg["start"] + len(seg["xyz"]) // 3)
                self.assertFalse(set(span) & set(range(r["start"], r["end"] + 1)), f"{st['id']} places part of the R region")
        # "residue 508 ... inside the first nucleotide-binding part" must be true of the data
        nbd1 = [d for d in DATA["domains"] if d["type"] == "Domain" and "transporter 1" in d["name"]][0]
        self.assertTrue(nbd1["start"] <= 508 <= nbd1["end"])
        # "runs from residue 507 straight to 509" must be true of both F508del structures, and 508 present in the other two
        for sid in ("8EIQ", "8EJ1"):
            placed = {n for seg in BY_ID[sid]["segments"] for n in range(seg["start"], seg["start"] + len(seg["xyz"]) // 3)}
            self.assertTrue({507, 509} <= placed and 508 not in placed, sid)
        for sid in ("6MSM", "5UAK"):
            self.assertIsNotNone(BY_ID[sid]["res508"])

    def test_the_caption_about_the_engineered_change_matches_the_files(self):
        text = norm(self.open().page.inner_text("main"))
        self.assertIn("an engineered change at position 1371 (glutamate replaced by glutamine), as do the files for 6MSM and 8EJ1; 5UAK does not have it.", text)
        carrying = {s["id"] for s in DATA["structures"] if any(d["position"] == 1371 and d["details"] == "engineered mutation" and d["file_residue"] == "GLN" and d["uniprot_residue"] == "GLU" for d in s["differences"])}
        self.assertEqual(carrying, {"6MSM", "8EIQ", "8EJ1"})
        self.assertEqual(DATA["uniprot"]["sequence"][1370], "E")
        # the table says the same, per structure, in the file's words
        rows = {r.query_selector("th").inner_text(): norm(r.inner_text()) for r in self.open().page.query_selector_all(".cftr-table tbody tr")}
        for sid in ("6MSM", "8EIQ", "8EJ1"):
            self.assertIn("1371: glutamine in the file, glutamate in UniProt (the file says: engineered mutation)", rows[sid])
        self.assertNotIn("1371", rows["5UAK"].split("1437 to 1480")[-1])
        self.assertIn("9 expression-tag residues beyond the protein, not drawn", rows["5UAK"])
        self.assertIn("508: phenylalanine is not in the file", rows["8EIQ"])

    def test_the_page_makes_no_claims_it_was_told_to_leave_out(self):
        text = norm(self.open().page.inner_text("main")).lower()
        allowed = "they do not say why, they do not say that any structure is open or closed, and they do not show what the medicines do."
        self.assertIn(allowed, text)
        rest = text.replace(allowed, " ")
        for imperative in ("open the entry on uniprot.org.", "open the pubmed record.", "open the file header and read those fields."):          # how-to-check instructions, not claims
            rest = rest.replace(imperative, " ")
        for banned in ("cystic fibrosis changes", "open", "closed", "broken", "cure", "pore is", "gate is", "% of", "percent", "miracle", "reaches the cell surface in"):
            self.assertNotIn(banned, rest, banned)
        self.assertIn("does not show where they sit or how they act", text)
        self.assertIn("they do not say why", text)

    # ---------------------------------------------------------------- motion, phones, no script, context loss
    def test_reduced_motion_draws_one_still_picture_and_then_stops(self):
        p = self.open(reduced_motion="reduce")
        p.wait_drawn()
        self.assertFalse(p.state()["playing"])
        self.assertEqual(p.page.get_attribute("#cftr-play", "aria-pressed"), "false")
        self.assertGreater(pngpixels.nonblank(p.shot()), 0.01, "a visitor who asked for less motion must still see the protein")
        n = p.frames()
        p.page.wait_for_timeout(800)
        self.assertEqual(p.frames(), n, "nothing may animate or redraw by itself")
        p.page.click("#cftr-play")                                     # they can still choose to turn it
        p.page.wait_for_timeout(600)
        self.assertGreater(p.frames(), n + 3)
        p.page.click("#cftr-play")
        p.page.wait_for_timeout(300)
        n = p.frames()
        p.page.wait_for_timeout(500)
        self.assertEqual(p.frames(), n)

    def test_normal_motion_turns_slowly_and_pause_stops_it(self):
        p = self.open()
        p.wait_drawn()
        p.show()
        self.assertTrue(p.state()["playing"])
        a = p.frames()
        t0 = p.page.input_value("#cftr-turn")
        p.page.wait_for_timeout(700)
        self.assertGreater(p.frames(), a + 3)
        self.assertNotEqual(p.page.input_value("#cftr-turn"), t0)
        p.page.click("#cftr-play")
        p.page.wait_for_timeout(300)
        b = p.frames()
        p.page.wait_for_timeout(500)
        self.assertEqual(p.frames(), b)

    def test_the_published_state_follows_the_pause_button(self):
        p = self.open()
        p.wait_drawn()
        self.assertTrue(p.state()["playing"])
        p.page.click("#cftr-play")
        self.assertFalse(p.state()["playing"])
        p.page.click("#cftr-play")
        self.assertTrue(p.state()["playing"])

    def test_a_hidden_tab_does_not_keep_drawing(self):
        p = self.open()
        p.wait_drawn()
        p.show()
        p.page.evaluate("() => { Object.defineProperty(document, 'hidden', {configurable: true, get: () => true}); document.dispatchEvent(new Event('visibilitychange')); }")
        p.page.wait_for_timeout(300)
        n = p.frames()
        p.page.wait_for_timeout(500)
        self.assertEqual(p.frames(), n)
        p.page.evaluate("() => { Object.defineProperty(document, 'hidden', {configurable: true, get: () => false}); document.dispatchEvent(new Event('visibilitychange')); }")
        p.page.wait_for_timeout(500)
        self.assertGreater(p.frames(), n)

    def test_an_off_screen_picture_stops_drawing_and_resumes(self):
        p = self.open()
        p.wait_drawn()
        self.assertTrue(p.state()["playing"], "the slow turn is on here, so the picture is redrawing while it is seen")
        p.show()
        self.assertIs(p.state()["offscreen"], False)
        p.page.evaluate("() => window.scrollTo({top: document.documentElement.scrollHeight, behavior: 'instant'})")
        p.page.wait_for_timeout(300)
        self.assertIs(p.state()["offscreen"], True, "scrolled past: the picture is off the screen")
        n = p.frames()
        p.page.wait_for_timeout(600)
        self.assertEqual(p.frames(), n, "nothing is drawn while the picture is off the screen")
        # a third of the picture in view counts as seen
        seen = p.page.evaluate("""() => { const c = document.getElementById('cftr-canvas'); const vis = () => { const r = c.getBoundingClientRect();
            return (Math.min(r.bottom, innerHeight) - Math.max(r.top, 0)) / r.height; };
          for (let y = scrollY; y >= 0; y -= 10) { window.scrollTo({top: y, behavior: 'instant'}); const v = vis(); if (v > 0.2 && v < 0.5) return v; }
          return -1; }""")
        self.assertGreater(seen, 0.2, "found a scroll position with part of the picture in view")
        p.page.wait_for_timeout(400)
        self.assertIs(p.state()["offscreen"], False, f"{seen:.2f} of the picture is in view, so it is on the screen")
        p.show()
        p.page.wait_for_timeout(200)
        self.assertIs(p.state()["offscreen"], False)
        self.assertGreater(p.frames(), n, "back on the screen, it draws again and the turn resumes")

    def test_a_page_can_start_with_the_slow_turn_off(self):
        ctx = self.browser.new_context(viewport={"width": 1100, "height": 800})
        self.addCleanup(ctx.close)
        page = ctx.new_page()
        marked = []

        def turn_off(route):
            if route.request.resource_type != "document":
                return route.continue_()
            r = route.fetch()
            body = r.text()
            if "text/html" in r.headers.get("content-type", "") and 'id="cftr-viewer"' in body:
                body = body.replace('id="cftr-viewer"', 'id="cftr-viewer" data-turn="off"', 1)
                marked.append(route.request.url)
            route.fulfill(response=r, body=body)
        page.route("**/*", turn_off)
        page.goto(self.base + PAGE)
        self.assertTrue(marked, "the page was served with data-turn=off")
        p = Page.__new__(Page)
        p.page = page
        p.wait_drawn()
        p.show()
        self.assertFalse(p.state()["playing"])
        self.assertEqual(page.get_attribute("#cftr-play", "aria-pressed"), "false")
        page.wait_for_timeout(300)
        n = p.frames()
        page.wait_for_timeout(600)
        self.assertEqual(p.frames(), n, "with the turn off, nothing redraws by itself")
        page.click("#cftr-play")
        self.assertTrue(p.state()["playing"])
        page.wait_for_timeout(600)
        self.assertGreater(p.frames(), n + 3)

    def test_each_arrow_key_moves_the_view_the_right_way(self):
        p = self.open()
        p.wait_drawn()
        p.hold()
        p.page.focus("#cftr-canvas")
        p.page.keyboard.press("r")
        t0, f0 = float(p.page.input_value("#cftr-turn")), float(p.page.input_value("#cftr-tilt"))
        turn = lambda: (float(p.page.input_value("#cftr-turn")) - t0) % 360
        tilt = lambda: float(p.page.input_value("#cftr-tilt")) - f0
        p.page.keyboard.press("ArrowRight")
        self.assertTrue(4 <= turn() <= 10, f"Right should turn by about 7 degrees, got {turn()}")
        p.page.keyboard.press("ArrowLeft")
        p.page.keyboard.press("ArrowLeft")
        self.assertTrue(turn() >= 350, f"Left twice should end up about 7 degrees the other way, got {turn()}")
        p.page.keyboard.press("r")
        p.page.keyboard.press("ArrowDown")
        self.assertGreater(tilt(), 0.05)
        p.page.keyboard.press("ArrowUp")
        p.page.keyboard.press("ArrowUp")
        self.assertLess(tilt(), -0.05)

    def test_the_overlay_really_draws_6msm_behind_the_other_structure(self):
        p = self.open()
        p.wait_drawn()
        p.hold()
        p.page.click("#cftr-picks button:nth-child(3)")
        p.settle()
        plain = p.look()
        p.page.click("#cftr-overlay")
        p.settle()
        with_ghost = p.look()
        self.assertNotEqual(plain, with_ghost, "the ghost of 6MSM must change the picture")
        self.assertGreater(with_ghost.changed_pixels(plain), 400)
        p.page.click("#cftr-overlay")
        p.settle()
        self.assertEqual(p.look(), plain, "switching the overlay off must give the plain picture back")

    def test_a_hidden_tab_draws_nothing_even_when_the_window_changes_size(self):
        p = self.open()
        p.wait_drawn()
        p.show()
        p.page.evaluate("() => { Object.defineProperty(document, 'hidden', {configurable: true, get: () => true}); document.dispatchEvent(new Event('visibilitychange')); }")
        p.page.wait_for_timeout(300)
        n = p.frames()
        p.page.set_viewport_size({"width": 1000, "height": 760})           # a resize asks for a redraw; a hidden page must not do one
        p.page.wait_for_timeout(500)
        self.assertEqual(p.frames(), n)
        p.page.evaluate("() => { Object.defineProperty(document, 'hidden', {configurable: true, get: () => false}); document.dispatchEvent(new Event('visibilitychange')); }")
        p.page.wait_for_timeout(500)
        self.assertGreater(p.frames(), n)

    def test_a_lost_webgl_context_is_rebuilt(self):
        p = self.open()
        p.wait_drawn()
        p.hold()
        p.page.evaluate("() => { const gl = document.getElementById('cftr-canvas').getContext('webgl2'); window.__ext = gl.getExtension('WEBGL_lose_context'); window.__ext.loseContext(); }")
        p.page.wait_for_timeout(300)
        n = p.frames()
        p.page.evaluate("() => window.__ext.restoreContext()")
        p.wait_until(lambda: p.frames() > n, "a frame after the context was restored", 10)
        p.page.wait_for_timeout(200)
        self.assertGreater(pngpixels.nonblank(p.shot()), 0.01)
        self.assertEqual([m for m in p.problems if "pageerror" in m], [])

    def test_a_phone_needs_no_sideways_scroll_and_has_big_enough_controls(self):
        p = self.open(viewport={"width": 360, "height": 740}, has_touch=True, is_mobile=True, device_scale_factor=2)
        p.wait_drawn()
        self.assertEqual(p.problems, [])
        self.assertLessEqual(p.page.evaluate("document.documentElement.scrollWidth"), 360)
        box = p.page.locator("#cftr-canvas").bounding_box()
        self.assertLessEqual(box["x"] + box["width"], 360 + 1)
        self.assertGreater(box["height"], 300)
        small = []
        for el in p.page.query_selector_all(".cftr-btn:not([hidden])"):
            b = el.bounding_box()
            if b and (b["height"] < 43.5 or b["width"] < 43.5):
                small.append((el.inner_text(), b))
        self.assertEqual(small, [])
        # the table scrolls inside its own box, not the page
        wraps = p.page.locator(".cftr-tablewrap")
        self.assertGreaterEqual(wraps.count(), 2)
        for i in range(wraps.count()):
            wrap = wraps.nth(i)
            self.assertTrue(wrap.evaluate("el => el.scrollWidth > el.clientWidth"), f"table {i} scrolls inside its own box on a phone")
            self.assertEqual(wrap.get_attribute("tabindex"), "0")

    def test_the_picture_is_pinned_below_the_sticky_header_not_under_it(self):
        p = self.open()
        p.wait_drawn()
        if p.page.query_selector('.site-header') is None:
            self.skipTest('this page has no sticky header to stay clear of')
        header = p.page.evaluate("() => document.querySelector('.site-header').getBoundingClientRect().height")
        pinned = p.page.evaluate("() => parseFloat(getComputedStyle(document.querySelector('.cftr-stage')).top)")
        self.assertEqual(p.page.evaluate("() => getComputedStyle(document.querySelector('.cftr-stage')).position"), "sticky")
        self.assertGreaterEqual(pinned, header, f"the picture pins at {pinned}px but the header reaches {header}px")
        # and it really sits clear of the header while it is pinned
        p.page.evaluate("() => window.scrollTo({top: 700, behavior: 'instant'})")
        p.page.wait_for_timeout(200)
        top = p.page.evaluate("() => document.querySelector('.cftr-stage').getBoundingClientRect().top")
        header_bottom = p.page.evaluate("() => document.querySelector('.site-header').getBoundingClientRect().bottom")
        self.assertGreaterEqual(top, header_bottom - 0.5, "while pinned, the picture starts below the header's bottom edge")

    def test_without_javascript_every_fact_is_still_on_the_page(self):
        p = self.open(java_script_enabled=False)
        self.assertTrue(p.page.is_visible("#cftr-fallback"))
        self.assertTrue(p.page.is_hidden("#cftr-viewer"))
        text = norm(p.page.inner_text("main"))
        for sid, st in BY_ID.items():
            self.assertIn(sid, text)
            self.assertIn(f"{st['residues_traced']} of {DATA['length']}", text)
        for st in DATA["structures"]:
            self.assertIn(f"PubMed {st['pubmed']}", text)
        self.assertIn("1 to 4", text)
        self.assertIn("1452 to 1480", text)
        self.assertEqual(len(p.page.query_selector_all("[data-step] button:visible")), 0, "step buttons stay hidden without a script")

    def test_if_the_data_cannot_load_the_page_says_so_and_keeps_the_text(self):
        ctx = self.browser.new_context(viewport={"width": 1100, "height": 800})
        self.addCleanup(ctx.close)
        page = ctx.new_page()
        page.route("**/cftr-data.js", lambda route: route.abort())
        page.goto(self.base + PAGE)
        page.wait_for_timeout(800)
        self.assertTrue(page.is_visible("#cftr-fallback"))
        self.assertTrue(page.is_hidden("#cftr-viewer"))
        self.assertIn("1452 to 1480", norm(page.inner_text("main")))

    def test_the_fallback_names_the_reason_when_webgl2_is_missing(self):
        ctx = self.browser.new_context(viewport={"width": 1100, "height": 800})
        self.addCleanup(ctx.close)
        page = ctx.new_page()
        page.add_init_script("const g = HTMLCanvasElement.prototype.getContext; HTMLCanvasElement.prototype.getContext = function (t, o) { return t === 'webgl2' ? null : g.call(this, t, o); };")
        page.goto(self.base + PAGE)
        page.wait_for_timeout(800)
        self.assertTrue(page.is_visible("#cftr-fallback"))
        self.assertIn("WebGL2", page.inner_text("#cftr-fallback"))
        self.assertEqual(page.inner_text("#cftr-fallback [data-why]"), "Your browser could not start WebGL2.")
        self.assertTrue(page.is_hidden("#cftr-viewer"))

    def test_focus_is_visible_on_the_picture_and_the_controls_have_names(self):
        p = self.open()
        p.wait_drawn()
        p.page.focus("#cftr-canvas")
        outline = p.page.evaluate("getComputedStyle(document.getElementById('cftr-canvas')).outlineStyle")
        self.assertNotEqual(outline, "none")
        for el in p.page.query_selector_all("#cftr-viewer button, #cftr-viewer input"):
            name = p.page.evaluate("""el => (el.getAttribute('aria-label') || el.innerText || (el.labels && el.labels[0] && el.labels[0].innerText) || '').trim()""", el)
            self.assertTrue(name, f"unnamed control {el.get_attribute('id')}")
        self.assertTrue(len(p.page.get_attribute("#cftr-canvas", "aria-label")) > 60)


if __name__ == "__main__":
    unittest.main()
