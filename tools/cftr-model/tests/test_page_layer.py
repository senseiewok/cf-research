"""The page layer: web/cftr-page.js (the model page: controls panel, variant box), web/cftr-variant.js (a pure parser for protein-level variant names),
web/cftr-controls.js (shared by both pages: the sheet, touch gestures, the walk by touch), web/cftr-home.js (the home page hero: a Controls button and
sheet, pinch and double tap, no tilt by touch) and their stylesheets. The lab's website uses them; the standalone index.html here does not load them yet,
so the sheets are tested on the website, where the markup for them is.

Three parts:
  StaticTests          standard library only: the files are here and listed, make no network call and use no storage, and the parser is exported.
  ParserTests          the parser's accepted and refused forms, run in a real Chromium from the module itself (needs Playwright, like test_viewer).
  SharedControlsTests  the page scripts loaded on the standalone page, which has none of their panel or variant markup: no error, and the shared
                       touch and walk code works there (needs Playwright).

Run from this folder:   python -m unittest test_page_layer -v"""
import re
import threading
import unittest
from pathlib import Path

PKG = Path(__file__).resolve().parents[1]
FILES = ("web/cftr-page.js", "web/cftr-variant.js", "web/cftr-page.css", "web/cftr-controls.js", "web/cftr-controls.css", "web/cftr-home.js", "web/cftr-home.css")
SCRIPTS = tuple(f for f in FILES if f.endswith(".js"))
STYLES = tuple(f for f in FILES if f.endswith(".css"))

# Anything that could send or keep what a visitor types. A dynamic import() of any kind counts, since it can load from a URL.
FORBIDDEN = {
    "fetch": r"\bfetch\s*\(",
    "XMLHttpRequest": r"\bXMLHttpRequest\b",
    "WebSocket": r"\bWebSocket\b",
    "EventSource": r"\bEventSource\b",
    "sendBeacon": r"\bsendBeacon\b",
    "dynamic import": r"\bimport\s*\(",
    "a URL": r"\b(?:https?|wss?):",
    "localStorage": r"\blocalStorage\b",
    "sessionStorage": r"\bsessionStorage\b",
    "indexedDB": r"\bindexedDB\b",
    "cookie": r"\bdocument\.cookie\b",
    "history": r"\bhistory\.(?:push|replace)State\b",
}
CSS_FORBIDDEN = {"url()": r"\burl\s*\(", "@import": r"@import\b"}


def text(rel):
    return (PKG / rel).read_bytes().decode("utf-8").replace("\r\n", "\n")


def code_only(src):
    """The script without its comments, so a comment that names what the file never does is not mistaken for doing it."""
    src = re.sub(r"/\*.*?\*/", "", src, flags=re.S)
    return re.sub(r"(?m)(^|[^:\\])//.*$", r"\1", src)


def problems(src, rules):
    return [name for name, pattern in rules.items() if re.search(pattern, src)]


class StaticTests(unittest.TestCase):
    def test_the_files_are_in_the_package(self):
        for rel in FILES:
            with self.subTest(rel=rel):
                self.assertTrue((PKG / rel).is_file(), rel)
                self.assertGreater((PKG / rel).stat().st_size, 0, rel)

    def test_the_readme_lists_each_file(self):
        readme = text("README.md")
        files_table = readme.split("## Files", 1)[1].split("\n## ", 1)[0]
        for rel in FILES:
            with self.subTest(rel=rel):
                self.assertIn(f"`{rel}`", files_table)

    def test_the_readme_carries_the_limit_and_the_cftr2_line(self):
        readme = text("README.md")
        self.assertIn("This shows where the position sits on one structure. It does not say what the change does or what it means for anyone.", readme)
        self.assertRegex(readme, r"CFTR2[^\n]*not included")

    def test_the_scripts_make_no_network_call_and_use_no_storage(self):
        for rel in SCRIPTS:
            with self.subTest(rel=rel):
                self.assertEqual(problems(code_only(text(rel)), FORBIDDEN), [])

    def test_the_stylesheets_load_nothing(self):
        for rel in STYLES:
            with self.subTest(rel=rel):
                self.assertEqual(problems(code_only(text(rel)), CSS_FORBIDDEN), [])

    def test_the_detector_would_catch_each_kind_of_call(self):
        """A check that can never fail proves nothing: each rule must fire on a line that does what it names."""
        samples = {"fetch": "fetch('/x')", "XMLHttpRequest": "new XMLHttpRequest()", "WebSocket": "new WebSocket(u)", "EventSource": "new EventSource(u)",
                   "sendBeacon": "navigator.sendBeacon(u, d)", "dynamic import": "await import('./x.js')", "a URL": "const u = 'https://example.org/'",
                   "localStorage": "localStorage.setItem('a', 1)", "sessionStorage": "sessionStorage.a = 1", "indexedDB": "indexedDB.open('a')",
                   "cookie": "document.cookie = 'a=1'", "history": "history.replaceState(null, '', u)"}
        self.assertEqual(set(samples), set(FORBIDDEN))
        for name, line in samples.items():
            with self.subTest(name=name):
                self.assertIn(name, problems(code_only(text("web/cftr-variant.js") + "\n" + line + "\n"), FORBIDDEN))
        self.assertEqual(problems("a { background: url(x.png); }", CSS_FORBIDDEN), ["url()"])

    def test_the_parser_is_exported_and_pure(self):
        src = code_only(text("web/cftr-variant.js"))
        self.assertRegex(src, r"(?m)^export function parseVariant\(text, max = 1480\)")
        self.assertNotRegex(src, r"\bimport\b", "the parser imports nothing")
        self.assertNotRegex(src, r"\b(?:document|window|navigator|globalThis)\b", "the parser touches no page")

    def test_the_page_scripts_import_only_their_neighbours(self):
        expected = {"web/cftr-page.js": ["./cftr-viewer.js", "./cftr-variant.js", "./cftr-controls.js"], "web/cftr-home.js": ["./cftr-viewer.js", "./cftr-controls.js"],
                    "web/cftr-controls.js": [], "web/cftr-variant.js": []}
        for rel, want in expected.items():
            with self.subTest(rel=rel):
                imports = re.findall(r"(?m)^import .* from '([^']+)';$", text(rel))
                self.assertEqual(imports, want)
                self.assertEqual(len(re.findall(r"(?m)^\s*import\b", text(rel))), len(want), "every import is on one line of its own")
                for dep in imports:
                    self.assertTrue((PKG / "web" / dep).is_file(), dep)

    def test_the_shared_controls_stop_quietly_without_their_elements(self):
        """A page that lacks the panel, the sheet buttons or the walk slider gets nothing built for them and no error; the browser test below loads them."""
        src = code_only(text("web/cftr-controls.js"))
        self.assertIn("if (!root || !panel || !open) return null;", src)
        self.assertIn("if (!walk) return null;", src)
        self.assertIn("close?.addEventListener", src)
        self.assertIn("title?.focus", src)
        self.assertIn("$('cftr-reset')?.click()", src)
        home = code_only(text("web/cftr-home.js"))
        self.assertIn("if (api && root && panel) {", home)
        self.assertIn("touch(api, { tilt: false })", home, "on the home page one finger only turns; a vertical swipe scrolls the page")
        page = code_only(text("web/cftr-page.js"))
        self.assertIn("form?.addEventListener('submit'", page)
        self.assertNotRegex(page, r"\$\('cftr-(?:panel|sheet-open|sheet-close|panel-title)'\)\.", "the model page script never dereferences the panel's elements itself")

    def test_the_walk_by_touch_is_built_by_the_script(self):
        """The step buttons, the bubble and the strip gestures are made here, with createElement and textContent, only when the walk slider exists;
        the website's Playwright tests (CftrPanelTests, HomeTouchTests) check how they behave on a phone."""
        src = code_only(text("web/cftr-controls.js"))
        start = "export function walkControls("
        self.assertIn(start, src, "the walk block is in the shared script")
        block = src.split(start, 1)[1]
        self.assertIn("if (!walk) return null;", block, "nothing runs without the slider")
        self.assertIn("for (const d of [-10, -1, 1, 10])", block)
        self.assertIn("mk('button', 'cftr-btn')", block)
        self.assertIn("type: 'button'", block)
        self.assertIn("ariaLabel: (d < 0 ? 'Back ' : 'Forward ')", block)
        self.assertIn("document.createElement(t)", src)
        self.assertNotRegex(src, r"\.(?:innerHTML|outerHTML)\b|insertAdjacentHTML", "no markup from strings")
        for event in ("pointerdown", "pointerup", "pointercancel", "pointerleave", "blur", "visibilitychange"):
            self.assertIn(event, block, f"the hold stops on {event}")
        self.assertIn("Math.max(1, Math.min(length, n))", block, "clamped to the chain")
        self.assertIn("new Event('input', { bubbles: true })", block, "the viewer hears the change as it hears the slider")
        self.assertIn("strip.setPointerCapture(e.pointerId)", block)
        self.assertIn("e.isPrimary", block)

    def test_the_stylesheet_makes_room_for_a_finger(self):
        css = code_only(text("web/cftr-controls.css"))
        self.assertIn("@media (pointer: coarse) {", css)
        coarse = css.split("@media (pointer: coarse) {", 1)[1].split("\n}", 1)[0]
        for thumb in ("::-webkit-slider-thumb", "::-moz-range-thumb"):
            with self.subTest(thumb=thumb):
                m = re.search(re.escape(".cftr-ruler-input" + thumb) + r"\s*\{([^}]*)\}", coarse)
                self.assertTrue(m, thumb)
                width = re.search(r"width:\s*([\d.]+)rem", m.group(1))
                self.assertGreaterEqual(float(width.group(1)) * 16, 28, "at least 28 px wide on a touch screen")
                self.assertRegex(m.group(1), r"border:\s*2px solid")
        self.assertRegex(css, r"#cftr-seq\s*\{[^}]*touch-action:\s*none", "the strip takes its own drag")
        self.assertRegex(css, r"\.cftr-nudge\s*\{[^}]*grid-template-columns:\s*repeat\(4, 1fr\)")
        self.assertNotIn(".cftr-steps", css, "the guided look already uses .cftr-steps; a second rule for it here broke the phone layout once")

    def test_every_reason_the_parser_gives_has_a_message(self):
        reasons = set(re.findall(r"no\('(\w+)'\)", text("web/cftr-variant.js")))
        self.assertEqual(reasons, {"empty", "range", "dna", "many", "fs", "unread"})
        block = text("web/cftr-page.js").split("export const MESSAGES = {", 1)[1].split("\n};", 1)[0]
        self.assertEqual(set(re.findall(r"(?m)^\s*(\w+):", block)), reasons)


try:
    import support as V
    from support import sync_playwright
except Exception:  # pragma: no cover - the data module is missing or Playwright is not importable
    V, sync_playwright = None, None


def module_path():
    """Where the parser is served from: web/ here, scripts/ on the website (CFTR_WEB_ROOT)."""
    for rel in ("cftr-variant.js", "scripts/cftr-variant.js"):
        if (V.WEB / rel).is_file():
            return "/" + rel
    return None


@unittest.skipIf(sync_playwright is None, "Playwright is not installed in this Python")
class ParserTests(unittest.TestCase):
    """The same cases as the website's CftrPanelTests, run from the module in the page (node is not needed)."""

    @classmethod
    def setUpClass(cls):
        chrome = V.find_chrome()
        if not chrome:
            raise unittest.SkipTest("no Chromium found")
        cls.module = module_path()
        if not cls.module:
            raise unittest.SkipTest("cftr-variant.js is not in the folder being served")
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

    def open(self):
        p = V.Page(self.browser, self.base, reduced_motion="reduce")
        self.addCleanup(p.close)
        return p

    def assertClean(self, p):
        """No page error and no warning, and every file the page asked for is in the folder being served. The standalone page has no icon, so Chromium's
        own request for /favicon.ico (which the page does not make, and so is not in p.requests) may be the one 404."""
        missing = [m for m in p.problems if "status of 404" in m]
        self.assertEqual([m for m in p.problems if m not in missing], [])
        self.assertLessEqual(len(missing), 1, p.problems)
        for u in p.requests:
            rel = u[len(self.base):].split("?", 1)[0].split("#", 1)[0]
            self.assertTrue((V.WEB / (rel or V.PAGE or "index.html")).exists(), u)

    def parse(self, p, value):
        return p.page.evaluate("async ([m, t]) => (await import(m)).parseVariant(t)", [self.module, value])

    def test_the_parser_reads_each_accepted_form(self):
        p = self.open()
        cases = {"F508del": (508, 508, ["F"], "del"), "delF508": (508, 508, ["F"], "del"), "G551D": (551, 551, ["G"], "sub"), "R117H": (117, 117, ["R"], "sub"),
                 "W1282X": (1282, 1282, ["W"], "stop"), "W1282*": (1282, 1282, ["W"], "stop"), "W1282Ter": (1282, 1282, ["W"], "stop"),
                 "p.Phe508del": (508, 508, ["F"], "del"), "p.(Gly551Asp)": (551, 551, ["G"], "sub"), "p.Trp1282Ter": (1282, 1282, ["W"], "stop"),
                 "p.Ile507_Phe508del": (507, 508, ["I", "F"], "del"), "I507_F508del": (507, 508, ["I", "F"], "del"), "  f508del ": (508, 508, ["F"], "del"),
                 "508": (508, 508, [], "number"), "1": (1, 1, [], "number"), "1480": (1480, 1480, [], "number")}
        for value, (start, end, ref, kind) in cases.items():
            with self.subTest(text=value):
                r = self.parse(p, value)
                self.assertTrue(r["ok"], r)
                self.assertEqual((r["start"], r["end"], r["ref"], r["kind"]), (start, end, ref, kind))
        self.assertClean(p)

    def test_the_parser_refuses_what_it_cannot_place_and_says_why(self):
        p = self.open()
        cases = {"c.1521_1523delCTT": "dna", "621+1G>T": "dna", "3849+10kbC>T": "dna", "1717-1G>A": "dna", "IVS8-5T": "dna", "3659delC": "dna",
                 "p.Phe508fs": "fs", "S466fs*5": "fs", "F508del/G551D": "many", "F508del G551D": "many", "p.[Phe508del];[Gly551Asp]": "many",
                 "0": "range", "1481": "range", "G1500D": "range", "p.Phe508_Ile507del": "range", "": "empty", "   ": "empty",
                 "hello": "unread", "X508del": "unread", "F508": "unread", "p.Xyz508del": "unread", "B508D": "unread", "F508Q2": "unread"}
        for value, why in cases.items():
            with self.subTest(text=value):
                self.assertEqual(self.parse(p, value), {"ok": False, "why": why})
        self.assertClean(p)

    def test_the_parser_respects_a_shorter_length(self):
        p = self.open()
        r = p.page.evaluate("async ([m]) => { const { parseVariant } = await import(m); return [parseVariant('508', 500), parseVariant('500', 500)]; }", [self.module])
        self.assertEqual(r[0], {"ok": False, "why": "range"})
        self.assertTrue(r[1]["ok"])


def served(name):
    """Where a page-layer script is served from: web/ here, scripts/ on the website (CFTR_WEB_ROOT)."""
    for rel in (name, "scripts/" + name):
        if (V.WEB / rel).is_file():
            return "/" + rel
    return None


@unittest.skipIf(sync_playwright is None, "Playwright is not installed in this Python")
class SharedControlsTests(unittest.TestCase):
    """The page scripts on the standalone page, which has the viewer, the walk slider and the strip but none of the panel, sheet or variant markup.
    Before cftr-controls.js, cftr-page.js stopped with an error on such a page. Now each script builds what its elements allow and nothing else."""

    @classmethod
    def setUpClass(cls):
        chrome = V.find_chrome()
        if not chrome:
            raise unittest.SkipTest("no Chromium found")
        cls.page_js, cls.home_js = served("cftr-page.js"), served("cftr-home.js")
        if not (cls.page_js and cls.home_js) or V.PAGE:
            raise unittest.SkipTest("these tests load the standalone page of this package")
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

    def load(self, **ctx):
        p = V.Page(self.browser, self.base, reduced_motion="reduce", **ctx)
        self.addCleanup(p.close)
        p.wait_drawn()
        self.assertEqual(p.page.locator("#cftr-panel, #cftr-sheet-open, #cftr-sheet-close, #cftr-panel-title, #cftr-variant").count(), 0, "the standalone page has none of them")
        p.page.evaluate("async (ms) => { for (const m of ms) await import(m); }", [self.page_js, self.home_js])
        p.page.wait_for_timeout(200)
        return p

    def view(self, p):
        return p.page.evaluate("async () => (await import('/cftr-viewer.js')).ready.then(a => a.state())")

    def test_both_page_scripts_load_without_their_markup_and_build_only_what_fits(self):
        p = self.load()
        self.assertEqual([m for m in p.problems if "status of 404" not in m], [], "no error from either script")
        self.assertEqual(p.page.locator(".cftr-nudge").count(), 1, "the model page script builds the step buttons once, under the ruler")
        self.assertEqual(p.page.locator(".cftr-nudge button").count(), 4)
        self.assertEqual(p.page.locator(".cftr-home-walkout").count(), 0, "the home page script builds nothing without its sheet")
        p.page.get_by_role("button", name="Forward 10 residues", exact=True).click()
        p.page.get_by_role("button", name="Back 1 residue", exact=True).click()
        self.assertEqual(self.view(p)["walk"], 507, "with no marker the first step marks 508; one back is 507")

    def test_the_shared_touch_code_pinches_and_resets_on_a_phone(self):
        p = self.load(viewport={"width": 390, "height": 844}, has_touch=True, is_mobile=True, device_scale_factor=2)
        page = p.page
        page.evaluate("() => document.getElementById('cftr-canvas').scrollIntoView({ block: 'center', behavior: 'instant' })")
        page.wait_for_timeout(200)
        b = page.locator("#cftr-canvas").bounding_box()
        x, y = b["x"] + b["width"] / 2, b["y"] + b["height"] / 2
        cdp = p.ctx.new_cdp_session(page)
        send = lambda kind, pts=(): cdp.send("Input.dispatchTouchEvent", {"type": kind, "touchPoints": [{"x": a, "y": c, "id": i} for i, (a, c) in enumerate(pts)]})
        z0 = self.view(p)["zoom"]
        send("touchStart", [(x - 30, y), (x + 30, y)])
        for k in range(1, 9):
            send("touchMove", [(x - 30 - 8 * k, y), (x + 30 + 8 * k, y)])
        send("touchEnd")
        self.assertGreater(self.view(p)["zoom"], z0 * 1.5, "spreading two fingers zooms in")
        page.wait_for_timeout(400)
        for _ in range(2):
            send("touchStart", [(x, y)])
            send("touchEnd")
        v = self.view(p)
        self.assertEqual((round(v["theta"], 5), round(v["phi"], 5), v["zoom"]), (0.6, 1.35, 1), "a double tap is the Reset view button")
        self.assertFalse(v["viewTouched"])
        self.assertEqual([m for m in p.problems if "status of 404" not in m], [])

    # ---------------------------------------------------------------- the home page's layout: the sheet inside the picture's parent
    # Found on an iPhone: the home page's sheet sits inside #cftr-viewer, the picture's parent, where touch() listened in the capture phase, so every touch in
    # the sheet was stopped there and captured to the picture: it turned the picture, a double tap on a button reset the view, a step button never saw its
    # own pointerdown, and the sheet reopened where it had been scrolled to. These tests build that layout here from the standalone page's own controls.
    HOME = """async ([css, home]) => {
      const el = (t, props) => Object.assign(document.createElement(t), props), loaded = [];
      for (const href of css) { const l = el('link', { rel: 'stylesheet', href }); loaded.push(new Promise(r => { l.onload = r; })); document.head.append(l); }
      const viewer = document.getElementById('cftr-viewer'), stage = document.getElementById('cftr-canvas').parentElement, turn = document.getElementById('cftr-turn');
      viewer.classList.add('cftr-hero-viewer');
      const panel = el('div', { id: 'cftr-panel', className: 'cftr-home-panel' }), head = el('div', { className: 'cftr-panel-head' });
      head.append(el('h2', { id: 'cftr-panel-title', tabIndex: -1, textContent: 'Controls' }), el('button', { id: 'cftr-sheet-close', className: 'cftr-btn', type: 'button', textContent: 'Close' }));
      panel.append(head, ...turn.parentElement.children);          // Turn, Tilt, their labels, the button row and the hint, moved with the viewer's listeners
      stage.append(el('button', { id: 'cftr-sheet-open', className: 'cftr-btn cftr-home-open', type: 'button', textContent: 'Controls' }), panel);
      await Promise.all(loaded);
      await import(home);
    }"""
    SHEET = "() => { const s = document.getElementById('cftr-panel'); return { top: s.scrollTop, y: scrollY }; }"

    def home_layout(self):
        if not all((V.WEB / f).is_file() for f in ("cftr-home.css", "cftr-controls.css")):
            self.skipTest("the package's stylesheets are not in the folder being served")
        p = V.Page(self.browser, self.base, reduced_motion="reduce", viewport={"width": 390, "height": 844}, has_touch=True, is_mobile=True, device_scale_factor=2)
        self.addCleanup(p.close)
        p.wait_drawn()
        p.page.evaluate(self.HOME, [["/cftr-controls.css", "/cftr-home.css"], self.home_js])
        p.page.wait_for_timeout(300)
        cdp = p.ctx.new_cdp_session(p.page)
        send = lambda kind, pts=(): cdp.send("Input.dispatchTouchEvent", {"type": kind, "touchPoints": [{"x": a, "y": c, "id": i} for i, (a, c) in enumerate(pts)]})
        p.page.locator("#cftr-sheet-open").tap()
        p.page.wait_for_timeout(300)
        self.assertTrue(p.page.evaluate("document.getElementById('cftr-viewer').classList.contains('is-open')"), "the sheet is open")
        self.assertTrue(p.page.evaluate("document.getElementById('cftr-panel').parentElement === document.getElementById('cftr-canvas').parentElement"), "inside the picture's parent")
        return p, send

    def paced(self, page, send, start, moves):
        send("touchStart", [start])
        for pt in moves:
            send("touchMove", [pt])
            page.wait_for_timeout(16)
        send("touchEnd")
        page.wait_for_timeout(150)

    def centre(self, page, sel):
        return page.evaluate("(s) => { const r = document.querySelector(s).getBoundingClientRect(); return [r.left + r.width / 2, r.top + r.height / 2]; }", sel)

    def test_on_the_home_layout_a_finger_drags_and_taps_the_sheet_sliders(self):
        p, send = self.home_layout()
        page = p.page
        for sel, key, to_view in (("#cftr-turn", "theta", lambda v: v * 3.141592653589793 / 180), ("#cftr-tilt", "phi", float)):
            with self.subTest(slider=sel):
                page.evaluate("(s) => document.querySelector(s).scrollIntoView({ block: 'nearest', behavior: 'instant' })", sel)
                page.wait_for_timeout(100)
                v0, sheet0, value0 = self.view(p), page.evaluate(self.SHEET), float(page.input_value(sel))
                b = page.locator(sel).bounding_box()
                f = (value0 - float(page.get_attribute(sel, "min"))) / (float(page.get_attribute(sel, "max")) - float(page.get_attribute(sel, "min")))
                x, y = b["x"] + 12 + f * (b["width"] - 24), b["y"] + b["height"] / 2
                self.paced(page, send, (x, y), [(x + 150 * k / 8, y) for k in range(1, 9)])
                value1, v1 = float(page.input_value(sel)), self.view(p)
                self.assertGreater(value1, value0, f"a sideways drag moves {sel}")
                self.assertAlmostEqual(v1[key], to_view(value1), places=4, msg="and the model follows it")
                others = [k for k in ("theta", "phi", "zoom") if k != key]
                self.assertEqual([round(v1[k], 6) for k in others], [round(v0[k], 6) for k in others], "the drag does not also turn or zoom the picture")
                self.assertEqual(page.evaluate(self.SHEET), sheet0, "nor scroll the sheet or the page")
                x += 150
                self.paced(page, send, (x, y), [(x - 150 * k / 8, y + 12 * k / 8) for k in range(1, 9)])
                value2 = float(page.input_value(sel))
                self.assertLess(value2, value1, "a slightly diagonal drag still moves it")
                self.assertEqual(page.evaluate(self.SHEET), sheet0, "and scrolls neither the sheet nor the page")
                send("touchStart", [(b["x"] + 0.85 * b["width"], y)])
                send("touchEnd")
                page.wait_for_timeout(150)
                lo, hi = (float(page.get_attribute(sel, a)) for a in ("min", "max"))
                self.assertAlmostEqual((float(page.input_value(sel)) - lo) / (hi - lo), 0.85, delta=0.05, msg="a tap on the track puts the thumb there")
                self.assertAlmostEqual(self.view(p)[key], to_view(float(page.input_value(sel))), places=4)
                page.wait_for_timeout(400)
        self.assertEqual([m for m in p.problems if "status of 404" not in m], [])

    def test_on_the_home_layout_a_touch_on_the_sheet_never_moves_the_picture(self):
        p, send = self.home_layout()
        page = p.page
        page.get_by_role("button", name="Forward 1 residue", exact=True).tap()
        page.locator("#cftr-turn").focus()
        page.keyboard.press("ArrowRight")
        page.wait_for_timeout(400)
        v0 = self.view(p)
        bx, by = self.centre(page, "#cftr-panel .cftr-nudge button:nth-child(3)")
        self.assertEqual(page.evaluate("document.querySelector('#cftr-panel .cftr-nudge button:nth-child(3)').getAttribute('aria-label')"), "Forward 1 residue")
        for _ in range(2):
            send("touchStart", [(bx, by)])
            send("touchEnd")
        v1 = self.view(p)
        self.assertEqual(v1["walk"], v0["walk"] + 2, "each tap is a step")
        self.assertTrue(v1["viewTouched"], "two quick taps on a button are not the picture's double tap (Reset view)")
        self.assertEqual([round(v1[k], 6) for k in ("theta", "phi", "zoom")], [round(v0[k], 6) for k in ("theta", "phi", "zoom")])
        page.wait_for_timeout(400)
        send("touchStart", [(bx, by)])
        page.wait_for_timeout(1000)
        send("touchEnd")
        self.assertGreater(self.view(p)["walk"], v1["walk"] + 3, "a held step button repeats: it sees its own pointerdown")
        page.wait_for_timeout(400)
        c = page.locator("#cftr-canvas").bounding_box()
        cx, cy = c["x"] + c["width"] / 2, max(c["y"], 0) + 60
        self.assertEqual(page.evaluate("([x, y]) => document.elementFromPoint(x, y).id", [cx, cy]), "cftr-canvas")
        t0 = self.view(p)["theta"]
        self.paced(page, send, (cx - 60, cy), [(cx - 60 + 12 * k, cy) for k in range(1, 11)])
        self.assertGreater(self.view(p)["theta"], t0 + 0.5, "a sideways drag on the picture still turns it")
        z0 = self.view(p)["zoom"]
        send("touchStart", [(cx - 30, cy), (cx + 30, cy)])
        for k in range(1, 9):
            send("touchMove", [(cx - 30 - 8 * k, cy), (cx + 30 + 8 * k, cy)])
        send("touchEnd")
        self.assertGreater(self.view(p)["zoom"], z0 * 1.5, "a pinch on the picture still zooms")
        page.wait_for_timeout(400)
        for _ in range(2):
            send("touchStart", [(cx, cy)])
            send("touchEnd")
        self.assertFalse(self.view(p)["viewTouched"], "a double tap on the picture still resets")
        for sel in ("#cftr-turn", "#cftr-tilt", "#cftr-panel-title", "#cftr-panel .cftr-home-walkout"):   # last: the page itself may zoom
            with self.subTest(start=sel):
                before = self.view(p)["zoom"]
                x, y = self.centre(page, sel)
                send("touchStart", [(x - 30, y), (x + 30, y)])
                for k in range(1, 9):
                    send("touchMove", [(x - 30 - 8 * k, y), (x + 30 + 8 * k, y)])
                send("touchEnd")
                page.wait_for_timeout(400)
                self.assertEqual(self.view(p)["zoom"], before, f"two fingers from {sel} do not zoom the picture")
        self.assertEqual([m for m in p.problems if "status of 404" not in m], [])

    def test_on_the_home_layout_the_sheet_scrolls_and_reopens_at_its_top(self):
        p, send = self.home_layout()
        page = p.page
        self.assertEqual(page.evaluate(self.SHEET)["top"], 0, "the sheet opens at its top")
        self.assertGreater(page.evaluate("(() => { const s = document.getElementById('cftr-panel'); return s.scrollHeight - s.clientHeight; })()"), 50, "it has more than a screen")
        x, y = self.centre(page, "#cftr-panel .cftr-home-walkout")
        self.paced(page, send, (x, y), [(x, y - 12 * k) for k in range(1, 16)])
        page.wait_for_timeout(300)
        self.assertGreater(page.evaluate(self.SHEET)["top"], 30, "a vertical swipe on the sheet scrolls it")
        page.locator("#cftr-sheet-close").tap()
        page.wait_for_timeout(100)
        page.locator("#cftr-sheet-open").tap()
        page.wait_for_timeout(300)
        self.assertEqual(page.evaluate(self.SHEET)["top"], 0, "opened again, it starts at its top")
        self.assertEqual([m for m in p.problems if "status of 404" not in m], [])

if __name__ == "__main__":
    unittest.main()
