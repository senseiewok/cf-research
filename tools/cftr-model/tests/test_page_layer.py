"""The page layer: web/cftr-page.js (controls panel, touch gestures, variant box), web/cftr-variant.js (a pure parser for protein-level variant names) and
web/cftr-page.css. The lab's website uses them on its model page; the standalone index.html here does not load them yet, so the panel and touch behaviour
are tested on the website, where the markup for them is.

Two parts:
  StaticTests   standard library only: the files are here and listed, make no network call and use no storage, and the parser is exported.
  ParserTests   the parser's accepted and refused forms, run in a real Chromium from the module itself (needs Playwright, like test_viewer).

Run from this folder:   python -m unittest test_page_layer -v"""
import re
import threading
import unittest
from pathlib import Path

PKG = Path(__file__).resolve().parents[1]
FILES = ("web/cftr-page.js", "web/cftr-variant.js", "web/cftr-page.css")

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
    def test_the_three_files_are_in_the_package(self):
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
        for rel in ("web/cftr-page.js", "web/cftr-variant.js"):
            with self.subTest(rel=rel):
                self.assertEqual(problems(code_only(text(rel)), FORBIDDEN), [])

    def test_the_stylesheet_loads_nothing(self):
        self.assertEqual(problems(code_only(text("web/cftr-page.css")), CSS_FORBIDDEN), [])

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

    def test_the_page_script_imports_only_its_neighbours(self):
        imports = re.findall(r"(?m)^import .* from '([^']+)';$", text("web/cftr-page.js"))
        self.assertEqual(imports, ["./cftr-viewer.js", "./cftr-variant.js"])
        for rel in imports:
            self.assertTrue((PKG / "web" / rel).is_file(), rel)

    def test_the_walk_by_touch_is_built_by_the_script(self):
        """The step buttons, the bubble and the strip gestures are made here, with createElement and textContent, only when the walk slider exists;
        the website's Playwright tests (CftrPanelTests) check how they behave on a phone."""
        src = code_only(text("web/cftr-page.js"))
        start = "const walk = $('cftr-walk')"
        self.assertIn(start, src, "the walk block is in the page script")
        block = src.split(start, 1)[1].split("// ---- the variant box", 1)[0]
        self.assertIn("if (walk) {", block, "nothing runs without the slider")
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
        css = code_only(text("web/cftr-page.css"))
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


if __name__ == "__main__":
    unittest.main()
