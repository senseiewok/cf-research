"""Shared helpers for the CFTR model's browser tests: which page to test, a local server, Chromium with a software GL, and a Page wrapper that records problems.

By default the tests run against the standalone page in ../web. To run the SAME tests against another page that has the same element ids (the lab's website does),
set these environment variables before starting Python:
    CFTR_WEB_ROOT   the folder to serve (default: ../web)
    CFTR_PAGE       the page to open, relative to that folder (default: "" = index.html)
    CFTR_DATA_JS    the data module the page imports (default: <web root>/cftr-data.js)
Playwright for Python is the only dependency, and only for the browser tests; it uses a Chromium already on disk (CHROME_PATH, or a Playwright browser folder)."""
import glob
import json
import os
import re
import sys
import time
import unittest
from pathlib import Path

try:
    from playwright.sync_api import sync_playwright
except ImportError:  # pragma: no cover
    sync_playwright = None

HERE = Path(__file__).resolve().parent
PKG = HERE.parent
sys.path.insert(0, str(PKG))
sys.path.insert(0, str(HERE))
import serve as model_server  # noqa: E402
import pngpixels  # noqa: E402

WEB = Path(os.environ.get("CFTR_WEB_ROOT") or PKG / "web")
SRC = WEB                                                           # the name the tests use for the folder being served
PAGE = os.environ.get("CFTR_PAGE", "")
DATA_JS = Path(os.environ.get("CFTR_DATA_JS") or WEB / "cftr-data.js")
DATA = json.loads(DATA_JS.read_text(encoding="utf-8").split("export default ", 1)[1].strip().rstrip(";"))
BY_ID = {s["id"]: s for s in DATA["structures"]}
ARGS = ["--use-angle=swiftshader", "--enable-unsafe-swiftshader", "--host-resolver-rules=MAP * ~NOTFOUND , EXCLUDE 127.0.0.1"]


def find_chrome():
    if os.environ.get("CHROME_PATH") and Path(os.environ["CHROME_PATH"]).is_file():
        return os.environ["CHROME_PATH"]
    for pattern in (Path.home() / "AppData/Local/ms-playwright/chromium-*/chrome-win64/chrome.exe", Path.home() / ".cache/ms-playwright/chromium-*/chrome-linux/chrome"):
        found = sorted(glob.glob(str(pattern)), reverse=True)
        if found:
            return found[0]
    return None


def norm(s):
    return re.sub(r"\s+", " ", s).strip()


class Look:
    """Decoded canvas pixels. Two looks are equal when at most a handful of pixels differ by more than 8 levels: a real change (a marker, a turn, a zoom)
    moves hundreds or thousands, and identical views measured here differ by 0."""
    MAX_STRAY = 20

    def __init__(self, png):
        self.w, self.h, self.px = pngpixels.pixels(png)

    def changed_pixels(self, other):
        if (self.w, self.h) != (other.w, other.h):
            return self.w * self.h
        return sum(1 for a, b in zip(self.px, other.px) if abs(a[0] - b[0]) > 8 or abs(a[1] - b[1]) > 8 or abs(a[2] - b[2]) > 8)

    def __eq__(self, other):
        return isinstance(other, Look) and self.changed_pixels(other) <= self.MAX_STRAY

    def __ne__(self, other):
        return not self == other

    __hash__ = None

    def __repr__(self):
        return f"<picture {self.w}x{self.h}>"


class Page:
    """One load of the page with its problems recorded."""

    def __init__(self, browser, base, path=PAGE, **ctx):
        self.ctx = browser.new_context(**ctx)
        self.page = self.ctx.new_page()
        self.problems, self.requests = [], []
        self.page.on("console", lambda m: self.problems.append(f"{m.type}: {m.text}") if m.type in ("error", "warning") else None)
        self.page.on("pageerror", lambda e: self.problems.append("pageerror: " + str(e)))
        self.page.on("requestfailed", lambda r: self.problems.append("requestfailed: " + r.url))
        self.page.on("request", lambda r: self.requests.append(r.url))
        self.page.goto(base + path)
        self.base = base

    def wait_until(self, condition, what, seconds=30):
        """Poll from here: Playwright's wait_for_function evaluates a string, which the page's Content-Security-Policy rightly refuses."""
        end = time.time() + seconds
        while time.time() < end:
            if condition():
                return
            self.page.wait_for_timeout(40)
        raise AssertionError("timed out waiting for " + what)

    def wait_drawn(self, minimum=1):
        self.wait_until(lambda: self.frames() >= minimum, f"{minimum} drawn frame(s)")

    def frames(self):
        return int(self.page.get_attribute("#cftr-canvas", "data-frames") or 0)

    def state(self):
        return json.loads(self.page.get_attribute("#cftr-canvas", "data-state"))

    def shot(self):
        """The canvas alone, always from the same scroll position so the pinned picture sits at the same sub-pixel place (and clear of the header)."""
        self.page.evaluate("() => { const v = document.getElementById('cftr-viewer').getBoundingClientRect(); window.scrollTo({top: window.scrollY + v.top + 250, behavior: 'instant'}); }")
        self.page.wait_for_timeout(150)
        return self.page.locator("#cftr-canvas").screenshot()

    def look(self):
        """The picture on the canvas right now, comparable with == and != (a few stray pixels of shading noise do not count as a change)."""
        return Look(self.shot())

    def show(self):
        """Bring the picture onto the screen: the viewer draws nothing while the canvas is off the screen, and on a page with a long introduction it starts below the fold."""
        self.page.evaluate("() => document.getElementById('cftr-canvas').scrollIntoView({block: 'center', behavior: 'instant'})")
        self.page.wait_for_timeout(300)

    def hold(self):
        """Pause the slow turn and let the last frame reach the screen, so a picture can be compared."""
        if self.state()["playing"]:
            self.page.click("#cftr-play")
        self.page.wait_for_timeout(300)

    def settle(self, ms=120):
        """Wait until drawing has gone quiet: the frame counter stops changing for one interval."""
        last = -1
        for _ in range(80):
            n = self.frames()
            if n == last:
                return
            last = n
            self.page.wait_for_timeout(ms)
        raise AssertionError("the picture never stopped redrawing")

    def close(self):
        self.ctx.close()
