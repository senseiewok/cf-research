"""Deterministic basics check for a one-file HTML page, for judging with-and-without trials of design guidance.

It decides only what a script can decide, in headless Chromium with all network access blocked and no dependency beyond the standard library:
  deps       the page loads nothing from outside itself (no external script, stylesheet, font, image, import or request)
  lang       <html lang> is set, <title> is not empty, there is exactly one <h1>, and there is a <main> landmark
  contrast   every visible text element has a measured contrast ratio of at least 4.5 against its background; text over a gradient or an image
             cannot be measured and fails (no evidence is no pass)
  mobile     at a 360 px wide window the page does not scroll sideways
  motion     no animation runs forever, and a page with any animation or transition has a prefers-reduced-motion rule
  copy       every sentence in the --copy file (one per line) appears word for word in the page text
  console    the page raised no uncaught error

It does not judge whether a page is well designed, calm, or good. A pass means the basics hold, nothing more.

Usage: python check_page_basics.py PAGE.html [--copy FILE] [--chrome PATH]
Prints one PASS or FAIL line per check, then `N/M checks passed`. Exit 0 when all pass, 1 when any fails, 2 for a missing page or browser."""
import argparse
import glob
import html
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from html.parser import HTMLParser
from pathlib import Path

FLAGS = ["--headless=new", "--disable-gpu", "--no-first-run", "--disable-extensions", "--disable-component-update", "--disable-background-networking",
         "--disable-sync", "--enable-logging=stderr", "--v=0", "--host-resolver-rules=MAP * ~NOTFOUND", "--virtual-time-budget=1500"]
FATAL_CONSOLE = re.compile(r"Uncaught|ReferenceError|TypeError|SyntaxError", re.I)

# Injected into a COPY of the page: measures what the page really renders, then publishes it as a data attribute on <html>.
MEASURE = r"""<script>addEventListener('load',function(){setTimeout(function(){try{
function rgba(s){var m=s.match(/rgba?\(([^)]+)\)/);if(!m)return null;var p=m[1].split(/[ ,\/]+/).filter(Boolean).map(parseFloat);return{r:p[0],g:p[1],b:p[2],a:p.length>3&&!isNaN(p[3])?p[3]:1}}
function lum(c){function f(v){v/=255;return v<=0.03928?v/12.92:Math.pow((v+0.055)/1.055,2.4)}return 0.2126*f(c.r)+0.7152*f(c.g)+0.0722*f(c.b)}
function over(top,bot){var a=top.a;return{r:top.r*a+bot.r*(1-a),g:top.g*a+bot.g*(1-a),b:top.b*a+bot.b*(1-a),a:1}}
function bgOf(el){var chain=[];for(var e=el;e;e=e.parentElement){var cs=getComputedStyle(e);if(cs.backgroundImage&&cs.backgroundImage!=='none')return null;var c=rgba(cs.backgroundColor);if(c&&c.a>0)chain.push(c);if(c&&c.a===1)break}
var base={r:255,g:255,b:255,a:1};for(var i=chain.length-1;i>=0;i--)base=over(chain[i],base);return base}
var min=null,minTag='',unmeasured=0,texts=0;
var all=document.body?document.body.querySelectorAll('*'):[];
for(var i=0;i<all.length;i++){var el=all[i];var own=false;for(var n=el.firstChild;n;n=n.nextSibling){if(n.nodeType===3&&n.nodeValue.trim())own=true}
if(!own)continue;var cs=getComputedStyle(el);if(cs.display==='none'||cs.visibility==='hidden'||parseFloat(cs.opacity)===0)continue;
if(['SCRIPT','STYLE','NOSCRIPT'].indexOf(el.tagName)>=0)continue;texts++;
var fg=rgba(cs.color);var bg=bgOf(el);if(!fg||!bg){unmeasured++;continue}
var f=over(fg,bg);var l1=lum(f),l2=lum(bg);var r=(Math.max(l1,l2)+0.05)/(Math.min(l1,l2)+0.05);
if(min===null||r<min){min=r;minTag=el.tagName.toLowerCase()}}
var anims=document.getAnimations?document.getAnimations():[];var inf=0;for(var j=0;j<anims.length;j++){try{if(anims[j].effect.getComputedTiming().iterations===Infinity)inf++}catch(e){}}
var h=document.documentElement;
h.setAttribute('data-m',JSON.stringify({lang:h.getAttribute('lang')||'',title:document.title||'',h1:document.querySelectorAll('h1').length,main:document.querySelectorAll('main,[role=main]').length,
innerWidth:innerWidth,overflow:h.scrollWidth-h.clientWidth,minContrast:min,minTag:minTag,unmeasured:unmeasured,texts:texts,animations:anims.length,infinite:inf}))}catch(e){document.documentElement.setAttribute('data-m',JSON.stringify({error:String(e)}))}},200)});</script>"""


def find_chrome(explicit=None):
    cands = [explicit, os.environ.get("CHROME_PATH")]
    home = Path.home()
    for pattern in (home / "AppData/Local/ms-playwright/chromium-*/chrome-win64/chrome.exe", home / ".cache/ms-playwright/chromium-*/chrome-linux/chrome",
                    home / "Library/Caches/ms-playwright/chromium-*/chrome-mac*/Chromium.app/Contents/MacOS/Chromium"):
        cands += sorted(glob.glob(str(pattern)), reverse=True)
    cands += [shutil.which(n) for n in ("chromium", "chromium-browser", "google-chrome", "chrome", "msedge")]
    cands += [r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe", r"C:\Program Files\Microsoft\Edge\Application\msedge.exe"]
    for c in cands:
        if c and Path(c).is_file():
            return str(c)
    return None


class _Refs(HTMLParser):
    """Everything the page would load from somewhere else."""

    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.refs, self.text, self._skip, self.h1 = [], [], 0, 0

    def handle_starttag(self, tag, attrs):
        a = {k: (v or "") for k, v in attrs}
        if tag in ("script", "style"):
            self._skip += 1
        for key, rel_ok in (("src", False), ("data", False), ("poster", False)):
            if key in a and a[key].strip() and not a[key].strip().lower().startswith("data:") and not (rel_ok and a[key].startswith("#")):
                if tag in ("script", "img", "source", "video", "audio", "iframe", "embed", "object", "track", "input"):
                    self.refs.append(f"<{tag} {key}={a[key][:60]}>")
        if tag == "link" and a.get("href", "").strip() and not a["href"].strip().lower().startswith("data:"):
            if (a.get("rel", "").lower() in ("stylesheet", "preload", "modulepreload", "prefetch", "icon", "shortcut icon", "manifest", "preconnect", "dns-prefetch")) or "rel" not in a:
                self.refs.append(f"<link rel={a.get('rel', '')} href={a['href'][:60]}>")
        for css in (a.get("style", ""),):
            for m in re.finditer(r"url\(\s*['\"]?\s*([^'\")\s]+)", css, re.I):
                if not m.group(1).lower().startswith("data:"):
                    self.refs.append(f"style url({m.group(1)[:60]})")

    def handle_endtag(self, tag):
        if tag in ("script", "style") and self._skip:
            self._skip -= 1

    def handle_data(self, data):
        if not self._skip:
            self.text.append(data)


def source_refs(text):
    p = _Refs()
    p.feed(text)
    refs = list(p.refs)
    # a url( inside an inline data: URI (for example an SVG filter reference such as url(%23n)) is part of that inline image, not a load
    text = re.sub(r"""url\(\s*(['"])\s*data:.*?\1\s*\)|url\(\s*data:[^)]*\)""", "url(data:)", text, flags=re.I | re.S)
    for block in re.findall(r"<style[^>]*>(.*?)</style>", text, re.I | re.S):
        for m in re.finditer(r"url\(\s*['\"]?\s*([^'\")\s]+)", block, re.I):
            if not m.group(1).lower().startswith("data:"):
                refs.append(f"css url({m.group(1)[:60]})")
        for m in re.finditer(r"@import\s+(?:url\()?\s*['\"]?([^'\")\s;]+)", block, re.I):
            refs.append(f"@import {m.group(1)[:60]}")
    for block in re.findall(r"<script(?![^>]*\bsrc=)[^>]*>(.*?)</script>", text, re.I | re.S):
        if re.search(r"\bfetch\s*\(|XMLHttpRequest|\bWebSocket\b|\bimport\s*\(\s*['\"]https?:|\bsendBeacon\b", block):
            refs.append("script makes a request")
    return refs, " ".join(html.unescape(" ".join(p.text)).split())


def inject(text):
    for pattern in (r"</head\s*>", r"</body\s*>"):
        m = re.search(pattern, text, re.I)
        if m:
            return text[:m.start()] + MEASURE + text[m.start():]
    return text + MEASURE


WRAPPER = """<!doctype html><html><head><meta charset="utf-8"></head><body style="margin:0"><iframe id="f" src="{src}" style="width:360px;height:640px;border:0"></iframe>
<script>document.getElementById('f').addEventListener('load',function(){setTimeout(function(){var d=document.getElementById('f').contentDocument;
document.documentElement.setAttribute('data-m',d?(d.documentElement.getAttribute('data-m')||''):'')},600)});</script></body></html>"""


def measure(chrome, page_text):
    """Measure the page inside a 360 px wide iframe. Headless Chromium will not make a window narrower than about 512 px (measured), so a
    window-size flag would measure the wrong width; an iframe gives the page a real 360 px viewport, media queries included."""
    base = Path(tempfile.gettempdir()) / f"check-page-basics-{os.getpid()}"
    page_copy, wrapper = Path(str(base) + "-page.html"), Path(str(base) + "-wrap.html")
    page_copy.write_text(inject(page_text), encoding="utf-8")
    wrapper.write_text(WRAPPER.replace("{src}", page_copy.as_uri()), encoding="utf-8")
    try:
        r = subprocess.run([chrome, *FLAGS, "--allow-file-access-from-files", "--dump-dom", wrapper.as_uri()], capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=90)
    finally:
        for f in (page_copy, wrapper):
            try:
                f.unlink()
            except OSError:
                pass
    m = re.search(r'<html[^>]*\sdata-m="([^"]*)"', r.stdout)
    data = {}
    if m:
        try:
            data = json.loads(html.unescape(m.group(1)))
        except ValueError:
            data = {}
    fatal = [re.sub(r"^\[[^\]]*\]\s*", "", ln.strip())[:140] for ln in r.stderr.splitlines() if "CONSOLE" in ln and FATAL_CONSOLE.search(ln)]
    return data, fatal


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("page")
    ap.add_argument("--copy", help="a file of sentences that must appear word for word, one per line")
    ap.add_argument("--chrome")
    a = ap.parse_args(argv)
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(errors="backslashreplace")
    page = Path(a.page)
    if not page.is_file():
        print(f"ERROR: no such page: {page}")
        return 2
    chrome = find_chrome(a.chrome)
    if not chrome:
        print("ERROR: no Chromium or Edge found; set CHROME_PATH or use --chrome")
        return 2
    text = page.read_text(encoding="utf-8", errors="replace")
    refs, plain = source_refs(text)
    m, fatal = measure(chrome, text)
    results = []

    def check(name, ok, detail=""):
        results.append((name, bool(ok), detail))

    check("deps", not refs, "loads from outside the page: " + "; ".join(refs[:3]) if refs else "")
    ok_struct = bool(m) and "error" not in m and m.get("lang", "") != "" and m.get("title", "").strip() != "" and m.get("h1") == 1 and (m.get("main") or 0) >= 1
    check("lang", ok_struct, f"lang={m.get('lang')!r} title={m.get('title')!r} h1={m.get('h1')} main={m.get('main')}" if m else "the page could not be measured")
    mc = m.get("minContrast") if m else None
    um = m.get("unmeasured", 0) if m else 0
    check("contrast", m and mc is not None and mc >= 4.5 and um == 0 and (m.get("texts") or 0) > 0,
          f"lowest measured ratio {mc if mc is None else round(mc, 2)} (on <{m.get('minTag')}>), {um} text element(s) over a gradient or image could not be measured" if m else "the page could not be measured")
    check("mobile", m and m.get("innerWidth") == 360 and (m.get("overflow") or 0) <= 1, f"window {m.get('innerWidth')} px wide, scrolls sideways by {m.get('overflow')} px" if m else "the page could not be measured")
    motion_src = bool(re.search(r"@keyframes|\banimation\s*:|\btransition\s*:|animation-name", text, re.I)) or (m and (m.get("animations") or 0) > 0)
    guard = "prefers-reduced-motion" in text
    check("motion", m and (m.get("infinite") or 0) == 0 and (not motion_src or guard),
          f"{m.get('infinite')} animation(s) run forever; motion present={bool(motion_src)}; prefers-reduced-motion rule present={guard}" if m else "the page could not be measured")
    if a.copy:
        wanted = [ln.strip() for ln in Path(a.copy).read_text(encoding="utf-8").splitlines() if ln.strip()]
        norm = " ".join(plain.split())
        missing = [w for w in wanted if " ".join(w.split()) not in norm]
        check("copy", not missing, f"{len(missing)} of {len(wanted)} sentence(s) not word for word: {missing[0][:70]!r}" if missing else "")
    check("console", not fatal, "; ".join(fatal[:2]) if fatal else "")
    for name, ok, detail in results:
        print(("PASS " if ok else "FAIL ") + name + ("" if ok or not detail else f": {detail}"))
    passed = sum(1 for _n, ok, _d in results if ok)
    print(f"{passed}/{len(results)} checks passed")
    return 0 if passed == len(results) else 1


if __name__ == "__main__":
    sys.exit(main())
