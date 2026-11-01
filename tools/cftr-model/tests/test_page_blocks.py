"""The generated parts of the standalone page match what the generator makes from the data. Run from this folder:   python -m unittest test_page_blocks -v"""
import subprocess
import sys
import unittest
from pathlib import Path

PKG = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PKG))
import make_page_blocks as gen  # noqa: E402


class PageBlockTests(unittest.TestCase):
    def test_the_standalone_page_is_up_to_date(self):
        r = subprocess.run([sys.executable, str(PKG / "make_page_blocks.py"), "--check"], capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)

    def test_a_page_that_differs_from_the_data_is_reported(self):
        import tempfile
        page = (PKG / "web" / "index.html").read_text(encoding="utf-8")
        data = gen.load_data()
        wrong = page.replace(f"{data['structures'][0]['residues_traced']} of {data['length']}", f"{data['structures'][0]['residues_traced'] + 1} of {data['length']}", 1)
        self.assertNotEqual(wrong, page)
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "page.html"
            p.write_text(wrong, encoding="utf-8")
            r = subprocess.run([sys.executable, str(PKG / "make_page_blocks.py"), "--page", str(p), "--check"], capture_output=True, text=True)
            self.assertEqual(r.returncode, 1, "the check must fail when a generated number is wrong")

    def test_every_generated_block_has_both_markers_in_the_standalone_page(self):
        page = (PKG / "web" / "index.html").read_text(encoding="utf-8")
        for name in ("table", "provenance", "cite"):
            self.assertEqual(page.count(f"<!--generated:cftr-{name}-->"), 1, name)
            self.assertEqual(page.count(f"<!--/generated:cftr-{name}-->"), 1, name)

    def test_the_caption_the_viewer_reproduces_is_the_one_the_generator_writes(self):
        data = gen.load_data()
        st = data["structures"][0]
        self.assertEqual(gen.caption_text(st, data["length"]),
                         "CFTR is a human protein, and variants of it are linked to cystic fibrosis. Structure 6MSM places 1,181 of its 1,480 amino acids. "
                         "Residue 508 is present in this structure. A model fitted to electron-microscope data, not a photograph.")

    def test_the_legend_block_is_the_key_the_viewer_writes(self):
        """The same seven rows as partsOf and legendRow in cftr-viewer.js (the website's hero test compares the two in a browser)."""
        import re
        data = gen.load_data()
        rows = re.findall(r'<li><span class="cftr-swatch ([^"]+)" aria-hidden="true"></span>([^<]+)</li>', gen.legend_block(data))
        self.assertEqual(len(rows), 7)
        parts = [d for d in data["domains"] if d["type"] == "Domain"]
        self.assertEqual(len(parts), 4)
        self.assertEqual(rows[0], ("sw-part-1", f"Membrane-spanning part 1, residues {parts[0]['start']} to {parts[0]['end']}"))
        self.assertEqual([c for c, _ in rows], ["sw-part-1", "sw-part-2", "sw-part-3", "sw-part-4", "sw-grey", "sw-grey is-dotted", "sw-band is-band"])
        self.assertEqual([t for _, t in rows[4:]], ["Gray: placed, outside the four marked parts", "Dotted gray: not placed in this structure",
                                                    "Flat band: helix or strand, as UniProt annotates"])
        self.assertEqual(sorted(t.split(",")[0] for _, t in rows[:4]), ["Membrane-spanning part 1", "Membrane-spanning part 2", "Nucleotide-binding part 1", "Nucleotide-binding part 2"])
        css = (PKG / "web" / "cftr-model.css").read_text(encoding="utf-8")
        for i, (r, g, b) in enumerate(gen.PART_RGB, start=1):
            self.assertIn(f".sw-part-{i} {{ background: rgb({r}, {g}, {b}); }}", css)


# Colour-vision simulation for the colour test below: Machado, Oliveira and Fernandes (2009), full severity, applied to linear sRGB; CIE76 difference in CIELAB (D65).
_MACHADO = {"deuteranopia": ((0.367322, 0.860646, -0.227968), (0.280085, 0.672501, 0.047413), (-0.011820, 0.042940, 0.968881)),
            "protanopia": ((0.152286, 1.052583, -0.204868), (0.114503, 0.786281, 0.099216), (-0.003882, -0.048116, 1.051998))}


def _lin(c):
    c /= 255
    return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4


def _enc(c):
    c = min(1.0, max(0.0, c))
    return 255 * (12.92 * c if c <= 0.0031308 else 1.055 * c ** (1 / 2.4) - 0.055)


def _simulate(rgb, kind):
    lin = [_lin(v) for v in rgb]
    return tuple(_enc(sum(row[j] * lin[j] for j in range(3))) for row in _MACHADO[kind])


def _lab(rgb):
    r, g, b = (_lin(v) for v in rgb)
    x, y, z = 0.4124 * r + 0.3576 * g + 0.1805 * b, 0.2126 * r + 0.7152 * g + 0.0722 * b, 0.0193 * r + 0.1192 * g + 0.9505 * b
    f = lambda t: t ** (1 / 3) if t > 0.008856 else 7.787 * t + 16 / 116  # noqa: E731
    fx, fy, fz = f(x / 0.95047), f(y), f(z / 1.08883)
    return 116 * fy - 16, 500 * (fx - fy), 200 * (fy - fz)


def _de76(a, b):
    return sum((p - q) ** 2 for p, q in zip(_lab(a), _lab(b))) ** 0.5


class ColourTests(unittest.TestCase):
    def test_the_viewer_and_the_generator_use_the_same_part_colours(self):
        import re
        js = (PKG / "web" / "cftr-viewer.js").read_text(encoding="utf-8")
        body = re.search(r"const PART_COLOURS = \[(.*?)\];", js).group(1)
        triples = [tuple(round(float(v) * 255) for v in t.split(",")) for t in re.findall(r"\[([^\[\]]+)\]", body)]
        self.assertEqual(triples, [tuple(c) for c in gen.PART_RGB])

    def test_the_four_parts_stay_apart_for_red_green_colour_blindness(self):
        """Every pair of part colours, and each against the gray of the rest of the chain, differs by CIE76 >= 20 under simulated deuteranopia and protanopia.
        Before 2026-10-31 nucleotide-binding part 1 (violet) and membrane-spanning part 2 (blue) were 10.3 apart under deuteranopia."""
        colours = [tuple(c) for c in gen.PART_RGB] + [gen.GREY_RGB]
        for kind in _MACHADO:
            for i in range(len(colours)):
                for j in range(i + 1, len(colours)):
                    d = _de76(_simulate(colours[i], kind), _simulate(colours[j], kind))
                    self.assertGreaterEqual(d, 20, f"{kind}: {colours[i]} and {colours[j]} are only {d:.1f} apart")
        self.assertAlmostEqual(_de76(_simulate((209, 92, 255), "deuteranopia"), _simulate((77, 140, 255), "deuteranopia")), 10.3, delta=0.1,
                               msg="the simulation reproduces the measured value for the old pair")


if __name__ == "__main__":
    unittest.main()
