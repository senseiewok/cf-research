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
        self.assertEqual([t for _, t in rows[4:]], ["Grey: placed, outside the four marked parts", "Dotted grey: not placed in this structure",
                                                    "Flat band: helix or strand, as UniProt annotates"])
        self.assertEqual(sorted(t.split(",")[0] for _, t in rows[:4]), ["Membrane-spanning part 1", "Membrane-spanning part 2", "Nucleotide-binding part 1", "Nucleotide-binding part 2"])
        css = (PKG / "web" / "cftr-model.css").read_text(encoding="utf-8")
        for i, (r, g, b) in enumerate(gen.PART_RGB, start=1):
            self.assertIn(f".sw-part-{i} {{ background: rgb({r}, {g}, {b}); }}", css)


if __name__ == "__main__":
    unittest.main()
