"""Write the text parts of a CFTR page from web/cftr-data.js, so the words and the picture read the same data.

Blocks, each between <!--generated:NAME--> and <!--/generated:NAME--> markers in the page:
  cftr-table       a table of the four structures: what each file calls itself, method, residues placed, residue 508, the paper, the stretches not placed, how the file differs from UniProt's sequence
  cftr-provenance  where each number on the page comes from and how to check it (values computed from the data)
  cftr-cite        the entries, papers and UniProt version to cite, as plain text
  cftr-poster      a still vector picture of the first structure in the viewer's starting view (shown without JavaScript or while the 3D view loads)
  cftr-caption     the caption the viewer's script reproduces word for word
  cftr-legend      the colour key, word for word what the viewer's script writes into the same list
Nothing is typed in: every number is computed from the data module.

Usage: python tools/cftr-model/make_page_blocks.py [--page PATH] [--data PATH] [--blocks table,provenance,cite] [--check]
       With no options it updates web/index.html. --check exits 1 if a page differs from what the data would produce.
The website imports this module (its pages use the same functions), so there is one generator.
Standard library only."""
import argparse
import html
import json
import math
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAGE = HERE / "web" / "index.html"
DATA = HERE / "web" / "cftr-data.js"
AMINO = {"ALA": "alanine", "ARG": "arginine", "ASN": "asparagine", "ASP": "aspartate", "CYS": "cysteine", "GLN": "glutamine", "GLU": "glutamate", "GLY": "glycine",
         "HIS": "histidine", "ILE": "isoleucine", "LEU": "leucine", "LYS": "lysine", "MET": "methionine", "PHE": "phenylalanine", "PRO": "proline", "SER": "serine",
         "THR": "threonine", "TRP": "tryptophan", "TYR": "tyrosine", "VAL": "valine"}       # the data names residues by their three-letter codes


def load_data(path=None):
    """The data module is one line of JSON after `export default`."""
    text = Path(path or DATA).read_text(encoding="utf-8")
    body = text[text.index("export default ") + len("export default "):].strip()
    return json.loads(body.rstrip(";"))


def unplaced(st, length):
    """The stretches of 1..length that the model does not place, as (first, last) pairs."""
    covered = set()
    for seg in st["segments"]:
        covered.update(range(seg["start"], seg["start"] + len(seg["xyz"]) // 3))
    out, start = [], None
    for n in range(1, length + 2):
        missing = n <= length and n not in covered
        if missing and start is None:
            start = n
        elif not missing and start is not None:
            out.append((start, n - 1))
            start = None
    return out


def stretch_text(pairs):
    return ", ".join(str(a) if a == b else f"{a} to {b}" for a, b in pairs)


def differences_text(st):
    """How the file's own sequence differs from UniProt's, in the file's words."""
    parts = []
    for d in st["differences"]:
        pos = d["position"]
        was = AMINO.get(d["uniprot_residue"], d["uniprot_residue"] or "?")
        now = AMINO.get(d["file_residue"], d["file_residue"]) if d["file_residue"] else None
        if now:
            parts.append(f"{pos}: {now} in the file, {was} in UniProt (the file says: {d['details']})")
        else:
            parts.append(f"{pos}: {was} is not in the file (the file says: {d['details']})")
    if st["expression_tag_residues"]:
        parts.append(f"{st['expression_tag_residues']} expression-tag residues beyond the protein, not drawn")
    return "; ".join(parts) if parts else "none recorded"


def table(data):
    esc = html.escape
    rows = []
    for st in data["structures"]:
        gaps = unplaced(st, data["length"])
        paper = f'{esc(st["journal"])} {esc(st["paper_year"])}, <a href="https://pubmed.ncbi.nlm.nih.gov/{esc(st["pubmed"])}/" rel="noopener noreferrer">PubMed {esc(st["pubmed"])}</a>'
        r508 = f'present ({esc(AMINO.get(st["res508"]["name"], st["res508"]["name"].lower()))})' if st.get("res508") else "absent: the chain runs from 507 to 509"
        rows.append(
            "      <tr>"
            f'<th scope="row"><a href="https://www.rcsb.org/structure/{esc(st["id"])}" rel="noopener noreferrer">{esc(st["id"])}</a></th>'
            f'<td>{esc(st["title"])}</td>'
            f'<td>Electron microscopy, {esc(st["resolution_A"])} &Aring;; deposited {esc(st["deposited"])}</td>'
            f'<td>{st["residues_traced"]} of {data["length"]}</td>'
            f"<td>{r508}</td>"
            f"<td>{paper}</td>"
            f"<td>{esc(stretch_text(gaps))}</td>"
            f"<td>{esc(differences_text(st))}</td>"
            "</tr>")
    head = ("      <tr><th scope=\"col\">Structure</th><th scope=\"col\">What the file calls it</th><th scope=\"col\">Method</th><th scope=\"col\">Residues placed</th>"
            "<th scope=\"col\">Residue 508</th><th scope=\"col\">Paper</th><th scope=\"col\">Residues not placed</th><th scope=\"col\">Differs from UniProt's sequence</th></tr>")
    return ('<div class="cftr-tablewrap" tabindex="0" role="region" aria-label="Table of the four structures, scrolls sideways on narrow screens">\n'
            '    <table class="cftr-table">\n'
            "    <caption>Four human CFTR structures from the Protein Data Bank (CC0). Residue numbers follow UniProt P13569.</caption>\n"
            f"    <thead>\n{head}\n    </thead>\n    <tbody>\n" + "\n".join(rows) + "\n    </tbody>\n    </table>\n  </div>")


HOME_THETA, HOME_PHI = 0.6, 1.35          # the viewer's starting view (cftr-viewer.js, HOME)
PART_RGB = [(0, 242, 255), (209, 92, 255), (77, 140, 255), (255, 173, 26)]
GREY_RGB, ROSE_RGB = (140, 153, 168), (255, 120, 150)


def rgb(c):
    return f"rgb({c[0]},{c[1]},{c[2]})"


def caption_text(st, length):
    r508 = "Residue 508 is present in this structure." if st.get("res508") else "Residue 508 is absent from this structure."
    return (f"CFTR is a human protein, and variants of it are linked to cystic fibrosis. Structure {st['id']} places {st['residues_traced']:,} of its {length:,} amino acids. "
            f"{r508} A model fitted to electron-microscope data, not a photograph.")


def poster_svg(data):
    """A still picture of the first structure in the viewer's starting view: an orthographic projection of the C-alpha trace, coloured by part."""
    st = data["structures"][0]
    parts = [d for d in data["domains"] if d["type"] == "Domain"]
    st_phi, st_theta = HOME_PHI, HOME_THETA
    e = (math.sin(st_phi) * math.cos(st_theta), math.sin(st_phi) * math.sin(st_theta), math.cos(st_phi))
    up = (0.0, 0.0, 1.0)
    cross = lambda a, b: (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])
    unit = lambda v: tuple(c / math.sqrt(sum(x * x for x in v)) for c in v)
    xa = unit(cross(up, e))
    ya = cross(e, xa)
    dot = lambda a, b: a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
    runs, allpts = [], []
    for seg in st["segments"]:
        pts = []
        for i in range(len(seg["xyz"]) // 3):
            p = seg["xyz"][3 * i:3 * i + 3]
            pts.append((seg["start"] + i, dot(p, xa), -dot(p, ya)))
        runs.append(pts)
        allpts += pts
    xs, ys = [q[1] for q in allpts], [q[2] for q in allpts]
    x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
    pad = 14.0
    scale = 600.0 / (y1 - y0)
    W, H = (x1 - x0) * scale + 2 * pad, (y1 - y0) * scale + 2 * pad
    sx = lambda x: round((x - x0) * scale + pad, 1)
    sy = lambda y: round((y - y0) * scale + pad, 1)
    colour = lambda n: next((PART_RGB[i % 4] for i, d in enumerate(parts) if d["start"] <= n <= d["end"]), GREY_RGB)
    paths = []
    for pts in runs:
        i = 0
        while i < len(pts):
            c = colour(pts[i][0])
            j = i
            while j + 1 < len(pts) and colour(pts[j + 1][0]) == c:
                j += 1
            seg_pts = pts[i:j + 2] if j + 1 < len(pts) else pts[i:j + 1]       # overlap one point so colour runs join
            if len(seg_pts) > 1:
                d = "M" + " L".join(f"{sx(p[1])} {sy(p[2])}" for p in seg_pts)
                paths.append(f'<path d="{d}" stroke="{rgb(c)}" />')
            i = j + 1
    bridges = []
    for a, b in zip(runs, runs[1:]):
        bridges.append(f'<path d="M{sx(a[-1][1])} {sy(a[-1][2])} L{sx(b[0][1])} {sy(b[0][2])}" />')
    ring = ""
    if st.get("res508"):
        p = st["res508"]["xyz"]
        cx, cy = sx(dot(p, xa)), sy(-dot(p, ya))
        ring = f'<circle cx="{cx}" cy="{cy}" r="11" fill="none" stroke="{rgb(ROSE_RGB)}" stroke-width="3" /><circle cx="{cx}" cy="{cy}" r="4.5" fill="{rgb(ROSE_RGB)}" />'
    label = (f"A still picture of CFTR structure {st['id']}, drawn as a line through the amino acids the model places and coloured by part. "
             "Residue 508 is ringed in rose. The interactive 3D view replaces it when the page can run it.")
    return (f'<svg class="cftr-poster-svg" viewBox="0 0 {round(W)} {round(H)}" role="img" aria-label="{html.escape(label)}" xmlns="http://www.w3.org/2000/svg">'
            f'<g fill="none" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round">{"".join(paths)}</g>'
            f'<g fill="none" stroke="{rgb(GREY_RGB)}" stroke-width="1.4" stroke-dasharray="1 4" stroke-linecap="round" opacity="0.7">{"".join(bridges)}</g>{ring}</svg>')


def legend_block(data):
    """The colour key as list items, word for word what cftr-viewer.js writes into #cftr-legend (partsOf and legendRow there), so a page can show it before and
    without the script and the script's rewrite changes nothing. Swatch colours are classes in cftr-model.css (sw-part-1 to 4 follow PART_RGB)."""
    kinds, rows = {"membrane": 0, "binding": 0}, []
    for i, d in enumerate(x for x in data["domains"] if x["type"] == "Domain"):
        kind = "membrane" if re.search("transmembrane", d["name"], re.I) else "binding"
        kinds[kind] += 1
        label = ("Membrane-spanning part " if kind == "membrane" else "Nucleotide-binding part ") + str(kinds[kind])
        rows.append((f"sw-part-{i % len(PART_RGB) + 1}", f"{label}, residues {d['start']} to {d['end']}"))
    rows += [("sw-grey", "Grey: placed, outside the four marked parts"), ("sw-grey is-dotted", "Dotted grey: not placed in this structure"),
             ("sw-band is-band", "Flat band: helix or strand, as UniProt annotates")]
    return "".join(f'<li><span class="cftr-swatch {c}" aria-hidden="true"></span>{html.escape(t)}</li>' for c, t in rows)


def provenance_block(data):
    """Where each number on the research page comes from and how to check it. The values are computed here from the data."""
    esc = html.escape
    sts = data["structures"]
    u = data["uniprot"]
    per = lambda f: "; ".join(f"{s['id']} {f(s)}" for s in sts)
    single = [v for v in u["variants"] if v["start"] == v["end"]]
    positions = sorted({v["start"] for v in single})
    spanning = len(u["variants"]) - len(single)
    parts = [d for d in data["domains"] if d["type"] == "Domain"]
    rregion = next(d for d in data["domains"] if d["name"] == "Disordered R region")
    checked = sum(s["residues_traced"] for s in sts)
    declared = sum(1 for s in sts for d in s["differences"] if d["file_residue"])
    rows = [
        ("Residues placed", per(lambda s: f"{s['residues_traced']}"), "The four .cif.gz files from the RCSB Protein Data Bank, chain A, model 1.",
         "Count the ATOM lines whose atom name is CA in chain A, model 1."),
        ("Method, resolution, date", per(lambda s: f"{s['resolution_A']} \u00c5, deposited {s['deposited']}"), "Fields _exptl.method, _em_3d_reconstruction.resolution and _pdbx_database_status.recvd_initial_deposition_date in each file.",
         "Open the file header and read those fields."),
        ("Paper", per(lambda s: f"{s['journal']} {s['paper_year']}, PubMed {s['pubmed']}"), "Fields _citation.journal_abbrev, _citation.year and _citation.pdbx_database_id_PubMed.", "Open the PubMed record."),
        ("Stretches not placed", "listed in the table above", "Computed from the residue numbers placed.", "List the residue numbers that have a CA atom and take the gaps between 1 and 1480."),
        ("How each file differs from UniProt\u2019s sequence", per(lambda s: ", ".join(str(d["position"]) for d in s["differences"]) or "none"), "The file\u2019s own _struct_ref_seq_dif list; expression-tag residues beyond the protein are counted, not drawn.", "Read that list in the file."),
        ("Length and part boundaries", f"{data['length']:,} residues; parts " + ", ".join(f"{d['start']} to {d['end']}" for d in parts) + f"; a region UniProt marks as disordered (the regulatory, or R, region) {rregion['start']} to {rregion['end']}",
         f"UniProt P13569, entry version {u['entry_version']}, annotation updated {u['last_annotation_update']}: the sequence, and the features of type Domain and Region.", "Open the entry on uniprot.org."),
        ("Residue numbering", f"{checked:,} placed residues checked; {declared} differ from UniProt, all declared by their files", "Every placed residue is compared with the amino acid UniProt has at that number; the build stops if one disagrees without its file saying so.",
         "Run tools/cftr-model/build_traces.py (in the cf-research repository) with the downloaded files; tools/cftr-model/tests/test_data.py also checks it."),
        ("Distance after fitting", per(lambda s: f"{s['superposition_rmsd_A']:.2f} \u00c5 over {s['superposition_residues']} residues" if s["id"] != sts[0]["id"] else "the reference"),
         "Computed for this page. Take the alpha carbons that 6MSM and the structure both place, with the same amino acid, inside UniProt residues 81 to 365 and 859 to 1155; fit one onto the other by a rigid rotation and shift (Horn\u2019s method); the number is the root-mean-square distance left.",
         "Run tools/cftr-model/build_traces.py (in the cf-research repository); tools/cftr-model/tests/test_data.py shows the fit recovers a known rotation."),
        ("Positions in the view", "angstrom, rounded to 0.1", "The alpha-carbon coordinates after that fit, centred on the placed residues of 6MSM.", "Save the CSV from the tools above."),
        ("Topology and helices", f"{len(u['transmembrane'])} transmembrane segments; {len(u['topology'])} topological domains", "UniProt P13569 features of type Transmembrane and Topological domain.", "Open the entry on uniprot.org."),
        ("Variant positions", f"{len(u['variants'])} entries: {len(single)} at a single residue ({len(positions)} different positions) and {spanning} spanning several residues", "UniProt P13569 features of type Natural variant (position, residues and wording as UniProt gives them). Only single-residue entries are drawn as beads; the others are listed in the residue look-up.", "Open the entry on uniprot.org."),
    ]
    head = "<tr><th scope=\"col\">What</th><th scope=\"col\">Value on this page</th><th scope=\"col\">Where it comes from</th><th scope=\"col\">How to check it</th></tr>"
    body = "\n".join("      <tr><th scope=\"row\">" + esc(a) + "</th><td>" + esc(b) + "</td><td>" + esc(c) + "</td><td>" + esc(d) + "</td></tr>" for a, b, c, d in rows)
    return ("\n    <div class=\"cftr-tablewrap\" tabindex=\"0\" role=\"region\" aria-label=\"Table of where each number comes from, scrolls sideways on narrow screens\">\n"
            "      <table class=\"cftr-table cftr-prov\">\n      <caption>Each number on this page, its source and a way to check it.</caption>\n      <thead>\n      " + head + "\n      </thead>\n      <tbody>\n" + body + "\n      </tbody>\n      </table>\n    </div>\n    ")


def cite_block(data):
    u = data["uniprot"]
    lines = ["Structures (RCSB Protein Data Bank, archive files under CC0 1.0; please cite the entry and its paper):"]
    for s in data["structures"]:
        lines.append(f"{s['id']}: {s['paper']} {s['journal']} {s['paper_year']}. PubMed {s['pubmed']}. Entry deposited {s['deposited']}.")
    lines.append(f"Annotations: UniProt P13569, entry version {u['entry_version']}, annotation last updated {u['last_annotation_update']} (CC BY 4.0; how to cite UniProt: https://www.uniprot.org/help/publications).")
    lines.append("Fit and numbers on this page: computed from those files with tools/cftr-model/build_traces.py in the cf-research repository (alpha-carbon atoms of chain A, model 1).")
    return "\n    <pre class=\"cftr-cite\" tabindex=\"0\">" + html.escape("\n".join(lines)) + "</pre>\n    "


def splice_named(text, name, block, where):
    pat = re.compile(rf"<!--generated:{name}-->.*?<!--/generated:{name}-->", re.S)
    if not pat.search(text):
        raise SystemExit(f"markers <!--generated:{name}--> not found in {where}")
    return pat.sub(lambda m: f"<!--generated:{name}-->{block}<!--/generated:{name}-->", text)


BLOCKS = {
    "table": lambda data: "\n  " + table(data) + "\n  ",
    "provenance": provenance_block,
    "cite": cite_block,
    "poster": poster_svg,
    "caption": lambda data: html.escape(caption_text(data["structures"][0], data["length"])),
    "legend": legend_block,
}


def run(jobs, data_path=None, check=False, label=lambda p: str(p)):
    """jobs: a list of (page path, [block names]). Writes the pages, or with check=True reports the ones that differ. Returns 0 or 1."""
    data = load_data(data_path)
    status = 0
    for path, names in jobs:
        path = Path(path)
        raw = path.read_bytes().decode("utf-8")
        crlf = "\r\n" in raw
        new = raw.replace("\r\n", "\n")
        for name in names:
            new = splice_named(new, "cftr-" + name, BLOCKS[name](data), path)
        if crlf:
            new = new.replace("\n", "\r\n")
        if check:
            if new != raw:
                print(f"{label(path)} differs from the data; run the generator")
                status = 1
            else:
                print(f"{label(path)} is up to date")
        elif new != raw:
            path.write_bytes(new.encode("utf-8"))
            print(f"wrote {label(path)}")
        else:
            print(f"{label(path)} already up to date")
    return status


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--page", default=str(PAGE))
    ap.add_argument("--data", default=str(DATA))
    ap.add_argument("--blocks", default="table,provenance,cite")
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args(argv)
    return run([(a.page, a.blocks.split(","))], a.data, a.check)


if __name__ == "__main__":
    sys.exit(main())
