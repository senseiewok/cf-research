#!/usr/bin/env python3
"""Independent Second-Method Verifier for Apo 5UAK Packing Compliance Benchmark (PRED-0013).

Zero shared code, subprocesses, or imports from `cryoem_packing_compliance_benchmark.py`,
`tools/`, `subprocess`, `re`, `shlex`, or `math`:
- AST self-check verifies zero forbidden imports before executing.
- Tokenizes `5UAK.cif.gz` character-by-character via a state-machine lexer.
- Parses canonical splice dinucleotides (`c.<digits>[+-][12]`) character-by-character.
- Parses decimal Cartesian strings directly into `fractions.Fraction` and evaluates
  r_{ij}^2 = dx^2 + dy^2 + dz^2, Pi_5UAK(i) = (361/25) * sum_{|j-i|>1} 1/r_{ij}^2,
  and s_flex(i) = 1 / (1 + Pi_5UAK(i)) in 100% EXACT rational Q arithmetic (zero square roots).
- Verifies exact Fraction(0) invariance across all 1,139 placed residues under SO(3, Q) |x Q^3.
- Computes exact rational ROC-AUC via `fractions.Fraction` Wilcoxon mid-rank sum:
      U = R_pos - |S_1|(|S_1| + 1) / 2.
"""
from __future__ import annotations

import argparse
import ast
import fractions
import gzip
import hashlib
import json
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
CIF_5UAK = REPO_ROOT / "tools" / "cftr-model" / "data" / "5UAK.cif.gz"
UNIPROT_JSON = REPO_ROOT / "tools" / "cftr-model" / "data" / "p13569.json"
ANNOTATIONS_JSON = REPO_ROOT / "benchmarks" / "theratype" / "theratype_open_annotations_2026-11-21.json"

EXPECTED_5UAK_SHA256 = "8bf00f2612dd961d3b1ccd2f2193781a534ad44f32666a6de8da948afb8d84cf"
EXPECTED_P13569_SHA256 = "c8feabfabe8690ef4a06af7de5b5f8dc3c8da554246c5ab3d403edf265d4a893"

R0_SQ_Q = Fraction(361, 25)  # (19/5)^2 = (3.8 A)^2

AA_3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

FORBIDDEN_IMPORTS = {
    "cryoem_packing_compliance_benchmark",
    "tools",
    "subprocess",
    "re",
    "shlex",
    "math",
}


def verify_ast_independence(self_path: Path) -> None:
    tree = ast.parse(self_path.read_text(encoding="utf-8"), filename=str(self_path))
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                root_mod = alias.name.split(".")[0]
                if root_mod in FORBIDDEN_IMPORTS:
                    raise AssertionError(f"Forbidden import in independent_audit.py: {alias.name}")
        elif isinstance(node, ast.ImportFrom) and node.module:
            root_mod = node.module.split(".")[0]
            if root_mod in FORBIDDEN_IMPORTS:
                raise AssertionError(f"Forbidden from-import in independent_audit.py: {node.module}")


def lex_mmcif_line(line: str) -> list[str]:
    lexemes: list[str] = []
    n = len(line)
    idx = 0
    while idx < n:
        ch = line[idx]
        if ch in (" ", "\t", "\r", "\n"):
            idx += 1
            continue
        if ch in ("'", '"'):
            q = ch
            idx += 1
            buf: list[str] = []
            while idx < n:
                if line[idx] == q and (idx + 1 == n or line[idx + 1] in (" ", "\t", "\r", "\n")):
                    idx += 1
                    break
                buf.append(line[idx])
                idx += 1
            lexemes.append("".join(buf))
        else:
            buf = []
            while idx < n and line[idx] not in (" ", "\t", "\r", "\n"):
                buf.append(line[idx])
                idx += 1
            lexemes.append("".join(buf))
    return lexemes


def is_canonical_splice_charwise(hgvs: str) -> bool:
    s = hgvs.strip()
    if not s.startswith("c.") or len(s) < 5:
        return False
    pos = 2
    while pos < len(s) and s[pos].isdigit():
        pos += 1
    if pos == 2 or pos >= len(s):
        return False
    if s[pos] not in ("+", "-"):
        return False
    pos += 1
    if pos >= len(s) or s[pos] not in ("1", "2"):
        return False
    pos += 1
    return pos == len(s) or not s[pos].isdigit()


def load_rational_5uak_ca() -> dict[int, tuple[Fraction, Fraction, Fraction]]:
    if hashlib.sha256(CIF_5UAK.read_bytes()).hexdigest() != EXPECTED_5UAK_SHA256:
        raise AssertionError("5UAK.cif.gz SHA-256 mismatch")
    if hashlib.sha256(UNIPROT_JSON.read_bytes()).hexdigest() != EXPECTED_P13569_SHA256:
        raise AssertionError("p13569.json SHA-256 mismatch")

    p13569_seq = str(json.loads(UNIPROT_JSON.read_text(encoding="utf-8"))["sequence"]["value"])
    if len(p13569_seq) != 1480:
        raise AssertionError("P13569 length != 1480")

    lines = gzip.open(CIF_5UAK, "rt", encoding="utf-8").read().splitlines()
    i = 0
    while i < len(lines) and not lines[i].startswith("_atom_site."):
        i += 1
    cols: list[str] = []
    while i < len(lines) and lines[i].startswith("_atom_site."):
        cols.append(lines[i].strip()[len("_atom_site."):])
        i += 1

    coords: dict[int, tuple[Fraction, Fraction, Fraction]] = {}
    while i < len(lines) and not lines[i].startswith(("#", "_", "loop_")):
        raw = lines[i].strip()
        i += 1
        if not raw:
            continue
        toks = lex_mmcif_line(raw)
        if len(toks) != len(cols):
            continue
        rec = dict(zip(cols, toks))
        if (
            rec["group_PDB"] == "ATOM"
            and rec.get("pdbx_PDB_model_num", "1") == "1"
            and rec["auth_asym_id"] == "A"
            and rec["label_asym_id"] == "A"
            and rec["label_atom_id"] == "CA"
            and rec.get("pdbx_PDB_ins_code", "?") in ("?", ".")
            and rec.get("label_alt_id", ".") in (".", "A")
        ):
            seq = int(rec["auth_seq_id"])
            if 1 <= seq <= 1480:
                aa1 = AA_3_TO_1[rec["label_comp_id"]]
                if aa1 != p13569_seq[seq - 1]:
                    raise AssertionError(f"Sequence mismatch at {seq}")
                if seq in coords:
                    raise AssertionError(f"Duplicate CA at {seq}")
                coords[seq] = (
                    Fraction(rec["Cartn_x"]),
                    Fraction(rec["Cartn_y"]),
                    Fraction(rec["Cartn_z"]),
                )
    return coords


def exact_s_flex_single(
    seq_i: int,
    coords: dict[int, tuple[Fraction, Fraction, Fraction]],
) -> Fraction:
    xi, yi, zi = coords[seq_i]
    wcn_q = Fraction(0)
    count = 0
    for seq_j, (xj, yj, zj) in coords.items():
        if abs(seq_j - seq_i) <= 1:
            continue
        dx = xi - xj
        dy = yi - yj
        dz = zi - zj
        d2 = dx * dx + dy * dy + dz * dz
        if d2 <= 0:
            raise AssertionError(f"Zero squared distance at ({seq_i}, {seq_j})")
        wcn_q += Fraction(1, 1) / d2
        count += 1
    if count == 0:
        raise AssertionError(f"Empty non-covalent set at {seq_i}")
    return Fraction(1, 1) / (Fraction(1, 1) + R0_SQ_Q * wcn_q)


def wilcoxon_midrank_auroc_exact(
    labels: list[int],
    scores: list[Fraction | None],
) -> Fraction:
    n_pos_all = sum(1 for y in labels if y == 1)
    n_neg_all = sum(1 for y in labels if y == 0)
    if n_pos_all == 0 or n_neg_all == 0:
        raise AssertionError("Empty positive or negative class")

    scored = [(scores[idx], labels[idx]) for idx in range(len(labels)) if scores[idx] is not None]
    s_pos = sum(1 for _, y in scored if y == 1)
    s_neg = sum(1 for _, y in scored if y == 0)
    if s_pos == 0 or s_neg == 0:
        return Fraction(1, 2)

    # Sort scored items ascending by exact Fraction value
    scored_sorted = sorted(scored, key=lambda pair: pair[0])  # type: ignore[arg-type]
    m = len(scored_sorted)
    rank_sum_pos = Fraction(0)
    k = 0
    while k < m:
        val = scored_sorted[k][0]
        end = k + 1
        while end < m and scored_sorted[end][0] == val:
            end += 1
        # 1-based ranks are (k + 1) .. end; mid-rank = (k + 1 + end) / 2
        mid_rank = Fraction(k + 1 + end, 2)
        pos_in_tie = sum(1 for idx_t in range(k, end) if scored_sorted[idx_t][1] == 1)
        rank_sum_pos += mid_rank * pos_in_tie
        k = end

    u_scored = rank_sum_pos - Fraction(s_pos * (s_pos + 1), 2)
    auc_cond = u_scored / Fraction(s_pos * s_neg)
    weight = Fraction(s_pos * s_neg, n_pos_all * n_neg_all)
    return Fraction(1, 2) + weight * (auc_cond - Fraction(1, 2))


def run_independent_audit() -> dict[str, object]:
    verify_ast_independence(Path(__file__).resolve())
    coords = load_rational_5uak_ca()
    if len(coords) != 1139:
        raise AssertionError(f"Expected 1139 placed CA, got {len(coords)}")
    if 507 not in coords or 508 not in coords:
        raise AssertionError("WT 5UAK must place codons 507 and 508")

    # Verify exact 4-point collinear 3.8 A chain in Q
    c4_coords = {
        k: (Fraction(19, 5) * (k - 1), Fraction(0), Fraction(0))
        for k in (1, 2, 3, 4)
    }
    c4_exact = [exact_s_flex_single(k, c4_coords) for k in (1, 2, 3, 4)]
    if c4_exact != [Fraction(36, 49), Fraction(4, 5), Fraction(4, 5), Fraction(36, 49)]:
        raise AssertionError(f"4-point exact Q control failed: {c4_exact}")

    # Proper Pythagorean rotation matrix R_+ in SO(3, Q) (det = +1, R_+^T R_+ = I)
    r_mat = (
        (Fraction(2, 3), Fraction(-2, 3), Fraction(1, 3)),
        (Fraction(2, 3), Fraction(1, 3), Fraction(-2, 3)),
        (Fraction(1, 3), Fraction(2, 3), Fraction(2, 3)),
    )
    det_r = (
        r_mat[0][0] * (r_mat[1][1] * r_mat[2][2] - r_mat[1][2] * r_mat[2][1])
        - r_mat[0][1] * (r_mat[1][0] * r_mat[2][2] - r_mat[1][2] * r_mat[2][0])
        + r_mat[0][2] * (r_mat[1][0] * r_mat[2][1] - r_mat[1][1] * r_mat[2][0])
    )
    if det_r != Fraction(1, 1):
        raise AssertionError(f"Expected det(R_+) == 1, got {det_r}")
    for ri in range(3):
        for ci in range(3):
            dot_rc = sum(r_mat[k][ri] * r_mat[k][ci] for k in range(3))
            if dot_rc != (Fraction(1, 1) if ri == ci else Fraction(0, 1)):
                raise AssertionError(f"Non-orthogonal R_+ at ({ri}, {ci}): {dot_rc}")

    rot_coords: dict[int, tuple[Fraction, Fraction, Fraction]] = {}
    t_x, t_y, t_z = Fraction(69, 4), Fraction(-83, 2), Fraction(73, 8)
    for k, (x, y, z) in coords.items():
        rx = r_mat[0][0] * x + r_mat[0][1] * y + r_mat[0][2] * z + t_x
        ry = r_mat[1][0] * x + r_mat[1][1] * y + r_mat[1][2] * z + t_y
        rz = r_mat[2][0] * x + r_mat[2][1] * y + r_mat[2][2] * z + t_z
        rot_coords[k] = (rx, ry, rz)

    # Exact SO(3, Q) |x Q^3 zero-drift certificate across ALL 1,139 residues (all 646,952 non-covalent pairs):
    # scaling by 24,000 clears denominators of 3-decimal mmCIF strings (1,000), 1/3 rotation, and 1/8 translation.
    scale_int = 24000
    sorted_keys = sorted(coords.keys())
    scaled_orig: list[tuple[int, int, int, int]] = []
    scaled_rot: list[tuple[int, int, int, int]] = []
    for k in sorted_keys:
        ox, oy, oz = (coords[k][0] * scale_int, coords[k][1] * scale_int, coords[k][2] * scale_int)
        rx, ry, rz = (rot_coords[k][0] * scale_int, rot_coords[k][1] * scale_int, rot_coords[k][2] * scale_int)
        if any(q.denominator != 1 for q in (ox, oy, oz, rx, ry, rz)):
            raise AssertionError(f"Unexpected rational denominator at codon {k}")
        scaled_orig.append((k, ox.numerator, oy.numerator, oz.numerator))
        scaled_rot.append((k, rx.numerator, ry.numerator, rz.numerator))

    noncovalent_pairs = 0
    n_all = len(scaled_orig)
    for idx_a in range(n_all):
        ka, ax, ay, az = scaled_orig[idx_a]
        _, rax, ray, raz = scaled_rot[idx_a]
        for idx_b in range(idx_a + 1, n_all):
            kb, bx, by, bz = scaled_orig[idx_b]
            if kb - ka <= 1:
                continue
            _, rbx, rby, rbz = scaled_rot[idx_b]
            d2_o = (ax - bx) * (ax - bx) + (ay - by) * (ay - by) + (az - bz) * (az - bz)
            d2_r = (rax - rbx) * (rax - rbx) + (ray - rby) * (ray - rby) + (raz - rbz) * (raz - rbz)
            if d2_o <= 0 or d2_o != d2_r:
                raise AssertionError(f"Non-zero SO(3, Q) pair drift at ({ka}, {kb}): {d2_o} vs {d2_r}")
            noncovalent_pairs += 1

    # Evaluate exact Q s_flex on all needed codons (controls + 30 benchmark codons) and verify c1 control bounds
    annotations = json.loads(ANNOTATIONS_JSON.read_text(encoding="utf-8"))
    variants = annotations["variants"]
    exposed_ctrls = (38, 1347, 1348, 1371)
    buried_ctrls = (195, 311, 361, 1030, 1102, 1293)
    needed_codons = set(exposed_ctrls) | set(buried_ctrls) | {178, 85}
    for v in variants:
        if v["variant_family"] == "missense":
            needed_codons.add(int(v["codon_position"]))

    s_flex_q: dict[int, Fraction] = {}
    for seq in sorted(needed_codons):
        if seq in coords:
            v_orig = exact_s_flex_single(seq, coords)
            v_rot = exact_s_flex_single(seq, rot_coords)
            if v_orig != v_rot:
                raise AssertionError(f"Non-zero SO(3, Q) drift at {seq}: {v_orig} != {v_rot}")
            s_flex_q[seq] = v_orig

    for ec in exposed_ctrls:
        if s_flex_q[ec] < Fraction(50936, 1000000):
            raise AssertionError(f"Exposed control {ec} below 0.050936: {float(s_flex_q[ec])}")
    for bc in buried_ctrls:
        if s_flex_q[bc] > Fraction(47689, 1000000):
            raise AssertionError(f"Buried control {bc} above 0.047689: {float(s_flex_q[bc])}")
    if not (s_flex_q[178] > s_flex_q[85]):
        raise AssertionError("Train control G178 > G85 failed in exact Q audit")

    test_vars = [v for v in variants if v["split"] == "test"]
    train_vars = [v for v in variants if v["split"] == "train"]

    def eval_q(sub_vars: list[dict[str, object]]) -> dict[str, object]:
        labels_all = [int(v["trikafta_v23_responsive"]) for v in sub_vars]  # type: ignore[arg-type]
        mis_vars = [v for v in sub_vars if v["variant_family"] == "missense"]
        labels_mis = [int(v["trikafta_v23_responsive"]) for v in mis_vars]  # type: ignore[arg-type]

        sf_mis: list[Fraction | None] = [
            s_flex_q.get(int(v["codon_position"]), None)  # type: ignore[arg-type]
            for v in mis_vars
        ]
        sam_mis: list[Fraction | None] = [
            (Fraction(1, 1) - Fraction(str(v["alphamissense_score"])))
            if v.get("alphamissense_score") is not None
            else None
            for v in mis_vars
        ]

        gf_all: list[Fraction | None] = []
        gam_all: list[Fraction | None] = []
        seq_all: list[Fraction | None] = []
        for v in sub_vars:
            fam = str(v["variant_family"])
            hgvs = str(v["hgvs_p_or_c"])
            if fam == "class_i_ptc_frameshift":
                gf_all.append(Fraction(0))
                gam_all.append(Fraction(0))
                seq_all.append(Fraction(0))
            elif fam == "splice_regulatory":
                sp_v = Fraction(0) if is_canonical_splice_charwise(hgvs) else Fraction(1, 2)
                gf_all.append(sp_v)
                gam_all.append(sp_v)
                seq_all.append(sp_v)
            elif fam == "inframe_indel":
                gf_all.append(Fraction(1, 2))
                gam_all.append(Fraction(1, 2))
                seq_all.append(Fraction(1, 2))
            elif fam == "missense":
                sf = s_flex_q.get(int(v["codon_position"]), None)  # type: ignore[arg-type]
                gf_all.append(Fraction(1, 2) + Fraction(1, 2) * sf if sf is not None else None)
                am = v.get("alphamissense_score")
                gam_all.append(
                    Fraction(1, 2) + Fraction(1, 2) * (Fraction(1, 1) - Fraction(str(am)))
                    if am is not None
                    else None
                )
                seq_all.append(Fraction(3, 4))
            else:
                raise AssertionError(f"Unexpected family {fam}")

        auc_cond_sf = wilcoxon_midrank_auroc_exact(labels_mis, sf_mis)
        auc_cond_sam = wilcoxon_midrank_auroc_exact(labels_mis, sam_mis)
        auc_abst_gf = wilcoxon_midrank_auroc_exact(labels_all, gf_all)
        auc_abst_gam = wilcoxon_midrank_auroc_exact(labels_all, gam_all)
        auc_seq = wilcoxon_midrank_auroc_exact(labels_all, seq_all)
        return {
            "auc_cond_s_flex_fraction": f"{auc_cond_sf.numerator}/{auc_cond_sf.denominator}",
            "auc_cond_s_flex": round(float(auc_cond_sf), 6),
            "auc_cond_s_am_fraction": f"{auc_cond_sam.numerator}/{auc_cond_sam.denominator}",
            "auc_cond_s_am": round(float(auc_cond_sam), 6),
            "auc_abst_g_flex_fraction": f"{auc_abst_gf.numerator}/{auc_abst_gf.denominator}",
            "auc_abst_g_flex": round(float(auc_abst_gf), 6),
            "auc_abst_g_am_fraction": f"{auc_abst_gam.numerator}/{auc_abst_gam.denominator}",
            "auc_abst_g_am": round(float(auc_abst_gam), 6),
            "auc_seq_fraction": f"{auc_seq.numerator}/{auc_seq.denominator}",
            "auc_seq": round(float(auc_seq), 6),
        }

    test_q = eval_q(test_vars)
    train_q = eval_q(train_vars)
    return {
        "placed_ca_5uak": len(coords),
        "noncovalent_pairs_verified": noncovalent_pairs,
        "so3_q_drift": "0/1",
        "train": train_q,
        "test": test_q,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Independent exact rational audit for PRED-0013.")
    parser.add_argument("--check", action="store_true", help="Run exact Q audit")
    parser.add_argument("--json", action="store_true", help="Output JSON")
    args = parser.parse_args()

    res = run_independent_audit()
    if args.json:
        print(json.dumps(res, indent=2))
    else:
        te = res["test"]  # type: ignore[index]
        print(
            f"PASS independent_audit (PRED-0013): 5UAK_CA={res['placed_ca_5uak']} "
            f"SO3_Q_drift={res['so3_q_drift']} "
            f"test_AUC_cond={te['auc_cond_s_flex']} ({te['auc_cond_s_flex_fraction']}) "  # type: ignore[index]
            f"test_AUC_abst={te['auc_abst_g_flex']} ({te['auc_abst_g_flex_fraction']})"  # type: ignore[index]
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
