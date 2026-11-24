#!/usr/bin/env python3
"""Independent Second-Method Verifier for deltaF508 Cryo-EM Rescue Benchmark (PRED-0014).

AST-isolated from `cryoem_f508del_rescue_benchmark.py`, `tools/`, `subprocess`, `re`, and `shlex`:
- Tokenizes `8EJ1.cif.gz` and `8EIQ.cif.gz` via a character-level mmCIF state-machine lexer.
- Parses canonical splice dinucleotides (`c.<digits>[+-][12]`) character-by-character without regex.
- Evaluates BOTH inclusion (`min(d_A, d_B) <= 15`) and Mariani et al. (2013) distance-difference
  thresholds (`|d_A - d_B| < t` for `t in {1/2, 1, 2, 4}`) in 100% EXACT rational `fractions.Fraction`
  polynomial arithmetic (`(A + B - t^2)^2 < 4AB` when `A + B >= t^2`) with ZERO square roots,
  cross-checked against 50-digit `decimal.Decimal` minimum cutoff margin.
- Verifies exact `Fraction(0)` drift under independent rational Pythagorean rigid motions
  `(R_1, t_1), (R_2, t_2) in SO(3, Q) |>< Q^3` across all non-covalent pairs of all 1,121 shared residues.
- Computes exact rational ROC-AUC via `fractions.Fraction` Wilcoxon mid-rank sum.
"""
from __future__ import annotations

import argparse
import ast
import decimal
import gzip
import hashlib
import json
from decimal import Decimal
from fractions import Fraction
from pathlib import Path
from typing import Any, Mapping, Sequence

decimal.getcontext().prec = 50

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
CIF_8EJ1 = REPO_ROOT / "tools" / "cftr-model" / "data" / "8EJ1.cif.gz"
CIF_8EIQ = REPO_ROOT / "tools" / "cftr-model" / "data" / "8EIQ.cif.gz"
UNIPROT_JSON = REPO_ROOT / "tools" / "cftr-model" / "data" / "p13569.json"
HOLDOUT_JSON = HERE / "holdout_open_annotations_2026-11-23.json"

EXPECTED_HASHES = {
    "8EJ1": "9582203b74c0add00a7917bed765cbc8a943fcc3a88ffbac8ae8dde7de2bc159",
    "8EIQ": "d13cbfc8ba6575ab6228f351578c6f9cc5b9dafe13e5913d9b27b3cacebc18d6",
    "P13569": "c8feabfabe8690ef4a06af7de5b5f8dc3c8da554246c5ab3d403edf265d4a893",
    "HOLDOUT_VARIANTS": "8d410145e5219bfa405dbf7a24f68cedb7a5f88c5035b0da6936b5d402024bb3",
    "PRE_REGISTERED_HOLDOUT_VARIANTS": "9b815b5bc581363fc8cc7c7c4b5c31fc346d60ee2047256bc264d544062eb703",
}

AA_3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

R0_SQ_Q = Fraction(225, 1)  # 15^2
THRESH_SQ_Q = (Fraction(1, 4), Fraction(1, 1), Fraction(4, 1), Fraction(16, 1))
R0_DEC = Decimal("15.0")
THRESH_DEC = (Decimal("0.5"), Decimal("1.0"), Decimal("2.0"), Decimal("4.0"))


def assert_ast_isolation() -> None:
    src = Path(__file__).read_text(encoding="utf-8")
    tree = ast.parse(src)
    forbidden = {"cryoem_f508del_rescue_benchmark", "tools", "subprocess", "re", "shlex", "math"}
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name.split(".")[0] in forbidden:
                    raise AssertionError(f"Forbidden import in independent_audit.py: {alias.name}")
        elif isinstance(node, ast.ImportFrom):
            if (node.module or "").split(".")[0] in forbidden:
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
    if not hgvs.startswith("c.") or len(hgvs) < 7:
        return False
    pos = 2
    while pos < len(hgvs) and hgvs[pos].isdigit():
        pos += 1
    if pos == 2 or pos >= len(hgvs) or hgvs[pos] not in ("+", "-"):
        return False
    pos += 1
    if pos >= len(hgvs) or hgvs[pos] not in ("1", "2"):
        return False
    pos += 1
    if pos + 2 >= len(hgvs) or hgvs[pos].isdigit():
        return False
    return (
        hgvs[pos] in ("A", "C", "G", "T")
        and hgvs[pos + 1] == ">"
        and hgvs[pos + 2] in ("A", "C", "G", "T")
        and pos + 3 == len(hgvs)
    )


def load_rational_ca(
    cif_path: Path,
    tag: str,
    uniprot_seq: str,
) -> dict[int, tuple[Fraction, Fraction, Fraction]]:
    digest = hashlib.sha256(cif_path.read_bytes()).hexdigest()
    if digest != EXPECTED_HASHES[tag]:
        raise AssertionError(f"{tag} SHA-256 mismatch: {digest}")

    with gzip.open(cif_path, "rt", encoding="utf-8") as fh:
        lines = fh.read().splitlines()
    mutations: dict[int, str] = {}
    deletions: set[int] = set()

    i = 0
    while i < len(lines) and not lines[i].startswith("_struct_ref_seq_dif."):
        i += 1
    if i < len(lines):
        prev = lines[i - 1].strip() if i > 0 else ""
        diff_records: list[dict[str, str]] = []
        if prev == "loop_":
            cols: list[str] = []
            while i < len(lines) and lines[i].startswith("_struct_ref_seq_dif."):
                cols.append(lines[i].strip()[len("_struct_ref_seq_dif."):])
                i += 1
            while i < len(lines) and not lines[i].startswith(("#", "_", "loop_")):
                if lines[i].strip():
                    toks = lex_mmcif_line(lines[i])
                    if len(toks) == len(cols):
                        diff_records.append(dict(zip(cols, toks)))
                i += 1
        for r in diff_records:
            if r.get("pdbx_pdb_strand_id") != "A" or r.get("pdbx_seq_db_accession_code") != "P13569":
                continue
            s_num = r.get("pdbx_seq_db_seq_num", "")
            if s_num.isdigit():
                pos = int(s_num)
                if 1 <= pos <= 1480:
                    mon = r.get("mon_id", "?")
                    if mon in ("?", "."):
                        deletions.add(pos)
                    else:
                        mutations[pos] = mon

    if 508 not in deletions or mutations.get(1371) != "GLN":
        raise AssertionError(f"{tag}: construct _struct_ref_seq_dif missing del508 or E1371Q")

    atom_cols: list[str] = []
    ca_map: dict[int, tuple[str, tuple[Fraction, Fraction, Fraction]]] = {}
    in_atom = False
    for ln in lines:
        if ln.startswith("_atom_site."):
            atom_cols.append(ln.strip()[len("_atom_site."):])
            in_atom = True
        elif in_atom and ln.startswith(("ATOM ", "HETATM ")):
            toks = lex_mmcif_line(ln)
            col = {name: toks[idx] for idx, name in enumerate(atom_cols)}
            if col["group_PDB"] != "ATOM" or col["pdbx_PDB_model_num"] != "1":
                continue
            if col["auth_asym_id"] != "A" or col["label_asym_id"] != "A":
                continue
            if col["pdbx_PDB_ins_code"] not in ("?", ".") or col["label_alt_id"] not in (".", "A"):
                continue
            if col["label_atom_id"] != "CA":
                continue
            raw_seq = col["auth_seq_id"]
            if not raw_seq.isdigit():
                continue
            seq_no = int(raw_seq)
            if not (1 <= seq_no <= 1480):
                continue
            if seq_no in ca_map:
                raise AssertionError(f"{tag}: duplicate CA at {seq_no}")
            ca_map[seq_no] = (
                col["label_comp_id"],
                (Fraction(col["Cartn_x"]), Fraction(col["Cartn_y"]), Fraction(col["Cartn_z"])),
            )
        elif in_atom and ca_map and not ln.startswith(("ATOM ", "HETATM ")):
            break

    if 508 in ca_map:
        raise AssertionError(f"{tag}: deleted 508 present")
    for pos, (comp, _) in sorted(ca_map.items()):
        if pos in mutations:
            if comp != mutations[pos]:
                raise AssertionError(f"{tag}: mutation mismatch at {pos}")
        elif AA_3_TO_1.get(comp) != uniprot_seq[pos - 1]:
            raise AssertionError(f"{tag}: sequence mismatch at {pos}")

    return {pos: xyz for pos, (_, xyz) in ca_map.items()}


def sq_dist_q(p: tuple[Fraction, Fraction, Fraction], q: tuple[Fraction, Fraction, Fraction]) -> Fraction:
    dx = p[0] - q[0]
    dy = p[1] - q[1]
    dz = p[2] - q[2]
    return dx * dx + dy * dy + dz * dz


def verify_so3_q_matrix(rot: Sequence[Sequence[Fraction]]) -> None:
    for r in range(3):
        for c in range(3):
            dot = sum(rot[k][r] * rot[k][c] for k in range(3))
            if dot != (Fraction(1) if r == c else Fraction(0)):
                raise AssertionError(f"Non-orthogonal rational matrix at ({r},{c})")
    det = (
        rot[0][0] * (rot[1][1] * rot[2][2] - rot[1][2] * rot[2][1])
        - rot[0][1] * (rot[1][0] * rot[2][2] - rot[1][2] * rot[2][0])
        + rot[0][2] * (rot[1][0] * rot[2][1] - rot[1][1] * rot[2][0])
    )
    if det != Fraction(1):
        raise AssertionError(f"Rational rotation determinant {det} != 1")


def apply_so3_q(
    coords: Mapping[int, tuple[Fraction, Fraction, Fraction]],
    rot: Sequence[Sequence[Fraction]],
    trans: tuple[Fraction, Fraction, Fraction],
) -> dict[int, tuple[Fraction, Fraction, Fraction]]:
    verify_so3_q_matrix(rot)
    out: dict[int, tuple[Fraction, Fraction, Fraction]] = {}
    for k, (x, y, z) in coords.items():
        out[k] = (
            rot[0][0] * x + rot[0][1] * y + rot[0][2] * z + trans[0],
            rot[1][0] * x + rot[1][1] * y + rot[1][2] * z + trans[1],
            rot[2][0] * x + rot[2][1] * y + rot[2][2] * z + trans[2],
        )
    return out


def abs_sqrt_diff_lt_t_exact(a_sq: Fraction, b_sq: Fraction, t_sq: Fraction) -> bool:
    """Exact rational test for |sqrt(A) - sqrt(B)| < t without square roots (A, B >= 0, t > 0)."""
    lhs = a_sq + b_sq - t_sq
    if lhs < 0:
        return True
    return lhs * lhs < 4 * a_sq * b_sq


def compute_exact_fraction_s_stab(
    ca_8ej1: Mapping[int, tuple[Fraction, Fraction, Fraction]],
    ca_8eiq: Mapping[int, tuple[Fraction, Fraction, Fraction]],
) -> tuple[dict[int, Fraction | None], int, Decimal]:
    common = sorted(set(ca_8ej1.keys()) & set(ca_8eiq.keys()))
    n = len(common)
    d1_sq: list[list[Fraction]] = [[Fraction(0)] * n for _ in range(n)]
    d2_sq: list[list[Fraction]] = [[Fraction(0)] * n for _ in range(n)]
    noncovalent_pairs = 0
    min_dec_margin = Decimal("999")

    for a in range(n):
        ra = common[a]
        p1 = ca_8ej1[ra]
        p2 = ca_8eiq[ra]
        for b in range(a + 1, n):
            rb = common[b]
            if abs(ra - rb) <= 1 or (ra == 507 and rb == 509):
                continue
            noncovalent_pairs += 1
            sq1 = sq_dist_q(p1, ca_8ej1[rb])
            sq2 = sq_dist_q(p2, ca_8eiq[rb])
            d1_sq[a][b] = d1_sq[b][a] = sq1
            d2_sq[a][b] = d2_sq[b][a] = sq2

            # Cross-check 50-digit Decimal cutoff margin
            dec1 = (Decimal(sq1.numerator) / Decimal(sq1.denominator)).sqrt()
            dec2 = (Decimal(sq2.numerator) / Decimal(sq2.denominator)).sqrt()
            m_dec = dec1 if dec1 <= dec2 else dec2
            mr = abs(m_dec - R0_DEC)
            if mr < min_dec_margin:
                min_dec_margin = mr
            if m_dec <= R0_DEC:
                dd_dec = abs(dec1 - dec2)
                for td in THRESH_DEC:
                    mt = abs(dd_dec - td)
                    if mt < min_dec_margin:
                        min_dec_margin = mt

    s_stab: dict[int, Fraction | None] = {}
    for a in range(n):
        ra = common[a]
        preserved = 0
        neighbors = 0
        for b in range(n):
            rb = common[b]
            if abs(ra - rb) <= 1 or (ra == 507 and rb == 509) or (ra == 509 and rb == 507):
                continue
            sq1 = d1_sq[a][b]
            sq2 = d2_sq[a][b]
            if (sq1 if sq1 <= sq2 else sq2) <= R0_SQ_Q:
                neighbors += 1
                for t_sq in THRESH_SQ_Q:
                    if abs_sqrt_diff_lt_t_exact(sq1, sq2, t_sq):
                        preserved += 1
        s_stab[ra] = Fraction(preserved, 4 * neighbors) if neighbors > 0 else None

    return s_stab, noncovalent_pairs, min_dec_margin


def verify_all_pairs_so3_q_zero_drift(
    ca_8ej1: Mapping[int, tuple[Fraction, Fraction, Fraction]],
    ca_8eiq: Mapping[int, tuple[Fraction, Fraction, Fraction]],
    common: Sequence[int],
) -> tuple[Fraction, int]:
    rot1 = (
        (Fraction(3, 5), Fraction(-4, 5), Fraction(0, 1)),
        (Fraction(4, 5), Fraction(3, 5), Fraction(0, 1)),
        (Fraction(0, 1), Fraction(0, 1), Fraction(1, 1)),
    )
    rot2 = (
        (Fraction(2, 3), Fraction(-2, 3), Fraction(1, 3)),
        (Fraction(2, 3), Fraction(1, 3), Fraction(-2, 3)),
        (Fraction(1, 3), Fraction(2, 3), Fraction(2, 3)),
    )
    m1 = apply_so3_q(ca_8ej1, rot1, (Fraction(35, 2), Fraction(-19, 2), Fraction(337, 8)))
    m2 = apply_so3_q(ca_8eiq, rot2, (Fraction(-31, 1), Fraction(59, 4), Fraction(-33, 4)))

    max_sq_drift = Fraction(0)
    checked = 0
    n = len(common)
    for a in range(n):
        ra = common[a]
        for b in range(a + 1, n):
            rb = common[b]
            if abs(ra - rb) <= 1 or (ra == 507 and rb == 509):
                continue
            checked += 1
            d1_orig = sq_dist_q(ca_8ej1[ra], ca_8ej1[rb])
            d1_mov = sq_dist_q(m1[ra], m1[rb])
            d2_orig = sq_dist_q(ca_8eiq[ra], ca_8eiq[rb])
            d2_mov = sq_dist_q(m2[ra], m2[rb])
            diff1 = abs(d1_orig - d1_mov)
            diff2 = abs(d2_orig - d2_mov)
            if diff1 > max_sq_drift:
                max_sq_drift = diff1
            if diff2 > max_sq_drift:
                max_sq_drift = diff2
    return max_sq_drift, checked


def wilcoxon_auroc_fraction(
    pairs: Sequence[tuple[Fraction | None, int]],
    *,
    abstention_aware: bool,
) -> Fraction:
    pts = list(pairs) if abstention_aware else [(s, y) for s, y in pairs if s is not None]
    pos = [s for s, y in pts if y == 1]
    neg = [s for s, y in pts if y == 0]
    if all(s is not None for s in pos) and all(s is not None for s in neg):
        # Exact Wilcoxon mid-rank sum formula U = R_pos - n1*(n1+1)/2
        combined = [(s, 1, idx) for idx, s in enumerate(pos) if s is not None] + [
            (s, 0, idx) for idx, s in enumerate(neg) if s is not None
        ]
        combined.sort(key=lambda t: t[0])
        m = len(combined)
        r_pos = Fraction(0)
        idx = 0
        while idx < m:
            j = idx
            while j < m and combined[j][0] == combined[idx][0]:
                j += 1
            mid_rank = Fraction((idx + 1) + j, 2)
            for k in range(idx, j):
                if combined[k][1] == 1:
                    r_pos += mid_rank
            idx = j
        n1 = len(pos)
        u_stat = r_pos - Fraction(n1 * (n1 + 1), 2)
        return u_stat / Fraction(len(pos) * len(neg), 1)

    total_u = Fraction(0)
    for sp in pos:
        for sn in neg:
            if sp is None or sn is None or sp == sn:
                total_u += Fraction(1, 2)
            elif sp > sn:
                total_u += Fraction(1, 1)
    return total_u / Fraction(len(pos) * len(neg), 1)


def dec_dist(a: tuple[Fraction, Fraction, Fraction], b: tuple[Fraction, Fraction, Fraction]) -> Decimal:
    sq = sq_dist_q(a, b)
    return (Decimal(sq.numerator) / Decimal(sq.denominator)).sqrt()


def run_independent_audit(*, verify_parity: bool = True) -> dict[str, Any]:
    assert_ast_isolation()
    if hashlib.sha256(UNIPROT_JSON.read_bytes()).hexdigest() != EXPECTED_HASHES["P13569"]:
        raise AssertionError("P13569 SHA-256 mismatch")
    uniprot_seq = json.loads(UNIPROT_JSON.read_text(encoding="utf-8"))["sequence"]["value"]

    ca_8ej1 = load_rational_ca(CIF_8EJ1, "8EJ1", uniprot_seq)
    ca_8eiq = load_rational_ca(CIF_8EIQ, "8EIQ", uniprot_seq)
    common = sorted(set(ca_8ej1.keys()) & set(ca_8eiq.keys()))

    d_junc_8ej1 = dec_dist(ca_8ej1[507], ca_8ej1[509])
    d_junc_8eiq = dec_dist(ca_8eiq[507], ca_8eiq[509])
    if not (Decimal("3.70") <= d_junc_8ej1 <= Decimal("3.90")) or not (Decimal("3.70") <= d_junc_8eiq <= Decimal("3.90")):
        raise AssertionError(f"deltaF508 (507, 509) peptide bond out of [3.70, 3.90] A: {d_junc_8ej1}, {d_junc_8eiq}")

    s_stab_q, noncovalent_pairs, min_dec_margin = compute_exact_fraction_s_stab(ca_8ej1, ca_8eiq)
    if min_dec_margin < Decimal("1e-9"):
        raise AssertionError(f"Minimum cutoff margin {min_dec_margin} < 1e-9 A")

    so3_drift, pairs_verified = verify_all_pairs_so3_q_zero_drift(ca_8ej1, ca_8eiq, common)
    if so3_drift != Fraction(0) or pairs_verified != noncovalent_pairs:
        raise AssertionError(f"SO(3, Q) drift check failed: drift={so3_drift}, pairs={pairs_verified}/{noncovalent_pairs}")

    # Positive controls in 50-digit Decimal & exact Fraction arithmetic
    d_site1_apo = dec_dist(ca_8ej1[464], ca_8ej1[1348])
    d_site1_eti = dec_dist(ca_8eiq[464], ca_8eiq[1348])
    d_site2_apo = dec_dist(ca_8ej1[582], ca_8ej1[1246])
    d_site2_eti = dec_dist(ca_8eiq[582], ca_8eiq[1246])
    if (d_site1_apo - d_site1_eti) <= Decimal("15.0") or d_site1_eti >= Decimal("12.0"):
        raise AssertionError("ATP Site 1 closure positive control failed in independent_audit")
    if (d_site2_apo - d_site2_eti) <= Decimal("15.0") or d_site2_eti >= Decimal("12.0"):
        raise AssertionError("ATP Site 2 closure positive control failed in independent_audit")

    for sc_codon in (38, 311, 361, 450, 1030, 1102):
        val = s_stab_q.get(sc_codon)
        if val is None or val <= Fraction(7, 10):
            raise AssertionError(f"Scaffold control {sc_codon} failed s_stab > 7/10: {val}")
    for atp_codon in (464, 582, 1246, 1347, 1348):
        val = s_stab_q.get(atp_codon)
        if val is None or val >= Fraction(14, 25):
            raise AssertionError(f"ATP-site control {atp_codon} failed s_stab < 14/25: {val}")

    holdout = json.loads(HOLDOUT_JSON.read_text(encoding="utf-8"))
    prov = holdout["provenance"]
    variants = holdout["variants"]
    v_sha = hashlib.sha256(json.dumps(variants, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()
    if v_sha != EXPECTED_HASHES["HOLDOUT_VARIANTS"]:
        raise AssertionError(f"Holdout variants_sha256 mismatch: {v_sha}")

    # Verify reversibility of CFTR-G970R auxiliary KALYDECO/SYMDEKO SPL erratum back to Commit 9a23e16 digest
    rev_variants = json.loads(json.dumps(variants))
    by_vid = {str(r["variant_id"]): r for r in rev_variants}
    for err in prov.get("post_audit_spl_errata", []):
        row_t = by_vid[str(err["variant_id"])]
        for fk, pair in err["fields_corrected"].items():
            old_v, new_v = pair
            if fk.startswith("label_status."):
                sk = fk.split(".", 1)[1]
                if row_t["label_status"][sk] != new_v:
                    raise AssertionError(f"Errata state mismatch on {fk}")
                row_t["label_status"][sk] = old_v
            else:
                if row_t[fk] != new_v:
                    raise AssertionError(f"Errata state mismatch on {fk}")
                row_t[fk] = old_v
    rev_sha = hashlib.sha256(json.dumps(rev_variants, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()
    if rev_sha != EXPECTED_HASHES["PRE_REGISTERED_HOLDOUT_VARIANTS"]:
        raise AssertionError(f"Pre-registered 9a23e16 digest reversal mismatch: {rev_sha}")

    miss_placed: list[tuple[Fraction | None, int]] = []
    miss_all14: list[tuple[Fraction | None, int]] = []
    miss_am: list[tuple[Fraction | None, int]] = []
    all24_g_stab: list[tuple[Fraction | None, int]] = []
    all24_g_am: list[tuple[Fraction | None, int]] = []
    all24_cls: list[tuple[Fraction | None, int]] = []

    for r in variants:
        y = int(r["trikafta_v23_responsive"])
        fam = r["variant_family"]
        codon = int(r["codon_position"])
        hgvs = str(r["hgvs_p_or_c"])
        if fam == "missense":
            sv = None if codon in (508, 1371) else s_stab_q.get(codon)
            gv = (Fraction(1, 2) + Fraction(1, 2) * sv) if sv is not None else None
            raw_am = r.get("alphamissense_score")
            s_am_q = (Fraction(1, 1) - Fraction(str(raw_am))) if raw_am is not None else None
            g_am_q = (Fraction(1, 2) + Fraction(1, 2) * s_am_q) if s_am_q is not None else None
            miss_all14.append((sv, y))
            miss_am.append((s_am_q, y))
            if sv is not None:
                miss_placed.append((sv, y))
            all24_g_stab.append((gv, y))
            all24_g_am.append((g_am_q, y))
            all24_cls.append((Fraction(3, 4), y))
        elif fam == "class_i_ptc_frameshift" or is_canonical_splice_charwise(hgvs):
            all24_g_stab.append((Fraction(0), y))
            all24_g_am.append((Fraction(0), y))
            all24_cls.append((Fraction(0), y))
        elif fam in ("inframe_indel", "splice_regulatory"):
            all24_g_stab.append((Fraction(1, 2), y))
            all24_g_am.append((Fraction(1, 2), y))
            all24_cls.append((Fraction(1, 2), y))
        else:
            raise AssertionError(f"Unsupported variant_family in independent_audit: {fam!r}")

    auc_placed = wilcoxon_auroc_fraction(miss_placed, abstention_aware=False)
    auc_all14 = wilcoxon_auroc_fraction(miss_all14, abstention_aware=True)
    auc_am = wilcoxon_auroc_fraction(miss_am, abstention_aware=True)
    auc_g_stab = wilcoxon_auroc_fraction(all24_g_stab, abstention_aware=True)
    auc_g_am = wilcoxon_auroc_fraction(all24_g_am, abstention_aware=True)
    auc_cls = wilcoxon_auroc_fraction(all24_cls, abstention_aware=True)

    if verify_parity:
        expected_counts = (1127, 1162, 1121, 626646)
        actual_counts = (len(ca_8ej1), len(ca_8eiq), len(common), pairs_verified)
        if actual_counts != expected_counts:
            raise AssertionError(f"Residue/pair count parity mismatch: {actual_counts} != {expected_counts}")
        expected_aucs = (
            Fraction(7, 10),
            Fraction(61, 98),
            Fraction(22, 49),
            Fraction(185, 256),
            Fraction(185, 256),
            Fraction(95, 128),
        )
        actual_aucs = (auc_placed, auc_all14, auc_am, auc_g_stab, auc_g_am, auc_cls)
        if actual_aucs != expected_aucs:
            raise AssertionError(f"Exact rational ROC-AUC parity mismatch: {actual_aucs} != {expected_aucs}")

    return {
        "placed_ca_8ej1": len(ca_8ej1),
        "placed_ca_8eiq": len(ca_8eiq),
        "shared_s_common": len(common),
        "noncovalent_pairs_verified": pairs_verified,
        "so3_q_drift": f"{so3_drift.numerator}/{so3_drift.denominator}",
        "min_decimal_cutoff_margin_angstrom": str(min_dec_margin),
        "placed_missense_auroc_rational": f"{auc_placed.numerator}/{auc_placed.denominator}",
        "placed_missense_auroc_float": round(float(auc_placed), 6),
        "missense_all14_abst_auroc_rational": f"{auc_all14.numerator}/{auc_all14.denominator}",
        "missense_all14_abst_auroc_float": round(float(auc_all14), 6),
        "standalone_am_auroc_rational": f"{auc_am.numerator}/{auc_am.denominator}",
        "standalone_am_auroc_float": round(float(auc_am), 6),
        "all24_g_stab_auroc_rational": f"{auc_g_stab.numerator}/{auc_g_stab.denominator}",
        "all24_g_stab_auroc_float": round(float(auc_g_stab), 6),
        "all24_g_am_auroc_rational": f"{auc_g_am.numerator}/{auc_g_am.denominator}",
        "all24_g_am_auroc_float": round(float(auc_g_am), 6),
        "all24_class_baseline_auroc_rational": f"{auc_cls.numerator}/{auc_cls.denominator}",
        "all24_class_baseline_auroc_float": round(float(auc_cls), 6),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Run exact rational audit and verify parity with Method 1")
    parser.add_argument("--json", action="store_true", help="Print JSON audit report")
    args = parser.parse_args(argv)

    rep = run_independent_audit(verify_parity=True if args.check else True)
    if args.json:
        print(json.dumps(rep, indent=2))
        return 0

    print(
        f"PASS independent_audit (PRED-0014): S_common={rep['shared_s_common']} | "
        f"pairs={rep['noncovalent_pairs_verified']} | SO(3,Q)_drift={rep['so3_q_drift']} | "
        f"placed_missense_AUC={rep['placed_missense_auroc_rational']} ({rep['placed_missense_auroc_float']:.6f}) | "
        f"all14_abst_AUC={rep['missense_all14_abst_auroc_rational']} ({rep['missense_all14_abst_auroc_float']:.6f}) | "
        f"all24_g_stab_AUC={rep['all24_g_stab_auroc_rational']} ({rep['all24_g_stab_auroc_float']:.6f})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
