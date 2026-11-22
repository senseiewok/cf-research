#!/usr/bin/env python3
"""Independent Second-Method Verifier for Cryo-EM Gating Kinematics Benchmark (PRED-0012).

Zero shared code, subprocesses, or imports from `cryoem_gating_kinematics_benchmark.py`,
`tools/`, `subprocess`, `re`, or `shlex`:
- Tokenizes `5UAK.cif.gz` and `6MSM.cif.gz` via a character-by-character mmCIF state-machine lexer.
- Parses canonical splice dinucleotides (`c.<digits>[+-][12]`) character-by-character without regex.
- Computes all-pairs Euclidean distances in 50-digit `decimal.Decimal` arithmetic
      d(p, q) = (dx^2 + dy^2 + dz^2).sqrt()
- Evaluates Mariani et al. (2013) conformational reorganization in EXACT rational arithmetic:
      s_kin(i) = Fraction(4 * |N_sym(i)| - sum_{j in N_sym(i)} c_ij, 4 * |N_sym(i)|)
- Computes exact rational ROC-AUC via `fractions.Fraction` Wilcoxon mid-rank sum:
      U = R_pos - |S_1|(|S_1| + 1) / 2
"""
from __future__ import annotations

import argparse
import decimal
import fractions
import gzip
import hashlib
import json
from decimal import Decimal
from fractions import Fraction
from pathlib import Path

decimal.getcontext().prec = 50

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
CIF_5UAK = REPO_ROOT / "tools" / "cftr-model" / "data" / "5UAK.cif.gz"
CIF_6MSM = REPO_ROOT / "tools" / "cftr-model" / "data" / "6MSM.cif.gz"
UNIPROT_JSON = REPO_ROOT / "tools" / "cftr-model" / "data" / "p13569.json"
ANNOTATIONS_JSON = REPO_ROOT / "benchmarks" / "theratype" / "theratype_open_annotations_2026-11-21.json"

EXPECTED_HASHES = {
    "5UAK": "8bf00f2612dd961d3b1ccd2f2193781a534ad44f32666a6de8da948afb8d84cf",
    "6MSM": "82933988abc3ed24a343819d410b39dc0d86913eaba4426dd945e3273e34fa59",
    "P13569": "c8feabfabe8690ef4a06af7de5b5f8dc3c8da554246c5ab3d403edf265d4a893",
}

AA_3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

R0_DEC = Decimal("15.0")
THRESH_DEC = (Decimal("0.5"), Decimal("1.0"), Decimal("2.0"), Decimal("4.0"))


def lex_mmcif_line(line: str) -> list[str]:
    """Character-by-character state-machine lexer for mmCIF lines (no shlex or re)."""
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
    if not hgvs.startswith("c.") or len(hgvs) < 5:
        return False
    pos = 2
    while pos < len(hgvs) and hgvs[pos].isdigit():
        pos += 1
    if pos == 2 or pos >= len(hgvs):
        return False
    if hgvs[pos] not in ("+", "-"):
        return False
    pos += 1
    if pos >= len(hgvs) or hgvs[pos] not in ("1", "2"):
        return False
    pos += 1
    return pos == len(hgvs) or not hgvs[pos].isdigit()


def dec_dist(p: tuple[Decimal, Decimal, Decimal], q: tuple[Decimal, Decimal, Decimal]) -> Decimal:
    dx = p[0] - q[0]
    dy = p[1] - q[1]
    dz = p[2] - q[2]
    return (dx * dx + dy * dy + dz * dz).sqrt()


def load_decimal_ca(
    cif_path: Path,
    tag: str,
    uniprot_seq: str,
) -> dict[int, tuple[Decimal, Decimal, Decimal]]:
    digest = hashlib.sha256(cif_path.read_bytes()).hexdigest()
    if digest != EXPECTED_HASHES[tag]:
        raise AssertionError(f"{tag} SHA-256 mismatch: {digest}")

    lines = gzip.open(cif_path, "rt", encoding="utf-8").read().splitlines()

    # Parse _struct_ref_seq_dif (loop or single record)
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
        else:
            rec: dict[str, str] = {}
            while i < len(lines) and lines[i].startswith("_struct_ref_seq_dif."):
                toks = lex_mmcif_line(lines[i])
                k = toks[0][len("_struct_ref_seq_dif."):]
                if len(toks) >= 2:
                    rec[k] = toks[1]
                i += 1
            if rec:
                diff_records.append(rec)

        for r in diff_records:
            if r.get("pdbx_pdb_strand_id") != "A" or r.get("pdbx_seq_db_accession_code") != "P13569":
                continue
            det = (r.get("details") or "").lower()
            if "expression tag" in det or "cloning" in det:
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

    atom_cols: list[str] = []
    ca_map: dict[int, tuple[str, tuple[Decimal, Decimal, Decimal]]] = {}
    in_atom = False
    for ln in lines:
        if ln.startswith("_atom_site."):
            atom_cols.append(ln.strip()[len("_atom_site."):])
            in_atom = True
        elif in_atom and ln.startswith(("ATOM ", "HETATM ")):
            toks = lex_mmcif_line(ln)
            col = {name: toks[idx] for idx, name in enumerate(atom_cols)}
            if col["group_PDB"] != "ATOM":
                continue
            if col["pdbx_PDB_model_num"] != "1":
                continue
            if col["auth_asym_id"] != "A" or col["label_asym_id"] != "A":
                continue
            if col["pdbx_PDB_ins_code"] not in ("?", "."):
                continue
            if col["label_alt_id"] not in (".", "A"):
                continue
            if col["label_atom_id"] != "CA":
                continue
            raw_seq = col["auth_seq_id"]
            if not (raw_seq.isdigit() or (raw_seq.startswith("-") and raw_seq[1:].isdigit())):
                continue
            seq_no = int(raw_seq)
            if seq_no < 1 or seq_no > 1480:
                continue
            if seq_no in ca_map:
                raise AssertionError(f"{tag}: duplicate CA at {seq_no}")
            ca_map[seq_no] = (
                col["label_comp_id"],
                (Decimal(col["Cartn_x"]), Decimal(col["Cartn_y"]), Decimal(col["Cartn_z"])),
            )
        elif in_atom and ca_map and not ln.startswith(("ATOM ", "HETATM ")):
            break

    for d_pos in deletions:
        if d_pos in ca_map:
            raise AssertionError(f"{tag}: deleted residue {d_pos} present")
    for pos, (comp, _) in sorted(ca_map.items()):
        if pos in mutations:
            if comp != mutations[pos]:
                raise AssertionError(f"{tag}: mutation mismatch at {pos}: {comp} != {mutations[pos]}")
        else:
            if AA_3_TO_1.get(comp) != uniprot_seq[pos - 1]:
                raise AssertionError(f"{tag}: P13569 mismatch at {pos}: {comp} vs {uniprot_seq[pos - 1]}")

    return {pos: xyz for pos, (_, xyz) in ca_map.items()}


def compute_exact_fraction_s_kin(
    ca_5uak: dict[int, tuple[Decimal, Decimal, Decimal]],
    ca_6msm: dict[int, tuple[Decimal, Decimal, Decimal]],
) -> tuple[dict[int, Fraction | None], dict[int, int]]:
    """Brute-force O(|S_common|^2) 50-digit Decimal distances -> exact Fraction s_kin(i)."""
    common = sorted(set(ca_5uak.keys()) & set(ca_6msm.keys()))
    n = len(common)
    # Precompute symmetric upper-triangular Decimal distance matrices on S_common
    d5: list[list[Decimal]] = [[Decimal(0)] * n for _ in range(n)]
    d6: list[list[Decimal]] = [[Decimal(0)] * n for _ in range(n)]
    for a in range(n):
        ra = common[a]
        p5 = ca_5uak[ra]
        p6 = ca_6msm[ra]
        for b in range(a + 1, n):
            rb = common[b]
            dist5 = dec_dist(p5, ca_5uak[rb])
            dist6 = dec_dist(p6, ca_6msm[rb])
            d5[a][b] = d5[b][a] = dist5
            d6[a][b] = d6[b][a] = dist6

    s_kin_frac: dict[int, Fraction | None] = {}
    deg_counts: dict[int, int] = {}
    for a in range(n):
        ri = common[a]
        preserved = 0
        count_neighbors = 0
        for b in range(n):
            rj = common[b]
            if abs(ri - rj) <= 1:
                continue
            da = d5[a][b]
            db = d6[a][b]
            m_dist = da if da <= db else db
            if m_dist <= R0_DEC:
                count_neighbors += 1
                diff = abs(da - db)
                for t in THRESH_DEC:
                    if diff < t:
                        preserved += 1
        deg_counts[ri] = count_neighbors
        if count_neighbors == 0:
            s_kin_frac[ri] = None
        else:
            denom = 4 * count_neighbors
            s_kin_frac[ri] = Fraction(denom - preserved, denom)

    return s_kin_frac, deg_counts


def wilcoxon_rank_sum_auroc(
    labels: list[int],
    scores: list[Fraction | Decimal | None],
) -> tuple[Fraction | None, Fraction]:
    if len(labels) != len(scores) or not labels:
        raise ValueError("Mismatched or empty labels/scores")
    for y in labels:
        if isinstance(y, bool) or y not in (0, 1):
            raise ValueError(f"Non-binary label {y!r}")

    n1 = sum(1 for y in labels if y == 1)
    n0 = sum(1 for y in labels if y == 0)
    if n1 == 0 or n0 == 0:
        raise ValueError("Both classes required")

    scored: list[tuple[Fraction, int]] = []
    for y, s in zip(labels, scores):
        if s is not None:
            f_val = s if isinstance(s, Fraction) else Fraction(str(s))
            scored.append((f_val, y))

    s1 = sum(1 for _, y in scored if y == 1)
    s0 = sum(1 for _, y in scored if y == 0)
    if s1 == 0 or s0 == 0:
        return None, Fraction(1, 2)

    scored.sort(key=lambda item: item[0])
    m = len(scored)
    r_pos = Fraction(0, 1)
    idx = 0
    while idx < m:
        end = idx + 1
        while end < m and scored[end][0] == scored[idx][0]:
            end += 1
        # 1-based mid-rank = ((idx + 1) + end) / 2
        mid_rank = Fraction((idx + 1) + end, 2)
        pos_in_tie = sum(1 for k in range(idx, end) if scored[k][1] == 1)
        r_pos += mid_rank * pos_in_tie
        idx = end

    u_stat = r_pos - Fraction(s1 * (s1 + 1), 2)
    auc_cond = u_stat / Fraction(s1 * s0, 1)
    auc_abst = Fraction(1, 2) + Fraction(s1 * s0, n1 * n0) * (auc_cond - Fraction(1, 2))
    return auc_cond, auc_abst


def class_baseline_fraction(v: dict[str, object]) -> Fraction:
    fam = str(v["variant_family"])
    hgvs = str(v["hgvs_p_or_c"])
    if fam == "class_i_ptc_frameshift":
        return Fraction(0, 1)
    if fam == "splice_regulatory":
        return Fraction(0, 1) if is_canonical_splice_charwise(hgvs) else Fraction(1, 2)
    if fam == "inframe_indel":
        return Fraction(1, 2)
    if fam == "missense":
        return Fraction(3, 4)
    raise ValueError(f"Unknown family {fam!r}")


def run_independent_audit() -> dict[str, object]:
    p_hash = hashlib.sha256(UNIPROT_JSON.read_bytes()).hexdigest()
    if p_hash != EXPECTED_HASHES["P13569"]:
        raise AssertionError("P13569 SHA-256 mismatch")
    uniprot_seq = str(json.loads(UNIPROT_JSON.read_text(encoding="utf-8"))["sequence"]["value"])

    ca_5uak = load_decimal_ca(CIF_5UAK, "5UAK", uniprot_seq)
    ca_6msm = load_decimal_ca(CIF_6MSM, "6MSM", uniprot_seq)
    s_kin_frac, deg_counts = compute_exact_fraction_s_kin(ca_5uak, ca_6msm)

    if len(ca_5uak) != 1139 or len(ca_6msm) != 1181 or len(s_kin_frac) != 1130:
        raise AssertionError(f"Unexpected counts: 5UAK={len(ca_5uak)}, 6MSM={len(ca_6msm)}, common={len(s_kin_frac)}")

    ann = json.loads(ANNOTATIONS_JSON.read_text(encoding="utf-8"))
    test_rows = [v for v in ann["variants"] if v["split"] == "test"]
    labels = [int(v["kalydeco_v36_responsive"]) for v in test_rows]

    missense_rows = [v for v in test_rows if v["variant_family"] == "missense"]
    m_labels = [int(v["kalydeco_v36_responsive"]) for v in missense_rows]
    m_skin = [s_kin_frac.get(int(v["codon_position"])) for v in missense_rows]
    m_sam = [
        (Fraction(1, 1) - Fraction(str(v["alphamissense_score"])))
        if v["alphamissense_score"] is not None
        else None
        for v in missense_rows
    ]

    g_kin: list[Fraction | None] = []
    g_am: list[Fraction | None] = []
    base_scores: list[Fraction | None] = []
    for v in test_rows:
        fam = str(v["variant_family"])
        b_val = class_baseline_fraction(v)
        base_scores.append(b_val)
        if fam != "missense":
            g_kin.append(b_val)
            g_am.append(b_val)
        else:
            sk = s_kin_frac.get(int(v["codon_position"]))
            g_kin.append(Fraction(1, 2) + Fraction(1, 2) * sk if sk is not None else None)
            am_raw = v["alphamissense_score"]
            g_am.append(
                Fraction(1, 2) + Fraction(1, 2) * (Fraction(1, 1) - Fraction(str(am_raw)))
                if am_raw is not None
                else None
            )

    cond_skin_q, _ = wilcoxon_rank_sum_auroc(m_labels, m_skin)
    cond_sam_q, _ = wilcoxon_rank_sum_auroc(m_labels, m_sam)
    _, abst_gkin_q = wilcoxon_rank_sum_auroc(labels, g_kin)
    _, abst_gam_q = wilcoxon_rank_sum_auroc(labels, g_am)
    _, abst_base_q = wilcoxon_rank_sum_auroc(labels, base_scores)

    assert cond_skin_q is not None and cond_sam_q is not None
    supported = bool(
        cond_skin_q >= Fraction(3, 4)
        and cond_skin_q > cond_sam_q
        and abst_gkin_q >= Fraction(3, 4)
        and abst_gkin_q > abst_base_q
        and abst_gkin_q > abst_gam_q
    )

    return {
        "n_5uak_ca": len(ca_5uak),
        "n_6msm_ca": len(ca_6msm),
        "n_common_ca": len(s_kin_frac),
        "test_missense_s_kin_fraction": f"{cond_skin_q.numerator}/{cond_skin_q.denominator}",
        "test_missense_s_kin_float": round(float(cond_skin_q), 6),
        "test_missense_s_am_fraction": f"{cond_sam_q.numerator}/{cond_sam_q.denominator}",
        "test_missense_s_am_float": round(float(cond_sam_q), 6),
        "test_gated_g_kin_fraction": f"{abst_gkin_q.numerator}/{abst_gkin_q.denominator}",
        "test_gated_g_kin_float": round(float(abst_gkin_q), 6),
        "test_gated_g_am_fraction": f"{abst_gam_q.numerator}/{abst_gam_q.denominator}",
        "test_gated_g_am_float": round(float(abst_gam_q), 6),
        "test_baseline_fraction": f"{abst_base_q.numerator}/{abst_base_q.denominator}",
        "test_baseline_float": round(float(abst_base_q), 6),
        "held_out_error_float": round(float(Fraction(1, 1) - abst_gkin_q), 6),
        "outcome": "supported" if supported else "refuted",
        "missense_test_details": [
            {
                "legacy_spl_name": v["legacy_spl_name"],
                "codon": int(v["codon_position"]),
                "kalydeco_v36": int(v["kalydeco_v36_responsive"]),
                "n_sym": deg_counts.get(int(v["codon_position"])),
                "s_kin_fraction": (
                    f"{s_kin_frac[int(v['codon_position'])].numerator}/{s_kin_frac[int(v['codon_position'])].denominator}"
                    if s_kin_frac.get(int(v["codon_position"])) is not None
                    else None
                ),
                "s_kin_float": (
                    round(float(s_kin_frac[int(v["codon_position"])]), 6)
                    if s_kin_frac.get(int(v["codon_position"])) is not None
                    else None
                ),
            }
            for v in missense_rows
        ],
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Independent Decimal/Fraction verifier for PRED-0012.")
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args(argv)

    rep = run_independent_audit()
    if args.json:
        print(json.dumps(rep, indent=2))
        return 0
    print(
        f"PASS independent_audit (PRED-0012 -> outcome={rep['outcome']}): "
        f"S_common={rep['n_common_ca']}; "
        f"s_kin AUC_cond={rep['test_missense_s_kin_float']:.6f} ({rep['test_missense_s_kin_fraction']}); "
        f"g_kin AUC_abst={rep['test_gated_g_kin_float']:.6f} ({rep['test_gated_g_kin_fraction']})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
