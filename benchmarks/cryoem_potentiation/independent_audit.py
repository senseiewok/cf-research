#!/usr/bin/env python3
"""Independent Second-Method Verifier for Cryo-EM Potentiation Benchmark (PRED-0011).

Zero shared code, subprocesses, or imports from `cryoem_potentiation_benchmark.py`,
`tools/`, `subprocess`, `re`, or `shlex`:
- Tokenizes `8EIQ.cif.gz` via a character-level state-machine mmCIF parser.
- Parses canonical splice dinucleotides (`c.<digits>[+-][12]`) character-by-character without regex.
- Evaluates 3D distances in 50-digit `decimal.Decimal` arithmetic (`(dx^2 + dy^2 + dz^2).sqrt()`).
- Computes exact rational ROC-AUC via `fractions.Fraction` Wilcoxon mid-rank sum
      U = R_pos - |S_1|(|S_1| + 1) / 2
  and the closed-form abstention identity
      AUC_abst = 1/2 + (|S_1| * |S_0| / (N_1 * N_0)) * (AUC_cond - 1/2).
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
CIF_8EIQ = REPO_ROOT / "tools" / "cftr-model" / "data" / "8EIQ.cif.gz"
UNIPROT_JSON = REPO_ROOT / "tools" / "cftr-model" / "data" / "p13569.json"
ANNOTATIONS_JSON = REPO_ROOT / "benchmarks" / "theratype" / "theratype_open_annotations_2026-11-21.json"

AA_MAP = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def tokenize_mmcif_line(line: str) -> list[str]:
    """Character-by-character state-machine parser for single-line mmCIF records (no shlex/regex)."""
    lexemes: list[str] = []
    n = len(line)
    idx = 0
    while idx < n:
        ch = line[idx]
        if ch in (" ", "\t", "\r", "\n"):
            idx += 1
            continue
        if ch in ("'", '"'):
            quote = ch
            idx += 1
            buf: list[str] = []
            while idx < n:
                if line[idx] == quote and (idx + 1 == n or line[idx + 1] in (" ", "\t", "\r", "\n")):
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
    """Charwise check for c.<digits>[+-][12] followed by non-digit or end-of-string (no regex)."""
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


def decimal_dist(p: tuple[Decimal, Decimal, Decimal], q: tuple[Decimal, Decimal, Decimal]) -> Decimal:
    dx = p[0] - q[0]
    dy = p[1] - q[1]
    dz = p[2] - q[2]
    return (dx * dx + dy * dy + dz * dz).sqrt()


def verify_ca_against_uniprot_strict(
    ca_decimal: dict[int, tuple[str, tuple[Decimal, Decimal, Decimal]]],
    uniprot_seq: str,
    declared_mutations: dict[int, str],
    declared_deletions: set[int],
) -> None:
    for del_pos in declared_deletions:
        if del_pos in ca_decimal:
            raise AssertionError(f"Deleted construct residue {del_pos} unexpectedly present in Chain A")
    for res_no, (comp, _) in sorted(ca_decimal.items()):
        if res_no < 1 or res_no > len(uniprot_seq):
            raise AssertionError(f"Residue index {res_no} out of P13569 range 1..{len(uniprot_seq)}")
        if res_no in declared_mutations:
            if comp != declared_mutations[res_no]:
                raise AssertionError(
                    f"Construct mutation mismatch at {res_no}: {comp} != declared {declared_mutations[res_no]}"
                )
            continue
        if AA_MAP.get(comp) != uniprot_seq[res_no - 1]:
            raise AssertionError(f"UniProt mismatch at {res_no}: {comp} vs {uniprot_seq[res_no - 1]}")


def wilcoxon_rank_sum_auroc(
    labels: list[int],
    scores: list[Decimal | Fraction | None],
) -> tuple[Fraction | None, Fraction, int]:
    """Exact rational ROC-AUC via Wilcoxon mid-rank sum and closed-form abstention identity."""
    if len(labels) != len(scores) or not labels:
        raise ValueError("Labels and scores must be non-empty and equal length")
    for y in labels:
        if isinstance(y, bool) or y not in (0, 1):
            raise ValueError(f"Labels must be binary integers in {{0, 1}}, got {y!r}")

    n_pos = sum(1 for y in labels if y == 1)
    n_neg = sum(1 for y in labels if y == 0)
    if n_pos == 0 or n_neg == 0:
        raise ValueError("Need positive and negative items")

    scored: list[tuple[Fraction, int]] = []
    for y, s in zip(labels, scores):
        if s is not None:
            if isinstance(s, bool) or (isinstance(s, Decimal) and not s.is_finite()):
                raise ValueError(f"Non-finite score {s!r}")
            f_val = s if isinstance(s, Fraction) else Fraction(str(s))
            scored.append((f_val, y))

    s_pos = sum(1 for _, y in scored if y == 1)
    s_neg = sum(1 for _, y in scored if y == 0)
    if s_pos == 0 or s_neg == 0:
        return None, Fraction(1, 2), len(scored)

    scored.sort(key=lambda item: item[0])
    m = len(scored)
    rank_sum_pos = Fraction(0, 1)
    idx = 0
    while idx < m:
        j = idx + 1
        while j < m and scored[j][0] == scored[idx][0]:
            j += 1
        mid_rank = Fraction(idx + 1 + j, 2)
        pos_in_tie = sum(1 for k in range(idx, j) if scored[k][1] == 1)
        rank_sum_pos += mid_rank * pos_in_tie
        idx = j

    u_stat = rank_sum_pos - Fraction(s_pos * (s_pos + 1), 2)
    auc_cond = u_stat / Fraction(s_pos * s_neg, 1)
    weight = Fraction(s_pos * s_neg, n_pos * n_neg)
    auc_abst = Fraction(1, 2) + weight * (auc_cond - Fraction(1, 2))
    return auc_cond, auc_abst, m


def run_independent_audit() -> dict:
    cif_sha = hashlib.sha256(CIF_8EIQ.read_bytes()).hexdigest()
    uni_sha = hashlib.sha256(UNIPROT_JSON.read_bytes()).hexdigest()
    if cif_sha != "d13cbfc8ba6575ab6228f351578c6f9cc5b9dafe13e5913d9b27b3cacebc18d6":
        raise AssertionError(f"Unexpected 8EIQ SHA-256: {cif_sha}")
    if uni_sha != "c8feabfabe8690ef4a06af7de5b5f8dc3c8da554246c5ab3d403edf265d4a893":
        raise AssertionError(f"Unexpected P13569 SHA-256: {uni_sha}")

    uniprot_seq = json.loads(UNIPROT_JSON.read_text(encoding="utf-8"))["sequence"]["value"]
    raw_lines = gzip.open(CIF_8EIQ, "rt", encoding="utf-8").read().splitlines()

    diff_cols: list[str] = []
    pos = 0
    while pos < len(raw_lines) and not raw_lines[pos].startswith("_struct_ref_seq_dif."):
        pos += 1
    while pos < len(raw_lines) and raw_lines[pos].startswith("_struct_ref_seq_dif."):
        diff_cols.append(raw_lines[pos].strip()[len("_struct_ref_seq_dif."):])
        pos += 1
    declared_mutations: dict[int, str] = {}
    declared_deletions: set[int] = set()
    while diff_cols and pos < len(raw_lines) and not raw_lines[pos].startswith(("#", "_", "loop_")):
        toks = tokenize_mmcif_line(raw_lines[pos])
        pos += 1
        if len(toks) == len(diff_cols):
            dmap = dict(zip(diff_cols, toks))
            if dmap.get("pdbx_pdb_strand_id") != "A" or dmap.get("pdbx_seq_db_accession_code") != "P13569":
                continue
            if dmap.get("details") == "expression tag":
                continue
            seq_str = dmap.get("pdbx_seq_db_seq_num", "")
            if seq_str.isdigit():
                p_num = int(seq_str)
                mon = dmap.get("mon_id", "?")
                if mon in ("?", "."):
                    declared_deletions.add(p_num)
                else:
                    declared_mutations[p_num] = mon

    atom_cols: list[str] = []
    ca_decimal: dict[int, tuple[str, tuple[Decimal, Decimal, Decimal]]] = {}
    a_pot: list[tuple[Decimal, Decimal, Decimal]] = []
    a_fold: list[tuple[Decimal, Decimal, Decimal]] = []

    in_atoms = False
    col_ix: dict[str, int] = {}
    for ln in raw_lines:
        if ln.startswith("_atom_site."):
            atom_cols.append(ln.strip()[len("_atom_site."):])
            in_atoms = True
            continue
        if in_atoms and not col_ix and atom_cols:
            col_ix = {name: k for k, name in enumerate(atom_cols)}
        if in_atoms and (ln.startswith("ATOM ") or ln.startswith("HETATM ")):
            t = tokenize_mmcif_line(ln)
            if t[col_ix["pdbx_PDB_model_num"]] != "1" or t[col_ix["auth_asym_id"]] != "A":
                continue
            if t[col_ix["label_alt_id"]] not in (".", "A"):
                continue
            if t[col_ix["type_symbol"]].upper() == "H":
                continue
            pt = (
                Decimal(t[col_ix["Cartn_x"]]),
                Decimal(t[col_ix["Cartn_y"]]),
                Decimal(t[col_ix["Cartn_z"]]),
            )
            grp = t[col_ix["group_PDB"]]
            comp = t[col_ix["label_comp_id"]]
            atom_name = t[col_ix["label_atom_id"]]
            s_raw = t[col_ix["auth_seq_id"]]
            seq_num = int(s_raw) if s_raw.lstrip("-").isdigit() else None
            if grp == "HETATM":
                if comp in ("ATP", "MG", "VX7"):
                    a_pot.append(pt)
                elif comp in ("CV6", "WJX"):
                    a_fold.append(pt)
            elif grp == "ATOM" and seq_num is not None:
                if atom_name == "CA" and seq_num not in ca_decimal:
                    ca_decimal[seq_num] = (comp, pt)
                if seq_num in (1063, 1070):
                    a_fold.append(pt)
        elif in_atoms and ca_decimal:
            break

    verify_ca_against_uniprot_strict(ca_decimal, uniprot_seq, declared_mutations, declared_deletions)
    if len(a_pot) != 93 or len(a_fold) != 103 or 508 in ca_decimal:
        raise AssertionError(f"Anchor/construct check failed: pot={len(a_pot)}, fold={len(a_fold)}")

    annot = json.loads(ANNOTATIONS_JSON.read_text(encoding="utf-8"))
    variants = annot["variants"]
    train_codons = [int(v["codon_position"]) for v in variants if v["split"] == "train"]
    test_codons = [int(v["codon_position"]) for v in variants if v["split"] == "test"]
    min_codon_gap = min(abs(a - b) for a in train_codons for b in test_codons)

    test_variants = [v for v in variants if v["split"] == "test"]
    y_test = [int(v["kalydeco_v36_responsive"]) for v in test_variants]

    s3d_all: list[Decimal | Fraction | None] = []
    sam_all: list[Decimal | Fraction | None] = []
    seq_all: list[Decimal | Fraction | None] = []
    g3d_all: list[Decimal | Fraction | None] = []
    gam_all: list[Decimal | Fraction | None] = []

    miss_y: list[int] = []
    miss_s3d: list[Decimal | Fraction | None] = []
    miss_sam: list[Decimal | Fraction | None] = []

    half = Decimal("0.5")
    for v in test_variants:
        fam = v["variant_family"]
        codon = int(v["codon_position"])
        hgvs = str(v["hgvs_p_or_c"])
        if fam == "class_i_ptc_frameshift":
            seq_val = Fraction(0, 1)
        elif fam == "splice_regulatory":
            seq_val = Fraction(0, 1) if is_canonical_splice_charwise(hgvs) else Fraction(1, 2)
        elif fam == "inframe_indel":
            seq_val = Fraction(1, 2)
        elif fam == "missense":
            seq_val = Fraction(3, 4)
        else:
            raise AssertionError(f"Unexpected family {fam}")
        seq_all.append(seq_val)

        if fam == "missense":
            miss_y.append(int(v["kalydeco_v36_responsive"]))
            am_raw = v.get("alphamissense_score")
            s_am = (Decimal("1") - Decimal(str(am_raw))) if am_raw is not None else None
            if codon in ca_decimal:
                ca_pt = ca_decimal[codon][1]
                dp = min(decimal_dist(ca_pt, q) for q in a_pot)
                df = min(decimal_dist(ca_pt, q) for q in a_fold)
                s_3d: Decimal | None = df / (dp + df)
            else:
                s_3d = None
            s3d_all.append(s_3d)
            sam_all.append(s_am)
            g3d_all.append((half + half * s_3d) if s_3d is not None else None)
            gam_all.append((half + half * s_am) if s_am is not None else None)
            miss_s3d.append(s_3d)
            miss_sam.append(s_am)
        else:
            s3d_all.append(None)
            sam_all.append(None)
            g3d_all.append(seq_val)
            gam_all.append(seq_val)

    cond_s3d, _, n_miss_scored = wilcoxon_rank_sum_auroc(miss_y, miss_s3d)
    cond_sam, _, _ = wilcoxon_rank_sum_auroc(miss_y, miss_sam)
    _, abst_s3d_ungated, _ = wilcoxon_rank_sum_auroc(y_test, s3d_all)
    _, abst_seq, _ = wilcoxon_rank_sum_auroc(y_test, seq_all)
    _, abst_g3d, _ = wilcoxon_rank_sum_auroc(y_test, g3d_all)
    _, abst_gam, _ = wilcoxon_rank_sum_auroc(y_test, gam_all)

    assert cond_s3d is not None and cond_sam is not None
    return {
        "placed_chain_a_residues": len(ca_decimal),
        "potentiation_heavy_atoms": len(a_pot),
        "folding_heavy_atoms": len(a_fold),
        "min_cross_split_codon_gap": min_codon_gap,
        "test_total": len(test_variants),
        "test_kalydeco_responsive": sum(y_test),
        "test_missense_scored": n_miss_scored,
        "exact_fractions": {
            "missense_conditional_auroc_s3d": f"{cond_s3d.numerator}/{cond_s3d.denominator}",
            "missense_conditional_auroc_sam": f"{cond_sam.numerator}/{cond_sam.denominator}",
            "all_ungated_abstention_auroc_s3d": f"{abst_s3d_ungated.numerator}/{abst_s3d_ungated.denominator}",
            "all_abstention_auroc_seq_baseline": f"{abst_seq.numerator}/{abst_seq.denominator}",
            "all_abstention_auroc_g3d": f"{abst_g3d.numerator}/{abst_g3d.denominator}",
            "all_abstention_auroc_gam": f"{abst_gam.numerator}/{abst_gam.denominator}",
        },
        "values": {
            "missense_conditional_auroc_s3d": round(float(cond_s3d), 6),
            "missense_conditional_auroc_sam": round(float(cond_sam), 6),
            "all_ungated_abstention_auroc_s3d": round(float(abst_s3d_ungated), 6),
            "all_abstention_auroc_seq_baseline": round(float(abst_seq), 6),
            "all_abstention_auroc_g3d": round(float(abst_g3d), 6),
            "all_abstention_auroc_gam": round(float(abst_gam), 6),
        },
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Independent second-method verifier for PRED-0011")
    ap.add_argument("--check", action="store_true", help="Verify exact rational invariants independently")
    ap.add_argument("--json", action="store_true", help="Print JSON report")
    args = ap.parse_args(argv)

    audit = run_independent_audit()
    if args.check:
        if audit["placed_chain_a_residues"] != 1162:
            raise AssertionError(f"Expected 1162 placed Chain A residues, got {audit['placed_chain_a_residues']}")
        if audit["potentiation_heavy_atoms"] != 93 or audit["folding_heavy_atoms"] != 103:
            raise AssertionError("Anchor heavy-atom count invariant failed")
        if audit["min_cross_split_codon_gap"] <= 5:
            raise AssertionError("Cross-split codon gap <= 5")
        if audit["test_kalydeco_responsive"] != 7 or audit["test_total"] != 20:
            raise AssertionError("Unexpected test split KALYDECO label distribution")

    if args.json:
        print(json.dumps(audit, indent=2))
        return 0

    v = audit["values"]
    ef = audit["exact_fractions"]
    print(
        f"PASS independent_audit: KALYDECO test responsive={audit['test_kalydeco_responsive']}/{audit['test_total']}, "
        f"missense s_3D={v['missense_conditional_auroc_s3d']:.4f} ({ef['missense_conditional_auroc_s3d']}) vs "
        f"s_AM={v['missense_conditional_auroc_sam']:.4f} ({ef['missense_conditional_auroc_sam']}); "
        f"g_3D={v['all_abstention_auroc_g3d']:.4f} ({ef['all_abstention_auroc_g3d']}) vs "
        f"seq={v['all_abstention_auroc_seq_baseline']:.4f} ({ef['all_abstention_auroc_seq_baseline']}), "
        f"g_AM={v['all_abstention_auroc_gam']:.4f} ({ef['all_abstention_auroc_gam']})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
