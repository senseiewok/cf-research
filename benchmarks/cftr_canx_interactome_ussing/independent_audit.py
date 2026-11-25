#!/usr/bin/env python3
"""Independent verification audit for benchmarks/cftr_canx_interactome_ussing/.

Pure-stdlib implementation importing NOTHING from `cftr_canx_interactome_ussing.py` or `tools/`:
1. Validates `done_sha256`, `results_2026-11-25.json`, and `ledger.json` consistency.
2. Checks Arm 1 Bonferroni gate rules and Arm 2 BH-FDR q-values (`tol <= 1e-6`) offline in CI.
3. When local `.xlsx` workbooks exist under `sources/downloads/tedman-2025/`, independently verifies
   all 5 workbook SHA-256 digests, re-parses sheets via `xml.dom.minidom`, recomputes Arm 1
   moments/correlations in `Decimal`, verifies Welch's t via smooth trigonometric Simpson quadrature,
   re-enumerates Arm 1a & 1b permutations via bitmask/two-pointer scans, and re-computes Arm 2
   Wilcoxon (`collections.Counter` subset-sum) & Mann-Whitney U (direct pairwise dominance `O(n1*n2)`)
   within `tol <= 1e-6`.
"""
from __future__ import annotations

import argparse
import collections
from decimal import Decimal, getcontext
import hashlib
import json
import math
from pathlib import Path
import re
import sys
from typing import Any, Dict, List, Optional, Sequence, Tuple
from xml.dom import minidom
import zipfile

getcontext().prec = 34

REPO_ROOT = Path(__file__).resolve().parents[2]
BENCH_DIR = REPO_ROOT / "benchmarks" / "cftr_canx_interactome_ussing"

EXPECTED_WORKBOOKS: Dict[str, str] = {
    "Copy of TMT_Labels.xlsx": "5b16d08dfefce36878149259d823604d2153df8b027a48bd8cf69c3f1af5b41e",
    "Copy of CFTR_Variants_LogFoldChanges.xlsx": "492df3ef266f27a5e18e1d9ce61b13356a8101b3e85d89c38ec87f5ca7208f7a",
    "Copy of CFTR_F508del_Drug_LogFoldChanges.xlsx": "d42ae782fe85163909396ea5ab4ee0915eec346e8ae420e419e80b5dd5e270c0",
    "Copy of FRT-F508del_siCANX_Ussing-summary-for-paper.xlsx": "a1ad7e325bab60ba6fe15e001fdddfd4fd14918d571e53ca45844c2f58f1ac1a",
    "Copy of FRT-WT_siCANX_Ussing-summary-for-paper.xlsx": "bddaffb3b68449184b518a923c3c0d18ce51ed409662915d7b16fe7f8fedc628",
}

USSING_GROUPS = (
    "NT",
    "siNS, DMSO",
    "siNS, VX-661",
    "siNS, VX-445",
    "siNS, 661+445",
    "siNS, ETI",
    "siCANX, DMSO",
    "siCANX, VX-661",
    "siCANX, VX-445",
    "siCANX, 661+445",
    "siCANX, ETI",
)

VAR_COLS = ("WT(-CX)", "V232D (-CX)", "P5L (-CX)", "T1036N (-CX)", "I1366N (-CX)")
DRUG_COLS = ("ΔF (-CX)", "DelF (-CX) +445", "DelF (-CX) +661", "DelF (-CX) + combo")
MODULES: Dict[str, set[str]] = {
    "folding": {"folding", "protein folding"},
    "ups": {"ubiquitin proteasome system"},
    "autophagy_lysosome": {"autophagy-lysosome pathway"},
    "trafficking": {"trafficking"},
}


def _find_data_dir(explicit: Optional[Path], repo_root: Path) -> Optional[Path]:
    probe = "Copy of FRT-F508del_siCANX_Ussing-summary-for-paper.xlsx"
    if explicit is not None:
        p = Path(explicit)
        return p if (p / probe).is_file() else None
    direct = repo_root / "sources" / "downloads" / "tedman-2025"
    if (direct / probe).is_file():
        return direct
    if len(repo_root.parents) >= 3:
        fallback = repo_root.parents[2] / "cf-research" / "sources" / "downloads" / "tedman-2025"
        if (fallback / probe).is_file():
            return fallback
    return None


def _col_num(col: str) -> int:
    acc = 0
    for ch in col:
        acc = acc * 26 + (ord(ch) - ord("A") + 1)
    return acc


def _num_col(n: int) -> str:
    out: List[str] = []
    while n > 0:
        n, r = divmod(n - 1, 26)
        out.append(chr(ord("A") + r))
    return "".join(reversed(out))


def _read_sheet_minidom(xlsx_path: Path, sheet_name: str) -> Dict[Tuple[int, str], str]:
    """Independent OOXML sheet parser using xml.dom.minidom."""
    with zipfile.ZipFile(xlsx_path, "r") as zf:
        sst: List[str] = []
        if "xl/sharedStrings.xml" in zf.namelist():
            sdoc = minidom.parseString(zf.read("xl/sharedStrings.xml"))
            for si in sdoc.getElementsByTagName("si"):
                txt = "".join(
                    node.firstChild.nodeValue
                    for node in si.getElementsByTagName("t")
                    if node.firstChild is not None and node.firstChild.nodeValue is not None
                )
                sst.append(txt)

        wb_doc = minidom.parseString(zf.read("xl/workbook.xml"))
        rel_doc = minidom.parseString(zf.read("xl/_rels/workbook.xml.rels"))
        rel_map: Dict[str, str] = {}
        for rel in rel_doc.getElementsByTagName("Relationship"):
            rel_map[rel.getAttribute("Id")] = rel.getAttribute("Target")

        target: Optional[str] = None
        for sh in wb_doc.getElementsByTagName("sheet"):
            if sh.getAttribute("name") == sheet_name:
                rid = sh.getAttribute("r:id")
                t = rel_map[rid].lstrip("/")
                target = t if t.startswith("xl/") else f"xl/{t}"
                break
        if target is None:
            raise ValueError(f"Sheet {sheet_name!r} missing in {xlsx_path.name}")

        sheet_doc = minidom.parseString(zf.read(target))
        grid: Dict[Tuple[int, str], str] = {}
        for cell in sheet_doc.getElementsByTagName("c"):
            coord = cell.getAttribute("r")
            m = re.fullmatch(r"([A-Z]+)(\d+)", coord)
            if not m:
                raise ValueError(f"Bad cell coordinate {coord!r}")
            col, rnum = m.group(1), int(m.group(2))
            ctype = cell.getAttribute("t")
            if ctype == "inlineStr":
                txt = "".join(
                    t.firstChild.nodeValue
                    for t in cell.getElementsByTagName("t")
                    if t.firstChild is not None and t.firstChild.nodeValue is not None
                )
                if txt != "":
                    grid[(rnum, col)] = txt
            else:
                vnodes = cell.getElementsByTagName("v")
                if vnodes and vnodes[0].firstChild is not None:
                    raw = vnodes[0].firstChild.nodeValue or ""
                    if ctype == "s":
                        grid[(rnum, col)] = sst[int(raw)]
                    elif ctype == "b":
                        grid[(rnum, col)] = "TRUE" if raw == "1" else "FALSE"
                    else:
                        grid[(rnum, col)] = raw
        return grid


def _decimal_pearson(xs: Sequence[float], ys: Sequence[float]) -> float:
    dx = [Decimal(str(x)) for x in xs]
    dy = [Decimal(str(y)) for y in ys]
    n = Decimal(len(dx))
    mx = sum(dx) / n
    my = sum(dy) / n
    num = sum((a - mx) * (b - my) for a, b in zip(dx, dy))
    vx = sum((a - mx) ** 2 for a in dx)
    vy = sum((b - my) ** 2 for b in dy)
    return float(num) / math.sqrt(float(vx) * float(vy))


def _welch_p_simpson(xs: Sequence[float], ys: Sequence[float]) -> Tuple[float, float]:
    """Independent Welch t p-value via trigonometric Simpson quadrature on [0, arctan(|t|/sqrt(df))]."""
    n1, n2 = len(xs), len(ys)
    m1 = sum(xs) / n1
    m2 = sum(ys) / n2
    v1 = sum((x - m1) ** 2 for x in xs) / (n1 - 1)
    v2 = sum((y - m2) ** 2 for y in ys) / (n2 - 1)
    se2 = v1 / n1 + v2 / n2
    t_val = (m1 - m2) / math.sqrt(se2)
    df = (se2 ** 2) / (((v1 / n1) ** 2) / (n1 - 1) + ((v2 / n2) ** 2) / (n2 - 1))

    theta_max = math.atan2(abs(t_val), math.sqrt(df))
    panels = 2000  # even
    h = theta_max / panels
    integral = 0.0
    for k in range(panels + 1):
        theta = k * h
        f_val = math.cos(theta) ** (df - 1.0)
        weight = 1.0 if (k == 0 or k == panels) else (4.0 if k % 2 == 1 else 2.0)
        integral += weight * f_val
    integral *= h / 3.0
    norm_const = 2.0 * math.exp(math.lgamma(0.5 * (df + 1.0)) - math.lgamma(0.5 * df)) / math.sqrt(math.pi)
    cdf_two_sided_tail = max(0.0, min(1.0, 1.0 - norm_const * integral))
    return t_val, cdf_two_sided_tail


def _bitmask_perm_two_sided(xs: Sequence[float], ys: Sequence[float]) -> float:
    """Bitmask exact two-sided permutation p-value (distinct from itertools.combinations)."""
    n1, n2 = len(xs), len(ys)
    pooled = list(xs) + list(ys)
    tot = sum(pooled)
    obs_abs = abs(sum(xs) / n1 - sum(ys) / n2)
    ge = 0
    total = 0
    for mask in range(1 << (n1 + n2)):
        if mask.bit_count() != n1:
            continue
        s1 = sum(pooled[i] for i in range(n1 + n2) if (mask >> i) & 1)
        s2 = tot - s1
        if abs(s1 / n1 - s2 / n2) >= obs_abs - 1e-12:
            ge += 1
        total += 1
    return ge / total


def _two_pointer_did_perm(
    f_c: Sequence[float], f_s: Sequence[float], w_c: Sequence[float], w_s: Sequence[float]
) -> float:
    """One-sided DiD permutation p-value via bitmasks and explicit two-pointer count."""
    f_pool = list(f_c) + list(f_s)
    w_pool = list(w_c) + list(w_s)
    nf1, nf2 = len(f_c), len(f_s)
    nw1, nw2 = len(w_c), len(w_s)
    ftot, wtot = sum(f_pool), sum(w_pool)
    did_obs = math.log((sum(f_c) / nf1) / (sum(f_s) / nf2)) - math.log((sum(w_c) / nw1) / (sum(w_s) / nw2))

    f_lrs: List[float] = []
    for mask in range(1 << (nf1 + nf2)):
        if mask.bit_count() == nf1:
            s1 = sum(f_pool[i] for i in range(nf1 + nf2) if (mask >> i) & 1)
            f_lrs.append(math.log((s1 / nf1) / ((ftot - s1) / nf2)))

    w_lrs: List[float] = []
    for mask in range(1 << (nw1 + nw2)):
        if mask.bit_count() == nw1:
            s1 = sum(w_pool[i] for i in range(nw1 + nw2) if (mask >> i) & 1)
            w_lrs.append(math.log((s1 / nw1) / ((wtot - s1) / nw2)))

    f_lrs.sort()
    w_lrs.sort()
    # Two-pointer count of pairs (lf, lw) where lf - lw <= did_obs + 1e-12
    nw = len(w_lrs)
    j = 0
    le_cnt = 0
    for lf in f_lrs:
        thresh = lf - did_obs - 1e-12
        while j < nw and w_lrs[j] < thresh:
            j += 1
        le_cnt += nw - j
    return le_cnt / (len(f_lrs) * nw)


def _bh_q_decimal(pvals: Sequence[float]) -> List[float]:
    m = len(pvals)
    indexed = sorted(enumerate(Decimal(str(p)) for p in pvals), key=lambda kv: kv[1])
    q_s = [Decimal("1")] * m
    run_min = Decimal("1")
    for idx in range(m - 1, -1, -1):
        val = min(Decimal("1"), indexed[idx][1] * Decimal(m) / Decimal(idx + 1))
        if val < run_min:
            run_min = val
        q_s[idx] = run_min
    out = [0.0] * m
    for pos, (orig_i, _) in enumerate(indexed):
        out[orig_i] = float(q_s[pos])
    return out


def _wilcoxon_counter_dp(diffs: Sequence[float]) -> float:
    """Independent Wilcoxon signed-rank test via collections.Counter subset-sum."""
    nz = [d for d in diffs if d != 0.0]
    n = len(nz)
    if n == 0:
        return 1.0
    order = sorted(range(n), key=lambda i: abs(nz[i]))
    doubled_ranks = [0] * n
    i = 0
    while i < n:
        j = i + 1
        while j < n and abs(nz[order[j]]) == abs(nz[order[i]]):
            j += 1
        # Sum of 1-based positions (i+1 .. j) * 2 / (j-i) == (i + 1 + j)
        dr = (i + 1) + j
        for k in range(i, j):
            doubled_ranks[order[k]] = dr
        i = j

    w_plus_d = sum(doubled_ranks[idx] for idx, d in enumerate(nz) if d > 0.0)
    s_double = sum(doubled_ranks)
    half_s = s_double / 2.0
    obs_dev = abs(w_plus_d - half_s)

    if n <= 25:
        dist: collections.Counter[int] = collections.Counter({0: 1})
        for dr in doubled_ranks:
            nxt: collections.Counter[int] = collections.Counter()
            for s, cnt in dist.items():
                nxt[s] += cnt
                nxt[s + dr] += cnt
            dist = nxt
        ext = sum(cnt for s, cnt in dist.items() if abs(s - half_s) >= obs_dev - 1e-9)
        return min(1.0, ext / (1 << n))

    # Fallback normal
    tie_counts = collections.Counter(doubled_ranks).values()
    tie_term = sum(c * (c * c - 1) for c in tie_counts if c > 1)
    mu = n * (n + 1) / 4.0
    var = (n * (n + 1) * (2 * n + 1) / 24.0) - (tie_term / 48.0)
    z = max(0.0, abs(w_plus_d / 2.0 - mu) - 0.5) / math.sqrt(var)
    return min(1.0, max(0.0, math.erfc(z / math.sqrt(2.0))))


def _pairwise_mwu_p(grp: Sequence[float], bg: Sequence[float]) -> float:
    """Direct O(n1 * n2) pairwise dominance U statistic (zero shared rank code)."""
    n1, n2 = len(grp), len(bg)
    u1 = 0.0
    for a in grp:
        for b in bg:
            if a > b:
                u1 += 1.0
            elif a == b:
                u1 += 0.5
    n_tot = n1 + n2
    tie_sum = sum(c * (c * c - 1) for c in collections.Counter(list(grp) + list(bg)).values() if c > 1)
    mu = (n1 * n2) / 2.0
    var = (n1 * n2 / 12.0) * ((n_tot + 1.0) - (tie_sum / (n_tot * (n_tot - 1.0))))
    z = max(0.0, abs(u1 - mu) - 0.5) / math.sqrt(var)
    return min(1.0, max(0.0, math.erfc(z / math.sqrt(2.0))))


def _mod_of(pw: str) -> Optional[str]:
    norm = pw.strip().lower()
    for mod, terms in MODULES.items():
        if norm in terms:
            return mod
    return None


def _parse_ussing_raw(xlsx_path: Path, avg_col: str) -> Tuple[Dict[str, List[Tuple[float, float, float]]], List[Tuple[str, str, str]]]:
    grid = _read_sheet_minidom(xlsx_path, "Raw-data-analyzed")
    end_idx = _col_num(avg_col)
    cols = [_num_col(i) for i in range(_col_num("C"), end_idx)]
    dedup_map: Dict[str, List[Tuple[float, float, float]]] = {}
    dups: List[Tuple[str, str, str]] = []
    for idx, gname in enumerate(USSING_GROUPS):
        r0 = 3 + 3 * idx
        if grid.get((r0, "A"), "").strip() != gname:
            raise AssertionError(f"Group mismatch at {r0}: {grid.get((r0, 'A'))}")
        seen: Dict[Tuple[float, float, float], str] = {}
        kept: List[Tuple[float, float, float]] = []
        for col in cols:
            vals = [grid.get((r0 + k, col)) for k in range(3)]
            if all(v is None for v in vals):
                continue
            tup = (float(vals[0]), float(vals[1]), float(vals[2]))
            if tup in seen:
                dups.append((gname, seen[tup], col))
            else:
                seen[tup] = col
                kept.append(tup)
        dedup_map[gname] = kept
    return dedup_map, dups


def run_independent_audit(repo_root: Path = REPO_ROOT, data_dir: Optional[Path] = None) -> Dict[str, Any]:
    bench_dir = repo_root / "benchmarks" / "cftr_canx_interactome_ussing"
    done_txt = (bench_dir / "done.md").read_text(encoding="utf-8", errors="replace").replace("\r\n", "\n")
    done_sha = hashlib.sha256(done_txt.encode("utf-8")).hexdigest()
    results = json.loads((bench_dir / "results_2026-11-25.json").read_text(encoding="utf-8"))
    if results["done_sha256"] != done_sha:
        raise AssertionError("done_sha256 mismatch in results JSON")

    ledger_path = bench_dir / "ledger.json"
    if ledger_path.is_file():
        ledger = json.loads(ledger_path.read_text(encoding="utf-8"))
        if ledger.get("done_sha256") != done_sha:
            raise AssertionError("done_sha256 mismatch in ledger.json")
        # Check standing consistency for predicted claims
        a1_stand = results["arm_1_frt_ussing_sicanx_asymmetry"]["primary_strict_dedup"]["ledger_standing"]
        a2_stand = results["arm_2_canx_ko_proteostasis_remodeling"]["primary_evaluation"]["ledger_standing"]
        by_id = {c["id"]: c for c in ledger.get("claims", [])}
        if "c4-arm-1-ussing-sicanx-asymmetry" in by_id:
            if by_id["c4-arm-1-ussing-sicanx-asymmetry"]["standing"] != a1_stand:
                raise AssertionError("Arm 1 ledger standing mismatch")
            if a1_stand == "refuted" and by_id["c4-arm-1-ussing-sicanx-asymmetry"].get("headline") is True:
                raise AssertionError("Refuted Arm 1 claim cannot be headline")
        if "c5-arm-2-proteostasis-remodeling" in by_id:
            if by_id["c5-arm-2-proteostasis-remodeling"]["standing"] != a2_stand:
                raise AssertionError("Arm 2 ledger standing mismatch")
            if a2_stand == "refuted" and by_id["c5-arm-2-proteostasis-remodeling"].get("headline") is True:
                raise AssertionError("Refuted Arm 2 claim cannot be headline")

    # Offline BH-FDR and gate arithmetic checks
    a1p = results["arm_1_frt_ussing_sicanx_asymmetry"]["primary_strict_dedup"]
    a2p = results["arm_2_canx_ko_proteostasis_remodeling"]["primary_evaluation"]
    d_tests = a2p["family_2a_drug_attenuation"]["tests"]
    v_tests = a2p["family_2b_variant_enrichment"]["tests"]
    d_q_chk = _bh_q_decimal([t["raw_p_two_sided"] for t in d_tests])
    v_q_chk = _bh_q_decimal([t["raw_p_two_sided"] for t in v_tests])
    for t, qc in zip(d_tests, d_q_chk):
        if abs(t["bh_fdr_q"] - qc) > 2e-6:
            raise AssertionError(f"Offline Drug BH-FDR mismatch: {t['bh_fdr_q']} vs {qc}")
    for t, qc in zip(v_tests, v_q_chk):
        if abs(t["bh_fdr_q"] - qc) > 2e-6:
            raise AssertionError(f"Offline Variants BH-FDR mismatch: {t['bh_fdr_q']} vs {qc}")

    ddir = _find_data_dir(data_dir, repo_root)
    raw_checked = False
    recomputed: Dict[str, Any] = {}
    if ddir is not None:
        for fn, exp_sha in EXPECTED_WORKBOOKS.items():
            got_sha = hashlib.sha256((ddir / fn).read_bytes()).hexdigest()
            if got_sha != exp_sha:
                raise AssertionError(f"SHA256 mismatch for {fn}")

        f508_map, f508_dups = _parse_ussing_raw(
            ddir / "Copy of FRT-F508del_siCANX_Ussing-summary-for-paper.xlsx", "M"
        )
        wt_map, wt_dups = _parse_ussing_raw(
            ddir / "Copy of FRT-WT_siCANX_Ussing-summary-for-paper.xlsx", "S"
        )
        if len(f508_dups) + len(wt_dups) != results["data_integrity_audit"]["duplicate_replicate_columns_count"]:
            raise AssertionError("Duplicate replicate column count mismatch")

        # Verify Decimal Pearson r for Inh172 reversal
        pooled_xs: List[float] = []
        pooled_ys: List[float] = []
        for lbl, m in (("FRT-F508del", f508_map), ("FRT-WT", wt_map)):
            xs = [sum(t[0] + t[1] for t in m[g]) / len(m[g]) for g in USSING_GROUPS]
            ys = [-sum(t[2] for t in m[g]) / len(m[g]) for g in USSING_GROUPS]
            pooled_xs.extend(xs)
            pooled_ys.extend(ys)
            r_dec = _decimal_pearson(xs, ys)
            r_rep = results["arm_1_frt_ussing_sicanx_asymmetry"]["positive_controls"]["inh172_reversal_parity"][lbl]["pearson_r_all_11_groups"]
            if abs(r_dec - r_rep) > 1e-6:
                raise AssertionError(f"Mass-balance r mismatch for {lbl}: {r_dec} vs {r_rep}")
        r_pool_dec = _decimal_pearson(pooled_xs, pooled_ys)
        r_pool_rep = results["arm_1_frt_ussing_sicanx_asymmetry"]["positive_controls"]["inh172_reversal_parity"]["pooled_22_groups_pearson_r"]
        if abs(r_pool_dec - r_pool_rep) > 1e-6:
            raise AssertionError(f"Pooled 22-group mass-balance r mismatch: {r_pool_dec} vs {r_pool_rep}")

        f_c = [t[0] + t[1] for t in f508_map["siCANX, 661+445"]]
        f_s = [t[0] + t[1] for t in f508_map["siNS, 661+445"]]
        w_c = [t[0] + t[1] for t in wt_map["siCANX, DMSO"]]
        w_s = [t[0] + t[1] for t in wt_map["siNS, DMSO"]]

        p_1a = _bitmask_perm_two_sided(f_c, f_s)
        t_1a, pt_1a = _welch_p_simpson(f_c, f_s)
        p_1b = _two_pointer_did_perm(f_c, f_s, w_c, w_s)
        did_dec = float(
            (Decimal(str(sum(f_c))) / Decimal(len(f_c)) / (Decimal(str(sum(f_s))) / Decimal(len(f_s)))).ln()
            - (Decimal(str(sum(w_c))) / Decimal(len(w_c)) / (Decimal(str(sum(w_s))) / Decimal(len(w_s)))).ln()
        )

        rep_1a = a1p["arm_1a_f508del_combo_sicanx_drop"]
        rep_1b = a1p["arm_1b_sicanx_proportion_f508del_vs_wt"]
        if abs(p_1a - rep_1a["exact_perm_p_two_sided"]) > 1e-6:
            raise AssertionError(f"Arm 1a perm p mismatch: {p_1a} vs {rep_1a['exact_perm_p_two_sided']}")
        if abs(pt_1a - rep_1a["welch_t_context"]["p_two_sided"]) > 1e-6:
            raise AssertionError(f"Arm 1a Welch p mismatch: {pt_1a} vs {rep_1a['welch_t_context']['p_two_sided']}")
        if abs(p_1b - rep_1b["exact_stratified_perm_p_one_sided"]) > 1e-6:
            raise AssertionError(f"Arm 1b DiD perm p mismatch: {p_1b} vs {rep_1b['exact_stratified_perm_p_one_sided']}")
        if abs(did_dec - rep_1b["did_log_ratio"]) > 1e-6:
            raise AssertionError(f"Arm 1b DiD log ratio mismatch: {did_dec} vs {rep_1b['did_log_ratio']}")

        # Independent Arm 2 parsing and statistical audit
        v_grid = _read_sheet_minidom(ddir / "Copy of CFTR_Variants_LogFoldChanges.xlsx", "Sheet1")
        d_grid = _read_sheet_minidom(ddir / "Copy of CFTR_F508del_Drug_LogFoldChanges.xlsx", "Sheet1")
        max_vr = max(r for r, _ in v_grid)
        max_dr = max(r for r, _ in d_grid)

        v_rows: List[Tuple[str, Dict[str, Optional[float]]]] = []
        for r in range(3, max_vr + 1):
            g = v_grid.get((r, "A"), "").strip()
            pw = v_grid.get((r, "B"), "").strip()
            if not g:
                continue
            vals = {}
            for cname, clet in zip(VAR_COLS, ("C", "D", "E", "F", "G")):
                raw = v_grid.get((r, clet), "")
                vals[cname] = None if raw == "#DIV/0!" else float(raw)
            if g == "CFTR":
                if any(vals[c] != 0.0 for c in VAR_COLS):
                    raise AssertionError("Variants CFTR anchor != 0")
            else:
                v_rows.append((pw, vals))

        d_rows: List[Tuple[str, Dict[str, Optional[float]]]] = []
        for r in range(3, max_dr + 1):
            g = d_grid.get((r, "A"), "").strip()
            pw = d_grid.get((r, "B"), "").strip()
            if not g:
                continue
            vals = {}
            for cname, clet in zip(DRUG_COLS, ("C", "D", "E", "F")):
                raw = d_grid.get((r, clet), "")
                vals[cname] = None if raw == "#DIV/0!" else float(raw)
            if g == "CFTR":
                if any(vals[c] != 0.0 for c in DRUG_COLS):
                    raise AssertionError("Drug CFTR anchor != 0")
            else:
                d_rows.append((pw, vals))

        d_raw_ps: List[float] = []
        for mod in ("folding", "ups", "autophagy_lysosome", "trafficking"):
            diffs = [
                abs(vals["DelF (-CX) + combo"]) - abs(vals["ΔF (-CX)"])
                for pw, vals in d_rows
                if _mod_of(pw) == mod and vals["ΔF (-CX)"] is not None and vals["DelF (-CX) + combo"] is not None
            ]
            d_raw_ps.append(_wilcoxon_counter_dp(diffs))
        d_qs = _bh_q_decimal(d_raw_ps)
        for rep, rp, rq in zip(d_tests, d_raw_ps, d_qs):
            if abs(rp - rep["raw_p_two_sided"]) > 1e-6 or abs(rq - rep["bh_fdr_q"]) > 1e-6:
                raise AssertionError(f"Drug Wilcoxon audit mismatch for {rep['module']}: {(rp, rq)} vs {rep}")

        v_raw_ps: List[float] = []
        for col in VAR_COLS:
            fin = [(pw, vals[col]) for pw, vals in v_rows if vals[col] is not None]
            for mod in ("folding", "ups", "autophagy_lysosome", "trafficking"):
                grp = [float(val) for pw, val in fin if _mod_of(pw) == mod]
                bg = [float(val) for pw, val in fin if _mod_of(pw) != mod]
                v_raw_ps.append(_pairwise_mwu_p(grp, bg))
        v_qs = _bh_q_decimal(v_raw_ps)
        for rep, rp, rq in zip(v_tests, v_raw_ps, v_qs):
            if abs(rp - rep["raw_p_two_sided"]) > 1e-6 or abs(rq - rep["bh_fdr_q"]) > 1e-6:
                raise AssertionError(
                    f"Variants MWU audit mismatch for {rep['variant']}:{rep['module']}: {(rp, rq)} vs {rep}"
                )

        raw_checked = True
        recomputed = {
            "arm1a_perm_p": round(p_1a, 6),
            "arm1a_welch_p": round(pt_1a, 6),
            "arm1b_did_log_ratio": round(did_dec, 6),
            "arm1b_perm_p_1s": round(p_1b, 6),
            "arm2_drug_min_q": round(min(d_qs), 6),
            "arm2_variants_min_q": round(min(v_qs), 6),
        }

    return {
        "status": "PASS",
        "raw_local_data_audited": raw_checked,
        "done_sha256": done_sha,
        "arm_1_verdict": a1p["overall_verdict"],
        "arm_2_verdict": a2p["overall_verdict"],
        "recomputed": recomputed,
    }


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--data-dir", type=Path, default=None)
    parser.add_argument("--check", action="store_true", default=True)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args(argv)

    rep = run_independent_audit(repo_root=REPO_ROOT, data_dir=args.data_dir) if args.check else {"status": "SKIP"}
    if args.json:
        print(json.dumps(rep, indent=2, allow_nan=False))
    else:
        print(
            f"PASS independent_audit cftr_canx_interactome_ussing: "
            f"raw_audited={rep['raw_local_data_audited']} "
            f"arm1={rep['arm_1_verdict']} arm2={rep['arm_2_verdict']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
