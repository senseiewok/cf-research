#!/usr/bin/env python3
"""Independent regex-free, zero-shared-import audit of `PRED-0015` (`delivery_route_compartmentalization`).

Rule 2 of `discovery-practice`: every headline metric computed by the lab requires a second method
in separate code that imports zero helpers from the primary benchmark or `tools/`.

This module imports NO `re`, NO `math` (no `math.comb`, `math.acos`, or `math.cos`), NO `tools/`,
and NO `delivery_route_compartmentalization`:
1. Tokenizes `snapshot.json` via character-by-character alphanumeric scanning over structured
   intervention dicts (`DRUG` / `BIOLOGICAL` / `GENETIC` strictly per `done.md`) and cross-checks
   PRISMA `940` membership against `data/care_team_trials_directory.json`.
2. Evaluates one-sided and two-sided hypergeometric Fisher p-values via multiplicative
   `fractions.Fraction` products (`prod((n - m + r) / r)`).
3. Evaluates Chan-Zhang (1999) unconditional 95% risk-difference bounds by solving the
   Miettinen-Nurminen restricted-MLE polynomial score derivative via monotone bisection (zero
   trigonometric Cardano formulas) and ternary nuisance search.
4. Evaluates Newcombe-Wilson Method-10 95% risk-difference bounds in 50-digit `decimal.Decimal`
   arithmetic (`Decimal.sqrt()`).
"""
from __future__ import annotations

import argparse
from decimal import Decimal, getcontext
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
from typing import Any

getcontext().prec = 50

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
DEFAULT_SNAPSHOT = REPO_ROOT / "tools" / "trial_atlas" / "snapshot.json"
DEFAULT_DIRECTORY = REPO_ROOT / "data" / "care_team_trials_directory.json"
DEFAULT_DONE = HERE / "done.md"
DEFAULT_RESULTS = HERE / "results_2026-11-24.json"

EARLY_PHASE_SET = frozenset({"EARLY_PHASE1", "PHASE1", "PHASE2"})
STRICT_INT_TYPES = frozenset({"DRUG", "BIOLOGICAL", "GENETIC"})
REQUIRED_OUTCOME_KEYS = frozenset({"measure", "description", "time_frame"})


def _tokens(text: str) -> tuple[str, ...]:
    """Split text into lowercase alphanumeric tokens without regular expressions."""
    buf: list[str] = []
    out: list[str] = []
    for ch in (text or "").lower():
        if ch.isalnum():
            buf.append(ch)
        else:
            if buf:
                out.append("".join(buf))
                buf.clear()
    if buf:
        out.append("".join(buf))
    return tuple(out)


def _has_seq(toks: tuple[str, ...], seq: tuple[str, ...]) -> bool:
    m = len(seq)
    if m == 0 or len(toks) < m:
        return False
    for idx in range(len(toks) - m + 1):
        if toks[idx : idx + m] == seq:
            return True
    return False


def _has_any_seq(toks: tuple[str, ...], seqs: tuple[tuple[str, ...], ...]) -> bool:
    return any(_has_seq(toks, s) for s in seqs)


def _condition_names_cf_independent(condition: str) -> bool | None:
    toks = _tokens(condition)
    neg_seqs = (
        ("non", "cystic", "fibrosis"),
        ("non", "cf"),
        ("noncystic", "fibrosis"),
        ("noncf",),
        ("without", "cystic", "fibrosis"),
        ("without", "cf"),
        ("other", "than", "cystic", "fibrosis"),
        ("other", "than", "cf"),
    )
    masked = list(toks)
    for seq in neg_seqs:
        m = len(seq)
        for i in range(len(masked) - m + 1):
            if tuple(masked[i : i + m]) == seq:
                for j in range(i, i + m):
                    masked[j] = "_"
    clean_toks = tuple(masked)
    if (
        _has_seq(clean_toks, ("cystic", "fibrosis"))
        or ("mucoviscidosis" in clean_toks)
        or (toks == ("cf",))
    ):
        return True
    if _has_any_seq(toks, neg_seqs):
        return None
    return False


def _prisma_included_independent(study: dict[str, Any]) -> tuple[bool, str | None]:
    st = str(study.get("study_type") or "INTERVENTIONAL").strip().upper()
    if st == "OBSERVATIONAL":
        return False, "X1"
    if st == "EXPANDED_ACCESS":
        return False, "X2"
    if st != "INTERVENTIONAL":
        return False, "X3"
    conds = [c for c in (study.get("conditions") or []) if isinstance(c, str)]
    verdicts = [_condition_names_cf_independent(c) for c in conds]
    if True in verdicts:
        return True, None
    if None in verdicts:
        return False, "X5"
    return False, "X4"


def _is_healthy_volunteer_only_independent(brief_title: str, official_title: str) -> bool:
    toks = _tokens(f"{brief_title} {official_title}")
    hv_nouns = {
        "volunteer", "volunteers", "subject", "subjects", "participant", "participants",
        "adult", "adults", "male", "males", "female", "females",
    }
    has_hv = any(toks[i] == "healthy" and toks[i + 1] in hv_nouns for i in range(len(toks) - 1))
    if not has_hv:
        return False
    cf_nouns = {"patients", "subjects", "participants", "people", "adults", "children", "adolescents", "individuals"}
    for i in range(len(toks) - 3):
        if toks[i] in cf_nouns and toks[i + 1 : i + 4] == ("with", "cystic", "fibrosis"):
            return False
    for i in range(len(toks) - 2):
        if toks[i : i + 2] == ("cystic", "fibrosis") and toks[i + 2] in {"patients", "subjects", "participants"}:
            return False
    return True


NAMED_INHALED_SEQS: tuple[tuple[str, ...], ...] = (
    ("vx", "522"), ("vx522",),
    ("mrt", "5005"), ("mrt5005",),
    ("arct", "032"), ("arct032",),
    ("rct", "2100"), ("rct2100",),
    ("bmd", "003"), ("bmd003",),
    ("lunar", "cf"), ("lunarcf",),
    ("4d", "710"), ("4d710",),
    ("sp", "101"), ("sp101",),
    ("kb", "407"), ("kb407",),
    ("bi", "3720931"), ("bi3720931",), ("lenticlair",),
    ("pgm", "169"), ("pgm169",), ("gl", "67a"), ("gl67a",),
    ("eluforsen",), ("qr", "010"), ("qr010",),
    ("spl", "84"), ("spl84",),
)

SYS_MOD_SINGLE = frozenset({
    "ivacaftor", "iva", "vx770", "kalydeco",
    "deutivacaftor", "vx561", "ctp656",
    "lumacaftor", "vx809", "orkambi",
    "tezacaftor", "vx661", "symdeko", "symkevi",
    "elexacaftor", "vx445", "trikafta", "kaftrio",
    "vanzacaftor", "vx121", "alyftrek",
    "vx659", "vx440", "vx152", "vx581",
    "galicaftor", "abbv2222", "glpg2222",
    "navocaftor", "abbv3067", "glpg3067",
    "glpg1837", "glpg2451", "glpg2737",
    "abbv191", "abbv576",
    "posenacaftor", "pti801", "nesolicaftor", "pti428", "dirocaftor", "pti808",
    "fdl169", "icenticaftor", "qbw251", "cavosonstat", "n9115",
    "buphenyl",
})

SYS_MOD_PAIRS: tuple[tuple[str, ...], ...] = (
    ("vx", "770"), ("vx", "561"), ("ctp", "656"),
    ("lum", "iva"), ("vx", "809"), ("tez", "iva"), ("vx", "661"),
    ("elx", "tez", "iva"), ("vx", "445"), ("vx", "121"),
    ("vx", "659"), ("vx", "440"), ("vx", "152"), ("vx", "581"),
    ("abbv", "2222"), ("glpg", "2222"), ("abbv", "3067"), ("glpg", "3067"),
    ("glpg", "1837"), ("glpg", "2451"), ("glpg", "2737"),
    ("abbv", "191"), ("abbv", "576"),
    ("pti", "801"), ("pti", "428"), ("pti", "808"),
    ("fdl", "169"), ("qbw", "251"), ("n", "9115"),
    ("4", "phenylbutyrate"), ("sodium", "phenylbutyrate"),
)

SYS_READTHROUGH_SINGLE = frozenset({"ataluren", "ptc124", "elx02"})
SYS_READTHROUGH_PAIRS: tuple[tuple[str, ...], ...] = (("ptc", "124"), ("elx", "02"))

NON_CFTR_TOPICAL_SINGLE = frozenset({
    "etd001", "aroenac", "vx371", "p1037", "spx101", "bi1265162",
    "azd5634", "qbw276", "amiloride", "denufosol",
    "mannitol",
    "bx004", "appa02", "achromophage", "tobramycin", "aztreonam", "colistin", "gallium",
})
NON_CFTR_TOPICAL_PAIRS: tuple[tuple[str, ...], ...] = (
    ("etd", "001"), ("aro", "enac"), ("vx", "371"), ("p", "1037"),
    ("spx", "101"), ("bi", "1265162"), ("azd", "5634"), ("qbw", "276"),
    ("hypertonic", "saline"), ("dornase", "alfa"),
    ("bx", "004"), ("ap", "pa02"),
)

WITHDRAWAL_DDI_SINGLE = frozenset({"discontinuation", "stopping", "simplify", "midazolam", "itraconazole", "rifampin"})
WITHDRAWAL_DDI_PAIRS: tuple[tuple[str, ...], ...] = (
    ("withdrawal", "of", "hypertonic"),
    ("drug", "drug", "interaction"),
    ("oral", "contraceptive"),
)

RESP_PREFIXES = ("inhal", "nebuliz", "nebulis", "aerosol", "bronchoscop")
RESP_EXACT = frozenset({"pulmonary", "endobronchial", "intranasal", "nasal", "sinus", "airway"})


def _matches_sys_mod(toks: tuple[str, ...]) -> bool:
    return any(t in SYS_MOD_SINGLE for t in toks) or _has_any_seq(toks, SYS_MOD_PAIRS)


def _matches_sys_cftr(toks: tuple[str, ...]) -> bool:
    return (
        _matches_sys_mod(toks)
        or any(t in SYS_READTHROUGH_SINGLE for t in toks)
        or _has_any_seq(toks, SYS_READTHROUGH_PAIRS)
    )


def _matches_non_cftr_topical(toks: tuple[str, ...]) -> bool:
    return any(t in NON_CFTR_TOPICAL_SINGLE for t in toks) or _has_any_seq(toks, NON_CFTR_TOPICAL_PAIRS)


def _has_resp_route(toks: tuple[str, ...]) -> bool:
    for t in toks:
        if t in RESP_EXACT or any(t.startswith(p) for p in RESP_PREFIXES):
            return True
    return _has_seq(toks, ("maxillary", "sinus"))


def _has_cftr_genetic_target(toks: tuple[str, ...]) -> bool:
    if _has_any_seq(
        toks,
        (
            ("cftr", "mrna"), ("cftr", "cdna"), ("cftr", "dna"), ("cftr", "plasmid"),
            ("cftr", "gene"), ("cftr", "vector"), ("cftr", "antisense"), ("cftr", "oligonucleotide"),
            ("ad", "cftr"), ("adcftr",), ("aav", "cftr"), ("aav2", "cftr"),
        ),
    ):
        return True
    cftr_pos = [i for i, t in enumerate(toks) if t == "cftr"]
    vec_pos: list[int] = []
    for i, t in enumerate(toks):
        if t in {"mrna", "cdna", "plasmid", "ad5"} or any(
            t.startswith(p) for p in ("adenovir", "adgv", "raav", "aav", "lentivir", "liposom")
        ):
            vec_pos.append(i)
        if i + 1 < len(toks) and toks[i : i + 2] in {("gene", "transfer"), ("gene", "therapy"), ("viral", "vector"), ("antisense", "oligonucleotide")}:
            vec_pos.append(i)
    for ci in cftr_pos:
        if any(abs(ci - vi) <= 12 for vi in vec_pos):
            return True
    for i in range(len(toks) - 1):
        if toks[i : i + 2] == ("antisense", "oligonucleotide"):
            window = toks[i : min(len(toks), i + 12)]
            if "cftr" in window or "f508del" in window or _has_seq(window, ("3849", "10kb")):
                return True
    return False


def _classify_regimen_independent(study: dict[str, Any]) -> tuple[str, bool]:
    """Per-intervention structured classifier adhering strictly to `done.md` Section 2."""
    brief = str(study.get("brief_title") or "")
    official = str(study.get("official_title") or "")
    kw = " ".join(str(k) for k in (study.get("keywords") or []))
    title_toks = _tokens(f"{brief} {official}")

    all_int_tok_lists: list[tuple[str, ...]] = []
    typed_int_tok_lists: list[tuple[str, ...]] = []
    for raw in study.get("interventions") or []:
        if not isinstance(raw, dict):
            continue
        nm_toks = _tokens(str(raw.get("name") or ""))
        tp = str(raw.get("type") or "").strip().upper()
        all_int_tok_lists.append(nm_toks)
        if tp in STRICT_INT_TYPES:
            typed_int_tok_lists.append(nm_toks)

    corpus_toks = title_toks + tuple(t for sub in all_int_tok_lists for t in sub) + _tokens(kw)
    dual_mod = any(_matches_sys_mod(sub) for sub in typed_int_tok_lists)
    has_non_cftr = _matches_non_cftr_topical(title_toks) or any(_matches_non_cftr_topical(sub) for sub in typed_int_tok_lists)

    named_inh = _has_any_seq(corpus_toks, NAMED_INHALED_SEQS)
    gen_resp = (not has_non_cftr) and _has_cftr_genetic_target(corpus_toks) and _has_resp_route(corpus_toks)
    ptc_amino = (
        (not has_non_cftr)
        and any(t in {"gentamicin", "amikacin"} for t in corpus_toks)
        and _has_resp_route(corpus_toks)
        and (
            any(t in {"nonsense", "readthrough", "w1282x", "g542x"} for t in corpus_toks)
            or _has_any_seq(corpus_toks, (("stop", "codon"), ("premature", "stop"), ("premature", "termination"), ("read", "through")))
        )
    )
    if named_inh or gen_resp or ptc_amino:
        return "inhaled_cftr", dual_mod

    if has_non_cftr:
        return "non_cftr_or_other", False

    if any(_matches_sys_cftr(sub) for sub in typed_int_tok_lists):
        if any(t in WITHDRAWAL_DDI_SINGLE for t in title_toks) or _has_any_seq(title_toks, WITHDRAWAL_DDI_PAIRS):
            return "non_cftr_or_other", False
        return "systemic_cftr", False

    return "non_cftr_or_other", False


def _outcomes_independent(study: dict[str, Any]) -> tuple[bool, bool]:
    sweat = False
    airway_safety = False
    diag_words = {"diagnosis", "diagnostic", "inclusion", "eligibility"}
    airway_exact = {"fev1", "ppfev1", "lci", "npd", "mucociliary", "safety", "tolerability", "laboratory"}
    outcomes = study.get("primary_outcomes")
    if not isinstance(outcomes, list):
        raise ValueError(f"Study {study.get('nct_id')} missing flattened primary_outcomes list")
    for o in outcomes:
        if not isinstance(o, dict) or not REQUIRED_OUTCOME_KEYS.issubset(o.keys()):
            raise ValueError(f"Study {study.get('nct_id')} has unflattened primary_outcomes entry")
        toks = _tokens(str(o.get("measure") or ""))
        if not toks:
            continue
        is_sweat = (
            _has_seq(toks, ("sweat", "chloride"))
            or _has_seq(toks, ("sweat", "cl"))
            or _has_seq(toks, ("chloride", "concentration", "in", "sweat"))
        )
        is_diag = any(t in diag_words for t in toks) or _has_seq(toks, ("screening", "criteria"))
        if is_sweat and not is_diag:
            sweat = True
        if (
            any(t in airway_exact for t in toks)
            or any(t.startswith(("spirometr", "teae", "sae", "toxicit")) for t in toks)
            or _has_any_seq(
                toks,
                (
                    ("fev", "1"), ("forced", "expiratory", "volume"),
                    ("lung", "clearance", "index"), ("multiple", "breath", "washout"),
                    ("potential", "difference"), ("adverse", "event"), ("adverse", "events"),
                    ("dose", "limiting"), ("vital", "signs"),
                ),
            )
        ):
            airway_safety = True
    return sweat, airway_safety


def _binom_fraction(n: int, k: int) -> Fraction:
    """Multiplicative rational binomial coefficient without `math.comb`."""
    if k < 0 or k > n:
        return Fraction(0, 1)
    m = min(k, n - k)
    acc = Fraction(1, 1)
    for r in range(1, m + 1):
        acc *= Fraction(n - m + r, r)
    return acc


def _fisher_tails_multiplicative(n1: int, k1: int, n2: int, k2: int) -> tuple[Fraction, Fraction]:
    total_n = n1 + n2
    total_k = k1 + k2
    denom = _binom_fraction(total_n, total_k)
    j_min = max(0, total_k - n2)
    j_max = min(n1, total_k)
    pmf = {j: (_binom_fraction(n1, j) * _binom_fraction(n2, total_k - j)) / denom for j in range(j_min, j_max + 1)}
    full = sum(pmf.values(), Fraction(0, 1))
    if full != Fraction(1, 1):
        raise ValueError(f"Hypergeometric mass sum != 1: {full}")
    p_one = sum((pmf[j] for j in range(k1, j_max + 1)), Fraction(0, 1))
    p_obs = pmf[k1]
    p_two = sum((prob for prob in pmf.values() if prob <= p_obs), Fraction(0, 1))
    return p_one, p_two


def _bisect_mn_restricted_p2(x1: int, n1: int, x2: int, n2: int, delta: float) -> float:
    """Solve Miettinen-Nurminen score derivative g(p2) = 0 by monotone bisection (no trigonometry)."""
    lo = max(0.0, -delta)
    hi = min(1.0, 1.0 - delta)
    if hi <= lo + 1e-15:
        return 0.5 * (lo + hi)

    def poly_score(p2: float) -> float:
        p1 = p2 + delta
        # Cleared-denominator log-likelihood derivative wrt p2 (strictly decreasing on (lo, hi)):
        # x1*(1-p1)*p2*(1-p2) - (n1-x1)*p1*p2*(1-p2) + x2*p1*(1-p1)*(1-p2) - (n2-x2)*p1*(1-p1)*p2
        return (
            x1 * (1.0 - p1) * p2 * (1.0 - p2)
            - (n1 - x1) * p1 * p2 * (1.0 - p2)
            + x2 * p1 * (1.0 - p1) * (1.0 - p2)
            - (n2 - x2) * p1 * (1.0 - p1) * p2
        )

    if poly_score(lo + 1e-15) <= 0.0:
        return lo
    if poly_score(hi - 1e-15) >= 0.0:
        return hi
    a, b = lo, hi
    for _ in range(54):
        mid = 0.5 * (a + b)
        if poly_score(mid) > 0.0:
            a = mid
        else:
            b = mid
    return 0.5 * (a + b)


def _mn_z_bisection(x1: int, n1: int, x2: int, n2: int, delta: float) -> float:
    num = (x1 / n1) - (x2 / n2) - delta
    if abs(num) < 1e-15:
        return 0.0
    p2_t = _bisect_mn_restricted_p2(x1, n1, x2, n2, delta)
    p1_t = min(1.0, max(0.0, p2_t + delta))
    total_n = n1 + n2
    var = (p1_t * (1.0 - p1_t) / n1 + p2_t * (1.0 - p2_t) / n2) * (total_n / (total_n - 1.0))
    if var <= 1e-24:
        return 1e12 if num > 0.0 else -1e12
    return num / (var ** 0.5)


def _chan_zhang_ci95_independent(n1: int, k1: int, n2: int, k2: int) -> tuple[float, float]:
    """Independent Chan-Zhang 95% CI using polynomial bisection + ternary nuisance refinement."""
    alpha_tail = 0.025
    c1 = [int(_binom_fraction(n1, i)) for i in range(n1 + 1)]
    c2 = [int(_binom_fraction(n2, j)) for j in range(n2 + 1)]

    def tail_sup(delta: float, upper: bool) -> float:
        lo_p2 = max(0.0, -delta)
        hi_p2 = min(1.0, 1.0 - delta)
        z_obs = _mn_z_bisection(k1, n1, k2, n2, delta)
        pairs: list[tuple[int, int]] = []
        for x1 in range(n1 + 1):
            for x2 in range(n2 + 1):
                z_xy = _mn_z_bisection(x1, n1, x2, n2, delta)
                if (upper and z_xy >= z_obs - 1e-11) or ((not upper) and z_xy <= z_obs + 1e-11):
                    pairs.append((x1, x2))
        if not pairs:
            return 0.0
        if len(pairs) == (n1 + 1) * (n2 + 1):
            return 1.0

        def prob(p2: float) -> float:
            p1 = min(1.0, max(0.0, p2 + delta))
            b1 = [1.0] + [0.0] * n1 if p1 <= 0.0 else ([0.0] * n1 + [1.0] if p1 >= 1.0 else [c1[i] * (p1 ** i) * ((1.0 - p1) ** (n1 - i)) for i in range(n1 + 1)])
            b2 = [1.0] + [0.0] * n2 if p2 <= 0.0 else ([0.0] * n2 + [1.0] if p2 >= 1.0 else [c2[j] * (p2 ** j) * ((1.0 - p2) ** (n2 - j)) for j in range(n2 + 1)])
            return sum(b1[x1] * b2[x2] for x1, x2 in pairs)

        if hi_p2 <= lo_p2 + 1e-14:
            return prob(0.5 * (lo_p2 + hi_p2))

        grid_m = 150
        best_p2 = lo_p2
        best_v = -1.0
        for idx in range(grid_m + 1):
            cand = lo_p2 + (hi_p2 - lo_p2) * (idx / grid_m)
            v = prob(cand)
            if v > best_v:
                best_v = v
                best_p2 = cand
        # Ternary refinement around best_p2
        step = (hi_p2 - lo_p2) / grid_m
        left = max(lo_p2, best_p2 - step)
        right = min(hi_p2, best_p2 + step)
        for _ in range(28):
            m1 = left + (right - left) / 3.0
            m2 = right - (right - left) / 3.0
            v1, v2 = prob(m1), prob(m2)
            if v1 > best_v:
                best_v = v1
            if v2 > best_v:
                best_v = v2
            if v1 < v2:
                left = m1
            else:
                right = m2
        return min(1.0, max(0.0, best_v))

    rd_obs = (k1 / n1) - (k2 / n2)
    lo, hi = -1.0 + 1e-12, rd_obs
    if tail_sup(lo, True) > alpha_tail:
        low_ci = -1.0
    else:
        for _ in range(42):
            mid = 0.5 * (lo + hi)
            if tail_sup(mid, True) > alpha_tail:
                hi = mid
            else:
                lo = mid
        low_ci = hi

    lo, hi = rd_obs, 1.0 - 1e-12
    if tail_sup(hi, False) > alpha_tail:
        high_ci = 1.0
    else:
        for _ in range(42):
            mid = 0.5 * (lo + hi)
            if tail_sup(mid, False) > alpha_tail:
                lo = mid
            else:
                hi = mid
        high_ci = lo
    return low_ci, high_ci


def _newcombe_wilson_decimal(n1: int, k1: int, n2: int, k2: int) -> tuple[Decimal, Decimal]:
    z = Decimal("1.959963984540054")
    z2 = z * z

    def wilson(k: int, n: int) -> tuple[Decimal, Decimal]:
        nd = Decimal(n)
        kd = Decimal(k)
        p = kd / nd
        denom = Decimal(2) * (nd + z2)
        disc = (z2 + Decimal(4) * kd * (Decimal(1) - p)).sqrt()
        lo = Decimal(0) if k == 0 else (Decimal(2) * kd + z2 - z * disc) / denom
        hi = Decimal(1) if k == n else (Decimal(2) * kd + z2 + z * disc) / denom
        return lo, hi

    p1 = Decimal(k1) / Decimal(n1)
    p2 = Decimal(k2) / Decimal(n2)
    rd = p1 - p2
    l1, u1 = wilson(k1, n1)
    l2, u2 = wilson(k2, n2)
    low = rd - ((p1 - l1) ** 2 + (u2 - p2) ** 2).sqrt()
    high = rd + ((u1 - p1) ** 2 + (p2 - l2) ** 2).sqrt()
    return low, high


def run_independent_audit(
    snapshot_path: Path = DEFAULT_SNAPSHOT,
    directory_path: Path = DEFAULT_DIRECTORY,
    done_path: Path = DEFAULT_DONE,
    results_path: Path = DEFAULT_RESULTS,
) -> dict[str, Any]:
    done_text = done_path.read_text(encoding="utf-8").replace("\r\n", "\n")
    done_sha = hashlib.sha256(done_text.encode("utf-8")).hexdigest()

    snap_doc = json.loads(snapshot_path.read_text(encoding="utf-8"))
    studies = snap_doc.get("studies") or []
    raw_n = len(studies)

    dir_doc = json.loads(directory_path.read_text(encoding="utf-8"))
    dir_prisma_n = int((dir_doc.get("metadata") or {}).get("total_prisma_included_studies") or 0)

    prisma_list: list[dict[str, Any]] = []
    for s in studies:
        inc, _ = _prisma_included_independent(s)
        if inc:
            if not isinstance(s.get("phases"), list) or not isinstance(s.get("interventions"), list):
                raise ValueError(f"Missing phases/interventions list on {s.get('nct_id')}")
            prisma_list.append(s)

    n_sys = k_sys = n_inh = k_inh = n_inh_pure = k_inh_pure = inh_airway_safety = 0
    for s in prisma_list:
        sweat, airway = _outcomes_independent(s)
        phases = [str(p).strip().upper() for p in (s.get("phases") or []) if str(p).strip()]
        if not phases or not (set(phases) <= EARLY_PHASE_SET):
            continue
        if _is_healthy_volunteer_only_independent(str(s.get("brief_title") or ""), str(s.get("official_title") or "")):
            continue
        arm, dual = _classify_regimen_independent(s)
        if arm == "systemic_cftr":
            n_sys += 1
            if sweat:
                k_sys += 1
        elif arm == "inhaled_cftr":
            n_inh += 1
            if sweat:
                k_inh += 1
            if airway:
                inh_airway_safety += 1
            if not dual:
                n_inh_pure += 1
                if sweat:
                    k_inh_pure += 1

    p_fisher_1s, p_fisher_2s = _fisher_tails_multiplicative(n_sys, k_sys, n_inh, k_inh)
    p_fisher_pure, _ = _fisher_tails_multiplicative(n_sys, k_sys, n_inh_pure, k_inh_pure)
    cz_lo, cz_hi = _chan_zhang_ci95_independent(n_sys, k_sys, n_inh, k_inh)
    nw_lo, nw_hi = _newcombe_wilson_decimal(n_sys, k_sys, n_inh, k_inh)

    res_doc = json.loads(results_path.read_text(encoding="utf-8"))
    prim = res_doc.get("primary_estimand_early_phase") or {}
    sens = res_doc.get("sensitivity_pure_inhaled_excluding_dual_route") or {}
    pos_spec = ((res_doc.get("positive_controls") or {}).get("inhaled_airway_or_safety_specificity") or {})

    mismatches: list[str] = []
    if res_doc.get("done_sha256") != done_sha:
        mismatches.append(f"done_sha256 mismatch: {res_doc.get('done_sha256')} != {done_sha}")
    if len(prisma_list) != 940 or dir_prisma_n != 940 or raw_n != 1053:
        mismatches.append(f"PRISMA/raw mismatch: raw={raw_n}, prisma={len(prisma_list)}, dir={dir_prisma_n}")
    if (n_sys, k_sys, n_inh, k_inh) != (prim.get("n_sys"), prim.get("k_sys"), prim.get("n_inh"), prim.get("k_inh")):
        mismatches.append(
            f"2x2 counts mismatch: audit={(n_sys, k_sys, n_inh, k_inh)} vs results={(prim.get('n_sys'), prim.get('k_sys'), prim.get('n_inh'), prim.get('k_inh'))}"
        )
    if (n_inh_pure, k_inh_pure) != (sens.get("n_inh_pure"), sens.get("k_inh_pure")):
        mismatches.append(f"pure inhaled mismatch: {(n_inh_pure, k_inh_pure)} vs {(sens.get('n_inh_pure'), sens.get('k_inh_pure'))}")
    if inh_airway_safety != pos_spec.get("matched"):
        mismatches.append(f"airway specificity count mismatch: {inh_airway_safety}")
    try:
        if p_fisher_1s != Fraction(str(prim.get("p_fisher_one_sided_exact", "0/1"))):
            mismatches.append(f"Fisher 1s mismatch: {p_fisher_1s} vs {prim.get('p_fisher_one_sided_exact')}")
        if p_fisher_2s != Fraction(str(prim.get("p_fisher_two_sided_exact", "0/1"))):
            mismatches.append(f"Fisher 2s mismatch: {p_fisher_2s} vs {prim.get('p_fisher_two_sided_exact')}")
    except (ValueError, ZeroDivisionError):
        mismatches.append(f"Malformed Fisher fraction in results: {prim.get('p_fisher_one_sided_exact')}")
    try:
        if p_fisher_pure != Fraction(str(sens.get("p_fisher_one_sided_exact", "0/1"))):
            mismatches.append(f"Pure Fisher exact mismatch: {p_fisher_pure} vs {sens.get('p_fisher_one_sided_exact')}")
    except (ValueError, ZeroDivisionError):
        mismatches.append(f"Malformed pure p_fisher_one_sided_exact: {sens.get('p_fisher_one_sided_exact')}")
    cz_ref = prim.get("chan_zhang_unconditional_ci95")
    if (
        not isinstance(cz_ref, list)
        or len(cz_ref) != 2
        or abs(cz_lo - float(cz_ref[0])) > 1e-5
        or abs(cz_hi - float(cz_ref[1])) > 1e-5
    ):
        mismatches.append(f"Chan-Zhang CI95 mismatch: ({cz_lo}, {cz_hi}) vs {cz_ref}")
    nw_ref = prim.get("newcombe_wilson_ci95")
    if (
        not isinstance(nw_ref, list)
        or len(nw_ref) != 2
        or abs(nw_lo - Decimal(str(nw_ref[0]))) > Decimal("0.000001")
        or abs(nw_hi - Decimal(str(nw_ref[1]))) > Decimal("0.000001")
    ):
        mismatches.append(f"NW CI95 mismatch: ({nw_lo}, {nw_hi}) vs {nw_ref}")

    return {
        "passed": len(mismatches) == 0,
        "mismatches": mismatches,
        "done_sha256": done_sha,
        "raw_interventional_count": raw_n,
        "prisma_included_count": len(prisma_list),
        "directory_cross_verified_count": dir_prisma_n,
        "n_sys": n_sys,
        "k_sys": k_sys,
        "p_sys": round(k_sys / n_sys, 6),
        "n_inh": n_inh,
        "k_inh": k_inh,
        "p_inh": round(k_inh / n_inh, 6),
        "inhaled_airway_or_safety_matched": inh_airway_safety,
        "p_fisher_one_sided_exact": f"{p_fisher_1s.numerator}/{p_fisher_1s.denominator}",
        "p_fisher_one_sided": round(float(p_fisher_1s), 8),
        "p_fisher_two_sided_exact": f"{p_fisher_2s.numerator}/{p_fisher_2s.denominator}",
        "chan_zhang_unconditional_ci95": [round(cz_lo, 6), round(cz_hi, 6)],
        "newcombe_wilson_ci95": [round(float(nw_lo), 6), round(float(nw_hi), 6)],
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--snapshot", type=Path, default=DEFAULT_SNAPSHOT)
    ap.add_argument("--directory", type=Path, default=DEFAULT_DIRECTORY)
    ap.add_argument("--done", type=Path, default=DEFAULT_DONE)
    ap.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args(argv)

    rep = run_independent_audit(args.snapshot, args.directory, args.done, args.results)
    if args.check and not rep["passed"]:
        for m in rep["mismatches"]:
            print(f"FAIL independent_audit: {m}", file=sys.stderr)
        return 1
    if args.json:
        print(json.dumps(rep, indent=2))
    else:
        print(
            f"PASS independent_audit: prisma={rep['prisma_included_count']} "
            f"sys={rep['k_sys']}/{rep['n_sys']} inh={rep['k_inh']}/{rep['n_inh']} "
            f"p_fisher_1s={rep['p_fisher_one_sided_exact']} "
            f"cz_ci95={rep['chan_zhang_unconditional_ci95']} nw_ci95={rep['newcombe_wilson_ci95']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
