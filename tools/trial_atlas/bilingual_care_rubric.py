#!/usr/bin/env python3
"""Bilingual (en & es-MX) Care-Team Rubric Evaluator & Real-Data Benchmark (T-0945 / PRED-0008).

Evaluates automated readability (`readability.py` / Flesch-Kincaid grade),
`cf-ai-safe-use` (`check_safe_output.py` + Spanish `es-MX` safety parity),
openFDA SPL Boxed Warning isolation (`AH-4`), Spanish (`es-MX`) cultural-equity
fidelity, and catalog license boundaries across CFF modulator-era whole-person
care domains and genetic-trial informed consent decision support.

Grounded strictly in verified public sources admitted in `sources/catalog.yaml`:
- `openfda-label`: TRIKAFTA (NDA212273/NDA217660, eff. 2026-06-22), ALYFTREK (NDA218730,
  eff. 2026-08-26), KALYDECO (NDA203188, eff. 2026-03-23), ORKAMBI (NDA206038, eff. 2026-03-24),
  SYMDEKO (NDA210491, eff. 2026-08-21).
- `cffpr-highlights-2025`: GI colonoscopy 48%, CFRD 24%, adult overweight/obesity 29%/14%,
  newly diagnosed Hispanic 13.0%, NTM 12%, depression 26% / anxiety 29%.
- `debray-2026-hepatolcommun` (`DOI 10.1097/HC9.0000000000001061`, `PMC13668289`, `CC BY-NC-ND 4.0`):
  strictly factual metadata only; zero derivative tables/flowcharts or verbatim translation.
- `mcgarry-2025-ijns` (`DOI 10.3390/ijns11020024`, `PMC12015897`, `CC BY 4.0`): NBS consensus equity.
- `kraft-2026-bmcmedresmethodol` (`DOI 10.1186/s12874-026-02885-2`, `PMC13721633`, `CC BY 4.0`):
  informed consent communication framework and decision support.
- `cff-research-priorities` (`kind: resource_pages`, `claim_label: unverified`):
  navigation pointer only; strictly excluded from `*-claims.json`.

100% pure Python 3.10+ standard library, zero network calls, zero LLM calls.

Usage:
    python tools/trial_atlas/bilingual_care_rubric.py --check
    python tools/trial_atlas/bilingual_care_rubric.py --evaluate-live --json
    python tools/trial_atlas/bilingual_care_rubric.py --benchmark --json
    python tools/trial_atlas/bilingual_care_rubric.py --falsify
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import re
import sys
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import care_team_trials as ctt  # noqa: E402
import whole_person_visit_prep as wpvp  # noqa: E402

NLM_DISCLAIMER_EN = "does not evaluate safety or efficacy"
NLM_DISCLAIMER_ES = "no significa que un estudio sea seguro, eficaz ni respaldado por el gobierno"

# Verified openFDA SPL Boxed Warning status (checked 2026-11-19 via openfda-label).
BOXED_WARNING_SPL_GROUND_TRUTH: dict[str, dict[str, Any]] = {
    "TRIKAFTA": {"nda": "NDA212273 / NDA217660", "has_boxed_warning": True, "effective_time": "2026-06-22", "boxed_warning_chars": 3039},
    "ALYFTREK": {"nda": "NDA218730", "has_boxed_warning": True, "effective_time": "2026-08-26", "boxed_warning_chars": 3254},
    "KALYDECO": {"nda": "NDA203188", "has_boxed_warning": False, "effective_time": "2026-03-23", "boxed_warning_chars": 0},
    "ORKAMBI": {"nda": "NDA206038", "has_boxed_warning": False, "effective_time": "2026-03-24", "boxed_warning_chars": 0},
    "SYMDEKO": {"nda": "NDA210491", "has_boxed_warning": False, "effective_time": "2026-08-21", "boxed_warning_chars": 0},
}

# Verbatim public-domain openFDA SPL & CC BY 4.0 Europe PMC excerpts retrieved 2026-11-20.
OPENFDA_SPL_EXCERPTS: dict[str, str] = {
    "TRIKAFTA": (
        "WARNING: DRUG-INDUCED LIVER INJURY AND LIVER FAILURE TRIKAFTA can cause serious and potentially fatal "
        "drug-induced liver injury. Cases of liver failure leading to transplantation and death have been reported "
        "in patients with and without a history of liver disease taking TRIKAFTA, in both clinical trials and the "
        "postmarketing setting. Assess liver function tests (ALT, AST, alkaline phosphatase, and bilirubin) in all "
        "patients prior to initiating TRIKAFTA."
    ),
    "ALYFTREK": (
        "WARNING: DRUG-INDUCED LIVER INJURY AND LIVER FAILURE Elevated transaminases have been observed in patients "
        "treated with ALYFTREK. Cases of serious and potentially fatal drug-induced liver injury and liver failure "
        "were reported in patients who were taking a fixed-dose combination drug containing elexacaftor, tezacaftor, "
        "and ivacaftor, which contains the same or similar active ingredients as ALYFTREK."
    ),
}

EUROPE_PMC_CC_BY_EXCERPTS: dict[str, str] = {
    "mcgarry-2025-ijns": (
        "Newborn screening for cystic fibrosis (CF) has been universal in the US since 2010; however, there is "
        "significant variation among newborn screening algorithms. Systematic reviews were used to develop seven "
        "recommendations for newborn screening program practices to improve timeliness, sensitivity, and equity in "
        "diagnosing infants with CF."
    ),
    "kraft-2026-bmcmedresmethodol": (
        "In clinical trials, the process of informed consent is intended to support prospective participants in "
        "making enrollment decisions that are well informed and aligned with their values and preferences. "
        "Despite the importance of relationship building during the consent process, there is little structured "
        "guidance for researchers."
    ),
}

# Bilingual safety and boundary patterns (en + es-MX parity).
DOSE_PATTERNS: tuple[re.Pattern[str], ...] = (
    re.compile(
        r"\b\d+(?:[.,]\d+)?\s*(?:mg|mcg|ug|ml|mL|milligrams?|miligramos?|micrograms?|microgramos?|"
        r"grams?|gramos?|units?|unidades?|tablets?|tabletas?|capsules?|cápsulas?|puffs?|inhalaciones?)\b",
        re.I,
    ),
    re.compile(
        r"\b(?:once|twice|three\s+times)\s+(?:a\s+day|daily)\b|"
        r"\b(?:una\s+vez|dos\s+veces|tres\s+veces)\s+al\s+día\b|"
        r"\b(?:every|cada)\s+\d+\s+(?:hours?|horas?|hrs?)\b|"
        r"\b(?:BID|TID|QID|PRN)\b",
        re.I,
    ),
)
DIRECTIVE_PATTERNS: tuple[re.Pattern[str], ...] = (
    re.compile(r"\b(?:you\s+should\s+take|stop\s+taking\s+your|you\s+are\s+eligible|you\s+qualify|may\s+be\s+eligible|could\s+qualify)\b", re.I),
    re.compile(r"\b(?:usted\s+(?:debe\s+tomar|califica|es\s+elegible)|puede\s+ser\s+elegible|suspenda\s+su\s+medicamento|le\s+recomiendo\s+tomar)\b", re.I),
)
PROTECTIVE_DIRECTIVE_BEFORE = re.compile(
    r"\b(?:no\s+podemos\s+decirle\s+si|no\s+evalúa\s+si|cannot\s+tell\s+you\s+if|whether|not\s+assess\s+if|nunca\s+indica\s+si)\b[^.!?\n]{0,60}$",
    re.I,
)
PROMISE_PATTERNS: tuple[re.Pattern[str], ...] = (
    re.compile(r"\b(?:proven\s+to\s+cure|guarantees\s+clinical\s+benefit|cures\s+cystic\s+fibrosis)\b", re.I),
    re.compile(r"\b(?:cura\s+la\s+fibrosis\s+quística|garantiza\s+(?:beneficio|cura)|100%\s+eficaz)\b", re.I),
)

# Identifier, URL, phone, and tier patterns for es-MX parity with check_safe_output.py.
ALLOWED_911_LITERAL = re.compile(r"\(911\s+in\s+the\s+US\)", re.I)
BARE_EMERGENCY_RE = re.compile(r"\b(?:911|988|112|999|000)\b(?!\s+(?:people|personas|studies|estudios)\b)", re.I)
ES_IDENTIFIER_PATTERNS: tuple[tuple[str, re.Pattern[str]], ...] = (
    ("date_of_birth", re.compile(r"\b(?:fecha\s+de\s+nacimiento|date\s+of\s+birth|d\.?o\.?b\.?)\b\W{0,3}\s*\S*\d", re.I)),
    ("slash_date", re.compile(r"\b\d{1,2}[/.]\d{1,2}[/.]\d{2,4}\b")),
    ("record_number", re.compile(r"\b(?:mrn|número\s+de\s+expediente|expediente\s+clínico|record\s+number)\b\W{0,3}\s*\S*\d", re.I)),
    ("long_digits", re.compile(r"(?<!NCT)(?<![\d.,])\b\d{7,}\b(?![.,]\d)")),
    ("email", re.compile(r"[\w.+-]+@[\w-]+\.[\w.-]+")),
)
URL_PATTERN = re.compile(r"\bhttps?://\S+|\bwww\.\S+|(?<![@\w.-])[\w-]+\.(?:org|com|gov|net|edu|nhs\.uk|org\.uk|gov\.uk)(?:/\S*)?\b", re.I)
ES_OFFICIAL_SITE_RE = re.compile(r"\b(?:sitio\s+web\s+oficial|página\s+oficial|check\s+(?:them\s+|it\s+)?on\s+the\s+official\s+(?:website|page|site))\b", re.I)
PERCENT_OR_STUDY_RE = re.compile(r"\b\d+(?:[.,]\d+)?\s?%|\bpor\s+ciento\b|\b(?:estudio|estudios|registro)\s+(?:muestra|muestran|encontró|reporta)\b", re.I)
TIER_LABEL_RE = re.compile(r"\bT[1-3]\b(?![\s-]*(?:diabetes|DM|D\b))|\bnot\s+stated\b|\bno\s+declarada\b", re.I)
NCT_TOKEN_RE = re.compile(r"\bnct(?:\d{8})?\b", re.I)

# AH-4 Boxed Warning patterns (order-independent, plurals, synonyms).
ISOLATED_BRAND_RE = re.compile(r"\b(?:kalydeco|orkambi|symdeko)\b", re.I)
BOXED_WARNING_TERM_RE = re.compile(r"\b(?:(?:boxed|black\s+box)\s+warnings?|advertencias?\s+de\s+recuadro)\b", re.I)
AH4_NEGATION_RE = re.compile(
    r"\b(?:"
    r"no\s+(?:fda\s+)?(?:boxed|black\s+box)\s+warnings?"
    r"|(?:boxed|black\s+box)\s+warnings?\s+(?:is|are|remains?)\s+(?:not\s+stated|absent|none)"
    r"|not\s+stated(?:\s+in\s+(?:the\s+)?spl)?"
    r"|no\s+tiene\s+(?:una\s+)?advertencias?\s+de\s+recuadro"
    r"|sin\s+advertencias?\s+de\s+recuadro"
    r"|advertencias?\s+de\s+recuadro\b[^.!?\n]{0,40}\b(?:no\s+declarada|ausente|ninguna)"
    r"|(?:figura\s+como\s+)?no\s+declarada"
    r")\b",
    re.I,
)
AH4_AFFIRMATIVE_BYPASS_RE = re.compile(
    r"(?<!\bno\s)(?<!\bnot\s)(?<!\bnever\s)\b(?:carries|has|with|applies\s+to|lleva|tiene|aplica\s+a)\b(?!\s+(?:no|ninguna|sin)\b)[^.!?\n]{0,45}\b(?:(?:boxed|black\s+box)\s+warnings?|advertencias?\s+de\s+recuadro)\b"
    r"|\b(?:(?:boxed|black\s+box)\s+warnings?|advertencias?\s+de\s+recuadro)\b[^.!?\n]{0,45}\b(?<!\bno\s)(?<!\bnot\s)(?:applies\s+to|aplica\s+a)\b"
    r"|\b(?:exceptions?|excepciones?)\s*:\s*(?:none|ninguna)\b",
    re.I,
)

# CC BY-NC-ND 4.0 Debray et al. (2026) patterns.
DEBRAY_ID_RE = re.compile(
    r"\bdebray\b|debray-2026-hepatolcommun|10\.1097/hc9\.0000000000001061|pmc13668289",
    re.I,
)
DEBRAY_ND_PHRASE_RE = re.compile(
    r"\b(?:"
    r"(?:reproduced|adapted)\s+tables?"
    r"|(?:adapted|diagnostic)\s+flowcharts?"
    r"|rubric\s+prose"
    r"|(?:verbatim\s+)?spanish\s+translations?"
    r"|bulk\s+xml\s+excerpts?"
    r"|tablas?\s+reproducidas?"
    r"|diagramas?\s+adaptados?"
    r"|traducci(?:ón|ones)\s+literal(?:es)?"
    r"|traducción\s+al\s+español"
    r")\b",
    re.I,
)

# Readability fallback regexes matching cf-plain-language-rewrite/scripts/readability.py.
WORD_RE = re.compile(r"[A-Za-z0-9À-ɏ]+(?:['’\-][A-Za-z0-9À-ɏ]+)*%?")
VOWEL_GROUP_RE = re.compile(r"[aeiouy]+")
SENTENCE_ABBREVIATIONS = ("et al.", "e.g.", "i.e.", "etc.", "Dr.", "Mr.", "Mrs.", "Ms.", "vs.", "approx.", "St.", "No.")

DOMAIN_PATTERNS: dict[str, re.Pattern[str]] = {
    "hepatobiliary_dili": re.compile(r"\b(?:alt|ast|bilirubin|bilirrubina|trikafta|alyftrek|hepatobiliar(?:y)?)\b", re.I),
    "newborn_screening_equity": re.compile(r"\b(?:irt|newborn\s+screen\w*|tamiz\s+neonatal|sweat\s+chloride|cloruro\s+en\s+sudor|mcgarry)\b", re.I),
    "non_sputum_infection_ntm": re.compile(r"\b(?:ntm|mnt|sputum|esputo|nct07312734|send-cf)\b", re.I),
    "gi_cfrd_nutrition": re.compile(r"\b(?:48%|24%|colonoscop\w*|cfrd|diabetes\s+relacionada\s+con\s+fq|nutrici[oó]n|nutrition)\b", re.I),
    "mental_pain_repro_telehealth": re.compile(r"\b(?:26%|depress\w*|depresi[oó]n|anxiet\w*|ansiedad|telehealth|telemedicina|reproductiv\w*|pain|dolor)\b", re.I),
    "genetic_trial_informed_consent": re.compile(r"\b(?:kraft|informed\s+consent|consentimiento\s+informado)\b|10\.1186/s12874-026-02885-2", re.I),
}


def _split_sentences_preserving_abbrev(text: str) -> list[str]:
    """Split text into logical sentences without breaking on 'et al.' or decimals."""
    protected = text
    for abbr in SENTENCE_ABBREVIATIONS:
        protected = protected.replace(abbr, abbr.replace(".", "\u0000"))
    protected = re.sub(r"(?<=\d)\.(?=\d)", "\u0000", protected)
    out: list[str] = []
    for line in protected.splitlines():
        for part in re.split(r"(?<=[.!?])[\"')\]]*\s+|(?<=[.!?])$", line):
            s = part.strip().replace("\u0000", ".")
            if s:
                out.append(s)
    return out


def _workspace_sibling_candidates(rel_parts: tuple[str, ...]) -> list[Path]:
    return [
        REPO_ROOT.parent / "cf-skills" / Path(*rel_parts),
        REPO_ROOT.parents[2] / "cf-skills" / Path(*rel_parts),
    ]


def _load_sibling_module(module_name: str, rel_parts: tuple[str, ...]) -> Any:
    for candidate in _workspace_sibling_candidates(rel_parts):
        if candidate.is_file():
            spec = importlib.util.spec_from_file_location(module_name, candidate)
            if spec and spec.loader:
                mod = importlib.util.module_from_spec(spec)
                sys.modules[module_name] = mod
                spec.loader.exec_module(mod)
                return mod
    return None


def _fallback_syllables(word: str) -> int:
    if any(c.isdigit() for c in word):
        return 1
    w = re.sub(r"[^a-z]", "", word.lower())
    if not w:
        return 1
    count = len(VOWEL_GROUP_RE.findall(w))
    if len(w) > 2 and w.endswith("e") and not w.endswith(("le", "ee", "ye")):
        count -= 1
    elif len(w) > 2 and w.endswith("le") and w[-3] in "aeiouy":
        count -= 1
    elif len(w) > 3 and w.endswith("es") and w[-3] not in "aeiouysxzgc" and not w.endswith(("ches", "shes")):
        count -= 1
    elif len(w) > 3 and w.endswith("ed") and w[-3] not in "aeiouytd":
        count -= 1
    return max(1, count)


def _fallback_measure_readability(text: str) -> dict[str, Any]:
    sents = _split_sentences_preserving_abbrev(text)
    per_sent = [WORD_RE.findall(re.sub(r"^\s*(?:[-*+#]|\d+[.)])\s*", "", s)) for s in sents]
    words = [w for ws in per_sent for w in ws]
    if not words or not sents:
        return {"words": 0, "sentences": 0, "syllables": 0, "average_sentence_length": 0.0, "flesch_kincaid_grade": 0.0}
    n_words, n_sents = len(words), len(sents)
    n_syll = sum(_fallback_syllables(w) for w in words)
    grade = 0.39 * (n_words / n_sents) + 11.8 * (n_syll / n_words) - 15.59
    return {
        "words": n_words,
        "sentences": n_sents,
        "syllables": n_syll,
        "average_sentence_length": round(n_words / n_sents, 1),
        "flesch_kincaid_grade": round(grade, 1),
    }


def measure_readability(text: str) -> dict[str, Any]:
    """Measure Flesch-Kincaid grade via cf-plain-language-rewrite/scripts/readability.py."""
    mod = _load_sibling_module(
        "cf_readability_mod",
        (".claude", "skills", "cf-plain-language-rewrite", "scripts", "readability.py"),
    )
    if mod is not None and hasattr(mod, "measure"):
        res = dict(mod.measure(text))
        res["source_tool"] = "cf-plain-language-rewrite/scripts/readability.py"
        return res
    res = _fallback_measure_readability(text)
    res["source_tool"] = "builtin-fallback"
    return res


def run_cf_ai_safe_use_check(
    text: str,
    lang: str = "en",
    require_visit_prep_boilerplate: bool = True,
) -> list[dict[str, Any]]:
    """Run check_safe_output.py on English text plus full bilingual es-MX parity checks."""
    flags: list[dict[str, Any]] = []
    is_es = lang.lower().startswith("es")

    if not is_es:
        mod = _load_sibling_module(
            "cf_check_safe_output_mod",
            (".claude", "skills", "cf-ai-safe-use", "scripts", "check_safe_output.py"),
        )
        if mod is not None and hasattr(mod, "check"):
            for fl in mod.check(text):
                if not require_visit_prep_boilerplate and fl.rule in ("care-team", "tier"):
                    continue
                flags.append({"rule": f"cf_ai_safe_use:{fl.rule}", "line": fl.line, "detail": fl.text, "why": fl.why})

    for pat in DOSE_PATTERNS:
        for m in pat.finditer(text):
            flags.append({"rule": "bilingual_safety:dose", "detail": m.group(0)})

    for pat in DIRECTIVE_PATTERNS:
        for m in pat.finditer(text):
            before = text[max(0, m.start() - 70):m.start()]
            if not PROTECTIVE_DIRECTIVE_BEFORE.search(before):
                flags.append({"rule": "bilingual_safety:directive", "detail": m.group(0)})

    for pat in PROMISE_PATTERNS:
        for m in pat.finditer(text):
            flags.append({"rule": "bilingual_safety:clinical_promise", "detail": m.group(0)})

    stripped_911 = ALLOWED_911_LITERAL.sub(" ", text)
    for m_emerg in BARE_EMERGENCY_RE.finditer(stripped_911):
        flags.append({"rule": "bilingual_safety:bare_emergency_number", "detail": m_emerg.group(0)})

    if is_es:
        for id_rule, pat in ES_IDENTIFIER_PATTERNS:
            for m in pat.finditer(text):
                flags.append({"rule": f"bilingual_safety:identifier:{id_rule}", "detail": m.group(0)})
        for para in re.split(r"\n\s*\n", text):
            if URL_PATTERN.search(para) and not ES_OFFICIAL_SITE_RE.search(para):
                flags.append({"rule": "bilingual_safety:unchecked_link", "detail": URL_PATTERN.search(para).group(0)})
            if require_visit_prep_boilerplate and PERCENT_OR_STUDY_RE.search(para) and not TIER_LABEL_RE.search(para):
                flags.append({"rule": "bilingual_safety:missing_tier_es", "detail": PERCENT_OR_STUDY_RE.search(para).group(0)})

    if require_visit_prep_boilerplate:
        if "(911 in the US)" not in text:
            flags.append({"rule": "bilingual_safety:missing_911_in_us", "detail": "Missing exact '(911 in the US)' notice."})
        lower = text.lower()
        if is_es:
            if "esto no es un diagnóstico" not in lower:
                flags.append({"rule": "bilingual_safety:missing_non_diagnostic_es", "detail": "Missing 'Esto no es un diagnóstico.'"})
            if "equipo de atención de fq" not in lower and "equipo de fq" not in lower:
                flags.append({"rule": "bilingual_safety:missing_care_team_es", "detail": "Missing Spanish CF care-team pointer."})
        else:
            if "not a diagnosis" not in lower and "not medical advice" not in lower:
                flags.append({"rule": "bilingual_safety:missing_non_diagnostic_en", "detail": "Missing non-clinical boundary."})
            if "care team" not in lower:
                flags.append({"rule": "bilingual_safety:missing_care_team_en", "detail": "Missing English CF care-team pointer."})

        if NCT_TOKEN_RE.search(text):
            if is_es:
                if NLM_DISCLAIMER_ES not in lower and NLM_DISCLAIMER_EN not in lower:
                    flags.append({"rule": "bilingual_safety:missing_nlm_disclaimer_es", "detail": "Missing NLM disclaimer for NCT listing."})
            elif NLM_DISCLAIMER_EN not in lower:
                flags.append({"rule": "bilingual_safety:missing_nlm_disclaimer_en", "detail": "Missing NLM disclaimer for NCT listing."})

    seen: set[tuple[str, str]] = set()
    deduped: list[dict[str, Any]] = []
    for item in flags:
        key = (str(item.get("rule", "")), str(item.get("detail", "")))
        if key not in seen:
            seen.add(key)
            deduped.append(item)
    return deduped


def check_ah4_boxed_warning_isolation(text: str) -> list[dict[str, str]]:
    """Verify AH-4 Boxed Warning isolation per logical sentence (order-independent)."""
    violations: list[dict[str, str]] = []
    for sent in _split_sentences_preserving_abbrev(text):
        if not ISOLATED_BRAND_RE.search(sent) or not BOXED_WARNING_TERM_RE.search(sent):
            continue
        is_negated = bool(AH4_NEGATION_RE.search(sent))
        has_affirmative_or_bypass = bool(AH4_AFFIRMATIVE_BYPASS_RE.search(sent))
        if (not is_negated) or has_affirmative_or_bypass:
            violations.append({"rule": "ah4:boxed_warning_extrapolation", "detail": sent[:180]})
    return violations


def check_es_mx_cultural_equity(text: str, lang: str) -> list[dict[str, str]]:
    """Verify Mexican-Spanish (es-MX) cultural-equity and interpreter guidance."""
    if not lang.lower().startswith("es"):
        return []
    lower = text.lower()
    issues: list[dict[str, str]] = []
    if "intérprete médico calificado" not in lower and "intérprete calificado" not in lower:
        issues.append(
            {
                "rule": "es_mx_equity:missing_qualified_interpreter",
                "detail": "Missing guidance to request a qualified medical interpreter ('intérprete médico calificado').",
            }
        )
    if "esto no es un diagnóstico" not in lower:
        issues.append({"rule": "es_mx_equity:missing_non_diagnostic_closing", "detail": "Missing 'Esto no es un diagnóstico.'"})
    if "(911 in the us)" not in lower:
        issues.append({"rule": "es_mx_equity:missing_911_in_us", "detail": "Missing exact '(911 in the US)' notice."})
    return issues


def check_license_and_pointer_boundaries(
    text: str,
    claims_json_payload: dict[str, Any] | list[Any] | None = None,
) -> list[dict[str, str]]:
    """Enforce CC BY-NC-ND 4.0 (debray-2026) and cff-research-priorities claims-JSON isolation."""
    violations: list[dict[str, str]] = []
    for sent in _split_sentences_preserving_abbrev(text):
        if DEBRAY_ID_RE.search(sent) and DEBRAY_ND_PHRASE_RE.search(sent):
            violations.append({"rule": "license:cc_by_nc_nd_derivative_debray_2026", "detail": sent[:180]})

    if claims_json_payload is not None:
        serialized = json.dumps(claims_json_payload, ensure_ascii=False)
        if "cff-research-priorities" in serialized.lower():
            violations.append(
                {
                    "rule": "provenance:cff_research_priorities_in_claims_json",
                    "detail": "cff-research-priorities must stay out of *-claims.json (unverified navigation pointer only).",
                }
            )
        if DEBRAY_ID_RE.search(serialized) and DEBRAY_ND_PHRASE_RE.search(serialized):
            violations.append({"rule": "license:cc_by_nc_nd_derivative_debray_2026", "detail": "Debray ND violation in claims payload."})
    return violations


def detect_domain_coverage(text: str) -> list[str]:
    """Detect which CFF whole-person & informed-consent domains are addressed with word boundaries."""
    return [dom for dom, pat in DOMAIN_PATTERNS.items() if pat.search(text)]


def evaluate_guide_text(
    text: str,
    lang: str = "en",
    max_grade: float = 8.0,
    enforce_readability: bool = False,
    require_es_interpreter: bool = True,
    require_visit_prep_boilerplate: bool = True,
    claims_json_payload: dict[str, Any] | list[Any] | None = None,
) -> dict[str, Any]:
    """Run all bilingual care-team rubric checks on a guide, SPL excerpt, or academic abstract."""
    safety_flags = run_cf_ai_safe_use_check(
        text, lang=lang, require_visit_prep_boilerplate=require_visit_prep_boilerplate
    )
    ah4_flags = check_ah4_boxed_warning_isolation(text)
    es_flags = check_es_mx_cultural_equity(text, lang=lang) if (require_es_interpreter and require_visit_prep_boilerplate) else []
    license_flags = check_license_and_pointer_boundaries(text, claims_json_payload=claims_json_payload)
    readability = measure_readability(text)
    grade = float(readability.get("flesch_kincaid_grade", 0.0))

    readability_flags: list[dict[str, str]] = []
    if enforce_readability and grade > max_grade:
        readability_flags.append(
            {"rule": "readability:above_max_grade", "detail": f"Flesch-Kincaid grade {grade} exceeds max_grade {max_grade}"}
        )

    all_flags = safety_flags + ah4_flags + es_flags + license_flags + readability_flags
    failed_rules = sorted({str(f["rule"]) for f in all_flags})

    return {
        "lang": lang,
        "rubric_pass": len(all_flags) == 0,
        "failed_rules": failed_rules,
        "flags": all_flags,
        "readability": {
            "words": readability["words"],
            "sentences": readability["sentences"],
            "syllables": readability["syllables"],
            "average_sentence_length": readability["average_sentence_length"],
            "flesch_kincaid_grade": grade,
            "within_grade_6": grade <= 6.0,
            "within_grade_8": grade <= 8.0,
            "source_tool": readability.get("source_tool", "unknown"),
            "calibration_note": (
                "Flesch-Kincaid is English-calibrated; Spanish estimates are surface length heuristics only."
                if lang.lower().startswith("es")
                else "English Flesch-Kincaid grade estimate (within ~1 grade; measures length, not clarity)."
            ),
        },
        "domains_covered": detect_domain_coverage(text),
    }


def _render_genetic_consent_guide(lang: str = "en") -> str:
    """Render a CC BY 4.0 genetic-trial informed consent decision-support guide (Kraft et al. 2026)."""
    if lang == "es-MX":
        return "\n".join(
            [
                "# Guía de Consentimiento Informado para Ensayos Genéticos de FQ [T1]",
                "",
                "Solo para preparación informativa de consulta (Tier T1). No es consejo médico ni diagnóstico. "
                "Si presenta un síntoma nuevo, grave, que empeora rápido o que le preocupa —especialmente en un bebé o niño—, "
                "no espere a la consulta: llame ahora a la línea de urgencias de su equipo de FQ o a los servicios locales de emergencia (911 in the US).",
                "",
                "[T1] Basado en el marco de comunicación de consentimiento informado de Kraft et al. (2026, "
                "DOI 10.1186/s12874-026-02885-2, licencia CC BY 4.0). La inscripción en un ensayo clínico siempre es voluntaria. "
                "Estar listado en el registro no significa que un estudio sea seguro, eficaz ni respaldado por el gobierno.",
                "",
                "## Preguntas Clave para Revisar con el Equipo de Estudio y su Equipo de Atención de FQ",
                "- ¿Qué implica este estudio genético (por ejemplo, vectores virales o ARNm inhalado en registros NCT) frente a nuestro cuidado actual?",
                "- ¿Qué riesgos conocidos, visitas presenciales y seguimientos a largo plazo contempla el protocolo antes de decidir?",
                "- ¿Podemos contar con un intérprete médico calificado y tiempo adicional para revisar el consentimiento en casa?",
                "",
                "Por favor consulte esto con su equipo de atención de FQ antes de cambiar cualquier tratamiento. "
                "Esto no es un diagnóstico.",
            ]
        )
    return "\n".join(
        [
            "# Genetic-Trial Informed Consent Decision-Support Questions [T1]",
            "",
            "Educational clinic-visit preparation only (Tier T1). Independent open-source research, "
            "not medical advice or diagnosis. If a symptom is new, severe, getting worse fast or worrying, "
            "above all in a baby or child: Do not wait for this sheet. Call your CF team's urgent or "
            "out-of-hours line, or local emergency services (911 in the US), now.",
            "",
            "[T1] Informed by Kraft et al. (2026, DOI 10.1186/s12874-026-02885-2, CC BY 4.0) communication framework "
            "for research enrollment decisions. Participation in any NCT study is voluntary. "
            "Please check them on the official website with your CF care team. "
            "ClinicalTrials.gov does not evaluate safety or efficacy; listing on the registry does not denote government endorsement.",
            "",
            "## Key Questions for Your CF Care Team and Study Investigators",
            "- What are the possible benefits, known risks, and daily visit commitments compared with standard CF care?",
            "- Could joining or waiting affect future eligibility for other genetic or modulator trials?",
            "- Can we take the consent form home and request a qualified medical interpreter before deciding?",
            "",
            "Please check this with your CF care team before changing anything. "
            "Do not wait for the visit. This is not a diagnosis.",
        ]
    )


def evaluate_live_lab_guides() -> dict[str, Any]:
    """Audit live whole_person_visit_prep.py, care_team_trials.py (limit=8 & limit=200), and claims JSONs."""
    guides: list[dict[str, Any]] = []

    for lang in ("en", "es-MX"):
        pkt = wpvp.build_whole_person_packet(domain="all", lang=lang)
        md = wpvp.render_whole_person_markdown(pkt)
        res = evaluate_guide_text(md, lang=lang, enforce_readability=False, require_es_interpreter=(lang == "es-MX"))
        res["guide_id"] = f"whole_person_visit_prep:all:{lang}"
        guides.append(res)

        for dom in wpvp.DOMAIN_ORDER:
            dom_pkt = wpvp.build_whole_person_packet(domain=dom, lang=lang)
            dom_md = wpvp.render_whole_person_markdown(dom_pkt)
            # Do NOT mask missing interpreter guidance on single-domain es-MX guides: report it honestly.
            dom_res = evaluate_guide_text(dom_md, lang=lang, enforce_readability=False, require_es_interpreter=(lang == "es-MX"))
            dom_res["guide_id"] = f"whole_person_visit_prep:{dom}:{lang}"
            guides.append(dom_res)

    meta, included = ctt.load_prisma_included_studies()
    snap_ts = str(meta.get("data_timestamp", "2026-10-14"))
    for limit in (8, 200):
        recruiting = ctt.filter_registered_trials(included, status="RECRUITING", limit=limit)
        for lang in ("en", "es-MX"):
            sheet = ctt.render_care_team_sheet(
                recruiting,
                {
                    "status": "RECRUITING",
                    "modality": "All modalities" if lang == "en" else "Todas las modalidades",
                    "phase": "All phases" if lang == "en" else "Todas las fases",
                    "query": "None" if lang == "en" else "Ninguno",
                },
                snap_ts,
                lang=lang,
            )
            res = evaluate_guide_text(sheet, lang=lang, enforce_readability=False, require_es_interpreter=(lang == "es-MX"))
            res["guide_id"] = f"care_team_trials:RECRUITING:limit={limit}:{lang}"
            guides.append(res)

    claims_dir = REPO_ROOT / "landscape"
    claims_leaks: list[str] = []
    claims_files_scanned = 0
    if claims_dir.is_dir():
        for cf in sorted(claims_dir.glob("*-claims.json")):
            claims_files_scanned += 1
            if "cff-research-priorities" in cf.read_text(encoding="utf-8").lower():
                claims_leaks.append(cf.name)

    return {
        "tool": "bilingual_care_rubric",
        "snapshot_prisma_included": len(included),
        "claims_files_scanned": claims_files_scanned,
        "cff_research_priorities_leaks": claims_leaks,
        "spl_boxed_warning_ground_truth": BOXED_WARNING_SPL_GROUND_TRUTH,
        "guides": guides,
    }


def _build_benchmark_cases() -> list[dict[str, Any]]:
    """Construct a 16-case benchmark anchored in measured public SPL/EuropePMC/CT.gov data + controls."""
    cases: list[dict[str, Any]] = []

    # 1-2: Measured public-domain openFDA SPL Boxed Warning texts (TRIKAFTA & ALYFTREK).
    for brand, spl_text in OPENFDA_SPL_EXCERPTS.items():
        cases.append(
            {
                "case_id": f"public_openfda_spl_{brand.lower()}_bw",
                "text": spl_text,
                "lang": "en",
                "enforce_readability": False,
                "require_es_interpreter": False,
                "require_visit_prep_boilerplate": False,
                "expected_pass": True,
                "expected_failed_rules": [],
            }
        )

    # 3-4: Measured Europe PMC CC BY 4.0 published abstracts (McGarry 2025 & Kraft 2026).
    for src_id, abstract_text in EUROPE_PMC_CC_BY_EXCERPTS.items():
        cases.append(
            {
                "case_id": f"public_epmc_cc_by_{src_id}",
                "text": abstract_text,
                "lang": "en",
                "enforce_readability": False,
                "require_es_interpreter": False,
                "require_visit_prep_boilerplate": False,
                "expected_pass": True,
                "expected_failed_rules": [],
            }
        )

    # 5-6: Real whole_person_visit_prep.py all-domain bilingual outputs (en & es-MX).
    for lang in ("en", "es-MX"):
        pkt = wpvp.build_whole_person_packet(domain="all", lang=lang)
        cases.append(
            {
                "case_id": f"live_wpvp_all_{lang}",
                "text": wpvp.render_whole_person_markdown(pkt),
                "lang": lang,
                "enforce_readability": False,
                "require_es_interpreter": (lang == "es-MX"),
                "require_visit_prep_boilerplate": True,
                "expected_pass": True,
                "expected_failed_rules": [],
            }
        )

    # 7-8: Real whole_person_visit_prep.py single-domain es-MX gap detection (lacks interpreter).
    pkt_es_dili = wpvp.build_whole_person_packet(domain="hepatobiliary_dili", lang="es-MX")
    cases.append(
        {
            "case_id": "live_wpvp_single_domain_dili_es_missing_interpreter",
            "text": wpvp.render_whole_person_markdown(pkt_es_dili),
            "lang": "es-MX",
            "enforce_readability": False,
            "require_es_interpreter": True,
            "require_visit_prep_boilerplate": True,
            "expected_pass": False,
            "expected_failed_rules": ["es_mx_equity:missing_qualified_interpreter"],
        }
    )

    pkt_en_nbs = wpvp.build_whole_person_packet(domain="newborn_screening_equity", lang="en")
    cases.append(
        {
            "case_id": "live_wpvp_nbs_en_exceeds_grade_8",
            "text": wpvp.render_whole_person_markdown(pkt_en_nbs),
            "lang": "en",
            "enforce_readability": True,
            "max_grade": 8.0,
            "require_es_interpreter": False,
            "require_visit_prep_boilerplate": True,
            "expected_pass": False,
            "expected_failed_rules": ["readability:above_max_grade"],
        }
    )

    # 9-10: Real care_team_trials.py sheets (limit=8, en & es-MX) from ClinicalTrials.gov snapshot.json.
    meta, included = ctt.load_prisma_included_studies()
    recruiting_8 = ctt.filter_registered_trials(included, status="RECRUITING", limit=8)
    snap_ts = str(meta.get("data_timestamp", "2026-10-14"))
    for lang in ("en", "es-MX"):
        cases.append(
            {
                "case_id": f"live_care_team_trials_limit8_{lang}",
                "text": ctt.render_care_team_sheet(
                    recruiting_8,
                    {"status": "RECRUITING", "modality": "All", "phase": "All", "query": "None"},
                    snap_ts,
                    lang=lang,
                ),
                "lang": lang,
                "enforce_readability": False,
                "require_es_interpreter": (lang == "es-MX"),
                "require_visit_prep_boilerplate": True,
                "expected_pass": True,
                "expected_failed_rules": [],
            }
        )

    # 11: True negative AH-4 phrasing (KALYDECO has no Boxed Warning).
    base_en = wpvp.render_whole_person_markdown(wpvp.build_whole_person_packet(domain="hepatobiliary_dili", lang="en"))
    cases.append(
        {
            "case_id": "true_neg_ah4_kalydeco_no_boxed_warning",
            "text": base_en + "\n\n[T1] KALYDECO has no Boxed Warning in its US FDA SPL.",
            "lang": "en",
            "enforce_readability": False,
            "require_es_interpreter": False,
            "require_visit_prep_boilerplate": True,
            "expected_pass": True,
            "expected_failed_rules": [],
        }
    )

    # 12: Genetic-trial informed consent decision support (Kraft et al. 2026, CC BY 4.0, es-MX).
    cases.append(
        {
            "case_id": "genetic_consent_kraft_2026_es",
            "text": _render_genetic_consent_guide(lang="es-MX"),
            "lang": "es-MX",
            "enforce_readability": False,
            "require_es_interpreter": True,
            "require_visit_prep_boilerplate": True,
            "expected_pass": True,
            "expected_failed_rules": [],
        }
    )

    # 13-16: Poisoned adversarial boundary controls (AH-4 order inversion + et al., Debray ND + et al., claims leak, (911) + es-MX dose).
    cases.append(
        {
            "case_id": "violation_ah4_inverted_order_and_et_al",
            "text": base_en + "\n\n[T1] A Black Box Warning applies to KALYDECO (see Debray et al. 2026), exceptions: none.",
            "lang": "en",
            "enforce_readability": False,
            "require_es_interpreter": False,
            "require_visit_prep_boilerplate": True,
            "expected_pass": False,
            "expected_failed_rules": ["ah4:boxed_warning_extrapolation"],
        }
    )

    base_es = wpvp.render_whole_person_markdown(wpvp.build_whole_person_packet(domain="newborn_screening_equity", lang="es-MX"))
    cases.append(
        {
            "case_id": "violation_cc_by_nc_nd_debray_et_al",
            "text": base_es + "\n\n[T1] Debray et al. (2026) adapted flowchart and reproduced tables.",
            "lang": "es-MX",
            "enforce_readability": False,
            "require_es_interpreter": True,
            "require_visit_prep_boilerplate": True,
            "expected_pass": False,
            "expected_failed_rules": ["license:cc_by_nc_nd_derivative_debray_2026"],
        }
    )

    cases.append(
        {
            "case_id": "violation_cff_research_priorities_in_claims_json",
            "text": base_en,
            "lang": "en",
            "enforce_readability": False,
            "require_es_interpreter": False,
            "require_visit_prep_boilerplate": True,
            "claims_json_payload": {"claims": [{"source_id": "cff-research-priorities", "claim": "unverified pointer"}]},
            "expected_pass": False,
            "expected_failed_rules": ["provenance:cff_research_priorities_in_claims_json"],
        }
    )

    cases.append(
        {
            "case_id": "violation_es_mx_dose_directive_and_bare_911",
            "text": base_es + "\n\n[T1] Usted debe tomar 290 microgramos una vez al día o llamar al (911).",
            "lang": "es-MX",
            "enforce_readability": False,
            "require_es_interpreter": True,
            "require_visit_prep_boilerplate": True,
            "expected_pass": False,
            "expected_failed_rules": [
                "bilingual_safety:bare_emergency_number",
                "bilingual_safety:directive",
                "bilingual_safety:dose",
            ],
        }
    )

    return cases


def run_benchmark() -> dict[str, Any]:
    """Run the 16-case PRED-0008 benchmark and compute exact accuracy, precision, and recall."""
    cases = _build_benchmark_cases()
    results: list[dict[str, Any]] = []
    tp = fp = tn = fn = 0

    for c in cases:
        ev = evaluate_guide_text(
            c["text"],
            lang=c["lang"],
            max_grade=float(c.get("max_grade", 8.0)),
            enforce_readability=bool(c.get("enforce_readability", False)),
            require_es_interpreter=bool(c.get("require_es_interpreter", False)),
            require_visit_prep_boilerplate=bool(c.get("require_visit_prep_boilerplate", True)),
            claims_json_payload=c.get("claims_json_payload"),
        )
        expected_pass = bool(c["expected_pass"])
        actual_pass = bool(ev["rubric_pass"])
        expected_rules = set(c.get("expected_failed_rules", []))
        actual_rules = set(ev["failed_rules"])
        rules_match = expected_rules.issubset(actual_rules) if not expected_pass else (actual_pass is True)
        correct = (actual_pass == expected_pass) and rules_match

        if not expected_pass and not actual_pass:
            tp += 1
        elif expected_pass and not actual_pass:
            fp += 1
        elif expected_pass and actual_pass:
            tn += 1
        else:
            fn += 1

        results.append(
            {
                "case_id": c["case_id"],
                "lang": c["lang"],
                "expected_pass": expected_pass,
                "actual_pass": actual_pass,
                "failed_rules": ev["failed_rules"],
                "flesch_kincaid_grade": ev["readability"]["flesch_kincaid_grade"],
                "correct": correct,
            }
        )

    total = len(results)
    exact_matches = sum(1 for r in results if r["correct"])
    accuracy = round(exact_matches / total, 4) if total else 0.0
    precision = round(tp / (tp + fp), 4) if (tp + fp) else 1.0
    recall = round(tp / (tp + fn), 4) if (tp + fn) else 1.0

    live_audit = evaluate_live_lab_guides()
    return {
        "prediction_id": "PRED-0008",
        "total_cases": total,
        "exact_matches": exact_matches,
        "accuracy": accuracy,
        "violation_precision": precision,
        "violation_recall": recall,
        "confusion": {"tp": tp, "fp": fp, "tn": tn, "fn": fn},
        "live_guide_summary": [
            {
                "guide_id": g["guide_id"],
                "rubric_pass": g["rubric_pass"],
                "failed_rules": g["failed_rules"],
                "flesch_kincaid_grade": g["readability"]["flesch_kincaid_grade"],
                "within_grade_8": g["readability"]["within_grade_8"],
                "within_grade_6": g["readability"]["within_grade_6"],
            }
            for g in live_audit["guides"]
        ],
        "cases": results,
    }


def run_falsification_suite() -> dict[str, Any]:
    """Run 4 deterministic negative controls (--falsify)."""
    controls: list[dict[str, Any]] = []
    base_en = wpvp.render_whole_person_markdown(wpvp.build_whole_person_packet(domain="hepatobiliary_dili", lang="en"))
    base_es = wpvp.render_whole_person_markdown(wpvp.build_whole_person_packet(domain="newborn_screening_equity", lang="es-MX"))

    # Control 1: AH-4 Boxed Warning extrapolation (inverted order + plural) to ORKAMBI/SYMDEKO fails.
    r1 = evaluate_guide_text(base_en + "\n\nBoxed Warnings apply to ORKAMBI and SYMDEKO for hepatic failure.", lang="en")
    controls.append(
        {
            "control": "ah4_orkambi_boxed_warning_blocked",
            "passed": (not r1["rubric_pass"]) and ("ah4:boxed_warning_extrapolation" in r1["failed_rules"]),
            "failed_rules": r1["failed_rules"],
        }
    )

    # Control 2: CC BY-NC-ND 4.0 derivative/translation violation on debray-2026-hepatolcommun with 'et al.' fails.
    r2 = evaluate_guide_text(
        base_es + "\n\nDebray et al. (2026) adapted flowchart and Spanish translation of Debray.",
        lang="es-MX",
        require_es_interpreter=True,
    )
    controls.append(
        {
            "control": "cc_by_nc_nd_debray_derivative_blocked",
            "passed": (not r2["rubric_pass"]) and ("license:cc_by_nc_nd_derivative_debray_2026" in r2["failed_rules"]),
            "failed_rules": r2["failed_rules"],
        }
    )

    # Control 3: cff-research-priorities leaked into *-claims.json fails.
    r3 = evaluate_guide_text(
        base_en,
        lang="en",
        claims_json_payload={"claims": [{"id": "C-1", "source_id": "cff-research-priorities"}]},
    )
    controls.append(
        {
            "control": "cff_research_priorities_claims_leak_blocked",
            "passed": (not r3["rubric_pass"]) and ("provenance:cff_research_priorities_in_claims_json" in r3["failed_rules"]),
            "failed_rules": r3["failed_rules"],
        }
    )

    # Control 4: Missing es-MX interpreter & parenthesized (911) fails.
    poisoned_es = base_es.replace("intérprete médico calificado", "personal").replace("(911 in the US)", "(911)")
    r4 = evaluate_guide_text(poisoned_es, lang="es-MX", require_es_interpreter=True)
    controls.append(
        {
            "control": "es_mx_cultural_equity_and_911_blocked",
            "passed": (
                (not r4["rubric_pass"])
                and ("es_mx_equity:missing_qualified_interpreter" in r4["failed_rules"])
                and ("bilingual_safety:bare_emergency_number" in r4["failed_rules"])
            ),
            "failed_rules": r4["failed_rules"],
        }
    )

    return {
        "tool": "bilingual_care_rubric",
        "falsify_passed": all(c["passed"] for c in controls),
        "negative_controls": controls,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Run fast deterministic self-check")
    parser.add_argument("--evaluate-live", action="store_true", help="Audit live whole_person_visit_prep and care_team_trials guides")
    parser.add_argument("--benchmark", action="store_true", help="Run 16-case PRED-0008 rubric benchmark")
    parser.add_argument("--falsify", action="store_true", help="Run 4-control negative falsification suite")
    parser.add_argument("--json", action="store_true", help="Emit JSON output")
    args = parser.parse_args(argv)

    if args.falsify:
        rep = run_falsification_suite()
        if args.json:
            print(json.dumps(rep, indent=2, ensure_ascii=False))
        else:
            print(f"PASS falsify bilingual_care_rubric ({len(rep['negative_controls'])}/4)" if rep["falsify_passed"] else "FAIL falsify bilingual_care_rubric")
        return 0 if rep["falsify_passed"] else 1

    if args.evaluate_live:
        live = evaluate_live_lab_guides()
        print(json.dumps(live, indent=2, ensure_ascii=False))
        return 0

    if args.benchmark:
        bm = run_benchmark()
        if args.json:
            print(json.dumps(bm, indent=2, ensure_ascii=False))
        else:
            print(
                f"PRED-0008 benchmark: {bm['exact_matches']}/{bm['total_cases']} "
                f"accuracy={bm['accuracy']} precision={bm['violation_precision']} recall={bm['violation_recall']}"
            )
        return 0 if bm["accuracy"] >= 0.90 else 1

    bm = run_benchmark()
    f_rep = run_falsification_suite()
    ok = bm["accuracy"] >= 0.90 and bm["violation_precision"] >= 0.90 and bm["violation_recall"] >= 0.90 and f_rep["falsify_passed"]
    if args.json:
        print(json.dumps({"benchmark": bm, "falsify": f_rep, "check_passed": ok}, indent=2, ensure_ascii=False))
    else:
        print(
            f"{'PASS' if ok else 'FAIL'} bilingual_care_rubric: "
            f"benchmark={bm['exact_matches']}/{bm['total_cases']} (acc={bm['accuracy']}, "
            f"prec={bm['violation_precision']}, rec={bm['violation_recall']}), "
            f"falsify={len(f_rep['negative_controls'])}/4"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
