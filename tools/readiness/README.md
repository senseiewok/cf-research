# CF Machine-Learning Paper Readiness Label Extractor (`tools/readiness`)

Deterministic sentence-grounded extraction of clinical readiness labels from Cystic Fibrosis machine-learning and artificial intelligence papers (board row T-0931).

## Motivation

As observed in our medical AI discovery audits (`.lab-files/scratch/ai-medicine-scan/discovery-round-2026-11-15.md`):
- Most machine learning applications in Cystic Fibrosis are retrospective, single-center, and lack independent external validation across geographic regions or patient registries.
- Secondary explainers, chatbots, and model summaries routinely overstate readiness—conflating in vitro cell-line models (e.g. HEK293T readthrough assays) with clinical trials, or presenting single-center retrospective scores as deployable clinical tools.
- TRIPOD+AI, CLAIM, and PRISMA-AI reporting standards mandate transparent disclosure of cohort size, split structure, external validation, and primary performance metrics.

`tools/readiness` extracts these four core parameters directly from paper abstracts and methods, anchors each extracted parameter with its exact verbatim source sentence, assigns a conservative epistemic readiness tier (`T0`–`T3`), and attaches anti-overclaim warnings.

## Extracted Fields

| Field | Description | Values / Types | Source Grounding |
| :--- | :--- | :--- | :--- |
| `cohort_size` | Number of human patients/individuals in study | `int` or `"not stated"` | Exact sentence mentioning patient count |
| `split_type` | Partitioning methodology | `patient_level`, `sample_or_record_level`, `unsupervised_or_zero_shot`, `not_stated` | Sentence describing split or cross-validation |
| `external_validation` | Independent external cohort testing | `external_multicenter_or_independent`, `absent_single_center_or_internal_only`, `not_stated` | Sentence stating single-center vs independent cohort |
| `headline_metric` | Primary performance metric | String (e.g., `AUCs of 0.75-0.80`, `AUC-ROC of 0.89`, `-1.40±1.53 points`) | Sentence reporting metric |
| `readiness_tier` | Standardized readiness level | `T0_PRECLINICAL_BENCHMARK`, `T1_RETROSPECTIVE_OR_UNVALIDATED`, `T2_INTERNALLY_VALIDATED`, `T3_EXTERNALLY_VALIDATED` | Deterministic rule from above |
| `anti_overclaim_flags` | Guardrail warnings | List of string flags (e.g. `SINGLE_CENTER_NO_EXTERNAL_VALIDATION`, `IN_VITRO_OR_VARIANT_ONLY_NO_HUMAN_COHORT`) | Rule-based flags |

### External-validation evidence rules

The extractor scans the supplied sentences before choosing a label. It recognizes
explicit completed validation/testing statements, not the mere availability of
independent data or a named external cohort. Negated validation does not qualify.
Future or hypothetical testing alone returns `not_stated`; it does not establish
whether some other external validation already occurred.
Statements denying evidence of completed validation also return `not_stated`.

A single-center development statement does not override an explicit completed
external-validation statement elsewhere. Calls for additional validation likewise
do not negate completed testing. Explicit positive and negative validation statements
produce `not_stated`, retaining a verbatim negative evidence sentence for review.
The existing single-sentence field cannot represent both sides of a contradiction.

These are bounded language rules, not a general natural-language inference system:
unrecognized phrasing may be missed, and model/cohort identity across statements is
not resolved. `not_stated` therefore includes unresolved evidence, not just silence.
Review the original text when that label appears. Under the unchanged tier rules,
an unresolved case can still receive T2; that is not confirmation of internal
validation, clinical readiness, or suitability for deployment.

## Epistemic Readiness Tiers

- **`T3_EXTERNALLY_VALIDATED`**: Large human patient cohort ($N \ge 50$), strict patient-level splitting, and independent external validation on a distinct registry or center (e.g., AutoPrognosis validated on the Canadian CF Registry).
- **`T2_INTERNALLY_VALIDATED`**: Human patient cohort ($N \ge 50$) with patient-level splitting and internal temporal or cross-validation, but lacking independent external center validation.
- **`T1_RETROSPECTIVE_OR_UNVALIDATED`**: Retrospective single-center cohorts, unsupervised morphological clustering, sample/record-level splits, or studies that explicitly state external validation is required before clinical use.
- **`T0_PRECLINICAL_BENCHMARK`**: Preclinical studies (in vitro cell lines, organoids, synthetic traces, or variant-level databases) with no human patient cohort.

## Benchmark (`cf_ml_benchmark.json`)

The extractor was evaluated on 6 published Cystic Fibrosis machine learning and computational papers:
1. **Jiang et al. 2026** (*MedComm*, `doi:10.1002/mco2.70894`): Retrospective single-center EHR AMR prediction in 209 pwCF; patient-level splitting; lacks external validation; `AUCs of 0.75-0.80` $\to$ `T1`.
2. **Campredon et al. 2022** (*Eur Respir J*, `doi:10.1183/13993003.01344-2021`): Chest CT radiomics in 283 patients; unsupervised k-means clustering; `-1.40±1.53 points` Bhalla decrease $\to$ `T1`.
3. **McDonald et al. 2024** (*PLoS ONE*, `doi:10.1371/journal.pone.0297560`): AlphaMissense on CFTR2 variants; variant database evaluation; zero patient cohort; fair classification $\to$ `T0`.
4. **Mangkalaphiban et al. 2024** (*Nat Commun*, `doi:10.1038/s41467-024-46703-z`): Extended stop codon machine learning in HEK293T cells; zero human cohort $\to$ `T0`.
5. **Bulcaen et al. 2024** (*Cell Reports Med*, `doi:10.1016/j.xcrm.2024.101560`): Prime editing and DETECTOR ML in rectal organoids; in vitro model; zero human cohort $\to$ `T0`.
6. **Alaa et al. 2018** (*Sci Rep*, `doi:10.1038/s41598-018-29229-5`): AutoPrognosis on 10,261 UK CF Registry patients; patient-level split; externally validated on Canadian CF Registry; `AUC-ROC of 0.89` $\to$ `T3`.

Benchmark Performance: **1.0000 exact-match accuracy** (30/30 fields) and **100% verbatim sentence grounding**.

This is agreement with the six curated fixtures, not an independent clinical
evaluation. Synthetic regression tests additionally exercise negation, prospective
language, cohort mentions, sentence order, and conflicting evidence; they do not
measure extraction accuracy on unseen papers.

## Evaluation Harness (`evaluate_readiness.py`)

`evaluate_readiness.py` separates factual extraction accuracy from heuristic `readiness_tier`
scoring and tracks statistical power explicitly:

- **Factual extraction metrics (per field)**: positive-assertion precision with
  95% Wilson score confidence intervals, recall on stated ground truth, answer coverage,
  abstention rate (`not stated` / `not_stated`), verbatim sentence grounding rate, and an
  error taxonomy (`false_positive_assertion`, `false_negative_miss`, `value_mismatch`,
  `ungrounded_quote`, `grounding_span_mismatch`). Achieving a 95% Wilson lower bound at zero
  errors requires $N \ge 73$ positive assertions; smaller samples are reported explicitly as
  `underpowered_n_lt_73`.
- **Heuristic tier isolation**: `readiness_tier` (`T0`–`T3`) is summarized in a separate
  `heuristic_readiness_tier` section tracking tier match rate, overclaims, and `T3` overclaims.
- **Version-pinned pre-PR #388 comparator**: `--compare-baseline` runs the frozen pre-PR #388
  `extract_external_validation` logic (`d8acab74d94c6d5bfdf1b4c2cbbf6769da4bb830`) on the same
  tokenized sentences to quantify negation/contradiction precision and `T3` overclaim reduction.
- **Synthetic stress fixtures (`stress_eval_synthetic.json`)**: Contains 18 explicitly labeled
  synthetic regression fixtures (`is_unseen_published_papers: false`) covering negation,
  prospective wording, uncertainty denial, section contradictions, multi-count text, in vitro
  models, silent abstracts, prompt injection, and novel metric phrasing. It does **not**
  measure accuracy on unseen published papers or clinical validity.
- **Dual-annotator protocol (`--compare-annotators A.json B.json`)**: Validates `readiness-eval-v1`
  annotation files, checks for `paper_id` / DOI / PMID / `study_group_id` leakage against the
  frozen 6-paper dev fixture, reports raw agreement over both all rows and jointly asserted rows,
  lists disagreements for adjudication, and computes nominal Cohen's kappa on closed categorical
  fields only (flagging `underpowered_n_lt_30` when $N < 30$).

## Usage

```bash
# Run frozen 6-paper fixture benchmark
python tools/readiness/extractor.py --benchmark --json

# Run labeled synthetic stress suite with pre-PR #388 baseline comparator
python tools/readiness/evaluate_readiness.py --stress --compare-baseline --json

# Compare two independent readiness-eval-v1 annotation files for unseen papers
python tools/readiness/evaluate_readiness.py --compare-annotators annotator_a.json annotator_b.json --json

# Parse a file
python tools/readiness/extractor.py --file abstract.txt

# Run unit tests
python -m unittest discover -s tools/readiness -t . -p 'test_*.py' -v
```

## Security & Medical Boundaries

- **Input Sanitization**: Untrusted input is quarantined against prompt injection attempts (`<system>`, `ignore previous instructions`, etc.).
- **Dignity & Medical Boundaries**: `tools/readiness` is an epistemic evaluation tool for scientific literature. It does not provide medical advice, diagnosis, treatment recommendations, or clinical validation.
