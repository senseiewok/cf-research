"""Unit tests for Round 106: Adenine Base Editor (ABE) Deaminase Window Bystander-vs-Target Kinetic Selectivity."""

from tools.discovery.round106_base_editor_bystander_markov_selectivity import (
    CITATIONS,
    HYPOTHESIS,
    ROUND_ID,
    TITLE,
    build_round106_summary,
    compute_theoretical_clean_ceiling,
    run_falsification,
    run_falsification_checks,
)


def test_round106_metadata_and_ah10_contract() -> None:
    assert ROUND_ID == "round106"
    assert "Adenine Base Editor" in TITLE
    assert "S_kin" in HYPOTHESIS
    assert len(CITATIONS) >= 4
    assert all(c.get("scite_verified") is True for c in CITATIONS)
    ah10 = run_falsification()
    assert ah10["rule"] == "AH-10"
    assert ah10["passed"] is True


def test_round106_falsification_checks_all_pass() -> None:
    res = run_falsification_checks()
    assert res["passed"] is True
    for name, ok in res["checks"].items():
        assert ok is True, f"Failed check: {name}"


def test_round106_markov_clean_ceiling_and_w1282x_comparison() -> None:
    summary = build_round106_summary()
    by_id = summary["by_locus"]

    r553x = by_id["r553x_narrow_abe"]
    w1282x_nar = by_id["w1282x_context_engineered_abe"]
    w1282x_wide = by_id["w1282x_wide_window_abe8e"]

    assert r553x["observed_clean_rescue_prob"] > 0.70
    assert w1282x_nar["observed_clean_rescue_prob"] > 0.50
    assert w1282x_wide["collateral_missense_prob"] > 0.70
    assert w1282x_wide["observed_clean_rescue_prob"] < 0.15
    assert compute_theoretical_clean_ceiling(10.0) > compute_theoretical_clean_ceiling(1.0)
