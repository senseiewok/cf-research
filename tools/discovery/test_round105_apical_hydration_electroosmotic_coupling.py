"""Unit tests for Round 105: Apical ASL Electro-Osmotic Hydration & Donnan Swelling Equilibrium."""

from tools.discovery.round105_apical_hydration_electroosmotic_coupling import (
    CITATIONS,
    HYPOTHESIS,
    ROUND_ID,
    TITLE,
    build_round105_summary,
    run_falsification,
    run_falsification_checks,
)


def test_round105_metadata_and_ah10_contract() -> None:
    assert ROUND_ID == "round105"
    assert "Donnan" in TITLE
    assert "7.0 um" in HYPOTHESIS
    assert len(CITATIONS) >= 4
    assert all(c.get("scite_verified") is True for c in CITATIONS)
    ah10 = run_falsification()
    assert ah10["rule"] == "AH-10"
    assert ah10["passed"] is True


def test_round105_falsification_checks_all_pass() -> None:
    res = run_falsification_checks()
    assert res["passed"] is True
    for name, ok in res["checks"].items():
        assert ok is True, f"Failed check: {name}"


def test_round105_donnan_swelling_and_pcl_clearance_ordering() -> None:
    summary = build_round105_summary()
    by_id = summary["by_genotype"]

    wt = by_id["wt_cftr"]
    eti = by_id["f508del_eti_rescued"]
    g551d = by_id["g551d_ivacaftor"]
    uncorr = by_id["f508del_uncorrected"]
    w1282x = by_id["w1282x_class_i_null"]

    assert wt["periciliary_height_um"] > eti["periciliary_height_um"] >= 7.0
    assert g551d["periciliary_height_um"] >= 7.0
    assert uncorr["periciliary_height_um"] < 3.5
    assert w1282x["periciliary_height_um"] < 3.0
    assert wt["clearance_restored"] is True
    assert eti["clearance_restored"] is True
    assert uncorr["clearance_restored"] is False
    assert w1282x["clearance_restored"] is False
