"""Unit tests for tools/discovery/science_model_benchmark.py (T-0265)."""
from __future__ import annotations

import pytest
from science_model_benchmark import (
    CANONICAL_VARIANT_BENCHMARK,
    disentangle_folding_vs_catalytic_defect,
    evaluate_science_model_concordance,
    spearman_rank_correlation,
)


def test_spearman_rank_correlation_exact():
    assert spearman_rank_correlation([1.0, 2.0, 3.0, 4.0], [10.0, 20.0, 30.0, 40.0]) == 1.0
    assert spearman_rank_correlation([1.0, 2.0, 3.0, 4.0], [40.0, 30.0, 20.0, 10.0]) == -1.0
    with pytest.raises(ValueError):
        spearman_rank_correlation([1.0], [1.0])


def test_disentangle_g551d_vs_f508del():
    """ESM-2 + ThermoMPNN separates Class III gating arrest (G551D) from Class II folding defect (F508del)."""
    g551d = disentangle_folding_vs_catalytic_defect(
        CANONICAL_VARIANT_BENCHMARK["G551D"].esm2_650m_llr,
        CANONICAL_VARIANT_BENCHMARK["G551D"].thermompnn_ddg_kcal_mol,
    )
    assert g551d["mechanistic_classification"] == "CLASS_III_CATALYTIC_GATING_ARREST"
    assert g551d["predicted_band_c_maturation_fraction"] > 0.75

    f508del = disentangle_folding_vs_catalytic_defect(
        CANONICAL_VARIANT_BENCHMARK["F508del"].esm2_650m_llr,
        CANONICAL_VARIANT_BENCHMARK["F508del"].thermompnn_ddg_kcal_mol,
    )
    assert f508del["mechanistic_classification"] == "CLASS_II_THERMODYNAMIC_DESTABILIZATION"
    assert f508del["predicted_band_c_maturation_fraction"] < 0.05


def test_evaluate_science_model_concordance_metrics():
    report = evaluate_science_model_concordance()
    assert report["task_id"] == "T-0265"
    metrics = report["concordance_metrics"]
    assert metrics["variant_count"] >= 10
    assert metrics["thermompnn_vs_mechanistic_ddg_spearman_rho"] >= 0.95
    assert metrics["thermompnn_vs_mechanistic_ddg_mae_kcal_mol"] < 0.50
    assert metrics["esm2_llr_vs_po_deficit_spearman_rho"] >= 0.85
