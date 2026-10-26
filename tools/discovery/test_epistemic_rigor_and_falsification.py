"""Unit tests for Epistemic Rigor (AH-9), Dynamic Falsification (AH-10), and Parametric Sensitivity (AH-11)."""
from __future__ import annotations

import pytest
from biophysical_models import calculate_mucus_penetration_index
from discovery_council import (
    DISCOVERY_MODULE_REGISTRY,
    check_discovery_index_drift,
    execute_discovery_council_round,
    run_dynamic_falsification_suite,
)
from epistemic_rigor import (
    CORE_BIOPHYSICAL_PARAMETERS,
    BiophysicalParameter,
    run_parametric_sensitivity,
    validate_physical_bounds,
)


def test_ah9_biophysical_parameter_provenance():
    """Verify AH-9 parameter provenance tagging and rejection of unanchored empirical parameters."""
    assert len(CORE_BIOPHYSICAL_PARAMETERS) >= 10
    categories_seen = {p.category for p in CORE_BIOPHYSICAL_PARAMETERS.values()}
    assert categories_seen == {"empirical", "derived", "phenomenological_heuristic"}

    for name, param in CORE_BIOPHYSICAL_PARAMETERS.items():
        assert param.name == name
        assert param.units
        if param.category in ("empirical", "derived"):
            assert len(param.doi_or_pdb) >= 4
            assert len(param.location_in_source) >= 4

    # Missing source citation on empirical parameter must raise ValueError
    with pytest.raises(ValueError, match="AH-9 violation"):
        BiophysicalParameter(
            name="UNANCHORED_EMPIRICAL",
            value=1.0,
            units="nM",
            category="empirical",
            doi_or_pdb="",
            location_in_source="Table 1",
        )

    # Value outside uncertainty_range must raise ValueError
    with pytest.raises(ValueError, match="uncertainty_range"):
        BiophysicalParameter(
            name="OUT_OF_RANGE",
            value=50.0,
            units="nM",
            category="phenomenological_heuristic",
            doi_or_pdb="internal",
            location_in_source="fit",
            uncertainty_range=(10.0, 20.0),
        )


def test_physical_boundary_validator():
    """Verify physical and thermodynamic boundary invariants."""
    assert validate_physical_bounds("p_open", 0.42, kind="probability") == 0.42
    assert validate_physical_bounds("wt_isc", 18.92, kind="percentage") == 18.92
    assert validate_physical_bounds("cassette", 3630.0, kind="aav_single_capsid_bp") == 3630.0

    with pytest.raises(ValueError, match="probability"):
        validate_physical_bounds("p_open_bad", 1.25, kind="probability")

    with pytest.raises(ValueError, match="4,700 bp"):
        validate_physical_bounds("oversized_aav", 5200.0, kind="aav_single_capsid_bp")

    with pytest.raises(ValueError, match="Kelvin"):
        validate_physical_bounds("temp_bad", 0.0, kind="positive_kelvin")


def test_ah10_dynamic_falsification_across_all_66_modules():
    """Verify Phase 2 falsification, non-empty mitigations, and negative controls across all 66 modules."""
    suite = run_dynamic_falsification_suite()
    assert suite["status"] == "PASS", f"Failed modules: {suite['failed_modules']}"
    assert suite["total_modules_tested"] == len(DISCOVERY_MODULE_REGISTRY)
    assert suite["total_modules_tested"] >= 66
    assert suite["passed_count"] == suite["total_modules_tested"]

    # Spot-check both thematic topics and individual registered modules
    for topic in ("aerosol_sort_lnp", "nonsense_trimodal", "alpha_evolve_pareto_regimen", "nbd1_thermodynamics"):
        res = execute_discovery_council_round(topic)
        p2 = res["phase_2_falsification"]
        assert p2["falsification_status"] == "STRESS_TESTED_AND_RESOLVED"
        assert len(p2["challenges_raised"]) >= 1
        assert len(p2["mitigations_resolved"]) == len(p2["challenges_raised"])
        assert p2["negative_control"]["passed"] is True


def test_ah11_parametric_sensitivity_25pct_sweep():
    """Verify AH-11 +/-25% perturbation on phenomenological heuristics preserves candidate rank order."""
    params = list(CORE_BIOPHYSICAL_PARAMETERS.values())

    def eval_lnp_candidates(p_map):
        alpha_elec = p_map["MUCUS_ELECTROSTATIC_PENALTY_COEFF"]
        weight_ratio = p_map["COMPOSITE_SAFETY_EFFICACY_WEIGHT_RATIO"]
        pore_nm = p_map["CF_SPUTUM_MUCIN_MESH_PORE_NM"]

        candidates = {
            "SORT_LNP_65nm_Neutral": (65.0, -1.2, 0.92),
            "Standard_LNP_95nm_MildAnionic": (95.0, -4.5, 0.78),
            "Unshielded_LNP_135nm_Cationic": (135.0, +8.0, 0.40),
        }
        scores = {}
        for cid, (diam_nm, zeta_mv, safety_score) in candidates.items():
            base = calculate_mucus_penetration_index(diam_nm, zeta_mv)
            elec_factor = 1.0 / (1.0 + alpha_elec * (zeta_mv ** 2))
            scores[cid] = (
                weight_ratio * base["effective_diffusivity_ratio"] * elec_factor + 0.25 * safety_score
            )
        return scores

    report = run_parametric_sensitivity(eval_lnp_candidates, params, perturbation_pct=0.25, only_heuristics=True)
    assert report["rule"] == "AH-11"
    assert report["parameters_tested"] == 3
    assert report["top_candidate_preserved"] is True
    assert report["full_ranking_preserved"] is True
    assert report["baseline_ranking"][0] == "SORT_LNP_65nm_Neutral"
    assert report["status"] == "PASS"


def test_cross_repo_discovery_index_drift_check():
    """Verify catalog.json, README.md, and biophysical-index.md are synchronized."""
    drift = check_discovery_index_drift()
    assert drift["status"] == "PASS", f"Drift check failed: {drift}"
    assert drift["catalog_json_in_sync"] is True
    assert drift["readme_md_in_sync"] is True
    assert drift["biophysical_index_in_sync"] is True
