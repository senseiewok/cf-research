"""Unit test suite for Inhaled CFTR Synthetic mRNA Translation Velocity Model (T-0272)."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

from mrna_translation_velocity import (
    CFTR_CDS_NUCLEOTIDES,
    CFTR_CODONS,
    CFTR_DOMAINS,
    DEFAULT_POLY_A_INITIAL_LENGTH,
    NLM_ATTRIBUTION,
    RESEARCH_BOUNDARY,
    PolyATailSpec,
    TranslationConfig,
    build_codon_elongation_profile,
    calculate_cap_factor,
    calculate_poly_a_synergy,
    evaluate_synthetic_mrna_design,
    get_cftr_domains,
    simulate_mrna_lifecycle,
    simulate_poly_a_length,
    simulate_tasep_ribosome_dynamics,
)


def test_cftr_sequence_and_domain_partition():
    """Verify CFTR coding sequence dimensions and contiguous domain partition."""
    assert CFTR_CDS_NUCLEOTIDES == 4443
    assert CFTR_CODONS == 1480
    # 1,480 sense amino acid codons (4,440 nt) plus 3 nt stop codon = 4,443 nt CDS
    assert CFTR_CDS_NUCLEOTIDES == CFTR_CODONS * 3 + 3

    domains = get_cftr_domains()
    assert len(domains) == 5

    # Check contiguous partition from codon 1 to 1480
    expected_start = 1
    total_aa = 0
    for domain in domains:
        assert domain.start_codon == expected_start
        assert domain.end_codon >= domain.start_codon
        assert domain.length_aa == domain.end_codon - domain.start_codon + 1
        total_aa += domain.length_aa
        expected_start = domain.end_codon + 1

    assert total_aa == 1480
    assert domains[-1].end_codon == 1480

    # Key mutation locations
    nbd1 = next(d for d in domains if d.name == "NBD1")
    assert nbd1.start_codon <= 508 <= nbd1.end_codon  # F508
    assert nbd1.start_codon <= 542 <= nbd1.end_codon  # G542X
    assert nbd1.start_codon <= 551 <= nbd1.end_codon  # G551D

    nbd2 = next(d for d in domains if d.name == "NBD2")
    assert nbd2.start_codon <= 1282 <= nbd2.end_codon  # W1282X
    assert nbd2.start_codon <= 1303 <= nbd2.end_codon  # N1303K


def test_codon_elongation_profiles():
    """Verify codon velocity profiles across optimization modes."""
    # Balanced folding mode
    balanced = build_codon_elongation_profile("balanced_folding")
    assert len(balanced) == CFTR_CODONS
    assert all(1.5 <= v <= 6.0 for v in balanced)
    # Domain boundary pause verification (e.g. MSD1->NBD1 transition at codon 385)
    assert balanced[384] == 2.0  # In pause window
    assert balanced[100] == 5.5  # Core globular region

    # Maximal speed mode
    fast = build_codon_elongation_profile("maximal_speed")
    assert len(fast) == CFTR_CODONS
    assert all(v == 5.5 for v in fast)

    # Wild-type mode
    wt = build_codon_elongation_profile("wild_type")
    assert len(wt) == CFTR_CODONS
    assert sum(wt) / len(wt) < sum(balanced) / len(balanced)

    with pytest.raises(ValueError):
        build_codon_elongation_profile("invalid_mode")


def test_poly_a_synergy_and_cap_factors():
    """Verify initiation synergy scaling with poly(A) length and cap structure."""
    assert calculate_poly_a_synergy(0.0) == 0.0
    assert calculate_poly_a_synergy(-10.0) == 0.0

    syn_30 = calculate_poly_a_synergy(30.0)
    syn_80 = calculate_poly_a_synergy(80.0)
    syn_110 = calculate_poly_a_synergy(110.0)

    assert 0.0 < syn_30 < syn_80 < syn_110 <= 1.0
    assert syn_110 > 0.90  # Optimal tail length achieves high synergy

    # Cap factors
    assert calculate_cap_factor("Cap1") == 1.00
    assert calculate_cap_factor("Cap0") == 0.65
    assert calculate_cap_factor("Cap2") == 1.08

    with pytest.raises(ValueError):
        calculate_cap_factor("CapInvalid")


def test_poly_a_deadenylation_kinetics():
    """Verify two-phase deadenylation kinetics from protected to unprotected state."""
    l0 = 110.0
    l_0h = simulate_poly_a_length(l0, 0.0)
    l_10h = simulate_poly_a_length(l0, 10.0)
    l_30h = simulate_poly_a_length(l0, 30.0)
    l_60h = simulate_poly_a_length(l0, 60.0)

    assert l_0h == 110.0
    assert l_0h > l_10h > l_30h > l_60h >= 0.0

    # Protected phase deadenylation: 10 hours = 600 min * 0.025 = 15 nt loss
    assert abs(l_10h - (110.0 - 15.0)) < 0.1


def test_tasep_ribosome_steady_state():
    """Verify TASEP steady-state occupancy, exclusion, and ribosome load."""
    rates = build_codon_elongation_profile("balanced_folding")
    profile = simulate_tasep_ribosome_dynamics(rates, initiation_rate=0.08)

    assert profile.codons_total == 1480
    assert 15.0 <= profile.ribosome_load <= 35.0
    assert 0.04 <= profile.steady_state_flux_codons_s <= 0.08
    assert all(0.0 <= rho <= 1.0 for rho in profile.steady_state_density)
    # Transit time: 1480 codons / ~5 codons/s ~ 300 s (~5 min)
    assert 200.0 <= profile.mean_transit_time_seconds <= 400.0
    assert profile.ngd_safety_score >= 85.0
    assert profile.disome_collision_count == 0


def test_ngd_collision_detection_under_extreme_load():
    """Pathologically elevated initiation rate induces disome collision and lowers NGD score."""
    rates = [5.5] * CFTR_CODONS
    # Artificially create a severe bottleneck at codon 500
    for i in range(495, 505):
        rates[i] = 0.2

    # High initiation pushing ribosomes into bottleneck
    profile_stress = simulate_tasep_ribosome_dynamics(rates, initiation_rate=0.60)
    assert profile_stress.disome_collision_count > 0
    assert profile_stress.ngd_safety_score < 80.0
    assert len(profile_stress.collision_hotspots) > 0


def test_mrna_lifecycle_nucleoside_modification():
    """100% m1Psi modification significantly enhances half-life and channel yield."""
    cfg_unmod = TranslationConfig(nucleoside_modification="unmodified")
    cfg_m1psi = TranslationConfig(nucleoside_modification="m1Psi")

    res_unmod = simulate_mrna_lifecycle(cfg_unmod, total_hours=48.0)
    res_m1psi = simulate_mrna_lifecycle(cfg_m1psi, total_hours=48.0)

    assert res_unmod.mrna_half_life_hours == 4.5
    assert res_m1psi.mrna_half_life_hours == 18.0
    # Yield should be substantially higher with m1Psi due to 4x longer half-life
    assert res_m1psi.total_cftr_yield_molecules_per_mrna > res_unmod.total_cftr_yield_molecules_per_mrna * 2.5


def test_full_lifecycle_simulation_and_verdict():
    """Verify nominal translation lifecycle result, folding fidelity, and verdict."""
    config = TranslationConfig(
        codon_optimization="balanced_folding",
        nucleoside_modification="m1Psi",
        cap_structure="Cap1",
    )
    res = simulate_mrna_lifecycle(config, total_hours=72.0)

    assert res.ngd_decay_risk_verdict == "OPTIMAL_COLLISION_FREE"
    assert res.co_translational_folding_fidelity == 0.91
    assert res.total_cftr_yield_molecules_per_mrna > 1000.0
    assert len(res.time_points) > 50

    # Compare against maximal speed (where folding fidelity is compromised)
    cfg_fast = TranslationConfig(codon_optimization="maximal_speed")
    res_fast = simulate_mrna_lifecycle(cfg_fast, total_hours=72.0)
    assert res.co_translational_folding_fidelity > res_fast.co_translational_folding_fidelity


def test_clinical_trial_surveillance_and_boundaries():
    """Verify trial surveillance anchors, NLM attribution, and absence of prohibited phrases."""
    evaluation = evaluate_synthetic_mrna_design()

    assert "NCT05668741" in evaluation["trial_surveillance_grounding"]
    assert "NCT03375047" in evaluation["trial_surveillance_grounding"]
    assert evaluation["trial_surveillance_grounding"]["nlm_attribution"] == NLM_ATTRIBUTION

    raw_json = json.dumps(evaluation)
    prohibited = [
        "proven to cure",
        "guarantees clinical benefit",
        "recommended treatment regimen",
        "safe and effective for patients",
        "cures cystic fibrosis",
    ]
    for phrase in prohibited:
        assert phrase not in raw_json.lower()


def test_cli_execution_and_verify(tmp_path):
    """Test CLI subprocess execution with --verify and --json flags."""
    script = Path(__file__).resolve().parent / "mrna_translation_velocity.py"

    # --verify flag
    res_verify = subprocess.run(
        [sys.executable, str(script), "--verify"],
        capture_output=True,
        text=True,
        check=False,
    )
    assert res_verify.returncode == 0
    assert "PASS [T-0272 self-verification]" in res_verify.stdout

    # --json flag
    res_json = subprocess.run(
        [sys.executable, str(script), "--json", "--hours", "12.0"],
        capture_output=True,
        text=True,
        check=False,
    )
    assert res_json.returncode == 0
    data = json.loads(res_json.stdout)
    assert "tasep_profile" in data
    assert "trial_surveillance_grounding" in data
