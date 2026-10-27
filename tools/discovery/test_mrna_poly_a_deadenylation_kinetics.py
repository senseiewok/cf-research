#!/usr/bin/env python3
"""Unit test suite for Inhaled CFTR mRNA Poly(A) Tail Deadenylation Kinetics & Ribosome Collision Engine (T-0449)."""

from __future__ import annotations

import sys
from pathlib import Path
import pytest

# Ensure tools/discovery is on path
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DISCOVERY_DIR = REPO_ROOT / "tools" / "discovery"
if str(DISCOVERY_DIR) not in sys.path:
    sys.path.insert(0, str(DISCOVERY_DIR))

from mrna_poly_a_deadenylation_kinetics import (
    MRNAPolyADeadenylationEngine,
    DeadenylationState,
    TASEPCollisionProfile,
    ChannelYieldResult,
    DOI_POLYA_PROFILING,
    DOI_CODON_OPTIMALITY,
    DOI_PSEUDOURIDINE_MRNA,
    NCT_VX522_MRNA_TRIAL,
    PDB_INITIATION_COMPLEX,
    PDB_ACTIVE_CFTR,
    CURATIVE_RESCUE_CHANNELS,
    CCR4_NOT_TRANSITION_THRESHOLD_NT,
    main,
)


@pytest.fixture
def engine() -> MRNAPolyADeadenylationEngine:
    return MRNAPolyADeadenylationEngine()


def test_metadata_and_primary_grounding(engine: MRNAPolyADeadenylationEngine) -> None:
    suite = engine.run_comprehensive_biophysical_suite()
    meta = suite["metadata"]

    assert meta["task"] == "T-0449"
    assert meta["track"] == "Track A"
    assert meta["discipline"] == "rna_gene_editing"

    # Primary literature grounding
    assert DOI_POLYA_PROFILING in meta["primary_literature"]
    assert DOI_CODON_OPTIMALITY in meta["primary_literature"]
    assert DOI_PSEUDOURIDINE_MRNA in meta["primary_literature"]

    # Clinical trial registry and structural grounding
    assert meta["clinical_registry"] == NCT_VX522_MRNA_TRIAL
    assert PDB_INITIATION_COMPLEX in meta["structural_grounding"]
    assert PDB_ACTIVE_CFTR in meta["structural_grounding"]


def test_biphasic_deadenylation_dynamics(engine: MRNAPolyADeadenylationEngine) -> None:
    # Length > 25 nt: PAN2/PAN3 protected slow phase
    rate_slow, phase_slow = engine.compute_deadenylation_rate(100.0, codon_optimized=False, m1psi_modified=False)
    assert phase_slow == "PAN2_PAN3_PROTECTED"
    assert rate_slow > 0.0

    # Length <= 25 nt: CCR4-NOT rapid terminal decay phase
    rate_fast, phase_fast = engine.compute_deadenylation_rate(15.0, codon_optimized=False, m1psi_modified=False)
    assert phase_fast == "CCR4_NOT_RAPID_DECAY"
    assert rate_fast > rate_slow * 3.5  # rapid decay is significantly faster than slow phase

    # Chemical modification and codon optimization stabilization
    rate_opt, _ = engine.compute_deadenylation_rate(100.0, codon_optimized=True, m1psi_modified=True)
    assert rate_opt < rate_slow  # stabilized construct decays much slower


def test_circularization_and_initiation_scaling(engine: MRNAPolyADeadenylationEngine) -> None:
    # Full-length poly(A) tail (120 nt)
    circ_120, k_init_120, pabp_120 = engine.compute_circularization_and_initiation(120.0)
    assert circ_120 > 0.90
    assert 0.04 <= k_init_120 <= 0.12  # within biophysical specification (0.04–0.12 Hz)
    assert pabp_120 >= 4  # 120 // 27 = 4 PABP monomers

    # Intermediate poly(A) tail (50 nt)
    circ_50, k_init_50, pabp_50 = engine.compute_circularization_and_initiation(50.0)
    assert circ_50 < circ_120
    assert k_init_50 < k_init_120
    assert pabp_50 == 1

    # Short poly(A) tail below threshold (10 nt)
    circ_10, k_init_10, pabp_10 = engine.compute_circularization_and_initiation(10.0)
    assert circ_10 < 0.10
    assert k_init_10 < 0.03
    assert pabp_10 == 0


def test_tasep_ribosome_collision_and_spacing(engine: MRNAPolyADeadenylationEngine) -> None:
    k_init = 0.080  # Hz
    tasep_wt = engine.evaluate_tasep_ribosome_dynamics(k_init, codon_optimized=False)
    tasep_opt = engine.evaluate_tasep_ribosome_dynamics(k_init, codon_optimized=True)

    # Elongation velocities
    assert tasep_opt.mean_elongation_velocity_codons_per_s > tasep_wt.mean_elongation_velocity_codons_per_s

    # Ribosome spacing should exceed 10-codon footprint
    assert tasep_wt.mean_ribosome_spacing_codons > 10.0
    assert tasep_opt.mean_ribosome_spacing_codons > tasep_wt.mean_ribosome_spacing_codons

    # Disome collision rate should be reduced in codon-optimized construct
    assert tasep_opt.disome_collision_rate_per_hour < tasep_wt.disome_collision_rate_per_hour
    assert tasep_opt.no_go_decay_probability_per_transcript < tasep_wt.no_go_decay_probability_per_transcript


def test_functional_half_life_and_trajectory(engine: MRNAPolyADeadenylationEngine) -> None:
    traj = engine.simulate_deadenylation_trajectory(
        initial_length_nt=120,
        codon_optimized=True,
        m1psi_modified=True,
    )
    assert len(traj) > 0
    assert traj[0].poly_a_length_nt == 120.0
    assert traj[-1].poly_a_length_nt <= CCR4_NOT_TRANSITION_THRESHOLD_NT

    # Half life calculation
    t_half_opt = engine.calculate_functional_half_life(120, codon_optimized=True, m1psi_modified=True)
    t_half_unmod = engine.calculate_functional_half_life(100, codon_optimized=False, m1psi_modified=False)

    assert t_half_opt > t_half_unmod * 2.0  # stabilized construct has >2x functional half-life


def test_channel_yield_and_curative_rescue(engine: MRNAPolyADeadenylationEngine) -> None:
    # Severely degraded tail (15 nt, 1 copy) should be sub-curative (< 150 channels/cell)
    yield_low = engine.simulate_channel_yield(
        initial_poly_a_length_nt=15,
        cytosolic_mrna_copies=1,
        codon_optimized=False,
        m1psi_modified=False,
    )
    assert not yield_low.curative_threshold_achieved
    assert yield_low.steady_state_surface_channels_per_cell < CURATIVE_RESCUE_CHANNELS

    # Therapeutic dose (15 copies, optimized 120 nt poly(A)) should achieve curative threshold (>= 150 channels/cell)
    yield_high = engine.simulate_channel_yield(
        initial_poly_a_length_nt=120,
        cytosolic_mrna_copies=15,
        codon_optimized=True,
        m1psi_modified=True,
    )
    assert yield_high.curative_threshold_achieved
    assert yield_high.steady_state_surface_channels_per_cell >= CURATIVE_RESCUE_CHANNELS
    assert yield_high.percent_wild_type_density >= 10.0


def test_negative_control_and_boundary_conditions(engine: MRNAPolyADeadenylationEngine) -> None:
    # 0 nt poly(A) length
    rate_zero, phase_zero = engine.compute_deadenylation_rate(0.0)
    assert rate_zero == 0.0
    assert phase_zero == "TERMINATED"

    circ_zero, k_init_zero, pabp_zero = engine.compute_circularization_and_initiation(0.0)
    assert circ_zero == 0.0
    assert pabp_zero == 0
    assert k_init_zero < 0.010

    # 0 cytosolic mRNA copies yields zero channels
    yield_zero = engine.simulate_channel_yield(
        initial_poly_a_length_nt=120,
        cytosolic_mrna_copies=0,
    )
    assert yield_zero.steady_state_surface_channels_per_cell == 0.0
    assert not yield_zero.curative_threshold_achieved


def test_cli_execution() -> None:
    assert main(["--run-all"]) == 0
    assert main(["--json", "--poly-a-len", "120", "--mrna-copies", "15"]) == 0
