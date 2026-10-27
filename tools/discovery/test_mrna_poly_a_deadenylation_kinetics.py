#!/usr/bin/env python3
"""Unit tests for Biophysical Modeling of Inhaled CFTR mRNA Poly(A) Tail Deadenylation Kinetics and Ribosomal Translocation Collision Protection (T-0449)."""

from __future__ import annotations

import json
from pathlib import Path
import sys
import pytest

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
for p in (str(HERE), str(REPO_ROOT)):
    if p not in sys.path:
        sys.path.insert(0, p)

from mrna_poly_a_deadenylation_kinetics import (
    CFTR_CDS_LENGTH_NT,
    CRITICAL_POLYA_THRESHOLD_NT,
    CURATIVE_CHANNEL_YIELD_THRESHOLD,
    DOI_CODON_OPTIMALITY,
    DOI_POLYA_PROFILING,
    DOI_PSEUDOURIDINE_MRNA,
    NCT_VX522_MRNA_TRIAL,
    NLM_ATTRIBUTION,
    OPTIMAL_POLYA_TAIL_NT,
    PABP1_FOOTPRINT_NT,
    PDB_ACTIVE_CFTR,
    PDB_CCR4_NOT_COMPLEX,
    RESEARCH_BOUNDARY,
    MRNAPolyADeadenylationEngine,
    main,
)


@pytest.fixture
def engine() -> MRNAPolyADeadenylationEngine:
    return MRNAPolyADeadenylationEngine()


def test_grounding_identifiers_and_disclaimers() -> None:
    """Verify primary literature DOIs, NCT registry IDs, PDB structures, and mandatory disclaimers."""
    assert NCT_VX522_MRNA_TRIAL == "NCT05668741"
    assert PDB_CCR4_NOT_COMPLEX == "6GQ1"
    assert PDB_ACTIVE_CFTR == "5UAK"
    assert DOI_POLYA_PROFILING == "10.1038/nature13007"
    assert DOI_CODON_OPTIMALITY == "10.1016/j.cell.2015.02.029"
    assert DOI_PSEUDOURIDINE_MRNA == "10.1038/mt.2008.200"
    assert CFTR_CDS_LENGTH_NT == 4443
    assert PABP1_FOOTPRINT_NT == 27
    assert CRITICAL_POLYA_THRESHOLD_NT == 25.0
    assert OPTIMAL_POLYA_TAIL_NT == 110.0
    assert "ClinicalTrials.gov does not evaluate safety or efficacy" in NLM_ATTRIBUTION
    assert "research, not medical advice" in RESEARCH_BOUNDARY


def test_biphasic_deadenylation_pan2_vs_ccr4_not(engine: MRNAPolyADeadenylationEngine) -> None:
    """Domain 1: Verify slow PAN2-PAN3 deadenylation above 25 nt and rapid CCR4-NOT decay below 25 nt."""
    rate_110, phase_110 = engine.calculate_biphasic_deadenylation_rate(
        polya_length_nt=110.0,
        codon_optimality_index=0.88,
        m1psi_modified=True,
        segmented_linker=True,
    )
    rate_18, phase_18 = engine.calculate_biphasic_deadenylation_rate(
        polya_length_nt=18.0,
        codon_optimality_index=0.88,
        m1psi_modified=True,
        segmented_linker=True,
    )
    assert phase_110 == "PAN2_PAN3_SLOW_PROTECTED"
    assert phase_18 == "CCR4_NOT_RAPID_TERMINAL"
    assert rate_110 < 2.5
    assert rate_18 > rate_110 * 5.0

    # Segmented poly(A) linker reduces PAN2-PAN3 deadenylation rate
    rate_unsegmented, _ = engine.calculate_biphasic_deadenylation_rate(
        polya_length_nt=110.0,
        codon_optimality_index=0.88,
        m1psi_modified=True,
        segmented_linker=False,
    )
    assert rate_110 < rate_unsegmented


def test_pabp_occupancy_and_closed_loop_circularization(engine: MRNAPolyADeadenylationEngine) -> None:
    """Domain 2: Verify PABP1 monomer stoichiometry and eIF4F/PABP closed-loop probability."""
    assert engine.calculate_pabp_occupancy(110.0) == 4  # 4 * 27 = 108 nt
    assert engine.calculate_pabp_occupancy(54.0) == 2
    assert engine.calculate_pabp_occupancy(10.0) == 0

    p_circ_110 = engine.calculate_closed_loop_probability(110.0, m1psi_modified=True)
    p_circ_20 = engine.calculate_closed_loop_probability(20.0, m1psi_modified=True)
    p_circ_0 = engine.calculate_closed_loop_probability(0.0, m1psi_modified=True)
    assert p_circ_110 >= 0.90
    assert p_circ_20 < 0.30
    assert p_circ_0 <= 0.05


def test_tasep_ribosome_translocation_and_ngd_collision(engine: MRNAPolyADeadenylationEngine) -> None:
    """Domain 3: Verify 1D TASEP spacing across 4,443 nt CDS and No-Go Decay (NGD) collision onset."""
    optimal_traffic = engine.evaluate_tasep_ribosome_traffic(
        initiation_rate_hz=0.08,
        codon_optimality_index=0.88,
        m1psi_modified=True,
    )
    stalled_traffic = engine.evaluate_tasep_ribosome_traffic(
        initiation_rate_hz=0.12,
        codon_optimality_index=0.50,
        m1psi_modified=False,
    )
    assert optimal_traffic.traffic_regime in ("FREE_FLOW_OPTIMAL", "HIGH_DENSITY_COOPERATIVE")
    assert optimal_traffic.collision_probability < 0.15
    assert optimal_traffic.mean_inter_ribosome_spacing_nt > 120.0

    assert stalled_traffic.traffic_regime == "QUEUING_COLLISION_NGD_RISK"
    assert stalled_traffic.collision_probability > 0.70
    assert stalled_traffic.ngd_cleavage_rate_hr > optimal_traffic.ngd_cleavage_rate_hr * 5.0


def test_24h_trajectory_and_construct_comparison(engine: MRNAPolyADeadenylationEngine) -> None:
    """Domain 4: Verify 24-hour CFTR channel yield across poly(A) tail length and m1Psi modification."""
    comp = engine.compare_therapeutic_constructs()
    opt = comp["optimized_construct"]
    short = comp["short_tail_control"]
    unmod = comp["unmodified_uridine_control"]

    assert opt["curative_yield_achieved"] is True
    assert opt["total_24h_cftr_channels"] >= CURATIVE_CHANNEL_YIELD_THRESHOLD
    assert opt["functional_half_life_hr"] > short["functional_half_life_hr"]
    assert opt["time_to_critical_threshold_hr"] > short["time_to_critical_threshold_hr"]
    assert comp["yield_fold_gain_vs_unmodified"] >= 2.0
    assert unmod["curative_yield_achieved"] is False


def test_validation_and_cli_execution(capsys: pytest.CaptureFixture) -> None:
    """Verify input validation errors and CLI flags (--check, --json, default)."""
    with pytest.raises(ValueError, match="CDS length must be positive"):
        MRNAPolyADeadenylationEngine(cds_length_nt=0)

    eng = MRNAPolyADeadenylationEngine()
    with pytest.raises(ValueError, match="Poly\\(A\\) tail length cannot be negative"):
        eng.calculate_pabp_occupancy(-5.0)
    with pytest.raises(ValueError, match="Codon optimality index must be in"):
        eng.calculate_biphasic_deadenylation_rate(100.0, codon_optimality_index=1.5)

    assert main(["--check"]) == 0
    out_check = capsys.readouterr().out
    assert "Self-check PASS" in out_check

    assert main(["--json", "--polya", "115.0"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["grounding"]["trials"][0] == "NCT05668741"
    assert "6GQ1" in payload["grounding"]["structures"]

    assert main([]) == 0
    out_std = capsys.readouterr().out
    assert "Sensei Ewok CF Lab - CFTR mRNA Poly(A) Deadenylation & TASEP Model (T-0449)" in out_std
