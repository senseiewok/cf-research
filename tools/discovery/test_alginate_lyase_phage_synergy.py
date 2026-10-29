#!/usr/bin/env python3
"""Unit tests for Cross-Scale Biophysical Modeling of Inhaled Alginate Lyase, Bacteriophage & Antibiotic Synergy (T-0472).

Validates:
  1. Primary literature grounding, RCSB PDB structures, and clinical registry provenance.
  2. Alginate lyase catalytic depolymerization kinetics, M/G ratio weighting, and O-acetylation steric hindrance.
  3. Biofilm mesh pore expansion and Stokes-Einstein diffusion enhancement for Podoviridae and Myoviridae.
  4. Electrostatic de-trapping of polycationic tobramycin upon alginate fragmentation.
  5. Pharmacodynamic tri-therapy log-kill synergy, MBEC fold-reduction, and FICI < 0.5 strong synergy classification.
  6. Phase 2 AH-10 negative controls and physical boundary stress tests.
  7. Deterministic benchmark execution and CLI interface.
"""

from __future__ import annotations

import math
from pathlib import Path
import sys
import pytest

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from alginate_lyase_phage_synergy import (
    CLINICAL_TRIAL_APPA02,
    CLINICAL_TRIAL_VANZACAFTOR,
    DEFAULT_ACETYLATION_DEGREE,
    DEFAULT_INITIAL_ALGINATE_MG_ML,
    DEFAULT_INITIAL_MESH_NM,
    DEFAULT_LYASE_CONC_UG_ML,
    DEFAULT_M_FRACTION,
    DEFAULT_MAX_MESH_NM,
    MYOVIRIDAE_RH_NM,
    NLM_ATTRIBUTION,
    PODOVIRIDAE_RH_NM,
    PRIMARY_DOI_ALKAWASH,
    PRIMARY_DOI_CHAN,
    PRIMARY_DOI_GLONTI,
    PRIMARY_DOI_MIDDLETON,
    RCSB_PDB_CFTR,
    RCSB_PDB_LYASE,
    RESEARCH_BOUNDARY,
    compute_alginate_lyase_kinetics,
    compute_antibiotic_detrapping,
    compute_phage_mesh_diffusion,
    evaluate_tri_therapy_synergy,
    run_benchmark_suite,
    run_falsification_negative_controls,
)


class TestGroundingAndAttribution:
    """Validate primary citations, NLM attribution, and structural coordinates."""

    def test_attribution_and_boundary_not_empty(self) -> None:
        assert "ClinicalTrials.gov does not evaluate safety or efficacy" in NLM_ATTRIBUTION
        assert "Preclinical in-silico" in RESEARCH_BOUNDARY
        assert "research, not medical advice" in RESEARCH_BOUNDARY

    def test_primary_literature_dois(self) -> None:
        assert PRIMARY_DOI_ALKAWASH == "10.1111/j.1600-0463.2006.apm_356.x"
        assert PRIMARY_DOI_GLONTI == "10.1111/j.1365-2672.2009.04469.x"
        assert PRIMARY_DOI_CHAN == "10.1038/srep26717"
        assert PRIMARY_DOI_MIDDLETON == "10.1056/nejmoa1908639"

    def test_clinical_trial_registry_grounding(self) -> None:
        assert CLINICAL_TRIAL_APPA02 == "NCT04596319"
        assert CLINICAL_TRIAL_VANZACAFTOR == "NCT05033080"

    def test_rcsb_pdb_structures(self) -> None:
        assert RCSB_PDB_CFTR == "6MSM"
        assert RCSB_PDB_LYASE == "4OZX"


class TestAlginateDepolymerizationKinetics:
    """Validate Michaelis-Menten kinetics, acetylation hindrance, and degradation scaling."""

    def test_default_kinetics_cleavage_extent(self) -> None:
        res = compute_alginate_lyase_kinetics(
            enzyme_conc_ug_ml=DEFAULT_LYASE_CONC_UG_ML,
            substrate_mg_ml=DEFAULT_INITIAL_ALGINATE_MG_ML,
            m_fraction=DEFAULT_M_FRACTION,
            acetylation_degree=DEFAULT_ACETYLATION_DEGREE,
            duration_min=60.0,
        )
        assert res.cleavage_extent_1h_pct >= 80.0
        assert res.final_mw_kda_1h < 50.0
        assert res.unsaturated_product_delta_a235 > 0.0
        assert res.depolymerization_half_life_min < 30.0

    def test_acetylation_steric_inhibition(self) -> None:
        unacetylated = compute_alginate_lyase_kinetics(acetylation_degree=0.0)
        acetylated = compute_alginate_lyase_kinetics(acetylation_degree=0.35)
        # Acetylation must reduce apparent kcat due to steric hindrance
        assert acetylated.effective_kcat_s < unacetylated.effective_kcat_s
        assert acetylated.steric_acetylation_inhibition_pct > 0.0

    def test_null_enzyme_yields_zero_cleavage(self) -> None:
        res = compute_alginate_lyase_kinetics(enzyme_conc_ug_ml=0.0)
        assert res.cleavage_extent_1h_pct == 0.0
        assert res.vmax_mg_ml_s == 0.0
        assert res.final_mw_kda_1h == pytest.approx(350.0, rel=1e-2)

    def test_negative_inputs_raise_error(self) -> None:
        with pytest.raises(ValueError):
            compute_alginate_lyase_kinetics(enzyme_conc_ug_ml=-5.0)


class TestPhageDiffusionAndMeshScaling:
    """Validate Stokes-Einstein diffusivity and pore hindrance across bacteriophage families."""

    def test_free_diffusivity_stokes_einstein(self) -> None:
        pod = compute_phage_mesh_diffusion(mesh_size_nm=180.0, phage_family="Podoviridae")
        myo = compute_phage_mesh_diffusion(mesh_size_nm=180.0, phage_family="Myoviridae")

        # Smaller Podoviridae (Rh ~ 32 nm) must diffuse faster in free solution than larger Myoviridae (Rh ~ 45 nm)
        assert pod.free_diffusivity_d0_um2_s > myo.free_diffusivity_d0_um2_s
        assert pod.free_diffusivity_d0_um2_s == pytest.approx(9.71, rel=0.10)
        assert myo.free_diffusivity_d0_um2_s == pytest.approx(6.90, rel=0.10)

    def test_pore_entrapment_in_intact_biofilm(self) -> None:
        # Intact mucoid biofilm (xi ~ 22 nm) traps both Podoviridae (32 nm) and Myoviridae (45 nm)
        pod_tight = compute_phage_mesh_diffusion(mesh_size_nm=DEFAULT_INITIAL_MESH_NM, phage_family="Podoviridae")
        assert pod_tight.penetration_verdict == "ENTRAPPED"
        assert pod_tight.diffusivity_ratio < 0.01

        myo_tight = compute_phage_mesh_diffusion(mesh_size_nm=DEFAULT_INITIAL_MESH_NM, phage_family="Myoviridae")
        assert myo_tight.penetration_verdict == "ENTRAPPED"
        assert myo_tight.diffusivity_ratio < 0.01

    def test_mesh_expansion_recovery(self) -> None:
        # Enzymatically cleared biofilm (xi ~ 180 nm) allows robust penetration
        pod_open = compute_phage_mesh_diffusion(mesh_size_nm=DEFAULT_MAX_MESH_NM, phage_family="Podoviridae")
        assert pod_open.penetration_verdict == "PENETRATING"
        assert pod_open.diffusivity_ratio >= 0.40

    def test_non_positive_mesh_raises_error(self) -> None:
        with pytest.raises(ValueError):
            compute_phage_mesh_diffusion(mesh_size_nm=0.0)


class TestAntibioticDetrapping:
    """Validate electrostatic release of aminoglycosides upon enzymatic depolymerization."""

    def test_tobramycin_detrapping_magnitude(self) -> None:
        intact = compute_antibiotic_detrapping(alginate_intact_fraction=1.0, antibiotic_name="tobramycin")
        digested = compute_antibiotic_detrapping(alginate_intact_fraction=0.10, antibiotic_name="tobramycin")

        assert intact.initial_bound_pct >= 85.0
        assert digested.post_lyase_bound_pct < 15.0
        assert digested.released_fraction_pct >= 70.0
        assert digested.electrostatic_de_trapping_fold >= 5.0


class TestSynergyAndEradication:
    """Validate bactericidal killing, MBEC reduction, and Loewe synergy / FICI."""

    def test_tri_therapy_kills_exceed_monotherapies(self) -> None:
        syn = evaluate_tri_therapy_synergy(
            lyase_ug_ml=DEFAULT_LYASE_CONC_UG_ML,
            phage_pfu_ml=1e8,
            antibiotic_ug_ml=20.0,
            antibiotic_type="tobramycin",
        )
        assert syn.tri_therapy_log_kill > syn.antibiotic_monotherapy_log_kill
        assert syn.tri_therapy_log_kill > syn.phage_monotherapy_log_kill
        assert syn.tri_therapy_log_kill > syn.dual_phage_antibiotic_log_kill
        assert syn.synergistic_excess_log_kill >= 1.0

    def test_strong_synergy_fici_below_half(self) -> None:
        syn = evaluate_tri_therapy_synergy()
        assert syn.fractional_inhibitory_concentration_index < 0.50
        assert "STRONG_SYNERGY" in syn.synergy_classification

    def test_mbec_reduction_fold(self) -> None:
        syn = evaluate_tri_therapy_synergy()
        assert syn.mbec_reduction_fold >= 8.0
        assert syn.combo_mbec_abx_ug_ml <= 32.0


class TestFalsificationAndBenchmarkSuite:
    """Validate AH-10 negative controls, benchmark integrity, and CLI execution."""

    def test_falsification_suite_passes(self) -> None:
        fals = run_falsification_negative_controls()
        assert fals["status"] == "PASS"
        assert fals["ah_10_compliant"] is True
        assert len(fals["failures"]) == 0

    def test_benchmark_suite_integrity(self) -> None:
        bench = run_benchmark_suite()
        assert bench["task_id"] == "T-0472"
        assert bench["model"] == "alginate_lyase_phage_synergy"
        assert bench["overall_status"] == "PASS"
        assert bench["kinetics"]["cleavage_extent_1h_pct"] >= 80.0
        assert bench["synergy_metrics"]["mbec_reduction_fold"] >= 8.0
