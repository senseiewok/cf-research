"""Unit tests for Boltzmann Spliceosome-ASO Steric-Blocking Kinetics Model (T-0254)."""
from __future__ import annotations

import sys
from pathlib import Path
import pytest

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from aso_splice_kinetics import (
    ALLELE_KINETICS_SPECS,
    CHEMISTRY_KINETICS,
    compute_boltzmann_splice_partition,
    resolve_allele_spec,
    simulate_allele_dose_response,
)


def test_allele_catalog_grounding():
    """All three canonical splicing alleles must carry verified NCT, PDB, and DOI provenance."""
    assert set(ALLELE_KINETICS_SPECS) == {"3849+10kbC->T", "2789+5G->A", "1811+1.6kbA->G"}
    for allele_id, spec in ALLELE_KINETICS_SPECS.items():
        assert spec.clinical_trial_nct in {"NCT06429176", "NCT04684641"}
        assert "6FF4" in spec.structural_pdbs
        assert "7ABI" in spec.structural_pdbs
        assert len(spec.primary_dois) >= 2


def test_boltzmann_partition_conservation_and_monotonicity():
    """Canonical + aberrant splice fractions must sum to 1.0 and increase monotonically with ASO dose."""
    for allele_id in ALLELE_KINETICS_SPECS:
        prev_canon = -1.0
        for dose in (0.0, 5.0, 25.0, 50.0, 100.0, 250.0, 500.0):
            res = compute_boltzmann_splice_partition(allele_id, dose, chemistry="2MOE")
            assert abs((res.canonical_splice_fraction + res.aberrant_splice_fraction) - 1.0) < 1e-3
            assert 0.0 <= res.aso_target_occupancy <= 1.0
            assert res.canonical_splice_fraction >= prev_canon
            prev_canon = res.canonical_splice_fraction


def test_3849_spl84_therapeutic_rescue():
    """3849+10kbC->T must be sub-therapeutic at 0 nM and exceed 40% WT current at >=100 nM 2'-O-MOE."""
    base = compute_boltzmann_splice_partition("3849+10kbC->T", 0.0, chemistry="2MOE")
    assert base.canonical_splice_fraction < 0.20
    assert base.curative_threshold_met is False

    treated = compute_boltzmann_splice_partition("3849+10kbC->T", 100.0, chemistry="2MOE")
    assert treated.canonical_splice_fraction > 0.70
    assert treated.restored_isc_pct_wt >= 40.0
    assert treated.curative_threshold_met is True


def test_2moe_exceeds_pmo_potency_at_sub_saturating_doses():
    """Higher-affinity 2'-O-MOE-PS must outperform uncharged PMO at 50 nM across all three alleles."""
    for allele_id in ALLELE_KINETICS_SPECS:
        moe = compute_boltzmann_splice_partition(allele_id, 50.0, chemistry="2MOE")
        pmo = compute_boltzmann_splice_partition(allele_id, 50.0, chemistry="PMO")
        assert moe.aso_target_occupancy > pmo.aso_target_occupancy
        assert moe.canonical_splice_fraction > pmo.canonical_splice_fraction
        assert moe.restored_isc_pct_wt > pmo.restored_isc_pct_wt


def test_eti_modulator_synergy_on_2789():
    """Co-administration of ETI modulators must amplify restored CFTR current for 2789+5G->A."""
    aso_only = compute_boltzmann_splice_partition("2789+5G->A", 80.0, chemistry="2MOE", eti_co_administered=False)
    combo = compute_boltzmann_splice_partition("2789+5G->A", 80.0, chemistry="2MOE", eti_co_administered=True)
    assert combo.restored_isc_pct_wt > aso_only.restored_isc_pct_wt


def test_dose_response_simulation_and_ec50():
    """Dose-response curve must return 8 points and a bounded apparent EC50."""
    sim = simulate_allele_dose_response("1811+1.6kbA->G", chemistry="2MOE")
    assert len(sim["dose_response"]) == 8
    assert 5.0 <= sim["apparent_ec50_nm"] <= 150.0


def test_invalid_inputs_rejected():
    """Negative ASO concentrations and unknown alleles/chemistries must raise ValueError."""
    with pytest.raises(ValueError):
        compute_boltzmann_splice_partition("3849+10kbC->T", -5.0)
    with pytest.raises(ValueError):
        compute_boltzmann_splice_partition("UNKNOWN_ALLELE", 50.0)
    with pytest.raises(ValueError):
        compute_boltzmann_splice_partition("3849+10kbC->T", 50.0, chemistry="LNA_UNKNOWN")
    assert resolve_allele_spec("3849").allele_id == "3849+10kbC->T"
