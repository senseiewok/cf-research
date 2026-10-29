#!/usr/bin/env python3
"""Comprehensive unit test suite for split_intein_cftr_reconstitution.py (T-0393)."""
from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import unittest

# Ensure tools/discovery is importable
HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from split_intein_cftr_reconstitution import (
    CFTR_MATURE_MW_KDA,
    CFTR_WT_RESIDUES,
    NLM_REGISTRY_DISCLAIMER,
    RESEARCH_BOUNDARY,
    THERAPEUTIC_CL_THRESHOLD_PCT_WT,
    build_full_reconstitution_report,
    compute_cotransduction_poisson,
    compute_intein_thermodynamics,
    map_split_junction_to_cryoem,
    self_test,
    simulate_trans_splicing_kinetics,
    simulate_vector_repeat_dosing,
)


class TestCoTransductionPoissonKinetics(unittest.TestCase):
    def test_equimolar_poisson_distribution(self):
        kinetics = compute_cotransduction_poisson(moi_n=3.0, moi_c=3.0, entry_efficiency=0.65)
        assert kinetics.optimal_ratio is True
        assert kinetics.percent_dual_positive > 70.0
        # Probabilities sum to 1.0
        total_p = (
            kinetics.p_uninfected
            + kinetics.p_n_only
            + kinetics.p_c_only
            + kinetics.p_dual_cotransduction
        )
        assert math.isclose(total_p, 1.0, abs_tol=1e-4)

    def test_zero_moi_boundary(self):
        kinetics = compute_cotransduction_poisson(moi_n=0.0, moi_c=0.0)
        assert kinetics.percent_dual_positive == 0.0
        assert math.isclose(kinetics.p_uninfected, 1.0, abs_tol=1e-4)
        assert kinetics.p_dual_cotransduction == 0.0

    def test_unbalanced_vector_ratio(self):
        # 1:5 ratio should flag as non-optimal
        kinetics = compute_cotransduction_poisson(moi_n=1.0, moi_c=5.0)
        assert kinetics.optimal_ratio is False
        assert kinetics.p_dual_cotransduction < compute_cotransduction_poisson(3.0, 3.0).p_dual_cotransduction

    def test_negative_moi_raises_value_error(self):
        with self.assertRaisesRegex(ValueError, "must be non-negative"):
            compute_cotransduction_poisson(-1.0, 3.0)


class TestInteinThermodynamics(unittest.TestCase):
    def test_cfa_intein_thermodynamics(self):
        thermo = compute_intein_thermodynamics("Cfa")
        assert thermo.intein_key == "Cfa"
        assert thermo.favorable_association is True
        assert thermo.delta_g_assoc_kj_mol < -25.0
        assert thermo.equilibrium_ka_m_inv > 1e4
        assert thermo.equilibrium_kd_nm < 200.0
        assert thermo.t_half_splice_seconds < 60.0

    def test_npu_dnae_thermodynamics(self):
        thermo = compute_intein_thermodynamics("Npu_DnaE")
        assert thermo.intein_key == "Npu_DnaE"
        assert thermo.favorable_association is True
        assert thermo.delta_g_assoc_kj_mol < -20.0
        assert thermo.t_half_splice_seconds > 100.0

    def test_gp41_1_fast_intein_thermodynamics(self):
        thermo = compute_intein_thermodynamics("gp41-1")
        assert thermo.intein_key == "gp41-1"
        assert thermo.favorable_association is True
        assert thermo.delta_g_assoc_kj_mol < -28.0
        assert thermo.t_half_splice_seconds < 35.0

    def test_unknown_intein_raises_key_error(self):
        with self.assertRaisesRegex(KeyError, "Unknown intein"):
            compute_intein_thermodynamics("UnknownInteinXYZ")


class TestCryoEMDomainMapping(unittest.TestCase):
    def test_ser744_in_disordered_r_domain(self):
        m = map_split_junction_to_cryoem(744, "Cfa")
        assert m["split_residue"] == 744
        assert m["assigned_domain"] == "R_DOMAIN"
        assert m["in_disordered_r_domain"] is True
        assert m["structural_perturbation_risk"] == "low"
        assert m["cryoem_pdb_active"] == "6MSM"
        assert m["cryoem_pdb_dephospho"] == "5UAK"
        assert m["both_fit_aav_limit"] is True
        assert m["vector_n_genome_nt"] <= 4700
        assert m["vector_c_genome_nt"] <= 4700

    def test_cys837_split_site(self):
        m = map_split_junction_to_cryoem(837, "Npu_DnaE")
        assert m["split_residue"] == 837
        assert m["assigned_domain"] == "R_DOMAIN"
        assert m["intein_junction_compatible"] is True
        assert m["both_fit_aav_limit"] is True

    def test_invalid_residue_raises_value_error(self):
        with self.assertRaisesRegex(ValueError, "split_residue must be between"):
            map_split_junction_to_cryoem(0)
        with self.assertRaisesRegex(ValueError, "split_residue must be between"):
            map_split_junction_to_cryoem(CFTR_WT_RESIDUES + 10)


class TestReconstitutionYieldAndODE(unittest.TestCase):
    def test_cfa_trans_splicing_yield(self):
        cotrans = compute_cotransduction_poisson(3.0, 3.0)
        thermo = compute_intein_thermodynamics("Cfa")
        recon = simulate_trans_splicing_kinetics(cotrans, thermo, split_site=744)

        assert recon.mature_cftr_yield_pct_wt > 20.0
        assert recon.restored_chloride_conductance_pct_wt > THERAPEUTIC_CL_THRESHOLD_PCT_WT
        assert recon.reaches_therapeutic_threshold is True
        assert recon.mature_mw_kda == CFTR_MATURE_MW_KDA

    def test_zero_transduction_yield(self):
        cotrans = compute_cotransduction_poisson(0.0, 0.0)
        thermo = compute_intein_thermodynamics("Cfa")
        recon = simulate_trans_splicing_kinetics(cotrans, thermo, split_site=744)

        assert recon.mature_cftr_yield_pct_wt == 0.0
        assert recon.restored_chloride_conductance_pct_wt == 0.0
        assert recon.reaches_therapeutic_threshold is False


class TestRepeatDosingHumoralAttenuation(unittest.TestCase):
    def test_redosing_nab_attenuation_trajectory(self):
        report = simulate_vector_repeat_dosing(
            moi_base=3.0,
            dosing_interval_days=90.0,
            total_doses=5,
            baseline_nab_titer=1.0,
            nab_boost_per_dose=80.0,
            nab_half_life_days=65.0,
        )
        assert len(report.doses) == 5
        # Dose 1 has lowest NAb, Dose 5 has accumulated NAb
        assert report.doses[0].accumulated_nab_titer < report.doses[-1].accumulated_nab_titer
        # Attenuation of entry and chloride transport
        assert report.doses[0].restored_cl_transport_pct_wt > report.doses[-1].restored_cl_transport_pct_wt
        assert report.redosing_attenuation_ratio < 1.0
        assert report.max_doses_above_threshold >= 1
        assert report.nct_id == "NCT05248230"

    def test_invalid_doses_raises_value_error(self):
        with self.assertRaisesRegex(ValueError, "total_doses must be >= 1"):
            simulate_vector_repeat_dosing(total_doses=0)


class TestFullReportAndProvenance(unittest.TestCase):
    def test_full_report_assembly(self):
        report = build_full_reconstitution_report()
        assert "Split-Intein Dual-AAV" in report.title
        assert report.nlm_attribution == NLM_REGISTRY_DISCLAIMER
        assert report.research_boundary == RESEARCH_BOUNDARY

        # Verify all primary citations present
        citation_ids = [c["id"] for c in report.primary_citations]
        assert "10.1126/scitranslmed.aav4523" in citation_ids
        assert "10.1126/scitranslmed.adw8886" in citation_ids
        assert "NCT05248230" in citation_ids
        assert "PDB:6MSM" in citation_ids
        assert "PDB:5UAK" in citation_ids

    def test_self_test_function(self):
        assert self_test() == 0

    def test_cli_execution_and_json(self):
        repo_root = HERE.parents[1]
        cmd = [sys.executable, str(HERE / "split_intein_cftr_reconstitution.py"), "--json"]
        res = subprocess.run(cmd, capture_output=True, text=True, check=True, cwd=str(repo_root))
        data = json.loads(res.stdout)
        assert "co_transduction" in data
        assert "intein_thermodynamics" in data
        assert "reconstitution_yield" in data
        assert "redosing_kinetics" in data
        assert "cryoem_domain_mapping" in data
        assert data["nlm_attribution"] == NLM_REGISTRY_DISCLAIMER


if __name__ == "__main__":
    unittest.main()
