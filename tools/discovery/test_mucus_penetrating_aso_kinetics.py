#!/usr/bin/env python3
"""Unit tests for Mucus-Penetrating PEGylated Peptide-ASO Splice & NMD Kinetics (T-0497)."""

from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import mucus_penetrating_aso_kinetics as mak


class TestMucusPenetratingAsoKinetics(unittest.TestCase):
    """Test suite verifying hydrodynamic sieving, cellular delivery, and functional CFTR rescue."""

    def test_grounding_identifiers_and_disclaimers(self):
        """Verify primary DOIs, clinical registry NCTs, PDB codes, and mandatory disclaimers."""
        self.assertEqual(mak.CLINICAL_TRIAL_SPL84, "NCT06429176")
        self.assertEqual(mak.CLINICAL_TRIAL_VX522, "NCT05668741")
        self.assertEqual(mak.CLINICAL_TRIAL_VANZACAFTOR, "NCT05033080")
        self.assertEqual(mak.PDB_DEPHOSPHO_CFTR, "6MSM")
        self.assertEqual(mak.PDB_TRIKAFTA_CFTR, "8EIQ")

        self.assertEqual(mak.DOI_HANES_MUCUS_PENETRATION, "10.1016/j.addr.2008.11.002")
        self.assertEqual(mak.DOI_SUK_PEGYLATION, "10.1016/j.addr.2015.09.012")
        self.assertEqual(mak.DOI_BUTTON_AIRWAY_MUCUS, "10.1152/physrev.00004.2021")
        self.assertEqual(mak.DOI_KEREM_SPL84_DEVELOPMENT, "10.1016/j.jcf.2021.06.003")
        self.assertEqual(mak.DOI_MANFREDI_W1282X_BARRIERS, "10.1021/acs.biochem.5c00590")

        self.assertIn("ClinicalTrials.gov does not evaluate safety or efficacy", mak.NLM_ATTRIBUTION)
        self.assertIn("research, not medical advice", mak.RESEARCH_BOUNDARY)
        self.assertEqual(mak.SYMPTOM_RELIEF_THRESHOLD_PCT_WT, 10.0)
        self.assertEqual(mak.CLINICAL_RESPONDER_THRESHOLD_PCT_WT, 15.0)
        self.assertEqual(mak.CURATIVE_BENCHMARK_PCT_WT, 50.0)

    def test_stokes_einstein_free_diffusivity(self):
        """Verify free Stokes-Einstein diffusivity scaling with particle size."""
        d0_small = mak.compute_stokes_einstein_d0(2.0)
        d0_med = mak.compute_stokes_einstein_d0(4.0)
        d0_large = mak.compute_stokes_einstein_d0(8.0)

        self.assertGreater(d0_small, d0_med)
        self.assertGreater(d0_med, d0_large)
        # Ratio for doubling radius should be exactly 0.5
        self.assertAlmostEqual(d0_med / d0_small, 0.5, places=2)
        self.assertGreater(d0_small, 50.0)  # > 50 um^2/s

    def test_ogston_renkin_hydrodynamic_sieving(self):
        """Verify steric pore and fiber obstruction factors."""
        prof_3nm = mak.compute_ogston_renkin_sieving(r_h_nm=3.0, surface_zeta_mv=-1.0, peg_mol_pct=3.0, xi_mesh_nm=100.0)
        prof_10nm = mak.compute_ogston_renkin_sieving(r_h_nm=10.0, surface_zeta_mv=-1.0, peg_mol_pct=3.0, xi_mesh_nm=100.0)

        self.assertGreater(prof_3nm.ogston_sieving_factor, prof_10nm.ogston_sieving_factor)
        self.assertGreater(prof_3nm.renkin_pore_factor, prof_10nm.renkin_pore_factor)
        self.assertGreater(prof_3nm.steric_diffusivity_ratio, prof_10nm.steric_diffusivity_ratio)
        self.assertLessEqual(prof_3nm.steric_diffusivity_ratio, 1.0)
        self.assertGreater(prof_3nm.steric_diffusivity_ratio, 0.0)

    def test_pegylation_and_zeta_potential_non_adhesive_boundary(self):
        """Verify non-adhesive transport for dense PEG coatings and entrapment of unshielded charges."""
        # Non-adhesive PEGylated conjugate
        stealth = mak.compute_ogston_renkin_sieving(r_h_nm=3.8, surface_zeta_mv=-1.5, peg_mol_pct=3.5)
        self.assertEqual(stealth.adhesion_retention_factor, 1.0)
        self.assertEqual(stealth.verdict, "MUCUS_PENETRATING")
        self.assertGreater(stealth.mucus_diffusion_ratio, 0.40)
        self.assertGreater(stealth.mucus_penetration_efficiency_pct, 60.0)

        # Unshielded cationic peptide conjugate
        cationic = mak.compute_ogston_renkin_sieving(r_h_nm=3.8, surface_zeta_mv=+18.0, peg_mol_pct=0.0)
        self.assertLess(cationic.adhesion_retention_factor, 0.05)
        self.assertEqual(cationic.verdict, "ENTRAPPED")
        self.assertLess(cationic.mucus_diffusion_ratio, 0.05)
        self.assertLess(cationic.mucus_penetration_efficiency_pct, 10.0)

    def test_mucus_transit_and_clearance_competition(self):
        """Verify transit time and clearance competition through varying mucus thicknesses."""
        thin_mucus = mak.compute_ogston_renkin_sieving(r_h_nm=3.8, surface_zeta_mv=-1.0, peg_mol_pct=3.0, mucus_thickness_um=10.0)
        thick_mucus = mak.compute_ogston_renkin_sieving(r_h_nm=3.8, surface_zeta_mv=-1.0, peg_mol_pct=3.0, mucus_thickness_um=30.0)

        self.assertLess(thin_mucus.transit_time_seconds, thick_mucus.transit_time_seconds)
        self.assertGreater(thin_mucus.mucus_penetration_efficiency_pct, thick_mucus.mucus_penetration_efficiency_pct)

    def test_cellular_uptake_and_endosomal_escape(self):
        """Verify cell-penetrating peptide enhancement of nuclear concentration."""
        bare = mak.compute_cellular_nuclear_delivery(apical_dose_nm=100.0, conjugate_type="bare_aso")
        peg_pep = mak.compute_cellular_nuclear_delivery(apical_dose_nm=100.0, conjugate_type="peg_peptide_aso")

        self.assertGreater(peg_pep.nuclear_concentration_nm, bare.nuclear_concentration_nm * 5.0)
        self.assertGreater(peg_pep.endosomal_escape_pct, 20.0)
        self.assertLess(bare.endosomal_escape_pct, 5.0)

    def test_splice_rescue_3849_dose_response(self):
        """Verify dose-dependent splicing correction for deep-intronic 3849+10kbC->T."""
        doses = [0.0, 10.0, 38.0, 100.0, 300.0]
        results = [mak.compute_splice_nmd_kinetics("3849+10kbC->T", d) for d in doses]

        # Target occupancy should increase monotonically
        for i in range(len(results) - 1):
            self.assertLess(results[i].target_occupancy, results[i + 1].target_occupancy)
            self.assertLess(results[i].restored_isc_pct_wt, results[i + 1].restored_isc_pct_wt)

        # High dose achieves clinical responder threshold
        self.assertTrue(results[-1].clinical_responder_met)
        self.assertGreater(results[-1].restored_isc_pct_wt, 15.0)

    def test_nonsense_rescue_w1282x_and_g542x(self):
        """Verify nonsense mutation rescue with Manfredi et al. (2025) ERAD chaperone and NMD grounding."""
        w1282x_low = mak.compute_splice_nmd_kinetics("W1282X", nuclear_aso_nm=5.0)
        w1282x_high = mak.compute_splice_nmd_kinetics("W1282X", nuclear_aso_nm=80.0, corrector_boost=1.8, potentiator_boost=2.0)

        self.assertGreater(w1282x_high.restored_isc_pct_wt, w1282x_low.restored_isc_pct_wt)
        self.assertTrue(w1282x_high.symptom_relief_met)

        g542x_low = mak.compute_splice_nmd_kinetics("G542X", nuclear_aso_nm=5.0)
        g542x_high = mak.compute_splice_nmd_kinetics("G542X", nuclear_aso_nm=80.0, corrector_boost=1.8, potentiator_boost=2.0)

        self.assertGreater(g542x_high.restored_isc_pct_wt, g542x_low.restored_isc_pct_wt)
        # G542X baseline NMD escape is lower than W1282X
        self.assertLess(g542x_low.canonical_or_rescued_mrna_pct_wt, w1282x_low.canonical_or_rescued_mrna_pct_wt)

    def test_clinical_threshold_crossings(self):
        """Verify boolean clinical milestone flags."""
        low_dose_sim = mak.simulate_mucus_penetrating_aso(variant="W1282X", apical_dose_nm=10.0, conjugate_type="bare_aso")
        high_dose_sim = mak.simulate_mucus_penetrating_aso(variant="W1282X", apical_dose_nm=200.0, conjugate_type="peg_peptide_aso")

        self.assertFalse(low_dose_sim.rescue.symptom_relief_met)
        self.assertTrue(high_dose_sim.rescue.symptom_relief_met)
        self.assertTrue(high_dose_sim.rescue.clinical_responder_met)

    def test_ah10_deterministic_negative_controls(self):
        """Verify AH-10 negative control suite passes."""
        res = mak.execute_negative_control_assertion()
        self.assertEqual(res["rule"], "AH-10")
        self.assertEqual(res["verdict"], "PASS")
        self.assertTrue(res["passed"])
        self.assertTrue(res["control_1_entrapped"])
        self.assertTrue(res["control_2_zero_dose_subtherapeutic"])

    def test_epistemic_rigor_parameters(self):
        """Verify parameter provenance under Rule AH-9."""
        params = mak.CORE_BIOPHYSICAL_PARAMETERS
        self.assertGreaterEqual(len(params), 10)
        for key, p in params.items():
            self.assertIsNotNone(p.doi_or_pdb)
            self.assertIsNotNone(p.location_in_source)
            self.assertIn(p.category, ("empirical", "derived", "phenomenological_heuristic"))
            self.assertIsInstance(p.uncertainty_range, tuple)
            self.assertEqual(len(p.uncertainty_range), 2)
            self.assertLessEqual(p.uncertainty_range[0], p.uncertainty_range[1])

    def test_cli_execution_and_json_serialization(self):
        """Verify CLI driver operates cleanly and produces valid JSON."""
        ret_check = mak.main(["--check", "--json"])
        self.assertEqual(ret_check, 0)

        ret_sim = mak.main(["--variant", "W1282X", "--dose-nm", "150.0", "--json"])
        self.assertEqual(ret_sim, 0)


if __name__ == "__main__":
    unittest.main()
