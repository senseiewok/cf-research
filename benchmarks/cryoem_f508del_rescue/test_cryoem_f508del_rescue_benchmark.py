#!/usr/bin/env python3
"""Unit tests for Modulator-Free vs Trikafta-Bound deltaF508 Cryo-EM Rescue Benchmark (PRED-0014)."""
from __future__ import annotations

import copy
from fractions import Fraction
from pathlib import Path
import sys
import unittest

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import cryoem_f508del_rescue_benchmark as m1  # noqa: E402
import independent_audit as m2  # noqa: E402


class TestCryoEMF508delRescueBenchmark(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.report_m1 = m1.run_benchmark(run_falsify=True, require_exposed_status=True)
        cls.report_m2 = m2.run_independent_audit()

    def test_01_construct_residue_counts_and_parity(self) -> None:
        pc = self.report_m1["positive_controls"]
        self.assertEqual(pc["placed_ca_8ej1"], 1127)
        self.assertEqual(pc["placed_ca_8eiq"], 1162)
        self.assertEqual(pc["shared_s_common"], 1121)
        self.assertEqual(pc["only_8eiq_count"], 41)
        self.assertEqual(pc["only_8ej1_count"], 6)
        self.assertEqual(self.report_m2["placed_ca_8ej1"], 1127)
        self.assertEqual(self.report_m2["placed_ca_8eiq"], 1162)
        self.assertEqual(self.report_m2["shared_s_common"], 1121)
        self.assertEqual(self.report_m2["noncovalent_pairs_verified"], 626646)

    def test_02_nbd_dimer_closure_and_scaffolding_positive_controls(self) -> None:
        pc = self.report_m1["positive_controls"]
        self.assertGreater(pc["atp_site1_contraction_angstrom"], 15.0)
        self.assertGreater(pc["atp_site2_contraction_angstrom"], 15.0)
        self.assertGreater(pc["min_scaffold_s_stab"], 0.70)
        self.assertLess(pc["max_atp_site_s_stab"], 0.56)
        self.assertGreater(pc["min_scaffold_s_stab"], pc["max_atp_site_s_stab"])
        self.assertGreater(pc["dev_exposed_a455_s_stab"], pc["dev_exposed_l467_s_stab"])
        self.assertGreater(pc["dev_exposed_cond_missense_auroc"], 0.50)

    def test_03_se3_x_se3_invariance_swap_symmetry_and_improper_reflection_rejection(self) -> None:
        pc = self.report_m1["positive_controls"]
        self.assertLess(pc["max_se3_x_se3_drift"], 1e-12)
        self.assertEqual(pc["max_swap_symmetry_diff"], 0.0)
        self.assertEqual(self.report_m2["so3_q_drift"], "0/1")
        improper_float = (
            (2.0 / 3.0, -2.0 / 3.0, 1.0 / 3.0),
            (1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0),
            (2.0 / 3.0, 1.0 / 3.0, -2.0 / 3.0),
        )
        with self.assertRaises(m1.CryoEMF508delRescueError):
            m1.verify_proper_rotation_matrix(improper_float)
        improper_q = (
            (Fraction(2, 3), Fraction(-2, 3), Fraction(1, 3)),
            (Fraction(1, 3), Fraction(2, 3), Fraction(2, 3)),
            (Fraction(2, 3), Fraction(1, 3), Fraction(-2, 3)),
        )
        with self.assertRaises(AssertionError):
            m2.verify_so3_q_matrix(improper_q)

    def test_04_delta508_peptide_bond_exclusion_control(self) -> None:
        chain_d508 = {506: (0.0, 0.0, 0.0), 507: (3.8, 0.0, 0.0), 509: (7.6, 0.0, 0.0), 510: (11.4, 0.0, 0.0)}
        _, deg_excl, _ = m1.compute_lddt_sym_conservation_map(chain_d508, chain_d508, exclude_507_509=True)
        _, deg_leak, _ = m1.compute_lddt_sym_conservation_map(chain_d508, chain_d508, exclude_507_509=False)
        self.assertEqual(deg_excl[507], 1)
        self.assertEqual(deg_excl[509], 1)
        self.assertEqual(deg_leak[507], 2)
        self.assertTrue(m1.is_covalent_deltaf508_pair(507, 509, exclude_507_509=True))
        self.assertFalse(m1.is_covalent_deltaf508_pair(507, 510, exclude_507_509=True))

    def test_05_unmodeled_holdout_missense_abstain_to_none(self) -> None:
        rows = {v["legacy_spl_name"]: v for v in self.report_m1["holdout_trikafta_v23"]["variants"]}
        for unmodeled_name in ("D110H", "D648V", "E822K"):
            self.assertIsNone(rows[unmodeled_name]["s_stab"])
            self.assertIsNone(rows[unmodeled_name]["g_stab"])
        self.assertEqual(rows["3199del6"]["g_stab"], 0.5)
        self.assertEqual(rows["1341+1G->A"]["g_stab"], 0.0)
        self.assertEqual(rows["3120+1G->A"]["g_stab"], 0.0)
        self.assertEqual(rows["3272-26A->G"]["g_stab"], 0.5)

    def test_06_all_seven_falsification_controls_pass(self) -> None:
        fc = self.report_m1["falsification_controls"]
        self.assertEqual(len(fc), 7)
        for name, verdict in fc.items():
            self.assertEqual(verdict, "PASS", msg=name)

    def test_07_manual_catalog_sources_blocked_fail_closed(self) -> None:
        cat = m1.load_catalog_permissions()
        for blocked_id in ("mavedb-dump-zenodo", "cftr-dms-mckee-2023", "cftr-dms-tedman-2025", "nonexistent-source"):
            with self.assertRaises(PermissionError):
                m1.check_manual_catalog_gate(blocked_id, catalog_entries=cat)

    def test_08_spl_provenance_and_exposure_guards(self) -> None:
        dev = m1.json.loads(m1.DEFAULT_DEV_JSON.read_text(encoding="utf-8"))
        hold = m1.json.loads(m1.DEFAULT_HOLDOUT_JSON.read_text(encoding="utf-8"))

        bad_spl = copy.deepcopy(hold)
        bad_spl["variants"][0]["label_status"]["trikafta_v23_spl"] = "listed_on_fda_spl"
        bad_spl["variants"][0]["trikafta_v23_responsive"] = 0
        with self.assertRaises(m1.CryoEMF508delRescueError):
            m1.validate_holdout_annotations(bad_spl, dev, expected_variants_sha256=None)

        # G970R is not listed on KALYDECO v36 Table 3 or SYMDEKO v22 Table 6
        g970r = next(r for r in hold["variants"] if r["variant_id"] == "CFTR-G970R")
        self.assertEqual(g970r["kalydeco_v36_responsive"], 0)
        self.assertEqual(g970r["symdeko_v22_responsive"], 0)

        # Manual catalog source inside provenance.catalog_sources fails closed
        bad_prov = copy.deepcopy(hold)
        bad_prov["provenance"]["catalog_sources"] = ["myvariant", "mavedb-dump-zenodo"]
        with self.assertRaises(PermissionError):
            m1.validate_holdout_annotations(bad_prov, dev)

        # Multi-codon indel endpoint within <=5 aa of reserved control codon 1030 fails closed
        bad_span = copy.deepcopy(hold)
        idx_indel = next(i for i, r in enumerate(bad_span["variants"]) if r["variant_id"] == "CFTR-Ile1023_Val1024del")
        bad_span["variants"][idx_indel]["codon_position"] = 1024
        bad_span["variants"][idx_indel]["hgvs_p_or_c"] = "p.Ile1024_Val1025del"
        with self.assertRaises(m1.CryoEMF508delRescueError):
            m1.validate_holdout_annotations(bad_span, dev, expected_variants_sha256=None)

        # Numeric imputation on ETI-only ECL1 codon 110 fails closed in evaluate_holdout_panel
        s_8ej1, s_8eiq, _ = m1.load_verified_structures()
        c_8ej1 = {k: pt for k, (_, pt) in s_8ej1.ca_by_seq.items()}
        c_8eiq = {k: pt for k, (_, pt) in s_8eiq.ca_by_seq.items()}
        common_set = frozenset(set(c_8ej1) & set(c_8eiq))
        s_map, _, _ = m1.compute_lddt_sym_conservation_map(c_8ej1, c_8eiq)
        with self.assertRaises(m1.CryoEMF508delRescueError):
            m1.evaluate_holdout_panel(hold["variants"], {**s_map, 110: 0.0}, valid_s_common=common_set)

        # Unknown variant_family and bool codon_position fail closed
        with self.assertRaises(m1.CryoEMF508delRescueError):
            m1.score_variant_baselines({"variant_family": "unknown_family", "hgvs_p_or_c": "p.Ser13Phe"})
        bad_bool = copy.deepcopy(hold)
        bad_bool["variants"][0]["codon_position"] = True
        with self.assertRaises(m1.CryoEMF508delRescueError):
            m1.validate_holdout_annotations(bad_bool, dev, expected_variants_sha256=None)

        with self.assertRaises(m1.CryoEMF508delRescueError):
            m1.compute_mann_whitney_auroc([(0.9, True), (0.1, False)], abstention_aware=True)  # type: ignore[list-item]
        with self.assertRaises(m1.CryoEMF508delRescueError):
            m1.compute_mann_whitney_auroc([(float("nan"), 1), (0.1, 0)], abstention_aware=True)

    def test_09_ast_isolation_and_exact_rational_method_parity(self) -> None:
        m2.assert_ast_isolation()
        ev1 = self.report_m1["holdout_trikafta_v23"]
        ev2 = self.report_m2
        self.assertEqual(ev2["placed_missense_auroc_rational"], "7/10")
        self.assertEqual(ev2["missense_all14_abst_auroc_rational"], "61/98")
        self.assertEqual(ev2["standalone_am_auroc_rational"], "22/49")
        self.assertEqual(ev2["all24_g_stab_auroc_rational"], "185/256")
        self.assertEqual(ev2["all24_g_am_auroc_rational"], "185/256")
        self.assertEqual(ev2["all24_class_baseline_auroc_rational"], "95/128")
        self.assertAlmostEqual(ev1["conditional_placed_missense_auroc_s_stab"], float(Fraction(7, 10)), places=6)
        self.assertAlmostEqual(ev1["missense_all14_abstention_auroc_s_stab"], float(Fraction(61, 98)), places=6)
        self.assertAlmostEqual(ev1["standalone_missense_auroc_s_am"], float(Fraction(22, 49)), places=6)
        self.assertAlmostEqual(ev1["gated_all24_abstention_auroc_g_stab"], float(Fraction(185, 256)), places=6)
        self.assertAlmostEqual(ev1["gated_all24_abstention_auroc_g_am"], float(Fraction(185, 256)), places=6)
        self.assertAlmostEqual(ev1["sequence_consequence_baseline_all24_auroc"], float(Fraction(95, 128)), places=6)
        self.assertEqual(ev1["prediction_outcome"], "refuted")


if __name__ == "__main__":
    unittest.main()
