#!/usr/bin/env python3
"""Unit tests for basal_stem_persistence.py (T-0220)."""
import math
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from basal_stem_persistence import (
    BasalStemPersistenceKinetics,
    CURATIVE_THRESHOLD_PCT,
    DEFAULT_MU_D_PER_DAY,
)


class TestBasalStemPersistence(unittest.TestCase):
    def test_homeostatic_steady_state_and_mass_balance(self):
        model = BasalStemPersistenceKinetics(b0=30.0, p_s=0.5)
        self.assertAlmostEqual(model.t_steady, 10.0, places=4)
        self.assertAlmostEqual(model.d_steady, 28.8539, places=4)
        db, dt, dd = model.derivatives(model.b0, model.t_steady, model.d_steady)
        self.assertAlmostEqual(db, 0.0, places=9)
        self.assertAlmostEqual(dt, 0.0, places=9)
        self.assertAlmostEqual(dd, 0.0, places=9)

        traj = model.simulate_hierarchy_rk4(days=365, dt=1.0)
        self.assertEqual(len(traj), 366)
        self.assertAlmostEqual(traj[0].total_density, traj[-1].total_density, places=6)
        self.assertAlmostEqual(traj[-1].basal_density, 30.0, places=6)

    def test_episomal_exponential_decay(self):
        model = BasalStemPersistenceKinetics()
        c0 = 0.40
        # Pure luminal episomal vector decays with half-life 100 days
        val_0 = model.episomal_retention(0.0, c0_luminal=c0, c_basal0=0.0)
        val_100 = model.episomal_retention(100.0, c0_luminal=c0, c_basal0=0.0)
        val_1000 = model.episomal_retention(1000.0, c0_luminal=c0, c_basal0=0.0)
        self.assertAlmostEqual(val_0, c0, places=6)
        self.assertAlmostEqual(val_100, 0.5 * c0, places=5)
        self.assertLess(val_1000, 0.001)

    def test_integrative_asymptotic_persistence(self):
        model = BasalStemPersistenceKinetics(b0=30.0)
        # 20% basal stem cell editing (b_edited = 6.0), 0% initial luminal editing
        f_0 = model.integrative_correction(0.0, b_edited=6.0, f_diff0=0.0)
        f_100 = model.integrative_correction(100.0, b_edited=6.0, f_diff0=0.0)
        f_3650 = model.integrative_correction(3650.0, b_edited=6.0, f_diff0=0.0)
        self.assertAlmostEqual(f_0, 0.0, places=6)
        self.assertAlmostEqual(f_100, 0.10, places=4)  # Halfway to 0.20 at t_1/2 = 100 d
        self.assertAlmostEqual(f_3650, 0.20, places=5)  # Asymptotically approaches B_edited / B_0 = 0.20

    def test_redosing_interval_calculation(self):
        model = BasalStemPersistenceKinetics()
        # With c0 = 0.40 and t_1/2 = 100 d, dropping to 0.10 (two half-lives) takes 200 days
        interval = model.calculate_redosing_interval_days(c0_luminal=0.40, c_basal0=0.0, threshold_fraction=0.10)
        self.assertIsNotNone(interval)
        self.assertAlmostEqual(interval, 200.0, delta=0.1)

    def test_clinical_modalities_evaluation(self):
        model = BasalStemPersistenceKinetics()
        summaries = model.evaluate_clinical_modalities()
        self.assertEqual(len(summaries), 5)
        by_id = {s.modality_id: s for s in summaries}
        self.assertEqual(by_id["MOD-MRNA-LNP"].durability_verdict, "REQUIRES_PERIODIC_REDOSING")
        self.assertLess(by_id["MOD-MRNA-LNP"].redosing_interval_days, 15.0)
        self.assertEqual(by_id["MOD-AAV-EPISOMAL"].durability_verdict, "REQUIRES_PERIODIC_REDOSING")
        self.assertGreater(by_id["MOD-AAV-EPISOMAL"].redosing_interval_days, 150.0)
        self.assertIn("PMID: 42529376", by_id["MOD-EV-BASAL-PE"].clinical_ref)
        self.assertEqual(by_id["MOD-BASAL-ENGRAFT"].durability_verdict, "DURABLE_HOMEOSTATIC_CURE")
        self.assertGreaterEqual(by_id["MOD-BASAL-ENGRAFT"].year_10_cftr_pct, CURATIVE_THRESHOLD_PCT)

    def test_invalid_parameters_raise(self):
        with self.assertRaises(ValueError):
            BasalStemPersistenceKinetics(b0=0.0)
        with self.assertRaises(ValueError):
            BasalStemPersistenceKinetics(p_s=1.5)
        model = BasalStemPersistenceKinetics()
        with self.assertRaises(ValueError):
            model.integrative_correction(day=10.0, b_edited=50.0)


if __name__ == "__main__":
    unittest.main()
