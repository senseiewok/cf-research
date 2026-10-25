#!/usr/bin/env python3
"""Unit test suite for AlphaEvolve-CF Pareto Regimen Optimizer, SMT Prover & Saltelli-Jansen Sobol Analyzer (T-0426)."""
from __future__ import annotations

from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))

from alpha_evolve_pareto_regimen import (
    PRIMARY_EVIDENCE_LEDGER,
    VERIFIED_TRIAL_ANCHORS,
    build_discovery_dossier,
    calculate_dnase_ionic_efficiency,
    calculate_thixotropic_pore_size,
    compute_pareto_hypervolume_2d,
    compute_sobol_decomposition,
    constrained_dominates,
    detect_math_backends,
    evolve_pareto_front,
    fast_non_dominated_sort,
    run_benchmark_comparison,
)


class TestAlphaEvolveParetoRegimen(unittest.TestCase):
    """Verify SMT unsat_core invariants, epistatic synergy (>6.5x), QMC Sobol variance, and NSGA-II Hypervolume."""

    def test_dnase_ionic_denaturation_guard(self) -> None:
        co_nebulized = calculate_dnase_ionic_efficiency(8.0, 8.0)
        lagged_51min = calculate_dnase_ionic_efficiency(7.0, 7.85)
        self.assertLess(co_nebulized, 0.01)
        self.assertGreater(lagged_51min, 0.90)

    def test_thixotropic_pore_window_and_shear_ejection(self) -> None:
        pore_during, aer_during = calculate_thixotropic_pore_size(15.0, -5.0, 0.95)
        self.assertEqual(aer_during, 0.12)
        self.assertGreater(pore_during, 100.0)

        pore_window, aer_window = calculate_thixotropic_pore_size(15.0, 12.0, 0.95)
        self.assertEqual(aer_window, 1.00)
        self.assertGreater(pore_window, 95.0)

        pore_late, _ = calculate_thixotropic_pore_size(15.0, 90.0, 0.95)
        self.assertLess(pore_late, 62.0)

    def test_benchmark_monotherapy_and_simultaneous_fail_with_unsat_core_while_pareto_passes(self) -> None:
        benchmarks = run_benchmark_comparison()
        mono = benchmarks["monotherapy"]
        simul = benchmarks["naive_simultaneous"]
        pareto = benchmarks["pareto_optimal"]

        self.assertFalse(mono.invariants_passed)
        self.assertGreaterEqual(len(mono.unsat_core), 2)
        self.assertLess(mono.delta_isc_pct_wt, 2.5)

        self.assertFalse(simul.invariants_passed)
        self.assertEqual(simul.violation_count, 4)
        self.assertEqual(len(simul.unsat_core), 4)
        self.assertTrue(any("UNSAT_INV1" in c for c in simul.unsat_core))
        self.assertLess(simul.delta_isc_pct_wt, 2.5)

        self.assertTrue(pareto.invariants_passed)
        self.assertEqual(pareto.violation_count, 0)
        self.assertEqual(len(pareto.unsat_core), 0)
        self.assertGreaterEqual(pareto.delta_isc_pct_wt, 16.8)
        self.assertGreaterEqual(pareto.asl_depth_um, 7.25)
        self.assertLessEqual(pareto.delta_lci_2p5, -1.85)
        self.assertGreaterEqual(pareto.epistatic_synergy_ratio, 6.50)
        self.assertTrue(constrained_dominates(pareto, simul))
        self.assertTrue(constrained_dominates(pareto, mono))

    def test_saltelli_jansen_qmc_sobol_variance_decomposition(self) -> None:
        sobol = compute_sobol_decomposition(n_samples=128)
        total = sobol["sum_first_order_additive_S"] + sobol["total_epistatic_variance_share"]
        self.assertAlmostEqual(total, 1.0, places=3)
        self.assertGreater(sobol["total_epistatic_variance_share"], 0.85)
        self.assertGreater(sobol["ST4_total_vector_infiltration"], sobol["S4_vector_infiltration"])

    def test_deb_nsga2_fast_non_dominated_sort_and_hypervolume(self) -> None:
        benchmarks = run_benchmark_comparison()
        pop = [benchmarks["pareto_optimal"], benchmarks["naive_simultaneous"], benchmarks["monotherapy"]]
        fronts = fast_non_dominated_sort(pop)
        self.assertEqual(fronts[0], [0])

        front = evolve_pareto_front(pop_size=30, generations=15)
        self.assertGreaterEqual(len(front), 1)
        for sol in front:
            self.assertTrue(sol.invariants_passed)
            self.assertEqual(sol.pareto_rank, 1)

        hv = compute_pareto_hypervolume_2d(front)
        self.assertGreater(hv, 300.0)

        dossier = build_discovery_dossier()
        self.assertEqual(dossier["status"], "PASS")
        self.assertEqual(dossier["task_id"], "T-0426")
        self.assertIn("moo_engine", detect_math_backends())
        self.assertGreaterEqual(len(PRIMARY_EVIDENCE_LEDGER), 6)
        self.assertIn("NCT05668741", VERIFIED_TRIAL_ANCHORS)


if __name__ == "__main__":
    unittest.main()
