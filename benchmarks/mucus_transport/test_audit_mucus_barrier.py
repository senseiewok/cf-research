#!/usr/bin/env python3
"""Unit and planted-defect tests for benchmarks/mucus_transport/audit_mucus_barrier.py."""
from __future__ import annotations

import json
import subprocess
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from audit_mucus_barrier import (  # noqa: E402
    DEFAULT_HEURISTIC,
    DEFAULT_MEASURED,
    REPRESENTATIVE_VECTORS,
    HeuristicModelParams,
    MeasuredLiteratureParams,
    VectorDesignParams,
    audit_dimensional_consistency,
    audit_discontinuities,
    audit_parity_and_invariants,
    independent_evaluate_mucus_barrier,
    prove_electrostatic_confounding,
    prove_single_size_non_identifiability,
    recover_pore_and_viscosity_from_two_neutral_probes,
    run_full_audit,
    run_planted_defect_checks,
    run_source_split_uncertainty_audit,
    run_two_probe_identifiability_proofs,
)

AUDIT_SCRIPT = HERE / "audit_mucus_barrier.py"


class TestNarrowCorrectnessAndDimensionality(unittest.TestCase):
    def test_dimensional_oracle_and_um_scaling(self) -> None:
        dim = audit_dimensional_consistency()
        self.assertTrue(dim["passed"])
        self.assertEqual(dim["d_water_dim"], [0, 2, -1, 0])
        self.assertEqual(dim["d_eff_dim"], [0, 2, -1, 0])
        self.assertEqual(dim["transit_time_dim"], [0, 0, 1, 0])
        self.assertAlmostEqual(dim["um_to_m_squared_scaling_ratio"], 1e12, delta=1.0)

    def test_parity_limits_monotonicity_and_active_inputs(self) -> None:
        par = audit_parity_and_invariants()
        self.assertTrue(par["passed"])
        self.assertTrue(par["cationic_stricter_than_anionic_pinned"])
        self.assertTrue(all(c["passed"] for c in par["parity_checks"]))
        self.assertTrue(all(par["limit_checks"].values()))
        self.assertTrue(all(par["monotonic_checks"].values()))
        self.assertTrue(all(par["active_input_checks"].values()))

    def test_discontinuities_detected_at_regime_and_zeta_kinks(self) -> None:
        disc = audit_discontinuities()
        self.assertTrue(disc["passed"])
        jumps = disc["c0_piecewise_shield_jumps"]
        self.assertEqual(len(jumps), 2)
        self.assertTrue(all(j["discontinuity_detected"] for j in jumps))
        self.assertTrue(disc["c1_zeta_threshold_kink"]["kink_detected"])

    def test_all_five_planted_defects_caught(self) -> None:
        res = run_planted_defect_checks()
        self.assertTrue(res["passed"])
        self.assertEqual(res["total"], 5)
        self.assertEqual(res["caught"], 5)

    def test_steric_singularity_and_vanishing_diffusivity_json_safe(self) -> None:
        blocked = independent_evaluate_mucus_barrier(
            VectorDesignParams("blocked", 140.0, 0.0, 0.0, 2.0),
            DEFAULT_MEASURED,
            DEFAULT_HEURISTIC,
        )
        self.assertTrue(blocked.steric_blocked)
        self.assertIsNone(blocked.transit_time_s)
        self.assertFalse(blocked.transit_time_is_finite)
        self.assertEqual(blocked.transit_time_reason, "sterically_blocked")

        # Strongly cationic bare particle (+300 mV) drives d_eff below 1e-24 m^2/s cutoff.
        vanishing = independent_evaluate_mucus_barrier(
            VectorDesignParams("vanish", 80.0, 300.0, 0.0, 2.0),
            DEFAULT_MEASURED,
            DEFAULT_HEURISTIC,
        )
        self.assertFalse(vanishing.steric_blocked)
        self.assertIsNone(vanishing.transit_time_s)
        self.assertFalse(vanishing.transit_time_is_finite)
        self.assertEqual(vanishing.transit_time_reason, "vanishing_diffusivity")


class TestSourceSplitUncertaintyAndIdentifiability(unittest.TestCase):
    def test_source_split_uncertainty_has_no_blended_bar_and_handles_boundary(self) -> None:
        unc = run_source_split_uncertainty_audit(REPRESENTATIVE_VECTORS)
        self.assertEqual(unc["epistemic_label"], "SYNTHETIC_KNOWN_ANSWER")
        self.assertFalse(unc["has_blended_uncertainty_bar"])

        for vec in REPRESENTATIVE_VECTORS:
            prof = unc["profiles"][vec.name]
            self.assertIn("measured_only", prof)
            self.assertIn("heuristic_only", prof)
            self.assertIn("joint_measured_plus_heuristic", prof)
            # Every numeric value across all three grids must serialize strictly without Inf/NaN.
            json.dumps(prof, allow_nan=False)

        large_prof = unc["profiles"]["large_lnp_120nm_neutral"]
        self.assertGreater(large_prof["measured_only"]["infinite_transit_samples"], 0)
        self.assertTrue(large_prof["joint_measured_plus_heuristic"]["verdict_flips_occurred"])

        peg_prof = unc["profiles"]["peg_lnp_80nm_neutral"]
        self.assertTrue(peg_prof["joint_measured_plus_heuristic"]["verdict_flips_occurred"])
        self.assertLess(peg_prof["joint_measured_plus_heuristic"]["verdict_stability"], 1.0)

    def test_single_size_non_identifiability_and_electrostatic_confounding(self) -> None:
        nid = prove_single_size_non_identifiability()
        self.assertEqual(nid["epistemic_label"], "SYNTHETIC_KNOWN_ANSWER")
        self.assertTrue(nid["not_empirical_validation"])
        self.assertTrue(nid["indistinguishable_proven"])
        self.assertLess(nid["relative_diff"], 1e-12)

        conf = prove_electrostatic_confounding()
        self.assertEqual(conf["epistemic_label"], "SYNTHETIC_KNOWN_ANSWER")
        self.assertTrue(conf["not_empirical_validation"])
        self.assertTrue(conf["confounded_proven"])
        self.assertLess(conf["relative_diff"], 1e-12)

    def test_two_probe_algebraic_recovery_and_closed_fail_guards(self) -> None:
        tp = run_two_probe_identifiability_proofs()
        self.assertEqual(tp["epistemic_label"], "SYNTHETIC_KNOWN_ANSWER")
        self.assertTrue(tp["not_empirical_validation"])
        self.assertTrue(tp["passed"])
        self.assertTrue(tp["steric_blocked_guard_triggered"])

        for bad_args in (
            (0.0, 80.0, 0.2, 0.1),
            (80.0, 40.0, 0.2, 0.1),
            (40.0, 80.0, 0.1, 0.1),
            (40.0, 80.0, 0.1, 0.15),
            (40.0, 80.0, 0.1, 0.0),
        ):
            with self.assertRaises(ValueError):
                recover_pore_and_viscosity_from_two_neutral_probes(*bad_args)


class TestFullReportAndCLI(unittest.TestCase):
    def test_full_audit_strict_json(self) -> None:
        rep = run_full_audit()
        self.assertEqual(rep["status"], "PASS")
        dumped = json.dumps(rep, allow_nan=False)
        self.assertIn("SYNTHETIC_KNOWN_ANSWER", dumped)

    def test_cli_check_and_json(self) -> None:
        res = subprocess.run(
            [sys.executable, str(AUDIT_SCRIPT), "--check", "--json"],
            capture_output=True,
            text=True,
        )
        self.assertEqual(res.returncode, 0, msg=res.stderr)
        payload = json.loads(res.stdout)
        self.assertEqual(payload["status"], "PASS")
        self.assertEqual(payload["planted_defect_checks"]["caught"], 5)


if __name__ == "__main__":
    unittest.main()
