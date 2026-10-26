#!/usr/bin/env python3
"""Unit tests for TMEM16A/ENaC dual synergy across ASL and CBF dynamics (T-0444)."""
from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from tmem16a_enac_asl_cbf_dynamics import (
    CLINICAL_TRIAL_4D710,
    CLINICAL_TRIAL_VX522,
    NLM_ATTRIBUTION,
    PDB_ACTIVE_CFTR,
    PDB_DEPHOSPHO_CFTR,
    PRIMARY_DOI_AIRWAY_MUCUS,
    PRIMARY_DOI_ASL_HYDRATION,
    PRIMARY_DOI_ENAC_STRUCT,
    PRIMARY_DOI_TMEM16A_POTENTIATION,
    PRIMARY_DOI_TMEM16A_STRUCT,
    RESEARCH_BOUNDARY,
    AirwayBiophysicalConstants,
    ENaCPharmacodynamics,
    TMEM16APharmacodynamics,
    calculate_ciliary_beat_frequency,
    calculate_mucociliary_transport_velocity,
    compute_osmotic_water_flux,
    compute_transepithelial_currents,
    evaluate_genotype_theratyping_matrix,
    simulate_dynamic_rehydration,
)


class TestTMEM16AENaCASLCBFDynamics(unittest.TestCase):
    """Test suite verifying biophysical equations, synergy mechanics, and clinical boundaries."""

    def setUp(self) -> None:
        self.constants = AirwayBiophysicalConstants()
        self.tmem_pd = TMEM16APharmacodynamics()
        self.enac_pd = ENaCPharmacodynamics()

    def test_tmem16a_calcium_sensitization(self) -> None:
        """Verify potentiator induces leftward shift in apparent calcium EC50."""
        ec50_basal = self.tmem_pd.calculate_apparent_calcium_ec50(0.0)
        self.assertAlmostEqual(ec50_basal, 1.20, places=3)

        ec50_half_max = self.tmem_pd.calculate_apparent_calcium_ec50(0.10)
        # At drug conc = Kd (0.10 uM), shift is 50% of max (0.85 * 0.5 = 42.5%)
        expected_half = 1.20 * (1.0 - 0.425)
        self.assertAlmostEqual(ec50_half_max, expected_half, places=3)

        # High potentiator concentration approaches 85% shift
        ec50_saturating = self.tmem_pd.calculate_apparent_calcium_ec50(5.0)
        self.assertLess(ec50_saturating, 0.25)

        # Negative drug concentration raises ValueError
        with self.assertRaises(ValueError):
            self.tmem_pd.calculate_apparent_calcium_ec50(-0.5)

    def test_tmem16a_open_probability(self) -> None:
        """Verify open probability responses under basal vs potentiated conditions."""
        # Basal resting Ca2+ = 0.15 uM without drug: Po should be very small (< 0.02)
        po_basal = self.tmem_pd.calculate_open_probability(0.15, 0.0)
        self.assertLess(po_basal, 0.02)

        # With 1.5 uM potentiator at resting Ca2+: Po increases significantly (> 0.20)
        po_pot = self.tmem_pd.calculate_open_probability(0.15, 1.50)
        self.assertGreater(po_pot, 0.20)

        # Zero or negative Ca2+ yields zero Po
        self.assertEqual(self.tmem_pd.calculate_open_probability(0.0, 1.0), 0.0)
        self.assertEqual(self.tmem_pd.calculate_open_probability(-0.1, 1.0), 0.0)

    def test_enac_pharmacodynamics(self) -> None:
        """Verify ENaC inhibition curves and CF disinhibition feedback."""
        # Zero blocker conc yields zero inhibition
        self.assertEqual(self.enac_pd.calculate_inhibition_fraction(0.0), 0.0)

        # Blocker conc at IC50 (0.10 uM) yields exactly 50% inhibition
        inhib_half = self.enac_pd.calculate_inhibition_fraction(0.10)
        self.assertAlmostEqual(inhib_half, 0.50, places=3)

        # High blocker conc approaches 100%
        self.assertGreater(self.enac_pd.calculate_inhibition_fraction(3.0), 0.95)

        # Negative blocker conc raises ValueError
        with self.assertRaises(ValueError):
            self.enac_pd.calculate_inhibition_fraction(-0.1)

        # Disinhibition: CFTR-null has higher baseline conductance than wild-type
        g_cf = self.enac_pd.calculate_active_conductance(4.20, cftr_function_fraction=0.0, inhibition_fraction=0.0)
        g_wt = self.enac_pd.calculate_active_conductance(4.20, cftr_function_fraction=1.0, inhibition_fraction=0.0)
        self.assertGreater(g_cf, g_wt * 2.0)

    def test_transepithelial_current_balance(self) -> None:
        """Verify current calculations and net secretory vs absorptive balance."""
        # Untreated CF: CFTR=0, TMEM16A Po low, ENaC high -> strongly negative (absorptive)
        curr_cf = compute_transepithelial_currents(
            constants=self.constants,
            cftr_fraction=0.0,
            tmem16a_po=0.015,
            enac_inhibition_fraction=0.0,
        )
        self.assertLess(curr_cf["i_sc_net_uA_cm2"], -200.0)
        self.assertGreater(curr_cf["i_na_uA_cm2"], 250.0)

        # Dual synergy: TMEM16A Po=0.40, ENaC inhib=0.85 -> positive secretory current
        curr_synergy = compute_transepithelial_currents(
            constants=self.constants,
            cftr_fraction=0.0,
            tmem16a_po=0.40,
            enac_inhibition_fraction=0.85,
        )
        self.assertGreater(curr_synergy["i_sc_net_uA_cm2"], 20.0)

    def test_osmotic_water_flux_and_compressive_stress(self) -> None:
        """Verify osmotic water flux directionality."""
        # Negative net current causes negative or zero water flux toward baseline
        dh_neg = compute_osmotic_water_flux(self.constants, -50.0, current_asl_height_um=5.0)
        self.assertLess(dh_neg, 0.0)

        # Positive net current causes positive water flux (secretion / hydration)
        dh_pos = compute_osmotic_water_flux(self.constants, +40.0, current_asl_height_um=5.0)
        self.assertGreater(dh_pos, 0.0)

    def test_cbf_and_mucociliary_velocity_transitions(self) -> None:
        """Verify CBF recovery (>= 8.0 Hz) and mucociliary transport velocity above 7.0 um."""
        # Collapsed ASL (2.5 um) yields suppressed CBF (< 3.0 Hz) and impaired velocity
        cbf_collapsed = calculate_ciliary_beat_frequency(self.constants, 2.50)
        mct_collapsed = calculate_mucociliary_transport_velocity(self.constants, cbf_collapsed, 2.50)
        self.assertLess(cbf_collapsed, 3.0)
        self.assertLess(mct_collapsed, 1.0)

        # Fully restored ASL (8.0 um) yields CBF >= 8.0 Hz and robust velocity (>= 5.0 mm/min)
        cbf_restored = calculate_ciliary_beat_frequency(self.constants, 8.00)
        mct_restored = calculate_mucociliary_transport_velocity(self.constants, cbf_restored, 8.00)
        self.assertGreaterEqual(cbf_restored, 8.0)
        self.assertGreaterEqual(mct_restored, 5.0)

        # Zero or negative ASL height returns 0.0 CBF
        self.assertEqual(calculate_ciliary_beat_frequency(self.constants, 0.0), 0.0)
        self.assertEqual(calculate_ciliary_beat_frequency(self.constants, -1.0), 0.0)

    def test_dual_synergy_demonstration(self) -> None:
        """Verify that monotherapies fail while combined therapy achieves ASL and CBF recovery."""
        # 1. Untreated CF
        sim_untreated = simulate_dynamic_rehydration(
            constants=self.constants,
            cftr_fraction=0.0,
            tmem16a_potentiator_um=0.0,
            enac_blocker_um=0.0,
            intracellular_ca_um=0.20,
            duration_min=45.0,
        )
        self.assertFalse(sim_untreated["synergistic_clearance_restored"])
        self.assertLess(sim_untreated["steady_state_asl_height_um"], 3.0)

        # 2. TMEM16A potentiator monotherapy at submaximal safe dose (1.50 uM)
        sim_tmem_alone = simulate_dynamic_rehydration(
            constants=self.constants,
            cftr_fraction=0.0,
            tmem16a_potentiator_um=1.50,
            enac_blocker_um=0.0,
            intracellular_ca_um=0.20,
            duration_min=45.0,
        )
        self.assertFalse(sim_tmem_alone["target_asl_restored"])
        self.assertLess(sim_tmem_alone["steady_state_asl_height_um"], 6.0)

        # 3. ENaC blocker monotherapy alone (0.50 uM)
        sim_enac_alone = simulate_dynamic_rehydration(
            constants=self.constants,
            cftr_fraction=0.0,
            tmem16a_potentiator_um=0.0,
            enac_blocker_um=0.50,
            intracellular_ca_um=0.20,
            duration_min=45.0,
        )
        self.assertFalse(sim_enac_alone["target_asl_restored"])
        self.assertLess(sim_enac_alone["steady_state_asl_height_um"], 6.0)

        # 4. Synergistic Combination (1.50 uM potentiator + 0.50 uM blocker)
        sim_synergy = simulate_dynamic_rehydration(
            constants=self.constants,
            cftr_fraction=0.0,
            tmem16a_potentiator_um=1.50,
            enac_blocker_um=0.50,
            intracellular_ca_um=0.20,
            duration_min=45.0,
        )
        self.assertTrue(sim_synergy["synergistic_clearance_restored"])
        self.assertGreaterEqual(sim_synergy["steady_state_asl_height_um"], 7.0)
        self.assertGreaterEqual(sim_synergy["steady_state_cbf_hz"], 8.0)
        self.assertGreaterEqual(sim_synergy["steady_state_mct_mm_min"], 5.0)

    def test_class_i_and_vii_theratyping_matrix(self) -> None:
        """Verify universal efficacy across modulator-ineligible genotypes."""
        matrix = evaluate_genotype_theratyping_matrix()
        for gid in ["W1282X", "G542X", "CFTRdele2,3"]:
            self.assertIn(gid, matrix)
            regimens = matrix[gid]["regimens"]
            self.assertFalse(regimens["Untreated"]["synergistic_clearance_restored"])
            self.assertFalse(regimens["TMEM16A_Monotherapy_Submaximal"]["synergistic_clearance_restored"])
            self.assertFalse(regimens["ENaC_Monotherapy_Submaximal"]["synergistic_clearance_restored"])
            self.assertTrue(regimens["Dual_Synergy_Regimen"]["synergistic_clearance_restored"])
            self.assertGreaterEqual(regimens["Dual_Synergy_Regimen"]["asl_height_um"], 7.0)
            self.assertGreaterEqual(regimens["Dual_Synergy_Regimen"]["cbf_hz"], 8.0)

    def test_clinical_registry_and_literature_provenance(self) -> None:
        """Verify clinical trial identifiers, structural coordinates, and primary DOIs."""
        self.assertEqual(CLINICAL_TRIAL_VX522, "NCT05668741")
        self.assertEqual(CLINICAL_TRIAL_4D710, "NCT05248230")
        self.assertEqual(PDB_DEPHOSPHO_CFTR, "6MSM")
        self.assertEqual(PDB_ACTIVE_CFTR, "5UAK")

        self.assertEqual(PRIMARY_DOI_TMEM16A_STRUCT, "10.1038/nature24652")
        self.assertEqual(PRIMARY_DOI_ENAC_STRUCT, "10.1016/j.cell.2022.04.024")
        self.assertEqual(PRIMARY_DOI_TMEM16A_POTENTIATION, "10.1113/JP286559")
        self.assertEqual(PRIMARY_DOI_ASL_HYDRATION, "10.1164/rccm.201908-1641OC")
        self.assertEqual(PRIMARY_DOI_AIRWAY_MUCUS, "10.1152/physrev.00004.2021")

        self.assertIn("ClinicalTrials.gov does not evaluate safety or efficacy", NLM_ATTRIBUTION)
        self.assertIn("research, not medical advice", RESEARCH_BOUNDARY)


if __name__ == "__main__":
    unittest.main()
