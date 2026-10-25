#!/usr/bin/env python3
"""Unit tests for TMEM16A / ENaC Synergy & ASL Rehydration Model (T-0430)."""
import math
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from tmem16a_enac_synergy_model import (
    CLINICAL_TRIAL_ARO_ENAC,
    CLINICAL_TRIAL_ENAC_CLEAR,
    CLINICAL_TRIAL_ETD002,
    CLINICAL_TRIAL_IONIS_ENAC,
    NLM_ATTRIBUTION,
    PDB_ACTIVE_CFTR,
    PDB_DEPHOSPHO_CFTR,
    PDB_GATING_CFTR,
    PDB_MODULATOR_CFTR,
    AirwayFluidHomeostasisCircuit,
    TMEM16APotentiationModel,
)


class TestTMEM16AENaCSynergyModel(unittest.TestCase):
    """Test suite verifying TMEM16A allostery, ENaC suppression, and ASL height resuscitation."""

    def setUp(self) -> None:
        self.tmem_model = TMEM16APotentiationModel()
        self.circuit = AirwayFluidHomeostasisCircuit()

    def test_tmem16a_apparent_ec50_shift(self) -> None:
        """Verify leftward shift in calcium EC50 with increasing potentiator concentration."""
        ec50_0 = self.tmem_model.apparent_calcium_ec50(0.0)
        self.assertAlmostEqual(ec50_0, 1.20, places=3)

        ec50_100nM = self.tmem_model.apparent_calcium_ec50(0.10) # at EC50_pot -> half max shift (42.5%)
        expected_ec50 = 1.20 * (1.0 - 0.425)
        self.assertAlmostEqual(ec50_100nM, expected_ec50, places=3)

        ec50_saturating = self.tmem_model.apparent_calcium_ec50(10.0)
        self.assertLess(ec50_saturating, 0.25)

    def test_tmem16a_open_probability_responses(self) -> None:
        """Verify TMEM16A open probability at basal versus stimulated/potentiated conditions."""
        # Basal Ca = 0.10 uM, no drug -> low Po (< 0.02)
        po_basal = self.tmem_model.open_probability(ca_conc_uM=0.10, potentiator_conc_uM=0.0)
        self.assertLess(po_basal, 0.02)

        # Basal Ca = 0.10 uM, with 1.0 uM potentiator -> significant increase
        po_pot = self.tmem_model.open_probability(ca_conc_uM=0.10, potentiator_conc_uM=1.0)
        self.assertGreater(po_pot, 0.10)

        # Stimulated Ca = 0.50 uM, with potentiator -> high open probability (> 0.70)
        po_stim = self.tmem_model.open_probability(ca_conc_uM=0.50, potentiator_conc_uM=1.0)
        self.assertGreater(po_stim, 0.70)

    def test_cf_untreated_state_is_absorptive(self) -> None:
        """Verify untreated CFTR-null epithelia exhibits negative secretory current and collapsed ASL."""
        res = self.circuit.evaluate_regimen(
            cftr_activity_fraction=0.0,
            potentiator_conc_uM=0.0,
            enac_inhibition_fraction=0.0,
            ca_conc_uM=0.10,
        )
        # Strongly negative net secretory current (hyperabsorptive)
        self.assertLess(res["i_sc_net_uA_cm2"], -50.0)
        # Collapsed ASL height (< 3.0 um, well below 7.0 um target)
        self.assertLess(res["asl_height_um"], 3.0)
        self.assertEqual(res["target_restored"], 0.0)

    def test_tmem16a_enac_synergistic_asl_resuscitation(self) -> None:
        """Verify combined TMEM16A potentiation and ENaC inhibition resuscitates ASL height >= 7.0 um."""
        # Potentiator 2.0 uM + ENaC inhibition 85%
        res = self.circuit.evaluate_regimen(
            cftr_activity_fraction=0.0,
            potentiator_conc_uM=2.0,
            enac_inhibition_fraction=0.85,
            ca_conc_uM=0.20,
        )
        self.assertGreater(res["i_sc_net_uA_cm2"], 25.0)
        self.assertGreaterEqual(res["asl_height_um"], 7.0)
        self.assertEqual(res["target_restored"], 1.0)

    def test_regulatory_and_pdb_anchors(self) -> None:
        """Verify all clinical trial tags and PDB structural identifiers."""
        self.assertEqual(CLINICAL_TRIAL_ETD002, "NCT04488705")
        self.assertEqual(CLINICAL_TRIAL_IONIS_ENAC, "NCT03647228")
        self.assertEqual(CLINICAL_TRIAL_ARO_ENAC, "NCT04375514")
        self.assertEqual(CLINICAL_TRIAL_ENAC_CLEAR, "NCT02343445")
        self.assertEqual(PDB_DEPHOSPHO_CFTR, "6MSM")
        self.assertEqual(PDB_ACTIVE_CFTR, "5UAK")
        self.assertEqual(PDB_MODULATOR_CFTR, "8EIQ")
        self.assertEqual(PDB_GATING_CFTR, "8EJ1")
        self.assertIn("ClinicalTrials.gov does not evaluate safety or efficacy", NLM_ATTRIBUTION)


if __name__ == "__main__":
    unittest.main()
