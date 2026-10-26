#!/usr/bin/env python3
"""Unit tests for Prime Editing & Suppressor tRNA Kinetic Model (T-0440)."""
import math
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from prime_editing_ptc_nmd_kinetics import (
    CLINICAL_TRIAL_ELX02,
    CLINICAL_TRIAL_VX522,
    NLM_ATTRIBUTION,
    PDB_DEPHOSPHO_CFTR,
    PDB_NUCLEIC_DUPLEX,
    NMDEvasionParameters,
    PrimeEditingParameters,
    RibosomalDecodingParameters,
    evaluate_ptc_genetic_rescue,
)


class TestPrimeEditingPTCNMDKinetics(unittest.TestCase):
    """Test suite verifying prime editing flap equilibria, decoding competition, NMD, and Delta Isc."""

    def test_prime_editing_mmr_evasion_gain(self) -> None:
        """Verify that evading MMR (MLH1dn phi=0.08) elevates prime editing efficiency over wild-type MMR."""
        # Uninhibited MMR (phi = 1.0)
        pe_nommr = PrimeEditingParameters(phi_mmr_evasion=1.0)
        eff_nommr = pe_nommr.compute_prime_editing_efficiency()

        # MMR evasion with MLH1dn (phi = 0.08)
        pe_mlh1dn = PrimeEditingParameters(phi_mmr_evasion=0.08)
        eff_mlh1dn = pe_mlh1dn.compute_prime_editing_efficiency()

        self.assertGreater(eff_mlh1dn, eff_nommr * 2.0)
        self.assertGreater(eff_mlh1dn, 0.50)

    def test_tevopreq1_pegrna_stabilization(self) -> None:
        """Verify tevopreQ1 3' structural motif improves editing efficiency."""
        pe_with_tevo = PrimeEditingParameters(has_tevopreq1=True)
        pe_without_tevo = PrimeEditingParameters(has_tevopreq1=False)

        self.assertGreater(
            pe_with_tevo.compute_prime_editing_efficiency(),
            pe_without_tevo.compute_prime_editing_efficiency(),
        )

    def test_ribosomal_decoding_flanking_context(self) -> None:
        """Verify W1282X (+4 A) achieves higher readthrough than G542X (+4 G) due to termination context."""
        dec_w1282x = RibosomalDecodingParameters(variant="W1282X", suppressor_trna_concentration_um=1.5)
        dec_g542x = RibosomalDecodingParameters(variant="G542X", suppressor_trna_concentration_um=1.5)

        rt_w1282x = dec_w1282x.compute_readthrough_efficiency()
        rt_g542x = dec_g542x.compute_readthrough_efficiency()

        self.assertGreater(rt_w1282x, rt_g542x)
        self.assertGreater(rt_w1282x, 0.40)

    def test_nmd_suppression_steady_state_mrna(self) -> None:
        """Verify SMG1 inhibition and coding rescue restores steady-state mRNA from depressed baseline."""
        # Basal unrescued W1282X without SMG1 inhibitor
        nmd_basal = NMDEvasionParameters(variant="W1282X", smg1_inhibition_pct=0.0)
        mrna_basal = nmd_basal.compute_steady_state_mrna(readthrough_or_edited_fraction=0.0)
        self.assertAlmostEqual(mrna_basal, 33.3, places=1)

        # Rescued W1282X with 85% SMG1 inhibition and 50% editing
        nmd_rescued = NMDEvasionParameters(variant="W1282X", smg1_inhibition_pct=85.0)
        mrna_rescued = nmd_rescued.compute_steady_state_mrna(readthrough_or_edited_fraction=0.50)
        self.assertGreater(mrna_rescued, 70.0)

    def test_multimodal_rescue_meets_clinical_threshold(self) -> None:
        """Verify combined prime editing, suppressor tRNA, and SMG1 inhibition exceeds 10% WT Delta Isc."""
        res_w1282x = evaluate_ptc_genetic_rescue(
            variant="W1282X",
            prime_editing_enabled=True,
            has_tevopreq1=True,
            has_mlh1dn=True,
            suppressor_trna_enabled=True,
            smg1_inhibitor_enabled=True,
        )
        self.assertTrue(res_w1282x["clinical_therapeutic_threshold_met"])
        self.assertGreaterEqual(res_w1282x["restored_delta_isc_pct_wt"], 10.0)
        self.assertEqual(res_w1282x["curative_verdict"], "THERAPEUTIC")

    def test_zero_intervention_negative_control(self) -> None:
        """Verify that zero genetic or pharmacological intervention fails therapeutic threshold."""
        res_null = evaluate_ptc_genetic_rescue(
            variant="G542X",
            prime_editing_enabled=False,
            suppressor_trna_enabled=False,
            smg1_inhibitor_enabled=False,
        )
        self.assertFalse(res_null["clinical_therapeutic_threshold_met"])
        self.assertLess(res_null["restored_delta_isc_pct_wt"], 1.0)
        self.assertEqual(res_null["curative_verdict"], "SUB_THERAPEUTIC")


if __name__ == "__main__":
    unittest.main()
