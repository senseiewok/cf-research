#!/usr/bin/env python3
"""Unit tests for Inhaled CFTR mRNA-LNP PEG Desorption Kinetics and Steric Mucus Penetration (T-0463)."""

from __future__ import annotations

import math
from pathlib import Path
import sys
import unittest

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
for p in (str(HERE), str(REPO_ROOT)):
    if p not in sys.path:
        sys.path.insert(0, p)

from lnp_peg_shedding_mucus_kinetics import (
    ANCHOR_DMG_C14,
    ANCHOR_DSPE_C18,
    DOI_AIRWAY_MUCUS,
    DOI_CODON_OPTIMALITY,
    DOI_PEGYLATION,
    DOI_PSEUDOURIDINE,
    DOI_SORT_LNP,
    LNPFormulationConfig,
    LNPPEGSheddingMucusEngine,
    MucusMeshEnvironment,
    NCT_4D710_AAV_TRIAL,
    NCT_VX522_MRNA_TRIAL,
    NLM_ATTRIBUTION,
    PDB_ACTIVE_CFTR,
    PDB_DEPHOSPHO_CFTR,
    RESEARCH_BOUNDARY,
)


class TestLNPPEGSheddingMucusKinetics(unittest.TestCase):
    """Test suite for biophysical PEG desorption and mucus transit modeling."""

    def setUp(self) -> None:
        self.default_config = LNPFormulationConfig()
        self.default_mucus = MucusMeshEnvironment()
        self.engine = LNPPEGSheddingMucusEngine(config=self.default_config, mucus_env=self.default_mucus)

    def test_anchor_desorption_kinetics(self) -> None:
        """Assert C14 DMG-PEG desorbs significantly faster than C18 DSPE-PEG."""
        k_c14 = self.engine.calculate_desorption_rate_constant(ANCHOR_DMG_C14)
        k_c18 = self.engine.calculate_desorption_rate_constant(ANCHOR_DSPE_C18)

        # C14 rate must be at least 10x faster than C18
        self.assertGreater(k_c14, 0.15)
        self.assertLess(k_c18, 0.02)
        self.assertGreater(k_c14 / k_c18, 10.0)

        # Half-lives must be biophysically consistent
        t_half_c14 = math.log(2.0) / k_c14
        t_half_c18 = math.log(2.0) / k_c18
        self.assertAlmostEqual(t_half_c14, 3.0, delta=0.5)
        self.assertGreater(t_half_c18, 50.0)

    def test_brush_extension_scaling(self) -> None:
        """Assert PEG brush thickness scales monotonically with PEG molar percentage."""
        l_low = self.engine.calculate_peg_brush_thickness(0.5)
        l_mid = self.engine.calculate_peg_brush_thickness(2.0)
        l_high = self.engine.calculate_peg_brush_thickness(4.0)

        self.assertGreater(l_mid, l_low)
        self.assertGreater(l_high, l_mid)
        # Bounded within de Gennes mushroom-brush range
        self.assertGreaterEqual(l_low, 1.0)
        self.assertLessEqual(l_high, 7.5)

    def test_mucus_diffusivity_and_charge_unmasking(self) -> None:
        """Assert shielded LNP penetrates mucus but unmasked cationic core causes severe entrapment."""
        state_0h = self.engine.evaluate_kinetic_state(0.0)
        state_6h = self.engine.evaluate_kinetic_state(6.0)

        # At t=0, near-neutral zeta potential and good mucus penetration
        self.assertLessEqual(state_0h.zeta_potential_mv, 0.0)
        self.assertGreaterEqual(state_0h.mucus_diffusivity_ratio, 0.15)
        self.assertEqual(state_0h.entrapment_status, "PENETRATING")

        # After 6 hours of C14 desorption (~2 half-lives), charge shifts positive and diffusivity drops
        self.assertGreater(state_6h.zeta_potential_mv, 10.0)
        self.assertLess(state_6h.mucus_diffusivity_ratio, state_0h.mucus_diffusivity_ratio)
        self.assertIn(state_6h.entrapment_status, ["RETARDED", "ENTRAPPED"])

    def test_cellular_uptake_and_competence_tradeoff(self) -> None:
        """Assert endosomal competence window balances mucus penetration against uptake inhibition."""
        uptake_0h = self.engine.evaluate_transit_and_uptake(0.0)
        uptake_2h = self.engine.evaluate_transit_and_uptake(2.0)

        # At t=0, penetration has not occurred yet and high PEG inhibits cellular uptake
        self.assertLess(uptake_0h.mucosal_penetration_fraction, 0.10)
        self.assertLess(uptake_0h.uptake_competence_score, 0.20)

        # At t=2h, mucus penetration is advanced and PEG shedding unlocks cellular uptake competence
        self.assertGreater(uptake_2h.mucosal_penetration_fraction, 0.50)
        self.assertGreater(uptake_2h.uptake_competence_score, uptake_0h.uptake_competence_score)

    def test_pareto_optimality_classification(self) -> None:
        """Assert C14 formulation at 2.0 mol% achieves OPTIMAL_DELIVERY while C18 is OVER_SHIELDED."""
        eval_c14 = self.engine.evaluate_formulation()
        self.assertEqual(eval_c14.pareto_optimality_verdict, "OPTIMAL_DELIVERY")
        self.assertGreater(eval_c14.effective_delivery_yield_pct, 5.0)

        # Over-shielded C18 formulation
        c18_cfg = LNPFormulationConfig(name="C18-Test", anchor=ANCHOR_DSPE_C18, peg_mol_pct=2.0)
        c18_engine = LNPPEGSheddingMucusEngine(config=c18_cfg, mucus_env=self.default_mucus)
        eval_c18 = c18_engine.evaluate_formulation()
        self.assertEqual(eval_c18.pareto_optimality_verdict, "OVER_SHIELDED")

        # Premature shedding (<1.0 mol%)
        premature_cfg = LNPFormulationConfig(name="Low-PEG-Test", peg_mol_pct=0.5)
        premature_engine = LNPPEGSheddingMucusEngine(config=premature_cfg, mucus_env=self.default_mucus)
        eval_premature = premature_engine.evaluate_formulation()
        self.assertEqual(eval_premature.pareto_optimality_verdict, "PREMATURE_SHEDDING")

    def test_provenance_and_boundaries(self) -> None:
        """Assert required regulatory, structural, and primary literature grounding."""
        evaluation = self.engine.evaluate_formulation()

        self.assertIn(NCT_VX522_MRNA_TRIAL, evaluation.provenance_nct)
        self.assertIn(NCT_4D710_AAV_TRIAL, evaluation.provenance_nct)
        self.assertIn(PDB_DEPHOSPHO_CFTR, evaluation.provenance_pdb)
        self.assertIn(PDB_ACTIVE_CFTR, evaluation.provenance_pdb)
        self.assertIn(DOI_AIRWAY_MUCUS, evaluation.provenance_doi)
        self.assertIn(DOI_CODON_OPTIMALITY, evaluation.provenance_doi)
        self.assertIn(DOI_PSEUDOURIDINE, evaluation.provenance_doi)
        self.assertIn(DOI_PEGYLATION, evaluation.provenance_doi)
        self.assertIn(DOI_SORT_LNP, evaluation.provenance_doi)

        self.assertIn("ClinicalTrials.gov does not evaluate safety or efficacy", NLM_ATTRIBUTION)
        self.assertIn("research, not medical advice", RESEARCH_BOUNDARY)


if __name__ == "__main__":
    unittest.main()
