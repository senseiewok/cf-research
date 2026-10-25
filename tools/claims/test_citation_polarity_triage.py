"""Unit tests for Section-Weighted Citation Polarity & Retraction Surveillance Model (T-0421)."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO_ROOT))

from citation_polarity_triage import (
    CitingSourceRecord,
    SectionCitationTally,
    audit_scite_polarity_provenance_text,
    compute_section_weighted_dispute_index,
    evaluate_citation_polarity_triage,
    filter_retracted_citing_sources,
)


class TestCitationPolarityTriage(unittest.TestCase):
    def test_section_weighted_dispute_index_empirical_vs_discursive(self):
        """Contrasting citations in Results/Methods must trigger EMPIRICAL_OR_METHODOLOGICAL_ATTENTION."""
        empirical_tallies = [
            SectionCitationTally(section="results", supporting=12, contrasting=3, mentioning=45),
            SectionCitationTally(section="discussion", supporting=8, contrasting=1, mentioning=30),
        ]
        res_emp = compute_section_weighted_dispute_index(empirical_tallies)
        self.assertEqual(res_emp["empirical_contrasting"], 3)
        self.assertEqual(res_emp["discursive_contrasting"], 1)
        self.assertEqual(res_emp["dispute_classification"], "EMPIRICAL_OR_METHODOLOGICAL_ATTENTION")
        self.assertGreater(res_emp["section_weighted_dispute_index"], 0.25)

        discursive_tallies = [
            SectionCitationTally(section="introduction", supporting=25, contrasting=1, mentioning=80),
            SectionCitationTally(section="discussion", supporting=15, contrasting=1, mentioning=40),
            SectionCitationTally(section="results", supporting=10, contrasting=0, mentioning=20),
        ]
        res_disc = compute_section_weighted_dispute_index(discursive_tallies)
        self.assertEqual(res_disc["empirical_contrasting"], 0)
        self.assertEqual(res_disc["discursive_contrasting"], 2)
        self.assertEqual(res_disc["dispute_classification"], "DISCURSIVE_DIFFERENCE")

        consensus_tallies = [
            SectionCitationTally(section="results", supporting=18, contrasting=0, mentioning=50),
        ]
        res_cons = compute_section_weighted_dispute_index(consensus_tallies)
        self.assertEqual(res_cons["dispute_classification"], "CONSENSUS_OR_UNCONTESTED")

    def test_filter_retracted_citing_sources(self):
        """Verify retracted citing publications are subtracted from supporting and contrasting tallies."""
        sources = [
            CitingSourceRecord(citing_doi="10.1000/retracted1", polarity="contrasting", section="results", is_retracted=True),
            CitingSourceRecord(citing_doi="10.1000/retracted2", polarity="supporting", section="discussion", editorial_notices=1),
            CitingSourceRecord(citing_doi="10.1000/clean1", polarity="supporting", section="results", is_retracted=False),
        ]
        filtered = filter_retracted_citing_sources(total_supporting=10, total_contrasting=2, citing_sources=sources)
        self.assertEqual(filtered["retracted_citing_sources_counted"], 2)
        self.assertEqual(filtered["adjusted_supporting"], 9)
        self.assertEqual(filtered["adjusted_contrasting"], 1)
        self.assertEqual(filtered["adjusted_net_support_ratio"], 0.8)

    def test_evaluate_citation_polarity_triage_end_to_end(self):
        """Verify end-to-end citation polarity evaluation with AH-8 provenance label."""
        tallies = [
            SectionCitationTally(section="methods", supporting=4, contrasting=2, mentioning=15),
            SectionCitationTally(section="discussion", supporting=16, contrasting=1, mentioning=40),
        ]
        sources = [
            CitingSourceRecord(citing_doi="10.1002/path.1597", polarity="contrasting", section="methods", is_retracted=True),
        ]
        eval_res = evaluate_citation_polarity_triage(
            "10.1056/NEJMoa1908639",
            section_tallies=tallies,
            citing_sources=sources,
            accessed_date="2026-10-24",
        )
        self.assertEqual(eval_res.target_doi, "10.1056/NEJMoa1908639")
        self.assertEqual(eval_res.dispute_classification, "EMPIRICAL_OR_METHODOLOGICAL_ATTENTION")
        self.assertEqual(eval_res.retracted_citing_sources_counted, 1)
        self.assertEqual(eval_res.adjusted_supporting, 20)
        self.assertEqual(eval_res.adjusted_contrasting, 2)
        self.assertIn("[observed: Scite ML citation classification, accessed 2026-10-24]", eval_res.provenance_label)

    def test_audit_scite_polarity_provenance_text(self):
        """Verify AH-8 provenance check catches unlabelled Scite polarity claims."""
        bad_text = "Scite reports 42 supporting citations and 2 contrasting statements for this trial."
        self.assertEqual(len(audit_scite_polarity_provenance_text(bad_text)), 1)

        good_text = (
            "Scite reports 42 supporting citations and 2 contrasting statements for this trial "
            "[observed: Scite ML citation classification, accessed 2026-10-24]."
        )
        self.assertEqual(audit_scite_polarity_provenance_text(good_text), [])

        with self.assertRaises(ValueError):
            compute_section_weighted_dispute_index([SectionCitationTally(section="results", supporting=-1)])


if __name__ == "__main__":
    unittest.main()
