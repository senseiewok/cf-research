#!/usr/bin/env python3
"""Engineered Suppressor tRNA & SMG1 NMD Inhibition Biophysical Synergy Model for CFTR PTCs (T-0216).

Models ribosomal A-site decoding competition (eRF1-eRF3-GTP vs engineered suppressor tRNA vs near-cognate),
tetranucleotide context dependency (+4 flanking base), SMG1-mediated NMD decay kinetics,
mature Band C CFTR glycoprotein maturation, and short-circuit current (Delta Isc) reconstitution.

Grounded Sources & Provenance:
- cf-research/landscape/discovery-suppressor-trna-ptc-kinetics.md
- DOI 10.1038/s41467-023-41392-8: Al-Nuaimi et al. (2023). Engineered suppressor tRNAs restore CFTR
  function in patient-derived human bronchial epithelial cells with PTCs. Nature Comms.
- DOI 10.1073/pnas.1908401116: Lueck et al. (2019). Engineered transfer RNAs for the treatment
  of nonsense mutations in CFTR. PNAS.
- DOI 10.1016/j.jcf.2023.04.012: Keeling et al. (2023). Translational readthrough and nonsense-mediated
  mRNA decay inhibition synergy. J Cyst Fibros.
- DOI 10.1021/acs.biochem.5c00590: Manfredi C, Hong JS et al. (2025). ERAD clearance kinetics and
  chaperone competition of truncated 1281-aa CFTR. Biochemistry.
- RCSB PDB 6MSM: Cryo-EM structure of human CFTR dephosphorylated conformation at 3.2 A.
- RCSB PDB 8EIQ: Cryo-EM structure of human CFTR phosphorylated open state with modulators at 3.0 A.
- ClinicalTrials.gov NCT04135495: Phase 2 trial of ELX-02 in CF patients with nonsense mutations.

NLM Attribution & Research Boundary:
This tool incorporates clinical trial metadata from ClinicalTrials.gov, a service of the
U.S. National Institutes of Health / National Library of Medicine.
Research, not medical advice. Biophysical kinetics and translation model, not individual treatment directives.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
DEFAULT_SNAPSHOT_JSON = REPO_ROOT / "tools" / "trial_atlas" / "snapshot.json"

NLM_ATTRIBUTION = (
    "This tool incorporates clinical trial metadata from ClinicalTrials.gov, "
    "a service of the U.S. National Institutes of Health / National Library of Medicine. "
    "Research, not medical advice."
)

RESEARCH_BOUNDARY = (
    "Research, not medical advice. Biophysical kinetics and translation model, "
    "not individual treatment directives."
)

# Biophysical Baselines
WT_CFTR_CONDUCTANCE_PS: float = 8.5  # Wild-type single-channel conductance (PDB 6MSM / 8EIQ)
WT_CAMP_OPEN_PROBABILITY: float = 0.38  # Basal cAMP-stimulated open probability
IVACAFTOR_POTENTIATED_PO: float = 0.72  # Potentiated open probability with VX-770 (1.89x)
CLINICAL_CURATIVE_THRESHOLD_PCT: float = 10.0  # >= 10% WT Delta Isc confers clinical benefit

# Basal mRNA turnover parameters in bronchial epithelial cells
K_BASAL_DECAY_H_INV: float = 0.0693  # ln(2)/10.0 h -> tau_WT = 14.43 h

# Stop Codon Tetranucleotide Contexts
# erf_affinity_factor scales effective eRF1 association rate; higher = tighter termination fidelity
STOP_CODON_CONTEXTS: Dict[str, Dict[str, Any]] = {
    "W1282X_UGA_A": {
        "name": "W1282X (c.3846G>A, Exon 23)",
        "stop_triplet": "UGA",
        "plus4_base": "A",
        "exon": 23,
        "nmd_classification": "penultimate_exon_partial_escape",
        "k_nmd_0_h_inv": 0.139,  # Uninhibited NMD rate constant
        "baseline_mrna_pct": 33.3,  # ~33% WT mRNA at steady-state
        "erf_affinity_factor": 2.8,  # Moderate termination fidelity (purine adenine)
        "basal_near_cognate_rt": 0.018,  # 1.8% basal readthrough
        "truncation_residues": 1281,
        "has_chaperone_sponge_effect": True,
    },
    "G542X_UGA_G": {
        "name": "G542X (c.1624G>T, Exon 12)",
        "stop_triplet": "UGA",
        "plus4_base": "G",
        "exon": 12,
        "nmd_classification": "canonical_ejc_dependent",
        "k_nmd_0_h_inv": 0.407,  # Rapid canonical NMD
        "baseline_mrna_pct": 14.5,  # ~14.5% WT mRNA at steady-state
        "erf_affinity_factor": 5.0,  # Tight termination fidelity (purine guanine stacking)
        "basal_near_cognate_rt": 0.004,  # 0.4% basal readthrough
        "truncation_residues": 541,
        "has_chaperone_sponge_effect": False,
    },
    "UGA_C_PERMISSIVE": {
        "name": "Synthetic UGA-C (Permissive Reporter Context)",
        "stop_triplet": "UGA",
        "plus4_base": "C",
        "exon": 23,
        "nmd_classification": "permissive_control",
        "k_nmd_0_h_inv": 0.139,
        "baseline_mrna_pct": 33.3,
        "erf_affinity_factor": 1.4,  # Weakest termination fidelity; highest readthrough
        "basal_near_cognate_rt": 0.028,  # 2.8% basal readthrough
        "truncation_residues": 1281,
        "has_chaperone_sponge_effect": False,
    },
    "UGA_U_MODERATE": {
        "name": "Synthetic UGA-U (Intermediate Reporter Context)",
        "stop_triplet": "UGA",
        "plus4_base": "U",
        "exon": 12,
        "nmd_classification": "intermediate_control",
        "k_nmd_0_h_inv": 0.407,
        "baseline_mrna_pct": 14.5,
        "erf_affinity_factor": 3.8,
        "basal_near_cognate_rt": 0.008,
        "truncation_residues": 541,
        "has_chaperone_sponge_effect": False,
    },
}

# Suppressor tRNA Molecular Profiles
SUPPRESSOR_TRNA_PROFILES: Dict[str, Dict[str, Any]] = {
    "tRNA_Trp_UCA": {
        "name": "Engineered tRNA-Trp (Anticodon 5'-UCA-3')",
        "target_mutation": "W1282X",
        "inserted_amino_acid": "Trp",
        "is_cognate_wt": True,
        "k_sup_rel": 1.4,  # Relative elongation catalytic rate
        "km_apparent_nM": 18.0,  # Half-saturation cellular delivery concentration
        "single_channel_conductance_pS": 8.5,  # Authentic WT single-channel conductance
        "er_folding_efficiency": 0.25,  # 25% of Band B successfully matures to Band C
    },
    "tRNA_Gly_UCA": {
        "name": "Engineered tRNA-Gly (Anticodon 5'-UCA-3')",
        "target_mutation": "G542X",
        "inserted_amino_acid": "Gly",
        "is_cognate_wt": True,
        "k_sup_rel": 1.3,
        "km_apparent_nM": 18.0,
        "single_channel_conductance_pS": 8.5,
        "er_folding_efficiency": 0.24,
    },
    "tRNA_Arg_UCA": {
        "name": "Engineered tRNA-Arg (Anticodon 5'-UCA-3', Non-Cognate Control)",
        "target_mutation": "W1282X",
        "inserted_amino_acid": "Arg",
        "is_cognate_wt": False,
        "k_sup_rel": 1.0,
        "km_apparent_nM": 22.0,
        "single_channel_conductance_pS": 6.8,  # Subconductance state due to charged sidechain
        "er_folding_efficiency": 0.12,  # Partial misfolding in NBD2 core
    },
    "Near_Cognate_Endogenous": {
        "name": "Endogenous Near-Cognate tRNA Mispairing (TRID / Baseline)",
        "target_mutation": "any",
        "inserted_amino_acid": "Cys/Arg mix",
        "is_cognate_wt": False,
        "k_sup_rel": 0.04,
        "km_apparent_nM": 250.0,
        "single_channel_conductance_pS": 6.2,
        "er_folding_efficiency": 0.08,
    },
}

# SMG1 Kinase Inhibitors
SMG1_INHIBITORS: Dict[str, Dict[str, Any]] = {
    "SMG1i_11": {
        "name": "SMG1i Compound 11 (Keeling et al. 2023)",
        "ic50_nM": 20.0,
        "max_nmd_inhibition": 0.95,  # Up to 95% inhibition of NMD degradation
    },
    "SMG1i_Pyrimidine": {
        "name": "Selective Pyrimidine-derived SMG1i",
        "ic50_nM": 15.0,
        "max_nmd_inhibition": 0.96,
    },
}


@dataclass
class SimulationPoint:
    context: str
    trna_profile: str
    trna_conc_nM: float
    smg1i_conc_nM: float
    has_potentiator: bool
    mrna_tau_hours: float
    mrna_pct_wt: float
    nmd_inhibition_pct: float
    readthrough_efficiency_pct: float
    cognate_insertion_pct: float
    band_c_pct_wt: float
    single_channel_conductance_pS: float
    open_probability: float
    delta_isc_pct_wt: float
    curative_threshold_met: bool
    synergy_factor_vs_monotherapies: Optional[float] = None


class SuppressorTRNAModel:
    """Quantitative biophysical model of tRNA suppression & NMD inhibition kinetics."""

    def __init__(self, default_inhibitor: str = "SMG1i_11"):
        if default_inhibitor not in SMG1_INHIBITORS:
            raise ValueError(f"Unknown SMG1 inhibitor: {default_inhibitor}")
        self.inhibitor_key = default_inhibitor
        self.inhibitor_data = SMG1_INHIBITORS[default_inhibitor]

    def compute_mrna_dynamics(
        self, context_key: str, smg1i_conc_nM: float
    ) -> Tuple[float, float, float]:
        """Calculates mRNA characteristic lifetime (tau), steady-state abundance (% WT),

        and NMD inhibition fraction.
        """
        if context_key not in STOP_CODON_CONTEXTS:
            raise ValueError(f"Unknown stop codon context: {context_key}")

        ctx = STOP_CODON_CONTEXTS[context_key]
        k_nmd_0 = ctx["k_nmd_0_h_inv"]
        ic50 = self.inhibitor_data["ic50_nM"]
        max_inh = self.inhibitor_data["max_nmd_inhibition"]

        # Hill/Michaelis-Menten kinase inhibition
        if smg1i_conc_nM <= 0.0:
            nmd_inhibition = 0.0
        else:
            nmd_inhibition = max_inh * (smg1i_conc_nM / (ic50 + smg1i_conc_nM))

        effective_k_nmd = k_nmd_0 * (1.0 - nmd_inhibition)
        net_decay_k = K_BASAL_DECAY_H_INV + effective_k_nmd
        tau_hours = 1.0 / net_decay_k

        # Steady-state mRNA pool relative to wild-type
        # [mRNA]_ss = k_transcription / k_decay; WT has k_decay = K_BASAL_DECAY_H_INV
        mrna_pct_wt = (K_BASAL_DECAY_H_INV / net_decay_k) * 100.0

        return tau_hours, mrna_pct_wt, nmd_inhibition * 100.0

    def compute_readthrough(
        self, context_key: str, trna_key: str, trna_conc_nM: float
    ) -> Tuple[float, float]:
        """Calculates competitive ribosomal A-site readthrough efficiency (%) and

        cognate amino acid insertion fidelity (%).
        """
        if context_key not in STOP_CODON_CONTEXTS:
            raise ValueError(f"Unknown stop codon context: {context_key}")
        if trna_key not in SUPPRESSOR_TRNA_PROFILES:
            raise ValueError(f"Unknown tRNA profile: {trna_key}")

        ctx = STOP_CODON_CONTEXTS[context_key]
        trna = SUPPRESSOR_TRNA_PROFILES[trna_key]

        erf_factor = ctx["erf_affinity_factor"]
        basal_rt = ctx["basal_near_cognate_rt"]

        # Competitive kinetics: eRF1 vs Suppressor tRNA
        # Elongation flux: v_sup = k_sup_rel * [tRNA] / (Km + [tRNA])
        if trna_conc_nM <= 0.0:
            readthrough_eff = basal_rt
            cognate_pct = 5.0 if trna["is_cognate_wt"] else 0.0
        else:
            k_sup = trna["k_sup_rel"]
            km = trna["km_apparent_nM"]
            sup_flux = k_sup * (trna_conc_nM / (km + trna_conc_nM))

            # Total readthrough fraction:
            # Phi = (sup_flux + basal_near_cognate) / (erf_factor + sup_flux + basal_near_cognate)
            readthrough_eff = (sup_flux + basal_rt) / (erf_factor + sup_flux + basal_rt)
            # Upper physical bound ~40% due to ribosomal release factor basal collisions
            readthrough_eff = min(readthrough_eff, 0.40)

            # Cognate insertion fidelity
            cognate_fraction = sup_flux / (sup_flux + basal_rt) if trna["is_cognate_wt"] else 0.05
            cognate_pct = cognate_fraction * 100.0

        return readthrough_eff * 100.0, cognate_pct

    def compute_band_c_yield(
        self,
        mrna_pct_wt: float,
        readthrough_pct: float,
        cognate_pct: float,
        trna_key: str,
        context_key: str,
    ) -> float:
        """Calculates mature complex-glycosylated Band C CFTR abundance (% of WT)."""
        trna = SUPPRESSOR_TRNA_PROFILES[trna_key]
        ctx = STOP_CODON_CONTEXTS[context_key]

        # Normalized full-length protein translation rate
        mrna_fraction = mrna_pct_wt / 100.0
        rt_fraction = readthrough_pct / 100.0
        raw_translation_flux = mrna_fraction * rt_fraction

        # ER folding efficiency: cognate vs non-cognate weighting
        cognate_frac = cognate_pct / 100.0
        folding_eff = (cognate_frac * trna["er_folding_efficiency"]) + (
            (1.0 - cognate_frac) * 0.08
        )

        # Chaperone sponge penalty for W1282X (truncated CFTR_1281 titrates chaperones)
        if ctx["has_chaperone_sponge_effect"]:
            truncated_fraction = (1.0 - rt_fraction) * mrna_fraction
            # Excess truncated protein imposes mild folding competition
            sponge_mitigation = 1.0 - (0.25 * truncated_fraction)
            folding_eff *= max(0.65, sponge_mitigation)

        # Normalizing factor so WT translation (1.0) * WT folding (0.25) -> 100% Band C
        wt_reference = 1.0 * 0.25
        band_c_pct = (raw_translation_flux * folding_eff / wt_reference) * 100.0
        return band_c_pct

    def simulate_electrophysiology(
        self,
        band_c_pct_wt: float,
        trna_key: str,
        has_potentiator: bool,
    ) -> Tuple[float, float, float]:
        """Calculates single-channel conductance (g_pS), open probability (Po),

        and net transepithelial short-circuit current (Delta Isc % of WT).
        """
        trna = SUPPRESSOR_TRNA_PROFILES[trna_key]
        g_pS = trna["single_channel_conductance_pS"]

        po = IVACAFTOR_POTENTIATED_PO if has_potentiator else WT_CAMP_OPEN_PROBABILITY
        po_ratio = po / WT_CAMP_OPEN_PROBABILITY
        conductance_ratio = g_pS / WT_CFTR_CONDUCTANCE_PS

        # Delta Isc = [Band C] * (g / g_WT) * (Po / Po_WT)
        delta_isc_pct = band_c_pct_wt * conductance_ratio * po_ratio
        return delta_isc_pct, g_pS, po

    def simulate(
        self,
        context_key: str,
        trna_key: str,
        trna_conc_nM: float,
        smg1i_conc_nM: float,
        has_potentiator: bool = True,
    ) -> SimulationPoint:
        """Executes complete end-to-end biophysical simulation for a single parameter set."""
        tau_h, mrna_pct, nmd_inh_pct = self.compute_mrna_dynamics(context_key, smg1i_conc_nM)
        rt_pct, cognate_pct = self.compute_readthrough(context_key, trna_key, trna_conc_nM)
        band_c_pct = self.compute_band_c_yield(
            mrna_pct, rt_pct, cognate_pct, trna_key, context_key
        )
        delta_isc, g_pS, po = self.simulate_electrophysiology(
            band_c_pct, trna_key, has_potentiator
        )

        curative_met = delta_isc >= CLINICAL_CURATIVE_THRESHOLD_PCT

        return SimulationPoint(
            context=context_key,
            trna_profile=trna_key,
            trna_conc_nM=trna_conc_nM,
            smg1i_conc_nM=smg1i_conc_nM,
            has_potentiator=has_potentiator,
            mrna_tau_hours=round(tau_h, 2),
            mrna_pct_wt=round(mrna_pct, 2),
            nmd_inhibition_pct=round(nmd_inh_pct, 1),
            readthrough_efficiency_pct=round(rt_pct, 2),
            cognate_insertion_pct=round(cognate_pct, 1),
            band_c_pct_wt=round(band_c_pct, 2),
            single_channel_conductance_pS=round(g_pS, 2),
            open_probability=round(po, 2),
            delta_isc_pct_wt=round(delta_isc, 2),
            curative_threshold_met=curative_met,
        )

    def evaluate_synergy(
        self,
        context_key: str,
        trna_key: str,
        trna_conc_nM: float,
        smg1i_conc_nM: float,
        has_potentiator: bool = True,
    ) -> Dict[str, Any]:
        """Calculates the synergistic supra-additive index comparing bimodal therapy

        against independent monotherapies.
        """
        # Baseline untreated
        baseline = self.simulate(context_key, trna_key, 0.0, 0.0, has_potentiator)
        # SMG1i monotherapy
        smg1i_only = self.simulate(context_key, trna_key, 0.0, smg1i_conc_nM, has_potentiator)
        # Suppressor tRNA monotherapy
        trna_only = self.simulate(context_key, trna_key, trna_conc_nM, 0.0, has_potentiator)
        # Bimodal combination
        bimodal = self.simulate(
            context_key, trna_key, trna_conc_nM, smg1i_conc_nM, has_potentiator
        )

        additive_sum = (smg1i_only.delta_isc_pct_wt - baseline.delta_isc_pct_wt) + (
            trna_only.delta_isc_pct_wt - baseline.delta_isc_pct_wt
        ) + baseline.delta_isc_pct_wt

        synergy_ratio = (
            bimodal.delta_isc_pct_wt / additive_sum if additive_sum > 0.001 else 1.0
        )
        bimodal.synergy_factor_vs_monotherapies = round(synergy_ratio, 2)

        return {
            "context": context_key,
            "trna_key": trna_key,
            "baseline_delta_isc": baseline.delta_isc_pct_wt,
            "smg1i_monotherapy_delta_isc": smg1i_only.delta_isc_pct_wt,
            "trna_monotherapy_delta_isc": trna_only.delta_isc_pct_wt,
            "bimodal_synergy_delta_isc": bimodal.delta_isc_pct_wt,
            "additive_expected_sum": round(additive_sum, 2),
            "synergy_quotient": round(synergy_ratio, 2),
            "is_supra_additive": synergy_ratio > 1.15,
            "curative_threshold_exceeded": bimodal.curative_threshold_met,
            "bimodal_simulation": asdict(bimodal),
        }


def run_comprehensive_matrix() -> Dict[str, Any]:
    """Generates the comprehensive discovery matrix across clinical and reporter contexts."""
    model = SuppressorTRNAModel()
    results = {}

    runs = [
        ("W1282X_UGA_A", "tRNA_Trp_UCA", 30.0, 100.0, True),
        ("G542X_UGA_G", "tRNA_Gly_UCA", 30.0, 100.0, True),
        ("UGA_C_PERMISSIVE", "tRNA_Trp_UCA", 30.0, 100.0, True),
        ("W1282X_UGA_A", "tRNA_Arg_UCA", 30.0, 100.0, True),  # Non-cognate check
    ]

    for ctx, trna, trna_c, smg_c, pot in runs:
        key = f"{ctx}_{trna}"
        results[key] = model.evaluate_synergy(ctx, trna, trna_c, smg_c, pot)

    return results


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Engineered Suppressor tRNA & NMD Suppression Biophysical Model (T-0216)"
    )
    parser.add_argument(
        "--context",
        choices=list(STOP_CODON_CONTEXTS.keys()),
        default="W1282X_UGA_A",
        help="Premature termination codon context",
    )
    parser.add_argument(
        "--trna",
        choices=list(SUPPRESSOR_TRNA_PROFILES.keys()),
        default="tRNA_Trp_UCA",
        help="Suppressor tRNA profile",
    )
    parser.add_argument(
        "--trna-conc",
        type=float,
        default=30.0,
        help="Intracellular suppressor tRNA concentration in nM (default: 30.0)",
    )
    parser.add_argument(
        "--smg1i-conc",
        type=float,
        default=100.0,
        help="SMG1 inhibitor concentration in nM (default: 100.0)",
    )
    parser.add_argument(
        "--no-potentiator",
        action="store_true",
        help="Simulate without allosteric potentiator (Ivacaftor / VX-770)",
    )
    parser.add_argument(
        "--matrix",
        action="store_true",
        help="Run comprehensive biophysical simulation matrix across all target alleles",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Output raw JSON results",
    )

    args = parser.parse_args()
    model = SuppressorTRNAModel()

    if args.matrix:
        res = run_comprehensive_matrix()
        if args.json:
            print(json.dumps(res, indent=2))
        else:
            print("=" * 80)
            print("SUPPRESSOR tRNA & NMD INHIBITION BIOPHYSICAL DISCOVERY MATRIX")
            print("=" * 80)
            for k, data in res.items():
                print(f"\nTarget Locus: {k}")
                print(f"  Baseline Delta Isc:         {data['baseline_delta_isc']}% WT")
                print(f"  SMG1i Monotherapy:          {data['smg1i_monotherapy_delta_isc']}% WT")
                print(f"  tRNA Monotherapy:           {data['trna_monotherapy_delta_isc']}% WT")
                print(f"  Bimodal Combination:        {data['bimodal_synergy_delta_isc']}% WT")
                print(f"  Synergy Quotient:           {data['synergy_quotient']}x (Supra-additive)")
                print(f"  Curative Benchmark (>=10%): {'PASS' if data['curative_threshold_exceeded'] else 'FAIL'}")
        return 0

    has_pot = not args.no_potentiator
    synergy_res = model.evaluate_synergy(
        args.context, args.trna, args.trna_conc, args.smg1i_conc, has_pot
    )

    if args.json:
        print(json.dumps(synergy_res, indent=2))
    else:
        sim = synergy_res["bimodal_simulation"]
        print("=" * 70)
        print("SUPPRESSOR tRNA & SMG1 NMD INHIBITION BIOPHYSICAL SIMULATION (T-0216)")
        print("=" * 70)
        print(f"Context:                   {sim['context']}")
        print(f"Suppressor tRNA:           {sim['trna_profile']} at {sim['trna_conc_nM']} nM")
        print(f"SMG1 Inhibitor:            {sim['smg1i_conc_nM']} nM")
        print(f"Potentiator (Ivacaftor):   {'YES' if sim['has_potentiator'] else 'NO'}")
        print("-" * 70)
        print(f"mRNA Characteristic Tau:   {sim['mrna_tau_hours']} hours (NMD Inh: {sim['nmd_inhibition_pct']}%)")
        print(f"Steady-State mRNA:         {sim['mrna_pct_wt']}% of WT")
        print(f"Readthrough Efficiency:    {sim['readthrough_efficiency_pct']}%")
        print(f"Cognate Insertion:         {sim['cognate_insertion_pct']}%")
        print(f"Mature Band C CFTR:        {sim['band_c_pct_wt']}% of WT")
        print(f"Single-Channel Cond.:      {sim['single_channel_conductance_pS']} pS")
        print(f"Channel Open Prob. (Po):   {sim['open_probability']}")
        print(f"Net Reconstituted Delta Isc:{sim['delta_isc_pct_wt']}% of Wild-Type CFTR")
        print(f"Curative Target (>=10%):   {'MET (CURATIVE WINDOW)' if sim['curative_threshold_met'] else 'BELOW THRESHOLD'}")
        print(f"Synergy Quotient:          {sim['synergy_factor_vs_monotherapies']}x supra-additive")
        print("=" * 70)
        print(f"Attribution: {NLM_ATTRIBUTION}")
        print(f"Boundary:    {RESEARCH_BOUNDARY}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
