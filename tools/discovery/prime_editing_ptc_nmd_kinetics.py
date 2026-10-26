#!/usr/bin/env python3
"""Engineered Prime Editing & Suppressor tRNA Kinetic Model for Refractory CFTR Nonsense Mutations (T-0440).

Models reverse transcription flap equilibria, mismatch repair (MMR) evasion via MLH1dn,
ribosomal A-site competitive decoding kinetics with +4 flanking base context,
SMG1-UPF1 nonsense-mediated decay (NMD) suppression, and transepithelial Delta Isc reconstitution.

Primary Evidence Grounding:
  - Nat Biotechnol (2021), DOI: 10.1038/s41587-021-01039-7
  - Cell (2021), DOI: 10.1016/j.cell.2021.09.018
  - Cell Rep Med (2024), DOI: 10.1016/j.xcrm.2024.101544
  - Nature (2023), DOI: 10.1038/s41586-023-06133-1
  - Nucleic Acids Res (2026), DOI: 10.1093/nar/gkag082
  - Am J Respir Cell Mol Biol (2019), DOI: 10.1165/rcmb.2018-0316oc
  - Biochemistry (2015), DOI: 10.1021/acs.biochem.5c00590

Clinical Trials:
  - NCT05668741 (VX-522 Phase 1/2 inhaled mRNA LNP in CF)
  - NCT04135495 (ELX-02 Phase 2 in CF patients with nonsense mutations)

Structures:
  - RCSB PDB 6MSM, 7UID
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Dict, List, Tuple

# Mandatory ClinicalTrials.gov Attribution
NLM_ATTRIBUTION = (
    "ClinicalTrials.gov does not evaluate safety or efficacy; "
    "listing on the registry does not denote government endorsement."
)

CLINICAL_TRIAL_VX522 = "NCT05668741"
CLINICAL_TRIAL_ELX02 = "NCT04135495"
PDB_DEPHOSPHO_CFTR = "6MSM"
PDB_NUCLEIC_DUPLEX = "7UID"

PRIMARY_DOI_ENG_PEGRNA = "10.1038/s41587-021-01039-7"
PRIMARY_DOI_PE_MMR_EVASION = "10.1016/j.cell.2021.09.018"
PRIMARY_DOI_PE_ORGANOIDS = "10.1016/j.xcrm.2024.101544"
PRIMARY_DOI_SUPPRESSOR_TRNA = "10.1038/s41586-023-06133-1"
PRIMARY_DOI_ACE_TRNA = "10.1093/nar/gkag082"
PRIMARY_DOI_NMD_INHIBITION = "10.1165/rcmb.2018-0316oc"
PRIMARY_DOI_W1282X_RESCUE = "10.1021/acs.biochem.5c00590"

R_GAS_KCAL_MOL_K = 0.0019872  # Gas constant in kcal/(mol*K)
T_KELVIN = 310.15              # 37 C in K
RT_KCAL = R_GAS_KCAL_MOL_K * T_KELVIN


@dataclass(frozen=True)
class PrimeEditingParameters:
    """Thermodynamic and enzymatic parameters for prime editing efficiency."""
    variant: str = "W1282X"              # Target mutation: W1282X (Exon 23) or G542X (Exon 12)
    delta_delta_g_flap_kcal_mol: float = -3.2 # Flap hybridization delta-delta-G (kcal/mol)
    has_tevopreq1: bool = True           # pegRNA 3' structural motif stabilizing reverse transcript
    phi_mmr_evasion: float = 0.08        # MMR active fraction (0.08 with MLH1dn; 1.0 without)
    k_mmr_antagonism: float = 3.5        # MMR repair penalty coefficient

    def compute_prime_editing_efficiency(self) -> float:
        """Compute net prime editing conversion efficiency in [0.0, 1.0]."""
        # Thermodynamic equilibrium of 3' flap vs unedited genomic flap
        # tevopreQ1 adds ~ -0.8 kcal/mol stabilization to pegRNA 3' extension
        ddg = self.delta_delta_g_flap_kcal_mol
        if self.has_tevopreq1:
            ddg -= 0.80

        k_flap = math.exp(-ddg / RT_KCAL)
        p_hybrid = k_flap / (1.0 + k_flap)

        # MMR antagonism factor: 1 / (1 + phi_mmr * K_mmr)
        mmr_factor = 1.0 / (1.0 + self.phi_mmr_evasion * self.k_mmr_antagonism)

        efficiency = p_hybrid * mmr_factor
        return max(0.0, min(1.0, efficiency))


@dataclass(frozen=True)
class RibosomalDecodingParameters:
    """Competitive Michaelis-Menten decoding parameters at premature termination codons."""
    variant: str = "W1282X"             # W1282X (+4 A, factor 2.8) vs G542X (+4 G, factor 5.0)
    stop_codon: str = "UGA"              # UGA stop codon
    suppressor_trna_concentration_um: float = 2.0 # Engineered ACE-tRNA concentration (uM)
    km_suppressor_um: float = 0.20       # Km for suppressor tRNA (uM)
    kcat_suppressor_s1: float = 22.0     # kcat for suppressor tRNA (1/s)
    erf_complex_concentration_um: float = 0.15 # eRF1-eRF3 complex concentration (uM)
    km_erf_um: float = 0.25              # Km for eRF1-eRF3 complex (uM)
    kcat_erf_s1: float = 12.0            # kcat for release factor termination (1/s)

    @property
    def flanking_context_factor(self) -> float:
        """Flanking nucleotide termination strength factor (+4 context)."""
        if self.variant == "W1282X":
            return 2.8  # UGA-A context: moderate termination efficiency
        return 5.0      # G542X UGA-G context: strong termination efficiency

    def compute_readthrough_efficiency(self) -> float:
        """Calculate catalytic readthrough fraction over termination in [0.0, 1.0]."""
        # Suppressor tRNA catalytic velocity
        v_sup = (self.kcat_suppressor_s1 * self.suppressor_trna_concentration_um) / (
            self.km_suppressor_um + self.suppressor_trna_concentration_um
        )
        # Termination velocity modulated by flanking nucleotide context
        v_term = (
            (self.kcat_erf_s1 * self.erf_complex_concentration_um)
            / (self.km_erf_um + self.erf_complex_concentration_um)
        ) * self.flanking_context_factor

        # Near-cognate basal decoding velocity (very low basal readthrough ~ 0.5-1%)
        v_near = 0.05

        v_total = v_sup + v_term + v_near
        if v_total <= 0.0:
            return 0.0
        return max(0.0, min(1.0, v_sup / v_total))


@dataclass(frozen=True)
class NMDEvasionParameters:
    """Nonsense-mediated decay and SMG1 phosphorylation kinetics."""
    variant: str = "W1282X"
    baseline_mrna_pct_wt: float = 33.3   # Depressed baseline (33.3% for W1282X; 14.5% for G542X)
    k_nmd_max_h1: float = 0.80           # Maximum NMD degradation rate (1/h)
    k_deg_basal_h1: float = 0.20         # Basal mRNA turnover rate (1/h)
    smg1_inhibition_pct: float = 85.0    # SMG1 small-molecule inhibition percentage

    def compute_steady_state_mrna(self, readthrough_or_edited_fraction: float) -> float:
        """Compute steady-state CFTR mRNA relative to wild-type (%) under editing/readthrough."""
        # Unreadthrough / unedited PTC transcripts undergo UPF1-SMG1 phosphorylation
        unrescued = max(0.0, 1.0 - readthrough_or_edited_fraction)
        smg1_activity = max(0.0, 1.0 - (self.smg1_inhibition_pct / 100.0))

        # Effective NMD degradation rate
        k_nmd_eff = self.k_nmd_max_h1 * unrescued * smg1_activity
        k_total = self.k_deg_basal_h1 + k_nmd_eff

        # Steady-state mRNA scales inversely with total turnover: mRNA_ss ~ k_tx / k_total
        # For WT, k_total_wt = k_deg_basal_h1 (no NMD)
        relative_mrna_pct = (self.k_deg_basal_h1 / k_total) * 100.0

        # Anchor to known mutant basal bounds when uninhibited
        min_baseline = self.baseline_mrna_pct_wt
        return max(min_baseline, min(100.0, relative_mrna_pct))


def evaluate_ptc_genetic_rescue(
    variant: str = "W1282X",
    prime_editing_enabled: bool = True,
    has_tevopreq1: bool = True,
    has_mlh1dn: bool = True,
    suppressor_trna_enabled: bool = True,
    smg1_inhibitor_enabled: bool = True,
) -> Dict[str, Any]:
    """Evaluate comprehensive genetic and transcriptional rescue for W1282X or G542X."""
    phi_mmr = 0.08 if has_mlh1dn else 1.0
    pe_params = PrimeEditingParameters(
        variant=variant,
        delta_delta_g_flap_kcal_mol=-3.2 if variant == "W1282X" else -2.9,
        has_tevopreq1=has_tevopreq1,
        phi_mmr_evasion=phi_mmr,
    )
    pe_efficiency = pe_params.compute_prime_editing_efficiency() if prime_editing_enabled else 0.0

    trna_conc = 1.5 if suppressor_trna_enabled else 0.0
    dec_params = RibosomalDecodingParameters(
        variant=variant,
        suppressor_trna_concentration_um=trna_conc,
    )
    readthrough_efficiency = dec_params.compute_readthrough_efficiency() if suppressor_trna_enabled else 0.005

    # Combined full-length transcript/protein rescue factor:
    # Alleles edited permanently (pe_efficiency) produce WT CFTR;
    # Remaining unedited alleles (1 - pe_efficiency) decoded by suppressor tRNA
    combined_coding_rescue = pe_efficiency + (1.0 - pe_efficiency) * readthrough_efficiency

    baseline_mrna = 33.3 if variant == "W1282X" else 14.5
    smg1_inhibition = 85.0 if smg1_inhibitor_enabled else 0.0
    nmd_params = NMDEvasionParameters(
        variant=variant,
        baseline_mrna_pct_wt=baseline_mrna,
        smg1_inhibition_pct=smg1_inhibition,
    )
    mrna_level_pct = nmd_params.compute_steady_state_mrna(combined_coding_rescue)

    # Post-translational Band C maturation factor:
    # W1282X full-length protein retains NBD2 and matures with ~85% WT efficiency
    # G542X restored with Glycine (GGA) is true wild-type sequence -> 95% maturation
    maturation_efficiency = 0.85 if variant == "W1282X" else 0.95

    # Transepithelial short-circuit current Delta Isc (% of wild-type)
    # Delta Isc = (mRNA / 100) * combined_coding_rescue * maturation_efficiency * 100
    delta_isc_pct_wt = (mrna_level_pct / 100.0) * combined_coding_rescue * maturation_efficiency * 100.0

    # Clinical therapeutic threshold: >= 10.0% of wild-type Delta Isc
    clinical_threshold_met = delta_isc_pct_wt >= 10.0

    return {
        "variant": variant,
        "prime_editing_efficiency": pe_efficiency,
        "suppressor_trna_readthrough": readthrough_efficiency,
        "combined_coding_rescue": combined_coding_rescue,
        "steady_state_mrna_pct_wt": mrna_level_pct,
        "restored_delta_isc_pct_wt": delta_isc_pct_wt,
        "clinical_therapeutic_threshold_met": clinical_threshold_met,
        "curative_verdict": "THERAPEUTIC" if clinical_threshold_met else "SUB_THERAPEUTIC",
    }
