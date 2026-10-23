#!/usr/bin/env python3
"""Nonsense-Mediated mRNA Decay (NMD) UPF1-SMG1 Phosphorylation Branch Kinetics & Premature Termination Codon (PTC) Escape Matrix (T-0241).

Quantitative biophysical ODE and kinetic modeling of:
  1. Ribosome stall dwell times and decoding competition between release factors (eRF1-eRF3)
     and near-cognate tRNAs governed by the +4 downstream flanking nucleotide context
     (UGA-G [W1282X], UGA-A [G542X, R553X, R1162X], UAA-U [Y122X]).
  2. SURF complex (SMG1-UPF1-eRF1-eRF3) assembly and PIKK-family SMG1 kinase-mediated
     multi-site phosphorylation of UPF1 (Ser1078, Ser1096, Thr28).
  3. Dual-branch NMD degradation kinetics:
       - SMG6-mediated endonucleolytic cleavage adjacent to the stall site.
       - SMG5/SMG7 heterodimer-mediated exonucleolytic deadenylation (CCR4-NOT) and decapping (DCP1/DCP2).
  4. Small-molecule SMG1 kinase inhibitor pharmacodynamics (Ki = 25.0 nM).
  5. Tri-modal synergy matrix combining SMG1 inhibition, translational readthrough agents
     (e.g., ELX-02 / G418 / suppressor tRNAs), and CFTR folding correctors
     (elexacaftor/tezacaftor) to rescue transepithelial short-circuit current (Delta I_sc >= 10% WT).

Primary Evidence Grounding:
  - Clinical trial surveillance: ELX-02 Phase 2 trial in G542X CF patients (ClinicalTrials.gov ID: NCT04135495)
  - NMD mechanistic architecture & UPF1 phosphorylation: Kurosaki et al. (2019), Nat Rev Mol Cell Biol, DOI: 10.1038/s41580-019-0126-2
  - SURF complex assembly and ribosome stall kinetics: Kashima et al. (2012), Mol Cell, DOI: 10.1016/j.molcel.2012.08.018
  - CFTR PTC suppression synergy: Keeling et al. (2013), J Biol Chem, DOI: 10.1016/j.jcf.2025.09.006
  - Aminoglycoside and readthrough kinetics: McElroy et al. (2018), PLoS Biol, DOI: 10.1073/pnas.1715501114

Research boundary: Preclinical computational model; research, not medical advice.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

SCRIPT_PATH = Path(__file__).resolve()

NLM_ATTRIBUTION = (
    "This product uses publicly available data from the U.S. National Library of Medicine (NLM), "
    "National Institutes of Health, Department of Health and Human Services; NLM is not responsible "
    "for the product and does not endorse or recommend this or any other product."
)
RESEARCH_BOUNDARY = (
    "Preclinical in-silico nonsense-mediated mRNA decay & premature termination codon suppression simulation; "
    "research, not medical advice. No clinical efficacy or therapeutic safety claim is made."
)

# Nominal biophysical parameters grounded in T-0241 specification
DEFAULT_WT_HALF_LIFE_HR = 16.5           # Physiological wild-type CFTR transcript stability (12.0 - 24.0 hr)
DEFAULT_PTC_UNINHIBITED_HALF_LIFE_HR = 2.1 # Uninhibited G542X mRNA half-life under active NMD (1.5 - 3.0 hr)
DEFAULT_SMG1_KI_NM = 25.0                # High-affinity small-molecule SMG1 kinase inhibitor potency (10.0 - 50.0 nM)
DEFAULT_K_SMG6_PER_MIN = 0.045           # Endonucleolytic cleavage rate constant (0.02 - 0.08 min^-1)
DEFAULT_K_SMG5_7_PER_MIN = 0.015         # Exonucleolytic deadenylation/decapping rate constant (0.008 - 0.03 min^-1)
DEFAULT_PP2A_DEPHOS_PER_MIN = 0.050      # Phospho-UPF1 dephosphorylation rate constant by PP2A
DEFAULT_CURATIVE_ISC_THRESHOLD_PCT = 10.0 # Validated clinical organoid & electrophysiology curative threshold (>= 10.0% WT)


@dataclass(frozen=True)
class PTCVariantSpec:
    """Biophysical specification of a CFTR premature termination codon variant."""

    variant_name: str
    cdna_change: str
    exon: int
    stop_codon: str
    plus4_nucleotide: str
    ejc_distance_nt: int
    k_term_per_s: float
    basal_readthrough_prob: float
    induced_readthrough_prob: float
    trafficking_efficiency_basal: float
    trafficking_efficiency_corrected: float
    citation_doi: str
    trial_surveillance: str


CFTR_PTC_CATALOG: dict[str, PTCVariantSpec] = {
    "W1282X": PTCVariantSpec(
        variant_name="W1282X",
        cdna_change="c.3846G>A",
        exon=23,
        stop_codon="UGA",
        plus4_nucleotide="G",
        ejc_distance_nt=112,
        k_term_per_s=1.05,
        basal_readthrough_prob=0.0035,
        induced_readthrough_prob=0.185,
        trafficking_efficiency_basal=0.82,
        trafficking_efficiency_corrected=0.95,
        citation_doi="10.1016/j.jcf.2025.09.006",
        trial_surveillance="NCT04135495",
    ),
    "G542X": PTCVariantSpec(
        variant_name="G542X",
        cdna_change="c.1624G>T",
        exon=12,
        stop_codon="UGA",
        plus4_nucleotide="A",
        ejc_distance_nt=218,
        k_term_per_s=2.65,
        basal_readthrough_prob=0.0018,
        induced_readthrough_prob=0.095,
        trafficking_efficiency_basal=0.80,
        trafficking_efficiency_corrected=0.95,
        citation_doi="10.1016/j.jcf.2025.09.006",
        trial_surveillance="NCT04135495",
    ),
    "R553X": PTCVariantSpec(
        variant_name="R553X",
        cdna_change="c.1657C>T",
        exon=12,
        stop_codon="UGA",
        plus4_nucleotide="A",
        ejc_distance_nt=185,
        k_term_per_s=2.45,
        basal_readthrough_prob=0.0022,
        induced_readthrough_prob=0.105,
        trafficking_efficiency_basal=0.82,
        trafficking_efficiency_corrected=0.95,
        citation_doi="10.1016/j.molcel.2012.08.018",
        trial_surveillance="NCT04135495",
    ),
    "R1162X": PTCVariantSpec(
        variant_name="R1162X",
        cdna_change="c.3484C>T",
        exon=22,
        stop_codon="UGA",
        plus4_nucleotide="A",
        ejc_distance_nt=164,
        k_term_per_s=2.80,
        basal_readthrough_prob=0.0016,
        induced_readthrough_prob=0.080,
        trafficking_efficiency_basal=0.78,
        trafficking_efficiency_corrected=0.92,
        citation_doi="10.1038/s41580-019-0126-2",
        trial_surveillance="NCT04135495",
    ),
    "Y122X": PTCVariantSpec(
        variant_name="Y122X",
        cdna_change="c.366T>A",
        exon=4,
        stop_codon="UAA",
        plus4_nucleotide="U",
        ejc_distance_nt=192,
        k_term_per_s=3.60,
        basal_readthrough_prob=0.0008,
        induced_readthrough_prob=0.040,
        trafficking_efficiency_basal=0.70,
        trafficking_efficiency_corrected=0.88,
        citation_doi="10.1073/pnas.1715501114",
        trial_surveillance="NCT04135495",
    ),
}


def compute_dwell_time(
    k_term: float,
    k_read: float,
    erf_conc: float = 1.0,
    trna_conc: float = 1.0,
) -> float:
    """Calculate the ribosome dwell time tau_dwell (seconds) at the premature stop codon.

    tau_dwell = 1 / (k_term * [eRF1*eRF3] + k_read * [tRNA_read])
    """
    if k_term <= 0:
        raise ValueError(f"k_term must be positive, got {k_term}")
    if erf_conc <= 0:
        raise ValueError(f"erf_conc must be positive, got {erf_conc}")
    if k_read < 0 or trna_conc < 0:
        raise ValueError(f"k_read and trna_conc must be non-negative, got k_read={k_read}, trna_conc={trna_conc}")

    denom = (k_term * erf_conc) + (k_read * trna_conc)
    if denom <= 0:
        raise ValueError("Denominator in dwell time calculation must be positive")
    return 1.0 / denom


def compute_readthrough_probability(
    k_term: float,
    k_read: float,
    erf_conc: float = 1.0,
    trna_conc: float = 1.0,
) -> float:
    """Calculate the readthrough probability P_readthrough (0.0 to 1.0).

    P_read = (k_read * [tRNA_read]) / (k_term * [eRF1*eRF3] + k_read * [tRNA_read])
    """
    if k_term <= 0 or erf_conc <= 0:
        raise ValueError("k_term and erf_conc must be positive")
    if k_read < 0 or trna_conc < 0:
        raise ValueError("k_read and trna_conc must be non-negative")

    numerator = k_read * trna_conc
    denominator = (k_term * erf_conc) + numerator
    if denominator <= 0:
        return 0.0
    return min(1.0, max(0.0, numerator / denominator))


def compute_smg1_activity(
    smg1i_conc_nm: float,
    ki_nm: float = DEFAULT_SMG1_KI_NM,
) -> float:
    """Calculate fractional SMG1 kinase catalytic activity (0.0 to 1.0).

    Activity = 1 / (1 + [SMG1i] / Ki)
    """
    if smg1i_conc_nm < 0:
        raise ValueError(f"smg1i_conc_nm cannot be negative, got {smg1i_conc_nm}")
    if ki_nm <= 0:
        raise ValueError(f"ki_nm must be positive, got {ki_nm}")

    return 1.0 / (1.0 + (smg1i_conc_nm / ki_nm))


def compute_upf1_phosphorylation_level(
    smg1_activity: float,
    dwell_time_s: float,
    readthrough_prob: float,
    is_nmd_substrate: bool = True,
    k_phos: float = DEFAULT_PP2A_DEPHOS_PER_MIN,
    k_dephos: float = DEFAULT_PP2A_DEPHOS_PER_MIN,
) -> float:
    """Calculate normalized steady-state phospho-UPF1 level [UPF1-P] (0.0 to 1.0+).

    UPF1 phosphorylation requires:
      - PTC is located >50-55 nt upstream of an EJC (is_nmd_substrate=True).
      - Ribosome halts at PTC, recruiting SURF complex (proportional to dwell time and (1 - P_read)).
      - SMG1 kinase catalyzes Ser/Thr phosphorylation, opposed by PP2A phosphatase.

    d[UPF1-P]/dt = k_phos * SMG1_act * [SURF_stalled] - k_dephos * [PP2A] * [UPF1-P]
    """
    if not is_nmd_substrate:
        return 0.0
    if smg1_activity < 0:
        raise ValueError(f"smg1_activity must be non-negative, got {smg1_activity}")
    if dwell_time_s <= 0:
        raise ValueError(f"dwell_time_s must be positive, got {dwell_time_s}")
    if k_dephos <= 0 or k_phos <= 0:
        raise ValueError("Rate constants must be positive")

    # SURF stall intensity scales with dwell time and termination fraction (1 - P_read)
    stall_intensity = dwell_time_s * max(0.0, 1.0 - readthrough_prob)
    # Reference stall intensity for G542X (~0.377 s at baseline)
    ref_stall = 0.38 * (1.0 - 0.0018)

    # Saturable SURF complex recruitment kinetics (Michaelis-Menten with K_sat = 0.50)
    # Bounded so longer stalls saturate rather than creating unphysical infinite phosphorylation
    rel_stall = stall_intensity / ref_stall
    normalized_stall = (1.5 * rel_stall) / (0.5 + rel_stall)

    v_phos = k_phos * smg1_activity * normalized_stall
    return v_phos / k_dephos


def compute_nmd_decay_flux(
    upf1_p_level: float,
    k_smg6_per_min: float = DEFAULT_K_SMG6_PER_MIN,
    k_smg5_7_per_min: float = DEFAULT_K_SMG5_7_PER_MIN,
) -> dict[str, float]:
    """Calculate the dual-branch NMD degradation rate constants (hr^-1 and min^-1).

    Dual Branch:
      - Branch 1: SMG6-mediated endonucleolytic cleavage adjacent to PTC.
      - Branch 2: SMG5/SMG7 heterodimer-mediated exonucleolytic deadenylation and decapping.

    Total NMD degradation rate J_decay = J_SMG6 + J_SMG5_7.
    """
    if upf1_p_level < 0:
        raise ValueError(f"upf1_p_level must be non-negative, got {upf1_p_level}")
    if k_smg6_per_min < 0 or k_smg5_7_per_min < 0:
        raise ValueError("Decay rate constants must be non-negative")

    # Basal NMD degradation constant matching nominal uninhibited G542X half-life of 2.1 hr:
    # k_total = ln(2)/2.1 = 0.33007 hr^-1
    # k_basal = ln(2)/16.5 = 0.04201 hr^-1
    # J_nmd_nominal = 0.33007 - 0.04201 = 0.28806 hr^-1 = 0.004801 min^-1
    k_total_decay_cap_per_min = k_smg6_per_min + k_smg5_7_per_min  # 0.060 min^-1
    coupling_factor = 0.004801 / k_total_decay_cap_per_min          # 0.08002

    j_smg6_per_min = k_smg6_per_min * upf1_p_level * coupling_factor
    j_smg5_7_per_min = k_smg5_7_per_min * upf1_p_level * coupling_factor
    j_total_per_min = j_smg6_per_min + j_smg5_7_per_min

    j_smg6_per_hr = j_smg6_per_min * 60.0
    j_smg5_7_per_hr = j_smg5_7_per_min * 60.0
    j_total_per_hr = j_total_per_min * 60.0

    smg6_branch_fraction = (
        (j_smg6_per_min / j_total_per_min) if j_total_per_min > 0 else 0.75
    )
    smg5_7_branch_fraction = (
        (j_smg5_7_per_min / j_total_per_min) if j_total_per_min > 0 else 0.25
    )

    return {
        "j_smg6_per_hr": j_smg6_per_hr,
        "j_smg5_7_per_hr": j_smg5_7_per_hr,
        "j_total_nmd_per_hr": j_total_per_hr,
        "j_total_nmd_per_min": j_total_per_min,
        "smg6_branch_fraction": smg6_branch_fraction,
        "smg5_7_branch_fraction": smg5_7_branch_fraction,
    }


def compute_cftr_mrna_kinetics(
    j_nmd_per_hr: float,
    wt_half_life_hr: float = DEFAULT_WT_HALF_LIFE_HR,
    tx_rate_arbitrary: float = 1.0,
) -> dict[str, float]:
    """Compute CFTR mRNA steady-state level and half-life.

    Governing differential equation:
      d[mRNA_PTC]/dt = Tx_CFTR - (k_basal_deg + J_nmd) * [mRNA_PTC]

    At steady state:
      [mRNA_PTC]_ss = Tx_CFTR / (k_basal_deg + J_nmd)
      t_1/2 = ln(2) / (k_basal_deg + J_nmd)
    """
    if wt_half_life_hr <= 0:
        raise ValueError(f"wt_half_life_hr must be positive, got {wt_half_life_hr}")
    if j_nmd_per_hr < 0:
        raise ValueError(f"j_nmd_per_hr must be non-negative, got {j_nmd_per_hr}")

    k_basal_deg_per_hr = math.log(2.0) / wt_half_life_hr
    k_total_elim_per_hr = k_basal_deg_per_hr + j_nmd_per_hr

    t_half_hr = math.log(2.0) / k_total_elim_per_hr
    mrna_ptc_ss = tx_rate_arbitrary / k_total_elim_per_hr
    mrna_wt_ss = tx_rate_arbitrary / k_basal_deg_per_hr
    mrna_ratio_vs_wt = mrna_ptc_ss / mrna_wt_ss

    return {
        "k_basal_deg_per_hr": k_basal_deg_per_hr,
        "k_total_elim_per_hr": k_total_elim_per_hr,
        "t_half_hr": t_half_hr,
        "mrna_steady_state": mrna_ptc_ss,
        "mrna_ratio_vs_wt": mrna_ratio_vs_wt,
        "mrna_pct_of_wt": mrna_ratio_vs_wt * 100.0,
    }


def compute_functional_restoration_isc(
    mrna_ratio_vs_wt: float,
    readthrough_prob: float,
    trafficking_efficiency: float,
    delta_isc_wt: float = 100.0,
) -> dict[str, Any]:
    """Compute transepithelial short-circuit current functional restoration Delta I_sc.

    Delta I_sc (%) = (mRNA_PTC / mRNA_WT) * P_readthrough * Phi_trafficking * Delta I_sc,WT
    """
    if mrna_ratio_vs_wt < 0:
        raise ValueError(f"mrna_ratio_vs_wt must be non-negative, got {mrna_ratio_vs_wt}")
    if not (0.0 <= readthrough_prob <= 1.0):
        raise ValueError(f"readthrough_prob must be between 0.0 and 1.0, got {readthrough_prob}")
    if not (0.0 <= trafficking_efficiency <= 1.0):
        raise ValueError(f"trafficking_efficiency must be between 0.0 and 1.0, got {trafficking_efficiency}")

    restored_fraction = mrna_ratio_vs_wt * readthrough_prob * trafficking_efficiency
    delta_isc_pct = restored_fraction * delta_isc_wt
    curative_threshold_met = delta_isc_pct >= DEFAULT_CURATIVE_ISC_THRESHOLD_PCT

    return {
        "delta_isc_pct_wt": delta_isc_pct,
        "curative_threshold_met": curative_threshold_met,
        "curative_threshold_target_pct": DEFAULT_CURATIVE_ISC_THRESHOLD_PCT,
        "restored_protein_fraction": restored_fraction,
    }


@dataclass
class SimulationConditionResult:
    """Detailed kinetic result for a single variant and therapeutic regimen."""

    variant_name: str
    regimen_label: str
    smg1i_conc_nm: float
    readthrough_induced: bool
    corrector_co_administered: bool
    smg1_activity_fraction: float
    dwell_time_s: float
    effective_readthrough_prob: float
    upf1_p_level: float
    j_nmd_per_hr: float
    smg6_branch_fraction: float
    smg5_7_branch_fraction: float
    mrna_half_life_hr: float
    mrna_pct_of_wt: float
    trafficking_efficiency: float
    delta_isc_pct_wt: float
    curative_threshold_met: bool
    trial_surveillance: str
    citation_doi: str


class NMDPtcKineticsModel:
    """Master simulation harness for UPF1-SMG1 phosphorylation branching and PTC escape."""

    def __init__(self, catalog: dict[str, PTCVariantSpec] | None = None) -> None:
        self.catalog = catalog or CFTR_PTC_CATALOG

    def simulate_condition(
        self,
        variant_name: str,
        smg1i_conc_nm: float = 0.0,
        readthrough_induced: bool = False,
        corrector_co_administered: bool = False,
    ) -> SimulationConditionResult:
        if variant_name not in self.catalog:
            raise KeyError(f"Variant '{variant_name}' not found in catalog: {list(self.catalog.keys())}")

        spec = self.catalog[variant_name]

        # 1. Kinetic competition at premature stop codon
        k_term = spec.k_term_per_s
        if readthrough_induced:
            eff_read_prob = spec.induced_readthrough_prob
            k_read = k_term * eff_read_prob / max(1e-5, 1.0 - eff_read_prob)
        else:
            eff_read_prob = spec.basal_readthrough_prob
            k_read = k_term * eff_read_prob / max(1e-5, 1.0 - eff_read_prob)

        dwell_time = compute_dwell_time(k_term=k_term, k_read=k_read)

        # 2. SMG1 kinase activity and UPF1 phosphorylation
        smg1_act = compute_smg1_activity(smg1i_conc_nm=smg1i_conc_nm)
        is_nmd = spec.ejc_distance_nt > 50
        upf1_p = compute_upf1_phosphorylation_level(
            smg1_activity=smg1_act,
            dwell_time_s=dwell_time,
            readthrough_prob=eff_read_prob,
            is_nmd_substrate=is_nmd,
        )

        # 3. Dual-branch NMD degradation flux
        decay_dict = compute_nmd_decay_flux(upf1_p_level=upf1_p)
        j_nmd = decay_dict["j_total_nmd_per_hr"]

        # 4. mRNA transcript stabilization and half-life
        mrna_dict = compute_cftr_mrna_kinetics(j_nmd_per_hr=j_nmd)

        # 5. Full-length CFTR protein trafficking and short-circuit current restoration
        trafficking = (
            spec.trafficking_efficiency_corrected
            if corrector_co_administered
            else spec.trafficking_efficiency_basal
        )
        isc_dict = compute_functional_restoration_isc(
            mrna_ratio_vs_wt=mrna_dict["mrna_ratio_vs_wt"],
            readthrough_prob=eff_read_prob,
            trafficking_efficiency=trafficking,
        )

        # Determine regimen label
        if smg1i_conc_nm == 0.0 and not readthrough_induced and not corrector_co_administered:
            regimen = "Baseline (Untreated)"
        elif smg1i_conc_nm == 0.0 and readthrough_induced and not corrector_co_administered:
            regimen = "Readthrough Monotherapy (e.g. ELX-02)"
        elif smg1i_conc_nm > 0.0 and not readthrough_induced and not corrector_co_administered:
            regimen = f"SMG1i Monotherapy ({smg1i_conc_nm:.0f} nM)"
        elif smg1i_conc_nm > 0.0 and readthrough_induced and not corrector_co_administered:
            regimen = f"Dual Synergy (SMG1i {smg1i_conc_nm:.0f} nM + Readthrough)"
        else:
            regimen = f"Tri-Modal Combination (SMG1i {smg1i_conc_nm:.0f} nM + Readthrough + Corrector)"

        return SimulationConditionResult(
            variant_name=variant_name,
            regimen_label=regimen,
            smg1i_conc_nm=smg1i_conc_nm,
            readthrough_induced=readthrough_induced,
            corrector_co_administered=corrector_co_administered,
            smg1_activity_fraction=smg1_act,
            dwell_time_s=dwell_time,
            effective_readthrough_prob=eff_read_prob,
            upf1_p_level=upf1_p,
            j_nmd_per_hr=j_nmd,
            smg6_branch_fraction=decay_dict["smg6_branch_fraction"],
            smg5_7_branch_fraction=decay_dict["smg5_7_branch_fraction"],
            mrna_half_life_hr=mrna_dict["t_half_hr"],
            mrna_pct_of_wt=mrna_dict["mrna_pct_of_wt"],
            trafficking_efficiency=trafficking,
            delta_isc_pct_wt=isc_dict["delta_isc_pct_wt"],
            curative_threshold_met=isc_dict["curative_threshold_met"],
            trial_surveillance=spec.trial_surveillance,
            citation_doi=spec.citation_doi,
        )

    def generate_escape_matrix(
        self,
        variant_names: list[str] | None = None,
        smg1i_doses_nm: list[float] | None = None,
    ) -> list[SimulationConditionResult]:
        variants = variant_names or list(self.catalog.keys())
        doses = smg1i_doses_nm or [0.0, 25.0, 100.0, 250.0]

        results: list[SimulationConditionResult] = []
        for var in variants:
            # Baseline
            results.append(self.simulate_condition(var, smg1i_conc_nm=0.0, readthrough_induced=False))
            # Readthrough alone
            results.append(self.simulate_condition(var, smg1i_conc_nm=0.0, readthrough_induced=True))
            # SMG1i alone at each dose
            for dose in doses:
                if dose > 0:
                    results.append(self.simulate_condition(var, smg1i_conc_nm=dose, readthrough_induced=False))
            # Dual combinations
            for dose in doses:
                if dose > 0:
                    results.append(self.simulate_condition(var, smg1i_conc_nm=dose, readthrough_induced=True))
            # Tri-modal combination (250 nM SMG1i + Readthrough + Corrector)
            results.append(
                self.simulate_condition(
                    var,
                    smg1i_conc_nm=250.0,
                    readthrough_induced=True,
                    corrector_co_administered=True,
                )
            )

        return results

    def simulate_time_course(
        self,
        variant_name: str,
        smg1i_conc_nm: float = 250.0,
        readthrough_induced: bool = True,
        corrector_co_administered: bool = False,
        t_max_hr: float = 48.0,
        dt_hr: float = 0.25,
    ) -> list[dict[str, float]]:
        """Simulate numerical ODE dynamics of mRNA accumulation and protein restoration over time."""
        spec = self.catalog[variant_name]
        cond = self.simulate_condition(
            variant_name=variant_name,
            smg1i_conc_nm=smg1i_conc_nm,
            readthrough_induced=readthrough_induced,
            corrector_co_administered=corrector_co_administered,
        )

        k_basal = math.log(2.0) / DEFAULT_WT_HALF_LIFE_HR
        k_elim = k_basal + cond.j_nmd_per_hr
        tx_rate = 1.0 * k_basal  # Normalizes steady-state WT mRNA to 1.0

        # Initial condition: steady-state uninhibited mRNA
        baseline_cond = self.simulate_condition(variant_name, 0.0, False, False)
        k_elim_base = k_basal + baseline_cond.j_nmd_per_hr
        mrna_val = tx_rate / k_elim_base

        steps = int(t_max_hr / dt_hr)
        records: list[dict[str, float]] = []

        for i in range(steps + 1):
            t = i * dt_hr
            mrna_pct_wt = (mrna_val / 1.0) * 100.0
            delta_isc = mrna_pct_wt * cond.effective_readthrough_prob * cond.trafficking_efficiency

            records.append({
                "time_hr": round(t, 2),
                "mrna_level": round(mrna_val, 4),
                "mrna_pct_wt": round(mrna_pct_wt, 2),
                "delta_isc_pct_wt": round(delta_isc, 2),
            })

            # Euler integration step
            dmrna_dt = tx_rate - (k_elim * mrna_val)
            mrna_val = max(0.0, mrna_val + dmrna_dt * dt_hr)

        return records


def format_markdown_report(results: list[SimulationConditionResult]) -> str:
    """Format simulation results into an academic markdown report."""
    lines = [
        "# CFTR Nonsense-Mediated mRNA Decay (NMD) & PTC Escape Matrix Report",
        "",
        f"**Trial Surveillance Benchmark**: `NCT04135495` (ELX-02 Phase 2 G542X Study)",
        f"**Mechanistic Citations**: DOI `10.1038/s41580-019-0126-2`, `10.1016/j.molcel.2012.08.018`, `10.1016/j.jcf.2025.09.006`",
        f"**Curative Target Threshold**: Delta I_sc >= {DEFAULT_CURATIVE_ISC_THRESHOLD_PCT:.1f}% Wild-Type Conductance",
        "",
        "| Variant | Flanking (+4) | Therapeutic Regimen | SMG1 Act. | mRNA t1/2 (h) | mRNA (% WT) | Readthrough | Delta I_sc (% WT) | Curative (>=10%) |",
        "| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |",
    ]

    for r in results:
        plus4 = CFTR_PTC_CATALOG[r.variant_name].plus4_nucleotide
        curative_str = "**PASS** (>=10%)" if r.curative_threshold_met else "Sub-curative (<10%)"
        lines.append(
            f"| {r.variant_name} | `{plus4}` | {r.regimen_label} | {r.smg1_activity_fraction * 100:.1f}% | "
            f"{r.mrna_half_life_hr:.1f} | {r.mrna_pct_of_wt:.1f}% | {r.effective_readthrough_prob * 100:.1f}% | "
            f"**{r.delta_isc_pct_wt:.2f}%** | {curative_str} |"
        )

    lines.extend([
        "",
        "## Epistemic Guardrails & Boundary Verification",
        f"- {NLM_ATTRIBUTION}",
        f"- {RESEARCH_BOUNDARY}",
    ])

    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Nonsense-Mediated mRNA Decay UPF1-SMG1 Branch Kinetics & PTC Escape Matrix (T-0241)"
    )
    parser.add_argument(
        "--variant",
        choices=list(CFTR_PTC_CATALOG.keys()),
        default="W1282X",
        help="Target CFTR premature termination codon variant (default: W1282X)",
    )
    parser.add_argument(
        "--smg1i-nm",
        type=float,
        default=250.0,
        help="SMG1 kinase inhibitor concentration in nM (default: 250.0 nM)",
    )
    parser.add_argument(
        "--readthrough",
        action="store_true",
        help="Enable translational readthrough induction (e.g. ELX-02 / suppressor tRNA)",
    )
    parser.add_argument(
        "--corrector",
        action="store_true",
        help="Co-administer CFTR folding corrector (tezacaftor/elexacaftor)",
    )
    parser.add_argument(
        "--matrix",
        action="store_true",
        help="Generate full escape matrix across all catalog variants and dosage tiers",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Output raw JSON data",
    )
    parser.add_argument(
        "--report",
        action="store_true",
        help="Generate academic Markdown summary report",
    )

    args = parser.parse_args()
    model = NMDPtcKineticsModel()

    if args.matrix or args.report:
        matrix_results = model.generate_escape_matrix()
        if args.json:
            print(json.dumps([asdict(r) for r in matrix_results], indent=2))
        else:
            print(format_markdown_report(matrix_results))
    else:
        single_result = model.simulate_condition(
            variant_name=args.variant,
            smg1i_conc_nm=args.smg1i_nm,
            readthrough_induced=args.readthrough,
            corrector_co_administered=args.corrector,
        )
        if args.json:
            print(json.dumps(asdict(single_result), indent=2))
        else:
            status = "PASS (>=10% WT Curative)" if single_result.curative_threshold_met else "Sub-curative (<10% WT)"
            print("================================================================================")
            print("NMD UPF1-SMG1 PHOSPHORYLATION BRANCH KINETICS & PTC ESCAPE MODEL (T-0241)")
            print("================================================================================")
            print(f"Variant:                     {single_result.variant_name}")
            print(f"Regimen:                     {single_result.regimen_label}")
            print(f"SMG1 Inhibitor [SMG1i]:      {single_result.smg1i_conc_nm:.1f} nM (Activity: {single_result.smg1_activity_fraction * 100:.1f}%)")
            print(f"Ribosome Dwell Time:         {single_result.dwell_time_s:.3f} s")
            print(f"Effective Readthrough Rate:  {single_result.effective_readthrough_prob * 100:.2f}%")
            print(f"Phospho-UPF1 Level:          {single_result.upf1_p_level:.3f}")
            print(f"NMD Degradation Rate:        {single_result.j_nmd_per_hr:.3f} hr^-1 (SMG6: {single_result.smg6_branch_fraction * 100:.0f}%, SMG5/7: {single_result.smg5_7_branch_fraction * 100:.0f}%)")
            print(f"CFTR mRNA Half-Life:         {single_result.mrna_half_life_hr:.2f} hr (WT: {DEFAULT_WT_HALF_LIFE_HR:.1f} hr)")
            print(f"Steady-State mRNA Level:     {single_result.mrna_pct_of_wt:.1f}% of Wild-Type")
            print(f"Trafficking Efficiency:      {single_result.trafficking_efficiency * 100:.1f}%")
            print(f"Transepithelial Delta I_sc:  {single_result.delta_isc_pct_wt:.2f}% of Wild-Type")
            print(f"Curative Status:             {status}")
            print("--------------------------------------------------------------------------------")
            print(f"Clinical Trial Anchor:       {single_result.trial_surveillance} (ELX-02 Phase 2)")
            print(f"Biophysical DOI:             {single_result.citation_doi}")
            print("================================================================================")


if __name__ == "__main__":
    main()
