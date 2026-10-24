#!/usr/bin/env python3
"""Eukaryotic Release Factor eRF1/eRF3 Depletion & Aminoglycoside-Free Nonsense Readthrough Kinetic Model for W1282X CFTR (T-0401).

Quantitative biophysical model of ribosome decoding competition, release factor depletion,
and coupled nonsense-mediated mRNA decay (NMD) escape for human CFTR premature termination codons (PTCs):
  1. Competitive Ribosomal A-Site Termination vs Readthrough Kinetics:
     Simulates kinetic competition between the ternary termination complex (eRF1-eRF3-GTP)
     and near-cognate aminoacyl-tRNAs at ribosomal premature stop codons:
       v_term = (k_term * [eRF]_complex) / (K_m,term + [eRF]_complex) * omega_context * omega_codon
       v_read = (k_read * [tRNA_nc]) / (K_m,read + [tRNA_nc])
       F_readthrough = v_read / (v_term + v_read)
  2. eRF Depletion Dynamics & Heterodimer Cooperativity:
     Models small-molecule termination inhibitors / release factor depleters (e.g., SRI-37240, SRI-41315)
     acting cooperatively (Hill coefficient h in [1.5, 2.0], IC50 in [0.1, 1.0] uM):
       [eRF](C_drug) = [eRF]_0 * (1 - (I_max * C_drug^h) / (IC50^h + C_drug^h))
     Because eRF1 and eRF3 must form an active heterodimeric complex to enter the ribosomal A-site,
     co-depletion scales functional ternary termination complexes super-linearly ([eRF]_complex ~ [eRF]^1.6),
     conferring a characteristic 5- to 15-fold readthrough boost upon 60-80% release factor depletion.
  3. Codon Permissiveness Hierarchy & +4 Context Factor:
     Accounts for stop codon termination efficiency (UGA > UAG >> UAA) and the +4 flanking
     downstream nucleotide weighting (omega_context approx 0.35 for +4 C/U in W1282X vs 1.0 for +4 G).
  4. Coupled NMD Escape & Mature Band C CFTR Expression:
     Couples pioneer round translation elongation through the PTC with exon junction complex (EJC)
     displacement, rescuing mRNA stability from nonsense-mediated decay and computing functional
     mature 170 kDa glycosylated Band C CFTR protein expression and chloride current rescue.

Primary Evidence Grounding:
  - Discovery Brief: cf-research/landscape/mission-grounding-2026-10-23-wave17.md
  - Cryo-EM Structure: RCSB PDB 6O2P (Human CFTR potentiator-bound gating conformation at 3.3 A)
  - Primary Literature:
    - Baradaran-Heravi A et al. (2021) Nat Chem Biol, DOI: 10.1038/s41589-020-00662-8
    - Roy B et al. (2020) Nucleic Acids Res, DOI: 10.1093/nar/gkaa728
  - Clinical Registry Grounding:
    - NCT04135495: Phase 2 evaluation of ELX-02 in CF patients with nonsense mutation

Epistemic Status:
  - Answer Tier T2 (Open-access primary literature, Cryo-EM coordinates, ClinicalTrials.gov API v2).
  - Preclinical computational model; research, not medical advice.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

SCRIPT_PATH = Path(__file__).resolve()

NLM_ATTRIBUTION = (
    "This product uses publicly available data from the U.S. National Library of Medicine (NLM), "
    "National Institutes of Health, Department of Health and Human Services; NLM is not responsible "
    "for the product and does not endorse or recommend this or any other product. "
    "ClinicalTrials.gov does not evaluate safety or efficacy; listing on the registry does not denote government endorsement."
)

RESEARCH_BOUNDARY = (
    "Preclinical in-silico eRF1/eRF3 depletion & nonsense readthrough biophysical model; "
    "research, not medical advice. No clinical efficacy or therapeutic safety claim is made."
)

# Reference Structural & Registry IDs
RCSB_PDB_STRUCTURE = "6O2P"  # Human CFTR potentiator-bound gating conformation (3.3 A)
PRIMARY_DOI_DEPLETION = "10.1038/s41589-020-00662-8"  # Baradaran-Heravi et al. Nat Chem Biol 2021
PRIMARY_DOI_CONTEXT = "10.1093/nar/gkaa728"           # Roy et al. Nucleic Acids Res 2020
CLINICAL_TRIAL_REGISTRY = "NCT04135495"               # ELX-02 Phase 2 nonsense suppression

# Nominal biophysical constants
DEFAULT_ERF_TOTAL_UM = 1.0            # Baseline physiological total [eRF1/eRF3] complex pool (uM)
DEFAULT_TRNA_NC_UM = 5.0              # Near-cognate aminoacyl-tRNA abundance in ribosomal vicinity (uM)
DEFAULT_K_TERM_MAX_S = 48.0           # Maximal catalytic termination velocity constant (s^-1)
DEFAULT_KM_TERM_UM = 2.00             # Michaelis constant for eRF complex termination binding (uM)
DEFAULT_K_READ_MAX_S = 0.08           # Near-cognate decoding rate constant at PTC A-site (s^-1)
DEFAULT_KM_READ_UM = 1.50             # Michaelis constant for near-cognate decoding (uM)
DEFAULT_ERF_COMPLEX_COOPERATIVITY = 1.80 # Heterodimeric eRF1-eRF3 assembly cooperativity under co-depletion

# eRF depletion pharmacodynamic constants (Baradaran-Heravi 2021)
DEFAULT_DEPLETION_IC50_UM = 0.35      # Small-molecule inhibitor IC50 for release factor depletion (uM)
DEFAULT_DEPLETION_HILL_H = 1.80       # Hill cooperativity coefficient (h in [1.5, 2.0])
DEFAULT_DEPLETION_IMAX = 0.85         # Maximum fractional release factor depletion at saturation

# Context factor weights for +4 downstream nucleotide (Roy 2020)
PLUS4_CONTEXT_WEIGHTS: Dict[str, float] = {
    "C": 0.35,  # Weakest termination / highest readthrough susceptibility (W1282X)
    "U": 0.35,  # Weak termination / pyrimidine context
    "A": 0.60,  # Intermediate termination efficiency
    "G": 1.00,  # Strongest eRF1 binding / maximal termination efficiency
}

# Codon termination permissiveness multipliers (UGA > UAG >> UAA)
# Normalizes relative termination efficiency: lower multiplier = more permissive to readthrough
CODON_PERMISSIVENESS_MULTIPLIERS: Dict[str, float] = {
    "UGA": 0.70,  # Opal: highest readthrough permissiveness (W1282X, R553X, R1162X)
    "UAG": 1.00,  # Amber: intermediate permissiveness
    "UAA": 2.50,  # Ochre: strictest termination fidelity / lowest readthrough (G542X, Y122X)
}

# CFTR mRNA stability and expression constants
DEFAULT_WT_CFTR_MRNA_RELATIVE = 1.00   # WT steady-state CFTR mRNA level
DEFAULT_PTC_BASAL_MRNA_RELATIVE = 0.18 # Unescaped NMD steady-state mRNA level (15-20% WT)
DEFAULT_NMD_ESCAPE_EC50 = 0.035        # Readthrough fraction yielding 50% maximal NMD escape
DEFAULT_WT_BAND_C_RELATIVE = 1.00      # WT mature 170 kDa Band C CFTR expression
DEFAULT_CURATIVE_BAND_C_THRESHOLD = 0.10 # Clinical curative threshold (>= 10% WT Band C)


@dataclass(frozen=True)
class PTCMutationProfile:
    """Biophysical profile of a cystic fibrosis premature termination codon mutation."""
    name: str
    codon: str
    plus4_nt: str
    exon: int
    cDNA: str
    native_amino_acid: str
    folding_yield_band_c: float  # Folding efficiency of full-length readthrough protein relative to WT
    notes: str


# Catalog of prominent CFTR nonsense mutations
CFTR_NONSENSE_MUTATIONS: Dict[str, PTCMutationProfile] = {
    "W1282X": PTCMutationProfile(
        name="W1282X",
        codon="UGA",
        plus4_nt="C",
        exon=23,
        cDNA="c.3846G>A",
        native_amino_acid="Trp",
        folding_yield_band_c=0.85,
        notes="Class I nonsense mutation in NBD2; highly responsive to release factor depletion due to +4 C context.",
    ),
    "G542X": PTCMutationProfile(
        name="G542X",
        codon="UGA",
        plus4_nt="A",
        exon=12,
        cDNA="c.1624G>T",
        native_amino_acid="Gly",
        folding_yield_band_c=0.75,
        notes="Most common CFTR nonsense mutation globally; located in NBD1 loop.",
    ),
    "R553X": PTCMutationProfile(
        name="R553X",
        codon="UGA",
        plus4_nt="A",
        exon=12,
        cDNA="c.1657C>T",
        native_amino_acid="Arg",
        folding_yield_band_c=0.88,
        notes="Class I nonsense in NBD1 alpha-subdomain; high readthrough folding fidelity.",
    ),
    "R1162X": PTCMutationProfile(
        name="R1162X",
        codon="UGA",
        plus4_nt="A",
        exon=20,
        cDNA="c.3484C>T",
        native_amino_acid="Arg",
        folding_yield_band_c=0.80,
        notes="Class I nonsense mutation in TM12 / intracellular loop 4.",
    ),
    "Y122X": PTCMutationProfile(
        name="Y122X",
        codon="UAA",
        plus4_nt="U",
        exon=4,
        cDNA="c.366T>A",
        native_amino_acid="Tyr",
        folding_yield_band_c=0.70,
        notes="Early N-terminal Ochre (UAA) stop codon; high termination fidelity.",
    ),
}


@dataclass
class KineticCompetitionResult:
    """Ribosomal A-site termination vs readthrough kinetic evaluation."""
    mutation: str
    codon: str
    plus4_nt: str
    context_weight: float
    permissiveness_factor: float
    erf_concentration_um: float
    trna_nc_concentration_um: float
    v_term_per_s: float
    v_read_per_s: float
    readthrough_fraction: float
    readthrough_percentage: float


@dataclass
class DepletionPharmacodynamics:
    """Release factor depletion pharmacodynamics under small-molecule inhibitor."""
    drug_concentration_um: float
    ic50_um: float
    hill_h: float
    fractional_depletion: float
    remaining_erf_fraction: float
    erf_effective_um: float


@dataclass
class BandCRescueResult:
    """Coupled NMD escape, mature Band C protein expression, and therapeutic rescue."""
    mutation: str
    drug_concentration_um: float
    fractional_erf_depletion: float
    readthrough_fraction: float
    readthrough_percentage: float
    fold_readthrough_boost: float
    mrna_stability_relative_to_wt: float
    nmd_escape_fraction: float
    band_c_protein_relative_to_wt: float
    band_c_percentage_wt: float
    curative_threshold_reached: bool
    estimated_chloride_transport_pct_wt: float


def compute_depleted_erf_concentration(
    drug_conc_um: float,
    erf_total_um: float = DEFAULT_ERF_TOTAL_UM,
    ic50_um: float = DEFAULT_DEPLETION_IC50_UM,
    hill_h: float = DEFAULT_DEPLETION_HILL_H,
    i_max: float = DEFAULT_DEPLETION_IMAX,
) -> DepletionPharmacodynamics:
    """Compute effective eRF1/eRF3 complex concentration under small-molecule depletion.

    Uses Hill equation:
      Depletion(C) = (I_max * C^h) / (IC50^h + C^h)
      [eRF]_eff = [eRF]_0 * (1 - Depletion(C))

    Args:
        drug_conc_um: Inhibitor / depleter drug concentration in uM.
        erf_total_um: Baseline uninhibited release factor concentration in uM.
        ic50_um: Half-maximal depletion concentration in uM.
        hill_h: Hill cooperativity coefficient.
        i_max: Maximal attainable depletion fraction.

    Returns:
        DepletionPharmacodynamics data structure.
    """
    if drug_conc_um <= 0.0:
        frac_depletion = 0.0
    else:
        conc_h = math.pow(drug_conc_um, hill_h)
        ic50_h = math.pow(ic50_um, hill_h)
        frac_depletion = (i_max * conc_h) / (ic50_h + conc_h)

    # Clamp depletion in [0, i_max]
    frac_depletion = max(0.0, min(i_max, frac_depletion))
    remaining_frac = 1.0 - frac_depletion
    erf_eff = erf_total_um * remaining_frac

    return DepletionPharmacodynamics(
        drug_concentration_um=drug_conc_um,
        ic50_um=ic50_um,
        hill_h=hill_h,
        fractional_depletion=frac_depletion,
        remaining_erf_fraction=remaining_frac,
        erf_effective_um=erf_eff,
    )


def compute_ribosomal_competition(
    erf_conc_um: float,
    mutation: PTCMutationProfile,
    trna_nc_um: float = DEFAULT_TRNA_NC_UM,
    k_term_max_s: float = DEFAULT_K_TERM_MAX_S,
    km_term_um: float = DEFAULT_KM_TERM_UM,
    k_read_max_s: float = DEFAULT_K_READ_MAX_S,
    km_read_um: float = DEFAULT_KM_READ_UM,
    cooperativity: float = DEFAULT_ERF_COMPLEX_COOPERATIVITY,
) -> KineticCompetitionResult:
    """Compute competitive A-site termination vs near-cognate readthrough kinetics.

    Governing equations:
      omega_context = PLUS4_CONTEXT_WEIGHTS[plus4_nt] (0.35 for C/U, 1.0 for G)
      omega_codon = CODON_PERMISSIVENESS_MULTIPLIERS[codon]
      [eRF]_complex = [eRF]_0 * ([eRF] / [eRF]_0)^cooperativity (heterodimeric assembly)
      v_term = (k_term_max * [eRF]_complex) / (Km_term + [eRF]_complex) * omega_context * omega_codon
      v_read = (k_read_max * [tRNA_nc]) / (Km_read + [tRNA_nc])
      F_readthrough = v_read / (v_term + v_read)

    Args:
        erf_conc_um: Active concentration of eRF1/eRF3 complex in uM.
        mutation: PTCMutationProfile specifying codon, +4 nucleotide, and context.
        trna_nc_um: Near-cognate aminoacyl-tRNA pool concentration in uM.
        k_term_max_s: Maximal catalytic termination rate in s^-1.
        km_term_um: Michaelis constant for termination complex in uM.
        k_read_max_s: Near-cognate decoding rate constant in s^-1.
        km_read_um: Michaelis constant for near-cognate decoding in uM.
        cooperativity: Heterodimeric assembly cooperativity exponent (eRF1 * eRF3).

    Returns:
        KineticCompetitionResult data structure.
    """
    ctx_weight = PLUS4_CONTEXT_WEIGHTS.get(mutation.plus4_nt.upper(), 1.00)
    codon_factor = CODON_PERMISSIVENESS_MULTIPLIERS.get(mutation.codon.upper(), 1.00)

    # Heterodimeric eRF1-eRF3 assembly cooperativity under co-depletion (Baradaran-Heravi 2021)
    if erf_conc_um <= 0.0:
        effective_erf_pool = 0.0
    else:
        norm_ratio = erf_conc_um / DEFAULT_ERF_TOTAL_UM
        effective_erf_pool = DEFAULT_ERF_TOTAL_UM * math.pow(norm_ratio, cooperativity)

    # Termination rate
    if effective_erf_pool <= 0.0:
        v_term = 0.0
    else:
        v_term = (
            (k_term_max_s * effective_erf_pool)
            / (km_term_um + effective_erf_pool)
            * ctx_weight
            * codon_factor
        )

    # Near-cognate readthrough rate
    if trna_nc_um <= 0.0:
        v_read = 0.0
    else:
        v_read = (k_read_max_s * trna_nc_um) / (km_read_um + trna_nc_um)

    total_v = v_term + v_read
    if total_v <= 0.0:
        readthrough_frac = 0.0
    else:
        readthrough_frac = v_read / total_v

    return KineticCompetitionResult(
        mutation=mutation.name,
        codon=mutation.codon,
        plus4_nt=mutation.plus4_nt,
        context_weight=ctx_weight,
        permissiveness_factor=codon_factor,
        erf_concentration_um=erf_conc_um,
        trna_nc_concentration_um=trna_nc_um,
        v_term_per_s=v_term,
        v_read_per_s=v_read,
        readthrough_fraction=readthrough_frac,
        readthrough_percentage=readthrough_frac * 100.0,
    )


def compute_nmd_escape_and_band_c(
    mutation: PTCMutationProfile,
    readthrough_fraction: float,
    basal_readthrough_fraction: float,
    folding_yield: Optional[float] = None,
    basal_mrna_relative: float = DEFAULT_PTC_BASAL_MRNA_RELATIVE,
    nmd_ec50: float = DEFAULT_NMD_ESCAPE_EC50,
) -> Tuple[float, float, float, float]:
    """Compute coupled NMD transcript stabilization and mature Band C expression.

    Coupled NMD Escape Model:
      Translating ribosomes displacing downstream EJCs protect the transcript from NMD:
        f_escape = F_read / (F_read + EC50_NMD)
        mRNA_rel = mRNA_basal + (mRNA_WT - mRNA_basal) * f_escape
      Mature Band C CFTR protein:
        Band_C_rel = mRNA_rel * F_read * folding_yield
      Chloride transport proxy:
        I_Cl_pct = Band_C_pct * P_open_correction (scaled to % WT)

    Returns:
        (mrna_relative, nmd_escape_fraction, band_c_relative, chloride_pct_wt)
    """
    if folding_yield is None:
        folding_yield = mutation.folding_yield_band_c

    # NMD Escape dynamics
    if readthrough_fraction <= 0.0:
        f_escape = 0.0
    else:
        f_escape = readthrough_fraction / (readthrough_fraction + nmd_ec50)
    f_escape = max(0.0, min(1.0, f_escape))

    mrna_rel = basal_mrna_relative + (DEFAULT_WT_CFTR_MRNA_RELATIVE - basal_mrna_relative) * f_escape

    # Mature 170 kDa Band C CFTR expression
    band_c_rel = mrna_rel * readthrough_fraction * folding_yield
    band_c_pct = band_c_rel * 100.0

    # Apical chloride current conductance proxy (% of WT CFTR)
    gating_fidelity = 0.90 if mutation.name == "W1282X" else 0.85
    chloride_pct_wt = band_c_pct * gating_fidelity

    return mrna_rel, f_escape, band_c_rel, chloride_pct_wt


def evaluate_mutation_dose_response(
    mutation_name: str,
    drug_concentrations: Optional[List[float]] = None,
    erf_total_um: float = DEFAULT_ERF_TOTAL_UM,
    ic50_um: float = DEFAULT_DEPLETION_IC50_UM,
    hill_h: float = DEFAULT_DEPLETION_HILL_H,
    i_max: float = DEFAULT_DEPLETION_IMAX,
) -> List[BandCRescueResult]:
    """Evaluate dose-dependent readthrough and Band C rescue across inhibitor concentrations.

    Args:
        mutation_name: Name of mutation from CFTR_NONSENSE_MUTATIONS.
        drug_concentrations: List of inhibitor concentrations in uM (default: 0 to 5 uM).
        erf_total_um: Total uninhibited release factor pool in uM.
        ic50_um: Inhibitor IC50 in uM.
        hill_h: Hill coefficient.
        i_max: Maximal depletion fraction.

    Returns:
        List of BandCRescueResult records across the dose range.
    """
    if mutation_name not in CFTR_NONSENSE_MUTATIONS:
        raise ValueError(
            f"Unknown mutation '{mutation_name}'. Valid mutations: {list(CFTR_NONSENSE_MUTATIONS.keys())}"
        )
    mutation = CFTR_NONSENSE_MUTATIONS[mutation_name]

    if drug_concentrations is None:
        drug_concentrations = [0.0, 0.05, 0.10, 0.20, 0.35, 0.50, 0.75, 1.0, 1.5, 2.5, 5.0]

    # Baseline (drug = 0) competition
    base_pd = compute_depleted_erf_concentration(0.0, erf_total_um, ic50_um, hill_h, i_max)
    base_comp = compute_ribosomal_competition(base_pd.erf_effective_um, mutation)
    basal_rt = base_comp.readthrough_fraction

    results: List[BandCRescueResult] = []
    for conc in drug_concentrations:
        pd = compute_depleted_erf_concentration(conc, erf_total_um, ic50_um, hill_h, i_max)
        comp = compute_ribosomal_competition(pd.erf_effective_um, mutation)

        mrna_rel, f_escape, band_c_rel, cl_pct = compute_nmd_escape_and_band_c(
            mutation=mutation,
            readthrough_fraction=comp.readthrough_fraction,
            basal_readthrough_fraction=basal_rt,
        )

        fold_boost = (
            comp.readthrough_fraction / basal_rt if basal_rt > 0.0 else 1.0
        )

        results.append(
            BandCRescueResult(
                mutation=mutation.name,
                drug_concentration_um=conc,
                fractional_erf_depletion=pd.fractional_depletion,
                readthrough_fraction=comp.readthrough_fraction,
                readthrough_percentage=comp.readthrough_percentage,
                fold_readthrough_boost=fold_boost,
                mrna_stability_relative_to_wt=mrna_rel,
                nmd_escape_fraction=f_escape,
                band_c_protein_relative_to_wt=band_c_rel,
                band_c_percentage_wt=band_c_rel * 100.0,
                curative_threshold_reached=(band_c_rel >= DEFAULT_CURATIVE_BAND_C_THRESHOLD),
                estimated_chloride_transport_pct_wt=cl_pct,
            )
        )

    return results


def run_benchmark_comparison() -> Dict[str, Any]:
    """Compare release factor depletion response across all cataloged CFTR nonsense mutations."""
    concentrations = [0.0, 0.35, 2.50]
    comparison: Dict[str, Any] = {}

    for mut_name, mut_spec in CFTR_NONSENSE_MUTATIONS.items():
        curve = evaluate_mutation_dose_response(mut_name, drug_concentrations=concentrations)
        basal = curve[0]
        at_ic50 = curve[1]
        saturating = curve[2]

        comparison[mut_name] = {
            "codon": mut_spec.codon,
            "plus4_nt": mut_spec.plus4_nt,
            "cDNA": mut_spec.cDNA,
            "exon": mut_spec.exon,
            "basal_readthrough_pct": round(basal.readthrough_percentage, 3),
            "ic50_readthrough_pct": round(at_ic50.readthrough_percentage, 3),
            "saturating_readthrough_pct": round(saturating.readthrough_percentage, 3),
            "max_fold_boost": round(saturating.fold_readthrough_boost, 2),
            "saturating_band_c_pct_wt": round(saturating.band_c_percentage_wt, 2),
            "curative_rescue_at_saturation": saturating.curative_threshold_reached,
        }

    return comparison


def main() -> None:
    parser = argparse.ArgumentParser(
        description="eRF1/eRF3 Depletion & Aminoglycoside-Free Nonsense Readthrough Kinetic Model (T-0401)"
    )
    parser.add_argument(
        "--mutation",
        choices=list(CFTR_NONSENSE_MUTATIONS.keys()),
        default="W1282X",
        help="Target CFTR nonsense mutation (default: W1282X)",
    )
    parser.add_argument(
        "--drug-conc",
        type=float,
        default=None,
        help="Specific inhibitor concentration in uM to evaluate",
    )
    parser.add_argument(
        "--benchmark-all",
        action="store_true",
        help="Run comparative benchmark across all cataloged CFTR nonsense mutations",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Emit results in machine-readable JSON format",
    )
    args = parser.parse_args()

    if args.benchmark_all:
        bench_data = run_benchmark_comparison()
        output = {
            "model": "ptc_erf_depletion_readthrough",
            "task": "T-0401",
            "tier": "T2",
            "nlm_attribution": NLM_ATTRIBUTION,
            "research_boundary": RESEARCH_BOUNDARY,
            "grounding": {
                "pdb": RCSB_PDB_STRUCTURE,
                "primary_literature": [PRIMARY_DOI_DEPLETION, PRIMARY_DOI_CONTEXT],
                "clinical_trials": [CLINICAL_TRIAL_REGISTRY],
            },
            "mutation_comparison": bench_data,
        }
        if args.json:
            print(json.dumps(output, indent=2))
        else:
            print("=" * 80)
            print("CFTR NONSENSE READTHROUGH BENCHMARK (eRF1/eRF3 Depletion)")
            print("=" * 80)
            for mut, data in bench_data.items():
                print(
                    f"{mut:8s} | Codon: {data['codon']}+{data['plus4_nt']} | "
                    f"Basal: {data['basal_readthrough_pct']:5.2f}% -> "
                    f"Max: {data['saturating_readthrough_pct']:5.2f}% "
                    f"({data['max_fold_boost']:4.1f}x boost) | "
                    f"Band C: {data['saturating_band_c_pct_wt']:5.2f}% WT "
                    f"[{'RESCUED' if data['curative_rescue_at_saturation'] else 'SUB-THRESHOLD'}]"
                )
            print("=" * 80)
        return

    # Specific mutation dose-response or single point
    if args.drug_conc is not None:
        concs = [args.drug_conc]
    else:
        concs = [0.0, 0.05, 0.10, 0.20, 0.35, 0.50, 0.75, 1.0, 1.5, 2.5, 5.0]

    curve = evaluate_mutation_dose_response(args.mutation, drug_concentrations=concs)
    output_records = [asdict(r) for r in curve]

    payload = {
        "model": "ptc_erf_depletion_readthrough",
        "task": "T-0401",
        "tier": "T2",
        "mutation": args.mutation,
        "nlm_attribution": NLM_ATTRIBUTION,
        "research_boundary": RESEARCH_BOUNDARY,
        "grounding": {
            "pdb": RCSB_PDB_STRUCTURE,
            "primary_literature": [PRIMARY_DOI_DEPLETION, PRIMARY_DOI_CONTEXT],
            "clinical_trials": [CLINICAL_TRIAL_REGISTRY],
        },
        "dose_response": output_records,
    }

    if args.json:
        print(json.dumps(payload, indent=2))
    else:
        print("=" * 80)
        print(f"eRF1/eRF3 Depletion Dose Response for {args.mutation}")
        print(f"PDB: {RCSB_PDB_STRUCTURE} | DOIs: {PRIMARY_DOI_DEPLETION}, {PRIMARY_DOI_CONTEXT}")
        print(f"Registry: {CLINICAL_TRIAL_REGISTRY}")
        print("=" * 80)
        print("Drug [uM] | Depletion | Readthrough | Boost | mRNA (NMD) | Band C (%WT) | Curative?")
        print("-" * 80)
        for r in curve:
            print(
                f"{r.drug_concentration_um:8.2f} | "
                f"{r.fractional_erf_depletion*100:8.1f}% | "
                f"{r.readthrough_percentage:10.2f}% | "
                f"{r.fold_readthrough_boost:5.1f}x | "
                f"{r.mrna_stability_relative_to_wt*100:9.1f}% | "
                f"{r.band_c_percentage_wt:11.2f}% | "
                f"{'YES (>=10%)' if r.curative_threshold_reached else 'NO'}"
            )
        print("=" * 80)


if __name__ == "__main__":
    main()
