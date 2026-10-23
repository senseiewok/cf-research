#!/usr/bin/env python3
"""Dual-AAV Split-Intein Trans-Splicing Biophysical Reconstitution Model for CFTR (T-0206).

Models vector genome packaging, split-intein junction engineering, Poisson
co-transduction probability, second-order protein trans-splicing kinetics, and
steady-state reconstituted mature Band C CFTR yield across airway epithelia.

Grounded Sources & Provenance:
- landscape/discovery-dual-aav-intein-packaging.md
- DOI 10.1126/scitranslmed.aav4523: Tornabene P et al. Intein-mediated dual AAV trans-splicing for large genes.
- DOI 10.3389/fgeed.2023.1271813: Dual AAV intein vectors for CFTR gene therapy.
- RCSB PDB 6MSM: Cryo-EM structure of human CFTR phosphorylated, ATP-bound active conductance states at 3.2 A.
- RCSB PDB 5UAK: Cryo-EM structure of human CFTR dephosphorylated conformation at 3.87 A.
- ClinicalTrials.gov NCT05248230:
  Phase 1/2 aerosolized CFTR gene therapy (4D-710).

NLM Attribution & Research Boundary:
This tool incorporates clinical trial metadata from ClinicalTrials.gov, a service of the
U.S. National Institutes of Health / National Library of Medicine.
Research, not medical advice. Vector biophysics and protein trans-splicing kinetics model,
not individual treatment directives.
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
    "Research, not medical advice. Vector biophysics and protein trans-splicing kinetics model, "
    "not individual treatment directives."
)

# Physical & Genomic Constants
AAV_PACKAGING_HARD_LIMIT_NT: int = 4700  # Strict encapsidation ceiling between ITRs
AAV_PACKAGING_MAX_STRETCH_NT: int = 5000  # Truncation and packaging failure above this
ITR_LENGTH_NT: int = 145  # Standard wild-type AAV2 inverted terminal repeat
CFTR_WT_RESIDUES: int = 1480  # Full-length human CFTR amino acids
CFTR_WT_CDS_BP: int = 4440  # 1,480 aa * 3 bp/codon (excluding stop)

# Cryo-EM Domain Boundaries (PDB 6MSM / 5UAK)
DOMAIN_BOUNDARIES = {
    "TMD1": (81, 365),
    "NBD1": (423, 646),
    "R_DOMAIN": (654, 831),  # Unplaced disordered region in cryo-EM models
    "TMD2": (859, 1155),
    "NBD2": (1210, 1443),
}

# Split Intein Biophysical Parameters (DOI 10.1126/scitranslmed.aav4523)
INTEIN_SYSTEMS: Dict[str, Dict[str, Any]] = {
    "Cfa": {
        "name": "Cfa (Consensus Fast DnaE Intein)",
        "int_n_aa": 101,
        "int_n_bp": 303,
        "int_c_aa": 36,
        "int_c_bp": 108,
        "catalytic_residues_plus1": ("C", "S"),  # Accepts Cys or Ser at extein +1
        "k_splice_s_inv": 0.015,  # ~1.5e-2 s^-1
        "t_half_seconds": 46.2,  # ln(2)/k
        "in_vitro_efficiency": 0.42,  # ~42% trans-splicing in primary airway ALI
        "airway_cl_transport_pct_wt": 42.0,  # 42% WT chloride current in CF epithelia
    },
    "Npu_DnaE": {
        "name": "Npu DnaE (Nostoc punctiforme Intein)",
        "int_n_aa": 102,
        "int_n_bp": 306,
        "int_c_aa": 36,
        "int_c_bp": 108,
        "catalytic_residues_plus1": ("C",),  # Strict requirement for Cys at +1
        "k_splice_s_inv": 0.005,  # ~5.0e-3 s^-1
        "t_half_seconds": 138.6,  # ln(2)/k
        "in_vitro_efficiency": 0.32,  # ~32% trans-splicing
        "airway_cl_transport_pct_wt": 34.0,  # 34% WT chloride current
    },
}


@dataclass(frozen=True)
class PackagingBudget:
    vector_name: str
    total_genome_nt: int
    hard_limit_nt: int
    fits_packaging_limit: bool
    payload_excess_nt: int
    elements: Dict[str, int]


@dataclass(frozen=True)
class SplitInteinDesign:
    intein_system: str
    split_residue: int
    junction_c_plus1_aa: str
    is_valid_junction: bool
    split_in_disordered_r_domain: bool
    n_vector: PackagingBudget
    c_vector: PackagingBudget
    int_n_size_bp: int
    int_c_size_bp: int
    k_splice_s_inv: float
    t_half_seconds: float
    predicted_splicing_efficiency: float


@dataclass(frozen=True)
class CoTransductionResult:
    moi_n: float
    moi_c: float
    p_transduction_n: float
    p_transduction_c: float
    p_dual_co_transduction: float
    percent_dual_positive_cells: float


@dataclass(frozen=True)
class DualAAVSimulationReport:
    title: str
    nlm_attribution: str
    research_boundary: str
    intein_design: SplitInteinDesign
    co_transduction: CoTransductionResult
    reconstituted_mature_cftr_yield_pct_wt: float
    restored_chloride_transport_pct_wt: float
    relevant_nct_ids: List[str]
    evidence_citations: List[str]


def check_single_aav_packaging_capacity(
    cds_bp: int = CFTR_WT_CDS_BP,
    promoter_bp: int = 600,  # e.g. CMV 600 bp
    poly_a_bp: int = 215,  # e.g. bGH 215 bp
    utr_bp: int = 80,  # 5' UTR / Kozak
    itr_count: int = 2,
) -> PackagingBudget:
    """Evaluate whether a single AAV vector configuration can package full-length CFTR."""
    itr_total = itr_count * ITR_LENGTH_NT
    elements = {
        "5_ITR": ITR_LENGTH_NT,
        "promoter": promoter_bp,
        "5_utr_kozak": utr_bp,
        "cds": cds_bp,
        "poly_a": poly_a_bp,
        "3_ITR": ITR_LENGTH_NT,
    }
    total = sum(elements.values())
    fits = total <= AAV_PACKAGING_HARD_LIMIT_NT
    excess = max(0, total - AAV_PACKAGING_HARD_LIMIT_NT)

    return PackagingBudget(
        vector_name="Single-AAV-Full-Length-CFTR",
        total_genome_nt=total,
        hard_limit_nt=AAV_PACKAGING_HARD_LIMIT_NT,
        fits_packaging_limit=fits,
        payload_excess_nt=excess,
        elements=elements,
    )


def design_split_cftr_vectors(
    split_residue: int = 744,
    intein_key: str = "Cfa",
    junction_c_plus1_aa: str = "S",
    promoter_n_bp: int = 600,
    promoter_c_bp: int = 250,  # Minimal promoter / splice acceptor in vector C
    poly_a_n_bp: int = 150,
    poly_a_c_bp: int = 215,
    utr_n_bp: int = 80,
    utr_c_bp: int = 40,
) -> SplitInteinDesign:
    """Design dual-AAV split-intein vectors partitioned at a given CFTR amino acid residue."""
    if intein_key not in INTEIN_SYSTEMS:
        raise ValueError(f"Unknown intein system: {intein_key}. Choices: {list(INTEIN_SYSTEMS.keys())}")

    sys_info = INTEIN_SYSTEMS[intein_key]
    allowed_plus1 = sys_info["catalytic_residues_plus1"]
    is_valid_junction = junction_c_plus1_aa in allowed_plus1

    # Check if split resides in disordered R-domain (654-831)
    r_start, r_end = DOMAIN_BOUNDARIES["R_DOMAIN"]
    in_r_domain = r_start <= split_residue <= r_end

    # Partition coding sequence
    n_cds_bp = split_residue * 3
    c_cds_bp = (CFTR_WT_RESIDUES - split_residue) * 3

    # N-terminal vector: ITR + promoter + 5' UTR + CFTR_N + IntN + polyA + ITR
    n_elements = {
        "5_ITR": ITR_LENGTH_NT,
        "promoter": promoter_n_bp,
        "5_utr": utr_n_bp,
        "cftr_n_cds": n_cds_bp,
        "int_n_cds": sys_info["int_n_bp"],
        "poly_a": poly_a_n_bp,
        "3_ITR": ITR_LENGTH_NT,
    }
    n_total = sum(n_elements.values())
    n_budget = PackagingBudget(
        vector_name="AAV-5'-CFTR-IntN",
        total_genome_nt=n_total,
        hard_limit_nt=AAV_PACKAGING_HARD_LIMIT_NT,
        fits_packaging_limit=n_total <= AAV_PACKAGING_HARD_LIMIT_NT,
        payload_excess_nt=max(0, n_total - AAV_PACKAGING_HARD_LIMIT_NT),
        elements=n_elements,
    )

    # C-terminal vector: ITR + promoter/acceptor + 5' UTR + IntC + CFTR_C + polyA + ITR
    c_elements = {
        "5_ITR": ITR_LENGTH_NT,
        "promoter_or_acceptor": promoter_c_bp,
        "5_utr": utr_c_bp,
        "int_c_cds": sys_info["int_c_bp"],
        "cftr_c_cds": c_cds_bp,
        "poly_a": poly_a_c_bp,
        "3_ITR": ITR_LENGTH_NT,
    }
    c_total = sum(c_elements.values())
    c_budget = PackagingBudget(
        vector_name="AAV-IntC-3'-CFTR",
        total_genome_nt=c_total,
        hard_limit_nt=AAV_PACKAGING_HARD_LIMIT_NT,
        fits_packaging_limit=c_total <= AAV_PACKAGING_HARD_LIMIT_NT,
        payload_excess_nt=max(0, c_total - AAV_PACKAGING_HARD_LIMIT_NT),
        elements=c_elements,
    )

    # Splicing efficiency calculation
    base_eff = sys_info["in_vitro_efficiency"]
    if not is_valid_junction:
        eff = 0.0
    elif not in_r_domain:
        # Penalize if split in a folded domain (misfolding / steric impediment)
        eff = base_eff * 0.15
    else:
        eff = base_eff

    return SplitInteinDesign(
        intein_system=intein_key,
        split_residue=split_residue,
        junction_c_plus1_aa=junction_c_plus1_aa,
        is_valid_junction=is_valid_junction,
        split_in_disordered_r_domain=in_r_domain,
        n_vector=n_budget,
        c_vector=c_budget,
        int_n_size_bp=sys_info["int_n_bp"],
        int_c_size_bp=sys_info["int_c_bp"],
        k_splice_s_inv=sys_info["k_splice_s_inv"],
        t_half_seconds=sys_info["t_half_seconds"],
        predicted_splicing_efficiency=round(eff, 4),
    )


def calculate_poisson_co_transduction(
    moi_n: float,
    moi_c: float,
) -> CoTransductionResult:
    """Calculate independent Poisson probabilities of cellular dual-vector entry."""
    if moi_n < 0 or moi_c < 0:
        raise ValueError("MOI must be non-negative.")

    # P(at least 1 vector entry) = 1 - e^(-MOI)
    p_n = 1.0 - math.exp(-moi_n)
    p_c = 1.0 - math.exp(-moi_c)
    p_dual = p_n * p_c

    return CoTransductionResult(
        moi_n=round(moi_n, 3),
        moi_c=round(moi_c, 3),
        p_transduction_n=round(p_n, 4),
        p_transduction_c=round(p_c, 4),
        p_dual_co_transduction=round(p_dual, 4),
        percent_dual_positive_cells=round(p_dual * 100.0, 2),
    )


def simulate_second_order_trans_splicing_kinetics(
    time_seconds: float,
    initial_n_um: float = 1.0,
    initial_c_um: float = 1.0,
    intein_key: str = "Cfa",
) -> Dict[str, float]:
    """Model second-order association and trans-splicing yield over time."""
    if intein_key not in INTEIN_SYSTEMS:
        raise ValueError(f"Unknown intein system: {intein_key}")

    k = INTEIN_SYSTEMS[intein_key]["k_splice_s_inv"]

    # When equal initial concentrations: [Product] = [A]0 * (k * [A]0 * t) / (1 + k * [A]0 * t)
    # Using normalized dimensionless fraction:
    frac = (k * time_seconds) / (1.0 + k * time_seconds) if time_seconds > 0 else 0.0
    product_um = min(initial_n_um, initial_c_um) * frac

    return {
        "time_seconds": round(time_seconds, 1),
        "fraction_spliced": round(frac, 4),
        "product_concentration_um": round(product_um, 4),
        "rate_constant_k": k,
    }


def simulate_dual_aav_cftr_reconstitution(
    moi: float = 3.0,
    split_residue: int = 744,
    intein_key: str = "Cfa",
    junction_c_plus1_aa: str = "S",
) -> DualAAVSimulationReport:
    """Simulate complete dual-AAV reconstitution pipeline."""
    design = design_split_cftr_vectors(
        split_residue=split_residue,
        intein_key=intein_key,
        junction_c_plus1_aa=junction_c_plus1_aa,
    )
    co_tx = calculate_poisson_co_transduction(moi_n=moi, moi_c=moi)

    # Mature Band C yield = P(dual) * splicing_efficiency * 100
    mature_yield = co_tx.p_dual_co_transduction * design.predicted_splicing_efficiency * 100.0

    # Chloride transport rescue benchmark (DOI 10.3389/fgeed.2023.1271813 reports ~42% WT at MOI 3)
    ref_cl = INTEIN_SYSTEMS[intein_key]["airway_cl_transport_pct_wt"]
    scaling = co_tx.p_dual_co_transduction / (1.0 - math.exp(-3.0)) ** 2 if moi > 0 else 0.0
    cl_transport = min(100.0, ref_cl * scaling * (design.predicted_splicing_efficiency / INTEIN_SYSTEMS[intein_key]["in_vitro_efficiency"]))

    return DualAAVSimulationReport(
        title="Dual-AAV Split-Intein Trans-Splicing CFTR Biophysical Reconstitution Model",
        nlm_attribution=NLM_ATTRIBUTION,
        research_boundary=RESEARCH_BOUNDARY,
        intein_design=design,
        co_transduction=co_tx,
        reconstituted_mature_cftr_yield_pct_wt=round(mature_yield, 2),
        restored_chloride_transport_pct_wt=round(cl_transport, 2),
        relevant_nct_ids=["NCT05248230"],
        evidence_citations=[
            "DOI: 10.1126/scitranslmed.aav4523",
            "DOI: 10.3389/fgeed.2023.1271813",
            "ClinicalTrials.gov: NCT05248230",
            "PDB: 6MSM",
            "PDB: 5UAK",
        ],
    )


def self_test() -> bool:
    """Deterministic internal self-test suite."""
    # 1. Single-AAV full-length CFTR must fail packaging capacity
    single_budget = check_single_aav_packaging_capacity()
    assert not single_budget.fits_packaging_limit, "Single AAV full-length CFTR must exceed 4,700 nt"
    assert single_budget.total_genome_nt > 5200, "Single AAV genome should exceed 5,200 nt"

    # 2. Split vectors at residue 744 with Cfa must both fit within 4,700 nt
    split_design = design_split_cftr_vectors(split_residue=744, intein_key="Cfa", junction_c_plus1_aa="S")
    assert split_design.n_vector.fits_packaging_limit, "Vector N must fit <= 4,700 nt"
    assert split_design.c_vector.fits_packaging_limit, "Vector C must fit <= 4,700 nt"
    assert split_design.split_in_disordered_r_domain, "Residue 744 must reside in disordered R-domain"
    assert split_design.is_valid_junction, "Serine at +1 must be valid for Cfa intein"

    # 3. Poisson co-transduction scaling
    co_1 = calculate_poisson_co_transduction(moi_n=1.0, moi_c=1.0)
    co_3 = calculate_poisson_co_transduction(moi_n=3.0, moi_c=3.0)
    co_10 = calculate_poisson_co_transduction(moi_n=10.0, moi_c=10.0)
    assert co_1.p_dual_co_transduction < co_3.p_dual_co_transduction < co_10.p_dual_co_transduction
    assert math.isclose(co_3.p_dual_co_transduction, (1.0 - math.exp(-3.0)) ** 2, rel_tol=1e-3)

    # 4. Compare Cfa vs Npu DnaE
    sim_cfa = simulate_dual_aav_cftr_reconstitution(moi=3.0, intein_key="Cfa", junction_c_plus1_aa="S")
    sim_npu = simulate_dual_aav_cftr_reconstitution(moi=3.0, intein_key="Npu_DnaE", junction_c_plus1_aa="C")
    assert sim_cfa.intein_design.k_splice_s_inv > sim_npu.intein_design.k_splice_s_inv
    assert sim_cfa.reconstituted_mature_cftr_yield_pct_wt > sim_npu.reconstituted_mature_cftr_yield_pct_wt

    # 5. Invalid junction residue check
    bad_junction = design_split_cftr_vectors(split_residue=744, intein_key="Npu_DnaE", junction_c_plus1_aa="P")
    assert not bad_junction.is_valid_junction, "Proline at +1 must be invalid for Npu DnaE"
    assert bad_junction.predicted_splicing_efficiency == 0.0

    return True


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Dual-AAV Split-Intein Trans-Splicing CFTR Biophysical Reconstitution Model (T-0206)"
    )
    parser.add_argument(
        "--moi",
        type=float,
        default=3.0,
        help="Multiplicity of Infection (MOI) per half-vector (default: 3.0)",
    )
    parser.add_argument(
        "--intein",
        type=str,
        default="Cfa",
        choices=list(INTEIN_SYSTEMS.keys()),
        help="Split intein architecture (default: Cfa)",
    )
    parser.add_argument(
        "--split-site",
        type=int,
        default=744,
        help="CFTR split residue (default: 744 in disordered R-domain)",
    )
    parser.add_argument(
        "--plus1-aa",
        type=str,
        default="S",
        help="Catalytic amino acid residue at C-extein +1 position (default: S)",
    )
    parser.add_argument("--json", action="store_true", help="Print structured JSON output")
    parser.add_argument("--export", type=str, help="Export JSON results to output path")
    parser.add_argument("--self-test", action="store_true", help="Run deterministic internal self-tests")
    args = parser.parse_args()

    if args.self_test:
        success = self_test()
        if success:
            print("PASS: dual_aav_intein_model self-test passed.")
            sys.exit(0)
        else:
            print("FAIL: dual_aav_intein_model self-test failed.")
            sys.exit(1)

    report = simulate_dual_aav_cftr_reconstitution(
        moi=args.moi,
        split_residue=args.split_site,
        intein_key=args.intein,
        junction_c_plus1_aa=args.plus1_aa,
    )
    report_dict = asdict(report)

    if args.export:
        out_p = Path(args.export)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        with open(out_p, "w", encoding="utf-8") as f:
            json.dump(report_dict, f, indent=2)
        print(f"Exported dual-AAV intein reconstitution model to {out_p}")

    if args.json:
        print(json.dumps(report_dict, indent=2))
    else:
        print("=" * 80)
        print(report.title.upper())
        print(RESEARCH_BOUNDARY)
        print("=" * 80)
        des = report.intein_design
        print(f"Intein Architecture: {des.intein_system} (Split at residue {des.split_residue}, +1 aa: '{des.junction_c_plus1_aa}')")
        print(f"  In Disordered R-Domain: {des.split_in_disordered_r_domain} | Valid Catalytic Junction: {des.is_valid_junction}")
        print(f"  Vector N Genome: {des.n_vector.total_genome_nt} nt (limit: {des.n_vector.hard_limit_nt} nt, fits: {des.n_vector.fits_packaging_limit})")
        print(f"  Vector C Genome: {des.c_vector.total_genome_nt} nt (limit: {des.c_vector.hard_limit_nt} nt, fits: {des.c_vector.fits_packaging_limit})")
        print(f"  Trans-Splicing Kinetics: k = {des.k_splice_s_inv} s^-1 (t1/2 = {des.t_half_seconds} s)")
        print(f"  Predicted Trans-Splicing Efficiency: {des.predicted_splicing_efficiency * 100:.1f}%")
        ctx = report.co_transduction
        print(f"\nPoisson Co-Transduction (MOI = {ctx.moi_n}):")
        print(f"  Single Vector Transduction: {ctx.p_transduction_n * 100:.1f}%")
        print(f"  Dual Co-Transduction: {ctx.percent_dual_positive_cells}% of airway cells receive both vectors")
        print(f"\nReconstitution Outcomes:")
        print(f"  Mature 170 kDa Band C CFTR Yield: {report.reconstituted_mature_cftr_yield_pct_wt}% of WT")
        print(f"  Restored Transepithelial Cl- Transport: {report.restored_chloride_transport_pct_wt}% of WT")
        print(f"  Clinical Reference Trial: {', '.join(report.relevant_nct_ids)}")
        print("=" * 80)


if __name__ == "__main__":
    main()
