#!/usr/bin/env python3
"""CFTR Delta-R Domain Minigene Packaging & Single-Channel Electrophysiology Simulator (T-0207).

Models single-vector AAV packaging limits, R-domain structural coordinates (PDB 6MSM/5UAK),
single-channel patch-clamp conductance (gamma = 8.5 pS), ATP-dependent gating open
probability (Po), macroscopic apical current, and ER Band C maturation.

Grounded Sources & Provenance:
- landscape/discovery-dual-aav-intein-packaging.md
- DOI 10.1073/pnas.1901174116: Ostedgaard LS, Welsh MJ et al. CFTR Delta-R domain minigene
  single-channel conductance and packaging in AAV. PNAS 2019.
- RCSB PDB 6MSM: Cryo-EM structure of human CFTR phosphorylated, ATP-bound active conductance states at 3.2 A.
- RCSB PDB 5UAK: Cryo-EM structure of human CFTR dephosphorylated conformation at 3.87 A.
- ClinicalTrials.gov NCT05248230 (referenced in programs under NCT05248022):
  Phase 1/2 aerosolized CFTR gene therapy (4D-710 encoding CFTR Delta-R).

NLM Attribution & Research Boundary:
This tool incorporates clinical trial metadata from ClinicalTrials.gov, a service of the
U.S. National Institutes of Health / National Library of Medicine.
Research, not medical advice. Minigene biophysics and electrophysiology simulation,
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
    "Research, not medical advice. Minigene biophysics and electrophysiology simulation, "
    "not individual treatment directives."
)

# Physical Constants & Packaging Limits
AAV_HARD_PACKAGING_LIMIT_NT: int = 4700  # Strict AAV2 encapsidation limit between ITRs
ITR_LENGTH_NT: int = 145  # 145 bp per wild-type AAV ITR
UNITARY_CONDUCTANCE_PS: float = 8.5  # 8.5 pS unitary single-channel conductance in symmetrical Cl-

# Structural Boundaries (PDB 6MSM / 5UAK)
R_DOMAIN_START_RESIDUE: int = 654
R_DOMAIN_END_RESIDUE: int = 831

# Regulatory Genetic Element Catalogs
PROMOTERS: Dict[str, Dict[str, Any]] = {
    "tg83": {
        "name": "tg83 (Synthetic Ultra-Compact Promoter)",
        "length_nt": 83,
        "description": "83-bp synthetic minimal core promoter designed for tight AAV packaging limits.",
    },
    "minCMV": {
        "name": "minCMV (Minimal Cytomegalovirus Promoter)",
        "length_nt": 200,
        "description": "200-bp minimal core CMV promoter without extended enhancer.",
    },
    "CBA_full": {
        "name": "CBA (Full Chicken Beta-Actin Promoter)",
        "length_nt": 800,
        "description": "800-bp robust ubiquitous promoter; too large for single AAV with CFTR cDNA.",
    },
}

POLYA_SIGNALS: Dict[str, Dict[str, Any]] = {
    "synPolyA_compact": {
        "name": "Synthetic Compact PolyA",
        "length_nt": 35,
        "description": "35-bp minimal synthetic polyadenylation signal (Levitt et al.).",
    },
    "synPolyA_standard": {
        "name": "Synthetic Standard PolyA",
        "length_nt": 95,
        "description": "95-bp standard synthetic polyA.",
    },
    "bGH_polyA": {
        "name": "bGH PolyA (Bovine Growth Hormone)",
        "length_nt": 215,
        "description": "215-bp high-efficiency polyA signal.",
    },
}

CONSTRUCTS: Dict[str, Dict[str, Any]] = {
    "CFTR_WT": {
        "name": "Full-Length Wild-Type CFTR",
        "deleted_residues": [],
        "deleted_residues_count": 0,
        "deleted_bp": 0,
        "total_amino_acids": 1480,
        "cds_length_bp": 4440,
        "apparent_mw_kda": 170.0,
        "fraction_band_c_mature": 0.28,  # ~28% wild-type ER exit efficiency
        "basal_po": 0.008,  # Tightly autoinhibited when unphosphorylated
        "stimulated_po": 0.40,  # Fully stimulated by PKA + ATP
        "structural_integrity": "Native cryo-EM structure with intact disordered R-domain (PDB 6MSM).",
    },
    "CFTR_DeltaR_708_759": {
        "name": "CFTR Delta-R(708-759) (Ostedgaard Minigene)",
        "deleted_residues": [708, 759],
        "deleted_residues_count": 52,
        "deleted_bp": 156,
        "total_amino_acids": 1428,
        "cds_length_bp": 4284,
        "apparent_mw_kda": 165.0,
        "fraction_band_c_mature": 0.30,  # Comparable or slightly superior ER export
        "basal_po": 0.075,  # Elevated basal activity due to deletion of autoinhibitory segment
        "stimulated_po": 0.42,  # Retains full response to PKA/cAMP + ATP
        "structural_integrity": "Excised 52 aa reside within cryo-EM unplaced disordered region; NBD1-TMD2 interfaces intact.",
    },
    "CFTR_DeltaR_708_835": {
        "name": "CFTR Delta-R(708-835) (Extended R-Deletion)",
        "deleted_residues": [708, 835],
        "deleted_residues_count": 128,
        "deleted_bp": 384,
        "total_amino_acids": 1352,
        "cds_length_bp": 4056,
        "apparent_mw_kda": 155.0,
        "fraction_band_c_mature": 0.25,
        "basal_po": 0.12,  # Higher basal leak
        "stimulated_po": 0.35,  # Slightly perturbed phosphorylation response
        "structural_integrity": "Removes majority of R-domain up to TMD2 boundary.",
    },
    "CFTR_Delta264": {
        "name": "CFTR Delta-264 (N-Terminal Truncation Control)",
        "deleted_residues": [1, 264],
        "deleted_residues_count": 264,
        "deleted_bp": 792,
        "total_amino_acids": 1216,
        "cds_length_bp": 3648,
        "apparent_mw_kda": 135.0,
        "fraction_band_c_mature": 0.02,  # Severe misfolding and ERAD degradation
        "basal_po": 0.0,
        "stimulated_po": 0.0,  # Fails to reach apical surface
        "structural_integrity": "Loss of N-terminal cytoplasmic hairpin and TM1/TM2 disrupts membrane insertion.",
    },
}


@dataclass(frozen=True)
class MinigenePackagingBudget:
    construct_key: str
    construct_name: str
    promoter_key: str
    promoter_bp: int
    polya_key: str
    polya_bp: int
    utr_bp: int
    itrs_bp: int
    cds_bp: int
    total_genome_nt: int
    hard_limit_nt: int
    fits_single_aav: bool
    packaging_headroom_nt: int  # Positive if fits, negative if over capacity


@dataclass(frozen=True)
class PatchClampSimulationResult:
    construct_key: str
    pka_stimulated: bool
    atp_concentration_mm: float
    open_probability_po: float
    unitary_conductance_ps: float
    single_channel_current_pa: float
    apical_channels_n: int
    fraction_band_c_mature: float
    macroscopic_current_pa: float
    v_membrane_mv: float
    e_cl_mv: float
    driving_force_mv: float


@dataclass(frozen=True)
class MinigeneSimulationReport:
    title: str
    nlm_attribution: str
    research_boundary: str
    packaging: MinigenePackagingBudget
    electrophysiology: PatchClampSimulationResult
    structural_notes: str
    relevant_nct_ids: List[str]
    evidence_citations: List[str]


def evaluate_minigene_packaging(
    construct_key: str = "CFTR_DeltaR_708_759",
    promoter_key: str = "tg83",
    polya_key: str = "synPolyA_compact",
    utr_bp: int = 6,
) -> MinigenePackagingBudget:
    """Calculate single-AAV genome packaging budget for a specified minigene design."""
    if construct_key not in CONSTRUCTS:
        raise ValueError(f"Unknown construct: {construct_key}. Choices: {list(CONSTRUCTS.keys())}")
    if promoter_key not in PROMOTERS:
        raise ValueError(f"Unknown promoter: {promoter_key}. Choices: {list(PROMOTERS.keys())}")
    if polya_key not in POLYA_SIGNALS:
        raise ValueError(f"Unknown polyA: {polya_key}. Choices: {list(POLYA_SIGNALS.keys())}")

    c_info = CONSTRUCTS[construct_key]
    p_info = PROMOTERS[promoter_key]
    a_info = POLYA_SIGNALS[polya_key]

    itrs_bp = 2 * ITR_LENGTH_NT  # 290 bp
    cds_bp = c_info["cds_length_bp"]
    total = itrs_bp + p_info["length_nt"] + utr_bp + cds_bp + a_info["length_nt"]
    fits = total <= AAV_HARD_PACKAGING_LIMIT_NT
    headroom = AAV_HARD_PACKAGING_LIMIT_NT - total

    return MinigenePackagingBudget(
        construct_key=construct_key,
        construct_name=c_info["name"],
        promoter_key=promoter_key,
        promoter_bp=p_info["length_nt"],
        polya_key=polya_key,
        polya_bp=a_info["length_nt"],
        utr_bp=utr_bp,
        itrs_bp=itrs_bp,
        cds_bp=cds_bp,
        total_genome_nt=total,
        hard_limit_nt=AAV_HARD_PACKAGING_LIMIT_NT,
        fits_single_aav=fits,
        packaging_headroom_nt=headroom,
    )


def simulate_patch_clamp_electrophysiology(
    construct_key: str = "CFTR_DeltaR_708_759",
    pka_stimulated: bool = True,
    atp_concentration_mm: float = 1.0,
    apical_channels_n: int = 1000,
    v_membrane_mv: float = -50.0,
    e_cl_mv: float = -30.0,
) -> PatchClampSimulationResult:
    """Simulate single-channel and macroscopic patch-clamp chloride current."""
    if construct_key not in CONSTRUCTS:
        raise ValueError(f"Unknown construct: {construct_key}")

    c_info = CONSTRUCTS[construct_key]

    # Zero ATP strictly prevents gating opening (DOI 10.1073/pnas.1901174116)
    if atp_concentration_mm <= 0.0 or construct_key == "CFTR_Delta264":
        po = 0.0
    else:
        # Hill equation for ATP binding dependence (Kd ~ 0.15 mM, n ~ 1.5)
        atp_factor = (atp_concentration_mm ** 1.5) / (0.15 ** 1.5 + atp_concentration_mm ** 1.5)
        base_po = c_info["stimulated_po"] if pka_stimulated else c_info["basal_po"]
        po = min(1.0, base_po * atp_factor)

    # Driving force: (Vm - E_Cl) in millivolts
    driving_force_mv = v_membrane_mv - e_cl_mv
    driving_force_v = driving_force_mv * 1e-3

    # Unitary current: i = gamma * (Vm - E_Cl)
    # gamma = 8.5 pS = 8.5e-12 Siemens
    gamma_siemens = UNITARY_CONDUCTANCE_PS * 1e-12
    single_channel_i_amp = gamma_siemens * driving_force_v
    single_channel_i_pa = single_channel_i_amp * 1e12  # in picoamperes (pA)

    # Effective active surface channels: N_surface = N_total * fraction_band_c
    f_band_c = c_info["fraction_band_c_mature"]
    n_active = apical_channels_n * f_band_c

    # Macroscopic apical current: I_Cl = N_active * i * Po
    macroscopic_i_pa = n_active * single_channel_i_pa * po

    return PatchClampSimulationResult(
        construct_key=construct_key,
        pka_stimulated=pka_stimulated,
        atp_concentration_mm=round(atp_concentration_mm, 3),
        open_probability_po=round(po, 4),
        unitary_conductance_ps=UNITARY_CONDUCTANCE_PS,
        single_channel_current_pa=round(single_channel_i_pa, 4),
        apical_channels_n=apical_channels_n,
        fraction_band_c_mature=round(f_band_c, 3),
        macroscopic_current_pa=round(macroscopic_i_pa, 2),
        v_membrane_mv=round(v_membrane_mv, 1),
        e_cl_mv=round(e_cl_mv, 1),
        driving_force_mv=round(driving_force_mv, 1),
    )


def simulate_cftr_minigene(
    construct_key: str = "CFTR_DeltaR_708_759",
    promoter_key: str = "tg83",
    polya_key: str = "synPolyA_compact",
    pka_stimulated: bool = True,
    atp_concentration_mm: float = 1.0,
) -> MinigeneSimulationReport:
    """Run end-to-end minigene packaging and patch-clamp electrophysiology audit."""
    packaging = evaluate_minigene_packaging(
        construct_key=construct_key,
        promoter_key=promoter_key,
        polya_key=polya_key,
    )
    ephys = simulate_patch_clamp_electrophysiology(
        construct_key=construct_key,
        pka_stimulated=pka_stimulated,
        atp_concentration_mm=atp_concentration_mm,
    )

    return MinigeneSimulationReport(
        title="CFTR Delta-R Domain Minigene Packaging & Single-Channel Electrophysiology Simulator",
        nlm_attribution=NLM_ATTRIBUTION,
        research_boundary=RESEARCH_BOUNDARY,
        packaging=packaging,
        electrophysiology=ephys,
        structural_notes=CONSTRUCTS[construct_key]["structural_integrity"],
        relevant_nct_ids=["NCT05248230"],
        evidence_citations=[
            "DOI: 10.1073/pnas.1901174116",
            "ClinicalTrials.gov: NCT05248230",
            "PDB: 6MSM",
            "PDB: 5UAK",
        ],
    )


def self_test() -> bool:
    """Deterministic internal self-test suite."""
    # 1. DeltaR(708-759) with tg83 and synPolyA_compact must fit within 4,700 nt
    pkg_opt = evaluate_minigene_packaging("CFTR_DeltaR_708_759", "tg83", "synPolyA_compact")
    assert pkg_opt.fits_single_aav, "DeltaR(708-759) with tg83 must fit within 4,700 nt"
    assert pkg_opt.total_genome_nt <= AAV_HARD_PACKAGING_LIMIT_NT

    # 2. Full-length WT CFTR even with tg83 must exceed 4,700 nt
    pkg_wt = evaluate_minigene_packaging("CFTR_WT", "tg83", "synPolyA_compact")
    assert not pkg_wt.fits_single_aav, "Full-length WT CFTR must exceed 4,700 nt"

    # 3. Unitary conductance must be invariant at 8.5 pS
    for ck in ("CFTR_WT", "CFTR_DeltaR_708_759", "CFTR_DeltaR_708_835"):
        ep = simulate_patch_clamp_electrophysiology(ck)
        assert math.isclose(ep.unitary_conductance_ps, 8.5), f"{ck} unitary conductance must be 8.5 pS"

    # 4. In absence of ATP, Po must be 0.0
    ep_no_atp = simulate_patch_clamp_electrophysiology("CFTR_DeltaR_708_759", atp_concentration_mm=0.0)
    assert ep_no_atp.open_probability_po == 0.0, "Zero ATP must yield Po = 0.0"

    # 5. Delta264 control fails ER Band C maturation
    ep_d264 = simulate_patch_clamp_electrophysiology("CFTR_Delta264")
    assert ep_d264.fraction_band_c_mature < 0.05, "Delta264 must fail ER export"
    assert ep_d264.macroscopic_current_pa == 0.0

    # 6. Basal Po is elevated in DeltaR(708-759) compared to WT
    ep_wt_basal = simulate_patch_clamp_electrophysiology("CFTR_WT", pka_stimulated=False)
    ep_dr_basal = simulate_patch_clamp_electrophysiology("CFTR_DeltaR_708_759", pka_stimulated=False)
    assert ep_dr_basal.open_probability_po > ep_wt_basal.open_probability_po

    return True


def main() -> None:
    parser = argparse.ArgumentParser(
        description="CFTR Delta-R Domain Minigene Packaging & Single-Channel Conductance Simulator (T-0207)"
    )
    parser.add_argument(
        "--construct",
        type=str,
        default="CFTR_DeltaR_708_759",
        choices=list(CONSTRUCTS.keys()),
        help="Target CFTR minigene construct (default: CFTR_DeltaR_708_759)",
    )
    parser.add_argument(
        "--promoter",
        type=str,
        default="tg83",
        choices=list(PROMOTERS.keys()),
        help="Promoter architecture (default: tg83)",
    )
    parser.add_argument(
        "--polya",
        type=str,
        default="synPolyA_compact",
        choices=list(POLYA_SIGNALS.keys()),
        help="Polyadenylation signal (default: synPolyA_compact)",
    )
    parser.add_argument(
        "--atp",
        type=float,
        default=1.0,
        help="Intracellular ATP concentration in mM (default: 1.0)",
    )
    parser.add_argument(
        "--no-pka",
        action="store_true",
        help="Simulate basal unphosphorylated state (default: PKA stimulated)",
    )
    parser.add_argument("--json", action="store_true", help="Print structured JSON output")
    parser.add_argument("--export", type=str, help="Export JSON results to output path")
    parser.add_argument("--self-test", action="store_true", help="Run deterministic internal self-tests")
    args = parser.parse_args()

    if args.self_test:
        success = self_test()
        if success:
            print("PASS: cftr_minigene_model self-test passed.")
            sys.exit(0)
        else:
            print("FAIL: cftr_minigene_model self-test failed.")
            sys.exit(1)

    report = simulate_cftr_minigene(
        construct_key=args.construct,
        promoter_key=args.promoter,
        polya_key=args.polya,
        pka_stimulated=not args.no_pka,
        atp_concentration_mm=args.atp,
    )
    report_dict = asdict(report)

    if args.export:
        out_p = Path(args.export)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        with open(out_p, "w", encoding="utf-8") as f:
            json.dump(report_dict, f, indent=2)
        print(f"Exported minigene simulation to {out_p}")

    if args.json:
        print(json.dumps(report_dict, indent=2))
    else:
        print("=" * 80)
        print(report.title.upper())
        print(RESEARCH_BOUNDARY)
        print("=" * 80)
        pkg = report.packaging
        print(f"Construct: {pkg.construct_name} ({pkg.construct_key})")
        print(f"  Regulatory Elements: {pkg.promoter_key} ({pkg.promoter_bp} nt) + {pkg.polya_key} ({pkg.polya_bp} nt)")
        print(f"  Vector Genome: {pkg.total_genome_nt} nt (limit: {pkg.hard_limit_nt} nt | headroom: {pkg.packaging_headroom_nt} nt)")
        print(f"  Single-AAV Packaging Feasibility: {'FITS CEILING' if pkg.fits_single_aav else 'OVER CAPACITY'}")
        print(f"  Structural Coordination: {report.structural_notes}")
        ep = report.electrophysiology
        print(f"\nPatch-Clamp Electrophysiology Simulation (ATP = {ep.atp_concentration_mm} mM, PKA = {ep.pka_stimulated}):")
        print(f"  Unitary Conductance (gamma): {ep.unitary_conductance_ps} pS")
        print(f"  Open Probability (Po): {ep.open_probability_po:.4f}")
        print(f"  Single-Channel Current (i): {ep.single_channel_current_pa:.3f} pA (at Vm = {ep.v_membrane_mv} mV, ECl = {ep.e_cl_mv} mV)")
        print(f"  Mature Band C ER Exit: {ep.fraction_band_c_mature * 100:.1f}%")
        print(f"  Macroscopic Apical Cl- Current (I_Cl): {ep.macroscopic_current_pa:.1f} pA per {ep.apical_channels_n} channels")
        print(f"  Clinical Reference Trial: {', '.join(report.relevant_nct_ids)}")
        print("=" * 80)


if __name__ == "__main__":
    main()
