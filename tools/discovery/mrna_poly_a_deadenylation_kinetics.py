#!/usr/bin/env python3
"""Biophysical Modeling of Inhaled CFTR mRNA Poly(A) Tail Deadenylation Kinetics and Ribosomal Translocation Collision Protection (T-0449).

Quantitative biophysical modeling of:
1. Poly(A) Tail Length & Biphasic Deadenylation Kinetics:
   - Formulates biphasic deadenylation dynamics: initial slow PAN2-PAN3 trimming
     with PABP protection (L_A > 25 nt) followed by rapid terminal decay by
     the CCR4-NOT complex (L_A <= 25 nt threshold).
   - Evaluates chemical modification stabilization (100% N1-methylpseudouridine / m1Psi)
     and codon-optimality suppression of deadenylase recruitment.
2. eIF4F / PABP Translation Initiation Circularization:
   - Quantifies closed-loop mRNP circularization stability via eIF4E-eIF4G-PABP bridging.
   - Formulates cooperative ribosome recruitment frequency (initiation rate: 0.04–0.12 Hz)
     as a function of instantaneous poly(A) tail length.
3. Ribosome Translocation & Collision Clearance:
   - Evaluates Totally Asymmetric Simple Exclusion Process (TASEP) ribosome spacing
     across the 4,443 nt (1,480 codons) CFTR coding sequence.
   - Quantifies disome collision rates and No-Go Decay (NGD) clearance comparing
     wild-type sequence vs synonymous codon-optimized synthetic mRNA.
4. Therapeutic mRNA Functional Half-Life & Channel Yield:
   - Simulates total mature Band C CFTR channel yield per transcript and per cell
     for aerosolized LNP delivery (VX-522, NCT05668741) across poly(A) tail lengths
     (100–120 nt) to determine clinical curative threshold rescue (>= 10% WT density).

Primary Literature Grounding:
  - Subtelny AO et al. (2014) Nature, DOI: 10.1038/nature13007
    (Poly(A)-tail profiling reveals an embryonic switch in translational control)
  - Presnyak V et al. (2015) Cell, DOI: 10.1016/j.cell.2015.02.029
    (Codon optimality is a major determinant of mRNA stability)
  - Kariko K et al. (2008) Mol Ther, DOI: 10.1038/mt.2008.200
    (Incorporation of pseudouridine into mRNA yields superior non-immunogenic vector with increased translational capacity)

Clinical Registry Grounding:
  - NCT05668741: Phase 1/2 Study of Inhaled mRNA LNP VX-522 in Participants With Cystic Fibrosis.

Structural Grounding:
  - RCSB PDB 6GQ1: Structure of the mammalian 48S pre-initiation complex with eIF4F.
  - RCSB PDB 5UAK: Phosphorylated active human CFTR dimer conformation.

Attribution & Research Boundary:
  ClinicalTrials.gov does not evaluate safety or efficacy; listing on the registry does not denote government endorsement.
  Preclinical computational biophysical simulation; research, not medical advice.
  No clinical efficacy or therapeutic safety claim is made.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
import sys
from typing import Any, Dict, List, Optional, Tuple

SCRIPT_PATH = Path(__file__).resolve()

NLM_ATTRIBUTION = (
    "This product uses publicly available data from the U.S. National Library of Medicine (NLM), "
    "National Institutes of Health, Department of Health and Human Services; NLM is not responsible "
    "for the product and does not endorse or recommend this or any other product."
)

RESEARCH_BOUNDARY = (
    "Preclinical in-silico biophysical modeling of inhaled CFTR mRNA poly(A) tail deadenylation kinetics, "
    "closed-loop circularization, and ribosome collision clearance; research, not medical advice. "
    "No clinical efficacy or therapeutic safety claim is made."
)

# Primary Literature & Registry Identifiers
DOI_POLYA_PROFILING = "10.1038/nature13007"
DOI_CODON_OPTIMALITY = "10.1016/j.cell.2015.02.029"
DOI_PSEUDOURIDINE_MRNA = "10.1038/mt.2008.200"

NCT_VX522_MRNA_TRIAL = "NCT05668741"

PDB_INITIATION_COMPLEX = "6GQ1"
PDB_DEPHOSPHO_CFTR = "6MSM"
PDB_ACTIVE_CFTR = "5UAK"

# Biophysical Constants
CFTR_CDS_NUCLEOTIDES: int = 4443
CFTR_CODONS: int = 1480
RIBOSOME_FOOTPRINT_CODONS: int = 10  # ~30 nucleotides occupied by 80S ribosome

# Poly(A) Tail & Deadenylation Constants
PABP_FOOTPRINT_NT: int = 27  # nucleotides covered per PABP monomer
CCR4_NOT_TRANSITION_THRESHOLD_NT: float = 25.0  # length below which PABP dissociates and CCR4-NOT accelerates

# Deadenylation rates in nucleotides / hour
K_DEADENYLATION_PAN2_WT: float = 3.20  # slow phase PAN2-PAN3 trimming rate (nt/h)
K_DEADENYLATION_CCR4_WT: float = 16.50  # rapid terminal decay CCR4-NOT rate (nt/h)

# Chemical modification & codon optimality stabilization multipliers (divisors of deadenylation rate)
M1PSI_STABILIZATION_FACTOR: float = 2.40  # 100% N1-methylpseudouridine stabilization
CODON_OPTIMALITY_STABILIZATION_FACTOR: float = 1.45  # high codon optimality (CAI ~0.95 vs 0.72)

# Translation Initiation Constants (Hz = ribosomes / second)
K_INIT_MAX_HZ: float = 0.100  # maximum initiation rate with intact closed-loop mRNP
K_INIT_BASAL_HZ: float = 0.012  # basal initiation rate of uncircularized/dead-tailed mRNA
CIRCULARIZATION_K_HALF_NT: float = 35.0  # poly(A) length for half-maximal closed-loop initiation
CIRCULARIZATION_HILL_N: float = 2.60  # cooperativity of PABP-eIF4G bridging

# Elongation Velocity Constants (codons / second)
V_ELONG_WILDTYPE_MEAN: float = 4.10  # mean elongation rate across wild-type CFTR CDS
V_ELONG_WILDTYPE_PAUSE: float = 1.70  # pause rate at rare codon clusters / domain boundaries
V_ELONG_OPTIMIZED_MEAN: float = 5.80  # mean elongation rate for synonymous codon-optimized CDS
V_ELONG_OPTIMIZED_PAUSE: float = 4.40  # pause rate for codon-optimized synthetic construct

# Folding, Maturation & Channel Yield Constants
CFTR_COTRANSLATIONAL_FOLDING_EFFICIENCY: float = 0.280  # mature Band C ER-to-surface fraction
CFTR_SURFACE_HALF_LIFE_HOURS: float = 20.00  # surface channel turnover half-life (~72,000 s)
WT_CFTR_CHANNELS_PER_CELL: float = 1500.0  # wild-type ciliated human epithelial channel density
CURATIVE_RESCUE_FRACTION_WT: float = 0.100  # >= 10% WT density clinical threshold (150 channels/cell)
CURATIVE_RESCUE_CHANNELS: float = WT_CFTR_CHANNELS_PER_CELL * CURATIVE_RESCUE_FRACTION_WT


@dataclass(frozen=True)
class DeadenylationState:
    """State of mRNA poly(A) tail at time t."""

    time_hours: float
    poly_a_length_nt: float
    pabp_molecules_bound: int
    deadenylation_phase: str
    instantaneous_decay_rate_nt_per_h: float
    circularization_probability: float
    effective_initiation_rate_hz: float


@dataclass(frozen=True)
class TASEPCollisionProfile:
    """Ribosomal translocation and disome collision profile."""

    mean_elongation_velocity_codons_per_s: float
    effective_initiation_rate_hz: float
    mean_ribosome_spacing_codons: float
    mean_ribosome_density_per_mrna: float
    disome_collision_rate_per_hour: float
    no_go_decay_probability_per_transcript: float
    transit_time_seconds: float


@dataclass(frozen=True)
class ChannelYieldResult:
    """Simulated CFTR channel synthesis and surface accumulation."""

    initial_poly_a_length_nt: int
    codon_optimized: bool
    m1psi_modified: bool
    functional_half_life_hours: float
    total_translations_per_mrna: float
    cytosolic_mrna_copies_per_cell: int
    steady_state_surface_channels_per_cell: float
    percent_wild_type_density: float
    curative_threshold_achieved: bool


class MRNAPolyADeadenylationEngine:
    """Quantitative biophysical model of CFTR mRNA deadenylation, initiation, and TASEP translocation."""

    def __init__(
        self,
        cftr_cds_nt: int = CFTR_CDS_NUCLEOTIDES,
        cftr_codons: int = CFTR_CODONS,
        ribosome_footprint: int = RIBOSOME_FOOTPRINT_CODONS,
    ) -> None:
        self.cftr_cds_nt = cftr_cds_nt
        self.cftr_codons = cftr_codons
        self.ribosome_footprint = ribosome_footprint

    def compute_deadenylation_rate(
        self,
        poly_a_length_nt: float,
        codon_optimized: bool = True,
        m1psi_modified: bool = True,
    ) -> Tuple[float, str]:
        """Compute instantaneous deadenylation rate (nt/h) and active phase."""
        if poly_a_length_nt <= 0.0:
            return 0.0, "TERMINATED"

        stabilization = 1.0
        if m1psi_modified:
            stabilization *= M1PSI_STABILIZATION_FACTOR
        if codon_optimized:
            stabilization *= CODON_OPTIMALITY_STABILIZATION_FACTOR

        if poly_a_length_nt > CCR4_NOT_TRANSITION_THRESHOLD_NT:
            # Slow PAN2-PAN3 trimming with PABP protection
            base_rate = K_DEADENYLATION_PAN2_WT
            phase = "PAN2_PAN3_PROTECTED"
        else:
            # Rapid terminal CCR4-NOT deadenylation below PABP binding threshold
            base_rate = K_DEADENYLATION_CCR4_WT
            phase = "CCR4_NOT_RAPID_DECAY"

        effective_rate = base_rate / stabilization
        return effective_rate, phase

    def compute_circularization_and_initiation(
        self,
        poly_a_length_nt: float,
    ) -> Tuple[float, float, int]:
        """Compute closed-loop circularization probability and effective initiation rate (Hz)."""
        if poly_a_length_nt <= 0.0:
            return 0.0, K_INIT_BASAL_HZ * 0.20, 0

        pabp_bound = int(poly_a_length_nt // PABP_FOOTPRINT_NT)

        # Cooperative circularization through eIF4G-PABP bridging (Hill function)
        n = CIRCULARIZATION_HILL_N
        k_half = CIRCULARIZATION_K_HALF_NT
        circ_prob = (poly_a_length_nt**n) / ((poly_a_length_nt**n) + (k_half**n))

        # Effective initiation rate interpolation between basal and maximum
        k_init = K_INIT_BASAL_HZ + (K_INIT_MAX_HZ - K_INIT_BASAL_HZ) * circ_prob
        return circ_prob, k_init, pabp_bound

    def simulate_deadenylation_trajectory(
        self,
        initial_length_nt: int = 120,
        codon_optimized: bool = True,
        m1psi_modified: bool = True,
        time_step_hours: float = 0.25,
        max_duration_hours: float = 160.0,
    ) -> List[DeadenylationState]:
        """Simulate time-course trajectory of poly(A) tail shortening and initiation capacity."""
        trajectory: List[DeadenylationState] = []
        current_len = float(initial_length_nt)
        t = 0.0

        while t <= max_duration_hours and current_len > 0.0:
            circ_prob, k_init, pabp = self.compute_circularization_and_initiation(current_len)
            rate, phase = self.compute_deadenylation_rate(
                current_len,
                codon_optimized=codon_optimized,
                m1psi_modified=m1psi_modified,
            )

            trajectory.append(
                DeadenylationState(
                    time_hours=round(t, 2),
                    poly_a_length_nt=round(current_len, 2),
                    pabp_molecules_bound=pabp,
                    deadenylation_phase=phase,
                    instantaneous_decay_rate_nt_per_h=round(rate, 3),
                    circularization_probability=round(circ_prob, 4),
                    effective_initiation_rate_hz=round(k_init, 4),
                )
            )

            # Advance Euler step
            delta_len = rate * time_step_hours
            current_len = max(0.0, current_len - delta_len)
            t += time_step_hours

        return trajectory

    def calculate_functional_half_life(
        self,
        initial_length_nt: int = 120,
        codon_optimized: bool = True,
        m1psi_modified: bool = True,
    ) -> float:
        """Calculate functional mRNA half-life (hours until poly(A) decays below 25 nt threshold)."""
        trajectory = self.simulate_deadenylation_trajectory(
            initial_length_nt=initial_length_nt,
            codon_optimized=codon_optimized,
            m1psi_modified=m1psi_modified,
            time_step_hours=0.10,
            max_duration_hours=160.0,
        )

        for state in trajectory:
            if state.poly_a_length_nt <= CCR4_NOT_TRANSITION_THRESHOLD_NT:
                return state.time_hours

        return trajectory[-1].time_hours if trajectory else 0.0

    def evaluate_tasep_ribosome_dynamics(
        self,
        k_init_hz: float,
        codon_optimized: bool = True,
    ) -> TASEPCollisionProfile:
        """Evaluate TASEP elongation, ribosome density, spacing, and disome collision probability."""
        if codon_optimized:
            v_mean = V_ELONG_OPTIMIZED_MEAN
            v_pause = V_ELONG_OPTIMIZED_PAUSE
        else:
            v_mean = V_ELONG_WILDTYPE_MEAN
            v_pause = V_ELONG_WILDTYPE_PAUSE

        # Average spacing between ribosomes (codons)
        spacing_codons = v_mean / max(0.001, k_init_hz)

        # Total transit time through 1,480 codons
        transit_time_s = float(self.cftr_codons) / v_mean

        # Steady-state number of ribosomes on the transcript
        mean_ribosomes = k_init_hz * transit_time_s

        # Collision probability at pause sites
        pause_ratio = v_pause / v_mean
        headway_ratio = spacing_codons / float(self.ribosome_footprint)

        # Disome formation rate per transcript per hour (exponential queuing approximation)
        if headway_ratio > 1.2:
            disome_rate_per_hour = 3600.0 * (k_init_hz**2) * ((1.0 - pause_ratio) ** 2) / (headway_ratio**2)
        else:
            disome_rate_per_hour = 3600.0 * k_init_hz * (1.0 - pause_ratio)

        # No-Go Decay (NGD) probability over transcript transit
        k_znf598_per_s = 0.0005
        ngd_prob = min(0.95, (disome_rate_per_hour / 3600.0) * k_znf598_per_s * transit_time_s)

        return TASEPCollisionProfile(
            mean_elongation_velocity_codons_per_s=round(v_mean, 2),
            effective_initiation_rate_hz=round(k_init_hz, 4),
            mean_ribosome_spacing_codons=round(spacing_codons, 1),
            mean_ribosome_density_per_mrna=round(mean_ribosomes, 2),
            disome_collision_rate_per_hour=round(disome_rate_per_hour, 4),
            no_go_decay_probability_per_transcript=round(ngd_prob, 5),
            transit_time_seconds=round(transit_time_s, 1),
        )

    def calculate_total_translations_per_mrna(
        self,
        initial_length_nt: int = 120,
        codon_optimized: bool = True,
        m1psi_modified: bool = True,
    ) -> float:
        """Integrate total polypeptide chains produced over the lifetime of a single mRNA transcript."""
        trajectory = self.simulate_deadenylation_trajectory(
            initial_length_nt=initial_length_nt,
            codon_optimized=codon_optimized,
            m1psi_modified=m1psi_modified,
            time_step_hours=0.20,
            max_duration_hours=160.0,
        )

        total_translations = 0.0
        dt_seconds = 0.20 * 3600.0

        for state in trajectory:
            tasep = self.evaluate_tasep_ribosome_dynamics(
                state.effective_initiation_rate_hz,
                codon_optimized=codon_optimized,
            )
            effective_synthesis_hz = state.effective_initiation_rate_hz * (
                1.0 - tasep.no_go_decay_probability_per_transcript
            )
            total_translations += effective_synthesis_hz * dt_seconds

        return total_translations

    def simulate_channel_yield(
        self,
        initial_poly_a_length_nt: int = 120,
        cytosolic_mrna_copies: int = 15,
        codon_optimized: bool = True,
        m1psi_modified: bool = True,
    ) -> ChannelYieldResult:
        """Simulate functional CFTR channel yield and surface density per cell."""
        half_life = self.calculate_functional_half_life(
            initial_length_nt=initial_poly_a_length_nt,
            codon_optimized=codon_optimized,
            m1psi_modified=m1psi_modified,
        )

        total_trans = self.calculate_total_translations_per_mrna(
            initial_length_nt=initial_poly_a_length_nt,
            codon_optimized=codon_optimized,
            m1psi_modified=m1psi_modified,
        )

        if cytosolic_mrna_copies <= 0:
            surface_channels = 0.0
        else:
            # 24-hour therapeutic dosing window synthesis and steady-state surface channels
            traj_24h = self.simulate_deadenylation_trajectory(
                initial_length_nt=initial_poly_a_length_nt,
                codon_optimized=codon_optimized,
                m1psi_modified=m1psi_modified,
                time_step_hours=0.20,
                max_duration_hours=24.0,
            )
            dt_seconds = 0.20 * 3600.0
            trans_24h = 0.0
            for state in traj_24h:
                if state.time_hours > 24.0:
                    break
                tasep = self.evaluate_tasep_ribosome_dynamics(
                    state.effective_initiation_rate_hz,
                    codon_optimized=codon_optimized,
                )
                effective_hz = state.effective_initiation_rate_hz * (
                    1.0 - tasep.no_go_decay_probability_per_transcript
                )
                trans_24h += effective_hz * dt_seconds

            mature_channels_24h = float(cytosolic_mrna_copies) * trans_24h * CFTR_COTRANSLATIONAL_FOLDING_EFFICIENCY
            r_synth_per_hour = mature_channels_24h / 24.0
            k_clearance_per_h = math.log(2.0) / CFTR_SURFACE_HALF_LIFE_HOURS
            surface_channels = r_synth_per_hour / k_clearance_per_h

        pct_wt = (surface_channels / WT_CFTR_CHANNELS_PER_CELL) * 100.0
        curative_achieved = surface_channels >= CURATIVE_RESCUE_CHANNELS

        return ChannelYieldResult(
            initial_poly_a_length_nt=initial_poly_a_length_nt,
            codon_optimized=codon_optimized,
            m1psi_modified=m1psi_modified,
            functional_half_life_hours=round(half_life, 2),
            total_translations_per_mrna=round(total_trans, 1),
            cytosolic_mrna_copies_per_cell=cytosolic_mrna_copies,
            steady_state_surface_channels_per_cell=round(surface_channels, 1),
            percent_wild_type_density=round(pct_wt, 2),
            curative_threshold_achieved=curative_achieved,
        )

    def run_comprehensive_biophysical_suite(self) -> Dict[str, Any]:
        """Execute full 4-domain biophysical evaluation suite across tail lengths and modifications."""
        # Domain 1 & 2: Poly(A) tail length sweep (0 to 140 nt)
        tail_lengths = [0, 20, 40, 60, 80, 100, 120, 140]
        initiation_profile = []
        for l_a in tail_lengths:
            circ, k_init, pabp = self.compute_circularization_and_initiation(float(l_a))
            rate_wt, _ = self.compute_deadenylation_rate(float(l_a), codon_optimized=False, m1psi_modified=False)
            rate_opt, _ = self.compute_deadenylation_rate(float(l_a), codon_optimized=True, m1psi_modified=True)
            initiation_profile.append({
                "poly_a_length_nt": l_a,
                "pabp_bound": pabp,
                "circularization_probability": round(circ, 4),
                "initiation_rate_hz": round(k_init, 4),
                "deadenylation_rate_wt_nt_h": round(rate_wt, 2),
                "deadenylation_rate_opt_nt_h": round(rate_opt, 2),
            })

        # Domain 3: TASEP ribosome spacing comparison
        tasep_wt = self.evaluate_tasep_ribosome_dynamics(0.080, codon_optimized=False)
        tasep_opt = self.evaluate_tasep_ribosome_dynamics(0.080, codon_optimized=True)

        # Domain 4: Formulation comparison & Dose-Response
        doses = [2, 5, 10, 15, 20, 30]
        dose_response_unmod = [asdict(self.simulate_channel_yield(100, d, codon_optimized=False, m1psi_modified=False)) for d in doses]
        dose_response_opt = [asdict(self.simulate_channel_yield(120, d, codon_optimized=True, m1psi_modified=True)) for d in doses]

        return {
            "metadata": {
                "task": "T-0449",
                "title": "Biophysical Modeling of Inhaled CFTR mRNA Poly(A) Tail Deadenylation Kinetics and Ribosomal Translocation Collision Protection",
                "track": "Track A",
                "discipline": "rna_gene_editing",
                "primary_literature": [
                    DOI_POLYA_PROFILING,
                    DOI_CODON_OPTIMALITY,
                    DOI_PSEUDOURIDINE_MRNA,
                ],
                "clinical_registry": NCT_VX522_MRNA_TRIAL,
                "structural_grounding": [PDB_INITIATION_COMPLEX, PDB_ACTIVE_CFTR],
                "nlm_attribution": NLM_ATTRIBUTION,
                "research_boundary": RESEARCH_BOUNDARY,
            },
            "domain_1_and_2_deadenylation_and_circularization": {
                "tail_length_sweep": initiation_profile,
                "ccr4_transition_threshold_nt": CCR4_NOT_TRANSITION_THRESHOLD_NT,
                "pan2_vs_ccr4_rate_ratio": round(K_DEADENYLATION_CCR4_WT / K_DEADENYLATION_PAN2_WT, 2),
            },
            "domain_3_tasep_ribosome_dynamics": {
                "unoptimized_wt": asdict(tasep_wt),
                "codon_optimized_m1psi": asdict(tasep_opt),
                "collision_rate_reduction_pct": round(
                    (1.0 - tasep_opt.disome_collision_rate_per_hour / max(0.001, tasep_wt.disome_collision_rate_per_hour)) * 100.0,
                    2,
                ),
            },
            "domain_4_channel_yield_and_curative_rescue": {
                "unmodified_100nt_dose_response": dose_response_unmod,
                "codon_optimized_m1psi_120nt_dose_response": dose_response_opt,
                "minimum_curative_dose_mRNAs_per_cell": next(
                    (d["cytosolic_mrna_copies_per_cell"] for d in dose_response_opt if d["curative_threshold_achieved"]),
                    None,
                ),
            },
        }


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Biophysical Modeling of Inhaled CFTR mRNA Poly(A) Tail Deadenylation Kinetics (T-0449)."
    )
    parser.add_argument(
        "--poly-a-len",
        type=int,
        default=None,
        help="Initial poly(A) tail length in nucleotides (e.g. 100, 120).",
    )
    parser.add_argument(
        "--mrna-copies",
        type=int,
        default=None,
        help="Cytosolic delivered CFTR mRNA copies per cell (e.g. 5 to 30).",
    )
    parser.add_argument(
        "--unmodified",
        action="store_true",
        help="Simulate unmodified/unoptimized mRNA.",
    )
    parser.add_argument(
        "--run-all",
        action="store_true",
        help="Execute comprehensive 4-domain biophysical evaluation suite.",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Output results in structured machine-readable JSON format.",
    )

    args = parser.parse_args(argv)
    engine = MRNAPolyADeadenylationEngine()

    if args.run_all or (args.poly_a_len is None and args.mrna_copies is None):
        results = engine.run_comprehensive_biophysical_suite()
        if args.json:
            print(json.dumps(results, indent=2))
        else:
            print("=" * 76)
            print("CFTR MRNA POLY(A) TAIL DEADENYLATION & RIBOSOME COLLISION ENGINE (T-0449)")
            print("=" * 76)
            print(f"Primary Registry: {NCT_VX522_MRNA_TRIAL} (Inhaled mRNA LNP VX-522)")
            print(f"Primary Literature: {DOI_POLYA_PROFILING}, {DOI_CODON_OPTIMALITY}, {DOI_PSEUDOURIDINE_MRNA}")
            print(f"Structural Grounding: PDB {PDB_INITIATION_COMPLEX}, {PDB_ACTIVE_CFTR}")

            d1 = results["domain_1_and_2_deadenylation_and_circularization"]["tail_length_sweep"]
            print("\n[Domain 1 & 2] Poly(A) Tail Length & Translation Initiation Circularization:")
            for p in d1[2:7]:
                print(f"  Length: {p['poly_a_length_nt']:3d} nt | PABP: {p['pabp_bound']} | Circ Prob: {p['circularization_probability']:.3f} | k_init: {p['initiation_rate_hz']:.4f} Hz | Decay: {p['deadenylation_rate_opt_nt_h']:.2f} nt/h")

            d3 = results["domain_3_tasep_ribosome_dynamics"]
            print(f"\n[Domain 3] Ribosome Collision Avoidance (TASEP Dynamics):")
            print(f"  Unoptimized WT Collision Rate:      {d3['unoptimized_wt']['disome_collision_rate_per_hour']:.4f} disomes/h")
            print(f"  Codon-Optimized m1Psi Rate:         {d3['codon_optimized_m1psi']['disome_collision_rate_per_hour']:.4f} disomes/h")
            print(f"  Collision Suppression:              {d3['collision_rate_reduction_pct']:.1f}% reduction")

            d4 = results["domain_4_channel_yield_and_curative_rescue"]["codon_optimized_m1psi_120nt_dose_response"]
            print(f"\n[Domain 4] Mature CFTR Channel Density per Cell (Optimized 120 nt Poly(A)):")
            for r in d4:
                status = "PASS (Curative)" if r["curative_threshold_achieved"] else "SUB-CURATIVE"
                print(f"  {r['cytosolic_mrna_copies_per_cell']:2d} mRNAs/cell | {r['steady_state_surface_channels_per_cell']:6.1f} channels ({r['percent_wild_type_density']:4.1f}% WT) -> {status}")
            print("=" * 76)
        return 0

    l_a = args.poly_a_len if args.poly_a_len is not None else 120
    copies = args.mrna_copies if args.mrna_copies is not None else 15
    res = engine.simulate_channel_yield(
        initial_poly_a_length_nt=l_a,
        cytosolic_mrna_copies=copies,
        codon_optimized=not args.unmodified,
        m1psi_modified=not args.unmodified,
    )
    print(json.dumps(asdict(res), indent=2) if args.json else res)
    return 0


if __name__ == "__main__":
    sys.exit(main())
