#!/usr/bin/env python3
"""Biophysical Modeling of Inhaled CFTR mRNA Poly(A) Tail Deadenylation Kinetics and Ribosomal Translocation Collision Protection (T-0449).

Quantitative biophysical dynamical modeling of synthetic CFTR mRNA longevity and
translational efficiency in human airway epithelial cells following lipid nanoparticle
(LNP) endosomal escape:

1. Poly(A) Tail Length & Biphasic Deadenylation Kinetics:
   - Formulates sequential two-phase deadenylation of the 3' poly(A) tail:
     Phase I (PAN2-PAN3 distributive trimming): When poly(A) tail length L_A(t) > L_crit (25 nt)
     and multiple poly(A)-binding protein 1 (PABP1) monomers (footprint 27 nt) occupy the tail,
     PAN2-PAN3 trims distal adenosines at a slow protected rate k_PAN2 (approx 2.2-3.5 nt/hr).
     Segmented poly(A) linkers (e.g., A60-GCAUAUGACU-A60) sterically impede exonuclease
     processivity, further slowing deadenylation.
     Phase II (CCR4-NOT processive terminal decay, PDB 6GQ1): Once L_A(t) drops below L_crit
     (25 nt, loss of the last PABP1 monomer), CCR4-NOT catalytic subunits (CNOT6L/CNOT7)
     rapidly deadenylate the residual tail (k_CCR4 approx 18.0 nt/hr), triggering DCP1/DCP2
     5' decapping and XRN1 5'->3' exonucleolytic decay.

2. eIF4F / PABP1 Translation Initiation Circularization:
   - Models cooperative PABP1 multimerization (n_PABP = floor(L_A / 27)) and high-affinity
     binding to the eIF4G scaffold of the 5' cap-binding eIF4F complex.
   - Closed-loop circularization probability P_circ(L_A) follows a Hill-cooperative partition
     that sustains ribosome 43S pre-initiation recruitment (k_init in [0.04, 0.12] Hz) while
     sterically shielding the 5' m7G cap from DCP2 decapping.

3. TASEP Ribosome Translocation & Collision Clearance:
   - Evaluates 1D Totally Asymmetric Simple Exclusion Process (TASEP) ribosome traffic across
     the 4,443 nt CFTR coding sequence (1,480 amino acids + stop codon; ribosome footprint
     L_ribo = 30 nt = 10 codons).
   - Quantifies ribosome density rho = alpha / (alpha + beta_eff), di-ribosome collision
     probability P_coll, and ZNF598/Hel2-mediated No-Go Decay (NGD) endonucleolytic cleavage
     when initiation rate k_init exceeds the critical queuing threshold relative to codon-optimized
     and N1-methylpseudouridine (m1Psi) elongation velocity (v_elong in [4.2, 6.8] codons/s).

4. Therapeutic mRNA Functional Half-Life & 24-Hour CFTR Channel Yield:
   - Integrates coupled deadenylation, NGD collision surveillance, and circularized polysome
     translation over a 24-hour therapeutic window to compute functional full-length CFTR
     channel yield (PDB 5UAK) per cytoplasmic mRNA transcript.

Primary Literature Grounding:
  - Subtelny AO et al. (2014) Nature, DOI: 10.1038/nature13007
    (Poly(A)-tail profiling reveals an embryonic switch in translational control)
  - Presnyak V et al. (2015) Cell, DOI: 10.1016/j.cell.2015.02.029
    (Codon optimality is a major determinant of mRNA stability)
  - Kariko K et al. (2008) Mol Ther, DOI: 10.1038/mt.2008.200
    (Incorporation of pseudouridine into mRNA yields superior nonimmunogenic vector with increased translational capacity)

Clinical Registry Grounding:
  - NCT05668741 (Phase 1/2 Study of Inhaled mRNA LNP VX-522 in Participants With Cystic Fibrosis)

Structural Grounding:
  - RCSB PDB 6GQ1 (Human CCR4-NOT deadenylase complex architecture)
  - RCSB PDB 5UAK (Phosphorylated ATP-bound human CFTR conformation)

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
import sys
from typing import Any, Dict, List, Optional, Tuple

# Mandatory ClinicalTrials.gov Attribution & Research Boundary
NLM_ATTRIBUTION: str = (
    "ClinicalTrials.gov does not evaluate safety or efficacy; "
    "listing on the registry does not denote government endorsement."
)

RESEARCH_BOUNDARY: str = (
    "Preclinical computational biophysical simulation; research, not medical advice. "
    "No clinical efficacy or therapeutic safety claim is made."
)

# Registry & Structural Anchors (AH-1, AH-3, AH-6)
NCT_VX522_MRNA_TRIAL: str = "NCT05668741"
PDB_CCR4_NOT_COMPLEX: str = "6GQ1"
PDB_ACTIVE_CFTR: str = "5UAK"

# Primary DOIs (AH-6)
DOI_POLYA_PROFILING: str = "10.1038/nature13007"
DOI_CODON_OPTIMALITY: str = "10.1016/j.cell.2015.02.029"
DOI_PSEUDOURIDINE_MRNA: str = "10.1038/mt.2008.200"

# CFTR Coding Sequence & Ribosome Biophysical Constants
CFTR_CDS_LENGTH_NT: int = 4443  # 1,480 amino acids + stop codon (4,443 nt)
CFTR_CDS_CODONS: int = 1480  # Functional CFTR polypeptide length
RIBOSOME_FOOTPRINT_NT: int = 30  # Eukaryotic 80S ribosome footprint (~30 nt = 10 codons)
RIBOSOME_FOOTPRINT_CODONS: int = 10

PABP1_FOOTPRINT_NT: int = 27  # Adenosines required per PABP1 monomer
CRITICAL_POLYA_THRESHOLD_NT: float = 25.0  # Threshold below which last PABP1 dissociates & CCR4-NOT accelerates
OPTIMAL_POLYA_TAIL_NT: float = 110.0  # Therapeutic synthetic mRNA poly(A) tail length (100-120 nt)

# Deadenylation Rate Constants (nt / hour)
BASE_PAN2_PAN3_RATE_NT_HR: float = 3.6  # Phase I distributive trimming rate (PAN2-PAN3)
CCR4_NOT_TERMINAL_RATE_NT_HR: float = 18.5  # Phase II rapid terminal decay rate (CCR4-NOT, PDB 6GQ1)
SEGMENTED_LINKER_PROTECTION_FACTOR: float = 0.62  # Reduction in deadenylation velocity via segmented poly(A) spacer

# Translation Initiation & Elongation Parameters
MIN_INITIATION_RATE_HZ: float = 0.04  # Baseline ribosome recruitment frequency (s^-1)
OPTIMAL_INITIATION_RATE_HZ: float = 0.085  # Balanced high-output initiation rate without collision stalling (s^-1)
MAX_INITIATION_RATE_HZ: float = 0.12  # Upper physiological initiation rate before severe queuing (s^-1)

UNMODIFIED_ELONGATION_CODONS_S: float = 4.1  # Elongation rate of unmodified uridine mRNA (PKR/OAS drag)
M1PSI_OPTIMIZED_ELONGATION_CODONS_S: float = 6.2  # Elongation rate with m1Psi + high codon optimality (Kariko 2008; Presnyak 2015)

CO_TRANSLATIONAL_FOLDING_EFFICIENCY: float = 0.72  # Fraction of completed nascent chains folding into active CFTR (PDB 5UAK)
CURATIVE_CHANNEL_YIELD_THRESHOLD: float = 850.0  # Target functional CFTR channels produced per mRNA over 24h


@dataclass(frozen=True)
class DeadenylationStatePoint:
    """Single time-step state along the 24-hour mRNA deadenylation and translation trajectory."""

    time_hr: float
    polya_length_nt: float
    bound_pabp_monomers: int
    closed_loop_probability: float
    effective_initiation_hz: float
    deadenylation_phase: str
    intact_transcript_fraction: float
    cumulative_cftr_channels: float


@dataclass(frozen=True)
class TASEPRibosomeCollisionMetrics:
    """1D TASEP ribosomal spacing, queuing density, and No-Go Decay (NGD) collision surveillance metrics."""

    initiation_rate_hz: float
    elongation_velocity_codons_s: float
    max_lattice_current_hz: float
    ribosome_coverage_density: float
    mean_inter_ribosome_spacing_nt: float
    active_polysome_ribosomes: float
    collision_probability: float
    ngd_cleavage_rate_hr: float
    steady_translation_flux_hz: float
    traffic_regime: str


@dataclass(frozen=True)
class PolyATrajectoryResult:
    """Complete 24-hour biophysical simulation summary for a synthetic CFTR mRNA design."""

    design_name: str
    initial_polya_nt: float
    m1psi_modified: bool
    segmented_linker: bool
    codon_optimality_index: float
    target_initiation_hz: float
    functional_half_life_hr: float
    time_to_critical_threshold_hr: float
    peak_polysome_size: float
    collision_probability: float
    total_24h_cftr_channels: float
    curative_yield_achieved: bool
    trajectory: List[DeadenylationStatePoint]


class MRNAPolyADeadenylationEngine:
    """Biophysical engine modeling CFTR mRNA poly(A) deadenylation, PABP circularization, and TASEP ribosome traffic."""

    def __init__(
        self,
        cds_length_nt: int = CFTR_CDS_LENGTH_NT,
        pabp_footprint_nt: int = PABP1_FOOTPRINT_NT,
        critical_polya_nt: float = CRITICAL_POLYA_THRESHOLD_NT,
    ) -> None:
        if cds_length_nt <= 0:
            raise ValueError("CDS length must be positive")
        if pabp_footprint_nt <= 0:
            raise ValueError("PABP1 footprint must be positive")
        if critical_polya_nt <= 0:
            raise ValueError("Critical poly(A) threshold must be positive")

        self.cds_length_nt = cds_length_nt
        self.cds_codons = cds_length_nt // 3
        self.pabp_footprint_nt = pabp_footprint_nt
        self.critical_polya_nt = critical_polya_nt

    def calculate_pabp_occupancy(self, polya_length_nt: float) -> int:
        """Compute number of PABP1 monomers bound to a poly(A) tail of length L_A."""
        if polya_length_nt < 0.0:
            raise ValueError("Poly(A) tail length cannot be negative")
        if polya_length_nt < 12.0:
            return 0
        # Minimum 12 nt for weak single monomer, 27 nt per cooperative PABP1 unit
        if polya_length_nt < self.pabp_footprint_nt:
            return 1
        return int(polya_length_nt // self.pabp_footprint_nt)

    def calculate_closed_loop_probability(
        self,
        polya_length_nt: float,
        m1psi_modified: bool = True,
    ) -> float:
        """Evaluate eIF4F-eIF4G-PABP1 closed-loop circularization probability P_circ(L_A).

        Follows cooperative Hill binding of PABP1 multimers to eIF4G (Subtelny et al. 2014).
        """
        if polya_length_nt < 0.0:
            raise ValueError("Poly(A) tail length cannot be negative")
        if polya_length_nt <= 5.0:
            return 0.02

        n_pabp = self.calculate_pabp_occupancy(polya_length_nt)
        # Cooperative stability scale: half-maximal circularization near L_half = 38 nt (~1.4 PABP units)
        l_half = 36.0 if m1psi_modified else 44.0
        hill_n = 2.4
        coop_bonus = 1.0 + 0.06 * min(n_pabp, 4)
        raw_prob = ((polya_length_nt ** hill_n) / ((l_half ** hill_n) + (polya_length_nt ** hill_n))) * coop_bonus
        return round(min(0.98, max(0.02, raw_prob)), 4)

    def calculate_biphasic_deadenylation_rate(
        self,
        polya_length_nt: float,
        codon_optimality_index: float = 0.85,
        m1psi_modified: bool = True,
        segmented_linker: bool = True,
    ) -> Tuple[float, str]:
        """Compute instantaneous poly(A) shortening rate (-dL_A/dt in nt/hr) and active phase.

        Phase I (L_A >= 25 nt): Slow PAN2-PAN3 distributive trimming, stabilized by PABP1
        multimer protection and high codon optimality (Presnyak et al. 2015).
        Phase II (L_A < 25 nt): Rapid CCR4-NOT processive terminal deadenylation (PDB 6GQ1).
        """
        if polya_length_nt < 0.0:
            raise ValueError("Poly(A) tail length cannot be negative")
        if not (0.0 <= codon_optimality_index <= 1.0):
            raise ValueError("Codon optimality index must be in [0.0, 1.0]")

        if polya_length_nt <= 0.0:
            return 0.0, "TERMINAL_DECAPPED"

        # Codon optimality slows Dhh1/CCR4-NOT recruitment (Presnyak et al. 2015)
        codon_protection = 1.35 - 0.55 * codon_optimality_index
        m1psi_factor = 0.82 if m1psi_modified else 1.25
        linker_factor = SEGMENTED_LINKER_PROTECTION_FACTOR if segmented_linker else 1.0

        if polya_length_nt >= self.critical_polya_nt:
            n_pabp = max(1, self.calculate_pabp_occupancy(polya_length_nt))
            pabp_shield = 1.0 / (1.0 + 0.14 * (n_pabp - 1))
            rate = BASE_PAN2_PAN3_RATE_NT_HR * codon_protection * m1psi_factor * linker_factor * pabp_shield
            return round(rate, 4), "PAN2_PAN3_SLOW_PROTECTED"
        else:
            # Smooth transition into rapid CCR4-NOT terminal catalysis below 25 nt
            terminal_accel = 1.0 + 0.45 * ((self.critical_polya_nt - polya_length_nt) / self.critical_polya_nt)
            rate = CCR4_NOT_TERMINAL_RATE_NT_HR * codon_protection * m1psi_factor * terminal_accel
            return round(rate, 4), "CCR4_NOT_RAPID_TERMINAL"

    def evaluate_tasep_ribosome_traffic(
        self,
        initiation_rate_hz: float,
        codon_optimality_index: float = 0.85,
        m1psi_modified: bool = True,
    ) -> TASEPRibosomeCollisionMetrics:
        """Evaluate 1D TASEP ribosome translocation, spacing, and No-Go Decay (NGD) collision risk.

        Models ribosomes with exclusion footprint L_ribo = 10 codons moving at velocity v_elong
        across the 1,480-codon CFTR transcript.
        """
        if initiation_rate_hz < 0.0:
            raise ValueError("Initiation rate cannot be negative")
        if not (0.0 <= codon_optimality_index <= 1.0):
            raise ValueError("Codon optimality index must be in [0.0, 1.0]")

        base_v = M1PSI_OPTIMIZED_ELONGATION_CODONS_S if m1psi_modified else UNMODIFIED_ELONGATION_CODONS_S
        # Elongation velocity scales with codon optimality (optimal tRNA pool matching)
        v_elong = base_v * (0.65 + 0.41 * codon_optimality_index)

        # Maximum steric TASEP current for extended particles of size l = 10 codons:
        # J_max = v_elong / (l * (1 + sqrt(l))^2) or continuum Lakatos-Chou limit;
        # practically, a ribosome clears its 10-codon footprint in t_clear = 10 / v_elong seconds.
        footprint_clearance_hz = v_elong / float(RIBOSOME_FOOTPRINT_CODONS)
        max_current_hz = footprint_clearance_hz * 0.25  # TASEP maximal current J_max = gamma / 4 in normalized units

        # Normalized queuing load ratio q = k_init / max_current_hz
        load_ratio = initiation_rate_hz / max(1e-6, max_current_hz)
        # Coverage density rho in [0, 1] (fraction of CDS nucleotides covered by 30-nt ribosomes)
        coverage_density = min(0.88, 0.52 * load_ratio / (1.0 + 0.45 * load_ratio))

        active_ribosomes = (coverage_density * self.cds_length_nt) / float(RIBOSOME_FOOTPRINT_NT)
        if active_ribosomes > 0.5:
            mean_spacing_nt = self.cds_length_nt / active_ribosomes
        else:
            mean_spacing_nt = float(self.cds_length_nt)

        # Di-ribosome collision probability rises sharply when load_ratio > 0.72 (spacing < 85 nt)
        collision_prob = 1.0 / (1.0 + math.exp(-8.5 * (load_ratio - 0.78)))
        if initiation_rate_hz == 0.0:
            collision_prob = 0.0

        # ZNF598/Hel2-mediated No-Go Decay (NGD) cleavage rate (hr^-1)
        ngd_cleavage_rate_hr = 0.38 * (collision_prob ** 1.6)

        # Effective steady-state protein completion flux (channels/sec before folding factor)
        collision_impedance = max(0.25, 1.0 - 0.65 * collision_prob)
        steady_flux_hz = min(initiation_rate_hz, max_current_hz) * collision_impedance

        if load_ratio < 0.68:
            regime = "FREE_FLOW_OPTIMAL"
        elif load_ratio <= 0.92:
            regime = "HIGH_DENSITY_COOPERATIVE"
        else:
            regime = "QUEUING_COLLISION_NGD_RISK"

        return TASEPRibosomeCollisionMetrics(
            initiation_rate_hz=round(initiation_rate_hz, 4),
            elongation_velocity_codons_s=round(v_elong, 3),
            max_lattice_current_hz=round(max_current_hz, 4),
            ribosome_coverage_density=round(coverage_density, 4),
            mean_inter_ribosome_spacing_nt=round(mean_spacing_nt, 1),
            active_polysome_ribosomes=round(active_ribosomes, 2),
            collision_probability=round(collision_prob, 4),
            ngd_cleavage_rate_hr=round(ngd_cleavage_rate_hr, 4),
            steady_translation_flux_hz=round(steady_flux_hz, 4),
            traffic_regime=regime,
        )

    def simulate_24h_trajectory(
        self,
        design_name: str = "VX522_Optimized_m1Psi_110A",
        initial_polya_nt: float = OPTIMAL_POLYA_TAIL_NT,
        target_initiation_hz: float = OPTIMAL_INITIATION_RATE_HZ,
        codon_optimality_index: float = 0.88,
        m1psi_modified: bool = True,
        segmented_linker: bool = True,
        duration_hr: float = 24.0,
        dt_hr: float = 0.5,
    ) -> PolyATrajectoryResult:
        """Integrate coupled poly(A) deadenylation, circularization, and CFTR translation over 24 hours."""
        if initial_polya_nt <= 0.0:
            raise ValueError("Initial poly(A) tail length must be positive")
        if duration_hr <= 0.0 or dt_hr <= 0.0:
            raise ValueError("Simulation duration and time step must be positive")

        tasep = self.evaluate_tasep_ribosome_traffic(
            initiation_rate_hz=target_initiation_hz,
            codon_optimality_index=codon_optimality_index,
            m1psi_modified=m1psi_modified,
        )

        polya_nt = float(initial_polya_nt)
        intact_fraction = 1.0
        cumulative_channels = 0.0
        half_life_hr: Optional[float] = None
        time_to_crit_hr: Optional[float] = None
        peak_polysome = 0.0

        steps = int(round(duration_hr / dt_hr))
        trajectory: List[DeadenylationStatePoint] = []

        for step in range(steps + 1):
            t_hr = round(step * dt_hr, 2)
            n_pabp = self.calculate_pabp_occupancy(polya_nt)
            p_circ = self.calculate_closed_loop_probability(polya_nt, m1psi_modified=m1psi_modified)
            deaden_rate, phase = self.calculate_biphasic_deadenylation_rate(
                polya_nt,
                codon_optimality_index=codon_optimality_index,
                m1psi_modified=m1psi_modified,
                segmented_linker=segmented_linker,
            )

            eff_init_hz = target_initiation_hz * p_circ
            step_tasep = self.evaluate_tasep_ribosome_traffic(
                initiation_rate_hz=eff_init_hz,
                codon_optimality_index=codon_optimality_index,
                m1psi_modified=m1psi_modified,
            )
            peak_polysome = max(peak_polysome, step_tasep.active_polysome_ribosomes)

            trajectory.append(
                DeadenylationStatePoint(
                    time_hr=t_hr,
                    polya_length_nt=round(polya_nt, 2),
                    bound_pabp_monomers=n_pabp,
                    closed_loop_probability=p_circ,
                    effective_initiation_hz=round(eff_init_hz, 4),
                    deadenylation_phase=phase,
                    intact_transcript_fraction=round(intact_fraction, 4),
                    cumulative_cftr_channels=round(cumulative_channels, 1),
                )
            )

            if step == steps:
                break

            # Functional CFTR channels synthesized during dt_hr (3600 s/hr)
            channels_step = (
                intact_fraction
                * step_tasep.steady_translation_flux_hz
                * CO_TRANSLATIONAL_FOLDING_EFFICIENCY
                * (dt_hr * 3600.0)
            )
            cumulative_channels += channels_step

            # Update poly(A) length
            next_polya = max(0.0, polya_nt - deaden_rate * dt_hr)
            if time_to_crit_hr is None and next_polya <= self.critical_polya_nt:
                time_to_crit_hr = t_hr + dt_hr

            # Transcript decay rate combines DCP2 decapping (accelerated when P_circ drops) + NGD + innate RNase L
            basal_innate_decay_hr = 0.022 if m1psi_modified else 0.085
            decapping_decay_hr = 0.018 * (1.0 - p_circ) + (0.32 if next_polya < self.critical_polya_nt else 0.0)
            total_decay_rate_hr = basal_innate_decay_hr + decapping_decay_hr + step_tasep.ngd_cleavage_rate_hr

            intact_fraction = intact_fraction * math.exp(-total_decay_rate_hr * dt_hr)
            if half_life_hr is None and intact_fraction <= 0.50:
                half_life_hr = t_hr + dt_hr

            polya_nt = next_polya

        if half_life_hr is None:
            half_life_hr = duration_hr
        if time_to_crit_hr is None:
            time_to_crit_hr = duration_hr

        return PolyATrajectoryResult(
            design_name=design_name,
            initial_polya_nt=initial_polya_nt,
            m1psi_modified=m1psi_modified,
            segmented_linker=segmented_linker,
            codon_optimality_index=codon_optimality_index,
            target_initiation_hz=target_initiation_hz,
            functional_half_life_hr=round(half_life_hr, 2),
            time_to_critical_threshold_hr=round(time_to_crit_hr, 2),
            peak_polysome_size=round(peak_polysome, 2),
            collision_probability=tasep.collision_probability,
            total_24h_cftr_channels=round(cumulative_channels, 1),
            curative_yield_achieved=cumulative_channels >= CURATIVE_CHANNEL_YIELD_THRESHOLD,
            trajectory=trajectory,
        )

    def compare_therapeutic_constructs(self) -> Dict[str, Any]:
        """Compare canonical synthetic CFTR mRNA constructs across poly(A), m1Psi, and initiation regimes."""
        optimized = self.simulate_24h_trajectory(
            design_name="VX522_m1Psi_110A_Segmented",
            initial_polya_nt=110.0,
            target_initiation_hz=0.085,
            codon_optimality_index=0.88,
            m1psi_modified=True,
            segmented_linker=True,
        )
        short_tail = self.simulate_24h_trajectory(
            design_name="Short_PolyA_45nt_Control",
            initial_polya_nt=45.0,
            target_initiation_hz=0.085,
            codon_optimality_index=0.88,
            m1psi_modified=True,
            segmented_linker=False,
        )
        unmodified_uridine = self.simulate_24h_trajectory(
            design_name="Unmodified_Uridine_100A",
            initial_polya_nt=100.0,
            target_initiation_hz=0.085,
            codon_optimality_index=0.70,
            m1psi_modified=False,
            segmented_linker=False,
        )
        collision_overdriven = self.simulate_24h_trajectory(
            design_name="Hyper_Initiation_Collision_Stalled",
            initial_polya_nt=110.0,
            target_initiation_hz=0.14,
            codon_optimality_index=0.65,
            m1psi_modified=False,
            segmented_linker=True,
        )

        yield_gain_vs_unmodified = round(
            optimized.total_24h_cftr_channels / max(1.0, unmodified_uridine.total_24h_cftr_channels),
            2,
        )
        return {
            "optimized_construct": asdict(optimized),
            "short_tail_control": asdict(short_tail),
            "unmodified_uridine_control": asdict(unmodified_uridine),
            "collision_overdriven_control": asdict(collision_overdriven),
            "yield_fold_gain_vs_unmodified": yield_gain_vs_unmodified,
        }


def main(argv: Optional[List[str]] = None) -> int:
    """CLI entry point for CFTR mRNA poly(A) deadenylation and TASEP collision modeling."""
    parser = argparse.ArgumentParser(
        description="Biophysical Modeling of Inhaled CFTR mRNA Poly(A) Tail Deadenylation Kinetics and Ribosomal Translocation Collision Protection (T-0449)."
    )
    parser.add_argument("--check", action="store_true", help="Run deterministic self-verification checks")
    parser.add_argument("--json", action="store_true", help="Print full JSON report")
    parser.add_argument("--polya", type=float, default=110.0, help="Initial poly(A) tail length in nt (default: 110.0)")
    parser.add_argument("--init-hz", type=float, default=0.085, help="Target ribosome initiation rate in Hz (default: 0.085)")
    args = parser.parse_args(argv)

    engine = MRNAPolyADeadenylationEngine()

    if args.check:
        comp = engine.compare_therapeutic_constructs()
        opt = comp["optimized_construct"]
        short = comp["short_tail_control"]
        unmod = comp["unmodified_uridine_control"]
        assert opt["total_24h_cftr_channels"] >= CURATIVE_CHANNEL_YIELD_THRESHOLD
        assert opt["functional_half_life_hr"] > short["functional_half_life_hr"]
        assert opt["total_24h_cftr_channels"] > unmod["total_24h_cftr_channels"] * 2.0
        print(
            f"Self-check PASS: 24h CFTR Yield={opt['total_24h_cftr_channels']} channels/mRNA, "
            f"t1/2={opt['functional_half_life_hr']}h, Gain vs Unmodified={comp['yield_fold_gain_vs_unmodified']}x"
        )
        return 0

    sim = engine.simulate_24h_trajectory(
        initial_polya_nt=args.polya,
        target_initiation_hz=args.init_hz,
    )

    if args.json:
        payload = {
            "simulation": asdict(sim),
            "grounding": {
                "trials": [NCT_VX522_MRNA_TRIAL],
                "structures": [PDB_CCR4_NOT_COMPLEX, PDB_ACTIVE_CFTR],
                "dois": [DOI_POLYA_PROFILING, DOI_CODON_OPTIMALITY, DOI_PSEUDOURIDINE_MRNA],
                "nlm_attribution": NLM_ATTRIBUTION,
                "research_boundary": RESEARCH_BOUNDARY,
            },
        }
        print(json.dumps(payload, indent=2))
        return 0

    print("==========================================================================")
    print("Sensei Ewok CF Lab - CFTR mRNA Poly(A) Deadenylation & TASEP Model (T-0449)")
    print(f"Construct: {sim.design_name} | Initial Poly(A): {sim.initial_polya_nt:.0f} nt | Init: {sim.target_initiation_hz:.3f} Hz")
    print(f"Grounded Registry: {NCT_VX522_MRNA_TRIAL} | Structures: PDB {PDB_CCR4_NOT_COMPLEX}, PDB {PDB_ACTIVE_CFTR}")
    print(f"Primary DOIs: {DOI_POLYA_PROFILING}, {DOI_CODON_OPTIMALITY}, {DOI_PSEUDOURIDINE_MRNA}")
    print("==========================================================================")
    print(f"Functional mRNA Half-Life:     {sim.functional_half_life_hr:.2f} hr")
    print(f"Time to 25-nt CCR4 Threshold:  {sim.time_to_critical_threshold_hr:.2f} hr")
    print(f"Peak Active Polysome Size:     {sim.peak_polysome_size:.1f} ribosomes / CDS")
    print(f"TASEP Collision Probability:   {sim.collision_probability:.4f}")
    print(f"24-Hour CFTR Channel Yield:    {sim.total_24h_cftr_channels:.1f} channels / mRNA")
    print(f"Curative Threshold (>=850):    {sim.curative_yield_achieved}")
    print("==========================================================================")
    print(f"Attribution: {NLM_ATTRIBUTION}")
    print(f"Boundary:    {RESEARCH_BOUNDARY}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
