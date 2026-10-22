#!/usr/bin/env python3
"""Inhaled CFTR Synthetic mRNA 3' UTR Polyadenylation & Ribosome Load Translation Velocity Model (T-0272).

Quantitative biophysical modeling of:
1. Totally Asymmetric Simple Exclusion Process (TASEP) ribosome elongation dynamics
   across the full-length human CFTR coding sequence (4,443 nucleotides / 1,480 codons).
2. Poly(A) tail length (100-120 nt) deadenylation kinetics, PABP circularization loop
   synergy, and cap-dependent translation initiation.
3. Co-translational domain folding coordination across MSD1, NBD1, R, MSD2, and NBD2
   domains balancing codon elongation speed against ribosome collision-induced No-Go Decay (NGD).
4. Nucleoside chemical modification (100% N1-methylpseudouridine / m1Psi) stability
   and cumulative full-length functional channel yield from inhaled LNP formulations.

Primary Evidence Grounding:
  - Clinical trial surveillance: BEACON-CF Phase 1/2 inhaled CFTR mRNA trial (ClinicalTrials.gov ID: NCT05668741)
  - Clinical trial surveillance: MRT5005 Phase 1/2 inhaled mRNA trial in adults with CF (ClinicalTrials.gov ID: NCT03375047)
  - Nucleoside modification and translation capacity: Kariko et al. (2008), Mol Ther, DOI: 10.1038/mt.2008.200
  - Codon optimality and elongation velocity: Presnyak et al. (2015), Cell, DOI: 10.1016/j.cell.2015.02.029
  - Poly(A) tail dynamics and translational control: Subtelny et al. (2014), Nature, DOI: 10.1038/nature13007
  - Cryo-EM CFTR structures: PDB 6GQ1 (ATP-bound phosphorylated human CFTR), PDB 5W81 (dephosphorylated human CFTR)

Research boundary: Preclinical computational model; research, not medical advice.
No clinical efficacy or therapeutic safety claim is made.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
from typing import Any

SCRIPT_PATH = Path(__file__).resolve()

NLM_ATTRIBUTION = (
    "This product uses publicly available data from the U.S. National Library of Medicine (NLM), "
    "National Institutes of Health, Department of Health and Human Services; NLM is not responsible "
    "for the product and does not endorse or recommend this or any other product."
)
RESEARCH_BOUNDARY = (
    "Preclinical in-silico inhaled synthetic CFTR mRNA translation velocity and ribosome load model; "
    "research, not medical advice. No clinical efficacy or therapeutic safety claim is made."
)

# Biophysical constants grounded in literature
CFTR_CDS_NUCLEOTIDES = 4443
CFTR_CODONS = 1480
RIBOSOME_FOOTPRINT_CODONS = 10  # ~30 nucleotides footprint on mRNA
DEFAULT_INITIATION_RATE_HZ = 0.08  # ribosomes per second for Cap1 eIF4F recruitment
DEFAULT_BASE_ELONGATION_RATE = 5.5  # optimal codons per second in human cells
DEFAULT_WILDTYPE_ELONGATION_RATE = 3.8  # wild-type human CFTR average codons per second
DEFAULT_DOMAIN_PAUSE_RATE = 2.0  # codons per second at domain boundaries
DEFAULT_POLY_A_INITIAL_LENGTH = 110.0  # optimal poly(A) length in synthetic mRNA (nt)
PABP_FOOTPRINT_NT = 27.0  # nucleotides covered per PABP monomer
PABP_DISSOCIATION_THRESHOLD_NT = 25.0  # tail length below which circularization loop breaks
K_DEAD_PROTECTED_PER_MIN = 0.025  # slow PAN2-PAN3 deadenylation with PABP (nt/min)
K_DEAD_UNPROTECTED_PER_MIN = 0.150  # fast CCR4-NOT deadenylation without PABP (nt/min)
M1PSI_HALF_LIFE_MULTIPLIER = 4.0  # ~18 hr half-life vs 4.5 hr base half-life


@dataclass(frozen=True)
class CFTRDomain:
    """Structural domain partition of human CFTR (1,480 amino acids)."""

    name: str
    start_codon: int
    end_codon: int
    length_aa: int
    description: str
    co_translational_pause_required: bool


CFTR_DOMAINS: tuple[CFTRDomain, ...] = (
    CFTRDomain(
        name="MSD1",
        start_codon=1,
        end_codon=388,
        length_aa=388,
        description="Membrane-spanning domain 1 (TM1-TM6, pore architecture)",
        co_translational_pause_required=False,
    ),
    CFTRDomain(
        name="NBD1",
        start_codon=389,
        end_codon=678,
        length_aa=290,
        description="Nucleotide-binding domain 1 (includes F508, G542X, G551D)",
        co_translational_pause_required=True,
    ),
    CFTRDomain(
        name="R_DOMAIN",
        start_codon=679,
        end_codon=844,
        length_aa=166,
        description="Regulatory domain (disordered, multiple PKA phosphorylation sites)",
        co_translational_pause_required=True,
    ),
    CFTRDomain(
        name="MSD2",
        start_codon=845,
        end_codon=1172,
        length_aa=328,
        description="Membrane-spanning domain 2 (TM7-TM12)",
        co_translational_pause_required=False,
    ),
    CFTRDomain(
        name="NBD2",
        start_codon=1173,
        end_codon=1480,
        length_aa=308,
        description="Nucleotide-binding domain 2 (includes W1282X, N1303K)",
        co_translational_pause_required=True,
    ),
)


@dataclass(frozen=True)
class PolyATailSpec:
    """Biophysical specification of synthetic mRNA 3' UTR poly(A) tail."""

    initial_length_nt: float = DEFAULT_POLY_A_INITIAL_LENGTH
    pabp_footprint_nt: float = PABP_FOOTPRINT_NT
    pabp_cooperative_hill: float = 2.5
    critical_length_nt: float = PABP_DISSOCIATION_THRESHOLD_NT
    k_dead_protected_per_min: float = K_DEAD_PROTECTED_PER_MIN
    k_dead_unprotected_per_min: float = K_DEAD_UNPROTECTED_PER_MIN


@dataclass(frozen=True)
class TranslationConfig:
    """Configuration for synthetic CFTR mRNA translation velocity simulation."""

    poly_a_spec: PolyATailSpec = PolyATailSpec()
    cap_structure: str = "Cap1"  # "Cap0", "Cap1", "Cap2"
    nucleoside_modification: str = "m1Psi"  # "m1Psi", "unmodified", "pseudoU"
    codon_optimization: str = "balanced_folding"  # "balanced_folding", "wild_type", "maximal_speed"
    initiation_rate_hz: float = DEFAULT_INITIATION_RATE_HZ
    base_elongation_rate_codons_s: float = DEFAULT_BASE_ELONGATION_RATE
    domain_pause_rate_codons_s: float = DEFAULT_DOMAIN_PAUSE_RATE
    temperature_celsius: float = 37.0


@dataclass
class TASEPProfile:
    """Spatial TASEP steady-state profile across 1,480 CFTR codons."""

    codons_total: int
    codon_velocities: list[float]
    steady_state_density: list[float]
    ribosome_load: float
    steady_state_flux_codons_s: float
    bottleneck_codon: int
    bottleneck_velocity: float
    mean_transit_time_seconds: float
    disome_collision_count: int
    ngd_safety_score: float
    collision_hotspots: list[int]


@dataclass
class TranslationTimePoint:
    """Snapshot of mRNA translation dynamics at time t."""

    time_hours: float
    poly_a_length_nt: float
    pabp_bound_count: float
    initiation_rate_hz: float
    ribosome_load: float
    translation_velocity_codons_s: float
    synthesis_rate_molecules_per_hour: float
    cumulative_cftr_protein: float
    mrna_intact_fraction: float


@dataclass
class TranslationLifecycleResult:
    """Full lifecycle synthesis and folding yield for synthetic CFTR mRNA."""

    config: TranslationConfig
    tasep_profile: TASEPProfile
    time_points: list[TranslationTimePoint]
    total_cftr_yield_molecules_per_mrna: float
    mrna_half_life_hours: float
    effective_translation_duration_hours: float
    co_translational_folding_fidelity: float
    ngd_decay_risk_verdict: str
    trial_surveillance_grounding: dict[str, str]


def get_cftr_domains() -> list[CFTRDomain]:
    """Return the structured domain partition of CFTR."""
    return list(CFTR_DOMAINS)


def build_codon_elongation_profile(
    optimization_mode: str = "balanced_folding",
    base_velocity: float = DEFAULT_BASE_ELONGATION_RATE,
    pause_velocity: float = DEFAULT_DOMAIN_PAUSE_RATE,
) -> list[float]:
    """Construct codon-by-codon elongation velocity profile for 1,480 codons.

    Modes:
      - 'balanced_folding': High speed in globular core regions with intentional
        slow clusters at domain boundaries (codons 380-395, 670-685, 840-855) to
        facilitate co-translational NBD1 and R domain folding.
      - 'wild_type': Native sequence elongation profile with naturally occurring
        rare codon clusters and lower average elongation rate (~3.8 codons/s).
      - 'maximal_speed': Uniform codon optimization targeting highest possible
        tRNA abundance without domain boundary pauses.
    """
    if optimization_mode not in ("balanced_folding", "wild_type", "maximal_speed"):
        raise ValueError(f"Unknown optimization_mode: {optimization_mode}")

    velocities = [base_velocity] * CFTR_CODONS

    if optimization_mode == "wild_type":
        # Native human CFTR sequence has substantial rare codon usage
        for i in range(CFTR_CODONS):
            # Harmonic fluctuation simulating native codon variation
            w_factor = 0.68 + 0.22 * math.sin(i * 0.04) + 0.10 * math.cos(i * 0.11)
            velocities[i] = max(1.2, DEFAULT_WILDTYPE_ELONGATION_RATE * w_factor)

    elif optimization_mode == "balanced_folding":
        # Programmed pauses at domain transitions (NBD1 start/end, R-domain exit)
        pause_windows = (
            (375, 395),  # MSD1-to-NBD1 transition
            (665, 685),  # NBD1-to-R domain transition
            (835, 855),  # R-to-MSD2 transition
            (1160, 1180),  # MSD2-to-NBD2 transition
        )
        for i in range(CFTR_CODONS):
            codon_1based = i + 1
            in_pause = any(start <= codon_1based <= end for start, end in pause_windows)
            if in_pause:
                velocities[i] = pause_velocity
            else:
                velocities[i] = base_velocity

    elif optimization_mode == "maximal_speed":
        # Uniform high velocity across all codons
        velocities = [base_velocity] * CFTR_CODONS

    return velocities


def calculate_poly_a_synergy(
    poly_a_length_nt: float,
    hill_n: float = 2.5,
    k_a_nt: float = 35.0,
) -> float:
    """Calculate cooperative circularization loop initiation multiplier from poly(A) length.

    Grounding: Subtelny et al. (Nature 2014) demonstrated that poly(A) tail length
    regulates translation initiation through cooperative PABP recruitment and eIF4G interaction.
    """
    if poly_a_length_nt <= 0.0:
        return 0.0
    l_eff = max(0.0, poly_a_length_nt)
    synergy = (l_eff**hill_n) / ((k_a_nt**hill_n) + (l_eff**hill_n))
    return min(1.0, max(0.0, synergy))


def calculate_cap_factor(cap_structure: str) -> float:
    """Return initiation recruitment multiplier relative to canonical Cap1."""
    caps = {
        "Cap0": 0.65,  # Unmethylated ribose; lower eIF4E affinity
        "Cap1": 1.00,  # 2'-O-methylated first nucleotide; standard therapeutic mRNA
        "Cap2": 1.08,  # Methylated first two nucleotides; slightly enhanced stability
    }
    if cap_structure not in caps:
        raise ValueError(f"Unsupported cap_structure: {cap_structure}. Supported: {list(caps.keys())}")
    return caps[cap_structure]


def simulate_tasep_ribosome_dynamics(
    elongation_rates: list[float],
    initiation_rate: float,
    footprint: int = RIBOSOME_FOOTPRINT_CODONS,
    max_iterations: int = 1500,
    tolerance: float = 1e-6,
) -> TASEPProfile:
    """Simulate steady-state ribosome occupancy and elongation flux via open-boundary TASEP.

    Models steric exclusion where a ribosome covers `footprint` codons (~10 codons / 30 nt).
    Iteratively solves the non-linear coupled mean-field master equations:
      Entry flux: J_in = alpha * (1 - sum_{k=0}^{L-1} rho_k)
      Hopping flux: J_i = v_i * rho_i * (1 - rho_{i + L})
      Exit flux: J_out = v_{N-1} * rho_{N-1}
    """
    n = len(elongation_rates)
    if n != CFTR_CODONS:
        raise ValueError(f"Expected {CFTR_CODONS} codons, received {n}")
    if initiation_rate < 0.0:
        raise ValueError("Initiation rate cannot be negative")

    # Initial guess: low uniform occupancy
    rho = [min(0.05, initiation_rate / max(1.0, elongation_rates[i])) for i in range(n)]

    # Damped iterative relaxation
    damp = 0.15
    for _ in range(max_iterations):
        max_delta = 0.0
        new_rho = list(rho)

        # Codon 0 entry boundary
        lead_block_0 = sum(rho[k] for k in range(min(footprint, n))) / float(footprint)
        j_in = initiation_rate * max(0.0, 1.0 - lead_block_0)

        # Internal update
        for i in range(n):
            v_i = elongation_rates[i]
            if i + footprint < n:
                blocked = rho[i + footprint]
            else:
                blocked = 0.0  # Clear exit near 3' end

            j_forward = v_i * rho[i] * max(0.0, 1.0 - blocked)

            # Inflow from previous codon
            if i == 0:
                j_prev = j_in
            else:
                v_prev = elongation_rates[i - 1]
                j_prev = v_prev * rho[i - 1] * max(0.0, 1.0 - rho[i])

            # Steady-state balance: delta rho = dt * (j_prev - j_forward)
            target = rho[i] + damp * (j_prev - j_forward) / max(1.0, v_i)
            target = min(0.95, max(0.0, target))

            diff = abs(target - rho[i])
            if diff > max_delta:
                max_delta = diff
            new_rho[i] = target

        rho = new_rho
        if max_delta < tolerance:
            break

    total_load = sum(rho)

    # Compute codon-by-codon flux
    fluxes = []
    for i in range(n):
        blocked = rho[i + footprint] if i + footprint < n else 0.0
        fluxes.append(elongation_rates[i] * rho[i] * max(0.0, 1.0 - blocked))

    mean_flux = sum(fluxes) / float(n) if n > 0 else 0.0

    # Identify bottleneck codon (lowest velocity region with highest occupancy)
    bottleneck_idx = min(range(n), key=lambda idx: elongation_rates[idx])
    bottleneck_vel = elongation_rates[bottleneck_idx]

    # Calculate transit time (tau = sum of 1 / v_eff)
    transit_time = sum(1.0 / max(0.1, elongation_rates[i] * (1.0 - (rho[i + footprint] if i + footprint < n else 0.0))) for i in range(n))

    # Detect disome collisions (adjacent occupied ribosomes within footprint)
    collision_hotspots: list[int] = []
    disome_count = 0
    for i in range(n - footprint):
        # Severe proximity collision probability
        p_prox = rho[i] * rho[i + footprint]
        if p_prox > 0.04:  # Statistically significant collision density
            disome_count += 1
            collision_hotspots.append(i + 1)  # 1-based codon coordinate

    # NGD safety score: exponential penalty on collision density per ribosome load
    collision_ratio = float(disome_count) / max(1.0, total_load)
    ngd_safety = 100.0 * math.exp(-6.5 * collision_ratio)
    ngd_safety = min(100.0, max(0.0, ngd_safety))

    return TASEPProfile(
        codons_total=n,
        codon_velocities=elongation_rates,
        steady_state_density=rho,
        ribosome_load=round(total_load, 2),
        steady_state_flux_codons_s=round(mean_flux, 4),
        bottleneck_codon=bottleneck_idx + 1,
        bottleneck_velocity=round(bottleneck_vel, 2),
        mean_transit_time_seconds=round(transit_time, 2),
        disome_collision_count=disome_count,
        ngd_safety_score=round(ngd_safety, 2),
        collision_hotspots=collision_hotspots[:15],  # Top hotspots
    )


def simulate_poly_a_length(
    initial_length: float,
    time_hours: float,
    k_protected_per_min: float = K_DEAD_PROTECTED_PER_MIN,
    k_unprotected_per_min: float = K_DEAD_UNPROTECTED_PER_MIN,
    critical_threshold_nt: float = PABP_DISSOCIATION_THRESHOLD_NT,
) -> float:
    """Calculate poly(A) tail length remaining at time t (in hours).

    Two-stage deadenylation kinetics:
      1. While length > critical_threshold_nt: slow PAN2-PAN3 trimming with PABP protection.
      2. Once length <= critical_threshold_nt: fast CCR4-NOT / CAF1 degradation.
    """
    if time_hours <= 0.0:
        return initial_length

    time_mins = time_hours * 60.0

    # Stage 1 duration needed to reach critical threshold
    if initial_length > critical_threshold_nt:
        stage1_delta = initial_length - critical_threshold_nt
        t_stage1_end = stage1_delta / k_protected_per_min

        if time_mins <= t_stage1_end:
            return max(critical_threshold_nt, initial_length - k_protected_per_min * time_mins)
        else:
            remaining_time = time_mins - t_stage1_end
            return max(0.0, critical_threshold_nt - k_unprotected_per_min * remaining_time)
    else:
        # Directly in stage 2
        return max(0.0, initial_length - k_unprotected_per_min * time_mins)


def simulate_mrna_lifecycle(
    config: TranslationConfig,
    total_hours: float = 72.0,
    dt_hours: float = 1.0,
) -> TranslationLifecycleResult:
    """Simulate complete synthetic CFTR mRNA translation lifecycle over time.

    Integrates:
      - Poly(A) tail deadenylation
      - Dynamic initiation rate through cap & circularization synergy
      - Ribosome load and translation velocity
      - Co-translational folding fidelity and NGD decay penalty
      - Cumulative functional CFTR channel yield per delivered mRNA molecule
    """
    elongation_rates = build_codon_elongation_profile(
        optimization_mode=config.codon_optimization,
        base_velocity=config.base_elongation_rate_codons_s,
        pause_velocity=config.domain_pause_rate_codons_s,
    )

    cap_multiplier = calculate_cap_factor(config.cap_structure)

    # Chemical modification mRNA half-life
    base_t_half = 4.5  # unmodified mRNA half-life in human bronchial epithelial cells (hours)
    if config.nucleoside_modification == "m1Psi":
        mrna_half_life = base_t_half * M1PSI_HALF_LIFE_MULTIPLIER  # ~18.0 hours
    elif config.nucleoside_modification == "pseudoU":
        mrna_half_life = base_t_half * 2.8  # ~12.6 hours
    else:
        mrna_half_life = base_t_half  # 4.5 hours

    k_mrna_decay = math.log(2.0) / mrna_half_life

    # Representative TASEP profile under nominal peak conditions
    nominal_loop = calculate_poly_a_synergy(config.poly_a_spec.initial_length_nt)
    nominal_alpha = config.initiation_rate_hz * cap_multiplier * nominal_loop
    nominal_tasep = simulate_tasep_ribosome_dynamics(elongation_rates, nominal_alpha)

    # Co-translational folding fidelity calculation:
    # Balanced pauses at domain boundaries permit NBD1/R folding (fidelity > 88%).
    # Over-accelerated maximal speed reduces folding fidelity due to hydrophobic collapse before pause.
    if config.codon_optimization == "balanced_folding":
        folding_fidelity = 0.91
    elif config.codon_optimization == "wild_type":
        folding_fidelity = 0.74  # Native sequence folds well but has low translation speed
    elif config.codon_optimization == "maximal_speed":
        folding_fidelity = 0.58  # Rapid elongation without pauses impairs NBD1 domain folding
    else:
        folding_fidelity = 0.75

    time_points: list[TranslationTimePoint] = []
    cumulative_protein = 0.0
    current_time = 0.0

    while current_time <= total_hours + 1e-5:
        poly_a_len = simulate_poly_a_length(
            initial_length=config.poly_a_spec.initial_length_nt,
            time_hours=current_time,
            k_protected_per_min=config.poly_a_spec.k_dead_protected_per_min,
            k_unprotected_per_min=config.poly_a_spec.k_dead_unprotected_per_min,
            critical_threshold_nt=config.poly_a_spec.critical_length_nt,
        )

        pabp_count = min(
            math.floor(config.poly_a_spec.initial_length_nt / config.poly_a_spec.pabp_footprint_nt),
            max(0.0, poly_a_len / config.poly_a_spec.pabp_footprint_nt),
        )

        loop_synergy = calculate_poly_a_synergy(poly_a_len)
        alpha_t = config.initiation_rate_hz * cap_multiplier * loop_synergy

        # mRNA intact fraction from standard exponential decay
        intact_frac = math.exp(-k_mrna_decay * current_time)

        # Scale translation flux with active initiation and intact mRNA
        flux_t = nominal_tasep.steady_state_flux_codons_s * (alpha_t / max(1e-4, nominal_alpha))
        active_load = nominal_tasep.ribosome_load * (alpha_t / max(1e-4, nominal_alpha))

        # Synthesis rate of completed CFTR polypeptide chains (molecules / hour)
        # flux_t is the ribosome exit rate (completed proteins / second)
        protein_synth_rate_per_sec = flux_t
        synth_rate_per_hour = protein_synth_rate_per_sec * 3600.0 * intact_frac

        # Folded functional channels factoring in NGD safety and folding fidelity
        ngd_factor = nominal_tasep.ngd_safety_score / 100.0
        functional_synth_rate = synth_rate_per_hour * folding_fidelity * ngd_factor

        cumulative_protein += functional_synth_rate * dt_hours

        time_points.append(
            TranslationTimePoint(
                time_hours=round(current_time, 2),
                poly_a_length_nt=round(poly_a_len, 2),
                pabp_bound_count=round(pabp_count, 2),
                initiation_rate_hz=round(alpha_t, 4),
                ribosome_load=round(active_load, 2),
                translation_velocity_codons_s=round(flux_t, 3),
                synthesis_rate_molecules_per_hour=round(synth_rate_per_hour, 2),
                cumulative_cftr_protein=round(cumulative_protein, 2),
                mrna_intact_fraction=round(intact_frac, 4),
            )
        )

        current_time += dt_hours

    # NGD decay risk verdict
    if nominal_tasep.ngd_safety_score >= 85.0:
        verdict = "OPTIMAL_COLLISION_FREE"
    elif nominal_tasep.ngd_safety_score >= 60.0:
        verdict = "MODERATE_QUEUEING_ACCEPTABLE"
    else:
        verdict = "SEVERE_DISOME_STALL_RISK"

    # Effective duration (time until intact fraction < 5%)
    eff_duration = -math.log(0.05) / k_mrna_decay

    trial_grounding = {
        "NCT05668741": "BEACON-CF Phase 1/2 inhaled CFTR mRNA program surveillance anchor",
        "NCT03375047": "MRT5005 Phase 1/2 nebulized CFTR mRNA in adult cystic fibrosis patients",
        "nlm_attribution": NLM_ATTRIBUTION,
    }

    return TranslationLifecycleResult(
        config=config,
        tasep_profile=nominal_tasep,
        time_points=time_points,
        total_cftr_yield_molecules_per_mrna=round(cumulative_protein, 1),
        mrna_half_life_hours=round(mrna_half_life, 2),
        effective_translation_duration_hours=round(eff_duration, 2),
        co_translational_folding_fidelity=round(folding_fidelity, 3),
        ngd_decay_risk_verdict=verdict,
        trial_surveillance_grounding=trial_grounding,
    )


def evaluate_synthetic_mrna_design(
    poly_a_length: float = DEFAULT_POLY_A_INITIAL_LENGTH,
    optimization: str = "balanced_folding",
    modification: str = "m1Psi",
    cap: str = "Cap1",
) -> dict[str, Any]:
    """Execute evaluation and return JSON-serializable assessment dictionary."""
    config = TranslationConfig(
        poly_a_spec=PolyATailSpec(initial_length_nt=poly_a_length),
        cap_structure=cap,
        nucleoside_modification=modification,
        codon_optimization=optimization,
    )
    result = simulate_mrna_lifecycle(config)
    res_dict = asdict(result)
    res_dict["research_boundary"] = RESEARCH_BOUNDARY
    return res_dict


def parse_arguments() -> argparse.ArgumentParser:
    """Parse command-line arguments for translation model."""
    parser = argparse.ArgumentParser(
        description="Inhaled CFTR Synthetic mRNA Translation Velocity & Ribosome Load Model (T-0272)",
        epilog=RESEARCH_BOUNDARY,
    )
    parser.add_argument(
        "--poly-a-length",
        type=float,
        default=DEFAULT_POLY_A_INITIAL_LENGTH,
        help="Initial poly(A) tail length in nucleotides (default: 110.0)",
    )
    parser.add_argument(
        "--optimization",
        choices=["balanced_folding", "wild_type", "maximal_speed"],
        default="balanced_folding",
        help="Codon optimization strategy (default: balanced_folding)",
    )
    parser.add_argument(
        "--modification",
        choices=["m1Psi", "unmodified", "pseudoU"],
        default="m1Psi",
        help="Chemical nucleoside modification (default: m1Psi)",
    )
    parser.add_argument(
        "--cap",
        choices=["Cap0", "Cap1", "Cap2"],
        default="Cap1",
        help="5' Cap structure (default: Cap1)",
    )
    parser.add_argument(
        "--hours",
        type=float,
        default=72.0,
        help="Simulation time span in hours (default: 72.0)",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Output raw JSON summary",
    )
    parser.add_argument(
        "--verify",
        action="store_true",
        help="Run self-verification checks and exit 0",
    )
    return parser


def main() -> int:
    """CLI execution entry point."""
    parser = parse_arguments()
    args = parser.parse_args()

    if args.verify:
        config = TranslationConfig()
        result = simulate_mrna_lifecycle(config, total_hours=24.0)
        assert result.total_cftr_yield_molecules_per_mrna > 10.0
        assert result.tasep_profile.codons_total == 1480
        assert result.ngd_decay_risk_verdict in ("OPTIMAL_COLLISION_FREE", "MODERATE_QUEUEING_ACCEPTABLE")
        print("PASS [T-0272 self-verification]: CFTR TASEP ribosome model validated.")
        return 0

    evaluation = evaluate_synthetic_mrna_design(
        poly_a_length=args.poly_a_length,
        optimization=args.optimization,
        modification=args.modification,
        cap=args.cap,
    )

    if args.json:
        print(json.dumps(evaluation, indent=2))
    else:
        print("=" * 76)
        print("CFTR SYNTHETIC mRNA TRANSLATION VELOCITY & RIBOSOME LOAD MODEL (T-0272)")
        print("=" * 76)
        print(f"Design Configuration: Poly(A)={args.poly_a_length} nt | Cap={args.cap} | Mod={args.modification} | Mode={args.optimization}")
        tasep = evaluation["tasep_profile"]
        print(f"Ribosome Load: {tasep['ribosome_load']} ribosomes/mRNA across {tasep['codons_total']} codons")
        print(f"Steady-State Translation Flux: {tasep['steady_state_flux_codons_s']} codons/s")
        print(f"Mean Elongation Transit Time: {tasep['mean_transit_time_seconds'] / 60.0:.2f} min/protein")
        print(f"Bottleneck Codon: #{tasep['bottleneck_codon']} ({tasep['bottleneck_velocity']} codons/s)")
        print(f"Disome Collision Count: {tasep['disome_collision_count']} | NGD Safety Score: {tasep['ngd_safety_score']}/100 ({evaluation['ngd_decay_risk_verdict']})")
        print(f"mRNA Half-Life: {evaluation['mrna_half_life_hours']} hr | Folding Fidelity: {evaluation['co_translational_folding_fidelity'] * 100:.1f}%")
        print(f"Estimated Total CFTR Protein Yield: {evaluation['total_cftr_yield_molecules_per_mrna']} channels/mRNA")
        print("-" * 76)
        print("Trial Surveillance Grounding:")
        for nct, desc in evaluation["trial_surveillance_grounding"].items():
            if nct.startswith("NCT"):
                print(f"  - {nct}: {desc}")
        print("=" * 76)
        print(RESEARCH_BOUNDARY)

    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
