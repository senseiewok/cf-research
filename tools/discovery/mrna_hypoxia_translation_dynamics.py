#!/usr/bin/env python3
"""Stochastic Single-Cell Model of Inhaled CFTR mRNA Translation Dynamics & Ribosome Collisions across Hypoxic Airway Epithelia (T-0452).

Quantitative biophysical modeling of:
1. Hypoxic Integrated Stress Response (ISR):
   - Formulates eIF2alpha Ser51 phosphorylation and ternary complex attenuation
     as a function of local oxygen tension (pO2: 0.5%–21.0% O2) in mucus-occluded airways.
   - Models initiation rate attenuation via cooperative Hill kinetics.
2. Ribosomal Collision Dynamics & No-Go Decay (NGD):
   - Computes elongation velocity across the 4,443 nt (1,480 codons) CFTR coding sequence.
   - Evaluates collision avoidance, disome/trisome accumulation, and ZNF598-mediated
     No-Go Decay (NGD) clearance comparing wild-type vs synonymous codon-optimized
     and 100% N1-methylpseudouridine (m1Psi) modified synthetic mRNA.
3. Translational Yield & Channel Output:
   - Predicts mature Band C CFTR channel synthesis rate per cell as a function of
     inhaled LNP cytosolic copy number (1–50 mRNAs/cell), co-translational folding efficiency,
     and peripheral membrane clearance kinetics.
4. Airway Heterogeneity & Curative Phenotype Rescue Thresholds:
   - Simulates single-cell distributions across conducting, intermediate, and severely
     occluded hypoxic ventilation zones.
   - Quantifies the fraction of epithelial cells achieving the clinical curative
     threshold (>= 10% wild-type CFTR channel density, >= 150 functional channels/cell).

Primary Literature Grounding:
  - Englehardt JF et al. (2023) Nature, DOI: 10.1038/s41586-023-06133-1
    (Engineered suppressor tRNAs restore CFTR function in patient-derived airway epithelia)
  - Manfredi C et al. (2025) Biochemistry, DOI: 10.1021/acs.biochem.5c00590
    (Barriers to the Pharmacologic Rescue of W1282X CFTR)
  - Middleton PG et al. (2019) N Engl J Med, DOI: 10.1056/NEJMoa1908639
    (Elexacaftor-Tezacaftor-Ivacaftor in Patients with Cystic Fibrosis with Phe508del-Minimal Function Genotypes)
  - Button B et al. (2021) Physiol Rev, DOI: 10.1152/physrev.00004.2021
    (Physiology and pathophysiology of human airway mucus)

Clinical Registries:
  - NCT05668741 (Phase 1/2 Study of Inhaled mRNA LNP VX-522 in Participants With Cystic Fibrosis)
  - NCT05248230 (Phase 1/2 Study of 4D-710 Aerosolized AAV Gene Therapy in Cystic Fibrosis)

Structural Coordinates:
  - RCSB PDB 6MSM (Dephosphorylated human CFTR inward-facing conformation)
  - RCSB PDB 5UAK (Phosphorylated active human CFTR dimer conformation)

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
import random
import sys
from typing import Any, Dict, List, Optional, Tuple

SCRIPT_PATH = Path(__file__).resolve()

NLM_ATTRIBUTION = (
    "This product uses publicly available data from the U.S. National Library of Medicine (NLM), "
    "National Institutes of Health, Department of Health and Human Services; NLM is not responsible "
    "for the product and does not endorse or recommend this or any other product."
)

RESEARCH_BOUNDARY = (
    "Preclinical in-silico inhaled synthetic CFTR mRNA translation dynamics and ribosome collision model "
    "under chronic airway hypoxia; research, not medical advice. No clinical efficacy or therapeutic safety claim is made."
)

# Registry & Literature Identifiers
DOI_SUPPRESSOR_TRNA = "10.1038/s41586-023-06133-1"
DOI_W1282X_RESCUE_BARRIER = "10.1021/acs.biochem.5c00590"
DOI_ETI_PHASE3_TRIAL = "10.1056/NEJMoa1908639"
DOI_AIRWAY_MUCUS_PHYSIOL = "10.1152/physrev.00004.2021"

NCT_VX522_MRNA_TRIAL = "NCT05668741"
NCT_4D710_AAV_TRIAL = "NCT05248230"

PDB_DEPHOSPHO_CFTR = "6MSM"
PDB_ACTIVE_CFTR = "5UAK"

# Biophysical Constants & Grounded Parameters
CFTR_CDS_NUCLEOTIDES: int = 4443
CFTR_CODONS: int = 1480
RIBOSOME_FOOTPRINT_CODONS: int = 10  # ~30 nucleotides footprint on mRNA

# Normoxic baseline translation parameters
DEFAULT_K_INIT_NORMOXIC_HZ: float = 0.080  # ribosomes/s initiated per mRNA under normoxia (21% O2)
MIN_K_INIT_FRACTION_HYPOXIA: float = 0.180  # residual initiation fraction under severe anoxia/hypoxia
HYPOXIA_KM_O2_PERCENT: float = 2.50  # pO2 percentage for half-maximal ISR initiation attenuation
HYPOXIA_HILL_COEFFICIENT: float = 2.00  # Hill cooperativity for eIF2alpha phosphorylation cascade

# Elongation velocity constants (codons / second)
V_ELONG_WILDTYPE_MEAN: float = 3.80  # Wild-type human CFTR mRNA average elongation rate
V_ELONG_WILDTYPE_PAUSE: float = 1.60  # Elongation rate at domain boundaries / rare codon clusters
V_ELONG_OPTIMIZED_MEAN: float = 5.60  # Codon-optimized synthetic mRNA average elongation rate
V_ELONG_OPTIMIZED_PAUSE: float = 4.20  # Codon-optimized rate at structural pause positions

# Chemical modification stability parameters
MRNA_HALF_LIFE_UNMODIFIED_HOURS: float = 4.50  # Inhaled unmodified mRNA half-life in airway epithelia
M1PSI_HALF_LIFE_MULTIPLIER: float = 4.00  # 100% N1-methylpseudouridine half-life extension (~18.0 h)
K_NGD_BASE_PER_HOUR: float = 0.350  # Base No-Go Decay mRNA cleavage rate per disome

# Folding, surface trafficking and clearance constants
CFTR_COTRANSLATIONAL_FOLDING_EFFICIENCY: float = 0.250  # ~25% Band C ER-to-surface maturation
CFTR_SURFACE_HALF_LIFE_HOURS: float = 20.00  # Mature CFTR surface half-life (~72,000 s)
WT_CFTR_CHANNELS_PER_CELL: float = 1500.0  # Healthy wild-type human ciliated epithelial channel density
CURATIVE_RESCUE_FRACTION_WT: float = 0.100  # >= 10% WT density (>= 150 channels/cell) clinical threshold
CURATIVE_RESCUE_CHANNELS_THRESHOLD: float = WT_CFTR_CHANNELS_PER_CELL * CURATIVE_RESCUE_FRACTION_WT


@dataclass(frozen=True)
class HypoxicISRState:
    """Hypoxic Integrated Stress Response and Initiation Attenuation State."""

    po2_percent: float
    effective_initiation_rate_hz: float
    normoxic_initiation_rate_hz: float
    initiation_attenuation_ratio: float
    eif2alpha_phosphorylation_index: float
    stress_regime: str
    description: str


@dataclass(frozen=True)
class CollisionNGDMetrics:
    """Ribosome Collision Dynamics, Disome Density, and No-Go Decay Risk."""

    is_codon_optimized: bool
    is_m1psi_modified: bool
    mean_elongation_velocity: float
    pause_elongation_velocity: float
    average_ribosome_occupancy: float
    collision_probability_per_transcript: float
    disome_density_per_mrna: float
    ngd_decay_rate_per_hour: float
    effective_mrna_half_life_hours: float
    ngd_escape_survival_fraction: float


@dataclass(frozen=True)
class SingleCellYieldMetrics:
    """Mature CFTR Channel Synthesis, Surface Density, and Curative Status."""

    mrna_copy_number: int
    po2_percent: float
    translation_run_time_seconds: float
    polypeptide_synthesis_rate_per_sec: float
    mature_cftr_synthesis_rate_per_sec: float
    steady_state_surface_channels: float
    percent_wild_type_cftr: float
    is_curative_threshold_met: bool


@dataclass(frozen=True)
class AirwayVentilationZone:
    """Microenvironmental airway zone with local oxygenation and epithelial density."""

    zone_name: str
    po2_percent: float
    population_fraction: float
    mucus_plug_status: str
    mean_lnp_mrna_uptake: float


class MRNAHypoxiaTranslationEngine:
    """Deterministic & Stochastic Biophysical Engine for CFTR mRNA Hypoxia Translation Dynamics."""

    def __init__(
        self,
        normoxic_initiation_hz: float = DEFAULT_K_INIT_NORMOXIC_HZ,
        km_o2: float = HYPOXIA_KM_O2_PERCENT,
        hill_coeff: float = HYPOXIA_HILL_COEFFICIENT,
    ) -> None:
        self.normoxic_initiation_hz = max(0.001, float(normoxic_initiation_hz))
        self.km_o2 = max(0.1, float(km_o2))
        self.hill_coeff = max(1.0, float(hill_coeff))

    def compute_hypoxic_isr_attenuation(self, po2_percent: float) -> HypoxicISRState:
        """Domain 1: Compute eIF2alpha phosphorylation and initiation rate attenuation under pO2."""
        po2 = max(0.0, float(po2_percent))

        # Cooperative Hill activation of ISR initiation attenuation
        # Fraction f_init = alpha_min + (1 - alpha_min) * (pO2^n / (pO2^n + Km^n))
        if po2 <= 0.0:
            frac = MIN_K_INIT_FRACTION_HYPOXIA
        else:
            po2_n = math.pow(po2, self.hill_coeff)
            km_n = math.pow(self.km_o2, self.hill_coeff)
            activation = po2_n / (po2_n + km_n)
            frac = MIN_K_INIT_FRACTION_HYPOXIA + (1.0 - MIN_K_INIT_FRACTION_HYPOXIA) * activation

        frac = max(MIN_K_INIT_FRACTION_HYPOXIA, min(1.0, frac))
        k_init = self.normoxic_initiation_hz * frac
        eif2_phospho = 1.0 - frac  # normalized phosphorylation index [0 = unphosphorylated, 1 = max]

        if po2 >= 15.0:
            regime = "NORMOXIA_CONDUCTING"
            desc = "Ventilated conducting airway with baseline unphosphorylated eIF2alpha and high initiation."
        elif po2 >= 4.0:
            regime = "INTERMEDIATE_HYPOXIA"
            desc = "Moderate hypoxia in mucociliary clearance-compromised zone; partial eIF2alpha phosphorylation."
        else:
            regime = "SEVERE_HYPOXIA_OCCLUSION"
            desc = "Severe mucus-occluded atelectatic plug zone; robust ISR activation and initiation shutdown."

        return HypoxicISRState(
            po2_percent=round(po2, 2),
            effective_initiation_rate_hz=round(k_init, 5),
            normoxic_initiation_rate_hz=round(self.normoxic_initiation_hz, 5),
            initiation_attenuation_ratio=round(frac, 4),
            eif2alpha_phosphorylation_index=round(eif2_phospho, 4),
            stress_regime=regime,
            description=desc,
        )

    def evaluate_collision_and_ngd(
        self,
        k_init_hz: float,
        codon_optimized: bool = True,
        m1psi_modified: bool = True,
    ) -> CollisionNGDMetrics:
        """Domain 2: Ribosome elongation, collision frequency, disome density, and No-Go Decay (NGD)."""
        k_init = max(0.0, float(k_init_hz))

        if codon_optimized:
            v_mean = V_ELONG_OPTIMIZED_MEAN
            v_pause = V_ELONG_OPTIMIZED_PAUSE
            beta_coll = 0.045  # reduced queuing susceptibility
        else:
            v_mean = V_ELONG_WILDTYPE_MEAN
            v_pause = V_ELONG_WILDTYPE_PAUSE
            beta_coll = 0.120  # higher susceptibility to pause-induced disomes

        # Average ribosome occupancy on the 1,480-codon transcript
        transit_time_sec = CFTR_CODONS / v_mean
        ribo_occupancy = k_init * transit_time_sec

        # Collision probability along transcript based on TASEP initiation vs pause bottleneck
        if k_init > 0.0 and v_pause > 0.0:
            flux_ratio = math.pow(k_init, 1.8) / (v_mean * v_pause)
            p_coll = 1.0 - math.exp(-beta_coll * flux_ratio * (CFTR_CODONS / 100.0))
            p_coll = max(0.0, min(0.999, p_coll))
        else:
            p_coll = 0.0

        # Disome density per mRNA transcript
        disomes = ribo_occupancy * p_coll

        # Base mRNA half-life and NGD accelerated turnover
        base_half_life = (
            MRNA_HALF_LIFE_UNMODIFIED_HOURS * M1PSI_HALF_LIFE_MULTIPLIER
            if m1psi_modified
            else MRNA_HALF_LIFE_UNMODIFIED_HOURS
        )
        base_decay_rate = math.log(2.0) / base_half_life

        # NGD decay enhancement from ZNF598 disome recognition
        k_ngd = K_NGD_BASE_PER_HOUR * (disomes / (disomes + 1.2))
        effective_decay_rate = base_decay_rate + k_ngd
        effective_half_life = math.log(2.0) / effective_decay_rate

        # NGD escape survival probability during transit
        transit_hours = transit_time_sec / 3600.0
        ngd_escape = math.exp(-k_ngd * transit_hours)

        return CollisionNGDMetrics(
            is_codon_optimized=codon_optimized,
            is_m1psi_modified=m1psi_modified,
            mean_elongation_velocity=round(v_mean, 2),
            pause_elongation_velocity=round(v_pause, 2),
            average_ribosome_occupancy=round(ribo_occupancy, 2),
            collision_probability_per_transcript=round(p_coll, 4),
            disome_density_per_mrna=round(disomes, 4),
            ngd_decay_rate_per_hour=round(k_ngd, 4),
            effective_mrna_half_life_hours=round(effective_half_life, 2),
            ngd_escape_survival_fraction=round(ngd_escape, 4),
        )

    def calculate_single_cell_yield(
        self,
        mrna_copy_number: int,
        po2_percent: float,
        codon_optimized: bool = True,
        m1psi_modified: bool = True,
        cotranslational_folding_efficiency: float = CFTR_COTRANSLATIONAL_FOLDING_EFFICIENCY,
        surface_half_life_hours: float = CFTR_SURFACE_HALF_LIFE_HOURS,
    ) -> SingleCellYieldMetrics:
        """Domain 3: Calculate mature CFTR channel synthesis rate and steady-state surface density."""
        copies = max(0, int(mrna_copy_number))
        isr = self.compute_hypoxic_isr_attenuation(po2_percent)
        coll = self.evaluate_collision_and_ngd(
            k_init_hz=isr.effective_initiation_rate_hz,
            codon_optimized=codon_optimized,
            m1psi_modified=m1psi_modified,
        )

        transit_time_sec = CFTR_CODONS / coll.mean_elongation_velocity

        if copies == 0 or isr.effective_initiation_rate_hz <= 0.0:
            polypeptide_rate = 0.0
            mature_cftr_rate = 0.0
            surface_channels = 0.0
        else:
            # Translation completion rate per cell accounting for mRNA stability and NGD clearance
            # Daily persistence factor over 24-hour therapeutic dosing interval
            persistence_factor = min(1.0, coll.effective_mrna_half_life_hours / 24.0)
            polypeptide_rate = (
                copies
                * persistence_factor
                * isr.effective_initiation_rate_hz
                * coll.ngd_escape_survival_fraction
            )
            mature_cftr_rate = polypeptide_rate * cotranslational_folding_efficiency
            # Steady-state surface channel density: N_ss = R_synth / k_deg
            surface_deg_rate_hz = math.log(2.0) / (surface_half_life_hours * 3600.0)
            surface_channels = mature_cftr_rate / surface_deg_rate_hz

        pct_wt = (surface_channels / WT_CFTR_CHANNELS_PER_CELL) * 100.0
        is_curative = surface_channels >= CURATIVE_RESCUE_CHANNELS_THRESHOLD

        return SingleCellYieldMetrics(
            mrna_copy_number=copies,
            po2_percent=round(po2_percent, 2),
            translation_run_time_seconds=round(transit_time_sec, 1),
            polypeptide_synthesis_rate_per_sec=round(polypeptide_rate, 4),
            mature_cftr_synthesis_rate_per_sec=round(mature_cftr_rate, 4),
            steady_state_surface_channels=round(surface_channels, 1),
            percent_wild_type_cftr=round(pct_wt, 2),
            is_curative_threshold_met=is_curative,
        )

    def simulate_airway_population_rescue(
        self,
        mean_lnp_dose: float = 20.0,
        cohort_size: int = 1000,
        codon_optimized: bool = True,
        m1psi_modified: bool = True,
        seed: int = 42,
    ) -> Dict[str, Any]:
        """Domain 4: Simulate epithelial population across heterogeneous airway ventilation zones."""
        rng = random.Random(seed)
        cohort_n = max(100, int(cohort_size))

        # Standard CF airway microenvironmental partition (Button et al., Physiol Rev 2021)
        zones = [
            AirwayVentilationZone("Conducting_Ventilated", 18.0, 0.40, "Mild_Mucus_Lining", mean_lnp_dose * 1.25),
            AirwayVentilationZone("Intermediate_Hypoxic", 5.0, 0.35, "Partial_Mucus_Plugging", mean_lnp_dose * 0.90),
            AirwayVentilationZone("Occluded_Atelectatic", 1.0, 0.25, "Dense_Mucus_Plug", mean_lnp_dose * 0.45),
        ]

        per_zone_results: Dict[str, Any] = {}
        total_rescued_count = 0
        all_channel_counts: List[float] = []

        for z in zones:
            zone_cells = int(cohort_n * z.population_fraction)
            rescued_in_zone = 0
            zone_channel_sum = 0.0

            for _ in range(zone_cells):
                # Poisson-distributed cellular LNP mRNA delivery
                lam = max(0.0, z.mean_lnp_mrna_uptake)
                # Knuth's Poisson generator
                l_val = math.exp(-lam) if lam < 30.0 else 0.0
                if l_val > 0.0:
                    k = 0
                    p = 1.0
                    while p > l_val:
                        k += 1
                        p *= rng.random()
                    mrna_copies = k - 1
                else:
                    mrna_copies = max(0, int(rng.gauss(lam, math.sqrt(lam))))

                cell_res = self.calculate_single_cell_yield(
                    mrna_copy_number=mrna_copies,
                    po2_percent=z.po2_percent,
                    codon_optimized=codon_optimized,
                    m1psi_modified=m1psi_modified,
                )
                zone_channel_sum += cell_res.steady_state_surface_channels
                all_channel_counts.append(cell_res.steady_state_surface_channels)
                if cell_res.is_curative_threshold_met:
                    rescued_in_zone += 1
                    total_rescued_count += 1

            mean_ch = zone_channel_sum / max(1, zone_cells)
            per_zone_results[z.zone_name] = {
                "po2_percent": z.po2_percent,
                "mucus_status": z.mucus_plug_status,
                "cells_in_zone": zone_cells,
                "mean_surface_channels_per_cell": round(mean_ch, 1),
                "percent_wt_cftr": round((mean_ch / WT_CFTR_CHANNELS_PER_CELL) * 100.0, 2),
                "cells_rescued_count": rescued_in_zone,
                "curative_rescue_percentage": round((rescued_in_zone / max(1, zone_cells)) * 100.0, 2),
            }

        overall_rescue_pct = (total_rescued_count / cohort_n) * 100.0
        overall_mean_channels = sum(all_channel_counts) / len(all_channel_counts)

        return {
            "cohort_size": cohort_n,
            "mean_lnp_delivery_dose": mean_lnp_dose,
            "is_codon_optimized": codon_optimized,
            "is_m1psi_modified": m1psi_modified,
            "curative_threshold_channels": CURATIVE_RESCUE_CHANNELS_THRESHOLD,
            "total_cells_rescued": total_rescued_count,
            "overall_curative_rescue_percentage": round(overall_rescue_pct, 2),
            "overall_mean_surface_channels": round(overall_mean_channels, 1),
            "overall_percent_wt_cftr": round((overall_mean_channels / WT_CFTR_CHANNELS_PER_CELL) * 100.0, 2),
            "per_zone_analysis": per_zone_results,
        }

    def run_comprehensive_translation_suite(self) -> Dict[str, Any]:
        """Execute full 4-domain biophysical evaluation of CFTR mRNA translation dynamics."""
        # Domain 1: pO2 gradient sweep
        po2_sweep = [21.0, 15.0, 10.0, 5.0, 2.5, 1.0, 0.5]
        isr_profile = [asdict(self.compute_hypoxic_isr_attenuation(p)) for p in po2_sweep]

        # Domain 2: Wild-type vs Codon-optimized & m1Psi collision comparison under normoxia and hypoxia
        coll_wt_norm = asdict(self.evaluate_collision_and_ngd(0.080, codon_optimized=False, m1psi_modified=False))
        coll_opt_norm = asdict(self.evaluate_collision_and_ngd(0.080, codon_optimized=True, m1psi_modified=True))
        coll_wt_hypox = asdict(self.evaluate_collision_and_ngd(0.025, codon_optimized=False, m1psi_modified=False))
        coll_opt_hypox = asdict(self.evaluate_collision_and_ngd(0.025, codon_optimized=True, m1psi_modified=True))

        # Domain 3: Dose-response curve across copy numbers (1 to 50 mRNAs/cell) under normoxia vs hypoxia
        doses = [0, 1, 5, 10, 20, 30, 50]
        dose_normoxia = [asdict(self.calculate_single_cell_yield(d, 21.0, True, True)) for d in doses]
        dose_hypoxia = [asdict(self.calculate_single_cell_yield(d, 2.0, True, True)) for d in doses]

        # Domain 4: Epithelial population simulation comparing Unmodified vs Optimized formulation
        pop_unmodified = self.simulate_airway_population_rescue(mean_lnp_dose=20.0, codon_optimized=False, m1psi_modified=False)
        pop_optimized = self.simulate_airway_population_rescue(mean_lnp_dose=20.0, codon_optimized=True, m1psi_modified=True)

        return {
            "metadata": {
                "task_id": "T-0452",
                "title": "Stochastic Single-Cell Model of Inhaled CFTR mRNA Translation Dynamics & Ribosome Collisions across Hypoxic Airway Epithelia",
                "answer_tier": "T2",
                "literature_grounding": [
                    DOI_SUPPRESSOR_TRNA,
                    DOI_W1282X_RESCUE_BARRIER,
                    DOI_ETI_PHASE3_TRIAL,
                    DOI_AIRWAY_MUCUS_PHYSIOL,
                ],
                "registry_grounding": [
                    NCT_VX522_MRNA_TRIAL,
                    NCT_4D710_AAV_TRIAL,
                ],
                "structural_grounding": [
                    PDB_DEPHOSPHO_CFTR,
                    PDB_ACTIVE_CFTR,
                ],
                "nlm_attribution": NLM_ATTRIBUTION,
                "research_boundary": RESEARCH_BOUNDARY,
            },
            "domain_1_hypoxic_isr_attenuation": {
                "po2_sweep": isr_profile,
                "attenuation_ratio_at_1pct_po2": round(isr_profile[-2]["initiation_attenuation_ratio"], 4),
            },
            "domain_2_ribosome_collision_and_ngd": {
                "normoxia_comparison": {
                    "unmodified_wildtype": coll_wt_norm,
                    "codon_optimized_m1psi": coll_opt_norm,
                    "disome_reduction_pct": round(
                        (1.0 - coll_opt_norm["disome_density_per_mrna"] / max(0.001, coll_wt_norm["disome_density_per_mrna"])) * 100.0,
                        2,
                    ),
                },
                "hypoxia_comparison": {
                    "unmodified_wildtype": coll_wt_hypox,
                    "codon_optimized_m1psi": coll_opt_hypox,
                },
            },
            "domain_3_translational_yield_and_channel_density": {
                "normoxia_dose_response": dose_normoxia,
                "hypoxia_dose_response": dose_hypoxia,
            },
            "domain_4_airway_population_curative_rescue": {
                "unmodified_formulation": pop_unmodified,
                "codon_optimized_m1psi_formulation": pop_optimized,
                "rescue_differential_pct": round(
                    pop_optimized["overall_curative_rescue_percentage"] - pop_unmodified["overall_curative_rescue_percentage"],
                    2,
                ),
            },
        }


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Stochastic Single-Cell Model of Inhaled CFTR mRNA Translation Dynamics & Ribosome Collisions across Hypoxic Airway Epithelia (T-0452)."
    )
    parser.add_argument(
        "--po2",
        type=float,
        default=None,
        help="Local oxygen tension percentage (pO2: 0.5% - 21.0% O2).",
    )
    parser.add_argument(
        "--mrna-copies",
        type=int,
        default=None,
        help="Cytosolic LNP-delivered CFTR mRNA copy number per cell (1 - 50).",
    )
    parser.add_argument(
        "--unmodified",
        action="store_true",
        help="Simulate unoptimized wild-type mRNA lacking m1Psi modification.",
    )
    parser.add_argument(
        "--simulate-population",
        action="store_true",
        help="Simulate 1,000-cell airway population across heterogeneous ventilation zones.",
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
    engine = MRNAHypoxiaTranslationEngine()

    if args.run_all or (args.po2 is None and args.mrna_copies is None and not args.simulate_population):
        results = engine.run_comprehensive_translation_suite()
        if args.json:
            print(json.dumps(results, indent=2))
        else:
            print("=" * 76)
            print("CFTR MRNA HYPOXIA TRANSLATION DYNAMICS & RIBOSOME COLLISION ENGINE (T-0452)")
            print("=" * 76)
            print(f"Primary Registries: {NCT_VX522_MRNA_TRIAL} (Inhaled mRNA LNP VX-522), {NCT_4D710_AAV_TRIAL}")
            print(f"Primary Literature: {DOI_SUPPRESSOR_TRNA}, {DOI_W1282X_RESCUE_BARRIER}, {DOI_ETI_PHASE3_TRIAL}")
            print(f"Structural Grounding: PDB {PDB_DEPHOSPHO_CFTR}, {PDB_ACTIVE_CFTR}")

            d1 = results["domain_1_hypoxic_isr_attenuation"]["po2_sweep"]
            print("\n[Domain 1] Hypoxic Integrated Stress Response (ISR) Initiation Attenuation:")
            for p in d1[:4]:
                print(f"  pO2 = {p['po2_percent']:4.1f}% O2 | k_init = {p['effective_initiation_rate_hz']:.4f} Hz ({p['initiation_attenuation_ratio'] * 100:.1f}% normoxic) | eIF2a-P: {p['eif2alpha_phosphorylation_index']:.2f}")

            d2 = results["domain_2_ribosome_collision_and_ngd"]["normoxia_comparison"]
            print(f"\n[Domain 2] Ribosome Collision Avoidance & No-Go Decay (NGD):")
            print(f"  Unmodified WT Disome Density:  {d2['unmodified_wildtype']['disome_density_per_mrna']:.4f} / mRNA")
            print(f"  Codon-Optimized m1Psi Density: {d2['codon_optimized_m1psi']['disome_density_per_mrna']:.4f} / mRNA")
            print(f"  Collision Reduction:           {d2['disome_reduction_pct']:.1f}% reduction via codon-optimization & m1Psi")

            d3_norm = results["domain_3_translational_yield_and_channel_density"]["normoxia_dose_response"]
            d3_hyp = results["domain_3_translational_yield_and_channel_density"]["hypoxia_dose_response"]
            print(f"\n[Domain 3] Mature CFTR Channel Density per Cell (Normoxia vs Hypoxia 2.0% O2):")
            for n, h in zip(d3_norm[2:6], d3_hyp[2:6]):
                print(f"  {n['mrna_copy_number']:2d} mRNAs/cell | Normoxia: {n['steady_state_surface_channels']:6.1f} channels ({n['percent_wild_type_cftr']:4.1f}% WT) | Hypoxia: {h['steady_state_surface_channels']:6.1f} channels ({h['percent_wild_type_cftr']:4.1f}% WT)")

            d4 = results["domain_4_airway_population_curative_rescue"]
            print(f"\n[Domain 4] Airway Population Curative Threshold Rescue (>= 10% WT / 150 channels):")
            print(f"  Unmodified WT Formulation:      {d4['unmodified_formulation']['overall_curative_rescue_percentage']:.1f}% cells rescued")
            print(f"  Codon-Optimized m1Psi Form.:    {d4['codon_optimized_m1psi_formulation']['overall_curative_rescue_percentage']:.1f}% cells rescued")
            print(f"  Formulation Rescue Benefit:     +{d4['rescue_differential_pct']:.1f}% absolute rescue increase")
            print("=" * 76)
        return 0

    if args.simulate_population:
        pop_res = engine.simulate_airway_population_rescue(
            mean_lnp_dose=20.0,
            codon_optimized=not args.unmodified,
            m1psi_modified=not args.unmodified,
        )
        print(json.dumps(pop_res, indent=2) if args.json else pop_res)
        return 0

    po2_val = args.po2 if args.po2 is not None else 21.0
    copies_val = args.mrna_copies if args.mrna_copies is not None else 10
    single_res = engine.calculate_single_cell_yield(
        mrna_copy_number=copies_val,
        po2_percent=po2_val,
        codon_optimized=not args.unmodified,
        m1psi_modified=not args.unmodified,
    )
    print(json.dumps(asdict(single_res), indent=2) if args.json else single_res)
    return 0


if __name__ == "__main__":
    sys.exit(main())
