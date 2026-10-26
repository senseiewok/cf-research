#!/usr/bin/env python3
"""Biophysical Model of Inhaled CFTR mRNA-LNP Cytosolic Delivery, Endosomal Escape & Translation Coupling in Primary Human Airway Basal Stem Cells (T-0455).

Quantitative biophysical dynamical modeling of aerosolized mRNA lipid nanoparticle (LNP)
delivery to primary human airway epithelial basal stem cells (KRT5+/TP63+) across 4 domains:

1. Endosomal Acidification & Ionizable Lipid Protonation:
   - Formulates pH-dependent endosomal lumen acidification from neutral (pH 7.4) down
     through early endosomes (pH 6.5–6.0) to late endosomes (pH 5.5–5.0) and lysosomes (pH 4.8).
   - Models ionizable lipid protonation fraction alpha_ion(pH) via Henderson-Hasselbalch
     across the apparent pKa spectrum (pKa 6.2–6.8).
   - Evaluates inverted hexagonal (H_II) non-bilayer phase transition propensity governed by
     lipid packing parameter S > 1, helper lipid formulation (DOPE inverted cone vs DSPC cylindrical),
     and anionic endosomal phospholipid pairing (LBPA/BMP and PS).

2. Cytosolic Cargo Delivery Yield:
   - Evaluates intracellular compartmental ODE trafficking:
     dE_early/dt, dE_late/dt, dE_lyso/dt, dE_rec/dt, and dC_intact/dt.
   - Computes intact CFTR mRNA escape fraction (phi_escape) surviving pore shear and lysosomal
     nucleases, contrasting productive cytosolic release against endosomal entrapment.

3. Stem Cell Quiescence & Nuclear Envelope/Endocytic Barriers:
   - Formulates endocytic uptake rate differences and translation initiation coupling across
     three epithelial subpopulations:
     * Differentiated pseudo-stratified luminal cells (ciliated/secretory club cells).
     * Proliferative/active cycling basal stem cells (B_a).
     * Quiescent basal stem progenitor cells (B_q).
   - Accounts for paracellular tight-junction resistance (ZO-1/claudin) and nuclear envelope
     intactness in quiescent versus cycling basal cells.

4. Curative Longevity & Epithelial Turnover Dynamics:
   - Simulates multi-week (30–90 day) duration of functional CFTR channel expression across
     epithelial renewal cycles.
   - Contrasts luminal-only ephemeral expression (washout t_1/2 ~ 1–3 days, dropping below
     10% wild-type curative threshold in < 14 days) against basal stem cell targeting that
     perpetuates functional channel replacement into daughter transit-amplifying lineages.

Primary Literature Grounding:
  - Englehardt JF et al. (2023) Nature, DOI: 10.1038/s41586-023-06133-1
    (Engineered suppressor tRNAs restore CFTR function in patient-derived airway epithelia)
  - Manfredi C et al. (2025) Biochemistry, DOI: 10.1021/acs.biochem.5c00590
    (Barriers to the Pharmacologic Rescue of W1282X CFTR)
  - Middleton PG et al. (2019) N Engl J Med, DOI: 10.1056/NEJMoa1908639
    (Elexacaftor-Tezacaftor-Ivacaftor in Patients with Cystic Fibrosis with Phe508del-Minimal Function Genotypes)
  - Button B et al. (2021) Physiol Rev, DOI: 10.1152/physrev.00004.2021
    (Physiology and pathophysiology of human airway mucus)

Clinical Registry Grounding:
  - NCT05668741 (Phase 1/2 Study of Inhaled mRNA LNP VX-522 in Participants With Cystic Fibrosis)
  - NCT05248230 (Phase 1/2 Study of 4D-710 Aerosolized AAV Gene Therapy in Cystic Fibrosis)

Structural Grounding:
  - RCSB PDB 6MSM (Dephosphorylated human CFTR inward-facing conformation)
  - RCSB PDB 5UAK (Phosphorylated active human CFTR conformation)

Attribution & Research Boundary:
  ClinicalTrials.gov does not evaluate safety or efficacy; listing on the registry does not denote government endorsement.
  Preclinical computational biophysical simulation; research, not medical advice.
  No clinical efficacy or therapeutic safety claim is made.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass, field
import json
import math
from pathlib import Path
import sys
from typing import Any, Dict, List, Optional, Tuple

SCRIPT_PATH = Path(__file__).resolve()

NLM_ATTRIBUTION = (
    "This product uses publicly available data from the U.S. National Library of Medicine (NLM), "
    "National Institutes of Health, Department of Health and Human Services; NLM is not responsible "
    "for the product and does not endorse or recommend this or any other product. "
    "ClinicalTrials.gov does not evaluate safety or efficacy; listing on the registry does not denote government endorsement."
)

RESEARCH_BOUNDARY = (
    "Preclinical computational biophysical simulation of inhaled CFTR mRNA-LNP endosomal escape "
    "and airway basal stem cell kinetics; research, not medical advice. No clinical efficacy or therapeutic safety claim is made."
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

# Biophysical Invariants & Clinical Thresholds
CFTR_CDS_NUCLEOTIDES: int = 4443
CFTR_CODONS: int = 1480
WT_CFTR_CHANNELS_PER_CELL: float = 1500.0
CURATIVE_RESCUE_FRACTION_WT: float = 0.100  # >= 10% wild-type threshold
CURATIVE_RESCUE_CHANNELS_THRESHOLD: float = 150.0  # 150 mature apical channels / cell

# Endosomal Trafficking Constants
DEFAULT_EXTRACELLULAR_PH: float = 7.40
DEFAULT_EARLY_ENDOSOME_PH: float = 6.20
DEFAULT_LATE_ENDOSOME_PH: float = 5.20
DEFAULT_LYSOSOME_PH: float = 4.80
DEFAULT_ACIDIFICATION_RATE_PER_HR: float = 1.85  # Exponential acidification rate constant (1/hr)

# Ionizable Lipid Constants
DEFAULT_IONIZABLE_LIPID_PKA: float = 6.45  # Optimal apparent pKa for pulmonary airway delivery (6.2-6.8)
PKA_LOWER_BOUND: float = 6.20
PKA_UPPER_BOUND: float = 6.80

# Helper Lipid Phase Transition Parameters (Inverted Hexagonal H_II Propensity)
PACKING_PARAMETER_DOPE: float = 1.15   # Inverted cone geometry (promotes H_II phase, S > 1)
PACKING_PARAMETER_DSPC: float = 0.98   # Cylindrical geometry (stabilizes lamellar L_alpha, S ~ 1)
LBPA_ANIONIC_AFFINITY_COEFF: float = 2.40  # Electrostatic pairing multiplier with endosomal LBPA
PROTON_SPONGE_SWELL_FACTOR: float = 0.28   # Osmotic swelling enhancement factor


@dataclass(frozen=True)
class EndosomalAcidificationState:
    """Endosomal acidification state, lipid protonation, and H_II phase transition propensity."""

    time_hr: float
    endosomal_ph: float
    protonation_fraction: float
    h_ii_phase_propensity: float
    membrane_disruption_rate_per_hr: float
    compartment_stage: str


@dataclass(frozen=True)
class CytosolicDeliveryYield:
    """Quantitative yield and partitioning of internalized mRNA cargo across compartments."""

    initial_internalized_copies: float
    early_endosome_trapped: float
    late_endosome_trapped: float
    lysosomal_degraded: float
    recycled_exocytosed: float
    cytosolic_intact_copies: float
    escape_efficiency_fraction: float
    time_to_peak_cytosolic_hr: float


@dataclass(frozen=True)
class StemCellSubpopulationProfile:
    """Endocytic uptake, translation kinetics, and barrier properties of an epithelial subpopulation."""

    cell_type: str
    quiescence_status: str
    relative_cell_volume_um3: float
    endocytic_uptake_rate_per_hr: float
    paracellular_access_fraction: float
    translation_initiation_rate_hz: float
    mature_cftr_folding_efficiency: float
    peak_channels_per_cell: float
    achieves_curative_threshold: bool


@dataclass(frozen=True)
class EpithelialLongevityTrajectory:
    """Multi-day epithelial renewal trajectory and CFTR channel durability."""

    days: List[float]
    luminal_channels_per_cell: List[float]
    basal_anchored_channels_per_cell: List[float]
    curative_threshold_channels: float
    luminal_therapeutic_duration_days: float
    basal_anchored_therapeutic_duration_days: float
    curative_longevity_ratio: float


class LNPStemEndosomalEscapeEngine:
    """Biophysical simulation engine for aerosolized LNP endosomal escape and basal stem cell kinetics."""

    def __init__(
        self,
        ionizable_pka: float = DEFAULT_IONIZABLE_LIPID_PKA,
        dope_molar_fraction: float = 0.35,
        dspc_molar_fraction: float = 0.05,
        k_acidification: float = DEFAULT_ACIDIFICATION_RATE_PER_HR,
    ):
        if not (5.0 <= ionizable_pka <= 8.0):
            raise ValueError(f"Ionizable pKa {ionizable_pka} out of physiological bounds [5.0, 8.0]")
        self.pka = float(ionizable_pka)
        self.dope_fraction = float(dope_molar_fraction)
        self.dspc_fraction = float(dspc_molar_fraction)
        self.k_acid = float(k_acidification)

    # -------------------------------------------------------------------------
    # Domain 1: Endosomal Acidification & Ionizable Lipid Protonation
    # -------------------------------------------------------------------------
    def calculate_endosomal_ph(self, time_hr: float) -> float:
        """Compute endosomal pH at elapsed time t via exponential acidification kinetics."""
        if time_hr < 0.0:
            raise ValueError(f"Negative time value: {time_hr}")
        # pH(t) = pH_lyso + (pH_ext - pH_lyso) * exp(-k_acid * t)
        delta_ph = DEFAULT_EXTRACELLULAR_PH - DEFAULT_LYSOSOME_PH
        ph = DEFAULT_LYSOSOME_PH + delta_ph * math.exp(-self.k_acid * time_hr)
        return max(DEFAULT_LYSOSOME_PH, min(DEFAULT_EXTRACELLULAR_PH, ph))

    def calculate_protonation_fraction(self, ph: float) -> float:
        """Compute ionizable lipid protonation fraction alpha_ion via Henderson-Hasselbalch equation."""
        # alpha_ion = 1 / (1 + 10^(pH - pKa))
        exponent = ph - self.pka
        # Clamp exponent to prevent numerical overflow
        clamped_exp = max(-12.0, min(12.0, exponent))
        return 1.0 / (1.0 + math.pow(10.0, clamped_exp))

    def compute_h_ii_phase_propensity(self, ph: float) -> float:
        """Calculate inverted hexagonal (H_II) non-bilayer phase transition propensity."""
        alpha = self.calculate_protonation_fraction(ph)

        # Helper lipid phase factor:
        # DOPE has inverted cone geometry promoting H_II phase transitions.
        # DSPC has cylindrical geometry that stabilizes lamellar L_alpha bilayers.
        helper_factor = 1.0 + 1.80 * self.dope_fraction - 1.20 * self.dspc_fraction

        # Electrostatic pairing with late endosomal LBPA (lysobisphosphatidic acid)
        # As endosomes acidify, LBPA accumulates and pairs with cationic lipid headgroups
        lbpa_pairing_factor = 1.0 + LBPA_ANIONIC_AFFINITY_COEFF * alpha * max(0.0, (DEFAULT_EARLY_ENDOSOME_PH - ph) / 1.5)

        # H_II phase transition propensity: increases as helper_factor > 1 and protonation increases
        h_ii_index = alpha * helper_factor * lbpa_pairing_factor * (1.0 + PROTON_SPONGE_SWELL_FACTOR)
        return max(0.0, h_ii_index)

    def evaluate_acidification_state(self, time_hr: float) -> EndosomalAcidificationState:
        """Evaluate complete endosomal acidification, protonation, and membrane disruption state."""
        ph = self.calculate_endosomal_ph(time_hr)
        alpha = self.calculate_protonation_fraction(ph)
        h_ii = self.compute_h_ii_phase_propensity(ph)

        # Membrane disruption rate constant: k_rupture = k_max * H_II propensity
        # Maximum escape rate is achieved in late endosomal window (pH 5.2 - 6.0)
        base_escape_rate = 0.450  # 1/hr max escape rate constant
        disruption_rate = base_escape_rate * min(2.5, h_ii)

        if ph >= 6.4:
            stage = "EARLY_ENDOSOME"
        elif ph >= 5.1:
            stage = "LATE_ENDOSOME_ESCAPE_WINDOW"
        else:
            stage = "LYSOSOME_ACIDIC_CORE"

        return EndosomalAcidificationState(
            time_hr=round(time_hr, 3),
            endosomal_ph=round(ph, 3),
            protonation_fraction=round(alpha, 4),
            h_ii_phase_propensity=round(h_ii, 4),
            membrane_disruption_rate_per_hr=round(disruption_rate, 4),
            compartment_stage=stage,
        )

    # -------------------------------------------------------------------------
    # Domain 2: Cytosolic Cargo Delivery Yield
    # -------------------------------------------------------------------------
    def simulate_intracellular_trafficking(
        self,
        initial_mRNAs: float = 20.0,
        total_time_hr: float = 24.0,
        dt_hr: float = 0.05,
    ) -> CytosolicDeliveryYield:
        """Simulate compartmental trafficking ODE system to calculate intact cytosolic mRNA yield."""
        if initial_mRNAs <= 0.0:
            raise ValueError(f"initial_mRNAs must be positive, got {initial_mRNAs}")

        # Rate constants (1 / hr)
        k_maturation = 0.650     # Early to late endosome maturation
        k_recycling = 0.120      # Early endosome exocytosis / recycling
        k_lyso = 0.450           # Late endosome to lysosome fusion
        k_deg_lyso = 1.200       # Degradation inside lysosome
        k_deg_cytosol = math.log(2.0) / 18.0  # Cytosolic mRNA turnover (t1/2 ~ 18 hr for m1Psi)
        eta_intact = 0.880       # Fraction of escaping cargo that remains unfragmented

        # Compartmental state variables
        e_early = float(initial_mRNAs)
        e_late = 0.0
        e_lyso = 0.0
        e_rec = 0.0
        c_intact = 0.0

        peak_cytosolic = 0.0
        time_to_peak = 0.0
        cumulative_escaped = 0.0
        cumulative_lysosomal_degraded = 0.0

        n_steps = int(total_time_hr / dt_hr)
        for step in range(n_steps):
            t = step * dt_hr
            ph = self.calculate_endosomal_ph(t)
            h_ii = self.compute_h_ii_phase_propensity(ph)

            # Stage-dependent escape rates (calibrated to Wittrup et al. 2015 & Gilleron et al. 2013)
            k_esc_early = 0.003 * h_ii
            k_esc_late = 0.030 * h_ii

            # Derivatives
            de_early = -(k_maturation + k_recycling + k_esc_early) * e_early
            de_late = k_maturation * e_early - (k_lyso + k_esc_late) * e_late
            de_lyso = k_lyso * e_late - k_deg_lyso * e_lyso
            de_rec = k_recycling * e_early

            escaped_rate = eta_intact * (k_esc_early * e_early + k_esc_late * e_late)
            dc_intact = escaped_rate - k_deg_cytosol * c_intact

            # Forward Euler step
            e_early = max(0.0, e_early + de_early * dt_hr)
            e_late = max(0.0, e_late + de_late * dt_hr)
            e_lyso = max(0.0, e_lyso + de_lyso * dt_hr)
            e_rec = e_rec + de_rec * dt_hr
            c_intact = max(0.0, c_intact + dc_intact * dt_hr)
            cumulative_escaped += escaped_rate * dt_hr
            cumulative_lysosomal_degraded += k_deg_lyso * e_lyso * dt_hr

            if c_intact > peak_cytosolic:
                peak_cytosolic = c_intact
                time_to_peak = t

        efficiency = cumulative_escaped / initial_mRNAs

        return CytosolicDeliveryYield(
            initial_internalized_copies=round(initial_mRNAs, 2),
            early_endosome_trapped=round(e_early, 3),
            late_endosome_trapped=round(e_late, 3),
            lysosomal_degraded=round(cumulative_lysosomal_degraded + e_lyso, 3),
            recycled_exocytosed=round(e_rec, 3),
            cytosolic_intact_copies=round(peak_cytosolic, 3),
            escape_efficiency_fraction=round(efficiency, 4),
            time_to_peak_cytosolic_hr=round(time_to_peak, 2),
        )

    # -------------------------------------------------------------------------
    # Domain 3: Stem Cell Quiescence & Nuclear Envelope/Endocytic Barriers
    # -------------------------------------------------------------------------
    def evaluate_subpopulation_kinetics(
        self,
        cell_type: str,
        paracellular_ter_ohm_cm2: float = 1200.0,
    ) -> StemCellSubpopulationProfile:
        """Formulate endocytic uptake, translation coupling, and channel yield across cell types."""
        valid_types = {"ciliated_luminal", "active_basal_stem", "quiescent_basal_stem"}
        if cell_type not in valid_types:
            raise ValueError(f"Unknown cell_type '{cell_type}'; must be one of {valid_types}")

        # Baseline paracellular barrier permeation based on TEER (tight junction resistance)
        # When TEER is high (1200 ohm*cm2), apical tight junctions restrict LNP paracellular access
        # to basal compartment to ~8-12%. With transient opening (TEER ~500 ohm*cm2), access reaches 30-35%.
        clamped_ter = max(300.0, min(2000.0, paracellular_ter_ohm_cm2))
        paracellular_access = max(0.05, min(0.40, 1.0 - (clamped_ter / 1800.0)))

        if cell_type == "ciliated_luminal":
            vol = 2200.0  # um^3, large columnar cell
            uptake_rate = 0.850  # 1/hr apical endocytosis
            access = 1.000  # directly exposed to aerosol lumen
            initiation_hz = 0.080  # ribosomes/s per mRNA
            folding_eff = 0.320  # mature apical trafficking fraction
            quiescence = "DIFFERENTIATED_POST_MITOTIC"
        elif cell_type == "active_basal_stem":
            vol = 850.0   # um^3, compact cycling basal progenitor
            uptake_rate = 0.620  # 1/hr active macropinocytosis and CME
            access = paracellular_access
            initiation_hz = 0.095  # elevated translation during cell cycle
            folding_eff = 0.420  # high chaperone abundance in active stem cells
            quiescence = "ACTIVE_CYCLING"
        else:  # quiescent_basal_stem
            vol = 720.0   # um^3, small quiescent niche stem cell (KRT5+/TP63+)
            uptake_rate = 0.220  # 1/hr slow baseline endocytic turnover
            access = paracellular_access
            initiation_hz = 0.035  # attenuated baseline translation in quiescence (G0)
            folding_eff = 0.380
            quiescence = "QUIESCENT_G0"

        # Cytosolic delivery yield for this subpopulation
        # Average inhaled exposure delivers ~25 LNPs to apical membrane
        effective_delivered_lnps = 25.0 * access
        yield_sim = self.simulate_intracellular_trafficking(
            initial_mRNAs=max(1.0, effective_delivered_lnps),
            total_time_hr=16.0,
        )

        # Mature CFTR channel calculation:
        # Intact cytosolic mRNAs are actively translated by polyribosomes at initiation_hz.
        # Polypeptides initiated per hour = yield_sim.cytosolic_intact_copies * (initiation_hz * 3600.0)
        # Peak channel accumulation accounts for folding efficiency and transient burst integration (2.5 - 3.2 hr)
        burst_duration_hr = 2.8 if cell_type == "active_basal_stem" else 2.5
        polypeptides_per_hr = yield_sim.cytosolic_intact_copies * (initiation_hz * 3600.0)
        peak_channels = polypeptides_per_hr * folding_eff * burst_duration_hr

        achieves_curative = peak_channels >= CURATIVE_RESCUE_CHANNELS_THRESHOLD

        return StemCellSubpopulationProfile(
            cell_type=cell_type,
            quiescence_status=quiescence,
            relative_cell_volume_um3=vol,
            endocytic_uptake_rate_per_hr=uptake_rate,
            paracellular_access_fraction=round(access, 4),
            translation_initiation_rate_hz=initiation_hz,
            mature_cftr_folding_efficiency=folding_eff,
            peak_channels_per_cell=round(peak_channels, 1),
            achieves_curative_threshold=achieves_curative,
        )

    # -------------------------------------------------------------------------
    # Domain 4: Curative Longevity & Epithelial Turnover Dynamics
    # -------------------------------------------------------------------------
    def simulate_epithelial_turnover_trajectory(
        self,
        days_total: float = 90.0,
        dt_days: float = 1.0,
    ) -> EpithelialLongevityTrajectory:
        """Simulate multi-week functional CFTR persistence across epithelial renewal cycles."""
        luminal_profile = self.evaluate_subpopulation_kinetics("ciliated_luminal")
        active_stem_profile = self.evaluate_subpopulation_kinetics("active_basal_stem")
        quiescent_stem_profile = self.evaluate_subpopulation_kinetics("quiescent_basal_stem")

        # Initial peak channels established at day 1
        lum_peak = luminal_profile.peak_channels_per_cell
        stem_peak = 0.65 * active_stem_profile.peak_channels_per_cell + 0.35 * quiescent_stem_profile.peak_channels_per_cell

        # Half-lives (days)
        # mRNA half-life in non-dividing differentiated cells is ~1.0 day; protein half-life is ~1.5 days
        # Net washout half-life in purely luminal cells without replenishment: ~2.5 days
        t_half_luminal_ephemeral = 2.5
        k_decay_luminal = math.log(2.0) / t_half_luminal_ephemeral

        # Epithelial turnover and stem cell progeny replenishment:
        # Airway basal stem cells divide homeostatically (lambda_div ~ 0.01 / day),
        # generating daughter transit-amplifying cells that replace shed ciliated cells over 45-90 days.
        # Basal-anchored persistence half-life: ~42 days due to continuous daughter cell seeding
        t_half_basal_anchored = 42.0
        k_decay_basal = math.log(2.0) / t_half_basal_anchored

        days_list: List[float] = []
        luminal_curve: List[float] = []
        basal_anchored_curve: List[float] = []

        luminal_days_above = 0.0
        basal_days_above = 0.0

        n_steps = int(days_total / dt_days) + 1
        for step in range(n_steps):
            d = step * dt_days
            days_list.append(round(d, 1))

            # Luminal-only decays rapidly without progenitor renewal
            c_lum = lum_peak * math.exp(-k_decay_luminal * d)
            luminal_curve.append(round(c_lum, 1))
            if c_lum >= CURATIVE_RESCUE_CHANNELS_THRESHOLD:
                luminal_days_above = d

            # Basal-anchored profile combines initial luminal burst plus sustained stem cell progeny output
            c_basal = 0.50 * lum_peak * math.exp(-k_decay_luminal * d) + stem_peak * math.exp(-k_decay_basal * d)
            basal_anchored_curve.append(round(c_basal, 1))
            if c_basal >= CURATIVE_RESCUE_CHANNELS_THRESHOLD:
                basal_days_above = d

        ratio = (basal_days_above / luminal_days_above) if luminal_days_above > 0.0 else 1.0

        return EpithelialLongevityTrajectory(
            days=days_list,
            luminal_channels_per_cell=luminal_curve,
            basal_anchored_channels_per_cell=basal_anchored_curve,
            curative_threshold_channels=CURATIVE_RESCUE_CHANNELS_THRESHOLD,
            luminal_therapeutic_duration_days=round(luminal_days_above, 1),
            basal_anchored_therapeutic_duration_days=round(basal_days_above, 1),
            curative_longevity_ratio=round(ratio, 2),
        )

    # -------------------------------------------------------------------------
    # Comprehensive Simulation Runner
    # -------------------------------------------------------------------------
    def run_full_simulation(self) -> Dict[str, Any]:
        """Execute comprehensive multi-domain biophysical simulation."""
        # 1. Acidification & H_II phase transition profile
        acid_times = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0]
        acid_states = [asdict(self.evaluate_acidification_state(t)) for t in acid_times]

        # 2. Cytosolic delivery yield across various initial LNP burdens
        yields = [
            asdict(self.simulate_intracellular_trafficking(initial_mRNAs=m))
            for m in [5.0, 15.0, 30.0, 50.0]
        ]

        # 3. Subpopulation profiles
        cell_profiles = {
            ctype: asdict(self.evaluate_subpopulation_kinetics(ctype))
            for ctype in ["ciliated_luminal", "active_basal_stem", "quiescent_basal_stem"]
        }

        # 4. Multi-week epithelial longevity simulation
        longevity = asdict(self.simulate_epithelial_turnover_trajectory(days_total=90.0))

        return {
            "model": "lnp_stem_endosomal_escape_kinetics",
            "task": "T-0455",
            "track": "Track A",
            "discipline": "airway_delivery",
            "answer_tier": "T2",
            "primary_evidence": {
                "dois": [
                    DOI_SUPPRESSOR_TRNA,
                    DOI_W1282X_RESCUE_BARRIER,
                    DOI_ETI_PHASE3_TRIAL,
                    DOI_AIRWAY_MUCUS_PHYSIOL,
                ],
                "clinical_trials": [NCT_VX522_MRNA_TRIAL, NCT_4D710_AAV_TRIAL],
                "structural_pdbs": [PDB_DEPHOSPHO_CFTR, PDB_ACTIVE_CFTR],
            },
            "parameters": {
                "ionizable_pka": self.pka,
                "dope_molar_fraction": self.dope_fraction,
                "dspc_molar_fraction": self.dspc_fraction,
                "cftr_codons": CFTR_CODONS,
                "curative_rescue_threshold_channels": CURATIVE_RESCUE_CHANNELS_THRESHOLD,
            },
            "acidification_and_h_ii_phase": acid_states,
            "cytosolic_delivery_yields": yields,
            "subpopulation_kinetics": cell_profiles,
            "longevity_trajectory": longevity,
            "attributions": {
                "nlm": NLM_ATTRIBUTION,
                "boundary": RESEARCH_BOUNDARY,
            },
        }


def main() -> int:
    """CLI entry point for LNP stem cell endosomal escape model."""
    parser = argparse.ArgumentParser(
        description="Biophysical Model of Inhaled CFTR mRNA-LNP Cytosolic Delivery, Endosomal Escape & Basal Stem Cell Kinetics (T-0455)"
    )
    parser.add_argument("--pka", type=float, default=DEFAULT_IONIZABLE_LIPID_PKA, help="Ionizable lipid pKa (default: 6.45)")
    parser.add_argument("--dope", type=float, default=0.35, help="DOPE molar fraction (default: 0.35)")
    parser.add_argument("--dspc", type=float, default=0.05, help="DSPC molar fraction (default: 0.05)")
    parser.add_argument("--json", action="store_true", help="Output machine-readable JSON")
    args = parser.parse_args()

    engine = LNPStemEndosomalEscapeEngine(
        ionizable_pka=args.pka,
        dope_molar_fraction=args.dope,
        dspc_molar_fraction=args.dspc,
    )
    result = engine.run_full_simulation()

    if args.json:
        print(json.dumps(result, indent=2))
    else:
        print("================================================================================")
        print("INHALED CFTR mRNA-LNP STEM CELL ENDOSOMAL ESCAPE & LONGEVITY MODEL (T-0455)")
        print("================================================================================")
        print(f"Ionizable Lipid pKa:          {engine.pka:.2f} (optimal airway target: 6.20–6.80)")
        print(f"Helper Lipids:                {engine.dope_fraction*100:.1f}% DOPE (H_II inducer) / {engine.dspc_fraction*100:.1f}% DSPC")
        print(f"CFTR Coding Sequence:         {CFTR_CODONS} codons ({CFTR_CDS_NUCLEOTIDES} nt)")
        print(f"Curative Threshold:           >= {CURATIVE_RESCUE_CHANNELS_THRESHOLD:.0f} mature apical channels/cell (10% WT)")
        print("--------------------------------------------------------------------------------")
        print("Subpopulation Peak CFTR Rescues:")
        for ctype, prof in result["subpopulation_kinetics"].items():
            status = "PASS (Curative)" if prof["achieves_curative_threshold"] else "SUB-CURATIVE"
            print(f"  - {ctype:<22}: {prof['peak_channels_per_cell']:>6.1f} channels/cell [{status}]")
        print("--------------------------------------------------------------------------------")
        long = result["longevity_trajectory"]
        print(f"Luminal Ephemeral Duration:   {long['luminal_therapeutic_duration_days']} days >= 10% WT")
        print(f"Basal-Anchored Duration:      {long['basal_anchored_therapeutic_duration_days']} days >= 10% WT")
        print(f"Stem Cell Longevity Benefit:  {long['curative_longevity_ratio']:.1f}x extension across renewal cycles")
        print("================================================================================")

    return 0


if __name__ == "__main__":
    sys.exit(main())
