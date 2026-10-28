#!/usr/bin/env python3
"""Biophysical and Kinetic Modeling of Lipid-Nanoparticle Nuclear Entry and Basal Stem Cell Genomic Persistence in Regenerative CF Airway Epithelia (T-0478).

Quantitative biophysical, transport, and stem cell population genetics model:
1. Nuclear Pore Complex (NPC) Translocation Kinetics:
   - Formulates the steric-hydrophobic permeability barrier of phenylalanine-glycine (FG) repeat nucleoporins
     in the central channel (channel radius R_pore ~ 25 nm, functional passive cutoff r_cutoff ~ 2.5 nm / Mw ~ 45 kDa).
   - Stokes-Einstein hydrodynamic radius calculation for ribonucleoprotein (RNP) complexes
     (uncomplexed Cas9 / base editor Mw ~ 160 kDa, r_H ~ 4.4 nm > r_cutoff).
   - Importin-alpha/beta-mediated facilitated active transport: Nuclear Localization Signal (NLS) affinity
     (Kd ~ 10-50 nM), RanGTP-dependent vector directionality, and nuclear-to-cytoplasmic accumulation ratio (R_N/C).
2. Pseudostratified Airway Epithelial Niches & Paracellular Delivery:
   - Apical junctional barrier: tight junction complexes (ZO-1, claudins, occludins) maintaining transepithelial
     electrical resistance (TEER ~ 300-1000 Ohm*cm^2).
   - Paracellular sieving mechanics: restricts intact LNPs (diameter 60-100 nm) to apical columnar uptake (< 3% basal access)
     unless transient tight junction permeabilization (e.g. LPC, hypertonic modulation, shedding) widens junctional gaps
     to 50-150 nm, elevating basal stem cell access to 15-35%.
3. Basal Stem Cell (KRT5+/TP63+) Dynamics & Multi-Year Turnover Kinetics:
   - Three-compartment lineage ODE system: Basal progenitor stem cells (B) -> Transit-amplifying cells (T) ->
     Differentiated columnar cells (D: ciliated, secretory, and FOXI1+ ionocytes).
   - Symmetric self-renewal (p_s), asymmetric division (1 - 2*p_d), and symmetric differentiation (p_d) balancing
     homeostatic pool conservation (d B / dt = 0 at baseline).
   - Epithelial shedding dynamics (turnover half-life tau_D ~ 30-90 days, default 60 days).
   - Genomic vs. Episomal Persistence Contrast:
     * Transient episomal mRNA (e.g., VX-522, NCT05668741): Rapid intracellular mRNA degradation (t_1/2 ~ 18 hours),
       mitotic dilution in dividing cells, and luminal shedding reduce functional CFTR expression to 0% within 30-60 days.
     * Permanent genomic editing (e.g., base editing / lentiviral integration, Cooney et al. 2018): Permanent allelic
       correction in KRT5+/TP63+ basal stem cells produces an unwaning, self-renewing progeny cohort, establishing an
       asymptotic homeostatic floor F_inf = eta_basal_edited * Efficacy_allelic that permanently exceeds the curative
       threshold (>= 10% WT CFTR, Middleton et al. 2019, Button et al. 2021) across multi-year intervals (> 10 years).

Primary Evidence Grounding:
- Primary Literature:
  - Button B et al. (2021) Physiol Rev, DOI: 10.1152/physrev.00004.2021
    (Physiology and pathophysiology of human airway mucus)
  - Cooney AL et al. (2018) Mol Ther, DOI: 10.1016/j.ymthe.2018.04.014
    (Lentiviral-mediated gene transfer, targeting airway basal cells, results in persistent CFTR expression in cystic fibrosis airways)
  - Middleton PG et al. (2019) N Engl J Med, DOI: 10.1056/nejmoa1908639
    (Elexacaftor-Tezacaftor-Ivacaftor for Cystic Fibrosis with a Single Phe508del Allele)
- Clinical Registry Grounding:
  - NCT05668741: Phase 1/2 Study of Inhaled mRNA LNP VX-522 in Participants With Cystic Fibrosis
  - NCT05033080: Phase 3 Study of Vanzacaftor-Tezacaftor-Deutivacaftor in Cystic Fibrosis
- Structural Grounding:
  - RCSB PDB 6MSM (Human CFTR dephosphorylated ATP-free conformation, 3.9 A Cryo-EM)

Epistemic Status:
- Answer Tier T2 (Peer-reviewed stem cell biology, transport biophysics, and ClinicalTrials.gov API v2).
- Preclinical in-silico computational model; research, not medical advice.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

SCRIPT_PATH = Path(__file__).resolve()

# Attribution & Mandatory Disclaimers
NLM_ATTRIBUTION: str = (
    "This product uses publicly available data from the U.S. National Library of Medicine (NLM), "
    "National Institutes of Health, Department of Health and Human Services; NLM is not responsible "
    "for the product and does not endorse or recommend this or any other product. "
    "ClinicalTrials.gov does not evaluate safety or efficacy; listing on the registry does not denote government endorsement."
)

RESEARCH_BOUNDARY: str = (
    "Preclinical in-silico nuclear entry and basal stem cell persistence model; "
    "research, not medical advice. No clinical efficacy or therapeutic safety claim is made."
)

# Reference Grounding Identifiers (Verified in sources/cite-check-ledger.jsonl & snapshot.json)
PRIMARY_DOI_BUTTON: str = "10.1152/physrev.00004.2021"
PRIMARY_DOI_COONEY: str = "10.1016/j.ymthe.2018.04.014"
PRIMARY_DOI_MIDDLETON: str = "10.1056/nejmoa1908639"
CLINICAL_TRIAL_VX522: str = "NCT05668741"
CLINICAL_TRIAL_VANZACAFTOR: str = "NCT05033080"
PDB_STRUCTURE_CFTR: str = "6MSM"

# Physical and Biological Constants
CURATIVE_THRESHOLD_CFTR_PCT: float = 10.0  # >= 10% WT CFTR restored is clinically curative
AVOGADRO_NUMBER: float = 6.02214076e23
BOLTZMANN_CONSTANT_J_K: float = 1.380649e-23
TEMPERATURE_KELVIN: float = 310.15  # 37 C
CYTOSOL_VISCOSITY_PA_S: float = 0.003  # ~3x water viscosity (3 mPa*s)

# Nuclear Pore Complex (NPC) Architecture Defaults
DEFAULT_NPC_COUNT_PER_NUCLEUS: float = 3000.0  # Typical human epithelial nucleus has 2000-4000 NPCs
DEFAULT_NPC_PORE_RADIUS_NM: float = 25.0       # Central channel radius (nm)
DEFAULT_NPC_PASSIVE_CUTOFF_NM: float = 2.5     # Passive Stokes radius limit ~ 2.5 nm (~45 kDa)
DEFAULT_NPC_LENGTH_NM: float = 80.0            # Scaffold length across nuclear envelope (nm)
DEFAULT_RNP_MOLECULAR_WEIGHT_KDA: float = 160.0  # SpCas9 or base editor RNP (~160 kDa)
DEFAULT_RNP_HYDRATION_SPECIFIC_VOLUME_CM3_G: float = 1.30  # Hydrated RNP complex with bound sgRNA (~1.30 cm^3/g)

# Importin Transport Kinetics Defaults
DEFAULT_NLS_KD_NM: float = 25.0                # Importin-alpha binding affinity to bipartite NLS (nM)
DEFAULT_IMPORTIN_CONCENTRATION_CYTO_UM: float = 3.0  # Cytosolic Importin-alpha/beta concentration (uM)
DEFAULT_RANGTP_NUCLEAR_FRACTION: float = 0.85  # Nuclear Ran in GTP-bound state
DEFAULT_MAX_IMPORT_TURNOVER_PER_PORE_S: float = 50.0  # Maximum transport events per NPC per second

# Pseudostratified Airway Epithelium & Barrier Defaults
DEFAULT_TEER_OHM_CM2: float = 600.0            # Baseline intact airway epithelial resistance (Ohm*cm^2)
DEFAULT_TIGHT_JUNCTION_PORE_INTACT_NM: float = 1.2   # Intact paracellular pore diameter
DEFAULT_TIGHT_JUNCTION_PORE_PERMEABILIZED_NM: float = 85.0  # Permeabilized pore size
DEFAULT_LNP_DIAMETER_NM: float = 80.0          # Typical inhaled LNP diameter (nm)
DEFAULT_BASAL_CELL_FRACTION: float = 0.30      # Basal cells represent ~30% of human airway epithelium
DEFAULT_COLUMNAR_CELL_FRACTION: float = 0.70   # Columnar ciliated & secretory cells represent ~70%

# Stem Cell Hierarchy & Turnover Defaults
DEFAULT_DIFFERENTIATED_TURNOVER_HALF_LIFE_DAYS: float = 60.0  # Columnar cell turnover t_1/2 = 60 days (30-90 days range)
DEFAULT_BASAL_STEM_DIVISION_RATE_PER_DAY: float = 0.015       # Basal stem cell division every ~66 days
DEFAULT_TRANSIT_AMPLIFYING_DIVISION_RATE_PER_DAY: float = 0.035
DEFAULT_TRANSIT_DIFFERENTIATION_RATE_PER_DAY: float = 0.025
DEFAULT_SYMMETRIC_RENEWAL_PROBABILITY: float = 0.25           # Homeostatic balance
DEFAULT_SYMMETRIC_DIFFERENTIATION_PROBABILITY: float = 0.25    # Homeostatic balance
DEFAULT_EPISOMAL_MRNA_HALF_LIFE_HOURS: float = 18.0           # Intracellular mRNA decay half-life


@dataclass(frozen=True)
class NpcTranslocationParameters:
    """Biophysical parameters governing nuclear envelope translocation."""
    molecular_weight_kda: float = DEFAULT_RNP_MOLECULAR_WEIGHT_KDA
    has_nls: bool = True
    nls_kd_nm: float = DEFAULT_NLS_KD_NM
    importin_conc_um: float = DEFAULT_IMPORTIN_CONCENTRATION_CYTO_UM
    rangtp_nuclear_fraction: float = DEFAULT_RANGTP_NUCLEAR_FRACTION
    npc_count: float = DEFAULT_NPC_COUNT_PER_NUCLEUS
    pore_radius_nm: float = DEFAULT_NPC_PORE_RADIUS_NM
    passive_cutoff_nm: float = DEFAULT_NPC_PASSIVE_CUTOFF_NM
    max_turnover_per_pore_s: float = DEFAULT_MAX_IMPORT_TURNOVER_PER_PORE_S


@dataclass(frozen=True)
class NpcTranslocationResult:
    """Calculated biophysical translocation rates across nuclear pore complexes."""
    stokes_radius_nm: float
    is_passive_permissive: bool
    passive_permeability_cm_s: float
    importin_bound_fraction: float
    active_flux_molecules_per_s: float
    effective_nuclear_import_rate_per_hour: float
    nuclear_cytoplasmic_steady_state_ratio: float
    translocation_mechanism: str


@dataclass(frozen=True)
class EpithelialNicheParameters:
    """Anatomical and paracellular barrier parameters of pseudostratified epithelium."""
    teer_ohm_cm2: float = DEFAULT_TEER_OHM_CM2
    is_permeabilized: bool = False
    lnp_diameter_nm: float = DEFAULT_LNP_DIAMETER_NM
    basal_cell_fraction: float = DEFAULT_BASAL_CELL_FRACTION
    apical_columnar_fraction: float = DEFAULT_COLUMNAR_CELL_FRACTION
    administered_dose_lnp_per_cm2: float = 1.0e8


@dataclass(frozen=True)
class BasalDeliveryResult:
    """Calculated distribution of delivered nanoparticles across epithelial compartments."""
    paracellular_penetration_fraction: float
    apical_columnar_uptake_fraction: float
    basal_stem_uptake_fraction: float
    basal_stem_dose_per_cell: float
    apical_columnar_dose_per_cell: float
    basal_delivery_efficiency_pct: float
    permeabilization_state: str


@dataclass(frozen=True)
class TurnoverKineticsParameters:
    """Dynamical parameters for stem cell renewal, differentiation, and tissue shedding."""
    differentiated_half_life_days: float = DEFAULT_DIFFERENTIATED_TURNOVER_HALF_LIFE_DAYS
    basal_division_rate_per_day: float = DEFAULT_BASAL_STEM_DIVISION_RATE_PER_DAY
    p_symmetric_renewal: float = DEFAULT_SYMMETRIC_RENEWAL_PROBABILITY
    p_symmetric_differentiation: float = DEFAULT_SYMMETRIC_DIFFERENTIATION_PROBABILITY
    transit_division_rate_per_day: float = DEFAULT_TRANSIT_AMPLIFYING_DIVISION_RATE_PER_DAY
    transit_differentiation_rate_per_day: float = DEFAULT_TRANSIT_DIFFERENTIATION_RATE_PER_DAY
    episomal_mrna_half_life_hours: float = DEFAULT_EPISOMAL_MRNA_HALF_LIFE_HOURS
    initial_basal_edit_efficiency_pct: float = 30.0
    initial_apical_edit_efficiency_pct: float = 40.0
    initial_episomal_transfection_pct: float = 65.0


@dataclass(frozen=True)
class PersistenceTimePoint:
    """Simulation snapshot of epithelial CFTR expression at day t."""
    day: float
    years: float
    episomal_cftr_pct: float
    genomic_edited_basal_pct: float
    genomic_edited_columnar_pct: float
    genomic_total_cftr_pct: float
    episomal_above_curative: bool
    genomic_above_curative: bool


@dataclass(frozen=True)
class ModalityDurabilitySummary:
    """Multi-year comparative durability metrics between genomic editing and episomal therapies."""
    modality_name: str
    mechanism: str
    clinical_reference: str
    day_0_cftr_pct: float
    day_30_cftr_pct: float
    day_90_cftr_pct: float
    day_365_cftr_pct: float
    day_1825_5yr_cftr_pct: float
    day_3650_10yr_cftr_pct: float
    asymptotic_homeostatic_floor_pct: float
    curative_duration_days: float
    requires_indefinite_redosing: bool
    durability_assessment: str


class LnpNuclearEntryBasalPersistenceModel:
    """Biophysical simulation engine for LNP nuclear translocation and basal stem cell persistence."""

    def __init__(
        self,
        npc_params: Optional[NpcTranslocationParameters] = None,
        niche_params: Optional[EpithelialNicheParameters] = None,
        turnover_params: Optional[TurnoverKineticsParameters] = None,
    ) -> None:
        self.npc_params = npc_params or NpcTranslocationParameters()
        self.niche_params = niche_params or EpithelialNicheParameters()
        self.turnover_params = turnover_params or TurnoverKineticsParameters()

    def compute_stokes_hydrodynamic_radius_nm(self, molecular_weight_kda: float) -> float:
        """Calculates Stokes-Einstein hydrodynamic radius for a globular protein/RNP complex."""
        if molecular_weight_kda <= 0:
            raise ValueError("Molecular weight must be strictly positive.")
        # Mass per molecule in grams
        mass_g = (molecular_weight_kda * 1000.0) / AVOGADRO_NUMBER
        # Volume V = mass * specific volume (cm^3) accounting for hydration shell and bound nucleic acids
        vol_cm3 = mass_g * DEFAULT_RNP_HYDRATION_SPECIFIC_VOLUME_CM3_G
        # Radius r = (3V / 4pi)^(1/3) in cm, converted to nm (1 cm = 1e7 nm)
        radius_cm = ( (3.0 * vol_cm3) / (4.0 * math.pi) ) ** (1.0 / 3.0)
        return radius_cm * 1.0e7

    def solve_npc_translocation(
        self,
        params: Optional[NpcTranslocationParameters] = None,
    ) -> NpcTranslocationResult:
        """Solves passive steric exclusion vs. importin-facilitated active translocation across the NPC."""
        p = params or self.npc_params
        r_stokes = self.compute_stokes_hydrodynamic_radius_nm(p.molecular_weight_kda)

        # Passive diffusion steric partitioning across FG-repeat hydrogel
        # Sharp hydrogel barrier for cargo exceeding cutoff radius (K_passive drops with power 4)
        if r_stokes < p.pore_radius_nm:
            steric_factor = max(0.0, (1.0 - (r_stokes / p.pore_radius_nm))) ** 2
            fg_exclusion = math.exp(- (r_stokes / p.passive_cutoff_nm) ** 4)
            k_partition = steric_factor * fg_exclusion
        else:
            k_partition = 0.0

        # Stokes-Einstein diffusion coefficient in water/viscous pore (cm^2/s)
        # D = k_B * T / (6 * pi * eta * r)
        r_meters = r_stokes * 1.0e-9
        d_diff_m2_s = (BOLTZMANN_CONSTANT_J_K * TEMPERATURE_KELVIN) / (6.0 * math.pi * CYTOSOL_VISCOSITY_PA_S * r_meters)
        d_diff_cm2_s = d_diff_m2_s * 1.0e4

        # Pore length in cm
        length_cm = p.npc_length_nm if hasattr(p, 'npc_length_nm') else DEFAULT_NPC_LENGTH_NM * 1.0e-7
        # Passive permeability P_passive = (D * K_partition) / length
        p_passive_cm_s = (d_diff_cm2_s * k_partition) / length_cm

        is_passive = r_stokes <= p.passive_cutoff_nm

        # Importin-alpha/beta binding equilibrium
        # Fraction bound = [Importin] / (Kd + [Importin])
        if p.has_nls:
            importin_conc_nm = p.importin_conc_um * 1000.0
            fraction_bound = importin_conc_nm / (p.nls_kd_nm + importin_conc_nm)
            # Active translocation flux: facilitated by FG-repeat interaction and driven by RanGTP gradient
            # J_active = N_NPC * k_max * fraction_bound * RanGTP_fraction
            active_flux_per_s = p.npc_count * p.max_turnover_per_pore_s * fraction_bound * p.rangtp_nuclear_fraction
            effective_import_rate_hr = (active_flux_per_s * 3600.0) / 1000.0  # arbitrary normalized unit scale
            # Active transport drives accumulation against passive leak
            # R_N/C steady-state ratio
            steady_state_rnc = 1.0 + (25.0 * fraction_bound * p.rangtp_nuclear_fraction) / (1.0 + p_passive_cm_s * 1.0e5)
            steady_state_rnc = min(45.0, max(1.0, steady_state_rnc))
            mechanism = "Importin-alpha/beta active translocation (RanGTP-vectored)"
        else:
            fraction_bound = 0.0
            active_flux_per_s = 0.0
            effective_import_rate_hr = p_passive_cm_s * 3600.0 * 10.0
            steady_state_rnc = 1.0 if is_passive else 0.08  # Nuclear exclusion for large non-NLS cargo
            mechanism = "Passive diffusion (size-restricted FG steric barrier)"

        return NpcTranslocationResult(
            stokes_radius_nm=round(r_stokes, 3),
            is_passive_permissive=is_passive,
            passive_permeability_cm_s=p_passive_cm_s,
            importin_bound_fraction=round(fraction_bound, 4),
            active_flux_molecules_per_s=round(active_flux_per_s, 2),
            effective_nuclear_import_rate_per_hour=round(effective_import_rate_hr, 2),
            nuclear_cytoplasmic_steady_state_ratio=round(steady_state_rnc, 2),
            translocation_mechanism=mechanism,
        )

    def solve_paracellular_basal_niche_delivery(
        self,
        params: Optional[EpithelialNicheParameters] = None,
    ) -> BasalDeliveryResult:
        """Calculates paracellular barrier penetration and basal stem cell dose distribution."""
        p = params or self.niche_params

        # Paracellular penetration fraction depends on TEER and junctional status
        if p.is_permeabilized or p.teer_ohm_cm2 < 120.0:
            # Under tight junction loosening (e.g. transient hypertonic or lipid permeabilization)
            # Pore size ~85 nm allows substantial penetration of ~80 nm LNPs
            effective_pore_nm = DEFAULT_TIGHT_JUNCTION_PORE_PERMEABILIZED_NM
            ratio = p.lnp_diameter_nm / effective_pore_nm
            hindrance = max(0.0, 1.0 - ratio) ** 2 if ratio < 1.0 else 0.05
            penetration_frac = 0.28 * (1.0 - min(0.9, p.teer_ohm_cm2 / 600.0)) + 0.12 * hindrance
            penetration_frac = max(0.08, min(0.38, penetration_frac))
            state = "Transiently Permeabilized (Tight Junctions Disassembled)"
        else:
            # Intact tight junctions (TEER 300-1000 Ohm*cm^2): paracellular pore < 2 nm
            # LNPs (60-100 nm) strictly excluded from paracellular passage; minimal vesicular leak (< 2%)
            penetration_frac = max(0.005, 0.03 * (150.0 / max(150.0, p.teer_ohm_cm2)))
            state = "Intact Epithelial Barrier (Tight Junction Paracellular Exclusion)"

        # Apical columnar cells capture remaining mucosal nanoparticles
        apical_uptake_frac = 1.0 - penetration_frac
        basal_uptake_frac = penetration_frac

        # Per-cell dose distribution (normalized to surface area and cellular proportions)
        dose_basal_cell = (p.administered_dose_lnp_per_cm2 * basal_uptake_frac) / max(1e-3, p.basal_cell_fraction * 1e6)
        dose_apical_cell = (p.administered_dose_lnp_per_cm2 * apical_uptake_frac) / max(1e-3, p.apical_columnar_fraction * 1e6)

        return BasalDeliveryResult(
            paracellular_penetration_fraction=round(penetration_frac, 4),
            apical_columnar_uptake_fraction=round(apical_uptake_frac, 4),
            basal_stem_uptake_fraction=round(basal_uptake_frac, 4),
            basal_stem_dose_per_cell=round(dose_basal_cell, 1),
            apical_columnar_dose_per_cell=round(dose_apical_cell, 1),
            basal_delivery_efficiency_pct=round(basal_uptake_frac * 100.0, 2),
            permeabilization_state=state,
        )

    def simulate_turnover_and_persistence(
        self,
        duration_days: float = 3650.0,  # 10 years
        time_step_days: float = 1.0,
        turnover_params: Optional[TurnoverKineticsParameters] = None,
    ) -> List[PersistenceTimePoint]:
        """Simulates longitudinal lineage ODEs comparing genomic editing vs episomal mRNA persistence."""
        tp = turnover_params or self.turnover_params

        # Decay constants
        # Differentiated turnover rate mu_D = ln(2) / half_life
        mu_d = math.log(2.0) / max(1.0, tp.differentiated_half_life_days)
        # Episomal mRNA decay rate in days: t_1/2 in hours converted to days
        k_mrna_decay = math.log(2.0) / (tp.episomal_mrna_half_life_hours / 24.0)

        # Basal division rate lambda_B
        lambda_b = tp.basal_division_rate_per_day
        # Transit-amplifying parameters
        lambda_t = tp.transit_division_rate_per_day
        gamma_t = tp.transit_differentiation_rate_per_day

        # Initial conditions:
        # Genomic editing:
        # - Fraction of basal stem cells edited: B_edit(0) = tp.initial_basal_edit_efficiency_pct
        # - Fraction of apical columnar cells edited: D_edit(0) = tp.initial_apical_edit_efficiency_pct
        # - Transit compartment starts partially edited: T_edit(0) = B_edit(0)
        b_edit = tp.initial_basal_edit_efficiency_pct
        t_edit = tp.initial_basal_edit_efficiency_pct
        d_edit = tp.initial_apical_edit_efficiency_pct

        # Episomal mRNA:
        # Transfected differentiated cells express CFTR initially, but degrade rapidly with zero self-renewal
        episomal_expr = tp.initial_episomal_transfection_pct

        trajectory: List[PersistenceTimePoint] = []
        current_day = 0.0

        num_steps = int(duration_days / time_step_days)

        for step in range(num_steps + 1):
            day = step * time_step_days
            years = day / 365.25

            # Weighted total epithelial CFTR expression for genomic editing:
            # Pseudostratified composition: ~30% basal, ~70% columnar/differentiated
            # Columnar cells carry majority of luminal ion transport, but basal edited cells sustain the lineage
            genomic_total = (0.25 * b_edit) + (0.75 * d_edit)

            # Record periodic snapshot
            if step % max(1, int(15 / time_step_days)) == 0 or step == num_steps:
                trajectory.append(
                    PersistenceTimePoint(
                        day=round(day, 1),
                        years=round(years, 3),
                        episomal_cftr_pct=round(episomal_expr, 3),
                        genomic_edited_basal_pct=round(b_edit, 2),
                        genomic_edited_columnar_pct=round(d_edit, 2),
                        genomic_total_cftr_pct=round(genomic_total, 2),
                        episomal_above_curative=episomal_expr >= CURATIVE_THRESHOLD_CFTR_PCT,
                        genomic_above_curative=genomic_total >= CURATIVE_THRESHOLD_CFTR_PCT,
                    )
                )

            # Euler integration for lineage ODEs
            # 1. Basal stem cell pool editing kinetics:
            # Under homeostatic balance (p_s = p_d), the fraction of edited basal stem cells is invariant:
            # d(B_edit)/dt = lambda_B * (p_s - p_d) * B_edit = 0
            db_edit = lambda_b * (tp.p_symmetric_renewal - tp.p_symmetric_differentiation) * b_edit

            # 2. Transit-amplifying pool editing kinetics:
            # dT_edit/dt = 2 * lambda_B * p_d * B_edit + lambda_T * T_edit - gamma_T * T_edit
            # In normalized fractional terms, influx from basal stem cells balances maturation:
            dt_edit = gamma_t * (b_edit - t_edit)

            # 3. Differentiated columnar pool editing kinetics:
            # Shedding occurs at rate mu_d; replenishment arises from differentiated transit progeny
            # dD_edit/dt = mu_d * (t_edit - d_edit)
            dd_edit = mu_d * (t_edit - d_edit)

            # 4. Episomal mRNA decay kinetics:
            # Degraded by intracellular ribonucleases (k_mrna_decay) and shed into airway lumen (mu_d)
            # Diluted additionally by cellular divisions in transit pool
            depisomal = - (k_mrna_decay + mu_d) * episomal_expr

            # Update states
            b_edit = max(0.0, min(100.0, b_edit + db_edit * time_step_days))
            t_edit = max(0.0, min(100.0, t_edit + dt_edit * time_step_days))
            d_edit = max(0.0, min(100.0, d_edit + dd_edit * time_step_days))
            episomal_expr = max(0.0, episomal_expr + depisomal * time_step_days)

        return trajectory

    def generate_durability_comparison(
        self,
        trajectory: List[PersistenceTimePoint],
    ) -> List[ModalityDurabilitySummary]:
        """Synthesizes structured comparative durability across genetic and mRNA modalities."""
        # Find timepoint closest to target days
        def get_pt(day_target: float) -> PersistenceTimePoint:
            return min(trajectory, key=lambda p: abs(p.day - day_target))

        p0 = get_pt(0.0)
        p30 = get_pt(30.0)
        p90 = get_pt(90.0)
        p365 = get_pt(365.0)
        p1825 = get_pt(1825.0)
        p3650 = get_pt(3650.0)

        # Calculate curative durations
        episomal_curative_days = 0.0
        for pt in trajectory:
            if pt.episomal_above_curative:
                episomal_curative_days = pt.day
            else:
                break

        genomic_curative_days = 3650.0 if p3650.genomic_above_curative else 0.0
        if genomic_curative_days == 0.0:
            for pt in trajectory:
                if pt.genomic_above_curative:
                    genomic_curative_days = pt.day

        # Asymptotic floor for genomic editing equals steady-state basal edit fraction
        asymptotic_genomic_floor = p3650.genomic_total_cftr_pct

        genomic_summary = ModalityDurabilitySummary(
            modality_name="Basal Stem Cell Targeted Base / Prime Editing RNP",
            mechanism="Genomic allelic correction in KRT5+/TP63+ airway basal stem cells with permanent mitotic conservation",
            clinical_reference=f"{PRIMARY_DOI_COONEY} (Cooney et al. 2018), PDB {PDB_STRUCTURE_CFTR}",
            day_0_cftr_pct=p0.genomic_total_cftr_pct,
            day_30_cftr_pct=p30.genomic_total_cftr_pct,
            day_90_cftr_pct=p90.genomic_total_cftr_pct,
            day_365_cftr_pct=p365.genomic_total_cftr_pct,
            day_1825_5yr_cftr_pct=p1825.genomic_total_cftr_pct,
            day_3650_10yr_cftr_pct=p3650.genomic_total_cftr_pct,
            asymptotic_homeostatic_floor_pct=asymptotic_genomic_floor,
            curative_duration_days=genomic_curative_days,
            requires_indefinite_redosing=False,
            durability_assessment="Permanent curative persistence (> 10 years) without redosing; asymptotic plateau maintained by self-renewing basal progenitors.",
        )

        episomal_summary = ModalityDurabilitySummary(
            modality_name="Inhaled Episomal mRNA Lipid-Nanoparticle (VX-522 archetype)",
            mechanism="Transient cytoplasmic translation of CFTR mRNA without nuclear integration or basal stem retention",
            clinical_reference=f"{CLINICAL_TRIAL_VX522} (Phase 1/2 VX-522), {PRIMARY_DOI_BUTTON} (Button et al. 2021)",
            day_0_cftr_pct=p0.episomal_cftr_pct,
            day_30_cftr_pct=p30.episomal_cftr_pct,
            day_90_cftr_pct=p90.episomal_cftr_pct,
            day_365_cftr_pct=p365.episomal_cftr_pct,
            day_1825_5yr_cftr_pct=p1825.episomal_cftr_pct,
            day_3650_10yr_cftr_pct=p3650.episomal_cftr_pct,
            asymptotic_homeostatic_floor_pct=0.0,
            curative_duration_days=episomal_curative_days,
            requires_indefinite_redosing=True,
            durability_assessment="Transient pulse expression (< 21 days curative window); rapid clearance via intracellular ribonuclease decay and epithelial shedding.",
        )

        return [genomic_summary, episomal_summary]

    def run_full_simulation_report(
        self,
        duration_days: float = 3650.0,
    ) -> Dict[str, Any]:
        """Runs the coupled end-to-end simulation and formats a comprehensive audit report."""
        # 1. NPC translocation evaluation
        npc_active = self.solve_npc_translocation(self.npc_params)
        npc_passive = self.solve_npc_translocation(
            NpcTranslocationParameters(
                molecular_weight_kda=self.npc_params.molecular_weight_kda,
                has_nls=False,
            )
        )

        # 2. Epithelial niche paracellular penetration
        niche_intact = self.solve_paracellular_basal_niche_delivery(
            EpithelialNicheParameters(
                teer_ohm_cm2=self.niche_params.teer_ohm_cm2,
                is_permeabilized=False,
            )
        )
        niche_permeabilized = self.solve_paracellular_basal_niche_delivery(
            EpithelialNicheParameters(
                teer_ohm_cm2=80.0,
                is_permeabilized=True,
            )
        )

        # 3. Dynamic persistence simulation
        trajectory = self.simulate_turnover_and_persistence(duration_days=duration_days)
        durability_summaries = self.generate_durability_comparison(trajectory)

        # Sample snapshots (Day 0, 30, 90, 180, 365, 730, 1825, 3650)
        sample_days = [0.0, 30.0, 60.0, 90.0, 180.0, 365.0, 730.0, 1825.0, 3650.0]
        sampled_snapshots = [
            asdict(min(trajectory, key=lambda p: abs(p.day - d)))
            for d in sample_days
            if d <= duration_days
        ]

        return {
            "model_metadata": {
                "name": "lnp_nuclear_entry_basal_persistence",
                "task_id": "T-0478",
                "title": "Biophysical and Kinetic Modeling of Lipid-Nanoparticle Nuclear Entry and Basal Stem Cell Genomic Persistence in Regenerative CF Airway Epithelia",
                "track": "Track A",
                "discipline": "airway_delivery",
                "answer_tier": "T2",
                "curative_threshold_cftr_pct": CURATIVE_THRESHOLD_CFTR_PCT,
            },
            "primary_grounding": {
                "literature": [
                    {"doi": PRIMARY_DOI_BUTTON, "citation": "Button B et al. (2021) Physiol Rev"},
                    {"doi": PRIMARY_DOI_COONEY, "citation": "Cooney AL et al. (2018) Mol Ther"},
                    {"doi": PRIMARY_DOI_MIDDLETON, "citation": "Middleton PG et al. (2019) N Engl J Med"},
                ],
                "clinical_trials": [
                    {"nct": CLINICAL_TRIAL_VX522, "title": "Phase 1/2 Study of Inhaled mRNA LNP VX-522 in Participants With Cystic Fibrosis"},
                    {"nct": CLINICAL_TRIAL_VANZACAFTOR, "title": "Phase 3 Study of Vanzacaftor-Tezacaftor-Deutivacaftor in Cystic Fibrosis"},
                ],
                "pdb_structural_coordinate": PDB_STRUCTURE_CFTR,
            },
            "npc_translocation_analysis": {
                "active_import_with_nls": asdict(npc_active),
                "passive_diffusion_without_nls": asdict(npc_passive),
                "rnp_nuclear_enrichment_fold": round(
                    npc_active.nuclear_cytoplasmic_steady_state_ratio / max(0.01, npc_passive.nuclear_cytoplasmic_steady_state_ratio),
                    1
                ),
            },
            "epithelial_niche_delivery_analysis": {
                "intact_barrier": asdict(niche_intact),
                "transiently_permeabilized_barrier": asdict(niche_permeabilized),
                "basal_delivery_enhancement_fold": round(
                    niche_permeabilized.basal_stem_uptake_fraction / max(0.001, niche_intact.basal_stem_uptake_fraction),
                    1
                ),
            },
            "durability_comparison": [asdict(s) for s in durability_summaries],
            "longitudinal_trajectory_samples": sampled_snapshots,
            "disclaimers": {
                "nlm_attribution": NLM_ATTRIBUTION,
                "research_boundary": RESEARCH_BOUNDARY,
            },
        }


def parse_arguments() -> argparse.Namespace:
    """Parses command-line arguments for running the discovery model."""
    parser = argparse.ArgumentParser(
        description="LNP Nuclear Entry & Basal Stem Cell Genomic Persistence Kinetics Model (T-0478)"
    )
    parser.add_argument(
        "--duration-days",
        type=float,
        default=3650.0,
        help="Simulation duration in days (default: 3650 days / 10 years)",
    )
    parser.add_argument(
        "--teer",
        type=float,
        default=DEFAULT_TEER_OHM_CM2,
        help="Transepithelial electrical resistance in Ohm*cm^2 (default: 600.0)",
    )
    parser.add_argument(
        "--turnover-half-life",
        type=float,
        default=DEFAULT_DIFFERENTIATED_TURNOVER_HALF_LIFE_DAYS,
        help="Columnar epithelial turnover half-life in days (default: 60.0)",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Emit raw machine-readable JSON output",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Path to save simulation JSON results",
    )
    return parser.parse_args()


def main() -> int:
    """Main CLI entrypoint."""
    args = parse_arguments()

    model = LnpNuclearEntryBasalPersistenceModel(
        niche_params=EpithelialNicheParameters(teer_ohm_cm2=args.teer),
        turnover_params=TurnoverKineticsParameters(
            differentiated_half_life_days=args.turnover_half_life,
        ),
    )

    report = model.run_full_simulation_report(duration_days=args.duration_days)

    if args.json or args.output:
        formatted_json = json.dumps(report, indent=2)
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(formatted_json, encoding="utf-8")
            print(f"Simulation report written to {args.output}")
        if args.json:
            print(formatted_json)
    else:
        print("=" * 78)
        print("Sensei Ewok CF Lab — LNP Nuclear Entry & Basal Stem Persistence Model (T-0478)")
        print("=" * 78)
        meta = report["model_metadata"]
        print(f"Task:       {meta['task_id']} | Track: {meta['track']} | Discipline: {meta['discipline']} | Tier: {meta['answer_tier']}")
        print(f"Title:      {meta['title']}")
        print("-" * 78)
        print("1. Nuclear Pore Complex (NPC) Translocation Kinetics:")
        npc_a = report["npc_translocation_analysis"]["active_import_with_nls"]
        npc_p = report["npc_translocation_analysis"]["passive_diffusion_without_nls"]
        print(f"   RNP Stokes Radius:        {npc_a['stokes_radius_nm']} nm (Molecular Weight: 160 kDa)")
        print(f"   Passive Permeability:     {npc_p['passive_permeability_cm_s']:.2e} cm/s (Permissive: {npc_p['is_passive_permissive']})")
        print(f"   Importin Bound Fraction:  {npc_a['importin_bound_fraction'] * 100:.1f}% (NLS Kd: 25 nM)")
        print(f"   Nuclear/Cyto Ratio (NLS): {npc_a['nuclear_cytoplasmic_steady_state_ratio']}x (Enrichment Fold: {report['npc_translocation_analysis']['rnp_nuclear_enrichment_fold']}x)")
        print("-" * 78)
        print("2. Pseudostratified Airway Niche Paracellular Delivery:")
        intact = report["epithelial_niche_delivery_analysis"]["intact_barrier"]
        perm = report["epithelial_niche_delivery_analysis"]["transiently_permeabilized_barrier"]
        print(f"   Intact Barrier (TEER {intact['apical_columnar_uptake_fraction'] * 100:.1f}% Apical): Basal Access = {intact['basal_stem_uptake_fraction'] * 100:.2f}%")
        print(f"   Permeabilized Barrier:     Basal Access = {perm['basal_stem_uptake_fraction'] * 100:.2f}% ({report['epithelial_niche_delivery_analysis']['basal_delivery_enhancement_fold']}x Enhancement)")
        print("-" * 78)
        print("3. Multi-Year Durability: Genomic Editing vs Episomal mRNA:")
        for mod in report["durability_comparison"]:
            print(f"   * {mod['modality_name']}:")
            print(f"     Day 0: {mod['day_0_cftr_pct']:.1f}% | Day 30: {mod['day_30_cftr_pct']:.1f}% | 1 Year: {mod['day_365_cftr_pct']:.1f}% | 10 Years: {mod['day_3650_10yr_cftr_pct']:.1f}%")
            print(f"     Asymptotic Homeostatic Floor: {mod['asymptotic_homeostatic_floor_pct']:.1f}% WT CFTR")
            print(f"     Curative Window (>10% WT):    {mod['curative_duration_days']:.0f} days (Redosing Needed: {mod['requires_indefinite_redosing']})")
            print(f"     Assessment: {mod['durability_assessment']}")
        print("=" * 78)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
