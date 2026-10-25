#!/usr/bin/env python3
"""Airway Basal Stem Cell (KRT5+/KRT14+) Paracellular Tight Junction Permeation & Clonal Persistence Model (T-0429).

Models paracellular transport of aerosolized vectors across tight junctions (ZO-1, claudins),
Renkin steric hindrance partitioning, and the long-term clonal durability of KRT5+/KRT14+
basal stem cells versus shedding ciliated cells.

Primary Evidence:
  - JCI Insight (2026), DOI: 10.1172/jci.insight.206375, PMID: 42611556
  - Front Bioeng Biotechnol (2026), DOI: 10.3389/fbioe.2026.1855907, PMID: 42529376
  - J Virol (2000), DOI: 10.1128/jvi.74.13.6050-6057.2000
  - J Virol (2004), DOI: 10.1128/JVI.78.24.13755-13768.2004
  - Nat Biotechnol (2021), DOI: 10.1038/s41587-021-00943-7
  - Nat Biomed Eng (2021), DOI: 10.1038/s41551-021-00724-5
Clinical Trials:
  - NCT05248230 (4D-710 aerosolized AAV in CF)
  - NCT05668741 (VX-522 aerosolized mRNA LNP in CF)
Structures:
  - RCSB PDB 6MSM, 5UAK
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Tuple

# Mandatory ClinicalTrials.gov Attribution
NLM_ATTRIBUTION = (
    "ClinicalTrials.gov does not evaluate safety or efficacy; "
    "listing on the registry does not denote government endorsement."
)

CLINICAL_TRIAL_4D710 = "NCT05248230"
CLINICAL_TRIAL_VX522 = "NCT05668741"
PDB_DEPHOSPHO_CFTR = "6MSM"
PDB_ACTIVE_CFTR = "5UAK"

PRIMARY_DOI_FERRET_BASAL = "10.1172/jci.insight.206375"
PRIMARY_DOI_EV_BASAL = "10.3389/fbioe.2026.1855907"
PRIMARY_DOI_CAR_RECEPTOR = "10.1128/jvi.74.13.6050-6057.2000"
PRIMARY_DOI_GLYCOCALYX = "10.1128/JVI.78.24.13755-13768.2004"
PRIMARY_DOI_BASE_EDITING = "10.1038/s41587-021-00943-7"
PRIMARY_DOI_MRNA_LNP = "10.1038/s41551-021-00724-5"


def calculate_renkin_partitioning(lambda_ratio: float) -> float:
    """Calculate the Renkin steric hindrance partitioning coefficient for cylindrical pores.

    Args:
        lambda_ratio: Ratio of particle hydrodynamic radius to pore radius (r_h / r_pore).

    Returns:
        Fractional partitioning factor F_Renkin in [0.0, 1.0].
    """
    if lambda_ratio <= 0.0:
        return 1.0
    if lambda_ratio >= 1.0:
        return 0.0
    lam = lambda_ratio
    return (1.0 - lam) ** 2 * (
        1.0 - 2.104 * lam + 2.09 * (lam ** 3) - 0.95 * (lam ** 5)
    )


@dataclass(frozen=True)
class TightJunctionKinetics:
    """Dynamic parameters for tight junction opening and paracellular vector transport."""
    r0_nm: float = 0.5            # Resting pore radius (nm)
    r_max_nm: float = 35.0        # Maximum dilated pore radius under transient permeabilization (nm)
    teer0_ohm_cm2: float = 500.0  # Baseline TEER (Ohm*cm^2)
    tau_reseal_hours: float = 2.0 # Spontaneous tight junction resealing half-life (hours)

    def teer_at_radius(self, r_pore_nm: float) -> float:
        """Calculate Transepithelial Electrical Resistance from effective pore radius."""
        if r_pore_nm <= self.r0_nm:
            return self.teer0_ohm_cm2
        return self.teer0_ohm_cm2 * (self.r0_nm / r_pore_nm) ** 2

    def dynamic_pore_radius(self, t_hours: float, duration_hours: float) -> float:
        """Compute effective pore radius at time t following transient opening."""
        if t_hours < 0.0:
            return self.r0_nm
        if t_hours <= duration_hours:
            # Saturated dilated phase
            return self.r_max_nm
        # Exponential resealing phase
        dt = t_hours - duration_hours
        decay = math.exp(-dt / self.tau_reseal_hours)
        return self.r0_nm + (self.r_max_nm - self.r0_nm) * decay

    def vector_paracellular_permeability(
        self,
        r_pore_nm: float,
        vector_radius_nm: float = 12.5,
    ) -> float:
        """Compute relative paracellular permeability factor in [0.0, 1.0]."""
        if r_pore_nm <= 0.0:
            return 0.0
        lam = vector_radius_nm / r_pore_nm
        return calculate_renkin_partitioning(lam)


@dataclass
class AirwayEpithelialLineageModel:
    """Epithelial hierarchy and clonal persistence ODE solver."""
    total_cells: float = 1000.0
    basal_fraction: float = 0.30       # 30% basal stem cells
    ta_fraction: float = 0.20          # 20% transit-amplifying cells
    ciliated_fraction: float = 0.50    # 50% differentiated ciliated/secretory cells
    ciliated_halflife_days: float = 60.0 # Shedding half-life (days)
    stem_division_rate_day: float = 0.012# Basal cell division rate (1/day)
    symmetric_self_renewal_prob: float = 0.50 # Neutral homeostatic self-renewal fraction (neutral drift)

    @property
    def k_shed(self) -> float:
        """Ciliated cell homeostatic shedding rate (1/day)."""
        return math.log(2.0) / self.ciliated_halflife_days

    @property
    def k_diff(self) -> float:
        """Transit amplifying differentiation rate into mature ciliated cells (1/day)."""
        return 0.05

    def simulate_persistence(
        self,
        days: float = 365.0,
        dt_days: float = 1.0,
        initial_transduced_basal: float = 50.0,
        initial_transduced_ciliated: float = 250.0,
    ) -> List[Dict[str, float]]:
        """Simulate longitudinal CFTR+ cell numbers across the epithelial hierarchy."""
        n_steps = int(days / dt_days) + 1
        history: List[Dict[str, float]] = []

        b_star = initial_transduced_basal
        ta_star = 0.0
        d_star = initial_transduced_ciliated

        for step in range(n_steps):
            t = step * dt_days
            total_transduced = b_star + ta_star + d_star
            ciliated_pct = (d_star / (self.total_cells * self.ciliated_fraction)) * 100.0

            history.append({
                "t_days": t,
                "transduced_basal": b_star,
                "transduced_ta": ta_star,
                "transduced_ciliated": d_star,
                "total_transduced": total_transduced,
                "ciliated_rescue_pct": ciliated_pct,
            })

            # ODE Derivatives
            # Basal stem cell self-renewal vs differentiation
            db_dt = self.stem_division_rate_day * (2.0 * self.symmetric_self_renewal_prob - 1.0) * b_star
            dta_dt = 2.0 * self.stem_division_rate_day * (1.0 - self.symmetric_self_renewal_prob) * b_star - self.k_diff * ta_star
            dd_dt = self.k_diff * ta_star - self.k_shed * d_star

            # Euler integration step
            b_star = max(0.0, b_star + db_dt * dt_days)
            ta_star = max(0.0, ta_star + dta_dt * dt_days)
            d_star = max(0.0, d_star + dd_dt * dt_days)

        return history
