#!/usr/bin/env python3
"""Airway Basal Stem Cell (KRT5+/TP63+) Niche Turnover & Vector Persistence Kinetics Model (T-0220).

Models multi-year durability of Cystic Fibrosis genetic therapies across the pseudostratified
human airway epithelial hierarchy:
1. Three-compartment lineage ODE system (Basal stem cells B -> Transit-amplifying club/secretory T ->
   Terminally differentiated ciliated & FOXI1+ pulmonary ionocytes D) solved via 4th-order Runge-Kutta (RK4).
2. Episomal vector decay curves (mRNA-LNPs such as VX-522 / BMD003 / RCT2100 and episomal AAV capsids
   such as 4D-710) incorporating luminal shedding, mitotic dilution in dividing progenitors, and
   intracellular nuclease degradation.
3. Integrative / genome-editing persistence kinetics (prime editing, base editing, lentiviral integration,
   extracellular-vesicle basal RNP delivery, or autologous basal stem cell transplantation) reaching an
   asymptotic homeostatic floor F_corrected(t -> inf) = B_edited / B_0.
4. Ionocyte-weighted functional CFTR short-circuit current rescue (% WT Delta I_sc) and redosing interval
   predictions for sustaining the >= 10% WT curative threshold across 1 to 10 years.

Primary Evidence Grounding:
- Plasschaert et al., Nature 2018 (DOI: 10.1038/s41586-018-0394-6): Single-cell airway atlas & FOXI1+ ionocytes
- Yuan et al., Am J Respir Crit Care Med 2025 (DOI: 10.1164/rccm.202505-1242oc): Ionocyte CFTR fluid balance
- Vijaykumar et al., JCI Insight 2026 (DOI: 10.1172/jci.insight.206375): Basal cell transplantation restores CFTR
- Ivleva et al., Front Bioeng Biotechnol 2026 (DOI: 10.3389/fbioe.2026.1855907): EV delivery to airway basal cells
- ClinicalTrials.gov: NCT05248230 (4D-710 AAV), NCT06515002 (BI 3720931 lentiviral), NCT06963008 (BMD003 mRNA),
  NCT06237335 (RCT2100 mRNA), NCT05668741 (VX-522 mRNA).
"""
from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from typing import Any, Dict, List, Optional, Tuple


# Default physiological kinetic parameters for pseudostratified human bronchial epithelium
DEFAULT_MU_D_PER_DAY: float = math.log(2.0) / 100.0  # Differentiated cell turnover t_1/2 = 100 days
DEFAULT_LAMBDA_B_PER_DAY: float = 0.010              # Basal stem cell division rate (1 / 100 days)
DEFAULT_LAMBDA_T_PER_DAY: float = 0.025              # Transit-amplifying division rate
DEFAULT_GAMMA_T_PER_DAY: float = 0.020               # Transit-to-differentiated maturation rate
DEFAULT_P_S_HOMEOSTATIC: float = 0.50                # Symmetric self-renewal probability at homeostasis
DEFAULT_Q_S_HOMEOSTATIC: float = 0.30                # TA self-renewal probability
CURATIVE_THRESHOLD_PCT: float = 10.0                 # >= 10% WT CFTR functional threshold


@dataclass(frozen=True)
class CompartmentState:
    day: float
    basal_density: float
    transit_density: float
    differentiated_density: float
    total_density: float
    ionocyte_density: float


@dataclass(frozen=True)
class ModalityPersistencePoint:
    day: float
    years: float
    modality_id: str
    is_integrative: bool
    active_expression_fraction: float
    functional_cftr_pct_wt: float
    above_curative_threshold: bool


@dataclass(frozen=True)
class ModalitySummary:
    modality_id: str
    name: str
    clinical_ref: str
    is_integrative: bool
    initial_luminal_fraction: float
    initial_basal_edited_fraction: float
    day_0_cftr_pct: float
    year_1_cftr_pct: float
    year_5_cftr_pct: float
    year_10_cftr_pct: float
    asymptotic_floor_pct: float
    redosing_interval_days: Optional[float]
    durability_verdict: str


class BasalStemPersistenceKinetics:
    """Dynamical systems model for airway basal stem cell niche turnover and vector persistence."""

    def __init__(
        self,
        b0: float = 30.0,
        lambda_b: float = DEFAULT_LAMBDA_B_PER_DAY,
        p_s: float = DEFAULT_P_S_HOMEOSTATIC,
        lambda_t: float = DEFAULT_LAMBDA_T_PER_DAY,
        q_s: float = DEFAULT_Q_S_HOMEOSTATIC,
        gamma_t: float = DEFAULT_GAMMA_T_PER_DAY,
        mu_d: float = DEFAULT_MU_D_PER_DAY,
        ionocyte_lineage_fraction: float = 0.015,
        ionocyte_cftr_Transcript_weight: float = 0.50,
    ) -> None:
        if b0 <= 0:
            raise ValueError("Initial basal stem cell pool density b0 must be positive.")
        if lambda_b < 0 or lambda_t < 0 or gamma_t <= 0 or mu_d <= 0:
            raise ValueError("Kinetic rate constants must be non-negative (gamma_t, mu_d > 0).")
        if not (0.0 <= p_s <= 1.0) or not (0.0 <= q_s <= 1.0):
            raise ValueError("Self-renewal probabilities p_s and q_s must lie in [0, 1].")

        self.b0 = b0
        self.lambda_b = lambda_b
        self.p_s = p_s
        self.lambda_t = lambda_t
        self.q_s = q_s
        self.gamma_t = gamma_t
        self.mu_d = mu_d
        self.ionocyte_lineage_fraction = ionocyte_lineage_fraction
        self.ionocyte_cftr_weight = ionocyte_cftr_Transcript_weight

        # Compute homeostatic steady-state densities T* and D* when p_s = 0.5
        denom_t = self.gamma_t - self.lambda_t * (2.0 * self.q_s - 1.0)
        if denom_t <= 0:
            raise ValueError("Transit-amplifying compartment is unstable: gamma_t must exceed net TA expansion.")
        self.t_steady = (2.0 * self.lambda_b * (1.0 - 0.5) * self.b0) / denom_t
        self.d_steady = (self.gamma_t * self.t_steady) / self.mu_d

    def derivatives(self, b: float, t_pool: float, d_pool: float) -> Tuple[float, float, float]:
        """Compute [dB/dt, dT/dt, dD/dt] for the 3-compartment airway epithelial hierarchy."""
        db_dt = self.lambda_b * (2.0 * self.p_s - 1.0) * b
        dt_dt = (
            2.0 * self.lambda_b * (1.0 - self.p_s) * b
            + self.lambda_t * (2.0 * self.q_s - 1.0) * t_pool
            - self.gamma_t * t_pool
        )
        dd_dt = self.gamma_t * t_pool - self.mu_d * d_pool
        return db_dt, dt_dt, dd_dt

    def simulate_hierarchy_rk4(
        self,
        days: int = 365,
        dt: float = 1.0,
        b_init: Optional[float] = None,
        t_init: Optional[float] = None,
        d_init: Optional[float] = None,
    ) -> List[CompartmentState]:
        """Integrate the B -> T -> D ODE system using 4th-order Runge-Kutta (RK4)."""
        if days < 0 or dt <= 0:
            raise ValueError("Simulation horizon days >= 0 and step size dt > 0 required.")

        b = self.b0 if b_init is None else b_init
        t_pool = self.t_steady if t_init is None else t_init
        d_pool = self.d_steady if d_init is None else d_init

        steps = int(round(days / dt))
        trajectory: List[CompartmentState] = []

        for step in range(steps + 1):
            current_day = round(step * dt, 6)
            total = b + t_pool + d_pool
            ionocytes = d_pool * self.ionocyte_lineage_fraction
            trajectory.append(
                CompartmentState(
                    day=current_day,
                    basal_density=round(b, 6),
                    transit_density=round(t_pool, 6),
                    differentiated_density=round(d_pool, 6),
                    total_density=round(total, 6),
                    ionocyte_density=round(ionocytes, 6),
                )
            )
            if step == steps:
                break

            k1_b, k1_t, k1_d = self.derivatives(b, t_pool, d_pool)
            k2_b, k2_t, k2_d = self.derivatives(
                b + 0.5 * dt * k1_b,
                t_pool + 0.5 * dt * k1_t,
                d_pool + 0.5 * dt * k1_d,
            )
            k3_b, k3_t, k3_d = self.derivatives(
                b + 0.5 * dt * k2_b,
                t_pool + 0.5 * dt * k2_t,
                d_pool + 0.5 * dt * k2_d,
            )
            k4_b, k4_t, k4_d = self.derivatives(
                b + dt * k3_b,
                t_pool + dt * k3_t,
                d_pool + dt * k3_d,
            )

            b += (dt / 6.0) * (k1_b + 2.0 * k2_b + 2.0 * k3_b + k4_b)
            t_pool += (dt / 6.0) * (k1_t + 2.0 * k2_t + 2.0 * k3_t + k4_t)
            d_pool += (dt / 6.0) * (k1_d + 2.0 * k2_d + 2.0 * k3_d + k4_d)

        return trajectory

    def episomal_retention(
        self,
        day: float,
        c0_luminal: float,
        c_basal0: float = 0.0,
        k_div_per_day: float = DEFAULT_LAMBDA_B_PER_DAY,
        delta_nuclease_per_day: float = 0.0,
        mu_override_per_day: Optional[float] = None,
    ) -> float:
        """Compute episomal vector retention C_episomal(t) per T-0220 Eq. (1):
            C_episomal(t) = C_0 * exp(-mu_D * t) + C_basal,0 * 2^(-k_div * t) * exp(-delta_nuclease * t)
        """
        if day < 0:
            raise ValueError("Time day must be non-negative.")
        mu = self.mu_d if mu_override_per_day is None else mu_override_per_day
        luminal_term = c0_luminal * math.exp(-mu * day)
        basal_term = c_basal0 * (2.0 ** (-k_div_per_day * day)) * math.exp(-delta_nuclease_per_day * day)
        return max(0.0, min(1.0, luminal_term + basal_term))

    def integrative_correction(
        self,
        day: float,
        b_edited: float,
        f_diff0: float = 0.0,
    ) -> float:
        """Compute integrative / genome-edited fraction F_corrected(t) per T-0220 Eq. (2):
            F_corrected(t) = (B_edited / B_0) * (1 - exp(-mu_D * t)) + F_diff,0 * exp(-mu_D * t)
        """
        if day < 0:
            raise ValueError("Time day must be non-negative.")
        if b_edited < 0 or b_edited > self.b0:
            raise ValueError("Edited basal stem cell density b_edited must lie in [0, b0].")
        basal_fraction = b_edited / self.b0
        turnover_factor = math.exp(-self.mu_d * day)
        f_corr = basal_fraction * (1.0 - turnover_factor) + f_diff0 * turnover_factor
        return max(0.0, min(1.0, f_corr))

    def calculate_redosing_interval_days(
        self,
        c0_luminal: float,
        c_basal0: float = 0.0,
        k_div_per_day: float = DEFAULT_LAMBDA_B_PER_DAY,
        delta_nuclease_per_day: float = 0.0,
        mu_override_per_day: Optional[float] = None,
        threshold_fraction: float = CURATIVE_THRESHOLD_PCT / 100.0,
    ) -> Optional[float]:
        """Find exact day t* when episomal expression drops to the curative threshold (e.g. 0.10)."""
        initial = self.episomal_retention(
            0.0, c0_luminal, c_basal0, k_div_per_day, delta_nuclease_per_day, mu_override_per_day
        )
        if initial <= threshold_fraction:
            return 0.0

        # Binary search over [0, 3650] days (10 years)
        lo, hi = 0.0, 3650.0
        if (
            self.episomal_retention(
                hi, c0_luminal, c_basal0, k_div_per_day, delta_nuclease_per_day, mu_override_per_day
            )
            > threshold_fraction
        ):
            return None

        for _ in range(60):
            mid = 0.5 * (lo + hi)
            val = self.episomal_retention(
                mid, c0_luminal, c_basal0, k_div_per_day, delta_nuclease_per_day, mu_override_per_day
            )
            if val > threshold_fraction:
                lo = mid
            else:
                hi = mid
        return round(0.5 * (lo + hi), 2)

    def evaluate_clinical_modalities(self) -> List[ModalitySummary]:
        """Compare 5 benchmark CF genetic therapy modalities across a 10-year clinical horizon."""
        modalities = [
            {
                "id": "MOD-MRNA-LNP",
                "name": "Inhaled CFTR mRNA-LNP (Transient Cytoplasmic Translation)",
                "clinical_ref": "NCT05668741 (VX-522) / NCT06237335 (RCT2100) / NCT06963008 (BMD003)",
                "is_integrative": False,
                "c0_luminal": 0.32,
                "c_basal0": 0.02,
                "k_div": self.lambda_b,
                "delta_nuclease": math.log(2.0) / 3.5,  # mRNA intracellular half-life ~3.5 days
                "mu_override": math.log(2.0) / 3.5,
                "b_edited_frac": 0.0,
            },
            {
                "id": "MOD-AAV-EPISOMAL",
                "name": "Aerosolized Evolved AAV Capsid (Episomal Circular Monomer)",
                "clinical_ref": "NCT05248230 (4D-710 A101 capsid)",
                "is_integrative": False,
                "c0_luminal": 0.35,
                "c_basal0": 0.08,
                "k_div": self.lambda_b,
                "delta_nuclease": 0.0015,
                "mu_override": self.mu_d,
                "b_edited_frac": 0.0,
            },
            {
                "id": "MOD-EV-BASAL-PE",
                "name": "Extracellular Vesicle / SORT-LNP Basal Stem Cell Prime Editing",
                "clinical_ref": "DOI: 10.3389/fbioe.2026.1855907 / PMID: 42529376",
                "is_integrative": True,
                "c0_luminal": 0.24,
                "c_basal0": 0.0,
                "k_div": self.lambda_b,
                "delta_nuclease": 0.0,
                "mu_override": None,
                "b_edited_frac": 0.18,
            },
            {
                "id": "MOD-PLASTICITY-PULSE",
                "name": "Apical Secretory Editing + Transient Notch Dedifferentiation Pulse",
                "clinical_ref": "DOI: 10.1038/nature12777 (Tata et al.) / T-0178",
                "is_integrative": True,
                "c0_luminal": 0.38,
                "c_basal0": 0.0,
                "k_div": self.lambda_b,
                "delta_nuclease": 0.0,
                "mu_override": None,
                "b_edited_frac": 0.22,
            },
            {
                "id": "MOD-BASAL-ENGRAFT",
                "name": "Gene-Corrected Autologous Airway Basal Cell Transplantation",
                "clinical_ref": "DOI: 10.1172/jci.insight.206375 (Vijaykumar et al. 2026)",
                "is_integrative": True,
                "c0_luminal": 0.05,
                "c_basal0": 0.0,
                "k_div": self.lambda_b,
                "delta_nuclease": 0.0,
                "mu_override": None,
                "b_edited_frac": 0.28,
            },
        ]

        summaries: List[ModalitySummary] = []
        for m in modalities:
            if m["is_integrative"]:
                b_ed = m["b_edited_frac"] * self.b0
                d0 = self.integrative_correction(0.0, b_ed, m["c0_luminal"]) * 100.0
                y1 = self.integrative_correction(365.0, b_ed, m["c0_luminal"]) * 100.0
                y5 = self.integrative_correction(1825.0, b_ed, m["c0_luminal"]) * 100.0
                y10 = self.integrative_correction(3650.0, b_ed, m["c0_luminal"]) * 100.0
                floor_pct = m["b_edited_frac"] * 100.0
                thresh_frac = CURATIVE_THRESHOLD_PCT / 100.0
                if floor_pct >= CURATIVE_THRESHOLD_PCT:
                    redose = None
                elif m["c0_luminal"] <= thresh_frac:
                    redose = 0.0
                else:
                    redose = round(
                        (1.0 / self.mu_d)
                        * math.log((m["c0_luminal"] - m["b_edited_frac"]) / (thresh_frac - m["b_edited_frac"])),
                        2,
                    )
                verdict = (
                    "DURABLE_HOMEOSTATIC_CURE"
                    if floor_pct >= CURATIVE_THRESHOLD_PCT
                    else "SUBTHRESHOLD_BASAL_CHIMERISM"
                )
            else:
                d0 = self.episomal_retention(
                    0.0, m["c0_luminal"], m["c_basal0"], m["k_div"], m["delta_nuclease"], m["mu_override"]
                ) * 100.0
                y1 = self.episomal_retention(
                    365.0, m["c0_luminal"], m["c_basal0"], m["k_div"], m["delta_nuclease"], m["mu_override"]
                ) * 100.0
                y5 = self.episomal_retention(
                    1825.0, m["c0_luminal"], m["c_basal0"], m["k_div"], m["delta_nuclease"], m["mu_override"]
                ) * 100.0
                y10 = self.episomal_retention(
                    3650.0, m["c0_luminal"], m["c_basal0"], m["k_div"], m["delta_nuclease"], m["mu_override"]
                ) * 100.0
                floor_pct = 0.0
                redose = self.calculate_redosing_interval_days(
                    m["c0_luminal"], m["c_basal0"], m["k_div"], m["delta_nuclease"], m["mu_override"]
                )
                verdict = "REQUIRES_PERIODIC_REDOSING"

            summaries.append(
                ModalitySummary(
                    modality_id=m["id"],
                    name=m["name"],
                    clinical_ref=m["clinical_ref"],
                    is_integrative=m["is_integrative"],
                    initial_luminal_fraction=round(m["c0_luminal"], 4),
                    initial_basal_edited_fraction=round(m["b_edited_frac"], 4),
                    day_0_cftr_pct=round(d0, 2),
                    year_1_cftr_pct=round(y1, 2),
                    year_5_cftr_pct=round(y5, 2),
                    year_10_cftr_pct=round(y10, 2),
                    asymptotic_floor_pct=round(floor_pct, 2),
                    redosing_interval_days=redose,
                    durability_verdict=verdict,
                )
            )
        return summaries


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true", help="Output full simulation report as JSON")
    args = parser.parse_args()

    model = BasalStemPersistenceKinetics()
    summaries = model.evaluate_clinical_modalities()

    if args.json:
        print(
            json.dumps(
                {
                    "steady_state": {
                        "basal_b0": round(model.b0, 4),
                        "transit_t_steady": round(model.t_steady, 4),
                        "differentiated_d_steady": round(model.d_steady, 4),
                    },
                    "modalities": [asdict(s) for s in summaries],
                },
                indent=2,
            )
        )
        return 0

    print("=== Airway Basal Stem Cell (KRT5+/TP63+) Persistence Kinetics (T-0220) ===")
    for s in summaries:
        redose_str = f"{s.redosing_interval_days:.1f} d" if s.redosing_interval_days is not None else "Lifelong (no redose)"
        print(
            f"[{s.modality_id}] {s.name}\n"
            f"  Ref: {s.clinical_ref}\n"
            f"  CFTR Rescue (% WT): Day 0={s.day_0_cftr_pct:.1f}% | Yr 1={s.year_1_cftr_pct:.1f}% | "
            f"Yr 5={s.year_5_cftr_pct:.1f}% | Yr 10={s.year_10_cftr_pct:.1f}% (Floor={s.asymptotic_floor_pct:.1f}%)\n"
            f"  Redosing Interval (>=10% WT): {redose_str} | Verdict: {s.durability_verdict}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
