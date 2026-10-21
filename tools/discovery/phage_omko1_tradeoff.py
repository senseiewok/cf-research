#!/usr/bin/env python3
"""Pseudomonas aeruginosa Phage OMKO1 Evolutionary Trade-off & MexAB-OprM Efflux Resensitization Dynamics (T-0246).

Biophysical modeling of:
  1. Lytic receptor engagement of bacteriophage OMKO1 attaching to the outer membrane
     porin OprM of the tripartite MexAB-OprM multidrug efflux system.
  2. Evolutionary divergence between:
       - B_S: Phage-sensitive, multidrug-resistant (OprM+, MIC_cipro approx 32.0 ug/mL)
       - B_R: Phage-resistant, antibiotic-resensitized (OprM-, MIC_cipro approx 1.0 ug/mL)
  3. Tripartite population dynamic differential equations with delayed burst kinetics:
       dB_S/dt = r_S * B_S * (1 - (B_S + B_R)/K) - k_ads * B_S * P - mu_A(C_A, MIC_S)*B_S - mu_mut*B_S + mu_rev*B_R
       dB_R/dt = r_R * B_R * (1 - (B_S + B_R)/K) - mu_A(C_A, MIC_R)*B_R + mu_mut*B_S - mu_rev*B_R
       dP/dt   = beta * k_ads * B_S(t - tau) * P(t - tau) * exp(-delta_B * tau) - k_ads * B_S * P - delta_P * P
     where r_R = r_S * (1 - c_fit) accounts for the outer membrane fitness penalty of OprM loss.
  4. Pharmacodynamic Hill antimicrobial kill rate:
       mu_A(C_A, MIC) = E_max * (C_A^H_A) / (MIC^H_A + C_A^H_A)
  5. The evolutionary trap & Phage-Antibiotic Synergy (PAS):
       - Phage OMKO1 alone lyses B_S and drives selection for B_R.
       - Ciprofloxacin alone (e.g. 10 ug/mL) fails to eradicate B_S (MIC 32 ug/mL).
       - Dual administration creates an inescapable evolutionary vice, driving total
         bacterial density B_tot -> 0 without treatment escape.

Primary Evidence Grounding:
  - Chan et al. (2018), Sci Rep, DOI: 10.1038/s41598-018-24541-1
    (Phage selection restores antibiotic sensitivity in MDR Pseudomonas aeruginosa).
  - Chan et al. (2021), J Cyst Fibros, DOI: 10.1016/j.jcf.2021.03.013
    (Phage therapy for refractory Pseudomonas aeruginosa in cystic fibrosis).
  - Clinical trial surveillance:
    * NCT04684641 (CYPHY Yale Cystic Fibrosis Bacteriophage Study).
    * NCT05010577 (Nebulized Bacteriophage Therapy in CF with Chronic P. aeruginosa).
    Both verified against tools/trial_atlas/snapshot.json.

Attribution & Disclaimers:
ClinicalTrials.gov data are used here with attribution to the U.S. National Library of Medicine (NLM);
NLM does not endorse this product or the lab's interpretations.
ClinicalTrials.gov does not evaluate safety or efficacy; listing on the registry does not denote government endorsement.
Research, not medical advice.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

NLM_REGISTRY_DISCLAIMER: str = (
    "ClinicalTrials.gov data are used here with attribution to the U.S. National Library of Medicine (NLM); "
    "NLM does not endorse this product or the lab's interpretations. "
    "ClinicalTrials.gov does not evaluate safety or efficacy; listing on the registry does not denote government endorsement."
)

RESEARCH_BOUNDARY: str = (
    "Preclinical biophysical modeling tool for research purposes only; not medical advice. "
    "No therapeutic safety or clinical efficacy claims are made."
)

# Nominal Biophysical Parameters (Grounded in Chan et al. 2018, 2021)
DEFAULT_R_S_PER_MIN: float = 0.018             # Wild-type growth rate (~38.5 min doubling time in sputum)
DEFAULT_C_FIT: float = 0.15                    # OprM loss outer membrane fitness cost (0.08 - 0.25)
DEFAULT_K_CFU_ML: float = 2.5e9                # Biofilm carrying capacity (1e8 - 1e10 CFU/mL)
DEFAULT_K_ADS_ML_PER_CFU_MIN: float = 1.8e-9   # OMKO1 adsorption affinity for OprM
DEFAULT_BURST_SIZE: float = 85.0               # Lytic burst size (virions / lysed cell)
DEFAULT_LATENT_PERIOD_MIN: float = 32.0        # Latent eclipse period before lysis (min)
DEFAULT_DELTA_P_PER_MIN: float = 0.004         # Phage decay / mucosal clearance rate
DEFAULT_DELTA_B_PER_MIN: float = 0.0005        # Baseline non-phage bacterial decay
DEFAULT_MU_MUT_PER_MIN: float = 1.2e-6         # Spontaneous oprM forward mutation rate (OprM+ -> OprM-)
DEFAULT_MU_REV_PER_MIN: float = 1.0e-8         # Reversion mutation rate (OprM- -> OprM+)

# Pharmacodynamic Hill Kill Parameters
DEFAULT_E_MAX_PER_MIN: float = 0.08            # Maximum ciprofloxacin kill rate (min^-1)
DEFAULT_HILL_COEFF: float = 2.5                # Hill steepness coefficient H_A
DEFAULT_MIC_S_UG_ML: float = 32.0              # Ciprofloxacin MIC for MexAB-OprM+ wild-type
DEFAULT_MIC_R_UG_ML: float = 1.0               # Ciprofloxacin MIC for OprM- resensitized mutant
DEFAULT_CIPRO_CONC_UG_ML: float = 10.0         # Clinically achievable epithelial lining fluid (ELF) concentration

# Initial Inoculum Reference
DEFAULT_B_S_0_CFU_ML: float = 1.0e7            # Starting chronic CFU/mL
DEFAULT_B_R_0_CFU_ML: float = 1.0e2            # Pre-existing low-frequency mutant pool
DEFAULT_PHAGE_0_PFU_ML: float = 1.0e8          # Therapeutic phage dose (PFU/mL)


@dataclass
class PhageTradeoffParams:
    """Biophysical and pharmacodynamic parameters for OMKO1-MexAB-OprM co-evolution."""

    r_s: float = DEFAULT_R_S_PER_MIN
    c_fit: float = DEFAULT_C_FIT
    carrying_capacity: float = DEFAULT_K_CFU_ML
    k_ads: float = DEFAULT_K_ADS_ML_PER_CFU_MIN
    burst_size: float = DEFAULT_BURST_SIZE
    latent_period_min: float = DEFAULT_LATENT_PERIOD_MIN
    delta_p: float = DEFAULT_DELTA_P_PER_MIN
    delta_b: float = DEFAULT_DELTA_B_PER_MIN
    mu_mut: float = DEFAULT_MU_MUT_PER_MIN
    mu_rev: float = DEFAULT_MU_REV_PER_MIN
    e_max: float = DEFAULT_E_MAX_PER_MIN
    hill_coeff: float = DEFAULT_HILL_COEFF
    mic_s: float = DEFAULT_MIC_S_UG_ML
    mic_r: float = DEFAULT_MIC_R_UG_ML
    cipro_conc_ug_ml: float = DEFAULT_CIPRO_CONC_UG_ML
    b_s_0: float = DEFAULT_B_S_0_CFU_ML
    b_r_0: float = DEFAULT_B_R_0_CFU_ML
    phage_0: float = DEFAULT_PHAGE_0_PFU_ML

    @property
    def r_r(self) -> float:
        """Growth rate of OprM-deficient mutant accounting for fitness cost."""
        return self.r_s * max(0.0, 1.0 - self.c_fit)

    @property
    def resensitization_fold_drop(self) -> float:
        """Fold reduction in antibiotic MIC conferred by OprM loss."""
        if self.mic_r <= 0.0:
            return float("inf")
        return self.mic_s / self.mic_r

    def validate(self) -> None:
        """Validate biophysical bounds."""
        if self.r_s <= 0.0:
            raise ValueError(f"r_s must be positive, got {self.r_s}")
        if not (0.0 <= self.c_fit < 1.0):
            raise ValueError(f"c_fit must be in [0, 1), got {self.c_fit}")
        if self.carrying_capacity <= 0.0:
            raise ValueError(f"carrying_capacity must be positive, got {self.carrying_capacity}")
        if self.k_ads <= 0.0:
            raise ValueError(f"k_ads must be positive, got {self.k_ads}")
        if self.burst_size <= 0.0:
            raise ValueError(f"burst_size must be positive, got {self.burst_size}")
        if self.latent_period_min < 0.0:
            raise ValueError(f"latent_period_min must be non-negative, got {self.latent_period_min}")
        if self.mic_s <= self.mic_r:
            raise ValueError(f"mic_s ({self.mic_s}) must exceed mic_r ({self.mic_r}) for trade-off")


def calculate_hill_kill_rate(
    cipro_conc_ug_ml: float,
    mic_ug_ml: float,
    e_max: float = DEFAULT_E_MAX_PER_MIN,
    hill_coeff: float = DEFAULT_HILL_COEFF,
) -> float:
    """Calculate antibiotic bactericidal kill rate using the Hill pharmacodynamic equation.

    Equation:
      mu_A(C_A, MIC) = E_max * (C_A^H_A) / (MIC^H_A + C_A^H_A)
    """
    if cipro_conc_ug_ml <= 0.0 or mic_ug_ml <= 0.0 or e_max <= 0.0:
        return 0.0
    c_pow = math.pow(cipro_conc_ug_ml, hill_coeff)
    mic_pow = math.pow(mic_ug_ml, hill_coeff)
    denominator = mic_pow + c_pow
    if denominator <= 0.0:
        return 0.0
    return e_max * (c_pow / denominator)


@dataclass
class SimulationPoint:
    """State vector at discrete time t."""

    t_min: float
    t_hours: float
    b_s: float
    b_r: float
    b_tot: float
    phage: float
    fraction_resistant: float
    kill_rate_s: float
    kill_rate_r: float


@dataclass
class SimulationResult:
    """Full trajectory and summary metrics of an OMKO1 co-evolution simulation."""

    scenario: str
    params: PhageTradeoffParams
    duration_hours: float
    dt_min: float
    trajectory: List[SimulationPoint]
    initial_b_tot: float
    final_b_tot: float
    final_b_s: float
    final_b_r: float
    final_phage: float
    peak_phage: float
    final_fraction_resistant: float
    log10_b_tot_reduction: float
    eradication_achieved: bool
    eradication_time_hours: Optional[float]
    evolutionary_vice_score: float

    def to_dict(self) -> Dict[str, Any]:
        """Convert simulation metrics to dictionary (omitting full trajectory for brevity if desired)."""
        return {
            "scenario": self.scenario,
            "duration_hours": self.duration_hours,
            "initial_b_tot_cfu_ml": self.initial_b_tot,
            "final_b_tot_cfu_ml": self.final_b_tot,
            "final_b_s_cfu_ml": self.final_b_s,
            "final_b_r_cfu_ml": self.final_b_r,
            "final_phage_pfu_ml": self.final_phage,
            "peak_phage_pfu_ml": self.peak_phage,
            "final_fraction_resistant": self.final_fraction_resistant,
            "log10_b_tot_reduction": self.log10_b_tot_reduction,
            "eradication_achieved": self.eradication_achieved,
            "eradication_time_hours": self.eradication_time_hours,
            "evolutionary_vice_score": self.evolutionary_vice_score,
            "mic_resensitization_fold": self.params.resensitization_fold_drop,
        }


def simulate_omko1_tradeoff_dynamics(
    params: Optional[PhageTradeoffParams] = None,
    duration_hours: float = 24.0,
    dt_min: float = 0.5,
    scenario_name: str = "custom",
    save_step_interval: int = 10,
    eradication_threshold_cfu: float = 1.0,
) -> SimulationResult:
    """Simulate tripartite delay-differential dynamics of OMKO1 predation and antibiotic kill.

    Integrates:
      dB_S/dt = r_S * B_S * (1 - (B_S + B_R)/K) - k_ads * B_S * P - mu_A(C_A, MIC_S)*B_S - mu_mut*B_S + mu_rev*B_R
      dB_R/dt = r_R * B_R * (1 - (B_S + B_R)/K) - mu_A(C_A, MIC_R)*B_R + mu_mut*B_S - mu_rev*B_R
      dP/dt   = beta * k_ads * B_S(t - tau) * P(t - tau) * exp(-delta_B * tau) - k_ads * B_S * P - delta_P * P
    """
    if params is None:
        params = PhageTradeoffParams()
    params.validate()

    total_steps = int(math.ceil((duration_hours * 60.0) / dt_min))
    delay_steps = max(1, int(round(params.latent_period_min / dt_min)))
    survival_factor = math.exp(-params.delta_b * params.latent_period_min)

    # Initial conditions
    b_s = float(params.b_s_0)
    b_r = float(params.b_r_0)
    phage = float(params.phage_0)

    # Pre-allocate history buffer for delayed burst term: B_S(t-tau) * P(t-tau)
    # For t <= tau, history represents pre-injection state (0 or baseline)
    history_adsorption = [b_s * phage] * (delay_steps + 1)

    trajectory: List[SimulationPoint] = []
    initial_b_tot = b_s + b_r
    peak_phage = phage
    eradication_achieved = False
    eradication_time_hours: Optional[float] = None

    mu_a_s = calculate_hill_kill_rate(
        params.cipro_conc_ug_ml, params.mic_s, params.e_max, params.hill_coeff
    )
    mu_a_r = calculate_hill_kill_rate(
        params.cipro_conc_ug_ml, params.mic_r, params.e_max, params.hill_coeff
    )

    r_r = params.r_r

    for step in range(total_steps + 1):
        t_min = step * dt_min
        t_hours = t_min / 60.0
        b_tot = b_s + b_r

        if phage > peak_phage:
            peak_phage = phage

        if b_tot < eradication_threshold_cfu and not eradication_achieved and step > 0:
            eradication_achieved = True
            eradication_time_hours = t_hours

        # Record trajectory snapshot
        if step % save_step_interval == 0 or step == total_steps:
            fraction_r = (b_r / b_tot) if b_tot > 0.0 else 1.0
            trajectory.append(
                SimulationPoint(
                    t_min=t_min,
                    t_hours=round(t_hours, 3),
                    b_s=max(0.0, b_s),
                    b_r=max(0.0, b_r),
                    b_tot=max(0.0, b_tot),
                    phage=max(0.0, phage),
                    fraction_resistant=fraction_r,
                    kill_rate_s=mu_a_s,
                    kill_rate_r=mu_a_r,
                )
            )

        if step == total_steps:
            break

        # Delayed burst rate from history buffer
        delayed_ads = history_adsorption[0]

        # Derivatives at step
        logistic_s = params.r_s * b_s * (1.0 - (b_tot / params.carrying_capacity))
        logistic_r = r_r * b_r * (1.0 - (b_tot / params.carrying_capacity))

        ads_rate = params.k_ads * b_s * phage
        mutation_rate = params.mu_mut * b_s
        reversion_rate = params.mu_rev * b_r

        db_s = logistic_s - ads_rate - (mu_a_s * b_s) - mutation_rate + reversion_rate
        db_r = logistic_r - (mu_a_r * b_r) + mutation_rate - reversion_rate

        burst_inflow = params.burst_size * delayed_ads * survival_factor
        d_phage = burst_inflow - ads_rate - (params.delta_p * phage)

        # Forward step (Euler with non-negativity constraint)
        b_s = max(0.0, b_s + db_s * dt_min)
        b_r = max(0.0, b_r + db_r * dt_min)
        phage = max(0.0, phage + d_phage * dt_min)

        # Update history buffer (FIFO queue)
        history_adsorption.pop(0)
        history_adsorption.append(params.k_ads * b_s * phage)

    final_b_tot = b_s + b_r

    if final_b_tot < eradication_threshold_cfu:
        # Complete sterile eradication: zero resistant or sensitive bacteria remain
        final_fraction_r = 0.0
        log10_reduction = (
            math.log10(initial_b_tot / eradication_threshold_cfu) if initial_b_tot > 0.0 else 12.0
        )
        vice_score = max(5.0, log10_reduction)
    else:
        final_fraction_r = (b_r / final_b_tot) if final_b_tot > 0.0 else 0.0
        if initial_b_tot > final_b_tot:
            log10_reduction = math.log10(initial_b_tot / final_b_tot)
        else:
            log10_reduction = 0.0
        # Evolutionary vice score: rewards overall kill while heavily penalizing resistant escape
        vice_score = max(0.0, log10_reduction * (1.0 - final_fraction_r))

    return SimulationResult(
        scenario=scenario_name,
        params=params,
        duration_hours=duration_hours,
        dt_min=dt_min,
        trajectory=trajectory,
        initial_b_tot=initial_b_tot,
        final_b_tot=final_b_tot,
        final_b_s=b_s,
        final_b_r=b_r,
        final_phage=phage,
        peak_phage=peak_phage,
        final_fraction_resistant=final_fraction_r,
        log10_b_tot_reduction=round(log10_reduction, 4),
        eradication_achieved=eradication_achieved,
        eradication_time_hours=round(eradication_time_hours, 3) if eradication_time_hours else None,
        evolutionary_vice_score=round(vice_score, 4),
    )


def compare_standard_regimens(
    duration_hours: float = 24.0,
    cipro_conc_ug_ml: float = DEFAULT_CIPRO_CONC_UG_ML,
    phage_dose_pfu_ml: float = DEFAULT_PHAGE_0_PFU_ML,
) -> Dict[str, SimulationResult]:
    """Execute the four standard therapeutic comparator arms to demonstrate the evolutionary trap:

    1. Untreated Control: C_A = 0, P_0 = 0
    2. Ciprofloxacin Monotherapy: C_A = 10 ug/mL, P_0 = 0 (fails due to high MIC_S = 32 ug/mL)
    3. Phage OMKO1 Monotherapy: C_A = 0, P_0 = 1e8 PFU/mL (suppresses B_S, selects for B_R outgrowth)
    4. Dual Synergy (Evolutionary Trap): C_A = 10 ug/mL, P_0 = 1e8 PFU/mL (simultaneous eradication)
    """
    # 1. Untreated Control
    params_ctrl = PhageTradeoffParams(
        cipro_conc_ug_ml=0.0,
        phage_0=0.0,
    )
    res_ctrl = simulate_omko1_tradeoff_dynamics(
        params_ctrl, duration_hours=duration_hours, scenario_name="Untreated Control"
    )

    # 2. Ciprofloxacin Monotherapy
    params_cipro = PhageTradeoffParams(
        cipro_conc_ug_ml=cipro_conc_ug_ml,
        phage_0=0.0,
    )
    res_cipro = simulate_omko1_tradeoff_dynamics(
        params_cipro, duration_hours=duration_hours, scenario_name="Ciprofloxacin Monotherapy"
    )

    # 3. Phage OMKO1 Monotherapy
    params_phage = PhageTradeoffParams(
        cipro_conc_ug_ml=0.0,
        phage_0=phage_dose_pfu_ml,
    )
    res_phage = simulate_omko1_tradeoff_dynamics(
        params_phage, duration_hours=duration_hours, scenario_name="Phage OMKO1 Monotherapy"
    )

    # 4. Dual Phage-Antibiotic Synergy (Evolutionary Trap)
    params_dual = PhageTradeoffParams(
        cipro_conc_ug_ml=cipro_conc_ug_ml,
        phage_0=phage_dose_pfu_ml,
    )
    res_dual = simulate_omko1_tradeoff_dynamics(
        params_dual, duration_hours=duration_hours, scenario_name="Dual Evolutionary Trap (OMKO1 + Cipro)"
    )

    return {
        "untreated": res_ctrl,
        "cipro_monotherapy": res_cipro,
        "phage_monotherapy": res_phage,
        "dual_evolutionary_trap": res_dual,
    }


def evaluate_synergy_dose_matrix(
    cipro_levels: Optional[List[float]] = None,
    phage_levels: Optional[List[float]] = None,
    duration_hours: float = 24.0,
) -> Dict[str, Any]:
    """Evaluate a grid of ciprofloxacin concentrations and phage doses.

    Quantifies fractional inhibitory synergy, log reduction, and B_R outgrowth suppression.
    """
    if cipro_levels is None:
        cipro_levels = [0.0, 1.0, 5.0, 10.0, 20.0]
    if phage_levels is None:
        phage_levels = [0.0, 1e6, 1e7, 1e8, 1e9]

    matrix_entries = []

    for c_cipro in cipro_levels:
        for p_dose in phage_levels:
            p = PhageTradeoffParams(cipro_conc_ug_ml=c_cipro, phage_0=p_dose)
            label = f"Cipro_{c_cipro}ug_Phage_{p_dose:.0e}"
            res = simulate_omko1_tradeoff_dynamics(p, duration_hours=duration_hours, scenario_name=label)
            matrix_entries.append({
                "cipro_conc_ug_ml": c_cipro,
                "phage_dose_pfu_ml": p_dose,
                "final_b_tot_cfu_ml": res.final_b_tot,
                "final_b_s_cfu_ml": res.final_b_s,
                "final_b_r_cfu_ml": res.final_b_r,
                "fraction_resistant": res.final_fraction_resistant,
                "log10_reduction": res.log10_b_tot_reduction,
                "eradication_achieved": res.eradication_achieved,
                "eradication_time_hours": res.eradication_time_hours,
                "vice_score": res.evolutionary_vice_score,
            })

    return {
        "cipro_levels_ug_ml": cipro_levels,
        "phage_levels_pfu_ml": phage_levels,
        "entries": matrix_entries,
        "total_simulations": len(matrix_entries),
    }


def format_markdown_report(comparisons: Dict[str, SimulationResult]) -> str:
    """Format comparative analysis into markdown."""
    lines = [
        "# Phage OMKO1 Evolutionary Trade-off & MexAB-OprM Resensitization Analysis",
        "",
        f"> **Regulatory Notice**: {NLM_REGISTRY_DISCLAIMER}",
        f"> **Research Scope**: {RESEARCH_BOUNDARY}",
        "",
        "## 1. Primary Evidence Grounding",
        "- **Chan et al. (2018)**, *Sci Rep*, DOI: `10.1038/s41598-018-24541-1`",
        "- **Chan et al. (2021)**, *J Cyst Fibros*, DOI: `10.1016/j.jcf.2021.03.013`",
        "- **CYPHY Trial**: `NCT04684641` (Yale CF Bacteriophage Study)",
        "- **Inhaled Bacteriophage Trial**: `NCT05010577`",
        "",
        "## 2. Comparative Regimen Outcomes (24-Hour Endpoint)",
        "",
        "| Regimen | Final CFU/mL | Log10 Reduction | % B_R Resistant | Eradication (<1 CFU) | Vice Score |",
        "| :--- | :--- | :--- | :--- | :--- | :--- |",
    ]

    for key, res in comparisons.items():
        pct_r = res.final_fraction_resistant * 100.0
        erad_str = f"Yes ({res.eradication_time_hours}h)" if res.eradication_achieved else "No"
        lines.append(
            f"| {res.scenario} | {res.final_b_tot:.2e} | {res.log10_b_tot_reduction:.2f} | "
            f"{pct_r:.1f}% | {erad_str} | {res.evolutionary_vice_score:.2f} |"
        )

    lines.extend([
        "",
        "## 3. Biophysical Insights & Evolutionary Vice Mechanism",
        "1. **Ciprofloxacin Alone**: Ineffective against MexAB-OprM+ wild-type (MIC 32 ug/mL), allowing biofilm expansion.",
        "2. **Phage OMKO1 Alone**: Successfully lyses B_S, but creates selective pressure driving mutant B_R (OprM-) outgrowth.",
        "3. **Dual Evolutionary Trap**: OMKO1 forces OprM deletion/alteration (32-fold MIC reduction to 1.0 ug/mL), rendering "
        "the bacterial population exquisitely vulnerable to ciprofloxacin (10 ug/mL) and driving sterile clearance.",
        "",
    ])
    return "\n".join(lines)


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Pseudomonas aeruginosa Phage OMKO1 Evolutionary Trade-off Model (T-0246)",
    )
    parser.add_argument(
        "--duration-hours",
        type=float,
        default=24.0,
        help="Simulation run time in hours (default: 24.0)",
    )
    parser.add_argument(
        "--cipro-ug-ml",
        type=float,
        default=DEFAULT_CIPRO_CONC_UG_ML,
        help=f"Ciprofloxacin concentration in ug/mL (default: {DEFAULT_CIPRO_CONC_UG_ML})",
    )
    parser.add_argument(
        "--phage-pfu-ml",
        type=float,
        default=DEFAULT_PHAGE_0_PFU_ML,
        help=f"Initial phage dose in PFU/mL (default: {DEFAULT_PHAGE_0_PFU_ML})",
    )
    parser.add_argument(
        "--dt-min",
        type=float,
        default=0.5,
        help="Integration time step in minutes (default: 0.5)",
    )
    parser.add_argument(
        "--matrix",
        action="store_true",
        help="Evaluate full 2D dose-synergy matrix",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Output raw JSON data",
    )
    return parser.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    """CLI Entry point."""
    args = parse_args(argv)

    if args.matrix:
        matrix_data = evaluate_synergy_dose_matrix(duration_hours=args.duration_hours)
        if args.json:
            print(json.dumps(matrix_data, indent=2))
        else:
            print(f"Evaluated {matrix_data['total_simulations']} dose combinations across cipro and OMKO1.")
        return 0

    comparisons = compare_standard_regimens(
        duration_hours=args.duration_hours,
        cipro_conc_ug_ml=args.cipro_ug_ml,
        phage_dose_pfu_ml=args.phage_pfu_ml,
    )

    if args.json:
        payload = {
            "disclaimer": NLM_REGISTRY_DISCLAIMER,
            "boundary": RESEARCH_BOUNDARY,
            "regimens": {k: v.to_dict() for k, v in comparisons.items()},
        }
        print(json.dumps(payload, indent=2))
    else:
        report = format_markdown_report(comparisons)
        print(report)

    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
