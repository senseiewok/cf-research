#!/usr/bin/env python3
"""Airway Secretory-to-Basal Dedifferentiation & Epithelial Plasticity Model (T-0178).

Simulates the 'dedifferentiation bypass' paradigm for Cystic Fibrosis curative therapy:
Proving that targeting directly accessible apical secretory/club cells (SCGB1A1+) paired
with a transient regenerative stimulus (such as a brief Notch/gamma-secretase inhibitor pulse)
permanently repopulates the basal stem cell niche (KRT5+/TP63+) without requiring tight junction disruption.

Primary literature references:
- Tata et al., Nature 2013 (DOI: 10.1038/nature12777)
- Pardo-Saganta et al., Nature 2015 (DOI: 10.1038/nature14548)
- Rawlins et al., Cell Stem Cell 2009 (DOI: 10.1016/j.stem.2009.04.002)
- Mou et al., Cell Stem Cell 2016 (DOI: 10.1016/j.stem.2016.11.001)
"""
from __future__ import annotations

import argparse
import json
import math
from typing import Dict, Any, List, Tuple


# Baseline human pseudostratified airway epithelial architecture
AIRWAY_EPITHELIAL_BASELINE = {
    "total_cells_per_cm2": 1.2e6,
    "cell_fractions": {
        "basal_stem_cells": 0.30,       # KRT5+, TP63+, NGFR+
        "secretory_club_cells": 0.25,   # SCGB1A1+, MUC5AC+
        "ciliated_cells": 0.40,         # FOXJ1+, acetylated alpha-tubulin+
        "ionocytes_and_rare": 0.05,     # FOXI1+, CFTR-high
    },
    "turnover_half_life_days": {
        "ciliated_cells": 60.0,
        "secretory_cells": 40.0,
        "basal_stem_cells": 180.0,      # Quiescent stem cell pool
    },
    "therapeutic_cftr_threshold_pct": 15.0, # Minimum functional CFTR expression for clinical benefit
}


def simulate_dedifferentiation_bypass(
    apical_editing_efficiency: float = 0.35,  # 35% of accessible secretory cells edited
    notch_inhibitor_pulse_days: float = 4.0,  # 4-day transient DAPT/DBZ pulse
    dediff_rate_per_day: float = 0.25,        # Rate of secretory-to-basal conversion during pulse (Tata et al. 2013)
    total_days: float = 180.0,
    dt_days: float = 0.5,
) -> Dict[str, Any]:
    """Simulate long-term epithelial dynamics and basal stem cell niche repopulation.
    
    Variables:
    - B_wt, B_ed: Basal stem cells (Wild-type vs Edited)
    - S_wt, S_ed: Secretory club cells (Wild-type vs Edited)
    - C_wt, C_ed: Ciliated cells (Wild-type vs Edited)
    """
    if not (0.0 <= apical_editing_efficiency <= 1.0):
        raise ValueError("Apical editing efficiency must be between 0.0 and 1.0")
    if notch_inhibitor_pulse_days < 0 or total_days <= 0:
        raise ValueError("Durations must be positive")

    # Initial state (normalized to total cell count 100.0)
    b_total = 100.0 * AIRWAY_EPITHELIAL_BASELINE["cell_fractions"]["basal_stem_cells"]
    s_total = 100.0 * AIRWAY_EPITHELIAL_BASELINE["cell_fractions"]["secretory_club_cells"]
    c_total = 100.0 * AIRWAY_EPITHELIAL_BASELINE["cell_fractions"]["ciliated_cells"]

    b_wt = b_total
    b_ed = 0.0

    # Apical editing applied specifically to exposed luminal secretory cells
    s_ed = s_total * apical_editing_efficiency
    s_wt = s_total - s_ed

    c_wt = c_total
    c_ed = 0.0

    # Basal niche carrying capacity
    b_capacity = b_total

    # Differentiation rates from basal stem cells (per day)
    k_diff_secretory = 0.015
    k_diff_ciliated = 0.018

    # Apical cell loss / turnover rates (per day)
    k_turnover_cilia = math.log(2.0) / AIRWAY_EPITHELIAL_BASELINE["turnover_half_life_days"]["ciliated_cells"]
    k_turnover_secretory = math.log(2.0) / AIRWAY_EPITHELIAL_BASELINE["turnover_half_life_days"]["secretory_cells"]

    steps = int(total_days / dt_days)
    trajectory_days = []
    basal_chimerism = []
    total_epithelial_cftr_pct = []

    for step in range(steps + 1):
        t = step * dt_days
        current_b_total = b_wt + b_ed
        chimerism_pct = (b_ed / max(0.001, current_b_total)) * 100.0
        
        current_total_cells = (b_wt + b_ed) + (s_wt + s_ed) + (c_wt + c_ed)
        current_edited_cells = b_ed + s_ed + c_ed
        cftr_functional_pct = (current_edited_cells / max(0.001, current_total_cells)) * 100.0

        trajectory_days.append(round(t, 1))
        basal_chimerism.append(round(chimerism_pct, 2))
        total_epithelial_cftr_pct.append(round(cftr_functional_pct, 2))

        # Check if currently within transient Notch inhibition pulse window
        in_notch_pulse = (notch_inhibitor_pulse_days > 0.0) and (t <= notch_inhibitor_pulse_days)

        if in_notch_pulse:
            # Notch inhibition triggers dedifferentiation of secretory cells to basal stem cells (Tata et al. 2013)
            dediff_flux_wt = s_wt * dediff_rate_per_day * dt_days
            dediff_flux_ed = s_ed * dediff_rate_per_day * dt_days
            
            s_wt = max(0.0, s_wt - dediff_flux_wt)
            s_ed = max(0.0, s_ed - dediff_flux_ed)
            
            # Competitively displace or enter vacant basal niche
            b_wt += dediff_flux_wt
            b_ed += dediff_flux_ed

        # Normalization / competition in basal niche to stay near carrying capacity
        b_sum = b_wt + b_ed
        if b_sum > b_capacity:
            excess_ratio = b_capacity / b_sum
            b_wt *= excess_ratio
            b_ed *= excess_ratio

        # Differentiation flux from basal stem cells according to current chimerism
        f_ed_basal = b_ed / max(0.001, b_wt + b_ed)
        f_wt_basal = b_wt / max(0.001, b_wt + b_ed)

        flux_to_s = b_capacity * k_diff_secretory * dt_days
        flux_to_c = b_capacity * k_diff_ciliated * dt_days

        s_wt += flux_to_s * f_wt_basal
        s_ed += flux_to_s * f_ed_basal

        c_wt += flux_to_c * f_wt_basal
        c_ed += flux_to_c * f_ed_basal

        # Turnover of differentiated apical cells
        s_wt = max(0.0, s_wt - s_wt * k_turnover_secretory * dt_days)
        s_ed = max(0.0, s_ed - s_ed * k_turnover_secretory * dt_days)

        c_wt = max(0.0, c_wt - c_wt * k_turnover_cilia * dt_days)
        c_ed = max(0.0, c_ed - c_ed * k_turnover_cilia * dt_days)

    final_basal_chimerism = basal_chimerism[-1]
    final_epithelial_cftr = total_epithelial_cftr_pct[-1]
    meets_therapeutic_threshold = final_epithelial_cftr >= AIRWAY_EPITHELIAL_BASELINE["therapeutic_cftr_threshold_pct"]

    return {
        "simulation_days": total_days,
        "apical_editing_efficiency_pct": round(apical_editing_efficiency * 100.0, 1),
        "notch_inhibitor_pulse_days": notch_inhibitor_pulse_days,
        "dedifferentiation_rate_per_day": dediff_rate_per_day,
        "final_basal_stem_cell_chimerism_pct": final_basal_chimerism,
        "final_epithelial_cftr_functional_pct": final_epithelial_cftr,
        "therapeutic_cftr_threshold_pct": AIRWAY_EPITHELIAL_BASELINE["therapeutic_cftr_threshold_pct"],
        "meets_therapeutic_threshold": meets_therapeutic_threshold,
        "verdict": (
            "PERMANENT_CURATIVE_REPOPULATION" if meets_therapeutic_threshold and final_basal_chimerism >= 15.0
            else "PARTIAL_REPOPULATION_SUBTHRESHOLD" if final_basal_chimerism > 0.0
            else "EPHEMERAL_EXPRESSION_NO_REPOPULATION"
        ),
        "trajectory_sample_day30": total_epithelial_cftr_pct[min(len(total_epithelial_cftr_pct)-1, int(30.0 / dt_days))],
        "trajectory_sample_day90": total_epithelial_cftr_pct[min(len(total_epithelial_cftr_pct)-1, int(90.0 / dt_days))],
        "trajectory_sample_day180": final_epithelial_cftr,
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Airway Secretory-to-Basal Dedifferentiation Plasticity Model")
    parser.add_argument("--edit-eff", type=float, default=0.35, help="Apical secretory editing efficiency (0.0 - 1.0)")
    parser.add_argument("--pulse-days", type=float, default=4.0, help="Notch inhibitor pulse duration in days")
    parser.add_argument("--total-days", type=float, default=180.0, help="Total simulation time in days")
    parser.add_argument("--json", action="store_true", help="Print JSON output")
    args = parser.parse_args(argv)

    res = simulate_dedifferentiation_bypass(
        apical_editing_efficiency=args.edit_eff,
        notch_inhibitor_pulse_days=args.pulse_days,
        total_days=args.total_days,
    )
    if args.json:
        print(json.dumps(res, indent=2))
    else:
        print("=== Airway Secretory-to-Basal Dedifferentiation Simulation ===")
        print(f"Apical Secretory Editing:      {res['apical_editing_efficiency_pct']}%")
        print(f"Notch Inhibitor Pulse:         {res['notch_inhibitor_pulse_days']} days")
        print(f"Final Basal Stem Cell Chimerism:{res['final_basal_stem_cell_chimerism_pct']}%")
        print(f"Final Epithelial CFTR Level:   {res['final_epithelial_cftr_functional_pct']}%")
        print(f"Therapeutic Threshold (15%):   {res['meets_therapeutic_threshold']}")
        print(f"Long-Term Persistence (Day 180): {res['trajectory_sample_day180']}%")
        print(f"Verdict:                       {res['verdict']}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
