#!/usr/bin/env python3
"""Apical-to-Basolateral Receptor Transcytosis Biophysical & Kinetic Model (T-0177).

Models the receptor-mediated transcytosis kinetics of aerosolized targeted nanoparticles
across polarized human airway epithelial cells to reach the basolateral basal stem cell compartment.
Simulates pH-dependent binding, endosomal sorting, lysosomal degradation, and basolateral exocytosis
for Neonatal Fc Receptor (FcRn), Transferrin Receptor (TfR), and Polymeric Immunoglobulin Receptor (pIgR).

Primary references:
- Yu et al., Sci Transl Med 2011 (DOI: 10.1126/scitranslmed.3002230)
- Spiekermann et al., J Exp Med 2002 (DOI: 10.1084/jem.20020400)
- Bitonti & Dumont, Adv Drug Deliv Rev 2006 (DOI: 10.1016/j.addr.2006.07.019)
- PDB IDs: 1CX8 (TfR), 6OKD (TfR-Ferritin), 4N0F (FcRn-Fc)
"""
from __future__ import annotations

import argparse
import json
import math
from typing import Dict, Any, List, Tuple


# Receptors and empirical kinetic parameters
RECEPTOR_PROFILES: Dict[str, Dict[str, Any]] = {
    "fcrn": {
        "name": "Neonatal Fc Receptor (FcRn)",
        "pdb_id": "4N0F",
        "description": "pH-dependent recycling receptor. High affinity at apical ASL pH 6.0-6.5; rapid release at basolateral pH 7.4.",
        "apical_surface_density_per_cell": 45000,
        "kd_apical_nm": 220.0,       # At pH 6.2 - 6.5
        "kd_basolateral_nm": 15000.0, # At pH 7.4 (virtually complete release)
        "k_int_per_min": 0.12,       # Internalization rate constant
        "k_sort_recycling": 0.08,    # Sorting to Rab11a exocytosis pathway
        "k_sort_lysosome": 0.015,    # Lysosomal routing
        "optimal_kd_window_nm": (150.0, 350.0),
    },
    "tfr": {
        "name": "Transferrin Receptor 1 (TfR1)",
        "pdb_id": "1CX8",
        "description": "Constitutive endocytic receptor. Low-affinity monovalent binders escape lysosomes; high-affinity binders degrade.",
        "apical_surface_density_per_cell": 28000,
        "kd_apical_nm": 250.0,       # Monovalent tuned affinity
        "kd_basolateral_nm": 250.0,
        "k_int_per_min": 0.18,
        "k_sort_recycling": 0.06,
        "k_sort_lysosome": 0.035,
        "optimal_kd_window_nm": (180.0, 450.0),
    },
    "pigr": {
        "name": "Polymeric Immunoglobulin Receptor (pIgR)",
        "pdb_id": "5D4K",
        "description": "Mucosal transcytosis receptor for dimeric IgA and secretory component targeting.",
        "apical_surface_density_per_cell": 18000,
        "kd_apical_nm": 310.0,
        "kd_basolateral_nm": 420.0,
        "k_int_per_min": 0.09,
        "k_sort_recycling": 0.05,
        "k_sort_lysosome": 0.02,
        "optimal_kd_window_nm": (200.0, 500.0),
    },
}


def calculate_fractional_occupancy(ligand_conc_nm: float, kd_nm: float) -> float:
    """Calculate steady-state receptor fractional occupancy [RL] / [R_total]."""
    if ligand_conc_nm < 0 or kd_nm <= 0:
        raise ValueError("Ligand concentration and Kd must be non-negative / positive")
    return ligand_conc_nm / (ligand_conc_nm + kd_nm)


def simulate_transcytosis_kinetics(
    receptor_id: str,
    ligand_conc_nm: float = 100.0,
    custom_kd_nm: float | None = None,
    time_minutes: float = 120.0,
    dt_minutes: float = 0.5,
) -> Dict[str, Any]:
    """Simulate transcytosis compartment dynamics over time.
    
    Compartments:
    - L_apical: Free apical ligand / nanoparticle pool
    - RL_surface: Receptor-bound complex on apical plasma membrane
    - Endosome_sort: Internalized sorting endosome pool
    - Lysosome: Degraded payload pool
    - L_basolateral: Successfully exocytosed cargo reaching basolateral basal stem cells
    """
    rec_id = receptor_id.lower()
    if rec_id not in RECEPTOR_PROFILES:
        raise ValueError(f"Unknown receptor '{receptor_id}'. Must be one of: {list(RECEPTOR_PROFILES.keys())}")
    
    profile = RECEPTOR_PROFILES[rec_id]
    kd_apical = custom_kd_nm if custom_kd_nm is not None else profile["kd_apical_nm"]
    
    # Internalization and sorting rates
    k_int = profile["k_int_per_min"]
    
    # High affinity (low Kd < 50 nM) increases lysosomal sorting due to failure of endosomal dissociation
    if kd_apical < 50.0:
        lysosome_bias = 2.5 * (50.0 / max(1.0, kd_apical))
        k_lys = profile["k_sort_lysosome"] * lysosome_bias
        k_recyc = profile["k_sort_recycling"] * 0.4
    else:
        k_lys = profile["k_sort_lysosome"]
        k_recyc = profile["k_sort_recycling"]
        
    k_exocytosis = 0.15  # Rate of exocytosis from recycling endosomes to basolateral space

    # Initial conditions (arbitrary units normalized to 100% total added dose)
    l_apical = 100.0
    rl_surface = 0.0
    endo_sort = 0.0
    rab11_recyc = 0.0
    lyso = 0.0
    l_basolateral = 0.0
    
    steps = int(time_minutes / dt_minutes)
    time_points = []
    basolateral_trajectory = []
    lysosomal_trajectory = []
    
    for step in range(steps + 1):
        current_time = step * dt_minutes
        time_points.append(current_time)
        basolateral_trajectory.append(round(l_basolateral, 3))
        lysosomal_trajectory.append(round(lyso, 3))
        
        # Equilibrium binding fraction on apical surface
        occupancy = calculate_fractional_occupancy(l_apical, kd_apical)
        d_bind = occupancy * 0.25 * l_apical * dt_minutes
        l_apical = max(0.0, l_apical - d_bind)
        rl_surface += d_bind
        
        # Internalization: RL_surface -> Endosome_sort
        d_int = rl_surface * k_int * dt_minutes
        rl_surface = max(0.0, rl_surface - d_int)
        endo_sort += d_int
        
        # Sorting: Endosome_sort -> Lysosome vs. Rab11a recycling pool
        d_lyso = endo_sort * k_lys * dt_minutes
        d_recyc = endo_sort * k_recyc * dt_minutes
        endo_sort = max(0.0, endo_sort - d_lyso - d_recyc)
        lyso += d_lyso
        rab11_recyc += d_recyc
        
        # Exocytosis to basolateral space from Rab11a recycling pool
        d_exo = rab11_recyc * k_exocytosis * dt_minutes
        rab11_recyc = max(0.0, rab11_recyc - d_exo)
        l_basolateral += d_exo

    total_internalized = l_basolateral + lyso + endo_sort + rab11_recyc
    transcytosis_efficiency_pct = (l_basolateral / max(0.001, total_internalized)) * 100.0
    lysosomal_escape_pct = 100.0 - ((lyso / max(0.001, total_internalized)) * 100.0)

    opt_min, opt_max = profile["optimal_kd_window_nm"]
    in_optimal_window = opt_min <= kd_apical <= opt_max

    return {
        "receptor": profile["name"],
        "receptor_id": rec_id,
        "pdb_id": profile["pdb_id"],
        "simulated_time_min": time_minutes,
        "tested_kd_nm": round(kd_apical, 2),
        "optimal_kd_window_nm": profile["optimal_kd_window_nm"],
        "in_optimal_affinity_window": in_optimal_window,
        "transcytosis_efficiency_pct": round(transcytosis_efficiency_pct, 2),
        "lysosomal_escape_pct": round(lysosomal_escape_pct, 2),
        "basolateral_delivery_yield_pct": round(l_basolateral, 2),
        "lysosomal_degradation_pct": round(lyso, 2),
        "apical_retention_pct": round(l_apical, 2),
        "verdict": (
            "OPTIMAL_TRANSCYTOSIS" if in_optimal_window and transcytosis_efficiency_pct >= 45.0
            else "HIGH_AFFINITY_LYSOSOMAL_TRAPPING" if kd_apical < opt_min
            else "LOW_AFFINITY_POOR_UPTAKE"
        ),
    }



class TranscytosisModel:
    """Apical-to-basolateral receptor transcytosis flux model."""

    DEFAULT_CONCENTRATIONS_NM: List[float] = [10.0, 25.0, 50.0, 100.0, 200.0, 300.0, 500.0, 1000.0]

    def __init__(self, receptors: List[str] | None = None):
        self.receptors = [r.lower() for r in (receptors or ["fcrn", "tfr", "pigr"])]
        for r in self.receptors:
            if r not in RECEPTOR_PROFILES:
                raise ValueError(f"Unknown receptor '{r}'. Must be one of: {list(RECEPTOR_PROFILES.keys())}")

    def simulate_flux(
        self,
        concentrations_nm: List[float] | None = None,
        time_minutes: float = 120.0,
    ) -> List[Dict[str, Any]]:
        """Simulate transcytosis flux across tested receptors and concentration range."""
        concs = concentrations_nm or self.DEFAULT_CONCENTRATIONS_NM
        results = []
        for rec_id in self.receptors:
            prof = RECEPTOR_PROFILES[rec_id]
            kd = prof["kd_apical_nm"]
            for conc in concs:
                occupancy = calculate_fractional_occupancy(conc, kd)
                surface_binding_pct = round(occupancy * 100.0, 2)
                sim = simulate_transcytosis_kinetics(
                    receptor_id=rec_id,
                    ligand_conc_nm=conc,
                    time_minutes=time_minutes,
                )
                internalized_pct = round(100.0 - sim["apical_retention_pct"], 2)
                basolateral_exocytosis_pct = sim["basolateral_delivery_yield_pct"]
                hours = max(0.1, time_minutes / 60.0)
                net_flux_pct_per_h = round(basolateral_exocytosis_pct / hours, 3)

                results.append({
                    "receptor_id": rec_id,
                    "receptor_name": prof["name"],
                    "pdb_id": prof["pdb_id"],
                    "concentration_nm": float(conc),
                    "surface_binding_pct": surface_binding_pct,
                    "internalized_pct": internalized_pct,
                    "basolateral_exocytosis_pct": basolateral_exocytosis_pct,
                    "net_flux_pct_per_h": net_flux_pct_per_h,
                    "transcytosis_efficiency_pct": sim["transcytosis_efficiency_pct"],
                    "lysosomal_escape_pct": sim["lysosomal_escape_pct"],
                    "verdict": sim["verdict"],
                })
        return results


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Apical-to-Basolateral Receptor Transcytosis Model")
    parser.add_argument("--receptor", default="fcrn", choices=list(RECEPTOR_PROFILES.keys()), help="Target receptor")
    parser.add_argument("--kd-nm", type=float, default=None, help="Custom ligand Kd in nM")
    parser.add_argument("--time-min", type=float, default=120.0, help="Simulation duration in minutes")
    parser.add_argument("--json", action="store_true", help="Print JSON output")
    args = parser.parse_args(argv)

    res = simulate_transcytosis_kinetics(
        receptor_id=args.receptor,
        custom_kd_nm=args.kd_nm,
        time_minutes=args.time_min
    )
    if args.json:
        print(json.dumps(res, indent=2))
    else:
        print(f"=== Receptor Transcytosis Simulation: {res['receptor']} (PDB: {res['pdb_id']}) ===")
        print(f"Tested Kd:                     {res['tested_kd_nm']} nM")
        print(f"Optimal Kd Window:             {res['optimal_kd_window_nm'][0]} - {res['optimal_kd_window_nm'][1]} nM")
        print(f"In Optimal Window:             {res['in_optimal_affinity_window']}")
        print(f"Transcytosis Efficiency:       {res['transcytosis_efficiency_pct']}%")
        print(f"Lysosomal Escape Rate:         {res['lysosomal_escape_pct']}%")
        print(f"Net Basolateral Delivery Yield:{res['basolateral_delivery_yield_pct']}% of dose")
        print(f"Verdict:                       {res['verdict']}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())

