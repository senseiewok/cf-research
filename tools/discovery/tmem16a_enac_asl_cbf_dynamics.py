#!/usr/bin/env python3
"""Theratyping Pharmacodynamic Model of Dual TMEM16A Potentiation & ENaC Blockade Synergy across CBF & ASL Height (T-0444).

Quantitative pharmacodynamic simulation of TMEM16A (ANO1) apical chloride/bicarbonate
conductance potentiation coupled with sub-maximal ENaC sodium hyperabsorption inhibition,
solving osmotic water flux to predict airway surface liquid (ASL) height recovery
(from collapsed 2.5 um to >= 7.0 um) and ciliary beat frequency (CBF >= 8.0 Hz)
across modulator-ineligible Class I and Class VII CFTR genotypes.

Primary Literature Grounding:
  - Paulino C et al. (2017) Nature, DOI: 10.1038/nature24652
    (Structural basis of TMEM16A channel architecture and activation)
  - Noreng S et al. (2022) Cell, DOI: 10.1016/j.cell.2022.04.024
    (Structure and regulation of human epithelial sodium channel ENaC)
  - Danahay H et al. (2024) J Physiol, DOI: 10.1113/JP286559
    (Functional potentiation of TMEM16A in human cystic fibrosis airway epithelia)
  - Boucher RC et al. (2019) Am J Respir Crit Care Med, DOI: 10.1164/rccm.201908-1641OC
    (Airway surface liquid hydration and mucociliary transport in cystic fibrosis)
  - Button B et al. (2021) Physiol Rev, DOI: 10.1152/physrev.00004.2021
    (Physiology and pathophysiology of human airway mucus)

Clinical Trials:
  - NCT05668741 (Phase 1/2 Study of Inhaled CFTR mRNA LNP VX-522 in Cystic Fibrosis)
  - NCT05248230 (Phase 1/2 Study of Aerosolized AAV Gene Therapy 4D-710 in Cystic Fibrosis)

Structural Coordinates:
  - RCSB PDB 6MSM (Dephosphorylated human CFTR)
  - RCSB PDB 5UAK (Phosphorylated active human CFTR)

Research Boundary:
  Preclinical computational biophysical simulation; research tool, not medical advice.
  No clinical efficacy or therapeutic safety claim is made.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
import sys
from typing import Any, Dict, List, Optional, Tuple

# Mandatory ClinicalTrials.gov Attribution
NLM_ATTRIBUTION = (
    "ClinicalTrials.gov does not evaluate safety or efficacy; "
    "listing on the registry does not denote government endorsement."
)

RESEARCH_BOUNDARY = (
    "Preclinical computational biophysical simulation; research, not medical advice. "
    "No clinical efficacy or therapeutic safety claim is made."
)

# Registry & Structural Anchors
CLINICAL_TRIAL_VX522 = "NCT05668741"
CLINICAL_TRIAL_4D710 = "NCT05248230"
PDB_DEPHOSPHO_CFTR = "6MSM"
PDB_ACTIVE_CFTR = "5UAK"

# Primary DOIs
PRIMARY_DOI_TMEM16A_STRUCT = "10.1038/nature24652"
PRIMARY_DOI_ENAC_STRUCT = "10.1016/j.cell.2022.04.024"
PRIMARY_DOI_TMEM16A_POTENTIATION = "10.1113/JP286559"
PRIMARY_DOI_ASL_HYDRATION = "10.1164/rccm.201908-1641OC"
PRIMARY_DOI_AIRWAY_MUCUS = "10.1152/physrev.00004.2021"

# Universal Physical Constants
GAS_CONSTANT_R = 8.314462618  # J / (mol * K)
FARADAY_CONSTANT_F = 96485.3321  # C / mol
BODY_TEMPERATURE_K = 310.15  # 37 C
WATER_MOLAR_VOLUME_VW = 1.8e-5  # m^3 / mol (18 cm^3 / mol)


@dataclass(frozen=True)
class AirwayBiophysicalConstants:
    """Standard physiological parameters of human bronchial epithelia."""
    # Membrane potentials (mV)
    v_apical_mv: float = -40.0
    e_cl_mv: float = +15.0  # Secretory chemical equilibrium potential across apical membrane
    e_hco3_mv: float = +10.0
    e_na_mv: float = +60.0

    # Baseline conductances (mS/cm^2)
    g_cftr_wt: float = 2.50
    g_tmem16a_max: float = 3.80
    g_enac_max: float = 4.20

    # Permeability ratios
    p_hco3_over_cl_tmem16a: float = 0.20
    p_hco3_over_cl_cftr: float = 0.25

    # Osmotic and fluid dynamics
    osmotic_water_permeability_pf_um_s: float = 420.0  # um/s
    evaporation_rate_um_min: float = 0.05  # um/min baseline airway fluid loss
    baseline_cf_asl_height_um: float = 2.50  # Collapsed CF ASL depth
    target_clearance_asl_height_um: float = 7.00  # Physiological mucociliary threshold
    wildtype_asl_height_um: float = 8.00  # Healthy reference depth
    tau_hydration_min: float = 12.0  # Characteristic ASL rehydration time constant (min)

    # Ciliary mechanics
    cilia_axoneme_length_um: float = 6.00
    cbf_quiescent_hz: float = 2.50  # Impaired CBF under collapsed mucus
    cbf_physiological_hz: float = 10.50  # Normal healthy CBF
    cbf_target_threshold_hz: float = 8.00  # Minimum therapeutic ciliary clearance frequency
    cbf_transition_asl_um: float = 5.50  # Midpoint of sigmoidal recovery curve
    cbf_transition_slope_um: float = 0.70  # Steepness of PCL uncoupling


@dataclass(frozen=True)
class TMEM16APharmacodynamics:
    """Allosteric calcium-sensitization model for TMEM16A potentiators."""
    ec50_calcium_basal_um: float = 1.20  # Resting Ca2+ EC50 without drug (uM)
    potentiator_kd_um: float = 0.10  # Potentiator binding affinity Kd (uM)
    max_ec50_shift_fraction: float = 0.85  # Max 85% leftward shift in calcium EC50
    calcium_hill_coefficient: float = 2.0  # Cooperative gating coefficient

    def calculate_apparent_calcium_ec50(self, drug_conc_um: float) -> float:
        """Compute apparent Ca2+ EC50 under potentiator concentration (uM)."""
        if drug_conc_um < 0.0:
            raise ValueError(f"Drug concentration cannot be negative: {drug_conc_um}")
        if drug_conc_um == 0.0:
            return self.ec50_calcium_basal_um
        fractional_shift = (self.max_ec50_shift_fraction * drug_conc_um) / (
            self.potentiator_kd_um + drug_conc_um
        )
        return self.ec50_calcium_basal_um * (1.0 - fractional_shift)

    def calculate_open_probability(
        self,
        intracellular_ca_um: float = 0.15,
        potentiator_conc_um: float = 0.0,
    ) -> float:
        """Calculate TMEM16A open probability Po in [0.0, 1.0]."""
        if intracellular_ca_um <= 0.0:
            return 0.0
        ec50_app = self.calculate_apparent_calcium_ec50(potentiator_conc_um)
        ca_pow = intracellular_ca_um ** self.calcium_hill_coefficient
        ec50_pow = ec50_app ** self.calcium_hill_coefficient
        po = ca_pow / (ca_pow + ec50_pow)
        return max(0.0, min(1.0, po))


@dataclass(frozen=True)
class ENaCPharmacodynamics:
    """ENaC sodium absorption and pharmacological blockade kinetics."""
    baseline_po_wt: float = 0.25  # Basal ENaC activity when CFTR is functional
    baseline_po_cf: float = 0.70  # Hyperabsorbed ENaC activity in CFTR-null epithelia
    blocker_ic50_um: float = 0.10  # Small molecule / peptide blocker IC50 (uM)
    hill_coefficient: float = 1.0  # Binding cooperativity

    def calculate_inhibition_fraction(self, blocker_conc_um: float) -> float:
        """Calculate fractional ENaC inhibition phi_inhib in [0.0, 1.0]."""
        if blocker_conc_um < 0.0:
            raise ValueError(f"Blocker concentration cannot be negative: {blocker_conc_um}")
        if blocker_conc_um == 0.0:
            return 0.0
        conc_pow = blocker_conc_um ** self.hill_coefficient
        ic50_pow = self.blocker_ic50_um ** self.hill_coefficient
        return conc_pow / (conc_pow + ic50_pow)

    def calculate_active_conductance(
        self,
        g_max: float,
        cftr_function_fraction: float,
        inhibition_fraction: float,
    ) -> float:
        """Compute net apical ENaC sodium conductance (mS/cm^2)."""
        po_basal = self.baseline_po_cf - (self.baseline_po_cf - self.baseline_po_wt) * min(
            1.0, max(0.0, cftr_function_fraction)
        )
        inhib = min(1.0, max(0.0, inhibition_fraction))
        return g_max * po_basal * (1.0 - inhib)


def compute_transepithelial_currents(
    constants: AirwayBiophysicalConstants,
    cftr_fraction: float,
    tmem16a_po: float,
    enac_inhibition_fraction: float,
) -> Dict[str, float]:
    """Calculate individual ion fluxes and net secretory short-circuit current (uA/cm^2)."""
    # 1. Apical Anion Secretion (Outward flow defined positive)
    g_cftr = constants.g_cftr_wt * max(0.0, min(1.0, cftr_fraction))
    g_tmem16a = constants.g_tmem16a_max * max(0.0, min(1.0, tmem16a_po))

    # Driving forces (V_apical - E_ion)
    df_cl = constants.v_apical_mv - constants.e_cl_mv  # typically -55 mV
    df_hco3 = constants.v_apical_mv - constants.e_hco3_mv  # typically -50 mV

    i_cl_cftr = -g_cftr * df_cl
    i_cl_tmem16a = -g_tmem16a * df_cl
    i_cl_total = i_cl_cftr + i_cl_tmem16a

    g_hco3_cftr = g_cftr * constants.p_hco3_over_cl_cftr
    g_hco3_tmem16a = g_tmem16a * constants.p_hco3_over_cl_tmem16a
    i_hco3_total = -(g_hco3_cftr + g_hco3_tmem16a) * df_hco3

    i_anion_total = i_cl_total + i_hco3_total

    # 2. Apical Sodium Influx (Absorptive inward current)
    enac_pd = ENaCPharmacodynamics()
    g_na = enac_pd.calculate_active_conductance(
        constants.g_enac_max, cftr_fraction, enac_inhibition_fraction
    )
    df_na = constants.v_apical_mv - constants.e_na_mv  # typically -100 mV
    i_na = g_na * (-df_na)  # absorptive magnitude

    # 3. Net Short-Circuit Secretory Current (uA/cm^2)
    i_sc_net = i_anion_total - i_na

    return {
        "g_cftr_mS_cm2": g_cftr,
        "g_tmem16a_mS_cm2": g_tmem16a,
        "g_na_mS_cm2": g_na,
        "i_cl_uA_cm2": i_cl_total,
        "i_hco3_uA_cm2": i_hco3_total,
        "i_anion_uA_cm2": i_anion_total,
        "i_na_uA_cm2": i_na,
        "i_sc_net_uA_cm2": i_sc_net,
    }


def compute_osmotic_water_flux(
    constants: AirwayBiophysicalConstants,
    i_sc_net_ua_cm2: float,
    current_asl_height_um: float,
) -> float:
    """Calculate rate of change of ASL height (dh_ASL/dt in um/min).

    Coupled water flux equation:
      dh_ASL/dt = (h_target(I_net) - h_ASL) / tau - J_evap
    """
    # Equilibrium target height from net current
    # Baseline CF height: 2.5 um; fluid gain: ~0.15 um per uA/cm^2 secretory current
    if i_sc_net_ua_cm2 <= 0.0:
        h_target = constants.baseline_cf_asl_height_um
    else:
        h_target = constants.baseline_cf_asl_height_um + i_sc_net_ua_cm2 * 0.15
        h_target = min(11.5, h_target)

    dh_dt = (h_target - current_asl_height_um) / constants.tau_hydration_min
    return dh_dt


def calculate_ciliary_beat_frequency(
    constants: AirwayBiophysicalConstants,
    asl_height_um: float,
) -> float:
    """Calculate ciliary beat frequency (CBF in Hz) as a function of ASL height."""
    if asl_height_um <= 0.0:
        return 0.0
    exponent = -(asl_height_um - constants.cbf_transition_asl_um) / constants.cbf_transition_slope_um
    clamped_exp = max(-30.0, min(30.0, exponent))
    logistic_factor = 1.0 / (1.0 + math.exp(clamped_exp))
    cbf = constants.cbf_quiescent_hz + (
        constants.cbf_physiological_hz - constants.cbf_quiescent_hz
    ) * logistic_factor
    return cbf


def calculate_mucociliary_transport_velocity(
    constants: AirwayBiophysicalConstants,
    cbf_hz: float,
    asl_height_um: float,
) -> float:
    """Calculate mucociliary clearance transport velocity (mm/min)."""
    if asl_height_um < constants.cilia_axoneme_length_um:
        penetration_ratio = max(0.0, (asl_height_um - 2.0) / (constants.cilia_axoneme_length_um - 2.0))
        coupling_efficiency = 0.15 * (penetration_ratio ** 2)
    else:
        excess_depth = min(4.0, asl_height_um - constants.cilia_axoneme_length_um)
        coupling_efficiency = 0.70 + 0.30 * (excess_depth / 4.0)

    base_velocity = 10.0 * (cbf_hz / constants.cbf_physiological_hz)
    return base_velocity * coupling_efficiency


def simulate_dynamic_rehydration(
    constants: AirwayBiophysicalConstants,
    cftr_fraction: float,
    tmem16a_potentiator_um: float,
    enac_blocker_um: float,
    intracellular_ca_um: float = 0.20,
    duration_min: float = 45.0,
    dt_min: float = 0.50,
) -> Dict[str, Any]:
    """Run dynamic time-course ODE integration of ASL height rehydration and CBF recovery."""
    tmem_pd = TMEM16APharmacodynamics()
    enac_pd = ENaCPharmacodynamics()

    po_tmem = tmem_pd.calculate_open_probability(
        intracellular_ca_um=intracellular_ca_um,
        potentiator_conc_um=tmem16a_potentiator_um,
    )
    enac_inhib = enac_pd.calculate_inhibition_fraction(enac_blocker_um)

    currents = compute_transepithelial_currents(
        constants=constants,
        cftr_fraction=cftr_fraction,
        tmem16a_po=po_tmem,
        enac_inhibition_fraction=enac_inhib,
    )

    steps = int(duration_min / dt_min)
    time_series_min: List[float] = [0.0]
    height_series_um: List[float] = [constants.baseline_cf_asl_height_um]
    cbf_series_hz: List[float] = [
        calculate_ciliary_beat_frequency(constants, constants.baseline_cf_asl_height_um)
    ]

    h_current = constants.baseline_cf_asl_height_um
    for step in range(steps):
        dh_dt = compute_osmotic_water_flux(constants, currents["i_sc_net_uA_cm2"], h_current)
        h_current += dh_dt * dt_min
        h_current = max(1.50, min(12.00, h_current))
        t_current = (step + 1) * dt_min
        cbf_current = calculate_ciliary_beat_frequency(constants, h_current)

        time_series_min.append(round(t_current, 2))
        height_series_um.append(round(h_current, 3))
        cbf_series_hz.append(round(cbf_current, 2))

    steady_state_h = height_series_um[-1]
    steady_state_cbf = cbf_series_hz[-1]
    steady_state_mct = calculate_mucociliary_transport_velocity(
        constants, steady_state_cbf, steady_state_h
    )

    target_asl_restored = steady_state_h >= constants.target_clearance_asl_height_um
    target_cbf_restored = steady_state_cbf >= constants.cbf_target_threshold_hz

    return {
        "cftr_fraction": cftr_fraction,
        "tmem16a_potentiator_conc_um": tmem16a_potentiator_um,
        "enac_blocker_conc_um": enac_blocker_um,
        "intracellular_ca_um": intracellular_ca_um,
        "tmem16a_po": round(po_tmem, 4),
        "enac_inhibition_fraction": round(enac_inhib, 4),
        "currents": {k: round(v, 3) for k, v in currents.items()},
        "steady_state_asl_height_um": steady_state_h,
        "steady_state_cbf_hz": steady_state_cbf,
        "steady_state_mct_mm_min": round(steady_state_mct, 2),
        "target_asl_restored": target_asl_restored,
        "target_cbf_restored": target_cbf_restored,
        "synergistic_clearance_restored": target_asl_restored and target_cbf_restored,
        "final_time_min": time_series_min[-1],
    }


def evaluate_genotype_theratyping_matrix() -> Dict[str, Any]:
    """Benchmark dual TMEM16A/ENaC synergy across Class I/VII modulator-ineligible CFTR genotypes."""
    constants = AirwayBiophysicalConstants()
    genotypes = {
        "W1282X": {"class": "Class I (Nonsense)", "cftr_fraction": 0.00, "description": "Exon 23 PTC; NBD2 deletion"},
        "G542X": {"class": "Class I (Nonsense)", "cftr_fraction": 0.00, "description": "Exon 12 PTC; early termination"},
        "CFTRdele2,3": {"class": "Class VII (Deletion)", "cftr_fraction": 0.00, "description": "Large deletion; complete transcript absence"},
        "F508del_Residual": {"class": "Class II (Residual)", "cftr_fraction": 0.03, "description": "Severe folding defect; minimal membrane leakage"},
    }

    regimens = {
        "Untreated": {"pot_um": 0.00, "block_um": 0.00},
        "TMEM16A_Monotherapy_Submaximal": {"pot_um": 1.50, "block_um": 0.00},
        "ENaC_Monotherapy_Submaximal": {"pot_um": 0.00, "block_um": 0.50},
        "Dual_Synergy_Regimen": {"pot_um": 1.50, "block_um": 0.50},
    }

    matrix: Dict[str, Dict[str, Any]] = {}
    for g_id, g_meta in genotypes.items():
        matrix[g_id] = {"meta": g_meta, "regimens": {}}
        for r_id, r_doses in regimens.items():
            sim = simulate_dynamic_rehydration(
                constants=constants,
                cftr_fraction=g_meta["cftr_fraction"],
                tmem16a_potentiator_um=r_doses["pot_um"],
                enac_blocker_um=r_doses["block_um"],
                intracellular_ca_um=0.20,
                duration_min=45.0,
            )
            matrix[g_id]["regimens"][r_id] = {
                "asl_height_um": sim["steady_state_asl_height_um"],
                "cbf_hz": sim["steady_state_cbf_hz"],
                "mct_mm_min": sim["steady_state_mct_mm_min"],
                "synergistic_clearance_restored": sim["synergistic_clearance_restored"],
            }

    return matrix


def run_theratyping_discovery_cli() -> int:
    """Command-line interface for the theratyping discovery engine."""
    parser = argparse.ArgumentParser(
        description="Dual TMEM16A/ENaC Alternative Conductance & ASL/CBF Synergy Simulator (T-0444)"
    )
    parser.add_argument("--potentiator", type=float, default=1.50, help="TMEM16A potentiator concentration (uM)")
    parser.add_argument("--blocker", type=float, default=0.50, help="ENaC blocker concentration (uM)")
    parser.add_argument("--cftr", type=float, default=0.00, help="Functional CFTR fraction [0.0, 1.0]")
    parser.add_argument("--ca", type=float, default=0.20, help="Resting intracellular calcium (uM)")
    parser.add_argument("--duration", type=float, default=45.0, help="Simulation duration (minutes)")
    parser.add_argument("--theratyping-matrix", action="store_true", help="Run full Class I/VII theratyping matrix")
    parser.add_argument("--json", action="store_true", help="Output raw JSON format")
    args = parser.parse_args()

    constants = AirwayBiophysicalConstants()

    if args.theratyping_matrix:
        matrix = evaluate_genotype_theratyping_matrix()
        if args.json:
            print(json.dumps(matrix, indent=2))
        else:
            print("=================================================================")
            print("TMEM16A / ENAC ALTERNATIVE CONDUCTANCE THERATYPING MATRIX")
            print("=================================================================")
            for gid, data in matrix.items():
                print(f"\nGenotype: {gid} ({data['meta']['class']})")
                for rname, rres in data["regimens"].items():
                    restored = "YES (REHYDRATED)" if rres["synergistic_clearance_restored"] else "NO (IMPAIRED)"
                    print(
                        f"  - {rname:32}: ASL={rres['asl_height_um']:5.2f} um | "
                        f"CBF={rres['cbf_hz']:5.2f} Hz | MCT={rres['mct_mm_min']:4.2f} mm/min | {restored}"
                    )
        return 0

    res = simulate_dynamic_rehydration(
        constants=constants,
        cftr_fraction=args.cftr,
        tmem16a_potentiator_um=args.potentiator,
        enac_blocker_um=args.blocker,
        intracellular_ca_um=args.ca,
        duration_min=args.duration,
    )

    if args.json:
        print(json.dumps(res, indent=2))
    else:
        print("=================================================================")
        print("TMEM16A POTENTIATION & ENAC BLOCKADE ASL / CBF DYNAMICS (T-0444)")
        print("=================================================================")
        print(f"CFTR Activity Fraction:           {res['cftr_fraction'] * 100:.1f}%")
        print(f"TMEM16A Potentiator Conc:         {res['tmem16a_potentiator_conc_um']:.3f} uM (Po={res['tmem16a_po']:.4f})")
        print(f"ENaC Blocker Conc:                {res['enac_blocker_conc_um']:.3f} uM (Inhib={res['enac_inhibition_fraction'] * 100:.1f}%)")
        print(f"Net Secretory Current:            {res['currents']['i_sc_net_uA_cm2']:.2f} uA/cm^2")
        print(f"Steady-State ASL Height:          {res['steady_state_asl_height_um']:.2f} um (Target >= 7.0 um: {res['target_asl_restored']})")
        print(f"Steady-State CBF:                 {res['steady_state_cbf_hz']:.2f} Hz (Target >= 8.0 Hz: {res['target_cbf_restored']})")
        print(f"Mucociliary Transport:            {res['steady_state_mct_mm_min']:.2f} mm/min")
        print(f"Synergistic Clearance Restored:   {res['synergistic_clearance_restored']}")
        print(f"\n{RESEARCH_BOUNDARY}")
        print(f"{NLM_ATTRIBUTION}")

    return 0


if __name__ == "__main__":
    sys.exit(run_theratyping_discovery_cli())
