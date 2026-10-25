#!/usr/bin/env python3
"""Universal Alternative Anion Conductance Resuscitation & Fluid Homeostasis Model (T-0430).

Models synergistic TMEM16A (ANO1) potentiation and ENaC inhibition in CFTR-null airway epithelia
to reverse osmotic water absorption and resuscitate Airway Surface Liquid (ASL) height >= 7.0 um.

Primary Evidence:
  - Am J Respir Crit Care Med (2020), DOI: 10.1164/rccm.201908-1641OC
  - J Clin Invest (2019), DOI: 10.1172/JCI128564
  - Eur Respir J (2020), DOI: 10.1183/13993003.congress-2020.4449
  - J Physiol (2024), DOI: 10.1113/JP286559
  - J Gen Physiol (2001), DOI: 10.1085/jgp.118.2.223
  - Science (2012), DOI: 10.1126/science.1223012
  - N Engl J Med (2019), DOI: 10.1056/NEJMra1813799
  - Physiol Rev (2022), DOI: 10.1152/physrev.00004.2021
Clinical Trials:
  - NCT04488705 (Inhaled ETD002 TMEM16A potentiator in CF)
  - NCT03647228 (IONIS-ENaCRx ASO in CF)
  - NCT04375514 (ARO-ENaC siRNA in CF)
  - NCT02343445 (ENaC inhibition in CF)
Structures:
  - RCSB PDB 6MSM, 5UAK, 8EIQ, 8EJ1
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Tuple

NLM_ATTRIBUTION = (
    "ClinicalTrials.gov does not evaluate safety or efficacy; "
    "listing on the registry does not denote government endorsement."
)

CLINICAL_TRIAL_ETD002 = "NCT04488705"
CLINICAL_TRIAL_IONIS_ENAC = "NCT03647228"
CLINICAL_TRIAL_ARO_ENAC = "NCT04375514"
CLINICAL_TRIAL_ENAC_CLEAR = "NCT02343445"

PDB_DEPHOSPHO_CFTR = "6MSM"
PDB_ACTIVE_CFTR = "5UAK"
PDB_MODULATOR_CFTR = "8EIQ"
PDB_GATING_CFTR = "8EJ1"

PRIMARY_DOI_TMEM16A_POT = "10.1164/rccm.201908-1641OC"
PRIMARY_DOI_TMEM16A_DIFF = "10.1172/JCI128564"
PRIMARY_DOI_TMEM16A_CONG = "10.1183/13993003.congress-2020.4449"
PRIMARY_DOI_ENAC_FURIN = "10.1113/JP286559"
PRIMARY_DOI_ASL_SURFACE = "10.1085/jgp.118.2.223"
PRIMARY_DOI_PERICILIARY = "10.1126/science.1223012"
PRIMARY_DOI_MUCO_OBSTRUCT = "10.1056/NEJMra1813799"
PRIMARY_DOI_AIRWAY_MUCUS = "10.1152/physrev.00004.2021"


@dataclass(frozen=True)
class TMEM16APotentiationModel:
    """Allosteric calcium and potentiator gating of TMEM16A (ANO1)."""
    ec50_basal_uM: float = 1.20    # Baseline calcium EC50 without drug (uM)
    ec50_pot_uM: float = 0.10      # Potentiator binding affinity EC50 (uM)
    max_shift_fraction: float = 0.85 # Maximum leftward shift in calcium EC50 (85%)
    hill_coefficient: float = 2.0   # Cooperative calcium binding Hill coefficient

    def apparent_calcium_ec50(self, potentiator_conc_uM: float) -> float:
        """Calculate apparent calcium EC50 under potentiator concentration."""
        if potentiator_conc_uM <= 0.0:
            return self.ec50_basal_uM
        pot_factor = (self.max_shift_fraction * potentiator_conc_uM) / (
            self.ec50_pot_uM + potentiator_conc_uM
        )
        return self.ec50_basal_uM * (1.0 - pot_factor)

    def open_probability(
        self,
        ca_conc_uM: float = 0.10,
        potentiator_conc_uM: float = 0.0,
    ) -> float:
        """Calculate TMEM16A open probability (Po) in [0.0, 1.0]."""
        if ca_conc_uM <= 0.0:
            return 0.0
        ec50_app = self.apparent_calcium_ec50(potentiator_conc_uM)
        ca_h = ca_conc_uM ** self.hill_coefficient
        ec50_h = ec50_app ** self.hill_coefficient
        return ca_h / (ca_h + ec50_h)


@dataclass
class AirwayFluidHomeostasisCircuit:
    """Coupled transepithelial equivalent circuit and ASL height dynamics."""
    # Conductances in mS/cm^2
    g_cftr_wt: float = 2.50
    g_tmem16a_max: float = 3.50
    g_enac_max: float = 4.00
    # Membrane potentials in mV
    v_apical_mV: float = -40.0
    e_cl_mV: float = 15.0
    e_na_mV: float = 60.0
    # ASL Homeostasis Parameters
    asl_target_um: float = 7.0     # Critical PCL clearance threshold
    evaporation_rate_um_min: float = 0.05 # Baseline evaporation rate

    def compute_short_circuit_current(
        self,
        cftr_activity_fraction: float = 0.0,
        tmem16a_po: float = 0.01,
        enac_inhibition_fraction: float = 0.0,
    ) -> Dict[str, float]:
        """Calculate apical ion currents and net short-circuit current (uA/cm^2)."""
        # Apical Chloride Conductance
        g_cl = (self.g_cftr_wt * cftr_activity_fraction) + (self.g_tmem16a_max * tmem16a_po)
        driving_force_cl = self.v_apical_mV - self.e_cl_mV # typically negative -> outward Cl- flow (secretory)
        i_cl = -g_cl * driving_force_cl # Secretory current defined positive

        # Apical Sodium Conductance (ENaC Po is elevated in CF when CFTR = 0)
        enac_baseline_po = 0.70 if cftr_activity_fraction < 0.05 else 0.25
        g_na = self.g_enac_max * enac_baseline_po * (1.0 - enac_inhibition_fraction)
        driving_force_na = self.v_apical_mV - self.e_na_mV # negative -> inward Na+ flow (absorptive)
        i_na = g_na * (-driving_force_na) # In absorptive direction

        # Net secretory short-circuit current (uA/cm^2)
        i_sc_net = i_cl - i_na

        return {
            "g_cl_mS_cm2": g_cl,
            "g_na_mS_cm2": g_na,
            "i_cl_uA_cm2": i_cl,
            "i_na_uA_cm2": i_na,
            "i_sc_net_uA_cm2": i_sc_net,
        }

    def steady_state_asl_height(self, i_sc_net_uA_cm2: float) -> float:
        """Estimate steady-state Airway Surface Liquid height (um) from net current."""
        # Water transport proportionality: 0.15 um height gain per uA/cm^2 net secretory current
        fluid_gain = i_sc_net_uA_cm2 * 0.15
        baseline_cf_height = 2.80 # Collapsed ASL height in CF
        h_asl = baseline_cf_height + fluid_gain
        # Clamp to realistic physical range [1.5, 12.0] um
        return max(1.5, min(12.0, h_asl))

    def evaluate_regimen(
        self,
        cftr_activity_fraction: float = 0.0,
        potentiator_conc_uM: float = 0.0,
        enac_inhibition_fraction: float = 0.0,
        ca_conc_uM: float = 0.10,
    ) -> Dict[str, float]:
        """Evaluate full alternative conductance regimen on ASL height restoration."""
        tmem_model = TMEM16APotentiationModel()
        po_tmem = tmem_model.open_probability(ca_conc_uM, potentiator_conc_uM)
        circuit_res = self.compute_short_circuit_current(
            cftr_activity_fraction=cftr_activity_fraction,
            tmem16a_po=po_tmem,
            enac_inhibition_fraction=enac_inhibition_fraction,
        )
        h_asl = self.steady_state_asl_height(circuit_res["i_sc_net_uA_cm2"])
        restored = h_asl >= self.asl_target_um

        return {
            "tmem16a_po": po_tmem,
            "i_sc_net_uA_cm2": circuit_res["i_sc_net_uA_cm2"],
            "asl_height_um": h_asl,
            "target_restored": 1.0 if restored else 0.0,
        }
